"""geosim の数値検証。

pytest は使わず、素の Python で走る（依存を増やさないため）。

    cd nea-simulation-pj
    .venv\\Scripts\\python tests\\test_geosim.py

各テストは「解析的に答えが分かる問題」を解かせて突き合わせる。
実装を書き換えたときに、数式レベルで壊れていないかを確かめるためのもの。

※ test.dxf の吸音率テーブル（absorption.csv）は Git 管理外なので、
  無い場合は既定の吸音率で走る（結果の判定には影響しない項目だけを見る）。
"""

import os
import sys
import itertools

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "geosim"))

import read_dxffile as rd            # noqa: E402
import mesh_method as mm             # noqa: E402
import sound_ray as sr               # noqa: E402
import loop_reflectionmesh as lr     # noqa: E402
import receiver_sphere as rs         # noqa: E402
import loop_deleteredundancy as ld   # noqa: E402
import loop_noredundancy as ln       # noqa: E402
import impulse as ip                 # noqa: E402
import reverberation as rv           # noqa: E402
import absorption as ab              # noqa: E402
import atmosphere as at              # noqa: E402

C = ln.SOUND_VELOCITY
TEST_DXF = os.path.join(ROOT, "test.dxf")
ABSORPTION = os.path.join(ROOT, "absorption.csv")
SAMPLE_ABSORPTION = os.path.join(ROOT, "data", "absorption_sample.csv")

_results = []


def check(name, ok, detail=""):
    _results.append((name, bool(ok)))
    mark = "OK " if ok else "NG "
    print(f"  [{mark}] {name}" + (f"  {detail}" if detail else ""))
    return ok


_room_cache = {}


def load_test_room(band_number=8):
    """test.dxf を読む（同じものを何度も読まないようにキャッシュする）。"""
    if band_number in _room_cache:
        return _room_cache[band_number]
    table = None
    if os.path.exists(ABSORPTION):
        library = ab.MaterialLibrary.from_csv(ABSORPTION, kind="normal")
        table = library.absorption_table(band_number=band_number, warn=False)
    model = rd.read_model(TEST_DXF, absorption_table=table,
                          band_number=band_number, verbose=False)
    _room_cache[band_number] = (model, np.asarray(model.source_points[0]),
                                np.asarray(model.receiver_points[0]))
    return _room_cache[band_number]


# 直方体 [0,2]x[0,3]x[0,1] の 6 面（法線は室内向き）
FACES = {"x=0": (np.array([1., 0, 0]), np.array([0., 0, 0])),
         "x=2": (np.array([-1., 0, 0]), np.array([2., 0, 0])),
         "y=0": (np.array([0., 1, 0]), np.array([0., 0, 0])),
         "y=3": (np.array([0., -1, 0]), np.array([0., 3, 0])),
         "z=0": (np.array([0., 0, 1]), np.array([0., 0, 0])),
         "z=1": (np.array([0., 0, -1]), np.array([0., 0, 1]))}


def mirror(point, face_name):
    normal, on_plane = FACES[face_name]
    return point - 2.0 * np.dot(normal, point - on_plane) * normal


# ---------------------------------------------------------------- ① DXF 読み込み
def test_read_dxf():
    print("\n[1] DXF 読み込み（test.dxf = 2x3x1 m の直方体）")
    model, src, rec = load_test_room()
    check("三角形 12 枚", len(model.mesh) == 12, f"{len(model.mesh)} 枚")
    check("mm → m の単位換算", np.allclose(model.extents[1] - model.extents[0],
                                            [2.0, 3.0, 1.0]),
          f"寸法 {(model.extents[1] - model.extents[0]).tolist()}")
    check("閉じた形状として認識", model.is_closed)
    check("体積 = 6.0 m^3", abs(abs(model.volume) - 6.0) < 1e-9,
          f"{model.volume:.6f}")
    check("法線が室内（空気側）を向いていると判定",
          model.shells and model.shells[0]["normals"] == "inward",
          model.shells[0]["normals"] if model.shells else "シェル情報なし")
    area = sum(0.5 * np.linalg.norm(np.cross(m.vertexes[1] - m.vertexes[0],
                                             m.vertexes[2] - m.vertexes[0]))
               for m in model.mesh)
    check("表面積 = 22.0 m^2", abs(area - 22.0) < 1e-9, f"{area:.6f}")


# ---------------------------------------------------------------- ② 音線生成
def test_soundray_generator():
    print("\n[2] 音線生成（Fibonacci 螺旋）")
    rays = sr.soundray_generator(5000)
    norms = np.linalg.norm(rays, axis=1)
    check("全て単位ベクトル", np.abs(norms - 1.0).max() < 1e-12,
          f"最大誤差 {np.abs(norms - 1.0).max():.2e}")
    check("重心がほぼ原点（偏りがない）", np.linalg.norm(rays.mean(axis=0)) < 1e-3,
          f"|重心| = {np.linalg.norm(rays.mean(axis=0)):.2e}")
    check("NaN が無い", not np.isnan(rays).any())


# ---------------------------------------------------------------- ③ エネルギー減衰
def test_energy_decay():
    print("\n[3] 斜入射のエネルギー反射率（書籍 式2.64）")
    normal = np.array([0.0, 0.0, 1.0])
    for alpha in (0.05, 0.3, 0.8):
        # 垂直入射では |R|^2 = 1 - α になるはず
        got = sr.energy_decay(np.array([0.0, 0.0, -1.0]), normal, alpha, 1.0)
        check(f"垂直入射 α={alpha} で 1-α", abs(got - (1.0 - alpha)) < 1e-12,
              f"{got:.12f} / 期待 {1 - alpha:.12f}")
    # かすめ入射（法線から 89.9°＝面すれすれ）では反射率 → 1（吸音が効かない）。
    # 局所反応性壁面の性質で、書籍 2.64 式が θ→90° で |R|→1 になることの確認
    theta = np.radians(89.9)
    grazing = np.array([np.sin(theta), 0.0, -np.cos(theta)])
    got = sr.energy_decay(grazing, normal, 0.8, 1.0)
    check("かすめ入射では反射率がほぼ 1（吸音が効かない）", got > 0.9,
          f"cosθ={abs(grazing[2]):.5f} で |R|^2 = {got:.6f}")


# ---------------------------------------------------------------- ④⑤ バックトレース
def test_backtrace():
    print("\n[4] 虚音源バックトレース（解析解と突き合わせ）")
    model, src, rec = load_test_room()
    mesh = model.mesh

    direct = float(np.linalg.norm(rec - src))
    first = {name: float(np.linalg.norm(rec - mirror(src, name))) for name in FACES}
    second = {}
    for a, b in itertools.permutations(FACES, 2):
        second[f"{a}->{b}"] = float(np.linalg.norm(rec - mirror(mirror(src, a), b)))

    rays = sr.soundray_generator(20000)
    history = ld.delete(lr.loop(src, rec, rays, 2, mesh, 0.15))
    pulses = ln.loop(src, rec, history, mesh, sound_velocity=C, verbose=False)

    # 直接音
    m0 = pulses.reflection_count == 0
    if check("直接音がちょうど 1 本", int(m0.sum()) == 1, f"{int(m0.sum())} 本"):
        check("直接音の距離が解析解と一致",
              abs(float(pulses.distance[m0][0]) - direct) < 1e-9,
              f"{float(pulses.distance[m0][0]):.9f} / {direct:.9f}")
        check("直接音のエネルギーが 1.0（無反射）",
              np.allclose(pulses.energy[m0][0], 1.0))
        check("直接音の到来方向が受音点→音源",
              np.allclose(pulses.direction[m0][0], (src - rec) / direct, atol=1e-12))

    # 1 次反射（凸な直方体なので 6 面ぶんちょうど 6 本）
    got1 = np.sort(pulses.distance[pulses.reflection_count == 1])
    want1 = np.sort(np.array(list(first.values())))
    if check("1 次反射がちょうど 6 本", len(got1) == 6, f"{len(got1)} 本"):
        check("1 次反射の距離が全て解析解と一致",
              np.abs(got1 - want1).max() < 1e-9,
              f"最大誤差 {np.abs(got1 - want1).max():.2e} m")

    # 2 次反射（本数は音線数次第だが、出たものは全て解析解のどれかに一致するはず）
    got2 = pulses.distance[pulses.reflection_count == 2]
    unmatched = [d for d in got2
                 if not any(abs(v - d) < 1e-9 for v in second.values())]
    check(f"2 次反射 {len(got2)} 本すべてが虚音源の解析距離と一致",
          len(unmatched) == 0, f"外れ {len(unmatched)} 本")

    # エネルギーの独立検算（床 z=0 の 1 次反射）
    image = mirror(src, "z=0")
    dist = float(np.linalg.norm(rec - image))
    vray = (image - rec) / dist
    node = rec + (-rec[2] / vray[2]) * vray
    floor_id = next(j for j, m in enumerate(mesh)
                    if np.allclose(m.normal, [0, 0, 1])
                    and abs(m.vertexes[0][2]) < 1e-9
                    and mm.collision_detection(node, m.vertexes))
    alpha = np.atleast_1d(mesh[floor_id].absorption_coefficient)
    cos_theta = abs(float(np.dot(vray, [0.0, 0.0, 1.0])))
    r0 = np.sqrt(1.0 - alpha)
    hand = np.abs(((1 + r0) * cos_theta - (1 - r0))
                  / ((1 + r0) * cos_theta + (1 - r0))) ** 2
    sel = np.where((pulses.reflection_count == 1)
                   & (np.abs(pulses.distance - dist) < 1e-9))[0]
    if len(sel):
        check("1 次反射のエネルギーが手計算の |R(θ)|^2 と一致",
              np.abs(pulses.energy[sel[0]] - hand).max() < 1e-12,
              f"最大誤差 {np.abs(pulses.energy[sel[0]] - hand).max():.2e}")


def test_image_sources():
    print("\n[5] 鏡像の生成（符号付き。元コードの abs() 版との違い）")
    model, src, _ = load_test_room()
    mesh = model.mesh
    # 面 0 に対して 2 回鏡像を取ると元に戻る（対合性）
    images = ln.image_sources(src, [0, 0], mesh)
    check("同じ面で 2 回鏡像を取ると元に戻る",
          np.allclose(images[2], src, atol=1e-12),
          f"誤差 {np.abs(images[2] - src).max():.2e}")
    # 面からの距離が鏡像前後で等しい
    normal = np.asarray(mesh[0].normal, dtype=float)
    d = mm.parameter_d(normal, mesh[0].vertexes[0])
    before = np.dot(normal, images[0]) + d
    after = np.dot(normal, images[1]) + d
    check("鏡像前後で面までの符号付き距離が符号反転",
          abs(before + after) < 1e-12, f"{before:.9f} / {after:.9f}")


# ---------------------------------------------------------------- ベクトル化
def test_vectorised_geometry():
    """F-1 の高速化。scalar 版と**ビット単位で一致**することを確かめる。

    速いだけでは意味がなく、結果が変わっていないことが前提。
    """
    print("\n[5.5] ベクトル化した交差判定（scalar 版との一致）")
    model, src, rec = load_test_room()
    mesh = model.mesh
    faces = mm.FaceArrays(mesh)

    rng = np.random.default_rng(1)
    origins = rng.uniform([0.1, 0.1, 0.1], [1.9, 2.9, 0.9], size=(500, 3))
    directions = rng.normal(size=(500, 3))
    directions /= np.linalg.norm(directions, axis=1)[:, None]

    hit_id, distance, node = faces.nearest_hit(origins, directions)

    scalar_id, scalar_distance = [], []
    for i in range(len(origins)):
        best_id, best = -1, np.inf
        for j in range(len(mesh)):
            collision, d = mm.collision_distance(directions[i], origins[i],
                                                 mesh[j].normal, mesh[j].vertexes)
            if collision and d < best:
                best, best_id = d, j
        scalar_id.append(best_id)
        scalar_distance.append(best)
    scalar_id = np.array(scalar_id)
    scalar_distance = np.array(scalar_distance)

    check("最寄り面の ID が scalar 版と完全一致",
          np.array_equal(hit_id, scalar_id),
          f"不一致 {int((hit_id != scalar_id).sum())} 件")
    finite = np.isfinite(scalar_distance)
    check("距離が scalar 版とビット単位で一致",
          np.array_equal(distance[finite], scalar_distance[finite]),
          f"最大差 {np.abs(distance[finite] - scalar_distance[finite]).max():.3e}")

    batch = rs.inside_sphere_batch(0.2, directions, origins, rec, distance)
    one = np.array([rs.inside_sphere(0.2, directions[i], origins[i], rec, distance[i])
                    for i in range(len(origins))])
    check("受音判定が scalar 版と完全一致", np.array_equal(batch, one),
          f"受音 {int(one.sum())} 本")

    # 分割して処理しても結果が変わらないこと
    small = mm.FaceArrays(mesh).nearest_hit(origins, directions, chunk_elements=37)
    check("音線を分割して処理しても結果が同じ",
          np.array_equal(small[0], hit_id) and np.array_equal(small[1], distance))


# ---------------------------------------------------------------- 大気
def test_atmosphere():
    print("\n[6] 大気（音速・空気吸収）")
    # 標準的な値との突き合わせ
    for t, rh, want, tol in [(0.0, 0.0, 331.5, 0.3), (20.0, 0.0, 343.2, 0.3),
                             (20.0, 40.0, 343.8, 0.5)]:
        got = at.sound_velocity(t, rh)
        check(f"音速 {t:.0f}℃/{rh:.0f}% が {want} m/s 付近", abs(got - want) < tol,
              f"{got:.2f} m/s")
    check("湿度が上がると音速が上がる",
          at.sound_velocity(20, 80) > at.sound_velocity(20, 0))
    check("飽和水蒸気圧 20℃ が 2.34 kPa 付近",
          abs(at.saturation_vapour_pressure(20.0) - 2.339) < 0.01,
          f"{at.saturation_vapour_pressure(20.0):.4f} kPa")
    check("元コードの 340 m/s はおよそ 14℃ 相当",
          abs(at.sound_velocity(14.0, 40.0) - 340.0) < 0.3,
          f"14℃/40% で {at.sound_velocity(14.0, 40.0):.2f} m/s")

    # 空気吸収: ISO 9613-1 の代表値（20℃ 70%RH）と桁が合うか
    got = at.absorption_db_per_metre(np.array([1000.0, 4000.0]), 20.0, 70.0) * 1000
    check("空気吸収 20℃70% 1kHz が 5 dB/km 付近", abs(got[0] - 4.8) < 1.5,
          f"{got[0]:.2f} dB/km")
    check("空気吸収 20℃70% 4kHz が 25 dB/km 付近", abs(got[1] - 25.0) < 8.0,
          f"{got[1]:.2f} dB/km")
    check("高音ほど強く減衰する",
          at.absorption_coefficient(8000.0) > at.absorption_coefficient(125.0))


# ---------------------------------------------------------------- 吸音率
def test_absorption():
    print("\n[7] 吸音率（Paris の式・材料ライブラリ）")
    z = np.array([0.5, 1.0, 2.0, 5.0, 20.0, 100.0])
    numeric = ab.statistical_absorption(z)
    closed = ab.statistical_absorption_closed_form(z)
    check("Paris の式：数値積分と解析解が一致",
          np.abs(numeric - closed).max() < 1e-8,
          f"最大差 {np.abs(numeric - closed).max():.2e}")

    check("α_s の最大値が 0.951（z≈1.57）",
          abs(ab.STATISTICAL_MAX - 0.951) < 0.002
          and abs(ab.STATISTICAL_MAX_IMPEDANCE - 1.567) < 0.01,
          f"{ab.STATISTICAL_MAX:.4f} (z={ab.STATISTICAL_MAX_IMPEDANCE:.3f})")

    alpha_n = np.array([0.02, 0.1, 0.3, 0.5, 0.8, 0.95])
    back = np.array([ab.random_to_normal(ab.normal_to_random(a), warn=False)
                     for a in alpha_n])
    check("垂直入射 → 残響室法 → 垂直入射 で元に戻る",
          np.abs(back - alpha_n).max() < 1e-9,
          f"最大誤差 {np.abs(back - alpha_n).max():.2e}")

    check("残響室法の値のほうが大きく出る（同じ材料なら α_s > α_n）",
          all(ab.normal_to_random(a) > a for a in alpha_n[:-1]))

    check("垂直入射吸音率と z の関係が整合",
          abs(ab.normal_absorption(ab.impedance_from_normal(0.3)) - 0.3) < 1e-12)

    # 1 を超える残響室法吸音率を丸める。
    # STATISTICAL_MAX_IMPEDANCE は格子探索の近似値、二分法の解はより正確なので
    # 完全一致はしない。0.951 付近に収まっていればよい
    got = ab.random_to_normal(1.15, warn=False)
    check("1 を超える残響室法吸音率を上限に丸める", abs(got - 0.9513) < 1e-3,
          f"α_n = {got:.4f}（上限 0.951 相当）")
    check("上限を超える値は全て同じ結果になる",
          abs(ab.random_to_normal(1.5, warn=False) - got) < 1e-9)

    # バンド定義
    check("8 バンドが 63〜8000 Hz",
          np.allclose(ab.octave_bands(8), [63, 125, 250, 500, 1000, 2000, 4000, 8000]))
    check("6 バンドが 125〜4000 Hz",
          np.allclose(ab.octave_bands(6), [125, 250, 500, 1000, 2000, 4000]))

    # 材料ライブラリ
    if os.path.exists(SAMPLE_ABSORPTION):
        lib = ab.MaterialLibrary.from_csv(SAMPLE_ABSORPTION, kind="normal")
        check("サンプル CSV を読める", len(lib) >= 7, f"{len(lib)} 種")
        lib.add("テスト材料", [0.2, 0.4, 0.7, 0.9, 0.9, 0.85], kind="random")
        check("材料を追加できる（GUI からの追加を想定）", "テスト材料" in lib)
        table = lib.absorption_table(band_number=8, warn=False)
        check("残響室法の材料が垂直入射に変換されて出てくる",
              table["テスト材料"].max() < 0.9,
              f"最大 {table['テスト材料'].max():.4f}（元は 0.9）")
        check("6 バンドの材料が 8 バンドに伸びる",
              len(table["カーペット"]) == 8)
        assignment = ab.LayerAssignment({"1": "カーペット"})
        table = lib.absorption_table(assignment, band_number=8, warn=False)
        check("LayerAssignment でレイヤに材料を割り当てられる",
              np.allclose(table["1"], table["カーペット"]))
        check("未割り当てのレイヤを検出できる",
              assignment.missing(["1", "存在しない"], lib) == ["存在しない"])


# ---------------------------------------------------------------- インパルス応答
def test_impulse():
    print("\n[8] インパルス応答の合成")
    # scipy 版が元コードの式と一致するか
    worst = 0.0
    for numtaps, lo, hi in [(1024, 0.05, 0.1), (4096, 0.001, 0.002)]:
        scipy_version = ip.filter_bandpass(numtaps, lo, hi)
        fortran_version = ip._fir1_bandpass_fortran(numtaps, lo, hi)
        worst = max(worst, np.abs(scipy_version - fortran_version).max()
                    / np.abs(fortran_version).max())
    check("scipy.signal.firwin が元コードの fir1_bandpass と一致", worst < 1e-12,
          f"最大相対誤差 {worst:.2e}")

    # 32 バンドの総和が平坦か
    nn = 8192
    fmax = ip.SAMPLING_FREQUENCY / 2
    mf = ip.third_octave_bands()
    lower, upper = ip.bandpass_edges(mf, fmax)
    total = sum(np.fft.fft(ip.filter_bandpass(nn, lower[j], upper[j]))
                for j in range(32))
    freq = np.fft.fftfreq(nn, 1 / ip.SAMPLING_FREQUENCY)
    mag = np.abs(total[(freq > 40) & (freq < 15000)])
    check("32 バンドの総和が平坦（40Hz〜15kHz）", 0.98 < mag.mean() < 1.02,
          f"平均振幅 {mag.mean():.4f}")

    # バンド割り当てが元コードの手書き表を再現するか
    fortran_table = np.array([0] * 11 + [1] * 3 + [2] * 3 + [3] * 3
                             + [4] * 3 + [5] * 9)
    got = ip.band_mapping(ab.octave_bands(6), ip.third_octave_bands())
    check("6 バンドの割り当てが元コードの手書き表（94〜125行）と完全一致",
          np.array_equal(got, fortran_table),
          f"差 {int((got != fortran_table).sum())} 個")

    got8 = ip.band_mapping(ab.octave_bands(8), ip.third_octave_bands())
    check("8 バンドでも 63Hz/8kHz が正しい 1/3 オクターブに対応",
          got8[6] == 0 and got8[27] == 7,
          f"63Hz→{got8[6]} / 8kHz→{got8[27]}")

    # 単一パルス（8 バンド）
    t0 = 0.01
    t, ir = ip.impulse_response(np.array([t0]), np.ones((1, 8)), verbose=False)
    dt = t[1] - t[0]
    peak = int(np.argmax(np.abs(ir)))
    check("単一パルスのピークが入力時刻に来る（遅れが無い）",
          abs(t[peak] - t0) < 2 * dt, f"{t[peak] * 1000:.4f} ms / 期待 {t0 * 1000:.4f} ms")
    check("出力長が max_time どおり", len(ir) == int(round(1.0 * 44100)), f"{len(ir)} 点")

    # 6 バンドでも動く
    t6, ir6 = ip.impulse_response(np.array([t0]), np.ones((1, 6)), verbose=False)
    check("6 バンドでも合成できる", len(ir6) == len(ir))

    # 距離減衰 1/r
    air = at.Atmosphere()
    e32 = ip.expand_to_third_octave(np.ones((1, 8)))
    hfp = ip.transfer_function(e32, np.array([t0]), 5, 1.0, air.sound_velocity)
    check("伝達関数の直流成分が 1/(t*c) と一致",
          abs(abs(hfp[0, 0]) - 1.0 / (t0 * air.sound_velocity)) < 1e-12,
          f"{abs(hfp[0, 0]):.9f} / {1.0 / (t0 * air.sound_velocity):.9f}")

    # 空気吸収が効いているか
    without = ip.expand_to_third_octave(np.ones((1, 8)))
    with_air = ip.apply_air_absorption(without, np.array([0.5]), air)
    check("空気吸収で高音ほど減る", with_air[-1, 0] < with_air[0, 0] < 1.0,
          f"最低バンド {with_air[0, 0]:.4f} / 最高バンド {with_air[-1, 0]:.6f}")


# ---------------------------------------------------------------- ⑦ 残響時間
def test_reverberation():
    print("\n[9] 残響時間（減衰率が既知の合成応答で検証）")
    fs = 44100.0
    dt = 1.0 / fs
    nn = 2 ** 17
    t = dt * np.arange(nn)
    rng = np.random.default_rng(0)

    for want_t60 in (0.3, 0.5, 1.0):
        tau = want_t60 / 6.907755          # 60 / (20 log10 e)
        ir = rng.standard_normal(nn) * np.exp(-t / tau)
        rt = rv.decay_curves(t, ir, verbose=False)["reverberation_time"]
        err = np.abs(rt - want_t60) / want_t60
        # 帯域幅×減衰時間が小さいと推定がばらつく（ISO 3382 でも既知）ので中央値で見る
        check(f"理論 T60 = {want_t60} s を再現（誤差の中央値）",
              np.nanmedian(err) < 0.03,
              f"中央値 {np.nanmedian(err) * 100:.2f} % / 最大 {np.nanmax(err) * 100:.2f} %")

    # 減衰曲線の傾き
    tau = 0.5 / 6.907755
    ir = rng.standard_normal(nn) * np.exp(-t / tau)
    result = rv.decay_curves(t, ir, verbose=False)
    # 全バンドの傾きから T60 を推定して平均で見る
    # （1 バンドだけだと雑音のばらつきで数 % ぶれる。BT 積の制約と同じ話）
    tt = result["time"]
    estimates = []
    for d in result["decay"]:
        mask = (d <= -5) & (d >= -35)
        if mask.sum() > 10:
            estimates.append(-60.0 / np.polyfit(tt[mask], d[mask], 1)[0])
    estimates = np.array(estimates)
    check("減衰曲線の傾きが理論と一致（全バンド平均）",
          abs(estimates.mean() - 0.5) / 0.5 < 0.03,
          f"平均 T60 = {estimates.mean():.4f} s / 理論 0.5000 s "
          f"（バンド別 {estimates.min():.3f}〜{estimates.max():.3f}）")
    d = result["decay"][3]
    check("フィルタの遅れが残っていない",
          tt[np.nonzero(d <= -1.0)[0][0]] < 0.05,
          f"-1 dB を切る時刻 {tt[np.nonzero(d <= -1.0)[0][0]] * 1000:.2f} ms")

    # FIR 方式でも同等に出るか
    rt_fir = rv.decay_curves(t, ir, method="fir", verbose=False)["reverberation_time"]
    err = np.abs(rt_fir - 0.5) / 0.5
    check("method='fir' でも同等の精度", np.nanmedian(err) < 0.05,
          f"誤差の中央値 {np.nanmedian(err) * 100:.2f} %")

    # EDT / T20 / T30 をまとめて出す
    measures = rv.decay_measures(t, ir, verbose=False)["measures"]
    check("EDT / T20 / T30 が揃って出る",
          set(measures) == {"EDT", "T20", "T30"}, str(sorted(measures)))
    for name in ("EDT", "T20", "T30"):
        err = np.abs(measures[name] - 0.5) / 0.5
        # 指数減衰なら EDT も T20 も T30 も同じ値になるはず
        check(f"{name} が理論 T60 = 0.5 s を再現", np.nanmedian(err) < 0.06,
              f"中央値 {np.nanmedian(err) * 100:.2f} % "
              f"（{np.nanmin(measures[name]):.3f}〜{np.nanmax(measures[name]):.3f} s）")
    check("直線減衰なら曲率がほぼ 0",
          np.nanmedian(np.abs(rv.decay_measures(t, ir, verbose=False)["curvature"])) < 10.0)


# ---------------------------------------------------------------- 統計残響式
def test_statistical_reverberation():
    print("\n[10] 統計残響式（Sabine / Eyring / Eyring-Knudsen）")
    model, _, _ = load_test_room()
    mesh = model.mesh

    areas = rv.triangle_areas(mesh)
    check("三角形の面積の合計が 22.0 m^2", abs(areas.sum() - 22.0) < 1e-9,
          f"{areas.sum():.9f}")

    # 全面を同じ吸音率にして解析解と突き合わせる
    alpha = 0.2
    import copy
    uniform = []
    for face in mesh:
        clone = copy.copy(face)
        clone.material = "uniform"
        clone.absorption_coefficient = np.full(8, alpha)
        uniform.append(clone)

    air = at.Atmosphere()
    result = rv.statistical_reverberation(uniform, 6.0, atmosphere=air,
                                          convert_to_random=False,
                                          include_air_absorption=False,
                                          verbose=False)
    constant = 24.0 * np.log(10.0) / air.sound_velocity
    want_sabine = constant * 6.0 / (22.0 * alpha)
    want_eyring = constant * 6.0 / (-22.0 * np.log(1.0 - alpha))
    check("Sabine が解析解と一致",
          np.abs(result["sabine"] - want_sabine).max() < 1e-12,
          f"{result['sabine'][0]:.9f} / 期待 {want_sabine:.9f}")
    check("Eyring が解析解と一致",
          np.abs(result["eyring"] - want_eyring).max() < 1e-12,
          f"{result['eyring'][0]:.9f} / 期待 {want_eyring:.9f}")
    check("空気吸収を切れば Eyring-Knudsen は Eyring と一致",
          np.abs(result["eyring_knudsen"] - result["eyring"]).max() < 1e-12)
    check("Eyring は Sabine より短く出る（吸音率が大きいほど差が開く）",
          np.all(result["eyring"] < result["sabine"]))
    check("平均自由行程が 4V/S", abs(result["mean_free_path"] - 4 * 6.0 / 22.0) < 1e-12,
          f"{result['mean_free_path']:.6f} m")

    # 乱入射への変換が効いているか（変換すると吸音率が上がるので残響が短くなる）
    converted = rv.statistical_reverberation(uniform, 6.0, atmosphere=air,
                                             convert_to_random=True,
                                             include_air_absorption=False,
                                             verbose=False)
    check("垂直入射→乱入射の変換で平均吸音率が上がる",
          np.all(converted["mean_absorption"] > result["mean_absorption"]),
          f"{result['mean_absorption'][0]:.4f} → {converted['mean_absorption'][0]:.4f}")
    check("そのぶん残響時間は短くなる",
          np.all(converted["sabine"] < result["sabine"]))

    # 空気吸収の項
    with_air = rv.statistical_reverberation(uniform, 6.0, atmosphere=air,
                                            convert_to_random=False,
                                            include_air_absorption=True,
                                            verbose=False)
    # ★Sabine と Eyring は**空気吸収を入れない素の形**。足したのが Eyring-Knudsen で、
    #   差がそのまま空気吸収の効きになる（ユーザー判断 2026-08-17）
    check("Sabine と Eyring は空気吸収の設定で変わらない",
          np.allclose(with_air["sabine"], result["sabine"])
          and np.allclose(with_air["eyring"], result["eyring"]))
    check("Eyring-Knudsen は Eyring より短い（空気吸収のぶん）",
          np.all(with_air["eyring_knudsen"] < with_air["eyring"]))
    check("空気吸収の効きは高音ほど大きい",
          (with_air["eyring"][-1] - with_air["eyring_knudsen"][-1])
          > (with_air["eyring"][0] - with_air["eyring_knudsen"][0]),
          f"63Hz {with_air['eyring'][0] - with_air['eyring_knudsen'][0]:.4f} s / "
          f"8kHz {with_air['eyring'][-1] - with_air['eyring_knudsen'][-1]:.4f} s")
    check("ミリントンは落とした", "millington" not in with_air)

    # 実際のモデルで面ごとに吸音率が違う場合
    real = rv.statistical_reverberation_from_model(model, atmosphere=air, verbose=False)
    check("DxfModel から直接計算できる", real is not None)
    check("Eyring-Knudsen < Eyring < Sabine の順になる",
          np.all(real["eyring_knudsen"] <= real["eyring"])
          and np.all(real["eyring"] < real["sabine"]))

    # 開いた形状では計算できない
    class OpenModel:
        is_closed = False
        open_edges = 11
    check("開いた形状では None を返す",
          rv.statistical_reverberation_from_model(OpenModel()) is None)


# ---------------------------------------------------------------- 可視化データ
def test_ray_log():
    """音線軌跡の記録と、可視化用の一括計算（G-1 / G-2 の土台）。"""
    print("\n[11] 音線軌跡（可視化用データ）")
    import tempfile
    import ray_recorder as rr
    import view_rays as vr

    model, src, rec = load_test_room()
    rays = sr.soundray_generator(400)
    recorder = rr.RayRecorder(total_rays=400, max_rays=100, sound_velocity=C,
                              band_number=8)
    lr.loop(src, rec, rays, 5, model.mesh, 0.2, recorder=recorder)

    check("軌跡が間引かれて記録される",
          len(recorder.trajectories) == 100 and recorder.stride == 4,
          f"{len(recorder.trajectories)} 本 / stride {recorder.stride}")
    check("音線番号が昇順に並ぶ（束処理でも順序が保たれる）",
          all(a.ray_index < b.ray_index for a, b in
              zip(recorder.trajectories, recorder.trajectories[1:])))
    check("先頭の節点が音源位置",
          np.allclose(recorder.trajectories[0].nodes[0], src))
    check("累積距離が単調増加",
          all(np.all(np.diff(t.distances) >= 0) for t in recorder.trajectories))

    with tempfile.TemporaryDirectory() as folder:
        path = os.path.join(folder, "log.npz")
        recorder.save_npz(path)
        log = vr.RayLog(path)

        check("保存して読み直せる", log.ray_count == 100, f"{log.ray_count} 本")
        check("受音の有無が一致",
              int(log.received.sum())
              == sum(1 for t in recorder.trajectories if len(t.receive_steps)),
              f"受音 {int(log.received.sum())} 本")

        # 一括計算した粒子位置が RayTrajectory.position_at と一致するか
        worst = 0.0
        for t in (0.002, 0.01, 0.03):
            position, rows = log.positions_at(t)
            for k, ray in enumerate(rows):
                want = recorder.trajectories[ray].position_at(t, C)
                worst = max(worst, float(np.abs(position[k] - want).max()))
        check("一括計算した粒子位置が 1 本ずつの計算と一致", worst < 1e-12,
              f"最大差 {worst:.2e}")

        # 折れ線の点数が節点数の合計と合う
        index = log.selection(max_rays=20)
        poly = log.line_polydata(index, colour="time")
        check("折れ線の点数が節点数の合計と一致",
              poly.n_points == int(log.node_counts[index].sum()),
              f"{poly.n_points} 点")

        truncated = log.line_polydata(index, colour="time", max_reflection=2)
        check("反射回数で打ち切ると点数が減る",
              truncated.n_points < poly.n_points,
              f"{truncated.n_points} 点（打ち切りなし {poly.n_points} 点）")

        check("受音した経路だけを選べる",
              np.all(log.received[log.selection(received_only=True)]))


# ---------------------------------------------------------------- 法線の向き
def test_normals():
    print("\n[13] 法線の向き（自動判定・手動反転）")
    triangles = [tuple(f.vertexes) for f in load_test_room()[0].mesh]

    # 閉じた直方体の内側は全方向で面に当たる。外側は当たらない
    inside = rd.encloses_point(triangles, [1.0, 1.5, 0.5])
    outside = rd.encloses_point(triangles, [10.0, 10.0, 10.0])
    check("閉じた室の内側は囲まれ度 1.0", abs(inside - 1.0) < 1e-9, f"{inside:.3f}")
    # 外の点は「その方向に箱が見える」ぶんだけ当たる（張る立体角ぶん）。
    # 0 にはならないが、閉じているかの判定に使う 0.95 には遠く及ばない
    check("室の外側は囲まれ度がしきい値をはるかに下回る",
          outside < 0.1, f"{outside:.3f}（しきい値 {rd.ENCLOSURE_THRESHOLD}）")

    # 'auto' は閉じた室では 'shells' を選ぶ
    model = rd.read_model(TEST_DXF, orient_normals="auto", verbose=False)
    check("閉じた室では 'auto' が 'shells' を選ぶ", model.orient_mode == "shells",
          model.orient_mode)
    centre = 0.5 * (model.extents[0] + model.extents[1])
    inward = sum(1 for f in model.mesh
                 if np.dot(f.normal, centre - f.vertexes.mean(axis=0)) > 0)
    check("'auto' の結果が全面内向き", inward == len(model.mesh),
          f"{inward}/{len(model.mesh)}")

    # 開いた形状（test2.dxf）は CAD のまま触らない
    test2 = os.path.join(ROOT, "test2.dxf")
    if os.path.exists(test2):
        opened = rd.read_model(test2, orient_normals="auto", verbose=False)
        plain = rd.read_model(test2, orient_normals="cad", verbose=False)
        same = all(np.allclose(a.normal, b.normal)
                   for a, b in zip(opened.mesh, plain.mesh))
        check("開いた形状では 'auto' が法線を触らない",
              opened.orient_mode == "cad" and same,
              f"mode={opened.orient_mode} 囲まれ度={opened.enclosure:.2f}")

    # レイの偶奇による内向き判定。わざと全反転した状態から復元できるか
    flipped_all = rd.read_model(TEST_DXF, orient_normals="flip", verbose=False)
    tri = [tuple(f.vertexes) for f in flipped_all.mesh]
    normals = np.array([f.normal for f in flipped_all.mesh])
    to_flip, ambiguous = rd.orient_inward(tri, normals)
    check("全反転した法線を偶奇判定が全部見つける",
          len(to_flip) == len(flipped_all.mesh) and ambiguous == 0,
          f"{len(to_flip)}/{len(flipped_all.mesh)} 判定割れ {ambiguous}")

    # flip_faces は「CAD の巻き順から反転する面」の絶対集合
    manual = rd.read_model(TEST_DXF, orient_normals="cad", flip_faces=[0, 3],
                           verbose=False)
    base = rd.read_model(TEST_DXF, orient_normals="cad", verbose=False)
    ok = all(np.allclose(m.normal, -b.normal if j in (0, 3) else b.normal)
             for j, (m, b) in enumerate(zip(manual.mesh, base.mesh)))
    check("flip_faces で指定した面だけが反転する", ok,
          f"反転 {sorted(manual.flipped_faces)}")

    # ★自動判定 → 保存 → 読み直しで**同じ法線に戻る**こと。
    #   ここが差分扱いだと、自動が反転した面を手動指定が二度反転して元に戻ってしまう
    #   （実際にそれで残響時間が半分近く変わった。2026-08-15）
    for dxf in (TEST_DXF, os.path.join(ROOT, "test2.dxf")):
        if not os.path.exists(dxf):
            continue
        auto = rd.read_model(dxf, orient_normals="auto", verbose=False)
        saved = sorted(auto.flipped_faces)
        again = rd.read_model(dxf, orient_normals="auto", flip_faces=saved,
                              verbose=False)
        same = all(np.allclose(a.normal, b.normal)
                   for a, b in zip(auto.mesh, again.mesh))
        check(f"自動判定を保存して読み直すと同じ法線になる（{os.path.basename(dxf)}）",
              same and sorted(again.flipped_faces) == saved,
              f"反転 {len(saved)}/{len(auto.mesh)} 枚")

    # **実際に反転が起きるケース**で往復を見る（test.dxf は自動反転が 0 枚なので
    # そのままだと二重反転の不具合を踏まない）。'flip' は全 12 枚を反転する
    flipped_model = rd.read_model(TEST_DXF, orient_normals="flip", verbose=False)
    saved = sorted(flipped_model.flipped_faces)
    restored = rd.read_model(TEST_DXF, orient_normals="flip", flip_faces=saved,
                             verbose=False)
    check("反転が起きるモデルでも保存→読み直しで同じ法線になる",
          len(saved) == len(flipped_model.mesh)
          and all(np.allclose(a.normal, b.normal)
                  for a, b in zip(flipped_model.mesh, restored.mesh)),
          f"反転 {len(saved)}/{len(flipped_model.mesh)} 枚")

    # 自動判定を人が上書きした場合も、そのとおりに再現される
    edited = sorted(set(sorted(rd.read_model(TEST_DXF, orient_normals="auto",
                                             verbose=False).flipped_faces)) ^ {2, 5})
    replayed = rd.read_model(TEST_DXF, orient_normals="auto", flip_faces=edited,
                             verbose=False)
    check("人が直した指定が自動判定より優先される",
          sorted(replayed.flipped_faces) == edited, f"反転 {edited}")

    # OCS→WCS。押し出し方向が Z なら恒等変換、水平なら鉛直面を張る
    same = rd.ocs_to_wcs([1.0, 2.0, 3.0], [0.0, 0.0, 1.0])
    check("押し出し方向 Z の OCS 変換は恒等", np.allclose(same, [1.0, 2.0, 3.0]),
          np.round(same, 6).tolist())
    axes = rd.ocs_axes([1.0, 0.0, 0.0])
    check("押し出し方向が水平なら OCS の 2 軸が鉛直面を張る",
          abs(np.dot(axes[0], [0, 0, 1])) < 1e-12 and abs(axes[1][2]) > 0.999,
          f"x軸={np.round(axes[0], 3).tolist()} y軸={np.round(axes[1], 3).tolist()}")
    check("OCS の 3 軸が正規直交",
          np.allclose(axes @ axes.T, np.eye(3), atol=1e-12))


# ---------------------------------------------------------------- 作図チェック
def test_check_model():
    print("\n[14] 作図ミスの自動チェック")
    model = rd.read_model(TEST_DXF, orient_normals="cad", verbose=False)
    issues = rd.check_model(model, verbose=False)
    errors = [i for i in issues if i["level"] == "error"]
    check("test.dxf にエラーは無い", not errors,
          f"{len(issues)} 件（{[i['level'] for i in issues]}）")

    # 受音点をわざと室外に置いたらエラーになる
    broken = rd.read_model(TEST_DXF, orient_normals="cad", verbose=False)
    broken.receiver_points = [np.array([100.0, 100.0, 100.0])]
    issues = rd.check_model(broken, verbose=False)
    check("室外の受音点をエラーにする",
          any(i["level"] == "error" and "外にあります" in i["message"] for i in issues))

    # 吸音率が引けないレイヤを指摘する（B-19 の状況）
    table = rd.AbsorptionTable()
    table["存在しない材料"] = np.zeros(8)
    issues = rd.check_model(model, absorption_table=table, verbose=False)
    check("吸音率が引けないレイヤを指摘する",
          any("吸音率が引けない" in i["message"] for i in issues))


# ---------------------------------------------------------------- 明瞭度の指標
def test_clarity():
    print("\n[15] 明瞭度の指標（C50 / C80 / D50 / Ts）")
    fs = 44100.0
    n = int(fs * 1.0)
    t = np.arange(n) / fs

    # 直接音のみ → 後続が無いので D50 = 1
    ir = np.zeros(n)
    ir[0] = 1.0
    result = rv.clarity_measures(t, ir, frequencies=[1000.0], verbose=False)
    check("直接音だけなら D50 = 1", abs(result["D50"][0] - 1.0) < 1e-6,
          f"D50={result['D50'][0]:.4f}")

    # 既知の指数減衰。D50 は解析的に 1 - exp(-2*50ms/tau)
    tau = 0.1
    noise = np.random.default_rng(0).normal(size=n)
    result = rv.clarity_measures(t, np.exp(-t / tau) * noise,
                                 frequencies=[1000.0], verbose=False)
    expected = 1.0 - np.exp(-2.0 * 0.050 / tau)
    check("指数減衰の D50 が解析値に近い",
          abs(result["D50"][0] - expected) < 0.10,
          f"計算 {result['D50'][0]:.3f} / 解析 {expected:.3f}")
    c50 = 10.0 * np.log10(result["D50"][0] / (1.0 - result["D50"][0]))
    check("C50 と D50 の関係 C50 = 10log10(D50/(1-D50))",
          abs(c50 - result["C50"][0]) < 0.05,
          f"{result['C50'][0]:.3f} vs {c50:.3f}")
    check("C80 は C50 より大きい", result["C80"][0] > result["C50"][0],
          f"C50={result['C50'][0]:.2f} C80={result['C80'][0]:.2f}")

    # 減衰が速いほど明瞭になる
    fast = rv.clarity_measures(t, np.exp(-t / 0.03) * noise,
                               frequencies=[1000.0], verbose=False)
    check("減衰が速いほど D50 が大きく Ts が小さい",
          fast["D50"][0] > result["D50"][0] and fast["Ts"][0] < result["Ts"][0],
          f"D50 {result['D50'][0]:.3f}→{fast['D50'][0]:.3f} / "
          f"Ts {result['Ts'][0] * 1000:.1f}→{fast['Ts'][0] * 1000:.1f} ms")


# ---------------------------------------------------------------- プロジェクト
def test_project():
    print("\n[16] プロジェクトの保存と読み込み")
    import tempfile
    import project as pj

    with tempfile.TemporaryDirectory() as folder:
        p = pj.Project(folder, dxf=TEST_DXF, absorption_kind="random",
                       rays=12345, radius=0.75, temperature=18.5, band_number=6,
                       volume=42.0)
        p.save()
        q = pj.Project.load(folder)
        same = all(getattr(p, k) == getattr(q, k)
                   for k in ("absorption_kind", "rays", "radius", "temperature",
                             "band_number", "volume"))
        check("条件が往復する", same, f"rays={q.rays} radius={q.radius}")
        check("解決したパスが元のファイルを指す",
              os.path.normcase(q.dxf_path) == os.path.normcase(os.path.abspath(TEST_DXF)),
              q.dxf_path)

        # 法線の反転指定。面数が合えば読め、合わなければ捨てる
        p.save_flipped_faces({1, 4, 7}, face_count=12, mode="shells")
        check("反転指定が往復する",
              pj.Project.load(folder).flipped_faces_for(12) == {1, 4, 7})
        check("面数が違えば反転指定を使わない",
              pj.Project.load(folder).flipped_faces_for(99) == set())

        # 面ごとの吸音材（materials.json）。normals.json と同じ約束にしてある
        p.save_face_materials({0: "カーペット", 3: "ガラス", 5: "カーペット"},
                              face_count=12)
        loaded = pj.Project.load(folder).face_materials_for(12)
        check("面ごとの吸音材が往復する",
              loaded == {0: "カーペット", 3: "ガラス", 5: "カーペット"}, str(loaded))
        check("面数が違えば吸音材の割り当ても使わない",
              pj.Project.load(folder).face_materials_for(99) == {})
        import json
        with open(os.path.join(folder, pj.MATERIALS_FILE), encoding="utf-8") as f:
            raw = json.load(f)
        check("materials.json は「材料 → 面番号」で持つ（人が読める）",
              raw["materials"]["カーペット"] == [0, 5], str(raw["materials"]))


# ------------------------------------------------------- 面グループと面ごとの吸音材
def test_face_groups():
    print("\n[25] 面グループ（同一平面パッチ）と面ごとの吸音材")

    # 単位立方体を三角形 12 枚で作る。同一平面で連結しているので **6 グループ**に
    # まとまるはず（3DSOLID を STL 経由で取り込んだときに、割られた壁を
    # 元の 1 枚に戻せるかどうかがこの機能の要）
    corner = np.array(list(itertools.product([0.0, 1.0], repeat=3)))
    # 立方体の 6 面。頂点番号は itertools.product の順（i → x=i>>2, y=(i>>1)&1, z=i&1）
    quads = [[0, 2, 6, 4], [1, 5, 7, 3],     # z=0 / z=1
             [0, 1, 3, 2], [4, 6, 7, 5],     # x=0 / x=1
             [0, 4, 5, 1], [2, 3, 7, 6]]     # y=0 / y=1
    triangles = []
    for indices in quads:
        a, b, c, d = (corner[i] for i in indices)
        triangles += [(a, b, c), (a, c, d)]
    # 法線は面から計算する（手で書くと取り違えるため）
    normals = np.array([rd.face_normal(*t) for t in triangles], dtype=float)
    group = rd.coplanar_groups(triangles, normals)
    check("立方体の三角形 12 枚 → 面 6 枚にまとまる",
          len(set(group.tolist())) == 6, f"{len(set(group.tolist()))} グループ")
    check("同じ面の 2 枚は同じグループ",
          all(group[2 * k] == group[2 * k + 1] for k in range(6)))
    check("向かい合う面は別グループ（平行でも連結していないので分かれる）",
          group[0] != group[2])

    # **法線を反転してもグループは変わらない**（幾何の性質。反転で選択単位が
    # 崩れてはいけない）
    flipped = rd.coplanar_groups(triangles, -normals)
    check("法線を反転してもグループ分けは同じ", np.array_equal(group, flipped))

    # 面ごとの吸音材がレイヤより優先されること
    table = {"レイヤ材": np.full(8, 0.10), "貼った材": np.full(8, 0.80)}
    plain = rd.read_model(TEST_DXF, absorption_table=table,
                          default_absorption=0.10, verbose=False)
    picked = {0: "貼った材", 2: "貼った材"}
    model = rd.read_model(TEST_DXF, absorption_table=table, default_absorption=0.10,
                          face_materials=picked, verbose=False)
    check("割り当てた面の吸音率が材料のものになる",
          np.allclose(model.mesh[0].absorption_coefficient, 0.80),
          str(np.round(model.mesh[0].absorption_coefficient[:3], 3).tolist()))
    check("割り当てた面の material が材料名になる（surface.csv が材料別になる）",
          model.mesh[0].material == "貼った材")
    check("割り当てていない面は元のまま",
          np.allclose(model.mesh[1].absorption_coefficient,
                      plain.mesh[1].absorption_coefficient))
    check("面ごとの割り当ての内訳を持っている",
          model.face_material_counts == {"貼った材": 2}, str(model.face_material_counts))
    check("レイヤ別の枚数は DXF のレイヤのまま（割り当てで書き換わらない）",
          model.layer_counts == plain.layer_counts)
    check("面ごとの DXF レイヤ名を別に持っている",
          model.face_layers == plain.face_layers
          and len(model.face_layers) == len(model.mesh))
    check("形状は変わらない（吸音率だけの差し替え）",
          np.allclose(model.mesh[0].vertexes, plain.mesh[0].vertexes))


# ------------------------------------------------------- 総表面積と容積（幾何だけ）
def test_area_and_volume():
    print("\n[27] 総表面積と容積（法線から＝発散定理）")

    # 2 x 3 x 4 m の直方体。表面積 2(6+8+12)=52 m2、体積 24 m3。
    # 法線は**内向き**（室として使うときの向き）にしておく
    size = np.array([2.0, 3.0, 4.0])
    corner = np.array(list(itertools.product([0.0, 1.0], repeat=3))) * size
    quads = [[0, 2, 6, 4], [1, 5, 7, 3], [0, 1, 3, 2],
             [4, 6, 7, 5], [0, 4, 5, 1], [2, 3, 7, 6]]
    box = []
    for indices in quads:
        a, b, c, d = (corner[i] for i in indices)
        box += [(a, b, c), (a, c, d)]
    centre = corner.mean(axis=0)
    inward = []
    for t in box:
        n = rd.face_normal(*t)
        inward.append(n if np.dot(n, centre - np.mean(t, axis=0)) > 0 else -n)
    inward = np.array(inward)

    check("総表面積", np.isclose(rd.surface_area(box), 52.0),
          f"{rd.surface_area(box):.3f} m2（正解 52）")
    check("容積（法線が内向き＝空気側）",
          np.isclose(rd.volume_from_normals(box, inward), 24.0),
          f"{rd.volume_from_normals(box, inward):.4f} m3（正解 24）")

    # **巻き順に依存しない**ことの確認。
    # ★全部まとめて裏返すと巻き順は「一貫したまま」なので符号が変わるだけ。
    #   一貫性を壊すには**一部だけ**裏返す必要がある（実モデルで起きるのはこちら）
    mixed = [(t[1], t[0], t[2]) if k % 2 == 0 else t for k, t in enumerate(box)]
    check("巻き順が一貫しなくなっている（前提の確認）",
          not rd.winding_is_consistent(mixed))
    check("★巻き順が崩れても容積は変わらない（法線を使うので）",
          np.isclose(rd.volume_from_normals(mixed, inward), 24.0),
          f"{rd.volume_from_normals(mixed, inward):.4f} m3")
    check("巻き順から出す signed_volume は崩れると合わない（だから法線を使う）",
          not np.isclose(abs(rd.signed_volume(mixed)), 24.0),
          f"{abs(rd.signed_volume(mixed)):.4f} m3（正解 24 にならない）")

    # 中に浮かせた 1 m 立方の物体。法線は物体の外＝空気側を向くので体積が引かれる
    inner_corner = np.array(list(itertools.product([0.0, 1.0], repeat=3))) + 0.5
    inner = []
    for indices in quads:
        a, b, c, d = (inner_corner[i] for i in indices)
        inner += [(a, b, c), (a, c, d)]
    inner_centre = inner_corner.mean(axis=0)
    outward = []
    for t in inner:
        n = rd.face_normal(*t)
        outward.append(n if np.dot(n, np.mean(t, axis=0) - inner_centre) > 0 else -n)
    both = box + inner
    normals = np.vstack([inward, np.array(outward)])
    check("★家具（宙に浮いた物体）の体積が自動で引かれる",
          np.isclose(rd.volume_from_normals(both, normals), 23.0),
          f"{rd.volume_from_normals(both, normals):.4f} m3（24 − 1 = 23）")

    # 開いた辺は本数だけでなく場所も返す（画面に赤で重ねて穴か T 字接合かを見るため）
    check("閉じた箱に開いた辺は無い", rd.open_edge_segments(box) == [])
    lid = box[:-2]      # 面を 1 枚（三角形 2 枚）欠かす
    segments = rd.open_edge_segments(lid)
    check("面を欠かすと開いた辺の場所が返る",
          len(segments) == rd.open_edge_count(lid) and len(segments) == 4,
          f"{len(segments)} 本")

    # 三角形の「形」の質。枚数は n 角形なら n-2 枚で最小なので減らせないが、
    # **細長い三角形（スリバー）はレイとの交差判定が丸めに左右されやすい**ので、
    # 耳刈り法は「最初に見つかった耳」ではなく**いちばん形の良い耳**を切る
    check("正三角形の最小角は 60 度",
          np.isclose(rd.triangle_min_angle([0, 0, 0], [1, 0, 0],
                                           [0.5, np.sqrt(3) / 2, 0]), 60.0))
    check("細長い三角形の最小角は小さい",
          rd.triangle_min_angle([0, 0, 0], [10, 0, 0], [10, 0.1, 0]) < 1.0,
          f"{rd.triangle_min_angle([0, 0, 0], [10, 0, 0], [10, 0.1, 0]):.3f}°")

    # **細長い面**で差が出る。10:1 の楕円上に 12 頂点を並べた輪。
    # 扇状に割ると端の三角形が潰れるが、質で選べばそうならない。
    # （正多角形ではどの耳も合同なので差が出ない。差が出るのは細長い形）
    ring = [np.array([10.0 * np.cos(a), 1.0 * np.sin(a), 0.0])
            for a in np.linspace(0.0, 2.0 * np.pi, 13)[:-1]]
    tris, info = rd.triangulate_polygon(ring)
    check("12 角形は 10 枚に分かれる（n-2 が最小で、これは減らせない）",
          len(tris) == 10, f"{len(tris)} 枚")
    worst = min(rd.triangle_min_angle(*t) for t in tris)
    fan = min(rd.triangle_min_angle(ring[0], ring[k], ring[k + 1])
              for k in range(1, len(ring) - 1))
    check("★扇状分割より形が良い（最小角が大きい）",
          worst > 5.0 * fan, f"耳刈り {worst:.2f}° / 扇状 {fan:.2f}°")
    check("面積は保たれる", info["area_error"] < 1.0e-9,
          f"相対誤差 {info['area_error']:.2e}")

    # read_model が総表面積とレイヤ別面積を持っていること
    m = rd.read_model(TEST_DXF, verbose=False)
    check("read_model が総表面積を持つ", m.surface_area > 0, f"{m.surface_area:.3f} m2")
    check("レイヤ別面積の合計が総表面積と一致",
          np.isclose(sum(m.layer_areas.values()), m.surface_area))
    check("容積が出ていて出し方も分かる",
          m.volume is not None and m.volume_source is not None,
          f"{m.volume:.3f} m3 / {m.volume_source}")
    check("三角形の形の統計を持っている",
          m.triangle_quality is not None and m.triangle_quality["median"] > 0,
          f"最小角の中央値 {m.triangle_quality['median']:.1f}°")


# ---------------------------------------------------------------- バンド数変換
def test_resample():
    print("\n[17] バンド数の変換（対数軸の補間）")
    six = [0.09, 0.36, 0.77, 0.97, 0.96, 0.95]
    m6 = ab.Material("6", six, kind="normal")
    up = m6.resample(8)
    check("6→8 は両端が隣のコピー（データが無いことの表れ）",
          abs(up[0] - six[0]) < 1e-12 and abs(up[-1] - six[-1]) < 1e-12,
          np.round(up, 3).tolist())
    check("6→8 で中身の 6 バンドは変わらない", np.allclose(up[1:7], six),
          np.round(up[1:7], 3).tolist())

    eight = [0.05, 0.09, 0.36, 0.77, 0.97, 0.96, 0.95, 0.90]
    down = ab.Material("8", eight, kind="normal").resample(6)
    check("8→6 は 63Hz と 8k を落とす", np.allclose(down, eight[1:7]),
          np.round(down, 3).tolist())
    check("同じバンド数なら素通し", np.allclose(m6.resample(6), six))

    # 中心周波数がずれた表からの内挿。log2 f 上で線形なら結果も線形になる
    source = np.array([125.0, 500.0, 2000.0])
    values = np.array([0.10, 0.30, 0.50])       # log2 f について線形
    target = np.array([250.0, 1000.0])
    got = np.interp(np.log2(target), np.log2(source), values)
    check("log2 f について線形なデータは中間点で線形補間される",
          np.allclose(got, [0.20, 0.40]), np.round(got, 4).tolist())


# ---------------------------------------------------------------- 伝搬方向
def test_direction():
    print("\n[18] 伝搬方向の集計（G-5）")
    import plots

    # 頭の正面を 0 とした相対方位に直せているか
    edges, totals = plots.direction_histogram([[1.0, 0, 0]], [1.0], head_azimuth=0.0)
    check("正面(+X)から来る音は 0° の区間に入る", int(np.argmax(totals)) == 0)

    edges, totals = plots.direction_histogram([[0, 1.0, 0]], [1.0], head_azimuth=0.0)
    check("左(+Y)から来る音は 90° の区間に入る",
          abs(np.rad2deg(edges[int(np.argmax(totals))]) - 90.0) < 1e-9,
          f"{np.rad2deg(edges[int(np.argmax(totals))]):.1f}°")

    edges, totals = plots.direction_histogram([[0, 1.0, 0]], [1.0], head_azimuth=90.0)
    check("頭を +Y に向けると、その音が正面扱いになる",
          int(np.argmax(totals)) == 0)

    # 上下は見ない（水平面へ投影する）
    slanted = [[1.0, 0.0, 5.0]]      # ほぼ真上から来るが方位は +X
    edges, totals = plots.direction_histogram(slanted, [1.0], head_azimuth=0.0)
    check("上下の傾きは方位に影響しない（水平面へ投影）",
          int(np.argmax(totals)) == 0)

    # 直接音（反射 0 回）の方位が、音源と受音点の位置から出る値と一致するか
    source, receiver = np.array([7.5, 7.0, 1.5]), np.array([4.5, 9.0, 1.2])
    toward_source = source - receiver
    want = np.rad2deg(np.arctan2(toward_source[1], toward_source[0])) % 360.0
    unit = toward_source / np.linalg.norm(toward_source)
    edges, totals = plots.direction_histogram([unit], [1.0], head_azimuth=0.0)
    got = np.rad2deg(edges[int(np.argmax(totals))])
    check("直接音の方位が幾何どおりの区間に入る",
          abs(((got + 5.0) - want + 180.0) % 360.0 - 180.0) < 5.0,
          f"区間 {got:.0f}° / 幾何 {want:.1f}°")

    # エネルギーが方位ごとに足し合わされるか
    edges, totals = plots.direction_histogram([[1.0, 0, 0], [1.0, 0, 0], [0, 1.0, 0]],
                                              [2.0, 3.0, 7.0], head_azimuth=0.0)
    check("同じ方位のエネルギーが合算される",
          abs(totals[0] - 5.0) < 1e-12 and abs(totals.sum() - 12.0) < 1e-12,
          f"正面 {totals[0]:.1f} / 合計 {totals.sum():.1f}")


# ---------------------------------------------------------------- モード分布
def test_modes():
    print("\n[19] モード分布（G-6）")
    import plots

    c = 343.0
    # 1 m 立方体の固有周波数を手計算と突き合わせる
    got, orders = plots.room_modes([1.0, 1.0, 1.0], 400.0, c)
    want = np.sort([0.5 * c * np.sqrt(a * a + b * b + d * d)
                    for a in range(4) for b in range(4) for d in range(4)
                    if (a or b or d) and 0.5 * c * np.sqrt(a * a + b * b + d * d) <= 400.0])
    check("1 m 立方体の固有周波数が手計算と一致",
          len(got) == len(want) and np.allclose(got, want),
          f"{len(got)} 個 / 最低次 {got[0]:.2f} Hz")
    check("最低次の軸モードが c/2L", abs(got[0] - 0.5 * c) < 1e-9,
          f"{got[0]:.3f} / 期待 {0.5 * c:.3f}")

    # (0,0,0) は音場ではないので含めない
    check("次数 (0,0,0) は含まれない", not np.any(np.all(orders == 0, axis=1)))
    check("すべて上限以下", got.max() <= 400.0 + 1e-9)

    # 細長い室では軸モードが低いほうへ寄る
    long_room, _ = plots.room_modes([10.0, 2.0, 2.0], 100.0, c)
    check("長い室ほど最低次が低い", long_room[0] < got[0],
          f"{long_room[0]:.2f} Hz（10 m の軸モード {0.5 * c / 10.0:.2f} Hz）")

    # スペクトルの山の検出：合成した山をそのまま拾えるか
    f = np.linspace(20.0, 200.0, 1800)
    level = np.full_like(f, -40.0)
    for centre in (40.0, 80.0, 150.0):
        level += 30.0 * np.exp(-0.5 * ((f - centre) / 2.0) ** 2)
    found = plots.spectrum_peaks(f, level)
    check("合成した 3 つの山を拾える", len(found) == 3, f"{len(found)} 個")
    check("拾った位置が合っている",
          np.allclose(np.round(f[found]), [40.0, 80.0, 150.0], atol=1.0),
          np.round(f[found], 1).tolist())

    flat = np.zeros_like(f)
    check("平坦なら山は拾わない", len(plots.spectrum_peaks(f, flat)) == 0)

    # シュレーダー周波数 2000√(T/V)
    check("シュレーダー周波数が 2000√(T/V)",
          abs(plots.schroeder_frequency(0.5, 100.0) - 2000.0 * np.sqrt(0.5 / 100.0)) < 1e-9,
          f"{plots.schroeder_frequency(0.5, 100.0):.1f} Hz")
    check("容積が無ければ None", plots.schroeder_frequency(0.5, None) is None)

    # 数え上げと近似式（体積・面積・辺長の 3 項）が近いこと
    lengths = np.array([8.5, 14.6, 3.0])
    counted, _ = plots.room_modes(lengths, 200.0, c)
    lx, ly, lz = lengths
    approx = ((4.0 * np.pi / 3.0) * (lx * ly * lz) * (200.0 / c) ** 3
              + (np.pi / 4.0) * 2.0 * (lx * ly + ly * lz + lz * lx) * (200.0 / c) ** 2
              + 4.0 * (lx + ly + lz) / 8.0 * (200.0 / c))
    check("累積モード数が近似式と 5% 以内で一致",
          abs(len(counted) - approx) / approx < 0.05,
          f"数え上げ {len(counted)} / 近似 {approx:.0f}")


# ---------------------------------------------------------------- 積み上げ
def test_mode_buildup():
    print("\n[20] 経路差からのモードの積み上げ（G-6b）")
    import plots

    c = 343.0

    # ---- 1 次元室（平行 2 面）で軸モードが出るか ----
    # 虚音源は x = 2nL ± xs。経路差が波長の整数倍になる周波数で強め合うはず
    L, xs, xr = 4.0, 1.0, 3.0
    n = np.arange(-400, 401)
    distance = np.abs(xr - np.concatenate([2 * n * L + xs, 2 * n * L - xs]))
    distance = distance[distance > 1e-9]
    time = distance / c

    f, H = plots.pulse_spectrum(time, np.ones_like(time),
                                max_time=time.max() * 1.05, bin_rate=16000.0)
    keep = (f > 5.0) & (f <= 200.0)
    f, magnitude = f[keep], np.abs(H[keep])
    level = 20.0 * np.log10(np.maximum(magnitude, 1e-9))
    found = plots.spectrum_peaks(f, level, prominence=6.0)
    found = found[magnitude[found] > np.sqrt(len(time))]
    peaks = f[found]

    # 軸モード f = m*c/(2L)。ただし音源・受音点が節に乗る次数は出ない（正しい）
    axial = np.array([m * c / (2.0 * L) for m in range(1, 6)])
    shape = np.array([np.cos(m * np.pi * xs / L) * np.cos(m * np.pi * xr / L)
                      for m in range(1, 6)])
    want = axial[np.abs(shape) > 1e-9]
    want = want[want <= 200.0]
    hit = [np.min(np.abs(peaks - w)) < 1.0 for w in want]
    check("平行 2 面の軸モードが山として出る",
          len(peaks) and all(hit),
          f"期待 {np.round(want, 1)} / 検出 {np.round(peaks, 1)}")

    # xs=1, xr=3, L=4 では 2 次（85.8 Hz）が両方とも節になるので出てはいけない
    node = axial[np.abs(shape) <= 1e-9]
    check("節に当たる次数は出ない（モード形状まで再現している）",
          all(np.min(np.abs(peaks - w)) > 1.0 for w in node),
          f"節 {np.round(node, 1)} Hz")

    # ---- 「重なった本数」がそのまま値になるか ----
    # 同じ時刻のパルスを k 本置けば、どの周波数でも |H| = k
    for k in (1, 2, 5):
        f2, H2 = plots.pulse_spectrum(np.full(k, 0.01), np.ones(k),
                                      max_time=0.1, bin_rate=8000.0)
        check(f"同位相に {k} 本重なれば値は {k}",
              np.allclose(np.abs(H2), float(k)),
              f"最大 {np.abs(H2).max():.3f} / 最小 {np.abs(H2).min():.3f}")

    # 半波長ずれた 2 本は打ち消し合う（位相を見ている証拠）
    half = 0.5 / 100.0                      # 100 Hz の半周期
    f3, H3 = plots.pulse_spectrum(np.array([0.01, 0.01 + half]), np.ones(2),
                                  max_time=0.2, bin_rate=64000.0)
    at100 = int(np.argmin(np.abs(f3 - 100.0)))
    check("半波長ずれた 2 本は打ち消し合う", abs(H3[at100]) < 0.02,
          f"|H(100 Hz)| = {abs(H3[at100]):.4f}")

    # ---- 重み（減衰）が効いているか ----
    # 完全反射は 1/d、吸音ありは sqrt(E)/d。1 本だけならその比がそのまま出る
    f4, H4 = plots.pulse_spectrum(np.array([0.02]), np.array([[1.0, 0.25]]),
                                  max_time=0.1, bin_rate=8000.0)
    check("重みの比がそのまま振幅の比になる",
          np.allclose(np.abs(H4[:, 1]) / np.abs(H4[:, 0]), 0.25),
          "0.25 倍")

    # ---- 周波数の刻みをまとめ直す（既定 1 Hz）----
    fine = np.arange(0.0, 20.0, 0.25)
    flat = np.full(len(fine), 3.0)
    coarse, value = plots.rebin_spectrum(fine, flat, 1.0)
    check("1 Hz 幅にまとめ直せる（刻みが 0.25 → 1 Hz）",
          abs(float(coarse[1] - coarse[0]) - 1.0) < 1e-9,
          f"{len(fine)} 点 → {len(coarse)} 点")
    check("一定値は値が変わらない（二乗平均なので）",
          np.allclose(value, 3.0))

    # **間引きではなく二乗平均**。ランダム位相の目安 √N が帯域幅に依らず保たれる、
    # というのが図の基準線を引き直さずに済む理由
    rng = np.random.default_rng(0)
    noise = rng.normal(size=len(fine)) + 1j * rng.normal(size=len(fine))
    _, wide = plots.rebin_spectrum(fine, noise, 1.0)
    check("ランダムな中身でも二乗平均の総量が保たれる",
          abs(np.sqrt(np.mean(wide ** 2)) - np.sqrt(np.mean(np.abs(noise) ** 2)))
          < 1e-9,
          f"まとめ前 {np.sqrt(np.mean(np.abs(noise)**2)):.4f} / "
          f"まとめ後 {np.sqrt(np.mean(wide**2)):.4f}")

    # 山を跨いで間引くと見落とすが、二乗平均なら残る
    spike = np.zeros(len(fine))
    spike[13] = 100.0                       # 3.25 Hz にだけ立つ鋭い山
    _, kept = plots.rebin_spectrum(fine, spike, 1.0)
    check("鋭い山は間引かれず残る", kept.max() > 40.0, f"最大 {kept.max():.1f}")

    check("刻みより細かい指定は何もしない",
          len(plots.rebin_spectrum(fine, flat, 0.1)[0]) == len(fine))

    # ---- バンドの割り当て ----
    bands = np.array([63.0, 125.0, 250.0, 500.0])
    index = plots._band_of(np.array([60.0, 90.0, 130.0, 400.0]), bands)
    check("周波数に最も近いオクターブバンドを当てる（対数距離）",
          list(index) == [0, 1, 1, 3], f"{list(index)}")

    # ---- 図が書けるか（吸音ありと完全反射で差が出るか）----
    import tempfile
    energy = np.column_stack([np.full(len(time), 0.5)] * 4)
    path = os.path.join(tempfile.mkdtemp(prefix="geosim_"), "mode_buildup.png")
    plots.mode_buildup(path, time, energy, distance=distance,
                       frequencies=bands, lengths=[L, 3.0, 2.5],
                       volume=L * 3.0 * 2.5, reverberation_time=0.5)
    check("mode_buildup.png が書ける", os.path.exists(path))


# ---------------------------------------------------------------- 図の描き直し
def test_redraw():
    print("\n[21] 保存済みの結果から図を描き直す")
    import shutil
    import tempfile

    import project as pj
    import run_project

    folder = tempfile.mkdtemp(prefix="geosim_redraw_")
    project = pj.Project(folder, dxf=TEST_DXF, band_number=6, rays=2000,
                         nref=8, radius=0.3, max_time=0.5, statistical=True,
                         volume=6.0)
    if os.path.exists(ABSORPTION):
        project.absorption_csv = ABSORPTION
    project.ensure_dirs()

    # まず本番と同じ手順で 1 回だけ計算し、結果 CSV を作る
    run_project.run(project, verbose=False, make_figures=True)
    before = sorted(os.listdir(project.path(pj.FIGURE_DIR)))
    check("計算すると図が書き出される", len(before) > 0, f"{len(before)} 枚")

    # 図を全部消してから、**計算し直さずに**描き直せるか
    shutil.rmtree(project.path(pj.FIGURE_DIR))
    stamps = {name: os.path.getmtime(project.path(pj.RESULT_DIR, name))
              for name in os.listdir(project.path(pj.RESULT_DIR))}
    written = run_project.redraw(project, verbose=False)
    after = sorted(os.listdir(project.path(pj.FIGURE_DIR)))
    check("描き直しで同じ図が揃う", after == before,
          f"{len(after)} 枚 / 足りない {sorted(set(before) - set(after))}")
    check("描き直しは結果 CSV を書き換えない",
          all(os.path.getmtime(project.path(pj.RESULT_DIR, n)) == t
              for n, t in stamps.items()),
          f"{len(stamps)} ファイル")
    check("書き出したファイルの一覧が返る", len(written) == len(after))

    # 結果が無いフォルダでは、黙って空の図を作らずに知らせる
    empty = pj.Project(tempfile.mkdtemp(prefix="geosim_empty_"), dxf=TEST_DXF)
    empty.ensure_dirs()
    try:
        run_project.redraw(empty, verbose=False)
        check("結果が無ければエラーにする", False, "例外が出なかった")
    except FileNotFoundError:
        check("結果が無ければエラーにする（黙って空の図を作らない）", True)

    shutil.rmtree(folder, ignore_errors=True)
    shutil.rmtree(empty.folder, ignore_errors=True)


# ---------------------------------------------------------------- 画面の保存
def test_capture():
    print("\n[22] 画面の保存とウィンドウのタイトル（G-12）")
    import shutil
    import tempfile

    import project as pj
    import view_model_gui as vg

    # ---- ウィンドウのタイトル ----
    check("ASCII ならそのまま通す", vg.ascii_tag("JR") == "JR")
    check("日本語は通さない（化けるので）", vg.ascii_tag("研修室") == "")
    check("空や空白は通さない",
          vg.ascii_tag("") == "" and vg.ascii_tag("   ") == "")

    # **物件名は入れない**（ユーザー判断）。画面の種類が分かればよい
    for screen, want in (("normals", "法線の確認"),
                         ("rays", "音線・音粒子"),
                         ("directions", "音線の飛び方"),
                         ("model", "モデルビューア")):
        bar, fallback = vg.window_titles(screen)
        check(f"{screen} の題が「{want}」だけ（物件名を入れない）",
              bar == want and "JR" not in fallback and "研修室" not in bar,
              f"{bar} / {fallback}")
    check("予備はすべて英字",
          all(vg.ascii_tag(f) == f for _, f in vg.WINDOW_TITLES.values()),
          str([f for _, f in vg.WINDOW_TITLES.values()]))
    check("知らない画面は渡した題を使う",
          vg.window_titles("なにか", "モデル") == ("モデル", "geosim"),
          str(vg.window_titles("なにか", "モデル")))
    check("知らない画面でも英字ならそれを予備に",
          vg.window_titles("なにか", "viewer") == ("viewer", "viewer"))

    # ---- キーの割り当て（★ここを間違えると画面が落ちる）----
    # `e` に数値入力を割り当てていたせいで、押すと入力ダイアログと同時に
    # **VTK の終了処理（ExitEvent）も走り**、値を入れた瞬間に画面が閉じていた
    # （ユーザー報告 2026-08-17）。同じ間違いを繰り返さないための番人。
    check("VTK の予約キーに e と q が入っている",
          "e" in vg.VTK_RESERVED_KEYS and "q" in vg.VTK_RESERVED_KEYS,
          str(sorted(vg.VTK_RESERVED_KEYS)))
    check("★数値入力のキーが VTK の予約キーでない",
          vg.VALUE_INPUT_KEY not in vg.VTK_RESERVED_KEYS,
          f"いまは {vg.VALUE_INPUT_KEY!r}")
    for key, what in (("g", "画像の保存"), ("b", "動画の保存"),
                      ("z", "視点(上)"), ("x", "視点(正面)"), ("c", "視点(横)"),
                      ("v", "視点(等角)"), ("n", "法線矢印"), ("o", "不透明度の対象"),
                      ("m", "モデル表示"), ("a", "自動判定"), ("d", "CAD の巻き順"),
                      ("i", "全反転")):
        check(f"{what}の {key!r} が予約キーでない", key not in vg.VTK_RESERVED_KEYS)

    # ---- 連番（撮るたびに増える。上書きしない）----
    folder = tempfile.mkdtemp(prefix="geosim_shot_")
    first = vg.next_free_path(folder, "法線")
    check("1 枚目は _01", os.path.basename(first) == "法線_01.png",
          os.path.basename(first))
    check("撮る前は同じ名前を返す", vg.next_free_path(folder, "法線") == first)
    open(first, "w").close()
    second = vg.next_free_path(folder, "法線")
    check("2 枚目は _02（上書きしない）",
          os.path.basename(second) == "法線_02.png", os.path.basename(second))
    check("種類が違えば別の連番",
          os.path.basename(vg.next_free_path(folder, "音線")) == "音線_01.png")

    # ---- プロジェクト側 ----
    project = pj.Project(tempfile.mkdtemp(prefix="geosim_proj_"), name="研修室")
    project.ensure_dirs()
    shots = project.screenshot_dir()
    check("画面の置き場が 図/画面/", os.path.isdir(shots)
          and os.path.basename(shots) == "画面"
          and os.path.basename(os.path.dirname(shots)) == "図", shots)
    # ★ここが肝。**計算し直しても手で撮った画像は消えない**こと。
    #   `clear_results()` は 図/ の PNG を消すので、同じ場所に置くと巻き添えになる
    manual = os.path.join(shots, "法線_01.png")
    open(manual, "w").close()
    auto = project.figure_path("modes.png")
    open(auto, "w").close()
    project.clear_results(verbose=False)
    check("計算し直すと自動生成の図は消える", not os.path.exists(auto))
    check("★手で撮った画像は消えない（図/画面/ に分けてあるから）",
          os.path.exists(manual))

    shutil.rmtree(folder, ignore_errors=True)
    shutil.rmtree(project.folder, ignore_errors=True)


# ------------------------------------------------------------ 反射（scalar 一致）
def test_reflection_vectorised():
    print("\n[24] ベクトル化した反射（scalar 版との一致）")

    # 音線追跡は速度のため反射を束ごと展開して書いている（loop_reflectionmesh 191行）。
    # scalar 版は `sound_ray.reflection_generator()`。**両者が一致することを押さえる**
    # （交差判定と同じ方針。CLAUDE.md「scalar 版は消さない」）。
    rng = np.random.default_rng(3)
    directions = rng.normal(size=(200, 3))
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    normals = rng.normal(size=(200, 3))
    normals /= np.linalg.norm(normals, axis=1)[:, None]

    # ベクトル化版（loop_reflectionmesh と同じ式）
    reflected = directions - 2.0 * np.einsum(
        "ij,ij->i", directions, normals)[:, None] * normals
    reflected /= np.linalg.norm(reflected, axis=1)[:, None]

    scalar = np.array([sr.reflection_generator(d, n)
                       for d, n in zip(directions, normals)])
    check("反射ベクトルが scalar 版と一致",
          np.allclose(reflected, scalar, atol=1e-12),
          f"最大差 {np.abs(reflected - scalar).max():.3e}")

    # 反射の性質そのものも確かめる（式を写し間違えていないか）
    cos_in = np.abs(np.einsum("ij,ij->i", directions, normals))
    cos_out = np.abs(np.einsum("ij,ij->i", reflected, normals))
    check("入射角と反射角が等しい", np.allclose(cos_in, cos_out, atol=1e-12),
          f"最大差 {np.abs(cos_in - cos_out).max():.3e}")
    check("長さが 1 に保たれる",
          np.allclose(np.linalg.norm(reflected, axis=1), 1.0, atol=1e-12))
    # 法線方向の成分だけが符号反転し、面に沿う成分は変わらない
    tangential_in = directions - np.einsum("ij,ij->i", directions, normals)[:, None] * normals
    tangential_out = reflected - np.einsum("ij,ij->i", reflected, normals)[:, None] * normals
    check("面に沿う成分は変わらない",
          np.allclose(tangential_in, tangential_out, atol=1e-12))


# ---------------------------------------------------------------- 表の向き
def test_table():
    print("\n[23] 表の並べ方（周波数は横）")
    import csv
    import shutil
    import tempfile

    import table as tb

    folder = tempfile.mkdtemp(prefix="geosim_table_")
    frequencies = [125.0, 250.0, 500.0, 1000.0]
    rows = {"EDT_s": [0.9, 0.6, 0.5, 0.3], "T30_s": [0.8, 0.5, 0.4, 0.2]}

    path = tb.write_frequency_table(os.path.join(folder, "rt.csv"),
                                    frequencies, rows)
    with open(path, encoding="utf-8-sig", newline="") as f:
        table = list(csv.reader(f))

    # ★ここが共通ルール。1 行目が周波数、2 行目以降が指標
    check("1 行目が周波数（横に並ぶ）",
          table[0] == ["項目", "125", "250", "500", "1000"], str(table[0]))
    check("2 行目以降が指標", [r[0] for r in table[1:]] == ["EDT_s", "T30_s"],
          str([r[0] for r in table[1:]]))
    check("Excel で開けるよう BOM 付き",
          open(path, "rb").read(3) == b"\xef\xbb\xbf")

    got_f, got_rows = tb.read_frequency_table(path)
    check("書いたものを読み戻せる",
          np.allclose(got_f, frequencies)
          and all(np.allclose(got_rows[k], v) for k, v in rows.items()))

    # NaN は空欄にする（Excel で「---」より扱いやすい）
    nan_path = tb.write_frequency_table(os.path.join(folder, "nan.csv"),
                                        frequencies,
                                        {"x": [1.0, np.nan, 3.0, 4.0]})
    with open(nan_path, encoding="utf-8-sig", newline="") as f:
        body = list(csv.reader(f))[1]
    check("NaN は空欄で書く", body[2] == "", repr(body))
    check("空欄は NaN で読み戻る",
          np.isnan(tb.read_frequency_table(nan_path)[1]["x"][1]))

    # **古い縦向きのファイルも読めること**（向きを変える前のプロジェクトが開けるように）
    old = os.path.join(folder, "old.csv")
    with open(old, "w", encoding="utf-8", newline="") as f:
        f.write("frequency_hz,EDT_s,T30_s\n125,0.9,0.8\n250,0.6,0.5\n")
    old_f, old_rows = tb.read_frequency_table(old)
    check("★古い縦向きの CSV も読める",
          np.allclose(old_f, [125, 250])
          and np.allclose(old_rows["EDT_s"], [0.9, 0.6])
          and np.allclose(old_rows["T30_s"], [0.8, 0.5]),
          f"{list(old_rows)}")

    check("無いファイルは (None, {}) を返す",
          tb.read_frequency_table(os.path.join(folder, "無い.csv")) == (None, {}))

    # 要素数が合わなければ黙って書かない
    try:
        tb.write_frequency_table(os.path.join(folder, "ng.csv"), frequencies,
                                 {"x": [1.0, 2.0]})
        check("要素数が合わなければエラー", False, "例外が出なかった")
    except ValueError:
        check("要素数が合わなければエラー（黙って壊れた表を書かない）", True)

    # 1 行が周波数以外のときの列名。**番号ではなく周波数そのものを入れる**
    check("バンドの列名が周波数入り",
          tb.band_column("alpha", 125.0) == "alpha_125Hz"
          and tb.band_column("decay", 1000.0, "db") == "decay_1000Hz_db",
          tb.band_column("decay", 1000.0, "db"))

    shutil.rmtree(folder, ignore_errors=True)


# ---------------------------------------------------------------- 音線の絞り込み
def test_ray_filter():
    print("\n[26] 音線の絞り込み（注目したい音線を残す）")
    import ray_filter as rfl

    # 解析的に答えが分かる軌跡を手で作る。
    # 音源 (0,0,0) から 3 本、それぞれ +x / +y / +z へ真っすぐ 10 m
    class FakeLog:
        pass

    log = FakeLog()
    log.ray_count = 3
    log.node_counts = np.array([2, 2, 2])
    log.pad_nodes = np.array([
        [[0.0, 0, 0], [10.0, 0, 0]],
        [[0.0, 0, 0], [0, 10.0, 0]],
        [[0.0, 0, 0], [0, 0, 10.0]],
    ])
    log.directions = np.array([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]])
    log.mesh_ids = np.array([5, 7, 5])
    log.mesh_offsets = np.array([0, 1, 2, 3])
    log.ray_indexes = np.array([0, 1, 2])
    log.total_distance = np.array([10.0, 10.0, 10.0])
    log.received = np.array([True, False, True])
    log.terminations = np.array(["受音", "nref", "受音"])
    log.sound_velocity = 343.0

    check("音源は 1 点目から取れる",
          np.allclose(rfl.source_point(log), [0, 0, 0]))

    # ---- 線分までの距離（節点までではない）----
    # (5, 1, 0) は 1 本目の**途中**の横 1 m。節点だけ見ると 5 m 以上に見える
    d = rfl.ray_distances(log, [5.0, 1.0, 0.0])
    check("線分までの距離で測る（節点までではない）",
          abs(d[0] - 1.0) < 1e-12, f"1 本目 {d[0]:.6f} m（期待 1.0）")
    # +y の音線は最近点が (0,1,0) なので 5.0。
    # +z の音線は最近点が原点（t が 0 に丸まる）なので hypot(5,1)
    check("他の音線の距離も解析解と一致",
          abs(d[1] - 5.0) < 1e-12 and abs(d[2] - np.hypot(5.0, 1.0)) < 1e-12,
          f"+y {d[1]:.4f}（期待 5.0）/ +z {d[2]:.4f}（期待 {np.hypot(5.0,1.0):.4f}）")

    # 線分の外側に落ちる点は端点までの距離になる（t を 0〜1 に丸めている）
    beyond = rfl.ray_distances(log, [13.0, 0.0, 0.0])
    check("線分の先にある点は端点までの距離",
          abs(beyond[0] - 3.0) < 1e-12, f"{beyond[0]:.6f} m（期待 3.0）")

    # ---- 近くを通る音線 ----
    near = rfl.near_point(log, [5.0, 1.0, 0.0], 1.5)
    check("半径 1.5 m で 1 本目だけ残る", list(near) == [0], str(list(near)))
    near = rfl.near_point(log, [5.0, 1.0, 0.0], 6.0)
    check("半径 6 m なら 3 本とも残る", list(near) == [0, 1, 2], str(list(near)))
    check("いちばん近い 1 本を選べる",
          rfl.nearest_ray(log, [5.0, 1.0, 0.0]) == 0)
    check("対象を絞ってから探せる（1 本目を除くと 2 本目か 3 本目）",
          rfl.nearest_ray(log, [5.0, 1.0, 0.0], index=[1, 2]) in (1, 2))

    # ---- 出射方向で絞る ----
    angles = rfl.launch_angles(log, [1.0, 0.0, 0.0])
    check("出射方向とのなす角が解析解と一致",
          np.allclose(angles, [0.0, 90.0, 90.0]), str(np.round(angles, 3)))
    check("半角 10°なら +x の 1 本だけ",
          list(rfl.in_direction(log, [1.0, 0, 0], 10.0)) == [0])
    check("半角 95°なら 3 本とも",
          list(rfl.in_direction(log, [1.0, 0, 0], 95.0)) == [0, 1, 2])
    check("長さ 0 の方向はエラーにする",
          _raises(ValueError, rfl.launch_angles, log, [0.0, 0.0, 0.0]))

    # 音源からクリック点への向き
    check("音源→点の向きが単位ベクトル",
          np.allclose(rfl.direction_to(log, [3.0, 0.0, 0.0]), [1, 0, 0]))
    check("音源と同じ位置はエラーにする",
          _raises(ValueError, rfl.direction_to, log, [0.0, 0.0, 0.0]))

    # ---- 方位角・仰角で方向を決める（クリックが難しいという指摘への対応）----
    # 約束は project.head_azimuth と同じ：0° = +X、反時計回り、仰角 0 が水平
    for azimuth, elevation, want in ((0, 0, [1, 0, 0]), (90, 0, [0, 1, 0]),
                                     (180, 0, [-1, 0, 0]), (270, 0, [0, -1, 0]),
                                     (0, 90, [0, 0, 1]), (0, -90, [0, 0, -1])):
        got = rfl.direction_from_angles(azimuth, elevation)
        check(f"方位 {azimuth:3d}° 仰角 {elevation:+3d}° → {want}",
              np.allclose(got, want, atol=1e-12), str(np.round(got, 3)))
    check("仰角 45°で水平成分と鉛直成分が等しい",
          abs(rfl.direction_from_angles(0, 45)[0]
              - rfl.direction_from_angles(0, 45)[2]) < 1e-12)
    check("いつでも単位ベクトル",
          all(abs(np.linalg.norm(rfl.direction_from_angles(a, e)) - 1.0) < 1e-12
              for a in (0, 37, 190, 359) for e in (-90, -30, 0, 30, 90)))
    # 逆変換（点を拾ったときにスライダの値を合わせるのに使う）
    for azimuth, elevation in ((37.0, 22.0), (300.0, -55.0), (0.0, 0.0)):
        back = rfl.angles_from_direction(
            rfl.direction_from_angles(azimuth, elevation))
        check(f"({azimuth}, {elevation}) を往復できる",
              abs(back[0] - azimuth) < 1e-9 and abs(back[1] - elevation) < 1e-9,
              f"{back[0]:.6f}, {back[1]:.6f}")
    check("方位角は 0〜360 に収める",
          0.0 <= rfl.angles_from_direction([-1.0, -1.0, 0.0])[0] < 360.0,
          f"{rfl.angles_from_direction([-1.0, -1.0, 0.0])[0]:.1f}°")
    check("長さ 0 は逆変換もエラーにする",
          _raises(ValueError, rfl.angles_from_direction, [0.0, 0.0, 0.0]))

    # ---- 反射回数で見る範囲を制限できる（絞り込みが効くようにするため）----
    # 1 本目を「途中で折れる」形にして、2 区間目だけが点の近くを通るようにする
    bent = FakeLog()
    bent.ray_count = 1
    bent.node_counts = np.array([3])
    bent.pad_nodes = np.array([[[0.0, 0, 0], [10.0, 0, 0], [10.0, 10.0, 0]]])
    bent.directions = np.array([[1.0, 0, 0]])
    check("全区間を見れば近い",
          abs(rfl.ray_distances(bent, [10.0, 5.0, 0.0])[0]) < 1e-12)
    check("★1 区間目までに限れば遠い（描いている範囲に合わせられる）",
          abs(rfl.ray_distances(bent, [10.0, 5.0, 0.0], max_reflection=1)[0] - 5.0)
          < 1e-12,
          f"{rfl.ray_distances(bent, [10.0, 5.0, 0.0], max_reflection=1)[0]:.3f} m")

    # ---- 面で絞る・説明文 ----
    check("その面で反射した音線を引ける",
          list(rfl.through_face(log, 5)) == [0, 2], str(list(rfl.through_face(log, 5))))
    text = rfl.describe_ray(log, 0)
    check("1 本の説明に反射回数・経路長・受音の有無が入る",
          "1 回反射" in text and "10.00 m" in text and "受音 した" in text,
          text.replace("\n", " / "))


def _raises(kind, func, *args):
    try:
        func(*args)
    except kind:
        return True
    except Exception:
        return False
    return False


def main():
    print("geosim 数値検証")
    print(f"  Python {sys.version.split()[0]} / numpy {np.__version__}")
    if not os.path.exists(ABSORPTION):
        print(f"  ※ {ABSORPTION} が無いので既定の吸音率で走ります")

    for fn in (test_read_dxf, test_soundray_generator, test_energy_decay,
               test_backtrace, test_image_sources, test_vectorised_geometry,
               test_atmosphere, test_absorption, test_impulse, test_reverberation,
               test_statistical_reverberation, test_ray_log,
               test_normals, test_check_model, test_clarity,
               test_project, test_resample, test_direction, test_modes,
               test_mode_buildup, test_redraw, test_capture, test_table,
               test_reflection_vectorised, test_face_groups,
               test_ray_filter, test_area_and_volume):
        fn()

    failed = [name for name, ok in _results if not ok]
    print(f"\n{'=' * 60}")
    print(f"  {len(_results) - len(failed)} / {len(_results)} 件 OK")
    for name in failed:
        print(f"  NG: {name}")
    print(f"{'=' * 60}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
