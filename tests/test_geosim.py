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


def main():
    print("geosim 数値検証")
    print(f"  Python {sys.version.split()[0]} / numpy {np.__version__}")
    if not os.path.exists(ABSORPTION):
        print(f"  ※ {ABSORPTION} が無いので既定の吸音率で走ります")

    for fn in (test_read_dxf, test_soundray_generator, test_energy_decay,
               test_backtrace, test_image_sources, test_vectorised_geometry,
               test_atmosphere, test_absorption, test_impulse, test_reverberation):
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
