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
    print("\n[10] 統計残響式（Sabine / Eyring / Millington）")
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
    check("吸音率が一様なら Millington と Eyring が一致",
          np.abs(result["millington"] - result["eyring"]).max() < 1e-12)
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
    check("空気吸収を入れると残響時間が短くなる（高音ほど顕著）",
          with_air["sabine"][-1] < with_air["sabine"][0] < result["sabine"][0],
          f"63Hz {with_air['sabine'][0]:.4f} / 8kHz {with_air['sabine'][-1]:.4f} s")

    # 実際のモデルで面ごとに吸音率が違う場合
    real = rv.statistical_reverberation_from_model(model, atmosphere=air, verbose=False)
    check("DxfModel から直接計算できる", real is not None)
    check("吸音率がばらつくと Millington < Eyring < Sabine",
          np.all(real["millington"] < real["eyring"])
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
               test_project, test_resample):
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
