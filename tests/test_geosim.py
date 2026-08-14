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
import loop_deleteredundancy as ld   # noqa: E402
import loop_noredundancy as ln       # noqa: E402
import impulse as ip                 # noqa: E402
import reverberation as rv           # noqa: E402

C = ln.SOUND_VELOCITY
TEST_DXF = os.path.join(ROOT, "test.dxf")
ABSORPTION = os.path.join(ROOT, "absorption.csv")

_results = []


def check(name, ok, detail=""):
    _results.append((name, bool(ok)))
    mark = "OK " if ok else "NG "
    print(f"  [{mark}] {name}" + (f"  {detail}" if detail else ""))
    return ok


def load_test_room():
    table = ABSORPTION if os.path.exists(ABSORPTION) else None
    model = rd.read_model(TEST_DXF, absorption_table=table, verbose=False)
    return (model, np.asarray(model.source_points[0]),
            np.asarray(model.receiver_points[0]))


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


# ---------------------------------------------------------------- ⑥ インパルス応答
def test_impulse():
    print("\n[6] インパルス応答の合成")
    # バンドパスが scipy と一致するか（E-6 の判断根拠）
    try:
        from scipy.signal import firwin
        worst = 0.0
        for numtaps, lo, hi in [(1024, 0.05, 0.1), (4096, 0.001, 0.002)]:
            mine = ip.filter_bandpass(numtaps, lo, hi)
            theirs = firwin(numtaps, [lo, hi], window="hamming",
                            pass_zero=False, scale=False)
            worst = max(worst, np.abs(mine - theirs).max() / np.abs(theirs).max())
        check("filter_bandpass が scipy.signal.firwin と一致", worst < 1e-12,
              f"最大相対誤差 {worst:.2e}")
    except ImportError:
        print("  [--] scipy が無いので firwin との比較は省略")

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

    # 6→32 展開表
    e32 = ip.expand_6_to_32(np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]]))
    want = {9: 1.0, 12: 2.0, 15: 3.0, 18: 4.0, 21: 5.0, 24: 6.0}
    check("6→32 バンド展開表が各オクターブ中心で正しい",
          all(e32[i, 0] == v for i, v in want.items()))

    # 単一パルス
    t0 = 0.01
    t, ir = ip.impulse_response(np.array([t0]), np.ones((1, 6)), verbose=False)
    dt = t[1] - t[0]
    delay = (len(ir) - 1) // 2
    peak = int(np.argmax(np.abs(ir)))
    check("単一パルスのピークが期待位置",
          abs(peak - (delay + int(round(t0 / dt)))) <= 1,
          f"{peak} / 期待 {delay + int(round(t0 / dt))}")

    t2, ir2 = ip.impulse_response(np.array([t0]), np.ones((1, 6)),
                                  compensate_filter_delay=True, verbose=False)
    peak2 = int(np.argmax(np.abs(ir2)))
    check("遅れ補正でピークが入力時刻に来る", abs(t2[peak2] - t0) < 2 * dt,
          f"{t2[peak2] * 1000:.4f} ms / 期待 {t0 * 1000:.4f} ms")

    # 距離減衰 1/r
    e32 = ip.apply_air_absorption(ip.expand_6_to_32(np.ones((1, 6))),
                                  np.array([t0]), C, np.zeros(32))
    hfp = ip.transfer_function(e32, np.array([t0]), 5, 1.0, C)
    check("伝達関数の直流成分が 1/(t*c) と一致",
          abs(abs(hfp[0, 0]) - 1.0 / (t0 * C)) < 1e-12,
          f"{abs(hfp[0, 0]):.9f} / {1.0 / (t0 * C):.9f}")


# ---------------------------------------------------------------- ⑦ 残響時間
def test_reverberation():
    print("\n[7] 残響時間（減衰率が既知の合成応答で検証）")
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
    d, tt = result["decay"][3], result["time"]
    mask = (d <= -5) & (d >= -35)
    slope = np.polyfit(tt[mask], d[mask], 1)[0]
    check("減衰曲線の傾きが理論と一致（1000 Hz）",
          abs(-60.0 / slope - 0.5) / 0.5 < 0.03,
          f"傾きから T60 = {-60.0 / slope:.4f} s / 理論 0.5000 s")
    check("フィルタの遅れが残っていない",
          tt[np.nonzero(d <= -1.0)[0][0]] < 0.05,
          f"-1 dB を切る時刻 {tt[np.nonzero(d <= -1.0)[0][0]] * 1000:.2f} ms")


def main():
    print("geosim 数値検証")
    print(f"  Python {sys.version.split()[0]} / numpy {np.__version__}")
    if not os.path.exists(ABSORPTION):
        print(f"  ※ {ABSORPTION} が無いので既定の吸音率で走ります")

    for fn in (test_read_dxf, test_soundray_generator, test_energy_decay,
               test_backtrace, test_image_sources, test_impulse,
               test_reverberation):
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
