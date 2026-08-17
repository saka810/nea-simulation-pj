"""インパルス応答 → 減衰曲線・残響時間。元コード ipls2rt_fortran.f90。

## 元コードからの変更方針（2026-08-14）

| 処理 | 元コード | 本実装 |
|---|---|---|
| オクターブフィルタ | 自作 FIR（タップ数 = データ点数） | `scipy.signal.butter` + `sosfilt`（**IEC 61260** 準拠の Butterworth） |
| 畳み込み | 巡回畳み込み | IIR の逐次フィルタなので回り込み自体が起きない |
| 遅れの補正 | nn/2 の巡回シフト | 不要（因果フィルタなので遅れが小さい） |
| バンド数 | 6 固定 | 可変（既定 8） |

**元コードのやり方は使えない。** バンドパスが左右対称な直線位相 FIR なので、
入力の立ち上がりより前の時刻から応答が始まり（プリリンギング）、
巡回畳み込みではそれが**バッファ末尾に回り込む**。しかもタップ数 = nn なので
回り込みはバッファの半分に及ぶ。結果、減衰曲線に -20 dB 程度の**床**ができて
-35 dB まで落ちなくなり、減衰率が既知の応答で試すと
**T30 が理論値の 5〜19 倍**になるか、そもそも算出できなかった。

Butterworth の因果 IIR にすると、この問題は構造的に起きない。
FIR（線形畳み込み）と精度は同等で（誤差の中央値 1.0〜1.9 %）、計算量は桁違いに少ない。
FIR で処理したい場合は `method='fir'` を指定できる。

## 流れ

    1. インパルス応答をオクターブバンドに分ける
    2. Schroeder 積分（後ろ向きの積分）で減衰曲線を作る
    3. -5 dB と -35 dB を横切る時刻の差から T30 を求める
"""

import numpy as np
from scipy.signal import butter, fftconvolve, firwin, sosfilt

import absorption as ab
from absorption import DEFAULT_OCTAVE_BANDS
from atmosphere import Atmosphere

# 元コード 9〜10 行。T30 を求めるための評価区間
DB_MAX = -5.0
DB_MIN = -35.0

# 実務で使う残響指標の評価区間（開始 dB, 終了 dB）。
# 60 dB 減を厳密に見ることは実際にはほとんど無く、
# 減衰の直線部分を測って 60 dB 相当に外挿するのが普通。
#   EDT … 初期減衰時間。0〜-10 dB。**聴感上の響きの短さに近い**とされる。
#          初期反射の密度を反映するので、後部残響が同じでも室形状で差が出る
#   T20 … -5〜-25 dB。暗騒音が高い実測でも取りやすい
#   T30 … -5〜-35 dB。最も一般的（ISO 3382 の標準）
DECAY_MEASURES = {
    "EDT": (0.0, -10.0),
    "T20": (-5.0, -25.0),
    "T30": (-5.0, -35.0),
}

# IEC 61260 のクラス 1 オクターブフィルタに相当する次数（Butterworth 6 次）
FILTER_ORDER = 6

# method='fir' のときのタップ数
FIR_NUMTAPS = 4096


# ------------------------------------------------------------------------------
# 統計的残響式（Sabine / Eyring-Knudsen / Millington）
#
# 音線を飛ばさず、**室容積 V と各面の面積・吸音率だけ**から残響時間を見積もる。
# 拡散音場（音がどの方向からも等確率に来る）を前提にした古典的な式で、
# シミュレーション結果の妥当性を確かめる物差しになる。
#
# ★閉じた室が前提★ 開いた形状（一面だけの壁など）では容積が定義できないので使えない。
# ------------------------------------------------------------------------------

def triangle_areas(mesh):
    """各三角形の面積 [m^2]。外積の大きさの半分。"""
    v = np.array([np.asarray(m.vertexes, dtype=float) for m in mesh])
    return 0.5 * np.linalg.norm(np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]), axis=1)


def surface_summary(mesh, convert_to_random=True, warn=True):
    """材料ごとの面積と吸音率をまとめる。

    引数:
        convert_to_random : Mesh が持つ**垂直入射**吸音率を、統計式が前提とする
            **乱入射**吸音率へ Paris の式で変換するか（既定 True）。
            ★変換しないと吸音を過大評価する（垂直入射 0.34 は乱入射 0.50 に相当）。

    戻り値: dict
        'materials' {材料名: {'area', 'absorption'(nf,)}}
        'total_area' float
        'areas' (材料数,) / 'absorption' (材料数, nf)  … 式に渡しやすい形
    """
    areas = triangle_areas(mesh)
    materials = {}
    for face, area in zip(mesh, areas):
        entry = materials.setdefault(face.material,
                                     {"area": 0.0, "absorption": None})
        entry["area"] += float(area)
        alpha = np.atleast_1d(np.asarray(face.absorption_coefficient, dtype=float))
        if entry["absorption"] is None:
            entry["absorption"] = alpha
        elif warn and not np.allclose(entry["absorption"], alpha):
            print(f"[reverberation] 警告: 材料 {face.material!r} に異なる吸音率の面が"
                  f"混ざっています。最初の値を使います")

    names = sorted(materials)
    area_array = np.array([materials[n]["area"] for n in names])
    alpha_array = np.array([materials[n]["absorption"] for n in names])

    if convert_to_random:
        # 面ごとに変換すると同じ値を何度も計算するので、重複を除いてから変換する
        flat = alpha_array.ravel()
        unique, inverse = np.unique(flat, return_inverse=True)
        converted = ab.normal_to_random(unique)
        alpha_array = converted[inverse].reshape(alpha_array.shape)
        for i, name in enumerate(names):
            materials[name]["absorption"] = alpha_array[i]

    return {"materials": materials, "total_area": float(area_array.sum()),
            "names": names, "areas": area_array, "absorption": alpha_array}


def statistical_reverberation(mesh, volume, frequencies=None, atmosphere=None,
                              convert_to_random=True, include_air_absorption=True,
                              verbose=True):
    """Sabine / Eyring / Millington-Sette の残響式で残響時間を見積もる。

    引数:
        mesh       list[Mesh]  室形状（`read_dxffile` の出力）
        volume     float       室容積 [m^3]。`DxfModel.volume` の絶対値
        atmosphere Atmosphere | None  音速と空気吸収をここから取る
        convert_to_random  垂直入射 → 乱入射の変換を行うか（既定 True。上記参照）
        include_air_absorption  空気吸収の項 4mV を入れるか（既定 True）

    **3 つの式の違い**

    共通して `T = 24 ln(10) V / (c A)`（A = 等価吸音面積）で、A の作り方が違う。
    係数 24 ln(10) / c は c = 343 m/s のとき 0.161 になる（よく見る 0.161V/A の形）。
    ここでは**大気条件から求めた音速**を使うので 0.161 に固定していない。

    | 式 | 等価吸音面積 A | 性質 |
    |---|---|---|
    | Sabine | `S·ᾱ + 4mV` | 最も古典的。**吸音率が小さいとき（〜0.2）に妥当**。ᾱ→1 でも T が 0 にならない欠点 |
    | Eyring-Knudsen | `-S·ln(1-ᾱ) + 4mV` | 反射のたびに (1-ᾱ) 倍になると考える。**吸音率が大きいときはこちら** |
    | Millington | `-Σ Sᵢ·ln(1-αᵢ) + 4mV` | 面ごとに個別に対数を取る。αᵢ→1 の面があると発散する |

    ᾱ は面積で重み付けした平均吸音率 `Σ Sᵢαᵢ / S`。

    ★**名前について。** アイリングの式そのものは `-S ln(1-ᾱ)` までで、
      **空気吸収の項 `4mV` を足した形はヌードセンの寄与**なので
      `アイリング・ヌードセンの式` と呼ぶのが正しい（ユーザー指摘 2026-08-17）。
      `include_air_absorption=False` にしたときだけ素の `Eyring`。
      表示名は `statistical_labels()` が切り替える。

    ★**ミリントンの式は参考値**。`ミリントン・セッテの式`（Millington 1932 /
      Sette 1933）で、平均してから対数を取る Eyring と違い**面ごとに対数を取る**。
      吸音率が面ごとに大きく違うときの理屈は通っているが、
      **αᵢ→1 の面が 1 枚でもあると A が発散して T→0 になる**（開口が 1 つあるだけで
      残響ゼロという結論になってしまう）ため、実務ではまず使われない。
      研修室（吸音面 α=0.951）では Eyring-Knudsen の半分近い値が出る。

    戻り値: dict
        'frequencies' / 'volume' / 'total_area' / 'mean_free_path'
        'mean_absorption' (nf,) / 'equivalent_area' (nf,)
        'sabine' / 'eyring' / 'millington' 各 (nf,) [s]
        'surface'  … `surface_summary()` の結果
    """
    if atmosphere is None:
        atmosphere = Atmosphere()
    surface = surface_summary(mesh, convert_to_random=convert_to_random, warn=verbose)

    areas = surface["areas"]
    alpha = surface["absorption"]                    # (材料数, nf)
    total_area = surface["total_area"]
    if frequencies is None:
        frequencies = ab.octave_bands(alpha.shape[1])
    frequencies = np.asarray(frequencies, dtype=float)

    sound_velocity = atmosphere.sound_velocity
    constant = 24.0 * np.log(10.0) / sound_velocity  # c=343 のとき 0.1611

    # 空気吸収の項 4mV（m はエネルギーの減衰係数 [1/m]）
    air = (4.0 * atmosphere.absorption_coefficient(frequencies) * volume
           if include_air_absorption else np.zeros(len(frequencies)))

    mean_absorption = (areas @ alpha) / total_area          # (nf,)

    # ln(1-α) は α → 1 で発散する。1 に張り付いた材料があると Eyring/Millington は
    # 意味を持たなくなるので、そこは NaN にして知らせる
    with np.errstate(divide="ignore", invalid="ignore"):
        eyring_area = -total_area * np.log(np.where(mean_absorption < 1.0,
                                                    1.0 - mean_absorption, np.nan))
        millington_area = -(areas @ np.log(np.where(alpha < 1.0, 1.0 - alpha, np.nan)))

    def to_time(equivalent_area):
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(equivalent_area > 0.0,
                            constant * volume / equivalent_area, np.nan)

    result = {
        "frequencies": frequencies,
        "volume": float(volume),
        "total_area": total_area,
        "mean_free_path": 4.0 * volume / total_area,
        "sound_velocity": sound_velocity,
        "mean_absorption": mean_absorption,
        "equivalent_area": total_area * mean_absorption,
        "air_absorption_area": air,
        "sabine": to_time(total_area * mean_absorption + air),
        "eyring": to_time(eyring_area + air),
        "millington": to_time(millington_area + air),
        "surface": surface,
    }

    if verbose:
        print(f"[統計残響] 容積 {volume:.3f} m3 / 表面積 {total_area:.3f} m2 / "
              f"平均自由行程 {result['mean_free_path']:.3f} m / "
              f"音速 {sound_velocity:.2f} m/s")
        kind = "乱入射（垂直入射から変換）" if convert_to_random else "垂直入射のまま"
        print(f"[統計残響] 吸音率は{kind}。空気吸収 "
              f"{'あり' if include_air_absorption else 'なし'}")
        for name in surface["names"]:
            entry = surface["materials"][name]
            print(f"[統計残響]   {name:<20} 面積 {entry['area']:8.3f} m2  "
                  f"α {np.array2string(entry['absorption'], precision=3)}")
        labels = statistical_labels(result)
        print(f"[統計残響] {'周波数':>10}{'平均α':>9}{labels['sabine']:>10}"
              f"{labels['eyring']:>16}{labels['millington']:>12}")
        for i, fc in enumerate(frequencies):
            def cell(key, width):
                value = result[key][i]
                return "---".rjust(width) if np.isnan(value) else f"{value:{width}.3f}"
            print(f"[統計残響] {fc:9.0f}Hz{mean_absorption[i]:9.3f}"
                  f"{cell('sabine', 10)}{cell('eyring', 16)}{cell('millington', 12)}")

    return result


def statistical_labels(result=None, include_air_absorption=None):
    """表示に使う式の名前。**空気吸収を入れているかで呼び名が変わる。**

    アイリングの式そのものは `-S ln(1-ᾱ)` までで、**空気吸収の項 `4mV` を
    足した形はヌードセンの寄与**なので `アイリング・ヌードセンの式` と呼ぶ
    （ユーザー指摘 2026-08-17。それまで空気吸収込みなのに `Eyring` と出していた）。

    `result` は `statistical_reverberation()` の戻り値。
    `air_absorption_area` が 0 でなければ空気吸収込みと判断する。
    """
    if include_air_absorption is None:
        air = np.asarray((result or {}).get("air_absorption_area", 0.0))
        include_air_absorption = bool(np.any(air > 0.0))
    return {
        "sabine": "Sabine",
        "eyring": "Eyring-Knudsen" if include_air_absorption else "Eyring",
        "millington": "Millington",
        "air": include_air_absorption,
    }


def statistical_reverberation_from_model(model, **kwargs):
    """`read_dxffile.read_model()` の結果からそのまま統計残響式を計算する。

    閉じていない形状では容積が定義できないので、その場合は None を返して警告する。
    """
    if not getattr(model, "is_closed", False):
        print(f"[統計残響] 形状が閉じていないので容積が決まりません"
              f"（開いた辺 {getattr(model, 'open_edges', '?')} 本）。"
              f"統計残響式は使えません")
        return None
    return statistical_reverberation(model.mesh, abs(model.volume), **kwargs)


def write_statistical_reverberation(filename, result):
    """統計残響式の結果を CSV に保存する。"""
    header = ["frequency_hz", "mean_absorption", "equivalent_area_m2",
              "sabine_s", "eyring_s", "millington_s"]
    rows = np.column_stack([result["frequencies"], result["mean_absorption"],
                            result["equivalent_area"], result["sabine"],
                            result["eyring"], result["millington"]])
    np.savetxt(filename, rows, delimiter=",", header=",".join(header),
               comments="", fmt="%.12g")
    return filename


def schroeder_integral(x):
    """後ろ向きに積分した残留エネルギー。元コード 114〜117 行。

        D(t) = Σ_{τ >= t} x(τ)^2

    「時刻 t 以降にまだ残っているエネルギー」。
    インパルス応答の 2 乗をそのまま dB 表示すると反射音の干渉で激しく暴れるが、
    積分すると滑らかになる（Schroeder の逆積分法）。
    """
    return np.cumsum(x[::-1] ** 2)[::-1]


def octave_bandpass(signal, centre_frequency, sampling_frequency,
                    method="butter", order=FILTER_ORDER, numtaps=FIR_NUMTAPS):
    """オクターブバンド（中心周波数の 1/√2 〜 √2 倍）を切り出す。

    method='butter' … `scipy.signal.butter` の SOS を `sosfilt` で適用（既定）。
                      IEC 61260 のオクターブフィルタは Butterworth 系なので実務に沿う。
    method='fir'    … `scipy.signal.firwin` + 線形畳み込み。遅れを切り落として返す。
    """
    nyquist = sampling_frequency / 2.0
    lower = centre_frequency / np.sqrt(2.0)
    upper = min(centre_frequency * np.sqrt(2.0), nyquist * 0.999)
    if lower >= upper:
        raise ValueError(f"{centre_frequency:.0f} Hz バンドがナイキスト周波数 "
                         f"{nyquist:.0f} Hz に収まりません")

    if method == "butter":
        sos = butter(order // 2, [lower, upper], btype="bandpass",
                     fs=sampling_frequency, output="sos")
        return sosfilt(sos, signal)
    if method == "fir":
        taps = firwin(numtaps, [lower / nyquist, upper / nyquist],
                      window="hamming", pass_zero=False, scale=False)
        filtered = fftconvolve(signal, taps)[:len(signal) + numtaps - 1]
        delay = (numtaps - 1) // 2
        return filtered[delay:delay + len(signal)]
    raise ValueError(f"method は 'butter' か 'fir' です: {method!r}")


def _crossing_index(decay_db, level):
    """減衰曲線が level [dB] を最初に下回る添字。届かなければ -1。"""
    below = np.nonzero(decay_db <= level)[0]
    return int(below[0]) if len(below) else -1


def _decay_time(decay_db, dt, db_start, db_end):
    """2 つのレベルを横切る時刻の差から 60 dB 減衰時間を外挿する。元コード 134 行。"""
    start = _crossing_index(decay_db, db_start)
    stop = _crossing_index(decay_db, db_end)
    if start < 0 or stop < 0 or stop <= start:
        return np.nan
    return (stop - start) * dt * 60.0 / (db_start - db_end)


def decay_curves(time, ir, frequencies=None, db_max=DB_MAX, db_min=DB_MIN,
                 method="butter", order=FILTER_ORDER, numtaps=FIR_NUMTAPS,
                 verbose=True):
    """オクターブバンドごとの減衰曲線と残響時間を求める。元コード 99〜135 行。

    引数:
        time        (n,) 時間ベクトル [s]（等間隔であること）
        ir          (n,) インパルス応答
        frequencies (nf,) | None  中心周波数。None なら 63〜8k の 8 バンド
        db_max / db_min  評価区間 [dB]。既定 -5 / -35（＝T30）
                    -5 / -25 なら T20、0 / -10 なら EDT 相当

    戻り値: dict
        'frequencies' / 'time' / 'decay' (nf, n) [dB] /
        'reverberation_time' (nf,) [s] / 'curvature' (nf,) [%]

    ※ 帯域幅×減衰時間（BT 積）が小さいと推定のばらつきが大きくなる
      （例: 125 Hz で T60 = 0.3 s）。これは手法の限界であって不具合ではない
      （ISO 3382 でも BT 積が小さい条件は不確かさが増すとされている）。
    """
    if frequencies is None:
        frequencies = DEFAULT_OCTAVE_BANDS
    frequencies = np.asarray(frequencies, dtype=float)

    ir = np.asarray(ir, dtype=float)
    time = np.asarray(time, dtype=float)
    dt = float(time[1] - time[0])
    sampling_frequency = 1.0 / dt

    decay = np.empty((len(frequencies), len(ir)))
    rt = np.empty(len(frequencies))
    curvature = np.empty(len(frequencies))

    for i, fc in enumerate(frequencies):
        band = octave_bandpass(ir, fc, sampling_frequency, method=method,
                               order=order, numtaps=numtaps)
        dc = schroeder_integral(band)
        if dc[0] <= 0.0:
            decay[i] = -np.inf
            rt[i] = curvature[i] = np.nan
            continue

        with np.errstate(divide="ignore"):
            decay[i] = 10.0 * np.log10(np.maximum(dc, 0.0) / dc[0])

        rt[i] = _decay_time(decay[i], dt, db_max, db_min)

        # 曲率（ISO 3382 の C）。評価区間を半分にした推定値との食い違い。
        # 減衰がまっすぐなら 0 % に近い。反射回数不足で音が途中で切れている場合や
        # 暗騒音がある場合は大きくずれるので、結果を信用してよいかの目安になる
        half = _decay_time(decay[i], dt, db_max, (db_max + db_min) / 2.0)
        curvature[i] = (np.nan if (np.isnan(half) or half == 0.0)
                        else (rt[i] / half - 1.0) * 100.0)

    if verbose:
        label = f"T{int(abs(db_min - db_max)):d}"
        for fc, value, c in zip(frequencies, rt, curvature):
            if np.isnan(value):
                print(f"[reverberation] {fc:7.0f} Hz : {label} = 算出不可"
                      f"（減衰が {db_min:.0f} dB に届いていません）")
                continue
            if not np.isnan(c) and abs(c) > 10.0:
                note = f"  ★曲率 {c:+.0f}% — 減衰が直線でないので信用しないこと"
            else:
                note = "" if np.isnan(c) else f"  (曲率 {c:+.0f}%)"
            print(f"[reverberation] {fc:7.0f} Hz : {label} = {value:.3f} s{note}")
        finite = curvature[~np.isnan(curvature)]
        if len(finite) and np.any(np.abs(finite) > 10.0):
            print("[reverberation] ※曲率が大きいときは、最大反射回数 nref が足りず"
                  "音が途中で途切れていることが多いです。"
                  "エネルギーが 35 dB 以上減衰するまで反射させる必要があります。")

    return {"frequencies": frequencies,
            "time": dt * np.arange(len(ir)),
            "decay": decay,
            "reverberation_time": rt,
            "curvature": curvature}


def decay_measures(time, ir, frequencies=None, measures=None, method="butter",
                   order=FILTER_ORDER, numtaps=FIR_NUMTAPS, verbose=True):
    """**EDT / T20 / T30 をまとめて求める。**

    減衰曲線は 1 回だけ作り、評価区間を変えて読み取るので `decay_curves` を
    3 回呼ぶより速い。60 dB 減を厳密に見ない実務の使い方に合わせた入口。

    引数:
        measures : {名前: (開始dB, 終了dB)} | None
                   None なら EDT / T20 / T30（`DECAY_MEASURES`）

    戻り値: dict
        'frequencies' / 'time' / 'decay' (nf, n) [dB]
        'measures'  {名前: (nf,) [s]}
        'curvature' (nf,) [%]   T30 と T20 の食い違い（ISO 3382 の C）
    """
    if measures is None:
        measures = DECAY_MEASURES
    base = decay_curves(time, ir, frequencies=frequencies, method=method,
                        order=order, numtaps=numtaps, verbose=False)
    dt = float(base["time"][1] - base["time"][0])

    values = {}
    for name, (db_start, db_end) in measures.items():
        values[name] = np.array([_decay_time(d, dt, db_start, db_end)
                                 for d in base["decay"]])

    if verbose:
        names = list(values)
        print("[reverberation] " + "周波数".rjust(8)
              + "".join(f"{n:>10}" for n in names) + f"{'曲率':>10}")
        for i, fc in enumerate(base["frequencies"]):
            cells = "".join(
                ("       ---" if np.isnan(values[n][i]) else f"{values[n][i]:10.3f}")
                for n in names)
            c = base["curvature"][i]
            mark = "" if np.isnan(c) else (" ★" if abs(c) > 10.0 else "")
            curvature = "      ---" if np.isnan(c) else f"{c:+9.0f}%"
            print(f"[reverberation] {fc:7.0f}Hz{cells}{curvature}{mark}")
        finite = base["curvature"][~np.isnan(base["curvature"])]
        if len(finite) and np.any(np.abs(finite) > 10.0):
            print("[reverberation] ★ 曲率が 10% を超えたバンドは減衰が直線でないので"
                  "そのまま信用しないこと。よくある原因:")
            print("[reverberation]   ・最大反射回数 nref の不足で後部残響が切れている"
                  "（procedure.py が別途エネルギーで判定して警告する）")
            print("[reverberation]   ・音場が拡散していない。小さい室・平行面・"
                  "面ごとに吸音率が大きく違う場合は、**減衰が本当に 2 段階になる**ので"
                  "曲率が出るのが正しい（EDT と T30 の差にも表れる）")

    return {"frequencies": base["frequencies"], "time": base["time"],
            "decay": base["decay"], "measures": values,
            "curvature": base["curvature"]}


def clarity_measures(time, ir, frequencies=None, method="butter",
                     order=FILTER_ORDER, numtaps=FIR_NUMTAPS, verbose=True):
    """**明瞭度系の指標**（C50 / C80 / D50 / Ts）をバンド別に求める（TODO G-4）。

    残響指標（EDT / T20 / T30）が「どれだけ長く響くか」を見るのに対し、
    こちらは「**直接音と初期反射が、後から来る音に対してどれだけ強いか**」を見る。
    値が大きいほど音が明瞭になる。会議室・教室のように「聞き取りやすさ」が
    問われる部屋では、残響時間よりこちらが効く。

    境目の時刻は用途で使い分ける（ISO 3382-1）。
    **C50 / D50 は音声**（50 ms までを有効な音とみなす）、**C80 は音楽**（80 ms）。

    | 指標 | 定義 | 単位 |
    |---|---|---|
    | C50 | `10 log10( ∫₀^50ms p² / ∫_50ms^∞ p² )` | dB |
    | C80 | 同上、境目 80 ms | dB |
    | D50 | `∫₀^50ms p² / ∫₀^∞ p²`（Deutlichkeit・明瞭度） | 0〜1 |
    | Ts  | `∫ t·p² / ∫ p²`（重心時刻。小さいほど明瞭） | s |

    C50 と D50 は同じものの言い換えで `C50 = 10 log10(D50 / (1 - D50))` の関係にある。
    両方出しているのは、資料によってどちらで書かれているかが違うため。

    **時刻の起点は直接音の到来時刻**にする（音源から受音点までの伝搬時間ぶん、
    インパルス応答の先頭には無音がある。そこを含めると 50 ms の窓がずれる）。
    ここでは各バンドで**エネルギーが最大になる時刻**を直接音とみなす。

    戻り値: dict  'frequencies' / 'C50' / 'C80' / 'D50' / 'Ts' 各 (nf,)
    """
    time = np.asarray(time, dtype=float)
    ir = np.asarray(ir, dtype=float)
    fs = 1.0 / float(time[1] - time[0])
    if frequencies is None:
        frequencies = np.asarray(DEFAULT_OCTAVE_BANDS, dtype=float)
    frequencies = np.atleast_1d(np.asarray(frequencies, dtype=float))

    result = {name: np.full(len(frequencies), np.nan)
              for name in ("C50", "C80", "D50", "Ts")}
    for i, fc in enumerate(frequencies):
        band = octave_bandpass(ir, fc, fs, method=method, order=order, numtaps=numtaps)
        energy = band ** 2
        total = energy.sum()
        if total <= 0.0:
            continue

        # 直接音の到来時刻を起点にする
        start = int(np.argmax(energy))
        energy = energy[start:]
        t = time[start:] - time[start]
        total = energy.sum()

        for name, limit in (("C50", 0.050), ("C80", 0.080)):
            early = energy[t < limit].sum()
            late = total - early
            if early > 0.0 and late > 0.0:
                result[name][i] = 10.0 * np.log10(early / late)
        result["D50"][i] = energy[t < 0.050].sum() / total
        result["Ts"][i] = float((t * energy).sum() / total)

    result["frequencies"] = frequencies
    if verbose:
        print("[reverberation] " + "周波数".rjust(8)
              + f"{'C50[dB]':>10}{'C80[dB]':>10}{'D50':>10}{'Ts[ms]':>10}")
        for i, fc in enumerate(frequencies):
            print(f"[reverberation] {fc:7.0f}Hz{result['C50'][i]:10.2f}"
                  f"{result['C80'][i]:10.2f}{result['D50'][i]:10.3f}"
                  f"{result['Ts'][i] * 1000.0:10.1f}")
    return result


def write_clarity_measures(filename, result):
    """明瞭度系の指標を CSV に保存する。"""
    header = ["frequency_hz", "C50_db", "C80_db", "D50", "Ts_s"]
    rows = np.column_stack([result["frequencies"], result["C50"], result["C80"],
                            result["D50"], result["Ts"]])
    np.savetxt(filename, rows, delimiter=",", header=",".join(header),
               comments="", fmt="%.12g")
    return filename


def write_decay_measures(filename, result):
    """EDT / T20 / T30 を CSV に保存する。"""
    names = list(result["measures"])
    header = ["frequency_hz"] + [f"{n}_s" for n in names] + ["curvature_percent"]
    rows = np.column_stack([result["frequencies"]]
                           + [result["measures"][n] for n in names]
                           + [result["curvature"]])
    np.savetxt(filename, rows, delimiter=",", header=",".join(header),
               comments="", fmt="%.12g")
    return filename


def write_reverberation_time(filename, result):
    """残響時間を CSV に保存する。元コード 141〜145 行。

    `decay_measures()` の結果（'measures' を持つ）でも
    `decay_curves()` の結果（'reverberation_time' を持つ）でも受け付ける。
    """
    if "measures" in result:
        return write_decay_measures(filename, result)
    np.savetxt(filename,
               np.column_stack([result["frequencies"], result["reverberation_time"],
                                result["curvature"]]),
               delimiter=",",
               header="frequency_hz,reverberation_time_s,curvature_percent",
               comments="", fmt="%.12g")
    return filename


# 減衰曲線を CSV に書き出すときの時間刻み [s]。1 ms（= 1 kHz 相当）。
DECAY_CSV_INTERVAL = 0.001


def write_decay_curve(filename, result, interval=DECAY_CSV_INTERVAL):
    """減衰曲線を CSV に保存する。元コード 146〜150 行。

    **サンプリング周波数のまま書かず、`interval` ごとに間引く。**
    Schroeder 逆積分は「時刻 t 以降に残っているエネルギー」なので
    **滑らかな単調減少**で、44.1 kHz の分解能で持つ意味がない。
    そのまま書くと 3 秒 × 6 バンドで 14 MB になり、Excel で開くのも一苦労になる。
    1 ms 刻みなら 300 KB 程度で、T30 の読み取り（0.1 秒単位の議論）には十分。

    間引きは**先頭から等間隔に抜く**だけ（平均は取らない）。単調減少なので
    抜いた点も曲線の上に乗っており、傾きは変わらない。

    `interval=None` を渡せば全点書ける（細かく見たいとき）。
    残響指標そのものは間引く前のデータから求めているので、この設定に依らない。
    """
    time = result["time"]
    decay = result["decay"]
    if interval:
        dt = float(time[1] - time[0])
        step = max(1, int(round(interval / dt)))
        time, decay = time[::step], decay[:, ::step]
    header = ["time_s"] + [f"decay_{f:.0f}Hz_db" for f in result["frequencies"]]
    rows = np.column_stack([time, decay.T])
    np.savetxt(filename, rows, delimiter=",", header=",".join(header),
               comments="", fmt="%.12g")
    return filename


def reverberation_time(time, ir, rt_filename=None, decay_filename=None,
                       frequencies=None, measures=None, method="butter",
                       verbose=True):
    """残響指標（EDT / T20 / T30）の算出と保存をまとめて行う。

    measures を渡せば評価区間を変えられる（既定は `DECAY_MEASURES`）。
    """
    result = decay_measures(time, ir, frequencies=frequencies, measures=measures,
                            method=method, verbose=verbose)
    if rt_filename is not None:
        write_reverberation_time(rt_filename, result)
        if verbose:
            print(f"[reverberation] 残響時間を書き出しました: {rt_filename}")
    if decay_filename is not None:
        write_decay_curve(decay_filename, result)
        if verbose:
            print(f"[reverberation] 減衰曲線を書き出しました: {decay_filename}")
    return result
