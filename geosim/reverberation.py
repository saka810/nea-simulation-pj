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

from absorption import DEFAULT_OCTAVE_BANDS

# 元コード 9〜10 行。T30 を求めるための評価区間
DB_MAX = -5.0
DB_MIN = -35.0

# IEC 61260 のクラス 1 オクターブフィルタに相当する次数（Butterworth 6 次）
FILTER_ORDER = 6

# method='fir' のときのタップ数
FIR_NUMTAPS = 4096


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


def write_reverberation_time(filename, result):
    """残響時間を CSV に保存する。元コード 141〜145 行。"""
    np.savetxt(filename,
               np.column_stack([result["frequencies"], result["reverberation_time"],
                                result["curvature"]]),
               delimiter=",",
               header="frequency_hz,reverberation_time_s,curvature_percent",
               comments="", fmt="%.12g")
    return filename


def write_decay_curve(filename, result):
    """減衰曲線を CSV に保存する。元コード 146〜150 行。"""
    header = ["time_s"] + [f"decay_{f:.0f}Hz_db" for f in result["frequencies"]]
    rows = np.column_stack([result["time"], result["decay"].T])
    np.savetxt(filename, rows, delimiter=",", header=",".join(header),
               comments="", fmt="%.12g")
    return filename


def reverberation_time(time, ir, rt_filename=None, decay_filename=None,
                       frequencies=None, db_max=DB_MAX, db_min=DB_MIN,
                       method="butter", verbose=True):
    """残響時間の算出と保存をまとめて行う。"""
    result = decay_curves(time, ir, frequencies=frequencies, db_max=db_max,
                          db_min=db_min, method=method, verbose=verbose)
    if rt_filename is not None:
        write_reverberation_time(rt_filename, result)
        if verbose:
            print(f"[reverberation] 残響時間を書き出しました: {rt_filename}")
    if decay_filename is not None:
        write_decay_curve(decay_filename, result)
        if verbose:
            print(f"[reverberation] 減衰曲線を書き出しました: {decay_filename}")
    return result
