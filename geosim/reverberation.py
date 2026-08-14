"""インパルス応答 → 残響時間。元コード ipls2rt_fortran.f90。

流れ（元コードの節と対応）:

    1. インパルス応答をオクターブバンドに分ける（元 87, 101〜105 行）
    2. フィルタの直線位相遅れを巡回シフトで打ち消す（元 107〜113 行）
    3. Schroeder 積分（後ろから積分した残留エネルギー）を取る（元 114〜117 行）
    4. dB 表示にして -5 dB と -35 dB を横切る時刻の差から T30 を求める（元 119〜134 行）

【元コードの nn/2 シフトについて】
    元コードは 108〜112 行で前半と後半を入れ替えている。これは
    タップ数 = nn の直線位相 FIR が持つ (nn-1)/2 サンプルの遅れを打ち消すため。
    つまり元コードも「バンドパスの遅れは邪魔だ」と分かっていて、
    残響時間の段では消している（インパルス応答の段では消していない）。
    impulse.py の `compensate_filter_delay` はこれと同じ発想。

【元コードとの相違点：巡回畳み込み → 線形畳み込み】★重要
    元コードは FFT で掛け算しているので**巡回畳み込み**になる。
    バンドパスは左右対称な直線位相 FIR なので、入力の立ち上がりに対する応答が
    **入力より前の時刻から始まる**（プリリンギング）。巡回畳み込みではその「負の時刻」ぶんが
    **バッファの末尾に回り込む**。しかも元コードはタップ数 = nn なので、
    回り込む量はバッファの半分に及ぶ。

    結果、減衰曲線の後半がまるごとプリリンギングで埋まり、
    Schroeder 積分に -20 dB 程度の**床**ができて -35 dB まで落ちなくなる。
    合成インパルス応答（減衰が既知の指数減衰雑音）で確かめたところ、
    元コードのやり方では T30 が理論値の 5〜19 倍になるか、そもそも算出できなかった。

    そこで本実装は**線形畳み込み**（ゼロ詰めして畳み込み、フィルタの遅れぶんを切り落とす）
    に変えている。同じ検証で誤差 0.4〜3.8 % に収まる。
    あわせてタップ数も既定 4096 にした（nn = 131072 だと遅れが 1.5 秒になり、
    切り落とす量が信号より長くなってしまうため）。

【残響時間が index の差で出る理由】
    T30 = (dbmin を切る時刻 - dbmax を切る時刻) * 60 / (dbmax - dbmin)
    は**時刻の差**なので、インパルス応答全体が一定量ずれていても影響を受けない。
    元コードのインパルス応答が 1.49 秒遅れたままでも残響時間だけは正しく出るのはこのため。
    ただし減衰曲線をグラフにするときはずれが見えるので、
    impulse.py 側で補正しておくほうが素直。
"""

import numpy as np

from impulse import OCTAVE_BAND_FREQUENCIES, filter_bandpass

# 元コード 9〜10 行。T30 を求めるための評価区間
DB_MAX = -5.0
DB_MIN = -35.0

# バンドパスのタップ数。元コードはデータ点数 nn（＝131072）を使っていたが、
# 遅れが 1.5 秒にもなって使いものにならないので既定を 4096 にした（上の docstring 参照）。
# 4096 タップだと 125 Hz バンド（幅 88 Hz）の遷移幅がおよそ 20 Hz で、
# オクターブバンドを切り出すのに十分。
NUMTAPS = 4096


def schroeder_integral(x):
    """後ろ向きに積分した残留エネルギー。元コード 114〜117 行。

        dc(j) = Σ_{m >= j} x(m)^2

    「時刻 j 以降に残っているエネルギー」を表す。これを dB にしたものが減衰曲線。
    ランダムな残響音をそのまま dB 表示すると激しく暴れるが、積分すると滑らかになる
    （Schroeder の逆積分法）。
    """
    return np.cumsum(x[::-1] ** 2)[::-1]


def bandpass(signal, numtaps, lower, upper):
    """線形畳み込みでバンドパスをかけ、フィルタの遅れを取り除いて返す。

    FFT で掛け算するのは元コードと同じだが、**ゼロ詰めしてから**掛けることで
    巡回畳み込みではなく線形畳み込みにしている（回り込みを起こさないため）。
    そのうえで直線位相 FIR の遅れ (numtaps-1)//2 サンプルを切り落とし、
    入力と同じ長さ・同じ時刻軸で返す。
    """
    h = filter_bandpass(numtaps, lower, upper)
    n_out = len(signal) + numtaps - 1
    nfft = 1 << (n_out - 1).bit_length()
    filtered = np.real(np.fft.ifft(np.fft.fft(signal, nfft) * np.fft.fft(h, nfft)))
    delay = (numtaps - 1) // 2
    return filtered[delay:delay + len(signal)]


def _crossing_index(decay_db, level):
    """減衰曲線が level [dB] を最初に下回る添字。届かなければ -1。"""
    below = np.nonzero(decay_db <= level)[0]
    return int(below[0]) if len(below) else -1


def _decay_time(decay_db, dt, db_start, db_end):
    """減衰曲線の 2 つのレベルを横切る時刻の差から 60 dB 減衰時間を外挿する。

    元コード 134 行と同じ式。db_start > db_end（例: -5 と -35）。
    """
    start = _crossing_index(decay_db, db_start)
    stop = _crossing_index(decay_db, db_end)
    if start < 0 or stop < 0 or stop <= start:
        return np.nan
    return (stop - start) * dt * 60.0 / (db_start - db_end)


def decay_curves(time, ir, frequencies=None, db_max=DB_MAX, db_min=DB_MIN,
                 numtaps=NUMTAPS, verbose=True):
    """オクターブバンドごとの減衰曲線と残響時間を求める。元コード 99〜135 行。

    引数:
        time        (n,) 時間ベクトル [s]（等間隔であること）
        ir          (n,) インパルス応答
        frequencies (nf,) | None  オクターブバンド中心周波数。
                    None なら 125〜4000 Hz の 6 バンド（元コード 72 行と同じ）
        db_max / db_min  評価区間 [dB]。既定 -5 / -35（＝T30）
                    -5 / -25 にすれば T20、0 / -10 にすれば EDT 相当になる
        numtaps     バンドパスのタップ数。既定 4096

    戻り値: dict
        'frequencies' (nf,)    中心周波数 [Hz]
        'time'        (n,)     時間ベクトル [s]
        'decay'       (nf, n)  減衰曲線 [dB]
        'reverberation_time' (nf,) 残響時間 [s]（求まらなければ NaN）
        'curvature'   (nf,)    減衰の曲率 [%]（ISO 3382 の C）。
                      0 に近いほど減衰が直線。10 % を超えたら結果を疑うこと

    ※ 低い周波数で減衰が短いと（例: 125 Hz で T60 = 0.3 s）、帯域幅×減衰時間が
      小さくなって推定のばらつきが大きくなる。これは手法の限界であって不具合ではない
      （ISO 3382 でも BT 積が小さい条件は不確かさが増すとされている）。
    """
    if frequencies is None:
        frequencies = OCTAVE_BAND_FREQUENCIES
    frequencies = np.asarray(frequencies, dtype=float)

    ir = np.asarray(ir, dtype=float)
    time = np.asarray(time, dtype=float)
    dt = float(time[1] - time[0])
    fmax = 1.0 / dt / 2.0

    if numtaps >= len(ir):
        raise ValueError(f"タップ数 {numtaps} がインパルス応答の長さ {len(ir)} 以上です")

    decay = np.empty((len(frequencies), len(ir)))
    rt = np.empty(len(frequencies))
    curvature = np.empty(len(frequencies))

    for i, fc in enumerate(frequencies):
        # オクターブバンド（中心周波数の 1/√2 〜 √2 倍）。元コード 87 行
        lower = fc / np.sqrt(2.0) / fmax
        upper = fc * np.sqrt(2.0) / fmax
        if upper >= 1.0:
            raise ValueError(
                f"{fc:.0f} Hz バンドの上端 {upper * fmax:.0f} Hz が"
                f"ナイキスト周波数 {fmax:.0f} Hz を超えます")

        band = bandpass(ir, numtaps, lower, upper)

        dc = schroeder_integral(band)
        if dc[0] <= 0.0:
            decay[i] = -np.inf
            rt[i] = np.nan
            curvature[i] = np.nan
            continue

        with np.errstate(divide="ignore"):
            decay[i] = 10.0 * np.log10(np.maximum(dc, 0.0) / dc[0])

        # 元コード 134 行。(db_max - db_min) は正（-5 - (-35) = 30）
        rt[i] = _decay_time(decay[i], dt, db_max, db_min)

        # 曲率（ISO 3382 の C）。評価区間を半分にした推定値とどれだけ食い違うかを見る。
        # 減衰がまっすぐなら 0 % に近い。反射回数が足りずに音が途中で切れている場合や、
        # 暗騒音・作り物の床がある場合は大きくずれるので、結果を信用してよいかの目安になる。
        half = _decay_time(decay[i], dt, db_max, (db_max + db_min) / 2.0)
        curvature[i] = np.nan if (np.isnan(half) or half == 0.0) else (rt[i] / half - 1.0) * 100.0

    if verbose:
        label = f"T{int(abs(db_min - db_max)):d}"
        for fc, value, c in zip(frequencies, rt, curvature):
            if np.isnan(value):
                print(f"[reverberation] {fc:7.0f} Hz : {label} = 算出不可"
                      f"（減衰が {db_min:.0f} dB に届いていません）")
                continue
            note = ""
            if not np.isnan(c) and abs(c) > 10.0:
                note = f"  ★曲率 {c:+.0f}% — 減衰が直線でないので信用しないこと"
            elif not np.isnan(c):
                note = f"  (曲率 {c:+.0f}%)"
            print(f"[reverberation] {fc:7.0f} Hz : {label} = {value:.3f} s{note}")
        if np.any(np.abs(curvature[~np.isnan(curvature)]) > 10.0):
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
               np.column_stack([result["frequencies"], result["reverberation_time"]]),
               delimiter=",", header="frequency_hz,reverberation_time_s",
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
                       numtaps=NUMTAPS, verbose=True):
    """残響時間の算出と保存をまとめて行う。"""
    result = decay_curves(time, ir, frequencies=frequencies, db_max=db_max,
                          db_min=db_min, numtaps=numtaps, verbose=verbose)
    if rt_filename is not None:
        write_reverberation_time(rt_filename, result)
        if verbose:
            print(f"[reverberation] 残響時間を書き出しました: {rt_filename}")
    if decay_filename is not None:
        write_decay_curve(decay_filename, result)
        if verbose:
            print(f"[reverberation] 減衰曲線を書き出しました: {decay_filename}")
    return result
