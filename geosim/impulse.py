"""パルス列 → インパルス応答の合成。元コード make_ipls_freq_monaural_fortran.f90。

## 元コードからの変更方針（2026-08-14）

元コードは Fortran で使えるライブラリが限られていたため、FFT もバンドパスフィルタも
自前で書いていた。Python には scipy があるので、**信号処理は既存ライブラリに任せる**。

| 処理 | 元コード | 本実装 |
|---|---|---|
| FFT | 自作 Cooley-Tukey | `scipy.fft`（実数信号なので `rfft`/`irfft`） |
| バンドパス FIR | 自作 `fir1_bandpass` | `scipy.signal.firwin` |
| 負の周波数への拡張 | 手動でエルミート対称に詰める | **不要**（`irfft` が実数信号を前提に処理する） |
| 畳み込み | 巡回畳み込み（長さ nn） | ゼロ詰めして**線形畳み込み** |
| フィルタ長 | データ点数 nn = 131072（遅れ 1.49 秒） | 既定 8192（遅れ 93 ms） |

結果として、**元コードにあった 1.49 秒の遅れは無くなった**。
出力されるインパルス応答は 0 秒が直接音より前（＝音の始まり）になる。

## 流れ

    1. オクターブバンドのエネルギーを 1/3 オクターブ 32 バンドに割り当てる
    2. 空気吸収を掛ける（温度・湿度から。`atmosphere.py`）
    3. パルス列を周波数領域の伝達関数にする（式2.67）
    4. バンドパスを掛けて足し合わせ、逆変換する

## 周波数バンド

エネルギーのバンド数は**可変**。既定は 8 バンド（63〜8k Hz）で、
63 Hz と 8 kHz を対象外にする場合は 6 バンド（125〜4k Hz）でもよい。
合成そのものは 1/3 オクターブ 32 バンド（15.625 Hz〜20 kHz）で行う。
オクターブ → 1/3 オクターブの割り当ては**中心周波数が対数軸で最も近いもの**を選ぶ。
この規則は元コードの手書き対応表（94〜125 行）を**完全に再現する**
（tests/test_geosim.py で確認済み）。
"""

import numpy as np
from scipy.fft import irfft, next_fast_len, rfft
from scipy.signal import firwin

from absorption import DEFAULT_OCTAVE_BANDS, octave_bands
from atmosphere import Atmosphere

# 元コード 9〜10 行
SAMPLING_FREQUENCY = 44100.0
MAX_TIME = 1.0

# 合成に使うバンドパスのタップ数。
# 元コードはデータ点数 nn（=131072）を使っていたが、遅れが 1.49 秒になるうえ
# 巡回畳み込みの回り込みを生むので、実用的な長さにした。
# 8192 タップだと遷移幅がおよそ 18 Hz で、オクターブ境界（最低 44 Hz）を切るのに十分。
NUMTAPS = 8192

# 後方互換。バンド数は絶対値ではなく `absorption.octave_bands()` で決める
OCTAVE_BAND_FREQUENCIES = DEFAULT_OCTAVE_BANDS


def third_octave_bands(count=32, start=15.625):
    """1/3 オクターブバンドの中心周波数 [Hz]。元コード 46 行。

    mf(i) = 15.625 * 2^((i-1)/3) → 15.625 Hz 〜 20 kHz
    """
    return start * 2.0 ** (np.arange(count) / 3.0)


def band_mapping(octave_frequencies, third_frequencies):
    """1/3 オクターブの各バンドに、どのオクターブバンドの値を使うかを返す。

    **対数軸で中心周波数が最も近いオクターブバンド**を選ぶ。

    元コードは 6 バンド用の対応表を手書きしていた（94〜125 行）が、
    バンド数を変えられるようにするためこの規則で置き換えた。
    6 バンド（125〜4k）の場合、この規則は元コードの表と**完全に一致する**。
    """
    octave = np.log2(np.asarray(octave_frequencies, dtype=float))
    third = np.log2(np.asarray(third_frequencies, dtype=float))
    return np.argmin(np.abs(third[:, None] - octave[None, :]), axis=1)


def expand_to_third_octave(energy, octave_frequencies=None, third_frequencies=None):
    """オクターブバンドのエネルギーを 1/3 オクターブに展開する。元コード 94〜125 行。

    引数:  energy (count, nband)  パルス列のバンド別エネルギー
    戻り値:        (32, count)    バンドが行（元コード enerred_air と同じ並び）
    """
    energy = np.asarray(energy, dtype=float)
    if energy.ndim != 2:
        raise ValueError(f"エネルギーは (count, バンド数) の 2 次元配列です: {energy.shape}")
    if octave_frequencies is None:
        octave_frequencies = octave_bands(energy.shape[1])
    if third_frequencies is None:
        third_frequencies = third_octave_bands()
    if len(octave_frequencies) != energy.shape[1]:
        raise ValueError(f"バンド数が合いません: エネルギー {energy.shape[1]} / "
                         f"中心周波数 {len(octave_frequencies)}")
    return energy[:, band_mapping(octave_frequencies, third_frequencies)].T


def apply_air_absorption(energy32, time, atmosphere=None, third_frequencies=None):
    """空気吸収による減衰を掛ける。元コード 127 行。

        E' = E * exp(-m * c * t)      c*t = 伝搬距離

    減衰係数 m は **ISO 9613-1**（`atmosphere.py`）で温度・湿度・気圧から求める。
    元コードは `1.81e-8 * f^1.57`（20℃・湿度 40% のフィット）で固定だった。
    低音側で元の近似は ISO の 1/3 程度しかなく、高音側は概ね一致する。

    引数:  energy32 (32, count)、time (count,)
    戻り値:         (32, count)
    """
    if atmosphere is None:
        atmosphere = Atmosphere()
    if third_frequencies is None:
        third_frequencies = third_octave_bands()
    m = atmosphere.absorption_coefficient(third_frequencies)
    distance = atmosphere.sound_velocity * np.asarray(time, dtype=float)
    return energy32 * np.exp(-m[:, None] * distance[None, :])


def transfer_function(energy32_air, time, nfreq, df, sound_velocity, chunk=4096):
    """パルス列を周波数領域の伝達関数に変換する。元コード 136〜143 行、書籍 式(2.67)。

        H_b(f) = Σ_n sqrt(E_{n,b}) / (c t_n) * exp(-i 2π f t_n)

    ・sqrt(エネルギー) → 振幅
    ・1/(c t_n) = 1/距離 → 球面波の距離減衰
    ・exp(-i 2π f t) → 到来時刻ぶんの位相回転（時間シフト定理）

    バンドによらない `exp` の部分とバンドごとの振幅を分けて行列積にしている
    （素直に 3 重ループすると exp を 32 回計算し直すことになる）。
    周波数方向に分割しているのはメモリ対策で、結果は変わらない。

    戻り値: (nfreq, 32) complex
    """
    time = np.asarray(time, dtype=float)
    amplitude = np.sqrt(energy32_air) / (time * sound_velocity)[None, :]   # (32, count)

    result = np.empty((nfreq, energy32_air.shape[0]), dtype=np.complex128)
    for start in range(0, nfreq, chunk):
        stop = min(start + chunk, nfreq)
        index = np.arange(start, stop, dtype=float)
        phase = np.exp(-2j * np.pi * df * index[:, None] * time[None, :])
        result[start:stop, :] = phase @ amplitude.T
    return result


def bandpass_edges(mf, fmax):
    """各 1/3 オクターブバンドの正規化遮断周波数。元コード 152〜154 行。

    下端 mf * 2^(-1/6)、上端 mf * 2^(1/6) をナイキスト周波数で割る。
    ナイキストを超えるバンドは 0.999 で頭打ちにする。
    """
    lower = mf * 2.0 ** (-1.0 / 6.0) / fmax
    upper = np.minimum(mf * 2.0 ** (1.0 / 6.0) / fmax, 0.999)
    return lower, upper


def filter_bandpass(numtaps, wmin, wmax):
    """ハミング窓付き FIR バンドパス。`scipy.signal.firwin` を使う。

    `scale=False` が要点。scipy は既定（`scale=True`）で通過域のゲインを 1 に
    正規化するが、それをすると**バンドを足し合わせたときに平坦にならない**。
    元コードの `fir1_bandpass` は正規化していないので `scale=False` が対応する。
    実際、両者は相対誤差 1e-16 で一致する（tests/test_geosim.py で確認）。
    """
    return firwin(numtaps, [wmin, wmax], window="hamming",
                  pass_zero=False, scale=False)


def _fir1_bandpass_fortran(numtaps, wmin, wmax):
    """元コード `fir1_bandpass`（241〜269 行）をそのまま移植したもの。

    **通常は使わない。** `filter_bandpass`（scipy 版）と一致することを
    テストで示すためだけに残している。
    """
    n = numtaps - 1
    m = n / 2.0
    shift = np.arange(numtaps, dtype=float) - m
    filter_bp = np.empty(numtaps)
    zero = (shift == 0.0)
    filter_bp[zero] = wmax - wmin
    s = shift[~zero]
    filter_bp[~zero] = (np.sin(np.pi * wmax * s) / (np.pi * s)
                        - np.sin(np.pi * wmin * s) / (np.pi * s))
    return filter_bp * (0.54 + 0.46 * np.cos(np.pi * 2.0 * shift / n))


def impulse_response(time, energy, octave_frequencies=None, atmosphere=None,
                     sampling_frequency=SAMPLING_FREQUENCY, max_time=MAX_TIME,
                     numtaps=NUMTAPS, verbose=True):
    """パルス列からインパルス応答を合成する。

    引数:
        time      (count,)        到来時刻 [s]
        energy    (count, nband)  オクターブバンド別エネルギー
        octave_frequencies (nband,) | None
            バンドの中心周波数。None ならバンド数から決める（8→63〜8k、6→125〜4k）
        atmosphere : Atmosphere | None
            温度・湿度・気圧。音速と空気吸収の両方をここから取る
        max_time  出力するインパルス応答の長さ [s]
        numtaps   バンドパスのタップ数。大きいほど帯域分割は鋭いが遅れが増える

    戻り値:
        (t, ir)  ともに (max_time * fs,)。t は 0 始まりの時間ベクトル [s]
    """
    time = np.asarray(time, dtype=float)
    energy = np.asarray(energy, dtype=float)
    if len(time) == 0:
        raise ValueError("パルス列が空です（受音に至った経路がありません）")
    if atmosphere is None:
        atmosphere = Atmosphere()
    if octave_frequencies is None:
        octave_frequencies = octave_bands(energy.shape[1])

    sound_velocity = atmosphere.sound_velocity
    n_out = int(round(max_time * sampling_frequency))

    # 出力より後に届くパルスは捨てる（そのぶん残響が切れることを伝える）
    inside = time < max_time
    if not np.all(inside):
        if verbose:
            print(f"[impulse] 注意: 到来時刻が max_time={max_time} s を超えるパルス "
                  f"{int((~inside).sum())} 本を捨てました"
                  f"（最大 {time.max():.3f} s）。max_time を伸ばしてください")
        time, energy = time[inside], energy[inside]
        if len(time) == 0:
            raise ValueError("max_time 内に届くパルスがありません")

    # ゼロ詰めして線形畳み込みにする（巡回畳み込みの回り込みを避ける）
    nfft = next_fast_len(n_out + numtaps)
    nfreq = nfft // 2 + 1
    df = sampling_frequency / nfft
    fmax = sampling_frequency / 2.0

    mf = third_octave_bands()

    if verbose:
        print(f"[impulse] {atmosphere.summary()}")
        print(f"[impulse] fs={sampling_frequency:.0f} Hz / 出力 {n_out} 点 "
              f"({max_time:.3f} s) / FFT長 {nfft} / パルス {len(time)} 本 / "
              f"バンド {len(octave_frequencies)}（{octave_frequencies[0]:.0f}〜"
              f"{octave_frequencies[-1]:.0f} Hz）")

    # ---- 1/3 オクターブへの展開と空気吸収 ----
    energy32 = expand_to_third_octave(energy, octave_frequencies, mf)
    energy32 = apply_air_absorption(energy32, time, atmosphere, mf)

    # ---- 伝達関数（式2.67）----
    spectrum = transfer_function(energy32, time, nfreq, df, sound_velocity)

    # ---- バンドごとにフィルタを掛けて足し合わせ、時間領域に戻す ----
    lower, upper = bandpass_edges(mf, fmax)
    delay = (numtaps - 1) // 2
    ir = np.zeros(nfft)
    for j in range(len(mf)):
        # 周波数領域の積 = 時間領域の畳み込み（ゼロ詰め済みなので線形畳み込み）
        band_spectrum = spectrum[:, j] * rfft(filter_bandpass(numtaps, lower[j],
                                                              upper[j]), nfft)
        ir += irfft(band_spectrum, nfft)

    # 直線位相 FIR の遅れを取り除いて 0 秒始まりにする
    ir = ir[delay:delay + n_out]
    t = np.arange(n_out) / sampling_frequency
    return t, ir


def write_impulseresponce(filename, t, ir):
    """インパルス応答を CSV に保存する。元コード 232〜235 行。"""
    np.savetxt(filename, np.column_stack([t, ir]), delimiter=",",
               header="time_s,ir", comments="", fmt="%.12g")
    return filename


def impulse_responce(filename, pulses, octave_frequencies=None, atmosphere=None,
                     sampling_frequency=SAMPLING_FREQUENCY, max_time=MAX_TIME,
                     numtaps=NUMTAPS, verbose=True):
    """パルス列 → インパルス応答 → CSV 保存までの一括処理。

    pulses : loop_noredundancy.PulseList（`.time` と `.energy` を使う）
    """
    t, ir = impulse_response(pulses.time, pulses.energy,
                             octave_frequencies=octave_frequencies,
                             atmosphere=atmosphere,
                             sampling_frequency=sampling_frequency,
                             max_time=max_time, numtaps=numtaps, verbose=verbose)
    write_impulseresponce(filename, t, ir)
    if verbose:
        print(f"[impulse] インパルス応答を書き出しました: {filename}")
    return t, ir
