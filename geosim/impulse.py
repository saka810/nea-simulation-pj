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


# 端数遅延（サブサンプルの位置にパルスを置く）に使う窓付き sinc のタップ数。
# ★**大きいほど厳密解（式(2)）に近づく**。実案件（パルス 8,825 本・応答 3 秒）で
#
#   タップ数   時間     実効値の比   波形の相関
#      64     0.29 s    0.99744     0.99911
#     128     0.47 s    0.99923     0.99970
#     256     0.89 s    0.99975     0.99988   ← 既定
#     512     2.62 s    0.99988     0.99996
#
# 単調に厳密解へ寄る（＝式(2) を切り詰めているだけだと確かめられる）。
# 256 でも厳密版（17.4 秒）の 20 倍速いので、**余裕を見て 256 を既定にした**
# （2026-08-23。差は 8 kHz バンドより上に集まる。`docs/技術説明書.md` 8.2.1）
FRACTIONAL_TAPS = 256


def impulse_train(energy32_air, time, nfft, sound_velocity,
                  sampling_frequency=SAMPLING_FREQUENCY,
                  taps=FRACTIONAL_TAPS):
    """パルス列を**時間領域の列**にする（1/3 オクターブバンドごと）。戻り値 (32, nfft)。

    ★**これが高速化の要**（2026-08-21。ユーザー要望「配列演算で高速化」）。

    `transfer_function` は周波数ごとに全パルスの位相を足すので
    **計算量が (周波数の数) × (パルスの数)** になる。実案件（研修室・受音点 1 点・
    パルス 8,825 本・応答 3 秒）で **15.7 秒**かかっていた（全体の 46%）。
    いっぽう時間領域に置いてから FFT すれば **(パルスの数) + (FFT)** で済み、
    同じ条件で **0.25 秒（63 倍）**になる。

    ただし到来時刻はサンプルの格子に乗らないので、**端数の遅延**を
    窓付き sinc（Hann 窓）で分配する。単純な最近傍（1 タップ）や
    線形補間（2 タップ）だと高域が落ちて実効値が 17% も減るが、
    64 タップなら 0.5% に収まる（`tests/test_geosim.py` の [39] で押さえてある）。

    厳密版（`transfer_function`）は**参照実装として残してある**。
    """
    energy32_air = np.asarray(energy32_air, dtype=float)
    time = np.asarray(time, dtype=float)
    amplitude = np.sqrt(energy32_air) / (time * sound_velocity)[None, :]

    half = int(taps) // 2
    # ★fs は**引数から**取る（モジュール定数を見ていて、fs を変えると
    #   パルスの置き場が全部ずれるバグがあった。2026-08-23 に気づいた）
    position = time * sampling_frequency
    index = np.floor(position).astype(np.int64)
    fraction = position - index

    offsets = np.arange(-half + 1, half + 1)

    # ★置き場は **nfft で折り返す**（2026-08-23 に直した）。
    #   窓付き sinc は前後に裾を持つので、始まりに近いパルスでは裾が 0 より前へ、
    #   終わりに近いパルスでは nfft より後へはみ出す。
    #   以前は**はみ出すパルスを丸ごと捨てていて、タップ数を増やすと
    #   直接音が消えていた**（taps=512 で実効値が 22% 落ちた）。
    #   厳密版は DFT なので周期的に折り返る。同じ扱いに揃えるのが正しい
    #   （はみ出した先は FIR の遅れを取り除くときに出力の外へ出る）。
    #
    # ★足し込みは **`np.bincount`**（2026-08-23。`np.add.at` から替えた）。
    #   結果は完全に同じで **256 タップで 1.9 倍・1024 タップで 3.7 倍**速い。
    #   `np.add.at` は要素ごとに回るので遅い。これでタップ数を増やす負担が減り、
    #   **厳密解に寄せるか速さを取るかで悩まなくなる**（ユーザー判断 2026-08-23:
    #   「複雑な部屋では微小な誤差より計算時間を優先したい」）。
    #   パルスを塊に分けているのはメモリ対策（パルス数 × タップ数の配列を作るため）。
    bands = energy32_air.shape[0]
    trains = np.zeros((bands, nfft))
    chunk = max(1, int(2_000_000 // max(1, len(offsets))))
    for start in range(0, len(index), chunk):
        stop = min(start + chunk, len(index))
        # 窓付き sinc。合計 1 に正規化して直流の重みを保つ
        x = fraction[start:stop, None] - offsets[None, :]
        kernel = np.sinc(x) * (0.5 + 0.5 * np.cos(np.pi * x / half))
        kernel /= kernel.sum(axis=1, keepdims=True)
        position = ((index[start:stop, None] + offsets[None, :]) % nfft).ravel()
        for band in range(bands):
            trains[band] += np.bincount(
                position,
                weights=(amplitude[band, start:stop, None] * kernel).ravel(),
                minlength=nfft)
    return trains


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
                     numtaps=NUMTAPS, verbose=True, method="fast",
                     taps=FRACTIONAL_TAPS):
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
        print(f"[impulse] 合成のやり方: "
              f"{'厳密（式2.67）' if method == 'exact' else '高速（時間領域→FFT）'}")
        print(f"[impulse] fs={sampling_frequency:.0f} Hz / 出力 {n_out} 点 "
              f"({max_time:.3f} s) / FFT長 {nfft} / パルス {len(time)} 本 / "
              f"バンド {len(octave_frequencies)}（{octave_frequencies[0]:.0f}〜"
              f"{octave_frequencies[-1]:.0f} Hz）")

    # ---- 1/3 オクターブへの展開と空気吸収 ----
    energy32 = expand_to_third_octave(energy, octave_frequencies, mf)
    energy32 = apply_air_absorption(energy32, time, atmosphere, mf)

    # ---- バンドごとの成分を周波数領域で作る ----
    #
    # `method='fast'`（既定）… 時間領域に置いてから FFT（`impulse_train`）。
    #   実案件で **63 倍**速い（2026-08-21）。端数遅延は窓付き sinc で分配する
    # `method='exact'`      … 式(2.67) をそのまま（`transfer_function`）。
    #   **参照実装**。周波数 × パルスの総当たりなので遅い
    lower, upper = bandpass_edges(mf, fmax)
    delay = (numtaps - 1) // 2
    if method == "exact":
        spectrum = transfer_function(energy32, time, nfreq, df, sound_velocity)
    else:
        spectrum = rfft(impulse_train(energy32, time, nfft, sound_velocity,
                                      sampling_frequency, taps),
                        n=nfft, axis=1).T

    # ★**足し合わせは周波数領域で行い、逆変換は 1 回だけ**にする（線形なので同じ）。
    #   以前はバンドごとに逆変換していて 32 回呼んでいた
    total = np.zeros(nfreq, dtype=np.complex128)
    for j in range(len(mf)):
        total += spectrum[:, j] * rfft(filter_bandpass(numtaps, lower[j],
                                                        upper[j]), nfft)
    ir = irfft(total, nfft)

    # 直線位相 FIR の遅れを取り除いて 0 秒始まりにする
    ir = ir[delay:delay + n_out]
    t = np.arange(n_out) / sampling_frequency
    return t, ir


def write_impulse_response(filename, t, ir):
    """インパルス応答を CSV に保存する。元コード 232〜235 行。"""
    np.savetxt(filename, np.column_stack([t, ir]), delimiter=",",
               header="time_s,ir", comments="", fmt="%.12g")
    return filename


def impulse_response_from_pulses(filename, pulses, octave_frequencies=None, atmosphere=None,
                     sampling_frequency=SAMPLING_FREQUENCY, max_time=MAX_TIME,
                     numtaps=NUMTAPS, verbose=True, method="fast",
                     taps=FRACTIONAL_TAPS):
    """パルス列 → インパルス応答 → CSV 保存までの一括処理。

    pulses : loop_noredundancy.PulseList（`.time` と `.energy` を使う）
    """
    t, ir = impulse_response(pulses.time, pulses.energy,
                             octave_frequencies=octave_frequencies,
                             atmosphere=atmosphere,
                             sampling_frequency=sampling_frequency,
                             max_time=max_time, numtaps=numtaps, verbose=verbose,
                             method=method, taps=taps)
    write_impulse_response(filename, t, ir)
    if verbose:
        print(f"[impulse] インパルス応答を書き出しました: {filename}")
    return t, ir
