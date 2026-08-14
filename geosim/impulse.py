"""パルス列 → インパルス応答の合成。元コード make_ipls_freq_monaural_fortran.f90。

流れ（元コードの節と対応）:

    1. 1/3 オクターブ 32 バンドの中心周波数と空気吸収係数を作る（元 45〜48 行）
    2. バックトレースの **6 バンド**エネルギーを **32 バンド**に展開する（元 94〜125 行）
    3. 空気吸収を掛ける（元 127 行）
    4. パルス列を伝達関数に変換する（元 136〜143 行）
    5. バンドパスフィルタ（FIR）を作って FFT する（元 151〜165 行）
    6. 正の周波数だけの伝達関数を負の周波数へ拡張する（元 173〜180 行）
    7. 掛け合わせて（＝畳み込み）逆 FFT し、32 バンドを足し合わせる（元 188〜216 行）

【周波数バンドについて】
    バックトレースが持つのは **6 オクターブバンド**。元コードの展開表（94〜125 行）から
    その中心周波数は **125 / 250 / 500 / 1000 / 2000 / 4000 Hz** と確定できる
    （例: 32 バンド側の 10 番目 = 125Hz が enerred(1) を参照している）。
    合成は 1/3 オクターブ 32 バンド（15.625 Hz〜20 kHz）で行う。

【フィルタの遅れに注意】
    元コードのバンドパスは **タップ数 = データ点数 nn** の直線位相 FIR なので、
    出力は (nn-1)/2 サンプル ≒ **1.49 秒**（fs=44100, nn=131072 のとき）遅れる。
    元コードはこれを補正していないので、出てくる CSV は
    「0 秒付近は無音、1.49 秒あたりから音が始まる」形になる。
    `compensate_filter_delay=True` を渡すと補正して 0 秒始まりにできる。
    残響時間や C50 を出す（G-4）ときは補正しないと使えないので注意。
"""

import numpy as np

# 元コード 9〜10 行
SAMPLING_FREQUENCY = 44100.0
MAX_TIME = 1.0

# 元コード 42 行。loop_noredundancy.SOUND_VELOCITY と揃えること
SOUND_VELOCITY = 340.0

# 6 オクターブバンド → 1/3 オクターブ 32 バンドの対応表。元コード 94〜125 行。
# 値は 6 バンド側の添字（0 始まり）。
#   32 バンド側:  1     2     3     4     5     6     7     8    (16〜80 Hz)
#                 9    10    11    12    13    14    15    16    (100〜500 Hz)
#                ...
BAND_EXPANSION_6_TO_32 = np.array([
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   # 16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160 Hz
    1, 1, 1,                           # 200, 250, 315 Hz
    2, 2, 2,                           # 400, 500, 630 Hz
    3, 3, 3,                           # 800, 1000, 1250 Hz
    4, 4, 4,                           # 1600, 2000, 2500 Hz
    5, 5, 5, 5, 5, 5, 5, 5, 5,         # 3150, 4000, 5000, 6300, 8000, 10k, 12.5k, 16k, 20k Hz
], dtype=int)

# 6 オクターブバンドの中心周波数 [Hz]。上の展開表から逆算した値。
OCTAVE_BAND_FREQUENCIES = np.array([125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0])


def third_octave_bands():
    """1/3 オクターブ 32 バンドの中心周波数 [Hz]。元コード 46 行。

    mf(i) = 15.625 * 2^((i-1)/3) → 15.625 Hz 〜 20 kHz
    """
    return 15.625 * 2.0 ** (np.arange(32) / 3.0)


def airdamping_coefficient(mf=None):
    """空気吸収の減衰係数 m [1/m]。元コード 47 行。

    mair(i) = 1.81e-8 * mf(i)^1.57
    20℃・湿度 40% の読み取り値をべき乗近似したもの（元コード 28 行のコメント）。

    2026-08-14 修正: 旧実装は `np.pow(15.625 * 2, 1/3 * i)` と括弧が崩れており、
    `(15.625*2)^(i/3)` を計算していた（正しくは `15.625 * 2^(i/3)`）。
    さらにループ内でスカラーを代入していたので配列にならず、最後の 1 バンド分しか残らなかった。
    """
    if mf is None:
        mf = third_octave_bands()
    return 1.81e-8 * mf ** 1.57


def expand_6_to_32(energy6):
    """6 オクターブバンドのエネルギーを 1/3 オクターブ 32 バンドに展開する。

    元コード 94〜125 行。単純に同じ値をコピーするだけ（帯域内は平坦とみなす）。

    引数:  energy6  (count, 6) … パルス列のバンド別エネルギー
    戻り値:         (32, count) … 元コード enerred_air と同じ並び（バンドが行）
    """
    energy6 = np.asarray(energy6, dtype=float)
    if energy6.ndim != 2 or energy6.shape[1] != len(OCTAVE_BAND_FREQUENCIES):
        raise ValueError(
            f"エネルギーの形状が (count, 6) ではありません: {energy6.shape}。"
            f"バックトレースは 6 オクターブバンド"
            f"（{OCTAVE_BAND_FREQUENCIES.tolist()} Hz）を前提にしています")
    return energy6[:, BAND_EXPANSION_6_TO_32].T


def apply_air_absorption(energy32, time, sound_velocity=SOUND_VELOCITY, mair=None):
    """空気吸収による減衰を掛ける。元コード 127 行。

        enerred_air(j,:) = enerred_air(j,:) * exp(-mair(j) * c0 * rtime(:))

    c0 * rtime は伝搬距離なので、距離に応じてバンドごとに減衰する。

    引数:  energy32 (32, count)、time (count,)
    戻り値:         (32, count)

    2026-08-14 修正: 旧実装は `reflection_timing[i]`（i はバンド添字）を使っており、
    パルスごとの到来時刻ではなくバンド番号で時刻を引いていた。正しくは全パルス分のベクトル。
    """
    if mair is None:
        mair = airdamping_coefficient()
    time = np.asarray(time, dtype=float)
    return energy32 * np.exp(-mair[:, None] * sound_velocity * time[None, :])


def transfer_function(energy32_air, time, nfreq, df, sound_velocity=SOUND_VELOCITY,
                      chunk=4096):
    """パルス列を伝達関数（周波数応答）に変換する。元コード 136〜143 行。

        hfp(j,i) = Σ_k sqrt(enerred_air(j,k))
                   * exp(-i 2π df (i-1) rtime(k)) / (rtime(k) * c0)

    ・sqrt(エネルギー) → 振幅
    ・exp(-i 2π f t)   → 到来時刻ぶんの位相回転（＝時間遅れ）
    ・1/(rtime * c0)   → 1/距離 の距離減衰

    2026-08-14 修正: 旧実装は `transfer[j,i] = transfer[j,i] * ...` と**掛け算**に
    なっていた。元コードは `hfp(j,i) + ...` の**足し込み**。初期値 0 に掛けると
    伝達関数が全部 0 になるので、そのままでは無音の応答しか出ない。
    あわせて `df`（周波数離散化幅）が引数に無く未定義参照になっていたのも直した（E-4）。

    引数:
        energy32_air (32, count) 空気吸収込みのバンド別エネルギー
        time         (count,)    到来時刻 [s]
        nfreq        int         周波数の点数（= nn/2 + 1）
        df           float       周波数離散化幅 [Hz]
        chunk        int         周波数方向の分割数（メモリ対策。結果は変わらない）
    戻り値:
        (nfreq, 32) complex
    """
    time = np.asarray(time, dtype=float)
    # バンドによらない部分（exp）とバンドごとの部分（振幅）を分けるのが要点。
    # 素直に 3 重ループすると exp を 32 回計算し直すことになる
    amplitude = np.sqrt(energy32_air) / (time * sound_velocity)[None, :]   # (32, count)

    hfp = np.empty((nfreq, 32), dtype=np.complex128)
    for start in range(0, nfreq, chunk):
        stop = min(start + chunk, nfreq)
        index = np.arange(start, stop, dtype=float)
        phase = np.exp(-2j * np.pi * df * index[:, None] * time[None, :])  # (chunk, count)
        hfp[start:stop, :] = phase @ amplitude.T
    return hfp


def filter_bandpass(numtaps, wmin, wmax):
    """ハミング窓をかけた理想バンドパスの FIR 係数。元コード 241〜269 行 fir1_bandpass。

    wmin / wmax はナイキスト周波数で正規化した遮断周波数（0〜1）。
    元コードの引数 n はタップ数 - 1 だったが、ここは**タップ数**を渡す形にした
    （`numtaps = n + 1`）。呼び出し側の off-by-one を防ぐため。

        h[k] = sinc 差分 * (0.54 + 0.46 cos(2π(k - m)/n))      m = n/2, n = numtaps-1

    2026-08-14 修正:
      ・ハミング窓が `0.54 * 0.46 * cos(...)` と**掛け算**になっていた（正しくは足し算）。
        掛け算だと窓が常に正負に振れる 0.248*cos となり、窓の役目をしない。
      ・`n_shift = i - 1.0 - m` が 0 始まりの Python 添字に対して 1 ずれていた。
        元コードは `dble(i) - 1 - m` で i は 1 始まり。Python では `k - m`。
    """
    n = numtaps - 1
    m = n / 2.0
    k = np.arange(numtaps, dtype=float)
    shift = k - m

    filter_bp = np.empty(numtaps)
    zero = (shift == 0.0)
    # sin(πwx)/(πx) は x→0 で w に収束する。中央のタップだけ別扱い（元コード 260〜261 行）
    filter_bp[zero] = wmax - wmin
    s = shift[~zero]
    filter_bp[~zero] = (np.sin(np.pi * wmax * s) / (np.pi * s)
                        - np.sin(np.pi * wmin * s) / (np.pi * s))

    window = 0.54 + 0.46 * np.cos(np.pi * 2.0 * shift / n)
    return filter_bp * window


def bandpass_edges(mf, fmax):
    """各 1/3 オクターブバンドの正規化遮断周波数を返す。元コード 152〜154 行。

    下端 mf * 2^(-1/6)、上端 mf * 2^(1/6) をナイキスト周波数で割る。
    最上バンド（20 kHz）だけは上端がナイキストを超えるので 0.99 で頭打ちにする。

    2026-08-14 修正: 旧実装は `(mf*2)^(-1/6)` と括弧が崩れていた
    （正しくは `mf * 2^(-1/6)`）。上端も `(mf*2)^(-1/6/fmax)` と指数の中に
    fmax が入り込んでいた。
    """
    lower = mf * 2.0 ** (-1.0 / 6.0) / fmax
    upper = mf * 2.0 ** (1.0 / 6.0) / fmax
    upper[-1] = fmax * 0.99 / fmax
    return lower, upper


def extend_to_negative(hfp, nn):
    """正の周波数だけの伝達関数を負の周波数へ折り返す。元コード 173〜180 行。

        hf(j,i)      = hfp(j,i)                   i = 1 .. nn/2
        hf(j,i) = conj(hfp(j, nn-i+2))            i = nn/2+1 .. nn

    実信号のスペクトルはエルミート対称（負の周波数は正の周波数の複素共役）なので、
    こうしてから逆 FFT すると実数の時間波形になる。

    2026-08-14 修正: 旧実装の添字は `transfer[j, nn-i+1]` で、1 始まり前提の式を
    そのまま 0 始まりに持ち込んでいた。0 始まりでは `hfp[nn - i]`（i は 0 始まり）。

    引数:  hfp (nfreq, 32)
    戻り値:     (nn, 32) complex
    """
    half = nn // 2
    hf = np.empty((nn, hfp.shape[1]), dtype=np.complex128)
    hf[:half] = hfp[:half]
    hf[half:] = np.conj(hfp[nn - np.arange(half, nn)])
    return hf


def impulse_response(time, energy6, sound_velocity=SOUND_VELOCITY,
                     sampling_frequency=SAMPLING_FREQUENCY, max_time=MAX_TIME,
                     compensate_filter_delay=False, verbose=True):
    """パルス列からインパルス応答を合成する。

    引数:
        time      (count,)    到来時刻 [s]
        energy6   (count, 6)  オクターブバンド別エネルギー
        compensate_filter_delay
            True にするとバンドパスフィルタの直線位相遅れ (nn-1)/2 サンプルを
            打ち消して 0 秒始まりにする。元コードは補正していない（既定 False）。

    戻り値:
        (t, ir)   ともに (nn,)。t は時間ベクトル [s]
    """
    time = np.asarray(time, dtype=float)
    energy6 = np.asarray(energy6, dtype=float)
    if len(time) == 0:
        raise ValueError("パルス列が空です（受音に至った経路がありません）")

    # ---- 定数の決定（元コード 51〜67 行）----
    # データ点数は 2*tmax*fs 以上で最小の 2 のべき乗
    ncount = 1
    while 2 ** ncount < int(round(2.0 * max_time * sampling_frequency)):
        ncount += 1
    nn = 2 ** ncount
    fmax = sampling_frequency / 2.0
    nfreq = nn // 2 + 1
    df = sampling_frequency / nn
    dt = 1.0 / sampling_frequency

    mf = third_octave_bands()
    mair = airdamping_coefficient(mf)

    if verbose:
        print(f"[impulse] fs={sampling_frequency:.0f} Hz / nn={nn} / df={df:.4f} Hz / "
              f"全長 {nn * dt:.3f} s / パルス {len(time)} 本")

    # ---- 6 → 32 バンド展開と空気吸収（元コード 94〜128 行）----
    energy32 = expand_6_to_32(energy6)
    energy32 = apply_air_absorption(energy32, time, sound_velocity, mair)

    # ---- 伝達関数（元コード 136〜143 行）----
    hfp = transfer_function(energy32, time, nfreq, df, sound_velocity)

    # ---- バンドごとにフィルタを掛けて逆 FFT し、足し合わせる ----
    # 元コードは (32, nn) の配列を 4 本同時に持つが、nn=131072 だと 1 本 67MB になる。
    # バンド単位で回せば結果は同じで数 MB で済むのでそうしている。
    lower, upper = bandpass_edges(mf, fmax)
    hf = extend_to_negative(hfp, nn)

    ir = np.zeros(nn)
    for j in range(32):
        # 時間領域フィルタ → その FFT（元コード 152, 163 行）
        # 元コードの fft(x, n, 1) は exp(-i...) なので np.fft.fft と同じ
        htf = np.fft.fft(filter_bandpass(nn, lower[j], upper[j]))
        # 畳み込み（周波数領域の積）→ 逆 FFT（元コード 190, 201 行）
        # 元コードの fft(x, n, -1) は exp(+i...) かつ 1/n 正規化なので np.fft.ifft と同じ
        ir += np.real(np.fft.ifft(hf[:, j] * htf))

    # ---- 時間ベクトル（元コード 225 行）----
    t = dt * np.arange(nn)

    if compensate_filter_delay:
        # 直線位相 FIR の群遅延 (nn-1)/2 サンプルを巻き戻す
        delay = (nn - 1) // 2
        ir = np.roll(ir, -delay)
        # ★巻き戻すと、フィルタの「負の時刻側の応答（プリリンギング）」が
        #   バッファの末尾に回り込む。これは合成が周波数領域の掛け算＝巡回畳み込みで
        #   あることによる作り物なので、有効な範囲だけ残して切り落とす。
        #   残す長さ nn - delay は「回り込みが混ざらない最後の時刻」。
        #   既定では 1.486 s あり、設計上の最大時間 tmax = 1 s より長いので実害はない。
        valid = nn - delay
        ir = ir[:valid]
        t = t[:valid]
        if verbose:
            print(f"[impulse] フィルタの遅れ {delay * dt:.3f} s を補正し、"
                  f"回り込みを除いた {valid * dt:.3f} s ぶんを返します")
    elif verbose:
        print(f"[impulse] 注意: バンドパスの直線位相により応答は "
              f"{((nn - 1) // 2) * dt:.3f} s 遅れています"
              f"（compensate_filter_delay=True で補正できます）")

    return t, ir


def write_impulseresponce(filename, t, ir):
    """インパルス応答を CSV に保存する。元コード 232〜235 行。"""
    np.savetxt(filename, np.column_stack([t, ir]), delimiter=",",
               header="time_s,ir", comments="", fmt="%.12g")
    return filename


def impulse_responce(filename, pulses, sound_velocity=SOUND_VELOCITY,
                     sampling_frequency=SAMPLING_FREQUENCY, max_time=MAX_TIME,
                     compensate_filter_delay=False, verbose=True):
    """パルス列 → インパルス応答 → CSV 保存までの一括処理。

    引数:
        pulses : loop_noredundancy.PulseList（`.time` と `.energy` を使う）

    ※ 旧シグネチャは (filename, sound_velocity, reflection_timing,
       soundsourceenergy_list, frequency_number, count, nn, mf, fmax, dt) だったが、
       nn / mf / fmax / dt / frequency_number / count は元コードでは
       fs と tmax から**導出される値**であって外から与えるものではないため、
       サンプリング周波数と最大時間だけを受け取る形に整理した。
    """
    t, ir = impulse_response(pulses.time, pulses.energy,
                             sound_velocity=sound_velocity,
                             sampling_frequency=sampling_frequency,
                             max_time=max_time,
                             compensate_filter_delay=compensate_filter_delay,
                             verbose=verbose)
    write_impulseresponce(filename, t, ir)
    if verbose:
        print(f"[impulse] インパルス応答を書き出しました: {filename}")
    return t, ir
