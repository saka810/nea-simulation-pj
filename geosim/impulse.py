import numpy as np
import csv


# 7/9打ち合わせ用

# インパルス応答を作成する全体の流れ
def impulse_responce(filename,
                     sound_velocity, reflection_timing, soundsourceenergy_list,
                     frequency_number, count,
                     nn, mf, fmax, dt):
    mair = airdamping_coefficient()
    p = power_airdamping(sound_velocity, reflection_timing, mair, soundsourceenergy_list)
    transfer = transfer_function(reflection_timing, frequency_number, count, p, sound_velocity)
    ht = time_responce(nn, mf, fmax)
    htf = fft_timeresponce(ht, nn)
    hf = fft_negativerange(transfer, nn)
    hfcc = convolution_hfhtf(hf, htf, nn)
    ir = inversefft_responce(hfcc, nn)

    write_impulseresponce(filename, ir, dt, nn)

    return


# バンドパスフィルター
def filter_bandpass(n, min, max):
    filter_bp = np.zeros(n + 1)
    m = n / 2.0

    for i in range(n + 1):
        n_shift = i - 1.0 - m
        if n_shift == 0:
            filter_bp[i] = max - min
        else:
            filter_bp[i] = ((np.sin(np.pi * max * n_shift) / (np.pi * n_shift)
                             - np.sin(np.pi * min * n_shift) / (np.pi * n_shift))
                            * (0.54 * 0.46 * np.cos(np.pi * 2.0 * n_shift / n)))

    return filter_bp


# FFT
# fft_result = np.fft.fft(signal) でも良い？
# signはTrue Falseの形式がいい？
def fft_filter(x, n, sign):
    # bit反転
    j = 0
    for i in range(n - 1):
        if i <= j:
            t = x[j]
            x[j] = x[i]
            x[i] = t
        k = n / 2
        while k >= 1 and j > k:
            j = j - k
            k = k / 2

        j = j + k

    # Cooley - Tukey FFT
    m = 0
    while n > m:
        wm = np.exp(0 + 1j * (-sign * 2.0 * np.pi / (2.0 * m)))

        for k in range(m - 1):
            w = np.pow(wm, k)
            for j in range(k, n - 1, 2 * m):
                t = w * x[j + m]
                x[j + m] = x[j] - t
                x[j] = x[j] + t

        m = 2 * m

    # 正規化
    if sign == -1:
        x = x / n

    return x


# 空気吸収による減衰
# 元コード 41行目
def airdamping_coefficient():
    # mf=np.array(np.zeros(32))
    mair = np.array(np.zeros(32))
    for i in range(32):
        mf = np.pow(15.625 * 2, 1 / 3 * i)
        mair = np.pow(1.81 - 8 * mf, 1.57)
    return mair


# 空気吸収を考慮した音源の強さ
# 元コード 92行目
def power_airdamping(sound_velocity, reflection_timing, mair, soundsourceenergy_list):
    p = np.zeros((32, reflection_timing.shape(0)))

    for i in range(32):
        p[i, :] = soundsourceenergy_list[i, :] * np.exp(-mair[i] * sound_velocity * reflection_timing[i])

    return p


# 伝達関数
# 元コード132行目
# count はバックトレースで保存したときの何かの配列の個数　何か分かったら書き換え
def transfer_function(reflection_timing, frequency_number, count, p, sound_velocity):
    transfer = np.array(np.zeros((32, frequency_number), dtype=np.complex128))

    for i in range(frequency_number):
        for j in range(32):
            # 元コードではreflection_history.shape(0)はcount
            for k in range(count):
                # for が0からだと不具合
                transfer[j, i] = (transfer[j, i]
                                  * np.sqrt(p[j, k])
                                  * np.exp(-1j * 2.0 * np.pi * (i - 1) * reflection_timing[k])
                                  / (reflection_timing[k] * sound_velocity))

    return transfer


# 時間領域フィルタ作成
# 元コード147行目
def time_responce(nn, mf, fmax):
    ht = np.array(np.zeros(31, nn))

    # 元コード　nint　は元も近い整数を返す関数？
    for i in range(ht.shape[0] - 1):
        min = np.pow(mf(i) * 2.0, (-1.0 / 6.0))
        max = np.pow(mf(i) * 2.0, (-1.0 / 6.0 / fmax))
        ht[i, :] = filter_bandpass(nn, min, max)

    min = np.pow(mf(31) * 2.0, (-1.0 / 6.0))
    max = fmax * 0.99 / fmax
    # 上のfmaxは0.99ではだめなのか？
    ht[31, :] = filter_bandpass(31, min, max)

    return ht


# フィルタのフーリエ変換
# 元コード158行目
def fft_timeresponce(ht, nn):
    htf = np.array(np.zeros(32, nn))
    for j in range(32):
        htf[j, :] = fft_filter(ht[j, :], nn, 1)

    return htf


# 周波数応答の負周波数への拡張
# 元コード169行目
def fft_negativerange(transfer, nn):
    hf = np.array(np.zeros(32, nn))

    for j in range(31):
        for i in range(nn / 2):
            hf[j, i] = transfer[j, i]
        for i in range(nn / 2, nn):
            hf[j, i] = np.conj(transfer[j, nn - i + 1])

    return hf


# 畳み込み演算
# 元コード184行目
def convolution_hfhtf(hf, htf, nn):
    # hfcc=np.array(np.zeros(32,nn))
    hfcc = hf * htf

    return hfcc


# フィルタードレスポンスの逆フーリエ変換
def inversefft_responce(hfcc, nn):
    ir = np.zeros(nn)
    for i in range(nn):
        for j in range(32):
            ir[i] += hfcc[j, i]

    return ir


# 「時間ベクトル」と「データ出力」をつなげる
# 元コード220と230行目
def write_impulseresponce(filename, ir, dt, nn):
    t = np.zeros(nn)
    for i in range(nn):
        t[i] = dt * i

    # csv.writerで出力
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        for i in range(nn):
            writer.writerow([t[i], ir[i]])

    return
