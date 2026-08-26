"""音圧レベル（帯域別）と STI（音声伝送指数）。

## 何を解決するモジュールか

バックトレースが出すパルス列（`loop_noredundancy.PulseList`）は
**虚音源の一覧**で、各パルスのエネルギーには

- 反射での吸音（`sound_ray.energy_decay` の累積）

だけが入っている。**距離減衰と空気吸収は入っていない**（インパルス応答を
合成する段階で `impulse.transfer_function` が `1/(c·t)` を掛け、
`impulse.apply_air_absorption` が `exp(-m·d)` を掛ける）。

そのため「受音点で実際にどれだけのエネルギーを受け取るか」は
パルス列だけを見ても分からない。ここで 1 か所にまとめる。

    受け取るエネルギー   E_i = A_i · exp(-m·d_i) / (4π d_i²)      … `received_energy()`
    音圧レベル           Lp = Lw + 10 log10(Σ E_i) + 10 log10(ρc/400)
    変調伝達関数         m(F) = |Σ E_i e^{-j2πF t_i}| / Σ E_i    … STI の材料

**この 2 つが同じ `received_energy()` を使うのが肝**。
音圧レベルは Σ E_i、STI は E_i の時間分布を見ているだけで、素は同じ。

## 音圧レベル（依頼 2026-08-21）

> 音圧レベルを知りたい。無論、帯域毎に必要。例えば、無響室の計算で
> 逆二乗がどれくらい成り立っているかを確認したりしたい。
> 音源に点音源の PWL を与えれば、絶対値もわかればよいですが、
> 未入力の場合、相対値でも良いです。

点音源が自由空間に置かれ、距離 d の 1 点だけに音が届くなら
`Σ E_i = 1/(4π d²)` なので

    Lp = Lw − 10 log10(4π d²) + 10 log10(ρc/400) = Lw − 20 log10 d − 11 + 0.1

となり、**教科書の `Lp = Lw − 20 log10 r − 11` にそのまま一致する**。
つまり逆二乗則は式の上では厳密に成り立つので、無響室のモデルで確認すると
**幾何音響の実装が正しいかどうかの物差しになる**（反射面が無ければ直接音だけ、
反射面があればそのぶん上に外れる）。`freefield_level()` がその理論値を返す。

★**出力（`spl.csv`・図）には自由音場の値と差を入れない**
（2026-08-21 ユーザー判断。「そういう場合もある、というだけ」）。
逆二乗則そのものは `tests/test_geosim.py` [31] で押さえてある。

PWL 未入力なら **Lw = 0 dB として計算する**（＝相対値。W = 1 pW 基準）。
帯域ごとの相対関係と距離依存性はそのまま読めるので、逆二乗の確認には足りる。

## STI（依頼 2026-08-21「STI も算出して」）

IEC 60268-16 の音声伝送指数。**変調伝達関数 m(F) から出す。**

    m(F) = |∫ h²(t) e^{-j2πFt} dt| / ∫ h²(t) dt

`h²(t)` はエネルギー的インパルス応答なので、**パルス列の (時刻, エネルギー) が
そのまま使える**（帯域フィルタを掛けたインパルス応答から積分するより素直で、
フィルタの遅れや裾の影響を受けない）。幾何音響では普通この形で出す。

    実効 SNR   SNR_eff = 10 log10( m / (1 − m) )   [−15, +15] dB に丸める
    伝送指数   TI = (SNR_eff + 15) / 30
    帯域指数   MTI = 14 個の変調周波数の TI の平均
    STI        Σ αₖ·MTIₖ − Σ βₖ·√(MTIₖ·MTIₖ₊₁)

α・β は IEC 60268-16 の男声／女声の重み。**63 Hz は使わない**（STI は
125 Hz〜8 kHz の 7 帯域）。6 バンド（125〜4k）で計算したときは 8 kHz が
足りないので、**その帯域を除いて重みを正規化し直して警告する**。

背景騒音と聴覚マスキングは、**音源 PWL と騒音レベルの両方が入っているときだけ**
効かせる（絶対値が分からないと SNR も受聴閾値も決まらない）。
"""

import numpy as np

import absorption as ab
import table as tb
from atmosphere import Atmosphere

# 音圧レベルの基準。p0 = 20 μPa、W0 = 1 pW。
# W0/p0² = 1e-12/4e-10 = 1/400 なので、定数項は 10 log10(ρc/400) になる
REFERENCE_PRESSURE = 2.0e-5
REFERENCE_POWER = 1.0e-12

# A 特性の補正値 [dB]（63 Hz〜8 kHz。IEC 61672-1）
A_WEIGHTING = {63.0: -26.2, 125.0: -16.1, 250.0: -8.6, 500.0: -3.2,
               1000.0: 0.0, 2000.0: 1.2, 4000.0: 1.0, 8000.0: -1.1}

# ---- STI（IEC 60268-16）------------------------------------------------------

# 変調周波数 [Hz]（1/3 オクターブ間隔で 14 個）
MODULATION_FREQUENCIES = np.array([0.63, 0.80, 1.00, 1.25, 1.60, 2.00, 2.50,
                                   3.15, 4.00, 5.00, 6.30, 8.00, 10.00, 12.50])

# STI が使うオクターブバンド [Hz]（**63 Hz は含まない**）
STI_BANDS = np.array([125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0])

# 帯域の重み。α は各帯域の寄与、β は隣り合う帯域の冗長性を差し引く分。
# 男声は 125 Hz から、女声は 125 Hz の寄与がゼロ（基本周波数が上にあるため）
STI_WEIGHTS = {
    "male":   {"alpha": np.array([0.085, 0.127, 0.230, 0.233, 0.309, 0.224, 0.173]),
               "beta":  np.array([0.085, 0.078, 0.065, 0.011, 0.047, 0.095])},
    "female": {"alpha": np.array([0.000, 0.117, 0.223, 0.216, 0.328, 0.250, 0.194]),
               "beta":  np.array([0.000, 0.099, 0.066, 0.062, 0.025, 0.076])},
}

# 音声の受聴閾値 [dB]（IEC 60268-16 表。125 Hz〜8 kHz）
SPEECH_RECEPTION_THRESHOLD = np.array([46.0, 27.0, 12.0, 6.5, 7.5, 8.0, 12.0])

# STI の評価区分（IEC 60268-16 / JIS）
STI_RATINGS = ((0.75, "優 (excellent)"), (0.60, "良 (good)"),
               (0.45, "可 (fair)"), (0.30, "不可 (poor)"), (0.0, "劣 (bad)"))


def sti_rating(value):
    """STI の値から評価区分の名前を返す。"""
    if value is None or np.isnan(value):
        return ""
    for limit, name in STI_RATINGS:
        if value >= limit:
            return name
    return STI_RATINGS[-1][1]


def a_weighting(frequencies):
    """A 特性の補正値 [dB]。表に無い周波数は対数補間する。"""
    known = np.array(sorted(A_WEIGHTING))
    values = np.array([A_WEIGHTING[f] for f in known])
    return np.interp(np.log10(np.asarray(frequencies, dtype=float)),
                     np.log10(known), values)


def level_constant(atmosphere=None):
    """`10 log10(ρc/400)` [dB]。20℃ で +0.12 dB 程度。

    音響インピーダンス ρc が 400 N·s/m³ からずれる分の補正。
    教科書の `Lp = Lw − 20 log10 r − 11` はこの項を落とした形。
    """
    if atmosphere is None:
        atmosphere = Atmosphere()
    rho_c = atmosphere.density * atmosphere.sound_velocity
    return 10.0 * np.log10(rho_c / 400.0)


def received_energy(times, energies, distances=None, atmosphere=None,
                    frequencies=None):
    """パルスごと・帯域ごとに「受音点で受け取るエネルギー」を返す (n, nf)。

        E_i = A_i · exp(-m·d_i) / (4π d_i²)

    A_i はパルス列が持っているエネルギー（反射の吸音だけ入っている）、
    `exp(-m·d)` が空気吸収、`1/(4π d²)` が点音源の距離減衰。
    **この形はインパルス応答の合成（`impulse.transfer_function` の
    `sqrt(E)/(c·t)`）と同じもの**を、エネルギーのまま書いたもの。
    4π はインパルス応答側では定数として省かれているが、
    音圧レベルの絶対値には要るのでここでは掛ける。

    引数:
        times     (n,)      到来時刻 [s]
        energies  (n, nf)   パルス列のバンド別エネルギー
        distances (n,) | None  経路長 [m]。None なら `時刻 × 音速`
    """
    if atmosphere is None:
        atmosphere = Atmosphere()
    times = np.atleast_1d(np.asarray(times, dtype=float))
    energies = np.atleast_2d(np.asarray(energies, dtype=float))
    if distances is None:
        distances = times * atmosphere.sound_velocity
    distances = np.atleast_1d(np.asarray(distances, dtype=float))
    if frequencies is None:
        frequencies = ab.octave_bands(energies.shape[1])
    frequencies = np.asarray(frequencies, dtype=float)

    m = atmosphere.absorption_coefficient(frequencies)          # (nf,)
    air = np.exp(-m[None, :] * distances[:, None])              # (n, nf)
    spreading = 1.0 / (4.0 * np.pi * np.maximum(distances, 1e-12) ** 2)
    return energies * air * spreading[:, None]


def freefield_level(distance, source_power_db=None, atmosphere=None,
                    band_number=None, frequencies=None, air_absorption=True):
    """点音源の自由音場での音圧レベル [dB]（`Lw − 20 log10 d − 11`）。

    **逆二乗則の物差し**。反射面が無いモデルなら計算結果がこれに一致する。
    空気吸収も入れる（入れないと高域で理論値のほうが高く出る）。
    """
    if atmosphere is None:
        atmosphere = Atmosphere()
    if frequencies is None:
        frequencies = ab.octave_bands(band_number or 8)
    frequencies = np.asarray(frequencies, dtype=float)
    power = _power_levels(source_power_db, len(frequencies))
    air = (10.0 * np.log10(np.exp(-atmosphere.absorption_coefficient(frequencies)
                                  * distance))
           if air_absorption else 0.0)
    return (power - 10.0 * np.log10(4.0 * np.pi * distance ** 2)
            + level_constant(atmosphere) + air)


def _power_levels(source_power_db, band_number):
    """音源パワーレベルを (nf,) の配列にする。None なら 0 dB（相対値）。"""
    if source_power_db is None:
        return np.zeros(band_number)
    values = np.atleast_1d(np.asarray(source_power_db, dtype=float)).ravel()
    if len(values) == 1:
        return np.full(band_number, float(values[0]))
    return np.array((list(values) + [values[-1]] * band_number)[:band_number])


def band_levels(times, energies, distances=None, atmosphere=None,
                frequencies=None, source_power_db=None, source_distance=None,
                early_limit=None, verbose=True):
    """帯域別の音圧レベルを求める。

    引数:
        source_power_db : float | (nf,) | None
            点音源のパワーレベル PWL [dB]。**None なら Lw = 0 dB（相対値）**
        source_distance : float | None
            音源から受音点までの直線距離 [m]。自由音場の理論値に使う。
            None ならパルス列のうち最短の経路長（＝直接音）を使う
        early_limit : float | None
            「初期」と「後期」を分ける時刻 [s]（直接音からの相対）。既定 50 ms

    戻り値: dict
        'frequencies' / 'levels'（(nf,) 音圧レベル [dB]）
        'direct' / 'reflected' / 'early' / 'late' … 成分ごとの (nf,)
        'freefield'（自由音場の理論値 (nf,)）/ 'excess'（差 = levels − freefield）
        'overall'（帯域を合成した値 [dB]）/ 'overall_a'（A 特性 [dB(A)]）
        'a_weighted'（帯域別の A 特性補正後 (nf,)）
        'source_power'（使った PWL (nf,)）/ 'relative'（PWL 未入力なら True）
        'source_distance' / 'direct_time'
    """
    if atmosphere is None:
        atmosphere = Atmosphere()
    times = np.atleast_1d(np.asarray(times, dtype=float))
    energies = np.atleast_2d(np.asarray(energies, dtype=float))
    if frequencies is None:
        frequencies = ab.octave_bands(energies.shape[1])
    frequencies = np.asarray(frequencies, dtype=float)
    if distances is None:
        distances = times * atmosphere.sound_velocity
    distances = np.atleast_1d(np.asarray(distances, dtype=float))

    received = received_energy(times, energies, distances, atmosphere, frequencies)
    power = _power_levels(source_power_db, len(frequencies))
    constant = level_constant(atmosphere)

    def to_level(total):
        with np.errstate(divide="ignore"):
            return np.where(total > 0.0, power + 10.0 * np.log10(total) + constant,
                            -np.inf)

    # 直接音＝いちばん早く届くパルス（反射回数 0 のものがあればそれ）
    direct_index = int(np.argmin(times))
    direct_time = float(times[direct_index])
    if source_distance is None:
        source_distance = float(distances[direct_index])
    limit = direct_time + (0.05 if early_limit is None else early_limit)

    direct = np.zeros(len(frequencies))
    direct[:] = received[direct_index]
    reflected = received.sum(axis=0) - direct
    early = received[times <= limit].sum(axis=0)
    late = received[times > limit].sum(axis=0)

    levels = to_level(received.sum(axis=0))
    freefield = freefield_level(source_distance, source_power_db, atmosphere,
                                frequencies=frequencies)
    weighted = levels + a_weighting(frequencies)

    def combine(values):
        finite = values[np.isfinite(values)]
        if not len(finite):
            return float("-inf")
        return float(10.0 * np.log10(np.sum(10.0 ** (finite / 10.0))))

    result = {
        "frequencies": frequencies,
        "levels": levels,
        "direct": to_level(direct),
        "reflected": to_level(reflected),
        "early": to_level(early),
        "late": to_level(late),
        "freefield": freefield,
        "excess": levels - freefield,
        "a_weighted": weighted,
        "overall": combine(levels),
        "overall_a": combine(weighted),
        "source_power": power,
        "relative": source_power_db is None,
        "source_distance": float(source_distance),
        "direct_time": direct_time,
    }
    if verbose:
        kind = ("相対値（PWL 未入力なので Lw = 0 dB として計算）"
                if result["relative"] else "絶対値（PWL 指定あり）")
        print(f"[音圧レベル] {kind} / 音源距離 {source_distance:.3f} m / "
              f"パルス {len(times)} 本")
        tb_rows = [("Lp[dB]", levels), ("直接音", result["direct"]),
                   ("反射音", result["reflected"])]
        print("[音圧レベル] " + " " * 12 + "".join(f"{f:>10.0f}" for f in frequencies))
        for name, values in tb_rows:
            print(f"[音圧レベル] {name:>12}"
                  + "".join("      -inf" if not np.isfinite(v) else f"{v:10.2f}"
                            for v in values))
        print(f"[音圧レベル] 帯域合成 {result['overall']:.2f} dB / "
              f"A 特性 {result['overall_a']:.2f} dB(A)")
    return result


# バンドの中で何点の周波数を計算して平均するか（複素和のとき）
COHERENT_LINES = 129


def coherent_band_levels(times, energies, distances=None, atmosphere=None,
                         frequencies=None, band_width="1/1",
                         lines=COHERENT_LINES, source_power_db=None):
    """**位相ごと重ねた**バンド別の音圧レベル [dB]（複素和）。

        p(f) = Σ_n √(A_n · e^{-m d_n} / 4π) / d_n · e^{-j 2π f d_n / c}
        Lp   = Lw + 10log10( <|p(f)|²> ) + 10log10(ρc/400)

    `< >` はバンドの中の周波数平均。`band_levels()`（エネルギー和）と
    **同じ量**を、足すときに位相を持たせただけのもの——干渉の山谷が出る。

    ★**反射の位相ずれは 0**（＝剛な面では厳密、吸音面では仮定）。
      材料の位相情報が無い段階では、山谷の位置ではなく**振れ幅の目安**として見る。
    """
    if atmosphere is None:
        atmosphere = Atmosphere()
    times = np.atleast_1d(np.asarray(times, dtype=float))
    energies = np.atleast_2d(np.asarray(energies, dtype=float))
    if distances is None:
        distances = times * atmosphere.sound_velocity
    distances = np.atleast_1d(np.asarray(distances, dtype=float))
    if frequencies is None:
        frequencies = ab.octave_bands(energies.shape[1])
    frequencies = np.asarray(frequencies, dtype=float)

    velocity = atmosphere.sound_velocity
    m = atmosphere.absorption_coefficient(frequencies)
    lower, upper = ab.band_edges(frequencies, band_width)
    power = 0.0 if source_power_db is None else float(np.mean(
        np.atleast_1d(np.asarray(source_power_db, dtype=float))))
    impedance = 10.0 * np.log10(atmosphere.density * velocity / 400.0)

    levels = np.empty(len(frequencies))
    for band in range(len(frequencies)):
        air = np.exp(-m[band] * distances)
        amplitude = np.sqrt(np.maximum(energies[:, band], 0.0) * air
                            / (4.0 * np.pi)) / distances
        line = np.linspace(lower[band], upper[band], lines)
        wave = 2.0 * np.pi * line / velocity
        pressure = np.exp(-1j * wave[:, None] * distances[None, :]) @ amplitude
        levels[band] = (power + 10.0 * np.log10(float(np.mean(
            np.abs(pressure) ** 2)) + 1.0e-300) + impedance)
    return levels


def write_levels(filename, result):
    """音圧レベルを CSV に保存する。**周波数は横**（区分付きの表）。"""
    # ★自由音場（逆二乗）の値と差は**出さない**（2026-08-21 ユーザー判断。
    #   「そういう場合もある、というだけ」）。計算はできるので
    #   `freefield_level()` は残してあり、逆二乗則の検算にテストで使っている
    rows = [
        ("音圧レベル", "Lp_dB", result["overall"], result["levels"]),
        ("音圧レベル", "Lp_A_dB", result["overall_a"], result["a_weighted"]),
        ("内訳", "直接音_dB", None, result["direct"]),
        ("内訳", "反射音_dB", None, result["reflected"]),
        ("内訳", "初期_50ms_dB", None, result["early"]),
        ("内訳", "後期_dB", None, result["late"]),
        ("参考", "音源パワーレベル_dB", result["source_distance"],
         result["source_power"]),
    ]
    # ★複素和（位相を含む）を出しているときは、区分を分けて一緒に入れる
    if result.get("coherent") is not None:
        rows.insert(2, ("複素和", "Lp_dB", None, result["coherent"]))
    return tb.write_sectioned_table(filename, result["frequencies"], rows,
                                    value_label="総合")


# ------------------------------------------------------------------------------
# STI（音声伝送指数）
# ------------------------------------------------------------------------------

def modulation_transfer(times, energies, modulation_frequencies=None):
    """変調伝達関数 m(F) を求める (nf, nF)。

        m(F) = |Σ E_i e^{-j2πF t_i}| / Σ E_i

    `energies` は**受音点で受け取るエネルギー**（`received_energy()` の出力）。
    エネルギー的インパルス応答のフーリエ変換そのものなので、
    パルス列から直接出せる（帯域フィルタを通す必要がない）。
    """
    times = np.atleast_1d(np.asarray(times, dtype=float))
    energies = np.atleast_2d(np.asarray(energies, dtype=float))
    if modulation_frequencies is None:
        modulation_frequencies = MODULATION_FREQUENCIES
    modulation_frequencies = np.asarray(modulation_frequencies, dtype=float)

    phase = np.exp(-2j * np.pi * modulation_frequencies[:, None] * times[None, :])
    numerator = np.abs(phase @ energies)                # (nF, nf)
    total = energies.sum(axis=0)                        # (nf,)
    with np.errstate(divide="ignore", invalid="ignore"):
        m = np.where(total > 0.0, numerator / total, np.nan)
    return m.T                                          # (nf, nF)


def _band_index(frequencies):
    """STI が使う 125 Hz〜8 kHz の帯域が `frequencies` の何番目かを返す。

    見つからない帯域は -1。**6 バンド（125〜4k）では 8 kHz が無い**ので、
    そこは重みを外して正規化し直す。
    """
    frequencies = np.asarray(frequencies, dtype=float)
    index = []
    for want in STI_BANDS:
        hit = np.nonzero(np.isclose(frequencies, want, rtol=0.05))[0]
        index.append(int(hit[0]) if len(hit) else -1)
    return np.array(index)


def speech_transmission_index(times, energies, distances=None, atmosphere=None,
                              frequencies=None, source_power_db=None,
                              noise_level_db=None, sex="male",
                              include_masking=True, verbose=True):
    """STI（音声伝送指数）を求める。IEC 60268-16。

    引数:
        times / energies / distances … パルス列（`PulseList` の中身）。
            エネルギーは**距離減衰の入っていない生の値**を渡す
            （ここで `received_energy()` を通す）
        source_power_db / noise_level_db :
            両方そろっているときだけ**背景騒音と聴覚マスキング**を効かせる。
            片方でも無ければ「騒音なし・理想の SNR」として計算する
        sex : 'male' | 'female'  帯域の重み（既定 男声）

    戻り値: dict
        'sti'（総合値）/ 'rating'（評価区分）/ 'sti_female'
        'bands'（使った帯域の中心周波数）/ 'mti'（帯域ごとの指数）
        'snr_effective'（帯域ごとの実効 SNR [dB]）/ 'mtf'（(帯域, 変調周波数)）
        'modulation_frequencies' / 'weights_note'（重みを直したときの説明）
    """
    if atmosphere is None:
        atmosphere = Atmosphere()
    energies = np.atleast_2d(np.asarray(energies, dtype=float))
    if frequencies is None:
        frequencies = ab.octave_bands(energies.shape[1])
    frequencies = np.asarray(frequencies, dtype=float)

    received = received_energy(times, energies, distances, atmosphere, frequencies)
    index = _band_index(frequencies)
    available = index >= 0
    note = ""
    if not np.all(available):
        missing = ", ".join(f"{f:.0f} Hz" for f in STI_BANDS[~available])
        note = (f"{missing} が計算対象のバンドに無いので、その帯域の重みを外して"
                f"正規化しました（STI は 125 Hz〜8 kHz の 7 帯域が本来の形）")
        if verbose:
            print(f"[STI] 注意: {note}")

    mtf = np.full((len(STI_BANDS), len(MODULATION_FREQUENCIES)), np.nan)
    if np.any(available):
        mtf[available] = modulation_transfer(times, received[:, index[available]])

    # ---- 背景騒音・聴覚マスキング（絶対値が分かるときだけ）----
    levels = None
    if source_power_db is not None and noise_level_db is not None:
        levels = band_levels(times, energies, distances, atmosphere, frequencies,
                             source_power_db=source_power_db,
                             verbose=False)["levels"]
        mtf = _apply_noise_and_masking(mtf, levels[index], noise_level_db,
                                       index, available, include_masking)

    # ---- m(F) → 実効 SNR → 伝送指数 ----
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = 10.0 * np.log10(np.clip(mtf, 1e-12, 1.0 - 1e-12)
                              / (1.0 - np.clip(mtf, 1e-12, 1.0 - 1e-12)))
    snr = np.clip(snr, -15.0, 15.0)
    ti = (snr + 15.0) / 30.0
    mti = np.nanmean(ti, axis=1)

    result = {"bands": STI_BANDS, "mti": mti, "mtf": mtf,
              "snr_effective": np.nanmean(snr, axis=1),
              "modulation_frequencies": MODULATION_FREQUENCIES,
              "weights_note": note, "levels": levels,
              "noise_used": levels is not None}
    for kind in ("male", "female"):
        result[f"sti_{kind}"] = _combine_bands(mti, available, kind)
    result["sti"] = result[f"sti_{sex}"]
    result["rating"] = sti_rating(result["sti"])
    result["sex"] = sex

    if verbose:
        print(f"[STI] STI = {result['sti']:.3f}（{result['rating']}）"
              f" / 男声 {result['sti_male']:.3f} 女声 {result['sti_female']:.3f}"
              f" / 騒音{'あり' if result['noise_used'] else 'なし（理想）'}")
        print("[STI] " + " " * 10 + "".join(f"{f:>9.0f}" for f in STI_BANDS))
        print("[STI] " + f"{'MTI':>10}"
              + "".join("      ---" if np.isnan(v) else f"{v:9.3f}" for v in mti))
    return result


def _apply_noise_and_masking(mtf, band_levels_db, noise_level_db, index,
                             available, include_masking):
    """背景騒音・聴覚マスキング・受聴閾値で m(F) を下げる。

        m' = m · I / (I + I_noise + I_masking + I_threshold)

    マスキングは**下の帯域から上の帯域へ**掛かる（IEC 60268-16）。
    マスキングの量は下の帯域のレベルで決まる（大きいほど広がる）。
    """
    levels = np.array(band_levels_db, dtype=float)
    noise = _power_levels(noise_level_db, len(STI_BANDS))
    signal = np.where(np.isfinite(levels), 10.0 ** (levels / 10.0), 0.0)
    noise_intensity = 10.0 ** (noise / 10.0)

    masking = np.zeros(len(STI_BANDS))
    if include_masking:
        for k in range(1, len(STI_BANDS)):
            level = levels[k - 1]
            if not np.isfinite(level):
                continue
            if level < 63.0:
                slope = 0.5 * level - 65.0
            elif level < 67.0:
                slope = 1.8 * level - 146.9
            elif level < 100.0:
                slope = 0.5 * level - 59.8
            else:
                slope = -10.0
            masking[k] = signal[k - 1] * 10.0 ** (slope / 10.0)

    threshold = 10.0 ** (SPEECH_RECEPTION_THRESHOLD / 10.0)
    total = signal + noise_intensity + masking + threshold
    with np.errstate(divide="ignore", invalid="ignore"):
        factor = np.where(total > 0.0, signal / total, 0.0)
    return mtf * factor[:, None]


def _combine_bands(mti, available, sex):
    """帯域ごとの MTI を重み付けして STI にする。

        STI = Σ αₖ·MTIₖ − Σ βₖ·√(MTIₖ·MTIₖ₊₁)

    使えない帯域（6 バンド計算の 8 kHz など）は**α を外して正規化し直す**。
    β は両側の帯域がそろっているものだけ残す。
    """
    alpha = STI_WEIGHTS[sex]["alpha"].copy()
    beta = STI_WEIGHTS[sex]["beta"].copy()
    usable = np.asarray(available) & ~np.isnan(mti)
    alpha[~usable] = 0.0
    for k in range(len(beta)):
        if not (usable[k] and usable[k + 1]):
            beta[k] = 0.0
    scale = alpha.sum() - beta.sum()
    if scale <= 0.0:
        return float("nan")
    values = np.nan_to_num(mti)
    total = (np.sum(alpha * values)
             - np.sum(beta * np.sqrt(np.maximum(values[:-1] * values[1:], 0.0))))
    return float(np.clip(total / scale, 0.0, 1.0))


def write_sti(filename, result):
    """STI を CSV に保存する。**周波数は横**（区分付きの表）。

    総合値（STI・評価）は帯域に分けられないので 3 列目に置く。
    変調伝達関数は 1 行が変調周波数、列がオクターブバンドになる。
    """
    bands = result["bands"]
    rows = [("総合", "STI", "%.3f" % result["sti"], None),
            ("総合", "評価", result["rating"], None),
            ("総合", "STI_男声", "%.3f" % result["sti_male"], None),
            ("総合", "STI_女声", "%.3f" % result["sti_female"], None),
            ("帯域別", "MTI", None, result["mti"]),
            ("帯域別", "実効SNR_dB", None, result["snr_effective"])]
    for j, f in enumerate(result["modulation_frequencies"]):
        rows.append(("変調伝達関数", f"m_{f:g}Hz", None, result["mtf"][:, j]))
    if result["weights_note"]:
        rows.append(("備考", "重みの調整", result["weights_note"], None))
    return tb.write_sectioned_table(filename, bands, rows, value_label="総合")
