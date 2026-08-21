"""大気の状態（温度・湿度・気圧）から音速と空気吸収を求める。

元コードは音速を `c0 = 340.0` の定数、空気吸収を `mair = 1.81e-8 * f^1.57` の
べき乗近似（「20℃・湿度 40%」の読み取り値をフィットしたもの）で固定していた。
本モジュールはどちらも**温度・湿度・気圧から計算する**。

  ・音速     … 湿り空気の状態方程式から（比熱比とモル質量が湿度で変わる）
  ・空気吸収 … **ISO 9613-1**（音の大気吸収の規格）

基準状態は `Atmosphere()` の既定値（20℃ / 湿度 40% / 101.325 kPa）。
40% は元コードの空気吸収近似が想定していた条件に合わせてある。

使い方:

    from atmosphere import Atmosphere
    air = Atmosphere()                          # 20℃ / 40% / 101.325 kPa
    air = Atmosphere(temperature=25.0, humidity=60.0)
    air.sound_velocity                          # -> 346.6 [m/s]
    air.absorption_coefficient([125, 1000])     # -> [1/m]（エネルギーの減衰係数）
"""

import numpy as np

# 基準状態。GUI から変えられるように、ここを既定値の置き場にしている
REFERENCE_TEMPERATURE = 20.0        # [℃]
REFERENCE_HUMIDITY = 40.0           # [%RH]
REFERENCE_PRESSURE = 101.325        # [kPa]

# 物理定数
GAS_CONSTANT = 8.31446261815324     # 気体定数 R [J/(mol K)]
MOLAR_MASS_DRY = 0.0289645          # 乾燥空気のモル質量 [kg/mol]
MOLAR_MASS_WATER = 0.01801528       # 水蒸気のモル質量 [kg/mol]
MOLAR_CP_DRY = 29.07                # 乾燥空気の定圧モル比熱 [J/(mol K)]
MOLAR_CP_WATER = 33.60              # 水蒸気の定圧モル比熱 [J/(mol K)]

KELVIN = 273.15
ISO_REFERENCE_TEMPERATURE = 293.15  # ISO 9613-1 の基準温度 T0 [K]
TRIPLE_POINT = 273.16               # 水の三重点 T01 [K]

# dB → ネーパ（エネルギーの指数減衰係数）への換算。
# E = E0 * 10^(-α_dB * d / 10) = E0 * exp(-m * d)  なので m = α_dB * ln(10) / 10
DB_TO_NEPER = np.log(10.0) / 10.0


def saturation_vapour_pressure(temperature, pressure=REFERENCE_PRESSURE):
    """飽和水蒸気圧 [kPa]。ISO 9613-1 の式を使う。

        p_sat / p_r = 10^(-6.8346 (T01/T)^1.261 + 4.6151)

    20℃ で 2.336 kPa（理科年表の 2.339 kPa とほぼ一致）。
    音速と空気吸収の両方で同じ式を使うことで、両者の整合を取っている。
    """
    t_kelvin = np.asarray(temperature, dtype=float) + KELVIN
    exponent = -6.8346 * (TRIPLE_POINT / t_kelvin) ** 1.261 + 4.6151
    return REFERENCE_PRESSURE * 10.0 ** exponent


def water_vapour_mole_fraction(temperature, humidity, pressure=REFERENCE_PRESSURE):
    """空気中の水蒸気のモル分率 x_w（0〜1）。

    相対湿度は「その温度での飽和水蒸気圧に対する割合」なので、
    実際の水蒸気分圧は p_w = RH/100 * p_sat。それを全圧で割ればモル分率になる。
    """
    p_water = humidity / 100.0 * saturation_vapour_pressure(temperature, pressure)
    return p_water / pressure


def sound_velocity(temperature=REFERENCE_TEMPERATURE, humidity=REFERENCE_HUMIDITY,
                   pressure=REFERENCE_PRESSURE):
    """湿り空気の音速 [m/s]。

        c = sqrt(γ R T / M)

    **なぜ湿度で音速が変わるのか**：水蒸気（M = 18.0 g/mol）は乾燥空気（28.96 g/mol）より
    軽いので、湿気を含むほど空気の平均モル質量 M が小さくなり、音速は**速くなる**。
    比熱比 γ も水蒸気の混入でわずかに下がるが、モル質量の効果のほうが大きい。

    温度の効果のほうがずっと大きく、20℃ 付近では 1℃ あたり約 0.6 m/s。

    | 条件 | 音速 |
    |---|---|
    | 0℃ / 0% | 331.5 m/s |
    | 20℃ / 0% | 343.3 m/s |
    | 20℃ / 40% | 343.8 m/s |
    | 20℃ / 80% | 344.4 m/s |
    | 30℃ / 60% | 350.2 m/s |

    元コードの 340.0 m/s はおよそ 14℃ 相当。
    """
    t_kelvin = np.asarray(temperature, dtype=float) + KELVIN
    x_water = water_vapour_mole_fraction(temperature, humidity, pressure)

    molar_mass = (1.0 - x_water) * MOLAR_MASS_DRY + x_water * MOLAR_MASS_WATER
    molar_cp = (1.0 - x_water) * MOLAR_CP_DRY + x_water * MOLAR_CP_WATER
    molar_cv = molar_cp - GAS_CONSTANT        # マイヤーの関係 Cp - Cv = R
    gamma = molar_cp / molar_cv

    return float(np.sqrt(gamma * GAS_CONSTANT * t_kelvin / molar_mass))


def density(temperature=REFERENCE_TEMPERATURE, humidity=REFERENCE_HUMIDITY,
            pressure=REFERENCE_PRESSURE):
    """湿り空気の密度 ρ [kg/m³]。

        ρ = p M / (R T)      M は湿り空気の平均モル質量

    音圧レベルの絶対値に要る（`sound_level.level_constant` の `10 log10(ρc/400)`）。
    20℃ / 湿度 40% / 101.325 kPa で 1.199 kg/m³。**湿るほど軽くなる**
    （水蒸気のモル質量 18.0 g/mol は乾燥空気 28.96 g/mol より小さい）。
    """
    t_kelvin = np.asarray(temperature, dtype=float) + KELVIN
    x_water = water_vapour_mole_fraction(temperature, humidity, pressure)
    molar_mass = (1.0 - x_water) * MOLAR_MASS_DRY + x_water * MOLAR_MASS_WATER
    return float(np.asarray(pressure, dtype=float) * 1000.0 * molar_mass
                 / (GAS_CONSTANT * t_kelvin))


def absorption_coefficient(frequency, temperature=REFERENCE_TEMPERATURE,
                           humidity=REFERENCE_HUMIDITY, pressure=REFERENCE_PRESSURE):
    """空気吸収によるエネルギー減衰係数 m [1/m]。**ISO 9613-1**。

    `E = E0 * exp(-m * 距離)` の m。dB/m が欲しいときは `absorption_db_per_metre()`。

    **中身の意味**：空気吸収は 2 つの機構の和になっている。

    1. **古典吸収＋回転緩和**（第 1 項）… 粘性と熱伝導による。周波数の 2 乗に比例
    2. **分子の振動緩和**（第 2・3 項）… 酸素と窒素の分子が振動励起されてエネルギーを
       食う。それぞれ緩和周波数 f_rO / f_rN を持ち、**湿度に強く依存する**
       （水分子が緩和を仲介するため）

    実務上の効き方:
      ・低音（〜250 Hz）ではほぼ無視できる
      ・高音（4k〜8k Hz）では 100 m で 10 dB を超えることもある
      ・**乾燥しているほど高音がよく減衰する**（意外に思われるが、緩和周波数が
        可聴域から外れて吸収のピークが可聴域に来るため）
    """
    frequency = np.asarray(frequency, dtype=float)
    t_kelvin = np.asarray(temperature, dtype=float) + KELVIN
    p_ratio = pressure / REFERENCE_PRESSURE
    t_ratio = t_kelvin / ISO_REFERENCE_TEMPERATURE

    # 水蒸気のモル濃度 [%]（ISO 9613-1 の h）
    h = humidity * (saturation_vapour_pressure(temperature, pressure)
                    / REFERENCE_PRESSURE) / p_ratio

    # 酸素・窒素の振動緩和周波数 [Hz]
    f_oxygen = p_ratio * (24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h))
    f_nitrogen = (p_ratio * t_ratio ** -0.5
                  * (9.0 + 280.0 * h * np.exp(-4.170 * (t_ratio ** (-1.0 / 3.0) - 1.0))))

    alpha_db = 8.686 * frequency ** 2 * (
        1.84e-11 / p_ratio * t_ratio ** 0.5
        + t_ratio ** -2.5 * (
            0.01275 * np.exp(-2239.1 / t_kelvin) / (f_oxygen + frequency ** 2 / f_oxygen)
            + 0.1068 * np.exp(-3352.0 / t_kelvin) / (f_nitrogen + frequency ** 2 / f_nitrogen)
        )
    )
    return alpha_db * DB_TO_NEPER


def absorption_db_per_metre(frequency, temperature=REFERENCE_TEMPERATURE,
                            humidity=REFERENCE_HUMIDITY, pressure=REFERENCE_PRESSURE):
    """空気吸収 [dB/m]。ISO 9613-1 の表と直接比べたいとき用。"""
    return absorption_coefficient(frequency, temperature, humidity,
                                  pressure) / DB_TO_NEPER


def legacy_absorption_coefficient(frequency):
    """元コードの空気吸収近似 `m = 1.81e-8 * f^1.57`（20℃・湿度 40% のフィット）。

    ISO 9613-1 との比較用に残してある。**新しいコードでは使わないこと。**
    温度・湿度を変えられないうえ、低音側で ISO より過大に出る。
    """
    return 1.81e-8 * np.asarray(frequency, dtype=float) ** 1.57


class Atmosphere:
    """大気の状態をひとまとめにして持ち回るためのクラス。

    音速と空気吸収は同じ温度・湿度から決まるので、別々に渡すと食い違いが起きる。
    このクラス 1 つを引き回せば整合が保たれる。GUI からは温度・湿度のスライダを
    このクラスに結ぶことになる。

        air = Atmosphere(temperature=25.0, humidity=60.0)
        air.sound_velocity                    # [m/s]
        air.absorption_coefficient(f)         # [1/m]
    """

    def __init__(self, temperature=REFERENCE_TEMPERATURE,
                 humidity=REFERENCE_HUMIDITY, pressure=REFERENCE_PRESSURE):
        if not -50.0 <= temperature <= 60.0:
            raise ValueError(f"温度 {temperature} ℃ は想定範囲外です（-50〜60℃）")
        if not 0.0 <= humidity <= 100.0:
            raise ValueError(f"相対湿度 {humidity} % は 0〜100 の範囲で指定してください")
        if pressure <= 0.0:
            raise ValueError(f"気圧 {pressure} kPa が不正です")
        self.temperature = float(temperature)
        self.humidity = float(humidity)
        self.pressure = float(pressure)

    @property
    def sound_velocity(self):
        return sound_velocity(self.temperature, self.humidity, self.pressure)

    @property
    def density(self):
        """湿り空気の密度 [kg/m³]。音圧レベルの絶対値に要る"""
        return density(self.temperature, self.humidity, self.pressure)

    @property
    def impedance(self):
        """空気の特性インピーダンス ρc [N·s/m³]。20℃ で約 412"""
        return self.density * self.sound_velocity

    def absorption_coefficient(self, frequency):
        return absorption_coefficient(frequency, self.temperature,
                                      self.humidity, self.pressure)

    def absorption_db_per_metre(self, frequency):
        return absorption_db_per_metre(frequency, self.temperature,
                                       self.humidity, self.pressure)

    def replace(self, **kwargs):
        """一部だけ変えた新しい Atmosphere を返す（GUI での値変更用）。"""
        return Atmosphere(kwargs.get("temperature", self.temperature),
                          kwargs.get("humidity", self.humidity),
                          kwargs.get("pressure", self.pressure))

    def summary(self):
        return (f"大気: {self.temperature:.1f}℃ / 湿度 {self.humidity:.0f}% / "
                f"{self.pressure:.3f} kPa → 音速 {self.sound_velocity:.2f} m/s")

    def __repr__(self):
        return (f"Atmosphere(temperature={self.temperature}, "
                f"humidity={self.humidity}, pressure={self.pressure})")
