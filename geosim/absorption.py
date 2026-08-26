"""吸音材のデータと、吸音率の種類の変換。

## 何を解決するモジュールか

1. **吸音率の「種類」を取り違えない**
   カタログに載っている吸音率には少なくとも 2 種類ある。

   | 種類 | 測り方 | 特徴 |
   |---|---|---|
   | **垂直入射吸音率** α_n | 音響管（インピーダンス管） | 正面から当てた場合だけ。1 を超えない |
   | **残響室法吸音率** α_s | 残響室でランダム入射 | **1 を超えることがある**（試料端部の回折等） |

   本 PJ の反射計算（`sound_ray.energy_decay`、書籍 式2.64）が要求するのは
   **垂直入射吸音率**。残響室法の値をそのまま入れると `sqrt(1 - α)` が NaN になる。
   そこで Paris の式を逆に解いて変換する（`random_to_normal`）。

2. **材料を CAD から切り離す**
   CAD のレイヤ名で材料を判別するが、条件を変えるたびに CAD を編集するのは大変。
   そこで「レイヤ → 材料名」の対応を **`LayerAssignment` として外に持つ**。
   GUI からはこの対応表と `MaterialLibrary` を編集すればよく、DXF は触らない。

3. **Excel に無い材料を後から足せる**
   `MaterialLibrary.add()` で追加、`to_csv()` で書き戻せる。

## 使い方

    library = MaterialLibrary.from_csv("absorption.csv", kind="random")
    library.add("自作グラスウール", [0.2, 0.4, 0.7, 0.9, 0.9, 0.85, 0.8, 0.75])
    assignment = LayerAssignment({"1": "コンクリート", "2": "石膏ボード"})
    table = library.absorption_table(assignment, bands)    # → read_dxffile に渡す
"""

import csv
import json
import os
import re

import numpy as np

# 吸音率としてありえない大きさ。これを超える行はヘッダ（列名に中心周波数を書いたもの）
# とみなして読み飛ばす。残響室法でも 1.3 程度までなので 2 で十分区別できる
MAX_PLAUSIBLE_ABSORPTION = 2.0

# オクターブバンドの中心周波数 [Hz]。
# **既定は 8 バンド（63〜8k Hz）**。63 と 8k を対象外にする運用もあるので
# 6 バンド（125〜4k）も使える。バンド数は CSV の列数から自動判別する。
OCTAVE_BAND_FREQUENCIES_8 = np.array([63.0, 125.0, 250.0, 500.0,
                                      1000.0, 2000.0, 4000.0, 8000.0])
OCTAVE_BAND_FREQUENCIES_6 = np.array([125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0])
DEFAULT_OCTAVE_BANDS = OCTAVE_BAND_FREQUENCIES_8

# 吸音率の種類
KIND_NORMAL = "normal"      # 垂直入射吸音率（音響管）
KIND_RANDOM = "random"      # 残響室法吸音率（乱入射）
KINDS = (KIND_NORMAL, KIND_RANDOM)


# 帯域の幅。`1/1` = オクターブ、`1/3` = 1/3 オクターブ（2026-08-26）
BAND_WIDTH_OCTAVE = "1/1"
BAND_WIDTH_THIRD = "1/3"
BAND_WIDTHS = (BAND_WIDTH_OCTAVE, BAND_WIDTH_THIRD)

# 1/3 オクターブの呼び中心周波数（IEC 61260 / JIS C 1513）。
# 計算には 10^(n/10) の厳密値ではなく、**表に載る呼び値**を使う
# （報告書や吸音率の表と突き合わせるときに数字がそろうため）
THIRD_OCTAVE_NOMINAL = (
    12.5, 16.0, 20.0, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0,
    160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0,
    1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0,
    8000.0, 10000.0, 12500.0, 16000.0, 20000.0)


def band_ratio(band_width=BAND_WIDTH_OCTAVE):
    """帯域幅をオクターブ数で返す（`1/1` → 1.0、`1/3` → 1/3）。"""
    text = str(band_width or BAND_WIDTH_OCTAVE).strip()
    if text in ("1/3", "1/3oct", "1/3 oct", "third", "3"):
        return 1.0 / 3.0
    return 1.0


def is_third_octave(band_width):
    return abs(band_ratio(band_width) - 1.0 / 3.0) < 1.0e-9


def band_edges(centres, band_width=BAND_WIDTH_OCTAVE):
    """各バンドの下端・上端 [Hz]。フィルタの遮断周波数に使う。

    幅 1 オクターブなら f/√2 〜 f√2、1/3 なら f·2^(-1/6) 〜 f·2^(1/6)。
    """
    centres = np.asarray(centres, dtype=float)
    half = band_ratio(band_width) / 2.0
    return centres * 2.0 ** (-half), centres * 2.0 ** half


def frequency_bands(band_number, band_width=BAND_WIDTH_OCTAVE, start=None):
    """計算に使う帯域の中心周波数。

    引数:
        band_number  バンド数
        band_width   `1/1`（オクターブ。既定）か `1/3`
        start        いちばん低いバンドの中心周波数 [Hz]。
                     None ならオクターブは従来どおり（8→63、6→125）、
                     1/3 は 100 Hz から

    ★1/3 のときは**呼び値の表**（`THIRD_OCTAVE_NOMINAL`）から連続して取る。
    """
    band_number = int(band_number)
    if not is_third_octave(band_width):
        if start is None:
            return octave_bands(band_number)
        return float(start) * 2.0 ** np.arange(band_number, dtype=float)

    table = np.array(THIRD_OCTAVE_NOMINAL, dtype=float)
    first = 100.0 if start is None else float(start)
    index = int(np.argmin(np.abs(np.log2(table / first))))
    chosen = table[index:index + band_number]
    if len(chosen) < band_number:
        # 表の外まで欲しいときは 2^(1/3) ずつ伸ばす（20 kHz より上など）
        extra = chosen[-1] * 2.0 ** (np.arange(1, band_number - len(chosen) + 1)
                                     / 3.0)
        chosen = np.concatenate([chosen, extra])
    return chosen


def octave_bands(band_number):
    """バンド数から中心周波数の並びを返す。

    8 → 63〜8k、6 → 125〜4k。それ以外は 1000 Hz を含む形で前後に伸ばす。
    """
    if band_number == 8:
        return OCTAVE_BAND_FREQUENCIES_8.copy()
    if band_number == 6:
        return OCTAVE_BAND_FREQUENCIES_6.copy()
    # 汎用: 1000 Hz を基準にオクターブで並べる（1000 Hz が必ず入るように中央寄せ）
    below = (band_number - 1) // 2
    return 1000.0 * 2.0 ** np.arange(-below, band_number - below, dtype=float)


# ------------------------------------------------------------------------------
# 吸音率の種類の変換（Paris の式）
# ------------------------------------------------------------------------------

def reflection_coefficient(impedance, cos_theta):
    """局所反応性壁面の圧力反射率 R(θ)。

        R = (z cosθ - 1) / (z cosθ + 1)

    z は規格化音響インピーダンス（実数と仮定）。書籍 2.2 節と同じモデル。
    """
    z_cos = np.asarray(impedance, dtype=float) * np.asarray(cos_theta, dtype=float)
    return (z_cos - 1.0) / (z_cos + 1.0)


def normal_absorption(impedance):
    """垂直入射吸音率 α_n = 1 - |R(0)|^2 = 4z / (1+z)^2。

    ※ z と 1/z が同じ値を与える（z → 1/z について対称）。
      つまり α_n だけからは z が一意に決まらない。逆変換に注意が要る理由。
    """
    z = np.asarray(impedance, dtype=float)
    return 4.0 * z / (1.0 + z) ** 2


def impedance_from_normal(alpha_normal):
    """垂直入射吸音率から規格化インピーダンス z を求める（z >= 1 の枝）。

        |R| = sqrt(1 - α_n),  z = (1 + |R|) / (1 - |R|)

    `sound_ray.energy_decay` が内部でやっているのと同じ計算。
    α_n = 1 のとき z = 1（完全吸音）。
    """
    alpha = np.asarray(alpha_normal, dtype=float)
    if np.any(alpha > 1.0) or np.any(alpha < 0.0):
        raise ValueError("垂直入射吸音率は 0〜1 の範囲でなければなりません")
    r = np.sqrt(1.0 - alpha)
    return (1.0 + r) / np.maximum(1.0 - r, 1e-15)


def statistical_absorption(impedance, method="closed", samples=20001):
    """Paris の式による乱入射（残響室法）吸音率 α_s。

        α_s = 2 ∫₀^{π/2} α(θ) sinθ cosθ dθ

    **意味**：あらゆる方向から等確率に音が来る（拡散音場）としたときの平均吸音率。
    重み `2 sinθ cosθ` は「立体角の重み sinθ」と「面に投影したときの割合 cosθ」の積で、
    ∫ の全体が 1 になるよう規格化されている。

    解析解もある。μ = cosθ と置くと 1 - R² = 4zμ/(zμ+1)² なので

        α_s = 8z ∫₀¹ μ²/(zμ+1)² dμ = (8/z²) [ (1+z) - 2 ln(1+z) - 1/(1+z) ]

    method='closed'（既定）はこの解析解、'numeric' は数値積分。
    最初は式の取り違えを避けるため数値積分だけで実装していたが、
    両者の一致を確認できた（tests/test_geosim.py）ので既定を解析解にした。
    逆解き（`impedance_from_statistical`）が二分法で何度も呼ぶため、
    数値積分のままだと重すぎる。
    """
    if method == "closed":
        return statistical_absorption_closed_form(impedance)
    if method != "numeric":
        raise ValueError(f"method は 'closed' か 'numeric' です: {method!r}")

    z = np.atleast_1d(np.asarray(impedance, dtype=float))
    theta = np.linspace(0.0, np.pi / 2, samples)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    weight = 2.0 * sin_t * cos_t

    alpha_theta = 1.0 - reflection_coefficient(z[:, None], cos_t[None, :]) ** 2
    result = np.trapezoid(alpha_theta * weight[None, :], theta, axis=1)
    return result if np.ndim(impedance) else float(result[0])


def statistical_absorption_closed_form(impedance):
    """Paris の式の解析解。`statistical_absorption` の検算用。

        α_s = (8/z²) [ (1+z) - 2 ln(1+z) - 1/(1+z) ]

    導出: μ = cosθ と置くと 1 - R² = 4zμ/(zμ+1)² なので
        α_s = 2∫₀¹ (1-R²) μ dμ = 8z ∫₀¹ μ²/(zμ+1)² dμ
    さらに u = zμ+1 と置換すると
        = (8/z²) [ u - 2 ln u - 1/u ]₁^{1+z}
    となり上式を得る。
    """
    z = np.asarray(impedance, dtype=float)
    return (8.0 / z ** 2) * ((1.0 + z) - 2.0 * np.log(1.0 + z) - 1.0 / (1.0 + z))


def _statistical_maximum(samples=200001):
    """α_s が取りうる最大値と、そのときの z を返す。

    α_s は z について単調ではなく、ある z で最大になってから両側で 0 に近づく。
    **この最大値を超える残響室法吸音率は、局所反応モデルでは表現できない。**
    """
    z_grid = np.geomspace(0.1, 10.0, samples)
    values = statistical_absorption(z_grid)
    index = int(np.argmax(values))
    return float(values[index]), float(z_grid[index])


STATISTICAL_MAX, STATISTICAL_MAX_IMPEDANCE = _statistical_maximum()


def impedance_from_statistical(alpha_random, branch="hard", tolerance=1e-12):
    """残響室法吸音率 α_s から規格化インピーダンス z を求める（Paris の式の逆解き）。

    **なぜ単純に解けないのか**：α_s は z について単調ではない。
    z → ∞（硬い壁）でも z → 0（完全に軟らかい壁）でも α_s → 0 になり、
    その間（z ≈ 1.567）で最大値 **0.951** を取る（局所反応・実数インピーダンスの古典的な結果）。
    したがって 1 つの α_s に対して z は 2 つある。

    - `branch='hard'`（既定）… z >= z_max の枝。建築材料はこちら
      （多孔質吸音材や剛壁は空気より音響インピーダンスが高い）
    - `branch='soft'` … z <= z_max の枝

    α_s が最大値を超えている場合（残響室法では 1 を超えることすらある）は、
    **最大値に丸めて警告を出す**。局所反応モデルで表せる限界だからである。
    """
    alpha = np.atleast_1d(np.asarray(alpha_random, dtype=float))
    if branch not in ("hard", "soft"):
        raise ValueError("branch は 'hard' か 'soft' を指定してください")

    clipped = np.minimum(alpha, STATISTICAL_MAX)
    result = np.empty_like(clipped)

    for i, target in enumerate(clipped):
        if target <= 0.0:
            result[i] = np.inf if branch == "hard" else 0.0
            continue
        # 単調な区間で二分法。hard 側は z が大きいほど α_s が小さい
        if branch == "hard":
            low, high = STATISTICAL_MAX_IMPEDANCE, 1e6
            decreasing = True
        else:
            low, high = 1e-6, STATISTICAL_MAX_IMPEDANCE
            decreasing = False
        for _ in range(200):
            mid = 0.5 * (low + high)
            value = statistical_absorption(mid)
            if (value > target) == decreasing:
                low = mid
            else:
                high = mid
            if high - low < tolerance * max(1.0, high):
                break
        result[i] = 0.5 * (low + high)

    return result if np.ndim(alpha_random) else float(result[0])


def random_to_normal(alpha_random, branch="hard", warn=True, label=""):
    """残響室法吸音率 → 垂直入射吸音率。

    Paris の式を逆に解いてインピーダンス z を求め、そこから垂直入射吸音率を出す。
    本 PJ の反射計算（書籍 式2.64）はこの値を要求する。

    残響室法吸音率が上限 0.951 を超える場合は丸めて警告する
    （測定値が 1 を超えるのは試料端部の回折などによるもので、
      局所反応モデルではそのまま扱えない）。

    変換の目安:

    | 残響室法 α_s | z | 垂直入射 α_n |
    |---|---|---|
    | 0.05 | 150.4 | 0.026 |
    | 0.20 | 32.6 | 0.116 |
    | 0.50 | 9.66 | 0.340 |
    | 0.80 | 3.88 | 0.652 |
    | 0.95 | 1.69 | 0.934 |

    **残響室法の値をそのまま垂直入射として使うと吸音を大幅に過大評価する**ことが分かる。
    """
    alpha = np.atleast_1d(np.asarray(alpha_random, dtype=float))
    over = alpha > STATISTICAL_MAX
    if warn and np.any(over):
        where = f"（{label}）" if label else ""
        print(f"[absorption] 警告{where}: 残響室法吸音率 "
              f"{np.array2string(alpha[over], precision=3)} が局所反応モデルの上限 "
              f"{STATISTICAL_MAX:.3f} を超えています。上限に丸めます。"
              f"（残響室法の値が 1 を超えるのは試料端部の回折などによるものです）")
    z = impedance_from_statistical(alpha, branch=branch)
    result = normal_absorption(z)
    return result if np.ndim(alpha_random) else float(result[0])


def normal_to_random(alpha_normal):
    """垂直入射吸音率 → 残響室法吸音率（Paris の式の順方向）。

    変換が往復で戻ることの確認や、カタログ値との突き合わせに使う。
    """
    return statistical_absorption(impedance_from_normal(alpha_normal))


# ------------------------------------------------------------------------------
# 材料
# ------------------------------------------------------------------------------

class Material:
    """吸音材 1 種類ぶんのデータ。

    属性:
        name         材料名
        coefficients (nband,) 吸音率。**kind が示す種類のままの生値**
        kind         'normal'（垂直入射）または 'random'（残響室法）
        note         出典や備考（GUI で表示する想定）
    """

    def __init__(self, name, coefficients, kind=KIND_NORMAL, note=""):
        if kind not in KINDS:
            raise ValueError(f"kind は {KINDS} のいずれかです: {kind!r}")
        self.name = str(name)
        self.coefficients = np.asarray(coefficients, dtype=float).ravel()
        if self.coefficients.size == 0:
            raise ValueError(f"材料 {name!r} の吸音率が空です")
        if np.any(self.coefficients < 0.0):
            raise ValueError(f"材料 {name!r} に負の吸音率があります")
        self.kind = kind
        self.note = str(note)

    @property
    def band_number(self):
        return int(self.coefficients.size)

    def normal_incidence(self, warn=True):
        """垂直入射吸音率にそろえて返す。反射計算に渡すのはこちら。"""
        if self.kind == KIND_NORMAL:
            if np.any(self.coefficients > 1.0):
                raise ValueError(
                    f"材料 {self.name!r} は垂直入射吸音率として登録されていますが "
                    f"1 を超える値があります。残響室法（kind='random'）ではありませんか？")
            return self.coefficients.copy()
        return random_to_normal(self.coefficients, warn=warn, label=self.name)

    def resample(self, band_number):
        """バンド数を合わせる（TODO E-14。2026-08-15 に対数軸の補間へ変更）。

        8 バンド（63〜8k Hz）と 6 バンド（125〜4k Hz）を混ぜて使えるようにするためのもの。

        **周波数軸（log₂ f）の上で線形補間する。** オクターブバンドは中心周波数が
        2 倍ずつ並ぶので、対数を取ると等間隔になり、素直に内挿できる。

        以前は「端の値を複製／切り詰める」だけだった。中心周波数が一致する場合
        （6↔8 バンドはまさにこれ）は結果が同じだが、素性の違う表
        （1/3 オクターブ表など）を混ぜると値がずれる。

        ⚠ **範囲の外は外挿せず、端の値をそのまま使う**（一定外挿）。
        吸音率は物理的に 0〜1 に収まるので、傾きを延長すると簡単に範囲外へ出てしまう。
        6 バンド表を 8 バンドで使えば、**63 Hz は 125 Hz のコピー、
        8 kHz は 4 kHz のコピー**になる。これは補間ではなく「その帯域のデータが無い」
        という事実の表れなので、値を鵜呑みにしないこと。
        """
        source = np.asarray(octave_bands(self.band_number), dtype=float)
        target = np.asarray(octave_bands(band_number), dtype=float)
        if len(source) == len(target) and np.allclose(source, target):
            return self.coefficients.copy()
        return np.interp(np.log2(target), np.log2(source), self.coefficients)

    def to_dict(self):
        return {"name": self.name, "coefficients": self.coefficients.tolist(),
                "kind": self.kind, "note": self.note}

    def __repr__(self):
        return (f"Material({self.name!r}, {np.array2string(self.coefficients, precision=3)}, "
                f"kind={self.kind!r})")


class MaterialLibrary:
    """材料の一覧。GUI から追加・変更・削除できるようにするための入れ物。

    キーは材料名。同じ材料に複数の別名（元コードの ID など）を付けられる。
    """

    def __init__(self, materials=None, aliases=None):
        self.materials = dict(materials or {})
        self.aliases = dict(aliases or {})    # 別名 -> 材料名

    # ---- 参照 ----
    def __contains__(self, key):
        return key in self.materials or key in self.aliases

    def __len__(self):
        return len(self.materials)

    def get(self, key, default=None):
        if key in self.materials:
            return self.materials[key]
        if key in self.aliases:
            return self.materials.get(self.aliases[key], default)
        return default

    def names(self):
        return sorted(self.materials)

    # ---- 編集（GUI 用） ----
    def add(self, name, coefficients, kind=KIND_NORMAL, note="", overwrite=False):
        """材料を追加する。Excel に無い材料を GUI から足すための入口。"""
        if name in self.materials and not overwrite:
            raise KeyError(f"材料 {name!r} は既にあります（overwrite=True で上書き）")
        self.materials[name] = Material(name, coefficients, kind, note)
        return self.materials[name]

    def update(self, name, coefficients=None, kind=None, note=None):
        """既存の材料を書き換える。"""
        if name not in self.materials:
            raise KeyError(f"材料 {name!r} がありません")
        current = self.materials[name]
        self.materials[name] = Material(
            name,
            current.coefficients if coefficients is None else coefficients,
            current.kind if kind is None else kind,
            current.note if note is None else note)
        return self.materials[name]

    def remove(self, name):
        self.materials.pop(name, None)
        for alias, target in list(self.aliases.items()):
            if target == name:
                self.aliases.pop(alias)

    def add_alias(self, alias, name):
        if name not in self.materials:
            raise KeyError(f"材料 {name!r} がありません")
        self.aliases[alias] = name

    # ---- 入出力 ----
    @classmethod
    def from_file(cls, file_name, kind=None, band_number=None, verbose=True):
        """**拡張子で読み方を選ぶ**（2026-08-21 ユーザー要望）。

        - `.xlsx` / `.xlsm` … 「吸音率」シートを読む（`condition_table` の形）。
          用意した吸音率表を設定画面で選べるようにするため
        - それ以外 … CSV として読む（`from_csv`）

        xlsx に「吸音率」シートが無ければ `ValueError`。
        """
        if str(file_name or "").lower().endswith((".xlsx", ".xlsm")):
            import condition_table as ct
            library = ct.library_from_book(file_name, kind=kind, verbose=verbose)
            if library is None:
                raise ValueError(
                    f"{os.path.basename(file_name)} に「{ct.ABSORPTION_SHEET}」"
                    f"シートが見つかりません（番号・材料名・吸音率の並びが要ります）")
            return library
        return cls.from_csv(file_name, kind=kind, band_number=band_number)

    @classmethod
    def from_csv(cls, file_name, kind=None, band_number=None):
        """吸音率 CSV を読み込む。

        対応する形式:

          (A) 元コードの absorption.csv 形式（ID + 材料名 + 吸音率）
                1,Concrete wall,0.01,0.02,...
              → ID も別名として登録するので、DXF のレイヤ名がどちらでも引ける
          (B) 材料名 + 吸音率
                コンクリート,0.02,0.02,...
          (C) 上記の後ろに種類や備考の列が付いたもの（`to_csv()` が書く形）
                コンクリート,0.02,0.02,...,normal,備考
              → **数値でないセルが出たらそこまでを吸音率として読む**

        **吸音率の種類**は次の順で決まる。

          1. 引数 `kind`
          2. ファイル内の宣言行 `# kind: random`
          3. 既定 'normal'（垂直入射）＋警告

        ★★ 種類が分からない表をそのまま使うのは危険です。★★
        残響室法の値を垂直入射として扱うと吸音を過大評価します（1 を超えると計算が壊れます）。
        カタログの出典を確認して `kind` を明示してください。
        """
        text = _read_text(file_name)

        # `# kind: random` のような宣言行を探す。
        # 説明文の中に kind という語があっても拾わないよう、書式を厳しく見る
        declared_kind = None
        for line in text.splitlines():
            match = re.match(r"^\s*#\s*kind\s*[:：]\s*(\w+)\s*$", line)
            if match and match.group(1).lower() in KINDS:
                declared_kind = match.group(1).lower()
                break

        resolved = kind or declared_kind
        if resolved is None:
            resolved = KIND_NORMAL
            print(f"[absorption] 警告: {os.path.basename(file_name)} の吸音率が "
                  f"垂直入射か残響室法か分かりません。垂直入射として扱います。"
                  f"（ファイル冒頭に `# kind: random` と書くか、kind= を指定してください）")

        library = cls()
        for row in csv.reader(text.splitlines()):
            if not row or not row[0].strip() or row[0].lstrip().startswith("#"):
                continue
            cells = [c.strip() for c in row]

            second_is_number = False
            if len(cells) >= 2:
                try:
                    float(cells[1])
                    second_is_number = True
                except ValueError:
                    second_is_number = False

            if second_is_number:
                keys, raw = [cells[0]], cells[1:]
            elif len(cells) >= 3:
                keys, raw = [cells[0], cells[1]], cells[2:]
            else:
                continue

            # ★数値でないセルが出たら**そこで打ち切る**（行は捨てない）。
            #   `to_csv()` は吸音率のあとに `kind` と `note` の列を書くので、
            #   捨てる作りだと**自分が書いた CSV を読み戻せなかった**
            #   （GUI で材料を足して保存 → 読み直しで消える。2026-08-21 に発見）
            values = []
            for cell in raw:
                if cell == "":
                    break
                try:
                    values.append(float(cell))
                except ValueError:
                    break
            if len(values) < 3:
                continue      # ヘッダ行や説明行（バンドは 6 か 8 なので 3 未満は無い）
            # ヘッダ行の見出しが数値（`material,63,125,250,...` のように
            # 中心周波数を列名にしている場合）だと材料として取り込まれてしまう。
            # 吸音率が 2 を超えることはないので、それで弾く
            if min(values) > MAX_PLAUSIBLE_ABSORPTION:
                continue
            if band_number is not None:
                values = (values + [values[-1]] * band_number)[:band_number]

            name = keys[-1]
            library.add(name, values, kind=resolved, overwrite=True)
            for alias in keys[:-1]:
                if alias:
                    library.aliases[alias] = name

        return library

    def to_csv(self, file_name, frequencies=None):
        """材料一覧を CSV に書き出す（GUI で追加したものを保存する用）。"""
        names = self.names()
        if not names:
            raise ValueError("材料が 1 つもありません")
        band_number = self.materials[names[0]].band_number
        if frequencies is None:
            frequencies = octave_bands(band_number)

        kinds = {self.materials[n].kind for n in names}
        with open(file_name, "w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            if len(kinds) == 1:
                f.write(f"# kind: {kinds.pop()}\n")
            writer.writerow(["material"] + [f"{v:g}" for v in frequencies]
                            + ["kind", "note"])
            for name in names:
                m = self.materials[name]
                writer.writerow([name] + [f"{c:g}" for c in m.coefficients]
                                + [m.kind, m.note])
        return file_name

    def to_json(self, file_name):
        payload = {"materials": [m.to_dict() for m in self.materials.values()],
                   "aliases": self.aliases}
        with open(file_name, "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)
        return file_name

    @classmethod
    def from_json(cls, file_name):
        with open(file_name, encoding="utf-8") as f:
            payload = json.load(f)
        library = cls(aliases=payload.get("aliases", {}))
        for item in payload.get("materials", []):
            library.add(item["name"], item["coefficients"],
                        item.get("kind", KIND_NORMAL), item.get("note", ""),
                        overwrite=True)
        return library

    # ---- 計算側へ渡す ----
    def absorption_table(self, assignment=None, band_number=None, warn=True):
        """`read_dxffile.read_model(absorption_table=...)` に渡す辞書を作る。

        戻り値は {レイヤ名またはキー: ndarray(band_number,)}。
        **値は必ず垂直入射吸音率に変換済み**なので、下流は種類を気にしなくてよい。

        assignment : LayerAssignment | dict | None
            レイヤ名 → 材料名の対応。None なら「レイヤ名 = 材料名」とみなす
            （＝これまでの挙動）。
        """
        if band_number is None:
            band_number = max((m.band_number for m in self.materials.values()),
                              default=len(DEFAULT_OCTAVE_BANDS))

        resampled = {m.band_number for m in self.materials.values()} - {band_number}
        if resampled and warn:
            print(f"[absorption] 注意: 材料表は {sorted(resampled)} バンドですが "
                  f"{band_number} バンドで計算します。"
                  f"{'両端を複製して伸ばします' if min(resampled) < band_number else '端を切り詰めます'}"
                  f"（8 バンドなら 63 Hz と 8 kHz が端の値のコピーになります）")

        table = {}
        for name, material in self.materials.items():
            values = Material(name, material.resample(band_number),
                              material.kind, material.note).normal_incidence(warn=warn)
            table[name] = values
        for alias, target in self.aliases.items():
            if target in table:
                table[alias] = table[target]

        if assignment is not None:
            mapping = (assignment.mapping if isinstance(assignment, LayerAssignment)
                       else dict(assignment))
            for layer, material_name in mapping.items():
                if material_name in table:
                    table[layer] = table[material_name]
                elif warn:
                    print(f"[absorption] 警告: レイヤ {layer!r} に割り当てられた材料 "
                          f"{material_name!r} が材料一覧にありません")
        return table

    def summary(self):
        kinds = {}
        for m in self.materials.values():
            kinds[m.kind] = kinds.get(m.kind, 0) + 1
        detail = " / ".join(f"{k}: {v}種" for k, v in sorted(kinds.items()))
        return f"材料 {len(self.materials)} 種（別名 {len(self.aliases)} 個） {detail}"


class LayerAssignment:
    """CAD のレイヤ名 → 材料名の対応表。

    **CAD を編集せずに材料を差し替えるための仕組み。**
    条件を変えるたびに DXF を書き直すのは大変なので、対応表だけを外に持つ。
    GUI からはこのオブジェクトを編集して保存/読み込みする。

        assignment = LayerAssignment({"壁": "石膏ボード"})
        assignment.assign("床", "カーペット")
        assignment.save("条件A.json")
    """

    def __init__(self, mapping=None):
        self.mapping = dict(mapping or {})

    def assign(self, layer, material_name):
        self.mapping[layer] = material_name

    def unassign(self, layer):
        self.mapping.pop(layer, None)

    def get(self, layer, default=None):
        return self.mapping.get(layer, default)

    def missing(self, layers, library):
        """割り当てが無い、または材料一覧に無いレイヤを返す（GUI の警告用）。"""
        result = []
        for layer in layers:
            name = self.mapping.get(layer, layer)
            if name not in library:
                result.append(layer)
        return result

    def save(self, file_name):
        with open(file_name, "w", encoding="utf-8") as f:
            json.dump(self.mapping, f, ensure_ascii=False, indent=2)
        return file_name

    @classmethod
    def load(cls, file_name):
        with open(file_name, encoding="utf-8") as f:
            return cls(json.load(f))

    def __repr__(self):
        return f"LayerAssignment({self.mapping!r})"


def _read_text(file_name):
    """元コード付属の absorption.csv は CP932 なので、順に試してデコードする。"""
    last_error = None
    for encoding in ("utf-8-sig", "cp932", "latin-1"):
        try:
            with open(file_name, encoding=encoding, newline="") as f:
                return f.read()
        except UnicodeDecodeError as e:
            last_error = e
    raise last_error
