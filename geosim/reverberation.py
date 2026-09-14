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
    3. **-5〜-35 dB の区間に直線を最小二乗で当て、傾きから T30 を外挿する**（ISO 3382）

## 読み取り方（2026-09-05 に ISO 準拠へ変更）

★**評価区間の全サンプルへの最小二乗回帰**が ISO 3382 の規定。
それ以前は「-5 dB と -35 dB を横切る**2 点**の時刻差」で出していた（元コード 134 行）。

- 減衰がまっすぐなら両者はほぼ一致する（合成した指数減衰で差 1% 未満）。
  **これまでの数字が的外れだったという意味ではない**
- 減衰が曲がっていると大きくずれる。`test.dxf` の実測で **125 Hz で 35.7%**
- 回帰にすると**適合の良し悪しが数値で出る**のが実利。
  ISO 3382-2 の非線形性 `ξ = 1000(1-r²)` が同時に手に入り、
  「このバンドは読んではいけない」を人の目でなく数字で言える

2 点法は `fit='crossing'` で残してある（元 Fortran との一致確認に使う）。
"""

import numpy as np
from scipy.signal import butter, fftconvolve, firwin, sosfilt

import absorption as ab
import table as tb
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

# 減衰曲線の読み取り方。
#   'least_squares' … 評価区間の全サンプルへ直線を最小二乗で当てる（**ISO 3382**。既定）
#   'crossing'      … 開始 dB と終了 dB を横切る 2 点の時刻差（元コード 134 行）。
#                     参照実装として残す。ISO 準拠の数字が要るときは使わないこと
DECAY_FIT_LEAST_SQUARES = "least_squares"
DECAY_FIT_CROSSING = "crossing"
DECAY_FITS = (DECAY_FIT_LEAST_SQUARES, DECAY_FIT_CROSSING)
DEFAULT_DECAY_FIT = DECAY_FIT_LEAST_SQUARES

# ISO 3382 の曲率 C を出すときの評価区間（T20 と T30）。
# ★C は「どの区間で残響時間を読んだか」に依らない**減衰曲線そのものの性質**なので、
#   db_max / db_min を変えてもこの 2 つで測る
CURVATURE_RANGES = (DECAY_MEASURES["T20"], DECAY_MEASURES["T30"])

# ISO 3382 が求める余裕。評価区間の下端より、さらにこれだけ下まで
# 減衰が見えていること（T30 なら -45 dB まで）。
# 実測では暗騒音の話だが、シミュレーションでは**反射回数の打ち切り**と
# **インパルス応答長の打ち切り**が同じ役目の「床」を作る
DECAY_MARGIN_DB = 10.0

# 曲率 C・非線形性 ξ の「ここから先は信用しない」の目安。
#   C  … ISO 3382-1 の曲率 [%]
#   ξ  … ISO 3382-2 の非線形性 1000(1-r²)。精密級 ξ<=5 / 工学級 ξ<=10 が目安
CURVATURE_LIMIT_PERCENT = 10.0
NONLINEARITY_LIMIT = 10.0

# IEC 61260 のオクターブバンドに合わせた Butterworth の次数（6 次）。
# ★「クラス 1」とは名乗らない：クラスは次数ではなく**減衰量マスクへの適合**で
#   決まるもので、その試験はこのリポジトリでは行っていない（2026-09-05）
FILTER_ORDER = 6

# method='fir' のときのタップ数
FIR_NUMTAPS = 4096


# ------------------------------------------------------------------------------
# 統計的残響式（Sabine / Eyring / Eyring-Knudsen）
#
# 音線を飛ばさず、**室容積 V と各面の面積・吸音率だけ**から残響時間を見積もる。
# 拡散音場（音がどの方向からも等確率に来る）を前提にした古典的な式で、
# シミュレーション結果の妥当性を確かめる物差しになる。
#
# ★閉じた室が前提★ 開いた形状（一面だけの壁など）では容積が定義できないので使えない。
# ------------------------------------------------------------------------------

def print_frequency_table(frequencies, rows, prefix="[reverberation]", width=10,
                          label_width=14):
    """**周波数を横に並べて**画面に出す（table.py の共通ルール）。

    CSV・図・画面で向きが揃っていないと、どこかで読み替えが要る。
    """
    print(f"{prefix} {'':>{label_width}}"
          + "".join(f"{f:>{width}.0f}" for f in frequencies))
    for name, values in rows:
        cells = "".join("       ---" if np.isnan(v) else f"{v:{width}.3f}"
                        for v in np.asarray(values, dtype=float))
        print(f"{prefix} {name:>{label_width}}{cells}")


def triangle_areas(mesh):
    """各三角形の面積 [m^2]。外積の大きさの半分。"""
    if not len(mesh):
        return np.zeros(0)      # 面が無くても落ちないように（呼ぶ側で判断する）
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
    """Sabine / Eyring / Eyring-Knudsen の残響式で残響時間を見積もる。

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
    | Sabine | `S·ᾱ` | 最も古典的。**吸音率が小さいとき（〜0.2）に妥当**。ᾱ→1 でも T が 0 にならない欠点 |
    | Eyring | `-S·ln(1-ᾱ)` | 反射のたびに (1-ᾱ) 倍になると考える。**吸音率が大きいときはこちら** |
    | Eyring-Knudsen | `-S·ln(1-ᾱ) + 4mV` | Eyring に**空気吸収**を足したもの。高域で効く |

    ᾱ は面積で重み付けした平均吸音率 `Σ Sᵢαᵢ / S`。

    ★**名前について**（ユーザー指摘 2026-08-17）。アイリングの式そのものは
      `-S ln(1-ᾱ)` までで、**空気吸収の項 `4mV` を足した形はヌードセンの寄与**。
      なので 2 つを別の列として並べる。差がそのまま**空気吸収の効き**になる。

    ★**ミリントンの式は落とした**（2026-08-17 ユーザー判断）。
      `ミリントン・セッテの式`（Millington 1932 / Sette 1933）は面ごとに対数を取る
      `-Σ Sᵢ ln(1-αᵢ)` で、**αᵢ→1 の面が 1 枚でもあると A が発散して T→0** になる
      （開口が 1 つあるだけで残響ゼロという結論になる）。実務では使われないため。
      研修室（吸音面 α=0.951）では Eyring-Knudsen の半分近い値が出ていた。

    戻り値: dict
        'frequencies' / 'volume' / 'total_area' / 'mean_free_path'
        'mean_absorption' (nf,) / 'equivalent_area' (nf,) / 'air_absorption_area' (nf,)
        'sabine' / 'eyring' / 'eyring_knudsen' 各 (nf,) [s]
        'surface'  … `surface_summary()` の結果
    """
    if atmosphere is None:
        atmosphere = Atmosphere()
    if not len(mesh):
        if verbose:
            print("[reverberation] 面が 1 枚も無いので統計残響式は計算できません")
        return None
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

    # ln(1-ᾱ) は ᾱ → 1 で発散する。1 に張り付いた材料ばかりだと Eyring は
    # 意味を持たなくなるので、そこは NaN にして知らせる
    with np.errstate(divide="ignore", invalid="ignore"):
        eyring_area = -total_area * np.log(np.where(mean_absorption < 1.0,
                                                    1.0 - mean_absorption, np.nan))

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
        # **Sabine と Eyring は空気吸収を入れない素の形**。
        # 空気吸収を足したのが Eyring-Knudsen で、差がそのまま空気吸収の効きになる
        "sabine": to_time(total_area * mean_absorption),
        "eyring": to_time(eyring_area),
        "eyring_knudsen": to_time(eyring_area + air),
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
        # 周波数を**横**に並べる（table.py の共通ルール）。
        # 画面の表とグラフと CSV で向きを揃えておく
        print(f"[統計残響] {'':>16}" + "".join(f"{f:>10.0f}" for f in frequencies))
        for key, label in ([("mean_absorption", "平均α")]
                           + list(STATISTICAL_LABELS.items())):
            def cell(value):
                return "     ---" if np.isnan(value) else f"{value:10.3f}"
            print(f"[統計残響] {label:>16}"
                  + "".join(cell(v) for v in result[key]))

    return result


# 統計残響式の表示名。**空気吸収の有無で名前が変わる**（ユーザー指摘 2026-08-17）。
# アイリングの式そのものは -S ln(1-ᾱ) までで、4mV を足した形はヌードセンの寄与
STATISTICAL_LABELS = {
    "sabine": "Sabine",
    "eyring": "Eyring",
    "eyring_knudsen": "Eyring-Knudsen",
}


def statistical_labels(result=None):
    """統計残響式の表示名 {キー: 名前}。`result` は互換のため受けるだけ。"""
    return dict(STATISTICAL_LABELS)


def statistical_reverberation_from_model(model, **kwargs):
    """`read_dxffile.read_model()` の結果からそのまま統計残響式を計算する。

    容積が決まらない形状（囲まれていない）では None を返して警告する。

    ★判定は `model.volume` があるかで行う。**`is_closed` で見てはいけない**。
    `is_closed` は「辺が 1 対 1 で閉じているか」なので、壁を高さの帯や開口で
    分割したモデル（T 字接合）では False になるが、面としては閉じていて容積は出せる。
    実際それで統計残響式が使えず、値が出ないモデルがあった（2026-08-19）。
    容積が怪しいかどうかは `model.free_edges` / `volume_note` が知らせる。
    """
    volume = getattr(model, "volume", None)
    if not volume:
        print(f"[統計残響] 容積が決まらないので統計残響式は使えません"
              f"（{getattr(model, 'volume_note', None) or '囲まれていない形状'}）")
        return None
    if getattr(model, "free_edges", None):
        print(f"[統計残響] 注意: 自由端が {len(model.free_edges)} 本あるので"
              f"容積 {abs(volume):.2f} m3 は目安です")
    return statistical_reverberation(model.mesh, abs(volume), **kwargs)


def write_statistical_reverberation(filename, result):
    """統計残響式の結果を CSV に保存する。**周波数は横**（table.py の共通ルール）。"""
    return tb.write_frequency_table(filename, result["frequencies"], {
        "mean_absorption": result["mean_absorption"],
        "equivalent_area_m2": result["equivalent_area"],
        "air_absorption_area_m2": result["air_absorption_area"],
        "sabine_s": result["sabine"],
        "eyring_s": result["eyring"],
        "eyring_knudsen_s": result["eyring_knudsen"],
    })


def schroeder_integral(x):
    """後ろ向きに積分した残留エネルギー。元コード 114〜117 行。

        D(t) = Σ_{τ >= t} x(τ)^2

    「時刻 t 以降にまだ残っているエネルギー」。
    インパルス応答の 2 乗をそのまま dB 表示すると反射音の干渉で激しく暴れるが、
    積分すると滑らかになる（Schroeder の逆積分法）。
    """
    return np.cumsum(x[::-1] ** 2)[::-1]


def octave_bandpass(signal, centre_frequency, sampling_frequency,
                    method="butter", order=FILTER_ORDER, numtaps=FIR_NUMTAPS,
                    band_width=ab.BAND_WIDTH_OCTAVE):
    """1 バンドぶんを切り出す（既定はオクターブ＝中心の 1/√2 〜 √2 倍）。

    band_width … `1/1`（オクターブ）か `1/3`。★2026-08-26 に 1/3 対応で追加。

    method='butter' … `scipy.signal.butter` の SOS を `sosfilt` で適用（既定）。
                      IEC 61260 のオクターブフィルタは Butterworth 系なので実務に沿う。
    method='fir'    … `scipy.signal.firwin` + 線形畳み込み。遅れを切り落として返す。
    """
    nyquist = sampling_frequency / 2.0
    # ★帯域幅は `band_width` で決める（`1/1` なら f/√2〜f√2、`1/3` なら
    #   f·2^(∓1/6)）。2026-08-26 に 1/3 オクターブ対応で追加
    half = ab.band_ratio(band_width) / 2.0
    lower = centre_frequency * 2.0 ** (-half)
    upper = min(centre_frequency * 2.0 ** half, nyquist * 0.999)
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


def _decay_time(decay_db, dt, db_start, db_end, fit=DEFAULT_DECAY_FIT,
                detail=False):
    """減衰曲線から 60 dB 減衰時間を求める。**既定は ISO 3382 の最小二乗回帰。**

    引数:
        decay_db  (n,) Schroeder 積分を dB にしたもの（先頭が 0 dB）
        dt        サンプル間隔 [s]
        db_start / db_end  評価区間 [dB]（T30 なら -5 / -35）
        fit       'least_squares'（ISO 3382。既定）| 'crossing'（2 点法。元コード 134 行）
        detail    True なら (値, 適合の情報) を返す

    戻り値:
        60 dB 減衰時間 [s]。求まらなければ `np.nan`
        （評価区間に届いていない／傾きが負でない＝減衰していない）

    `detail=True` の情報:
        'slope_db_per_s' 傾き / 'intercept_db' 切片 / 'r2' 決定係数 /
        'xi' ISO 3382-2 の非線形性 1000(1-r²) / 'range_db' 使った区間 [dB] /
        'range_s' 使った区間 [s] / 'samples' 使った点の数 / 'fit' 読み取り方
    """
    if fit not in DECAY_FITS:
        raise ValueError(f"fit は {DECAY_FITS} のいずれかです: {fit!r}")

    start = _crossing_index(decay_db, db_start)
    stop = _crossing_index(decay_db, db_end)
    if start < 0 or stop < 0 or stop <= start:
        return (np.nan, None) if detail else np.nan

    if fit == DECAY_FIT_CROSSING:
        # **参照実装**（元コード 134 行）。区間の 2 点しか見ない。
        # ISO 3382 が求めているのは区間全体への回帰なので、規格準拠の数字には使わない
        value = (stop - start) * dt * 60.0 / (db_start - db_end)
        if not detail:
            return value
        return value, {"slope_db_per_s": (db_end - db_start) / ((stop - start) * dt),
                       "intercept_db": np.nan, "r2": np.nan, "xi": np.nan,
                       "range_db": (float(db_start), float(db_end)),
                       "range_s": (start * dt, stop * dt),
                       "samples": 2, "fit": fit}

    # ---- ISO 3382：評価区間の全サンプルへ直線を最小二乗で当てる ----
    time = np.arange(start, stop + 1) * dt
    level = np.asarray(decay_db[start:stop + 1], dtype=float)
    slope, intercept = np.polyfit(time, level, 1)
    if not np.isfinite(slope) or slope >= 0.0:
        # 減衰していない（適合不良）。黙って変な値を返さない
        return (np.nan, None) if detail else np.nan
    value = -60.0 / slope
    if not detail:
        return value

    # 決定係数。直線からのずれ＝減衰が曲がっている／途中で床に当たっている度合い。
    # ISO 3382-2 はこれを ξ = 1000(1-r²) という形で使う（小さいほど直線的）
    residual = level - (slope * time + intercept)
    variance = float(np.sum((level - level.mean()) ** 2))
    r2 = 1.0 - float(np.sum(residual ** 2)) / variance if variance > 0.0 else np.nan
    return value, {"slope_db_per_s": float(slope), "intercept_db": float(intercept),
                   "r2": r2, "xi": np.nan if np.isnan(r2) else 1000.0 * (1.0 - r2),
                   "range_db": (float(db_start), float(db_end)),
                   "range_s": (float(time[0]), float(time[-1])),
                   "samples": int(len(time)), "fit": fit}


def curvature_percent(decay_db, dt, fit=DEFAULT_DECAY_FIT):
    """**ISO 3382 の曲率 C** [%]。

        C = (T30 / T20 - 1) × 100

    減衰がまっすぐなら 0 に近い。反射回数不足で音が途中で切れている場合や、
    音場が拡散していない（小さい室・平行面・面ごとに吸音率が大きく違う）場合に大きくなる。

    ★**残響時間をどの区間で読んだかに依らない**。C は減衰曲線そのものの性質なので、
      `decay_curves` に db_max / db_min を渡していても T20・T30 の区間で測る。
      2026-09-05 より前は「評価区間を半分にした値との食い違い」（-5〜-20 dB）を
      曲率と呼んでいたが、これは ISO の C ではなく、実測で偽警報を出していた
      （`test.dxf` の 2 kHz で ISO の C が +0.4% なのに -11.8% と出た）。
    """
    (t20_start, t20_end), (t30_start, t30_end) = CURVATURE_RANGES
    t20 = _decay_time(decay_db, dt, t20_start, t20_end, fit=fit)
    t30 = _decay_time(decay_db, dt, t30_start, t30_end, fit=fit)
    if np.isnan(t20) or np.isnan(t30) or t20 == 0.0:
        return np.nan
    return (t30 / t20 - 1.0) * 100.0


def _decay_floor_db(decay_db):
    """減衰曲線がどこまで下がったか [dB]（見えている範囲）。届いていなければ大きい値。"""
    finite = np.asarray(decay_db, dtype=float)
    finite = finite[np.isfinite(finite)]
    return float(finite.min()) if len(finite) else np.nan


def decay_curves(time, ir, frequencies=None, db_max=DB_MAX, db_min=DB_MIN,
                 band_width=ab.BAND_WIDTH_OCTAVE,
                 method="butter", order=FILTER_ORDER, numtaps=FIR_NUMTAPS,
                 fit=DEFAULT_DECAY_FIT, verbose=True):
    """オクターブバンドごとの減衰曲線と残響時間を求める。元コード 99〜135 行。

    引数:
        time        (n,) 時間ベクトル [s]（等間隔であること）
        ir          (n,) インパルス応答
        frequencies (nf,) | None  中心周波数。None なら 63〜8k の 8 バンド
        db_max / db_min  評価区間 [dB]。既定 -5 / -35（＝T30）
                    -5 / -25 なら T20、0 / -10 なら EDT 相当
        fit         'least_squares'（ISO 3382。既定）| 'crossing'（2 点法・参照実装）

    戻り値: dict
        'frequencies' / 'time' / 'decay' (nf, n) [dB] /
        'reverberation_time' (nf,) [s] / 'curvature' (nf,) [%]（**ISO 3382 の C**）/
        'nonlinearity' (nf,)（**ISO 3382-2 の ξ = 1000(1-r²)**）/
        'r2' (nf,) / 'floor_db' (nf,)（減衰曲線が下がりきった深さ）/
        'fit' (nf,) 適合の情報 dict のリスト / 'fit_method'

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
    nonlinearity = np.empty(len(frequencies))
    r_squared = np.empty(len(frequencies))
    floor_db = np.empty(len(frequencies))
    fits = [None] * len(frequencies)

    for i, fc in enumerate(frequencies):
        band = octave_bandpass(ir, fc, sampling_frequency, method=method,
                               band_width=band_width,
                               order=order, numtaps=numtaps)
        dc = schroeder_integral(band)
        if dc[0] <= 0.0:
            decay[i] = -np.inf
            rt[i] = curvature[i] = np.nan
            nonlinearity[i] = r_squared[i] = floor_db[i] = np.nan
            continue

        with np.errstate(divide="ignore"):
            decay[i] = 10.0 * np.log10(np.maximum(dc, 0.0) / dc[0])

        # ★ISO 3382：評価区間の全サンプルへの最小二乗回帰（`fit` で 2 点法にも切替可）
        rt[i], info = _decay_time(decay[i], dt, db_max, db_min, fit=fit, detail=True)
        fits[i] = info
        nonlinearity[i] = np.nan if info is None else info["xi"]
        r_squared[i] = np.nan if info is None else info["r2"]

        # 曲率（ISO 3382 の C = (T30/T20 - 1)×100）。減衰がまっすぐなら 0 % に近い。
        # 反射回数不足で音が途中で切れている場合や、音場が拡散していない場合に大きくなる
        curvature[i] = curvature_percent(decay[i], dt, fit=fit)
        floor_db[i] = _decay_floor_db(decay[i])

    if verbose:
        label = f"T{int(abs(db_min - db_max)):d}"
        for fc, value, c, xi in zip(frequencies, rt, curvature, nonlinearity):
            if np.isnan(value):
                print(f"[reverberation] {fc:7.0f} Hz : {label} = 算出不可"
                      f"（減衰が {db_min:.0f} dB に届いていません）")
                continue
            marks = []
            if not np.isnan(c) and abs(c) > CURVATURE_LIMIT_PERCENT:
                marks.append(f"★曲率 {c:+.0f}%")
            elif not np.isnan(c):
                marks.append(f"曲率 {c:+.0f}%")
            if not np.isnan(xi) and xi > NONLINEARITY_LIMIT:
                marks.append(f"★ξ {xi:.0f}")
            elif not np.isnan(xi):
                marks.append(f"ξ {xi:.1f}")
            note = f"  ({' / '.join(marks)})" if marks else ""
            print(f"[reverberation] {fc:7.0f} Hz : {label} = {value:.3f} s{note}")
        _print_quality_note(frequencies, curvature, nonlinearity, floor_db, db_min)

    return {"frequencies": frequencies,
            "time": dt * np.arange(len(ir)),
            "decay": decay,
            "reverberation_time": rt,
            "curvature": curvature,
            "nonlinearity": nonlinearity,
            "r2": r_squared,
            "floor_db": floor_db,
            "fit": fits,
            "fit_method": fit}


def _print_quality_note(frequencies, curvature, nonlinearity, floor_db, db_min):
    """曲率・非線形性・減衰の余裕について、読むときの注意を出す。

    ★**ISO 3382 は評価区間の下端よりさらに 10 dB 下まで減衰が見えていることを求める**
      （T30 なら -45 dB）。実測では暗騒音の話だが、シミュレーションでは
      反射回数の打ち切りとインパルス応答長の打ち切りが同じ役目の「床」を作る。
    """
    bad_c = [f"{f:.0f}Hz" for f, c in zip(frequencies, curvature)
             if not np.isnan(c) and abs(c) > CURVATURE_LIMIT_PERCENT]
    bad_xi = [f"{f:.0f}Hz" for f, x in zip(frequencies, nonlinearity)
              if not np.isnan(x) and x > NONLINEARITY_LIMIT]
    if bad_c or bad_xi:
        if bad_c:
            print(f"[reverberation] ★ 曲率が {CURVATURE_LIMIT_PERCENT:.0f}% を超えたバンド"
                  f"（{' / '.join(bad_c)}）は減衰が直線でないのでそのまま信用しないこと。")
        if bad_xi:
            print(f"[reverberation] ★ 非線形性 ξ（ISO 3382-2）が "
                  f"{NONLINEARITY_LIMIT:.0f} を超えたバンド（{' / '.join(bad_xi)}）は"
                  f"回帰直線に乗っていません。")
        print("[reverberation]   よくある原因:")
        print("[reverberation]   ・最大反射回数 nref の不足で後部残響が切れている"
              "（procedure.py が別途エネルギーで判定して警告する）")
        print("[reverberation]   ・音場が拡散していない。小さい室・平行面・"
              "面ごとに吸音率が大きく違う場合は、**減衰が本当に 2 段階になる**ので"
              "曲率が出るのが正しい（EDT と T30 の差にも表れる）")

    # ISO 3382 が求める余裕（評価区間の下端 + 10 dB）まで減衰が見えているか
    needed = db_min - DECAY_MARGIN_DB
    short = [f"{f:.0f}Hz" for f, floor in zip(frequencies, floor_db)
             if not np.isnan(floor) and floor > needed]
    if short:
        print(f"[reverberation] ★ ISO 3382 は評価区間の下端よりさらに "
              f"{DECAY_MARGIN_DB:.0f} dB 下（{needed:.0f} dB）まで減衰が見えていることを"
              f"求めます。届いていないバンド: {' / '.join(short)}")
        print("[reverberation]   最大反射回数 nref を増やすか、"
              "インパルス応答の長さ（max_time）を延ばしてください。")


def decay_measures(time, ir, frequencies=None, measures=None, method="butter",
                   band_width=ab.BAND_WIDTH_OCTAVE,
                   order=FILTER_ORDER, numtaps=FIR_NUMTAPS,
                   fit=DEFAULT_DECAY_FIT, verbose=True):
    """**EDT / T20 / T30 をまとめて求める。**

    減衰曲線は 1 回だけ作り、評価区間を変えて読み取るので `decay_curves` を
    3 回呼ぶより速い。60 dB 減を厳密に見ない実務の使い方に合わせた入口。

    引数:
        measures : {名前: (開始dB, 終了dB)} | None
                   None なら EDT / T20 / T30（`DECAY_MEASURES`）
        fit      : 'least_squares'（ISO 3382。既定）| 'crossing'（2 点法・参照実装）

    戻り値: dict
        'frequencies' / 'time' / 'decay' (nf, n) [dB]
        'measures'      {名前: (nf,) [s]}
        'nonlinearity'  {名前: (nf,)}   ISO 3382-2 の ξ = 1000(1-r²)
        'fit'           {名前: [適合の情報 dict]}
        'curvature'     (nf,) [%]   **ISO 3382 の C = (T30/T20 - 1)×100**
        'floor_db'      (nf,)       減衰曲線が下がりきった深さ
        'fit_method'
    """
    if measures is None:
        measures = DECAY_MEASURES
    base = decay_curves(time, ir, frequencies=frequencies, method=method,
                        order=order, numtaps=numtaps, verbose=False,
                        band_width=band_width, fit=fit)
    dt = float(base["time"][1] - base["time"][0])

    values, nonlinearity, fits = {}, {}, {}
    for name, (db_start, db_end) in measures.items():
        pairs = [_decay_time(d, dt, db_start, db_end, fit=fit, detail=True)
                 for d in base["decay"]]
        values[name] = np.array([v for v, _ in pairs])
        nonlinearity[name] = np.array([np.nan if i is None else i["xi"]
                                       for _, i in pairs])
        fits[name] = [i for _, i in pairs]

    if verbose:
        # **周波数は横**（table.py の共通ルール）
        rows = list(values.items())
        rows += [(f"ξ({name})", nonlinearity[name]) for name in values]
        rows.append(("曲率%", base["curvature"]))
        print_frequency_table(base["frequencies"], rows)
        # 横向きの表にすると行末に ★ を付けられないので、
        # どのバンドが該当するかを注意書きの側に書く
        deepest = min((db_end for _, db_end in measures.values()), default=DB_MIN)
        _print_quality_note(base["frequencies"], base["curvature"],
                            base["nonlinearity"], base["floor_db"], deepest)

    return {"frequencies": base["frequencies"], "time": base["time"],
            "decay": base["decay"], "measures": values,
            "nonlinearity": nonlinearity, "fit": fits,
            "curvature": base["curvature"], "floor_db": base["floor_db"],
            "fit_method": fit}


# インパルス応答の立ち上がりとみなすレベル（ピークからの差 [dB]）。ISO 3382-1
ONSET_THRESHOLD_DB = -20.0


def _onset_index(energy, threshold_db=ONSET_THRESHOLD_DB):
    """**インパルス応答の立ち上がり**の添字（ISO 3382-1）。

    ピークより `threshold_db` 低いレベルを**最初に超える**点。
    探すのはピークまでの範囲だけ（後ろの残響に同じレベルの点があっても拾わない）。
    """
    energy = np.asarray(energy, dtype=float)
    peak = int(np.argmax(energy))
    if energy[peak] <= 0.0:
        return 0
    threshold = energy[peak] * 10.0 ** (threshold_db / 10.0)
    above = np.nonzero(energy[:peak + 1] >= threshold)[0]
    return int(above[0]) if len(above) else peak


def clarity_measures(time, ir, frequencies=None, method="butter",
                     band_width=ab.BAND_WIDTH_OCTAVE,
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

    ★**起点は「ピークより 20 dB 低いレベルを最初に超える点」**（ISO 3382-1。
      2026-09-05 に変更）。以前は各バンドで**エネルギーが最大になる時刻**（`argmax`）を
      直接音としていたが、これには 2 つの問題があった。

      1. 強い初期反射が直接音より大きいと、起点がその反射まで飛ぶ
      2. **帯域フィルタを通すと直接音のパルスが時間的に広がる**ので、
         低域では `argmax` が真の立ち上がりより**あと**に来る。63 Hz では
         フィルタの応答自体が数十 ms に及び、50 ms の窓が構造的にずれて
         **C50 / D50 が系統的に大きく出ていた**

      探すのはピークまでの範囲だけ（「最初に超える」の意味）。

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
        band = octave_bandpass(ir, fc, fs, method=method, order=order,
                               numtaps=numtaps, band_width=band_width)
        energy = band ** 2
        total = energy.sum()
        if total <= 0.0:
            continue

        # 直接音の到来時刻を起点にする（ISO 3382-1：ピークより 20 dB 低い点）
        start = _onset_index(energy)
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
        # **周波数は横**（table.py の共通ルール）
        print_frequency_table(frequencies, [
            ("C50[dB]", result["C50"]), ("C80[dB]", result["C80"]),
            ("D50", result["D50"]), ("Ts[ms]", result["Ts"] * 1000.0)])
    return result


def write_clarity_measures(filename, result):
    """明瞭度系の指標を CSV に保存する。**周波数は横**（table.py の共通ルール）。"""
    return tb.write_frequency_table(filename, result["frequencies"], {
        "C50_db": result["C50"], "C80_db": result["C80"],
        "D50": result["D50"], "Ts_s": result["Ts"],
    })


def write_decay_measures(filename, result):
    """EDT / T20 / T30 を CSV に保存する。**周波数は横**（table.py の共通ルール）。

    ★**適合の質も一緒に書く**（2026-09-05）。ISO 3382 の回帰にすると
    決定係数が手に入るので、`ξ = 1000(1-r²)`（ISO 3382-2 の非線形性）を
    指標ごとに並べ、曲率 C と「減衰が実際にどこまで下がったか」も添える。
    **どのバンドを読んでよいかを、警告文ではなく数値で残す**ため。

    区分付きの表（`table.write_sectioned_table`）にはしない。
    `summary.py` が行の名前で拾う形になっており、行を足すだけなら
    まとめ表の側に手を入れずに済むため（読み手にとっても 1 行 1 指標のまま）。
    """
    rows = {f"{name}_s": values for name, values in result["measures"].items()}
    for name, values in result.get("nonlinearity", {}).items():
        rows[f"{name}_xi"] = values
    rows["curvature_percent"] = result["curvature"]
    if result.get("floor_db") is not None:
        rows["decay_floor_db"] = result["floor_db"]
    return tb.write_frequency_table(filename, result["frequencies"], rows)


def write_reverberation_time(filename, result):
    """残響時間を CSV に保存する。元コード 141〜145 行。

    `decay_measures()` の結果（'measures' を持つ）でも
    `decay_curves()` の結果（'reverberation_time' を持つ）でも受け付ける。
    """
    if "measures" in result:
        return write_decay_measures(filename, result)
    return tb.write_frequency_table(filename, result["frequencies"], {
        "reverberation_time_s": result["reverberation_time"],
        "curvature_percent": result["curvature"],
    })


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
    # 1 行が時刻なので、周波数は**列**（table.py の共通ルール）
    header = ["time_s"] + [tb.band_column("decay", f, "db")
                           for f in result["frequencies"]]
    rows = np.column_stack([time, decay.T])
    np.savetxt(filename, rows, delimiter=",", header=",".join(header),
               comments="", fmt="%.12g")
    return filename


def reverberation_time(time, ir, rt_filename=None, decay_filename=None,
                       frequencies=None, measures=None, method="butter",
                       band_width=ab.BAND_WIDTH_OCTAVE,
                       fit=DEFAULT_DECAY_FIT, verbose=True):
    """残響指標（EDT / T20 / T30）の算出と保存をまとめて行う。

    measures を渡せば評価区間を変えられる（既定は `DECAY_MEASURES`）。
    fit を渡せば読み取り方を変えられる（既定は ISO 3382 の最小二乗回帰）。
    """
    result = decay_measures(time, ir, frequencies=frequencies, measures=measures,
                            method=method, verbose=verbose, band_width=band_width,
                            fit=fit)
    if rt_filename is not None:
        write_reverberation_time(rt_filename, result)
        if verbose:
            print(f"[reverberation] 残響時間を書き出しました: {rt_filename}")
    if decay_filename is not None:
        write_decay_curve(decay_filename, result)
        if verbose:
            print(f"[reverberation] 減衰曲線を書き出しました: {decay_filename}")
    return result
