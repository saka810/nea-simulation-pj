"""計算結果を図にして PNG で保存する（TODO G-3 / G-4）。

プロジェクトフォルダの `図/` に書き出す。**インパルス応答は正規化して保存する**
（絶対振幅の校正は未了なので、そのまま出しても値に意味が無い。TODO E-11）。

    impulse_response.png   インパルス応答（正規化）。広帯域とバンド別
    decay.png              Schroeder 逆積分の減衰曲線（バンド別）
    reverberation.png      EDT / T20 / T30 と統計残響式の比較
    clarity.png            C50 / C80 / D50 / Ts（明瞭度系）
    absorption.png         レイヤ別の吸音率と面積
    pulses.png             パルス列（到来時刻とエネルギー）
    direction.png          伝搬方向（真上から見た人と、音の来る向き）
    modes.png              受音点のスペクトルと直方体の固有周波数
    mode_buildup.png       経路差から見たモードの積み上げと、吸音の効果

画面には出さず**ファイルに書くだけ**なので Agg バックエンドを使う
（GUI のイベントループと取り合わないようにするため）。
"""

import os

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

import view_model_gui as vg

BACKGROUND = "#12151c"
PANEL = "#1a1f29"
TEXT = "#d6dae2"
GRID = "#2a3140"
ACCENT = "#4cc9f0"

# バンドごとの色（低域 → 高域で寒色 → 暖色）
BAND_COLORS = ["#4361ee", "#4cc9f0", "#4cc38a", "#f7b801", "#f18701", "#e5484d",
               "#b5179e", "#7209b7"]

_FONT_READY = False


def use_japanese_font():
    """matplotlib に日本語フォントを登録する（無ければ黙って諦める）。"""
    global _FONT_READY
    if _FONT_READY:
        return
    _FONT_READY = True
    path = vg.japanese_font()
    if not path or not os.path.exists(path):
        return
    try:
        from matplotlib import font_manager
        font_manager.fontManager.addfont(path)
        name = font_manager.FontProperties(fname=path).get_name()
        plt.rcParams["font.family"] = name
    except Exception as e:      # フォントが読めなくても図は出したい
        print(f"[plots] 日本語フォントを登録できませんでした（英字で描きます）: {e}")


def _style(ax, title=None, xlabel=None, ylabel=None):
    ax.set_facecolor(PANEL)
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color(GRID)
    ax.tick_params(colors=TEXT, labelsize=9)
    if title:
        ax.set_title(title, color=TEXT, fontsize=11, pad=10)
    if xlabel:
        ax.set_xlabel(xlabel, color=TEXT, fontsize=9)
    if ylabel:
        ax.set_ylabel(ylabel, color=TEXT, fontsize=9)
    return ax


def _figure(nrows=1, ncols=1, figsize=(10, 6)):
    use_japanese_font()
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, facecolor=BACKGROUND)
    return fig, axes


def _save(fig, path, legend_axes=()):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    for ax in np.atleast_1d(legend_axes).ravel():
        legend = ax.get_legend()
        if legend is not None:
            legend.get_frame().set_facecolor(PANEL)
            legend.get_frame().set_edgecolor(GRID)
            for text in legend.get_texts():
                text.set_color(TEXT)
    fig.tight_layout()
    fig.savefig(path, dpi=140, facecolor=BACKGROUND)
    plt.close(fig)
    return path


def _band_color(i):
    return BAND_COLORS[i % len(BAND_COLORS)]


# ------------------------------------------------------------------------------
# ③ インパルス応答
# ------------------------------------------------------------------------------

def impulse_response(path, time, ir, frequencies=None, bands=None, max_time=None):
    """インパルス応答を**正規化して**描く。

    上段が波形（最大値を 1 に正規化）、下段がエネルギー時間曲線（dB, 最大値を 0 dB に）。
    正規化するのは、絶対振幅の校正が未了で値そのものに意味が無いため（TODO E-11）。
    どちらも**形と減衰の様子を読むための図**として使う。
    """
    time = np.asarray(time, dtype=float)
    ir = np.asarray(ir, dtype=float)
    peak = float(np.max(np.abs(ir))) or 1.0
    normalised = ir / peak

    if max_time is None:
        # 音が -70 dB まで落ちたところで切る（後ろの無音を延々と描かない）
        energy = normalised ** 2
        loud = np.nonzero(energy > 10 ** (-7.0))[0]
        max_time = float(time[loud[-1]]) * 1.05 if len(loud) else float(time[-1])
    keep = time <= max_time

    fig, axes = _figure(2, 1, figsize=(11, 7))
    _style(axes[0], "インパルス応答（最大値で正規化）", None, "振幅 [-]")
    axes[0].plot(time[keep], normalised[keep], color=ACCENT, linewidth=0.6)
    axes[0].set_xlim(0, max_time)
    axes[0].set_ylim(-1.05, 1.05)

    _style(axes[1], "エネルギー時間曲線（最大値を 0 dB）", "時間 [s]", "レベル [dB]")
    with np.errstate(divide="ignore"):
        db = 10.0 * np.log10(np.maximum(normalised ** 2, 1e-12))
    # バンド別を重ねるので、広帯域の生波形は背景として薄く敷く
    axes[1].plot(time[keep], db[keep], color=ACCENT, linewidth=0.4, alpha=0.35)
    axes[1].set_xlim(0, max_time)
    axes[1].set_ylim(-80, 5)

    if bands is None and frequencies is not None:
        # バンド別のエネルギー時間曲線を重ねる。低域ほど長く尾を引くのが見える
        import reverberation as rv
        fs = 1.0 / float(time[1] - time[0])
        bands = [rv.octave_bandpass(ir, fc, fs) for fc in frequencies]

    if bands is not None and frequencies is not None:
        # 生のままだと反射音の干渉で激しく暴れて重ねられないので、
        # **5 ms の移動平均で包絡線にしてから**描く（減衰の傾きを読むための図なので
        # 細かい山谷は落として構わない）
        window = max(1, int(round(0.005 / float(time[1] - time[0]))))
        kernel = np.ones(window) / window
        for i, (fc, band) in enumerate(zip(frequencies, bands)):
            envelope = np.convolve((band / peak) ** 2, kernel, mode="same")
            with np.errstate(divide="ignore"):
                band_db = 10.0 * np.log10(np.maximum(envelope, 1e-12))
            axes[1].plot(time[keep], band_db[keep], color=_band_color(i),
                         linewidth=1.1, alpha=0.9, label=f"{fc:.0f} Hz")
        axes[1].plot([], [], color=ACCENT, linewidth=0.8, alpha=0.35, label="広帯域")
        axes[1].legend(loc="upper right", fontsize=8, ncol=4)
        axes[1].set_title("エネルギー時間曲線（最大値を 0 dB / バンド別は 5 ms 移動平均）",
                          color=TEXT, fontsize=11, pad=10)
    return _save(fig, path, axes[1])


def pulses(path, time, energy, frequencies, distance=None):
    """バックトレースが出したパルス列を、到来時刻とエネルギーで描く。

    インパルス応答を作る前の**素の反射音の並び**。
    エコーの粗密や、後期に経路を取りこぼしていないかがここで分かる。
    """
    time = np.asarray(time, dtype=float)
    energy = np.atleast_2d(np.asarray(energy, dtype=float))
    if distance is not None:
        # 受音点に届くエネルギー（距離減衰を入れる）
        energy = energy / np.asarray(distance, dtype=float)[:, None] ** 2

    fig, axes = _figure(2, 1, figsize=(11, 7))
    _style(axes[0], "パルス列（受音点に届くエネルギー、最大を 0 dB）", None, "レベル [dB]")
    reference = float(energy.max()) or 1.0
    for i, fc in enumerate(frequencies):
        with np.errstate(divide="ignore"):
            db = 10.0 * np.log10(np.maximum(energy[:, i] / reference, 1e-12))
        axes[0].scatter(time, db, s=6, color=_band_color(i), alpha=0.6,
                        label=f"{fc:.0f} Hz", edgecolors="none")
    axes[0].set_ylim(-100, 5)
    axes[0].legend(loc="upper right", fontsize=8, ncol=3)

    # 経路数の時間分布。**拡散音場なら時刻の 2 乗で増える**はずなので、
    # 減っていく形なら後期の経路を取りこぼしている（受音球が小さすぎる）
    _style(axes[1], "経路数の時間分布（拡散音場なら時刻とともに増える）",
           "到来時刻 [s]", "経路数 / 50 ms")
    if len(time):
        edges = np.arange(0.0, float(time.max()) + 0.05, 0.05)
        axes[1].hist(time, bins=edges, color=ACCENT, alpha=0.85)
    return _save(fig, path, axes[0])


# ------------------------------------------------------------------------------
# ④ 音響指標
# ------------------------------------------------------------------------------

def decay_curves(path, time, decay, frequencies):
    """Schroeder 逆積分の減衰曲線をバンド別に描く。"""
    fig, ax = _figure(figsize=(10, 6))
    _style(ax, "減衰曲線（Schroeder 逆積分）", "時間 [s]", "残留エネルギー [dB]")
    for i, fc in enumerate(frequencies):
        ax.plot(time, decay[i], color=_band_color(i), linewidth=1.2,
                label=f"{fc:.0f} Hz")
    for level, style in ((-5, ":"), (-25, ":"), (-35, ":")):
        ax.axhline(level, color=GRID, linestyle=style, linewidth=1.0)
    ax.set_ylim(-70, 2)
    ax.legend(loc="upper right", fontsize=9, ncol=2)
    return _save(fig, path, ax)


def reverberation(path, frequencies, measures, statistical=None, curvature=None):
    """EDT / T20 / T30 と、統計残響式の値を並べて描く。

    **どちらが正しいという話ではない**。統計式は拡散音場を前提にした平均像、
    シミュレーションは特定の受音点での実際の減衰。
    大きく食い違うときは音場が拡散していないか設定に問題があるかの手がかりになる。
    """
    fig, axes = _figure(1, 2 if curvature is not None else 1,
                        figsize=(13, 6) if curvature is not None else (10, 6))
    axes = np.atleast_1d(axes)
    ax = axes[0]
    _style(ax, "残響時間", "周波数 [Hz]", "残響時間 [s]")

    x = np.arange(len(frequencies))
    styles = {"EDT": ("o-", "#4cc9f0"), "T20": ("s-", "#4cc38a"), "T30": ("^-", "#f7b801")}
    for name, values in measures.items():
        marker, colour = styles.get(name, ("d-", ACCENT))
        ax.plot(x, values, marker, color=colour, linewidth=1.8, markersize=6, label=name)

    if statistical:
        import reverberation as rv
        dashes = {"sabine": (2, 2), "eyring": (1, 1), "eyring_knudsen": (5, 2)}
        # Eyring と Eyring-Knudsen の差がそのまま空気吸収の効き
        labels = rv.STATISTICAL_LABELS
        for key, dash in dashes.items():
            if key in statistical:
                ax.plot(x, statistical[key], color="#8b929e", linewidth=1.2,
                        dashes=dash, label=labels[key])
    ax.set_xticks(x)
    ax.set_xticklabels([f"{f:.0f}" for f in frequencies])
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", fontsize=9, ncol=2)

    if curvature is not None:
        _style(axes[1], "曲率 C（T30 と T20 の食い違い）", "周波数 [Hz]", "曲率 [%]")
        colours = ["#e5484d" if abs(c) > 10 else "#4cc38a" for c in curvature]
        axes[1].bar(x, curvature, color=colours, alpha=0.9)
        axes[1].axhline(0, color=GRID, linewidth=1.0)
        for limit in (-10, 10):
            axes[1].axhline(limit, color="#e5484d", linestyle=":", linewidth=1.0)
        axes[1].set_xticks(x)
        axes[1].set_xticklabels([f"{f:.0f}" for f in frequencies])
    return _save(fig, path, axes)


def clarity(path, result):
    """C50 / C80 / D50 / Ts を描く（`reverberation.clarity_measures()` の結果）。"""
    frequencies = result["frequencies"]
    x = np.arange(len(frequencies))
    fig, axes = _figure(1, 3, figsize=(14, 5))

    _style(axes[0], "明瞭度 C（大きいほど明瞭）", "周波数 [Hz]", "C [dB]")
    axes[0].bar(x - 0.2, result["C50"], width=0.4, color="#4cc9f0", label="C50（音声）")
    axes[0].bar(x + 0.2, result["C80"], width=0.4, color="#f7b801", label="C80（音楽）")
    axes[0].axhline(0, color=GRID, linewidth=1.0)
    axes[0].legend(loc="upper left", fontsize=9)

    _style(axes[1], "D50（明瞭度・0〜1）", "周波数 [Hz]", "D50 [-]")
    axes[1].bar(x, result["D50"], color="#4cc38a")
    axes[1].axhline(0.5, color="#e5484d", linestyle=":", linewidth=1.0)
    axes[1].set_ylim(0, 1)

    _style(axes[2], "重心時刻 Ts（小さいほど明瞭）", "周波数 [Hz]", "Ts [ms]")
    axes[2].bar(x, np.asarray(result["Ts"]) * 1000.0, color="#b5179e")

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels([f"{f:.0f}" for f in frequencies])
    return _save(fig, path, axes[0])


def absorption(path, surface, frequencies):
    """レイヤ別の吸音率と面積（`reverberation.surface_summary()` の結果）。"""
    names = surface["names"]
    x = np.arange(len(frequencies))
    fig, axes = _figure(1, 2, figsize=(14, 6))

    _style(axes[0], "レイヤ別の吸音率（乱入射）", "周波数 [Hz]", "吸音率 α [-]")
    for i, (name, alpha) in enumerate(zip(names, surface["absorption"])):
        axes[0].plot(x, alpha, "o-", color=_band_color(i), linewidth=1.5,
                     markersize=5, label=name)
    axes[0].set_ylim(0, 1.05)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([f"{f:.0f}" for f in frequencies])
    axes[0].legend(loc="upper left", fontsize=8)

    _style(axes[1], "材料別の面積", "面積 [m²]", None)
    y = np.arange(len(names))
    axes[1].barh(y, surface["areas"], color=ACCENT, alpha=0.9)
    axes[1].set_yticks(y)
    axes[1].set_yticklabels(names, fontsize=8)
    axes[1].invert_yaxis()
    for j, area in enumerate(surface["areas"]):
        axes[1].text(area, j, f" {area:.1f}", va="center", color=TEXT, fontsize=8)
    return _save(fig, path, axes[0])


# ------------------------------------------------------------------------------
# ⑤ 伝搬方向（G-5）
# ------------------------------------------------------------------------------

# 方向を何分割して集計するか（10° 刻み）
DIRECTION_SECTORS = 36

# 初期と後期の境目 [s]。C50 と同じ 50 ms にしてある
EARLY_LIMIT = 0.050


def _head_patch(ax, radius, colour="#d6dae2"):
    """真上から見た人の絵を原点に描く（正面は +X 方向）。

    円（頭）＋鼻の三角＋両耳。**どちらが前か**が一目で分かればよいので、
    細かい形は追わない。図の主役は方向分布のほうなので控えめな色にしてある。
    """
    from matplotlib.patches import Circle, Ellipse, Polygon

    ax.add_patch(Circle((0, 0), radius, facecolor="#2a3140",
                        edgecolor=colour, linewidth=1.2, zorder=5))
    # 鼻（正面 = +X）
    nose = radius * 0.42
    ax.add_patch(Polygon([[radius * 0.92, -nose * 0.5],
                          [radius * 1.38, 0.0],
                          [radius * 0.92, nose * 0.5]],
                         closed=True, facecolor=colour, edgecolor=colour, zorder=6))
    # 両耳（±Y）
    for sign in (-1.0, 1.0):
        ax.add_patch(Ellipse((0.0, sign * radius * 1.02),
                             width=radius * 0.30, height=radius * 0.52,
                             facecolor="#2a3140", edgecolor=colour,
                             linewidth=1.0, zorder=6))


def _head_azimuth(project, results=None):
    """図に使う顔の向き [度]。**受音点ごとの値が結果に入っていればそれを使う。**

    `project.head_azimuth` は受音点ごとのリストになりうる（2026-08-20）。
    そのまま渡すと図が壊れるので、ここで 1 つの数値に落とす。
    """
    if isinstance(results, dict) and results.get("head_azimuth") is not None:
        return float(results["head_azimuth"])
    value = getattr(project, "head_azimuth", 0.0)
    if value is None:
        return 0.0
    return float(value) if np.isscalar(value) else float(value[0] if len(value) else 0.0)


def direction_histogram(direction, energy, head_azimuth=0.0,
                        sectors=DIRECTION_SECTORS):
    """到来方向を水平面で集計する。

    引数:
        direction (N,3) 到来方向の単位ベクトル（受音点から音が来る向き）
        energy    (N,)  その経路のエネルギー
        head_azimuth    人の正面方向 [度]（真上から見て +X から反時計回り）

    戻り値 (角度の境目 [rad], 区間ごとのエネルギー)。
    角度は**頭の正面を 0 とした相対方位**。0=正面 / 90=左 / 180=後ろ / 270=右。

    ★上下は見ない（水平面へ投影する）。実務では水平面で足りるため（ユーザー判断）。
    """
    direction = np.atleast_2d(np.asarray(direction, dtype=float))
    energy = np.asarray(energy, dtype=float).ravel()

    azimuth = np.arctan2(direction[:, 1], direction[:, 0])
    azimuth = azimuth - np.deg2rad(head_azimuth)        # 頭の正面を 0 にする
    azimuth = np.mod(azimuth, 2.0 * np.pi)

    edges = np.linspace(0.0, 2.0 * np.pi, sectors + 1)
    totals, _ = np.histogram(azimuth, bins=edges, weights=energy)
    return edges, totals


def propagation_direction(path, direction, energy, distance=None, time=None,
                          head_azimuth=0.0, frequencies=None, band=None,
                          dynamic_range=30.0, reflection_count=None):
    """受音点にどの方向から音が来ているかを、真上から見た図にする（G-5）。

    左が**全体**（初期と後期を重ねる）、右が**バンド別**。
    真ん中に人の絵を置き、正面を上に向けて描く。

    半径は**最大を 0 dB としたときのレベル**で、外側ほど強い。
    `dynamic_range` dB より弱い方向は中心に潰れる（見たいのは強い方向の偏りなので）。

    `reflection_count` を渡すと、**直接音（反射 0 回）の方向を矢印で示す**。
    どちらから直接届いているかは真っ先に知りたい情報なので、
    分布の線に埋もれないよう別に描く（遮蔽されていて直接音が無ければ描かない）。
    """
    direction = np.atleast_2d(np.asarray(direction, dtype=float))
    energy = np.atleast_2d(np.asarray(energy, dtype=float))
    if distance is not None:
        # 受音点に届くエネルギー（距離減衰を入れる）
        energy = energy / np.asarray(distance, dtype=float)[:, None] ** 2

    # 直接音（反射 0 回）の方位。頭の正面を 0 とした相対方位で持つ
    direct_azimuth = None
    if reflection_count is not None:
        direct = np.nonzero(np.asarray(reflection_count) == 0)[0]
        if len(direct):
            # 複数あることは無いはずだが、あれば最も強いものを採る
            pick = direct[int(np.argmax(energy[direct].sum(axis=1)))]
            angle = np.arctan2(direction[pick, 1], direction[pick, 0])
            direct_azimuth = float(np.mod(angle - np.deg2rad(head_azimuth),
                                          2.0 * np.pi))

    fig, axes = _figure(1, 2, figsize=(13, 6.5))

    def draw(ax, series, title):
        """series = [(ラベル, 色, エネルギー (N,))]"""
        ax.set_facecolor(PANEL)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(title, color=TEXT, fontsize=11, pad=12)

        # 目盛りの円（dB）と方位の線
        for level in np.arange(0.25, 1.01, 0.25):
            ax.add_patch(plt.Circle((0, 0), level, fill=False, color=GRID,
                                    linewidth=0.7, zorder=1))
            decibel = -dynamic_range * (1.0 - level)
            ax.text(0.0, level, f"{decibel + 0.0:.0f} dB".replace("-0 dB", "0 dB"),
                    color="#7f8794", fontsize=7, ha="center", va="bottom", zorder=2)
        for angle, name in ((0, "正面"), (90, "左"), (180, "後ろ"), (270, "右")):
            a = np.deg2rad(angle)
            # 図では正面を上に描く（真上から見た人の絵に合わせる）
            x, y = -np.sin(a), np.cos(a)
            ax.plot([0, x * 1.12], [0, y * 1.12], color=GRID, linewidth=0.7, zorder=1)
            ax.text(x * 1.22, y * 1.22, name, color=TEXT, fontsize=9,
                    ha="center", va="center")

        for label, colour, values in series:
            edges, totals = direction_histogram(direction, values,
                                                head_azimuth=head_azimuth)
            peak = totals.max()
            if peak <= 0.0:
                continue
            with np.errstate(divide="ignore"):
                level = 10.0 * np.log10(np.maximum(totals / peak, 1e-12))
            radius = np.clip(1.0 + level / dynamic_range, 0.0, 1.0)
            centres = 0.5 * (edges[:-1] + edges[1:])
            # 閉じた折れ線にする
            centres = np.append(centres, centres[0])
            radius = np.append(radius, radius[0])
            x = -radius * np.sin(centres)
            y = radius * np.cos(centres)
            ax.plot(x, y, color=colour, linewidth=1.8, label=label, zorder=4)
            ax.fill(x, y, color=colour, alpha=0.18, zorder=3)

        # 直接音の方向。分布の線に埋もれないよう、外から中心へ向かう矢印で示す
        if direct_azimuth is not None:
            x, y = -np.sin(direct_azimuth), np.cos(direct_azimuth)
            ax.annotate("", xy=(x * 0.20, y * 0.20), xytext=(x * 1.30, y * 1.30),
                        arrowprops=dict(arrowstyle="-|>", color="#ff5f5f",
                                        linewidth=2.4, shrinkA=0, shrinkB=0),
                        zorder=7)
            ax.text(x * 1.36, y * 1.36, "直接音", color="#ff5f5f", fontsize=9,
                    ha="center", va="center", zorder=7)

        _head_patch(ax, 0.15)
        ax.set_xlim(-1.55, 1.55)
        ax.set_ylim(-1.45, 1.45)
        legend = ax.legend(loc="lower right", fontsize=8)
        if legend is not None:
            legend.get_frame().set_facecolor(PANEL)
            legend.get_frame().set_edgecolor(GRID)
            for text in legend.get_texts():
                text.set_color(TEXT)

    broadband = energy.sum(axis=1)
    series = [("すべて", ACCENT, broadband)]
    if time is not None:
        t = np.asarray(time, dtype=float)
        early = np.where(t < t.min() + EARLY_LIMIT, broadband, 0.0)
        late = np.where(t < t.min() + EARLY_LIMIT, 0.0, broadband)
        series = [("すべて", "#8b929e", broadband),
                  (f"初期（〜{EARLY_LIMIT * 1000:.0f} ms）", "#4cc9f0", early),
                  ("後期", "#f7b801", late)]
    draw(axes[0], series, "伝搬方向（真上から見た図）")

    if frequencies is not None and energy.shape[1] == len(frequencies):
        bands = [(f"{f:.0f} Hz", _band_color(i), energy[:, i])
                 for i, f in enumerate(frequencies)]
        draw(axes[1], bands, "バンド別")
    else:
        axes[1].axis("off")

    fig.text(0.5, 0.02, f"人の正面 = 方位 {head_azimuth:.0f}°"
                        f"（真上から見て +X から反時計回り）",
             color="#7f8794", fontsize=8, ha="center")
    return _save(fig, path)


# ------------------------------------------------------------------------------
# ⑥ モード分布（G-6）
# ------------------------------------------------------------------------------

def room_modes(lengths, limit, sound_velocity=343.0):
    """直方体の固有周波数を `limit` [Hz] まで列挙する。

        f = (c/2) √((nx/Lx)² + (ny/Ly)² + (nz/Lz)²)

    書籍『建築音響物理学』2.1 節（波動音響理論）の直方体音場。
    実際の室は直方体でないことが多いが、**低域でどのあたりに固有周波数が
    並ぶか**の目安としては外形寸法から出したもので十分役に立つ。

    戻り値 (周波数 (M,), 次数 (M,3))。周波数の昇順。
    """
    lengths = np.asarray(lengths, dtype=float)
    if np.any(lengths <= 0.0):
        return np.array([]), np.zeros((0, 3), dtype=int)
    # 各軸で必要な最大次数（他の軸が 0 のときに limit に届く数）
    top = np.maximum(1, np.ceil(2.0 * limit * lengths / sound_velocity).astype(int))
    grid = np.stack(np.meshgrid(*[np.arange(n + 1) for n in top], indexing="ij"),
                    axis=-1).reshape(-1, 3)
    grid = grid[np.any(grid > 0, axis=1)]           # (0,0,0) は音場ではない
    frequency = 0.5 * sound_velocity * np.sqrt(
        np.sum((grid / lengths[None, :]) ** 2, axis=1))
    keep = frequency <= limit
    frequency, grid = frequency[keep], grid[keep]
    order = np.argsort(frequency)
    return frequency[order], grid[order]


def schroeder_frequency(reverberation_time, volume):
    """シュレーダー周波数 `f = 2000 √(T/V)` [Hz]。

    **これより低い帯域は個々の固有振動が分離して見える**（モードの世界）、
    高い帯域はモードが重なり合って統計的に扱える（幾何音響の世界）、という境目。
    幾何音響シミュレーションの結果を低域でどこまで信用してよいかの目安になる。
    """
    if not volume or reverberation_time is None or not np.isfinite(reverberation_time):
        return None
    return 2000.0 * np.sqrt(float(reverberation_time) / float(volume))


def spectrum_peaks(frequency, level, prominence=3.0, distance_hz=2.0):
    """スペクトルの山を拾う。

    戻り値は山の添字。`prominence` [dB] より高く盛り上がっていて、
    互いに `distance_hz` 以上離れているものだけを残す
    （隣り合う細かい揺れを全部拾うと印だらけになって逆に読めない）。
    """
    from scipy.signal import find_peaks

    step = float(frequency[1] - frequency[0])
    indices, _ = find_peaks(level, prominence=prominence,
                            distance=max(1, int(round(distance_hz / step))))
    return indices


def mode_distribution(path, time, ir, lengths=None, volume=None,
                      reverberation_time=None, sound_velocity=343.0,
                      max_frequency=200.0, label_peaks=8):
    """受音点のモード分布（G-6）。

    上段：インパルス応答の**スペクトル**（低域）。山が立っているところが
          その受音点で強く出ている固有振動。**山には印を付ける**
          （減衰が大きい室ではモードの幅が広く、目で山を拾いにくいため）。
    下段：外形寸法から求めた**直方体の固有周波数**と、その累積数。

    ★寸法は外形（バウンディングボックス）なので、室が直方体でなければ
      **目安**にしかならない。それでも「低域にどれくらい隙間があるか」は分かる。

    ★**山が鋭くないのは正しい**。モードの半値幅はおよそ `2.2 / T` [Hz] で、
      T = 0.5 s なら 4.4 Hz にもなる。残響が短い室ほど山はなだらかになる。
      図にもこの値を書き添えてある。
    """
    time = np.asarray(time, dtype=float)
    ir = np.asarray(ir, dtype=float)
    fs = 1.0 / float(time[1] - time[0])

    spectrum = np.abs(np.fft.rfft(ir))
    frequency = np.fft.rfftfreq(len(ir), d=1.0 / fs)
    keep = (frequency > 0) & (frequency <= max_frequency)
    frequency, spectrum = frequency[keep], spectrum[keep]
    peak = spectrum.max() or 1.0
    with np.errstate(divide="ignore"):
        level = 20.0 * np.log10(np.maximum(spectrum / peak, 1e-6))

    modes, orders = (room_modes(lengths, max_frequency, sound_velocity)
                     if lengths is not None else (np.array([]), None))
    f_schroeder = schroeder_frequency(reverberation_time, volume)

    # 個々の固有周波数を線で描くのは、**分離して見える帯域だけ**にする。
    # 全部引くと灰色のベタ塗りになって情報が消える（実際にそうなった）。
    # 境目はシュレーダー周波数そのもの（これより上はモードが重なり合う、という定義）
    separable = min(max_frequency,
                    f_schroeder if f_schroeder else max_frequency * 0.25)

    fig, axes = _figure(2, 1, figsize=(12, 8))
    _style(axes[0], "受音点のスペクトル（最大を 0 dB）", None, "レベル [dB]")
    axes[0].axvspan(separable, max_frequency, color="#2a3140", alpha=0.35, zorder=0)
    for f in modes[modes <= separable]:
        axes[0].axvline(f, color="#8b929e", linewidth=0.7, alpha=0.65, zorder=1)
    axes[0].plot(frequency, level, color=ACCENT, linewidth=1.2, zorder=3)

    # 山に印を付ける。減衰が大きい室ではモードの幅が広く、
    # 線を眺めただけでは「どこが山か」を拾いにくい（ユーザー指摘）
    found = spectrum_peaks(frequency, level)
    if len(found):
        axes[0].plot(frequency[found], level[found], "v", color="#ffd166",
                     markersize=6, linestyle="none", zorder=5,
                     label=f"スペクトルの山（{len(found)} 個）")
        # 強い順にいくつかだけ周波数を書く（全部書くと読めない）
        strongest = found[np.argsort(level[found])[::-1][:label_peaks]]
        for i in strongest:
            axes[0].annotate(f"{frequency[i]:.0f}", (frequency[i], level[i]),
                             textcoords="offset points", xytext=(0, 9),
                             ha="center", color="#ffd166", fontsize=8, zorder=6)

    axes[0].set_xlim(0, max_frequency)
    axes[0].set_ylim(-60, 8)
    if len(modes):
        axes[0].plot([], [], color="#8b929e", linewidth=0.7,
                     label="直方体の固有周波数（外形寸法から）")
    axes[0].axvspan(np.nan, np.nan, color="#2a3140", alpha=0.35,
                    label="モードが重なり合う帯域")

    if f_schroeder is not None:
        for ax in axes:
            ax.axvline(f_schroeder, color="#e5484d", linestyle="--", linewidth=1.2,
                       zorder=2)
        axes[0].plot([], [], color="#e5484d", linestyle="--", linewidth=1.2,
                     label=f"シュレーダー周波数 {f_schroeder:.0f} Hz")
    axes[0].legend(loc="lower right", fontsize=8)

    # **山がなだらかなのは正しい**ことを図の中で断っておく。
    # モードの半値幅は 2.2/T で決まり、残響が短い室ほど広がる
    if reverberation_time:
        axes[0].set_title(
            f"受音点のスペクトル（最大を 0 dB / モードの半値幅は "
            f"約 2.2÷T = {2.2 / reverberation_time:.1f} Hz）",
            color=TEXT, fontsize=11, pad=10)

    _style(axes[1], "固有周波数の累積数（外形寸法から）", "周波数 [Hz]", "累積モード数")
    if len(modes):
        axes[1].step(modes, np.arange(1, len(modes) + 1), where="post",
                     color="#4cc38a", linewidth=1.5, label="数え上げ")
        if lengths is not None:
            # 累積モード数の近似式（体積・面積・辺長の 3 項）。
            # 体積項だけだと数え上げより下に出て比較にならない
            lx, ly, lz = np.asarray(lengths, dtype=float)
            box_volume = lx * ly * lz
            surface = 2.0 * (lx * ly + ly * lz + lz * lx)
            edges = 4.0 * (lx + ly + lz)
            f = np.linspace(1.0, max_frequency, 400)
            approx = ((4.0 * np.pi / 3.0) * box_volume * (f / sound_velocity) ** 3
                      + (np.pi / 4.0) * surface * (f / sound_velocity) ** 2
                      + (edges / 8.0) * (f / sound_velocity))
            axes[1].plot(f, approx, color="#8b929e", linestyle=":", linewidth=1.4,
                         label="近似 4πVf³/3c³ + πSf²/4c² + Lf/8c")
        axes[1].legend(loc="upper left", fontsize=8)
        axes[1].set_xlim(0, max_frequency)
    else:
        axes[1].text(0.5, 0.5, "寸法が分からないので固有周波数を出せません",
                     color="#7f8794", ha="center", transform=axes[1].transAxes)
    return _save(fig, path, axes)


# ------------------------------------------------------------------------------
# ⑥-2 モードの積み上げ（経路差からの重なり）
# ------------------------------------------------------------------------------

def pulse_spectrum(time, amplitude, max_time=None, bin_rate=8000.0):
    """パルス列を**位相込みで足し合わせて**周波数特性にする。

        H(f) = Σ_n a_n exp(-i 2π f t_n)

    到来時刻の差がそのまま経路差なので、`exp` の中身は経路差ぶんの位相。
    **同じ位相で重なった経路だけが足し算になる**（ずれていれば打ち消し合う）。
    だから `a_n = 1`（減衰を考えない）で計算すると、`|H(f)|` は
    「その周波数で同位相に重なった経路の本数」そのものになる。
    1 本なら 1、2 本重なれば 2。3 本がばらばらの向きなら 1 前後にしかならない。

    **経路差から固有周波数が出る理由**：往復の経路差 Δ で戻ってきた音は、
    Δ が波長の整数倍のときに元の音と強め合う。つまり `f = m·c/Δ`。
    直方体の x 方向の往復は Δ = 2Lx なので `f = m·c/(2Lx)` となり、
    これは軸モード `f = (c/2)·(m/Lx)` に一致する。`room_modes` の式そのもの。

    **速さ**：素直に書くと「パルス本数 × 周波数点数」の掛け算になって重いが、
    到来時刻を細かい格子に**投げ込んで（ヒストグラムにして）から FFT** すれば
    O(n log n) で済む。要はインパルス応答を作るのと同じ手順なので、
    パルスが 100 万本あっても一瞬で終わる（重み違いで 3 回やっても同じ）。

    引数:
        time      (n,)  到来時刻 [s]
        amplitude (n,) | (n,k)  重み（**振幅**。エネルギーなら sqrt を取ってから渡す）
        bin_rate  時刻を丸める格子の細かさ [Hz]。8 kHz なら 0.125 ms = 4.3 cm 刻みで、
                  200 Hz での位相誤差は最大 4.5°。低域の図には十分

    戻り値: (周波数 (m,), H (m,) または (m,k) complex)
    """
    time = np.asarray(time, dtype=float)
    amplitude = np.asarray(amplitude, dtype=float)
    single = amplitude.ndim == 1
    amplitude = amplitude.reshape(len(time), -1)

    if max_time is None:
        max_time = float(time.max()) * 1.05 + 1e-3
    count = max(2, int(round(max_time * bin_rate)))

    index = np.rint(time * bin_rate).astype(int)
    keep = (index >= 0) & (index < count)

    binned = np.empty((count, amplitude.shape[1]))
    for k in range(amplitude.shape[1]):
        binned[:, k] = np.bincount(index[keep], weights=amplitude[keep, k],
                                   minlength=count)

    spectrum = np.fft.rfft(binned, axis=0)
    frequency = np.fft.rfftfreq(count, d=1.0 / bin_rate)
    return frequency, (spectrum[:, 0] if single else spectrum)


def rebin_spectrum(frequency, spectrum, step):
    """スペクトルを `step` [Hz] 幅にまとめ直す（＝分析帯域幅 `step` の狭帯域分析）。

    FFT の刻みは応答長で決まってしまう（3 秒なら 0.333 Hz）。細かすぎると
    線がぎざぎざして読みにくいので、実務で見たい幅（1 Hz 程度）にまとめる。

    **点を間引くのではなく二乗平均する**。`pulse_spectrum` の①（減衰なし）は
    減衰が無いぶん櫛の歯が鋭く、**幅は応答長でしか決まらない**（3 秒なら 0.333 Hz）。
    間引くと山を跨いで見落とす。二乗平均なら「その 1 Hz の中にどれだけ入っているか」
    になるので取りこぼさない。

    位相がばらばらのときの目安 √N も**そのまま保たれる**（帯域幅に依らない）ので、
    図に引いた基準線を引き直さずに済む。

    引数の `spectrum` は複素でも振幅でもよい。戻り値は (帯域の中心周波数, 振幅)。
    """
    frequency = np.asarray(frequency, dtype=float)
    magnitude = np.abs(np.asarray(spectrum))
    single = magnitude.ndim == 1
    magnitude = magnitude.reshape(len(frequency), -1)

    native = float(frequency[1] - frequency[0])
    if not step or step <= native:
        return frequency, (magnitude[:, 0] if single else magnitude)

    index = np.floor((frequency - frequency[0]) / step + 1e-9).astype(int)
    count = np.bincount(index)
    centre = np.bincount(index, weights=frequency) / np.maximum(count, 1)
    power = np.column_stack([np.bincount(index, weights=magnitude[:, k] ** 2)
                             for k in range(magnitude.shape[1])])
    keep = count > 0
    result = np.sqrt(power[keep] / count[keep, None])
    return centre[keep], (result[:, 0] if single else result)


def _band_of(frequency, band_frequencies):
    """各周波数点に、いちばん近いオクターブバンドの添字を割り当てる（対数距離）。"""
    band_frequencies = np.asarray(band_frequencies, dtype=float)
    ratio = np.log2(np.maximum(frequency, 1e-6)[:, None] / band_frequencies[None, :])
    return np.argmin(np.abs(ratio), axis=1)


def _spread_labels(index, frequency, strength, count, min_gap):
    """強い順に拾いつつ、**すでに拾ったものと近すぎるものは飛ばす**。

    山の印は全部付けてよいが、文字まで全部書くと重なって読めなくなる
    （実際に「125 Hz」と「128 Hz」が重なった）。
    """
    picked = []
    for i in index[np.argsort(strength[index])[::-1]]:
        if len(picked) >= count:
            break
        if all(abs(frequency[i] - frequency[j]) >= min_gap for j in picked):
            picked.append(i)
    return picked


def mode_buildup(path, time, energy, distance=None, frequencies=None,
                 lengths=None, volume=None, reverberation_time=None,
                 sound_velocity=343.0, max_frequency=200.0,
                 min_frequency=None, frequency_step=1.0, max_time=None,
                 bin_rate=8000.0, label_peaks=6):
    """経路差から見たモードの積み上げと、**吸音の効き方**（G-6b）。

    上段：**減衰をいっさい考えずに**、経路差から決まる位相だけで経路を積み上げる。
          その周波数で同位相に重なった経路の本数がそのまま縦軸になる。
          室形状だけで決まる図で、吸音材を貼っても変わらない。
          ばらばらに重なったときの目安 √N を引いてあり、**これを超えている山が
          「本当に強め合っている周波数」**。固有周波数の線とだいたい合う。

    下段：同じ足し算に減衰を入れる。2 本の線の差が**吸音の効果**。
          ・完全反射：距離減衰 1/d だけ。壁は 100% 跳ね返る
          ・設計の吸音：距離減衰 + 各面の吸音率（斜入射のエネルギー反射率）
          どちらも直接音（1/d₀）を 0 dB にしてあるので、**0 dB より上に出た分が
          室の反響ぶん**。塗りつぶした帯が、その周波数で吸音が削った量。

    ★空気吸収は入れていない。この帯域では効かないため
      （125 Hz で 100 m 進んで 0.05 dB。ISO 9613-1、20℃ 湿度 40%）。

    ★吸音率はオクターブバンドでしか持っていないので、各周波数点には**いちばん
      近いバンドの値**を当てる。6 バンド（125 Hz 始まり）だと 88 Hz より下は
      すべて 125 Hz の吸音率になる。低域を細かく見たいときは 8 バンドにすること。

    ★パルス列は受音球が拾えた経路の**標本**なので、本数の絶対値は音線数と
      受音球の半径で変わる。**山と谷の位置**と**2 本の線の差**を読むための図。

    引数:
        time     (n,)      到来時刻 [s]
        energy   (n,b)     バンド別エネルギー（**面の吸音だけ**が入っている。
                           距離減衰と空気吸収は impulse.py 側で掛ける約束）
        distance (n,)      経路長 [m]。省略すると time × 音速
        frequency_step     図の周波数刻み [Hz]。既定 1 Hz。
                           FFT の刻み（応答長の逆数。3 秒なら 0.333 Hz）は細かすぎて
                           線がぎざぎざするので、**間引かず二乗平均**でまとめ直す
                           （`rebin_spectrum`）。None なら FFT の刻みのまま
    """
    time = np.asarray(time, dtype=float)
    energy = np.atleast_2d(np.asarray(energy, dtype=float))
    if distance is None:
        distance = time * sound_velocity
    distance = np.maximum(np.asarray(distance, dtype=float), 1e-9)
    if frequencies is None:
        from absorption import octave_bands
        frequencies = octave_bands(energy.shape[1])
    frequencies = np.asarray(frequencies, dtype=float)

    direct = float(distance.min())          # 直接音（いちばん短い経路）を基準にする

    # 重み 3 種。**振幅**なのでエネルギーは sqrt を取る
    weights = np.column_stack([
        np.ones_like(distance),             # ①減衰なし＝本数の積み上げ
        1.0 / distance,                     # ②完全反射（距離減衰のみ）
        np.sqrt(energy) / distance[:, None],  # ③設計の吸音（バンドごと）
    ])
    f, H = pulse_spectrum(time, weights, max_time=max_time, bin_rate=bin_rate)

    modes, _ = (room_modes(lengths, max_frequency, sound_velocity)
                if lengths is not None else (np.array([]), None))

    # **0 Hz では全部の経路が同位相になる**（位相差が消えるので当たり前）。
    # そのぶんが図の左端に山として立ってしまい、肝心の帯域が潰れる。
    # いちばん低い固有周波数より下には見るものが無いので、そこから描く
    if min_frequency is None:
        min_frequency = float(modes[0]) if len(modes) else max_frequency / 40.0

    keep = (f >= min_frequency) & (f <= max_frequency)
    f, H = f[keep], H[keep]

    # FFT の刻みは応答長で決まる（3 秒なら 0.333 Hz）。細かすぎて線がぎざぎざするので
    # 実務で見る幅にまとめる。**間引かず二乗平均**（理由は rebin_spectrum の説明）
    width = max(frequency_step or 0.0, float(f[1] - f[0]))
    f, magnitude = rebin_spectrum(f, H, frequency_step)

    counted = magnitude[:, 0]
    free = magnitude[:, 1]
    absorbed = magnitude[:, 2:][np.arange(len(f)), _band_of(f, frequencies)]

    with np.errstate(divide="ignore"):
        free_db = 20.0 * np.log10(np.maximum(free * direct, 1e-6))
        absorbed_db = 20.0 * np.log10(np.maximum(absorbed * direct, 1e-6))

    f_schroeder = schroeder_frequency(reverberation_time, volume)
    separable = min(max_frequency,
                    f_schroeder if f_schroeder else max_frequency * 0.25)

    fig, axes = _figure(2, 1, figsize=(12, 8.5))

    # ---- 上段：本数の積み上げ ----
    _style(axes[0], f"経路の重なり（減衰を考えない＝室形状だけで決まる "
                    f"/ 分析帯域幅 {width:.2g} Hz）",
           None, "同位相で重なった経路の本数")
    axes[0].axvspan(separable, max_frequency, color="#2a3140", alpha=0.35, zorder=0)
    for fm in modes[modes <= separable]:
        axes[0].axvline(fm, color="#8b929e", linewidth=0.7, alpha=0.65, zorder=1)
    axes[0].plot(f, counted, color=ACCENT, linewidth=1.2, zorder=3)

    # 位相がばらばらなら和はランダムウォークになって √N 程度にしかならない。
    # **この線を超えている山だけが、本当に強め合っている周波数**
    random_walk = np.sqrt(len(time))
    axes[0].axhline(random_walk, color="#7f8794", linestyle="--", linewidth=1.0,
                    zorder=2, label=f"ばらばらに重なった場合の目安 √N = {random_walk:.0f}")
    found = spectrum_peaks(f, 20.0 * np.log10(np.maximum(counted, 1e-9)))
    found = found[counted[found] > random_walk]
    if len(found):
        axes[0].plot(f[found], counted[found], "v", color="#ffd166", markersize=6,
                     linestyle="none", zorder=5,
                     label=f"強め合っている周波数（{len(found)} 個）")
        for i in _spread_labels(found, f, counted, label_peaks,
                                max_frequency / 12.0):
            axes[0].annotate(f"{f[i]:.0f} Hz / {counted[i]:.0f} 本",
                             (f[i], counted[i]), textcoords="offset points",
                             xytext=(0, 8), ha="center", color="#ffd166",
                             fontsize=8, zorder=6)
    if len(modes):
        axes[0].plot([], [], color="#8b929e", linewidth=0.7,
                     label="直方体の固有周波数（外形寸法から）")
    axes[0].set_yscale("log")
    axes[0].set_xlim(min_frequency, max_frequency)
    # 深い谷まで入れると肝心の帯域が潰れるので下は √N の 1/20 で切る。
    # 上は凡例と山の文字がぶつからないよう余白を取る
    axes[0].set_ylim(max(counted.min() * 0.7, random_walk / 20.0),
                     counted.max() * 5.0)
    axes[0].legend(loc="upper right", fontsize=8)

    # ---- 下段：吸音の効果 ----
    _style(axes[1], "吸音の効果（直接音を 0 dB とした重なり）",
           "周波数 [Hz]", "レベル [dB]")
    axes[1].axvspan(separable, max_frequency, color="#2a3140", alpha=0.35, zorder=0)
    axes[1].fill_between(f, absorbed_db, free_db, color="#e5484d", alpha=0.22,
                         zorder=2, label="吸音が削っている分")
    axes[1].plot(f, free_db, color="#f7b801", linewidth=1.2, zorder=3,
                 label="完全反射（距離減衰のみ）")
    axes[1].plot(f, absorbed_db, color=ACCENT, linewidth=1.2, zorder=4,
                 label="設計の吸音あり")
    axes[1].axhline(0.0, color="#7f8794", linestyle=":", linewidth=1.0, zorder=1)
    axes[1].set_xlim(min_frequency, max_frequency)

    gap = float(np.mean(free_db - absorbed_db))
    axes[1].set_title(f"吸音の効果（直接音を 0 dB とした重なり / "
                      f"平均 {gap:.1f} dB 下がっている）", color=TEXT,
                      fontsize=11, pad=10)

    if f_schroeder is not None:
        for ax in axes:
            ax.axvline(f_schroeder, color="#e5484d", linestyle="--",
                       linewidth=1.2, zorder=2)
        axes[1].plot([], [], color="#e5484d", linestyle="--", linewidth=1.2,
                     label=f"シュレーダー周波数 {f_schroeder:.0f} Hz")
    axes[1].legend(loc="upper right", fontsize=8)
    return _save(fig, path, axes)


# ------------------------------------------------------------------------------
# まとめて書き出す
# ------------------------------------------------------------------------------

def save_all(project, results, verbose=True):
    """計算結果 dict から図を一式書き出す。書けたファイルのリストを返す。

    `results` は `procedure.process()` の戻り値と同じ形。
    **無い項目は飛ばす**（インパルス応答を作っていない場合など）。
    """
    written = []

    def emit(name, func, *args, **kwargs):
        path = project.figure_path(name)
        try:
            func(path, *args, **kwargs)
        except Exception as e:
            print(f"[plots] {name} を描けませんでした: {type(e).__name__}: {e}")
            return
        written.append(path)
        if verbose:
            print(f"[plots] {path}")

    impulse = results.get("impulse")
    if impulse is not None:
        time, ir = impulse[0], impulse[1]
        bands = impulse[2] if len(impulse) > 2 else None
        emit("impulse_response.png", impulse_response, time, ir,
             frequencies=results.get("frequencies"), bands=bands)

    pulse_list = results.get("pulses")
    if pulse_list is not None and len(pulse_list):
        emit("pulses.png", pulses, pulse_list.time, pulse_list.energy,
             results.get("frequencies"), distance=pulse_list.distance)

    rt = results.get("reverberation")
    if rt is not None:
        stat = results.get("statistical")
        emit("reverberation.png", reverberation, rt["frequencies"], rt["measures"],
             statistical=stat, curvature=rt.get("curvature"))
        emit("decay.png", decay_curves, rt["time"], rt["decay"], rt["frequencies"])

    clarity_result = results.get("clarity")
    if clarity_result is not None:
        emit("clarity.png", clarity, clarity_result)

    stat = results.get("statistical")
    if stat is not None and "surface" in stat:
        emit("absorption.png", absorption, stat["surface"], stat["frequencies"])

    # ⑤ 伝搬方向。人の正面方向はプロジェクトが持つ（GUI で決める）
    if pulse_list is not None and len(pulse_list):
        emit("direction.png", propagation_direction,
             pulse_list.direction, pulse_list.energy,
             distance=pulse_list.distance, time=pulse_list.time,
             head_azimuth=_head_azimuth(project, results),
             frequencies=results.get("frequencies"),
             reflection_count=pulse_list.reflection_count)

    # ⑥ モード分布。寸法は外形（バウンディングボックス）から
    model = results.get("model")
    if model is not None and model.extents is not None:
        lengths = np.asarray(model.extents[1]) - np.asarray(model.extents[0])
        rt = results.get("reverberation")
        # 代表の残響時間として中域（500 Hz に最も近いバンド）の T30 を使う
        t_mid = None
        if rt is not None:
            frequencies = np.asarray(rt["frequencies"], dtype=float)
            t_mid = float(rt["measures"]["T30"][
                int(np.argmin(np.abs(frequencies - 500.0)))])
        volume = (project.volume if project.volume
                  else (abs(model.volume) if model.volume else None))
        atmosphere = results.get("atmosphere")
        velocity = atmosphere.sound_velocity if atmosphere else 343.0
        if impulse is not None:
            emit("modes.png", mode_distribution, impulse[0], impulse[1],
                 lengths=lengths, volume=volume, reverberation_time=t_mid,
                 sound_velocity=velocity)

        # ⑥-2 経路差からの積み上げと吸音の効果。パルス列だけで描ける
        if pulse_list is not None and len(pulse_list):
            emit("mode_buildup.png", mode_buildup, pulse_list.time,
                 pulse_list.energy, distance=pulse_list.distance,
                 frequencies=results.get("frequencies"), lengths=lengths,
                 volume=volume, reverberation_time=t_mid,
                 sound_velocity=velocity)

    return written
