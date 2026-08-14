"""計算結果を図にして PNG で保存する（TODO G-3 / G-4）。

プロジェクトフォルダの `図/` に書き出す。**インパルス応答は正規化して保存する**
（絶対振幅の校正は未了なので、そのまま出しても値に意味が無い。TODO E-11）。

    impulse_response.png   インパルス応答（正規化）。広帯域とバンド別
    decay.png              Schroeder 逆積分の減衰曲線（バンド別）
    reverberation.png      EDT / T20 / T30 と統計残響式の比較
    clarity.png            C50 / C80 / D50 / Ts（明瞭度系）
    absorption.png         レイヤ別の吸音率と面積
    pulses.png             パルス列（到来時刻とエネルギー）

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
        dashes = {"sabine": (2, 2), "eyring": (5, 2), "millington": (1, 1)}
        labels = {"sabine": "Sabine", "eyring": "Eyring", "millington": "Millington"}
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

    _style(axes[1], "レイヤ別の面積", "面積 [m²]", None)
    y = np.arange(len(names))
    axes[1].barh(y, surface["areas"], color=ACCENT, alpha=0.9)
    axes[1].set_yticks(y)
    axes[1].set_yticklabels(names, fontsize=8)
    axes[1].invert_yaxis()
    for j, area in enumerate(surface["areas"]):
        axes[1].text(area, j, f" {area:.1f}", va="center", color=TEXT, fontsize=8)
    return _save(fig, path, axes[0])


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

    return written
