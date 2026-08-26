# -*- coding: utf-8 -*-
"""**逆二乗則からのずれ**を測線ごとに出す（別ルートの出力）。

★2026-08-24 ユーザー要望

> 逆二乗用にテンプレートを常備しておく必要はないですが、
> これは別途別ルート？として逆二乗だしてほしいです。
> 受音点群も 3 つの方向別にレイヤー分けしています。

半無響室・無響室の**適合性の確認**（ISO 3745 / ISO 26101 の考え方）に使うもの。
音源から測線に沿って離れながら音圧レベルを測り、
**理想の 1/r²（距離 2 倍で −6 dB）からどれだけずれるか**を見る。

決めごと:
  ・**ふだんの出力には入れない**（2026-08-21 のユーザー判断で、音圧レベルの表に
    自由音場の値と差は入れないことにした）。**必要なときにこれを回す**
  ・測線は**受音点のレイヤ**で分ける（`rec1` `rec2` `rec3`）。
    `結果/<室>_測定点.csv` の「レイヤ」の列を見る
  ・基準の高さは**最小二乗**で決める（傾きは −20log10 r に固定）。
    どれか 1 点を基準にすると、その点の誤差が全部に乗るため。ISO 3745 と同じ考え方
  ・許容値（参考）は ISO 3745 の値を添える。**計算値なので合否判定ではない**

読むのは計算済みの結果だけ（`まとめ_音圧レベル.csv` と `測定点.csv`）。
**計算はやり直さない。**
"""
import csv
import io
import math
import os

import numpy as np

import project as pj
import table as tb

FILE_NAME = "逆二乗.csv"

# ★**許容偏差は JIS Z 8732**（2026-08-26 ユーザー指定。1/3 オクターブバンド）
#
#     1/3 オクターブバンド中心周波数 [Hz]   許容偏差 [dB]
#     ≦ 630                                ±2.5
#     800 〜 5000                           ±2.0
#     ≧ 6300                                ±3.0
#
# 以前は ISO 3745 の値（±1.5 / ±1.0）を参考として添えていた。
# **報告書に載せるのは JIS Z 8732 の値**なのでこちらに合わせる
TOLERANCE_STANDARD = "JIS Z 8732"
TOLERANCE = ((630.0, 2.5), (5000.0, 2.0), (float("inf"), 3.0))

# 受音点のレイヤ名 → 測線の呼び名（★2026-08-26 ユーザー指定）。
# 報告書とグラフの見出しに使う
TRACE_LABELS = {"rec1": "真上方向", "rec2": "短辺稜線方向", "rec3": "長辺稜線方向"}


def trace_label(layer):
    """レイヤ名を測線の呼び名にする（決めていないものはそのまま）。"""
    return TRACE_LABELS.get(str(layer), str(layer))

# バンドの色（図。8 バンドまではこの順で回す）
BAND_COLORS = ("#4c78a8", "#f58518", "#54a24b", "#e45756", "#72b7b2",
               "#b279a2", "#9d755d", "#eeca3b")
TRACE_COLORS = BAND_COLORS      # 昔の名前（参照しているところがあれば拾えるように）


def tolerance_of(frequency):
    """その帯域の許容偏差 [dB]（JIS Z 8732）。"""
    for limit, value in TOLERANCE:
        if float(frequency) <= limit:
            return value
    return TOLERANCE[-1][1]


# ---- 入力（計算済みの結果を読む）--------------------------------------------

def read_points(project):
    """`測定点.csv` から {名前: {"距離": …, "レイヤ": …, "座標": …}}。"""
    path = project.existing_result_path("points")
    if path is None or not os.path.exists(path):
        return {}
    points = {}
    with io.open(path, encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        return {}
    header = rows[0]
    distance_column = next((k for k, name in enumerate(header)
                            if name.endswith("からの距離_m")), None)
    layer_column = next((k for k, name in enumerate(header)
                         if name == "レイヤ"), None)
    for row in rows[1:]:
        if not row or row[0] != "受音点":
            continue
        name = row[1]
        try:
            distance = float(row[distance_column]) if distance_column else float("nan")
        except (TypeError, ValueError):
            distance = float("nan")
        points[name] = {
            "距離": distance,
            "レイヤ": (row[layer_column] if layer_column is not None
                       and layer_column < len(row) else ""),
            "座標": [float(v) for v in row[2:5]] if len(row) >= 5 else None,
        }
    return points


def read_levels(project):
    """`まとめ_音圧レベル.csv` から (周波数, {受音点名: Lp[バンド]})。"""
    import summary as sm

    path = os.path.join(project.folder, pj.RESULT_DIR,
                        project.prefixed(sm.LEVEL_FILE))
    if not os.path.exists(path):
        return None, {}
    table = tb.read_sectioned_table(path)
    if table is None:
        return None, {}
    levels = {}
    # まとめ表は **区分 = 受音点の名前（rec1…）／項目 = Lp_dB** の形
    for section, items in table["sections"].items():
        if section == "平均" or "Lp_dB" not in items:
            continue
        levels[section] = np.asarray(items["Lp_dB"], dtype=float)
    return table["frequencies"], levels


# ---- 評価 ------------------------------------------------------------------

def fit_reference(distances, levels):
    """傾きを −20log10 r に固定した最小二乗の高さ A [dB]。

    L(r) = A − 20 log10 r とみなしたときの A。**平均のずれが 0 になる**ので、
    どれか 1 点を基準にするより素直（ISO 3745 と同じ考え方）。
    """
    distances = np.asarray(distances, dtype=float)
    levels = np.asarray(levels, dtype=float)
    ok = np.isfinite(distances) & np.isfinite(levels) & (distances > 0)
    if not np.any(ok):
        return float("nan")
    return float(np.mean(levels[ok] + 20.0 * np.log10(distances[ok])))


def deviations(distances, levels):
    """理想の逆二乗からのずれ ΔL [dB]（と、当てはめた高さ）。"""
    reference = fit_reference(distances, levels)
    distances = np.asarray(distances, dtype=float)
    ideal = reference - 20.0 * np.log10(np.where(distances > 0, distances, np.nan))
    return np.asarray(levels, dtype=float) - ideal, reference


def free_field_levels(distances, source_power_db=None, atmosphere=None,
                      hemisphere=True):
    """**半自由音場の理論上の距離減衰** Lp(r) [dB]。

        Lp = Lw + 10log10(Q / (4π r²)) + 10log10(ρc/400)

    半空間（床の上に置いた音源）は指向係数 Q = 2。
    `source_power_db` が無ければ Lw = 0 dB として**相対値**を返す
    （計算側も PWL 未入力なら相対値なので、そのまま重ねられる）。
    """
    import atmosphere as at

    distances = np.asarray(distances, dtype=float)
    if atmosphere is None:
        atmosphere = at.Atmosphere()
    power = 0.0 if source_power_db is None else float(np.mean(
        np.atleast_1d(np.asarray(source_power_db, dtype=float))))
    directivity = 2.0 if hemisphere else 1.0
    impedance = atmosphere.density * atmosphere.sound_velocity
    return (power + 10.0 * np.log10(directivity
                                    / (4.0 * np.pi * distances ** 2))
            + 10.0 * np.log10(impedance / 400.0))


# バンドの中で何点の周波数を計算して平均するか（参考案件の「周波数平均」に相当）
COHERENT_LINES = 129


def coherent_band_levels(pulses, frequencies, atmosphere, band_width="1/3",
                         lines=None, source_power_db=None):
    """パルス列を**位相ごと**重ねたバンド別の音圧レベル [dB]（B の中身）。

    ★中身は `sound_level.coherent_band_levels()`（正式な置き場はあちら。
    音圧レベルの出力でも同じものを使うので、2026-08-26 に一本化した）。
    """
    import sound_level as sl

    return sl.coherent_band_levels(
        pulses.time, pulses.energy, pulses.distance, atmosphere, frequencies,
        band_width=band_width,
        lines=lines or sl.COHERENT_LINES, source_power_db=source_power_db)


def read_coherent_levels(project, verbose=True):
    """受音点ごとのパルス列から B の音圧レベルを出す。→ {受音点名: Lp[バンド]}

    パルス列（`結果/recN/…_pulses.csv`）を読むだけで、**計算はやり直さない**。
    """
    import absorption as ab
    import atmosphere as at
    import loop_noredundancy as ln

    frequencies = ab.frequency_bands(project.band_number,
                                     getattr(project, "band_width", "1/1"),
                                     getattr(project, "band_start", None))
    air = at.Atmosphere(temperature=project.temperature,
                        humidity=project.humidity, pressure=project.pressure)
    levels, missing = {}, 0
    for index in range(1, 1000):
        sub = pj.Project(project.folder,
                         **{k: getattr(project, k) for k in pj.DEFAULTS})
        sub.receiver_index = index
        path = sub.existing_result_path("pulses")
        if path is None or not os.path.exists(path):
            missing += 1
            if missing > 3 and levels:
                break
            continue
        missing = 0
        try:
            pulses = ln.PulseList.from_csv(path)
        except Exception as error:
            if verbose:
                print(f"[逆二乗] rec{index} のパルス列が読めません: "
                      f"{type(error).__name__}: {error}")
            continue
        levels[f"rec{index}"] = coherent_band_levels(
            pulses, frequencies, air,
            band_width=getattr(project, "band_width", "1/1"),
            source_power_db=project.source_power_db)
    if verbose and levels:
        print(f"[逆二乗] 複素和（参考）を {len(levels)} 点ぶん計算しました"
              f"（反射の位相ずれは 0 と仮定）")
    return levels


def evaluate(project, verbose=True, coherent=True):
    """測線ごとに逆二乗からのずれを出す。→ 結果の辞書（無ければ None）"""
    frequencies, levels = read_levels(project)
    if frequencies is None or not levels:
        if verbose:
            print("[逆二乗] 音圧レベルのまとめが見つかりません（先に計算してください）")
        return None
    points = read_points(project)
    if not points:
        if verbose:
            print("[逆二乗] 測定点の一覧が見つかりません")
        return None

    # ★B（位相を含む複素和）。**参考値**として一緒に持つ（2026-08-26）
    coherent_levels = read_coherent_levels(project, verbose=verbose) \
        if coherent else {}

    traces = {}
    for name, level in levels.items():
        info = points.get(name)
        if info is None or not np.isfinite(info["距離"]):
            continue
        traces.setdefault(trace_label(info["レイヤ"] or "測線"), []).append(
            (info["距離"], name, np.asarray(level, dtype=float),
             coherent_levels.get(name)))

    result = {"frequencies": np.asarray(frequencies, dtype=float), "traces": {}}
    for trace, entries in traces.items():
        entries.sort(key=lambda item: item[0])
        distances = np.array([d for d, _, _, _ in entries], dtype=float)
        names = [n for _, n, _, _ in entries]
        block = np.array([values for _, _, values, _ in entries], dtype=float)
        delta = np.empty_like(block)
        reference = np.empty(block.shape[1], dtype=float)
        for band in range(block.shape[1]):
            delta[:, band], reference[band] = deviations(distances, block[:, band])
        entry = {
            "names": names, "distances": distances, "levels": block,
            "deviation": delta, "reference": reference,
            "worst": np.nanmax(np.abs(delta), axis=0),
        }
        # B（複素和）も同じ形で持つ。全点そろっているときだけ
        pack = [values for _, _, _, values in entries]
        if all(values is not None for values in pack) and pack:
            block_c = np.array(pack, dtype=float)
            delta_c = np.empty_like(block_c)
            reference_c = np.empty(block_c.shape[1], dtype=float)
            for band in range(block_c.shape[1]):
                delta_c[:, band], reference_c[band] = deviations(
                    distances, block_c[:, band])
            entry["coherent"] = {
                "levels": block_c, "deviation": delta_c,
                "reference": reference_c,
                "worst": np.nanmax(np.abs(delta_c), axis=0)}
        result["traces"][trace] = entry
    if verbose:
        _report(result)
    return result


def _report(result):
    frequencies = result["frequencies"]
    print("[逆二乗] 測線ごとの最大のずれ [dB]（理想の 1/r² との差）")
    print("[逆二乗]   " + "測線".ljust(10)
          + "".join(f"{v:>9.0f}Hz" for v in frequencies))
    for trace, data in result["traces"].items():
        print("[逆二乗]   " + str(trace).ljust(10)
              + "".join(f"{v:>11.2f}" for v in data["worst"]))
        if "coherent" in data:
            print("[逆二乗]   " + (str(trace) + "(参考:複素和)").ljust(10)
                  + "".join(f"{v:>11.2f}"
                            for v in data["coherent"]["worst"]))
    print("[逆二乗]   " + "許容(参考)".ljust(10)
          + "".join(f"{tolerance_of(v):>11.1f}" for v in frequencies))


# ---- 出力 ------------------------------------------------------------------

def write_csv(project, result, verbose=True):
    """`結果/<室>_<条件>_逆二乗.csv`（区分付きの表。周波数は横）。"""
    frequencies = result["frequencies"]
    rows = []
    for trace, data in result["traces"].items():
        for name, distance, level, delta in zip(data["names"], data["distances"],
                                                data["levels"], data["deviation"]):
            rows.append((f"{trace} 音圧レベル", name, distance, level))
        for name, distance, delta in zip(data["names"], data["distances"],
                                         data["deviation"]):
            rows.append((f"{trace} 逆二乗からのずれ", name, distance, delta))
        rows.append((f"{trace} まとめ", "最大のずれ_dB", None, data["worst"]))
        rows.append((f"{trace} まとめ", "当てはめた高さ_dB", None, data["reference"]))
        coherent = data.get("coherent")
        if coherent is not None:
            # ★B（参考）。**位相ずれ 0 と仮定した複素和**
            for name, distance, level in zip(data["names"], data["distances"],
                                             coherent["levels"]):
                rows.append((f"{trace} 参考:複素和 音圧レベル", name, distance,
                             level))
            for name, distance, delta in zip(data["names"], data["distances"],
                                             coherent["deviation"]):
                rows.append((f"{trace} 参考:複素和 ずれ", name, distance, delta))
            rows.append((f"{trace} 参考:複素和 まとめ", "最大のずれ_dB", None,
                         coherent["worst"]))
    rows.append(("参考", "許容偏差_dB(ISO 3745)", None,
                 [tolerance_of(v) for v in frequencies]))

    path = os.path.join(project.folder, pj.RESULT_DIR,
                        project.prefixed(FILE_NAME))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tb.write_sectioned_table(path, frequencies, rows, value_label="音源距離_m")
    if verbose:
        print(f"[逆二乗] 表を書き出しました: {path}")
    return path


def bands_in_range(frequencies, low=None, high=None):
    """使うバンドの番号（`low` 〜 `high` [Hz]）。指定が無ければ全部。"""
    numbers = [k for k, v in enumerate(frequencies)
               if (low is None or v >= low - 1e-6)
               and (high is None or v <= high + 1e-6)]
    return numbers or list(range(len(frequencies)))


def write_figures(project, result, bands=None, low=None, high=500.0,
                  show_coherent=True, verbose=True):
    """★**測線ごとに 1 枚**の図（2026-08-26 ユーザー要望）。

    1 枚に 2 段:
      上段 … 音圧レベルの偏移。**半自由音場の理論上の距離減衰**も重ねる
      下段 … 逆二乗からのずれ。**JIS Z 8732 の許容偏差**（上限・下限）を引く

    実線が **A（エネルギー和）**、細い破線が **B（位相を含む複素和・参考）**。

    → 書き出したファイルの一覧
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    import plots as pl

    pl.use_japanese_font()
    frequencies = result["frequencies"]
    if bands is None:
        bands = bands_in_range(frequencies, low, high)

    written = []
    for trace, data in result["traces"].items():
        figure, (top, bottom) = plt.subplots(2, 1, figsize=(7.2, 8.4),
                                             sharex=True)
        distances = data["distances"]
        coherent = data.get("coherent") if show_coherent else None
        for order, band in enumerate(bands):
            colour = BAND_COLORS[order % len(BAND_COLORS)]
            label = f"{frequencies[band]:.0f} Hz"
            top.plot(distances, data["levels"][:, band], "o-", color=colour,
                     markersize=3.5, linewidth=1.2, label=label)
            bottom.plot(distances, data["deviation"][:, band], "o-",
                        color=colour, markersize=3.5, linewidth=1.2, label=label)
            if coherent is not None:
                # ★B（参考）は**細い破線**。同じ色で「同じ帯域の別の足し方」と分かる
                top.plot(distances, coherent["levels"][:, band], "--",
                         color=colour, linewidth=0.9, alpha=0.75)
                bottom.plot(distances, coherent["deviation"][:, band], "--",
                            color=colour, linewidth=0.9, alpha=0.75)

        # ★半自由音場（Q = 2）の理論線。計算側と同じ「Lw 未入力なら相対値」
        theory = free_field_levels(distances,
                                   source_power_db=getattr(project,
                                                           "source_power_db", None))
        top.plot(distances, theory, "--", color="#333333", linewidth=1.4,
                 label="半自由音場の理論（Q=2）")

        # ★許容偏差（JIS Z 8732）。バンドで値が変わるので、使ったバンドの分だけ引く
        allowed = sorted({tolerance_of(frequencies[b]) for b in bands})
        for value in allowed:
            bottom.axhline(+value, color="#c0392b", linewidth=1.2,
                           linestyle="--")
            bottom.axhline(-value, color="#c0392b", linewidth=1.2,
                           linestyle="--")
        # 線の少し上に置く（線に重ねると読めない）
        bottom.text(distances[0], allowed[-1] + 0.12,
                    f"許容偏差の上限値 ±{allowed[-1]:.1f} dB"
                    f"（{TOLERANCE_STANDARD}）", color="#c0392b", fontsize=8,
                    va="bottom")
        bottom.text(distances[0], -allowed[-1] - 0.12,
                    f"許容偏差の下限値 −{allowed[-1]:.1f} dB",
                    color="#c0392b", fontsize=8, va="top")
        bottom.axhline(0.0, color="#888888", linewidth=0.8)

        top.set_ylabel("音圧レベル [dB]")
        top.set_title(f"{project.room_label}／{project.condition_label}\n"
                      f"{trace}　音圧レベルの偏移", fontsize=11)
        bottom.set_title(f"{trace}　逆二乗測からの偏差", fontsize=11)
        bottom.set_ylabel("逆二乗測からの偏差 [dB]")
        bottom.set_xlabel("音源からの距離 [m]")
        # 参考（複素和）が入るときは、その振れ幅も収まるように広げる
        span = max(4.0, allowed[-1] * 1.6)
        if coherent is not None:
            span = max(span, float(np.nanmax(np.abs(
                coherent["deviation"][:, bands]))) * 1.1)
        bottom.set_ylim(-span, span)
        for axis in (top, bottom):
            axis.grid(True, which="both", alpha=0.25)
        if coherent is not None:
            top.plot([], [], "k--", linewidth=0.9,
                     label="参考：位相を含む複素和")
        top.legend(fontsize=7, ncol=2)
        figure.tight_layout()
        path = project.figure_path(f"逆二乗_{trace}.png")
        figure.savefig(path, dpi=140)
        plt.close(figure)
        written.append(path)
        if verbose:
            print(f"[逆二乗] 図（{trace}）: {path}")
    return written


def write_figure(project, result, bands=None, verbose=True):
    """**参照実装**：測線を横に並べた 1 枚もの（2026-08-26 に測線ごとへ移行）。"""
    return write_figures(project, result, bands=bands, verbose=verbose)


def run(project, verbose=True):
    """1 条件ぶん（表と図）。→ (表, 図)"""
    result = evaluate(project, verbose=verbose)
    if result is None:
        return None, None
    return (write_csv(project, result, verbose=verbose),
            write_figures(project, result, verbose=verbose))


def main(argv=None):
    import argparse

    import condition_table as ct

    parser = argparse.ArgumentParser(
        description="逆二乗則からのずれを測線ごとに出す（計算はやり直さない）")
    parser.add_argument("folder", help="プロジェクトのフォルダ")
    parser.add_argument("--sheet", default=None,
                        help="条件（シート名）。省略すると入っている条件を全部")
    args = parser.parse_args(argv)

    base = pj.Project.load(args.folder)
    sheets = ([args.sheet] if args.sheet
              else (ct.sheets(ct.path(base)) if ct.exists(base) else [None]))
    done = 0
    for sheet in sheets:
        project = pj.Project(base.folder,
                             **{k: getattr(base, k) for k in pj.DEFAULTS})
        if sheet:
            project.condition_sheet = sheet
        print(f"[逆二乗] 条件『{project.condition_label}』")
        table, figure = run(project)
        done += 1 if table else 0
    print(f"[逆二乗] {done} 条件ぶん書き出しました")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
