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

# ISO 3745 の許容値（参考。バンド中心周波数 [Hz] → 許容偏差 [dB]）
TOLERANCE = ((630.0, 1.5), (5000.0, 1.0), (float("inf"), 1.5))

# 測線の色（図。3 本まではこの順）
TRACE_COLORS = ("#4c78a8", "#f58518", "#54a24b", "#e45756", "#72b7b2")


def tolerance_of(frequency):
    """その帯域の許容偏差 [dB]（ISO 3745。参考値）。"""
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


def evaluate(project, verbose=True):
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

    traces = {}
    for name, level in levels.items():
        info = points.get(name)
        if info is None or not np.isfinite(info["距離"]):
            continue
        traces.setdefault(info["レイヤ"] or "測線", []).append(
            (info["距離"], name, np.asarray(level, dtype=float)))

    result = {"frequencies": np.asarray(frequencies, dtype=float), "traces": {}}
    for trace, entries in traces.items():
        entries.sort(key=lambda item: item[0])
        distances = np.array([d for d, _, _ in entries], dtype=float)
        names = [n for _, n, _ in entries]
        block = np.array([values for _, _, values in entries], dtype=float)
        delta = np.empty_like(block)
        reference = np.empty(block.shape[1], dtype=float)
        for band in range(block.shape[1]):
            delta[:, band], reference[band] = deviations(distances, block[:, band])
        result["traces"][trace] = {
            "names": names, "distances": distances, "levels": block,
            "deviation": delta, "reference": reference,
            "worst": np.nanmax(np.abs(delta), axis=0),
        }
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
    rows.append(("参考", "許容偏差_dB(ISO 3745)", None,
                 [tolerance_of(v) for v in frequencies]))

    path = os.path.join(project.folder, pj.RESULT_DIR,
                        project.prefixed(FILE_NAME))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tb.write_sectioned_table(path, frequencies, rows, value_label="音源距離_m")
    if verbose:
        print(f"[逆二乗] 表を書き出しました: {path}")
    return path


def write_figure(project, result, bands=None, verbose=True):
    """測線ごとの図（上：音圧レベルと理想の線／下：ずれ）。"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    import plots as pl

    pl.use_japanese_font()
    frequencies = result["frequencies"]
    if bands is None:
        # ★要望は「63〜500 Hz」。データがそれより広ければ 500 Hz までに絞る
        bands = [k for k, v in enumerate(frequencies) if v <= 500.0] or \
            list(range(len(frequencies)))

    traces = list(result["traces"])
    if not traces:
        return None
    figure, axes = plt.subplots(2, len(traces), figsize=(5.2 * len(traces), 7.4),
                                squeeze=False)
    for column, trace in enumerate(traces):
        data = result["traces"][trace]
        top, bottom = axes[0][column], axes[1][column]
        for order, band in enumerate(bands):
            colour = TRACE_COLORS[order % len(TRACE_COLORS)]
            label = f"{frequencies[band]:.0f} Hz"
            top.plot(data["distances"], data["levels"][:, band], "o-",
                     color=colour, markersize=3, linewidth=1.2, label=label)
            bottom.plot(data["distances"], data["deviation"][:, band], "o-",
                        color=colour, markersize=3, linewidth=1.2, label=label)
        # ★理想の線は**どのバンドに合わせたか**を書く（当てはめた高さは
        #   バンドごとに違うので、1 本だけ引くと他のバンドがずれて見える）
        ideal = data["reference"][bands[0]] - 20.0 * np.log10(data["distances"])
        top.plot(data["distances"], ideal, "--", color="#888888", linewidth=1.0,
                 label=f"理想（1/r²・{frequencies[bands[0]]:.0f} Hz に合わせた線）")
        allowed = tolerance_of(frequencies[bands[-1]])
        bottom.axhspan(-allowed, allowed, color="#4c78a8", alpha=0.10)
        bottom.axhline(0.0, color="#888888", linewidth=0.8)
        top.set_xscale("log")
        bottom.set_xscale("log")
        top.set_title(f"{trace}")
        top.set_ylabel("音圧レベル [dB]")
        bottom.set_ylabel("逆二乗からのずれ [dB]")
        bottom.set_xlabel("音源からの距離 [m]")
        top.grid(True, which="both", alpha=0.25)
        bottom.grid(True, which="both", alpha=0.25)
        if column == 0:
            top.legend(fontsize=7)
    figure.suptitle(f"{project.room_label} / {project.condition_label}"
                    f"  逆二乗則からのずれ")
    figure.tight_layout()
    path = project.figure_path("逆二乗.png")
    figure.savefig(path, dpi=140)
    plt.close(figure)
    if verbose:
        print(f"[逆二乗] 図を書き出しました: {path}")
    return path


def run(project, verbose=True):
    """1 条件ぶん（表と図）。→ (表, 図)"""
    result = evaluate(project, verbose=verbose)
    if result is None:
        return None, None
    return (write_csv(project, result, verbose=verbose),
            write_figure(project, result, verbose=verbose))


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
