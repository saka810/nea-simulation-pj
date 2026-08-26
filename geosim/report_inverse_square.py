# -*- coding: utf-8 -*-
"""逆二乗の**検討報告書**（Markdown ＋ PDF）を作る。

★2026-08-26 ユーザー要望

> 直下に検討報告書のサンプルを入れました。もし可能であれば、このような体裁に
> 合わせて欲しいです。逆二乗は 500Hz までで良いです。
> ちなみに、グラフの雰囲気をまねしなくて良いです。このようなグラフが欲しいと
> いうだけです。また、追加で、逆二乗の様子（音圧レベルの偏移）も載せてください。

サンプル（`無響室シミュレーション検討書.docx`）の並びに合わせてある:

  概要 → 目標音響性能（JIS Z 8732 の許容偏差の表）→ 計算条件 →
  吸音材の吸音率 → 測定方向（図）→ 方向ごとの結果（音圧レベルの偏移＋偏差）→ まとめ

決めごと:
  ・**計算はやり直さない**。`inverse_square.evaluate()` が読んだ結果を並べるだけ
  ・PDF 化は `docs/build_pdf.py` を使い回す（**追加のライブラリを入れない**方針のまま。
    数式は画像に焼き、印刷は Edge / Chrome のヘッドレス）
  ・図は**測線ごとに 1 枚**（`inverse_square.write_figures`）
  ・**判定は書くが、合否の断定はしない**（計算値なので「目標を満たす／外れる」まで）
"""
import io
import os
import sys

import numpy as np

import inverse_square as iq
import project as pj

FILE_NAME = "逆二乗_検討報告書"          # 拡張子は .md / .pdf

# 報告書に載せる帯域の範囲 [Hz]（★ユーザー指定：100〜500）
LOW, HIGH = 100.0, 500.0

# 目標を判定するときに参考として外すバンド（サンプルでは 100 Hz が対象外）
REFERENCE_ONLY = (100.0,)


def _fmt(value, digits=2):
    return "—" if value is None or not np.isfinite(value) else f"{value:.{digits}f}"


def _tolerance_table():
    lines = ["| 1/3 オクターブバンド中心周波数 [Hz] | 許容偏差 [dB] |",
             "|---|---|"]
    for limit, value in iq.TOLERANCE:
        if limit == 630.0:
            lines.append(f"| ≦ 630 | ±{value:.1f} |")
        elif limit == 5000.0:
            lines.append(f"| 800 〜 5000 | ±{value:.1f} |")
        else:
            lines.append(f"| ≧ 6300 | ±{value:.1f} |")
    lines.append("")
    lines.append(f"{iq.TOLERANCE_STANDARD} に準ずる")
    return "\n".join(lines)


def _absorption_table(project):
    """条件表から、この条件で使っている材料の吸音率を拾う。"""
    import condition_table as ct
    import run_project as rp

    try:
        library = rp._library_for(project)
        assignment = ct.assignment_for(project, verbose=False)
    except Exception:
        return None
    if not library or not assignment:
        return None
    frequencies = None
    rows = []
    for layer, material in assignment.items():
        entry = library.get(material) if hasattr(library, "get") else None
        if entry is None:
            continue
        # `MaterialLibrary` は `Material`（名前・吸音率・種類）を返す
        values = getattr(entry, "coefficients", entry)
        kind = getattr(entry, "kind", None)
        name = getattr(entry, "name", material)
        values = np.atleast_1d(np.asarray(values, dtype=float))
        rows.append((layer, f"{name}"
                     + (f"（{ct.KIND_TEXT.get(kind, kind)}）" if kind else ""),
                     values))
    if not rows:
        return None
    frequencies = getattr(library, "frequencies", None)
    return frequencies, rows


def _judgement(result, bands, frequencies):
    """測線ごとに「許容偏差に収まっているか」を見る。→ [(測線, 収まったか, 最大, バンド)]"""
    verdict = []
    for trace, data in result["traces"].items():
        worst, worst_band, inside = 0.0, None, True
        for band in bands:
            if frequencies[band] in REFERENCE_ONLY:
                continue
            value = float(data["worst"][band])
            allowed = iq.tolerance_of(frequencies[band])
            if value > worst:
                worst, worst_band = value, frequencies[band]
            if value > allowed:
                inside = False
        verdict.append((trace, inside, worst, worst_band))
    return verdict


def build_markdown(project, result, figures, bands=None, verbose=True):
    """報告書の Markdown を組み立てる。→ 文字列"""
    frequencies = result["frequencies"]
    if bands is None:
        bands = iq.bands_in_range(frequencies, LOW, HIGH)
    used = [frequencies[b] for b in bands]
    figure_dir = os.path.dirname(figures[0]) if figures else project.figure_dir()

    def relative(path):
        return os.path.relpath(path, os.path.join(project.folder, pj.RESULT_DIR)) \
            .replace("\\", "/")

    lines = [f"# {project.room_label}　逆二乗測特性シミュレーション報告書", ""]
    lines += [f"条件：**{project.condition_label}**", ""]

    lines += ["## 概要", "",
              f"{project.room_label}（半無響室）に設置する吸音材について、"
              f"室内の逆二乗測特性を計算で推定した。"
              f"音源を床面中央に置いたと仮定し、"
              f"{len(result['traces'])} 方向の測線について"
              f"逆二乗測からの偏差を算出している。", ""]

    lines += ["## 目標音響性能", "",
              f"- 半自由音場　{used[0]:.0f}〜{used[-1]:.0f} Hz（1/3 oct.）",
              "- 逆二乗測からの偏差", "", _tolerance_table(), ""]
    if REFERENCE_ONLY:
        lines += [f"（{'・'.join(f'{v:.0f}' for v in REFERENCE_ONLY)} Hz の"
                  f"計算結果については参考として記載）", ""]

    lines += ["## 計算条件", "", "| 項目 | 値 |", "|---|---|"]
    model = _model_note(project)
    for key, value in model:
        lines.append(f"| {key} | {value} |")
    lines.append("")

    absorption = _absorption_table(project)
    if absorption is not None:
        centres, rows = absorption
        centres = frequencies if centres is None else np.asarray(centres)
        lines += ["## 吸音材の吸音率（計算に使った値）", "",
                  "| 面 | 材料 | " + " | ".join(f"{v:.0f} Hz" for v in used) + " |",
                  "|---|---|" + "---|" * len(used)]
        for layer, material, values in rows:
            cells = []
            for band in bands:
                cells.append(f"{values[band]:.2f}" if band < len(values) else "—")
            lines.append(f"| {layer} | {material} | " + " | ".join(cells) + " |")
        lines.append("")

    points = project.existing_result_path("points")
    figure_points = project.figure_path("points.png", shared=True)
    if os.path.exists(figure_points):
        lines += ["## 測定方向", "",
                  f"![測定方向]({relative(figure_points)})", "",
                  "図-1　室寸法と逆二乗測の計算方向"
                  "（左：平面図／中・右：立面図）", ""]

    lines += ["## 方向ごとの結果", ""]
    for index, (trace, data) in enumerate(result["traces"].items(), start=1):
        figure = next((f for f in figures if trace in os.path.basename(f)), None)
        lines += [f"### {trace}", ""]
        if figure:
            lines += [f"![{trace}]({relative(figure)})", "",
                      f"図-{index + 1}　{trace}　"
                      f"上：音圧レベルの偏移（半自由音場の理論線つき）／"
                      f"下：逆二乗測からの偏差（{iq.TOLERANCE_STANDARD} の許容偏差つき）",
                      ""]
        lines += ["| 帯域 [Hz] | " + " | ".join(f"{v:.0f}" for v in used) + " |",
                  "|---|" + "---|" * len(used),
                  "| 最大の偏差 [dB] | "
                  + " | ".join(_fmt(data["worst"][b]) for b in bands) + " |",
                  "| 許容偏差 [dB] | "
                  + " | ".join(f"±{iq.tolerance_of(frequencies[b]):.1f}"
                               for b in bands) + " |", ""]

    verdict = _judgement(result, bands, frequencies)
    lines += ["## まとめ", "", "| 方向 | 最大の偏差 [dB] | 判定 |", "|---|---|---|"]
    for trace, inside, worst, band in verdict:
        note = "許容偏差の内側" if inside else "許容偏差を外れる帯域あり"
        where = f"{worst:.2f}（{band:.0f} Hz）" if band else _fmt(worst)
        lines.append(f"| {trace} | {where} | {note} |")
    lines.append("")
    if all(inside for _, inside, _, _ in verdict):
        lines.append("いずれの方向についても、目標音響性能を満たす結果となった。")
    else:
        outside = [t for t, inside, _, _ in verdict if not inside]
        lines.append(f"**{'・'.join(outside)}** で許容偏差を外れる帯域がある。")
    lines += ["",
              f"（{'・'.join(f'{v:.0f}' for v in REFERENCE_ONLY)} Hz は"
              f"目標値の対象周波数外として判定から除いている）" if REFERENCE_ONLY else "",
              "",
              "※ 本書はシミュレーションの結果であり、実測による確認に代わるものではない。",
              ""]
    return "\n".join(lines)


def _model_note(project):
    """計算条件の表に載せる項目。"""
    import read_dxffile as rd

    rows = []
    try:
        model = rd.read_model(project.dxf_path, unit=project.unit,
                              band_number=project.band_number, verbose=False)
        low, high = model.extents
        size = np.asarray(high) - np.asarray(low)
        rows.append(("室の寸法", f"{size[0]:.2f} × {size[1]:.2f} × {size[2]:.2f} m"
                                 f"（{model.volume:.1f} m³）"))
        if model.source_points:
            point = model.source_points[0]
            rows.append(("音源", f"[{point[0]:.2f}, {point[1]:.2f}, {point[2]:.2f}] m"
                                 f"（床面上）"))
    except Exception:
        pass
    import absorption as ab
    centres = ab.frequency_bands(project.band_number,
                                 getattr(project, "band_width", "1/1"),
                                 getattr(project, "band_start", None))
    rows.append(("帯域", f"{getattr(project, 'band_width', '1/1')} オクターブ　"
                         f"{centres[0]:.0f}〜{centres[-1]:.0f} Hz"))
    rows.append(("音線", f"{project.rays:,} 本／最大反射 {project.nref} 回／"
                         f"受音球 {project.radius:.2f} m"))
    rows.append(("大気", f"{project.temperature:.1f} ℃／湿度 {project.humidity:.0f} %／"
                         f"{project.pressure:.3f} kPa"))
    rows.append(("音源のパワーレベル",
                 "未入力（相対値。Lw = 0 dB として計算）"
                 if project.source_power_db is None else f"{project.source_power_db}"))
    return rows


def build(project, verbose=True, pdf=True):
    """報告書（Markdown と PDF）を作る。→ (md, pdf)"""
    result = iq.evaluate(project, verbose=verbose)
    if result is None:
        return None, None
    bands = iq.bands_in_range(result["frequencies"], LOW, HIGH)
    figures = iq.write_figures(project, result, bands=bands, verbose=verbose)
    iq.write_csv(project, result, verbose=verbose)

    text = build_markdown(project, result, figures, bands=bands, verbose=verbose)
    md_path = os.path.join(project.folder, pj.RESULT_DIR,
                           project.prefixed(FILE_NAME + ".md"))
    os.makedirs(os.path.dirname(md_path), exist_ok=True)
    with io.open(md_path, "w", encoding="utf-8") as handle:
        handle.write(text)
    if verbose:
        print(f"[報告書] {md_path}")

    pdf_path = None
    if pdf:
        try:
            pdf_path = _to_pdf(md_path, verbose=verbose)
        except Exception as error:
            print(f"[報告書] PDF にはできませんでした（{type(error).__name__}: {error}）。"
                  f"Markdown はできています")
    return md_path, pdf_path


def _to_pdf(md_path, verbose=True):
    """`docs/build_pdf.py` を使い回して PDF にする。"""
    here = os.path.dirname(os.path.abspath(__file__))
    docs = os.path.join(os.path.dirname(here), "docs")
    if docs not in sys.path:
        sys.path.insert(0, docs)
    import build_pdf

    html_path = os.path.splitext(md_path)[0] + ".html"
    with io.open(html_path, "w", encoding="utf-8") as handle:
        handle.write(build_pdf.build_html(md_path))
    pdf_path = os.path.splitext(md_path)[0] + ".pdf"
    build_pdf.html_to_pdf(html_path, pdf_path, verbose=verbose)
    os.remove(html_path)
    if verbose:
        print(f"[報告書] {pdf_path}")
    return pdf_path


def main(argv=None):
    import argparse

    import condition_table as ct

    parser = argparse.ArgumentParser(
        description="逆二乗の検討報告書を作る（計算はやり直さない）")
    parser.add_argument("folder", help="プロジェクトのフォルダ")
    parser.add_argument("--sheet", default=None, help="条件（シート名）")
    parser.add_argument("--no-pdf", action="store_true")
    args = parser.parse_args(argv)

    base = pj.Project.load(args.folder)
    sheets = ([args.sheet] if args.sheet
              else (ct.sheets(ct.path(base)) if ct.exists(base) else [None]))
    for sheet in sheets:
        project = pj.Project(base.folder,
                             **{k: getattr(base, k) for k in pj.DEFAULTS})
        if sheet:
            project.condition_sheet = sheet
        print(f"[報告書] 条件『{project.condition_label}』")
        build(project, pdf=not args.no_pdf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
