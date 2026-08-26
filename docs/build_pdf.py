# -*- coding: utf-8 -*-
"""`docs/` の Markdown 資料を **PDF** にする（2026-08-23 ユーザー要望）。

    py -3.10 docs\build_pdf.py                 # docs/*.md を全部 → docs/pdf/*.pdf
    py -3.10 docs\build_pdf.py 技術説明書.md   # 1 つだけ

**新しいライブラリは入れない。**手元にあるもので済ませてある。

- Markdown → HTML … このファイルの中の簡易変換（`markdown_to_html`）。
  docs で実際に使っている記法（見出し・表・箇条書き・コード・引用・数式）だけ扱う
- 数式 … **matplotlib の mathtext** で PNG に焼いて埋め込む（`requirements.txt`
  に元から入っている）。`$…$` は行の中、`$$…$$` は中央寄せの別行。
  ★mathtext が読めない書き方（`\\xrightarrow` など）は**そのまま等幅で出す**。
  黙って消えると式が抜けたことに気づけないので、枠を付けて分かるようにする
- HTML → PDF … **Edge か Chrome のヘッドレス印刷**（`--print-to-pdf`）。
  Windows なら標準で入っているので追加インストールが要らない

★**画像は data URI で埋め込む**（PDF 化のときに別ファイルを探しに行かせない）。
★mermaid の図はブラウザだけでは描けないので、**元のテキストを枠で囲んで出す**
  （GitHub では図として見えるので、PDF では読める形にしておくだけにする）。
"""

import base64
import glob
import html as html_module
import io
import os
import re
import subprocess
import sys

DOCS_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(DOCS_DIR, "pdf")

# ブラウザの探し場所（ヘッドレス印刷に使う）。見つかった順に使う
BROWSERS = [
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
]

# 数式を焼くときの解像度と色。本文の字と並べて浮かないように黒に近い色にする
MATH_DPI = 220
MATH_COLOR = "#111111"
# 数式の字の大きさ [pt]。本文 10.5pt に対し、mathtext の字は少し小さく見えるので
# 行の中は 11、別行は 14 にしてある。
# ★**表示の大きさは焼いた PNG の実寸から決める**（`height:1.15em` のように
#   決め打ちすると分数や総和のある式が縦につぶれる。2026-08-23）
MATH_SIZE_INLINE = 11.0
MATH_SIZE_DISPLAY = 14.0

# ★mathtext が知らない書き方の言い換え（2026-08-23）。
#   そのまま渡すと式が焼けず生の LaTeX が出てしまうので、意味の変わらない
#   範囲で置き換える。**置き換えられないもの（`pmatrix` など）は生で出す**
MATH_REPLACEMENTS = [
    (r"\lVert", r"\|"), (r"\rVert", r"\|"),
    (r"\lvert", "|"), (r"\rvert", "|"),
    (r"\textrm", r"\mathrm"), (r"\text", r"\mathrm"),
    (r"\xrightarrow", r"\rightarrow"),
    (r"\xleftarrow", r"\leftarrow"),
    (r"\bigl", r"\left"), (r"\bigr", r"\right"),
    (r"\Bigl", r"\left"), (r"\Bigr", r"\right"),
    (r"\nonumber", ""), (r"\!", ""), (r"\;", r"\,"),
]

STYLE = """
@page { size: A4; margin: 16mm 14mm 18mm 14mm; }
body { font-family: "Yu Gothic", "Meiryo", "MS Gothic", sans-serif;
       font-size: 10.5pt; line-height: 1.75; color: #111; }
h1 { font-size: 19pt; border-bottom: 2px solid #333; padding-bottom: 6px;
     margin: 0 0 14px 0; }
h2 { font-size: 15pt; border-bottom: 1px solid #bbb; padding-bottom: 4px;
     margin: 26px 0 10px 0; page-break-after: avoid; }
h3 { font-size: 12.5pt; margin: 20px 0 8px 0; page-break-after: avoid; }
h4 { font-size: 11pt; margin: 16px 0 6px 0; page-break-after: avoid; }
p { margin: 8px 0; }
ul, ol { margin: 8px 0 8px 0; padding-left: 22px; }
li { margin: 3px 0; }
code { font-family: "Consolas", "MS Gothic", monospace; font-size: 9.5pt;
       background: #f2f3f5; padding: 1px 4px; border-radius: 3px; }
pre { font-family: "Consolas", "MS Gothic", monospace; font-size: 9pt;
      background: #f6f7f9; border: 1px solid #dcdfe4; border-radius: 4px;
      padding: 9px 11px; overflow-x: auto; line-height: 1.5;
      white-space: pre-wrap; word-break: break-all; page-break-inside: avoid; }
pre code { background: none; padding: 0; font-size: 9pt; }
table { border-collapse: collapse; margin: 10px 0; font-size: 9.5pt;
        width: 100%; page-break-inside: avoid; }
th, td { border: 1px solid #c9ced6; padding: 4px 7px; text-align: left;
         vertical-align: top; }
th { background: #eceff3; font-weight: bold; }
blockquote { margin: 10px 0; padding: 6px 12px; border-left: 4px solid #9aa3b0;
             background: #f5f6f8; color: #333; }
blockquote p { margin: 4px 0; }
hr { border: none; border-top: 1px solid #ccc; margin: 18px 0; }
a { color: #1a4f8a; text-decoration: none; }
img.math-inline { vertical-align: -0.28em; }
div.math-block { text-align: center; margin: 12px 0; page-break-inside: avoid; }
span.math-raw { font-family: "Consolas", monospace; font-size: 9.5pt;
                background: #fff3d6; border: 1px solid #e0c070;
                padding: 0 4px; border-radius: 3px; }
div.mermaid-box { border: 1px dashed #9aa3b0; background: #fbfbfc;
                  padding: 8px 11px; margin: 10px 0; page-break-inside: avoid; }
div.mermaid-box .caption { font-size: 9pt; color: #666; margin-bottom: 4px; }
div.mermaid-box pre { border: none; background: none; padding: 0; }
"""


# ------------------------------------------------------------------------------
# 数式（matplotlib の mathtext で PNG にする）
# ------------------------------------------------------------------------------

def _png_size(png):
    """PNG のバイト列から (幅, 高さ) を読む（IHDR は先頭の決まった位置にある）。"""
    return (int.from_bytes(png[16:20], "big"), int.from_bytes(png[20:24], "big"))


def _math_png(latex, display=False):
    """LaTeX を PNG のバイト列にする。読めなければ None。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for before, after in MATH_REPLACEMENTS:
        latex = latex.replace(before, after)
    size = MATH_SIZE_DISPLAY if display else MATH_SIZE_INLINE
    figure = plt.figure(figsize=(0.01, 0.01))
    try:
        figure.text(0, 0, f"${latex}$", fontsize=size, color=MATH_COLOR)
        buffer = io.BytesIO()
        figure.savefig(buffer, format="png", dpi=MATH_DPI, transparent=True,
                       bbox_inches="tight", pad_inches=0.02)
        return buffer.getvalue()
    except Exception:
        return None                 # mathtext が読めない書き方
    finally:
        plt.close(figure)


def _math_tag(latex, display=False):
    """数式を <img>（焼けた場合）か等幅の生テキスト（焼けなかった場合）にする。"""
    png = _math_png(latex, display=display)
    if png is None:
        # ★黙って消さない。式が抜けたことが分かるように枠を付けて生で出す
        raw = html_module.escape(latex)
        if display:
            return f'<div class="math-block"><span class="math-raw">{raw}</span></div>'
        return f'<span class="math-raw">{raw}</span>'
    data = base64.b64encode(png).decode("ascii")
    # ★大きさは**焼いた PNG の実寸から**決める（MATH_DPI で焼いたので 72/dpi 倍が pt）
    height = _png_size(png)[1] * 72.0 / MATH_DPI
    if display:
        return (f'<div class="math-block">'
                f'<img src="data:image/png;base64,{data}" '
                f'style="height:{height:.1f}pt;max-width:100%"></div>')
    return (f'<img class="math-inline" src="data:image/png;base64,{data}" '
            f'style="height:{height:.1f}pt">')


# ------------------------------------------------------------------------------
# Markdown → HTML（docs で実際に使っている記法だけ）
# ------------------------------------------------------------------------------

INLINE_PATTERNS = [
    (re.compile(r"\*\*(.+?)\*\*", re.S), r"<strong>\1</strong>"),
    (re.compile(r"(?<![\w*])\*([^*\n]+?)\*(?![\w*])"), r"<em>\1</em>"),
    (re.compile(r"\[([^\]]+)\]\(([^)]+)\)"), r'<a href="\2">\1</a>'),
]


# 画像の base（`![図](…)` の相対パスをここから解決する）。
# `markdown_to_html(base_dir=…)` が入れ替える
IMAGE_BASE = [None]

# 画像の横幅（PDF の紙幅に収める）
IMAGE_STYLE = "max-width:100%;height:auto;display:block;margin:0.6em auto;"


def _image_tag(alt, source):
    """`![alt](path)` を **中身を埋め込んだ** <img> にする（1 枚の HTML にするため）。"""
    path = source.strip()
    if not os.path.isabs(path) and IMAGE_BASE[0]:
        path = os.path.join(IMAGE_BASE[0], path)
    try:
        with io.open(path, "rb") as handle:
            data = base64.b64encode(handle.read()).decode("ascii")
    except OSError:
        return (f'<span class="warn">画像が読めません: '
                f'{html_module.escape(source)}</span>')
    kind = "png" if path.lower().endswith(".png") else "jpeg"
    return (f'<img alt="{html_module.escape(alt)}" style="{IMAGE_STYLE}" '
            f'src="data:image/{kind};base64,{data}"/>')


def _inline(text):
    """段落の中の記法を HTML にする。**コードと数式は先に取り分ける。**"""
    kept = []

    def keep(fragment):
        kept.append(fragment)
        return f"\x00{len(kept) - 1}\x00"

    # ① `code` を守る（中の * や _ を記法として解釈させない）
    text = re.sub(r"`([^`]+)`",
                  lambda m: keep(f"<code>{html_module.escape(m.group(1))}</code>"),
                  text)
    # ★画像 `![alt](path)` を先に取り分ける（中身を埋め込む）
    # ★丸括弧を含むファイル名（`…(GW10K+32K+48K)_逆二乗_真上方向.png`）も拾えるように、
    #   **かっこの入れ子を 1 段だけ許す**（`[^)]+` だと最初の `)` で切れる）
    text = re.sub(r"!\[([^\]]*)\]\(((?:[^()]|\([^()]*\))+)\)",
                  lambda m: keep(_image_tag(m.group(1), m.group(2))), text)
    # ② $…$ の数式を守る（$ が 1 つだけの行は数式ではないので触らない）
    text = re.sub(r"(?<!\$)\$([^$\n]+?)\$(?!\$)",
                  lambda m: keep(_math_tag(m.group(1), display=False)), text)

    text = html_module.escape(text)
    for pattern, replacement in INLINE_PATTERNS:
        text = pattern.sub(replacement, text)
    text = text.replace("&lt;br/&gt;", "<br/>").replace("&lt;br&gt;", "<br/>")
    for index, fragment in enumerate(kept):
        text = text.replace(f"\x00{index}\x00", fragment)
    return text


def _table(lines):
    """`| a | b |` の並びを <table> にする（2 行目の区切りは読み飛ばす）。"""
    rows = [[cell.strip() for cell in line.strip().strip("|").split("|")]
            for line in lines]
    body = rows[2:] if len(rows) >= 2 and set("-: |") >= set("".join(rows[1])) \
        else rows[1:]
    out = ["<table>", "<thead><tr>"]
    out += [f"<th>{_inline(cell)}</th>" for cell in rows[0]]
    out.append("</tr></thead><tbody>")
    for row in body:
        out.append("<tr>" + "".join(f"<td>{_inline(cell)}</td>" for cell in row)
                   + "</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def _list_block(lines):
    """箇条書き（`- ` / `1. `。インデントで入れ子）を <ul>/<ol> にする。"""
    out, stack = [], []          # stack は (インデント幅, タグ)
    for line in lines:
        match = re.match(r"^(\s*)([-*+]|\d+\.)\s+(.*)$", line)
        if match is None:                       # 続きの行は前の項目にくっつける
            if out:
                out.append("<br/>" + _inline(line.strip()))
            continue
        indent, marker, text = len(match.group(1)), match.group(2), match.group(3)
        tag = "ol" if marker[:1].isdigit() else "ul"
        while stack and indent < stack[-1][0]:
            out.append(f"</{stack.pop()[1]}>")
        if not stack or indent > stack[-1][0]:
            stack.append((indent, tag))
            out.append(f"<{tag}>")
        out.append(f"<li>{_inline(text)}</li>")
    while stack:
        out.append(f"</{stack.pop()[1]}>")
    return "".join(out)


def markdown_to_html(text, base_dir=None):
    """docs の Markdown を HTML の本文にする。

    `base_dir` を渡すと `![図](相対パス)` をそこから解決する（報告書で使う）。
    """
    IMAGE_BASE[0] = base_dir
    lines = text.replace("\r\n", "\n").split("\n")
    out, index = [], 0
    while index < len(lines):
        line = lines[index]

        # ---- コードブロック（``` で囲む）----
        fence = re.match(r"^\s*```+\s*(\S*)\s*$", line)
        if fence:
            language, index = fence.group(1), index + 1
            block = []
            while index < len(lines) and not re.match(r"^\s*```+\s*$", lines[index]):
                block.append(lines[index])
                index += 1
            index += 1
            body = html_module.escape("\n".join(block))
            if language.lower() == "mermaid":
                # ★ブラウザだけでは図にできないので、読める形で残す
                out.append('<div class="mermaid-box">'
                           '<div class="caption">図（mermaid の記述。'
                           'GitHub 上では図として表示されます）</div>'
                           f'<pre><code>{body}</code></pre></div>')
            else:
                out.append(f"<pre><code>{body}</code></pre>")
            continue

        # ---- 別行の数式（$$ … $$）----
        if line.strip().startswith("$$"):
            block = [line.strip()[2:]]
            if not line.strip().endswith("$$") or len(line.strip()) <= 4:
                index += 1
                while index < len(lines) and "$$" not in lines[index]:
                    block.append(lines[index])
                    index += 1
                if index < len(lines):
                    block.append(lines[index].replace("$$", ""))
            else:
                block = [line.strip()[2:-2]]
            index += 1
            out.append(_math_tag(" ".join(b.strip() for b in block).strip(),
                                 display=True))
            continue

        # ---- 見出し ----
        heading = re.match(r"^(#{1,6})\s+(.*)$", line)
        if heading:
            level = len(heading.group(1))
            out.append(f"<h{level}>{_inline(heading.group(2))}</h{level}>")
            index += 1
            continue

        # ---- 水平線 ----
        if re.match(r"^\s*(-{3,}|\*{3,}|_{3,})\s*$", line):
            out.append("<hr/>")
            index += 1
            continue

        # ---- 表 ----
        if line.lstrip().startswith("|"):
            block = []
            while index < len(lines) and lines[index].lstrip().startswith("|"):
                block.append(lines[index])
                index += 1
            out.append(_table(block))
            continue

        # ---- 引用 ----
        if line.lstrip().startswith(">"):
            block = []
            while index < len(lines) and (lines[index].lstrip().startswith(">")
                                          or lines[index].strip()):
                if not lines[index].lstrip().startswith(">") and block:
                    break
                block.append(re.sub(r"^\s*>\s?", "", lines[index]))
                index += 1
            paragraphs = "".join(f"<p>{_inline(b)}</p>" for b in block if b.strip())
            out.append(f"<blockquote>{paragraphs}</blockquote>")
            continue

        # ---- 箇条書き ----
        if re.match(r"^\s*([-*+]|\d+\.)\s+", line):
            block = []
            while index < len(lines) and lines[index].strip() \
                    and not lines[index].lstrip().startswith("|") \
                    and not re.match(r"^\s*```", lines[index]) \
                    and not re.match(r"^#{1,6}\s", lines[index]):
                block.append(lines[index])
                index += 1
            out.append(_list_block(block))
            continue

        # ---- 段落 ----
        if not line.strip():
            index += 1
            continue
        block = []
        while index < len(lines) and lines[index].strip() \
                and not re.match(r"^\s*(```|\||>|#{1,6}\s|[-*+]\s|\d+\.\s)",
                                 lines[index]) \
                and not lines[index].strip().startswith("$$"):
            block.append(lines[index].strip())
            index += 1
        if block:
            out.append(f"<p>{_inline(' '.join(block))}</p>")
        else:
            index += 1
    return "\n".join(out)


def build_html(md_path):
    """Markdown 1 本を、画像まで埋め込んだ 1 枚の HTML にする。"""
    with io.open(md_path, encoding="utf-8") as f:
        text = f.read()
    title = os.path.splitext(os.path.basename(md_path))[0]
    body = markdown_to_html(
        text, base_dir=os.path.dirname(os.path.abspath(md_path)))
    return (f"<!doctype html><html lang=\"ja\"><head><meta charset=\"utf-8\">"
            f"<title>{html_module.escape(title)}</title>"
            f"<style>{STYLE}</style></head><body>{body}</body></html>")


# ------------------------------------------------------------------------------
# HTML → PDF（Edge / Chrome のヘッドレス印刷）
# ------------------------------------------------------------------------------

def find_browser():
    for path in BROWSERS:
        if os.path.exists(path):
            return path
    return None


def html_to_pdf(html_path, pdf_path, browser=None, verbose=True):
    """ヘッドレスのブラウザで印刷する。うまくいけば PDF のパスを返す。"""
    browser = browser or find_browser()
    if browser is None:
        raise RuntimeError("Edge か Chrome が見つかりません（PDF 化に使います）。"
                           "HTML はできているので、それをブラウザで開いて"
                           "「印刷 → PDF として保存」でも同じものが作れます")
    command = [browser, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
               f"--print-to-pdf={pdf_path}", f"file:///{html_path}"]
    result = subprocess.run(command, capture_output=True, timeout=180)
    if not os.path.exists(pdf_path):
        # 古い Edge は `--headless=new` を知らないので、素の --headless で試す
        command[1] = "--headless"
        result = subprocess.run(command, capture_output=True, timeout=180)
    if not os.path.exists(pdf_path):
        raise RuntimeError(f"PDF を作れませんでした: "
                           f"{result.stderr.decode('utf-8', 'replace')[:400]}")
    return pdf_path


def build(names=None, out_dir=OUT_DIR, keep_html=False, verbose=True):
    """`docs/*.md` を PDF にする。作ったパスの一覧を返す。"""
    if names:
        paths = [n if os.path.isabs(n) else os.path.join(DOCS_DIR, n) for n in names]
    else:
        paths = sorted(glob.glob(os.path.join(DOCS_DIR, "*.md")))
    os.makedirs(out_dir, exist_ok=True)
    browser = find_browser()
    if verbose:
        print(f"[pdf] ブラウザ: {browser}")

    written = []
    for md_path in paths:
        stem = os.path.splitext(os.path.basename(md_path))[0]
        html_path = os.path.join(out_dir, stem + ".html")
        pdf_path = os.path.join(out_dir, stem + ".pdf")
        with io.open(html_path, "w", encoding="utf-8", newline="") as f:
            f.write(build_html(md_path))
        try:
            html_to_pdf(os.path.abspath(html_path).replace("\\", "/"), pdf_path,
                        browser=browser, verbose=verbose)
            size = os.path.getsize(pdf_path) / 1024.0
            if verbose:
                print(f"[pdf] {stem}.pdf（{size:.0f} KB）")
            written.append(pdf_path)
        except Exception as error:
            print(f"[pdf] {stem} は PDF にできませんでした: {error}")
            written.append(html_path)
        finally:
            if not keep_html and os.path.exists(html_path) \
                    and os.path.exists(pdf_path):
                os.remove(html_path)
    if verbose:
        print(f"[pdf] {len(written)} 件 → {out_dir}")
    return written


if __name__ == "__main__":
    arguments = [a for a in sys.argv[1:] if not a.startswith("--")]
    build(arguments or None, keep_html="--keep-html" in sys.argv)
