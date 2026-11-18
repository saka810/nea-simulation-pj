"""社内説明資料_標準テンプレート.pptx から A4 横・提出用書体のひな形を作る。

    py -3.13 資料\\テンプレート\\powerpoint\\make_a4_template.py

出力は同じフォルダの 提出資料_A4テンプレート.pptx（見本 4 ページ）。
用紙 A4 横（297×210mm）、和文 ＭＳ ゴシック・英数字 Times New Roman（テーマにも設定）。
帯・カード・背景は標準テンプレートの図形をそのまま複製するので、影や角丸も引き継ぐ。
"""
import copy
from pathlib import Path

from lxml import etree
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Emu, Inches, Pt

HERE = Path(__file__).parent
SOURCE = HERE / "社内説明資料_標準テンプレート.pptx"
OUT = HERE / "提出資料_A4テンプレート.pptx"

LATIN = "Times New Roman"
EA = "ＭＳ ゴシック"

TEAL = "0F766E"
TEAL_DARK = "0F4C45"
HEAD_SUB = "D5EEEA"
INK = "2B2F36"
SUB = "6B7280"
GRID = "D0D5DD"

prs = Presentation(SOURCE)
prs.slide_width, prs.slide_height = Emu(297 * 36000), Emu(210 * 36000)
SW = prs.slide_width / 914400  # 11.69 in

# 標準テンプレートのページから部品の見本を取り、元のページは外す
src = list(prs.slides)
PROTO = {key: copy.deepcopy(shape._element) for key, shape in {
    "cover_band": src[0].shapes[0], "head_band": src[1].shapes[0],
    "text": src[1].shapes[2], "card": src[1].shapes[4], "tint": src[1].shapes[6],
}.items()}
BG = copy.deepcopy(src[1]._element.cSld.bg)
for sid in list(prs.slides._sldIdLst):
    prs.part.drop_rel(sid.rId)
    prs.slides._sldIdLst.remove(sid)

# テーマの書体（後から足した文字にも効くように）
for part in prs.part.package.iter_parts():
    if "/theme/" in str(part.partname):
        root = etree.fromstring(part.blob)
        for font in root.iter(qn("a:majorFont"), qn("a:minorFont")):
            font.find(qn("a:latin")).set("typeface", LATIN)
            font.find(qn("a:ea")).set("typeface", EA)
            for script in font.findall(qn("a:font")):
                if script.get("script") == "Jpan":
                    script.set("typeface", EA)
        part._blob = etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)

BLANK = prs.slide_layouts[6]
TOTAL = 4


def set_font(rpr, size, bold, color):
    rpr.set("sz", str(int(size * 100)))
    rpr.set("b", "1" if bold else "0")
    for tag in ("a:solidFill", "a:latin", "a:ea", "a:cs"):
        for e in rpr.findall(qn(tag)):
            rpr.remove(e)
    fill = etree.SubElement(rpr, qn("a:solidFill"))
    etree.SubElement(fill, qn("a:srgbClr")).set("val", color)
    etree.SubElement(rpr, qn("a:latin")).set("typeface", LATIN)
    etree.SubElement(rpr, qn("a:ea")).set("typeface", EA)
    etree.SubElement(rpr, qn("a:cs")).set("typeface", LATIN)


def write(tf, paras, size=12, color=INK, bold=False, align=PP_ALIGN.LEFT, space_after=0):
    for p in list(tf.paragraphs)[1:]:
        p._p.getparent().remove(p._p)
    for i, line in enumerate(paras):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        for r in list(p.runs):
            r._r.getparent().remove(r._r)
        p.alignment = align
        if space_after:
            p.space_after = Pt(space_after)
        r = p.add_run()
        r.text = line
        set_font(r._r.get_or_add_rPr(), size, bold, color)


def clone(slide, key, x, y, w, h):
    el = copy.deepcopy(PROTO[key])
    nv = el.find(".//" + qn("p:cNvPr"))
    for ext in nv.findall(qn("a:extLst")):
        nv.remove(ext)
    ids = [int(e.get("id")) for e in slide.shapes._spTree.iter(qn("p:cNvPr"))]
    nv.set("id", str(max(ids + [1]) + 1))
    slide.shapes._spTree.append(el)
    shape = slide.shapes[-1]
    shape.left, shape.top, shape.width, shape.height = Inches(x), Inches(y), Inches(w), Inches(h)
    return shape


def text(slide, x, y, w, h, paras, anchor=MSO_ANCHOR.TOP, **kw):
    shape = clone(slide, "text", x, y, w, h)
    tf = shape.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    write(tf, paras, **kw)
    return shape


def footer(slide, no):
    text(slide, 0.4, 7.88, 6.5, 0.24, ["案件名・資料名"], size=9, color=SUB)
    text(slide, SW - 3.4, 7.88, 3.0, 0.24, [f"会社名　{no:02d} / {TOTAL:02d}"], size=9, color=SUB, align=PP_ALIGN.RIGHT)


def page(no, label, title, subtitle):
    slide = prs.slides.add_slide(BLANK)
    slide._element.cSld.insert(0, copy.deepcopy(BG))
    clone(slide, "head_band", 0, 0, SW, 0.82)
    text(slide, 0.4, 0.13, 1.0, 0.28, [label], size=12, color="FFFFFF", bold=True)
    text(slide, 1.45, 0.07, SW - 1.85, 0.40, [title], size=19, color="FFFFFF", bold=True, anchor=MSO_ANCHOR.MIDDLE)
    text(slide, 1.45, 0.48, SW - 1.85, 0.26, [subtitle], size=11, color=HEAD_SUB)
    footer(slide, no)
    return slide


def table(slide, x, y, rows, col_w, row_h=0.42, size=12):
    shape = slide.shapes.add_table(len(rows), len(rows[0]), Inches(x), Inches(y),
                                   Inches(sum(col_w)), Inches(row_h * len(rows)))
    tbl = shape.table
    tbl_pr = tbl._tbl.tblPr
    for k in list(tbl_pr.attrib):
        del tbl_pr.attrib[k]  # 既定の表スタイル（縞・太字）を外す
    for e in list(tbl_pr):
        tbl_pr.remove(e)
    for i, w in enumerate(col_w):
        tbl.columns[i].width = Inches(w)
    for r, row in enumerate(rows):
        tbl.rows[r].height = Inches(row_h)
        head = r == 0
        for c, value in enumerate(row):
            cell = tbl.cell(r, c)
            cell.fill.solid()
            cell.fill.fore_color.rgb = RGBColor.from_string(TEAL if head else "FFFFFF")
            cell.margin_left = cell.margin_right = Inches(0.06)
            cell.vertical_anchor = MSO_ANCHOR.MIDDLE
            write(cell.text_frame, [value], size=size, bold=head,
                  color="FFFFFF" if head else INK, align=PP_ALIGN.CENTER)
            tc_pr = cell._tc.get_or_add_tcPr()
            for i, tag in enumerate(("a:lnL", "a:lnR", "a:lnT", "a:lnB")):
                ln = etree.Element(qn(tag), w="6350")
                fill = etree.SubElement(ln, qn("a:solidFill"))
                etree.SubElement(fill, qn("a:srgbClr")).set("val", TEAL if head else GRID)
                tc_pr.insert(i, ln)  # 線は塗りより前（スキーマの順）


# 1. 表紙
s = prs.slides.add_slide(BLANK)
s._element.cSld.insert(0, copy.deepcopy(BG))
clone(s, "cover_band", 0, 0, SW, 3.6)
text(s, 0.75, 1.0, 10.2, 0.3, ["資料の種類（例：騒音測定　比較資料）"], size=14, color=HEAD_SUB, bold=True)
text(s, 0.75, 1.42, 10.2, 1.2, ["案件名", "資料のタイトル"], size=28, color="FFFFFF", bold=True)
text(s, 0.75, 2.75, 10.2, 0.35, ["サブタイトル・この資料で示すこと"], size=15, color=HEAD_SUB)
text(s, 0.75, 4.15, 10.2, 1.3, ["条件・出典など（例：測定日、計算条件）", "2 行目"], size=13, space_after=6)
text(s, 0.75, 6.3, 10.2, 0.75, ["2026年○月", "会社名"], size=14, color="334155", space_after=4)
footer(s, 1)

# 2. 2 列のカード＋まとめ
s = page(2, "分類", "このページで伝えたいこと", "結論を補足する説明")
for x in (0.4, 5.94):
    clone(s, "card", x, 1.05, 5.35, 4.75)
    text(s, x + 0.25, 1.25, 4.85, 0.35, ["カードの見出し"], size=14, color=TEAL_DARK, bold=True)
    text(s, x + 0.25, 1.8, 4.85, 3.8, ["本文 12〜13pt。根拠・条件・具体例を書く", "2 段落目"], size=12.5, space_after=5)
clone(s, "tint", 0.4, 6.0, 10.89, 1.65)
text(s, 0.65, 6.18, 10.4, 1.35, ["・ページのまとめ・注記（12pt）", "・2 行目"], size=12, color=TEAL_DARK, space_after=4)

# 3. 図＋説明
s = page(3, "結果", "図・グラフと説明のページ", "左に図（ネイティブのグラフ・図・写真）、右に要点")
clone(s, "card", 0.4, 1.05, 7.65, 6.6)
text(s, 0.6, 4.1, 7.25, 0.4, ["ここに図・グラフを挿入"], size=16, color=SUB, bold=True, align=PP_ALIGN.CENTER)
clone(s, "card", 8.25, 1.05, 3.04, 6.6)
text(s, 8.45, 1.25, 2.7, 0.35, ["要点"], size=14, color=TEAL_DARK, bold=True)
text(s, 8.45, 1.75, 2.7, 5.7, ["数値・読み方・注意点", "2 段落目"], size=11.5, space_after=5)

# 4. 表＋まとめ
s = page(4, "比較", "表で比べるページ", "表は PowerPoint 上で編集できる形式")
clone(s, "card", 0.4, 1.05, 10.89, 4.75)
table(s, 0.65, 1.3, [["項目", "案A", "案B", "案C"]] + [[f"項目{i}", "内容", "内容", "内容"] for i in (1, 2, 3, 4)],
      [2.4, 2.66, 2.66, 2.66])
clone(s, "tint", 0.4, 6.0, 10.89, 1.65)
text(s, 0.65, 6.18, 10.4, 1.35, ["・単位・凡例・注記（表の数値は 10〜12pt）"], size=12, color=TEAL_DARK)

prs.save(OUT)
print("保存:", OUT)
