# 7 枚版のコピーに「社内説明資料_標準テンプレート」の見た目（緑の帯・灰色の地・角丸カード・BIZ UDPGothic）を当てる
# 使い方：PowerPoint で対象のコピーを開いた状態で  powershell -File apply_template.ps1 -Name "<ファイル名の一部>"
param([string]$Name = "標準テンプレート")

function RGBc([string]$h) { [int]("0x" + $h.Substring(0, 2)) + [int]("0x" + $h.Substring(2, 2)) * 256 + [int]("0x" + $h.Substring(4, 2)) * 65536 }
$FONT = "BIZ UDPGothic"
# 旧配色 → テンプレートの配色（塗り）
$fillMap = @{
    (RGBc "1F2A36") = (RGBc "0F766E")   # 紺 → 帯の緑
    (RGBc "2B3A4A") = (RGBc "0F4C45")   # 紺2 → 濃い緑
    (RGBc "EEF2F5") = (RGBc "E3F2F0")   # 薄灰 → 薄緑
    (RGBc "FBEBE1") = (RGBc "E3F2F0")   # 薄橙の帯 → まとめカード
    (RGBc "E07A3F") = (RGBc "0F4C45")   # 橙の印 → 濃い緑
}
# 旧配色 → テンプレートの配色（文字）
$textMap = @{
    (RGBc "1F2A36") = (RGBc "2B2F36")
    (RGBc "2B3A4A") = (RGBc "2B2F36")
    (RGBc "E07A3F") = (RGBc "0F766E")
    (RGBc "8A99A8") = (RGBc "6B7280")
}
$lineMap = @{ (RGBc "D5DCE3") = (RGBc "D0D5DD"); (RGBc "1F2A36") = (RGBc "0F766E"); (RGBc "E07A3F") = (RGBc "0F766E") }

function Fix-Text($tr) {
    $tr.Font.Name = $FONT; $tr.Font.NameFarEast = $FONT
    foreach ($run in $tr.Runs()) {
        $c = $run.Font.Color.RGB
        if ($textMap.ContainsKey($c)) { $run.Font.Color.RGB = $textMap[$c] }
    }
}
function Fix-Shape($sh) {
    if ($sh.Type -eq 6) { foreach ($g in $sh.GroupItems) { Fix-Shape $g }; return }   # グループ
    if ($sh.HasTable) {
        $tb = $sh.Table
        for ($r = 1; $r -le $tb.Rows.Count; $r++) { for ($c = 1; $c -le $tb.Columns.Count; $c++) {
            $cs = $tb.Cell($r, $c).Shape
            try { $f = $cs.Fill.ForeColor.RGB; if ($cs.Fill.Visible -and $fillMap.ContainsKey($f)) { $cs.Fill.ForeColor.RGB = $fillMap[$f] } } catch {}
            if ($cs.TextFrame.HasText) { Fix-Text $cs.TextFrame.TextRange; $cs.TextFrame.TextRange.Font.Size = [Math]::Round($cs.TextFrame.TextRange.Font.Size * 0.92 * 2) / 2 }   # BIZ UDPGothic は幅が広いので少し小さく
        } }
        return
    }
    if ($sh.HasChart) {
        $ch = $sh.Chart
        # ChartArea にまとめて書体を当てると文字の大きさまで揃ってしまうので、要素ごとに名前だけ替える
        function Set-TF2($tf2) { try { $f = $tf2.TextRange.Font; $f.Name = $FONT; $f.NameFarEast = $FONT } catch {} }
        try { if ($ch.HasTitle) { Set-TF2 $ch.ChartTitle.Format.TextFrame2 } } catch {}
        try { if ($ch.HasLegend) { $ch.Legend.Font.Name = $FONT } } catch {}
        foreach ($ai in 1, 2) {
            try {
                if ($ch.HasAxis($ai)) {
                    $ax = $ch.Axes($ai); $ax.TickLabels.Font.Name = $FONT
                    if ($ax.HasTitle) { Set-TF2 $ax.AxisTitle.Format.TextFrame2 }
                }
            } catch {}
        }
        foreach ($se in $ch.SeriesCollection()) {
            try { $lc = $se.Format.Line.ForeColor.RGB; if ($lineMap.ContainsKey($lc)) { $se.Format.Line.ForeColor.RGB = $lineMap[$lc]; $se.MarkerBackgroundColor = $lineMap[$lc]; $se.MarkerForegroundColor = $lineMap[$lc] } } catch {}
            try { $fc = $se.Format.Fill.ForeColor.RGB; if ($fillMap.ContainsKey($fc)) { $se.Format.Fill.ForeColor.RGB = $fillMap[$fc] } } catch {}
        }
        return
    }
    try { if ($sh.Fill.Visible) { $f = $sh.Fill.ForeColor.RGB; if ($fillMap.ContainsKey($f)) { $sh.Fill.ForeColor.RGB = $fillMap[$f] } } } catch {}
    try { if ($sh.Line.Visible) { $l = $sh.Line.ForeColor.RGB; if ($lineMap.ContainsKey($l)) { $sh.Line.ForeColor.RGB = $lineMap[$l] } } } catch {}
    if ($sh.HasTextFrame -and $sh.TextFrame.HasText) { Fix-Text $sh.TextFrame.TextRange }
    # 塗りのある四角（細い印・帯は除く）はテンプレートのカードと同じ角丸＋影に
    if ($sh.Type -eq 1 -and $sh.AutoShapeType -eq 1 -and $sh.Fill.Visible -and $sh.Width -ge 40 -and $sh.Height -ge 30) {
        $sh.AutoShapeType = 5
        $sh.Adjustments.Item(1) = [Math]::Min(0.5, 14 / [Math]::Min($sh.Width, $sh.Height))
        $sh.Shadow.Visible = -1; $sh.Shadow.ForeColor.RGB = 0; $sh.Shadow.Blur = 3.15
        $sh.Shadow.OffsetX = 0; $sh.Shadow.OffsetY = 1.81; $sh.Shadow.Transparency = 0.65
    }
}

$pp = [Runtime.InteropServices.Marshal]::GetActiveObject("PowerPoint.Application")
$p = $pp.Presentations | Where-Object { $_.Name -like "*$Name*" }
if (-not $p) { throw "開いている資料に '$Name' がありません" }
$n = $p.Slides.Count
$SHIFT = 22   # 帯が細くなったぶん、本文を上へ
foreach ($s in $p.Slides) {
    $s.FollowMasterBackground = 0; $s.Background.Fill.Solid(); $s.Background.Fill.ForeColor.RGB = (RGBc "F2F4F7")
    $s.DisplayMasterShapes = 0
    $items = @(); foreach ($sh in $s.Shapes) { $items += $sh }
    foreach ($sh in $items) {
        $isPh = ($sh.Type -eq 14)
        $txt = ""; if ($sh.HasTextFrame -and $sh.TextFrame.HasText) { $txt = $sh.TextFrame.TextRange.Text }
        if (-not $isPh -and $sh.Top -ge 495 -and $txt -match '^\s*\d+\s*$') { $sh.Delete(); continue }   # 旧ページ番号
        if (-not $isPh -and $sh.Top -ge 495) {   # 下端の注記はページ番号と重ならないよう左下へ
            $sh.Left = 32; $sh.Width = 720; $sh.TextFrame.TextRange.ParagraphFormat.Alignment = 1
        }
        if (-not $isPh -and $sh.Top -lt 495) { $sh.Top = $sh.Top - $SHIFT }
    }
    # 角丸にするカードに貼り付いた細い印（左端・上端の線）は角と合わないので消す
    $cards = @(); foreach ($sh in $s.Shapes) { if ($sh.Type -eq 1 -and $sh.AutoShapeType -eq 1 -and $sh.Fill.Visible -and $sh.Width -ge 40 -and $sh.Height -ge 30) { $cards += $sh } }
    $thin = @(); foreach ($sh in $s.Shapes) { if ($sh.Type -eq 1 -and $sh.AutoShapeType -eq 1 -and ($sh.Width -le 10 -or $sh.Height -le 10)) { $thin += $sh } }
    foreach ($t in $thin) {
        foreach ($c in $cards) {
            if ($t.Left -lt $c.Left + $c.Width + 2 -and $t.Left + $t.Width -gt $c.Left - 2 -and $t.Top -lt $c.Top + $c.Height + 2 -and $t.Top + $t.Height -gt $c.Top - 2) { $t.Delete(); break }
        }
    }
    foreach ($sh in $s.Shapes) { Fix-Shape $sh }
    # 見出しの帯
    $band = $s.Shapes.AddShape(1, 0, 0, 960, 66); $band.Name = "TplBand"; $band.Line.Visible = 0
    $band.Fill.ForeColor.RGB = (RGBc "0F766E")
    $band.Shadow.Visible = -1; $band.Shadow.ForeColor.RGB = 0; $band.Shadow.Blur = 15; $band.Shadow.OffsetX = 0; $band.Shadow.OffsetY = 3; $band.Shadow.Transparency = 0.86
    $band.ZOrder(1)   # 最背面
    foreach ($sh in $s.Shapes) {
        if ($sh.Type -ne 14) { continue }
        $tf = $sh.TextFrame; $tf.AutoSize = 0; $tf.WordWrap = -1; $tf.MarginLeft = 0; $tf.MarginTop = 0; $tf.MarginBottom = 0
        $tr = $tf.TextRange; $tr.Font.Name = $FONT; $tr.Font.NameFarEast = $FONT
        if ($sh.Top -lt 40) {   # 区分ラベル
            $sh.Left = 32; $sh.Top = 8; $sh.Width = 880; $sh.Height = 18
            $tr.Font.Size = 12; $tr.Font.Bold = -1; $tr.Font.Color.RGB = (RGBc "D5EEEA")
        } else {                # 見出し
            $sh.Left = 32; $sh.Top = 27; $sh.Width = 900; $sh.Height = 34
            $tr.Font.Size = 22; $tr.Font.Bold = -1; $tr.Font.Color.RGB = (RGBc "FFFFFF"); $tf.VerticalAnchor = 1
        }
    }
    # ページ番号（テンプレートと同じ書き方）
    $pg = $s.Shapes.AddTextbox(1, 773, 508, 158, 19); $pg.Name = "TplPage"
    $pg.TextFrame.TextRange.Text = ("社内向け　{0:D2} / {1:D2}" -f $s.SlideIndex, $n)
    $pg.TextFrame.TextRange.Font.Name = $FONT; $pg.TextFrame.TextRange.Font.NameFarEast = $FONT
    $pg.TextFrame.TextRange.Font.Size = 9.75; $pg.TextFrame.TextRange.Font.Color.RGB = (RGBc "6B7280"); $pg.TextFrame.TextRange.ParagraphFormat.Alignment = 3
}
$p.Save()
"done: " + $p.FullName
