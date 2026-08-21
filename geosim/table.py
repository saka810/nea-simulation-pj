"""表（CSV）の並べ方の共通ルール。

# ★共通ルール：周波数は「横」に並べる（2026-08-17 ユーザー判断）

**グラフにしたときの横軸が周波数だから、表も横方向を周波数に揃える。**
縦に並んでいると、Excel でグラフを作るたびに行と列を選び直す手間が要る。

    ○ 正しい向き
        項目,125,250,500,1000,2000,4000
        EDT_s,0.945,0.571,0.459,0.295,0.267,0.284
        T20_s,0.633,0.569,0.404,0.255,0.212,0.207

    × 直す前
        frequency_hz,EDT_s,T20_s
        125,0.945,0.633
        250,0.571,0.569

**これから新しい出力を書くときも、この向きにすること。**
周波数が「行」になるのは、1 行が周波数以外の何か（時刻・面・経路）で、
周波数がそのバリエーションになっている場合だけ。その場合も**列**を周波数にする。

| ファイル | 1 行の意味 | 周波数の位置 |
|---|---|---|
| `rt.csv` / `clarity.csv` | 指標 | **列**（このモジュールの形） |
| `decay.csv` | 時刻 | 列（`decay_125Hz_db` …） |
| `吸音率と理論値.csv` | 材料／指標（1 列目が区分） | 列（`project.write_room_csv`） |
| `pulses.csv` | 経路 | 列（`energy_125Hz` …） |
| `ir.csv` | 時刻 | 周波数の軸を持たない |

**どうしても縦にしたい事情が出たら、勝手に決めずに確認すること。**
"""

import csv
import os

import numpy as np

# 1 列目の見出し。ここが `frequency_hz` なら**古い縦向き**と判断する
LABEL_HEADER = "項目"
OLD_LABEL_HEADER = "frequency_hz"


def band_column(prefix, frequency, unit=""):
    """`alpha_125Hz` のような列名を作る。

    1 行が周波数以外（材料・経路・時刻）のとき、周波数は列で表す。
    番号（`energy_1`）ではなく**周波数そのもの**を入れるのは、
    バンド数が 6 と 8 で意味が変わってしまうため。
    """
    return f"{prefix}_{float(frequency):.0f}Hz" + (f"_{unit}" if unit else "")


def write_frequency_table(filename, frequencies, rows, label=LABEL_HEADER,
                          fmt="%.12g"):
    """周波数を**横**に並べた CSV を書く。

    引数:
        frequencies (nf,)   列になる周波数 [Hz]
        rows                {行の名前: (nf,) の値} または [(名前, 値), ...]
        label               1 列目の見出し

    Excel でそのまま開けるよう **BOM 付き UTF-8** で書く。
    """
    os.makedirs(os.path.dirname(os.path.abspath(filename)) or ".", exist_ok=True)
    frequencies = np.asarray(frequencies, dtype=float).ravel()
    items = list(rows.items()) if isinstance(rows, dict) else list(rows)

    with open(filename, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([label] + [f"{v:.0f}" for v in frequencies])
        for name, values in items:
            values = np.asarray(values, dtype=float).ravel()
            if len(values) != len(frequencies):
                raise ValueError(
                    f"{name!r} の要素数 {len(values)} が周波数 {len(frequencies)} と違います")
            writer.writerow([name] + ["" if np.isnan(v) else fmt % v
                                      for v in values])
    return filename


def read_frequency_table(path):
    """`write_frequency_table` が書いた CSV を読む。

    **古い縦向き（1 列目が `frequency_hz`）もそのまま読める。**
    向きを変えた前に作ったプロジェクトを開いても壊れないようにするため。

    戻り値 (周波数 (nf,), {行の名前: (nf,)})。読めなければ (None, {})。
    """
    if not os.path.exists(path):
        return None, {}
    with open(path, encoding="utf-8-sig", newline="") as f:
        table = [row for row in csv.reader(f) if row and any(c.strip() for c in row)]
    if len(table) < 2:
        return None, {}

    header, body = table[0], table[1:]

    def number(text):
        text = (text or "").strip()
        if not text:
            return np.nan
        try:
            return float(text)
        except ValueError:
            return np.nan

    if header[0].strip().lower() == OLD_LABEL_HEADER:
        # 古い縦向き：1 列目が周波数、2 列目以降が指標
        frequencies = np.array([number(row[0]) for row in body])
        rows = {name: np.array([number(row[i + 1]) if i + 1 < len(row) else np.nan
                                for row in body])
                for i, name in enumerate(header[1:])}
        return frequencies, rows

    frequencies = np.array([number(c) for c in header[1:]])
    rows = {row[0].strip(): np.array([number(c) for c in row[1:1 + len(frequencies)]])
            for row in body}
    return frequencies, rows
