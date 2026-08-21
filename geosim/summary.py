"""受音点をまたいだ結果のまとめ表（G-21）。

受音点ごとの CSV は受音点ごとのフォルダにあるが、**全測定点を 1 つのファイルで
見たい**という要望（2026-08-21）。設計の検討では点ごとの値を並べて眺めるので、
受音点ごとに開き直すのが手間になる。

作るもの（プロジェクト直下の `結果/`。名前の頭に対象室＋条件名が付く）:

    <対象室_条件>_まとめ_残響時間.csv   EDT / T20 / T30 を受音点ごと ＋ 平均 ＋ **理論値**
                                       （Sabine / Eyring / Eyring-Knudsen）
    <対象室_条件>_まとめ_明瞭度.csv     C50 / C80 / D50 / Ts を受音点ごと ＋ 平均

**周波数は横**（`table.py` の共通ルール）。1 列目が「受音点」、2 列目が「項目」で、
3 列目以降が周波数。受音点ごとの CSV（`rt.csv` など）はそのまま残す。

理論値を残響時間のほうに混ぜているのは、
「音線法の値が統計残響式とどれだけ離れているか」をその場で見たいため
（実案件では 2 倍前後離れる。減衰が 2 段になるモデルでは離れるのが正しい）。
"""

import csv
import os

import numpy as np

import project as pj

# まとめ表のファイル名
REVERBERATION_FILE = "まとめ_残響時間.csv"
CLARITY_FILE = "まとめ_明瞭度.csv"

# `rt.csv` / `clarity.csv` から拾う行と、まとめ表での並び順
REVERBERATION_ROWS = ["EDT_s", "T20_s", "T30_s"]
STATISTICAL_ROWS = ["sabine_s", "eyring_s", "eyring_knudsen_s"]
CLARITY_ROWS = ["C50_db", "C80_db", "D50", "Ts_s"]


def _find(project, folder, filename):
    """受音点フォルダの中からその CSV を探す。

    ファイル名の頭には対象室＋条件名が付くが、**頭を付ける前に計算した
    プロジェクトも読める**ように、付いていない名前も探す（2026-08-21）。
    """
    for name in (project.prefixed(filename), filename):
        path = os.path.join(folder, name)
        if os.path.exists(path):
            return path
    return os.path.join(folder, project.prefixed(filename))


def receiver_folders(project):
    """受音点ごとの結果フォルダを順に返す [(表示名, フォルダ), …]。

    新しい置き方（2026-08-21 以降）は `結果/rec1/` `結果/rec2/` …。
    古い置き方（1 点目が `結果/` 直下、2 点目以降が `rec2/結果/`）も読める
    ようにしてある（作り直す前の結果をそのまま見たいことがあるため）。
    """
    folders = []
    index = 1
    while True:
        folder = os.path.join(project.folder, pj.RESULT_DIR,
                              pj.RECEIVER_DIR % index)
        if not os.path.isdir(folder):
            break
        folders.append((pj.RECEIVER_DIR % index, folder))
        index += 1
    if folders:
        return folders

    # ---- 古い置き方 ----
    folders = [("rec1", project.path(pj.RESULT_DIR))]
    index = 2
    while True:
        folder = os.path.join(project.folder, f"rec{index}", pj.RESULT_DIR)
        if not os.path.isdir(folder):
            break
        folders.append((f"rec{index}", folder))
        index += 1
    return folders


def _read(path):
    """周波数を横に並べた CSV を (周波数, {行名: 値}) で読む。無ければ (None, {})。"""
    if not os.path.exists(path):
        return None, {}
    with open(path, encoding="utf-8-sig", newline="") as f:
        rows = list(csv.reader(f))
    if len(rows) < 2:
        return None, {}
    frequencies = [float(v) for v in rows[0][1:]]
    values = {r[0]: [float(v) if v else np.nan for v in r[1:]] for r in rows[1:]}
    return frequencies, values


def _write(filename, frequencies, records):
    """1 列目「受音点」・2 列目「項目」・3 列目以降が周波数の CSV を書く。

    `table.write_frequency_table` は 1 列目が 1 つだけなので、受音点と項目の
    2 段になるここでは使わずに書く（周波数を横に並べる約束は同じ）。
    """
    os.makedirs(os.path.dirname(os.path.abspath(filename)) or ".", exist_ok=True)
    with open(filename, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["受音点", "項目"] + [f"{v:.0f}" for v in frequencies])
        for who, what, values in records:
            writer.writerow([who, what] + ["" if v is None or np.isnan(v)
                                           else "%.12g" % v for v in values])
    return filename


def _statistical_rows(project):
    """残響時間の理論値 {項目: (nf,)} を `結果/` 直下から読む。

    いまは『吸音率と理論値.csv』（材料別の吸音率 → 平均吸音率 → 理論値）に
    入っている。1 枚にまとめる前の `rt_statistical.csv` も読めるようにしてある。
    """
    # 「room」は受音点に依らない（`pj.SHARED_RESULTS`）ので `結果/` 直下を見る
    room = pj.read_room_csv(project.existing_result_path("room"))
    if room is not None:
        return room["rows"]
    _, values = _read(os.path.join(project.folder, pj.RESULT_DIR,
                                   "rt_statistical.csv"))
    return values


def write_reverberation_summary(project, verbose=True):
    """残響時間のまとめ表。受音点ごと ＋ 平均 ＋ 理論値。"""
    folders = receiver_folders(project)
    frequencies, records, gathered = None, [], {k: [] for k in REVERBERATION_ROWS}

    for name, folder in folders:
        freq, values = _read(_find(project, folder, "rt.csv"))
        if freq is None:
            continue
        frequencies = frequencies or freq
        for key in REVERBERATION_ROWS:
            if key in values:
                records.append((name, key, values[key]))
                gathered[key].append(values[key])
    if frequencies is None:
        return None

    # 全測定点の平均。JIS / ISO の運用でも点の平均を代表値にするので並べておく
    for key in REVERBERATION_ROWS:
        if gathered[key]:
            records.append(("平均", key, np.nanmean(np.array(gathered[key]), axis=0)))
    # ばらつき（最大 − 最小）。点による差が見えると設計の判断に使える
    for key in REVERBERATION_ROWS:
        if len(gathered[key]) > 1:
            block = np.array(gathered[key])
            records.append(("ばらつき", key,
                            np.nanmax(block, axis=0) - np.nanmin(block, axis=0)))

    # 理論値（統計残響式）。受音点に依らないので `結果/` 直下のものを 1 回だけ。
    # 『吸音率と理論値.csv』にまとめる前の `rt_statistical.csv` も読める
    statistical = _statistical_rows(project)
    for key in STATISTICAL_ROWS:
        if key in statistical:
            records.append(("理論値", key, statistical[key]))

    path = os.path.join(project.folder, pj.RESULT_DIR,
                        project.prefixed(REVERBERATION_FILE))
    _write(path, frequencies, records)
    if verbose:
        print(f"[まとめ] 残響時間（{len(folders)} 点 ＋ 平均 ＋ 理論値）: {path}")
    return path


def write_clarity_summary(project, verbose=True):
    """明瞭度のまとめ表。受音点ごと ＋ 平均。"""
    folders = receiver_folders(project)
    frequencies, records, gathered = None, [], {k: [] for k in CLARITY_ROWS}

    for name, folder in folders:
        freq, values = _read(_find(project, folder, "clarity.csv"))
        if freq is None:
            continue
        frequencies = frequencies or freq
        for key in CLARITY_ROWS:
            if key in values:
                records.append((name, key, values[key]))
                gathered[key].append(values[key])
    if frequencies is None:
        return None

    for key in CLARITY_ROWS:
        if gathered[key]:
            records.append(("平均", key, np.nanmean(np.array(gathered[key]), axis=0)))
    for key in CLARITY_ROWS:
        if len(gathered[key]) > 1:
            block = np.array(gathered[key])
            records.append(("ばらつき", key,
                            np.nanmax(block, axis=0) - np.nanmin(block, axis=0)))

    path = os.path.join(project.folder, pj.RESULT_DIR,
                        project.prefixed(CLARITY_FILE))
    _write(path, frequencies, records)
    if verbose:
        print(f"[まとめ] 明瞭度（{len(folders)} 点 ＋ 平均）: {path}")
    return path


def write_all(project, verbose=True):
    """まとめ表をすべて書く。書けたファイルのパスのリストを返す。"""
    written = [write_reverberation_summary(project, verbose=verbose),
               write_clarity_summary(project, verbose=verbose)]
    return [p for p in written if p]


def main():
    import argparse

    p = argparse.ArgumentParser(
        description="受音点をまたいだ結果のまとめ表を作る（計算はやり直さない）")
    p.add_argument("folder", help="プロジェクトフォルダ")
    a = p.parse_args()

    project = pj.Project.load(a.folder)
    written = write_all(project)
    if not written:
        raise SystemExit(f"{a.folder} に結果が見つかりません（先に計算してください）")


if __name__ == "__main__":
    main()
