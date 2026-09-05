"""受音点をまたいだ結果のまとめ表（G-21）。

受音点ごとの CSV は受音点ごとのフォルダにあるが、**全測定点を 1 つのファイルで
見たい**という要望（2026-08-21）。設計の検討では点ごとの値を並べて眺めるので、
受音点ごとに開き直すのが手間になる。

作るもの（プロジェクト直下の `結果/`。名前の頭に対象室＋条件名が付く）:

    <対象室_条件>_まとめ_残響時間.csv   EDT / T20 / T30 を受音点ごと ＋ 平均 ＋ **理論値**
                                       （Sabine / Eyring / Eyring-Knudsen）
    <対象室_条件>_まとめ_明瞭度.csv     C50 / C80 / D50 / Ts を受音点ごと ＋ 平均
    <対象室_条件>_まとめ_音圧レベル.csv  帯域別の Lp（合計・A 特性・直接音・反射音）
    <対象室_条件>_まとめ_STI.csv         STI と帯域別 MTI
    <対象室>_まとめ_条件比較.csv         **条件を横に並べた比較**（一括計算のとき）

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
LEVEL_FILE = "まとめ_音圧レベル.csv"
STI_FILE = "まとめ_STI.csv"
# 条件（材料条件表）を横に並べた比較表。**条件名を頭に付けない**（条件をまたぐので）
CONDITION_FILE = "まとめ_条件比較.csv"

# 音圧レベルのまとめに載せる行（`spl.csv` の項目名）。
# 自由音場（逆二乗）との比較は**入れない**（2026-08-21 ユーザー判断）
LEVEL_ROWS = ["Lp_dB", "Lp_A_dB", "直接音_dB", "反射音_dB"]

# レベルの平均は**エネルギー平均**で取る（dB をそのまま平均してはいけない）
LEVEL_ENERGY_AVERAGE = True

# `rt.csv` / `clarity.csv` から拾う行と、まとめ表での並び順
REVERBERATION_ROWS = ["EDT_s", "T20_s", "T30_s"]
STATISTICAL_ROWS = ["sabine_s", "eyring_s", "eyring_knudsen_s"]
CLARITY_ROWS = ["C50_db", "C80_db", "D50", "Ts_s"]


def _find(project, folder, filename):
    """受音点フォルダの中からその CSV を探す。

    ファイル名の頭には対象室＋条件名が付くが、**頭を付ける前に計算した
    プロジェクトも読める**ように、付いていない名前も探す（2026-08-21）。
    """
    for name in project.name_candidates(filename):
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


def _write(filename, frequencies, records, extra_label=None,
           first_label="受音点"):
    """1 列目「受音点」・2 列目「項目」・3 列目以降が周波数の CSV を書く。

    `table.write_frequency_table` は 1 列目が 1 つだけなので、受音点と項目の
    2 段になるここでは使わずに書く（周波数を横に並べる約束は同じ）。
    """
    os.makedirs(os.path.dirname(os.path.abspath(filename)) or ".", exist_ok=True)

    def cell(value):
        if value is None:
            return ""
        if isinstance(value, str):
            return value
        return "" if np.isnan(value) else "%.12g" % value

    with open(filename, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.writer(f)
        head = [first_label, "項目"] + ([extra_label] if extra_label else [])
        writer.writerow(head + [f"{v:.0f}" for v in frequencies])
        for record in records:
            if extra_label:
                who, what, extra, values = record
                lead = [who, what, cell(extra)]
            else:
                who, what, values = record
                lead = [who, what]
            band = ([""] * len(frequencies) if values is None
                    else [cell(v) for v in values])
            writer.writerow(lead + band)
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


def _read_sectioned(project, folder, filename):
    """区分付きの表（`spl.csv` / `sti.csv`）を読む。無ければ None。"""
    import table as tb
    return tb.read_sectioned_table(_find(project, folder, filename))


def _level_average(values):
    """レベルの平均。**dB をそのまま平均せずエネルギーで平均する。**

    3 dB 違う 2 点の平均は 1.5 dB ではなく約 1.8 dB。
    測定点をまたいだ代表値を出すときはこちらが正しい。
    """
    block = np.array(values, dtype=float)
    with np.errstate(invalid="ignore"):
        return 10.0 * np.log10(np.nanmean(10.0 ** (block / 10.0), axis=0))


def write_level_summary(project, verbose=True):
    """音圧レベルのまとめ表。受音点ごと ＋ 平均（エネルギー平均）。

    3 列目に**音源距離**を入れてある。距離と Lp を並べれば
    そのまま逆二乗（6 dB/倍距離）の確認に使える（依頼 2026-08-21）。
    """
    folders = receiver_folders(project)
    frequencies, records, gathered = None, [], {k: [] for k in LEVEL_ROWS}

    for name, folder in folders:
        table = _read_sectioned(project, folder, "spl.csv")
        if table is None:
            continue
        frequencies = frequencies if frequencies is not None else table["frequencies"]
        # ★音源距離は独立した行（2026-09-05）。**昔のファイルも読める**ように、
        #   以前の置き場（音源パワーレベルの行の 3 列目）も見る
        distance = (table["values"].get("音源距離_m")
                    or table["values"].get("音源パワーレベル_dB", ""))
        for key in LEVEL_ROWS:
            if key in table["rows"]:
                records.append((name, key, distance if key == "Lp_dB" else None,
                                table["rows"][key]))
                gathered[key].append(table["rows"][key])
    if frequencies is None:
        return None

    for key in LEVEL_ROWS:
        if gathered[key]:
            average = (_level_average(gathered[key]) if LEVEL_ENERGY_AVERAGE
                       else np.nanmean(np.array(gathered[key]), axis=0))
            records.append(("平均", key, None, average))

    path = os.path.join(project.folder, pj.RESULT_DIR,
                        project.prefixed(LEVEL_FILE))
    _write(path, frequencies, records, extra_label="音源距離_m")
    if verbose:
        print(f"[まとめ] 音圧レベル（{len(folders)} 点 ＋ 平均）: {path}")
    return path


def write_sti_summary(project, verbose=True):
    """STI のまとめ表。受音点ごと（総合値 ＋ 帯域別 MTI）＋ 平均。"""
    folders = receiver_folders(project)
    frequencies, records, values, band_values = None, [], [], []

    for name, folder in folders:
        table = _read_sectioned(project, folder, "sti.csv")
        if table is None:
            continue
        frequencies = frequencies if frequencies is not None else table["frequencies"]
        sti = table["values"].get("STI", "")
        records.append((name, "STI", sti, table["rows"].get("MTI")))
        records.append((name, "評価", table["values"].get("評価", ""), None))
        try:
            values.append(float(sti))
        except ValueError:
            pass
        if "MTI" in table["rows"]:
            band_values.append(table["rows"]["MTI"])
    if frequencies is None:
        return None

    if values:
        import sound_level as sl
        mean = float(np.mean(values))
        records.append(("平均", "STI", "%.3f" % mean,
                        np.nanmean(np.array(band_values), axis=0)
                        if band_values else None))
        records.append(("平均", "評価", sl.sti_rating(mean), None))
        if len(values) > 1:
            records.append(("ばらつき", "STI",
                            "%.3f" % (max(values) - min(values)), None))

    path = os.path.join(project.folder, pj.RESULT_DIR, project.prefixed(STI_FILE))
    _write(path, frequencies, records, extra_label="総合")
    if verbose:
        print(f"[まとめ] STI（{len(folders)} 点 ＋ 平均）: {path}")
    return path


def write_condition_summary(project, conditions=None, verbose=True):
    """**条件（材料条件表）を横に並べた比較表**（依頼 2026-08-21 の一括計算用）。

    1 行が「条件 × 指標」、列が周波数。どの条件がいちばん良いかを 1 枚で見る。
    載せるのは全受音点の**平均**（点ごとの値は条件別のまとめ表にある）。

        条件,項目,総合,63,125,…
        現状,T30_s,,1.282,1.027,…
        現状,STI,0.852,0.651,…
        吸音追加案,T30_s,,0.951,0.804,…

    `conditions` は条件表のパスのリスト。None ならフォルダの条件表を全部。
    結果が無い条件は飛ばす。
    """
    import condition_table as ct

    if conditions is None:
        conditions = ct.discover(project.folder)
    records, frequencies = [], None

    for item in conditions:
        file_name, sheet = item if isinstance(item, (tuple, list)) else (item, None)
        sub = pj.Project(project.folder,
                         **{k: getattr(project, k) for k in pj.DEFAULTS})
        sub.condition_csv = file_name
        sub.condition_sheet = sheet or ""
        label = sub.condition_label or "（既定）"

        for filename, keys, skip in (
                (REVERBERATION_FILE, ("EDT_s", "T30_s"), 2),
                (CLARITY_FILE, ("C50_db", "D50"), 2),
                (LEVEL_FILE, ("Lp_dB",), 3)):
            rows = _read_summary(sub, filename, skip)
            if rows is None:
                continue
            frequencies = frequencies if frequencies is not None else rows[0]
            for key in keys:
                if ("平均", key) in rows[1]:
                    records.append((label, key, None, rows[1][("平均", key)]))

        sti = _read_summary(sub, STI_FILE, 3)
        if sti is not None and ("平均", "STI") in sti[1]:
            records.append((label, "STI", sti[2].get(("平均", "STI"), ""),
                            sti[1][("平均", "STI")]))
            records.append((label, "評価", sti[2].get(("平均", "評価"), ""), None))

    if frequencies is None or not records:
        if verbose:
            print("[まとめ] 条件の比較表は作れません（条件ごとの結果がまだありません）")
        return None

    room = project.room_label
    name = f"{room}_{CONDITION_FILE}" if room else CONDITION_FILE
    path = os.path.join(project.folder, pj.RESULT_DIR, name)
    _write(path, frequencies, records, extra_label="総合", first_label="条件")
    if verbose:
        conditions_found = len({r[0] for r in records})
        print(f"[まとめ] 条件の比較（{conditions_found} 条件）: {path}")
    return path


def _read_summary(project, filename, skip):
    """まとめ表を読む。戻り値 (周波数, {(受音点, 項目): 値}, {(受音点, 項目): 3 列目})。

    まとめ表は 1 列目が受音点、2 列目が項目で、`skip` が 3 なら 3 列目に
    周波数に依らない値が入る（`_write` と対応）。読めなければ None。
    """
    path = _find(project, os.path.join(project.folder, pj.RESULT_DIR), filename)
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8-sig", newline="") as f:
        table = [row for row in csv.reader(f) if row and any(c.strip() for c in row)]
    if len(table) < 2:
        return None

    def number(text):
        try:
            return float((text or "").strip())
        except ValueError:
            return np.nan

    frequencies = np.array([number(c) for c in table[0][skip:]])
    values, extras = {}, {}
    for row in table[1:]:
        key = (row[0].strip(), row[1].strip())
        values[key] = np.array([number(c) for c in row[skip:skip + len(frequencies)]])
        extras[key] = row[2].strip() if skip >= 3 and len(row) > 2 else ""
    return frequencies, values, extras


def write_all(project, verbose=True):
    """まとめ表をすべて書く。書けたファイルのパスのリストを返す。"""
    written = [write_reverberation_summary(project, verbose=verbose),
               write_clarity_summary(project, verbose=verbose),
               write_level_summary(project, verbose=verbose),
               write_sti_summary(project, verbose=verbose)]
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
