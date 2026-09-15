"""**音源が複数あるときの結果の見方**（合成）。

2026-09-15 ユーザー要望（不具合報告 ⑨ の直し方の相談）:

> 複数音源の場合は、音源 1 つずつ計算して。その後、S1 の場合の結果を見るのか、
> S2 の結果を見るのか、S1 と S2 の結果の平均を見るのか、S1 と S2 の
> インパルス応答を重ね合わせてそこから結果を出すのか、重ね合わせる場合
> 時間遅れを考慮するのか、結果の見方は選択できるようにして欲しい。

**計算そのものは音源 1 点ずつ**（`run_project.run`）。ここはその結果を
組み合わせるだけで、★**音線追跡もバックトレースもやり直さない**
（`view_images.py` や `summary.py` と同じ流儀）。

置き場（`結果/` の下の棚）:

    結果/src1/rec1/   … 音源 1 だけの結果（S1 を見る）
    結果/src2/rec1/   … 音源 2 だけの結果（S2 を見る）
    結果/合成_平均/rec1/                   … 指標の平均（ISO 3382 の作法）
    結果/合成_重ね合わせ/rec1/             … パルス列を**時間差ありのまま**重ねる
    結果/合成_重ね合わせ_時間差なし/rec1/  … 音源ごとに直接音を 0 s に揃えてから重ねる

3 つの違い（★ここが判断の分かれ目）:

  ・**平均** … 音源位置ごとに測って平均する ISO 3382 の作法。
    「1 つずつ鳴らしたときの代表値」。干渉は起きない（別々の測定の平均なので）
  ・**重ね合わせ（時間差あり）** … **2 台を同時に鳴らした**ときの音。
    音源から受音点までの伝搬遅れがそのまま入るので、**2 つの音が干渉する**。
    PA のように実際に同時に鳴らす使い方はこれ
  ・**重ね合わせ（時間差なし）** … 音源ごとに**自分の直接音を 0 s に揃えて**から
    重ねる。伝搬距離の違いによる到達時刻のずれを消した形。
    ★**大きさ（距離減衰・空気吸収）はそのまま**で、時刻だけ揃える

★★**重ね合わせは「エネルギーを足す」のではなく、パルス列を 1 本にまとめてから
  いつもの合成に通す**。インパルス応答は位相ごと重ねる複素和（式(2)）なので、
  パルス列を並べた時点で干渉は正しく入る（`CLAUDE.md`「プログラムの足し方」）。

★**音源のパワーレベル（PWL）は音源ごとに変えられない**。パルスのエネルギーは
  「音源パワー 1 あたり」の値で、`Lw` は最後に 1 回足すためである
  （重ね合わせは**同じ PWL の音源が複数**という扱いになる）。
  音源ごとに出力を変えたいときは別途重みが要る（TODO F-7b）。

使い方（計算のあとで見方を足したくなったとき。計算はやり直さない）:

    cd geosim
    python source_mix.py "C:\\...\\プロジェクト"
"""

import os

import numpy as np

import absorption as ab
import loop_noredundancy as ln
import project as pj
import reverberation as rv
from atmosphere import Atmosphere

# 見方（`project.json` の `source_combination`）
MIX_NONE = "個別のみ"
MIX_AVERAGE = "平均"
MIX_SUM = "重ね合わせ"
MIX_SUM_ALIGNED = "重ね合わせ_時間差なし"
MIX_ALL = "すべて"

# 設定画面に出す選択肢（表示文, 値）。**個別（音源ごと）は必ず作る**ので、
# ここで選ぶのは「そのうえでどの合成を作るか」
MIX_CHOICES = [
    ("すべて（個別＋平均＋重ね合わせ 2 通り）", MIX_ALL),
    ("個別のみ（音源ごとの結果だけ）", MIX_NONE),
    ("個別＋平均（音源位置ごとの平均。ISO 3382 の作法）", MIX_AVERAGE),
    ("個別＋重ね合わせ（同時に鳴らす。時間差あり）", MIX_SUM),
    ("個別＋重ね合わせ（直接音を揃える。時間差なし）", MIX_SUM_ALIGNED),
]

# 棚の名前（`結果/` の下のフォルダ名）
FOLDERS = {
    MIX_AVERAGE: pj.MIX_PREFIX + "平均",
    MIX_SUM: pj.MIX_PREFIX + "重ね合わせ",
    MIX_SUM_ALIGNED: pj.MIX_PREFIX + "重ね合わせ_時間差なし",
}

# 平均のとり方。**音圧レベルだけエネルギー平均**（dB をそのまま平均してはいけない）。
# 残響時間・明瞭度・STI は算術平均（ISO 3382-1 / IEC 60268-16 の扱い）
ENERGY_AVERAGE_FILES = {"spl"}
AVERAGE_KEYS = ("rt", "clarity", "spl", "sti")


def modes_for(value):
    """設定の値 → 作る合成のリスト。知らない値は「すべて」とみなす。"""
    value = (str(value) if value is not None else MIX_ALL).strip()
    if value == MIX_NONE:
        return []
    if value == MIX_ALL:
        return [MIX_AVERAGE, MIX_SUM, MIX_SUM_ALIGNED]
    if value in FOLDERS:
        return [value]
    return [MIX_AVERAGE, MIX_SUM, MIX_SUM_ALIGNED]


# ------------------------------------------------------------------------------
# 棚を指す Project
# ------------------------------------------------------------------------------

def tagged(project, tag=None, index=None, receiver_index=None):
    """その棚・その受音点を指す `Project` を作る（フォルダは同じ）。

    `tag` は合成の棚（`合成_平均` など）、`index` は音源の番号（1 始まり）。
    ★**中身の条件は親と同じ**にする（条件名が結果ファイル名に入るため）。
    """
    sub = pj.Project(project.folder,
                     **{k: getattr(project, k) for k in pj.DEFAULTS})
    sub.source_tag = tag
    sub.source_index = index
    sub.receiver_index = receiver_index
    return sub


def source_tags(project):
    """`結果/` にある**音源ごと**の棚（`src1` `src2` …）。合成の棚は含まない。"""
    return [name for name in project.source_folders()
            if not name.startswith(pj.MIX_PREFIX)]


def default_shelf(project, verbose=True):
    """**結果を見るときに開く棚**を返す（音源が 1 点なら `project` のまま）。

    音源が複数あると結果は `結果/src1/` … に分かれるので、`結果/` 直下を見る
    画面（3D の可視化・保存済み結果の読み出し）はそのままでは何も見つけられない。
    ★ここでは **1 番目の音源**を既定にして、そう言う（黙って空にしない）。
    他の棚を見たいときは `結果/srcM/` を直接開く。
    """
    if project.source_folder:
        return project
    tags = source_tags(project)
    if len(tags) < 2:
        return project
    if verbose:
        print(f"[合成] 音源が {len(tags)} 点あります。"
              f"ここでは『{tags[0]}』を開きます"
              f"（他は 結果/{tags[1]}/ などを直接見てください）")
    return tagged(project, tag=tags[0])


def receiver_count(project, tag):
    """その棚にある受音点フォルダの数。"""
    root = os.path.join(project.path(pj.RESULT_DIR), tag)
    index = 0
    while os.path.isdir(os.path.join(root, pj.RECEIVER_DIR % (index + 1))):
        index += 1
    return index


# ------------------------------------------------------------------------------
# ① パルス列を重ねる
# ------------------------------------------------------------------------------

def merge_pulses(lists, aligned=False):
    """パルス列を 1 本にまとめる。→ `PulseList`

    `aligned=True` なら**音源ごとに自分の最初の到来を 0 s に揃えて**から重ねる
    （＝伝搬遅れの違いを消す）。★**時刻だけ動かし、距離とエネルギーは触らない**
    （距離減衰と空気吸収は距離で決まるので、揃えても大きさは変わらない）。

    ★最初の到来 = 最短経路 = 直接音（遮られていればいちばん早い反射）。
    `argmin` で拾うので、直接音が無い配置でも破綻しない。
    """
    lists = [p for p in lists if p is not None and len(p)]
    if not lists:
        return None
    bands = {p.energy.shape[1] for p in lists}
    if len(bands) != 1:
        raise ValueError(f"音源によってバンド数が違います: {sorted(bands)}")

    times = []
    for p in lists:
        shift = float(p.time.min()) if aligned else 0.0
        times.append(p.time - shift)

    out = ln.PulseList(lists[0].energy.shape[1], lists[0].sound_velocity)
    out.reflection_count = np.concatenate([p.reflection_count for p in lists])
    out.time = np.concatenate(times)
    out.distance = np.concatenate([p.distance for p in lists])
    out.direction = np.vstack([p.direction for p in lists])
    out.energy = np.vstack([p.energy for p in lists])
    return out.sort_by_time()


def read_pulses(project, tags, receiver_index):
    """音源ごとのパルス列を読む。1 つでも無ければ None。"""
    lists = []
    for order, tag in enumerate(tags):
        sub = tagged(project, tag=tag, receiver_index=receiver_index)
        path = sub.existing_result_path("pulses")
        if not os.path.exists(path):
            return None
        lists.append(ln.PulseList.from_csv(path))
    return lists


def aligned_impulse(project, lists, frequencies, air, band_width):
    """音源ごとにインパルス応答を作り、**それぞれの直接音を 0 s に寄せて**足す。

    ★★**時刻をずらしたパルス列を 1 本にまとめて合成してはいけない。**
      `impulse.impulse_response` は**到来時刻から伝搬距離を出す**
      （`振幅 = √E /(t·c)`）ので、時刻をずらすと**大きさまで変わる**。
      0 s に寄せた直接音は距離 0 になって発散し、応答が丸ごと NaN になる
      （2026-09-15 に実際に踏んだ。残響時間の表が全部空になった）。
      音源ごとに**正しい時刻のまま**合成してから、**波形をずらして足す**。

    ずらす量は「その音源の最初の到来」＝直接音の時刻。サンプルに丸めるので
    端数（最大 1/2 サンプル ≒ 11 μs）のずれは残る。★**時間差を無視する**
    という指定なので、ここは端数まで合わせる意味がない。
    """
    import impulse as ir

    fs = ir.SAMPLING_FREQUENCY
    method = getattr(project, "impulse_method", "fast")
    time_axis, total = None, None
    for one in lists:
        t, y = ir.impulse_response(one.time, one.energy,
                                   octave_frequencies=frequencies,
                                   atmosphere=air, band_width=band_width,
                                   max_time=project.max_time, method=method,
                                   verbose=False)
        shift = int(round(float(one.time.min()) * fs))
        if shift > 0:
            y = np.roll(y, -shift)
            y[-shift:] = 0.0
        time_axis = t
        total = y if total is None else total + y
    return time_axis, total


def _atmosphere(project):
    return Atmosphere(temperature=project.temperature,
                      humidity=project.humidity,
                      pressure=project.pressure)


def write_superposed(project, tags, mode, verbose=True):
    """パルス列を重ねて、その 1 本から指標を出し直す。→ 書いたフォルダのリスト

    やることは `procedure.process()` の後半と同じ（インパルス応答 → 残響指標 →
    明瞭度、パルス列 → 音圧レベル・STI）。**前半（音線追跡・バックトレース）は
    やり直さない**。
    """
    import impulse as ir
    import sound_level as sl

    aligned = (mode == MIX_SUM_ALIGNED)
    folder_tag = FOLDERS[mode]
    frequencies = ab.frequency_bands(project.band_number,
                                     getattr(project, "band_width", "1/1"),
                                     getattr(project, "band_start", None))
    air = _atmosphere(project)
    band_width = getattr(project, "band_width", "1/1")
    written = []

    for index in range(1, receiver_count(project, tags[0]) + 1):
        lists = read_pulses(project, tags, index)
        if not lists:
            if verbose:
                print(f"[合成] 受音点 {index}: パルス列がそろっていないので飛ばします")
            continue
        pulses = merge_pulses(lists, aligned=aligned)
        if pulses is None or not len(pulses):
            continue

        sub = tagged(project, tag=folder_tag, receiver_index=index)
        sub.ensure_dirs()
        # ★昔の名前が残ると結果フォルダに並んで、どちらが新しいか分からなくなる。
        #   音源に依らないもの（`結果/` 直下）は消さない（`clear_results` が守る）
        sub.clear_results(verbose=False)
        pulses.save_csv(sub.result_path("pulses"))

        if aligned:
            # ★音源ごとに合成してから波形をずらして足す（`aligned_impulse` 参照）
            impulse = aligned_impulse(project, lists, frequencies, air, band_width)
            ir.write_impulse_response(sub.result_path("ir"), impulse[0], impulse[1])
        else:
            impulse = ir.impulse_response_from_pulses(
                sub.result_path("ir"), pulses, octave_frequencies=frequencies,
                atmosphere=air, band_width=band_width, max_time=project.max_time,
                method=getattr(project, "impulse_method", "fast"), verbose=False)
        rt = rv.reverberation_time(impulse[0], impulse[1],
                                   rt_filename=sub.result_path("rt"),
                                   decay_filename=sub.result_path("decay"),
                                   frequencies=frequencies,
                                   band_width=band_width, verbose=False)
        clarity = rv.clarity_measures(impulse[0], impulse[1],
                                      frequencies=frequencies,
                                      band_width=band_width, verbose=False)
        rv.write_clarity_measures(sub.clarity_path(), clarity)

        # ★音圧レベルの音源距離は**いちばん近い音源まで**（自由音場の目安に使うだけ）。
        #   None を渡すと最短経路長＝最も近い音源の直接音になる
        level = sl.band_levels(pulses.time, pulses.energy, pulses.distance, air,
                               frequencies,
                               source_power_db=project.source_power_db,
                               verbose=False)
        if str(getattr(project, "level_method", "both")).lower() in ("coherent", "both"):
            coherent = sl.coherent_levels(pulses.time, pulses.energy,
                                          pulses.distance, air, frequencies,
                                          band_width=band_width,
                                          source_power_db=project.source_power_db)
            level["coherent"] = coherent["levels"]
        sl.write_levels(sub.result_path("spl"), level)

        if not ab.is_third_octave(band_width):
            sti = sl.speech_transmission_index(
                pulses.time, pulses.energy, pulses.distance, air, frequencies,
                source_power_db=project.source_power_db,
                noise_level_db=project.noise_level_db, verbose=False)
            sl.write_sti(sub.result_path("sti"), sti)

        written.append(sub.result_dir())
        if verbose:
            print(f"[合成] {folder_tag} / rec{index}: "
                  f"パルス {len(pulses)} 本（{' + '.join(str(len(p)) for p in lists)}）")
    return written


# ------------------------------------------------------------------------------
# ② 指標を平均する
# ------------------------------------------------------------------------------

def _average_rows(tables, energy=False):
    """同じ形の表（行のリスト）を平均する。文字の欄は 1 つ目のものを残す。

    `energy=True` の欄は**エネルギー平均**（`10log10(mean(10^(L/10)))`）。
    ★dB をそのまま平均してはいけない（`summary.py` と同じ約束）。
    """
    base = [list(row) for row in tables[0]]
    for r, row in enumerate(base):
        for c, cell in enumerate(row):
            values = []
            for table in tables:
                if r >= len(table) or c >= len(table[r]):
                    values = []
                    break
                try:
                    values.append(float(str(table[r][c]).strip()))
                except (TypeError, ValueError):
                    values = []
                    break
            if not values:
                continue        # 文字（項目名・評価など）はそのまま
            if energy:
                mean = 10.0 * np.log10(np.mean(np.power(10.0, np.array(values) / 10.0)))
            else:
                mean = float(np.mean(values))
            base[r][c] = "%.12g" % mean
    return base


def _read_csv(path):
    import csv
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return [row for row in csv.reader(handle)]


def _write_csv(path, rows):
    import csv
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as handle:
        csv.writer(handle).writerows(rows)
    return path


def _energy_columns(key, rows):
    """音圧レベルの表で**エネルギー平均にしてよい行**か（dB の行だけ）。

    `spl.csv` は区分付きの表で、`音源距離_m` のように dB でない行も混じる。
    項目名（2 列目）が `_dB` で終わる行だけをエネルギー平均にする。
    """
    if key not in ENERGY_AVERAGE_FILES:
        return [False] * len(rows)
    flags = []
    for row in rows:
        label = row[1] if len(row) > 1 else ""
        flags.append(str(label).endswith("_dB") or str(label).endswith("_db"))
    return flags


def write_average(project, tags, verbose=True):
    """音源ごとの指標を平均して `結果/合成_平均/recN/` に書く。

    ★**表（rt / clarity / spl / sti）だけ**を作る。インパルス応答や減衰曲線は
    「平均」には意味が無い（別々の測定を平均した値なので波形は存在しない）。
    """
    folder_tag = FOLDERS[MIX_AVERAGE]
    written = []
    for index in range(1, receiver_count(project, tags[0]) + 1):
        sub = tagged(project, tag=folder_tag, receiver_index=index)
        sub.ensure_dirs()
        sub.clear_results(verbose=False)
        made = 0
        for key in AVERAGE_KEYS:
            paths = []
            for tag in tags:
                one = tagged(project, tag=tag, receiver_index=index)
                path = one.existing_result_path(key)
                if os.path.exists(path):
                    paths.append(path)
            if len(paths) != len(tags):
                continue
            tables = [_read_csv(p) for p in paths]
            if len({len(t) for t in tables}) != 1:
                if verbose:
                    print(f"[合成] 平均: {key} の表の形が音源ごとに違うので飛ばします")
                continue
            flags = _energy_columns(key, tables[0])
            if any(flags):
                # 行ごとに平均のとり方を変える（dB の行だけエネルギー平均）
                rows = _average_rows(tables, energy=False)
                energy_rows = _average_rows(tables, energy=True)
                for r, use_energy in enumerate(flags):
                    if use_energy:
                        # ★3 列目（総合）も dB なのでまとめてエネルギー平均にする
                        rows[r] = energy_rows[r]
            else:
                rows = _average_rows(tables, energy=False)
            target = (sub.clarity_path() if key == "clarity"
                      else sub.result_path(key))
            _write_csv(target, rows)
            made += 1
        if made:
            written.append(sub.result_dir())
            if verbose:
                print(f"[合成] {folder_tag} / rec{index}: {made} 種の表を平均しました")
    return written


# ------------------------------------------------------------------------------
# まとめ
# ------------------------------------------------------------------------------

def write_all(project, verbose=True):
    """設定（`source_combination`）に従って合成を作る。→ 書いたフォルダのリスト

    音源ごとの棚が 2 つ以上そろっていないときは何もしない。
    """
    tags = source_tags(project)
    if len(tags) < 2:
        return []
    modes = modes_for(getattr(project, "source_combination", MIX_ALL))
    if not modes:
        if verbose:
            print("[合成] 設定が『個別のみ』なので合成は作りません")
        return []

    written = []
    for mode in modes:
        if verbose:
            print(f"\n[合成] ── {FOLDERS[mode]}（音源 {len(tags)} 点）")
        try:
            if mode == MIX_AVERAGE:
                folders = write_average(project, tags, verbose=verbose)
            else:
                folders = write_superposed(project, tags, mode, verbose=verbose)
        except Exception as error:      # 1 つ失敗しても他の見方は残す
            print(f"[合成] {FOLDERS[mode]} を作れませんでした: "
                  f"{type(error).__name__}: {error}")
            continue
        written.extend(folders)
        if folders:
            _summaries(project, FOLDERS[mode], verbose=verbose)
    return written


def _summaries(project, tag, verbose=True):
    """その棚のまとめ表と Excel を作る（`summary` / `workbook` を棚ごとに呼ぶ）。"""
    shelf = tagged(project, tag=tag)
    try:
        import summary as sm
        sm.write_all(shelf, verbose=verbose)
    except Exception as error:
        print(f"[合成] まとめ表を作れませんでした: {type(error).__name__}: {error}")
    try:
        import workbook as wb
        wb.write(shelf, verbose=verbose)
    except Exception as error:
        print(f"[合成] 結果一式（Excel）を作れませんでした: "
              f"{type(error).__name__}: {error}")


def main():
    import argparse

    p = argparse.ArgumentParser(
        description="音源ごとの結果を組み合わせる（計算はやり直さない）")
    p.add_argument("folder", help="プロジェクトフォルダ")
    p.add_argument("--mode", choices=[MIX_ALL, MIX_NONE, MIX_AVERAGE,
                                      MIX_SUM, MIX_SUM_ALIGNED],
                   help="作る見方（省略すると project.json の設定）")
    a = p.parse_args()

    project = pj.Project.load(a.folder)
    if a.mode:
        project.source_combination = a.mode
    tags = source_tags(project)
    if len(tags) < 2:
        raise SystemExit(f"{a.folder} に音源ごとの結果（結果/src1, src2 …）が"
                         f"ありません。先に計算してください")
    written = write_all(project)
    print(f"[合成] {len(written)} フォルダぶん書き出しました")


if __name__ == "__main__":
    main()
