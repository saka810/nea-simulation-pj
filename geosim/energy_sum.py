"""パルス列を「エネルギー和」で読む（位相を見ない読み方）。**計算はやり直さない。**

本体の残響時間・明瞭度は、パルス列を**位相ごと重ねた波形**（インパルス応答）を
オクターブフィルタに通して読んでいる。低音の固有モードのうなりや、直接音と
床反射の干渉（くし形の谷）がそのまま入る。ここでは同じパルス列を
**エネルギーのまま足して**減衰曲線を作り、同じ物差し（ISO 3382 の最小二乗回帰。
`reverberation._decay_time`）で読む。

2026-09-24 から案件フォルダの `エネルギー和の読み取り.py` で手で回していたものを、
結果を見比べる画面（`result_viewer.py`）から使えるように本体へ移した（2026-09-30）。
★**やり方は同じ**（数字も同じになる）：
  - パルスごとのエネルギー … `sound_level.received_energy`
    （空気吸収 ＋ 1/(4πd²)。音圧レベル・STI と同じ足し方）
  - 時間の格子 … 44.1 kHz（本体のインパルス応答と同じ刻み）に `bincount` で置く
  - 明瞭度 … 直接音の到来を時刻 0 として 50 / 80 ms で分ける
★**既存の結果（rt.csv / clarity.csv / まとめ表）は変えない**。読み方を 1 つ足すだけ。
"""

import numpy as np

import reverberation as rv
import sound_level as sl

# 本体のインパルス応答と同じ刻み
FS = 44100.0
# 評価区間（本体と同じ）
RANGES = {"EDT": (0.0, -10.0), "T20": (-5.0, -25.0), "T30": (-5.0, -35.0)}


def read_pulses(path):
    """`pulses.csv` を読む。戻り値 dict（time / distance / direction / energy / frequencies）。"""
    with open(path, encoding="utf-8-sig") as f:
        head = f.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2, encoding="utf-8-sig")
    columns = [i for i, name in enumerate(head) if name.startswith("energy_")]
    frequencies = np.array([float(head[i][len("energy_"):].rstrip("Hz")) for i in columns])
    out = {
        "time": data[:, head.index("time_s")],
        "energy": data[:, columns],
        "frequencies": frequencies,
    }
    out["distance"] = data[:, head.index("distance_m")] if "distance_m" in head else None
    if all(k in head for k in ("dir_x", "dir_y", "dir_z")):
        out["direction"] = data[:, [head.index("dir_x"), head.index("dir_y"),
                                    head.index("dir_z")]]
    else:
        out["direction"] = None
    if "reflection_count" in head:
        out["reflections"] = data[:, head.index("reflection_count")]
    return out


def received(pulses, atmosphere):
    """受音点で受け取るエネルギー (n, nf)。音圧レベル・STI と同じ足し方。"""
    return sl.received_energy(pulses["time"], pulses["energy"],
                              distances=pulses["distance"], atmosphere=atmosphere,
                              frequencies=pulses["frequencies"])


def energy_grid(time, energy, fs=FS):
    """パルスを時間の格子に置く (nf, n)。"""
    time = np.asarray(time, float)
    n = int(np.ceil(time.max() * fs)) + 2 if len(time) else 2
    index = np.minimum(np.round(time * fs).astype(int), n - 1)
    return np.array([np.bincount(index, weights=energy[:, j], minlength=n)
                     for j in range(energy.shape[1])])


def decay_curves(time, energy, fs=FS):
    """Schroeder の逆積分 [dB] (nf, n)。先頭が 0 dB。"""
    grid = energy_grid(time, energy, fs)
    tail = np.cumsum(grid[:, ::-1], axis=1)[:, ::-1]
    with np.errstate(divide="ignore", invalid="ignore"):
        return 10.0 * np.log10(np.maximum(tail / np.maximum(tail[:, :1], 1e-300), 1e-300))


def decay_time(curve, dt, start_db, end_db):
    """(T [s], 回帰の情報 | None)。`reverberation._decay_time` そのもの。"""
    return rv._decay_time(curve, dt, start_db, end_db, detail=True)


def clarity(time, energy):
    """C50 / C80 [dB]・D50 [-]・Ts [s]（帯域ごと）。直接音の到来を 0 とする。"""
    time = np.asarray(time, float)
    onset = time.min() if len(time) else 0.0
    total = energy.sum(axis=0)
    out = {}
    for key, limit in (("C50_db", 0.050), ("C80_db", 0.080)):
        early = energy[time < onset + limit].sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            out[key] = 10.0 * np.log10(early / (total - early))
    early50 = energy[time < onset + 0.050].sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        out["D50"] = early50 / total
        out["Ts_s"] = ((time - onset)[:, None] * energy).sum(axis=0) / total
    return out


def time_curve(time, energy, step=0.001):
    """エネルギー時間曲線（`step` ごとに足したもの）[dB]、帯域ごとに最大を 0 dB。"""
    grid = energy_grid(time, energy, 1.0 / step)
    peak = np.maximum(grid.max(axis=1, keepdims=True), 1e-300)
    with np.errstate(divide="ignore"):
        return np.arange(grid.shape[1]) * step, 10.0 * np.log10(np.maximum(grid / peak, 1e-30))
