"""結果を見比べる ― 計算した結果を、窓を並べて条件・足し方ごとに見比べる画面。

    cd geosim
    python result_viewer.py <プロジェクトフォルダ>
    python result_viewer.py <プロジェクトフォルダ> --port 8770 --no-browser

2026-09-30 ユーザー要望：
> 結果の図たちを GUI 上でも見たい。条件で比較できるようにしたい。
> GUI 上でも減衰曲線をエネルギー和や位相の切り替え、同時表示ができるようにしたい。
> 例えば、デフォルトで 4 つ窓があって、何を比較するか任意で選べるとか。
> ウィンドウ毎に条件やエネルギーか位相かを選択できるようにして。
> 結果の図は png をそのまま持ってくるのはやめて、今まで結果フォルダの図に
> 出してくれていたような描画を各ウィンドウで選べるようにして欲しい。
> ウィンドウの縦軸横軸を任意に設定できるようにして。ホールド機能も付けて。
> RTany を任意に設定できる（sti_measurement を参考に）。

★作りは可聴化（`auralize.py`）と同じ：**画面はブラウザ（Edge のアプリ窓）**、
  この Python は `結果/` を読んで渡すだけの小さなサーバ（**127.0.0.1 だけ**）。
  **追加のライブラリは入れていない**。画面は `result_viewer.html`。
★**計算はやり直さない**。位相考慮（波形）は本体が書いた `decay.csv` / `rt.csv` /
  `clarity.csv` / `spl.csv` / `ir.csv` を読む。エネルギー和は `pulses.csv` から
  その場で読む（`energy_sum.py`。案件フォルダの `エネルギー和の読み取り.py` と同じやり方）。
★**開いたときは並べるだけ、中身は使うときに読む**（実案件は OneDrive のクラウドにしか
  無いファイルがある。可聴化と同じ約束）。

「足し方」の対応：
  減衰曲線・残響時間・明瞭度 … 位相考慮＝本体（波形をフィルタに通して読む）／
                                エネルギー和＝パルス列のエネルギーをそのまま足す
  音圧レベル                  … 位相考慮＝`spl.csv` の複素和／エネルギー和＝同じ表の音圧レベル
  インパルス応答              … 位相考慮＝波形（`ir.csv`）／エネルギー和＝パルス列（反射音の並び）
  STI・到来方向・吸音率       … 足し方に依らない（STI と到来方向はエネルギーで定義される）
"""

import argparse
import csv
import json
import os
import re
import sys
import threading
import time
import urllib.parse
import warnings
from collections import OrderedDict
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import numpy as np

import auralize as au
import energy_sum as es
import project as pj
import reverberation as rv
import table as tb
from atmosphere import Atmosphere

HERE = os.path.dirname(os.path.abspath(__file__))
PAGE = os.path.join(HERE, "result_viewer.html")
# 窓の並び・RTany の区間を残すファイル（プロジェクト直下。**入力**なので頭を付けない。
# `視点.json` と同じ扱い）
LAYOUT_FILE = "見比べ.json"
# 画像で保存した窓の置き場（`図/` の下）
FIGURE_SUBDIR = "見比べ"
AVERAGE = "平均"
# 到来方向の分割（10° 刻み。`plots.DIRECTION_SECTORS` と同じ）
SECTORS = 36
EARLY = 0.050
# 包絡線の移動平均（`plots.impulse_response` と同じ 5 ms）
ENVELOPE_S = 0.005
# 画面へ渡す曲線の刻み（本体の decay.csv と同じ 1 ms）
DISPLAY_STEP = 0.001


def _round(values, digits=4):
    """JSON に載せる（NaN・inf は null）。"""
    array = np.asarray(values, dtype=float)
    out = np.round(array, digits).tolist()
    if np.all(np.isfinite(array)):
        return out

    def clean(v):
        if isinstance(v, list):
            return [clean(x) for x in v]
        return v if v is not None and np.isfinite(v) else None
    return clean(out)


def _band_key(fc):
    return str(int(round(float(fc))))


class Results:
    """`結果/` の中身を並べ、頼まれたものだけ読む。"""

    KEEP_PULSES = 16      # 読んだパルス列を覚えておく数（1 本 数 MB）
    KEEP_IR = 12

    def __init__(self, folder):
        self.folder = os.path.abspath(folder)
        self.lock = threading.Lock()
        self.pulse_cache = OrderedDict()
        self.ir_cache = OrderedDict()
        self.index_cache = OrderedDict()
        self.decay_cache = OrderedDict()
        self._bands = None
        self.scan()

    # ---- 並べる --------------------------------------------------------------

    def scan(self):
        project = self._project()
        self.project = project
        self.room = project.room_label
        self.title = project.display_name
        try:
            self.atmosphere = Atmosphere(project.temperature, project.humidity, project.pressure)
        except Exception:
            self.atmosphere = Atmosphere()
        self.rt_any = (getattr(project, "rt_any_start_db", -5.0),
                       getattr(project, "rt_any_end_db", -35.0))
        root = os.path.join(self.folder, pj.RESULT_DIR)
        entries = {}          # (shelf, rec, cond) -> {key: path}
        suffixes = ("decay.csv", "rt.csv", "clarity.csv", "spl.csv", "sti.csv",
                    "pulses.csv", "ir.csv")
        for directory, dirs, files in os.walk(root):
            parts = [p for p in os.path.relpath(directory, root).replace("\\", "/").split("/")
                     if p and p != "."]
            if not parts:
                dirs[:] = [d for d in dirs if au.is_source_folder(d) or au.is_receiver_folder(d)]
            elif len(parts) == 1 and au.is_source_folder(parts[0]):
                dirs[:] = [d for d in dirs if au.is_receiver_folder(d)]
            else:
                dirs[:] = []
            receiver = next((p for p in parts if au.is_receiver_folder(p)), "")
            if not receiver:
                continue
            shelf = next((p for p in parts if p != receiver), "")
            for name in files:
                for suffix in suffixes:
                    if name.endswith("_" + suffix) or name == suffix:
                        stem = name[:-len(suffix)].rstrip("_")
                        cond = self._condition(stem)
                        entries.setdefault((shelf, receiver, cond), {})[suffix[:-4]] = \
                            os.path.join(directory, name)
        self.entries = entries
        self.shelves = sorted({k[0] for k in entries}, key=au._source_order)
        self.conditions = sorted({k[2] for k in entries}, key=au.natural_key)
        self.receivers = {s: sorted({k[1] for k in entries if k[0] == s}, key=au._receiver_order)
                          for s in self.shelves}
        with self.lock:
            self.pulse_cache.clear()
            self.ir_cache.clear()
            self.index_cache.clear()
            self.decay_cache.clear()
        self._bands = None
        return entries

    def _project(self):
        try:
            return pj.Project.load(self.folder) if pj.Project.exists(self.folder) \
                else pj.Project(self.folder)
        except Exception as error:
            print(f"[見比べ] project.json を読めませんでした（{error}）")
            return pj.Project(self.folder)

    def _condition(self, stem):
        if self.room and stem.startswith(self.room + "_"):
            return stem[len(self.room) + 1:]
        if stem == self.room:
            return ""
        return stem

    def listing(self):
        if self._bands is None:
            self._bands = self._find_bands()
        return {
            "project": self.title, "room": self.room, "folder": self.folder,
            "shelves": self.shelves, "conditions": self.conditions,
            "receivers": self.receivers, "bands": self._bands,
            "rt_any": list(self.rt_any),
            "head_azimuth": getattr(self.project, "head_azimuth", 0.0),
        }

    def _find_bands(self):
        bands = []
        for files in self.entries.values():
            if "rt" in files:
                frequencies, _ = tb.read_frequency_table(files["rt"])
                if frequencies is not None:
                    bands = [float(f) for f in frequencies]
                    break
        return bands

    # ---- 読む ---------------------------------------------------------------

    def _files(self, shelf, rec, cond):
        files = self.entries.get((shelf, rec, cond))
        if files is None:
            raise KeyError(f"{shelf or ''}/{rec}/{cond}")
        return files

    def _cached(self, cache, key, loader, keep):
        with self.lock:
            if key in cache:
                cache.move_to_end(key)
                return cache[key]
        value = loader()
        with self.lock:
            cache[key] = value
            while len(cache) > keep:
                cache.popitem(last=False)
        return value

    def pulses(self, shelf, rec, cond):
        path = self._files(shelf, rec, cond).get("pulses")
        if not path:
            return None

        def load():
            p = es.read_pulses(path)
            p["received"] = es.received(p, self.atmosphere)
            return p
        return self._cached(self.pulse_cache, path, load, self.KEEP_PULSES)

    def impulse(self, shelf, rec, cond):
        path = self._files(shelf, rec, cond).get("ir")
        if not path:
            return None
        return self._cached(self.ir_cache, path, lambda: au.read_ir(path), self.KEEP_IR)

    def _room_path(self, cond):
        stem = f"{self.room}_{cond}" if cond else self.room
        return os.path.join(self.folder, pj.RESULT_DIR, f"{stem}_吸音率と理論値.csv")

    def room_table(self, cond):
        return tb.read_sectioned_table(self._room_path(cond))

    # ---- 減衰曲線 -----------------------------------------------------------

    def _phase_decay(self, files):
        """位相考慮の減衰曲線 (時刻の刻み, {帯域: 曲線})。

        ★**`ir.csv` から本体と同じ手順で作り直す**（`reverberation.decay_curves`。
          フィルタ・帯域幅・読み取り方も本体に合わせる）。`decay.csv` は 1 ms に
          間引いてあるので、そこで RTany を読むと本体の T20 と 0.5% ほど食い違う
          （−5〜−25 dB を入れても T20 と同じにならない）。`ir.csv` が無いときだけ
          `decay.csv` を使う。
        """
        ir_path = files.get("ir")
        if ir_path:
            def load():
                fs, ir = au.read_ir(ir_path)
                time_axis = np.arange(len(ir)) / fs
                base = rv.decay_curves(time_axis, ir.astype(float), frequencies=self.bands(),
                                       band_width=getattr(self.project, "band_width", "1/1"),
                                       fit=self.fit_method(), verbose=False)
                return 1.0 / fs, {_band_key(fc): base["decay"][j].astype(np.float32)
                                  for j, fc in enumerate(base["frequencies"])}
            return self._cached(self.decay_cache, ir_path, load, self.KEEP_IR)
        path = files.get("decay")
        if not path:
            return None
        data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2, encoding="utf-8-sig")
        with open(path, encoding="utf-8-sig") as f:
            head = f.readline().strip().split(",")
        curves = {}
        for i, name in enumerate(head):
            match = re.fullmatch(r"decay_([\d.]+)Hz_db", name)
            if match:
                curves[_band_key(match.group(1))] = data[:, i]
        dt = float(data[1, 0] - data[0, 0]) if len(data) > 1 else DISPLAY_STEP
        return dt, curves

    def bands(self):
        return self.listing()["bands"] or None

    def fit_method(self):
        return getattr(self.project, "decay_fit", None) or rv.DEFAULT_DECAY_FIT

    @staticmethod
    def _fits(curve, dt, ranges):
        out = {}
        for key, (start, end) in ranges.items():
            value, info = es.decay_time(curve, dt, start, end)
            out[key] = None if info is None else {
                "T": float(value), "a": info["slope_db_per_s"], "b": info["intercept_db"],
                "xi": info["xi"], "start": start, "end": end}
        return out

    def decay(self, shelf, rec, cond, rt_any):
        """両方の足し方の減衰曲線（1 ms 刻み）と、EDT/T20/T30/RTany の回帰。"""
        files = self._files(shelf, rec, cond)
        out = {}
        phase = self._phase_decay(files)
        if phase is not None:
            dt, curves = phase
            step = max(1, int(round(DISPLAY_STEP / dt)))
            fits = {}
            for fc, curve in curves.items():
                ranges = dict(es.RANGES)
                ranges["RTany"] = tuple(rt_any.get(fc, rt_any.get("*", self.rt_any)))
                fits[fc] = self._fits(curve.astype(float), dt, ranges)
            out["phase"] = {"t0": 0.0, "dt": dt * step,
                            "curves": {fc: _round(np.maximum(c[::step], -200.0), 2)
                                       for fc, c in curves.items()},
                            "fits": fits}
        p = self.pulses(shelf, rec, cond)
        if p is not None and len(p["time"]):
            full = es.decay_curves(p["time"], p["received"])
            step = int(round(es.FS * DISPLAY_STEP))
            fits, curves = {}, {}
            for j, fc in enumerate(p["frequencies"]):
                key = _band_key(fc)
                ranges = dict(es.RANGES)
                ranges["RTany"] = tuple(rt_any.get(key, rt_any.get("*", self.rt_any)))
                fits[key] = self._fits(full[j], 1.0 / es.FS, ranges)
                curves[key] = _round(np.maximum(full[j, ::step], -200.0), 2)
            out["energy"] = {"t0": 0.0, "dt": step / es.FS, "curves": curves, "fits": fits}
        return out

    # ---- 周波数特性（残響時間・明瞭度・音圧レベル・STI）---------------------------

    def indices(self, shelf, rec, cond, rt_any):
        """帯域ごとの指標。{bands, phase: {指標: [...]}, energy: {…}, theory, sti}。"""
        key = (shelf, rec, cond, json.dumps(rt_any, sort_keys=True))
        return self._cached(self.index_cache, key,
                            lambda: self._average(shelf, cond, rt_any) if rec == AVERAGE
                            else self._indices(shelf, rec, cond, rt_any), 256)

    def _indices(self, shelf, rec, cond, rt_any):
        files = self._files(shelf, rec, cond)
        decay = self.decay(shelf, rec, cond, rt_any)
        bands = None
        out = {"phase": {}, "energy": {}, "only": {}}
        for method in ("phase", "energy"):
            if method not in decay:
                continue
            fits = decay[method]["fits"]
            keys = list(fits)
            bands = bands or [float(k) for k in keys]
            for index in ("EDT", "T20", "T30", "RTany"):
                out[method][index] = [fits[k][index]["T"] if fits[k].get(index) else None
                                      for k in keys]
        if bands is None:
            bands = self.listing()["bands"]
        # 明瞭度：位相考慮＝clarity.csv／エネルギー和＝パルス列から
        if files.get("clarity"):
            freq, rows = tb.read_frequency_table(files["clarity"])
            for key in ("C50_db", "C80_db", "D50", "Ts_s"):
                if key in rows:
                    out["phase"][key] = self._on_bands(freq, rows[key], bands)
        p = self.pulses(shelf, rec, cond)
        if p is not None and len(p["time"]):
            for key, values in es.clarity(p["time"], p["received"]).items():
                out["energy"][key] = self._on_bands(p["frequencies"], values, bands)
        # 音圧レベル：spl.csv の「音圧レベル」がエネルギー和、「複素和」が位相考慮
        table = tb.read_sectioned_table(files.get("spl")) if files.get("spl") else None
        if table:
            freq = table["frequencies"]
            sections = table["sections"]
            if "Lp_dB" in sections.get("音圧レベル", {}):
                out["energy"]["Lp_dB"] = self._on_bands(freq, sections["音圧レベル"]["Lp_dB"], bands)
            if "Lp_dB" in sections.get("複素和", {}):
                out["phase"]["Lp_dB"] = self._on_bands(freq, sections["複素和"]["Lp_dB"], bands)
            for name, values in sections.get("内訳", {}).items():
                out["only"][name] = self._on_bands(freq, values, bands)
            if "Lp_A_dB" in sections.get("音圧レベル", {}):
                out["only"]["Lp_A_dB"] = self._on_bands(freq, sections["音圧レベル"]["Lp_A_dB"], bands)
        # STI（足し方に依らない。エネルギーで定義される）
        table = tb.read_sectioned_table(files.get("sti")) if files.get("sti") else None
        if table:
            mti = table["sections"].get("帯域別", {}).get("MTI")
            if mti is not None:
                out["only"]["MTI"] = self._on_bands(table["frequencies"], mti, bands)
            try:
                out["sti"] = float(table["values"].get("STI"))
            except (TypeError, ValueError):
                pass
        out["bands"] = bands
        out["theory"] = self._theory(cond, bands)
        return self._clean(out)

    @staticmethod
    def _on_bands(frequencies, values, bands):
        lookup = {_band_key(f): v for f, v in zip(frequencies, values)}
        return [lookup.get(_band_key(b)) for b in bands]

    def _theory(self, cond, bands):
        table = self.room_table(cond)
        if not table:
            return {}
        rows = table["sections"].get(pj.ROOM_SECTION_STATISTICAL, {})
        out = {}
        for key, label in rv.STATISTICAL_LABELS.items():
            if f"{key}_s" in rows:
                out[label] = self._on_bands(table["frequencies"], rows[f"{key}_s"], bands)
        return out

    @staticmethod
    def _clean(value):
        if isinstance(value, dict):
            return {k: Results._clean(v) for k, v in value.items()}
        if isinstance(value, np.ndarray):
            value = value.tolist()
        if isinstance(value, (list, tuple)):
            return [Results._clean(v) for v in value]
        if isinstance(value, (float, np.floating)):
            return None if not np.isfinite(value) else round(float(value), 4)
        if isinstance(value, np.integer):
            return int(value)
        return value

    def _average(self, shelf, cond, rt_any):
        """受音点の平均。音圧レベルは**エネルギー平均**、ほかは算術平均（`summary.py` と同じ）。"""
        recs = [r for r in self.receivers.get(shelf, []) if (shelf, r, cond) in self.entries]
        each = [self.indices(shelf, r, cond, rt_any) for r in recs]
        if not each:
            raise KeyError(f"{shelf}/{AVERAGE}/{cond}")
        out = {"bands": each[0]["bands"], "theory": each[0]["theory"], "count": len(each)}
        for group in ("phase", "energy", "only"):
            out[group] = {}
            keys = set().union(*(e[group].keys() for e in each))
            for key in keys:
                stack = np.array([[np.nan if v is None else v for v in e[group].get(key, [])]
                                  for e in each if e[group].get(key)], dtype=float)
                if not len(stack):
                    continue
                with np.errstate(all="ignore"), warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    if key.endswith("_dB") and key.startswith("L") or key in ("直接音_dB", "反射音_dB",
                                                                             "初期_50ms_dB", "後期_dB"):
                        mean = 10.0 * np.log10(np.nanmean(10.0 ** (stack / 10.0), axis=0))
                    else:
                        mean = np.nanmean(stack, axis=0)
                out[group][key] = mean
        stis = [e["sti"] for e in each if e.get("sti") is not None]
        if stis:
            out["sti"] = float(np.mean(stis))
        return self._clean(out)

    # ---- 時間の図（インパルス応答・エネルギー時間曲線）-------------------------

    def impulse_view(self, shelf, rec, cond, band):
        """位相考慮＝波形（帯域指定ならフィルタに通す）／エネルギー和＝パルス列。"""
        out = {}
        loaded = self.impulse(shelf, rec, cond)
        if loaded is not None:
            fs, ir = loaded
            ir = ir.astype(float)
            if band and band != "all":
                ir = rv.octave_bandpass(ir, float(band), fs,
                                        band_width=getattr(self.project, "band_width", "1/1"))
            # ★**帯域を通したあとの最大で割る**（2026-09-30 ユーザー指摘「正規化されてなくない？」）。
            #   以前は広帯域の最大で割っていたので、帯域を選ぶと最大が 1 にならなかった
            peak = float(np.max(np.abs(ir))) or 1.0
            out["phase"] = {"t0": 0.0, "dt": 1.0 / fs, "y": _round(ir / peak, 5)}
        p = self.pulses(shelf, rec, cond)
        if p is not None and len(p["time"]):
            energy = p["received"]
            e = energy.sum(axis=1) if not band or band == "all" else \
                energy[:, self._band_column(p, band)]
            ref = float(e.max()) or 1.0
            with np.errstate(divide="ignore"):
                db = 10.0 * np.log10(np.maximum(e / ref, 1e-12))
            order = np.argsort(p["time"])
            out["energy"] = {"t": _round(p["time"][order], 6), "db": _round(db[order], 2)}
        return out

    @staticmethod
    def _band_column(p, band):
        keys = [_band_key(f) for f in p["frequencies"]]
        return keys.index(_band_key(band))

    def time_curve(self, shelf, rec, cond, band):
        """エネルギー時間曲線（1 ms 刻み・最大を 0 dB）。"""
        out = {}
        loaded = self.impulse(shelf, rec, cond)
        if loaded is not None:
            fs, ir = loaded
            ir = ir.astype(float)
            if band and band != "all":
                ir = rv.octave_bandpass(ir, float(band), fs,
                                        band_width=getattr(self.project, "band_width", "1/1"))
            window = max(1, int(round(ENVELOPE_S * fs)))
            envelope = np.convolve(ir ** 2, np.ones(window) / window, mode="same")
            step = max(1, int(round(DISPLAY_STEP * fs)))
            # 刻みの中の最大（間引きで山を落とさない）
            n = len(envelope) // step
            env = envelope[:n * step].reshape(n, step).max(axis=1)
            with np.errstate(divide="ignore"):
                db = 10.0 * np.log10(np.maximum(env / (env.max() or 1.0), 1e-12))
            out["phase"] = {"t0": 0.0, "dt": step / fs, "y": _round(db, 2)}
        p = self.pulses(shelf, rec, cond)
        if p is not None and len(p["time"]):
            energy = p["received"]
            e = energy.sum(axis=1, keepdims=True) if not band or band == "all" else \
                energy[:, [self._band_column(p, band)]]
            # ★位相考慮と同じ 5 ms の移動平均をかける（見比べられるように）
            grid = es.energy_grid(p["time"], e, 1.0 / DISPLAY_STEP)[0]
            window = max(1, int(round(ENVELOPE_S / DISPLAY_STEP)))
            smooth = np.convolve(grid, np.ones(window) / window, mode="same")
            with np.errstate(divide="ignore"):
                db = 10.0 * np.log10(np.maximum(smooth / (smooth.max() or 1.0), 1e-12))
            out["energy"] = {"t0": 0.0, "dt": DISPLAY_STEP, "y": _round(np.maximum(db, -120.0), 2)}
        return out

    # ---- 到来方向 ------------------------------------------------------------

    def direction(self, shelf, rec, cond, band, part):
        p = self.pulses(shelf, rec, cond)
        if p is None or p.get("direction") is None:
            return {}
        energy = p["received"]
        e = energy.sum(axis=1) if not band or band == "all" else energy[:, self._band_column(p, band)]
        onset = float(p["time"].min())
        if part == "early":
            e = np.where(p["time"] < onset + EARLY, e, 0.0)
        elif part == "late":
            e = np.where(p["time"] >= onset + EARLY, e, 0.0)
        head = self._head_azimuth(rec)
        direction = p["direction"]
        azimuth = np.mod(np.arctan2(direction[:, 1], direction[:, 0]) - np.deg2rad(head),
                         2.0 * np.pi)
        edges = np.linspace(0.0, 2.0 * np.pi, SECTORS + 1)
        totals, _ = np.histogram(azimuth, bins=edges, weights=e)
        # 直接音の向き（最初のパルス）
        first = int(np.argmin(p["time"]))
        direct = float(np.rad2deg(azimuth[first]))
        total_all = float(energy.sum()) if energy.size else 1.0
        return {"sectors": SECTORS, "energy": _round(totals / (total_all or 1.0), 8),
                "direct_deg": round(direct, 1), "head_azimuth": head}

    def _head_azimuth(self, rec):
        value = getattr(self.project, "head_azimuth", 0.0)
        if value is None:
            return 0.0
        if np.isscalar(value):
            return float(value)
        digits = re.sub(r"\D", "", rec)
        index = int(digits) - 1 if digits else 0
        return float(value[index]) if 0 <= index < len(value) else float(value[0])

    # ---- 吸音率 --------------------------------------------------------------

    def absorption(self, cond):
        table = self.room_table(cond)
        if not table:
            return {}
        sections = table["sections"]
        # ★面積（3 列目）はファイルから直に読む。`read_sectioned_table` の 'values' は
        #   項目名で引くので、同じ名前が後の区分（垂直入射など。面積は空欄）に
        #   出てくると上書きされて空になる
        areas = {}
        with open(self._room_path(cond), encoding="utf-8-sig", newline="") as f:
            for row in csv.reader(f):
                if len(row) > 2 and row[0].strip() == pj.ROOM_SECTION_MATERIALS:
                    try:
                        areas[row[1].strip()] = float(row[2])
                    except ValueError:
                        pass
        materials = []
        for section, item in table["order"]:
            if section == pj.ROOM_SECTION_MATERIALS:
                materials.append({"name": item, "area": areas.get(item),
                                  "alpha": _round(sections[section][item], 4)})
        mean = sections.get(pj.ROOM_SECTION_MEAN, {}).get("平均吸音率")
        return {"bands": [float(f) for f in table["frequencies"]], "materials": materials,
                "mean": None if mean is None else _round(mean, 4)}

    # ---- 窓の並び・画像 ------------------------------------------------------

    def read_layout(self):
        path = os.path.join(self.folder, LAYOUT_FILE)
        if not os.path.exists(path):
            return {}
        try:
            with open(path, encoding="utf-8") as f:
                return json.load(f)
        except (OSError, ValueError) as error:
            print(f"[見比べ] {LAYOUT_FILE} を読めませんでした（{error}）")
            return {}

    def write_layout(self, data):
        path = os.path.join(self.folder, LAYOUT_FILE)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=1)
        return path

    def save_figure(self, name, data):
        name = re.sub(r'[\\/:*?"<>|]', "_", name).strip() or "見比べ"
        folder = os.path.join(self.folder, pj.FIGURE_DIR if hasattr(pj, "FIGURE_DIR") else "図",
                              FIGURE_SUBDIR)
        os.makedirs(folder, exist_ok=True)
        path = os.path.join(folder, name + ".png")
        stem, number = path[:-4], 2
        while os.path.exists(path):          # ★上書きしない
            path = f"{stem}_{number}.png"
            number += 1
        with open(path, "wb") as f:
            f.write(data)
        return path


def make_handler(results, state):

    class Handler(BaseHTTPRequestHandler):

        def log_message(self, fmt, *args):
            pass

        def _send(self, body, kind="application/json; charset=utf-8", status=200):
            if isinstance(body, (dict, list)):
                body = json.dumps(body, ensure_ascii=False, allow_nan=False).encode("utf-8")
            elif isinstance(body, str):
                body = body.encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", kind)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            try:
                self.wfile.write(body)
            except (BrokenPipeError, ConnectionResetError):
                pass

        def do_GET(self):
            state["last_seen"] = time.time()
            url = urllib.parse.urlparse(self.path)
            q = {k: v[0] for k, v in urllib.parse.parse_qs(url.query).items()}
            try:
                if url.path in ("/", "/index.html"):
                    with open(PAGE, "rb") as f:
                        self._send(f.read(), "text/html; charset=utf-8")
                    return
                if url.path == "/api/catalog":
                    if "rescan" in q:
                        results.scan()
                    self._send(results.listing())
                    return
                if url.path == "/api/layout":
                    self._send(results.read_layout())
                    return
                if url.path == "/api/ping":
                    self._send({"ok": True})
                    return
                shelf, rec, cond = q.get("shelf", ""), q.get("rec", ""), q.get("cond", "")
                rt_any = json.loads(q.get("rtany", "{}") or "{}")
                band = q.get("band", "all")
                if url.path == "/api/decay":
                    self._send(results.decay(shelf, rec, cond, rt_any))
                elif url.path == "/api/indices":
                    self._send(results.indices(shelf, rec, cond, rt_any))
                elif url.path == "/api/impulse":
                    self._send(results.impulse_view(shelf, rec, cond, band))
                elif url.path == "/api/etc":
                    self._send(results.time_curve(shelf, rec, cond, band))
                elif url.path == "/api/direction":
                    self._send(results.direction(shelf, rec, cond, band, q.get("part", "all")))
                elif url.path == "/api/absorption":
                    self._send(results.absorption(cond))
                else:
                    self._send({"error": "not found"}, status=404)
            except KeyError as error:
                self._send({"error": f"結果が見つかりません: {error}"}, status=404)
            except Exception as error:
                print(f"[見比べ] {url.path}: {error}")
                self._send({"error": str(error)}, status=500)

        def do_POST(self):
            url = urllib.parse.urlparse(self.path)
            if url.path == "/api/bye":
                state["last_seen"] = time.time() - au.IDLE_SECONDS + au.BYE_GRACE
                self._send({"ok": True})
                return
            state["last_seen"] = time.time()
            length = int(self.headers.get("Content-Length", "0"))
            body = self.rfile.read(length)
            try:
                if url.path == "/api/layout":
                    path = results.write_layout(json.loads(body.decode("utf-8")))
                    self._send({"path": path})
                elif url.path == "/api/figure":
                    name = urllib.parse.parse_qs(url.query).get("name", ["見比べ"])[0]
                    self._send({"path": results.save_figure(name, body)})
                else:
                    self._send({"error": "not found"}, status=404)
            except Exception as error:
                self._send({"error": str(error)}, status=400)

    return Handler


def serve(folder, port=0, browser=True, stay=False):
    print(f"[見比べ] 結果を並べています（中身は使うときに読みます）: {folder}")
    results = Results(folder)
    print(f"[見比べ] 条件 {len(results.conditions)} ／ 置き場 {len(results.entries)} か所")
    if not results.entries:
        print("[見比べ] ★結果がありません。先に計算してください（画面は開きます）")
    state = {"last_seen": time.time()}
    server = ThreadingHTTPServer(("127.0.0.1", port), make_handler(results, state))
    url = f"http://127.0.0.1:{server.server_address[1]}/"
    print(f"[見比べ] {url}  （止めるときは Ctrl+C か、画面の窓を閉じる）")
    if browser:
        au.open_window(url)
    if not stay:
        def watchdog():
            while True:
                time.sleep(2.0)
                if time.time() - state["last_seen"] > au.IDLE_SECONDS:
                    print("[見比べ] 画面が閉じられたので終わります")
                    server.shutdown()
                    return
        threading.Thread(target=watchdog, daemon=True).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


def main(argv=None):
    p = argparse.ArgumentParser(description="計算結果を窓を並べて見比べる")
    p.add_argument("folder", nargs="?", help="プロジェクトフォルダ（結果/ があるところ）")
    p.add_argument("--port", type=int, default=0, help="待ち受けるポート（既定は空きを自動）")
    p.add_argument("--no-browser", action="store_true", help="画面を開かない")
    p.add_argument("--stay", action="store_true",
                   help="画面を閉じてもサーバを止めない（Ctrl+C で止める）")
    args = p.parse_args(argv)
    folder = args.folder
    if not folder:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        import dialog_dirs as dd
        import setup_window as sw
        # ★前回のプロジェクトのフォルダから開く（全アプリ共通の約束。`dialog_dirs.py`）
        dd.enable_persistence(sw.DIALOG_APP)
        folder = filedialog.askdirectory(title="見比べるプロジェクトフォルダ",
                                         initialdir=dd.base_dir(dd.PROJECT) or None)
        root.destroy()
        if not folder:
            return 1
        dd.remember(dd.PROJECT, folder)
    if not os.path.isdir(folder):
        print(f"[見比べ] フォルダがありません: {folder}")
        return 1
    serve(folder, port=args.port, browser=not args.no_browser,
          stay=args.stay or args.no_browser)
    return 0


if __name__ == "__main__":
    sys.exit(main())
