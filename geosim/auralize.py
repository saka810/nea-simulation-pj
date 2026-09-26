"""可聴化 ― 計算したインパルス応答をドライソースに畳み込んで聴く。

    cd geosim
    python auralize.py <プロジェクトフォルダ>
    python auralize.py <プロジェクトフォルダ> --port 8765 --no-browser

2026-09-26 ユーザー要望：
> 計算したインパルス応答をドライソースに畳み込みをして，可聴化したい。
> GUI 上で，選択した音源に，選択した条件の結果を畳み込むイメージです。
> 再生しながらも，畳み込む ON/OFF の切り替えや，条件の切り替えをしても，
> 再生を続けながらも即座に切り替わるようにして欲しい。
> もし可能であれば，マイク入力を畳み込んでスピーカ等へ出力できる機能も欲しい。

★**画面はブラウザ（Edge / Chrome）、音は Web Audio の `ConvolverNode`**。
  リアルタイムの畳み込み（分割 FFT）・マイク入力・波形の描画がブラウザに
  そろっているので、**追加のライブラリを入れずに**作れる（`docs/build_pdf.py`
  と同じ方針）。この Python はプロジェクトの結果を読んで渡すだけの小さな
  サーバで、**127.0.0.1 だけ**で待ち受ける（外からは見えない）。
  マイクは `getUserMedia` が要るが、localhost は安全な場所として扱われるので使える。

★★**条件どうしの音量の比は保つ**（`common_gain`）。インパルス応答は
  `振幅 = √E/(t·c)` で**絶対値を持っている**ので、条件ごとに正規化すると
  「吸音を増やすと静かになる」が聞こえなくなる。全結果の中で
  いちばんエネルギーの大きい応答を 1 とする**共通の 1 つの係数**を掛ける。
  ブラウザ側も `ConvolverNode.normalize = false`（既定の true だと応答ごとに
  正規化されて同じことが起きる）。

★サンプリング周波数はブラウザの出力に合わせて**ここで変換してから**渡す
  （`ConvolverNode` は再生側と同じ周波数の応答しか受け付けない）。
  ドライソースはブラウザが読み込むときに自分で合わせる。

ドライソースはリポジトリ直下の `ドライソース/` と、
プロジェクトフォルダ直下の `ドライソース/`（あれば）から探す。
"""

import argparse
import concurrent.futures
import json
import os
import re
import subprocess
import sys
import threading
import time
import urllib.parse
import webbrowser
from fractions import Fraction
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import numpy as np

import project as pj

HERE = os.path.dirname(os.path.abspath(__file__))
REPOSITORY = os.path.dirname(HERE)
PAGE = os.path.join(HERE, "auralize.html")
DRY_DIR = "ドライソース"
AUDIO_EXTENSIONS = (".wav", ".mp3", ".flac", ".ogg", ".m4a")
IR_SUFFIX = "ir.csv"
# 画面の減衰曲線（カード）の点の数と下限
DECAY_POINTS = 96
DECAY_FLOOR_DB = -60.0
# 残響時間を代表させる帯域（中音域の平均）
MID_BANDS = (500.0, 1000.0)
# 画面から便りが途絶えてサーバを止めるまでの秒数（`--stay` で止めない）。
# ★短くしない：窓を最小化するとブラウザは定期処理を**1 分に 1 回**まで間引く。
#   窓を閉じたときは画面が `/api/bye` を送るので、すぐ（数秒で）止まる
IDLE_SECONDS = 150.0
BYE_GRACE = 5.0
# 画面をアプリ風の窓で開くブラウザ（見つかった順に使う。`docs/build_pdf.py` と同じ）
BROWSERS = [
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Google\Chrome\Application\chrome.exe",
    r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
]


# ------------------------------------------------------------------------------
# インパルス応答を読む
# ------------------------------------------------------------------------------

def read_ir(path):
    """`impulse.write_impulse_response` の CSV → (サンプリング周波数, 応答 float64)。"""
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    t, ir = data[:, 0], data[:, 1]
    fs = 1.0 / np.median(np.diff(t[:64])) if len(t) > 1 else 44100.0
    return float(round(fs)), ir


def schroeder_db(ir, points=DECAY_POINTS, floor=DECAY_FLOOR_DB):
    """全帯域のシュレーダー積分 [dB]（画面の小さな減衰曲線に使う。0 dB 始まり）。"""
    energy = np.cumsum((ir.astype(np.float64) ** 2)[::-1])[::-1]
    if energy[0] <= 0:
        return [floor] * points
    level = 10.0 * np.log10(np.maximum(energy / energy[0], 1e-30))
    index = np.linspace(0, len(level) - 1, points).astype(int)
    return [round(max(float(v), floor), 2) for v in level[index]]


def mid_value(frequencies, values):
    """中音域（500 / 1k Hz）の平均。無ければ None。"""
    if frequencies is None or values is None:
        return None
    picked = [v for f, v in zip(frequencies, values)
              if any(abs(f - m) < 1.0 for m in MID_BANDS) and np.isfinite(v)]
    return round(float(np.mean(picked)), 3) if picked else None


def resample(ir, fs_in, fs_out):
    """`fs_in` → `fs_out`。同じならそのまま（有理数の比で多相フィルタ）。

    ★**`fs_in / fs_out` を掛ける**。畳み込みは「サンプルの和」なので、
    同じ応答でも 1 秒あたりのサンプルが増えるとそのぶん大きく出る
    （44.1 → 48 kHz で +0.74 dB）。`resample_poly` は波形の高さを保つので、
    和（＝低域の利得）を保つにはこの係数が要る
    """
    if int(fs_in) == int(fs_out):
        return ir
    from scipy.signal import resample_poly
    ratio = Fraction(int(fs_out), int(fs_in)).limit_denominator(1000)
    return (resample_poly(ir, ratio.numerator, ratio.denominator)
            * (float(fs_in) / float(fs_out)))


# ------------------------------------------------------------------------------
# プロジェクトの結果を並べる
# ------------------------------------------------------------------------------

def _receiver_order(name):
    digits = re.sub(r"\D", "", name)
    return (0, int(digits)) if digits else (1, name)


def _source_order(name):
    if not name:
        return (0, 0, "")
    if re.fullmatch(r"src\d+", name):
        return (1, int(name[3:]), name)
    return (2, 0, name)


class Catalog:
    """プロジェクトの `結果/` を見て、畳み込める応答を並べる。

    ★**計算はやり直さない**（`ir.csv` を読むだけ）。置き場は
    `結果/[srcM/]recN/<室>_<条件>_ir.csv`。条件はファイル名から、
    音源位置と受音点はフォルダから決める。
    """

    def __init__(self, folder):
        self.folder = os.path.abspath(folder)
        self.lock = threading.Lock()
        self.cache = {}          # id -> (fs, ir float32)
        self.resampled = {}      # (id, rate) -> bytes
        self.results = []
        self.gain = 1.0
        self.room = ""
        self.scan()

    # ---- 探す -------------------------------------------------------------

    def _project(self):
        try:
            return pj.Project.load(self.folder) if pj.Project.exists(self.folder) \
                else pj.Project(self.folder)
        except Exception as error:            # 読めなくても結果は並べられる
            print(f"[可聴化] project.json を読めませんでした（{error}）。"
                  "フォルダ名を室名として扱います")
            return pj.Project(self.folder)

    def _condition(self, stem):
        """`研修室_条件A` → `条件A`（室名の頭を外す）。"""
        if self.room and stem.startswith(self.room + "_"):
            return stem[len(self.room) + 1:]
        if stem == self.room:
            return ""
        return stem

    def scan(self):
        project = self._project()
        self.room = project.room_label
        self.title = project.display_name
        root = os.path.join(self.folder, pj.RESULT_DIR)
        found = []
        for directory, _dirs, files in os.walk(root):
            for name in files:
                if not name.endswith(IR_SUFFIX):
                    continue
                stem = name[:-len(IR_SUFFIX)].rstrip("_")
                parts = os.path.relpath(directory, root).replace("\\", "/").split("/")
                parts = [p for p in parts if p and p != "."]
                receiver = next((p for p in parts if re.fullmatch(r"rec\d+", p)), "")
                source = next((p for p in parts if p != receiver), "")
                path = os.path.join(directory, name)
                found.append({
                    "id": os.path.relpath(path, self.folder).replace("\\", "/"),
                    "condition": self._condition(stem),
                    "source": source,
                    "receiver": receiver,
                    "path": path,
                })
        # 条件は**条件表のシートの順**に並べたいが、結果だけからは分からないので
        # 見つかった名前の順（五十音）にする。受音点・音源は番号順
        found.sort(key=lambda r: (r["condition"], _source_order(r["source"]),
                                  _receiver_order(r["receiver"])))
        self._measure(found)
        self.results = found
        return found

    def _measure(self, results):
        """応答ごとの大きさ・減衰曲線・残響時間を読む（並列）。"""
        import table

        def one(item):
            fs, ir = read_ir(item["path"])
            with self.lock:
                self.cache[item["id"]] = (fs, ir.astype(np.float32))
            item["fs"] = fs
            item["seconds"] = round(len(ir) / fs, 3)
            item["norm"] = float(np.sqrt(np.sum(ir.astype(np.float64) ** 2)))
            item["decay"] = schroeder_db(ir)
            rt_path = item["path"][:-len(IR_SUFFIX)] + "rt.csv"
            frequencies, rows = table.read_frequency_table(rt_path)
            for key, row in (("t30", "T30_s"), ("edt", "EDT_s")):
                item[key] = mid_value(frequencies, rows.get(row))
            return item

        with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
            list(pool.map(one, results))
        norms = [r["norm"] for r in results if r["norm"] > 0]
        # ★**共通の 1 つの係数**。条件ごとに正規化すると音量の差が聞こえなくなる
        reference = max(norms) if norms else 1.0
        self.gain = 1.0 / reference
        for r in results:
            r["level_db"] = (round(20 * np.log10(r["norm"] / reference), 1)
                             if r["norm"] > 0 else None)
        self.resampled.clear()

    # ---- 渡す -------------------------------------------------------------

    def listing(self):
        public = [{k: v for k, v in r.items() if k not in ("path", "norm")}
                  for r in self.results]
        return {
            "project": self.title,
            "room": self.room,
            "folder": self.folder,
            "results": public,
            "sources": dry_sources(self.folder),
        }

    def ir_bytes(self, identifier, rate):
        key = (identifier, int(rate))
        with self.lock:
            if key in self.resampled:
                return self.resampled[key]
            if identifier not in self.cache:
                raise KeyError(identifier)
            fs, ir = self.cache[identifier]
        out = (resample(ir.astype(np.float64), fs, rate) * self.gain).astype("<f4")
        data = out.tobytes()
        with self.lock:
            self.resampled[key] = data
        return data


def dry_sources(project_folder):
    """ドライソースを探す（プロジェクトの分 → リポジトリの分）。"""
    places = [("project", os.path.join(project_folder, DRY_DIR)),
              ("common", os.path.join(REPOSITORY, DRY_DIR))]
    found = []
    for where, folder in places:
        if not os.path.isdir(folder):
            continue
        for name in sorted(os.listdir(folder)):
            if name.lower().endswith(AUDIO_EXTENSIONS):
                found.append({"id": f"{where}/{name}",
                              "name": os.path.splitext(name)[0],
                              "where": where,
                              "bytes": os.path.getsize(os.path.join(folder, name))})
    return found


def dry_path(project_folder, identifier):
    where, _, name = identifier.partition("/")
    base = {"project": os.path.join(project_folder, DRY_DIR),
            "common": os.path.join(REPOSITORY, DRY_DIR)}.get(where)
    # ★フォルダの外は渡さない（`..` を含む名前を弾く）
    if base is None or os.path.basename(name) != name:
        raise KeyError(identifier)
    path = os.path.join(base, name)
    if not os.path.isfile(path):
        raise KeyError(identifier)
    return path


# ------------------------------------------------------------------------------
# サーバ
# ------------------------------------------------------------------------------

AUDIO_TYPES = {".wav": "audio/wav", ".mp3": "audio/mpeg", ".flac": "audio/flac",
               ".ogg": "audio/ogg", ".m4a": "audio/mp4"}


def make_handler(catalog, state):

    class Handler(BaseHTTPRequestHandler):

        def log_message(self, fmt, *args):      # 端末を埋めない
            pass

        def _send(self, body, kind="application/json; charset=utf-8", status=200):
            if isinstance(body, (dict, list)):
                body = json.dumps(body, ensure_ascii=False).encode("utf-8")
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
            query = urllib.parse.parse_qs(url.query)
            try:
                if url.path in ("/", "/index.html"):
                    with open(PAGE, "rb") as f:
                        self._send(f.read(), "text/html; charset=utf-8")
                elif url.path == "/api/catalog":
                    if "rescan" in query:
                        catalog.scan()
                    self._send(catalog.listing())
                elif url.path == "/api/ir":
                    rate = int(float(query.get("rate", ["48000"])[0]))
                    data = catalog.ir_bytes(query["id"][0], rate)
                    self._send(data, "application/octet-stream")
                elif url.path == "/api/dry":
                    path = dry_path(catalog.folder, query["id"][0])
                    with open(path, "rb") as f:
                        kind = AUDIO_TYPES.get(os.path.splitext(path)[1].lower(),
                                               "application/octet-stream")
                        self._send(f.read(), kind)
                elif url.path == "/api/ping":
                    self._send({"ok": True})
                else:
                    self._send({"error": "not found"}, status=404)
            except KeyError as error:
                self._send({"error": f"見つかりません: {error}"}, status=404)
            except Exception as error:
                print(f"[可聴化] {url.path}: {error}")
                self._send({"error": str(error)}, status=500)

        def do_POST(self):
            # 窓を閉じたときの知らせ（`navigator.sendBeacon`）。★読み込み直しでも
            #   届くので、すぐには止めず猶予を置く（次の要求が来れば続ける）
            if urllib.parse.urlparse(self.path).path == "/api/bye":
                state["last_seen"] = time.time() - IDLE_SECONDS + BYE_GRACE
            self._send({"ok": True})

    return Handler


def find_browser():
    for path in BROWSERS:
        if os.path.exists(path):
            return path
    return None


def open_window(url):
    """アプリ風の窓（アドレス欄なし）で開く。ブラウザが見つからなければ既定の開き方。"""
    browser = find_browser()
    if browser:
        try:
            subprocess.Popen([browser, f"--app={url}", "--window-size=1440,920",
                              "--autoplay-policy=no-user-gesture-required"])
            return
        except OSError:
            pass
    webbrowser.open(url)


def serve(folder, port=0, browser=True, stay=False):
    """サーバを立てて画面を開く。窓が閉じたら（応答が途絶えたら）止まる。"""
    print(f"[可聴化] 結果を読んでいます: {folder}")
    catalog = Catalog(folder)
    print(f"[可聴化] インパルス応答 {len(catalog.results)} 本 / "
          f"ドライソース {len(dry_sources(catalog.folder))} 個")
    if not catalog.results:
        print("[可聴化] ★結果にインパルス応答（*_ir.csv）がありません。"
              "先に計算してください（画面は開きます）")
    state = {"last_seen": time.time()}
    server = ThreadingHTTPServer(("127.0.0.1", port), make_handler(catalog, state))
    url = f"http://127.0.0.1:{server.server_address[1]}/"
    print(f"[可聴化] {url}  （止めるときは Ctrl+C か、画面の窓を閉じる）")
    if browser:
        open_window(url)
    if not stay:
        def watchdog():
            while True:
                time.sleep(2.0)
                if time.time() - state["last_seen"] > IDLE_SECONDS:
                    print("[可聴化] 画面が閉じられたので終わります")
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
    p = argparse.ArgumentParser(description="インパルス応答をドライソースに畳み込んで聴く")
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
        folder = filedialog.askdirectory(title="可聴化するプロジェクトフォルダ")
        root.destroy()
        if not folder:
            return 1
    if not os.path.isdir(folder):
        print(f"[可聴化] フォルダがありません: {folder}")
        return 1
    serve(folder, port=args.port, browser=not args.no_browser,
          stay=args.stay or args.no_browser)
    return 0


if __name__ == "__main__":
    sys.exit(main())
