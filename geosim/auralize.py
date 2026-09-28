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

ドライソースの置き場（★**ここに置けば一覧に並ぶ**。サブフォルダも見る）：
  リポジトリ直下の `ドライソース/`        … みんなで使う音源（**画面へ落とした音もここに保存**）
  プロジェクトフォルダ直下の `ドライソース/` … その案件だけで使う音源
★音声本体は Git に入れない（public のため。`.gitignore` で除外）。
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


def natural_key(text):
    """`条件2-1` が `条件10-1` より前に来る並べ方（数字は数として比べる）。"""
    return [(0, int(part), "") if part.isdigit() else (1, 0, part)
            for part in re.split(r"(\d+)", text or "")]


def is_source_folder(name):
    """音源ごとの棚（`src1` / `合成_平均`）か。★`旧_…` などの控えは拾わない。"""
    return bool(re.fullmatch(pj.SOURCE_DIR.replace("%d", r"\d+"), name)
                or name.startswith(pj.MIX_PREFIX))


def is_receiver_folder(name):
    return bool(re.fullmatch(r"rec\d+", name))


class Catalog:
    """プロジェクトの `結果/` を見て、畳み込める応答を並べる。

    ★**計算はやり直さない**（`ir.csv` を読むだけ）。置き場は
    `結果/[srcM/]recN/<室>_<条件>_ir.csv`。条件はファイル名から、
    音源位置と受音点はフォルダから決める。

    ★★**開いたときは並べるだけで、中身は使うときに読む**（2026-09-26。
    実案件の階段教室は応答が **1,045 本 × 7.4 MB ＝ 7.7 GB**、しかも OneDrive の
    共有フォルダで**クラウドにしか無い**ファイル。全部を先に読むと全部を落とすことになる）。
      カードの T30 / EDT … 画面に出ている分の `rt.csv`（小さい）だけ読む（`info`）
      インパルス応答     … 畳み込みに使うもの（温めておく最大 12 本）だけ読む
    """

    # 読んだ応答を覚えておく本数（1 本 = 3 秒・44.1 kHz で 1.2 MB 前後）
    KEEP = 48

    def __init__(self, folder):
        self.folder = os.path.abspath(folder)
        self.lock = threading.Lock()
        self.cache = {}          # id -> (fs, ir float32)。使った順に古いものを捨てる
        self.resampled = {}      # (id, rate) -> (bytes, level_db)
        self.infos = {}          # id -> {t30, edt}
        self.results = []
        self.by_id = {}
        self.gain = None
        self.reference = None
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
        """ファイルを並べるだけ（**中身は読まない**）。"""
        project = self._project()
        self.room = project.room_label
        self.title = project.display_name
        root = os.path.join(self.folder, pj.RESULT_DIR)
        found = []
        for directory, dirs, files in os.walk(root):
            parts = [p for p in os.path.relpath(directory, root).replace("\\", "/").split("/")
                     if p and p != "."]
            # ★棚（src / 合成）と受音点（recN）以外のフォルダには降りない
            #   （`旧_単一音源S1のみ_260915` のような控えを別の音源位置として拾わない）
            if not parts:
                dirs[:] = [d for d in dirs if is_source_folder(d) or is_receiver_folder(d)]
            elif len(parts) == 1 and is_source_folder(parts[0]):
                dirs[:] = [d for d in dirs if is_receiver_folder(d)]
            else:
                dirs[:] = []
            receiver = next((p for p in parts if is_receiver_folder(p)), "")
            source = next((p for p in parts if p != receiver), "")
            for name in files:
                if not name.endswith(IR_SUFFIX):
                    continue
                stem = name[:-len(IR_SUFFIX)].rstrip("_")
                path = os.path.join(directory, name)
                found.append({
                    "id": os.path.relpath(path, self.folder).replace("\\", "/"),
                    "condition": self._condition(stem),
                    "source": source,
                    "receiver": receiver,
                    "path": path,
                })
        # 条件は**数を数として**並べる（条件2 → 条件10）。受音点・音源は番号順
        found.sort(key=lambda r: (natural_key(r["condition"]), _source_order(r["source"]),
                                  _receiver_order(r["receiver"])))
        with self.lock:
            self.results = found
            self.by_id = {r["id"]: r for r in found}
            self.cache.clear()
            self.resampled.clear()
            self.infos.clear()
            self.gain = None
            self.reference = None
        return found

    # ---- 読む（使うときに）---------------------------------------------------

    def _load(self, identifier):
        with self.lock:
            if identifier in self.cache:
                item = self.cache.pop(identifier)
                self.cache[identifier] = item            # 使った順の末尾へ
                return item
            if identifier not in self.by_id:
                raise KeyError(identifier)
            path = self.by_id[identifier]["path"]
        fs, ir = read_ir(path)
        item = (fs, ir.astype(np.float32))
        with self.lock:
            self.cache[identifier] = item
            while len(self.cache) > self.KEEP:
                old = next(iter(self.cache))
                del self.cache[old]
                for key in [k for k in self.resampled if k[0] == old]:
                    del self.resampled[key]
        return item

    def _gain(self):
        """★**共通の 1 つの係数**。条件ごとに正規化すると音量の差が聞こえなくなる。

        基準は**並びの先頭の応答**（最初の条件・src1・rec1）。全部を読まずに
        決まり、開き直しても同じになる。応答は絶対値（√E/(t·c)）を持つので、
        どれを基準にしても**比は変わらない**。
        """
        if self.gain is None:
            if not self.results:
                return 1.0
            first = self.results[0]
            _fs, ir = self._load(first["id"])
            norm = float(np.sqrt(np.sum(ir.astype(np.float64) ** 2))) or 1.0
            self.gain = 1.0 / norm
            self.reference = first
        return self.gain

    def level_db(self, identifier):
        """基準の応答に対する大きさ [dB]（エネルギーの比）。"""
        gain = self._gain()
        _fs, ir = self._load(identifier)
        norm = float(np.sqrt(np.sum(ir.astype(np.float64) ** 2)))
        return round(float(20 * np.log10(norm * gain)), 1) if norm > 0 else None

    def info(self, identifiers):
        """カードに出す T20 / T30 / EDT（中音域の平均）。`rt.csv` だけ読む（小さい）。"""
        import table

        def one(identifier):
            with self.lock:
                if identifier in self.infos:
                    return identifier, self.infos[identifier]
                item = self.by_id.get(identifier)
            if item is None:
                return identifier, None
            rt_path = item["path"][:-len(IR_SUFFIX)] + "rt.csv"
            frequencies, rows = table.read_frequency_table(rt_path)
            value = {key: mid_value(frequencies, rows.get(row))
                     for key, row in (("t20", "T20_s"), ("t30", "T30_s"),
                                     ("edt", "EDT_s"))}
            with self.lock:
                self.infos[identifier] = value
            return identifier, value

        with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
            return dict(pool.map(one, identifiers))

    def reference_info(self):
        self._gain()
        ref = self.reference
        return ({"condition": ref["condition"], "source": ref["source"],
                 "receiver": ref["receiver"]} if ref else {})

    # ---- 渡す -------------------------------------------------------------

    def listing(self):
        public = [{k: v for k, v in r.items() if k != "path"} for r in self.results]
        return {
            "project": self.title,
            "room": self.room,
            "folder": self.folder,
            "results": public,
            "sources": dry_sources(self.folder),
        }

    def ir_bytes(self, identifier, rate):
        """ブラウザの周波数に直し、共通の係数を掛けた応答（float32）と大きさ [dB]。"""
        key = (identifier, int(rate))
        gain = self._gain()
        with self.lock:
            if key in self.resampled:
                return self.resampled[key]
        fs, ir = self._load(identifier)
        out = (resample(ir.astype(np.float64), fs, rate) * gain).astype("<f4")
        result = (out.tobytes(), self.level_db(identifier))
        with self.lock:
            self.resampled[key] = result
        return result


def dry_folders(project_folder):
    """ドライソースの置き場（プロジェクトの分 → リポジトリの分）。"""
    return [("project", os.path.join(project_folder, DRY_DIR)),
            ("common", os.path.join(REPOSITORY, DRY_DIR))]


def dry_sources(project_folder):
    """ドライソースを探す（プロジェクトの分 → リポジトリの分）。

    ★**呼ぶたびにフォルダを見直す**（控えを持たない）。利用者が自分の音源を
    置いたら、画面の「一覧を更新」だけで並ぶ（結果を読み直さなくてよい）。
    ★**サブフォルダも見る**（`ドライソース/坂吉/声.wav` → 名前は `坂吉/声`）。
    人ごと・種類ごとに分けて置けるように。
    """
    found = []
    for where, folder in dry_folders(project_folder):
        if not os.path.isdir(folder):
            continue
        for directory, dirs, files in os.walk(folder):
            dirs.sort()
            for name in sorted(files):
                if not name.lower().endswith(AUDIO_EXTENSIONS):
                    continue
                relative = os.path.relpath(os.path.join(directory, name), folder)
                relative = relative.replace("\\", "/")
                found.append({"id": f"{where}/{relative}",
                              "name": os.path.splitext(relative)[0],
                              "where": where,
                              "bytes": os.path.getsize(os.path.join(directory, name))})
    return found


def dry_path(project_folder, identifier):
    where, _, name = identifier.partition("/")
    base = dict(dry_folders(project_folder)).get(where)
    if base is None or not name:
        raise KeyError(identifier)
    # ★フォルダの外は渡さない（`..` で抜ける名前を弾く）
    root = os.path.realpath(base)
    path = os.path.realpath(os.path.join(root, name))
    if os.path.commonpath([root, path]) != root or not os.path.isfile(path):
        raise KeyError(identifier)
    return path


def save_dry(name, data):
    """画面へ落とされた音声を**共通の `ドライソース/`** に保存する。戻り値は保存した名前。

    ★**上書きしない**。同じ名前で中身が同じならそのまま、違えば `名前_2.wav` にする。
    """
    name = os.path.basename(name.replace("\\", "/")).strip()
    if not name.lower().endswith(AUDIO_EXTENSIONS) or name.startswith("."):
        raise ValueError(f"音声ファイルではありません: {name}")
    folder = os.path.join(REPOSITORY, DRY_DIR)
    os.makedirs(folder, exist_ok=True)
    stem, ext = os.path.splitext(name)
    candidate, number = name, 1
    while os.path.exists(os.path.join(folder, candidate)):
        path = os.path.join(folder, candidate)
        if os.path.getsize(path) == len(data):
            with open(path, "rb") as f:
                if f.read() == data:
                    return candidate
        number += 1
        candidate = f"{stem}_{number}{ext}"
    with open(os.path.join(folder, candidate), "wb") as f:
        f.write(data)
    print(f"[可聴化] ドライソースに加えました: {os.path.join(folder, candidate)}")
    return candidate


# ------------------------------------------------------------------------------
# サーバ
# ------------------------------------------------------------------------------

AUDIO_TYPES = {".wav": "audio/wav", ".mp3": "audio/mpeg", ".flac": "audio/flac",
               ".ogg": "audio/ogg", ".m4a": "audio/mp4"}


def make_handler(catalog, state):

    class Handler(BaseHTTPRequestHandler):

        def log_message(self, fmt, *args):      # 端末を埋めない
            pass

        def _send(self, body, kind="application/json; charset=utf-8", status=200,
                  headers=None):
            if isinstance(body, (dict, list)):
                body = json.dumps(body, ensure_ascii=False).encode("utf-8")
            elif isinstance(body, str):
                body = body.encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", kind)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            for name, value in (headers or {}).items():
                self.send_header(name, value)
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
                    data, level = catalog.ir_bytes(query["id"][0], rate)
                    self._send(data, "application/octet-stream",
                               headers={"X-Level-dB": "" if level is None else str(level)})
                elif url.path == "/api/info":
                    ids = [i for i in query.get("ids", [""])[0].split("\n") if i]
                    self._send(catalog.info(ids))
                elif url.path == "/api/reference":
                    self._send(catalog.reference_info())
                elif url.path == "/api/dry":
                    path = dry_path(catalog.folder, query["id"][0])
                    with open(path, "rb") as f:
                        kind = AUDIO_TYPES.get(os.path.splitext(path)[1].lower(),
                                               "application/octet-stream")
                        self._send(f.read(), kind)
                elif url.path == "/api/sources":
                    self._send({"sources": dry_sources(catalog.folder),
                                "folder": os.path.join(REPOSITORY, DRY_DIR)})
                elif url.path == "/api/open-dry":
                    # 置き場をエクスプローラで開く（無ければ作る）
                    folder = os.path.join(REPOSITORY, DRY_DIR)
                    os.makedirs(folder, exist_ok=True)
                    os.startfile(folder)
                    self._send({"ok": True})
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
            url = urllib.parse.urlparse(self.path)
            # 窓を閉じたときの知らせ（`navigator.sendBeacon`）。★読み込み直しでも
            #   届くので、すぐには止めず猶予を置く（次の要求が来れば続ける）
            if url.path == "/api/bye":
                state["last_seen"] = time.time() - IDLE_SECONDS + BYE_GRACE
                self._send({"ok": True})
                return
            state["last_seen"] = time.time()
            if url.path == "/api/dry-upload":
                try:
                    name = urllib.parse.parse_qs(url.query)["name"][0]
                    length = int(self.headers.get("Content-Length", "0"))
                    saved = save_dry(name, self.rfile.read(length))
                    self._send({"id": f"common/{saved}", "name": os.path.splitext(saved)[0]})
                except (KeyError, ValueError, OSError) as error:
                    self._send({"error": str(error)}, status=400)
                return
            self._send({"error": "not found"}, status=404)

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
    print(f"[可聴化] 結果を並べています（中身は使うときに読みます）: {folder}")
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
