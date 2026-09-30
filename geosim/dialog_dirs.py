"""ファイルダイアログの初期フォルダを覚える（Qt 非依存）。全アプリ共通の約束。

**`Desktop/Claude/python/meas_project/common/dialog_dirs.py` の写し**（2026-09-30 ユーザー指示
「Claude フォルダのルールに則って、ファイルダイアログの初期フォルダのルールに合わせて」）。
このリポジトリは別の venv・別の Git なので写して持っている。挙動は元と同じに保つこと
（直すときは元も揃える。`Desktop/Claude/CLAUDE.md`「ファイルダイアログの初期フォルダ」）。
geosim では ``import dialog_dirs as dd``。記録は ``%LOCALAPPDATA%/nea-simulation-pj/last_dirs.json``。

以下は元の説明。

DA 分析の `projects/da_analyzer/paths.py` と同じ考え方を、どのアプリからも使えるように
したもの（2026-09-26。ルールはワークスペースの CLAUDE.md「ファイルダイアログの初期フォルダ」）。
DA 分析は自前の `paths.py` をそのまま使っている（中身は同じ約束）。

- **用途ごとに最後に使ったフォルダ**を覚える（CSV を読む・出力先・プロジェクト …）。
  毎回「ドキュメント」や EXE の直下から辿り直さずに済む
- 覚えていない用途は**プロジェクトのフォルダ**から開く（全体の拠り所）
- **プロジェクトを開いた・保存したら用途ごとの記憶を捨て**、そのプロジェクトに書かれて
  いるフォルダを入れ直す（`reset_to_project()`）。捨てないと、案件 A の出力先で
  案件 B のダイアログが開く（DA 分析 2026-09-24 指摘）
- **最後のプロジェクトのフォルダはファイルに残す**（`enable_persistence()`）。起動し直しても
  プロジェクトを開く画面が前回のフォルダから始まる（DA 分析 2026-09-26 指摘）。
  プロジェクトの無いアプリは `kinds=ALL` で全用途を残す

使い方:

    from common import dialog_dirs as dd

    dd.enable_persistence("TubeAnalyzer")          # 起動時に 1 回
    path, _ = QFileDialog.getOpenFileName(self, "CSV", dd.base_dir(dd.DATA), "CSV (*.csv)")
    if path:
        dd.remember(dd.DATA, path)
    out, _ = QFileDialog.getSaveFileName(self, "保存", dd.suggest(dd.OUTPUT, "a.xlsx"))
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

# 用途。ダイアログごとに分けておくと、出力先を選んだあとにプロジェクトの保存が
# 出力フォルダで開く、といった迷子にならない。アプリが独自の用途名を足してもよい
PROJECT = "project"      # プロジェクトファイル（全体のフォールバックも兼ねる）
OUTPUT = "output"        # 書き出し先
DATA = "data"            # 測定データ（CSV・WAV・フォルダ）
IMPORT = "import"        # 読み込む結果・設定ファイル

ALL = None               # enable_persistence(kinds=ALL): 全用途をファイルに残す

# 確認用のスクリプト・スクリーンショット・テストでユーザーの記録を上書きしないための目印
NO_PERSIST_ENV = ("DIALOG_DIRS_NO_PERSIST", "DA_NO_PERSIST")

_dirs: dict[str, Path] = {}
_persist_file: Path | None = None
_persist_kinds: tuple[str, ...] | None = (PROJECT,)


def default_persist_file(app_name: str) -> Path:
    base = Path(os.environ.get("LOCALAPPDATA") or tempfile.gettempdir())
    return base / app_name / "last_dirs.json"


def enable_persistence(app_name: str, kinds=(PROJECT,), path=None) -> None:
    """覚えたフォルダを設定ファイルに残し、今ある分を読み込む。

    ``kinds`` は残す用途。既定はプロジェクトのフォルダだけ（案件をまたいで出力先が
    残ると迷子になるため）。プロジェクトの無いアプリは ``ALL``。
    """
    global _persist_file, _persist_kinds
    if any(os.environ.get(k) for k in NO_PERSIST_ENV):
        return
    _persist_file = Path(path) if path else default_persist_file(app_name)
    _persist_kinds = None if kinds is ALL else tuple(kinds)
    try:
        data = json.loads(_persist_file.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return
    if not isinstance(data, dict):
        return
    for kind, value in data.items():
        if _persist_kinds is not None and kind not in _persist_kinds:
            continue
        d = Path(value or "")
        if str(value or "") and d.is_dir():
            _dirs.setdefault(kind, d)


def disable_persistence() -> None:
    """ファイルに残すのをやめる（テスト用）。"""
    global _persist_file
    _persist_file = None


def _save_persist() -> None:
    if _persist_file is None:
        return
    keep = {k: str(v) for k, v in _dirs.items()
            if _persist_kinds is None or k in _persist_kinds}
    if not keep:
        return
    try:
        _persist_file.parent.mkdir(parents=True, exist_ok=True)
        _persist_file.write_text(json.dumps(keep, ensure_ascii=False), encoding="utf-8")
    except OSError:
        pass


def remember(kind: str, path) -> None:
    """選ばれたファイル（またはフォルダ）の場所を覚える。

    ファイルなら親フォルダ、フォルダならそのフォルダ。存在しないパスは覚えない。
    複数選んだときはリストを渡してよい（先頭を使う）。
    """
    if isinstance(path, (list, tuple)):
        path = path[0] if path else None
    if not path:
        return
    p = Path(path)
    d = p if p.is_dir() else p.parent
    if d.is_dir():
        _dirs[kind] = d
        if _persist_kinds is None or kind in _persist_kinds:
            _save_persist()


def base_dir(kind: str = PROJECT) -> str:
    """その用途の初期フォルダ。無ければプロジェクトの場所、それも無ければ空。

    `getOpenFileName()` / `getExistingDirectory()` の初期フォルダに渡す。
    """
    d = _dirs.get(kind) or _dirs.get(PROJECT)
    return str(d) if d else ""


def suggest(kind: str, name: str = "") -> str:
    """保存ダイアログに渡す初期パス（フォルダ ＋ ``name``）。覚えが無ければ ``name`` だけ。"""
    d = _dirs.get(kind) or _dirs.get(PROJECT)
    if d is None:
        return name
    return str(d / name) if name else str(d)


def forget() -> None:
    """全部忘れる（ファイルは消さない）。"""
    _dirs.clear()


def reset_to_project(project_path, seeds: dict | None = None) -> None:
    """**プロジェクトを開いた／保存したときに呼ぶ。**

    用途ごとの記憶を捨て、プロジェクトの場所を覚え直す。``seeds`` は
    ``{用途: [候補のパス, …]}`` で、用途ごとに**先に見つかった実在のフォルダ**を入れる
    （プロジェクトに書かれているデータ・出力先など）。
    """
    forget()
    remember(PROJECT, project_path)
    for kind, candidates in (seeds or {}).items():
        for cand in candidates or ():
            remember(kind, cand)
            if kind in _dirs:
                break
