"""計算条件を入力するウィンドウ（TODO G-7 の入口）。

**最終的には 1 つのウィンドウで完結する GUI にする**予定だが、
どんな情報が要るかがまだ固まっていないので、いまは
「必要なウィンドウをそのとき出す」形にしてある（`app.py` が順に開く）。

ここで決めるもの:

    モデル        DXF ファイル
    吸音率        CSV ファイルと、その値が垂直入射か残響室法か
    プロジェクト  保存先フォルダ（既存なら条件を読み込む）
    計算条件      音線数 / 最大反射回数 / 受音球の半径 / インパルス応答の長さ / バンド数
    大気条件      温度 / 湿度 / 気圧（**音速と空気吸収の両方が連動する**）
    法線          自動 / CAD のまま など

**tkinter（標準ライブラリ）で書いている。** 依存を増やさないため。
3D 表示が要るところ（法線の確認・可視化）は PyVista 側の別ウィンドウに任せる。
"""

import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import project as pj

BAND_CHOICES = [("8 バンド（63〜8k Hz）", 8), ("6 バンド（125〜4k Hz）", 6)]
# ★インパルス応答の合成のやり方（2026-08-23 ユーザー指摘）。
# 査読論文 Toyoda & Sakayoshi (2021) 式(2) の厳密解をそのまま解くか、
# 時間領域に置いてから FFT するか。**波動解と突き合わせるときは厳密版**
IMPULSE_CHOICES = [("高速（時間領域→FFT。指標を出すならこれで足りる）", "fast"),
                   ("厳密（式(2) をそのまま。波動解と比べるとき）", "exact")]
NORMAL_CHOICES = [("自動（閉じていれば内向き、開いていれば CAD のまま）", "auto"),
                  ("CAD の巻き順をそのまま使う", "cad"),
                  ("面ごとに室内側へ揃える", "inward"),
                  ("シェル単位で空気側へ揃える", "shells"),
                  ("全反転", "flip")]

# ★**区分ごとに分けてある**（2026-08-21 ユーザー要望）
GEOMETRY_FIELDS = [
    ("rays", "音線数", int, "多いほど後期の経路を拾える。20 万本が目安"),
    ("nref", "最大反射回数", int, "35 dB 減衰するまで必要。足りないと残響が短く出る"),
    ("radius", "受音球の半径 [m]", float,
     "★経路を見つけるための網。小さすぎると残響が短く出る"),
    ("max_time", "インパルス応答の長さ [s]", float, "残響時間より十分長くとる"),
]
ROOM_FIELDS = [
    ("volume", "室容積 [m³]", float, "空欄なら閉じた形状から自動算出"),
]
SOURCE_FIELDS = [
    # 音圧レベルの絶対値・STI の SNR に要る（2026-08-21 ユーザー要望）
    ("source_power_db", "音源 PWL [dB]", float,
     "空欄なら音圧レベルは相対値（Lw=0 dB）"),
    ("noise_level_db", "背景騒音 [dB]", float,
     "STI の SNR に使う。PWL と両方要る"),
]
# 大気は説明なしで 1 行に収める（ユーザー判断）
ATMOSPHERE_FIELDS = [
    ("temperature", "温度 [℃]", float, ""),
    ("humidity", "相対湿度 [%]", float, ""),
    ("pressure", "気圧 [kPa]", float, ""),
]
NUMBER_FIELDS = GEOMETRY_FIELDS + ROOM_FIELDS + SOURCE_FIELDS + ATMOSPHERE_FIELDS

# インパルス応答の長さを理論残響時間の何倍にするか（1 秒単位に切り上げる）。
# T30 は 35 dB 減るまで見るので 0.58 倍あれば測れるが、余裕を見て 1.5 倍
MAX_TIME_MARGIN = 1.5


def estimate_volume(dxf_path, unit=None):
    """DXF から室容積の目安を出す（統計残響式の入力を埋めるため）。

    閉じていれば符号付き体積、閉じていなければ**床面積 × 天井高**で見積もる。
    板の寄せ集めモデルでは容積が自動で決まらず、空欄のままだと統計残響式が
    黙って飛ばされる（実際にそれで比較できない結果が出た）ので、
    ここで目安を出して入力欄に入れられるようにした。

    **総表面積も一緒に返す**。統計残響式に効くのは容積だけでなく面積もで、
    「拾えている面が想定どおりか」を確かめる目安になるため（2026-08-19 ユーザー要望）。

    戻り値 (容積 [m³], 総表面積 [m²], 求め方の説明)。
    容積が求められなければ (None, 総表面積, 理由)。
    """
    import numpy as np
    import read_dxffile as rd

    # ★法線を揃えてから読む。CAD のままだと**天井の法線も上を向いている**ことがあり、
    #   「法線が真上＝床」で拾うと天井まで床に数えて容積が倍になる（実際になった）
    model = rd.read_model(dxf_path, unit=unit, orient_normals="auto", verbose=False)
    if not model.mesh:
        return None, 0.0, "面が読めません"
    area = float(model.surface_area)
    if model.volume:
        # `read_model` が法線から発散定理で出した値。辺が 1 対 1 で閉じていなくても
        # （T 字接合でも）正しく、家具・反射板の体積も引かれている
        return float(model.volume), area, f"囲まれた形状の空気容積（{model.volume_source}）"

    # 床（法線がほぼ真上＝室内を向いている水平面）の面積 × 高さ。
    # 壁が鉛直な部屋ならこれで足りる
    floor = 0.0
    for face in model.mesh:
        v = face.vertexes
        if face.normal[2] > 0.9:
            floor += 0.5 * float(np.linalg.norm(
                np.cross(v[1] - v[0], v[2] - v[0])))
    height = float(model.extents[1][2] - model.extents[0][2])
    if floor <= 0.0 or height <= 0.0:
        return None, area, "床が見つからないので見積もれません"
    return floor * height, area, f"床 {floor:.1f} m² × 高さ {height:.2f} m の目安"


class SetupWindow:
    """条件入力のウィンドウ。`run()` が (Project, 押されたボタン) を返す。

    ボタンは 'run'（計算する）/ 'run_all'（**全条件を一括**）/
    'normals'（面を確認する。法線と吸音材）/ 'view'（前回の結果を見る）/ None（閉じた）。
    """

    def __init__(self, project=None, folder=None):
        self.project = project or pj.Project.load(folder or os.getcwd())
        self.action = None
        self.root = None

    # ---- 組み立て ------------------------------------------------------

    def run(self):
        self.root = tk.Tk()
        self.root.title("幾何音響シミュレーション — 計算条件")
        # ★**ボタンが隠れない大きさで開く**（2026-08-24 ユーザー指摘。
        #   900 px だと「面を確認」のボタンが右へ押し出されて見えなかった）。
        #   `minsize` も入れて、**縮めても隠れない**ようにする
        self.root.geometry("1180x840")
        self.root.minsize(1040, 560)

        # ★**縦にスクロールできるようにする**（2026-08-21 ユーザー指摘。
        #   画面を広げないと下のボタンが見切れていた）。
        #   ボタンは下に固定し、条件の並びだけを送る
        buttons = ttk.Frame(self.root, padding=(14, 0, 14, 12))
        buttons.pack(side="bottom", fill="x")
        canvas = tk.Canvas(self.root, highlightthickness=0, borderwidth=0)
        bar = ttk.Scrollbar(self.root, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=bar.set)
        bar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        outer = ttk.Frame(canvas, padding=14)
        window = canvas.create_window((0, 0), window=outer, anchor="nw")

        def fit(_event=None):
            canvas.configure(scrollregion=canvas.bbox("all"))
            canvas.itemconfigure(window, width=canvas.winfo_width())

        outer.bind("<Configure>", fit)
        canvas.bind("<Configure>", fit)
        # マウスホイールで送る（Windows は <MouseWheel>）
        canvas.bind_all("<MouseWheel>",
                        lambda e: canvas.yview_scroll(-e.delta // 120, "units"))

        self._build_files(outer)
        self._build_numbers(outer)
        self._build_options(outer)
        self._build_buttons(buttons)
        self._load_into_widgets()

        self.root.mainloop()
        return self.project, self.action

    def _section(self, parent, title):
        frame = ttk.LabelFrame(parent, text=title, padding=10)
        frame.pack(fill="x", pady=(0, 12))
        frame.columnconfigure(1, weight=1)
        return frame

    def _build_files(self, parent):
        frame = self._section(parent, "ファイル")
        self.vars = {}

        rows = [
            ("folder", "プロジェクトフォルダ", "dir",
             "結果と図をここに保存します。既存フォルダを選ぶと前回の条件を読み込みます"),
            # ★対象室の名前。**空欄なら DXF のファイル名を使う**（2026-08-21 ユーザー要望。
            #   「DXF のファイル名は物件名や部屋名」）。結果ファイル名の頭になる
            ("name", "対象室名（任意）", "text",
             "結果ファイル名の頭に付きます。**空欄なら DXF のファイル名**を使います"),
            ("dxf", "モデル（DXF）", "file", "室形状。音源・受音点も src / rec レイヤから読みます"),
            # ★条件はファイルで選ぶ。用意していなければ「条件表を作成」で作る。
            #   xlsx なら**シート 1 枚が条件 1 つ**で、シート名が条件名になる
            ("condition_csv", "条件表（xlsx）", "file",
             "レイヤー名 → 材料番号の対応表。**シート 1 枚が条件 1 つ**です"
             "（用意が無ければ下の「条件表を作成」で作れます）"),
        ]
        for row, (key, label, kind, hint) in enumerate(rows):
            ttk.Label(frame, text=label).grid(row=row * 2, column=0, sticky="w", pady=2)
            var = tk.StringVar()
            self.vars[key] = var
            entry = ttk.Entry(frame, textvariable=var)
            entry.grid(row=row * 2, column=1, sticky="ew", padx=8)
            if kind != "text":      # 名前は手で書くので「参照…」は付けない
                ttk.Button(frame, text="参照…", width=8,
                           command=lambda k=key, t=kind: self._browse(k, t)
                           ).grid(row=row * 2, column=2)
            if key == "condition_csv":
                # ★条件（シート）を選ぶ。xlsx はシート 1 枚が条件 1 つ
                cell = ttk.Frame(frame)
                cell.grid(row=row * 2 + 1, column=1, columnspan=2, sticky="w",
                          padx=8)
                ttk.Label(cell, text="条件（シート）").pack(side="left")
                self.vars["condition_sheet"] = tk.StringVar()
                self.sheet_box = ttk.Combobox(
                    cell, textvariable=self.vars["condition_sheet"],
                    width=24, state="readonly", values=[])
                self.sheet_box.pack(side="left", padx=6)
                ttk.Button(cell, text="条件表を作成", width=14,
                           command=self._create_condition_table).pack(side="left",
                                                                      padx=4)
                ttk.Button(cell, text="開く", width=6,
                           command=self._open_condition_table).pack(side="left")
                continue        # 補足はこの行で使ったので下の Label は出さない
            ttk.Label(frame, text=hint, foreground="#666").grid(
                row=row * 2 + 1, column=1, columnspan=2, sticky="w", padx=8)

    def _number_row(self, frame, row, key, label, hint, width=14):
        """数値の入力欄 1 行（見出し・欄・補足）。戻り値は補足を置く枠。"""
        ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w", pady=2)
        var = tk.StringVar()
        self.vars[key] = var
        ttk.Entry(frame, textvariable=var, width=width).grid(row=row, column=1,
                                                            sticky="w", padx=8)
        cell = ttk.Frame(frame)
        cell.grid(row=row, column=2, sticky="w")
        if hint:
            ttk.Label(cell, text=hint, foreground="#666").pack(side="right")
        return cell

    def _build_numbers(self, parent):
        """計算条件。**区分ごとに枠を分ける**（2026-08-21 ユーザー要望）。

        幾何音響の設定（音線数・反射回数・受音球・応答の長さ）と、
        室の情報（容積）、音源の情報（PWL・騒音）、大気（1 行で済ませる）を分ける。
        """
        frame = self._section(parent, "計算条件（幾何音響）")
        for row, (key, label, _type, hint) in enumerate(GEOMETRY_FIELDS):
            cell = self._number_row(frame, row, key, label, hint)
            if key == "rays":
                # 音線がどう飛ぶかは**室形状と関係ない**ので、ここで先に見られる
                ttk.Button(cell, text="音線の飛び方を見る", width=18,
                           command=self._show_directions).pack(side="left",
                                                               padx=(0, 8))
            if key == "max_time":
                # ★容積・表面積・平均吸音率から理論上の残響時間を出して決める
                ttk.Button(cell, text="残響から推定", width=14,
                           command=self._estimate_max_time).pack(side="left",
                                                                 padx=(0, 8))

        frame = self._section(parent, "室（統計残響式に使う）")
        cell = self._number_row(frame, 0, "volume", "室容積 [m³]",
                                "空欄なら閉じた形状から自動算出")
        ttk.Button(cell, text="モデルから見積", width=16,
                   command=self._estimate_volume).pack(side="left", padx=(0, 8))
        self.room_note = ttk.Label(frame, text="", foreground="#3a6ea5",
                                   justify="left")
        self.room_note.grid(row=1, column=0, columnspan=3, sticky="w", pady=(4, 0))

        frame = self._section(parent, "音源（音圧レベル・STI に使う）")
        for row, (key, label, _type, hint) in enumerate(SOURCE_FIELDS):
            self._number_row(frame, row, key, label, hint)

        # ★大気は 1 行で済ませる（説明は要らないというユーザー判断 2026-08-21）
        frame = self._section(parent, "大気（音速と空気吸収が決まる）")
        line = ttk.Frame(frame)
        line.grid(row=0, column=0, sticky="w")
        for key, label, _type, _hint in ATMOSPHERE_FIELDS:
            ttk.Label(line, text=label).pack(side="left")
            var = tk.StringVar()
            self.vars[key] = var
            ttk.Entry(line, textvariable=var, width=8).pack(side="left",
                                                            padx=(4, 14))

    def _build_options(self, parent):
        frame = self._section(parent, "設定")

        def combo(row, key, label, choices):
            ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w", pady=3)
            var = tk.StringVar()
            self.vars[key] = var
            box = ttk.Combobox(frame, textvariable=var, state="readonly", width=52,
                               values=[text for text, _ in choices])
            box.grid(row=row, column=1, sticky="w", padx=8)
            return box

        # ★吸音率の種類は**条件表の「吸音率」シート**に材料ごとに書くので、
        #   ここには置かない（2026-08-21 ユーザー指摘）
        combo(0, "band_number", "周波数バンド", BAND_CHOICES)
        combo(1, "orient_normals", "法線の向き", NORMAL_CHOICES)
        combo(2, "impulse_method", "インパルス応答の合成", IMPULSE_CHOICES)
        self.vars["statistical"] = tk.BooleanVar()
        ttk.Checkbutton(frame, text="統計残響式（Sabine / Eyring / Eyring-Knudsen）も計算する",
                        variable=self.vars["statistical"]).grid(row=3, column=1,
                                                                sticky="w", padx=8)

    def _build_buttons(self, parent):
        """下端のボタン。**状態の文とボタンで行を分ける**（2026-08-24 ユーザー指摘）。

        以前は 1 行に「状態の文（左）＋ボタン（右）」を詰めていたので、
        窓が狭いと**ボタンが右へ押し出されて見えなかった**
        （「面を確認」が消えていた）。行を分け、ボタンは `grid` で並べて
        **足りなければ 2 段目に折り返す**ようにした。
        """
        # ① 状態の文（左のログっぽいもの）は独立した行
        line = ttk.Frame(parent)
        line.pack(fill="x", pady=(4, 2))
        self.status = ttk.Label(line, text="", foreground="#333", anchor="w",
                                justify="left")
        self.status.pack(side="left", fill="x", expand=True)

        # ② ボタンの行。**右端から並べる**が、幅が足りなければ折り返す
        frame = ttk.Frame(parent)
        frame.pack(fill="x", pady=(2, 0))
        self._button_row = frame

        # （左寄りの補助的なもの → 右寄りの本命）の順に並べる
        items = [
            ("条件だけ保存", self._on_save),
            ("面を確認…（法線・吸音材）", self._on_normals),
            ("前回の結果を見る", self._on_view),
            # ★同じフォルダの条件表を全部回す（2026-08-21 ユーザー要望）。
            #   経路は吸音に依らないので 2 件目以降は一瞬で終わる（F-9）
            ("全条件を一括 ▶▶", self._on_run_all),
            ("計算する ▶", self._on_run),
            ("閉じる", self._close),
        ]
        self._buttons = []
        for text, command in items:
            button = ttk.Button(frame, text=text, command=command)
            self._buttons.append(button)
            if text == "前回の結果を見る":
                self.view_button = button
        self._layout_buttons()
        frame.bind("<Configure>", lambda _e: self._layout_buttons())

    def _layout_buttons(self):
        """ボタンを幅に合わせて並べ直す。**入りきらなければ 2 段にする。**

        ★隠さないのが約束（ユーザー指摘）。`pack(side="right")` は入りきらない
        ぶんを黙って画面外へ出してしまうので、`grid` で段を作る。
        """
        frame = getattr(self, "_button_row", None)
        if frame is None:
            return
        width = frame.winfo_width() or 1
        needed = [b.winfo_reqwidth() + 8 for b in self._buttons]
        # 1 行に何個入るか（最低 1 個）
        per_row, total = 0, 0
        for w in needed:
            if total + w > width and per_row:
                break
            total += w
            per_row += 1
        per_row = max(1, per_row)
        if per_row == getattr(self, "_button_per_row", None):
            return                          # 並びが変わらないなら触らない
        self._button_per_row = per_row
        for b in self._buttons:
            b.grid_forget()
        for index, button in enumerate(self._buttons):
            row, column = divmod(index, per_row)
            button.grid(row=row, column=column, padx=4, pady=2, sticky="ew")
        for column in range(per_row):
            frame.columnconfigure(column, weight=1)

    # ---- 値のやりとり --------------------------------------------------

    def _load_into_widgets(self):
        p = self.project
        self.vars["folder"].set(p.folder)
        self.vars["name"].set(p.name or "")
        self.vars["dxf"].set(p.dxf_path or "")
        self.vars["condition_csv"].set(p.condition_csv and
                                       (p.resolve(p.condition_csv) or "") or "")
        self._refresh_sheets(p.condition_sheet or "")
        for key, _label, _type, _hint in NUMBER_FIELDS:
            value = getattr(p, key)
            self.vars[key].set("" if value is None else str(value))
        self._set_combo("band_number", BAND_CHOICES, p.band_number)
        self._set_combo("orient_normals", NORMAL_CHOICES, p.orient_normals)
        self._set_combo("impulse_method", IMPULSE_CHOICES,
                        getattr(p, "impulse_method", "fast"))
        self.vars["statistical"].set(bool(p.statistical))
        self._show_status()

    def _set_combo(self, key, choices, value):
        for text, item in choices:
            if item == value:
                self.vars[key].set(text)
                return
        self.vars[key].set(choices[0][0])

    def _combo_value(self, key, choices):
        text = self.vars[key].get()
        for label, value in choices:
            if label == text:
                return value
        return choices[0][1]

    def _collect(self):
        """ウィジェットの値を Project に取り込む。問題があれば文字列で返す。"""
        folder = self.vars["folder"].get().strip()
        if not folder:
            return "プロジェクトフォルダを指定してください"
        dxf = self.vars["dxf"].get().strip()
        if not dxf or not os.path.exists(dxf):
            return "DXF ファイルが見つかりません"

        if os.path.abspath(folder) != self.project.folder:
            # フォルダを変えたら、そのフォルダの既存条件を土台にする
            self.project = pj.Project.load(folder)

        self.project.dxf = dxf
        # 対象室名。**空欄のままにする**（空なら DXF のファイル名が使われる）
        self.project.name = self.vars["name"].get().strip()
        # 条件表と条件シート。**シート名が条件名**になり、結果ファイル名にも入る
        self.project.condition_csv = self.vars["condition_csv"].get().strip()
        self.project.condition_sheet = self.vars["condition_sheet"].get().strip()
        for key, _label, cast, _hint in NUMBER_FIELDS:
            text = self.vars[key].get().strip()
            if not text:
                # 空欄を「未入力（None）」として扱う条件。既定値に戻してはいけない
                # （PWL は None のとき相対値、0 dB のときは絶対値 0 dB を意味する）
                setattr(self.project, key,
                        None if pj.DEFAULTS[key] is None else pj.DEFAULTS[key])
                continue
            try:
                setattr(self.project, key, cast(float(text)) if cast is int
                        else cast(text))
            except ValueError:
                return f"{_label} の値が数値になっていません: {text!r}"

        self.project.band_number = self._combo_value("band_number", BAND_CHOICES)
        self.project.orient_normals = self._combo_value("orient_normals",
                                                        NORMAL_CHOICES)
        self.project.impulse_method = self._combo_value("impulse_method",
                                                        IMPULSE_CHOICES)
        self.project.statistical = bool(self.vars["statistical"].get())
        return None

    def _show_status(self):
        parts = []
        if pj.Project.exists(self.project.folder):
            parts.append("既存プロジェクトを読み込みました")
        if os.path.exists(self.project.path(pj.NORMALS_FILE)):
            flipped, _ = self.project.load_flipped_faces()
            parts.append(f"法線の指定あり（{len(flipped)} 枚反転）")
        if pj.has_results(self.project):
            parts.append("計算結果あり")
        self.status.config(text="   ".join(parts))

    # ---- ボタン --------------------------------------------------------

    def _browse(self, key, kind):
        current = self.vars[key].get().strip()
        start = current if os.path.isdir(current) else os.path.dirname(current)
        if kind == "dir":
            path = filedialog.askdirectory(title="プロジェクトフォルダ",
                                           initialdir=start or os.getcwd())
            if path:
                self.vars[key].set(os.path.normpath(path))
                # そのフォルダに project.json があれば読み込んで反映する
                if pj.Project.exists(path):
                    self.project = pj.Project.load(path)
                    self._load_into_widgets()
            return
        if key == "dxf":
            patterns = [("DXF", "*.dxf"), ("すべて", "*.*")]
        elif key == "condition_csv":
            # 条件表は xlsx が本命。昔の CSV も選べるようにしておく
            patterns = [("条件表", "*.xlsx *.xlsm *.csv"), ("すべて", "*.*")]
        else:
            patterns = [("CSV", "*.csv"), ("すべて", "*.*")]
        path = filedialog.askopenfilename(title="ファイルを選ぶ", filetypes=patterns,
                                          initialdir=start or os.getcwd())
        if path:
            self.vars[key].set(os.path.normpath(path))
            if key == "condition_csv":
                self._refresh_sheets()      # 選んだ表の条件シートを並べ直す

    def _show_directions(self):
        """音線がどの向きへ飛ぶかを見る（室形状は関係ないので単体で開ける）。"""
        text = self.vars["rays"].get().strip()
        try:
            total = int(float(text))
        except ValueError:
            messagebox.showerror("音線数を確認してください",
                                 f"数値になっていません: {text!r}")
            return
        if total < 1:
            messagebox.showerror("音線数を確認してください", "1 以上にしてください")
            return
        import view_directions

        # 保存先はプロジェクトフォルダが決まっているときだけ（`図/画面/`）
        folder = self.vars["folder"].get().strip()
        save_dir = None
        if folder and os.path.isdir(folder):
            import project as pj
            save_dir = os.path.join(folder, pj.SCREENSHOT_DIR)
        try:
            view_directions.show(total=total, save_dir=save_dir)
        except Exception as e:
            messagebox.showerror("表示できませんでした", f"{type(e).__name__}: {e}")

    def _room_summary(self):
        """いまの条件で**容積・総表面積・平均吸音率・理論残響時間**を出す。

        ★「モデルから見積」と「残響から推定」で共通に使う（2026-08-21 ユーザー要望）。
        平均吸音率は**上で選んだ条件（シート）に則った**値にする
        （条件表の材料番号 → 吸音率シート → 安全率まで通した結果）。

        戻り値 dict … 'volume' / 'area' / 'note' / 'alpha'（(nf,) or None）/
                      'frequencies' / 'reverberation'（統計残響式の結果 or None）
        """
        problem = self._collect()
        if problem:
            raise ValueError(problem)
        self.project.save()

        import reverberation as rv
        import run_project

        model = run_project._model_for(self.project)
        area = float(model.surface_area)
        volume = self.project.volume or (abs(model.volume) if model.volume else None)
        note = ("指定した室容積" if self.project.volume
                else (f"囲まれた形状から（{model.volume_source}）"
                      if model.volume else "容積が決まりません"))
        if volume is None:
            # 閉じていないモデルは床面積 × 高さで見積もる（従来と同じ手）
            volume, area_guess, note = estimate_volume(self.project.dxf_path,
                                                       unit=self.project.unit)
            area = area or area_guess
        statistical = None
        if volume:
            statistical = rv.statistical_reverberation(
                model.mesh, abs(volume), atmosphere=self._atmosphere(),
                verbose=False)
        return {"volume": volume, "area": area, "note": note,
                "model": model, "statistical": statistical}

    def _atmosphere(self):
        from atmosphere import Atmosphere
        return Atmosphere(temperature=self.project.temperature,
                          humidity=self.project.humidity,
                          pressure=self.project.pressure)

    def _estimate_volume(self):
        """モデルから容積を見積もり、**表面積と平均吸音率も併せて出す**。

        ユーザー要望（2026-08-21）：入力するのは容積だけで、
        表面積と平均吸音率は**確認のための表示**（入力欄は作らない）。
        """
        try:
            found = self._room_summary()
        except ValueError as error:
            messagebox.showwarning("条件を確認してください", str(error))
            return
        except Exception as e:
            messagebox.showerror("見積もれませんでした", f"{type(e).__name__}: {e}")
            return

        volume, area = found["volume"], found["area"]
        if volume:
            self.vars["volume"].set(f"{volume:.2f}")
        lines = [f"総表面積 {area:.1f} m²",
                 (f"室容積 {volume:.2f} m³（{found['note']}）" if volume
                  else f"室容積は決まりません（{found['note']}）")]
        statistical = found["statistical"]
        if statistical is not None:
            alpha = statistical["mean_absorption"]
            frequencies = statistical["frequencies"]
            lines.append("平均吸音率 ᾱ（条件『"
                         + (self.project.condition_label or "既定") + "』）")
            lines.append("  " + " / ".join(f"{f:.0f}Hz {a:.2f}"
                                           for f, a in zip(frequencies, alpha)))
            lines.append(f"  平均自由行程 4V/S = "
                         f"{4.0 * volume / area:.2f} m" if area else "")
        else:
            lines.append("平均吸音率は出せません（容積か材料が決まっていません）")
        self.room_note.configure(text=chr(10).join(l for l in lines if l))
        self._show_status()

    def _estimate_max_time(self):
        """理論上の残響時間から**インパルス応答の長さ**を決める（1 秒刻み）。

        ★容積・表面積・平均吸音率がそろえば残響時間は式で出る、というユーザー指摘。
        いちばん長いバンドの理論値（Eyring-Knudsen）の **1.5 倍**を取って
        1 秒単位に切り上げる。T30 は 35 dB 減るまで見るので `0.58 × T60` あれば
        足りるが、初期の遅れと余裕を見て 1.5 倍にしてある。
        入力欄はそのまま残すので、任意の値も打ち込める。
        """
        try:
            found = self._room_summary()
        except ValueError as error:
            messagebox.showwarning("条件を確認してください", str(error))
            return
        except Exception as e:
            messagebox.showerror("推定できませんでした", f"{type(e).__name__}: {e}")
            return

        statistical = found["statistical"]
        if statistical is None:
            messagebox.showwarning(
                "推定できませんでした",
                "室容積か材料が決まっていないので理論残響時間が出せません。"
                "「モデルから見積」で容積を入れてから試してください。")
            return
        import numpy as np

        times = np.asarray(statistical["eyring_knudsen"], dtype=float)
        if not np.any(np.isfinite(times)):
            times = np.asarray(statistical["sabine"], dtype=float)
        longest = float(np.nanmax(times))
        want = max(1.0, float(np.ceil(longest * MAX_TIME_MARGIN)))
        self.vars["max_time"].set(f"{want:.0f}")
        band = statistical["frequencies"][int(np.nanargmax(times))]
        messagebox.showinfo(
            "インパルス応答の長さ",
            f"{want:.0f} 秒を入れました。{chr(10)}{chr(10)}"
            f"理論上いちばん長いのは {band:.0f} Hz の {longest:.2f} 秒"
            f"（Eyring-Knudsen）。{chr(10)}"
            f"その {MAX_TIME_MARGIN:g} 倍を 1 秒単位に切り上げています"
            f"（T30 は 35 dB 減るまで見るので、"
            f"理論値の 0.58 倍あれば測れます）。{chr(10)}"
            f"任意の値を打ち込んでも構いません。")

    def _on_save(self):
        error = self._collect()
        if error:
            messagebox.showerror("入力を確認してください", error)
            return
        path = self.project.save()
        self._show_status()
        messagebox.showinfo("保存しました", f"条件を保存しました:\n{path}")

    def _on_normals(self):
        error = self._collect()
        if error:
            messagebox.showerror("入力を確認してください", error)
            return
        self.project.save()
        self.action = "normals"
        self.root.destroy()

    def _on_run_all(self):
        """同じフォルダの材料条件表を**全部**回す（一括計算）。"""
        import condition_table as ct

        problem = self._collect()
        if problem:
            messagebox.showwarning("条件を確認してください", problem)
            return
        found = ct.discover(self.project.folder)
        if not found:
            messagebox.showinfo(
                "条件表がありません",
                f"{self.project.folder} に材料条件表が見つかりません。"
                f"「材料条件表を開く」で作ってから、条件ごとに名前を付けて"
                f"保存してください。")
            return
        names = chr(10).join("  ・" + os.path.basename(f) for f in found)
        if not messagebox.askokcancel(
                "全条件を一括計算",
                f"{len(found)} 条件を続けて計算します。{chr(10)}{chr(10)}{names}{chr(10)}{chr(10)}"
                f"音線追跡は 1 回だけで、2 件目以降は吸音率を当て直すだけなので"
                f"すぐ終わります。"):
            return
        self.project.save()
        self.action = "run_all"
        self.root.destroy()

    def _refresh_sheets(self, chosen=None):
        """条件表の中の条件シートを読んで、コンボボックスに並べる。"""
        import condition_table as ct

        file_name = self.vars["condition_csv"].get().strip()
        if not file_name:
            file_name = self.project.condition_path
        names = [n for n in ct.sheets(file_name) if n]
        # ★先頭に空欄を置く。**選ばなければ結果ファイル名に条件名を付けない**
        #   （条件表を作った瞬間に名前が変わって、前の結果と揃わなくなるのを防ぐ。
        #     読むときは最初のシートの割り当てを使う）
        self.sheet_box.configure(values=[""] + names)
        current = (chosen if chosen is not None
                   else self.vars["condition_sheet"].get()).strip()
        self.vars["condition_sheet"].set(current if current in names else "")

    def _create_condition_table(self):
        """**DXF のレイヤから条件表（xlsx）を作る**（ユーザー要望 2026-08-21）。

        「吸音率」シートに材料一覧、条件シートにレイヤーと材料番号の欄が入る。
        既にあるファイルは上書きせず、面数・面積だけ更新する。
        """
        import condition_table as ct
        import run_project

        problem = self._collect()
        if problem:
            messagebox.showwarning("条件を確認してください", problem)
            return
        try:
            self.project.save()
            model = run_project._model_for(self.project)
            library = run_project._library_for(self.project)
            if library is None:
                messagebox.showwarning(
                    "吸音率の一覧がありません",
                    "吸音率（CSV）を指定すると、その材料一覧を「吸音率」シートに"
                    "書き込みます。指定しないまま作ると番号を手で書くことになります。")
            path = ct.create(self.project, model, library)
        except Exception as e:
            messagebox.showerror("条件表を作れませんでした",
                                 f"{type(e).__name__}: {e}")
            return
        self.vars["condition_csv"].set(path)
        self._refresh_sheets()
        layers = len(getattr(model, "layer_counts", {}) or {})
        messagebox.showinfo(
            "条件表を作りました",
            f"{path}{chr(10)}{chr(10)}レイヤー {layers} 件を並べました。"
            f"Excel で「材料番号」の列に番号を入れてください"
            f"（材料名は隣に自動で出ます）。"
            f"{chr(10)}条件を増やすときはシートを複製して名前を変えます。")
        self._open_condition_table()

    def _open_condition_table(self):
        """条件表を既定のアプリ（Excel など）で開く。無ければ作る。"""
        import condition_table as ct

        problem = self._collect()
        if problem:
            messagebox.showwarning("条件を確認してください", problem)
            return
        path = ct.path(self.project)
        if not os.path.exists(path):
            self._create_condition_table()
            return
        try:
            os.startfile(path)
        except Exception:
            messagebox.showinfo("条件表", f"ここにあります: {path}")

    def _on_view(self):
        """計算し直さずに、保存済みの結果を開く。"""
        if not pj.has_results(self.project):
            messagebox.showinfo("結果がありません",
                                f"{self.project.folder} にまだ計算結果がありません。\n"
                                f"「計算する ▶」を先に実行してください。")
            return
        self.action = "view"
        self.root.destroy()

    def _on_run(self):
        error = self._collect()
        if error:
            messagebox.showerror("入力を確認してください", error)
            return
        self.project.save()
        self.action = "run"
        self.root.destroy()

    def _close(self):
        self.action = None
        self.root.destroy()


def ask(project=None, folder=None):
    """条件入力ウィンドウを開く。(Project, 'run'|'normals'|None) を返す。

    'normals' は**面の確認ウィンドウ**（`face_editor.py`。法線と吸音材）。
    キーの名前は project.json 側との互換のためそのままにしてある。
    """
    return SetupWindow(project=project, folder=folder).run()


if __name__ == "__main__":
    import sys

    result, action = ask(folder=sys.argv[1] if len(sys.argv) > 1 else None)
    print(f"action={action}\n{result.summary()}")
