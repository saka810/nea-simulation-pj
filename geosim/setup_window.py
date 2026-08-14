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
KIND_CHOICES = [("残響室法（乱入射）", "random"),
                ("垂直入射", "normal"),
                ("CSV の # kind: 宣言に従う", "")]
NORMAL_CHOICES = [("自動（閉じていれば内向き、開いていれば CAD のまま）", "auto"),
                  ("CAD の巻き順をそのまま使う", "cad"),
                  ("面ごとに室内側へ揃える", "inward"),
                  ("シェル単位で空気側へ揃える", "shells"),
                  ("全反転", "flip")]

# (属性名, ラベル, 型, 補足)
NUMBER_FIELDS = [
    ("rays", "音線数", int, "多いほど後期の経路を拾える。20 万本が目安"),
    ("nref", "最大反射回数", int, "35 dB 減衰するまで必要。足りないと残響が短く出る"),
    ("radius", "受音球の半径 [m]", float,
     "★経路を見つけるための網。小さすぎると残響が短く出る（半径を変えて値が動かないか確認）"),
    ("max_time", "インパルス応答の長さ [s]", float, "残響時間より十分長くとる"),
    ("temperature", "温度 [℃]", float, "音速と空気吸収が変わる"),
    ("humidity", "相対湿度 [%]", float, "高域の空気吸収に効く"),
    ("pressure", "気圧 [kPa]", float, "既定 101.325"),
    ("volume", "室容積 [m³]", float,
     "統計残響式に使う。空欄なら閉じた形状から自動算出（開いた形状では要指定）"),
]


class SetupWindow:
    """条件入力のウィンドウ。`run()` が (Project, 押されたボタン) を返す。

    ボタンは 'run'（計算する）/ 'normals'（法線を確認する）/ None（閉じた）。
    """

    def __init__(self, project=None, folder=None):
        self.project = project or pj.Project.load(folder or os.getcwd())
        self.action = None
        self.root = None

    # ---- 組み立て ------------------------------------------------------

    def run(self):
        self.root = tk.Tk()
        self.root.title("幾何音響シミュレーション — 計算条件")
        self.root.geometry("880x760")

        outer = ttk.Frame(self.root, padding=14)
        outer.pack(fill="both", expand=True)

        self._build_files(outer)
        self._build_numbers(outer)
        self._build_options(outer)
        self._build_buttons(outer)
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
            ("dxf", "モデル（DXF）", "file", "室形状。音源・受音点も src / rec レイヤから読みます"),
            ("absorption_csv", "吸音率（CSV）", "file",
             "1列目=材料名または ID。レイヤ名の先頭の数字でも引けます"),
        ]
        for row, (key, label, kind, hint) in enumerate(rows):
            ttk.Label(frame, text=label).grid(row=row * 2, column=0, sticky="w", pady=2)
            var = tk.StringVar()
            self.vars[key] = var
            entry = ttk.Entry(frame, textvariable=var)
            entry.grid(row=row * 2, column=1, sticky="ew", padx=8)
            ttk.Button(frame, text="参照…", width=8,
                       command=lambda k=key, t=kind: self._browse(k, t)
                       ).grid(row=row * 2, column=2)
            ttk.Label(frame, text=hint, foreground="#666").grid(
                row=row * 2 + 1, column=1, columnspan=2, sticky="w", padx=8)

    def _build_numbers(self, parent):
        frame = self._section(parent, "計算条件")
        for row, (key, label, _type, hint) in enumerate(NUMBER_FIELDS):
            ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w", pady=2)
            var = tk.StringVar()
            self.vars[key] = var
            ttk.Entry(frame, textvariable=var, width=14).grid(row=row, column=1,
                                                              sticky="w", padx=8)
            ttk.Label(frame, text=hint, foreground="#666").grid(row=row, column=2,
                                                               sticky="w")

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

        combo(0, "absorption_kind", "吸音率の種類", KIND_CHOICES)
        ttk.Label(frame, text="★取り違えると吸音を大きく誤ります", foreground="#a33"
                  ).grid(row=0, column=2, sticky="w")
        combo(1, "band_number", "周波数バンド", BAND_CHOICES)
        combo(2, "orient_normals", "法線の向き", NORMAL_CHOICES)

        self.vars["two_sided"] = tk.BooleanVar()
        ttk.Checkbutton(frame, text="面の裏からの入射も当てる（法線がまちまちなモデルの代替手段）",
                        variable=self.vars["two_sided"]).grid(row=3, column=1,
                                                              sticky="w", padx=8, pady=3)
        self.vars["statistical"] = tk.BooleanVar()
        ttk.Checkbutton(frame, text="統計残響式（Sabine / Eyring / Millington）も計算する",
                        variable=self.vars["statistical"]).grid(row=4, column=1,
                                                                sticky="w", padx=8)

    def _build_buttons(self, parent):
        frame = ttk.Frame(parent)
        frame.pack(fill="x", pady=(6, 0))

        self.status = ttk.Label(frame, text="", foreground="#333")
        self.status.pack(side="left")

        ttk.Button(frame, text="閉じる", command=self._close).pack(side="right", padx=4)
        ttk.Button(frame, text="計算する ▶", command=self._on_run).pack(side="right", padx=4)
        ttk.Button(frame, text="法線を確認…", command=self._on_normals).pack(side="right",
                                                                           padx=4)
        ttk.Button(frame, text="条件だけ保存", command=self._on_save).pack(side="right",
                                                                        padx=4)

    # ---- 値のやりとり --------------------------------------------------

    def _load_into_widgets(self):
        p = self.project
        self.vars["folder"].set(p.folder)
        self.vars["dxf"].set(p.dxf_path or "")
        self.vars["absorption_csv"].set(p.absorption_path or "")
        for key, _label, _type, _hint in NUMBER_FIELDS:
            value = getattr(p, key)
            self.vars[key].set("" if value is None else str(value))
        self._set_combo("absorption_kind", KIND_CHOICES, p.absorption_kind or "")
        self._set_combo("band_number", BAND_CHOICES, p.band_number)
        self._set_combo("orient_normals", NORMAL_CHOICES, p.orient_normals)
        self.vars["two_sided"].set(bool(p.two_sided))
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
        absorption = self.vars["absorption_csv"].get().strip()
        if absorption and not os.path.exists(absorption):
            return "吸音率 CSV が見つかりません"

        if os.path.abspath(folder) != self.project.folder:
            # フォルダを変えたら、そのフォルダの既存条件を土台にする
            self.project = pj.Project.load(folder)

        self.project.dxf = dxf
        self.project.absorption_csv = absorption
        for key, _label, cast, _hint in NUMBER_FIELDS:
            text = self.vars[key].get().strip()
            if not text:
                setattr(self.project, key, None if key == "volume" else
                        pj.DEFAULTS[key])
                continue
            try:
                setattr(self.project, key, cast(float(text)) if cast is int
                        else cast(text))
            except ValueError:
                return f"{_label} の値が数値になっていません: {text!r}"

        self.project.absorption_kind = self._combo_value("absorption_kind",
                                                         KIND_CHOICES) or None
        self.project.band_number = self._combo_value("band_number", BAND_CHOICES)
        self.project.orient_normals = self._combo_value("orient_normals",
                                                        NORMAL_CHOICES)
        self.project.two_sided = bool(self.vars["two_sided"].get())
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
        patterns = ([("DXF", "*.dxf"), ("すべて", "*.*")] if key == "dxf"
                    else [("CSV", "*.csv"), ("すべて", "*.*")])
        path = filedialog.askopenfilename(title="ファイルを選ぶ", filetypes=patterns,
                                          initialdir=start or os.getcwd())
        if path:
            self.vars[key].set(os.path.normpath(path))

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
    """条件入力ウィンドウを開く。(Project, 'run'|'normals'|None) を返す。"""
    return SetupWindow(project=project, folder=folder).run()


if __name__ == "__main__":
    import sys

    result, action = ask(folder=sys.argv[1] if len(sys.argv) > 1 else None)
    print(f"action={action}\n{result.summary()}")
