"""プロジェクト（計算条件・法線の修正・計算結果）の保存と読み込み。

**ねらい**：一度作った条件と結果を、あとから開き直せるようにする。
条件を JSON、結果を CSV、図を PNG で**プロジェクトフォルダに全部置く**。
CSV にしてあるのは Excel でそのまま開けるようにするため。

    プロジェクトフォルダ/
      project.json          条件（DXF・吸音率・音線数・受音球・温度湿度…）
      normals.json          法線の反転指定（面ごと。face_editor.py が書く）
      materials.json        面ごとの吸音材の割り当て（同上。レイヤで分けられないモデル用）
      結果/
        まとめ_残響時間.csv   全受音点 ＋ 平均 ＋ 理論値（summary.py が書く）
        まとめ_明瞭度.csv     全受音点 ＋ 平均
        rt_statistical.csv  統計残響式 Sabine / Eyring / Eyring-Knudsen
        surface.csv         材料別の面積と吸音率（面ごとの割り当てが無ければレイヤ別）
        raylog.npz          音線軌跡（可視化用。可変長なので npz）
        rec1/               ← **受音点ごと**
          pulses.csv        パルス列（反射回数・到来時刻・到来方向・バンド別エネルギー）
          ir.csv            インパルス応答
          rt.csv            残響指標 EDT / T20 / T30
          decay.csv         減衰曲線
          clarity.csv       明瞭度 C50 / C80 / D50 / Ts
        rec2/ …
      図/
        rec1/ *.png         正規化したインパルス応答・減衰曲線・残響時間 ほか
        rec2/ …
        画面/               画面から手で撮った画像・動画

★**受音点ごとのものは `結果/recN/`、受音点に依らないものは `結果/` 直下**
（2026-08-21 にこの形へ。それまでは 1 点目だけ `結果/` 直下、2 点目以降が
`rec2/結果/` という不揃いな置き方で、ユーザー指摘で直した）。
統計残響式・材料別面積・音線軌跡は受音点に依らないので 1 つだけ持つ。

DXF や吸音率 CSV は**プロジェクトフォルダからの相対パスで持つ**（フォルダごと
別の端末へ移してもそのまま開ける）。フォルダの外にある場合は絶対パスのままにする。
"""

import json
import os

import numpy as np

PROJECT_FILE = "project.json"
NORMALS_FILE = "normals.json"
MATERIALS_FILE = "materials.json"
RESULT_DIR = "結果"
FIGURE_DIR = "図"
# 受音点ごとのフォルダ名（`結果/rec1/` `図/rec1/`）
RECEIVER_DIR = "rec%d"
# 画面から手で撮った画像・動画の置き場。**`図/` の直下ではなく子フォルダにする。**
# `clear_results()` が `図/` の PNG を消してしまうので、
# 同じ所に置くと計算し直すたびに撮った画像が巻き添えで消える
SCREENSHOT_DIR = os.path.join(FIGURE_DIR, "画面")

# 結果ファイルの名前。**キーはコード側の呼び名**で、値が実ファイル名
RESULT_FILES = {
    "pulses": "pulses.csv",
    "ir": "ir.csv",
    "rt": "rt.csv",
    "statistical": "rt_statistical.csv",
    "decay": "decay.csv",
    "surface": "surface.csv",
    "raylog": "raylog.npz",
}

# **受音点に依らない**結果。受音点ごとのフォルダではなく `結果/` 直下に置く。
#   統計残響式・材料別面積 … 室形状と材料だけで決まる
#   音線軌跡             … 音源から出た音線の形。受音点をまたいで共有している（F-6）
SHARED_RESULTS = {"statistical", "surface", "raylog"}

# project.json に書き出す条件と既定値。
# ここに無いキーは保存されないので、**新しい計算条件を足したらここにも足すこと**
DEFAULTS = {
    "name": "",
    "dxf": "",
    "absorption_csv": "",
    "absorption_kind": None,       # 'normal' | 'random' | None（CSV の # kind: を見る）
    "assignment": None,            # レイヤ → 材料の対応（辞書）
    "band_number": 6,
    "unit": None,                  # None なら DXF の $INSUNITS
    "orient_normals": "auto",
    "two_sided": False,
    "rays": 200000,
    "nref": 120,
    "radius": 1.0,
    "max_time": 3.0,
    "temperature": 20.0,
    "humidity": 40.0,
    "pressure": 101.325,
    "volume": None,                # 統計残響式に使う容積。None なら閉形状から自動
    "source": None,                # None なら DXF の src レイヤ
    "receiver": None,              # None なら DXF の rec レイヤ
    # 受音点に置く「人」の正面方向。真上から見た方位角 [度]。
    # **0° = +X 方向、反時計回り**（真上から見て）。
    # 伝搬方向の図（G-5）で「前・後ろ・左・右」を決めるのに使う。
    # 上下の向きは扱わない（実務では水平面で足りるため。ユーザー判断）
    #
    # ★**受音点ごとに持てる**（2026-08-20 ユーザー要望）。値は数値でもリストでもよい。
    #   数値なら全受音点に同じ向きを使う（従来の project.json をそのまま読める）。
    #   リストなら k 番目の受音点に k 番目の向きを使う。`head_azimuth_for()` を通すこと
    "head_azimuth": 0.0,
    "raylog_max_rays": 2000,
    "statistical": True,
}


class Project:
    """計算条件の入れ物。`folder` に紐づく。

    属性は `DEFAULTS` のキーがそのまま生える（`project.rays` のように読める）。
    """

    def __init__(self, folder, **values):
        self.folder = os.path.abspath(folder)
        for key, default in DEFAULTS.items():
            setattr(self, key, values.get(key, default))
        if not self.name:
            self.name = os.path.basename(self.folder)
        # いま何番目の受音点を扱っているか（1 始まり）。**保存する条件ではない**ので
        # DEFAULTS には入れない。`結果/recN/` `図/recN/` の振り分けにだけ使う
        self.receiver_index = values.get("receiver_index")

    # ---- パス ---------------------------------------------------------

    def path(self, *parts):
        return os.path.join(self.folder, *parts)

    def result_dir(self, shared=False):
        """結果の置き場。受音点が決まっていれば `結果/recN/`。

        `shared=True` は受音点に依らないもの（統計残響式など）で、
        受音点を扱っていても `結果/` 直下を返す。
        """
        if shared or self.receiver_index is None:
            return self.path(RESULT_DIR)
        return self.path(RESULT_DIR, RECEIVER_DIR % self.receiver_index)

    def result_path(self, key):
        """結果ファイルのパス。`key` は RESULT_FILES のキー。

        受音点に依らないもの（`SHARED_RESULTS`）は `結果/` 直下、
        それ以外は `結果/recN/` に置く。
        """
        return os.path.join(self.result_dir(shared=key in SHARED_RESULTS),
                            RESULT_FILES[key])

    def clarity_path(self):
        """明瞭度の CSV。`RESULT_FILES` に入れていないので別に持つ。"""
        return os.path.join(self.result_dir(), "clarity.csv")

    def figure_dir(self):
        """図の置き場。受音点が決まっていれば `図/recN/`。"""
        if self.receiver_index is None:
            return self.path(FIGURE_DIR)
        return self.path(FIGURE_DIR, RECEIVER_DIR % self.receiver_index)

    def figure_path(self, name):
        return os.path.join(self.figure_dir(), name)

    def screenshot_dir(self):
        """画面から手で撮った画像・動画の置き場（`図/画面/`）。"""
        folder = self.path(SCREENSHOT_DIR)
        os.makedirs(folder, exist_ok=True)
        return folder

    def resolve(self, value):
        """project.json に入っている相対パスを実際のパスに直す。"""
        if not value:
            return None
        return value if os.path.isabs(value) else self.path(value)

    def _relative(self, value):
        """フォルダの中にあるファイルは相対パスで持つ（持ち運べるように）。"""
        if not value:
            return ""
        value = os.path.abspath(value)
        try:
            rel = os.path.relpath(value, self.folder)
        except ValueError:      # ドライブが違うと相対化できない
            return value
        return value if rel.startswith("..") else rel.replace("\\", "/")

    @property
    def dxf_path(self):
        return self.resolve(self.dxf)

    @property
    def absorption_path(self):
        return self.resolve(self.absorption_csv)

    # ---- 保存・読み込み ------------------------------------------------

    def ensure_dirs(self):
        for folder in (self.folder, self.path(RESULT_DIR), self.path(FIGURE_DIR),
                       self.result_dir(), self.figure_dir()):
            os.makedirs(folder, exist_ok=True)
        return self

    def clear_results(self, verbose=True):
        """前回の計算結果と図を消す。

        条件を変えて回し直したとき、**前回の条件で作ったファイルが残っていると
        新しい結果と混ざる**。たとえば容積を空にして統計残響式が飛ばされると、
        前回の `rt_statistical.csv` だけが古いまま残り、いまの条件の値だと
        思って読んでしまう（実際に起きた）。

        消すのは**このプログラムが作るファイルだけ**。フォルダごと消すと
        利用者が置いた資料まで巻き添えにするので、名前を決め打ちにしてある。
        """
        removed = 0
        # 受音点ごとのものと、受音点に依らないものの両方（result_path が振り分ける）
        for key in RESULT_FILES:
            path = self.result_path(key)
            if os.path.exists(path):
                os.remove(path)
                removed += 1
        for path in (self.clarity_path(),):
            if os.path.exists(path):
                os.remove(path)
                removed += 1
        figures = self.figure_dir()
        if os.path.isdir(figures):
            for name in os.listdir(figures):
                if name.lower().endswith(".png"):
                    os.remove(os.path.join(figures, name))
                    removed += 1
        if removed and verbose:
            print(f"[project] 前回の結果と図を {removed} 個消しました"
                  f"（古い条件のファイルが混ざらないように）")
        return removed

    def to_dict(self):
        data = {key: getattr(self, key) for key in DEFAULTS}
        data["dxf"] = self._relative(self.dxf_path)
        data["absorption_csv"] = self._relative(self.absorption_path)
        for key in ("source", "receiver"):
            value = data[key]
            if value is not None:
                data[key] = [float(v) for v in np.asarray(value).ravel()]
        # 受音点ごとの向きはリストで持つ（1 点しかなければ数値のままでよい）
        if data["head_azimuth"] is not None and not np.isscalar(data["head_azimuth"]):
            data["head_azimuth"] = [float(v) for v in data["head_azimuth"]]
        return data

    def save(self):
        self.ensure_dirs()
        with open(self.path(PROJECT_FILE), "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, ensure_ascii=False, indent=2)
        return self.path(PROJECT_FILE)

    @classmethod
    def load(cls, folder):
        """フォルダから読む。project.json が無ければ既定値の新規プロジェクト。"""
        folder = os.path.abspath(folder)
        path = os.path.join(folder, PROJECT_FILE)
        if not os.path.exists(path):
            return cls(folder)
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        unknown = sorted(set(data) - set(DEFAULTS) - {"folder"})
        if unknown:
            print(f"[project] 注意: 知らない項目があります（無視します）: {unknown}")
        return cls(folder, **{k: v for k, v in data.items() if k in DEFAULTS})

    @classmethod
    def exists(cls, folder):
        return os.path.exists(os.path.join(os.path.abspath(folder), PROJECT_FILE))

    # ---- 法線の反転指定 ------------------------------------------------

    def load_flipped_faces(self):
        """`normals.json` から、反転する面のインデックス集合を読む。

        DXF が変わると面の番号がずれるので、**書いたときの DXF と面数を控えてある**。
        食い違ったら使わずに知らせる（黙って間違った面を反転させないため）。
        """
        path = self.path(NORMALS_FILE)
        if not os.path.exists(path):
            return set(), {}
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        return set(int(i) for i in data.get("flipped_faces", [])), data

    def save_flipped_faces(self, flipped, face_count, mode="auto"):
        self.ensure_dirs()
        data = {
            "dxf": self._relative(self.dxf_path),
            "face_count": int(face_count),
            "orient_normals": mode,
            "flipped_faces": sorted(int(i) for i in flipped),
        }
        with open(self.path(NORMALS_FILE), "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        return self.path(NORMALS_FILE)

    def flipped_faces_for(self, face_count):
        """面数が合っている場合だけ反転指定を返す。合わなければ空集合＋警告。"""
        flipped, data = self.load_flipped_faces()
        if not data:
            return set()
        saved = int(data.get("face_count", -1))
        if saved != face_count:
            print(f"[project] 警告: {NORMALS_FILE} は面 {saved} 枚のときの指定ですが、"
                  f"いま読んだ DXF は {face_count} 枚です。"
                  f"面の番号が対応しないので**使いません**。法線の確認をやり直してください")
            return set()
        return flipped

    # ---- 受音点ごとの顔の向き --------------------------------------------

    def head_azimuth_list(self, count):
        """受音点 `count` 個ぶんの正面方向 [度] のリストを返す。

        `head_azimuth` は数値でもリストでもよい。数値なら全点に同じ値を使う
        （従来の project.json との互換）。リストが短ければ足りない分は 0° で埋める。
        """
        value = self.head_azimuth
        if value is None:
            return [0.0] * count
        if np.isscalar(value):
            return [float(value)] * count
        values = [float(v) for v in value]
        return (values + [0.0] * count)[:count]

    def head_azimuth_for(self, index):
        """`index` 番目の受音点の正面方向 [度]。"""
        return self.head_azimuth_list(index + 1)[index]

    # ---- 面ごとの吸音材（materials.json） --------------------------------
    #
    # レイヤで吸音材を分けられないモデル（1 つの 3DSOLID で出来ていて面ごとの
    # レイヤが無いなど）のための逃げ道。**normals.json と同じ作り**にしてある
    # ＝ 面の番号で持ち、DXF と面数を控えて食い違いを検出する。
    #
    # 中身は「材料名 → 面番号の配列」。逆向き（面番号 → 材料名）より圧倒的に短く、
    # 人が開いたときも「どの材料をどこに貼ったか」が読み取れる。

    def load_face_materials(self):
        """`materials.json` から {面インデックス: 材料名} を読む。"""
        path = self.path(MATERIALS_FILE)
        if not os.path.exists(path):
            return {}, {}
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        assignment = {}
        for name, faces in (data.get("materials") or {}).items():
            for index in faces:
                assignment[int(index)] = name
        return assignment, data

    def save_face_materials(self, face_materials, face_count):
        """{面インデックス: 材料名} を `materials.json` に保存する。"""
        self.ensure_dirs()
        grouped = {}
        for index, name in sorted(face_materials.items()):
            if name:
                grouped.setdefault(name, []).append(int(index))
        data = {
            "dxf": self._relative(self.dxf_path),
            "face_count": int(face_count),
            "absorption_csv": self._relative(self.absorption_path),
            "materials": grouped,
        }
        with open(self.path(MATERIALS_FILE), "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        return self.path(MATERIALS_FILE)

    def face_materials_for(self, face_count):
        """面数が合っている場合だけ割り当てを返す。合わなければ空＋警告。"""
        assignment, data = self.load_face_materials()
        if not data:
            return {}
        saved = int(data.get("face_count", -1))
        if saved != face_count:
            print(f"[project] 警告: {MATERIALS_FILE} は面 {saved} 枚のときの割り当てですが、"
                  f"いま読んだ DXF は {face_count} 枚です。"
                  f"面の番号が対応しないので**使いません**。吸音材の割り当てをやり直してください")
            return {}
        return assignment

    # ---- 表示 ----------------------------------------------------------

    def summary(self):
        return (f"プロジェクト『{self.name}』\n"
                f"  フォルダ  {self.folder}\n"
                f"  モデル    {self.dxf}\n"
                f"  吸音率    {self.absorption_csv}"
                f"（{ {'normal': '垂直入射', 'random': '残響室法'}.get(self.absorption_kind, '未指定')}）\n"
                f"  音線 {self.rays} 本 / 最大反射 {self.nref} 回 / 受音球 {self.radius} m\n"
                f"  {self.temperature}℃ / 湿度 {self.humidity}% / {self.pressure} kPa"
                f" / {self.band_number} バンド")


# ------------------------------------------------------------------------------
# 結果の読み書き
# ------------------------------------------------------------------------------

def write_surface_csv(filename, surface, frequencies):
    """材料別の面積と吸音率を CSV にする（`reverberation.surface_summary` の結果）。

    区分は `Mesh.material`。面ごとの吸音材を割り当てていなければ DXF のレイヤ名、
    割り当てていればその材料名になる（`read_dxffile.read_model(face_materials=...)`）。
    """
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    import table as tb
    # 1 行が材料なので、周波数は**列**（table.py の共通ルール）
    header = ["material", "area_m2"] + [tb.band_column("alpha", f)
                                        for f in frequencies]
    lines = [",".join(header)]
    for name, area, alpha in zip(surface["names"], surface["areas"],
                                 surface["absorption"]):
        cells = [str(name), f"{area:.6f}"] + [f"{a:.6f}" for a in alpha]
        lines.append(",".join(cells))
    with open(filename, "w", encoding="utf-8-sig", newline="") as f:
        f.write("\n".join(lines) + "\n")
    return filename


def _read_csv(path):
    if not os.path.exists(path):
        return None
    return np.genfromtxt(path, delimiter=",", names=True, encoding="utf-8-sig")


# 周波数を横に並べてある CSV（`table.write_frequency_table` が書いたもの）。
# 1 行が指標なので、普通の genfromtxt では読めない
FREQUENCY_TABLES = ("rt", "statistical", "clarity")


def load_results(project):
    """プロジェクトフォルダの結果を読み戻す。

    戻り値の dict。**無いものは None**（途中まで計算した状態も開けるように）。

    `rt` / `statistical` / `clarity` は**周波数が横**に並んだ表なので、
    `{'frequencies': (nf,), 'rows': {行の名前: (nf,)}}` の形で返す
    （`table.read_frequency_table`。古い縦向きのファイルもそのまま読める）。
    それ以外は `np.genfromtxt` の構造化配列。
    """
    import table as tb

    result = {}
    for key in ("pulses", "ir", "decay", "surface"):
        result[key] = _read_csv(project.result_path(key))
    for key in FREQUENCY_TABLES:
        path = (project.path(RESULT_DIR, "clarity.csv") if key == "clarity"
                else project.result_path(key))
        frequencies, rows = tb.read_frequency_table(path)
        result[key] = None if frequencies is None else {"frequencies": frequencies,
                                                        "rows": rows}
    raylog = project.result_path("raylog")
    result["raylog"] = raylog if os.path.exists(raylog) else None
    return result


def has_results(project):
    return os.path.exists(project.result_path("rt")) or \
           os.path.exists(project.result_path("pulses"))
