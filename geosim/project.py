"""プロジェクト（計算条件・法線の修正・計算結果）の保存と読み込み。

**ねらい**：一度作った条件と結果を、あとから開き直せるようにする。
条件を JSON、結果を CSV、図を PNG で**プロジェクトフォルダに全部置く**。
CSV にしてあるのは Excel でそのまま開けるようにするため。

    プロジェクトフォルダ/
      project.json          条件（DXF・吸音率・音線数・受音球・温度湿度…）
      normals.json          法線の反転指定（面ごと。normal_editor.py が書く）
      結果/
        pulses.csv          パルス列（反射回数・到来時刻・到来方向・バンド別エネルギー）
        ir.csv              インパルス応答
        rt.csv              残響指標 EDT / T20 / T30
        rt_statistical.csv  統計残響式 Sabine / Eyring / Millington
        decay.csv           減衰曲線
        surface.csv         レイヤ別の面積と吸音率
        raylog.npz          音線軌跡（可視化用。可変長なので npz）
      図/
        *.png               正規化したインパルス応答・減衰曲線・残響時間 ほか

DXF や吸音率 CSV は**プロジェクトフォルダからの相対パスで持つ**（フォルダごと
別の端末へ移してもそのまま開ける）。フォルダの外にある場合は絶対パスのままにする。
"""

import json
import os

import numpy as np

PROJECT_FILE = "project.json"
NORMALS_FILE = "normals.json"
RESULT_DIR = "結果"
FIGURE_DIR = "図"
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

    # ---- パス ---------------------------------------------------------

    def path(self, *parts):
        return os.path.join(self.folder, *parts)

    def result_path(self, key):
        """結果ファイルのパス。`key` は RESULT_FILES のキー。"""
        return self.path(RESULT_DIR, RESULT_FILES[key])

    def figure_path(self, name):
        return self.path(FIGURE_DIR, name)

    def screenshot_dir(self):
        """画面から手で撮った画像・動画の置き場（`図/画面/`）。"""
        folder = self.path(SCREENSHOT_DIR)
        os.makedirs(folder, exist_ok=True)
        return folder

    def ascii_tag(self):
        """ウィンドウのタイトルの予備に使う短い名前（ASCII で書けるときだけ）。

        VTK のウィンドウは日本語のタイトルが化けることがあり、その場合に
        英字の題へ差し替える（`view_model_gui.set_window_title`）。
        プロジェクト名は日本語のことが多いので**フォルダ名のほうを先に見る**。
        """
        for candidate in (os.path.basename(self.folder), self.name):
            candidate = (candidate or "").strip()
            if candidate and all(ord(c) < 128 for c in candidate):
                return candidate
        return ""

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
        for sub in ("", RESULT_DIR, FIGURE_DIR):
            os.makedirs(self.path(sub) if sub else self.folder, exist_ok=True)
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
        for key in RESULT_FILES:
            path = self.result_path(key)
            if os.path.exists(path):
                os.remove(path)
                removed += 1
        for name in ("clarity.csv",):
            path = self.path(RESULT_DIR, name)
            if os.path.exists(path):
                os.remove(path)
                removed += 1
        figures = self.path(FIGURE_DIR)
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
    """レイヤ別の面積と吸音率を CSV にする（`reverberation.surface_summary` の結果）。"""
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    header = ["material", "area_m2"] + [f"alpha_{int(f)}Hz" for f in frequencies]
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


def load_results(project):
    """プロジェクトフォルダの結果を読み戻す。

    戻り値の dict。**無いものは None**（途中まで計算した状態も開けるように）。
    """
    result = {}
    for key in ("pulses", "ir", "rt", "statistical", "decay", "surface"):
        result[key] = _read_csv(project.result_path(key))
    raylog = project.result_path("raylog")
    result["raylog"] = raylog if os.path.exists(raylog) else None
    return result


def has_results(project):
    return os.path.exists(project.result_path("rt")) or \
           os.path.exists(project.result_path("pulses"))
