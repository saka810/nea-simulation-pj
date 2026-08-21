"""プロジェクト（計算条件・法線の修正・計算結果）の保存と読み込み。

**ねらい**：一度作った条件と結果を、あとから開き直せるようにする。
条件を JSON、結果を CSV、図を PNG で**プロジェクトフォルダに全部置く**。
CSV にしてあるのは Excel でそのまま開けるようにするため。

    プロジェクトフォルダ/
      project.json          条件（DXF・吸音率・音線数・受音球・温度湿度…）
      normals.json          法線の反転指定（面ごと。face_editor.py が書く）
      materials.json        面ごとの吸音材の割り当て（同上。レイヤで分けられないモデル用）
      結果/
        研修室_条件A_まとめ_残響時間.csv    全受音点 ＋ 平均 ＋ 理論値（summary.py が書く）
        研修室_条件A_まとめ_明瞭度.csv      全受音点 ＋ 平均
        研修室_条件A_吸音率と理論値.csv     材料別の吸音率 → 平均吸音率 → 残響時間理論値
        研修室_条件A_raylog.npz             音線軌跡（可視化用。可変長なので npz）
        rec1/               ← **受音点ごと**
          研修室_条件A_pulses.csv   パルス列（反射回数・到来時刻・到来方向・エネルギー）
          研修室_条件A_ir.csv       インパルス応答
          研修室_条件A_rt.csv       残響指標 EDT / T20 / T30
          研修室_条件A_decay.csv    減衰曲線
          研修室_条件A_clarity.csv  明瞭度 C50 / C80 / D50 / Ts
        rec2/ …
      図/
        rec1/ 研修室_条件A_*.png   正規化したインパルス応答・減衰曲線・残響時間 ほか
        rec2/ …
        画面/               画面から手で撮った画像・動画

★**受音点ごとのものは `結果/recN/`、受音点に依らないものは `結果/` 直下**
（2026-08-21 にこの形へ。それまでは 1 点目だけ `結果/` 直下、2 点目以降が
`rec2/結果/` という不揃いな置き方で、ユーザー指摘で直した）。
統計残響式・材料別面積・音線軌跡は受音点に依らないので 1 つだけ持つ。

★**ファイル名の頭に「対象室＋条件名」を付ける**（2026-08-21 ユーザー要望）。
プロジェクト名（`name`。既定はフォルダ名）をそのまま頭に付ける。
結果ファイルは報告書やメールでフォルダの外へ出ることが多く、
`rt.csv` のままだとどの室・どの条件のものか分からなくなるため。
図（PNG）にも同じ頭を付ける（貼り付けたあとで見分けが付くように）。
**読むときは頭が付いていないファイルも探す**ので、
名前を変える前に計算したプロジェクトもそのまま開ける。

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

# 材料条件表の既定名（`condition_table.CONDITION_FILE` と同じ）。
# この名前のときは条件名を付けない（条件を分けていないということなので）
DEFAULT_CONDITION_FILE = "材料条件表.csv"
DEFAULT_CONDITION_STEM = "材料条件表"

# 経路の幾何（吸音材を変えた再計算に使う）。**条件名は付けない**
PATHS_FILE = "経路.npz"

# 結果ファイルの名前。**キーはコード側の呼び名**で、値が実ファイル名。
# 実際のファイル名にはプロジェクト名（対象室＋条件名）が頭に付く（`Project.prefixed`）
RESULT_FILES = {
    "pulses": "pulses.csv",
    "ir": "ir.csv",
    "rt": "rt.csv",
    "decay": "decay.csv",
    "clarity": "clarity.csv",
    "spl": "spl.csv",
    "sti": "sti.csv",
    "room": "吸音率と理論値.csv",
    "raylog": "raylog.npz",
    # 経路の幾何（反射面の並びと入射角）。**吸音材を変えた再計算に使う**（F-9）
    "paths": PATHS_FILE,
}

# **受音点に依らない**結果。受音点ごとのフォルダではなく `結果/` 直下に置く。
#   室の吸音と理論値 … 室形状と材料だけで決まる
#   音線軌跡         … 音源から出た音線の形。受音点をまたいで共有している（F-6）
SHARED_RESULTS = {"room", "raylog"}

# **条件（吸音材）に依らない**結果。ファイル名に条件名を付けず、対象室名だけにする。
#   経路の幾何 … 吸音に依らない（それを使い回すのがこの仕組みの目的）
#   音線軌跡   … 形は吸音に依らない（色分けに使うエネルギーだけ条件に依る）
ROOM_SCOPED_RESULTS = {"paths", "raylog"}

# `clear_results()` で**消さない**結果。作り直すのが高くつき、
# かつ中身が古いかどうかを自分で判定できるもの（経路は指紋を突き合わせる）
KEEP_ON_CLEAR = {"paths"}

# 昔の名前。**読むときだけ**探す（作り直す前のプロジェクトを開けるように）。
# `clear_results()` はこちらも消す（古い条件のファイルが残って混ざらないように）
LEGACY_RESULT_FILES = {
    "room": ["rt_statistical.csv", "surface.csv"],
}

# ファイル名に使えない文字（Windows）。対象室・条件名から作るので置き換える
UNSAFE_CHARACTERS = '\\/:*?"<>|'


def _stem(path_or_name):
    """パスから拡張子なしのファイル名を取り出す（空なら ""）。"""
    if not path_or_name:
        return ""
    name = str(path_or_name).replace(chr(92), "/")
    return os.path.splitext(os.path.basename(name))[0]


def safe_name(text):
    """ファイル名の一部に使える形に直す。使えない文字と空白は `_` にする。"""
    text = (text or "").strip()
    for character in UNSAFE_CHARACTERS:
        text = text.replace(character, "_")
    return "_".join(text.split())


# project.json に書き出す条件と既定値。
# ここに無いキーは保存されないので、**新しい計算条件を足したらここにも足すこと**
DEFAULTS = {
    "name": "",
    "dxf": "",
    "absorption_csv": "",
    "absorption_kind": None,       # 'normal' | 'random' | None（CSV の # kind: を見る）
    # 材料条件表（レイヤー名 → 吸音材）の CSV。**条件名はこのファイル名から取る**。
    # 同じフォルダに条件ごとのファイルを置いて選び替える使い方（2026-08-21 ユーザー要望）。
    # 空なら既定名 `材料条件表.csv` を探す
    "condition_csv": "",
    "assignment": None,            # レイヤ → 材料の対応（辞書。条件表が優先）
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
    # 点音源のパワーレベル PWL [dB]。数値なら全帯域同じ、リストなら帯域ごと。
    # **None なら相対値**（Lw = 0 dB として計算し、出力にそう書く）。
    # 絶対値の音圧レベルを出すのに要る（`sound_level.py`）
    "source_power_db": None,
    # 背景騒音の音圧レベル [dB]（帯域ごと。数値なら全帯域同じ）。
    # STI の SNR に使う。**PWL と両方そろっていないと使えない**（絶対値が要る）
    "noise_level_db": None,
}


class Project:
    """計算条件の入れ物。`folder` に紐づく。

    属性は `DEFAULTS` のキーがそのまま生える（`project.rays` のように読める）。
    """

    def __init__(self, folder, **values):
        self.folder = os.path.abspath(folder)
        for key, default in DEFAULTS.items():
            setattr(self, key, values.get(key, default))
        # ★`name` は**空のままにしておく**（既定はフォルダ名ではない）。
        #   空なら DXF のファイル名を対象室の名前として使う（2026-08-21 ユーザー要望。
        #   「DXF のファイル名は物件名や部屋名、条件表のファイル名は条件名にする」）
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

    # ---- ファイル名の頭（対象室＋条件名）------------------------------
    #
    # 結果は報告書やメールでフォルダの外へ出るので、ファイル名だけで
    # どの室・どの条件のものか分かるようにする（2026-08-21 ユーザー要望）

    @property
    def display_name(self):
        """画面やログに出す名前。`name` が空なら DXF 名、それも無ければフォルダ名。"""
        return (self.name or _stem(self.dxf) or os.path.basename(self.folder))

    @property
    def room_label(self):
        """**対象室の名前**。`name` の指定が最優先、無ければ DXF のファイル名。

        条件をまたいで共通なので、条件に依らない置き場（経路のキャッシュ）に使う。
        """
        return safe_name(self.name or _stem(self.dxf)
                         or os.path.basename(self.folder))

    @property
    def condition_label(self):
        """**条件の名前**。材料条件表のファイル名から取る。既定名なら空。"""
        stem = _stem(self.condition_csv)
        return "" if not stem or stem == DEFAULT_CONDITION_STEM else safe_name(stem)

    @property
    def file_prefix(self):
        """結果ファイル名の頭 ＝ 対象室名 ＋ 条件名。

        ★**DXF のファイル名（物件名・部屋名）と条件表のファイル名（条件名）を
        そのまま使う**（2026-08-21 ユーザー要望）。同じフォルダに条件ごとの
        結果を並べても混ざらない。
        """
        room, condition = self.room_label, self.condition_label
        if room and condition:
            return f"{room}_{condition}"
        return room or condition

    def prefixed(self, filename):
        """`rt.csv` → `研修室_条件A_rt.csv`。名前が空なら元のまま。"""
        prefix = self.file_prefix
        return f"{prefix}_{filename}" if prefix else filename

    def result_path(self, key):
        """結果ファイルの**書き出し先**。`key` は RESULT_FILES のキー。

        受音点に依らないもの（`SHARED_RESULTS`）は `結果/` 直下、
        それ以外は `結果/recN/` に置く。名前には対象室＋条件名が頭に付く。
        """
        return os.path.join(self.result_dir(shared=key in SHARED_RESULTS),
                            self._named(key, RESULT_FILES[key]))

    def result_candidates(self, key):
        """その結果として**読める名前**を、探す順に返す。

        ① 対象室＋条件名の付いた今の名前
        ② 頭の付いていない名前（頭を付ける前に計算したプロジェクト）
        ③ さらに古い名前（`rt_statistical.csv` など。`LEGACY_RESULT_FILES`）
        """
        folder = self.result_dir(shared=key in SHARED_RESULTS)
        names = [RESULT_FILES[key]] + LEGACY_RESULT_FILES.get(key, [])
        paths = []
        for name in names:
            paths.append(os.path.join(folder, self._named(key, name)))
            if self.file_prefix:
                paths.append(os.path.join(folder, name))
        return paths

    def _named(self, key, filename):
        """結果ファイル名。**条件に依らないものは対象室名だけ**を頭に付ける。"""
        if key in ROOM_SCOPED_RESULTS:
            room = self.room_label
            return f"{room}_{filename}" if room else filename
        return self.prefixed(filename)

    def existing_result_path(self, key):
        """**読むため**のパス。無ければ昔の名前も探す。全部無ければ今の名前を返す。"""
        for path in self.result_candidates(key):
            if os.path.exists(path):
                return path
        return self.result_path(key)

    def clarity_path(self):
        """明瞭度の CSV（`結果/recN/…clarity.csv`）。"""
        return self.result_path("clarity")

    def figure_dir(self):
        """図の置き場。受音点が決まっていれば `図/recN/`。"""
        if self.receiver_index is None:
            return self.path(FIGURE_DIR)
        return self.path(FIGURE_DIR, RECEIVER_DIR % self.receiver_index)

    def figure_path(self, name):
        """図のパス。**図にも対象室＋条件名を付ける**（貼ってから見分けが付くように）。"""
        return os.path.join(self.figure_dir(), self.prefixed(name))

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

    @property
    def condition_path(self):
        """材料条件表のパス。指定が無ければ既定名（`材料条件表.csv`）。"""
        return self.resolve(self.condition_csv) or self.path(DEFAULT_CONDITION_FILE)

    def paths_cache(self):
        """経路の幾何（`path_cache`）の置き場（`結果/recN/<室>_経路.npz`）。

        ★**条件名を付けない。**条件（吸音材）が違っても経路は同じなので、
        条件をまたいで共有する（それがこの仕組みの目的）。
        """
        return self.result_path("paths")

    # ---- 保存・読み込み ------------------------------------------------

    def ensure_dirs(self):
        for folder in (self.folder, self.path(RESULT_DIR), self.path(FIGURE_DIR),
                       self.result_dir(), self.figure_dir()):
            os.makedirs(folder, exist_ok=True)
        return self

    def clear_results(self, verbose=True, keep=()):
        """前回の計算結果と図を消す。

        条件を変えて回し直したとき、**前回の条件で作ったファイルが残っていると
        新しい結果と混ざる**。たとえば容積を空にして統計残響式が飛ばされると、
        前回の `rt_statistical.csv` だけが古いまま残り、いまの条件の値だと
        思って読んでしまう（実際に起きた）。

        消すのは**このプログラムが作るファイルだけ**。フォルダごと消すと
        利用者が置いた資料まで巻き添えにするので、名前を決め打ちにしてある。

        ★**昔の名前のファイルも消す**（`result_candidates`）。名前を変える前の
        `rt.csv` や `rt_statistical.csv` が残っていると、今回の結果と並んでしまう。

        ★**経路の幾何（`経路.npz`）は消さない**（`KEEP_ON_CLEAR`）。
        作り直すのに音線追跡が丸ごと要るうえ、古いかどうかは指紋で判定できるため。
        `keep` に鍵を足せば他のものも残せる（経路を使い回すときの音線軌跡など）。
        """
        removed = 0
        skip = set(KEEP_ON_CLEAR) | set(keep)
        # 受音点ごとのものと、受音点に依らないものの両方（result_path が振り分ける）
        for key in RESULT_FILES:
            if key in skip:
                continue
            for path in self.result_candidates(key):
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
        if self.condition_csv:
            data["condition_csv"] = self._relative(self.resolve(self.condition_csv))
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

    # ---- 帯域ごとの値（音源パワーレベル・背景騒音）--------------------

    def band_values(self, key, band_number=None):
        """`source_power_db` のような「数値でもリストでもよい」条件を配列で返す。

        数値なら全帯域に同じ値、リストなら帯域ごと（足りない分は最後の値で伸ばす）。
        未設定なら None。**None と 0 は意味が違う**（未入力か 0 dB か）ので、
        呼ぶ側は None を「未入力」として扱うこと。
        """
        value = getattr(self, key)
        if value is None:
            return None
        count = band_number or self.band_number
        if np.isscalar(value):
            return np.full(count, float(value))
        values = [float(v) for v in value]
        if not values:
            return None
        return np.array((values + [values[-1]] * count)[:count])

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
        # ★受音点の番号は `receiver_index` から出す。**`name` に足してはいけない**
        #   （`name` は結果ファイル名の頭に付くので、受音点ごとに変わると
        #     ファイル名が受音点ごとに違ってしまう）
        who = "" if self.receiver_index is None else f"（受音点 {self.receiver_index}）"
        condition = f" / 条件『{self.condition_label}』" if self.condition_label else ""
        return (f"プロジェクト『{self.display_name}』{condition}{who}\n"
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

# 「吸音率と理論値」の CSV の区分（この順に並べる。2026-08-21 ユーザー指定）
ROOM_SECTIONS = ("材料別の吸音率", "平均吸音率", "残響時間理論値")

# 「平均吸音率」の区分に入れる行。**平均吸音率だけでなく等価吸音面積も入れる**。
# 統計残響式は A = Sᾱ（Sabine）／-S ln(1-ᾱ)（Eyring）から T を出すので、
# ᾱ と A を並べておくと理論値がどこから来たのか表の中で追える
ROOM_MEAN_ROWS = (("平均吸音率", "mean_absorption"),
                  ("等価吸音面積_m2", "equivalent_area"),
                  ("空気吸収_4mV_m2", "air_absorption_area"))

# 「残響時間理論値」の区分に入れる行（`reverberation.STATISTICAL_LABELS` と同じ 3 式）
ROOM_STATISTICAL_ROWS = (("sabine_s", "sabine"),
                         ("eyring_s", "eyring"),
                         ("eyring_knudsen_s", "eyring_knudsen"))


def write_room_csv(filename, statistical, frequencies=None):
    """室の吸音と残響時間理論値を**1 つの CSV**にする。

    元は `rt_statistical.csv`（統計残響式）と `surface.csv`（材料別の面積・吸音率）
    に分かれていたが、**どちらも受音点に依らない室の性質**で、見るときは
    「材料の吸音率 → 平均吸音率 → 理論値」と続けて追うので 1 枚にした
    （2026-08-21 ユーザー要望。並び順もユーザー指定）。

        区分,項目,面積_m2,63,125,…
        材料別の吸音率,コンクリート,120.5,0.02,…      ← 材料ごとの α（乱入射）
        材料別の吸音率,吸音板,45.0,0.15,…
        平均吸音率,平均吸音率,340.2,0.153,…          ← 面積で重み付けした ᾱ
        平均吸音率,等価吸音面積_m2,,52.1,…
        平均吸音率,空気吸収_4mV_m2,,0.0,…
        残響時間理論値,sabine_s,,1.91,…
        残響時間理論値,eyring_s,,…
        残響時間理論値,eyring_knudsen_s,,…

    **周波数は横**（`table.py` の共通ルール。「区分付きの表」の形）。

    引数:
        statistical … `reverberation.statistical_reverberation()` の戻り値
                      （`['surface']` に材料別の面積・吸音率が入っている）
    """
    import table as tb

    if frequencies is None:
        frequencies = statistical["frequencies"]
    surface = statistical["surface"]

    rows = [(ROOM_SECTIONS[0], name, area, alpha)
            for name, area, alpha in zip(surface["names"], surface["areas"],
                                         surface["absorption"])]
    for label, key in ROOM_MEAN_ROWS:
        # 面積の欄は平均吸音率の行だけ埋める（総表面積）。他は帯域の値だけ
        area = surface["total_area"] if key == "mean_absorption" else None
        rows.append((ROOM_SECTIONS[1], label, area, statistical[key]))
    for label, key in ROOM_STATISTICAL_ROWS:
        rows.append((ROOM_SECTIONS[2], label, None, statistical[key]))
    return tb.write_sectioned_table(filename, frequencies, rows,
                                    value_label="面積_m2")


def read_room_csv(path):
    """`write_room_csv` が書いた CSV を読む。読めなければ None。

    戻り値: dict
        'frequencies' (nf,)
        'surface'  … `reverberation.surface_summary()` と同じ形
                     （'names' / 'areas' / 'absorption' / 'total_area'）
        'rows'     … 区分を除いた {項目: (nf,)}（`平均吸音率` / `sabine_s` …）
        'sections' … {区分: {項目: (nf,)}}
    """
    import table as tb

    table = tb.read_sectioned_table(path)
    if table is None or table["labels"][0] != "区分":
        return None

    names, areas, alphas = [], [], []
    for section, item in table["order"]:
        if section != ROOM_SECTIONS[0]:
            continue
        names.append(item)
        try:
            areas.append(float(table["values"][item]))
        except ValueError:
            areas.append(np.nan)
        alphas.append(table["sections"][section][item])

    surface = _surface_dict(names, areas, alphas, len(table["frequencies"]))
    return {"frequencies": table["frequencies"], "surface": surface,
            "rows": table["rows"], "sections": table["sections"]}


def _surface_dict(names, areas, alphas, band_number):
    """`reverberation.surface_summary()` と同じ形の dict を作る。"""
    return {"names": list(names), "areas": np.array(areas, dtype=float),
            "absorption": (np.array(alphas, dtype=float) if len(alphas)
                           else np.zeros((0, band_number))),
            "total_area": float(np.nansum(areas)) if len(areas) else 0.0,
            "materials": {n: {"area": a, "absorption": np.asarray(v)}
                          for n, a, v in zip(names, areas, alphas)}}


def _read_csv(path):
    if not os.path.exists(path):
        return None
    return np.genfromtxt(path, delimiter=",", names=True, encoding="utf-8-sig")


# 周波数を横に並べてある CSV（`table.write_frequency_table` が書いたもの）。
# 1 行が指標なので、普通の genfromtxt では読めない
FREQUENCY_TABLES = ("rt", "clarity")


def load_results(project):
    """プロジェクトフォルダの結果を読み戻す。

    戻り値の dict。**無いものは None**（途中まで計算した状態も開けるように）。

    `rt` / `clarity` は**周波数が横**に並んだ表なので、
    `{'frequencies': (nf,), 'rows': {行の名前: (nf,)}}` の形で返す
    （`table.read_frequency_table`。古い縦向きのファイルもそのまま読める）。
    `room` は `read_room_csv` の形。`statistical` / `surface` はその中身を
    昔と同じ形で取り出したもの（呼ぶ側を変えずに済むように）。
    それ以外は `np.genfromtxt` の構造化配列。
    """
    import table as tb

    result = {}
    for key in ("pulses", "ir", "decay"):
        result[key] = _read_csv(project.existing_result_path(key))
    for key in FREQUENCY_TABLES:
        frequencies, rows = tb.read_frequency_table(
            project.existing_result_path(key))
        result[key] = None if frequencies is None else {"frequencies": frequencies,
                                                        "rows": rows}

    # 室の吸音と理論値。**1 枚にまとめる前の 2 ファイルもそのまま読める**
    room = read_room_csv(project.existing_result_path("room"))
    if room is None:
        room = _load_legacy_room(project)
    result["room"] = room
    result["statistical"] = (None if room is None else
                             {"frequencies": room["frequencies"],
                              "rows": room["rows"]})
    result["surface"] = None if room is None else room["surface"]

    raylog = project.existing_result_path("raylog")
    result["raylog"] = raylog if os.path.exists(raylog) else None
    return result


def _load_legacy_room(project):
    """`rt_statistical.csv` ＋ `surface.csv` に分かれていた頃の結果を読む。

    1 枚にまとめる前（2026-08-21 より前）に計算したプロジェクトを
    開き直したときのため。**見え方を今の形に揃えて返す**。
    """
    import csv as _csv
    import table as tb

    folder = project.result_dir(shared=True)
    frequencies, rows = tb.read_frequency_table(
        os.path.join(folder, "rt_statistical.csv"))
    surface_path = os.path.join(folder, "surface.csv")
    if frequencies is None and not os.path.exists(surface_path):
        return None

    names, areas, alphas = [], [], []
    if os.path.exists(surface_path):
        # ★`np.genfromtxt` では読めない（材料名の列が float 扱いで nan になる）
        with open(surface_path, encoding="utf-8-sig", newline="") as f:
            table = [row for row in _csv.reader(f) if row and row[0].strip()]
        header, body = table[0], table[1:]
        columns = [i for i, name in enumerate(header) if name.startswith("alpha_")]
        for row in body:
            names.append(row[0].strip())
            areas.append(float(row[1]))
            alphas.append([float(row[i]) for i in columns])
        if frequencies is None:
            frequencies = np.array([float(header[i].split("_")[1].rstrip("Hz"))
                                    for i in columns])
    surface = _surface_dict(names, areas, alphas, len(frequencies))
    return {"frequencies": frequencies, "surface": surface,
            "rows": rows, "sections": {}}


def has_results(project):
    """計算結果があるか（条件入力画面の「計算結果あり」表示に使う）。

    ★**受音点ごとのフォルダ（`結果/rec1/`）も見る。**
    置き場を `結果/recN/` に変えたとき、`結果/` 直下しか見ていなかったので
    結果があるのに「ありません」と言われていた。
    """
    saved = project.receiver_index
    try:
        for index in ((None, 1) if saved is None else (saved,)):
            project.receiver_index = index
            if any(os.path.exists(project.existing_result_path(key))
                   for key in ("rt", "pulses")):
                return True
    finally:
        project.receiver_index = saved
    return False
