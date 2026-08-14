"""虚音源法バックトレース。元コード backtrace.f90 876〜1134 行。

目的:
    音線法（3）で見つかり重複削除（4）を通った「反射面の並び」ごとに、
    **虚音源を鏡像で組み立て、受音点から逆向きに経路を検算**して、
    到来時刻・到来方向・バンド別エネルギーのパルス列を作る。

なぜ 2 段階なのか:
    音線法は「どの壁をどの順に反射する経路がありそうか」を**発見**するだけの手段で、
    受音球に入ったかどうかで判定しているぶん位置に誤差がある。
    虚音源法は経路が分かっていれば**厳密**な到来時刻とエネルギーを出せる。
    そこで「探索は音線法、計算は虚音源法」という役割分担になっている。

このモジュールが検算していること（＝経路が却下される条件）:
    ・各反射で本当にその壁に当たるか（`jtmp .ne. tracred(i,k)` → 却下）
    ・途中で別の壁に遮られないか（最寄りの壁が経路上の壁でなければ却下）
    ・最後に音源が見えているか（音源より手前に壁があれば却下）
"""

import numpy as np

import mesh_method as mm
import sound_ray as sr
from atmosphere import Atmosphere

# 音速 [m/s]。**基準大気（20℃ / 湿度 40% / 101.325 kPa）から計算した値**で 343.8 m/s。
# 元コードは `c0 = 340.0d0` の定数だったが、これは約 14℃ 相当にあたる。
# 温度・湿度を変えたい場合は `atmosphere.Atmosphere` を作って
# その `sound_velocity` を渡すこと（空気吸収も同じ大気条件から計算される）。
SOUND_VELOCITY = Atmosphere().sound_velocity


class PulseList:
    """バックトレースが出すパルス列。元コード 1080 / 1124 行の write に対応。

    元コードの出力は 11 列:

        ktmp, rtime, -vtgt(1:3), enertmp(1:6)
         ↑      ↑        ↑           ↑
        反射回数 到来時刻  到来方向    バンド別エネルギー

    Python 版も中身は同じだが、**到来方向は単位ベクトルにして距離を別に持つ**。
    元コードの `-vtgt` は正規化されておらず、その長さが経路長そのものなので、
    「方向」と「距離」が 1 つの列に混ざっていた。分けておくほうが
    出力⑤（伝搬方向の可視化）でそのまま使える。
    `direction * distance` が元コードの `-vtgt` に一致する。

    属性:
        reflection_count : (n,)  int   反射回数（元 ktmp）
        time             : (n,)  float 到来時刻 [s]（元 rtime）
        distance         : (n,)  float 経路長 [m]（= time * 音速）
        direction        : (n,3) float 到来方向の単位ベクトル（受音点 → 虚音源）
        energy           : (n,b) float バンド別エネルギー（元 enertmp）
    """

    def __init__(self, band_number, sound_velocity=SOUND_VELOCITY):
        self.band_number = band_number
        self.sound_velocity = sound_velocity
        self.reflection_count = np.zeros(0, dtype=int)
        self.time = np.zeros(0)
        self.distance = np.zeros(0)
        self.direction = np.zeros((0, 3))
        self.energy = np.zeros((0, band_number))

    def __len__(self):
        return len(self.time)

    @classmethod
    def from_records(cls, records, band_number, sound_velocity=SOUND_VELOCITY):
        self = cls(band_number, sound_velocity)
        if not records:
            return self
        self.reflection_count = np.array([r["reflection_count"] for r in records], dtype=int)
        self.time = np.array([r["time"] for r in records], dtype=float)
        self.distance = np.array([r["distance"] for r in records], dtype=float)
        self.direction = np.array([r["direction"] for r in records], dtype=float)
        self.energy = np.array([r["energy"] for r in records], dtype=float)
        return self

    def sort_by_time(self):
        """到来時刻の昇順に並べ替える（元コードは音線の発見順のまま）。

        インパルス応答の合成は順序に依存しないが、目視・デバッグでは
        時刻順のほうが圧倒的に読みやすい。直接音が先頭に来る。
        """
        order = np.argsort(self.time, kind="stable")
        self.reflection_count = self.reflection_count[order]
        self.time = self.time[order]
        self.distance = self.distance[order]
        self.direction = self.direction[order]
        self.energy = self.energy[order]
        return self

    def save_csv(self, filename):
        """パルス列を CSV に保存する。

        列: reflection_count, time_s, distance_m, dir_x, dir_y, dir_z, energy_1..energy_b
        元コードの 11 列との対応は本クラスの docstring を参照。
        """
        header = ["reflection_count", "time_s", "distance_m", "dir_x", "dir_y", "dir_z"]
        header += [f"energy_{b + 1}" for b in range(self.band_number)]
        rows = np.column_stack([self.reflection_count, self.time, self.distance,
                                self.direction, self.energy])
        np.savetxt(filename, rows, delimiter=",", header=",".join(header),
                   comments="", fmt="%.12g")
        return filename

    def summary(self):
        if len(self) == 0:
            return "パルス列: 0 本（受音に至った経路なし）"
        total = self.energy.sum(axis=0)
        return (f"パルス列: {len(self)} 本 / "
                f"到来時刻 {self.time.min() * 1000.0:.2f}〜{self.time.max() * 1000.0:.2f} ms / "
                f"反射回数 0〜{int(self.reflection_count.max())} 回 / "
                f"バンド別エネルギー合計 {np.array2string(total, precision=4)}")


def image_sources(soundsource_point, wall_ids, mesh):
    """反射面の並びから虚音源の列を作る。元コード 900〜926 行。

    戻り値は (len(wall_ids)+1, 3)。先頭が実音源、以降が 1 回・2 回…反射の虚音源。

    【元コードとの相違点】元コードは面までの距離に **abs() を掛けている**（911行）:

        temp = |n・p + d| / |n|
        isrc(k) = isrc(k-1) - 2 * temp * n

    鏡像の正しい式は符号付きで

        isrc(k) = isrc(k-1) - 2 * (n・p + d) / |n|^2 * n

    で、両者が一致するのは p が面の**表側**（n・p + d > 0）にあるときだけ。
    裏側にある虚音源に abs() を使うと、面を跨いで鏡像を作るかわりに
    **面から遠ざかる方向へ動いてしまう**。凸な部屋なら表側にしかならないので
    元コードでも実害は出にくいが、凹んだ形状や法線を反転させたモデルでは崩れる。
    ここでは符号付き（数学的に正しい形）で実装している。
    """
    wall_ids = list(wall_ids)
    images = np.empty((len(wall_ids) + 1, 3))
    images[0] = np.asarray(soundsource_point, dtype=float)

    for k, wall_id in enumerate(wall_ids):
        normal = np.asarray(mesh[wall_id].normal, dtype=float)
        norm = np.linalg.norm(normal)
        if norm == 0.0:
            raise ValueError(f"面 {wall_id} の法線が零ベクトルです")
        normal = normal / norm
        d = mm.parameter_d(normal, mesh[wall_id].vertexes[0])
        signed_distance = np.dot(normal, images[k]) + d
        images[k + 1] = images[k] - 2.0 * signed_distance * normal

    return images


def backtrace_path(soundsource_point, reciever_point, wall_ids, mesh,
                   band_number, sound_velocity=SOUND_VELOCITY, faces=None):
    """経路 1 本分のバックトレース。元コード 948〜1132 行。

    引数:
        wall_ids : list[int] 反射面 ID の並び（先頭の番兵は含めないこと）
                   空リストなら直接音の経路。

    戻り値:
        受音に至れば dict、途中で却下されれば None。
    """
    if faces is None:
        faces = mm.FaceArrays(mesh)
    images = image_sources(soundsource_point, wall_ids, mesh)
    ktmp = len(wall_ids)

    energy = np.ones(band_number)

    # 最初の区間は「受音点から最後の虚音源を見込む向き」（元コード 951〜969 行）。
    # 音は逆向きに進むが、経路の折れ線をたどるだけなので向きはこれでよい。
    vini = np.asarray(reciever_point, dtype=float)
    vray = sr.noramlized_soundray(images[ktmp] - vini)

    # ■バックトレースループ■ 元コード 948行 do k = ktmp, 0, -1
    for k in range(ktmp, -1, -1):

        # 基点から音線方向で最も手前に当たる面を探す（元コード 978〜1049 行）。
        # 音線追跡と同じ `FaceArrays` を使うので、全面との判定が 1 回の配列演算で済む
        hit_id_array, hit_distance_array, node_array = faces.nearest_hit(
            vini[None, :], vray[None, :])
        hit_id = int(hit_id_array[0])
        hit_distance = float(hit_distance_array[0])

        if hit_id >= 0:
            # ---- 壁面に当たる場合（元コード 1052 行）----
            if k == 0:
                # 音源より手前に壁があれば、直接音は遮蔽されている（元コード 1066〜1072 行）
                source_distance = np.linalg.norm(
                    np.asarray(soundsource_point, dtype=float) - vini)
                if hit_distance < source_distance:
                    return None
                break

            # 経路が言う壁と違う壁に当たったら却下（元コード 1086 行）
            if hit_id != wall_ids[k - 1]:
                return None

            node = node_array[0]

            # エネルギーを減衰させる（元コード 1091〜1094 行）。
            # 斜入射のエネルギー反射率 |R(θ)|^2 を掛け込む。式の根拠は sound_ray.energy_decay。
            absorption = np.atleast_1d(mesh[hit_id].absorption_coefficient)
            for b in range(band_number):
                energy[b] = sr.energy_decay(vray, mesh[hit_id].normal,
                                            absorption[b], energy[b])

            # 次の区間へ（元コード 1096〜1110 行）
            vini = sr.soundraycomesfrom_renew(node)
            vray = sr.noramlized_soundray(sr.soundray_renew(images[k - 1], vini))

        else:
            # ---- 壁面に当たらない場合（元コード 1115 行）----
            # k == 0 なら音源までさえぎるものが無いということなので受音成立。
            # 開いた形状（一面だけの壁など）ではこちらを通る。
            if k != 0:
                return None
            break

    # ■受音リストへの書き込み■ 元コード 1074〜1080 行
    vtgt = np.asarray(reciever_point, dtype=float) - images[ktmp]
    distance = float(np.linalg.norm(vtgt))
    if distance == 0.0:
        # 虚音源が受音点に一致（音源＝受音点）。時刻 0 で割れないので捨てる
        return None

    return {
        "reflection_count": ktmp,
        "time": distance / sound_velocity,
        "distance": distance,
        # 元コードの -vtgt。受音点から虚音源を見込む向き＝音が到来する方向
        "direction": (-vtgt) / distance,
        "energy": energy.copy(),
        "wall_ids": list(wall_ids),
    }


def loop(soundsource_point, reciever_point, reflectionmeshid_history, mesh,
         sound_velocity=SOUND_VELOCITY, band_number=None, filename=None,
         verbose=True):
    """非重複経路ループ。元コード 878 行 `do i = 1, countred`。

    引数:
        soundsource_point       : (3,)   音源座標
        reciever_point          : (3,)   受音点座標
        reflectionmeshid_history: list[list[int]]
            重複削除後の反射面 ID 履歴。`loop_deleteredundancy.delete()` の出力。
        mesh                    : list[Mesh] 室形状
        sound_velocity          : float  音速 [m/s]
        band_number             : int | None
            周波数バンド数。None なら mesh の吸音率の長さから決める
        filename                : str | None  パルス列を書き出す CSV のパス
        verbose                 : bool   進捗を表示するか

    戻り値:
        PulseList

    ※ 元コードは `absper`（面ごとの吸音率）を大域配列で参照していたので
      吸音率リストを引数で受け取る形になっていたが、Python 版では
      `Mesh.absorption_coefficient` が面ごとに持っているため引数から外した。
    """
    if band_number is None:
        if not mesh:
            raise ValueError("メッシュが空です")
        band_number = len(np.atleast_1d(mesh[0].absorption_coefficient))

    faces = mm.FaceArrays(mesh)
    records = []
    rejected = 0

    for wall_ids in reflectionmeshid_history:
        record = backtrace_path(soundsource_point, reciever_point, wall_ids, mesh,
                                band_number, sound_velocity, faces=faces)
        if record is None:
            rejected += 1
        else:
            records.append(record)

    pulses = PulseList.from_records(records, band_number, sound_velocity).sort_by_time()

    if verbose:
        print(f"[loop_noredundancy] 経路 {len(reflectionmeshid_history)} 本 → "
              f"受音 {len(pulses)} 本 / 却下 {rejected} 本")
        print(f"[loop_noredundancy] {pulses.summary()}")

    if filename is not None:
        pulses.save_csv(filename)
        if verbose:
            print(f"[loop_noredundancy] パルス列を書き出しました: {filename}")

    return pulses
