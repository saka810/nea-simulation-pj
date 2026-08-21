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
        walls            : (n,K) int   反射面の並び（-1 詰め）。**吸音率の差し替え用**
        cos_theta        : (n,K) float 各反射の |cosθ|（同上）
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
        # 反射の幾何。**吸音率を変えるだけなら再計算に使える**（`energy_from_geometry`）。
        # 古い CSV から読み戻した場合は空のまま（幾何は CSV に入らない）
        self.walls = np.zeros((0, 0), dtype=np.int32)
        self.cos_theta = np.zeros((0, 0))

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
        self.walls, self.cos_theta = _geometry_arrays(records)
        return self

    def has_geometry(self):
        """反射の幾何を持っているか（吸音率の差し替えができるか）。"""
        return len(self.walls) == len(self.time) and self.walls.size >= 0             and len(self.walls) == len(self.time)

    def recompute_energy(self, mesh, band_number=None):
        """**吸音率だけ差し替えてエネルギーを計算し直す。**

        音線追跡もバックトレースもやり直さない（`energy_from_geometry`）。
        材料条件表を変えて回すときの本体（F-9）。
        """
        absorption = np.array([np.atleast_1d(m.absorption_coefficient)
                               for m in mesh], dtype=float)
        if band_number is not None and absorption.shape[1] != band_number:
            raise ValueError(f"吸音率のバンド数 {absorption.shape[1]} が "
                             f"{band_number} と違います")
        self.band_number = absorption.shape[1]
        self.energy = energy_from_geometry(self.walls, self.cos_theta, absorption)
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
        if len(self.walls) == len(order):
            self.walls = self.walls[order]
            self.cos_theta = self.cos_theta[order]
        return self

    def save_csv(self, filename):
        """パルス列を CSV に保存する。

        列: reflection_count, time_s, distance_m, dir_x, dir_y, dir_z, energy_1..energy_b
        元コードの 11 列との対応は本クラスの docstring を参照。
        """
        # バンドは**周波数そのもの**を列名にする（`energy_1` だと 6 バンドと
        # 8 バンドで意味が変わって読み違える）。table.py の共通ルール
        import absorption as ab
        import table as tb
        bands = ab.octave_bands(self.band_number)
        header = ["reflection_count", "time_s", "distance_m", "dir_x", "dir_y", "dir_z"]
        header += [tb.band_column("energy", f) for f in bands]
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


def _geometry_arrays(records):
    """記録から (反射面 (n,K), cosθ (n,K)) を作る。K は最大反射回数。"""
    if not records or "cos_theta" not in records[0]:
        return np.zeros((len(records), 0), dtype=np.int32),                np.zeros((len(records), 0))
    kmax = max(len(r["wall_ids"]) for r in records)
    walls = np.full((len(records), kmax), -1, dtype=np.int32)
    cosines = np.zeros((len(records), kmax))
    for i, r in enumerate(records):
        n = len(r["wall_ids"])
        if n:
            walls[i, :n] = r["wall_ids"]
            cosines[i, :n] = r["cos_theta"]
    return walls, cosines


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


def backtrace_path(soundsource_point, receiver_point, wall_ids, mesh,
                   band_number, sound_velocity=SOUND_VELOCITY, faces=None):
    """経路 1 本分のバックトレース。元コード 948〜1132 行。

    引数:
        wall_ids : list[int] 反射面 ID の並び（先頭の番兵は含めないこと）
                   空リストなら直接音の経路。

    戻り値:
        受音に至れば dict、途中で却下されれば None。
    """
    if faces is None:
        faces = mm.collision_arrays(mesh)
    images = image_sources(soundsource_point, wall_ids, mesh)
    ktmp = len(wall_ids)

    energy = np.ones(band_number)

    # 最初の区間は「受音点から最後の虚音源を見込む向き」（元コード 951〜969 行）。
    # 音は逆向きに進むが、経路の折れ線をたどるだけなので向きはこれでよい。
    vini = np.asarray(receiver_point, dtype=float)
    vray = sr.normalized_soundray(images[ktmp] - vini)

    # 直前に当たった面。両面判定のとき同じ面に当たり直すのを防ぐ（音線追跡側と同じ理由）
    last_face = -1

    # ■バックトレースループ■ 元コード 948行 do k = ktmp, 0, -1
    for k in range(ktmp, -1, -1):

        # 基点から音線方向で最も手前に当たる面を探す（元コード 978〜1049 行）。
        # 音線追跡と**同じ入れ物**を使う（`mesh_method.collision_arrays`）。
        # 反射面の番号を突き合わせるので、ここが食い違うと経路が全部却下される
        hit_id_array, hit_distance_array, node_array = faces.nearest_hit(
            vini[None, :], vray[None, :],
            ignore=[last_face] if faces.two_sided else None)
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
            last_face = hit_id
            vini = sr.soundraycomesfrom_renew(node)
            vray = sr.normalized_soundray(sr.soundray_renew(images[k - 1], vini))

        else:
            # ---- 壁面に当たらない場合（元コード 1115 行）----
            # k == 0 なら音源までさえぎるものが無いということなので受音成立。
            # 開いた形状（一面だけの壁など）ではこちらを通る。
            if k != 0:
                return None
            break

    # ■受音リストへの書き込み■ 元コード 1074〜1080 行
    vtgt = np.asarray(receiver_point, dtype=float) - images[ktmp]
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


# ------------------------------------------------------------------------------
# バックトレースの配列演算化（F-5。2026-08-21 追加）
#
# `backtrace_path()` は経路 1 本ずつ Python のループで回る。1 経路につき反射回数ぶん
# `nearest_hit` を**単発で**呼ぶので、研修室（経路 37,592 本・最大反射 160 回）では
# 実測で全体の 65% を占めていた（音線追跡は 34%）。
#
# ここでは**全経路を同時に、反射回数 k を大きいほうから 1 段ずつ**下ろす。
# 経路の長さはまちまちだが、長さ L の経路は k = L のときに入ってくると考えれば、
# **どの時点でも動いている経路はすべて同じ k にいる**ので、素直に束で処理できる。
#
#   k = Kmax        L = Kmax の経路が入る
#   k = Kmax-1      L = Kmax-1 の経路が入り、上の経路は 1 段下りている
#   …
#   k = 0           直接音の区間。ここまで残った経路が受音成立
#
# 虚音源も同じく束で作る（面での鏡像は p - 2(n·p + d)n の 1 行で、経路ごとに面が
# 違うだけなので配列で引ける）。
#
# ★`backtrace_path()` は**参照実装として残す**。同じ結果になることを
#   tests/test_geosim.py の「バックトレースの配列演算（1 本ずつとの一致）」で確かめる。
# ------------------------------------------------------------------------------

def _image_source_table(soundsource_point, walls, lengths, faces):
    """経路ごとの虚音源の列 (P, Kmax+1, 3) をまとめて作る。

    `image_sources()` を束にしたもの。k 段目は「k 回反射ぶんの虚音源」で、
    長さ L の経路では k = 0..L だけが意味を持つ（それ以外は使わない）。
    """
    n_path, kmax = len(lengths), walls.shape[1]
    images = np.empty((n_path, kmax + 1, 3))
    images[:, 0, :] = np.asarray(soundsource_point, dtype=float)
    for k in range(kmax):
        active = np.nonzero(lengths > k)[0]
        if not len(active):
            break
        wall = walls[active, k]
        # 1 本ずつの版（image_sources）と同じく、念のため正規化してから使う
        normal = faces.normal[wall]
        normal = normal / np.linalg.norm(normal, axis=1)[:, None]
        offset = -np.einsum("ij,ij->i", normal, faces.plane_point[wall])
        signed = np.einsum("ij,ij->i", normal, images[active, k, :]) + offset
        images[active, k + 1, :] = images[active, k, :] - 2.0 * signed[:, None] * normal
    return images


def _energy_decay_batch(vray, normal, absorption, energy):
    """`sound_ray.energy_decay` を束で。式は 1 本ずつの版とまったく同じ順序で書く。

    vray (A,3) / normal (A,3) / absorption (A,b) / energy (A,b) → (A,b)
    """
    coefficient = np.abs(np.einsum("ij,ij->i", vray, normal))[:, None]
    root = np.sqrt(1.0 - absorption)
    reflection = (1.0 + root) * coefficient - (1.0 - root)
    reflection = reflection / ((1.0 + root) * coefficient + (1.0 - root))
    reflection = np.abs(reflection)
    return energy * reflection * reflection


def energy_from_geometry(walls, cosines, absorption, energy=None):
    """**反射の幾何（面と入射角）からエネルギーを計算し直す。**

    バックトレースで求まる経路の形（どの面に、どの角度で当たったか）は
    **吸音率に依らない**。だから幾何だけ残しておけば、吸音材を変えたときに
    音線追跡もバックトレースもやり直さず、ここだけ回せば済む（F-9）。

        E = Π_k |R(cosθ_k, α_{面_k})|²

    `R` は書籍 式(2.64) の斜入射反射係数で、`sound_ray.energy_decay` と同じ式。
    掛ける順番が逆（経路の先頭から）になるが、積なので結果は同じ（丸めのみ差）。

    引数:
        walls    (P, K) int   反射面のインデックス。**-1 は「反射なし」**の詰め物
        cosines  (P, K) float 各反射での |cosθ|
        absorption (面数, b)  面ごとの垂直入射吸音率
        energy   (P, b) | None  初期エネルギー（既定は 1）

    戻り値: (P, b)
    """
    walls = np.atleast_2d(np.asarray(walls))
    cosines = np.atleast_2d(np.asarray(cosines, dtype=float))
    absorption = np.asarray(absorption, dtype=float)
    n_path, kmax = walls.shape
    band_number = absorption.shape[1]
    result = (np.ones((n_path, band_number)) if energy is None
              else np.array(energy, dtype=float))

    for k in range(kmax):
        active = walls[:, k] >= 0
        if not np.any(active):
            continue
        wall = walls[active, k]
        coefficient = cosines[active, k][:, None]
        root = np.sqrt(1.0 - absorption[wall])
        reflection = ((1.0 + root) * coefficient - (1.0 - root))
        reflection = reflection / ((1.0 + root) * coefficient + (1.0 - root))
        reflection = np.abs(reflection)
        result[active] *= reflection * reflection
    return result


def backtrace_batch(soundsource_point, receiver_point, histories, mesh, faces,
                    band_number, sound_velocity=SOUND_VELOCITY):
    """経路の束をまとめてバックトレースする。戻り値は `backtrace_path` と同じ dict のリスト。

    却下された経路は入らない（1 本ずつの版が None を返すのと同じ扱い）。
    """
    n_path = len(histories)
    if n_path == 0:
        return []
    lengths = np.array([len(h) for h in histories], dtype=np.int64)
    kmax = int(lengths.max())
    walls = np.full((n_path, max(kmax, 1)), -1, dtype=np.int64)
    for i, h in enumerate(histories):
        if len(h):
            walls[i, :len(h)] = h

    images = _image_source_table(soundsource_point, walls, lengths, faces)

    receiver = np.asarray(receiver_point, dtype=float)
    source = np.asarray(soundsource_point, dtype=float)
    absorption = np.array([np.atleast_1d(m.absorption_coefficient) for m in mesh])

    origin = np.repeat(receiver[None, :], n_path, axis=0)
    ray = np.zeros((n_path, 3))
    energy = np.ones((n_path, band_number))
    # ★反射ごとの入射角の余弦を残す。**吸音率を変えるだけならここから
    #   エネルギーを計算し直せる**（`energy_from_geometry`）ので、
    #   条件を変えて回すときに音線追跡とバックトレースを省ける（F-9）
    cosines = np.zeros((n_path, max(kmax, 1)))
    last_face = np.full(n_path, -1, dtype=np.int64)
    alive = np.zeros(n_path, dtype=bool)        # すでに入ってきて、まだ却下されていない
    done = np.zeros(n_path, dtype=bool)         # 受音成立
    source_distance = np.linalg.norm(receiver - source)

    for k in range(kmax, -1, -1):
        entering = np.nonzero(lengths == k)[0]
        if len(entering):
            # 受音点から「最後の虚音源」を見込む向きで始める（元コード 951〜969 行）
            vector = images[entering, k, :] - receiver
            length = np.linalg.norm(vector, axis=1)
            ok = length > 0.0
            ray[entering[ok]] = vector[ok] / length[ok][:, None]
            alive[entering[ok]] = True
        index = np.nonzero(alive)[0]
        if not len(index):
            continue

        hit_id, hit_distance, node = faces.nearest_hit(
            origin[index], ray[index],
            ignore=last_face[index] if faces.two_sided else None)
        hit = hit_id >= 0

        if k == 0:
            # 音源より手前に壁があれば直接音は遮蔽（元コード 1066〜1072 行）。
            # 壁に当たらないなら遮るものが無いということで受音成立
            blocked = np.zeros(len(index), dtype=bool)
            if np.any(hit):
                reach = np.linalg.norm(origin[index] - source, axis=1)
                blocked = hit & (hit_distance < reach)
            done[index[~blocked]] = True
            alive[index] = False
            continue

        # 経路が言う壁と違う壁に当たったら却下（元コード 1086 行）
        want = walls[index, k - 1]
        good = hit & (hit_id == want)
        alive[index[~good]] = False
        keep = index[good]
        if not len(keep):
            continue
        local = np.nonzero(good)[0]

        wall = want[good]
        cosines[keep, k - 1] = np.abs(np.einsum("ij,ij->i", ray[keep],
                                                faces.normal[wall]))
        energy[keep] = _energy_decay_batch(ray[keep], faces.normal[wall],
                                           absorption[wall], energy[keep])
        last_face[keep] = wall
        origin[keep] = node[local]
        vector = images[keep, k - 1, :] - origin[keep]
        length = np.linalg.norm(vector, axis=1)
        alive[keep[length <= 0.0]] = False
        good_length = length > 0.0
        ray[keep[good_length]] = vector[good_length] / length[good_length][:, None]

    records = []
    for i in np.nonzero(done)[0]:
        vtgt = receiver - images[i, lengths[i], :]
        distance = float(np.linalg.norm(vtgt))
        if distance == 0.0:
            continue          # 虚音源が受音点に一致。時刻 0 で割れないので捨てる
        records.append({
            "reflection_count": int(lengths[i]),
            "time": distance / sound_velocity,
            "distance": distance,
            "direction": (-vtgt) / distance,
            "energy": energy[i].copy(),
            "wall_ids": list(histories[i]),
            # 吸音率を変えたときに再利用する材料（`energy_from_geometry`）
            "cos_theta": cosines[i, :lengths[i]].copy(),
        })
    return records


def loop(soundsource_point, receiver_point, reflectionmeshid_history, mesh,
         sound_velocity=SOUND_VELOCITY, band_number=None, filename=None,
         verbose=True, two_sided=False, progress=None, path_chunk=20000):
    """非重複経路ループ。元コード 878 行 `do i = 1, countred`。

    引数:
        soundsource_point       : (3,)   音源座標
        receiver_point          : (3,)   受音点座標
        reflectionmeshid_history: list[list[int]]
            重複削除後の反射面 ID 履歴。`loop_deleteredundancy.delete()` の出力。
        mesh                    : list[Mesh] 室形状
        sound_velocity          : float  音速 [m/s]
        band_number             : int | None
            周波数バンド数。None なら mesh の吸音率の長さから決める
        filename                : str | None  パルス列を書き出す CSV のパス
        verbose                 : bool   進捗を表示するか
        two_sided               : bool
            面の裏からの入射も当てるか。**音線追跡と必ず揃えること**。
            食い違うと、追跡側が通した経路をバックトレース側が全部却下する。
            詳細は mesh_method.FaceArrays
        path_chunk              : int
            一度に束で処理する経路の本数。虚音源の表 (経路 × 反射回数 × 3) の
            メモリ量を抑えるためのもので、結果には影響しない

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

    faces = mm.collision_arrays(mesh, two_sided=two_sided)
    records = []

    total = len(reflectionmeshid_history)
    # 経路を塊に分けて束で処理する（`backtrace_batch`）。1 本ずつの `backtrace_path` は
    # 参照実装として残してあり、結果が一致することをテストで確かめている。
    # 塊に分けるのは虚音源の表 (経路 × 反射回数 × 3) がメモリに乗る大きさに保つため
    step = max(1, int(path_chunk))
    for start in range(0, total, step):
        stop = min(start + step, total)
        records.extend(backtrace_batch(
            soundsource_point, receiver_point,
            reflectionmeshid_history[start:stop], mesh, faces,
            band_number, sound_velocity))
        if progress is not None:
            progress(stop / total)
    rejected = total - len(records)

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
