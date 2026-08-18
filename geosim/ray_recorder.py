"""可視化用の音線軌跡レコーダ（出力①②のためのデータ）。

**本線の計算とは目的が違うので、別チャンネルにしてある。**
音線追跡（`loop_reflectionmesh`）が残すのは受音に至った経路の反射面 ID 列だが、
可視化で見たいのは「音がどう広がるか」なので、**受音しなかった音線も含めた軌跡**が要る。

`recorder=None` を渡せば何も記録しない。本線の速度には影響しない。
間引き（`max_rays`）も等間隔なので、密度の見え方が偏らない。

役割と保存形式の詳細は下のコメントに書いてある。
"""


import numpy as np

import sound_ray as sr


# 可視化用の音線軌跡レコーダ
#
# 【役割】
# loop_reflectionmesh.loop() の「本線」の出力は受音経路の反射面ID列（元コード traceff）だが、
# 可視化で見たいのは「音がどう広がるか」なので、受音しなかった音線も含めた軌跡が要る。
# 目的が違うので同じ配列に混ぜず、このレコーダを副チャンネルとして分けている。
#
# 【この軌跡から作れる出力】（docs/出力・可視化方針.md 参照）
#   ① 音線の可視化      … nodes をそのまま折れ線として描く
#   ② 音粒子の可視化（動画）… times を使って任意時刻の粒子位置を線形補間する
#   ⑤ 伝搬方向確認      … 受音イベント時の入射方向（本命は虚音源バックトレース側の出力）
#
# 【間引きについて】
# 100 万本規模では全音線の記録はメモリに乗らないので間引く。
# 既定は「絶対本数の上限（max_rays）＋ 等間隔ストライド」。理由は docs/出力・可視化方針.md に記載。

SOUND_VELOCITY_DEFAULT = 343.0
MAX_RAYS_DEFAULT = 2000


class RayTrajectory:
    """間引いて記録した音線 1 本分の軌跡。

    属性（finalize() 後は ndarray）:
        ray_index    int         間引き前の元の音線番号
        direction    (3,)        音源から出たときの初期方向ベクトル
        nodes        (m+1, 3)    節点座標。先頭は音源位置、以降は各反射点
        distances    (m+1,)      音源からの累積距離 [m]
        mesh_ids     (m,)        各セグメントの終端で反射した面 ID
        energies     (m+1, nb)   各節点でのバンド別エネルギー（音源で 1.0）。無効時は None
        receive_steps(r,)        受音したときのセグメント番号 k
        termination  str         'no_hit'（当たる壁がなくなった）/ 'nref'（反射回数上限）
    """

    __slots__ = ("ray_index", "direction", "nodes", "distances", "mesh_ids",
                 "energies", "receive_steps", "termination")

    def __init__(self, ray_index, direction, origin, initial_energy):
        self.ray_index = int(ray_index)
        self.direction = np.array(direction, dtype=float)
        self.nodes = [np.array(origin, dtype=float)]
        self.distances = [0.0]
        self.mesh_ids = []
        self.energies = None if initial_energy is None else [np.array(initial_energy, dtype=float)]
        self.receive_steps = []
        self.termination = None

    def times(self, sound_velocity=SOUND_VELOCITY_DEFAULT):
        """各節点への到達時刻 [s]。音粒子アニメーションの時間軸に使う。"""
        return np.asarray(self.distances, dtype=float) / sound_velocity

    def position_at(self, t, sound_velocity=SOUND_VELOCITY_DEFAULT):
        """時刻 t [s] における音粒子の位置。まだ出発前／既に消滅後は None。

        音粒子アニメーション（出力②）はこれを全軌跡について毎フレーム呼べばよい。
        """
        d = t * sound_velocity
        dist = np.asarray(self.distances, dtype=float)
        if d < 0.0 or d > dist[-1]:
            return None
        k = int(np.searchsorted(dist, d, side="right")) - 1
        k = min(max(k, 0), len(dist) - 2) if len(dist) >= 2 else 0
        if len(dist) < 2:
            return np.asarray(self.nodes[0], dtype=float)
        span = dist[k + 1] - dist[k]
        w = 0.0 if span <= 0.0 else (d - dist[k]) / span
        n0 = np.asarray(self.nodes[k], dtype=float)
        n1 = np.asarray(self.nodes[k + 1], dtype=float)
        return n0 + w * (n1 - n0)

    def finalize(self):
        self.nodes = np.asarray(self.nodes, dtype=float)
        self.distances = np.asarray(self.distances, dtype=float)
        self.mesh_ids = np.asarray(self.mesh_ids, dtype=np.int32)
        self.receive_steps = np.asarray(self.receive_steps, dtype=np.int32)
        if self.energies is not None:
            self.energies = np.asarray(self.energies, dtype=float)
        return self


class RayRecorder:
    """音線追跡ループに渡すと、間引いた音線の軌跡を記録する。

    使い方:
        recorder = RayRecorder(total_rays=soundray_number)
        history = lr.loop(..., recorder=recorder)
        recorder.save_npz("rays.npz")
    """

    def __init__(self, total_rays, max_rays=MAX_RAYS_DEFAULT,
                 sound_velocity=SOUND_VELOCITY_DEFAULT,
                 band_number=8, record_energy=True):
        """
        total_rays     : 飛ばす音線の総数
        max_rays       : 記録する音線の本数の上限。None なら全記録
        sound_velocity : 音速 [m/s]。時刻の算出に使う
        band_number    : エネルギーの周波数バンド数
        record_energy  : 反射ごとのエネルギー減衰も記録するか
        """
        self.total_rays = int(total_rays)
        self.max_rays = max_rays
        self.sound_velocity = float(sound_velocity)
        self.band_number = int(band_number)
        self.record_energy = bool(record_energy)

        # 等間隔ストライドで間引く。max_rays=None なら全記録（stride=1）
        if max_rays is None or max_rays <= 0 or self.total_rays <= max_rays:
            self.stride = 1
        else:
            self.stride = max(1, self.total_rays // int(max_rays))

        self.trajectories = []
        # 追跡中の軌跡を音線番号で持つ。
        # 2026-08-14: 音線追跡を束処理にしたので、音線が順番に完結しなくなった。
        # 「今の 1 本」だけを持つ形（_current）から辞書に変えてある。
        self._open = {}
        self._current_index = None

    # ---- 記録対象かの判定 ----------------------------------------------

    def is_recording(self, ray_index):
        return ray_index % self.stride == 0

    # ---- 音線追跡ループから呼ばれるフック --------------------------------
    #
    # ray_index を省略すると「直前に start_ray した音線」に対して働く（逐次処理向け）。
    # 束処理では ray_index を明示して呼ぶ。

    def start_ray(self, ray_index, direction, origin):
        """音線 1 本の追跡開始時に呼ぶ。"""
        ray_index = int(ray_index)
        self._current_index = ray_index
        if not self.is_recording(ray_index):
            return
        initial_energy = np.ones(self.band_number) if self.record_energy else None
        self._open[ray_index] = RayTrajectory(ray_index, direction, origin,
                                              initial_energy)

    def _trajectory(self, ray_index):
        key = self._current_index if ray_index is None else int(ray_index)
        return self._open.get(key)

    def add_reflection(self, node, mesh_id, incident_ray, mesh_object, ray_index=None):
        """反射が確定したときに呼ぶ。node は反射点（＝交点）。"""
        traj = self._trajectory(ray_index)
        if traj is None:
            return
        node = np.asarray(node, dtype=float)
        step = np.linalg.norm(node - traj.nodes[-1], ord=2)
        traj.nodes.append(node)
        traj.distances.append(traj.distances[-1] + float(step))
        traj.mesh_ids.append(int(mesh_id))

        if traj.energies is not None:
            absorption = getattr(mesh_object, "absorption_coefficient", None)
            if absorption is None:
                traj.energies.append(np.array(traj.energies[-1], dtype=float))
            else:
                decayed = sr.energy_decay(incident_ray, mesh_object.normal,
                                          absorption, traj.energies[-1])
                traj.energies.append(np.broadcast_to(
                    np.asarray(decayed, dtype=float), (self.band_number,)).copy())

    def mark_receive(self, step, ray_index=None):
        """受音球を通過したときに呼ぶ。step は反射回数ループの k。"""
        traj = self._trajectory(ray_index)
        if traj is not None:
            traj.receive_steps.append(int(step))

    def end_ray(self, termination, ray_index=None):
        """音線 1 本の追跡終了時に呼ぶ。termination は 'no_hit' か 'nref'。"""
        key = self._current_index if ray_index is None else int(ray_index)
        traj = self._open.pop(key, None)
        if traj is None:
            return
        traj.termination = termination
        self.trajectories.append(traj.finalize())

    # ---- 保存・要約 ------------------------------------------------------

    def summary(self):
        counts = [len(t.mesh_ids) for t in self.trajectories]
        return {
            "total_rays": self.total_rays,
            "recorded_rays": len(self.trajectories),
            "stride": self.stride,
            "sound_velocity": self.sound_velocity,
            "max_reflection": max(counts) if counts else 0,
            "max_time": (max(float(t.distances[-1]) for t in self.trajectories)
                         / self.sound_velocity) if self.trajectories else 0.0,
        }

    def save_npz(self, filename):
        """音線ごとに長さが違う可変長データなので、連結＋オフセット方式で保存する。

        読み出し側は offsets[i]:offsets[i+1] でスライスすれば i 本目の軌跡が取れる。
        """
        if not self.trajectories:
            raise ValueError("記録された軌跡がありません")

        node_offsets = np.zeros(len(self.trajectories) + 1, dtype=np.int64)
        mesh_offsets = np.zeros(len(self.trajectories) + 1, dtype=np.int64)
        recv_offsets = np.zeros(len(self.trajectories) + 1, dtype=np.int64)
        for i, t in enumerate(self.trajectories):
            node_offsets[i + 1] = node_offsets[i] + len(t.nodes)
            mesh_offsets[i + 1] = mesh_offsets[i] + len(t.mesh_ids)
            recv_offsets[i + 1] = recv_offsets[i] + len(t.receive_steps)

        payload = {
            "nodes": np.concatenate([t.nodes for t in self.trajectories]),
            "distances": np.concatenate([t.distances for t in self.trajectories]),
            "mesh_ids": np.concatenate([t.mesh_ids for t in self.trajectories])
                        if mesh_offsets[-1] else np.zeros(0, dtype=np.int32),
            "receive_steps": np.concatenate([t.receive_steps for t in self.trajectories])
                             if recv_offsets[-1] else np.zeros(0, dtype=np.int32),
            "node_offsets": node_offsets,
            "mesh_offsets": mesh_offsets,
            "recv_offsets": recv_offsets,
            "ray_indexes": np.array([t.ray_index for t in self.trajectories], dtype=np.int64),
            "directions": np.array([t.direction for t in self.trajectories]),
            "terminations": np.array([t.termination for t in self.trajectories]),
            "total_rays": np.array(self.total_rays),
            "stride": np.array(self.stride),
            "sound_velocity": np.array(self.sound_velocity),
        }
        if self.trajectories[0].energies is not None:
            payload["energies"] = np.concatenate([t.energies for t in self.trajectories])

        np.savez_compressed(filename, **payload)
