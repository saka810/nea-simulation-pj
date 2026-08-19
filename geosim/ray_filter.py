"""音線の絞り込み（注目したい音線・音粒子だけを選び出す）。

**ねらい**：可視化で「この経路をもっと見たい」「この辺りの粒子を見たい」
「この方向に飛ぶ音線に注目したい」を叶えるための道具（ユーザー要望 2026-08-19）。

音線を 60 本だけ描いても、室内では線が重なって**どれがどれだか追えない**。
本数を増やすともっと読めなくなる。そこで**数を減らす方向ではなく、
見たいものだけを残す方向**で絞る。

ここは**純粋な計算だけ**（画面には触らない）。
`RayLog` と条件を渡すと、条件に合う音線の添字を返す。
`view_rays` 側はその添字を `RayDisplay` / `ParticleAnimation` に渡し直すだけで済む
（どちらも「表示する添字」を差し替えられる作りになっている）。

絞り込みは 3 通り。**どれも「基準点をクリックで決める」**のが入口。

| 絞り込み | 関数 | 何を残すか |
|---|---|---|
| 近く | `near_point()` | 折れ線が基準点の近くを通った音線 |
| 方向 | `in_direction()` | **音源から出たときの向き**が基準の向きに近い音線 |
| 1 本 | `nearest_ray()` | 基準点にいちばん近い 1 本だけ |

おまけで「この面で反射した経路」（`through_face()`）も引ける。

★**方向は「出射方向」で見る**（`RayLog.directions`）。途中の向きで見ると、
  反射のたびに向きが変わるのでどの音線も条件に引っかかってしまい、絞り込みにならない。
"""

import numpy as np


def _pool(raylog, index):
    """対象の添字。`None` なら全部。"""
    if index is None:
        return np.arange(raylog.ray_count)
    return np.asarray(index, dtype=int)


def source_point(raylog):
    """音源の位置。**どの音線も 1 点目が音源**なのでそこから取る。"""
    return np.asarray(raylog.pad_nodes[0, 0], dtype=float)


def ray_distances(raylog, point, index=None, max_reflection=None):
    """各音線の**折れ線**から `point` までの最短距離 [m]。

    節点までの距離ではなく**線分までの距離**を測る。
    節点だけで見ると、長い区間の途中を通っている音線を拾い落とす。

    音線 × 区間をまとめて配列で計算する（`RayLog.pad_nodes` が
    行ごとに揃えた形で持っているのでそのまま使える）。

    ★`max_reflection` は**何回目の反射までを見るか**。
      画面に描いている範囲（`RayDisplay.max_reflection`）と合わせるためにある。
      反射を 50 回も追うと、どの音線も室内を回って基準点の近くを通ってしまい、
      「近く」で絞ったつもりが絞れていない状態になる
      （研修室で半径 0.5 m に 2000 本中 741 本が該当した）。
      初期反射だけに限ると本当に近くを通ったものだけが残る。
    """
    pool = _pool(raylog, index)
    point = np.asarray(point, dtype=float)

    start = raylog.pad_nodes[pool, :-1, :]          # (m, width-1, 3)
    end = raylog.pad_nodes[pool, 1:, :]
    segment = end - start
    length2 = np.einsum("ijk,ijk->ij", segment, segment)

    # 区間上のいちばん近い点。長さ 0 の区間（詰め物）は始点そのもの
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(length2 > 0.0,
                     np.einsum("ijk,ijk->ij", point[None, None, :] - start,
                               segment) / np.where(length2 > 0.0, length2, 1.0),
                     0.0)
    t = np.clip(t, 0.0, 1.0)
    nearest = start + t[:, :, None] * segment
    distance = np.linalg.norm(nearest - point[None, None, :], axis=2)

    # 詰め物の区間（節点数を超えた分）は見ない。
    # 最後の節点が延々と並んでいるので、そのままだと「端点までの距離」を
    # 何度も数えるだけだが、区間数が 0 の音線で min が取れなくなるのを避ける
    counts = raylog.node_counts[pool]
    limit = counts - 1
    if max_reflection is not None:
        limit = np.minimum(limit, max(1, int(max_reflection)))
    valid = np.arange(start.shape[1])[None, :] < limit[:, None]
    distance = np.where(valid, distance, np.inf)
    # 節点が 1 つしかない音線（区間なし）は、その点までの距離にする
    lonely = counts < 2
    if np.any(lonely):
        only = np.linalg.norm(raylog.pad_nodes[pool[lonely], 0, :] - point, axis=1)
        distance[lonely, 0] = only
    return distance.min(axis=1)


def near_point(raylog, point, radius, index=None, max_reflection=None):
    """折れ線が `point` から `radius` [m] 以内を通った音線の添字。"""
    pool = _pool(raylog, index)
    distance = ray_distances(raylog, point, pool, max_reflection)
    return pool[distance <= float(radius)]


def nearest_ray(raylog, point, index=None, max_reflection=None):
    """`point` にいちばん近い音線の添字（1 つ）。対象が空なら None。"""
    pool = _pool(raylog, index)
    if len(pool) == 0:
        return None
    return int(pool[np.argmin(ray_distances(raylog, point, pool, max_reflection))])


def launch_angles(raylog, direction, index=None):
    """各音線の**出射方向**と `direction` のなす角 [度]。

    途中の向きではなく音源から出たときの向きで見る（このファイルの冒頭参照）。
    """
    pool = _pool(raylog, index)
    reference = np.asarray(direction, dtype=float)
    norm = np.linalg.norm(reference)
    if norm == 0.0:
        raise ValueError("方向ベクトルの長さが 0 です")
    reference = reference / norm

    launched = raylog.directions[pool]
    launched = launched / np.linalg.norm(launched, axis=1)[:, None]
    cosine = np.clip(launched @ reference, -1.0, 1.0)
    return np.rad2deg(np.arccos(cosine))


def in_direction(raylog, direction, half_angle, index=None):
    """出射方向が `direction` から `half_angle` [度] 以内の音線の添字。

    音源を頂点とする円錐の中に入った音線だけが残る。
    """
    pool = _pool(raylog, index)
    return pool[launch_angles(raylog, direction, pool) <= float(half_angle)]


def direction_to(raylog, point):
    """音源から `point` へ向かう単位ベクトル。クリックした先を「方向」にするため。"""
    vector = np.asarray(point, dtype=float) - source_point(raylog)
    norm = np.linalg.norm(vector)
    if norm == 0.0:
        raise ValueError("音源と同じ位置なので方向が決まりません")
    return vector / norm


def through_face(raylog, face_id, index=None):
    """指定した面で反射した音線の添字。「この壁を経由した経路」を見るのに使う。"""
    pool = _pool(raylog, index)
    keep = []
    for i in pool:
        start, stop = raylog.mesh_offsets[i], raylog.mesh_offsets[i + 1]
        if int(face_id) in raylog.mesh_ids[start:stop]:
            keep.append(i)
    return np.asarray(keep, dtype=int)


def describe_ray(raylog, i):
    """音線 1 本の中身を文字にする（1 本だけ表示したときに何が見えているか出す）。"""
    i = int(i)
    start, stop = raylog.mesh_offsets[i], raylog.mesh_offsets[i + 1]
    faces = [int(v) for v in raylog.mesh_ids[start:stop]]
    total = float(raylog.total_distance[i])
    received = bool(raylog.received[i])
    lines = [f"音線 {int(raylog.ray_indexes[i])}（{len(faces)} 回反射）",
             f"経路長 {total:.2f} m / 到達 {total / raylog.sound_velocity * 1000:.1f} ms",
             f"受音 {'した' if received else 'しなかった'}"
             f" / 打ち切り {raylog.terminations[i]}"]
    if faces:
        shown = ", ".join(str(f) for f in faces[:8])
        lines.append(f"反射面 {shown}" + (" …" if len(faces) > 8 else ""))
    return "\n".join(lines)
