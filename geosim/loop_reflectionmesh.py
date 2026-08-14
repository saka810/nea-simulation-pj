import numpy as np

import sound_ray as sr
import mesh_method as mm
import receiver_sphere as rs


# 目的：音線を追跡し、受音に至った経路が「どの壁面をどの順に反射したか」を記録する
# 元コード backtrace.f90 524〜717 行
#
# ------------------------------------------------------------------------------
# 【traceff の扱い】2026-08-12 に元コードを読み解いて確定（12/04 の疑問1〜4への回答）
#
# 元コードは 2 つの配列を使い分けている。
#   tractmp(0:nref+1)   … 一時反射履歴。音線 1 本ごとに初期化し、反射のたびに壁 ID を追記
#   traceff(count,0:nref) … 受音した瞬間のスナップショットを貯める 2 次元配列
#
# 鍵は「受音判定（663行）が壁 ID の追記（677行）より前にある」こと。
# 受音した時点の履歴は "そこまでに済んだ反射" であるべきで、これから当たる壁は含めない。
# したがって受音と衝突は独立した条件で、
#   ・壁 ID の追記      → 衝突するたびに無条件
#   ・スナップショット保存 → 受音するたび
#
#   [疑問1] do j のループを k ループの外に出せるか
#       → 出せない。1 本の音線が別々の k で複数回受音しうるので、受音の瞬間ごとに要る。
#   [疑問2] nray*nref*nref 回まわるのか
#       → まわらない。do j = 0, nref のコピーは受音時のみ（count 回）。
#         固定長配列の行まるごとコピーなだけで、Python なら list() のコピー 1 回で等価。
#   [疑問3] do j = 0, nref を do j = 0, k にしてよいか
#       → よい。k+1 以降は必ず 0。0 埋めが要るのは Fortran が固定長配列＋traceffn（有効反射回数）
#         で行長を管理しているため。Python の可変長リストなら 0 埋め自体が不要。
#   [疑問4] inside かつ collision を条件にしてよいか
#       → よくない。両者は独立。上記のとおり追記条件が別。
# ------------------------------------------------------------------------------
def loop(soundsource_point, reciever_point, soundray_list, nref, mesh, sphere_radius,
         recorder=None, ray_chunk=20000):
    """音線追跡。

    引数:
        soundsource_point : (3,)     音源座標
        reciever_point    : (3,)     受音点座標
        soundray_list     : (nray,3) 音源から出る音線の単位ベクトル群
        nref              : int      最大反射回数
        mesh              : list[Mesh] 室形状
        sphere_radius     : float    受音球の半径
        recorder          : RayRecorder | None
            渡すと可視化用の軌跡を副チャンネルで記録する（本線の結果には影響しない）。
            詳細は ray_recorder.py と docs/出力・可視化方針.md を参照。
        ray_chunk         : int  一度にまとめて追跡する音線の本数（メモリと速度の兼ね合い）

    【2026-08-14 高速化（TODO F-1）】
    元コードは「音線 1 本ずつ × 面 1 枚ずつ」の二重ループだったが、Python では
    ループ 1 回あたりの手間が大きすぎて実用速度に届かなかった。
    そこで **音線を束にして、全面との交差判定を 1 回の配列演算で行う**形に書き直した。
    処理の順序（受音判定 → 打ち切り → 壁 ID 追記 → 反射）は元コードのまま。
    計算式・比較・同点の扱いも同じなので、結果はビット単位で一致する。

    戻り値:
        list[list[int]] 受音した経路ごとの反射面 ID 履歴（元コード traceff 相当）
        反射した順に面 ID が並ぶ。**反射回数 = len(history)**。
        直接音（無反射）で受音した経路は空リストになる。

        【番兵 -1 について】2026-08-14 に廃止（TODO A-9）。
        元コードは `tractmp(0) = -1`、`tractmp(k+1) = jtmp` と 1 つずらして詰めており、
        バックトレース側は `tracred(i,k)` を **k = 1 から** 読む（元コード 903行）。
        つまり添字 0 の -1 は**一度も参照されない固定長配列のパディング**でしかない。
        Python は可変長リストなので持つ意味がなく、持っていると
        `mesh[history[0]]` が -1 で最後の面を指す事故のもとになるため落とした。
    """
    faces = mm.FaceArrays(mesh)
    soundray_list = np.atleast_2d(np.asarray(soundray_list, dtype=float))
    results = []

    # 音線を塊に分けて処理する。1 塊のなかは配列演算なので Python のループが回らない
    step = max(1, int(ray_chunk))
    for start in range(0, len(soundray_list), step):
        stop = min(start + step, len(soundray_list))
        _trace_chunk(soundsource_point, reciever_point,
                     soundray_list[start:stop], np.arange(start, stop),
                     nref, mesh, faces, sphere_radius, results, recorder)
    return results


def _trace_chunk(soundsource_point, reciever_point, soundray_list, ray_ids,
                 nref, mesh, faces, sphere_radius, results, recorder):
    """音線の塊 1 つぶんを追跡する。処理の順序は元コードと同じ。"""
    n_ray = len(soundray_list)

    origins = np.repeat(np.asarray(soundsource_point, dtype=float)[None, :], n_ray, axis=0)
    directions = soundray_list / np.linalg.norm(soundray_list, axis=1)[:, None]

    # 反射履歴。元コードの固定長配列 tractmp と同じ持ち方にした。
    # 可変長リストを音線ごとに持つより、行をスライスするだけで済むので速い
    history = np.full((n_ray, nref + 1), -1, dtype=np.int64)
    depth = np.zeros(n_ray, dtype=np.int64)      # 各音線がこれまでに反射した回数

    # まだ追跡中の音線（元コードで exit していないもの）の添字
    alive = np.arange(n_ray)

    # 記録対象の音線だけを相手にする。
    # ここを全音線で回すと、せっかく配列演算にしても Python のループが復活してしまう
    recorded = np.zeros(n_ray, dtype=bool)
    if recorder is not None:
        recorded = np.array([recorder.is_recording(int(r)) for r in ray_ids])
        for i in np.nonzero(recorded)[0]:
            recorder.start_ray(int(ray_ids[i]), directions[i], origins[i])
    termination = np.full(n_ray, "nref", dtype=object)

    # ■反射回数ループ■ 元コード 545行 do k = 0, nref
    # k = 0 は直接音の区間。k = nref まで回すので反射回数は 0〜nref 回になる。
    for k in range(nref + 1):
        if len(alive) == 0:
            break

        origin_alive = origins[alive]
        direction_alive = directions[alive]

        # ■壁面ループ■ 元コード 576行 do j = 1, sfcount をまとめて処理
        # 向き判定・t>0 判定・交点算出・三角形内部判定・距離算出・最寄り面の決定を
        # 1 回の配列演算で行う（元コード 579〜636行）
        hit_id, min_distance, node = faces.nearest_hit(origin_alive, direction_alive)

        # 受音判定（元コード 649〜663行）
        # ★ 壁 ID の追記より前に行う。この時点の履歴が「そこまでに済んだ反射」になる
        inside = rs.inside_sphere_batch(sphere_radius, direction_alive, origin_alive,
                                        reciever_point, min_distance)
        for local in np.nonzero(inside)[0]:
            ray = alive[local]
            # 受音経路として保存（元コード 666〜669行 traceff への書き込み）
            results.append(history[ray, :depth[ray]].tolist())
            if recorded[ray]:
                recorder.mark_receive(k, ray_index=int(ray_ids[ray]))

        # 当たる壁がなければその音線は終了（元コード 709〜712行）
        collision = hit_id >= 0
        termination[alive[~collision]] = "no_hit"

        alive_next = alive[collision]
        if len(alive_next) == 0:
            alive = alive_next
            break

        hit_id = hit_id[collision]
        node = node[collision]
        direction_hit = direction_alive[collision]

        # 一時反射履歴に壁面番号を書き込み（元コード 677行）
        # ★ 受音したかどうかに関わらず、衝突のたびに無条件で追記する
        history[alive_next, depth[alive_next]] = hit_id
        depth[alive_next] += 1

        if recorder is not None:
            for local in np.nonzero(recorded[alive_next])[0]:
                ray = alive_next[local]
                recorder.add_reflection(node[local], int(hit_id[local]),
                                        direction_hit[local], mesh[hit_id[local]],
                                        ray_index=int(ray_ids[ray]))

        # 次回の基点 = 最寄り面との交点（元コード 680〜694行）
        origins[alive_next] = node

        # 反射音線のベクトル（元コード 697〜706行）
        # v' = v - 2(v·n)n を計算して正規化する
        normal_hit = faces.normal[hit_id]
        reflected = direction_hit - 2.0 * np.einsum(
            "ij,ij->i", direction_hit, normal_hit)[:, None] * normal_hit
        directions[alive_next] = reflected / np.linalg.norm(reflected, axis=1)[:, None]

        alive = alive_next

    if recorder is not None:
        for ray in np.nonzero(recorded)[0]:
            recorder.end_ray(termination[ray], ray_index=int(ray_ids[ray]))
