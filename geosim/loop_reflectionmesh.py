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
         recorder=None):
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

    戻り値:
        list[list[int]] 受音した経路ごとの反射面 ID 履歴（元コード traceff 相当）

        ⚠ 各履歴の先頭要素は元コード tractmp(0) = -1 に対応する番兵で、面 ID ではない。
          反射回数 = len(history) - 1。下流で mesh[history[k]] のように使うと
          -1 が「最後の面」を指してしまうので、必ず先頭を除いてから扱うこと。
    """
    reflectionmeshid_history_2dim = []

    # ■音線ループ■ 元コード 527行 do i = 1, nray
    for soundray_i in range(len(soundray_list)):

        # 一時反射履歴の初期化。音線 1 本ごとにリセットする（元コード 535〜542行）
        # 先頭の -1 は元コード tractmp(0) = -1 に対応する番兵
        reflectionmeshid_history = [-1]

        soundray_comesfrom = np.asarray(soundsource_point, dtype=float)
        sound_ray = soundray_list[soundray_i, :]

        if recorder is not None:
            recorder.start_ray(soundray_i, sound_ray, soundray_comesfrom)

        termination = "nref"

        # ■反射回数ループ■ 元コード 545行 do k = 0, nref
        # k = 0 は直接音の区間。k = nref まで回すので反射回数は 0〜nref 回になる。
        for k in range(nref + 1):

            sound_ray = sr.noramlized_soundray(sound_ray)

            # 「壁面に当たる」フラグと最寄り距離の初期化（元コード 570, 573行）
            # ※ 壁面ループの中ではなく外で初期化する。中で初期化すると
            #    最後に調べた面の結果しか残らない
            collision = False
            mesh_nearestid = -1
            min_distance = np.inf

            # ■壁面ループ■ 元コード 576行 do j = 1, sfcount
            for j in range(len(mesh)):
                # 向き判定・t>0 判定・交点算出・三角形内部判定・距離算出をまとめて行う
                # （元コード 579〜630行）
                # ※ ここでは soundray_comesfrom を書き換えない。元コードも交点は一時変数 node に
                #    置くだけで、基点の更新は最寄り面が確定してから（692〜694行）行う
                collision_j, distance_j = mm.collision_distance(
                    sound_ray, soundray_comesfrom, mesh[j].normal, mesh[j].vertexes)

                if collision_j:
                    collision = True
                    # 最も手前の交点を持つ面を採用（元コード 633〜636行）
                    if distance_j < min_distance:
                        min_distance = distance_j
                        mesh_nearestid = j

            # 受音判定（元コード 649〜663行）
            # ★ 壁 ID の追記より前に行う。この時点の履歴が「そこまでに済んだ反射」になる
            inside = rs.inside_sphere(sphere_radius, sound_ray, soundray_comesfrom,
                                      reciever_point, min_distance)

            if inside:
                # 受音経路として保存（元コード 666〜669行 traceff への書き込み）
                # ★ list() でコピーする。同じリストを参照のまま append すると
                #    以降の反射で全要素が書き換わってしまう
                reflectionmeshid_history_2dim.append(list(reflectionmeshid_history))
                if recorder is not None:
                    recorder.mark_receive(k)

            # 当たる壁がなければこの音線は終了（元コード 709〜712行）
            if not collision:
                termination = "no_hit"
                break

            # 一時反射履歴に壁面番号を書き込み（元コード 677行）
            # ★ 受音したかどうかに関わらず、衝突のたびに無条件で追記する
            reflectionmeshid_history.append(mesh_nearestid)

            # 次回の基点 = 最寄り面との交点（元コード 680〜694行）
            t = mm.parameter_t(sound_ray, soundray_comesfrom,
                               mesh[mesh_nearestid].normal, mesh[mesh_nearestid].vertexes)
            node = mm.node_renew(sound_ray, soundray_comesfrom, t)

            if recorder is not None:
                recorder.add_reflection(node, mesh_nearestid, sound_ray, mesh[mesh_nearestid])

            soundray_comesfrom = sr.soundraycomesfrom_renew(node)

            # 反射音線のベクトル（元コード 697〜706行）
            sound_ray = sr.reflection_generator(sound_ray, mesh[mesh_nearestid].normal)

        if recorder is not None:
            recorder.end_ray(termination)

    return reflectionmeshid_history_2dim
