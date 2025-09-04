import numpy as np
import sound_ray as sr
import mesh_method as mm
from mesh import Mesh
import receiver_sphere as rs


# 9月 打ち合わせ用




# 受音球の半径
# sphere_radius = 0.1
# 音線数
# soundray_number = 100000
# 反射回数
# nref = 10000
# soundray_list = np.zeros((soundray_number, 3))
# soundsource_point = np.zeros(3)
# reciever_point = np.zeros(3)
# mesh_count = 1000  # メッシュの個数（ダミー）
# mesh = []  # メッシュの個数分Meshオブジェクトを準備
# for i in range(mesh_count):
#     mesh.append(Mesh(np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3), "name", 0.1))


# 目的：音線を追跡した際、反射する面のインデックスを記録する
# 音源位置、受音点、音源から出る音線全て持つ配列、反射回数、室形状のメッシュ全ての情報、受音球の半径
# を用いて、
# 面に反射する音線を追跡し、反射した面のインデックス（番号）を
# 反射履歴　reflectionmeshid_history として記録する
# このreflectionmeshid_historyを出力として返す
# 元コード524
def loop(soundsource_point, reciever_point, soundray_list, nref, mesh, sphere_radius):
    reflectionmeshid_history = []
    reflectionmeshid_history.append(-1)  # 0番目は-1を入れてる 後で追加するときにインデックスをずらす必要があるみたい

    for soundray_i in range(len(soundray_list)):

        soundray_comesfrom = soundsource_point

        for k in range(nref):
            sound_ray = soundray_list[soundray_i, :]
            sound_ray = sr.noramlized_soundray(sound_ray)

            # min_distance = 1000000000.0  # とりあえず大きい数字
            min_distance = np.inf  # とりあえず大きい数字

            # for j in range(mesh_count):
            for j in range(len(mesh)):
                collision = False

                if np.dot(sound_ray, mesh[j].normal) < 0:
                    t = mm.parameter_t(sound_ray, soundray_comesfrom, mesh[j].normal, mesh[j].vertexes)

                    if t > 0:
                        soundray_comesfrom = mm.node_renew(sound_ray, soundray_comesfrom, t)
                        collision, distance = mm.collision_distance(sound_ray, soundray_comesfrom, mesh[j].normal,
                                                                    mesh[j].vertexes)

                        if collision:
                            if distance < min_distance:
                                min_distance = distance
                                mesh_nearestid = j

            # ここから元コード649  音線基点・受音点間ベクトル 受音球に入っているかの判定
            inside = rs.inside_sphere(sphere_radius, sound_ray, soundray_comesfrom, reciever_point, min_distance)

            # if inside:
            # 元コードのこれはtractmp(0か-1)で塗りつぶし？なぜ？ 反射履歴に書き込み if inside and collisionの後＊の位置ではないのか
            # count = count + 1
            # do j = 0, nref
            #     traceff(count, j) = tractmp(j)
            # end  do


            # 662行目以降　再度打ち合わせ
            # 662~671は if inside内で行うもの ????
            if inside:
                # 元コード↓
                # ! 一時反射履歴に壁面番号を書き込み
                # tractmp(k + 1) = jtmp
                # これはtractmpはnref個の配列だが、inside判定がない場合は0が入って、
                # 入っている場合は最も近い交点になる壁番号が入る？
                # 0の配列は必要なのか？必要ないなら以下のコードで書く　必要ならelseで0を足す　＊＊部分
                reflectionmeshid_history.append(mesh_nearestid)
                # ここ＊
                t = mm.parameter_t(sound_ray, soundray_comesfrom, mesh[mesh_nearestid].normal,
                                   mesh[mesh_nearestid].vertexes)
                node = mm.node_renew(sound_ray, soundray_comesfrom, t)
                soundray_comesfrom = sr.soundraycomesfrom_renew(node)
                sound_ray = sr.reflection_generator(sound_ray, mesh[mesh_nearestid].normal)

            else:
                # ＊＊0配列が必要な場合
                reflectionmeshid_history.append(mesh_nearestid)


    return reflectionmeshid_history
