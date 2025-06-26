import numpy as np
import sound_ray as sr
import mesh_method as mm
from geosim.sound_ray import reflection_generator
from mesh import Mesh
import receiver_sphere as rs


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
# 元コード524
def loop(soundsource_point, reciever_point, soundray_list, nref, mesh, sphere_radius):
    reflectionmeshid_history = []
    reflectionmeshid_history.append(-1)  # 0番目は-1を入れてる 後で追加するときにインデックスをずらす必要があるみたい

    for soundray_i in range(len(soundray_list)):

        for k in range(nref):
            sound_ray = soundray_list[soundray_i, :]
            sound_ray = sr.noramlized_soundray(sound_ray)
            soundray_comesfrom = soundsource_point

            min_distance = 1000000000.0  # とりあえず大きい数字

            # for j in range(mesh_count):
            for j in range(len(mesh)):
                collision = False

                if np.dot(sound_ray, mesh[j].normal) < 0:
                    t = mm.parameter_t(sound_ray, soundray_comesfrom, mesh[j].normal, mesh[j].vertexes)

                    if t > 0:
                        node = mm.node_renew(sound_ray, soundray_comesfrom, t)
                        collision, distance = mm.collision_distance(sound_ray, soundray_comesfrom, mesh[j].normal,
                                                                    mesh[j].vertexes)

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

            if inside and collision:
                reflectionmeshid_history.append(k)
                # ここ＊
                t = mm.parameter_t(sound_ray, soundray_comesfrom, mesh[mesh_nearestid].normal,
                                   mesh[mesh_nearestid].vertexes)
                node = mm.node_renew(sound_ray, soundray_comesfrom, t)
                soundray_comesfrom = sr.soundraycomesfrom_renew(node)
                sound_ray = sr.reflection_generator(sound_ray, mesh[mesh_nearestid].normal)

    return reflectionmeshid_history
