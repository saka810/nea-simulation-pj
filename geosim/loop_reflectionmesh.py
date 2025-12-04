import numpy as np
from numpy.ma.core import count

import sound_ray as sr
import mesh_method as mm
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
# 音源位置、受音点、音源から出る音線全て持つ配列、反射回数、室形状のメッシュ全ての情報、受音球の半径
# を用いて、
# 面に反射する音線を追跡し、反射した面のインデックス（番号）を
# 反射履歴　reflectionmeshid_history として記録する
# このreflectionmeshid_historyを出力として返す
# 元コード524
def loop(soundsource_point, reciever_point, soundray_list, nref, mesh, sphere_radius):
    reflectionmeshid_history_2dim=[]
    collision = False

    for soundray_i in range(len(soundray_list)):
        # reflectionmeshid_historyは今のコードだと一時反射履歴なので，ここでいちど初期化する必要あり
        # 「loop_deleteredundancy.py」に持って行くための二次元配列[受音数，受音するまでの反射経路（受音毎に変わる）]が別途必要
        reflectionmeshid_history = []
        reflectionmeshid_history.append(-1)  # 0番目は-1を入れてる 後で追加するときにインデックスをずらす必要があるみたい

        soundray_comesfrom = soundsource_point

        # 一回目のみなので，外に出しました
        sound_ray = soundray_list[soundray_i, :]
        for k in range(nref):

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

            # 12/04 追記メモ
            # 確認したい点はこの if insideで保存する反射経路リストについて
            # if inside:
                #受音した場合は，反射経路リストに保存する（重複経路削除に持って行く形式）
                #ここの配列は2次元になるはず

                # 12/04追記メモ
                # 元コード662行目～677行目が理解ができていない部分です。申し訳ございません。
                # ここのif insideは元コードtraceffに代入する部分と想像しています。
                # 以下は662~671部分
                # count = count + 1
                # do j = 0, nref
                #     traceff(count, j) = tractmp(j)
                # end  do

                # [疑問1]
                # do jのループは
                # このループ１個外側の
                # !■反射回数ループ■　
                # pythonコード for k in range(nref):　（元コード545行目do k = 0, nref）
                # のループと同じ範囲なので else:breakの後にfor k in range(nref):を抜けてから記載しても良いでしょうか

                # [疑問2]
                # traceffに代入するtractmpのが元コード668行目
                # tractmpに値を代入するのが元コード677行目
                # 代入の順序がこの順番になっていて、668行目で毎回nrefまで代入をするのはどのような手順を想定されていますでしょうか。
                # この部分は一つ外側　元コード545行目 do k = 0, nref　と　
                # さらに外側の元コード527行目do i = 1, nrayと
                # 元コード667行目 do j = 0, nref でnray*nref*nref回ループが回り
                # 代入されるtraceffは
                # 配列の行数がinside判定があるか無いかが関わるので countの数＜＝nray*nref　、
                # 配列の列数がnref個のように見えたので、
                # 元コード545行目do k = 0, nrefを抜けてからnref回のループでも良いものでしょうか。
                # 今のループだと　nray*nref*nrefのループが行われていますでしょうか？
                # 必要なループの個数はnray*nref個でしょうか。

                # [疑問3]
                # 545行目do k = 0, nrefを考えると
                # do j = 0, nrefは do j = 0, kとしても良いものでしょうか。
                # k+1~nref では traceff[k+1:nref]=0 ?

                # [疑問4]
                # 最終的にinside も collisionの条件を満たすものが残れば良いでしょうか。
                # その場合、reflectionmeshid_history（元コードtractmp）追加条件をinside も collisionの条件を満たすに変更し
                # reflectionmeshid_historyが一つのkに対して完成してからtraceffに追加しても良いか
                # この場合の変更想定箇所　** 、***、****

            # 662行目以降　再度打ち合わせ
            # 662~671は if inside内で行うもの ????

            # 元コード↓
            # ! 一時反射履歴に壁面番号を書き込み
            # tractmp(k + 1) = jtmp
            # これはtractmpはnref個の配列だが、inside判定がない場合は0が入って、
            # 入っている場合は最も近い交点になる壁番号が入る？
            # 0の配列は必要なのか？必要ないなら以下のコードで書く　必要ならelseで0を足す　＊＊部分

            #　受音したかにかかわらず，次の壁に当たる場合に，衝突する壁・交点等を求める。
            if collision:
                # 12/04追記メモ **
                # reflectionmeshid_historyはここで追加はなしにして***で追加
                # reflectionmeshid_history.append(mesh_nearestid)
                #　↑これは一時的な反射経路（1次元）


                # 12/04追記メモ
                # ***
                if inside:
                    reflectionmeshid_history.append(mesh_nearestid)

                t = mm.parameter_t(sound_ray, soundray_comesfrom, mesh[mesh_nearestid].normal,
                                   mesh[mesh_nearestid].vertexes)
                node = mm.node_renew(sound_ray, soundray_comesfrom, t)
                soundray_comesfrom = sr.soundraycomesfrom_renew(node)
                sound_ray = sr.reflection_generator(sound_ray, mesh[mesh_nearestid].normal)
            else:
                break
                #衝突する壁が無ければループを抜ける


            # 12/04追記メモ
            # ****
            # ここでtraceffに追加
            reflectionmeshid_history_2dim.append(reflectionmeshid_history)



            #受音しない場合は，履歴に保存する必要はない
            # else:
            #     # ＊＊0配列が必要な場合
            #     reflectionmeshid_history.append(mesh_nearestid)

    # 12/04追記メモ
    # returnはtractmp(reflectionmeshid_history)ではなくtraceff(reflectionmeshid_history_2dim)でしょうか。
    # return reflectionmeshid_history
    return reflectionmeshid_history_2dim
