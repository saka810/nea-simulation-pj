import numpy as np


# 7/9打ち合わせ用



# 音源、虚音源に関するものをこちらに集めています

# 音源から出る音線を作成します
def soundray_generator(ray_number):
    sound_rays = np.array([np.zeros(3), 3])
    dt = np.pi * (3.0 - np.sqrt(5.0))
    dz = 2.0 / ray_number

    for i in range(ray_number):
        sound_rays[i, 2] = dz * i - 1.0 + 1.0 / ray_number
        sound_rays[i, 0] = np.sqrt(1.0 - sound_rays(i, 2) ** 2.0) * np.cos(dt * i)
        sound_rays[i, 1] = np.sqrt(1.0 - sound_rays(i, 2) ** 2.0) * np.sin(dt * i)

    return sound_rays


# 音線の正規化
def noramlized_soundray(sound_ray):
    # normalized_ray = np.array(np.zeros(3))
    distance = np.linalg.norm(sound_ray, ord=2)
    normalized_ray = sound_ray / distance
    return normalized_ray


# 反射音線ベクトルの作成と正規化 (バックトレースではこれはつかわない)
def reflection_generator(sound_ray, normal):
    # reflection = np.array(np.zeros(3))
    t = np.dot(sound_ray, normal)
    reflection = sound_ray - 2.0 * t * normal
    reflection = noramlized_soundray(reflection)
    return reflection


# 反射音の音線ベクトルの基点を更新する　交点と入れ替えるだけのもの
def soundraycomesfrom_renew(node):
    return node


# 虚音源の座標と反射音の音線ベクトルの基点を結んで反射音の音線ベクトルを作成する
def soundray_renew(imaginarysound_point, soundray_comesfrom):
    new_soundray = np.array(np.zeros(3))
    new_soundray = imaginarysound_point - soundray_comesfrom
    return new_soundray


# 虚音源になったときのエネルギーを減衰させるメソッドをここに追加予定
# 元コード1091
# 想定されている吸音率は垂直入射？ＮＥＡは残響室吸音率が多い？
def energy_decay(sound_ray, normal, absorption):
    coefficient = np.linalg.norm(sound_ray * normal, ord=2)
    energy = 1 + np.sqrt(1 - absorption) * coefficient - (1 + np.sqrt(1 - absorption))
    energy = energy / (1 + np.sqrt(1 - absorption) * coefficient + (1 + np.sqrt(1 - absorption)))
    energy = abs(energy)
    energy = energy * energy * energy
    # 元コードが2条の後自身を掛ける形になっているが3条？
    return energy
