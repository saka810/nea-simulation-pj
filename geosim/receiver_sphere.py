import numpy as np



# 7/9打ち合わせ用


# 受音点に関わるものをこちらに記述します

# 受音球の中を音線が通過したかを判定する
# OK
def inside_sphere(sphere_radius, sound_ray, soundray_comesfrom, receiver_point, min_distance):
    inside = False

    vector = receiver_point - soundray_comesfrom
    inner_product = np.dot(sound_ray, vector)
    distance = np.linalg.norm(soundray_comesfrom + sound_ray * inner_product - receiver_point, ord=2)

    if distance <= sphere_radius:
        if distance <= min_distance:
            if inner_product >= 0:
                inside = True

    return inside
