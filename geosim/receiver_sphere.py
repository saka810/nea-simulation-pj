import numpy as np



# 7/9打ち合わせ用


# 受音点に関わるものをこちらに記述します

# 受音球の中を音線が通過したかを判定する
# OK
def inside_sphere(sphere_radius, sound_ray, soundray_comesfrom, receiver_point, min_distance):
    inside = False

    # 元コード649〜663
    vector = receiver_point - soundray_comesfrom            # vtgt: 基点→受音点
    inner_product = np.dot(sound_ray, vector)               # distd: 音線上への射影距離（垂線の足まで）
    distance = np.linalg.norm(soundray_comesfrom + sound_ray * inner_product - receiver_point, ord=2)
                                                            # distr: 受音点から音線への垂線距離

    if distance <= sphere_radius:                           # distr <= rcvr
        # 2026-08-12 修正: 元コード663行は distd .le. disttmp（射影距離と壁までの距離の比較）。
        # 旧実装は distance（垂線距離）と比較しており、壁の向こう側の受音球でも
        # 受音と誤判定しうる状態だった。
        if inner_product <= min_distance:                   # distd <= disttmp
            if inner_product >= 0:                          # distd >= 0（前方）
                inside = True

    return inside
