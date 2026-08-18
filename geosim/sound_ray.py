"""音線ベクトルそのものを扱う関数群（生成・正規化・反射・エネルギー減衰）。

元コード `backtrace.f90` の音源まわり。面や受音点には触らず、
**ベクトルの計算だけ**をここに集めてある。

| この関数の引数 | 意味 | 元コードの変数 |
|---|---|---|
| `sound_ray` | 音線の方向ベクトル | `vray` |
| `soundray_comesfrom` | 音線の基点 | `vini` |
| `normal` | 面の法線 | `vn` |
| `imaginarysound_point` | 虚音源の座標 | `isrc` |

**注意：`sound_ray` は必ずしも音源から出ていない。** 壁で反射したあとの音線も
`sound_ray` と呼んでいる（むしろそちらが多い）。一度も反射していないときだけ
`soundray_comesfrom` が音源の位置になる。

各関数の詳しい説明・元コードの行番号・数式の根拠は、関数の直前のコメントに書いてある。
"""


import numpy as np


# 7/9打ち合わせ用


# 音源、虚音源に関するものをこちらに集めています

# 音源から出る音線を作成します
# OK
def soundray_generator(ray_number):
    # 2026-08-12 修正: np.zeros は形状をタプルで渡す（np.zeros(ray_number, 3) は TypeError）
    sound_rays = np.zeros((ray_number, 3))
    dt = np.pi * (3.0 - np.sqrt(5.0))
    dz = 2.0 / ray_number

    # 添字は 0 始まり（元コードは do i = 1, nray）。
    # 元コードは i = nray のとき z = 1 + 1/nray > 1 となり sqrt(1 - z^2) が NaN になるが、
    # 0 始まりだと z は -1 + 1/nray 〜 1 - 1/nray に収まり、球面上に正しく均等分布する。
    for i in range(ray_number):
        # 2026-08-12 修正: sound_rays(i, 2) → sound_rays[i, 2]（丸括弧では呼び出し扱いになる）
        sound_rays[i, 2] = dz * i - 1.0 + 1.0 / ray_number
        sound_rays[i, 0] = np.sqrt(1.0 - sound_rays[i, 2] ** 2.0) * np.cos(dt * i)
        sound_rays[i, 1] = np.sqrt(1.0 - sound_rays[i, 2] ** 2.0) * np.sin(dt * i)

    return sound_rays


# 音線の正規化
# OK
def normalized_soundray(sound_ray):
    # normalized_ray = np.array(np.zeros(3))
    distance = np.linalg.norm(sound_ray, ord=2)
    normalized_ray = sound_ray / distance
    return normalized_ray


# 反射音線ベクトルの作成と正規化 (バックトレースではこれはつかわない)
# OK
def reflection_generator(sound_ray, normal):
    # r reflection = np.array(np.zeros(3))
    t = np.dot(sound_ray, normal)
    reflection = sound_ray - 2.0 * t * normal
    reflection = normalized_soundray(reflection)
    return reflection


# 反射音の音線ベクトルの基点を更新する　交点と入れ替えるだけのもの
# OK
def soundraycomesfrom_renew(node):
    return node


# 虚音源の座標と反射音の音線ベクトルの基点を結んで反射音の音線ベクトルを作成する
# 元コード1101
# OK
def soundray_renew(imaginarysound_point, soundray_comesfrom):
    # new_soundray = np.array(np.zeros(3))
    new_soundray = imaginarysound_point - soundray_comesfrom
    return new_soundray


# 虚音源になったときのエネルギーを減衰させるメソッド
# 元コード1091〜1094
#
# 【元の数式の根拠】局所反応性壁面の斜入射エネルギー反射率。
#   法線入射吸音率 α0 → 法線入射の圧力反射率 |R0| = sqrt(1 - α0)
#   → 規格化音響インピーダンス（実数と仮定） z = (1 + |R0|) / (1 - |R0|)
#   → 斜入射の圧力反射率 R(θ) = (z*cosθ - 1) / (z*cosθ + 1)
#   分母分子に (1 - |R0|) を掛けて整理すると
#     R(θ) = ((1 + sqrt(1-α0)) * cosθ - (1 - sqrt(1-α0)))
#          / ((1 + sqrt(1-α0)) * cosθ + (1 - sqrt(1-α0)))
#   エネルギー反射率は |R(θ)|^2。元コードはこの形で、数式として正しい。
#   （θ→90° で |R|→1 = 吸音なし、α0=1 でも斜入射では反射が残る。局所反応モデルとして妥当な挙動）
#
# 想定されている吸音率は垂直入射（法線入射）吸音率。
# ＮＥＡは残響室吸音率が多い？ 残響室吸音率の場合は書き換えが必要
#   （残響室法 α は 1 を超えることがあり sqrt(1 - α) が NaN になる。
#     Paris の式を逆解きして α0 または z を求める前処理が要る）
def energy_decay(sound_ray, normal, absorption, initial_energy):
    coefficient = np.abs(np.dot(sound_ray, normal))  # cosθ（音線・法線とも単位ベクトル前提）
    # 2026-08-12 修正: 第2項の符号。元コードは (1 - sqrt(1-α))、旧実装は 2 箇所とも (1 + sqrt(1-α)) だった。
    # 旧実装のままだと R = (cosθ - 1)/(cosθ + 1) となり吸音率に一切依存しなくなる。
    reflection_energy = (1 + np.sqrt(1 - absorption)) * coefficient - (1 - np.sqrt(1 - absorption))
    reflection_energy = reflection_energy / (
                (1 + np.sqrt(1 - absorption)) * coefficient + (1 - np.sqrt(1 - absorption)))
    reflection_energy = abs(reflection_energy)
    # 元コードは abs(...)**2.0d0 * enertmp(l)。2 乗しているのは圧力反射率 R で、
    # 最後に掛けるのは累積エネルギー enertmp。3 乗ではない。
    reflection_energy = reflection_energy * reflection_energy
    reflection_energy = initial_energy * reflection_energy
    return reflection_energy
