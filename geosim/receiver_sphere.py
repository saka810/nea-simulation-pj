"""受音判定。**受音点に半径 `sphere_radius` の球を置き、音線が通ったかを見る。**

元コード `backtrace.f90` 649〜663 行。

音線は線なので、点である受音点にはまず当たらない。そこで球を「経路を見つけるための網」
として置く。**半径は精度のつまみ**で、小さすぎると後期の経路を取りこぼし、
残響時間が短く出る（`t_max = r√N / (2c)` が拾える時刻の目安）。
半径を変えて値が動かないことを確認してから数値を使うこと。

★受音球はあくまで**経路を見つけるための網**であって、エネルギーの係数ではない。
  エネルギーは虚音源から `1/(4πd²)` で出し直す（`loop_noredundancy`）ので、
  音線法で普通に問題になる**受音球の断面積による正規化バイアスは原理的に無い**。

| 引数 | 意味 | 元コードの変数 |
|---|---|---|
| `sphere_radius` | 受音球の半径 | `rcvr` |
| `sound_ray` / `soundray_comesfrom` | 音線の方向 / 基点 | `vray` / `vini` |
| `receiver_point` | 受音点の座標 | `rcv` |
| `min_distance` | その音線が壁に当たるまでの距離 | `disttmp` |

`inside_sphere`（1 本ずつ）と `inside_sphere_batch`（束をまとめて）は同じ式・同じ比較で、
結果は一致する。本番は後者を使う。

## 判定の中身（2026-09-05 に修正）

★**「基点から壁までの線分」と球が交わるか**を見る（`method='segment'`。既定）。

    t* = clip(音線方向・(受音点 - 基点), 0, 壁までの距離)     ← 線分上の最近点
    |基点 + 方向 t* - 受音点| <= 半径

それまでは**垂線の足**が球の中にあり、かつ足までの距離が壁までの距離以内、
という条件だった（`method='foot'`。元コードのやり方）。これだと

- **足が壁より先にある**とき、線分自身は球を通っていても落とす
  （壁まで 1.0 m・球の中心 1.1 m・半径 0.2 m なら、線分は 0.9〜1.0 m で球の中にいる）
- **足が基点より後ろにある**とき（`t < 0`）も同様に落とす。
  基点が球の中にいても拾えない

どちらも**受音点が壁から受音球の半径以内**にあるときに起きる。
線分と球の判定にすると、落としていた経路が拾えるようになる（判定は必ず**増える**方向で、
以前拾えていたものが拾えなくなることはない）。

`method='foot'` は元コードのやり方を再現するために残してある（**参照実装**）。
"""


import numpy as np


# 判定の仕方。
#   'segment' … 基点から壁までの**線分**と球の交差（既定。2026-09-05）
#   'foot'    … 垂線の足が球の中かどうか（元コード 649〜663 行。**参照実装**）
METHOD_SEGMENT = "segment"
METHOD_FOOT = "foot"
METHODS = (METHOD_SEGMENT, METHOD_FOOT)
DEFAULT_METHOD = METHOD_SEGMENT


# 7/9打ち合わせ用


# 受音点に関わるものをこちらに記述します

# 受音球の中を音線が通過したかを判定する
# OK
def inside_sphere(sphere_radius, sound_ray, soundray_comesfrom, receiver_point,
                  min_distance, method=DEFAULT_METHOD):
    """音線 1 本が受音球を通過したか。**参照実装**（本番は `inside_sphere_batch`）。

    `method='segment'`（既定）… 基点から壁までの**線分**が球と交わるか。
      線分上で受音点にいちばん近い点を求め、その距離を半径と比べる。

    `method='foot'`（元コード 649〜663 行）… 次の 3 つを満たしたら受音。
      ・受音点から音線への**垂線距離**が半径以内
      ・垂線の足までの**射影距離**が壁に当たるまでの距離以内（＝壁の手前）
      ・射影距離が 0 以上（＝前方。後ろの受音点を拾わない）
    """
    if method not in METHODS:
        raise ValueError(f"method は {METHODS} のいずれかです: {method!r}")

    # 元コード649〜663
    vector = receiver_point - soundray_comesfrom            # vtgt: 基点→受音点
    inner_product = np.dot(sound_ray, vector)               # distd: 音線上への射影距離（垂線の足まで）

    if method == METHOD_FOOT:
        distance = np.linalg.norm(
            soundray_comesfrom + sound_ray * inner_product - receiver_point, ord=2)
                                                            # distr: 受音点から音線への垂線距離
        if distance <= sphere_radius:                       # distr <= rcvr
            # 2026-08-12 修正: 元コード663行は distd .le. disttmp（射影距離と壁までの距離の比較）。
            # 旧実装は distance（垂線距離）と比較しており、壁の向こう側の受音球でも
            # 受音と誤判定しうる状態だった。
            if inner_product <= min_distance:               # distd <= disttmp
                if inner_product >= 0:                      # distd >= 0（前方）
                    return True
        return False

    # ★線分と球の交差（既定）。垂線の足を線分の中へ押し込んでから距離を測る
    nearest = min(max(inner_product, 0.0), float(min_distance))
    closest = soundray_comesfrom + sound_ray * nearest
    return bool(np.linalg.norm(closest - receiver_point, ord=2) <= sphere_radius)


def inside_sphere_batch(sphere_radius, directions, origins, receiver_point,
                        min_distance, method=DEFAULT_METHOD):
    """`inside_sphere` を音線の束についてまとめて判定する（F-1 高速化）。

    引数:
        directions (A,3) 音線の方向（単位ベクトル前提）
        origins    (A,3) 音線の基点
        receiver_point (3,)
        min_distance (A,) 各音線が壁に当たるまでの距離（当たらなければ inf）
    戻り値:
        (A,) の bool

    計算は scalar 版とまったく同じ式・同じ比較なので結果は一致する。
    """
    if method not in METHODS:
        raise ValueError(f"method は {METHODS} のいずれかです: {method!r}")

    receiver = np.asarray(receiver_point, dtype=float)
    vector = receiver[None, :] - origins
    inner_product = np.einsum("ij,ij->i", directions, vector)

    if method == METHOD_FOOT:
        foot = origins + directions * inner_product[:, None]  # 音線への垂線の足
        distance = np.linalg.norm(foot - receiver, axis=1)
        return ((distance <= sphere_radius)
                & (inner_product <= min_distance)
                & (inner_product >= 0.0))

    # ★線分上の最近点。`min_distance` は壁に当たらなければ inf なので clip の上限に使える
    nearest = np.clip(inner_product, 0.0, min_distance)
    closest = origins + directions * nearest[:, None]
    return np.linalg.norm(closest - receiver, axis=1) <= sphere_radius
