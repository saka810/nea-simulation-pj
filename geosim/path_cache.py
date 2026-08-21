"""経路の幾何を保存して、吸音材だけ差し替えた再計算に使う（F-9。依頼 2026-08-21）。

> 音線ループからバックトレースまでは、とくに吸音の違いは関係なく、
> 最後の虚像法による計算の時にようやく吸音の効果が反映されるのかなと思っています。
> つまり、わざわざ最初から計算し直さなくても、最後の部分だけデータを残しておけば、
> データの流用および、再計算時に、そこから再開できませんか？

## そのとおり。ただし正確には「幾何」と「吸音」の境目はもう一段あとにある

| 段階 | 吸音に依存するか |
|---|---|
| ① 音線追跡 | **しない**（受音しても打ち切らない。反射は鏡面反射だけ） |
| ② 重複削除 | しない |
| ③ バックトレース・経路の検証 | **しない**（虚音源が成立するかは純粋に幾何） |
| ③' バックトレース・エネルギー | **する**（各反射で `|R(θ,α)|²` を掛ける） |
| ④ インパルス応答以降 | する（③' の結果を使う） |

つまり吸音が効くのは③'だけで、そこは

    E = Π_k |R(cosθ_k, α_{面_k})|²

という**掛け算のかたまり**にすぎない。だから

- 各経路の**反射面の並び**と**各反射の入射角 cosθ**
- 到来時刻・距離・到来方向（これも幾何）

を残しておけば、吸音材を変えたときの再計算は**配列の掛け算だけ**で済む。
実測（研修室・受音経路 8,825 本・最大 300 反射）で **①〜③ の 10 分 → 0.1 秒未満**。

## ★ただし「吸音に依らない」には 1 つ落とし穴がある

交差判定は**同一平面パッチ単位**で、パッチは「同一平面＋辺で連結＋**同じ材料**＋
同じ法線の向き」でまとめている（材料で割らないと吸音率が引けないため）。
つまり**材料の割り当てを変えるとパッチの切れ目が変わり、見つかる経路も変わる**。

- 材料の**値**（吸音率）を変えるだけ → パッチは同じ。**再利用できる**
- 材料の**割り当て方**を変えて、隣り合う同一平面の面が別材料になったり
  同じ材料になったりする → パッチが変わる。**再利用できない**

そこで保存時に**パッチの分け方の指紋**を残し、読むときに突き合わせる。
食い違ったら黙って使わずに、理由を告げて①から計算し直す。

## ファイル

`結果/recN/<室>_経路.npz`（**条件名は付けない**。条件をまたいで共有するため）。
大きさは経路 × 最大反射回数で決まる（研修室 5 点で約 40 MB。圧縮あり）。
"""

import hashlib
import os

import numpy as np

# 指紋に入れる項目（ここが 1 つでも違えば経路が変わる）
FINGERPRINT_KEYS = ("geometry", "patches", "source", "receiver",
                    "rays", "nref", "radius", "two_sided")

CACHE_VERSION = 1


def _digest(*arrays):
    """配列の並びから 16 進のダイジェストを作る（浮動小数のビットをそのまま見る）。"""
    h = hashlib.sha1()
    for array in arrays:
        a = np.ascontiguousarray(array)
        h.update(str(a.dtype).encode())
        h.update(str(a.shape).encode())
        h.update(a.tobytes())
    return h.hexdigest()


def fingerprint(mesh, faces, source, receiver, rays, nref, radius, two_sided):
    """経路が決まる条件の指紋。**吸音率の値は入れない**（変えても経路は同じ）。

    入れるもの:
        geometry … 三角形の頂点と法線（DXF が変わった／法線を反転した）
        patches  … パッチの分け方（**材料の割り当てで変わる**。上記の落とし穴）
        source / receiver / rays / nref / radius / two_sided
    """
    vertices = np.array([np.asarray(m.vertexes, dtype=float) for m in mesh])
    normals = np.array([np.asarray(m.normal, dtype=float) for m in mesh])
    patches = np.asarray(getattr(faces, "patch_of_face",
                                 np.arange(len(mesh))), dtype=np.int64)
    return {
        "geometry": _digest(vertices, normals),
        "patches": _digest(patches),
        "source": np.round(np.asarray(source, dtype=float), 9).tolist(),
        "receiver": np.round(np.asarray(receiver, dtype=float), 9).tolist(),
        "rays": int(rays),
        "nref": int(nref),
        "radius": float(radius),
        "two_sided": bool(two_sided),
    }


def compare(saved, current):
    """指紋を突き合わせる。合っていれば None、違えば理由の文字列を返す。"""
    reasons = {
        "geometry": "モデルの形か法線が変わっています",
        "patches": "材料の割り当て方が変わってパッチの切れ目が動いています"
                   "（吸音率の値を変えるだけなら再利用できます）",
        "source": "音源の位置が違います",
        "receiver": "受音点の位置が違います",
        "rays": "音線の本数が違います",
        "nref": "最大反射回数が違います",
        "radius": "受音球の半径が違います",
        "two_sided": "面の裏からの入射の扱いが違います",
    }
    for key in FINGERPRINT_KEYS:
        a, b = saved.get(key), current.get(key)
        if isinstance(b, float):
            same = a is not None and np.isclose(float(a), b)
        elif isinstance(b, list):
            same = (a is not None and len(a) == len(b)
                    and np.allclose(np.asarray(a, dtype=float),
                                    np.asarray(b, dtype=float)))
        else:
            same = a == b
        if not same:
            return f"{reasons[key]}（{key}: 保存 {a!r} / いま {b!r}）"
    return None


def save(filename, pulses, mark, verbose=True):
    """経路の幾何を npz に保存する。エネルギーは**入れない**（吸音に依るので）。"""
    if not len(pulses) or not len(pulses.walls):
        return None
    os.makedirs(os.path.dirname(os.path.abspath(filename)) or ".", exist_ok=True)
    payload = {
        "version": CACHE_VERSION,
        "reflection_count": np.asarray(pulses.reflection_count, dtype=np.int32),
        "distance": np.asarray(pulses.distance, dtype=np.float64),
        "direction": np.asarray(pulses.direction, dtype=np.float64),
        "walls": np.asarray(pulses.walls, dtype=np.int32),
        "cos_theta": np.asarray(pulses.cos_theta, dtype=np.float64),
    }
    for key in FINGERPRINT_KEYS:
        payload[f"fp_{key}"] = np.asarray(mark[key])
    np.savez_compressed(filename, **payload)
    if verbose:
        size = os.path.getsize(filename) / 1024.0 / 1024.0
        print(f"[経路] 保存しました（{len(pulses)} 本 × 最大 {pulses.walls.shape[1]} 反射 / "
              f"{size:.1f} MB）: {filename}")
    return filename


def load(filename, mark=None, sound_velocity=None, verbose=True):
    """保存した経路を `PulseList` に戻す。使えなければ None。

    `mark` を渡すと指紋を突き合わせ、食い違えば理由を告げて None を返す。
    **エネルギーは入っていない**ので、戻したあと `recompute_energy(mesh)` を呼ぶこと。
    """
    import loop_noredundancy as ln

    if not filename or not os.path.exists(filename):
        return None
    try:
        data = np.load(filename, allow_pickle=False)
    except Exception as error:
        print(f"[経路] {filename} を読めませんでした: {type(error).__name__}: {error}")
        return None
    if int(data.get("version", 0)) != CACHE_VERSION:
        print(f"[経路] 形式が古いので使いません: {filename}")
        return None

    if mark is not None:
        saved = {}
        for key in FINGERPRINT_KEYS:
            value = data[f"fp_{key}"]
            if value.ndim == 0:
                saved[key] = value.item()
            else:
                saved[key] = value.tolist()
        reason = compare(saved, mark)
        if reason:
            if verbose:
                print(f"[経路] 保存した経路は使えません: {reason}")
                print("[経路] 音線追跡からやり直します")
            return None

    pulses = ln.PulseList(0, sound_velocity or ln.SOUND_VELOCITY)
    pulses.reflection_count = data["reflection_count"].astype(int)
    pulses.distance = data["distance"].astype(float)
    pulses.direction = data["direction"].astype(float)
    pulses.walls = data["walls"].astype(np.int32)
    pulses.cos_theta = data["cos_theta"].astype(float)
    pulses.time = pulses.distance / (sound_velocity or ln.SOUND_VELOCITY)
    pulses.energy = np.zeros((len(pulses.distance), 0))
    if verbose:
        print(f"[経路] 保存した経路を使います（{len(pulses.distance)} 本）: {filename}")
    return pulses


def reuse(filename, mesh, mark, band_number, sound_velocity, verbose=True):
    """経路を読み、**吸音率だけ差し替えてエネルギーを作り直した** PulseList を返す。

    使えなければ None（呼び側は従来どおり音線追跡から回す）。
    """
    pulses = load(filename, mark, sound_velocity, verbose=verbose)
    if pulses is None:
        return None
    pulses.recompute_energy(mesh, band_number=band_number)
    pulses.sort_by_time()
    if verbose:
        print(f"[経路] 吸音率を当てて作り直しました: {pulses.summary()}")
    return pulses
