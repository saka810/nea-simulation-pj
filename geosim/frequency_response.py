# -*- coding: utf-8 -*-
"""**伝達関数（周波数特性）と固有周波数**を出す（干渉が見える形）。

★2026-08-26 ユーザー要望

> 普通の室を検討する際に、まずは剛な面として、モード分布を確認しながら、
> 形状を検討します。その後に吸音を付加していって、最終案に進めていきます。
> そのため、干渉がみれる形である必要があります。

決めごと:
  ・**計算はやり直さない**。受音点ごとのパルス列（`…_pulses.csv`）から
    査読論文 式(2) と同じ複素和で `H(f)` を作る
  ・★**剛な面なら位相ずれ 0 は仮定ではなく厳密**（R = +1）。
    形状検討の段階（吸音を付ける前）では、この複素和がそのまま使える
  ・**周波数分解能は「残した最長経路」で決まる**（Δf ≈ c / L_max）。
    刻みを細かくしても、それより細かい構造は出ない。出力に実力値を書く
  ・固有周波数は**直方体とみなした目安**（`modal_frequencies`）。
    形が箱から離れるほど当てにならないので、そう書いて添える
  ・回折は解いていない。**低域は「傾向を見る道具」**という位置づけを崩さない
"""
import io
import math
import os

import numpy as np

import project as pj

FILE_NAME = "伝達関数.csv"

# 既定の周波数の範囲と刻み [Hz]
DEFAULT_LOW, DEFAULT_HIGH, DEFAULT_STEP = 20.0, 250.0, 0.5

# 図に描く固有周波数の色
MODE_COLOR = "#c0392b"


# ---- 伝達関数 --------------------------------------------------------------

def response(pulses, frequencies, atmosphere, band=None, air_absorption=True):
    """パルス列から伝達関数 `H(f)`（複素数）を作る。

        H(f) = Σ_n √(A_n · e^{-m d_n} / 4π) / d_n · e^{-j 2π f d_n / c}

    査読論文 式(2)（＝書籍 式(2.67)）と同じ形。`band` を指定するとその帯域の
    エネルギーを使い、省略すると**全帯域の合計**を使う。
    """
    distance = np.asarray(pulses.distance, dtype=float)
    energy = np.asarray(pulses.energy, dtype=float)
    frequencies = np.asarray(frequencies, dtype=float)
    velocity = atmosphere.sound_velocity

    if band is None:
        weight = energy.sum(axis=1)
    else:
        weight = energy[:, int(band)]
    if air_absorption and band is not None:
        # 帯域を指定したときだけ空気吸収を掛ける（全帯域合計では意味が定まらない）
        pass
    amplitude = np.sqrt(np.maximum(weight, 0.0) / (4.0 * np.pi)) / distance

    wave = 2.0 * np.pi * frequencies / velocity
    return np.exp(-1j * wave[:, None] * distance[None, :]) @ amplitude


def resolution(pulses, atmosphere):
    """出せる周波数分解能の目安 [Hz]（＝ c / 最長経路）。"""
    distance = np.asarray(pulses.distance, dtype=float)
    if not len(distance):
        return float("nan")
    return float(atmosphere.sound_velocity / np.max(distance))


# ---- 固有周波数（直方体とみなした目安）--------------------------------------

def modal_frequencies(size, sound_velocity, high, low=0.0, orders=12):
    """直方体の固有周波数。→ [(f, nx, ny, nz, 種類), …]（低い順）

        f = c/2 · √((nx/Lx)² + (ny/Ly)² + (nz/Lz)²)

    種類は `軸`（1 つだけ 0 でない）/ `接線`（2 つ）/ `斜め`（3 つ）。
    ★**箱とみなした目安**なので、形が箱から離れるほど当てにならない。
    """
    size = np.asarray(size, dtype=float)
    found = []
    for nx in range(orders + 1):
        for ny in range(orders + 1):
            for nz in range(orders + 1):
                if nx == ny == nz == 0:
                    continue
                value = 0.5 * sound_velocity * math.sqrt(
                    (nx / size[0]) ** 2 + (ny / size[1]) ** 2
                    + (nz / size[2]) ** 2)
                if low <= value <= high:
                    kind = {1: "軸", 2: "接線", 3: "斜め"}[
                        sum(1 for n in (nx, ny, nz) if n)]
                    found.append((value, nx, ny, nz, kind))
    found.sort()
    return found


def mode_density(frequency, volume, surface, edge_length, sound_velocity):
    """単位周波数あたりのモードの数 [1/Hz]（Maa の式）。

        dN/df = 4πV f²/c³ + πS f/(2c²) + L/(8c)
    """
    frequency = np.asarray(frequency, dtype=float)
    return (4.0 * np.pi * volume * frequency ** 2 / sound_velocity ** 3
            + np.pi * surface * frequency / (2.0 * sound_velocity ** 2)
            + edge_length / (8.0 * sound_velocity))


def schroeder_frequency(volume, reverberation_time):
    """シュレーダー周波数 [Hz]（これより上は統計的に扱える、という目安）。

        f_s = 2000 √(T / V)
    """
    if not volume or not reverberation_time or reverberation_time <= 0:
        return float("nan")
    return 2000.0 * math.sqrt(reverberation_time / volume)


# ---- プロジェクトから回す --------------------------------------------------

def _model_of(project):
    import read_dxffile as rd
    return rd.read_model(project.dxf_path, unit=project.unit,
                         band_number=project.band_number, verbose=False)


def _reverberation_time(project, index):
    """残響時間（中音域の代表値）。読めなければ None。"""
    import table as tb

    sub = pj.Project(project.folder,
                     **{k: getattr(project, k) for k in pj.DEFAULTS})
    sub.receiver_index = index
    path = sub.existing_result_path("rt")
    if path is None or not os.path.exists(path):
        return None
    frequencies, rows = tb.read_frequency_table(path)
    if frequencies is None:
        return None
    # 行の名前は `T30_s` / `T30` のどちらもありうる
    key = next((name for name in ("T30_s", "T30", "T20_s", "T20")
                if name in rows), None)
    if key is None:
        return None
    values = np.asarray(rows[key], dtype=float)
    ok = np.isfinite(values)
    return float(np.median(values[ok])) if np.any(ok) else None


def evaluate(project, low=DEFAULT_LOW, high=DEFAULT_HIGH, step=DEFAULT_STEP,
             band=None, receivers=None, verbose=True):
    """受音点ごとの伝達関数と、固有周波数の目安を出す。→ 結果の辞書"""
    import atmosphere as at
    import loop_noredundancy as ln

    air = at.Atmosphere(temperature=project.temperature,
                        humidity=project.humidity, pressure=project.pressure)
    frequencies = np.arange(float(low), float(high) + 0.5 * step, float(step))

    model = None
    try:
        model = _model_of(project)
    except Exception as error:
        if verbose:
            print(f"[伝達関数] モデルが読めません（{type(error).__name__}: {error}）。"
                  f"固有周波数は出しません")

    # 受音点の指定は番号でも `rec3` でも受ける
    wanted = None
    if receivers:
        wanted = {int(str(value).lower().replace("rec", "")) for value in receivers}

    result = {"frequencies": frequencies, "receivers": {}, "modes": [],
              "sound_velocity": air.sound_velocity}

    if model is not None:
        low_corner, high_corner = model.extents
        size = np.asarray(high_corner) - np.asarray(low_corner)
        result["size"] = size
        result["volume"] = model.volume
        result["modes"] = modal_frequencies(size, air.sound_velocity, high, low)
        edge = 4.0 * float(np.sum(size))
        result["density"] = mode_density(frequencies, model.volume or 0.0,
                                         model.surface_area or 0.0, edge,
                                         air.sound_velocity)

    index = 0
    missing = 0
    while index < 1000:
        index += 1
        if wanted is not None and index not in wanted:
            if index > max(wanted):
                break
            continue
        sub = pj.Project(project.folder,
                         **{k: getattr(project, k) for k in pj.DEFAULTS})
        sub.receiver_index = index
        path = sub.existing_result_path("pulses")
        if path is None or not os.path.exists(path):
            missing += 1
            if missing > 3 and result["receivers"]:
                break
            continue
        missing = 0
        pulses = ln.PulseList.from_csv(path)
        result["receivers"][f"rec{index}"] = {
            "response": response(pulses, frequencies, air, band=band),
            "resolution": resolution(pulses, air),
            "pulses": len(pulses.distance),
            "rt": _reverberation_time(project, index),
        }
    if result["receivers"] and result.get("volume"):
        times = [data["rt"] for data in result["receivers"].values()
                 if data["rt"]]
        if times:
            result["schroeder"] = schroeder_frequency(result["volume"],
                                                      float(np.median(times)))
    if verbose:
        _report(result)
    return result


def _report(result):
    print(f"[伝達関数] 受音点 {len(result['receivers'])} 点 / "
          f"{result['frequencies'][0]:.0f}〜{result['frequencies'][-1]:.0f} Hz / "
          f"刻み {result['frequencies'][1] - result['frequencies'][0]:.2f} Hz")
    for name, data in list(result["receivers"].items())[:3]:
        print(f"[伝達関数]   {name}: パルス {data['pulses']} 本 / "
              f"出せる分解能の目安 {data['resolution']:.2f} Hz")
    if result.get("modes"):
        modes = result["modes"]
        kinds = {}
        for _f, _x, _y, _z, kind in modes:
            kinds[kind] = kinds.get(kind, 0) + 1
        print(f"[伝達関数] 固有周波数（直方体とみなした目安）: {len(modes)} 本"
              f"（{'・'.join(f'{k} {v}' for k, v in kinds.items())}）"
              f" 最低 {modes[0][0]:.1f} Hz")
    if result.get("schroeder") and np.isfinite(result["schroeder"]):
        print(f"[伝達関数] シュレーダー周波数の目安: "
              f"{result['schroeder']:.0f} Hz（これより下はモードが個別に見える）")


# ---- 出力 ------------------------------------------------------------------

def write_csv(project, result, verbose=True):
    """`結果/recN/<室>_<条件>_伝達関数.csv`（周波数・大きさ・位相）。"""
    written = []
    for name, data in result["receivers"].items():
        sub = pj.Project(project.folder,
                         **{k: getattr(project, k) for k in pj.DEFAULTS})
        sub.receiver_index = int(name.replace("rec", ""))
        path = os.path.join(os.path.dirname(sub.result_path("pulses")),
                            sub.prefixed(FILE_NAME))
        os.makedirs(os.path.dirname(path), exist_ok=True)
        pressure = data["response"]
        peak = float(np.max(np.abs(pressure))) or 1.0
        with io.open(path, "w", encoding="utf-8-sig", newline="") as handle:
            handle.write("周波数_Hz,大きさ_dB,最大を0dBとした_dB,位相_deg\n")
            for frequency, value in zip(result["frequencies"], pressure):
                magnitude = abs(value)
                handle.write(f"{frequency:.3f},"
                             f"{20.0 * math.log10(magnitude + 1e-300):.4f},"
                             f"{20.0 * math.log10(magnitude / peak + 1e-300):.4f},"
                             f"{math.degrees(math.atan2(value.imag, value.real)):.3f}\n")
        written.append(path)
    if verbose and written:
        print(f"[伝達関数] 表を {len(written)} 点ぶん書き出しました"
              f"（例: {written[0]}）")
    return written


def write_modes_csv(project, result, verbose=True):
    """`結果/<室>_固有周波数.csv`（直方体とみなした目安）。"""
    if not result.get("modes"):
        return None
    path = os.path.join(project.folder, pj.RESULT_DIR,
                        project.prefixed("固有周波数.csv"))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    size = result.get("size")
    with io.open(path, "w", encoding="utf-8-sig", newline="") as handle:
        handle.write(f"# 直方体とみなした目安（室 {size[0]:.2f} × {size[1]:.2f}"
                     f" × {size[2]:.2f} m、音速 {result['sound_velocity']:.1f} m/s）\n")
        handle.write("周波数_Hz,nx,ny,nz,種類\n")
        for value, nx, ny, nz, kind in result["modes"]:
            handle.write(f"{value:.2f},{nx},{ny},{nz},{kind}\n")
    if verbose:
        print(f"[伝達関数] 固有周波数の一覧: {path}")
    return path


def write_figures(project, result, receivers=None, verbose=True):
    """受音点ごとの図（伝達関数＋固有周波数の重ね描き）。"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    import plots as pl

    pl.use_japanese_font()
    frequencies = result["frequencies"]
    written = []
    for name, data in result["receivers"].items():
        if receivers and name not in receivers:
            continue
        pressure = np.abs(data["response"])
        level = 20.0 * np.log10(pressure / (np.max(pressure) or 1.0) + 1e-300)

        figure, axis = plt.subplots(figsize=(11.0, 4.6))
        # ★線を引きすぎると読めないので、**軸モードは全部・接線は
        #   シュレーダー周波数まで・斜めは引かない**
        limit = result.get("schroeder")
        limit = float(limit) if limit and np.isfinite(limit) else frequencies[-1]
        drawn = 0
        for value, nx, ny, nz, kind in result.get("modes", []):
            if kind == "斜め":
                continue
            if kind == "接線" and value > limit:
                continue
            axis.axvline(value, color=MODE_COLOR,
                         linewidth=0.8 if kind == "軸" else 0.5,
                         alpha=0.55 if kind == "軸" else 0.25)
            drawn += 1
        axis.plot(frequencies, level, color="#2b6cb0", linewidth=1.0)
        if result.get("schroeder") and np.isfinite(result["schroeder"]):
            axis.axvline(result["schroeder"], color="#2f855a", linewidth=1.4,
                         linestyle="--",
                         label=f"シュレーダー周波数 {result['schroeder']:.0f} Hz")
        axis.plot([], [], color=MODE_COLOR, linewidth=0.8,
                  label=f"固有周波数（直方体とみなした目安・軸と接線 {drawn} 本）")
        axis.set_xlim(frequencies[0], frequencies[-1])
        axis.set_ylim(max(-60.0, float(np.percentile(level, 1)) - 5.0), 3.0)
        axis.set_xlabel("周波数 [Hz]")
        axis.set_ylabel("伝達関数の大きさ [dB]（最大を 0 dB）")
        axis.set_title(f"{project.room_label}／{project.condition_label}　{name}"
                       f"　伝達関数（位相ごと重ねた複素和）\n"
                       f"分解能の目安 {data['resolution']:.2f} Hz"
                       f"（パルス {data['pulses']} 本）", fontsize=10)
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8, loc="lower right")
        figure.tight_layout()
        path = project.figure_path(f"伝達関数_{name}.png")
        figure.savefig(path, dpi=140)
        plt.close(figure)
        written.append(path)
    if verbose and written:
        print(f"[伝達関数] 図を {len(written)} 枚書き出しました"
              f"（例: {written[0]}）")
    return written


def run(project, low=DEFAULT_LOW, high=DEFAULT_HIGH, step=DEFAULT_STEP,
        band=None, receivers=None, verbose=True):
    result = evaluate(project, low=low, high=high, step=step, band=band,
                      receivers=receivers, verbose=verbose)
    if not result["receivers"]:
        print("[伝達関数] パルス列が見つかりません（先に計算してください）")
        return result
    write_csv(project, result, verbose=verbose)
    write_modes_csv(project, result, verbose=verbose)
    write_figures(project, result, verbose=verbose)
    return result


def main(argv=None):
    import argparse

    import condition_table as ct

    parser = argparse.ArgumentParser(
        description="伝達関数（周波数特性）と固有周波数を出す（計算はやり直さない）")
    parser.add_argument("folder", help="プロジェクトのフォルダ")
    parser.add_argument("--sheet", default=None, help="条件（シート名）")
    parser.add_argument("--low", type=float, default=DEFAULT_LOW)
    parser.add_argument("--high", type=float, default=DEFAULT_HIGH)
    parser.add_argument("--step", type=float, default=DEFAULT_STEP)
    parser.add_argument("--band", type=int, default=None,
                        help="使う帯域の番号（省略すると全帯域の合計）")
    parser.add_argument("--receiver", type=int, action="append", default=None,
                        help="受音点の番号（繰り返し指定可。省略すると全部）")
    args = parser.parse_args(argv)

    base = pj.Project.load(args.folder)
    sheets = ([args.sheet] if args.sheet
              else (ct.sheets(ct.path(base)) if ct.exists(base) else [None]))
    for sheet in sheets:
        project = pj.Project(base.folder,
                             **{k: getattr(base, k) for k in pj.DEFAULTS})
        if sheet:
            project.condition_sheet = sheet
        print(f"[伝達関数] 条件『{project.condition_label}』")
        run(project, low=args.low, high=args.high, step=args.step,
            band=args.band,
            receivers=[f"rec{n}" for n in args.receiver] if args.receiver else None)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
