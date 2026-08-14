import numpy as np
from mesh import Mesh
import read_dxffile as rd
import sound_ray as sr
import loop_reflectionmesh as lr
import loop_deleteredundancy as ld
import loop_noredundancy as ln
import impulse as ir
import reverberation as rv
import absorption as ab
from atmosphere import Atmosphere
from ray_recorder import RayRecorder


###
###全体の流れをここに記述
###
### 目的
### 音源位置、受音点、室形状、音線法に使う受音球の半径、反射回数、音線の数
### を元に（def process()内の順に記載）
### パルス列、または、インパルス応答（wavファイル）を保存する
def process(soundsource_point, reciever_point, dxf_filename, sphere_radius, nref, soundray_number,
            absorption_csv=None, absorption_kind=None, layer_assignment=None,
            band_number=rd.DEFAULT_BAND_NUMBER, material_library=None,
            unit=None, orient_normals="cad",
            atmosphere=None, raylog_filename=None, raylog_max_rays=2000,
            pulse_filename=None, impulse_filename=None,
            sampling_frequency=ir.SAMPLING_FREQUENCY, max_time=ir.MAX_TIME,
            reverberation_filename=None, decay_filename=None):
    """
    閉じた室でも、一面だけの壁のような**開いた形状**でも計算できる
    （当たる壁がなくなった音線はそこで打ち切られる）。

    soundsource_point / reciever_point : (3,) | None
        None を渡すと DXF の src / rec レイヤの POINT から取る。
        CAD 側で音源・受音点まで作図しておけば、ここでの指定は不要。
    absorption_csv : str | None
        吸音率テーブルの CSV パス（1列目=レイヤ名または ID、以降=バンド別吸音率）。
        None ならレイヤ全部が既定値になり警告が出る。data/absorption_sample.csv を参照。
    absorption_kind : 'normal' | 'random' | None
        CSV の吸音率が**垂直入射**か**残響室法（乱入射）**か。
        残響室法なら Paris の式で垂直入射に変換してから使う。
        None ならファイル冒頭の `# kind:` 宣言を見る。それも無ければ垂直入射として扱い警告する。
        ★取り違えると吸音を大きく誤るので、出典を確認して明示すること。
    layer_assignment : absorption.LayerAssignment | dict | None
        CAD のレイヤ名 → 材料名の対応。**CAD を編集せずに材料を差し替えるための仕組み。**
        None なら「レイヤ名 = 材料名」とみなす。
    band_number : int
        周波数バンド数。既定 8（63〜8k Hz）。63 と 8k を外すなら 6（125〜4k Hz）。
    material_library : absorption.MaterialLibrary | None
        材料一覧を直接渡す場合（GUI から編集したものなど）。
        指定すると absorption_csv より優先される。
    atmosphere : atmosphere.Atmosphere | None
        温度・湿度・気圧。**音速と空気吸収の両方がここから決まる。**
        None なら基準状態（20℃ / 湿度 40% / 101.325 kPa → 音速 343.8 m/s）。
        元コードの 340.0 m/s を再現したい場合は Atmosphere(temperature=14.0)。
    unit : None | str | float
        None なら DXF ヘッダの $INSUNITS から自動判定。'mm' / 'm' などで明示指定も可。
    orient_normals : str
        法線の向きの扱い。既定 'cad' は CAD の巻き順をそのまま使う
        （法線が空気側を向くようにモデルを作るのは CAD 側の責任）。
        'flip' で全反転、'shells' でシェル単位の自動補正。
        詳細は read_dxffile.read_model() を参照。
    raylog_filename : str | None
        可視化用の音線軌跡を npz で保存するパス。None なら記録しない。
        出力①音線の可視化・②音粒子の可視化（動画）に使う。docs/出力・可視化方針.md 参照。
    raylog_max_rays : int
        記録する音線の本数の上限。総数がこれを超えたら等間隔に間引く。
    pulse_filename : str | None
        バックトレースが出すパルス列（到来時刻・到来方向・バンド別エネルギー）の
        CSV 出力先。None なら書かない。出力⑤伝搬方向の素材になる。
    impulse_filename : str | None
        インパルス応答の CSV 出力先。None ならインパルス応答の合成自体を行わない
        （音線追跡だけ見たいときは None にすると速い）。
    sampling_frequency / max_time : float
        インパルス応答のサンプリング周波数と長さ（元コード fs / tmax）。
    reverberation_filename / decay_filename : str | None
        残響指標（EDT / T20 / T30）・減衰曲線の CSV 出力先。
        どちらかを指定すると算出する（インパルス応答の合成が前提）。

    戻り値:
        dict … 'model' / 'pulses' / 'impulse' / 'reverberation'
               （計算しなかったものは None）
    """
    # オクターブバンドの中心周波数。既定は 8 バンド（63〜8k Hz）。
    frequencies = ab.octave_bands(band_number)

    if atmosphere is None:
        atmosphere = Atmosphere()
    sound_velocity = atmosphere.sound_velocity
    print(f"[procedure] {atmosphere.summary()}")

    # 吸音率テーブルを作る。
    # ・残響室法の値なら Paris の式で垂直入射に変換してから渡す
    # ・レイヤ → 材料の対応は layer_assignment で差し替えられる（CAD を触らずに済む）
    absorption_table = None
    if material_library is None and absorption_csv is not None:
        material_library = ab.MaterialLibrary.from_csv(absorption_csv,
                                                       kind=absorption_kind)
    if material_library is not None:
        print(f"[procedure] {material_library.summary()}")
        absorption_table = material_library.absorption_table(layer_assignment,
                                                             band_number=band_number)

    # 室形状・吸音率・音源・受音点をまとめて DXF から読む
    # 元コード132〜283行目に対応
    model = rd.read_model(dxf_filename, unit=unit, absorption_table=absorption_table,
                          orient_normals=orient_normals, band_number=band_number)
    mesh = model.mesh

    if soundsource_point is None:
        if not model.source_points:
            raise ValueError("音源が指定されておらず、DXF に src レイヤの POINT もありません")
        soundsource_point = model.source_points[0]
    if reciever_point is None:
        if not model.receiver_points:
            raise ValueError("受音点が指定されておらず、DXF に rec レイヤの POINT もありません")
        reciever_point = model.receiver_points[0]

    # 音線ベクトルを作成
    soundray_list = sr.soundray_generator(soundray_number)

    # 可視化用の軌跡レコーダ（本線の計算には影響しない副チャンネル）
    recorder = None
    if raylog_filename is not None:
        recorder = RayRecorder(total_rays=soundray_number, max_rays=raylog_max_rays,
                               sound_velocity=sound_velocity,
                               band_number=len(frequencies))

    # 音線ループで反射面のIDを履歴として記録します
    # 元コード524行目に対応
    reflection_history = lr.loop(soundsource_point, reciever_point, soundray_list, nref, mesh,
                                 sphere_radius, recorder=recorder)

    if recorder is not None:
        print("音線軌跡:", recorder.summary())
        recorder.save_npz(raylog_filename)

    # 重複経路の削除
    # 元コード721行目に対応
    reflection_history = ld.delete(reflection_history)

    # 非重複経路　バックトレース（虚音源法）
    # 元コード876行目に対応。
    # 吸音率は Mesh が面ごとに持っているので、ここで別途渡す必要はない。
    pulses = ln.loop(soundsource_point, reciever_point, reflection_history, mesh,
                     sound_velocity=sound_velocity, band_number=len(frequencies),
                     filename=pulse_filename)

    # 後部残響が nref で切れていないかの確認。
    # 残響時間は「エネルギーが 35 dB 減衰するまで」を見るので、
    # そこに達する前に反射回数の上限で打ち切られていると過小評価になる。
    #
    # ただし「nref に達した経路がある」だけでは警告にならない。そこまで反射した音は
    # たいてい十分弱まっているからである。そこで**受音点に届くエネルギー**
    # （バンド別エネルギー / 距離^2。インパルス応答の振幅の 2 乗に比例する量）で見て、
    # 打ち切り時点の音がまだ大きいときだけ警告する。
    if len(pulses):
        top = int(pulses.reflection_count.max())
        received = pulses.energy.max(axis=1) / pulses.distance ** 2
        at_limit = received[pulses.reflection_count == top].max()
        drop_db = 10.0 * np.log10(at_limit / received.max())
        if top >= nref and drop_db > -35.0:
            print(f"[procedure] 警告: 最大反射回数 nref={nref} で打ち切られた音が"
                  f"まだ {drop_db:.1f} dB までしか減衰していません。"
                  f"後部残響が途中で切れているので残響時間は過小評価になります。"
                  f"（35 dB 以上減衰するまで nref を増やしてください。"
                  f"吸音率が小さい部屋ほど大きな nref が要ります）")
        elif top >= nref:
            print(f"[procedure] 最大反射回数 nref={nref} に達した経路がありますが、"
                  f"その時点で {drop_db:.1f} dB まで減衰しているので残響時間には影響しません")

    # 離散的なパルス列を、バンドパスフィルタを通して自然なインパルス応答に変換する
    # 元コード make_ipls_freq_monaural_fortran.f90 に対応
    impulse = None
    if impulse_filename is not None:
        if len(pulses) == 0:
            print("[procedure] 受音した経路が無いのでインパルス応答は作れません")
        else:
            impulse = ir.impulse_responce(
                impulse_filename, pulses, octave_frequencies=frequencies,
                atmosphere=atmosphere, sampling_frequency=sampling_frequency,
                max_time=max_time)

    # 残響時間の算出。元コード ipls2rt_fortran.f90 に対応
    reverberation = None
    if reverberation_filename is not None or decay_filename is not None:
        if impulse is None:
            print("[procedure] インパルス応答が無いので残響時間は算出できません"
                  "（impulse_filename を指定してください）")
        else:
            reverberation = rv.reverberation_time(
                impulse[0], impulse[1], rt_filename=reverberation_filename,
                decay_filename=decay_filename, frequencies=frequencies)

    return {"model": model, "pulses": pulses, "impulse": impulse,
            "reverberation": reverberation}


def main():
    """コマンドラインから通し実行する。

        cd geosim
        python procedure.py ..\\test.dxf --absorption ..\\absorption.csv --out ..\\結果
    """
    import argparse
    import os

    p = argparse.ArgumentParser(description="幾何音響シミュレーションの通し実行")
    p.add_argument("dxf", help="室形状の DXF ファイル")
    p.add_argument("--absorption", help="吸音率 CSV")
    p.add_argument("--absorption-kind", choices=["normal", "random"],
                   help="吸音率が垂直入射(normal)か残響室法(random)か。"
                        "省略時はファイルの # kind: 宣言を見る")
    p.add_argument("--assignment", help="レイヤ→材料の対応 JSON（CAD を触らずに差し替える）")
    p.add_argument("--bands", type=int, default=8, choices=[6, 8],
                   help="周波数バンド数（8=63〜8kHz / 6=125〜4kHz）")
    p.add_argument("--out", default=".", help="出力先フォルダ（既定: カレント）")
    p.add_argument("--rays", type=int, default=20000, help="音線の本数")
    p.add_argument("--nref", type=int, default=8, help="最大反射回数")
    p.add_argument("--radius", type=float, default=0.15, help="受音球の半径 [m]")
    p.add_argument("--temperature", type=float, default=None, help="気温 [℃]（既定 20）")
    p.add_argument("--humidity", type=float, default=None, help="相対湿度 [%%]（既定 40）")
    p.add_argument("--pressure", type=float, default=None, help="気圧 [kPa]（既定 101.325）")
    p.add_argument("--max-time", type=float, default=ir.MAX_TIME,
                   help="インパルス応答の長さ [s]")
    p.add_argument("--unit", help="'mm' / 'm' など。省略すると $INSUNITS から自動判定")
    p.add_argument("--orient-normals", default="cad", choices=["cad", "flip", "shells"])
    p.add_argument("--source", type=float, nargs=3, metavar=("X", "Y", "Z"),
                   help="音源座標 [m]。省略すると DXF の src レイヤから取る")
    p.add_argument("--receiver", type=float, nargs=3, metavar=("X", "Y", "Z"),
                   help="受音点座標 [m]。省略すると DXF の rec レイヤから取る")
    p.add_argument("--no-impulse", action="store_true",
                   help="インパルス応答・残響時間を計算しない（音線追跡だけ見たいとき）")
    a = p.parse_args()

    os.makedirs(a.out, exist_ok=True)
    base = os.path.splitext(os.path.basename(a.dxf))[0]

    def out(suffix):
        return os.path.join(a.out, f"{base}_{suffix}")

    air_kwargs = {k: v for k, v in (("temperature", a.temperature),
                                    ("humidity", a.humidity),
                                    ("pressure", a.pressure)) if v is not None}
    assignment = (ab.LayerAssignment.load(a.assignment) if a.assignment else None)

    process(soundsource_point=a.source, reciever_point=a.receiver,
            dxf_filename=a.dxf, sphere_radius=a.radius, nref=a.nref,
            soundray_number=a.rays, absorption_csv=a.absorption,
            absorption_kind=a.absorption_kind, layer_assignment=assignment,
            band_number=a.bands, unit=a.unit, orient_normals=a.orient_normals,
            atmosphere=Atmosphere(**air_kwargs),
            raylog_filename=out("raylog.npz"),
            pulse_filename=out("pulses.csv"),
            impulse_filename=None if a.no_impulse else out("ir.csv"),
            max_time=a.max_time,
            reverberation_filename=None if a.no_impulse else out("rt.csv"),
            decay_filename=None if a.no_impulse else out("decay.csv"))


if __name__ == "__main__":
    main()
