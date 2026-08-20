"""パイプライン全体の通し実行。**DXF を読んで残響時間まで一気に出す。**

`geosim` の各モジュールを順に呼ぶだけの層で、計算そのものは各モジュールにある。
GUI（`app.py` → `run_project.py`）もコマンドラインもここを通る。

    ① DXF 読込          read_dxffile.read_model()      元 backtrace.f90 132〜283 行
    ② 音線生成          sound_ray.soundray_generator() 元 318〜326 行
    ③ 音線追跡          loop_reflectionmesh.loop()     元 524〜717 行
    ④ 重複経路の削除    loop_deleteredundancy.loop()   元 721〜841 行
    ⑤ バックトレース    loop_noredundancy.loop()       元 876〜1134 行
    ⑥ インパルス応答    impulse.impulse_response_from_pulses()     元 make_ipls_freq_monaural_fortran.f90
    ⑦ 残響時間          reverberation.reverberation_time()  元 ipls2rt_fortran.f90

統計残響式（Sabine / Eyring / Eyring-Knudsen）は音線を飛ばす前に出せるので、
②より前に計算して結果と並べて表示する。**どちらが正しいという話ではなく**、
食い違いが大きいときに設定を疑う手がかりにする。

出力先を渡さない引数（`*_filename=None`）はその段階を飛ばす。
"""


import numpy as np
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
def process(soundsource_point, receiver_point, dxf_filename, sphere_radius, nref, soundray_number,
            absorption_csv=None, absorption_kind=None, layer_assignment=None,
            band_number=rd.DEFAULT_BAND_NUMBER, material_library=None,
            unit=None, orient_normals="cad", two_sided=False, volume=None,
            atmosphere=None, raylog_filename=None, raylog_max_rays=2000,
            pulse_filename=None, impulse_filename=None,
            sampling_frequency=ir.SAMPLING_FREQUENCY, max_time=ir.MAX_TIME,
            reverberation_filename=None, decay_filename=None,
            statistical_filename=None, statistical=True,
            clarity=True, clarity_filename=None, surface_filename=None,
            flip_faces=None, face_materials=None, traced_history=None,
            progress=None):
    """
    閉じた室でも、一面だけの壁のような**開いた形状**でも計算できる
    （当たる壁がなくなった音線はそこで打ち切られる）。

    soundsource_point / receiver_point : (3,) | None
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
        **'inward' は面ごとにレイの偶奇で室内側へ揃える**。面のつながりを要求しないので、
        床・壁・天井を 1 枚ずつ描いた「板の寄せ集め」モデルはこれで直る。
        詳細は read_dxffile.read_model() を参照。
    two_sided : bool
        面の**裏からの入射も当てる**か（既定 False）。
        既定では法線の側から来た音線しか当たらないので、法線が逆を向いた壁はすり抜ける。
        CAD で床・壁・天井を 1 枚ずつ描いた「板の寄せ集め」モデルは巻き順や押し出し方向で
        法線がまちまちになりがちで、そのままだと壁が抜けて音線が室外へ逃げる。
        そういうモデルは True にする。詳細は mesh_method.FaceArrays。
    volume : float | None
        統計残響式に使う室容積 [m³]。None なら DXF が閉じている場合に自動算出する。
        **閉じていないモデル（板の寄せ集めなど）では自動算出できない**ので、
        統計式と比べたければここで与える（例: 床面積 × 天井高）。
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
    statistical : bool
        Sabine / Eyring / Eyring-Knudsen の統計残響式でも残響時間を出すか（既定 True）。
        **室が閉じている場合のみ**（容積が決まらないと計算できない）。
        音線を飛ばさず面積と吸音率だけから出るので、シミュレーション結果の物差しになる。
    statistical_filename : str | None
        統計残響式の結果の CSV 出力先。
    clarity : bool
        明瞭度系の指標（C50 / C80 / D50 / Ts）も出すか（既定 True）。
        残響時間が「どれだけ長く響くか」なのに対し、こちらは
        「初期の音が後から来る音に対してどれだけ強いか」を見る。会議室・教室で効く。
    clarity_filename / surface_filename : str | None
        明瞭度の指標／材料別の面積と吸音率の CSV 出力先。
    flip_faces : iterable[int] | None
        法線を反転する面インデックス（`face_editor.py` で目で見て直したぶん）。
        自動判定のあとに重ねて適用される。
    face_materials : dict[int, str] | None
        **面ごとの吸音材の割り当て** {面インデックス: 材料名}（同じく `face_editor.py`）。
        レイヤで吸音材を分けられないモデル（1 つの 3DSOLID で出来ていて
        面ごとのレイヤが無いなど）で使う。指定した面はレイヤより優先される。
    traced_history : list[list[int]] | None
        **すでに済ませた音線追跡の結果**（この受音点ぶんの反射面 ID 履歴）。
        渡すと音線の生成と追跡を飛ばす。受音点が複数あるとき、追跡は受音点に
        依らないので `run_project` が 1 回だけ回して各受音点に配る（F-6）。
    progress : callable(段階名: str, 割合: float|None) | None
        **進み具合の通知先**（GUI の進捗表示用）。割合は 0〜1、分からない段階は None。
        重い段階（音線追跡・バックトレース）は途中でも何度か呼ばれる。
        渡さなければ何もしないので、本線の計算には影響しない。

    戻り値:
        dict … 'model' / 'pulses' / 'impulse' / 'reverberation' / 'statistical'
               （計算しなかったものは None）
    """
    def report(stage, fraction=None):
        if progress is not None:
            progress(stage, fraction)

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
    report("モデルを読み込み中")
    model = rd.read_model(dxf_filename, unit=unit, absorption_table=absorption_table,
                          orient_normals=orient_normals, band_number=band_number,
                          flip_faces=flip_faces, face_materials=face_materials)
    mesh = model.mesh

    # 作図ミスの洗い出し（TODO B-10）。計算に入る前に指摘するほうが早い
    issues = rd.check_model(model, absorption_table=absorption_table)
    if any(i["level"] == "error" for i in issues):
        print("[procedure] ★ エラーがあります。結果は当てになりません")

    if soundsource_point is None:
        if not model.source_points:
            raise ValueError("音源が指定されておらず、DXF に src レイヤの POINT もありません")
        soundsource_point = model.source_points[0]
    if receiver_point is None:
        if not model.receiver_points:
            raise ValueError("受音点が指定されておらず、DXF に rec レイヤの POINT もありません")
        receiver_point = model.receiver_points[0]

    # 統計残響式（Sabine / Eyring / Eyring-Knudsen）。音線を飛ばす前に出せる。
    # 面積と吸音率だけから決まるので、あとの計算結果と突き合わせる物差しになる
    report("統計残響式")
    statistical_result = None
    if statistical:
        if volume is not None:
            print(f"[統計残響] 容積は指定値 {volume:.2f} m3 を使います"
                  f"（DXF からの自動算出はしません）")
            statistical_result = rv.statistical_reverberation(
                mesh, volume, frequencies=frequencies, atmosphere=atmosphere)
        else:
            statistical_result = rv.statistical_reverberation_from_model(
                model, frequencies=frequencies, atmosphere=atmosphere)
        if statistical_result is not None and statistical_filename is not None:
            rv.write_statistical_reverberation(statistical_filename, statistical_result)
            print(f"[統計残響] 書き出しました: {statistical_filename}")
        if statistical_result is not None and surface_filename is not None:
            import project as pj
            pj.write_surface_csv(surface_filename, statistical_result["surface"],
                                 frequencies)
            print(f"[統計残響] 材料別の面積・吸音率: {surface_filename}")

    if traced_history is None:
        # 音線ベクトルを作成
        report("音線を生成中")
        soundray_list = sr.soundray_generator(soundray_number)

        # 可視化用の軌跡レコーダ（本線の計算には影響しない副チャンネル）
        recorder = None
        if raylog_filename is not None:
            recorder = RayRecorder(total_rays=soundray_number, max_rays=raylog_max_rays,
                                   sound_velocity=sound_velocity,
                                   band_number=len(frequencies))

        # 音線ループで反射面のIDを履歴として記録します
        # 元コード524行目に対応
        report("音線追跡", 0.0)
        reflection_history = lr.loop(soundsource_point, receiver_point, soundray_list,
                                     nref, mesh, sphere_radius, recorder=recorder,
                                     two_sided=two_sided,
                                     progress=lambda f: report("音線追跡", f))

        if recorder is not None:
            print("音線軌跡:", recorder.summary())
            recorder.save_npz(raylog_filename)
    else:
        # ★受音点が複数あるとき、追跡は**呼び出し側で 1 回だけ**済ませてある
        #   （追跡は受音点に依らない。`loop_reflectionmesh.loop` の受音判定を参照）
        report("音線追跡（済み・受音点間で共有）", 1.0)
        reflection_history = traced_history

    # 重複経路の削除
    # 元コード721行目に対応
    report("重複経路の削除")
    reflection_history = ld.delete(reflection_history)

    # 非重複経路　バックトレース（虚音源法）
    # 元コード876行目に対応。
    # 吸音率は Mesh が面ごとに持っているので、ここで別途渡す必要はない。
    pulses = ln.loop(soundsource_point, receiver_point, reflection_history, mesh,
                     sound_velocity=sound_velocity, band_number=len(frequencies),
                     filename=pulse_filename, two_sided=two_sided,
                     progress=lambda f: report("バックトレース", f))

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
            report("インパルス応答の合成")
            impulse = ir.impulse_response_from_pulses(
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
            report("残響時間の算出")
            reverberation = rv.reverberation_time(
                impulse[0], impulse[1], rt_filename=reverberation_filename,
                decay_filename=decay_filename, frequencies=frequencies)

    # 明瞭度系の指標（C50 / C80 / D50 / Ts）。残響時間とは見ている中身が違う。
    # 「どれだけ長く響くか」ではなく「初期の音が後から来る音に対してどれだけ強いか」
    clarity_result = None
    if clarity and impulse is not None:
        report("明瞭度の指標")
        clarity_result = rv.clarity_measures(impulse[0], impulse[1],
                                            frequencies=frequencies)
        if clarity_filename is not None:
            rv.write_clarity_measures(clarity_filename, clarity_result)
            print(f"[procedure] 明瞭度の指標を書き出しました: {clarity_filename}")

    # 統計残響式との突き合わせ。**どちらが正しいという話ではない**。
    # 統計式は拡散音場を前提にした平均像、シミュレーションは特定の受音点での実際の減衰。
    # 大きく食い違うときは、音場が拡散していないか設定に問題があるかの手がかりになる
    if statistical_result is not None and reverberation is not None:
        # **周波数は横**（table.py の共通ルール）
        measured = reverberation["measures"]["T30"]
        knudsen = statistical_result["eyring_knudsen"]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(knudsen > 0, measured / knudsen, np.nan)
        rows = [("T30(計算)", measured)]
        rows += [(label, statistical_result[key])
                 for key, label in rv.STATISTICAL_LABELS.items()]
        rows += [("比 T30/E-K", ratio)]
        print("[procedure] " + " " * 14
              + "".join(f"{f:>10.0f}" for f in frequencies))
        for name, values in rows:
            print(f"[procedure] {name:>14}"
                  + "".join("       ---" if np.isnan(v) else f"{v:10.3f}"
                            for v in values))

    return {"model": model, "pulses": pulses, "impulse": impulse,
            "reverberation": reverberation, "statistical": statistical_result,
            "clarity": clarity_result, "frequencies": frequencies,
            "atmosphere": atmosphere,
            "soundsource_point": np.asarray(soundsource_point, dtype=float),
            "receiver_point": np.asarray(receiver_point, dtype=float)}


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
    p.add_argument("--orient-normals", default="cad",
                   choices=["cad", "flip", "shells", "inward"],
                   help="法線の扱い。cad=そのまま / flip=全反転 / shells=シェル単位 / "
                        "inward=面ごとにレイの偶奇で室内側へ揃える"
                        "（CAD で法線を意識せずに描いたモデル用）")
    p.add_argument("--two-sided", action="store_true",
                   help="面の裏からの入射も当てる。CAD の法線がまちまちな"
                        "「板の寄せ集め」モデル用（既定は表からのみ）")
    p.add_argument("--volume", type=float, default=None,
                   help="統計残響式に使う室容積 [m3]。"
                        "閉じていないモデルで統計式と比べたいときに指定する")
    p.add_argument("--source", type=float, nargs=3, metavar=("X", "Y", "Z"),
                   help="音源座標 [m]。省略すると DXF の src レイヤから取る")
    p.add_argument("--receiver", type=float, nargs=3, metavar=("X", "Y", "Z"),
                   help="受音点座標 [m]。省略すると DXF の rec レイヤから取る")
    p.add_argument("--no-impulse", action="store_true",
                   help="インパルス応答・残響時間を計算しない（音線追跡だけ見たいとき）")
    p.add_argument("--no-statistical", action="store_true",
                   help="Sabine / Eyring-Knudsen の統計残響式を計算しない")
    a = p.parse_args()

    os.makedirs(a.out, exist_ok=True)
    base = os.path.splitext(os.path.basename(a.dxf))[0]

    def out(suffix):
        return os.path.join(a.out, f"{base}_{suffix}")

    air_kwargs = {k: v for k, v in (("temperature", a.temperature),
                                    ("humidity", a.humidity),
                                    ("pressure", a.pressure)) if v is not None}
    assignment = (ab.LayerAssignment.load(a.assignment) if a.assignment else None)

    process(soundsource_point=a.source, receiver_point=a.receiver,
            dxf_filename=a.dxf, sphere_radius=a.radius, nref=a.nref,
            soundray_number=a.rays, absorption_csv=a.absorption,
            absorption_kind=a.absorption_kind, layer_assignment=assignment,
            band_number=a.bands, unit=a.unit, orient_normals=a.orient_normals,
            two_sided=a.two_sided, volume=a.volume,
            atmosphere=Atmosphere(**air_kwargs),
            raylog_filename=out("raylog.npz"),
            pulse_filename=out("pulses.csv"),
            impulse_filename=None if a.no_impulse else out("ir.csv"),
            max_time=a.max_time,
            reverberation_filename=None if a.no_impulse else out("rt.csv"),
            decay_filename=None if a.no_impulse else out("decay.csv"),
            statistical=not a.no_statistical,
            statistical_filename=None if a.no_statistical else out("rt_statistical.csv"))


if __name__ == "__main__":
    main()
