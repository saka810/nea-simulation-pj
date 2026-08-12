import numpy as np
from mesh import Mesh
import read_dxffile as rd
import sound_ray as sr
import loop_reflectionmesh as lr
import loop_deleteredundancy as ld
import loop_noredundancy as ln
import impulse as ir
from ray_recorder import RayRecorder


###
###全体の流れをここに記述
###
### 目的
### 音源位置、受音点、室形状、音線法に使う受音球の半径、反射回数、音線の数
### を元に（def process()内の順に記載）
### パルス列、または、インパルス応答（wavファイル）を保存する
def process(soundsource_point, reciever_point, dxf_filename, sphere_radius, nref, soundray_number,
            raylog_filename=None, raylog_max_rays=2000):
    """
    raylog_filename : str | None
        可視化用の音線軌跡を npz で保存するパス。None なら記録しない。
        出力①音線の可視化・②音粒子の可視化（動画）に使う。docs/出力・可視化方針.md 参照。
    raylog_max_rays : int
        記録する音線の本数の上限。総数がこれを超えたら等間隔に間引く。
    """
    frequencies = np.array([63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0])
    print(frequencies)

    # コマンドプロンプトで入力、ＵＩから読み込み、別ｃｓｖファイルから読み込み
    # このファイルに直接書き込みなど
    # 音源、受音点、形状、「吸音率のリスト」が必要
    # soundsource_point=
    # reciever_point=
    mesh = rd.read(dxf_filename)
    # absorption_list

    # 音線ベクトルを作成
    soundray_list = sr.soundray_generator(soundray_number)

    # 可視化用の軌跡レコーダ（本線の計算には影響しない副チャンネル）
    recorder = None
    if raylog_filename is not None:
        recorder = RayRecorder(total_rays=soundray_number, max_rays=raylog_max_rays)

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

    # 非重複経路　バックトレース
    # 元コード876行目に対応
    # 吸音率のリストはいつ読み込むか
    # ここでcsvファイルに計算結果を書き込んだら終了
    # ln.loop(soundsource_point,reciever_point,absorption_list,reflection_history,mesh)

    # この後おそらく離散的なパルス列をより自然なインパルス応答に変換しているものと推測
    # ir.impulse_responce(filename, sound_velocity, reflection_timing, soundsourceenergy_list,
    #                     frequency_number, count,
    #                     nn, mf, fmax, dt)

    # インパルス応答に変換した後の残響時間を算出するプログラムは他ソフトでもよい？
