import numpy as np
import sound_ray as sr
import mesh_method as mm
from geosim.sound_ray import reflection_generator
from mesh import Mesh


# 元コード721 重複経路の削除に該当
# 何をしているか分からない
def delete():
    # 反射履歴traceffn=0の部分以外はnref+1としてインデックス範囲外の値を代入する意図？
    # 初めからtraceffnを作成するようにはできないのか
    # do i = 1, count
    #     do j = 0, nref
    #         if (traceff(i, j).eq. 0) then
    #             traceffn(i) = j
    #             exit
    #         end if
    #         traceffn(i) = nref + 1
    #     end do
    # end do


    #もし重複だけ除くのであれば以下の記述を採用します
    list_dummy = [[0, 1, 2], [0, 1, 1], [0, 1, 2], [0, 2, 3]]
    print("重複あり")
    print(list_dummy)
    tuple_dummy = [tuple(data) for data in list_dummy]
    tuple_dummy = set(tuple_dummy)
    list_dummy = [list(data) for data in tuple_dummy]
    print("重複なし")
    print(list_dummy)

    return

if __name__ == '__main__':
    delete()
