import numpy as np
import sound_ray as sr
import mesh_method as mm
from mesh import Mesh


# 7/9打ち合わせ用




# 元コード721 重複経路の削除に該当
# 何をしているか分からないでいます
# 想像として、反射履歴の配列を重複のない配列にできれば問題ないように見えるが
# 実際にどのようなことをしているか分からない
# 以下のコードは
# 重複ありの反射履歴　reflectionhistory_redundancy　を
# 重複なしの反射履歴　reflection_history　に変更するデモ
def delete(reflectionhistory_redundancy):
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

    # もし重複だけ除くのであれば以下の記述を採用します
    # reflectionhistory_redundancy
    # reflectionhistory_redundancy = [[0, 1, 2], [0, 1, 1], [0, 1, 2], [0, 2, 3]]
    print("重複あり")
    print(reflectionhistory_redundancy)
    tuple_dummy = [tuple(data) for data in reflectionhistory_redundancy]
    tuple_dummy = set(tuple_dummy)
    reflection_history = [list(data) for data in tuple_dummy]
    print("重複なし")
    print(reflection_history)

    return reflection_history


if __name__ == '__main__':
    test = [[0, 1, 2], [0, 1, 1], [0, 1, 2], [0, 2, 3]]
    delete(test)
