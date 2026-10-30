! ifort backtrace_rev.f90 /O3 /Qopenmp   （BACKTRACE\build_ifort.bat でビルド）
!
! backtrace.f90 の修正版（2026-10-19）。旧版からの変更点：
!  1. 音線追跡を1回にした。受音した反射履歴は各スレッドが自分の作業領域に積み、最後に1つにまとめる。
!     （旧版は、件数を数える並列の追跡と、履歴を書き込む逐次の追跡を2回行っており、
!       2回の件数が食い違うと配列の範囲外に書き込むおそれがあった）
!  2. 反射履歴に三角形の番号ではなく平面の番号を記録する。同一平面上の三角形は同じ虚音源を作るので、
!     逆追跡で隣の三角形に当たっても同じ平面なら有効とする。反射率は実際に当たった三角形の材料で計算する。
!     （旧版は、受音球の中心からの経路が記録と別の三角形に当たると、正しい虚音源を捨てていた）
!  3. 音線数は、計算する時間長 tmax と、時刻 tmax の虚音源1つあたりに受音球へ届く音線数の期待値 nhitmin から、
!     原著論文（Vorlander 1989）の式(9)で決める。届かずに見落とす確率は exp(-nhitmin)。
!     到達時間が tmax を超える経路は追跡も出力もしない。
!  4. 重複経路の削除は、反射履歴を辞書式に整列（マージソート）して隣どうしを比べる方法に変えた。
!     （旧版は反射回数ごとの総当たりで、計算量が件数の2乗に比例した）
!  5. 出力は1スレッドで書き、行の順序を実行ごとに変わらないようにした（旧版は並列に書き込んでいた。順序は 8. を参照）。
!  6. 音線方向のフィボナッチ格子の式を z = (2i-1)/N - 1 に直した（旧版は (2i+1)/N - 1 で、最後の1本が |z|>1 となりNaNだった）。
!  7. 反射回数の上限 nref を自動で決める（2026-10-21）。方向を均等に散らした nsample 本の音線を、伝搬距離が c0*tmax に
!     達するまで反射回数の制限なしで追跡し、その最大反射回数に nrefmrg を掛けた値を上限とする。
!     本番の追跡で上限に達した音線があれば警告する（旧版は nref を定数で与えており、小さい室では tmax まで届かなかった）。
!  8. 音線を組に分け、組ごとに「音線追跡→重複経路の削除→バックトレース」を行って有効な経路だけを残す（2026-10-21）。
!     必要メモリは音線の総数ではなく1組の大きさで決まる。組の数は、試し追跡で数えた1本あたりの受音件数から、
!     1組の処理が memgb [GB] 程度に収まるように決める。組 b には b, b+組の数, b+2*組の数, … 番目の音線を割り当てる
!     （どの組も全方向をまんべんなく含む）。組をまたいで同じ経路が見つかった場合は、結果（11列）が完全に一致するので、
!     最後に到達時間の順に並べ替えて一致する行を除く。出力は到達時間の順になる。

program main

!$ use omp_lib

implicit none

real*8,			parameter	:: tmax		= 3.0d0				! 計算する時間長[s]
real*8,			parameter	:: nhitmin	= 30.0d0			! 時刻tmaxの虚音源1つあたりに受音球へ届く音線数の期待値（見落とす確率はexp(-nhitmin)）
real*8,			parameter	:: srcx		= -1.2d0			! 音源x座標[m]
real*8,			parameter	:: srcy		=  0.0d0			! 音源y座標[m]
real*8,			parameter	:: srcz		=  1.1d0			! 音源z座標[m]
real*8,			parameter	:: rcvx		= -1.25d0			! 受音点x座標[m]
real*8,			parameter	:: rcvy		=  0.0d0			! 受音点y座標[m]
real*8,			parameter	:: rcvz		=  1.4d0			! 受音点z座標[m]
real*8,			parameter	:: rcvr		=  1.0d0			! 受音球半径[m]
real*8,			parameter	:: c0		= 340.0d0			! 音速
integer*4,		parameter	:: nsample	= 100000			! 反射回数の上限を決める試し追跡の音線数
real*8,			parameter	:: nrefmrg	= 10.0d0			! 試し追跡で得た最大反射回数に掛ける余裕
integer*4,		parameter	:: nrefcap	= 1000000			! 試し追跡で1本の音線を追う反射回数の打ち切り（無限ループよけ）
real*8,			parameter	:: memgb	= 32.0d0			! 1組の処理に使うメモリの目安[GB]（組の数を決める）
real*8,			parameter	:: tolnm	= 1.0d-10			! 同一平面の判定：1-(単位法線ベクトルの内積)の許容値
real*8,			parameter	:: told		= 1.0d-4			! 同一平面の判定：平面の方程式ax+by+cz+d=0のdの差の許容値[m]
character*128,	parameter	:: fnabs	= 'absorption_mag_normal.csv'	! 吸音率リストファイル名
character*128,	parameter	:: fndxf	= '00_mag.dxf'		! 形状ファイル名
character*128,	parameter	:: fnout	= '00_mag_stage.txt'	! 出力ファイル名
character*1,	parameter	:: ynshade	= 'y'				! shadeで書き出したかどうか（y: shade, n: shade以外）
character*1,	parameter	:: ynnmrev	= 'n'				! 法線が外向きかどうか（y: 外向き（反転する）, n: 内向き（反転しない））
character*1,	parameter	:: ynmeter	= 'y'				! dxfの単位がメートルかどうか（y: メートル, n: ミリメートル）

type histbuf
	integer*4,	dimension(:,:),	allocatable	:: h	! 反射履歴(反射順,通し番号)
	integer*4								:: n	! 格納数
end type histbuf

integer*4,		dimension(:,:),	allocatable	:: surf	! 面情報(使用する吸音材番号,頂点番号)
real*8,			dimension(:,:),	allocatable	:: vwpt	! 頂点座標(頂点番号,xyz)
real*8,			dimension(:,:),	allocatable	:: vwnm	! 全壁面の単位法線ベクトル(壁面番号,xyz)

integer*4,		dimension(:),	allocatable	:: num0						! 素材番号（読み込み用）
real*8,			dimension(:),	allocatable	:: a1, a2, a3, a4, a5, a6	! 各周波数の垂直入射吸音率（読み込み用）
character*128,	dimension(:),	allocatable	:: sozai					! 素材名
real*8,			dimension(:,:),	allocatable	:: absper					! 吸音率

character*128,	dimension(:),	allocatable	:: grpno	! dxfファイル読み込みグループナンバー
character*128,	dimension(:),	allocatable	:: enti		! dxfファイル読み込みエンティティナンバー

integer*4,		dimension(:),	allocatable	:: pid		! 各壁面が属する平面の番号(壁面番号)
real*8,			dimension(:,:),	allocatable	:: pnm		! 平面の単位法線ベクトル(平面番号,xyz)
real*8,			dimension(:),	allocatable	:: pd		! 平面の方程式ax+by+cz+d=0のd(平面番号)

type(histbuf),	dimension(:),	allocatable	:: buf		! スレッドごとの受音反射履歴(スレッド番号)
integer*4,		dimension(:,:),	allocatable	:: traceff	! 組の受音反射履歴(反射順,経路番号)
integer*4,		dimension(:),	allocatable	:: idx		! 整列用の経路番号

real*8,			dimension(:,:),	allocatable	:: res		! 組の受音結果(反射回数・受音時間・虚音源方向xyz・6帯域のエネルギー,経路番号)
logical,		dimension(:),	allocatable	:: valid	! 組の受音結果が有効かどうか(経路番号)
real*8,			dimension(:,:),	allocatable	:: resall	! 全組の有効な受音結果(11列,通し番号)
real*8,			dimension(:,:),	allocatable	:: restmp	! resallを広げるときの作業用

integer*4,		dimension(:),	allocatable	:: tractmp	! 各音線用の一時反射履歴（平面番号）(反射順)
real*8,			dimension(:,:),	allocatable	:: isrc		! 虚音源リスト(反射順,xyz)
real*8		:: enertmp(6)		! 各経路用の一時エネルギー

integer*4	:: i, j, k, l, m	! ループ変数
integer*4	:: nray				! 音線数
integer*4	:: nref				! 反射回数の上限（試し追跡で決める）
integer*4	:: nsampl			! 試し追跡の音線数
integer*4	:: kmax, nesc0		! 試し追跡の最大反射回数、室外へ抜けた音線数
integer*4	:: ncap				! 試し追跡で反射回数の打ち切りに達した音線数
real*8		:: ksum				! 試し追跡の反射回数の合計
real*8		:: nrayr			! 音線数（実数）
real*8		:: dt, dz			! 方位角（角度）・仰角（z座標）それぞれの離散化幅
real*8		:: pi				! 円周率
real*8		:: distr			! 音線と受音点の距離
real*8		:: distd			! 基点と受音点の距離
real*8		:: temp				! 作業用
real*8		:: temp2			! 平面の方程式の定数項
integer*4	:: flag				! 「壁面に当たる」フラグ
integer*4	:: sfcount			! 壁面数
integer*4	:: ptcount			! 壁面の頂点数
integer*4	:: plcount			! 平面数
integer*4	:: nthr, tid		! スレッド数、スレッド番号
integer*4	:: ntrunc			! 反射回数の上限で打ち切った音線数
integer*4	:: nescape			! 壁面に当たらず室外へ抜けた音線数
integer*4	:: nbatch, ib		! 組の数、組の番号
integer*4	:: nrb, jr, ic		! 組の音線数、組の中の音線番号、候補の経路番号
integer*4	:: ncap0			! 組の作業配列の大きさ（受音件数）
integer*4	:: nall				! 全組の有効な受音結果の件数（組をまたいだ重複を含む）
integer*8	:: ndet, nuniq		! 受音した経路数（重複を含む）、組ごとの非重複経路数の合計
integer*8	:: nraydone			! 処理を終えた音線数
real*8		:: ndet0			! 試し追跡で受音球を通った回数
real*8		:: rate, bpd		! 音線1本あたりの受音件数、受音1件あたりのメモリ[バイト]
integer*4	:: nvalid			! 有効な経路数
real*8		:: lmax				! tmaxに対応する伝搬距離[m]
real*8		:: plen				! 音源からの伝搬距離[m]
real*8		:: dhit, thit		! 最も手前の壁面までの距離、そのときの直線のパラメータt

! 一時変数
integer*4		:: count, countred, jtmp, ktmp
real*8			:: node(3), disttmp, costmp
real*8			:: vray(3), vini(3), vtgt(3)
character*128	:: entitmp, strtmp
integer*4		:: surfid, pointid, absid, vertexflag, ptindex(3), pointidtmp
real*8			:: pointtmp(3), vectmp(3,3), revtmp
character*8		:: date
character*10	:: time
character*5		:: zone
integer			:: values(8)

!■■■プログラム開始時刻■■■
call date_and_time(date, time, zone, values)
write(*,'(i4,"/",i2,"/",i2," ",i2,":",i2,":",i2,".",i3)') values(1), values(2), values(3), values(5), values(6), values(7), values(8)


!■■■吸音率リストの読み込み■■■
write(*,'(a,$)') '吸音率リストの読み込み中'
open(1, file = trim(fnabs), status = 'old')
count	= 0
do
	read(1,*,end=100)
	count	= count + 1
end do
100 continue
allocate(num0(count))
allocate(sozai(count))
allocate(a1(count))
allocate(a2(count))
allocate(a3(count))
allocate(a4(count))
allocate(a5(count))
allocate(a6(count))
allocate(absper(count,6))
rewind(1)
do i = 1, count
	read(1,*) num0(i), sozai(i), a1(i), a2(i), a3(i), a4(i), a5(i), a6(i)
	absper(num0(i),1)	= a1(i)
	absper(num0(i),2)	= a2(i)
	absper(num0(i),3)	= a3(i)
	absper(num0(i),4)	= a4(i)
	absper(num0(i),5)	= a5(i)
	absper(num0(i),6)	= a6(i)
end do
close(1)
write(*,*) '...完了'


!■■■室形状の読み込み■■■
write(*,'(a,$)') '室形状の読み込み'

!■形状ファイルオープン■
open(1, file = trim(fndxf), status = 'old')

!■行数カウント■
count	= 0
do
	read(1,*,end=200)
	count	= count + 1
end do
200 continue

!■配列の確保■
allocate(grpno(count/2))
allocate(enti(count/2))

!■頂点数、面数のカウント■
rewind(1)
do i = 1, count / 2
	read(1,'(a)') strtmp
	grpno(i)	= trim(adjustl(strtmp))
	read(1,'(a)') strtmp
	enti(i)		= trim(adjustl(strtmp))
end do
close(1)
ptcount	= 0
sfcount	= 0
do i = 1, count / 2
	if((grpno(i) .eq. '70') .and. (enti(i) .eq. '192')) ptcount	= ptcount + 1
	if((grpno(i) .eq. '70') .and. (enti(i) .eq. '128')) sfcount	= sfcount + 1
end do

!■配列の確保■
allocate(vwpt(ptcount,3))
allocate(surf(sfcount,4))
allocate(vwnm(sfcount,3))

!■形状データ読み込み■
i		= 1
surfid	= 0
pointid	= 0
do
	if((grpno(i) .eq. '0') .and. (enti(i) .eq. 'EOF')) then
		exit
	elseif((grpno(i) .eq. '0') .and. (enti(i) .eq. 'POLYLINE')) then
		do
			i	= i + 1
			if(grpno(i) .eq. '8') then
				entitmp	= enti(i)
				if(ynshade .eq. 'y') then
					read(entitmp(index(entitmp,'_')+1:len(entitmp)),*) absid
				else
					read(entitmp,*) absid
				end if
			elseif(grpno(i) .eq. '70') then
				if(enti(i) .ne. '64') then
					exit
				else
					pointidtmp	= pointid
					i	= i + 1
					do
						if(grpno(i) .eq. '0') then
							if(enti(i) .eq. 'SEQEND') then
								exit
							elseif(enti(i) .eq. 'VERTEX') then
								do
									i	= i + 1
									if(grpno(i) .eq. '10') then
										read(enti(i),*) pointtmp(1)
									elseif(grpno(i) .eq. '20') then
										read(enti(i),*) pointtmp(2)
									elseif(grpno(i) .eq. '30') then
										read(enti(i),*) pointtmp(3)
									elseif(grpno(i) .eq. '70') then
										read(enti(i),*) vertexflag
									elseif(grpno(i) .eq. '71') then
										read(enti(i),*) ptindex(1)
									elseif(grpno(i) .eq. '72') then
										read(enti(i),*) ptindex(2)
									elseif(grpno(i) .eq. '73') then
										read(enti(i),*) ptindex(3)
									elseif(grpno(i) .eq. '0') then
										if(vertexflag .eq. 128) then
											surfid	= surfid + 1
											surf(surfid,1)	= absid
											surf(surfid,2)	= abs(ptindex(1)) + pointidtmp
											surf(surfid,3)	= abs(ptindex(2)) + pointidtmp
											surf(surfid,4)	= abs(ptindex(3)) + pointidtmp
											do l = 1, 3
												vectmp(1,l)	= vwpt(surf(surfid,3),l) - vwpt(surf(surfid,2),l)
												vectmp(2,l)	= vwpt(surf(surfid,4),l) - vwpt(surf(surfid,2),l)
											end do
											if((vectmp(1,2) * vectmp(2,3) - vectmp(1,3) * vectmp(2,2)) ** 2.0d0 &
											 + (vectmp(1,3) * vectmp(2,1) - vectmp(1,1) * vectmp(2,3)) ** 2.0d0 &
											 + (vectmp(1,1) * vectmp(2,2) - vectmp(1,2) * vectmp(2,1)) ** 2.0d0 .eq. 0.0d0) then
												surfid	= surfid - 1
												sfcount	= sfcount - 1
											end if
										elseif(vertexflag .eq. 192) then
											pointid	= pointid + 1
											if(ynmeter .eq. 'y') then
												vwpt(pointid,1)	= pointtmp(1)
												vwpt(pointid,2)	= pointtmp(2)
												vwpt(pointid,3)	= pointtmp(3)
											else
												vwpt(pointid,1)	= pointtmp(1) / 1000.0d0
												vwpt(pointid,2)	= pointtmp(2) / 1000.0d0
												vwpt(pointid,3)	= pointtmp(3) / 1000.0d0
											end if
										end if
										exit
									end if
								end do
							end if
						else
							i	= i + 1
						end if
					end do
					exit
				end if
			end if
		end do
		i	= i + 1
	else
		i	= i + 1
	end if
end do

!■法線の算出■
do j = 1, sfcount
	do l = 1, 3
		vectmp(1,l)	= vwpt(surf(j,3),l) - vwpt(surf(j,2),l)
		vectmp(2,l)	= vwpt(surf(j,4),l) - vwpt(surf(j,2),l)
	end do
	vwnm(j,1)	= vectmp(1,2) * vectmp(2,3) - vectmp(1,3) * vectmp(2,2)
	vwnm(j,2)	= vectmp(1,3) * vectmp(2,1) - vectmp(1,1) * vectmp(2,3)
	vwnm(j,3)	= vectmp(1,1) * vectmp(2,2) - vectmp(1,2) * vectmp(2,1)
	disttmp = sqrt(vwnm(j,1) ** 2.0d0 + vwnm(j,2) ** 2.0d0 + vwnm(j,3) ** 2.0d0)
	do l = 1, 3
		vwnm(j,l)	= vwnm(j,l) / disttmp
	end do
end do
if(ynnmrev .eq. 'y') then
	do j = 1, sfcount
		do l = 1, 3
			vwnm(j,l)	= -vwnm(j,l)
		end do
	end do
end if
write(*,*) '...完了'


!■■■同一平面の判定■■■
! 同一平面上の三角形は同じ虚音源を作るので、反射履歴には平面番号を記録する
write(*,'(a,$)') '同一平面の判定'
allocate(pid(sfcount))
allocate(pnm(sfcount,3))
allocate(pd(sfcount))
plcount	= 0
do j = 1, sfcount

	! 平面の方程式ax+by+cz+d=0のdを算出
	temp2	= -(vwnm(j,1) * vwpt(surf(j,2),1) + vwnm(j,2) * vwpt(surf(j,2),2) + vwnm(j,3) * vwpt(surf(j,2),3))

	! 既出の平面と法線・dが一致すれば同じ平面番号
	pid(j)	= 0
	do m = 1, plcount
		if((1.0d0 - (vwnm(j,1) * pnm(m,1) + vwnm(j,2) * pnm(m,2) + vwnm(j,3) * pnm(m,3)) .le. tolnm) &
			.and. (abs(temp2 - pd(m)) .le. told)) then
			pid(j)	= m
			exit
		end if
	end do

	! 一致しなければ新しい平面
	if(pid(j) .eq. 0) then
		plcount	= plcount + 1
		pid(j)	= plcount
		do l = 1, 3
			pnm(plcount,l)	= vwnm(j,l)
		end do
		pd(plcount)	= temp2
	end if

end do
write(*,*) '...完了'
write(*,'(a,i8,a,i8)') '  壁面数', sfcount, '  平面数', plcount


!■■■音線数の決定■■■
! 原著論文の式(8)：距離 c*t の虚音源から半径rcvrの受音球に届く音線数の期待値は nray*rcvr**2/(4*(c*t)**2)
! これが時刻tmaxでnhitminになるように音線数を決める（式(9)）
pi		= acos(-1.0d0)
lmax	= c0 * tmax
nrayr	= 4.0d0 * nhitmin * lmax ** 2.0d0 / rcvr ** 2.0d0
if(nrayr .gt. dble(huge(nray))) then
	write(*,*) '音線数が多すぎます。tmax、nhitmin、rcvrを見直してください'
	stop
end if
nray	= ceiling(nrayr)
write(*,'(a,i12)') '  音線数', nray
write(*,'(a,f10.3,a)') '  届く音線数の期待値が1本になる時刻', sqrt(dble(nray) * rcvr ** 2.0d0 / 4.0d0) / c0, ' s'


!■■■反射回数の上限の決定■■■
! 方向を均等に散らしたnsampl本の音線を、伝搬距離がc0*tmax（受音球の半径の分の余裕を含む）に達するまで、
! 反射回数の制限なしで追跡する（受音の記録はしない）。その最大反射回数にnrefmrgを掛けた値を上限とする。
! あわせて受音球を通った回数を数え、組の数を決めるのに使う
write(*,'(a,$)') '反射回数の上限の決定'
nsampl	= min(nsample, nray)
dt		= pi * (3.0d0 - sqrt(5.0d0))
dz		= 2.0d0 / dble(nsampl)
kmax	= 0
ksum	= 0.0d0
nesc0	= 0
ncap	= 0
ndet0	= 0.0d0
!$omp parallel do schedule(dynamic,100) private(i,k,l,vray,disttmp,vini,plen,flag,jtmp,dhit,thit,node,revtmp,vtgt,distd,distr) reduction(max: kmax) reduction(+: ksum, nesc0, ncap, ndet0)
do i = 1, nsampl

	!■音線ベクトルの作成（本番と同じ球面上のフィボナッチ格子）■
	vray(3)	= dz * dble(i) - 1.0d0 - 1.0d0 / dble(nsampl)
	vray(1)	= sqrt(1.0d0 - vray(3) ** 2.0d0) * cos(dt * dble(i))
	vray(2)	= sqrt(1.0d0 - vray(3) ** 2.0d0) * sin(dt * dble(i))
	disttmp	= sqrt(vray(1) ** 2.0d0 + vray(2) ** 2.0d0 + vray(3) ** 2.0d0)
	do l = 1, 3
		vray(l)	= vray(l) / disttmp
	end do
	vini(1)	= srcx
	vini(2)	= srcy
	vini(3)	= srcz
	plen	= 0.0d0

	!■反射回数ループ（本番の追跡と同じ条件で反射を数える）■
	k	= 0
	do
		call nearest_wall(vini, vray, flag, jtmp, dhit, thit)

		! 受音球を通るかどうか（本番の追跡と同じ判定）
		vtgt(1)	= rcvx - vini(1)
		vtgt(2)	= rcvy - vini(2)
		vtgt(3)	= rcvz - vini(3)
		distd	= vray(1) * vtgt(1) + vray(2) * vtgt(2) + vray(3) * vtgt(3)
		distr	= sqrt((vini(1) + vray(1) * distd - rcvx) ** 2.0d0 &
					 + (vini(2) + vray(2) * distd - rcvy) ** 2.0d0 &
					 + (vini(3) + vray(3) * distd - rcvz) ** 2.0d0)
		if ((distr .le. rcvr) .and. (distd .le. dhit) .and. (distd .ge. 0.0d0) .and. (plen + distd - rcvr .le. lmax)) then
			ndet0	= ndet0 + 1.0d0
		end if

		if(flag .eq. 0) then
			nesc0	= nesc0 + 1
			exit
		end if
		plen	= plen + dhit
		if(plen - rcvr .gt. lmax) exit
		if(k .ge. nrefcap) then
			ncap	= ncap + 1
			exit
		end if
		k	= k + 1
		do l = 1, 3
			node(l)	= vini(l) + thit * vray(l)
		end do
		vini	= node
		revtmp	= vray(1) * vwnm(jtmp,1) + vray(2) * vwnm(jtmp,2) + vray(3) * vwnm(jtmp,3)
		do l = 1, 3
			vray(l)	= vray(l) - 2.0d0 * revtmp * vwnm(jtmp,l)
		end do
		disttmp	= sqrt(vray(1) ** 2.0d0 + vray(2) ** 2.0d0 + vray(3) ** 2.0d0)
		do l = 1, 3
			vray(l)	= vray(l) / disttmp
		end do
	end do

	kmax	= max(kmax, k)
	ksum	= ksum + dble(k)

end do
nref	= max(1, ceiling(nrefmrg * dble(kmax)))
write(*,*) '...完了'
write(*,'(a,i12)') '  試し追跡の音線数', nsampl
write(*,'(a,f12.1,a,i10)') '  反射回数の平均', ksum / dble(nsampl), '  最大', kmax
write(*,'(a,i12,a,i10,a)') '  反射回数の上限', nref, '  （反射履歴1件の大きさ', 4 * nref, ' バイト）'
if(ncap .gt. 0) then
	write(*,'(a,i10,a,i10,a)') '  警告：試し追跡で', ncap, ' 本の音線が反射回数', nrefcap, ' 回に達しました（音線が閉じ込められている可能性）'
end if
if(nesc0 .gt. 0) then
	write(*,'(a,i10,a)') '  試し追跡で室外へ抜けた音線数', nesc0, '（形状データのすき間）'
end if

!■組の数の決定（1組の処理が memgb [GB] 程度に収まるようにする）■
! 受音1件あたりのメモリ：反射履歴（4*nref バイト）がスレッドの作業領域と組の配列の2か所にあり、作業領域は最大2倍まで広がる
rate	= ndet0 / dble(nsampl)
bpd		= 3.0d0 * 4.0d0 * dble(nref) + 120.0d0
nbatch	= max(1, min(nray, ceiling(dble(nray) * rate * bpd / (memgb * 1.0d9))))
write(*,'(a,f12.2)') '  音線1本あたりの受音件数', rate
write(*,'(a,es12.3)') '  受音件数の見込み（重複を含む）', dble(nray) * rate
write(*,'(a,i12,a,i12,a,f8.1,a)') '  組の数', nbatch, '  1組の音線数', (nray - 1) / nbatch + 1, &
	'  1組のメモリの見込み', dble((nray - 1) / nbatch + 1) * rate * bpd / 1.0d9, ' GB'


!■■■音線追跡・重複経路の削除・バックトレース（組ごと）■■■
write(*,'(a)') '音線追跡・重複経路の削除・バックトレース（組ごとに行う）'
allocate(tractmp(nref))
allocate(isrc(0:nref,3))

!■スレッドごとの受音反射履歴の確保（1組分の見込みの1.2倍。足りなければ広げる）■
nthr	= 1
!$ nthr	= omp_get_max_threads()
allocate(buf(0:nthr-1))
do m = 0, nthr - 1
	allocate(buf(m)%h(nref, max(1024, ceiling(1.2d0 * rate * dble((nray - 1) / nbatch + 1) / dble(nthr)))))
	buf(m)%n	= 0
end do

!■組の作業配列（受音件数がこれを超えたら広げる）と、全組の有効な受音結果■
ncap0	= max(1024, ceiling(1.2d0 * rate * dble((nray - 1) / nbatch + 1)))
allocate(traceff(nref,ncap0))
allocate(idx(ncap0))
allocate(res(11,ncap0))
allocate(valid(ncap0))
nall	= 0
allocate(resall(11,1024))

!■球面上のフィボナッチ格子の離散化幅■
dt	= pi * (3.0d0 - sqrt(5.0d0))
dz	= 2.0d0 / dble(nray)

ntrunc		= 0
nescape		= 0
ndet		= 0
nuniq		= 0
nraydone	= 0
do ib = 1, nbatch

	!■■音線追跡（組 ib の音線 i = ib, ib+nbatch, ib+2*nbatch, …）■■
	nrb	= (nray - ib) / nbatch + 1
	do m = 0, nthr - 1
		buf(m)%n	= 0
	end do
	!$omp parallel do schedule(dynamic,1000) private(jr,i,k,l,tid,tractmp,vray,disttmp,vini,plen,flag,jtmp,dhit,thit,node,vtgt,distd,distr,revtmp) reduction(+: ntrunc, nescape)
	do jr = 1, nrb

	i	= ib + (jr - 1) * nbatch

	tid	= 0
	!$ tid	= omp_get_thread_num()

	!■一時反射履歴の初期化■
	do k = 1, nref
		tractmp(k)	= 0
	end do

	!■直接音音線ベクトルの作成■
	vray(3)	= dz * dble(i) - 1.0d0 - 1.0d0 / dble(nray)
	vray(1)	= sqrt(1.0d0 - vray(3) ** 2.0d0) * cos(dt * dble(i))
	vray(2)	= sqrt(1.0d0 - vray(3) ** 2.0d0) * sin(dt * dble(i))

	!■直接音音線ベクトルの単位ベクトル化■
	disttmp	= sqrt(vray(1) ** 2.0d0 + vray(2) ** 2.0d0 + vray(3) ** 2.0d0)
	do l = 1, 3
		vray(l)	= vray(l) / disttmp
	end do

	!■直接音音線ベクトルの基点（音源位置）の初期化■
	vini(1)	= srcx
	vini(2)	= srcy
	vini(3)	= srcz
	plen	= 0.0d0

	!■反射回数ループ■
	do k = 0, nref

		! 最も手前の壁面を探す
		call nearest_wall(vini, vray, flag, jtmp, dhit, thit)

		! 音線基点・受音点間ベクトル
		vtgt(1)	= rcvx - vini(1)
		vtgt(2)	= rcvy - vini(2)
		vtgt(3)	= rcvz - vini(3)

		! 「基点」と「受音点から音線ベクトルにおろした垂線の足」との距離
		distd	= vray(1) * vtgt(1) + vray(2) * vtgt(2) + vray(3) * vtgt(3)	! 内積

		! 「受音点」と「受音点から音線ベクトルにおろした垂線の足」との距離
		distr	= sqrt((vini(1) + vray(1) * distd - rcvx) ** 2.0d0 &
					 + (vini(2) + vray(2) * distd - rcvy) ** 2.0d0 &
					 + (vini(3) + vray(3) * distd - rcvz) ** 2.0d0)

		! 音線と受音点の距離が受音半径以内 且つ 壁面までの距離以内 且つ 音線ベクトルと基点・受音点間ベクトルの内積が正
		! 且つ tmax以内に届き得るとき（受音）
		if ((distr .le. rcvr) .and. (distd .le. dhit) .and. (distd .ge. 0.0d0) .and. (plen + distd - rcvr .le. lmax)) then

			! スレッドの受音反射履歴に追加
			call append_hist(buf(tid), tractmp)

		end if

		! 壁面に当たらない場合は反射回数ループを抜ける
		if(flag .eq. 0) then
			nescape	= nescape + 1
			exit
		end if

		! 次の区間の始点がtmaxに対応する距離を超えたら打ち切る
		plen	= plen + dhit
		if(plen - rcvr .gt. lmax) exit

		! tmax以内なのに反射回数の上限に達したら打ち切る
		if(k .eq. nref) then
			ntrunc	= ntrunc + 1
			exit
		end if

		! 一時反射履歴に平面番号を書き込み
		tractmp(k+1)	= pid(jtmp)

		! 交点の座標を算出
		do l = 1, 3
			node(l)	= vini(l) + thit * vray(l)
		end do

		! 次回の基点 = 交点座標
		vini(1)	= node(1)
		vini(2)	= node(2)
		vini(3)	= node(3)

		! 反射音線のベクトル
		revtmp	= vray(1) * vwnm(jtmp,1) + vray(2) * vwnm(jtmp,2) + vray(3) * vwnm(jtmp,3)
		vray(1)	= vray(1) - 2.0d0 * revtmp * vwnm(jtmp,1)
		vray(2)	= vray(2) - 2.0d0 * revtmp * vwnm(jtmp,2)
		vray(3)	= vray(3) - 2.0d0 * revtmp * vwnm(jtmp,3)

		! 反射音線ベクトルの単位ベクトル化
		disttmp	= sqrt(vray(1) ** 2.0d0 + vray(2) ** 2.0d0 + vray(3) ** 2.0d0)
		do l = 1, 3
			vray(l)	= vray(l) / disttmp
		end do

	end do

	end do

	!■スレッドごとの受音反射履歴を組の配列にまとめる■
	count	= 0
	do m = 0, nthr - 1
		count	= count + buf(m)%n
	end do
	if(count .gt. ncap0) then
		ncap0	= ceiling(1.2d0 * dble(count))
		deallocate(traceff, idx, res, valid)
		allocate(traceff(nref,ncap0))
		allocate(idx(ncap0))
		allocate(res(11,ncap0))
		allocate(valid(ncap0))
	end if
	count	= 0
	do m = 0, nthr - 1
		do i = 1, buf(m)%n
			do k = 1, nref
				traceff(k,count+i)	= buf(m)%h(k,i)
			end do
		end do
		count	= count + buf(m)%n
	end do
	ndet	= ndet + count

	!■■重複経路の削除（反射履歴を辞書式に整列し、隣どうしを比べる）■■
	do i = 1, count
		idx(i)	= i
	end do
	call merge_sort(traceff, idx, count, .false.)
	countred	= 0
	do i = 1, count
		if(countred .eq. 0) then
			countred	= 1
			idx(countred)	= idx(i)
		elseif(histcmp(traceff(:,idx(i)), traceff(:,idx(countred))) .ne. 0) then
			countred	= countred + 1
			idx(countred)	= idx(i)
		end if
	end do
	nuniq	= nuniq + countred

	!■■バックトレース計算（候補 i の反射履歴は traceff(:,idx(i))）■■
	!$omp parallel do schedule(dynamic,100) private(i,ic,k,l,ktmp,isrc,enertmp,temp,vray,disttmp,vini,flag,jtmp,dhit,thit,node,distr,vtgt,costmp)
	do i = 1, countred

	ic	= idx(i)

	valid(i)	= .false.

	!■一時エネルギーの初期化■
	do k = 1, 6
		enertmp(k)	= 1.0d0
	end do

	!■虚音源リストの作成■

	! 音源位置
	isrc(0,1)	= srcx
	isrc(0,2)	= srcy
	isrc(0,3)	= srcz

	! 反射回数ループ（平面に対して前回の（虚）音源を鏡像に移す）
	ktmp	= 0
	do k = 1, nref
		if(traceff(k,ic) .eq. 0) exit
		ktmp	= k
		temp	= pnm(traceff(k,ic),1) * isrc(k-1,1) + pnm(traceff(k,ic),2) * isrc(k-1,2) &
				+ pnm(traceff(k,ic),3) * isrc(k-1,3) + pd(traceff(k,ic))
		do l = 1, 3
			isrc(k,l)	= isrc(k-1,l) - 2.0d0 * temp * pnm(traceff(k,ic),l)
		end do
	end do

	!■最後の反射：受音点から最後の虚音源へ向かう音線■
	vray(1)	= isrc(ktmp,1) - rcvx
	vray(2)	= isrc(ktmp,2) - rcvy
	vray(3)	= isrc(ktmp,3) - rcvz
	disttmp	= sqrt(vray(1) ** 2.0d0 + vray(2) ** 2.0d0 + vray(3) ** 2.0d0)
	do l = 1, 3
		vray(l)	= vray(l) / disttmp
	end do
	vini(1)	= rcvx
	vini(2)	= rcvy
	vini(3)	= rcvz

	!■バックトレースループ■
	do k = ktmp, 0, -1

		! 最も手前の壁面を探す
		call nearest_wall(vini, vray, flag, jtmp, dhit, thit)

		if(k .eq. 0) then

			! 音源より手前に壁面があればバックトレースループを抜ける
			distr	= sqrt((srcx - vini(1)) ** 2.0d0 + (srcy - vini(2)) ** 2.0d0 + (srcz - vini(3)) ** 2.0d0)
			if((flag .eq. 1) .and. (dhit .lt. distr)) exit

			! 受音結果を格納
			vtgt(1)	= rcvx - isrc(ktmp,1)
			vtgt(2)	= rcvy - isrc(ktmp,2)
			vtgt(3)	= rcvz - isrc(ktmp,3)
			disttmp	= sqrt(vtgt(1) ** 2.0d0 + vtgt(2) ** 2.0d0 + vtgt(3) ** 2.0d0)
			res(1,i)	= dble(ktmp)
			res(2,i)	= disttmp / c0
			do l = 1, 3
				res(2+l,i)	= -vtgt(l)
			end do
			do l = 1, 6
				res(5+l,i)	= enertmp(l)
			end do
			valid(i)	= (res(2,i) .le. tmax)

		else

			! 壁面に当たらない、または当たった壁面が記録の平面上になければバックトレースループを抜ける
			if(flag .eq. 0) exit
			if(pid(jtmp) .ne. traceff(k,ic)) exit

			! エネルギーを減衰させる（実際に当たった三角形の材料を使う）
			costmp	= abs(vwnm(jtmp,1) * vray(1) + vwnm(jtmp,2) * vray(2) + vwnm(jtmp,3) * vray(3))
			do l = 1, 6
				enertmp(l)	= abs(((1.0d0 + sqrt(1.0d0 - absper(surf(jtmp,1),l))) * costmp - (1.0d0 - sqrt(1.0d0 - absper(surf(jtmp,1),l)))) / ((1.0d0 + sqrt(1.0d0 - absper(surf(jtmp,1),l))) * costmp + (1.0d0 - sqrt(1.0d0 - absper(surf(jtmp,1),l))))) ** 2.0d0 * enertmp(l)	! 斜入射エネルギー反射率
			end do

			! 次回の基点 = 交点座標
			do l = 1, 3
				node(l)	= vini(l) + thit * vray(l)
			end do
			vini(1)	= node(1)
			vini(2)	= node(2)
			vini(3)	= node(3)

			! 反射音線のベクトル（交点から1つ前の虚音源へ）
			vray(1)	= isrc(k-1,1) - node(1)
			vray(2)	= isrc(k-1,2) - node(2)
			vray(3)	= isrc(k-1,3) - node(3)
			disttmp	= sqrt(vray(1) ** 2.0d0 + vray(2) ** 2.0d0 + vray(3) ** 2.0d0)
			do l = 1, 3
				vray(l)	= vray(l) / disttmp
			end do

		end if

	end do

	end do

	!■有効な受音結果を全組の配列に加える（足りなければ2倍に広げる）■
	do i = 1, countred
		if(valid(i)) then
			if(nall .eq. size(resall, 2)) then
				allocate(restmp(11, 2 * size(resall, 2)))
				restmp(:,1:nall)	= resall(:,1:nall)
				call move_alloc(restmp, resall)
			end if
			nall	= nall + 1
			resall(:,nall)	= res(:,i)
		end if
	end do

	!■進捗表示■
	nraydone	= nraydone + nrb
	if(ib .lt. nbatch) then
		write(*,'(f10.2,a,i12,a,i12,a,i8,a,i8,a,a,$)') dble(nraydone) / dble(nray) * 100.0d0, '% ', nraydone, '/', nray, '  （組', ib, '/', nbatch, '）', char(13)
	else
		write(*,'(f10.2,a,i12,a,i12,a,i8,a,i8,a)') 100.0d0, '% ', nraydone, '/', nray, '  （組', ib, '/', nbatch, '）'
	end if

end do
write(*,*) '...完了'
write(*,'(a,i14)') '  受音した経路数（重複を含む）', ndet
write(*,'(a,i14)') '  組ごとの非重複経路数の合計', nuniq
write(*,'(a,i14)') '  反射回数の上限で打ち切った音線数', ntrunc
if(ntrunc .gt. 0) then
	write(*,'(a,i10,a)') '  警告：tmaxに達する前に反射回数の上限', nref, ' 回で打ち切った音線があります。nrefmrgを大きくしてください'
end if
write(*,'(a,i14)') '  壁面に当たらず室外へ抜けた音線数', nescape
do m = 0, nthr - 1
	deallocate(buf(m)%h)
end do
deallocate(buf)
deallocate(traceff, idx, res, valid)
deallocate(tractmp)
deallocate(isrc)


!■■■組をまたいだ重複の削除と書き出し■■■
! 同じ反射履歴のバックトレース結果は11列すべてが一致するので、到達時間の順に並べ替えて一致する行を除く
write(*,'(a,$)') '組をまたいだ重複の削除と書き出し'
allocate(idx(nall))
do i = 1, nall
	idx(i)	= i
end do
call sort_res(resall, idx, nall)
open(1, file = trim(fnout), status = 'replace')
nvalid	= 0
do i = 1, nall
	if(i .gt. 1) then
		if(all(resall(:,idx(i)) .eq. resall(:,idx(i-1)))) cycle
	end if
	nvalid	= nvalid + 1
	write(1,'(11e30.20e3)') (resall(l,idx(i)), l=1,11)
end do
close(1)
deallocate(idx)
deallocate(resall)
write(*,*) '...完了'
write(*,'(a,i12)') '  有効な経路数（可視な虚音源数）', nvalid


!■■■メモリの解放■■■
write(*,'(a,$)') 'メモリの解放'
deallocate(num0)
deallocate(sozai)
deallocate(a1)
deallocate(a2)
deallocate(a3)
deallocate(a4)
deallocate(a5)
deallocate(a6)
deallocate(absper)
deallocate(grpno)
deallocate(enti)
deallocate(vwpt)
deallocate(surf)
deallocate(vwnm)
deallocate(pid)
deallocate(pnm)
deallocate(pd)
write(*,*) '...完了'


!■■■プログラム終了時刻■■■
call date_and_time(date, time, zone, values)
write(*,'(i4,"/",i2,"/",i2," ",i2,":",i2,":",i2,".",i3)') values(1), values(2), values(3), values(5), values(6), values(7), values(8)


contains

!■■■最も手前の壁面を探すサブルーチン■■■
! OpenMPの並列領域から呼ぶので、recursiveにして局所変数をスタックに置く
! （ホストの変数と同名の局所変数はすべてここで宣言し、スレッド間で共有しないようにする）
recursive subroutine nearest_wall(vini, vray, flag, jhit, dhit, thit)

	implicit none

	real*8,		intent(in)	:: vini(3)	! 音線の基点
	real*8,		intent(in)	:: vray(3)	! 音線の単位方向ベクトル
	integer*4,	intent(out)	:: flag		! 「壁面に当たる」フラグ
	integer*4,	intent(out)	:: jhit		! 最も手前の壁面番号
	real*8,		intent(out)	:: dhit		! 最も手前の壁面までの距離
	real*8,		intent(out)	:: thit		! そのときの直線のパラメータt
	integer*4	:: j, l
	real*8		:: temp, temp2, tpar, node(3), vectmp(3,3), nmtmp(2,3), inpro(2)

	flag	= 0
	jhit	= 0
	dhit	= 1.0d50
	thit	= 1.0d50

	! 壁面ループ
	do j = 1, sfcount

		! 音線ベクトルと壁面の法線ベクトルの内積が0より小さいとき（音線が壁面に向かっているとき）
		if (vray(1) * vwnm(j,1) + vray(2) * vwnm(j,2) + vray(3) * vwnm(j,3) .lt. 0.0d0) then

			! 平面の方程式ax+by+cz+d=0のdを算出
			temp2	= -(vwnm(j,1) * vwpt(surf(j,2),1) + vwnm(j,2) * vwpt(surf(j,2),2) + vwnm(j,3) * vwpt(surf(j,2),3))

			! 直線と壁面が交わるときのパラメータtを算出
			tpar	= -(vwnm(j,1) * vini(1) + vwnm(j,2) * vini(2) + vwnm(j,3) * vini(3) + temp2) &
					/ (vwnm(j,1) * vray(1) + vwnm(j,2) * vray(2) + vwnm(j,3) * vray(3))

			! 基点から音線方向を向いて正の方向にある面を検出
			if(tpar .gt. 0.0d0) then

				! 交点の座標を算出
				do l = 1, 3
					node(l)	= vini(l) + tpar * vray(l)
				end do

				! 頂点2を基準とした面を張るベクトルとと交点までのベクトルの外積の内積を算出
				do l = 1, 3
					vectmp(1,l)	= vwpt(surf(j,3),l) - vwpt(surf(j,2),l)
					vectmp(2,l)	= vwpt(surf(j,4),l) - vwpt(surf(j,2),l)
					vectmp(3,l)	= node(l) - vwpt(surf(j,2),l)
				end do
				nmtmp(1,1)	= vectmp(3,2) * vectmp(1,3) - vectmp(3,3) * vectmp(1,2)
				nmtmp(1,2)	= vectmp(3,3) * vectmp(1,1) - vectmp(3,1) * vectmp(1,3)
				nmtmp(1,3)	= vectmp(3,1) * vectmp(1,2) - vectmp(3,2) * vectmp(1,1)
				nmtmp(2,1)	= vectmp(3,2) * vectmp(2,3) - vectmp(3,3) * vectmp(2,2)
				nmtmp(2,2)	= vectmp(3,3) * vectmp(2,1) - vectmp(3,1) * vectmp(2,3)
				nmtmp(2,3)	= vectmp(3,1) * vectmp(2,2) - vectmp(3,2) * vectmp(2,1)
				inpro(1)	= nmtmp(1,1) * nmtmp(2,1) + nmtmp(1,2) * nmtmp(2,2) + nmtmp(1,3) * nmtmp(2,3)

				! 頂点3を基準とした面を張るベクトルとと交点までのベクトルの外積の内積を算出
				do l = 1, 3
					vectmp(1,l)	= vwpt(surf(j,2),l) - vwpt(surf(j,3),l)
					vectmp(2,l)	= vwpt(surf(j,4),l) - vwpt(surf(j,3),l)
					vectmp(3,l)	= node(l) - vwpt(surf(j,3),l)
				end do
				nmtmp(1,1)	= vectmp(3,2) * vectmp(1,3) - vectmp(3,3) * vectmp(1,2)
				nmtmp(1,2)	= vectmp(3,3) * vectmp(1,1) - vectmp(3,1) * vectmp(1,3)
				nmtmp(1,3)	= vectmp(3,1) * vectmp(1,2) - vectmp(3,2) * vectmp(1,1)
				nmtmp(2,1)	= vectmp(3,2) * vectmp(2,3) - vectmp(3,3) * vectmp(2,2)
				nmtmp(2,2)	= vectmp(3,3) * vectmp(2,1) - vectmp(3,1) * vectmp(2,3)
				nmtmp(2,3)	= vectmp(3,1) * vectmp(2,2) - vectmp(3,2) * vectmp(2,1)
				inpro(2)	= nmtmp(1,1) * nmtmp(2,1) + nmtmp(1,2) * nmtmp(2,2) + nmtmp(1,3) * nmtmp(2,3)

				! 交点が面内かどうかを判定
				if((inpro(1) .le. 0.0d0) .and. (inpro(2) .le. 0.0d0)) then

					! 基点からの距離を算出
					temp	= sqrt((vini(1) - node(1)) ** 2.0d0 &
								 + (vini(2) - node(2)) ** 2.0d0 &
								 + (vini(3) - node(3)) ** 2.0d0)

					! 最も手前の交点を持つ面を算出
					if(temp .lt. dhit) then
						dhit	= temp
						thit	= tpar
						jhit	= j
					end if

					! 「壁面に当たる」フラグを立てる
					flag	= 1

				end if

			end if

		end if

	end do

end subroutine nearest_wall

!■■■反射履歴をスレッドの作業領域に追加するサブルーチン■■■
! 作業領域が一杯になったら2倍に広げる
recursive subroutine append_hist(b, h)

	implicit none

	type(histbuf),	intent(inout)	:: b
	integer*4,		intent(in)		:: h(:)
	integer*4,	dimension(:,:),	allocatable	:: tmp
	integer*4	:: ii

	if(b%n .eq. size(b%h,2)) then
		allocate(tmp(size(b%h,1),2*size(b%h,2)))
		do ii = 1, b%n
			tmp(:,ii)	= b%h(:,ii)
		end do
		call move_alloc(tmp, b%h)
	end if
	b%n	= b%n + 1
	b%h(:,b%n)	= h(:)

end subroutine append_hist

!■■■反射履歴を辞書式に比べる関数■■■
! a<bなら-1、a=bなら0、a>bなら1。履歴の末尾は0で埋まっているので、両方0になった時点で一致
integer*4 function histcmp(a, b)

	implicit none

	integer*4,	intent(in)	:: a(:), b(:)
	integer*4	:: kk

	histcmp	= 0
	do kk = 1, size(a)
		if(a(kk) .lt. b(kk)) then
			histcmp	= -1
			return
		elseif(a(kk) .gt. b(kk)) then
			histcmp	= 1
			return
		elseif(a(kk) .eq. 0) then
			return
		end if
	end do

end function histcmp

!■■■反射履歴の経路番号を辞書式に整列するサブルーチン（ボトムアップのマージソート）■■■
subroutine merge_sort(dat, idx, n, verbose)

	implicit none

	integer*4,	intent(in)		:: dat(:,:)
	integer*4,	intent(in)		:: n
	integer*4,	intent(inout)	:: idx(:)
	logical,	intent(in)		:: verbose	! 進捗を表示するかどうか
	integer*4,	dimension(:),	allocatable	:: tmp
	integer*4	:: width, lo, mid, hi, ia, ib, ic
	integer*4	:: npass					! 統合の段階数
	integer*8	:: ndone8, ntot8, nstep8, nnext8	! 進捗表示：処理した件数、全件数（段階数×件数）、表示の間隔、次に表示する件数

	!■進捗表示の準備（各段階で全件を1回ずつ処理する）■
	npass	= 0
	width	= 1
	do while(width .lt. n)
		npass	= npass + 1
		width	= 2 * width
	end do
	ntot8	= int(npass, 8) * int(n, 8)
	nstep8	= max(1_8, ntot8 / 1000_8)
	nnext8	= nstep8
	ndone8	= 0

	allocate(tmp(n))
	width	= 1
	do while(width .lt. n)
		lo	= 1
		do while(lo .le. n)
			mid	= min(lo + width - 1, n)
			hi	= min(lo + 2 * width - 1, n)
			ia	= lo
			ib	= mid + 1
			do ic = lo, hi
				!■進捗表示（1件ごとに数え、nstep8件ごとに表示する）■
				ndone8	= ndone8 + 1
				if(verbose .and. ndone8 .ge. nnext8 .and. ndone8 .lt. ntot8) then
					write(*,'(f10.2,a,i12,a,i12,a,$)') dble(ndone8) / dble(ntot8) * 100.0d0, '% ', ndone8, '/', ntot8, char(13)
					nnext8	= nnext8 + nstep8
				end if
				if(ib .gt. hi) then
					tmp(ic)	= idx(ia)
					ia		= ia + 1
				elseif(ia .gt. mid) then
					tmp(ic)	= idx(ib)
					ib		= ib + 1
				elseif(histcmp(dat(:,idx(ia)), dat(:,idx(ib))) .le. 0) then
					tmp(ic)	= idx(ia)
					ia		= ia + 1
				else
					tmp(ic)	= idx(ib)
					ib		= ib + 1
				end if
			end do
			lo	= lo + 2 * width
		end do
		idx(1:n)	= tmp(1:n)
		width	= 2 * width
	end do
	deallocate(tmp)
	if(verbose .and. npass .gt. 0) then
		write(*,'(f10.2,a,i12,a,i12)') 100.0d0, '% ', ntot8, '/', ntot8
	end if

end subroutine merge_sort

!■■■受音結果を比べる関数（到達時間、虚音源方向xyz、反射回数、エネルギーの順に比べる）■■■
integer*4 function rescmp(a, b)

	implicit none

	real*8,	intent(in)	:: a(:), b(:)
	integer*4	:: kk
	integer*4,	parameter	:: ord(11) = (/2, 3, 4, 5, 1, 6, 7, 8, 9, 10, 11/)

	rescmp	= 0
	do kk = 1, 11
		if(a(ord(kk)) .lt. b(ord(kk))) then
			rescmp	= -1
			return
		elseif(a(ord(kk)) .gt. b(ord(kk))) then
			rescmp	= 1
			return
		end if
	end do

end function rescmp

!■■■受音結果の番号を到達時間の順に整列するサブルーチン（ボトムアップのマージソート）■■■
subroutine sort_res(dat, idx, n)

	implicit none

	real*8,		intent(in)		:: dat(:,:)
	integer*4,	intent(in)		:: n
	integer*4,	intent(inout)	:: idx(:)
	integer*4,	dimension(:),	allocatable	:: tmp
	integer*4	:: width, lo, mid, hi, ia, ib, ic

	allocate(tmp(n))
	width	= 1
	do while(width .lt. n)
		lo	= 1
		do while(lo .le. n)
			mid	= min(lo + width - 1, n)
			hi	= min(lo + 2 * width - 1, n)
			ia	= lo
			ib	= mid + 1
			do ic = lo, hi
				if(ib .gt. hi) then
					tmp(ic)	= idx(ia)
					ia		= ia + 1
				elseif(ia .gt. mid) then
					tmp(ic)	= idx(ib)
					ib		= ib + 1
				elseif(rescmp(dat(:,idx(ia)), dat(:,idx(ib))) .le. 0) then
					tmp(ic)	= idx(ia)
					ia		= ia + 1
				else
					tmp(ic)	= idx(ib)
					ib		= ib + 1
				end if
			end do
			lo	= lo + 2 * width
		end do
		idx(1:n)	= tmp(1:n)
		width	= 2 * width
	end do
	deallocate(tmp)

end subroutine sort_res

end program main
