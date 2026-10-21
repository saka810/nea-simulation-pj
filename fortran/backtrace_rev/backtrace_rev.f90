! ifort backtrace_rev.f90 /O3 /Qopenmp
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
!  5. 出力は反射履歴の辞書式順に1スレッドで書く（旧版は並列に書き込み、行の順序が実行ごとに変わった）。
!  6. 音線方向のフィボナッチ格子の式を z = (2i-1)/N - 1 に直した（旧版は (2i+1)/N - 1 で、最後の1本が |z|>1 となりNaNだった）。

program main

!$ use omp_lib

implicit none

real*8,			parameter	:: tmax		= 3.0d0				! 計算する時間長[s]
real*8,			parameter	:: nhitmin	= 3.0d0				! 時刻tmaxの虚音源1つあたりに受音球へ届く音線数の期待値（見落とす確率はexp(-nhitmin)）
real*8,			parameter	:: srcx		=  5.00d0			! 音源x座標[m]
real*8,			parameter	:: srcy		= 10.00d0			! 音源y座標[m]
real*8,			parameter	:: srcz		=  4.20d0			! 音源z座標[m]
real*8,			parameter	:: rcvx		= 20.00d0			! 受音点x座標[m]
real*8,			parameter	:: rcvy		= 10.00d0			! 受音点y座標[m]
real*8,			parameter	:: rcvz		=  1.20d0			! 受音点z座標[m]
real*8,			parameter	:: rcvr		=  1.0d0			! 受音球半径[m]
real*8,			parameter	:: c0		= 340.0d0			! 音速
integer*4,		parameter	:: nref		= 100				! 反射回数の上限
real*8,			parameter	:: tolnm	= 1.0d-10			! 同一平面の判定：1-(単位法線ベクトルの内積)の許容値
real*8,			parameter	:: told		= 1.0d-4			! 同一平面の判定：平面の方程式ax+by+cz+d=0のdの差の許容値[m]
character*128,	parameter	:: fnabs	= 'absorption.csv'	! 吸音率リストファイル名
character*128,	parameter	:: fndxf	= 'hall1.dxf'		! 形状ファイル名
character*128,	parameter	:: fnout	= 'hall1-r1.txt'		! 出力ファイル名
character*1,	parameter	:: ynshade	= 'y'				! shadeで書き出したかどうか（y: shade, n: shade以外）
character*1,	parameter	:: ynnmrev	= 'y'				! 法線が外向きかどうか（y: 外向き（反転する）, n: 内向き（反転しない））
character*1,	parameter	:: ynmeter	= 'n'				! dxfの単位がメートルかどうか（y: メートル, n: ミリメートル）

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
integer*4,		dimension(:,:),	allocatable	:: traceff	! 有効反射経路履歴(反射順,経路番号)
integer*4,		dimension(:,:),	allocatable	:: tracred	! 非重複反射履歴(反射順,経路番号)
integer*4,		dimension(:),	allocatable	:: idx		! 整列用の経路番号

real*8,			dimension(:,:),	allocatable	:: res		! 受音結果(反射回数・受音時間・虚音源方向xyz・6帯域のエネルギー,経路番号)
logical,		dimension(:),	allocatable	:: valid	! 受音結果が有効かどうか(経路番号)

integer*4	:: tractmp(nref)	! 各音線用の一時反射履歴（平面番号）
real*8		:: enertmp(6)		! 各経路用の一時エネルギー
real*8		:: isrc(0:nref,3)	! 虚音源リスト

integer*4	:: i, j, k, l, m	! ループ変数
integer*4	:: nray				! 音線数
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
integer*4	:: ndone, nd		! 進捗表示：追跡を終えた音線数、各スレッドが数えた時点の値
integer*4	:: nshown, nstep	! 進捗表示：表示済みの音線数、表示の間隔
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


!■■■音線ループ■■■
write(*,'(a)') '音線追跡'

!■スレッドごとの受音反射履歴の確保■
nthr	= 1
!$ nthr	= omp_get_max_threads()
allocate(buf(0:nthr-1))
do m = 0, nthr - 1
	allocate(buf(m)%h(nref,1024))
	buf(m)%n	= 0
end do

!■球面上のフィボナッチ格子の離散化幅■
dt	= pi * (3.0d0 - sqrt(5.0d0))
dz	= 2.0d0 / dble(nray)

ntrunc	= 0
nescape	= 0
ndone	= 0
nshown	= 0
nstep	= max(1, nray / 1000)
!$omp parallel do schedule(dynamic,1000) private(i,k,l,tid,tractmp,vray,disttmp,vini,plen,flag,jtmp,dhit,thit,node,vtgt,distd,distr,revtmp,nd) reduction(+: ntrunc, nescape)
do i = 1, nray

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

	!■進捗表示（追跡を終えた音線数を全スレッドで数え、nstep本ごとに表示する。表示済みより小さい値は表示しない）■
	!$omp atomic capture
	ndone	= ndone + 1
	nd		= ndone
	!$omp end atomic
	if(mod(nd, nstep) .eq. 0 .and. nd .lt. nray) then
		!$omp critical (progress)
		if(nd .gt. nshown) then
			nshown	= nd
			write(*,'(f10.2,a,i10,a,i10,a,$)') dble(nd) / dble(nray) * 100.0d0, '% ', nd, '/', nray, char(13)
		end if
		!$omp end critical (progress)
	end if

end do
write(*,'(f10.2,a,i10,a,i10)') 100.0d0, '% ', nray, '/', nray
write(*,*) '...完了'

!■スレッドごとの受音反射履歴を1つにまとめる■
count	= 0
do m = 0, nthr - 1
	count	= count + buf(m)%n
end do
allocate(traceff(nref,count))
count	= 0
do m = 0, nthr - 1
	do i = 1, buf(m)%n
		do k = 1, nref
			traceff(k,count+i)	= buf(m)%h(k,i)
		end do
	end do
	count	= count + buf(m)%n
	deallocate(buf(m)%h)
end do
deallocate(buf)
write(*,'(a,i12)') '  受音した経路数（重複を含む）', count
write(*,'(a,i12)') '  反射回数の上限で打ち切った音線数', ntrunc
write(*,'(a,i12)') '  壁面に当たらず室外へ抜けた音線数', nescape


!■■■重複経路の削除■■■
! 進捗は時間のかかる整列（マージソート）の処理件数で表示する
write(*,'(a)') '重複経路の削除'

!■反射履歴を辞書式に整列■
allocate(idx(count))
do i = 1, count
	idx(i)	= i
end do
call merge_sort(traceff, idx, count)

!■隣どうしを比べ、直前に残した履歴と異なるものだけ残す■
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

!■非重複経路を格納■
allocate(tracred(nref,countred))
do i = 1, countred
	do k = 1, nref
		tracred(k,i)	= traceff(k,idx(i))
	end do
end do
deallocate(traceff)
deallocate(idx)
write(*,*) '...完了'
write(*,'(a,i12)') '  非重複経路数', countred


!■■■バックトレース計算■■■
write(*,'(a)') 'バックトレース計算'
allocate(res(11,countred))
allocate(valid(countred))

!■非重複経路ループ■
ndone	= 0
nshown	= 0
nstep	= max(1, countred / 1000)
!$omp parallel do schedule(dynamic,100) private(i,k,l,ktmp,isrc,enertmp,temp,vray,disttmp,vini,flag,jtmp,dhit,thit,node,distr,vtgt,costmp,nd)
do i = 1, countred

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
		if(tracred(k,i) .eq. 0) exit
		ktmp	= k
		temp	= pnm(tracred(k,i),1) * isrc(k-1,1) + pnm(tracred(k,i),2) * isrc(k-1,2) &
				+ pnm(tracred(k,i),3) * isrc(k-1,3) + pd(tracred(k,i))
		do l = 1, 3
			isrc(k,l)	= isrc(k-1,l) - 2.0d0 * temp * pnm(tracred(k,i),l)
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
			if(pid(jtmp) .ne. tracred(k,i)) exit

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

	!■進捗表示（計算を終えた経路数を全スレッドで数え、nstep本ごとに表示する。表示済みより小さい値は表示しない）■
	!$omp atomic capture
	ndone	= ndone + 1
	nd		= ndone
	!$omp end atomic
	if(mod(nd, nstep) .eq. 0 .and. nd .lt. countred) then
		!$omp critical (progress)
		if(nd .gt. nshown) then
			nshown	= nd
			write(*,'(f10.2,a,i10,a,i10,a,$)') dble(nd) / dble(countred) * 100.0d0, '% ', nd, '/', countred, char(13)
		end if
		!$omp end critical (progress)
	end if

end do
write(*,'(f10.2,a,i10,a,i10)') 100.0d0, '% ', countred, '/', countred
write(*,*) '...完了'

!■受音結果を書き出す■
open(1, file = trim(fnout), status = 'replace')
nvalid	= 0
do i = 1, countred
	if(valid(i)) then
		nvalid	= nvalid + 1
		write(1,'(11e30.20e3)') (res(l,i), l=1,11)
	end if
end do
close(1)
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
deallocate(tracred)
deallocate(res)
deallocate(valid)
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
subroutine merge_sort(dat, idx, n)

	implicit none

	integer*4,	intent(in)		:: dat(:,:)
	integer*4,	intent(in)		:: n
	integer*4,	intent(inout)	:: idx(:)
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
				if(ndone8 .ge. nnext8 .and. ndone8 .lt. ntot8) then
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
	if(npass .gt. 0) then
		write(*,'(f10.2,a,i12,a,i12)') 100.0d0, '% ', ntot8, '/', ntot8
	end if

end subroutine merge_sort

end program main
