;;; ============================================================================
;;;  閉じ判定 — 面が閉じているかを AutoCAD の中で確かめる
;;;  日本環境アメニティ株式会社／幾何音響シミュレーション PJ
;;;
;;;  ★2026-08-28 ユーザー要望
;;;    「この閉じた判定ができるシステムだけ、他の人に渡して確認できるようにしたい。
;;;      モデル作成者に確認してもらいながら、作成して欲しいので。
;;;      ただし、元々のデータは dxf で、リージョンのデータを持っています。」
;;;
;;;  使い方（詳しくは 使い方.md）
;;;    CHECKCLOSED … 閉じ判定を実行する。自由端が赤い線で図面に出る
;;;    CHECKNEXT   … 自由端のかたまりを 1 つずつ拡大して見る
;;;    CHECKCLEAR  … 判定で描いた線を消す
;;;
;;;  何をしているか
;;;    面をぐるりと囲む辺を全部集め、**何枚の面に共有されているか**を数える。
;;;    閉じた立体では、どの辺もちょうど 2 枚に共有される。
;;;    1 枚にしか属さない辺が「開いた辺」。
;;;
;;;    ★ただし開いた辺には**性質のまったく違う 2 種類**が混ざっている。
;;;
;;;      自由端   … 他の辺に覆われていない。宙に浮いた片面の板の外周・面の抜け。
;;;                 **閉じた室のつもりならここが作図ミス**
;;;      T字接合 … 同じ直線上の他の開いた辺に覆われている。壁を帯で分割した
;;;                 継ぎ目など。**面としては閉じている**（辺の分け方が違うだけ）
;;;
;;;    この 2 つを混ぜて報告すると「まだ穴だらけだ」と誤解するので必ず分ける。
;;;
;;;  ★REGION / 3DSOLID はその場で分解して辺を拾い、**UNDO で元に戻す**。
;;;    図面は変わらない（描き足す赤い線だけが残る）。
;;;
;;;  ★判定の細かさは**拾った面の大きさ**から決める（図面の EXTMIN/EXTMAX では
;;;    ない）。図面には通芯や文字が遠くに散らばっていることがあり、そこから
;;;    決めると許容が緩くなって**離れた点が同じ点とみなされ、開いた辺を
;;;    見落とす**（2026-08-28 の検算で、27 本あるはずが 7 本になった）。
;;; ============================================================================

(setq *cc-version* "1.3")

;; 判定で描く線の置き場（CHECKCLEAR はこの 2 つを消す）
(setq *cc-layer-free* "_閉じ判定_自由端")
(setq *cc-layer-tee*  "_閉じ判定_T字接合")

;; ---------------------------------------------------------------- 小物 ----

(defun cc-fix (v)
  (fix (+ v (if (minusp v) -0.5 0.5)))
)

(defun cc-key (p res)
  (strcat (itoa (cc-fix (* (car p) res))) ","
          (itoa (cc-fix (* (cadr p) res))) ","
          (itoa (cc-fix (* (caddr p) res))))
)

;; 辺の呼び名。**向きに依らない**（a→b と b→a を同じものとして数える）
(defun cc-edge-key (a b / ka kb)
  (setq ka (cc-key a *cc-res*) kb (cc-key b *cc-res*))
  (if (< ka kb) (strcat ka "|" kb) (strcat kb "|" ka))
)

(defun cc-note-layer (name)
  (if (and name (not (member name *cc-face-layers*)))
    (setq *cc-face-layers* (cons name *cc-face-layers*)))
)

(defun cc-3d (p)
  (list (float (car p)) (float (cadr p))
        (if (caddr p) (float (caddr p)) 0.0))
)

;; ★ここでは**点をそのまま溜めるだけ**。判定の細かさは全部集めてから決める
(defun cc-add-edge (a b)
  (setq a (cc-3d a) b (cc-3d b))
  (setq *cc-raw* (cons (list a b) *cc-raw*))
  (setq *cc-lo* (mapcar 'min *cc-lo* a b))
  (setq *cc-hi* (mapcar 'max *cc-hi* a b))
)

;; 点の並び（輪）から辺を作る
(defun cc-ring (pts / n i)
  (setq n (length pts) i 0)
  (while (< i n)
    (cc-add-edge (nth i pts) (nth (rem (1+ i) n) pts))
    (setq i (1+ i))
  )
  ;; ★面の数は**輪を 1 つ拾ったら 1 枚**と数える（数え方を 1 か所にまとめる）
  (setq *cc-faces* (1+ *cc-faces*))
)

;; --------------------------------------------------- 面から辺を拾う ----

;; 分解しなくても辺が分かるもの（3DFACE / 閉じたポリライン）
(defun cc-plain (e / d ty pts p flag z)
  (setq d (entget e) ty (cdr (assoc 0 d)))
  (cc-note-layer (cdr (assoc 8 d)))
  (cond
    ((= ty "3DFACE")
     (setq pts (list (cdr (assoc 10 d)) (cdr (assoc 11 d))
                     (cdr (assoc 12 d)) (cdr (assoc 13 d))))
     ;; 三角形は 4 点目が 3 点目と同じ
     (if (equal (nth 2 pts) (nth 3 pts) 1e-9)
       (setq pts (list (nth 0 pts) (nth 1 pts) (nth 2 pts))))
     (cc-ring pts)
     T)
    ((= ty "LWPOLYLINE")
     (setq flag (cdr (assoc 70 d)))
     (if (= 1 (logand 1 flag))
       (progn
         (setq pts nil z (cdr (assoc 38 d)))
         (if (null z) (setq z 0.0))
         (foreach x d
           (if (= 10 (car x))
             (setq pts (cons (list (cadr x) (caddr x) z) pts))))
         (if (> (length pts) 2) (cc-ring (reverse pts)))
         T)
       nil))
    ((= ty "POLYLINE")
     (setq flag (cdr (assoc 70 d)))
     (cond
       ;; ポリフェイスメッシュ … 頂点が並んだあとに「面レコード」が続く
       ((/= 0 (logand 64 flag)) (cc-polyface e))
       ;; 閉じたポリライン … 輪郭をそのまま面とみなす
       ((= 1 (logand 1 flag))
        (setq pts nil p (entnext e))
        (while (and p (= "VERTEX" (cdr (assoc 0 (entget p)))))
          (setq pts (cons (cdr (assoc 10 (entget p))) pts))
          (setq p (entnext p)))
        (if (> (length pts) 2) (cc-ring (reverse pts)))
        T)
       (T nil)))
    (T nil)
  )
)

;; ポリフェイスメッシュ（`POLYLINE` の 70 に 64）
;;   頂点 VERTEX（70 に 192）が並び、そのあとに面レコード VERTEX（70 に 128）が続く。
;;   面レコードの 71〜74 が頂点番号（1 始まり。負値は「辺を表示しない」という意味だけ）
(defun cc-polyface (e / p d vflag coords faces rec idx pts k n)
  (setq coords nil faces nil p (entnext e))
  (while (and p (= "VERTEX" (cdr (assoc 0 (entget p)))))
    (setq d (entget p) vflag (cdr (assoc 70 d)))
    (if (and (/= 0 (logand 128 vflag)) (= 0 (logand 64 vflag)))
      (setq faces (cons (list (cdr (assoc 71 d)) (cdr (assoc 72 d))
                              (cdr (assoc 73 d)) (cdr (assoc 74 d))) faces))
      (setq coords (cons (cdr (assoc 10 d)) coords)))
    (setq p (entnext p))
  )
  (setq coords (reverse coords) n (length coords))
  (foreach rec (reverse faces)
    (setq pts nil)
    (foreach k rec
      (if (and k (/= 0 k) (<= (abs k) n))
        (setq pts (cons (nth (1- (abs k)) coords) pts))))
    (if (> (length pts) 2)
      (cc-ring (reverse pts)))
  )
  T
)

;; ACIS（REGION / 3DSOLID）は分解して拾う。
;; ★`EXPLODE` は選択を**空行で閉じる**こと。閉じないとコマンドが開いたままになり、
;;   次の行を選択として食べて止まる（2026-08-28 に実際に 2 回止めた）
(defun cc-burst (e / marker c out)
  (setq marker (entlast))
  (command "_.EXPLODE" e "")
  (setq c marker out nil)
  (while (setq c (entnext c)) (setq out (cons c out)))
  (reverse out)
)

(defun cc-acis (queue / e d ty kids guard c)
  (setq guard 0)
  (while (and queue (< guard 20000))
    (setq guard (1+ guard))
    (setq e (car queue) queue (cdr queue))
    ;; ★種別は**分解する前に**控える（分解すると元の図形は消える）
    ;; ★種別は**分解する前に**控える（分解すると元の図形は消える）。
    ;;   REGION は 1 個が 1 枚の面。3DSOLID は面の REGION に割れるので数えない
    (setq ty (cdr (assoc 0 (entget e))))
    (cc-note-layer (cdr (assoc 8 (entget e))))
    (if (= ty "REGION") (setq *cc-faces* (1+ *cc-faces*)))
    (setq kids (cc-burst e))
    (foreach c kids
      (setq d (entget c) ty (cdr (assoc 0 d)))
      (cond
        ((= ty "LINE")
         (cc-add-edge (cdr (assoc 10 d)) (cdr (assoc 11 d))))
        ((member ty '("REGION" "3DSOLID"))
         (setq queue (cons c queue)))
        ((member ty '("ARC" "CIRCLE" "ELLIPSE" "SPLINE"))
         ;; ★曲線の辺は端点だけ拾う。数え方は合うが、覆い判定は当てにならない
         (setq *cc-curved* (1+ *cc-curved*))
         (cc-add-edge (vlax-curve-getStartPoint c)
                      (vlax-curve-getEndPoint c)))
        (T nil)
      )
    )
  )
)

;; ------------------------------------------------- 開いた辺を見つける ----

;; 1 枚にしか属さない辺（＝開いた辺）を返す
(defun cc-open-edges (/ sorted out run key item)
  (setq sorted (vl-sort *cc-edges* '(lambda (a b) (< (car a) (car b)))))
  (setq out nil run nil key nil)
  (foreach item sorted
    (if (and key (= key (car item)))
      (setq run (cons item run))
      (progn
        (if (and run (= 1 (length run))) (setq out (cons (car run) out)))
        (setq key (car item) run (list item))
      )
    )
  )
  (if (and run (= 1 (length run))) (setq out (cons (car run) out)))
  out
)

;; ---------------------------- 自由端か T字接合か（同じ直線上で覆われているか）

(defun cc-unit (a b / d L)
  (setq L (distance a b))
  (if (> L *cc-gap*)
    (mapcar '(lambda (v) (/ v L)) (mapcar '- b a))
    nil)
)

(defun cc-dot (a b) (apply '+ (mapcar '* a b)))

(defun cc-vkey (u)
  (strcat (itoa (cc-fix (* (car u) 1000000.0))) ","
          (itoa (cc-fix (* (cadr u) 1000000.0))) ","
          (itoa (cc-fix (* (caddr u) 1000000.0))))
)

;; 直線の呼び名（向きを一意にしてから、原点からの足で表す）
(defun cc-line-key (a b / u foot k1 k2)
  (setq u (cc-unit a b))
  (if (null u)
    nil
    (progn
      ;; 向きを一意にする（u と -u を同じ直線とみなす）
      (setq k1 (cc-vkey u) k2 (cc-vkey (mapcar '- '(0.0 0.0 0.0) u)))
      (if (> k1 k2) (setq u (mapcar '- '(0.0 0.0 0.0) u)))
      (setq foot (mapcar '- a (mapcar '(lambda (v) (* v (cc-dot a u))) u)))
      ;; ★足の丸めは点の丸めより**粗く**する（同じ直線に乗る辺を取り逃さない）
      (list (strcat (cc-vkey u) "/" (cc-key foot *cc-fres*)) u)
    )
  )
)

;; その区間を覆っている辺の本数（自分も 1 本に数える）
(defun cc-depth (lo hi spans / c s)
  (setq c 0)
  (foreach s spans
    (if (and (>= lo (- (car s) *cc-gap*)) (<= hi (+ (cadr s) *cc-gap*)))
      (setq c (1+ c))))
  c
)

;; 区間ごとに見て、1 つでも「自分しか覆っていない」ところがあれば自由端
(defun cc-uncovered (lo hi bounds spans / i n a b bad)
  (setq bad nil i 0 n (length bounds))
  (while (and (< i (1- n)) (not bad))
    (setq a (nth i bounds) b (nth (1+ i) bounds))
    (if (and (>= a (- lo *cc-gap*)) (<= b (+ hi *cc-gap*))
             (> (- b a) *cc-gap*)
             (< (cc-depth a b spans) 2))
      (setq bad T))
    (setq i (1+ i))
  )
  bad
)

;; 開いた辺を「自由端」と「T字接合」に分ける
(defun cc-split (open / groups lk item found free tee spans bounds
                        u s e bnds members k)
  (setq groups nil)
  (foreach item open
    (setq lk (cc-line-key (nth 1 item) (nth 2 item)))
    (if lk
      (progn
        (setq found (assoc (car lk) groups))
        (if found
          (setq groups (subst (list (car lk) (cadr lk)
                                    (cons item (caddr found)))
                              found groups))
          (setq groups (cons (list (car lk) (cadr lk) (list item)) groups))
        )
      )
      (setq free (cons item free))          ; 長さ 0（起きないはず）
    )
  )

  (setq free nil tee nil)
  (foreach found groups
    (setq u (nth 1 found) members (nth 2 found))
    ;; 直線上の位置（スカラー）に直す
    (setq spans nil bnds nil)
    (foreach item members
      (setq s (cc-dot (nth 1 item) u) e (cc-dot (nth 2 item) u))
      (setq spans (cons (list (min s e) (max s e)) spans))
      (setq bnds (cons (min s e) (cons (max s e) bnds)))
    )
    (setq spans (reverse spans))
    (setq bounds (vl-sort bnds '(lambda (a b) (< a b))))
    ;; 辺ごとに、区間の覆われ具合を見る
    (setq k 0)
    (foreach item members
      (setq s (nth k spans) k (1+ k))
      (if (cc-uncovered (car s) (cadr s) bounds spans)
        (setq free (cons item free))
        (setq tee (cons item tee)))
    )
  )
  (list free tee)
)

;; ------------------------------------------ かたまり（連なり）にまとめる ----
;; ★端点を共有する辺どうしを**ひとつながり**としてまとめる（連結成分）。
;;   1 本ずつ並べても場所が分からないので、かたまりで数える
;;
;; ★★取り込んだ辺だけを `rest` から外すこと。端点の一覧で一括して削ると、
;;   **かたまりに入れていない辺まで消える**（実際に長さが 3 割減った）

(defun cc-drop (items lst / out keep x y)
  (setq out nil)
  (foreach x lst
    (setq keep T)
    (foreach y items (if (eq x y) (setq keep nil)))
    (if keep (setq out (cons x out)))
  )
  (reverse out)
)

(defun cc-chains (edges / rest group keys moved out item x ka kb hits)
  (setq rest edges out nil)
  (while rest
    (setq item (car rest) rest (cdr rest))
    (setq group (list item))
    (setq keys (list (cc-key (nth 1 item) *cc-res*)
                     (cc-key (nth 2 item) *cc-res*)))
    (setq moved T)
    (while moved
      (setq moved nil hits nil)
      (foreach x rest
        (setq ka (cc-key (nth 1 x) *cc-res*) kb (cc-key (nth 2 x) *cc-res*))
        (if (or (member ka keys) (member kb keys))
          (progn
            (setq hits (cons x hits))
            (if (not (member ka keys)) (setq keys (cons ka keys)))
            (if (not (member kb keys)) (setq keys (cons kb keys)))
            (setq moved T)
          )
        )
      )
      (if hits
        (progn
          (foreach x hits (setq group (cons x group)))
          (setq rest (cc-drop hits rest))
        )
      )
    )
    (setq out (cons group out))
  )
  (reverse out)
)

(defun cc-length (chain / L item)
  (setq L 0.0)
  (foreach item chain (setq L (+ L (distance (nth 1 item) (nth 2 item)))))
  L
)

;; ------------------------------------------------------------ 描く ----

(defun cc-draw (edges layer colour / item)
  (if edges
    (progn
      (if (not (tblsearch "LAYER" layer))
        (command "_.-LAYER" "_Make" layer "_Color" colour layer ""))
      (foreach item edges
        (entmake (list '(0 . "LINE") (cons 8 layer)
                       (cons 10 (nth 1 item)) (cons 11 (nth 2 item)))))
    )
  )
)

(defun c:CHECKCLEAR ( / ss name)
  (setvar "CMDECHO" 0)
  (foreach name (list *cc-layer-free* *cc-layer-tee*)
    (if (tblsearch "LAYER" name)
      (progn
        (if (= name (getvar "CLAYER")) (setvar "CLAYER" "0"))
        (command "_.-LAYER" "_Thaw" name "_On" name "_Unlock" name "")
        (setq ss (ssget "_X" (list (cons 8 name))))
        (if ss (command "_.ERASE" ss ""))
      )
    )
  )
  (setvar "CMDECHO" 1)
  (princ "\n[閉じ判定] 判定で描いた線を消しました。")
  (princ)
)

;; ------------------------------------------------ 画面を見やすくする ----
;; ★何が邪魔かを決め打ちせず、**面を拾えた画層かどうか**で決める。
;;   図面ごとに通芯・文字の画層名は違うので、名前で当てにいかない

(defun cc-join (names / out)
  (setq out "")
  (foreach n names
    (setq out (if (= out "") n (strcat out "," n))))
  out
)

(defun cc-all-layers ( / tb out)
  (setq out nil tb (tblnext "LAYER" T))
  (while tb
    (setq out (cons (cdr (assoc 2 tb)) out))
    (setq tb (tblnext "LAYER")))
  (reverse out)
)

(defun cc-saved-p (name / found item)
  (setq found nil)
  (foreach item *cc-view-saved*
    (if (= name (cdr (assoc 2 item))) (setq found T)))
  found
)

(defun c:CHECKVIEW ( / name others judged)
  (if (null *cc-face-layers*)
    (progn (princ "\n[閉じ判定] 先に CHECKCLOSED を実行してください。") (princ))
    (progn
      (setvar "CMDECHO" 0)
      ;; ---- 元の状態を控える（entget を丸ごと。オン/オフ・色・透過が戻る）----
      ;; ★★控えるのは**画層ごとに 1 回目だけ**。2 回目も控え直すと
      ;;   「整えたあとの状態」を元の状態として覚えてしまい、
      ;;   CHECKVIEWOFF で戻らなくなる。
      ;;   ★作りながら何度も回す使い方なので、**あとから増えた画層も
      ;;     そのとき控える**（増えた画層だけ戻らない、を防ぐ）
      (if (null *cc-view-saved*) (setq *cc-view-lw* (getvar "LWDISPLAY")))
      (foreach name (cc-all-layers)
        (if (not (cc-saved-p name))
          (setq *cc-view-saved*
                (cons (entget (tblobjname "LAYER" name)) *cc-view-saved*))))

      ;; 判定の画層を現在層にしておく（消す画層が現在層だと確認を求められる）
      (if (tblsearch "LAYER" *cc-layer-free*)
        (setvar "CLAYER" *cc-layer-free*)
        (setvar "CLAYER" *cc-layer-tee*))

      (setq judged (list *cc-layer-free* *cc-layer-tee*))
      (setq others nil)
      (foreach name (cc-all-layers)
        (if (and (not (member name *cc-face-layers*))
                 (not (member name judged)))
          (setq others (cons name others))))

      ;; ---- 面の無い画層を消す（通芯・文字・ハッチなど）----
      (if others (command "_.-LAYER" "_Off" (cc-join others) ""))
      ;; ---- 面の画層を透かす（中の赤い線が見えるように）----
      (command "_.-LAYER" "_TRansparency" "70" (cc-join *cc-face-layers*) "")
      ;; ---- 判定の線は透かさず、自由端を太く ----
      (command "_.-LAYER" "_TRansparency" "0" (cc-join judged) "")
      (if (tblsearch "LAYER" *cc-layer-free*)
        (command "_.-LAYER" "_LWeight" "0.50" *cc-layer-free* ""))
      (setvar "LWDISPLAY" 1)

      ;; ---- 表示スタイルは X線（面が透ける）----
      ;; ★カメラ（視線の向き）は触らない。合わせた向きを崩さないため
      (command "_.VSCURRENT" "_X")
      (setvar "VSFACEOPACITY" 30)
      (setvar "CMDECHO" 1)
      (princ (strcat "\n[閉じ判定] 画面を整えました（面の無い画層 "
                     (itoa (length others)) " 枚を非表示、面を透過 70%、"
                     "自由端を太線、表示スタイルを X線）。"))
      (princ "\n           CHECKVIEWOFF で元に戻せます。")
      (princ)
    )
  )
)

(defun c:CHECKVIEWOFF ( / item)
  (if (null *cc-view-saved*)
    (princ "\n[閉じ判定] 戻す状態を控えていません（CHECKVIEW を実行していません）。")
    (progn
      (setvar "CMDECHO" 0)
      (setvar "CLAYER" "0")
      (foreach item *cc-view-saved* (entmod item))
      (setvar "LWDISPLAY" *cc-view-lw*)
      (command "_.VSCURRENT" "_2")
      (setq *cc-view-saved* nil)      ; 次に整えるときは今の状態を控え直す
      (setvar "CMDECHO" 1)
      (princ "\n[閉じ判定] 画層の表示・透過・線の太さを元に戻しました。")
      (princ "\n           ★表示スタイルだけは元の設定を読み取れないので、")
      (princ "\n             2D ワイヤフレームにしてあります。必要なら選び直してください。")
    )
  )
  (princ)
)

;; -------------------------------------------------------- 本体 ----

(defun cc-report (free tee groups / total)
  (setq total 0.0)
  (foreach g groups (setq total (+ total (cc-length g))))
  (princ (strcat "\n\n=== 閉じ判定 " *cc-version* " ==="))
  (princ (strcat "\n  面 " (itoa *cc-faces*)
                 " 枚 → 辺 " (itoa (length *cc-edges*)) " 本"))
  (if (= 0 (length free))
    (princ "\n  ★自由端はありません。面は閉じています。")
    (progn
      (princ (strcat "\n  ★自由端 " (itoa (length free)) " 本／"
                     (itoa (length groups)) " かたまり／計 "
                     (rtos total 2 2)))
      (princ "\n      ＝ 宙に浮いた片面の板の外周、または面の抜け。")
      (princ "\n        閉じた室のつもりなら、ここが作図ミスです。")
    )
  )
  (princ (strcat "\n    T字接合 " (itoa (length tee)) " 本"))
  (princ "\n      ＝ 同じ直線上の他の辺に覆われている。**面としては閉じています**")
  (princ "\n        （壁を帯で分けた継ぎ目など。直す必要はありません）")
  (if (> *cc-curved* 0)
    (princ (strcat "\n  ※ 円弧・スプラインの辺が " (itoa *cc-curved*)
                   " 本あります。端点だけで見ているので、その周りは目安です。")))
  (princ (strcat "\n\n  赤い線＝自由端（レイヤ " *cc-layer-free* "）"))
  (princ (strcat "\n  橙の線＝T字接合（レイヤ " *cc-layer-tee* "）"))
  (if (> (length groups) 0)
    (princ "\n  CHECKNEXT で自由端のかたまりを 1 つずつ拡大できます。"))
  (princ "\n  CHECKCLEAR で判定の線を消せます。")
  (princ "\n  CHECKVIEWOFF で画面の見え方を元に戻せます。\n")
)

(defun c:CHECKCLOSED ( / ss n e acis plain diag open split free tee before)
  (vl-load-com)
  (setvar "CMDECHO" 0)
  (c:CHECKCLEAR)

  (setq *cc-raw* nil *cc-curved* 0 *cc-faces* 0 *cc-face-layers* nil)
  (setq *cc-lo* (list 1e20 1e20 1e20) *cc-hi* (list -1e20 -1e20 -1e20))

  ;; ---- ① 分解せずに読めるもの（3DFACE / 閉じたポリライン）----
  (setq plain 0)
  (setq ss (ssget "_X" '((0 . "3DFACE,LWPOLYLINE,POLYLINE") (410 . "Model"))))
  (if ss
    (progn
      (setq n 0)
      (repeat (sslength ss)
        (if (cc-plain (ssname ss n)) (setq plain (1+ plain)))
        (setq n (1+ n))
      )
    )
  )

  ;; ---- ② ACIS（REGION / 3DSOLID）は分解して読み、**UNDO で元に戻す** ----
  (setq acis nil)
  (setq ss (ssget "_X" '((0 . "REGION,3DSOLID") (410 . "Model"))))
  (if ss
    (progn
      (setq n 0)
      (repeat (sslength ss) (setq acis (cons (ssname ss n) acis) n (1+ n)))
      (setq before (sslength ss))
      (princ (strcat "\n[閉じ判定] ACIS " (itoa before)
                     " 個を分解して辺を拾います（すぐ元に戻します）…"))
      ;; ★★戻すのは **`UNDO マーク` → `UNDO 後退`**。
      ;;   `UNDO 開始/終了` ＋ `U` だと「グループの開始点が見つかりました」で
      ;;   止まり、**分解したまま戻らない**（2026-08-28 に実際に踏んだ。
      ;;   1 回目から戻っておらず、2 回目の判定が「面が 1 枚も見つかりません」に
      ;;   なって気づいた）
      (command "_.UNDO" "_Mark")
      (cc-acis (reverse acis))
      (command "_.UNDO" "_Back")
      ;; ★戻ったことを必ず確かめる（黙って図面を壊さない）
      (setq ss (ssget "_X" '((0 . "REGION,3DSOLID") (410 . "Model"))))
      (if (or (null ss) (/= before (sslength ss)))
        (princ (strcat "\n[閉じ判定] ★注意: 元に戻せていないかもしれません"
                       "（分解前 " (itoa before) " 個 → いま "
                       (if ss (itoa (sslength ss)) "0")
                       " 個）。保存せずに閉じ、開き直してください。"))
      )
    )
  )

  (if (null *cc-raw*)
    (progn
      (princ "\n[閉じ判定] 面が 1 枚も見つかりません。")
      (princ "\n           読めるのは REGION / 3DSOLID / 3DFACE / 閉じたポリラインです。")
      (setvar "CMDECHO" 1)
    )
    (progn
      ;; ---- ③ 判定の細かさを**拾った面の大きさ**から決める ----
      (setq diag (distance *cc-lo* *cc-hi*))
      (if (or (null diag) (<= diag 0.0)) (setq diag 1000.0))
      (setq *cc-gap* (/ diag 10000000.0))       ; 同じ点とみなす距離
      (setq *cc-res* (/ 1.0 *cc-gap*))
      (setq *cc-fres* (/ 100000.0 diag))        ; 直線の足の丸め（少し粗く）

      ;; ---- ④ 辺に呼び名を付けて、1 枚にしか属さない辺を探す ----
      (setq *cc-edges* nil)
      (foreach item *cc-raw*
        (if (> (distance (car item) (cadr item)) *cc-gap*)
          (setq *cc-edges* (cons (list (cc-edge-key (car item) (cadr item))
                                       (car item) (cadr item)) *cc-edges*))))
      (princ (strcat "\n[閉じ判定] 辺 " (itoa (length *cc-edges*))
                     " 本を調べています…"))
      (setq open (cc-open-edges))
      (setq split (cc-split open))
      (setq free (car split) tee (cadr split))

      ;; ---- ⑤ 描いて知らせる ----
      (cc-draw free *cc-layer-free* 1)      ; 赤
      (cc-draw tee  *cc-layer-tee*  40)     ; 橙
      (setq *cc-groups*
            (vl-sort (cc-chains free)
                     '(lambda (a b) (> (cc-length a) (cc-length b)))))
      (setq *cc-index* 0)
      (setvar "CMDECHO" 1)
      ;; ★判定のあとに**自動で見やすくする**（2026-08-28 ユーザー要望）
      (c:CHECKVIEW)
      (cc-report free tee *cc-groups*)
    )
  )
  (princ)
)

;; ---- 自由端のかたまりを 1 つずつ見る ----

(defun c:CHECKNEXT ( / g lo hi pad item)
  (if (null *cc-groups*)
    (princ "\n[閉じ判定] 先に CHECKCLOSED を実行してください。")
    (if (>= *cc-index* (length *cc-groups*))
      (progn
        (setq *cc-index* 0)
        (princ "\n[閉じ判定] 最後まで見ました。もう一度 CHECKNEXT で先頭に戻ります。"))
      (progn
        (setq g (nth *cc-index* *cc-groups*))
        (setq lo (nth 1 (car g)) hi (nth 1 (car g)))
        (foreach item g
          (setq lo (mapcar 'min lo (nth 1 item) (nth 2 item)))
          (setq hi (mapcar 'max hi (nth 1 item) (nth 2 item)))
        )
        (setq pad (* 0.6 (max 1.0 (distance lo hi))))
        (setvar "CMDECHO" 0)
        (command "_.ZOOM" "_Window"
                 (mapcar '(lambda (v) (- v pad)) lo)
                 (mapcar '(lambda (v) (+ v pad)) hi))
        (setvar "CMDECHO" 1)
        (princ (strcat "\n[閉じ判定] 自由端 " (itoa (1+ *cc-index*)) "/"
                       (itoa (length *cc-groups*)) "　長さ "
                       (rtos (cc-length g) 2 2) "　辺 " (itoa (length g)) " 本"))
        (setq *cc-index* (1+ *cc-index*))
      )
    )
  )
  (princ)
)

(princ (strcat "\n閉じ判定 " *cc-version*
               " を読み込みました。CHECKCLOSED で実行します。\n"))
(princ)
