"""計算中の進み具合を出すウィンドウ（tkinter）。

条件入力ウィンドウを閉じてから結果が出るまで**何も出ない時間**があり、
動いているのか止まっているのか分からなかった（ユーザー指摘）ので用意した。

  ・いまどの段階か（モデル読み込み／音線追跡／バックトレース …）
  ・重い段階は**その中の進み具合**もバーで出す
  ・コンソールに出ていたログをそのまま画面にも流す

**計算は別スレッドで走らせる。** tkinter のイベントループは主スレッドから
動かす決まりなので、計算を主スレッドでやると画面が固まって（応答なしになって）
進捗どころではなくなる。

スレッド間の受け渡しは `queue.Queue` 1 本だけにしてある。
計算側は queue に積むだけ、画面側は `after()` で定期的に取り出して描く。
ウィジェットを計算スレッドから触らないので競合しない。
"""

import queue
import sys
import threading
import tkinter as tk
import traceback
from tkinter import ttk


class _Tee:
    """標準出力を、元の出力先と queue の両方へ流す。

    計算の進み具合は既存の `print` に詳しく出ているので、
    それを捨てずに画面へ回す。元の出力先にも流すのは、
    ウィンドウを閉じたあともコンソールに残すため。
    """

    def __init__(self, original, sink):
        self.original = original
        self.sink = sink

    def write(self, text):
        if self.original is not None:
            self.original.write(text)
        if text.strip():
            self.sink(("log", text.rstrip("\n")))
        return len(text)

    def flush(self):
        if self.original is not None:
            self.original.flush()


class ProgressWindow:
    """計算の進み具合を出しながら、渡された処理を別スレッドで走らせる。

        window = ProgressWindow("研修室 を計算中")
        result = window.run(lambda progress: run_project.run(project, progress=progress))

    `run()` は処理の戻り値を返す。処理が例外を投げたらそれを送出し直す
    （画面には内容を出したうえで）。
    """

    POLL_MS = 80        # queue を覗く間隔

    def __init__(self, title="計算中", subtitle="", width=760, height=460):
        self.title = title
        self.subtitle = subtitle
        self.size = (width, height)
        self.queue = queue.Queue()
        self.result = None
        self.error = None
        self.finished = False
        self.cancelled = False

    # ---- 計算側から呼ばれる ------------------------------------------

    def _progress(self, stage, fraction=None):
        self.queue.put(("stage", (stage, fraction)))

    # ---- 画面 ----------------------------------------------------------

    def _build(self):
        self.root = tk.Tk()
        self.root.title(self.title)
        self.root.geometry(f"{self.size[0]}x{self.size[1]}")
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)

        outer = ttk.Frame(self.root, padding=14)
        outer.pack(fill="both", expand=True)

        ttk.Label(outer, text=self.title, font=("", 12, "bold")).pack(anchor="w")
        if self.subtitle:
            ttk.Label(outer, text=self.subtitle, foreground="#666").pack(anchor="w")

        self.stage_label = ttk.Label(outer, text="準備中…")
        self.stage_label.pack(anchor="w", pady=(12, 4))

        self.bar = ttk.Progressbar(outer, mode="indeterminate", length=100)
        self.bar.pack(fill="x")
        self.bar.start(12)

        ttk.Label(outer, text="ログ", foreground="#666").pack(anchor="w", pady=(12, 2))
        frame = ttk.Frame(outer)
        frame.pack(fill="both", expand=True)
        self.log = tk.Text(frame, height=14, wrap="none", font=("Consolas", 9))
        scroll = ttk.Scrollbar(frame, command=self.log.yview)
        self.log.configure(yscrollcommand=scroll.set, state="disabled")
        self.log.pack(side="left", fill="both", expand=True)
        scroll.pack(side="right", fill="y")

        self.button = ttk.Button(outer, text="閉じる", command=self._on_close,
                                 state="disabled")
        self.button.pack(anchor="e", pady=(10, 0))

    def _append(self, text):
        self.log.configure(state="normal")
        self.log.insert("end", text + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    def _set_stage(self, stage, fraction):
        if fraction is None:
            # 進み具合が分からない段階は流れるバーにする
            if self.bar["mode"] != "indeterminate":
                self.bar.configure(mode="indeterminate")
                self.bar.start(12)
            self.stage_label.configure(text=stage)
        else:
            if self.bar["mode"] != "determinate":
                self.bar.stop()
                self.bar.configure(mode="determinate", maximum=100.0)
            self.bar["value"] = max(0.0, min(1.0, fraction)) * 100.0
            self.stage_label.configure(text=f"{stage}   {fraction * 100:.0f} %")

    def _pump(self):
        """queue にたまったものを画面へ反映する（主スレッドから定期的に呼ばれる）。"""
        try:
            while True:
                kind, payload = self.queue.get_nowait()
                if kind == "log":
                    self._append(payload)
                elif kind == "stage":
                    self._set_stage(*payload)
                elif kind == "done":
                    self._on_done()
        except queue.Empty:
            pass
        if not self.finished:
            self.root.after(self.POLL_MS, self._pump)

    def _on_done(self):
        self.finished = True
        self.bar.stop()
        self.bar.configure(mode="determinate", maximum=100.0)
        self.bar["value"] = 100.0
        if self.error is None:
            self.stage_label.configure(text="完了しました")
            self.button.configure(state="normal")
            # うまくいったときは待たせずに閉じる。失敗したときだけ読ませる
            self.root.after(400, self.root.destroy)
        else:
            self.stage_label.configure(text="エラーで止まりました")
            self.button.configure(state="normal")

    def _on_close(self):
        if self.finished:
            self.root.destroy()
        else:
            # 計算スレッドは止められないので、閉じるボタンは効かないことを伝える
            self._append("[進捗] 計算中は閉じられません（終わるまでお待ちください）")

    # ---- 入口 ----------------------------------------------------------

    def run(self, work):
        """`work(progress)` を別スレッドで走らせ、終わるまで画面を出す。"""
        self._build()

        def worker():
            stdout = sys.stdout
            sys.stdout = _Tee(stdout, self.queue.put)
            try:
                self.result = work(self._progress)
            except BaseException as e:      # 何が起きても画面に出して伝える
                self.error = e
                self.queue.put(("log", "".join(traceback.format_exc())))
            finally:
                sys.stdout = stdout
                self.queue.put(("done", None))

        thread = threading.Thread(target=worker, daemon=True)
        thread.start()
        self.root.after(self.POLL_MS, self._pump)
        self.root.mainloop()
        thread.join(timeout=1.0)

        if self.error is not None:
            raise self.error
        return self.result


def run_with_progress(title, work, subtitle=""):
    """進捗ウィンドウを出しながら `work(progress)` を実行する。"""
    return ProgressWindow(title, subtitle=subtitle).run(work)
