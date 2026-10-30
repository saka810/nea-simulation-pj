"""音線の本数を「計算したい時間の長さ」から決める（2026-10-21）。

豊田氏の修正版 backtrace_rev.f90（2026-10-19）の改善 3 に合わせたもの。

距離 c·t の虚音源から、半径 r の受音球に届く音線の本数の期待値は
（Vorländer 1989 の式(8)）

    k(t) = N r² / (4 (c t)²)

1 本も届かなければその虚音源は見つからない。その確率は e^{-k}。
そこで「応答の長さ t の時刻で k 本は届いてほしい」から本数を決める（同 式(9)）：

    N = 4 k (c t)² / r²

修正版の既定は t = 3 s・k = 3（見落とす確率 e^{-3} ≈ 5%）で約 1,250 万本。

★**計算時間を優先して本数を指定することもできる**（2026-10-21 ユーザー指示
「計算時間優先で、任意で本数を決めたい場合もあります。その場合、足りない時は、
精度怪しい、という表示を出すようにしてください」）。そのときは
応答の長さの時刻で期待値が 1 本を切ったら**精度が怪しい**と知らせ、
期待値が 1 本になる時刻（それより後は見落としが多い）を必ず添える。
"""
import math

AUTO = "auto"           # 応答の長さから決める
MANUAL = "manual"       # 本数を指定する
MODES = (AUTO, MANUAL)

DEFAULT_HITS = 3.0      # 応答の長さの時刻で届いてほしい本数（修正版の既定）
# ★これを切ったら「精度が怪しい」。期待値 1 本＝見落とす確率 e^{-1} ≈ 37%
DOUBTFUL_HITS = 1.0

# ★★**反射回数の上限も応答の長さから決める**（2026-10-30。豊田氏の修正版 10-23 版の改善 7）。
#   それまでは `nref` を定数（既定 120）で与えていたので、**小さい室では応答の長さまで
#   届かなかった**（JR 研修室：平均自由行程 3.8 m × 120 回 ≒ 460 m ＝ 1.5 s で 500 万本の
#   ほぼ全部が上限で止まり、T20 が 3〜4 s の条件は後半が欠けていた）。
#   方向を均等に散らした `NREF_SAMPLE` 本を、伝搬距離が「音速 × 応答の長さ」に
#   達するまで反射回数の制限なしで追い（`loop_reflectionmesh.sample_reflections`）、
#   その最大反射回数に `nref_margin` を掛けた値を上限にする
#   - 余裕の既定は 1.25（豊田氏のメールの「AI の推奨値」。豊田氏自身は正確さ重視で 10）
#   - ★昔の project.json（`nref_mode` が無い）と反射回数を渡したプロジェクトは `manual`
#     （`rays_mode` と同じ約束。結果が黙って変わらないように）
NREF_AUTO = "auto"
NREF_MANUAL = "manual"
NREF_MODES = (NREF_AUTO, NREF_MANUAL)
DEFAULT_NREF_MARGIN = 1.25
NREF_SAMPLE = 100_000       # 試し追跡の本数（修正版の nsample）
NREF_CAP = 1_000_000        # 試し追跡で 1 本を追う反射回数の打ち切り（修正版の nrefcap）


def reflection_limit(max_reflections, margin=DEFAULT_NREF_MARGIN):
    """試し追跡の最大反射回数から上限を決める（修正版の `nref = ceiling(nrefmrg * kmax)`）。"""
    return max(1, int(math.ceil(float(margin) * int(max_reflections))))


def rays_for(time, radius, hits=DEFAULT_HITS, sound_velocity=343.0):
    """時刻 `time` [s] で受音球に `hits` 本届く音線の本数。"""
    if time <= 0 or radius <= 0 or hits <= 0:
        raise ValueError("応答の長さ・受音球の半径・届いてほしい本数は正の値にしてください")
    return int(math.ceil(4.0 * hits * (sound_velocity * time) ** 2 / radius ** 2))


def expected_hits(rays, time, radius, sound_velocity=343.0):
    """音線 `rays` 本のとき、時刻 `time` [s] で受音球に届く本数の期待値。"""
    return rays * radius ** 2 / (4.0 * (sound_velocity * time) ** 2)


def time_for_hits(rays, radius, hits=1.0, sound_velocity=343.0):
    """届く本数の期待値が `hits` 本になる時刻 [s]。それより後は見落としが増える。"""
    return math.sqrt(rays * radius ** 2 / (4.0 * hits)) / sound_velocity


def miss_probability(hits):
    """期待値 `hits` 本のとき、その虚音源を見落とす確率 e^{-k}。"""
    return math.exp(-hits)


class Plan:
    """使う音線の本数と、その見込み（期待値・見落とし・注意）。"""

    def __init__(self, mode, rays, time, radius, hits_wanted, sound_velocity):
        self.mode = mode
        self.rays = int(rays)
        self.time = float(time)
        self.radius = float(radius)
        self.hits_wanted = float(hits_wanted)
        self.sound_velocity = float(sound_velocity)
        self.hits = expected_hits(self.rays, self.time, self.radius, self.sound_velocity)
        self.one_hit_time = time_for_hits(self.rays, self.radius, 1.0, self.sound_velocity)
        self.miss = miss_probability(self.hits)

    @property
    def doubtful(self):
        """★精度が怪しいか（応答の長さの時刻で期待値が 1 本を切る）。"""
        return self.hits < DOUBTFUL_HITS

    @property
    def mode_label(self):
        return "自動（応答の長さから）" if self.mode == AUTO else "指定"

    def summary(self):
        """1 行で：本数・期待値・見落とし。"""
        return (f"音線 {self.rays:,} 本（{self.mode_label}）／応答の長さ {self.time:g} s の時刻で"
                f"届く本数の期待値 {self.hits:.2f} 本（見落とす確率 {self.miss * 100:.0f}%）／"
                f"期待値が 1 本になる時刻 {self.one_hit_time:.2f} s")

    def warning(self):
        """精度が怪しいときの知らせ。怪しくなければ None。"""
        if not self.doubtful:
            return None
        need = rays_for(self.time, self.radius, DOUBTFUL_HITS, self.sound_velocity)
        return (f"★精度が怪しい：音線 {self.rays:,} 本では、{self.one_hit_time:.2f} s より後の"
                f"反射音の多くを見落とします（応答の長さ {self.time:g} s の時刻で期待値 "
                f"{self.hits:.2f} 本・見落とす確率 {self.miss * 100:.0f}%）。"
                f"後半の残響（T20 / T30）は当てにしないでください。"
                f"期待値 1 本なら {need:,} 本、{DEFAULT_HITS:g} 本なら "
                f"{rays_for(self.time, self.radius, DEFAULT_HITS, self.sound_velocity):,} 本です")

    def rows(self):
        """計算情報の表に並べる行（項目, 値）。"""
        rows = [("音線数の決め方", self.mode_label),
                ("音線数", self.rays),
                ("応答の長さの時刻で届く本数の期待値", round(self.hits, 3)),
                ("その時刻で見落とす確率 [%]", round(self.miss * 100, 1)),
                ("届く本数の期待値が 1 本になる時刻 [s]", round(self.one_hit_time, 3))]
        if self.mode == AUTO:
            rows.insert(2, ("応答の長さの時刻で届いてほしい本数", self.hits_wanted))
        if self.doubtful:
            rows.append(("精度の注意", self.warning()))
        return rows


def plan(project, sound_velocity=None):
    """プロジェクトの設定から使う本数を決める。→ `Plan`"""
    if sound_velocity is None:
        from atmosphere import Atmosphere
        sound_velocity = Atmosphere(temperature=project.temperature,
                                    humidity=project.humidity,
                                    pressure=project.pressure).sound_velocity
    mode = getattr(project, "rays_mode", MANUAL) or MANUAL
    if mode not in MODES:
        raise ValueError(f"音線数の決め方が分かりません: {mode!r}（{' / '.join(MODES)}）")
    hits_wanted = float(getattr(project, "rays_hits", None) or DEFAULT_HITS)
    time, radius = float(project.max_time), float(project.radius)
    if mode == AUTO:
        rays = rays_for(time, radius, hits_wanted, sound_velocity)
    else:
        rays = int(project.rays)
    return Plan(mode, rays, time, radius, hits_wanted, sound_velocity)
