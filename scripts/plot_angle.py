#!/usr/bin/env python
"""把 `ferro traj angle` 的产物画成发表级的键角分布图。

布局与 `plot_gr.py` 相同：一个 csv 一张子图，子图里一条轨迹一条曲线。

**单轴，不像 gr 那样上双轴**：产物里的 `count` 与 `p` 是同一条曲线的两种刻度
（`p` 是 `count` 归一化的结果），不是两个量，画成双轴会得到两条完全重合的线。
`count` 的用处是跟 dump2analysis 逐 bin 对拍（见 `compare_angle.py`），不是出图。

用法：
    python plot_angle.py angle_O-P-O.csv angle_O-Al-O.csv --outdir figs
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

import ferroplot as fp
from plot_gr import expand

# ── 配置：常改的量都在这里 ───────────────────────────────────────────────────

CFG = {
    "panel_w": 3.3,
    "panel_h": 2.5,
    "ncols": None,          # None = 一行排开
    "wspace": 0.35,
    "hspace": 0.45,

    "lw": 1.0,
    "smooth": 0,            # >0 则对 P(θ) 做该窗口的滑动平均（bin 很细时有用）

    "xlabel": r"$\theta$ (deg)",
    "ylabel": r"$P(\theta)$",
    "title_fmt": "{a}--{b}--{c}",
    "legend_loc": "upper left",
    "legend_ncol": 1,

    "xlim": (0, 180),
    "ylim": (0, None),
}


# ── 绘图 ─────────────────────────────────────────────────────────────────────

def triplet_of(df):
    """子图标题取自数据列。多于一个三元组就拒绝——一张子图说不清是哪个角。"""
    cols = ["end_a", "center", "end_c"]
    tri = df[cols].drop_duplicates()
    if len(tri) != 1:
        raise SystemExit(
            f"这份产物含 {len(tri)} 个三元组，一张子图画不下。"
            f"用 `ferro traj angle -a END_A -b CENTRE -c END_C` 重跑。"
        )
    row = tri.iloc[0]
    return str(row["end_a"]), str(row["center"]), str(row["end_c"])


def build_panel(df, ax, colors):
    a, b, c = triplet_of(df)
    for name in fp.files_in(df):
        sub = df[df["file"] == name].sort_values("angle")
        y = sub["p"]
        if CFG["smooth"] > 1:
            y = y.rolling(CFG["smooth"], center=True, min_periods=1).mean()
        ax.plot(sub["angle"], y, color=colors[name], lw=CFG["lw"], label=name)

    ax.set_xlabel(CFG["xlabel"])
    ax.set_ylabel(CFG["ylabel"])
    ax.set_title(CFG["title_fmt"].format(a=a, b=b, c=c))
    ax.set_xlim(*CFG["xlim"])
    ax.set_ylim(*CFG["ylim"])


def build_figure(frames):
    n = len(frames)
    ncols = CFG["ncols"] or n
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(CFG["panel_w"] * ncols, CFG["panel_h"] * nrows),
        squeeze=False,
    )
    fig.subplots_adjust(wspace=CFG["wspace"], hspace=CFG["hspace"])

    names = list(dict.fromkeys(sum((fp.files_in(df) for _, df in frames), [])))
    colors = fp.color_map(names)

    flat = axes.ravel()
    for ax, (_, df) in zip(flat, frames):
        build_panel(df, ax, colors)
    for ax in flat[n:]:
        ax.set_visible(False)

    handles, labels = flat[0].get_legend_handles_labels()
    flat[0].legend(handles, labels, loc=CFG["legend_loc"], ncol=CFG["legend_ncol"])
    return fig


def main():
    ap = argparse.ArgumentParser(description="键角分布 P(theta) 发表级绘图")
    ap.add_argument("inputs", nargs="+", help="ferro traj angle 的 csv（可用 glob）")
    ap.add_argument("-o", "--output", default="angle", help="产物文件名 stem")
    fp.add_outdir_arg(ap)
    args = ap.parse_args()

    fp.use_style()
    frames = [(p, fp.read(p)) for p in expand(args.inputs)]
    print(f"Inputs: {len(frames)} file(s)")

    fig = build_figure(frames)
    fp.save(fig, args.output, args.outdir)


if __name__ == "__main__":
    main()
