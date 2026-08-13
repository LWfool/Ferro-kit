#!/usr/bin/env python
"""把 `ferro traj gr` 的产物画成发表级的 g(r) / CN(r) 图。

一个 csv 一张子图，子图里一条轨迹一条曲线：`gr_P-O.csv gr_Al-O.csv` → 两张子图，
每张里若干条成分曲线叠着比。左轴 g(r) 实线，右轴 CN(r) 虚线，同一条轨迹同色。

用法：
    python plot_gr.py gr_P-O.csv gr_Al-O.csv --outdir figs
    python plot_gr.py 'gr_*.csv' --outdir figs -o rdf
"""

import argparse
import glob
from pathlib import Path

import matplotlib.pyplot as plt

import ferroplot as fp

# ── 配置：常改的量都在这里 ───────────────────────────────────────────────────

CFG = {
    # 一张子图的尺寸（英寸）；总图宽 = PANEL_W * 列数
    "panel_w": 3.3,
    "panel_h": 2.5,
    "ncols": None,          # None = 一行排开；给整数则换行
    "wspace": 0.55,         # 子图横向间距（右轴标签要占地方，别调太小）
    "hspace": 0.45,

    "gr_lw": 1.0,
    "cn_lw": 0.9,
    "cn_dash": (3, 1.5),    # CN 的虚线样式

    "xlabel": r"$r$ (\AA)",
    "ylabel_gr": r"$g(r)$",
    "ylabel_cn": r"$CN(r)$",
    "title_fmt": "{center}--{neighbor}",   # 子图标题
    "legend_loc": "upper right",
    "legend_ncol": 1,

    "xlim": (0, None),      # None = 自动
    "gr_ylim": (0, None),
    "cn_ylim": (0, None),
}


# ── 绘图 ─────────────────────────────────────────────────────────────────────

def pair_of(df):
    """子图标题取自**数据列**而不是文件名。

    文件名可能被改过，`center`/`neighbor` 列不会。长表里一个 gr 产物通常只含一个
    配对（`-a P -b O` 时），不带筛选跑出来的 `gr_all.csv` 则含全部有序对——那种
    情况下标题说不清是哪一对，直接拒绝，让调用者去指定配对重跑。
    """
    pairs = df[["center", "neighbor"]].drop_duplicates()
    if len(pairs) != 1:
        raise SystemExit(
            f"这份产物含 {len(pairs)} 个配对，一张子图画不下。"
            f"用 `ferro traj gr -a CENTRE -b NEIGHBOUR` 重跑，或先在 pandas 里筛。"
        )
    row = pairs.iloc[0]
    return str(row["center"]), str(row["neighbor"])


def build_panel(df, ax, colors):
    """在 ax 上画一份产物：左轴 g(r) 实线，右轴 CN(r) 虚线。"""
    center, neighbor = pair_of(df)

    # science 样式默认在右侧画主轴的镜像刻度，会跟 CN 轴的刻度叠住
    ax.tick_params(right=False)
    ax_cn = ax.twinx()
    ax_cn.tick_params(right=True)

    for name in fp.files_in(df):
        sub = df[df["file"] == name].sort_values("r")
        c = colors[name]
        ax.plot(sub["r"], sub["gr"], color=c, lw=CFG["gr_lw"], label=name)
        ax_cn.plot(sub["r"], sub["cn"], color=c, lw=CFG["cn_lw"],
                   dashes=CFG["cn_dash"])

    ax.set_xlabel(CFG["xlabel"])
    ax.set_ylabel(CFG["ylabel_gr"])
    ax_cn.set_ylabel(CFG["ylabel_cn"])
    ax.set_title(CFG["title_fmt"].format(center=center, neighbor=neighbor))
    ax.set_xlim(*CFG["xlim"])
    ax.set_ylim(*CFG["gr_ylim"])
    ax_cn.set_ylim(*CFG["cn_ylim"])
    return ax_cn


def build_figure(frames):
    """frames: [(path, df)]，返回 fig。"""
    n = len(frames)
    ncols = CFG["ncols"] or n
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(CFG["panel_w"] * ncols, CFG["panel_h"] * nrows),
        squeeze=False,
    )
    fig.subplots_adjust(wspace=CFG["wspace"], hspace=CFG["hspace"])

    # 颜色按全部产物里出现过的 file 统一分配，保证跨子图一致
    names = list(dict.fromkeys(sum((fp.files_in(df) for _, df in frames), [])))
    colors = fp.color_map(names)

    flat = axes.ravel()
    for ax, (_, df) in zip(flat, frames):
        build_panel(df, ax, colors)
    for ax in flat[n:]:
        ax.set_visible(False)

    # 图例只画一次：曲线的颜色语义在所有子图里相同
    handles, labels = flat[0].get_legend_handles_labels()
    flat[0].legend(handles, labels, loc=CFG["legend_loc"], ncol=CFG["legend_ncol"])
    return fig


# ── 入口 ─────────────────────────────────────────────────────────────────────

def expand(patterns):
    """自展开 glob，与 ferro 的 `-i` 同义（省得依赖 shell）。"""
    out = []
    for p in patterns:
        hits = sorted(glob.glob(p)) or ([p] if Path(p).is_file() else [])
        if not hits:
            raise SystemExit(f"没有文件匹配 '{p}'")
        out.extend(Path(h) for h in hits)
    return out


def main():
    ap = argparse.ArgumentParser(description="g(r) / CN(r) 发表级绘图")
    ap.add_argument("inputs", nargs="+", help="ferro traj gr 的 csv（可用 glob）")
    ap.add_argument("-o", "--output", default="gr", help="产物文件名 stem")
    fp.add_outdir_arg(ap)
    args = ap.parse_args()

    fp.use_style()
    paths = expand(args.inputs)
    frames = [(p, fp.read(p)) for p in paths]
    print(f"Inputs: {len(frames)} file(s)")

    fig = build_figure(frames)
    fp.save(fig, args.output, args.outdir)


if __name__ == "__main__":
    main()
