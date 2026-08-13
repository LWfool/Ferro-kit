#!/usr/bin/env python
"""把 `ferro traj sq` 的产物画成发表级的结构因子图。

**一个 csv 一行两格**：左 $S^{N}(Q)$（中子），右 $S^{X}(Q)$（X 射线）。两条是不同
探针下的量，不共享刻度，所以分格而不是叠在一起。给多个 csv 就是 N 行 × 2 格。

画的是两条 total。加权 partial（每对三列）是能加回 total 的诊断分解，不在这里画——
`ferro traj sq` 恒定把它们写在同一份 csv 里，要看在 pandas 里选列即可。

用法：
    python plot_sq.py sq.csv --outdir figs
"""

import argparse

import matplotlib.pyplot as plt

import ferroplot as fp
from plot_gr import expand

# ── 配置：常改的量都在这里 ───────────────────────────────────────────────────

CFG = {
    "panel_w": 3.3,
    "panel_h": 2.5,
    "wspace": 0.35,
    "hspace": 0.45,

    "lw": 1.0,

    "xlabel": r"$Q$ (\AA$^{-1}$)",
    # 左格在前：列顺序即绘制顺序
    "panels": [
        ("total_neutron", r"$S^{N}(Q)$"),
        ("total_xrd",     r"$S^{X}(Q)$"),
    ],
    "legend_loc": "upper right",
    "legend_ncol": 1,

    "xlim": (0, None),
    "ylim": (None, None),
    "hline": 1.0,           # S(Q) 大 Q 处趋于 1；None 则不画这条参考线
    "hline_kw": {"color": "0.6", "lw": 0.5, "dashes": (2, 2), "zorder": 0},
}


# ── 绘图 ─────────────────────────────────────────────────────────────────────

def build_row(df, axes, colors, path):
    """在一行的两个格子上画一份产物。"""
    for ax, (col, ylabel) in zip(axes, CFG["panels"]):
        if col not in df.columns:
            raise SystemExit(
                f"{path} 里没有 '{col}' 列——那一跑的 --weighting 排除了它。"
                f"用 --weighting both 重跑。"
            )
        if CFG["hline"] is not None:
            ax.axhline(CFG["hline"], **CFG["hline_kw"])
        for name in fp.files_in(df):
            sub = df[df["file"] == name].sort_values("q")
            ax.plot(sub["q"], sub[col], color=colors[name], lw=CFG["lw"], label=name)
        ax.set_xlabel(CFG["xlabel"])
        ax.set_ylabel(ylabel)
        ax.set_xlim(*CFG["xlim"])
        ax.set_ylim(*CFG["ylim"])


def build_figure(frames):
    nrows = len(frames)
    ncols = len(CFG["panels"])
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(CFG["panel_w"] * ncols, CFG["panel_h"] * nrows),
        squeeze=False,
    )
    fig.subplots_adjust(wspace=CFG["wspace"], hspace=CFG["hspace"])

    names = list(dict.fromkeys(sum((fp.files_in(df) for _, df in frames), [])))
    colors = fp.color_map(names)

    for row, (path, df) in zip(axes, frames):
        build_row(df, row, colors, path)

    handles, labels = axes[0][0].get_legend_handles_labels()
    axes[0][0].legend(handles, labels, loc=CFG["legend_loc"], ncol=CFG["legend_ncol"])
    return fig


def main():
    ap = argparse.ArgumentParser(description="结构因子 S(Q) 发表级绘图")
    ap.add_argument("inputs", nargs="+", help="ferro traj sq 的 csv（可用 glob）")
    ap.add_argument("-o", "--output", default="sq", help="产物文件名 stem")
    fp.add_outdir_arg(ap)
    args = ap.parse_args()

    fp.use_style()
    frames = [(p, fp.read(p)) for p in expand(args.inputs)]
    print(f"Inputs: {len(frames)} file(s)")

    fig = build_figure(frames)
    fp.save(fig, args.output, args.outdir)


if __name__ == "__main__":
    main()
