#!/usr/bin/env python
r"""把 `ferro traj msd` 的产物画成发表级的 MSD(t) 图。

一个 csv 一张子图（`msd_all.csv msd_Zn.csv` → 两张），子图里一条轨迹一条曲线。
给了 `--fit-range` 的产物再叠一条同色虚线：拟合直线，**只画在拟合窗口内**，
图例带 $D \pm$ 误差（cm²/s）。

拟合不在这里做。斜率、截距、D、误差全部读 ferro 写在 `[inputs]` 清单里的值——
这份脚本只画，算法只有 ferro 一份。头部共享区只有全批一致的参数（窗口比例、
元素选择），逐文件的数一律在 `[inputs]`。

`--loglog` 是 MDAnalysis 推荐的诊断：双对数下扩散区斜率为 1，拟合窗口应落在
这一段里。参考线锚在第一条轨迹拟合窗口的起点。

用法：
    python plot_msd.py msd_all.csv --outdir figs
    python plot_msd.py msd_all.csv msd_Zn.csv --components --outdir figs   # x/y/z 分量
    python plot_msd.py msd_all.csv --loglog --outdir figs -o msd_loglog
"""

import argparse
import math
from pathlib import Path

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

import ferroplot as fp
from plot_gr import expand

# ── 配置：常改的量都在这里 ───────────────────────────────────────────────────

CFG = {
    "panel_w": 3.3,
    "panel_h": 2.5,
    "ncols": None,          # None = 一行排开；给整数则换行
    "wspace": 0.35,
    "hspace": 0.7,          # 图例挂在每格下方，行距要给它留地方

    # ferro 的时间恒为 fs；这里只换显示单位
    "time_unit": "ps",      # "ps" 或 "fs"
    "lw": 1.0,
    "fit_lw": 1.2,
    "fit_dash": (4, 2),
    "comp_lw": 0.7,
    "comp_styles": {"msd_x": ":", "msd_y": "-.", "msd_z": (0, (1, 3))},
    "ref_color": "0.5",

    "ylabel": r"MSD (\AA$^2$)",
    # 图例挂在子图下方：带 D ± 误差的标签比 3.3 英寸的面板还宽，放在框内哪个角都会
    # 压住曲线（扩散体系占对角线，玻璃态占上沿）。每格 D 不同，不能并成一份 fig 图例
    "legend_anchor": (0.5, -0.22),
    "legend_fontsize": 6,

    "xlim": (0, None),      # 仅线性坐标生效
    "ylim": (0, None),
}

TIME_SCALE = {"fs": 1.0, "ps": 1e-3}


# ── 读取 ─────────────────────────────────────────────────────────────────────

def fit_of(inputs, name):
    """一条轨迹的拟合参数；没给 `--fit-range` 或该行失败时返回 None。"""
    if "slope" not in inputs.columns or name not in inputs.index:
        return None
    row = inputs.loc[name]
    keys = ("t_lo", "t_hi", "slope", "intercept", "d_ang2_per_fs", "d_err", "r2")
    vals = {k: float(row[k]) for k in keys}
    return None if math.isnan(vals["slope"]) else vals


# ── 绘图 ─────────────────────────────────────────────────────────────────────

def d_label(fit):
    r"""`$D = (1.75 \pm 1.45)\times10^{-8}$ cm$^2$/s`，指数取自 D 本身。"""
    d = fit["d_ang2_per_fs"] * 0.1          # Å²/fs → cm²/s
    err = fit["d_err"] * 0.1
    exp = math.floor(math.log10(abs(d))) if d != 0 else 0
    m = d / 10 ** exp
    body = f"{m:.2f}" if math.isnan(err) else f"({m:.2f} \\pm {err / 10 ** exp:.2f})"
    return f"$D = {body}\\times10^{{{exp}}}$ cm$^2$/s"


def build_panel(path, df, ax, ctx, args):
    shared, inputs = fp.read_header(path)
    scale = TIME_SCALE[CFG["time_unit"]]
    first_fit = None

    for name in ctx.shown(df):
        sub = df[df["file"] == name].sort_values("time")
        if args.loglog:
            # t=0 处 MSD=0，双对数下是 -inf
            sub = sub[sub["time"] > 0]
        c = ctx.colors[name]
        fit = fit_of(inputs, name)
        label = ctx.label(name) if fit is None else f"{ctx.label(name)}, {d_label(fit)}"
        ax.plot(sub["time"] * scale, sub["msd"], color=c, lw=CFG["lw"], label=label)
        ctx.record_curve(path.stem, name, "msd", sub["time"] * scale, sub["msd"])

        if args.components:
            for col, ls in CFG["comp_styles"].items():
                ax.plot(sub["time"] * scale, sub[col], color=c,
                        lw=CFG["comp_lw"], ls=ls, alpha=0.8)
                ctx.record_curve(path.stem, name, col, sub["time"] * scale, sub[col])

        if fit is not None:
            t = np.linspace(fit["t_lo"], fit["t_hi"], 50)
            # 压在含噪曲线上面，否则同色虚线会被淹掉
            ax.plot(t * scale, fit["slope"] * t + fit["intercept"], color=c,
                    lw=CFG["fit_lw"], dashes=CFG["fit_dash"], zorder=3)
            ctx.record_curve(path.stem, name, "fit", t * scale,
                             fit["slope"] * t + fit["intercept"])
            if first_fit is None:
                first_fit = (fit, sub)

    if args.loglog:
        ax.set_xscale("log")
        ax.set_yscale("log")
        add_slope_one(ax, first_fit, df, scale)
    else:
        ax.set_xlim(*CFG["xlim"])
        ax.set_ylim(*CFG["ylim"])

    sel = shared.get("elements", "all")
    ax.set_title("all atoms" if sel == "all" else sel)
    ax.set_xlabel(f"$t$ ({CFG['time_unit']})")
    ax.set_ylabel(CFG["ylabel"])
    handles, _ = ax.get_legend_handles_labels()
    if args.components:
        # 分量的线型条目并进同一份图例（第二份 legend 会顶掉第一份）：
        # 颜色仍是轨迹，灰色线型只说明方向
        handles += [mlines.Line2D([], [], color="0.3", lw=CFG["comp_lw"], ls=ls,
                                  label=f"{col.split('_')[1]} component")
                    for col, ls in CFG["comp_styles"].items()]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=CFG["legend_anchor"],
              fontsize=CFG["legend_fontsize"], frameon=False)


def add_slope_one(ax, first_fit, df, scale):
    """斜率 1 参考线，最多跨一个数量级，不越过数据的时间范围。

    锚点：第一条轨迹拟合窗口的起点——窗口若选在扩散区，曲线应沿着这条线走。
    没有拟合时退回第一条轨迹时间轴的几何中点。锚点离末端不足 2 倍时向前挪，
    否则参考线短得看不出斜率。
    """
    if first_fit is not None:
        fit, sub = first_fit
        t0 = fit["t_lo"]
    else:
        sub = df[(df["file"] == fp.files_in(df)[0]) & (df["time"] > 0)].sort_values("time")
        t0 = math.sqrt(sub["time"].iloc[0] * sub["time"].iloc[-1])
    t_min, t_max = sub["time"].iloc[0], sub["time"].iloc[-1]
    if min(10 * t0, t_max) < 2 * t0:
        t0 = max(t_max / 10, t_min)
    t1 = min(10 * t0, t_max)
    y0 = float(np.interp(t0, sub["time"], sub["msd"]))
    if not (t0 > 0 and y0 > 0 and t1 > t0):
        return
    ax.plot(np.array([t0, t1]) * scale, [y0, y0 * t1 / t0],
            color=CFG["ref_color"], lw=0.8, ls="--")
    ax.annotate("slope 1", (t1 * scale, y0 * t1 / t0), fontsize=6,
                color=CFG["ref_color"], ha="right", va="bottom")


def build_figure(frames, ctx, args):
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

    flat = axes.ravel()
    for ax, (path, df) in zip(flat, frames):
        build_panel(path, df, ax, ctx, args)
    for ax in flat[n:]:
        ax.set_visible(False)
    return fig


# ── 入口 ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="MSD(t) 发表级绘图")
    ap.add_argument("inputs", nargs="+", help="ferro traj msd 的 csv（可用 glob）")
    ap.add_argument("-o", "--output", default="msd", help="产物文件名 stem")
    ap.add_argument("--components", action="store_true",
                    help="叠画笛卡尔 x/y/z 分量（看各向异性；三者之和等于总量）")
    ap.add_argument("--loglog", action="store_true",
                    help="双对数坐标 + 斜率 1 参考线（检查拟合窗口是否在扩散区）")
    fp.add_outdir_arg(ap)
    args = ap.parse_args()

    fp.use_style()
    paths = expand(args.inputs)
    frames = [(p, fp.read(p)) for p in paths]
    print(f"Inputs: {len(frames)} file(s)")

    ctx = fp.context(frames, args.outdir)
    fig = build_figure(frames, ctx, args)
    fp.save(fig, args.output, args.outdir, ctx)


if __name__ == "__main__":
    main()
