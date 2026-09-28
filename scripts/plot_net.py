#!/usr/bin/env python
r"""把 `ferro net` 的产物画成发表级的网络拓扑图。

两类图，形状相同：**x = 成分（`file` 列），100 % 堆积柱，条带 = 位点类型**。
问的是同一个问题——某个位点的分布随成分怎么变——所以用同一套模板。

    qn   network_qn.csv           一个形成子一张图，条带 = Q0…Q4
    cn   network_coordination.csv 一个元素一张图，条带 = 配位数

堆积柱相对折线的好处是它把「和为 1」画成了图形约束：读者一眼看到此消彼长，不必在
脑子里再加一次。代价是**小分量看不出趋势**——占 2 % 的条带薄得只剩一条线。真要追
小分量的变化，把那一列单独拉出来画折线（或用导出的 `<stem>_data.csv` 重绘）。

**口径**：`qn` 只数**同核**连接（P–O–P），即文献 $Q^n_m$ 的 $n$。异核分解
（`network_qn_partner.csv` 的 `m_<X>` 列）不画：色相套明度的两级堆积柱在实物上
难以读懂，2026-09-28 按用户决定删去 `--partner`。那张表仍由 `ferro net` 照常产出，
要看 $Q^n(m\mathrm{Al})$ 在 pandas 里筛。

## 多 csv

`-o` 的 suffix 落在文件名里，所以 `network_qn_CMD.csv` 与 `network_qn_MLMD.csv` 天然
区分。给多个 csv 就是**同一成分刻度下并排多根柱**，组标签取自 suffix（也可在
`labels.csv` 里改名或隐藏）。

用法：
    python plot_net.py qn network_qn.csv --outdir figs
    python plot_net.py cn network_coordination.csv --element Al,Zn --outdir figs
    python plot_net.py qn network_qn_CMD.csv network_qn_MLMD.csv --outdir figs
"""

import argparse

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import ferroplot as fp
from plot_gr import expand

# ── 配置：常改的量都在这里 ───────────────────────────────────────────────────

CFG = {
    "panel_w": 3.3,
    "panel_h": 2.5,
    "ncols": None,          # None = 一行排开
    "wspace": 0.35,
    "hspace": 0.55,

    "bar_width": 0.55,      # 一个成分刻度上所有柱子合计占的宽度
    "bar_gap": 0.08,        # 并排多组时组间空隙（占 bar_width 的比例）
    "edge_lw": 0.4,         # 每个条带的描边
    "edge_color": "white",

    "xlabel": "",           # 成分名已在刻度上，默认不重复
    "ylabel": r"Fraction (\%)",
    "title_qn": "{former}",
    "title_cn": "{element}",
    "xtick_rotation": 30,
    # 多 csv 并排时组标签（CMD / MLMD）竖排在柱顶，标题得让开
    "title_pad_grouped": 22,
    "legend_loc": "upper center",
    # 图例挂在**整张图**下方而不是某个子图下方：成分名是旋转的长标签，
    # 挂在 ax 上会被它们顶穿。y 要留够旋转标签的高度
    "legend_bbox": (0.5, -0.16),
    "legend_ncol": 6,

    "ylim": (0, 100),
}


# ── 数据整形 ─────────────────────────────────────────────────────────────────

def series_qn(df, former):
    """返回 [(段标签, 类别值, 该段的行)]，按 qn 升序（堆叠自下而上）。"""
    sub = df[df["former"] == former]
    return [(f"Q{int(q)}", int(q), sub[sub["qn"] == q])
            for q in sorted(sub["qn"].unique())]


def series_cn(df, element):
    """段标签只写配位数,不带元素符号。

    颜色是**全图统一按类别值分配**的,所以图例是共享的一份;元素名在子图标题上。
    若标签写成 `Al_4`,Zn 那一格里同一个橙色就会挂着 `Al` 的名字。
    """
    sub = df[df["element"] == element]
    return [(f"CN = {int(c)}", int(c), sub[sub["cn"] == c])
            for c in sorted(sub["cn"].unique())]


def nontrivial_elements(df):
    """只有一个配位档的元素跳过——那是一条平的 100 %，是噪音不是信息。

    与 `ferro net` 对空 Qn 表的处理同一条原则：不产出「测了但没内容」的东西。
    """
    keep, skipped = [], []
    for e in dict.fromkeys(df["element"]):
        (keep if df[df["element"] == e]["cn"].nunique() > 1 else skipped).append(e)
    return keep, skipped


# ── 绘图 ─────────────────────────────────────────────────────────────────────

def build_panel(ax, groups, title, colors, handles, seen, ctx):
    """groups: [(组名, segs)]，segs 来自 series_*。

    一个 x 刻度 = 一个成分；刻度上每个组一根柱；柱内自下而上堆叠各段。

    `colors` 是**跨全图**的 {类别值: 颜色}，`handles`/`seen` 也跨全图累加：
    颜色的含义必须在每一格里相同，否则一份共享图例就是在撒谎。
    """
    files = ctx.names
    n_grp = len(groups)
    total_w = CFG["bar_width"]
    w = total_w / n_grp * (1 - CFG["bar_gap"])
    offsets = [(-total_w / 2 + total_w / n_grp * (i + 0.5)) for i in range(n_grp)]

    for gi, (gname, segs) in enumerate(groups):
        xs = [i + offsets[gi] for i in range(len(files))]
        bottom = [0.0] * len(files)
        for label, key, rows in segs:
            vals = [float(rows[rows["file"] == f]["fraction"].sum()) * 100 for f in files]
            ax.bar(xs, vals, width=w, bottom=bottom, color=colors[key],
                   edgecolor=CFG["edge_color"], linewidth=CFG["edge_lw"], zorder=2)
            bottom = [b + v for b, v in zip(bottom, vals)]
            for f, v in zip(files, vals):
                if (rows["file"] == f).any():
                    ctx.record(panel=title, group=gname, file=f, category=label,
                               percent=round(v, 6), count=int(rows[rows["file"] == f]["count"].sum()))
            if key not in seen:
                seen.add(key)
                handles.append((key, Patch(facecolor=colors[key], label=label)))

        if n_grp > 1:
            ax.text(offsets[gi], 101, ctx.label(gname), ha="center", va="bottom",
                    fontsize=plt.rcParams["font.size"] - 1, rotation=90)

    ax.set_xticks(range(len(files)))
    ax.set_xticklabels([ctx.label(f) for f in files],
                       rotation=CFG["xtick_rotation"], ha="right")
    ax.set_xlabel(CFG["xlabel"])
    ax.set_ylabel(CFG["ylabel"])
    ax.set_ylim(*CFG["ylim"])
    ax.set_title(title, pad=CFG["title_pad_grouped"] if n_grp > 1 else None)
    ax.tick_params(top=False)          # 分类轴上的镜像刻度没有意义


def build_figure(panels, groups_of, ctx):
    """panels: [(title, key)]，groups_of(key) -> [(组名, segs)]。"""
    n = len(panels)
    ncols = CFG["ncols"] or n
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(CFG["panel_w"] * ncols, CFG["panel_h"] * nrows),
        squeeze=False,
    )
    fig.subplots_adjust(wspace=CFG["wspace"], hspace=CFG["hspace"])

    per_panel = [groups_of(key) for _, key in panels]
    keys = sorted({k for groups in per_panel for _, segs in groups for _, k, _ in segs})
    colors = fp.color_map(keys)

    flat = axes.ravel()
    handles, seen = [], set()
    for ax, (title, _), groups in zip(flat, panels, per_panel):
        build_panel(ax, groups, title, colors, handles, seen, ctx)
    for ax in flat[n:]:
        ax.set_visible(False)
    # 图例按类别值排序，与堆叠自下而上的次序一致
    handles = [h for _, h in sorted(handles, key=lambda kh: kh[0])]

    fig.legend(handles=handles, loc=CFG["legend_loc"],
               bbox_to_anchor=CFG["legend_bbox"], ncol=CFG["legend_ncol"],
               bbox_transform=fig.transFigure)
    return fig


# ── 入口 ─────────────────────────────────────────────────────────────────────

def run_qn(frames, ctx, args):
    formers = list(dict.fromkeys(f for _, df in frames for f in df["former"]))

    def groups_of(former):
        return [(fp.suffix_of(p), series_qn(df, former)) for p, df in frames]

    return build_figure([(CFG["title_qn"].format(former=f), f) for f in formers],
                        groups_of, ctx)


def run_cn(frames, ctx, args):
    if args.element:
        elements = [e.strip() for e in args.element.split(",") if e.strip()]
    else:
        elements, skipped = [], []
        for _, df in frames:
            k, s = nontrivial_elements(df)
            elements += [e for e in k if e not in elements]
            skipped += [e for e in s if e not in skipped]
        for e in skipped:
            if e not in elements:
                print(f"Note  : {e} 只有一个配位档，跳过（画出来是一条平的 100 %）")
    if not elements:
        raise SystemExit("没有分布不平凡的元素可画；用 --element 强制指定")

    def groups_of(element):
        return [(fp.suffix_of(p), series_cn(df, element)) for p, df in frames]

    return build_figure([(CFG["title_cn"].format(element=e), e) for e in elements],
                        groups_of, ctx)


RUN = {"qn": run_qn, "cn": run_cn}


def main():
    ap = argparse.ArgumentParser(description="ferro net 发表级绘图")
    ap.add_argument("kind", choices=list(RUN), help="qn = Qn 分布；cn = 配位数分布")
    ap.add_argument("inputs", nargs="+", help="ferro net 的 csv（可用 glob）")
    ap.add_argument("--element", default=None,
                    help="[cn] 只画这些元素，逗号分隔（默认：分布不平凡的全部）")
    ap.add_argument("-o", "--output", default=None, help="产物文件名 stem")
    fp.add_outdir_arg(ap)
    args = ap.parse_args()

    fp.use_style()
    frames = [(p, fp.read(p)) for p in expand(args.inputs)]
    print(f"Inputs: {len(frames)} file(s)")

    # 组名只在多 csv 并排时显示，单 csv 时它退回整个 stem，登记进去只会弄脏 labels.csv
    groups = [fp.suffix_of(p) for p, _ in frames] if len(frames) > 1 else []
    ctx = fp.context(frames, args.outdir, groups=groups)
    # labels.csv 里 show=0 的组整份不画
    frames = [(p, df) for p, df in frames if not groups or ctx.group_shown(fp.suffix_of(p))]
    if not frames:
        raise SystemExit(f"{ctx.path} 把本次全部组都设成了 show=0，没有可画的")

    fig = RUN[args.kind](frames, ctx, args)
    fp.save(fig, args.output or f"net_{args.kind}", args.outdir, ctx)


if __name__ == "__main__":
    main()
