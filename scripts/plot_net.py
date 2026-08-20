#!/usr/bin/env python
r"""把 `ferro net` 的产物画成发表级的网络拓扑图。

两类图，形状相同：**x = 成分（`file` 列），100 % 堆积柱，条带 = 位点类型**。
问的是同一个问题——某个位点的分布随成分怎么变——所以用同一套模板。

    qn   network_qn.csv           一个形成子一张图，条带 = Q0…Q4
    cn   network_coordination.csv 一个元素一张图，条带 = 配位数

堆积柱相对折线的好处是它把「和为 1」画成了图形约束：读者一眼看到此消彼长，不必在
脑子里再加一次。代价是**小分量看不出趋势**——占 2 % 的条带薄得只剩一条线。真要追
小分量的变化，把那一列单独拉出来画折线。

## `--partner`（仅 qn）

展开成 `network_qn_partner.csv` 的两级结构：色相仍是 Qn，同色系内按 m_<X> 分明度，
每个 Qn 组外加一圈框线标出粗结构边界。这样第一眼仍落在 Qn 的色块上，伙伴分解只是
补充。段标签即文献记号 $Q^n(m\mathrm{X})$——`Q1(2Al)` = 1 个 P–O–P + 2 个 P–O–Al。

**口径（2026-08-20 起）**：`qn` 只数**同核**连接（P–O–P），`m_<X>` 是**异核**连接
（P–O–Al…），总连接数 = `qn + Σm`。形成子自身那一列已从 csv 里去掉（它恒等于 `qn`），
所以剩下的每一列都携带 `qn` 没有的信息：

- **退化不再发生**。旧口径下单形成子体系 `m_P ≡ qn`，展开只画出一条假的两级结构，
  这是 `--partner` 当初被设成开关而非默认的唯一理由；那个理由现在没有了。
  是否改回默认属出图形式问题，另行决定。
- **无异核形成子的体系没有 m_<X> 列**（如 Zn–P–O），此时 `--partner` 自动退化为普通
  Qn 图并打印说明——不是错误，那个体系确实没有伙伴维可展开。
- **三簇配体按连接数计入**：它把本位点连上两个伙伴，故贡献 2。于是 `qn + Σm` 可能
  **超过**该位点的桥氧个数（差额即三簇桥），两者在 `ferro net` 里分别是
  `mean_qn`/`m_` 与 `mean_n_bo`。

## 多 csv

`-o` 的 suffix 落在文件名里，所以 `network_qn_CMD.csv` 与 `network_qn_MLMD.csv` 天然
区分。给多个 csv 就是**同一成分刻度下并排多根柱**，组标签取自 suffix。

用法：
    python plot_net.py qn network_qn.csv --outdir figs
    python plot_net.py qn network_qn_partner.csv --partner --outdir figs
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
    "group_lw": 0.7,        # --partner 下 Qn 组的框线
    "group_color": "black",

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
    # --partner 下明度维的说明。图例只画 Qn 色相（第一眼要落在粗结构上），
    # 但同色系的深浅档携带 m_<X>，不说明读者无从得知它编码什么
    "partner_note": r"Shading within each hue: {m} high (dark) $\to$ low (light)",
    "note_bbox": (0.5, -0.30),

    "ylim": (0, 100),
}


# ── 数据整形 ─────────────────────────────────────────────────────────────────

def series_qn(df, former, partner):
    """返回 [(段标签, 所属 Qn, 该段的行)]，按堆叠自下而上排序。

    **一个段是一个物种，不是一行**。`partner` 模式下同一个 `(qn, m_<X>…)` 组合在不同
    成分里各有一行，它们必须归入同一段——否则明度档会按行号分配，同一个 m 组合在两个
    成分里拿到不同深浅，颜色就不再编码 m 而是编码"第几行"。

    段序：`partner=False` 按 qn 升序；`True` 时先 qn 升序，组内按第一个 m_ 列**降序**，
    于是同色系里最深的一档 m 最高。

    `m_<X>` 列在新口径下全是**异核**伙伴（形成子自身那一列已从 csv 去掉），故没有
    m_ 列就是真的没有异核伙伴，此时退化为普通 Qn 分布是正确行为而非降级。
    """
    sub = df[df["former"] == former]
    m_cols = [c for c in sub.columns if c.startswith("m_")]

    if not partner or not m_cols:
        return [(f"Q{int(q)}", int(q), sub[sub["qn"] == q])
                for q in sorted(sub["qn"].unique())]

    # 列并集下某个成分缺整列 m_<X>（那个体系根本没有 X），空字段读回来是 NaN。
    # 这里的 NaN 确实等于 0：没有 X 就没有连向 X 的桥。不填的话 `== key` 匹配不上，
    # 整根柱子会静默消失——一般不该把 NaN 当 0，这是有明确语义的例外。
    sub = sub.copy()
    sub[m_cols] = sub[m_cols].fillna(0)

    combos = sub[["qn", *m_cols]].drop_duplicates()
    combos = combos.sort_values(["qn", m_cols[0]], ascending=[True, False])

    segs = []
    for _, key in combos.iterrows():
        mask = sub["qn"] == key["qn"]
        for c in m_cols:
            mask &= sub[c] == key[c]
        # 数字在前：文献写 Q^1(2Al) 而不是 Q^1(Al2)
        tag = ",".join(f"{int(key[c])}{c[2:]}" for c in m_cols if key[c] > 0)
        label = f"Q{int(key['qn'])}" + (f"({tag})" if tag else "")
        segs.append((label, int(key["qn"]), sub[mask]))
    return segs


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

def build_panel(ax, groups, files, title, partner, palette, handles, seen):
    """groups: [(组标签, segs)]，segs 来自 series_*。

    一个 x 刻度 = 一个成分；刻度上每个组一根柱；柱内自下而上堆叠各段。

    `palette` 是**跨全图**的 {类别值: [明度档]}，`handles`/`seen` 也跨全图累加：
    颜色的含义必须在每一格里相同，否则一份共享图例就是在撒谎。
    """
    n_grp = len(groups)
    total_w = CFG["bar_width"]
    w = total_w / n_grp * (1 - CFG["bar_gap"])
    offsets = [(-total_w / 2 + total_w / n_grp * (i + 0.5)) for i in range(n_grp)]
    cycle = {q: shades[0] for q, shades in palette.items()}
    shade_of = palette

    for gi, (gname, segs) in enumerate(groups):
        xs = [i + offsets[gi] for i in range(len(files))]
        bottom = [0.0] * len(files)
        used = {}
        # 每个 Qn 组的上下沿，供框线用
        span = {}
        for label, q, rows in segs:
            k = used.get(q, 0)
            used[q] = k + 1
            color = shade_of[q][k]
            vals = [float(rows[rows["file"] == f]["fraction"].sum()) * 100 for f in files]
            ax.bar(xs, vals, width=w, bottom=bottom, color=color,
                   edgecolor=CFG["edge_color"], linewidth=CFG["edge_lw"], zorder=2)
            lo = list(bottom)
            bottom = [b + v for b, v in zip(bottom, vals)]
            for i in range(len(files)):
                s = span.setdefault((q, i), [lo[i], bottom[i]])
                s[1] = bottom[i]
            key = q if partner else label
            if key not in seen:
                seen.add(key)
                handles.append(Patch(facecolor=cycle[q], label=f"Q{q}" if partner else label))

        if partner:
            # Qn 的粗结构边界：一圈框线，把同一 Qn 的若干明度档圈起来
            for (q, i), (lo, hi) in span.items():
                ax.bar(xs[i], hi - lo, width=w, bottom=lo, facecolor="none",
                       edgecolor=CFG["group_color"], linewidth=CFG["group_lw"], zorder=3)

        if n_grp > 1:
            ax.text(offsets[gi], 101, gname, ha="center", va="bottom",
                    fontsize=plt.rcParams["font.size"] - 1, rotation=90)

    ax.set_xticks(range(len(files)))
    ax.set_xticklabels(files, rotation=CFG["xtick_rotation"], ha="right")
    ax.set_xlabel(CFG["xlabel"])
    ax.set_ylabel(CFG["ylabel"])
    ax.set_ylim(*CFG["ylim"])
    ax.set_title(title, pad=CFG["title_pad_grouped"] if n_grp > 1 else None)
    ax.tick_params(top=False)          # 分类轴上的镜像刻度没有意义


def build_palette(all_groups):
    """跨全图的 {类别值: [明度档]}。

    明度档数取**所有格、所有组**里该类别的最大段数：不同成分下伙伴分解的行数不同
    （某个成分一个 `Q2(m_Al=2)` 都没有），只按第一组数会在别的组上索引越界。
    """
    keys = sorted({q for groups in all_groups for _, segs in groups for _, q, _ in segs})
    cycle = fp.color_map(keys)
    per = {
        q: max((sum(1 for _, qq, _ in segs if qq == q)
                for groups in all_groups for _, segs in groups), default=1) or 1
        for q in keys
    }
    return {q: fp.shades(cycle[q], per[q]) for q in keys}


def build_figure(panels, groups_of, files, partner, note=None):
    """panels: [(title, key)]，groups_of(key) -> [(组标签, segs)]。

    `note` 是图例下方的一行小字，用于解释图例本身表达不了的维度（明度档）。
    """
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
    palette = build_palette(per_panel)

    flat = axes.ravel()
    handles, seen = [], set()
    for ax, (title, _), groups in zip(flat, panels, per_panel):
        build_panel(ax, groups, files, title, partner, palette, handles, seen)
    for ax in flat[n:]:
        ax.set_visible(False)
    # 图例按类别值排序，与堆叠自下而上的次序一致
    handles.sort(key=lambda h: h.get_label())

    fig.legend(handles=handles, loc=CFG["legend_loc"],
               bbox_to_anchor=CFG["legend_bbox"], ncol=CFG["legend_ncol"],
               bbox_transform=fig.transFigure)
    if note:
        fig.text(*CFG["note_bbox"], note, ha="center", va="top",
                 fontsize=plt.rcParams["font.size"] - 1,
                 transform=fig.transFigure)
    return fig


# ── 入口 ─────────────────────────────────────────────────────────────────────

def run_qn(frames, args):
    partner = args.partner
    formers = list(dict.fromkeys(f for _, df in frames for f in df["former"]))
    files = list(dict.fromkeys(f for _, df in frames for f in fp.files_in(df)))

    # 无 m_<X> 列有两种可能：给的是 network_qn.csv，或者这个体系根本没有异核形成子
    # （Zn-P-O 的 qn_partner 与 qn 列结构完全相同）。两者无法从列上区分，而后者是
    # 正常情形，所以退化而不是报错 —— series_qn 已有对应分支
    if partner and not any(c.startswith("m_") for _, df in frames for c in df.columns):
        print("Note  : 没有 m_<X> 列，按普通 Qn 分布绘制。"
              "该体系没有异核形成子，或给的是 network_qn.csv")
        partner = False

    def groups_of(former):
        return [(fp.suffix_of(p), series_qn(df, former, partner)) for p, df in frames]

    note = None
    if partner:
        m_elems = sorted({c[2:] for _, df in frames for c in df.columns
                          if c.startswith("m_")})
        m_tex = ", ".join(rf"$m_{{\mathrm{{{e}}}}}$" for e in m_elems)
        note = CFG["partner_note"].format(m=m_tex)

    return build_figure([(CFG["title_qn"].format(former=f), f) for f in formers],
                        groups_of, files, partner, note)


def run_cn(frames, args):
    files = list(dict.fromkeys(f for _, df in frames for f in fp.files_in(df)))

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
                        groups_of, files, False)


def main():
    ap = argparse.ArgumentParser(description="ferro net 发表级绘图")
    ap.add_argument("kind", choices=["qn", "cn"], help="qn = Qn 分布；cn = 配位数分布")
    ap.add_argument("inputs", nargs="+", help="ferro net 的 csv（可用 glob）")
    ap.add_argument("--partner", action="store_true",
                    help="[qn] 展开 m_<X> 异核伙伴分解，即文献的 Q^n(mX)"
                         "（需 network_qn_partner.csv）")
    ap.add_argument("--element", default=None,
                    help="[cn] 只画这些元素，逗号分隔（默认：分布不平凡的全部）")
    ap.add_argument("-o", "--output", default=None, help="产物文件名 stem")
    fp.add_outdir_arg(ap)
    args = ap.parse_args()

    fp.use_style()
    frames = [(p, fp.read(p)) for p in expand(args.inputs)]
    print(f"Inputs: {len(frames)} file(s)")

    fig = run_qn(frames, args) if args.kind == "qn" else run_cn(frames, args)
    stem = args.output or (f"net_{args.kind}" + ("_partner" if args.partner else ""))
    fp.save(fig, stem, args.outdir)


if __name__ == "__main__":
    main()
