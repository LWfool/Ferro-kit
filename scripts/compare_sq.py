#!/usr/bin/env python
"""比较 `ferro traj sq` 与 dump2sq 计算的结构因子 S(q)（XRD 与中子加权总曲线）。

除下述那个已被实测证实为纯常数的定义差之外，两侧输出**不做任何缩放**，
差异如实呈现，由使用者自行判断。

读图前需要知道的两件事（均为源码事实，不是推断）：

1. 两者对 S(q) 的定义相差一个常数项。dump2sq 的 `CalcSq` 计算的是
   `4πρ/q · Σ_r r[g(r)−1] sin(qr) Δr`，**不含 +1**；ferro 的
   `calc_sq_from_gr` 写作 `1.0 + prefactor * integral`。

   完整轨迹实测：残差在 q>15 处稳定为 1.0003（XRD）/ 1.0002（中子），确认这
   就是一个纯常数偏移。**默认把 dump2sq 的曲线整体 +1** 补成标准 S(q)，两条
   曲线才落在同一基准上、残差图才能反映真实差异而不是被常数淹没。选择加到
   dump2sq 而不是从 ferro 减，是因为 S(q) 的定义要求大 q 极限趋于 1，
   ferro 已是该标准形式。加 `--raw` 可关掉这个补偿、看两侧的原始输出。

2. **本脚本默认那条轨迹上，dump2sq 的输出不可信**，差异不构成对 ferro 的质疑。

   成因是 dump2sq 的一个源码缺陷（早先版本的本文件把它描述成"伪峰"，那只是
   症状不是成因，已更正）：`InitializeType` 只在第 0 帧写入 `data[i].type_new`，
   而 `ReadDataLammpstrj` 的 sscanf 此后只更新 id/type/elem/x/y/z，**从不更新
   type_new**，于是原子类型按**数组下标**而非实际元素归类。examples 下这条
   43Z 轨迹没加 `dump_modify sort id`，帧间同位置同 id 仅 210/2004，第 1 帧起
   partial 归类即失效 —— 全部 10 条 partial 相互相关系数 0.959，XRD 与中子
   总曲线相关 0.9986，Faber-Ziman 加权被完全冲掉。两条总曲线都退化成总 g(r)
   的傅里叶变换（与 FT[dump2sq 自身 total g(r)] 相关 0.9987 / 0.9999，而
   `.gr` 的 Total 列不经过 type_new，正是全文件里唯一干净的一列）。

   反证：换成按 id 排序的 70Z30P00A_NVT（5003/5003），同样两个程序、同样参数，
   ferro 与 dump2sq 的总曲线 max|Δ| 降到 8e-4 ~ 2.7e-3。
   见 scripts/compare_sq_experiment.py。

   脚本每次运行都做这项自检并把结论标在图上。第三行并排的三方 O-O g(r) 是同一
   缺陷在实空间的表现：该处本应严格为 0（O-O 最近邻在 2.48 Å），ferro 与
   dump2analysis 都给 0，只有 dump2sq 给出约 4.5。dump2analysis 每帧按元素重新
   选原子，不受影响。

**0.2.0 的产物形状**（详见 ferrocmp.py 的模块 docstring）：`-o` 是文件名后缀不是
路径，故以 outdir 为工作目录调用；sq 是宽表（`total_xrd` / `total_neutron` /
每对三列），gr 是长表（按 center/neighbor 选行）。

用法:
    python compare_sq.py [--traj FILE] [--outdir DIR] [--raw]
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import ferrocmp as fc
from trajcheck import partial_independence, report_order_stability, scan_traj

# ── 统一计算参数 ──────────────────────────────────────────────────────────────
R_MIN = 0.005   # 两侧显式对齐（ferro 默认 0.001）
R_MAX = 15.0    # 小于 ferro 的最小镜像 clamp 上界（本轨迹 15.16）
DR = 0.01
Q_MIN = 0.1
Q_MAX = 25.0
DQ = 0.05

# dump2sq 的 CalcSq 不含 S(q) 的 +1 项，补上后两侧才在同一基准（见模块 docstring）
SQ_SHIFT_REF = 1.0

REPO = Path(__file__).resolve().parent.parent
DEFAULT_TRAJ = REPO / "examples" / "43Z43P15A_NPT.lammpstrj"


def compute(traj, outdir, ferro_bin, d2sq_bin, d2a_bin):
    print("[S(q)]")
    # ferro：不给配对参数则输出全部 partial + 两条 total
    fe_sq = fc.run_ferro(
        ferro_bin,
        ["traj", "sq", "-i", traj, "-o", "cmp",
         "--q-min", Q_MIN, "--q-max", Q_MAX, "--dq", DQ,
         "--r-min", R_MIN, "--r-max", R_MAX, "--dr", DR,
         "--weighting", "both"],
        outdir, fc.product_name("sq", suffix="cmp"))   # sq 无类型选择,故无 label 段

    # dump2sq：无元素选择参数，一次算出全部 partial 与两条 total
    d2sq_base = outdir / "d2sq"
    d2sq_sq = Path(str(d2sq_base) + ".sq")
    fc.run([d2sq_bin, "-o", str(d2sq_base),
            f"--GRmin={R_MIN}", f"--GRmax={R_MAX}", f"--GRdel={DR}",
            f"--SQmin={Q_MIN}", f"--SQmax={Q_MAX}", f"--SQdel={DQ}",
            str(traj)], d2sq_sq)

    print("[O-O g(r) 归因]")
    # 配对已由 label 段进文件名，-o 与 S(q) 共用同一个批次标记
    fe_gr = fc.run_ferro(
        ferro_bin,
        ["traj", "gr", "-a", "O", "-b", "O", "-i", traj, "-o", "cmp",
         "--r-min", R_MIN, "--r-max", R_MAX, "--dr", DR],
        outdir,
        fc.product_name("gr", label=fc.file_label("O", "O"), suffix="cmp"))
    d2a_gr = outdir / "d2a_O-O.gr"
    fc.run([d2a_bin, "-m", "gr", "-x", "O", "-y", "O",
            "--rmin", str(R_MIN), "--rmax", str(R_MAX), "--dr", str(DR),
            "-i", str(traj), "-o", str(d2a_gr)], d2a_gr)

    # ── S(q)：ferro 是宽表，直接按列名取 ────────────────────────────────────
    fe = fc.one_file(fc.read(fe_sq), fe_sq)
    ref = fc.load_legacy(d2sq_sq)
    ref_names = fc.legacy_header(d2sq_sq, "Sq-XRD")

    q_fe = fc.column(fe, fe_sq, "q")
    n = fc.same_grid(q_fe, ref[:, 0], "S(q) q 网格")

    sq = {
        "q": q_fe[:n],
        "XRD": (fc.column(fe, fe_sq, "total_xrd")[:n],
                ref[:n, ref_names.index("Sq-XRD")]),
        "Neutron": (fc.column(fe, fe_sq, "total_neutron")[:n],
                    ref[:n, ref_names.index("Sq-ND")]),
        # dump2sq 的 partial 列（Sq-XRD 之前的全部），用于归类失效体检
        "partials": ref[:n, 1:ref_names.index("Sq-XRD")],
        "version": fc.version(fe_sq),
        "n_fe": len(fe),
        "n_ref": len(ref),
    }

    # ── O-O g(r)：ferro 是长表，按 center/neighbor 选行 ──────────────────────
    g_fe = fc.pick(fc.one_file(fc.read(fe_gr), fe_gr), fe_gr, "r",
                   center="O", neighbor="O")
    g_d2a = fc.load_legacy(d2a_gr)
    g_d2a = g_d2a[g_d2a[:, 0] > 0]          # 去掉 dump2analysis 固定输出的 r=0 行
    d2sq_gr = Path(str(d2sq_base) + ".gr")
    g_ref = fc.load_legacy(d2sq_gr)
    g_ref_names = fc.legacy_header(d2sq_gr, "Total")

    m = min(len(g_fe), len(g_d2a), len(g_ref))
    gr = {
        "r": g_fe["r"].to_numpy()[:m],
        "fe": g_fe["gr"].to_numpy()[:m],
        "d2a": g_d2a[:m, 1],
        "d2sq": g_ref[:m, g_ref_names.index("O-O")],
    }
    return sq, gr


def plot(sq, gr, traj, png, shift, order_ok):
    fig, axes = plt.subplots(3, 2, figsize=(13, 12))
    q = sq["q"]
    ref_label = "dump2sq{}".format(f" + {shift:g}" if shift else "")

    for col, key in enumerate(["XRD", "Neutron"]):
        fe, ref = sq[key]
        ref = ref + shift
        diff = fe - ref

        ax = axes[0][col]
        ax.plot(q, fe, lw=1.0, label=f"ferro  total_{key.lower()}")
        ax.plot(q, ref, lw=1.0, ls="--",
                label=f"{ref_label}  Sq-{'XRD' if key == 'XRD' else 'ND'}")
        ax.axhline(1.0, color="gray", lw=0.5, ls=":")
        note = f"dump2sq shifted by +{shift:g}" if shift else "raw output, no shift"
        ax.set_title(f"{key}-weighted total S(q)   ({note})", fontsize=11)
        ax.set_ylabel("S(q)")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

        axd = axes[1][col]
        axd.plot(q, diff, lw=0.9, color="crimson")
        axd.axhline(0.0, color="k", lw=0.6)
        tail = diff[q > 15]
        axd.set_title(f"ferro − {ref_label}    mean over q>15 = {tail.mean():.5f}   "
                      f"max|Δ| = {np.abs(diff).max():.4f}", fontsize=9)
        axd.set_ylabel("Δ S(q)")
        axd.set_xlabel("q [Å⁻¹]")
        axd.grid(alpha=0.3)

    # ── 第三行：上游 O-O g(r) 三方对照，用于归因 ─────────────────────────────
    r = gr["r"]
    for col, (lo, hi, title) in enumerate([
        (0.0, R_MAX, "O-O g(r), full range"),
        (1.0, 3.2, "O-O g(r), zoom on the 1.48 Å region"),
    ]):
        ax = axes[2][col]
        sel = (r >= lo) & (r <= hi)
        ax.plot(r[sel], gr["fe"][sel], lw=1.2, label="ferro")
        ax.plot(r[sel], gr["d2a"][sel], lw=1.0, ls="--", label="dump2analysis")
        ax.plot(r[sel], gr["d2sq"][sel], lw=1.0, ls="-.", label="dump2sq (internal)")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("r [Å]")
        ax.set_ylabel("g(r)")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        if col == 1:
            k = int(np.argmin(np.abs(r - 1.48)))
            ax.annotate(
                f"at r=1.48 Å\nferro {gr['fe'][k]:.3f} / d2a {gr['d2a'][k]:.3f} / "
                f"dump2sq {gr['d2sq'][k]:.3f}",
                xy=(1.48, gr["d2sq"][k]), xytext=(1.9, gr["d2sq"][k] * 0.9),
                fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.8))

    shift_note = (f"dump2sq shifted by +{shift:g} (its CalcSq omits the +1 of S(q))"
                  if shift else "raw values, no shift")
    fig.suptitle(
        f"S(q):  {sq['version']}  vs  dump2sq\n"
        f"{Path(traj).name}   |   r: {R_MIN}–{R_MAX} step {DR}   "
        f"q: {Q_MIN}–{Q_MAX} step {DQ}   |   {shift_note}",
        fontsize=12)
    fig.text(0.5, 0.005,
             "Only the constant offset is compensated. The remaining deviation comes from "
             "dump2sq's own g(r) (bottom row), not from the transform.",
             ha="center", fontsize=8, style="italic")

    # 未排序的 dump 会让 dump2sq 的归类从第 1 帧起失效，此时整张图的"差异"都不
    # 应读作 ferro 的问题 —— 直接印在图上，避免脱离终端输出后被误读
    if not order_ok:
        fig.text(0.5, 0.955,
                 "⚠ This dump is not written in a stable atom order: dump2sq classifies atoms by "
                 "array index (type_new is set on frame 0 only),\nso its partials and both weighted "
                 "totals are meaningless here. The differences below are NOT evidence against ferro.",
                 ha="center", va="top", fontsize=9, color="crimson", weight="bold")
    fig.tight_layout(rect=(0, 0.02, 1, 0.94 if order_ok else 0.90))
    fig.savefig(png, dpi=140)
    print(f"\n图已写入 {png}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--traj", default=str(DEFAULT_TRAJ))
    p.add_argument("--outdir", default="./cmp_out")
    p.add_argument("--ferro", default="ferro")
    p.add_argument("--dump2sq", default="dump2sq")
    p.add_argument("--dump2analysis", default="dump2analysis")
    p.add_argument("--raw", action="store_true",
                   help=f"不给 dump2sq 补 +{SQ_SHIFT_REF:g}，显示两侧的原始输出")
    args = p.parse_args()
    shift = 0.0 if args.raw else SQ_SHIFT_REF

    traj = Path(args.traj).resolve()
    if not traj.exists():
        sys.exit(f"轨迹不存在: {traj}")
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"轨迹   : {traj}")
    print(f"输出   : {outdir}")
    print(f"参数   : r {R_MIN}–{R_MAX} step {DR}   q {Q_MIN}–{Q_MAX} step {DQ}")
    print(f"对齐   : {'无（--raw）' if args.raw else f'dump2sq + {SQ_SHIFT_REF:g}'}\n")

    _, _, _, stable = scan_traj(traj)
    order_ok = report_order_stability(stable)
    print()

    sq, gr = compute(traj, outdir, args.ferro, args.dump2sq, args.dump2analysis)

    if "partials" in sq:
        c = partial_independence(sq["partials"][sq["q"] > 2.0])
        print(f"dump2sq partial 平均非对角相关 (q>2) = {c:+.3f}"
              "    （接近 1 说明归类已失效）")

    q = sq["q"]
    # 先报告未补偿的残差尾部 —— 它是「偏移确为常数 1」这一判断的依据
    print("\n常数偏移核验（未补偿时的 fe − dump2sq）")
    for key in ["XRD", "Neutron"]:
        fe, ref = sq[key]
        d = fe - ref
        print(f"  {key:<8} q>15 均值 = {d[q > 15].mean():.5f}   "
              f"q>15 标准差 = {d[q > 15].std():.5f}")

    tag = "原始值" if args.raw else f"dump2sq 已 +{SQ_SHIFT_REF:g}"
    print(f"\nS(q) 数值摘要（{tag}）")
    print(f"{'weighting':>10}{'fe 均值':>11}{'ref 均值':>11}{'Δ 全域均值':>13}"
          f"{'Δ (q>15)':>11}{'max|Δ|':>11}")
    for key in ["XRD", "Neutron"]:
        fe, ref = sq[key]
        ref = ref + shift
        d = fe - ref
        print(f"{key:>10}{fe.mean():>11.4f}{ref.mean():>11.4f}{d.mean():>13.4f}"
              f"{d[q > 15].mean():>11.5f}{np.abs(d).max():>11.4f}")

    r, k = gr["r"], int(np.argmin(np.abs(gr["r"] - 1.48)))
    print(f"\nO-O g(r) 在 r={r[k]:.2f} Å（O-O 最近邻在 2.48 Å，此处应为 0）")
    print(f"  ferro         {gr['fe'][k]:.4f}")
    print(f"  dump2analysis {gr['d2a'][k]:.4f}")
    print(f"  dump2sq       {gr['d2sq'][k]:.4f}")

    plot(sq, gr, traj, outdir / "compare_sq.png", shift, order_ok)


if __name__ == "__main__":
    main()
