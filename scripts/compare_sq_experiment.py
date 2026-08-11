#!/usr/bin/env python
"""把 `ferro traj sq` 与 dump2sq 算出的 S(q) 和实验数据画在一起。

轨迹默认用 70ZnO-30P2O5（无 Al）的 NVT 轨迹：优先 examples/70Z30P00A_NVT.lammpstrj
（50 帧，8 MB，未入库），缺席时退到 tests/70Z30P00A_NVT_5.lammpstrj（5 帧，已入库）。
实验数据取 tests/Experiment.xlsx 的 Sheet1 —— 该表的四列
nq_70 / SQn_70 / xq_70 / SQx_70 正是同一组分的中子与 X 射线总结构因子。

除下面这一处已被实测证实为纯常数的定义差之外，三方曲线**不做任何缩放或平移**，
差异如实呈现。

**dump2sq 的 +1**：其 `CalcSq` 算的是 `4πρ/q · Σ_r r[g(r)−1] sin(qr) Δr`，不含
S(q) 定义里的 +1；ferro 与实验数据都是 q→∞ 趋于 1 的标准形式。默认给 dump2sq
整条曲线 +1 补到同一基准，`--raw` 可关掉。

**这条轨迹为什么值得单独跑**：dump2sq 只在第 0 帧写入原子类型，dump 未按 id 排序时
其 partial 归类会从第 1 帧起失效 —— 成因、实测数据与反证见 scripts/trajcheck.py。
本轨迹按 id 排序且帧间顺序完全稳定（5003/5003），bug 不触发，因此是一次干净的三方
对照。脚本每次运行都会做这项自检并把结论打印出来。

**0.2.0 的产物形状**（详见 ferrocmp.py 的模块 docstring）：`-o` 是文件名后缀不是
路径，故以 outdir 为工作目录调用；sq 是宽表，按 `total_xrd` / `total_neutron` 取列。

用法:
    python compare_sq_experiment.py [--traj FILE] [--outdir DIR] [--raw]
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl

import ferrocmp as fc
from trajcheck import partial_independence, report_order_stability, scan_traj

# ── 统一计算参数 ──────────────────────────────────────────────────────────────
R_MIN = 0.005   # 两侧显式对齐（ferro 默认 0.001）
DR = 0.01
Q_MIN = 0.1
Q_MAX = 45.0    # 覆盖实验中子数据的完整范围（0.5–45）
DQ = 0.05       # 与实验 q 网格同步长，且 0.1+0.05k 恰好命中实验的 0.5+0.05k

# dump2sq 的 CalcSq 不含 S(q) 的 +1 项（见模块 docstring）
SQ_SHIFT_REF = 1.0

REPO = Path(__file__).resolve().parent.parent
DEFAULT_XLSX = REPO / "tests" / "Experiment.xlsx"

# 优先用完整的 50 帧轨迹（8 MB，未入库）；没有就退到已入库的 5 帧子集。
# 两者结论一致：50 帧 rms(fe−exp)=0.0194/0.0277，5 帧 0.0197/0.0282，
# fe 与 dump2sq 的一致性同样在 1e-3 量级，FSDP 峰位相同。
TRAJ_CANDIDATES = [
    REPO / "examples" / "70Z30P00A_NVT.lammpstrj",
    REPO / "tests" / "70Z30P00A_NVT_5.lammpstrj",
]


def default_traj():
    for p in TRAJ_CANDIDATES:
        if p.exists():
            return p
    return TRAJ_CANDIDATES[-1]

# Experiment.xlsx / Sheet1 的列（0 起）: nq_70 SQn_70 xq_70 SQx_70
EXP_SHEET = "Sheet1"
EXP_COLS = {"neutron": (0, 1), "xrd": (2, 3)}


def load_experiment(xlsx):
    ws = openpyxl.load_workbook(xlsx, data_only=True)[EXP_SHEET]
    grid = np.array([[np.nan if c is None else float(c) for c in row]
                     for row in ws.iter_rows(min_row=2, values_only=True)])
    out = {}
    for key, (cq, cs) in EXP_COLS.items():
        q, s = grid[:, cq], grid[:, cs]
        m = ~(np.isnan(q) | np.isnan(s))
        out[key] = (q[m], s[m])
    return out


def compute(traj, outdir, r_max, ferro_bin, d2sq_bin):
    fe_out = fc.run_ferro(
        ferro_bin,
        ["traj", "sq", "-i", traj, "-o", "exp",
         "--q-min", Q_MIN, "--q-max", Q_MAX, "--dq", DQ,
         "--r-min", R_MIN, "--r-max", r_max, "--dr", DR,
         "--weighting", "both"],
        outdir, "sq_exp.csv")

    # dump2sq 把 argv 的最后一项当输入文件；输出为 <outbase>.{gr,cn,sq,coeff}
    base = outdir / "d2sq"
    ref_out = base.with_suffix(".sq")
    fc.run([d2sq_bin, "-o", str(base),
            f"--GRmin={R_MIN}", f"--GRmax={r_max}", f"--GRdel={DR}",
            f"--SQmin={Q_MIN}", f"--SQmax={Q_MAX}", f"--SQdel={DQ}",
            str(traj)], ref_out)

    fe = fc.one_file(fc.read(fe_out), fe_out)
    ref = fc.load_legacy(ref_out)
    ref_names = fc.legacy_header(ref_out, "Sq-XRD")

    q = fc.column(fe, fe_out, "q")
    n = fc.same_grid(q, ref[:, 0], "S(q) q 网格")

    return {
        "q": q[:n],
        "fe": {"neutron": fc.column(fe, fe_out, "total_neutron")[:n],
               "xrd": fc.column(fe, fe_out, "total_xrd")[:n]},
        "ref": {"neutron": ref[:n, ref_names.index("Sq-ND")],
                "xrd": ref[:n, ref_names.index("Sq-XRD")]},
        # partial 之间是否相互独立 —— type_new bug 的直接体检项
        "ref_partials": ref[:n, 1:ref_names.index("Sq-XRD")],
        "version": fc.version(fe_out),
    }


def first_peak(q, y, q_lo=1.0, q_hi=4.0):
    """q_lo 之后的第一个**局部**极大（FSDP）。

    不能用区间内的 argmax：中子曲线在 4 Å⁻¹ 附近仍在爬升，argmax 会落在窗口
    右边界上，报出一个并不存在的"峰"。
    """
    m = (q >= q_lo) & (q <= q_hi) & ~np.isnan(y)
    qi, yi = q[m], y[m]
    for k in range(1, len(qi) - 1):
        if yi[k] >= yi[k - 1] and yi[k] > yi[k + 1]:
            return qi[k]
    return float("nan")


def resample(q_src, y_src, q_dst):
    """把实验曲线取到模拟的 q 网格上。

    两者步长同为 0.05 且相位一致，绝大多数点是精确命中；仍走 interp 以容纳
    实验数据缺点的情况。超出实验范围的位置返回 NaN，不外推。
    """
    y = np.interp(q_dst, q_src, y_src, left=np.nan, right=np.nan)
    y[(q_dst < q_src.min()) | (q_dst > q_src.max())] = np.nan
    return y


def plot(data, exp, traj, png, shift, info):
    q = data["q"]
    fig, axes = plt.subplots(3, 2, figsize=(15, 12))
    styles = {
        "exp": dict(color="k", lw=1.6, label="Experiment", zorder=3),
        # ferro 画粗、dump2sq 细虚线叠在其上 —— 两者若重合，蓝色会从红虚线的
        # 缝隙里透出来；若只看到一条实线或一条虚线，那才是真的有差异。
        "fe": dict(color="tab:blue", lw=2.6, alpha=0.85, label="ferro", zorder=2),
        "ref": dict(color="tab:red", lw=1.0, ls="--",
                    label=f"dump2sq (+{shift:g})" if shift else "dump2sq (raw)", zorder=2),
    }

    for col, key in enumerate(["neutron", "xrd"]):
        qe, se = exp[key]
        fe, ref = data["fe"][key], data["ref"][key] + shift
        se_on_q = resample(qe, se, q)
        title = "Neutron" if key == "neutron" else "X-ray"
        xmax = min(q.max(), qe.max())

        for row, (lo, hi, tag) in enumerate([(0.0, xmax, "full range"),
                                             (0.0, 8.0, "low-q zoom")]):
            ax = axes[row][col]
            ax.plot(qe, se, **styles["exp"])
            ax.plot(q, fe, **styles["fe"])
            ax.plot(q, ref, **styles["ref"])
            ax.axhline(1.0, color="grey", lw=0.5, ls=":")
            ax.set_xlim(lo, hi)
            ax.set_ylabel("S(q)")
            ax.set_title(f"{title}  —  {tag}", fontsize=11)
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)
            if row == 1:
                m = (q >= lo) & (q <= hi)
                ys = np.concatenate([fe[m], ref[m], se[(qe >= lo) & (qe <= hi)]])
                pad = 0.08 * (np.nanmax(ys) - np.nanmin(ys))
                ax.set_ylim(np.nanmin(ys) - pad, np.nanmax(ys) + pad)

        ax = axes[2][col]
        for name, y in (("ferro", fe), ("dump2sq", ref)):
            d = y - se_on_q
            c = styles["fe"]["color"] if name == "ferro" else styles["ref"]["color"]
            ax.plot(q, d, lw=0.9, color=c,
                    ls="-" if name == "ferro" else "--",
                    label=f"{name} − exp   rms={np.sqrt(np.nanmean(d ** 2)):.4f}")
        ax.axhline(0.0, color="k", lw=0.6)
        ax.set_xlim(0.0, xmax)
        ax.set_xlabel("q [Å⁻¹]")
        ax.set_ylabel("Δ S(q)")
        ax.set_title(f"{title}  —  simulation − experiment", fontsize=11)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    comp = "  ".join(f"{k}:{v}" for k, v in sorted(info["elems"].items()))
    fig.suptitle(
        f"S(q):  experiment  vs  {data['version']}  vs  dump2sq\n"
        f"{Path(traj).name}   {info['natoms']} atoms ({comp})   "
        f"box {info['box'][0]:.3f} Å   |   r: {R_MIN}–{info['r_max']} step {DR}   "
        f"q: {Q_MIN}–{Q_MAX} step {DQ}",
        fontsize=12)
    fig.text(0.5, 0.005,
             "Raw values, no rescaling — apart from the constant "
             f"{'+%g added to dump2sq (its CalcSq omits the +1 of S(q))' % shift if shift else 'shift (disabled by --raw)'}.",
             ha="center", fontsize=8, style="italic")
    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    fig.savefig(png, dpi=140)
    print(f"\n图已写入 {png}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--traj", default=str(default_traj()))
    p.add_argument("--xlsx", default=str(DEFAULT_XLSX))
    p.add_argument("--outdir", default="./cmp_out/exp")
    p.add_argument("--ferro", default="ferro")
    p.add_argument("--dump2sq", default="dump2sq")
    p.add_argument("--raw", action="store_true",
                   help=f"不给 dump2sq 补 +{SQ_SHIFT_REF:g}，显示两侧原始输出")
    args = p.parse_args()

    traj, xlsx = Path(args.traj).resolve(), Path(args.xlsx).resolve()
    for f in (traj, xlsx):
        if not f.exists():
            sys.exit(f"文件不存在: {f}")
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    shift = 0.0 if args.raw else SQ_SHIFT_REF

    box, natoms, elems, stable = scan_traj(traj)
    # 两侧都必须小于最小盒长的一半：ferro 会 clamp 到最小面间距的一半，
    # dump2sq 则以 GRmax < lmax*0.5 决定走最小镜像分支。留一格余量。
    r_max = round((min(box) * 0.5 - DR) // DR * DR, 6)

    print(f"轨迹   : {traj}")
    print(f"实验   : {xlsx}  [{EXP_SHEET}]")
    print(f"输出   : {outdir}")
    print(f"体系   : {natoms} atoms  " + "  ".join(f"{k}:{v}" for k, v in sorted(elems.items()))
          + f"   box = {box[0]:.4f} × {box[1]:.4f} × {box[2]:.4f} Å")
    print(f"参数   : r {R_MIN}–{r_max} step {DR}   q {Q_MIN}–{Q_MAX} step {DQ}")
    print(f"对齐   : {'无（--raw）' if args.raw else f'dump2sq + {SQ_SHIFT_REF:g}'}\n")

    report_order_stability(stable)
    print()

    data = compute(traj, outdir, r_max, args.ferro, args.dump2sq)
    exp = load_experiment(xlsx)

    # partial 相互独立性：全部雷同即为归类失效的直接证据
    c = partial_independence(data["ref_partials"][data["q"] > 2.0])
    print(f"dump2sq partial 平均非对角相关 (q>2) = {c:+.3f}"
          "    （接近 1 说明归类已失效）\n")

    q = data["q"]
    print("S(q) 数值摘要" + ("（原始值）" if args.raw else f"（dump2sq 已 +{SQ_SHIFT_REF:g}）"))
    print(f"{'':>9}{'rms(fe−exp)':>13}{'rms(ref−exp)':>14}{'rms(fe−ref)':>13}"
          f"{'max|fe−ref|':>13}{'FSDP q_exp':>12}{'q_fe':>8}{'q_ref':>8}")
    for key in ("neutron", "xrd"):
        qe, se = exp[key]
        fe, ref = data["fe"][key], data["ref"][key] + shift
        se_q = resample(qe, se, q)
        d_fe, d_ref, d_fr = fe - se_q, ref - se_q, fe - ref
        print(f"{key:>9}{np.sqrt(np.nanmean(d_fe ** 2)):>13.4f}"
              f"{np.sqrt(np.nanmean(d_ref ** 2)):>14.4f}"
              f"{np.sqrt(np.nanmean(d_fr ** 2)):>13.4f}"
              f"{np.nanmax(np.abs(d_fr)):>13.4f}"
              f"{first_peak(qe, se):>12.2f}"
              f"{first_peak(q, fe):>8.2f}"
              f"{first_peak(q, ref):>8.2f}")

    info = {"box": box, "natoms": natoms, "elems": elems, "r_max": r_max}
    plot(data, exp, traj, outdir / "compare_sq_experiment.png", shift, info)
    plot(data, exp, traj, outdir / "compare_sq_experiment.pdf", shift, info)


if __name__ == "__main__":
    main()
