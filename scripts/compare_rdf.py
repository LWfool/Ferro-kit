#!/usr/bin/env python
"""比较 `ferro traj gr` 与 dump2analysis 计算的 g(r) 与 CN(r)。

对比 P-O 与 Al-O 两个对。两侧输出的**原始数值一律不做任何缩放或平移**，
差异如实呈现，由使用者自行判断。

参数选择上的一个不对称（无法回避，只能靠调另一侧对齐）：
  * dump2analysis 的 --rmax 默认是 20.005（其 -h 里写的 10.005 已过时），
    且它不做最小镜像截断；ferro 会把 r_max clamp 到最小面间距的一半。
    取 15.0 是因为它同时小于 ferro 的 clamp 上界（本轨迹 15.16）。

`--r-min` 自 0.1.15 起已上 CLI（默认 0.001），本脚本显式传 0.005 与 dump2analysis
对齐 —— ferro 的 bin 中心是 `r_min + (i+0.5)·dr`，两侧同 r_min 同 dr 即逐点重合。

另一处需要留意：dump2analysis 的 --rcut 只用于其输出头部那行
`# distance: mean +/- std`，不进入 g(r)/CN(r)。ferro 在 0.1.11 已删除对应的
PairStats，故此参数无对应物，这里显式传 2.3 但不使用其结果。

**0.2.0 的产物形状**（详见 ferrocmp.py 的模块 docstring）：
  * `-o` 是文件名后缀不是路径 → 以 outdir 为工作目录调用
  * gr 是长表 `file,r,center,neighbor,gr,cn` → 按 center/neighbor 选行，不按列名取列
  * 产物名带 label 段（2026-08-13 起）：`gr_P-O_cmp.csv` → 用 `fc.product_name()` 拼

用法:
    python compare_rdf.py [--traj FILE] [--outdir DIR]
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import ferrocmp as fc

# ── 统一计算参数 ──────────────────────────────────────────────────────────────
R_MIN = 0.005   # 两侧显式对齐（ferro 默认 0.001，dump2analysis 默认 0.0）
R_MAX = 15.0
DR = 0.01
R_CUT = 2.3     # 仅 dump2analysis 的键长统计用，不影响曲线

PAIRS = [("P", "O"), ("Al", "O")]
SUFFIX = "cmp"  # -o 的批次标记；配对由文件名的 label 段区分，不必再进后缀

REPO = Path(__file__).resolve().parent.parent
DEFAULT_TRAJ = REPO / "examples" / "43Z43P15A_NPT.lammpstrj"


def compute(traj, outdir, ferro_bin, d2a_bin):
    """跑两侧程序，返回 {pair: dict} 。"""
    results = {}
    for a, b in PAIRS:
        tag = f"{a}-{b}"
        print(f"[{tag}]")

        # ferro：-a 是中心元素、-b 是近邻元素；-o 给的是后缀，产物落在 outdir 下。
        # 配对已经由 label 段进了文件名，故 -o 只留一个批次标记，不再重复 tag
        fe_out = fc.run_ferro(
            ferro_bin,
            ["traj", "gr", "-a", a, "-b", b, "-i", traj, "-o", SUFFIX,
             "--r-min", R_MIN, "--r-max", R_MAX, "--dr", DR],
            outdir,
            fc.product_name("gr", label=fc.file_label(a, b), suffix=SUFFIX))

        # dump2analysis：-x/-y 才是【元素】，-a/-b 是原子 id
        d2a_out = outdir / f"d2a_{tag}.gr"
        fc.run([d2a_bin, "-m", "gr", "-x", a, "-y", b,
                "--rmin", str(R_MIN), "--rmax", str(R_MAX), "--dr", str(DR),
                "--rcut", str(R_CUT),
                "-i", str(traj), "-o", str(d2a_out)], d2a_out)

        # 长表：按 center/neighbor 选出这一对。给了 -a/-b 时表里本就只有这一对，
        # 仍显式选一次 —— 元素名拼错时 pick() 会把表里实际有的组合列出来
        fe = fc.pick(fc.one_file(fc.read(fe_out), fe_out), fe_out, "r",
                     center=a, neighbor=b)
        d2a = fc.load_legacy(d2a_out)

        # dump2analysis 额外写了一行 r=0（gr.c 的 OutputGr 固定输出 "0 0 0"），
        # 去掉后两侧的 bin 中心逐点一致
        d2a = d2a[d2a[:, 0] > 0]

        r_fe = fe["r"].to_numpy()
        n = fc.same_grid(r_fe, d2a[:, 0], tag)

        results[tag] = {
            "r": r_fe[:n],
            "fe_gr": fe["gr"].to_numpy()[:n],
            "fe_cn": fe["cn"].to_numpy()[:n],
            "d2a_gr": d2a[:n, 1],
            "d2a_cn": d2a[:n, 2],
            "version": fc.version(fe_out),
            "n_fe": len(fe),
            "n_d2a": len(d2a),
        }
    return results


def plot(results, traj, png):
    fig, axes = plt.subplots(4, 2, figsize=(13, 14), sharex="col")
    ver = next(iter(results.values()))["version"]

    for col, (a, b) in enumerate(PAIRS):
        tag = f"{a}-{b}"
        d = results[tag]
        r = d["r"]

        for row, (key, label) in enumerate([("gr", "g(r)"), ("cn", "CN(r)")]):
            fe, ref = d[f"fe_{key}"], d[f"d2a_{key}"]
            diff = fe - ref

            ax = axes[2 * row][col]
            ax.plot(r, fe, lw=1.2, label=f"ferro  {tag}")
            ax.plot(r, ref, lw=1.0, ls="--", label=f"dump2analysis  {tag}")
            ax.set_ylabel(label)
            ax.set_title(f"{tag}   {label}", fontsize=11)
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)

            # CN(r) 到 15 Å 已累积到数百，第一配位壳（~4）在全域纵轴上被压成一条线。
            # 加一个放大视窗；数据未做任何改动，只是换个视角看同一条曲线。
            if key == "cn":
                sel = (r >= 1.0) & (r <= 3.5)
                axin = ax.inset_axes([0.12, 0.45, 0.45, 0.5])
                axin.plot(r[sel], fe[sel], lw=1.2)
                axin.plot(r[sel], ref[sel], lw=1.0, ls="--")
                axin.set_title("first shell (zoom)", fontsize=7)
                axin.tick_params(labelsize=6)
                axin.grid(alpha=0.3)

            axd = axes[2 * row + 1][col]
            axd.plot(r, diff, lw=0.9, color="crimson")
            axd.axhline(0.0, color="k", lw=0.6)
            axd.set_ylabel(f"Δ {label}")
            mx = np.abs(diff).max()
            rms = float(np.sqrt(np.mean(diff ** 2)))
            axd.set_title(f"ferro − dump2analysis    max|Δ|={mx:.3e}   rms={rms:.3e}",
                          fontsize=9)
            axd.grid(alpha=0.3)

    for ax in axes[-1]:
        ax.set_xlabel("r [Å]")

    fig.suptitle(
        f"g(r) / CN(r):  {ver}  vs  dump2analysis\n"
        f"{Path(traj).name}   |   r_min={R_MIN}, r_max={R_MAX}, dr={DR}   |   "
        f"raw values, no rescaling",
        fontsize=12)
    # 残差呈现为离散台阶，是 dump2analysis 以 %12g（6 位有效数字）写出的量化误差 ——
    # 台阶高度随曲线自身量级增大，正是有效数字固定的特征。
    fig.text(0.5, 0.005,
             "Residual steps are the quantisation of dump2analysis' %12g output "
             "(6 significant digits), not a difference in the computation.",
             ha="center", fontsize=8, style="italic")
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    fig.savefig(png, dpi=140)
    print(f"\n图已写入 {png}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--traj", default=str(DEFAULT_TRAJ), help="LAMMPS dump 轨迹")
    p.add_argument("--outdir", default="./cmp_out", help="输出目录（保留全部中间文件）")
    p.add_argument("--ferro", default="ferro")
    p.add_argument("--dump2analysis", default="dump2analysis")
    args = p.parse_args()

    traj = Path(args.traj).resolve()
    if not traj.exists():
        sys.exit(f"轨迹不存在: {traj}")
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"轨迹   : {traj}")
    print(f"输出   : {outdir}")
    print(f"参数   : r_min={R_MIN}  r_max={R_MAX}  dr={DR}\n")

    results = compute(traj, outdir, args.ferro, args.dump2analysis)

    print("\n数值摘要（原始值，未做任何修正）")
    print(f"{'pair':>6}{'量':>6}{'fe 峰值':>12}{'ref 峰值':>12}{'峰位差[Å]':>12}"
          f"{'max|Δ|':>12}{'rms Δ':>12}")
    for tag, d in results.items():
        for key, label in [("gr", "g(r)"), ("cn", "CN(r)")]:
            fe, ref, r = d[f"fe_{key}"], d[f"d2a_{key}"], d["r"]
            diff = fe - ref
            print(f"{tag:>6}{label:>6}{fe.max():>12.4f}{ref.max():>12.4f}"
                  f"{r[np.argmax(fe)] - r[np.argmax(ref)]:>12.4f}"
                  f"{np.abs(diff).max():>12.3e}{np.sqrt(np.mean(diff**2)):>12.3e}")

    plot(results, traj, outdir / "compare_rdf.png")


if __name__ == "__main__":
    main()
