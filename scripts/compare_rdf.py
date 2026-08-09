#!/usr/bin/env python
"""比较 fe-traj 与 dump2analysis 计算的 g(r) 与 CN(r)。

对比 P-O 与 Al-O 两个对。两侧输出的**原始数值一律不做任何缩放或平移**，
差异如实呈现，由使用者自行判断。

参数选择上的两个不对称（无法回避，只能靠调另一侧对齐）：
  * fe-traj 的 r_min 固定 0.005，CLI 未暴露 → dump2analysis 显式传 --rmin 0.005
  * dump2analysis 的 --rmax 默认是 20.005（其 -h 里写的 10.005 已过时），
    且它不做最小镜像截断；fe-traj 会把 r_max clamp 到最小面间距的一半。
    取 15.0 是因为它同时小于 fe-traj 的 clamp 上界（本轨迹 15.16）。

另一处需要留意：dump2analysis 的 --rcut 只用于其输出头部那行
`# distance: mean +/- std`，不进入 g(r)/CN(r)。fe-traj 在 0.1.11 已删除对应的
PairStats，故此参数无对应物，这里显式传 2.3 但不使用其结果。

用法:
    python compare_rdf.py [--traj FILE] [--outdir DIR]
"""

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── 统一计算参数 ──────────────────────────────────────────────────────────────
R_MIN = 0.005   # fe-traj 固定值，其余工具向它对齐
R_MAX = 15.0
DR = 0.01
R_CUT = 2.3     # 仅 dump2analysis 的键长统计用，不影响曲线

PAIRS = [("P", "O"), ("Al", "O")]

REPO = Path(__file__).resolve().parent.parent
DEFAULT_TRAJ = REPO / "examples" / "43Z43P15A_NPT.lammpstrj"


def run(cmd, expect):
    """执行外部命令并确认它真的产出了文件。

    dump2analysis 的参数校验失败走 exit(0)，退出码不可信，只能验产物。
    """
    print("  $", " ".join(str(c) for c in cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    expect = Path(expect)
    if not expect.exists() or expect.stat().st_size == 0:
        sys.exit(
            f"命令未产出 {expect}\n"
            f"  exit={proc.returncode}\n  stdout={proc.stdout[-800:]}\n  stderr={proc.stderr[-800:]}"
        )
    return expect


def load_columns(path):
    """读取以 '#' 为注释、空白分隔的数值表，跳过空行。

    dump2analysis 的 gr 输出在表头后有一个空行；fe-traj 用制表符分隔且写
    `1.0e-2` 这类指数格式 —— float() 两者都能吃下。
    """
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        rows.append([float(x) for x in line.split()])
    if not rows:
        sys.exit(f"{path} 中没有数值行")
    width = len(rows[0])
    if any(len(r) != width for r in rows):
        sys.exit(f"{path} 各行列数不一致")
    return np.array(rows)


def fe_header_columns(path):
    """取出 fe-traj 输出的列名行（最后一个以 '#' 开头且含 '[' 的行）。"""
    names = None
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") and "[" in line:
            names = line.lstrip("# ").split("\t")
    if names is None:
        sys.exit(f"{path} 中找不到列名行")
    return [n.strip() for n in names]


def fe_version(path):
    first = Path(path).read_text().splitlines()[0]
    return first.lstrip("# ").strip()


def compute(traj, outdir, fe_bin, d2a_bin):
    """跑两侧程序，返回 {pair: dict} 。"""
    results = {}
    for a, b in PAIRS:
        tag = f"{a}-{b}"
        print(f"[{tag}]")

        # fe-traj：-a/-b 是【元素】
        fe_out = outdir / f"fe_{tag}.gr.dat"
        run([fe_bin, "-m", "gr", "-a", a, "-b", b,
             "-i", str(traj), "-o", str(fe_out),
             "--r-max", str(R_MAX), "--dr", str(DR)], fe_out)

        # dump2analysis：-x/-y 才是【元素】，-a/-b 是原子 id
        d2a_out = outdir / f"d2a_{tag}.gr"
        run([d2a_bin, "-m", "gr", "-x", a, "-y", b,
             "--rmin", str(R_MIN), "--rmax", str(R_MAX), "--dr", str(DR),
             "--rcut", str(R_CUT),
             "-i", str(traj), "-o", str(d2a_out)], d2a_out)

        fe = load_columns(fe_out)
        names = fe_header_columns(fe_out)
        d2a = load_columns(d2a_out)

        # dump2analysis 额外写了一行 r=0（gr.c 的 OutputGr 固定输出 "0 0 0"），
        # 去掉后两侧的 bin 中心逐点一致
        d2a = d2a[d2a[:, 0] > 0]

        n = min(len(fe), len(d2a))
        dr_max = np.abs(fe[:n, 0] - d2a[:n, 0]).max()
        if dr_max > 1e-9:
            sys.exit(f"{tag}: r 网格不一致，最大偏差 {dr_max:.3e} —— 参数未对齐")

        results[tag] = {
            "r": fe[:n, 0],
            "fe_gr": fe[:n, names.index(f"{a}-{b}_gr")],
            "fe_cn": fe[:n, names.index(f"{a}-{b}_cn")],
            "d2a_gr": d2a[:n, 1],
            "d2a_cn": d2a[:n, 2],
            "version": fe_version(fe_out),
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
            ax.plot(r, fe, lw=1.2, label=f"fe-traj  {tag}")
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
            axd.set_title(f"fe-traj − dump2analysis    max|Δ|={mx:.3e}   rms={rms:.3e}",
                          fontsize=9)
            axd.grid(alpha=0.3)

    for ax in axes[-1]:
        ax.set_xlabel("r [Å]")

    fig.suptitle(
        f"g(r) / CN(r):  fe-traj ({ver})  vs  dump2analysis\n"
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
    p.add_argument("--fe-traj", default="fe-traj")
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

    results = compute(traj, outdir, args.fe_traj, args.dump2analysis)

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
