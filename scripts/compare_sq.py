#!/usr/bin/env python
"""比较 fe-traj 与 dump2sq 计算的结构因子 S(q)（XRD 与中子加权总曲线）。

除下述那个已被实测证实为纯常数的定义差之外，两侧输出**不做任何缩放**，
差异如实呈现，由使用者自行判断。

读图前需要知道的两件事（均为源码事实，不是推断）：

1. 两者对 S(q) 的定义相差一个常数项。dump2sq 的 `CalcSq` 计算的是
   `4πρ/q · Σ_r r[g(r)−1] sin(qr) Δr`，**不含 +1**；fe-traj 的
   `calc_sq_from_gr` 写作 `1.0 + prefactor * integral`。

   完整轨迹实测：残差在 q>15 处稳定为 1.0003（XRD）/ 1.0002（中子），确认这
   就是一个纯常数偏移。**默认把 dump2sq 的曲线整体 +1** 补成标准 S(q)，两条
   曲线才落在同一基准上、残差图才能反映真实差异而不是被常数淹没。选择加到
   dump2sq 而不是从 fe-traj 减，是因为 S(q) 的定义要求大 q 极限趋于 1，
   fe-traj 已是该标准形式。加 `--raw` 可关掉这个补偿、看两侧的原始输出。

2. dump2sq 自带的那份 g(r) 与 dump2analysis / fe-traj 的不一致：它的每一条
   partial g(r) 都叠加了一个位于 P-O 键长（约 1.48 Å）处的伪峰 —— 全部 10 条
   partial 的全域最大值都落在 1.48~1.53 Å，而各自 r>1.8 之后的峰位却完全正常。
   最干净的证据是 O-O：该处本应严格为 0（O-O 最近邻在 2.48 Å），fe-traj 与
   dump2analysis 都给 0，只有 dump2sq 给出约 4.5。

   S(q) 是 g(r) 的傅里叶变换，上游 g(r) 不同则 S(q) 必然不同。因此本图第三行
   并排画出三方的 O-O g(r)，使 S(q) 的差异可以直接归因到上游，而不至于被误读
   成 fe-traj 的变换有问题。

用法:
    python compare_sq.py [--traj FILE] [--outdir DIR]
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
R_MAX = 15.0    # 小于 fe-traj 的最小镜像 clamp 上界（本轨迹 15.16）
DR = 0.01
Q_MIN = 0.1     # fe-traj 固定值，dump2sq 向它对齐
Q_MAX = 25.0
DQ = 0.05

# dump2sq 的 CalcSq 不含 S(q) 的 +1 项，补上后两侧才在同一基准（见模块 docstring）
SQ_SHIFT_REF = 1.0

REPO = Path(__file__).resolve().parent.parent
DEFAULT_TRAJ = REPO / "examples" / "43Z43P15A_NPT.lammpstrj"


def run(cmd, expect):
    """执行外部命令并确认它真的产出了文件（这些 C 程序报错时也返回 0）。"""
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
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        rows.append([float(x) for x in line.split()])
    if not rows:
        sys.exit(f"{path} 中没有数值行")
    return np.array(rows)


def header_columns(path, marker, sep=None):
    """取含 `marker` 的最后一个注释行作为列名行。

    fe-traj 用制表符分隔列名，dump2sq 用空格对齐 —— 用 sep 区分。
    """
    names = None
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") and marker in line:
            names = line.lstrip("# ").split(sep)
    if names is None:
        sys.exit(f"{path} 中找不到含 '{marker}' 的列名行")
    return [n.strip() for n in names if n.strip()]


def fe_version(path):
    return Path(path).read_text().splitlines()[0].lstrip("# ").strip()


def compute(traj, outdir, fe_bin, d2sq_bin, d2a_bin):
    print("[S(q)]")
    # fe-traj：不给配对参数则输出宽表，含 total_xrd / total_neutron
    fe_sq = outdir / "fe.sq.dat"
    run([fe_bin, "-m", "sq", "-i", str(traj), "-o", str(fe_sq),
         "--q-max", str(Q_MAX), "--dq", str(DQ),
         "--r-max", str(R_MAX), "--dr", str(DR),
         "--weighting", "both"], fe_sq)

    # dump2sq：无元素选择参数，一次算出全部 partial 与两条 total
    d2sq_base = outdir / "d2sq"
    run([d2sq_bin, "-o", str(d2sq_base),
         f"--GRmin={R_MIN}", f"--GRmax={R_MAX}", f"--GRdel={DR}",
         f"--SQmin={Q_MIN}", f"--SQmax={Q_MAX}", f"--SQdel={DQ}",
         str(traj)], Path(str(d2sq_base) + ".sq"))

    print("[O-O g(r) 归因]")
    fe_gr = outdir / "fe_O-O.gr.dat"
    run([fe_bin, "-m", "gr", "-a", "O", "-b", "O", "-i", str(traj), "-o", str(fe_gr),
         "--r-max", str(R_MAX), "--dr", str(DR)], fe_gr)
    d2a_gr = outdir / "d2a_O-O.gr"
    run([d2a_bin, "-m", "gr", "-x", "O", "-y", "O",
         "--rmin", str(R_MIN), "--rmax", str(R_MAX), "--dr", str(DR),
         "-i", str(traj), "-o", str(d2a_gr)], d2a_gr)

    # ── S(q) ────────────────────────────────────────────────────────────────
    fe = load_columns(fe_sq)
    fe_names = header_columns(fe_sq, "q[Ang", sep="\t")
    ref = load_columns(Path(str(d2sq_base) + ".sq"))
    ref_names = header_columns(Path(str(d2sq_base) + ".sq"), "Sq-XRD")

    n = min(len(fe), len(ref))
    dq_max = np.abs(fe[:n, 0] - ref[:n, 0]).max()
    if dq_max > 1e-9:
        sys.exit(f"q 网格不一致，最大偏差 {dq_max:.3e} —— 参数未对齐")

    sq = {
        "q": fe[:n, 0],
        "XRD": (fe[:n, fe_names.index("total_xrd")], ref[:n, ref_names.index("Sq-XRD")]),
        "Neutron": (fe[:n, fe_names.index("total_neutron")], ref[:n, ref_names.index("Sq-ND")]),
        "version": fe_version(fe_sq),
        "n_fe": len(fe),
        "n_ref": len(ref),
    }

    # ── O-O g(r) ────────────────────────────────────────────────────────────
    g_fe = load_columns(fe_gr)
    g_fe_names = header_columns(fe_gr, "r[Ang", sep="\t")
    g_d2a = load_columns(d2a_gr)
    g_d2a = g_d2a[g_d2a[:, 0] > 0]          # 去掉 dump2analysis 固定输出的 r=0 行
    g_ref = load_columns(Path(str(d2sq_base) + ".gr"))
    g_ref_names = header_columns(Path(str(d2sq_base) + ".gr"), "Total")

    m = min(len(g_fe), len(g_d2a), len(g_ref))
    gr = {
        "r": g_fe[:m, 0],
        "fe": g_fe[:m, g_fe_names.index("O-O_gr")],
        "d2a": g_d2a[:m, 1],
        "d2sq": g_ref[:m, g_ref_names.index("O-O")],
    }
    return sq, gr


def plot(sq, gr, traj, png, shift):
    fig, axes = plt.subplots(3, 2, figsize=(13, 12))
    q = sq["q"]
    ref_label = "dump2sq{}".format(f" + {shift:g}" if shift else "")

    for col, key in enumerate(["XRD", "Neutron"]):
        fe, ref = sq[key]
        ref = ref + shift
        diff = fe - ref

        ax = axes[0][col]
        ax.plot(q, fe, lw=1.0, label=f"fe-traj  total_{key.lower()}")
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
        axd.set_title(f"fe-traj − {ref_label}    mean over q>15 = {tail.mean():.5f}   "
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
        ax.plot(r[sel], gr["fe"][sel], lw=1.2, label="fe-traj")
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
                f"at r=1.48 Å\nfe-traj {gr['fe'][k]:.3f} / d2a {gr['d2a'][k]:.3f} / "
                f"dump2sq {gr['d2sq'][k]:.3f}",
                xy=(1.48, gr["d2sq"][k]), xytext=(1.9, gr["d2sq"][k] * 0.9),
                fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.8))

    shift_note = (f"dump2sq shifted by +{shift:g} (its CalcSq omits the +1 of S(q))"
                  if shift else "raw values, no shift")
    fig.suptitle(
        f"S(q):  fe-traj ({sq['version']})  vs  dump2sq\n"
        f"{Path(traj).name}   |   r: {R_MIN}–{R_MAX} step {DR}   "
        f"q: {Q_MIN}–{Q_MAX} step {DQ}   |   {shift_note}",
        fontsize=12)
    fig.text(0.5, 0.005,
             "Only the constant offset is compensated. The remaining low-q deviation comes from "
             "dump2sq's own g(r), which carries a spurious peak at the P-O bond length (bottom row).",
             ha="center", fontsize=8, style="italic")
    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    fig.savefig(png, dpi=140)
    print(f"\n图已写入 {png}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--traj", default=str(DEFAULT_TRAJ))
    p.add_argument("--outdir", default="./cmp_out")
    p.add_argument("--fe-traj", default="fe-traj")
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

    sq, gr = compute(traj, outdir, args.fe_traj, args.dump2sq, args.dump2analysis)

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
    print(f"  fe-traj       {gr['fe'][k]:.4f}")
    print(f"  dump2analysis {gr['d2a'][k]:.4f}")
    print(f"  dump2sq       {gr['d2sq'][k]:.4f}")

    plot(sq, gr, traj, outdir / "compare_sq.png", shift)


if __name__ == "__main__":
    main()
