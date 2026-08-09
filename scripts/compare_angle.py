#!/usr/bin/env python
"""比较 fe-traj 与 dump2analysis 计算的键角分布 P(θ)。

对比 O-P-O 与 O-Al-O 两个三元组。两侧输出的**原始计数一律不做归一化**，
差异如实呈现，由使用者自行判断。

两处必须知道的约定差异（都不修改数据，只影响读图）：

1. 计数相差约 2 倍。dump2analysis 对 A-B-C 枚举时把 O₁-P-O₂ 与 O₂-P-O₁ 各数一次，
   fe-traj 用 canonical triplet 只数一次。两者的 mean/std 完全一致，说明只是
   计数约定不同、分布形状相同。这里**不做归一化**，倍率直接显示在叠加图上。

2. 横轴标注约定不同：dump2analysis 标 bin 的**右边界**（0.1, 0.2 … 180.0），
   fe-traj 标 bin 的**中心**（0.05, 0.15 … 179.95）。同一索引的 bin 覆盖同一个
   角度区间。叠加图按各自的原始横轴绘制（0.05° 偏移如实可见）；残差则按
   **bin 索引**对齐 —— 否则在半高宽仅数度、bin 宽 0.1° 的峰上，相邻 bin 相减会
   产生纯属标注错位的巨大伪残差，把真实差异完全淹没。

dump2analysis 的角度 bin 宽度固定 0.1°（无对应参数），故 fe-traj 显式传
--d-angle 0.1 与之对齐。

用法:
    python compare_angle.py [--traj FILE] [--outdir DIR]
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
R_CUT_AB = 2.3
R_CUT_BC = 2.3
D_ANGLE = 0.1   # dump2analysis 固定 0.1°，fe-traj 显式对齐

TRIPLETS = [("O", "P", "O"), ("O", "Al", "O")]

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
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        rows.append([float(x) for x in line.split()])
    if not rows:
        sys.exit(f"{path} 中没有数值行")
    return np.array(rows)


def fe_version(path):
    return Path(path).read_text().splitlines()[0].lstrip("# ").strip()


def fe_stats(path):
    """抽出 fe-traj 头部的 `# O-P-O : mean=... std=... count=...` 一行。"""
    for line in Path(path).read_text().splitlines():
        if not line.startswith("#") or "mean=" not in line:
            continue
        body = line.lstrip("# ")
        try:
            mean = float(body.split("mean=")[1].split()[0])
            std = float(body.split("std=")[1].split()[0])
            count = float(body.split("count=")[1].split()[0])
            return mean, std, count
        except (IndexError, ValueError):
            continue
    return None


def d2a_stats(path):
    """抽出 dump2analysis 头部的 `# angle: mean +/- std` 与 `# Total pairs:`。"""
    mean = std = count = None
    for line in Path(path).read_text().splitlines():
        if not line.startswith("#"):
            break
        if line.startswith("# angle:"):
            parts = line.split(":")[1].split("+/-")
            mean, std = float(parts[0]), float(parts[1])
        elif "Total pairs:" in line:
            count = float(line.split(":")[1])
    return mean, std, count


def compute(traj, outdir, fe_bin, d2a_bin):
    results = {}
    for a, b, c in TRIPLETS:
        tag = f"{a}-{b}-{c}"
        print(f"[{tag}]")

        # fe-traj：-a/-b/-c 是【元素】，-b 为中心原子
        fe_out = outdir / f"fe_{tag}.angle.dat"
        run([fe_bin, "-m", "angle", "-a", a, "-b", b, "-c", c,
             "-i", str(traj), "-o", str(fe_out),
             "--r-cut-ab", str(R_CUT_AB), "--r-cut-bc", str(R_CUT_BC),
             "--d-angle", str(D_ANGLE)], fe_out)

        # dump2analysis：-x/-y/-z 才是【元素】，-a/-b/-c 是原子 id
        d2a_out = outdir / f"d2a_{tag}.angle"
        run([d2a_bin, "-m", "angle", "-x", a, "-y", b, "-z", c,
             "--rcut_ab", str(R_CUT_AB), "--rcut_bc", str(R_CUT_BC),
             "-i", str(traj), "-o", str(d2a_out)], d2a_out)

        fe = load_columns(fe_out)
        d2a = load_columns(d2a_out)

        n = min(len(fe), len(d2a))
        # 按 bin 索引对齐（见模块 docstring 第 2 条）；先确认两侧确实是同宽度的 bin
        step_fe = float(np.median(np.diff(fe[:, 0])))
        step_ref = float(np.median(np.diff(d2a[:, 0])))
        if abs(step_fe - step_ref) > 1e-9 or abs(step_fe - D_ANGLE) > 1e-9:
            sys.exit(f"{tag}: bin 宽度不一致 fe={step_fe} ref={step_ref}")

        results[tag] = {
            "x_fe": fe[:n, 0],
            "x_ref": d2a[:n, 0],
            "y_fe": fe[:n, 1],
            "y_ref": d2a[:n, 1],
            "fe_stats": fe_stats(fe_out),
            "ref_stats": d2a_stats(d2a_out),
            "version": fe_version(fe_out),
        }
    return results


def plot(results, traj, png):
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    ver = next(iter(results.values()))["version"]

    for col, (a, b, c) in enumerate(TRIPLETS):
        tag = f"{a}-{b}-{c}"
        d = results[tag]
        fe, ref = d["y_fe"], d["y_ref"]
        ratio = fe.sum() / ref.sum() if ref.sum() else float("nan")

        # 叠加：各用自己的原始横轴，原始计数
        ax = axes[0][col]
        ax.plot(d["x_fe"], fe, lw=1.0, label="fe-traj (bin centre)")
        ax.plot(d["x_ref"], ref, lw=1.0, ls="--", label="dump2analysis (bin right edge)")
        ax.set_title(f"{tag}   raw counts   Σfe/Σref = {ratio:.4f}", fontsize=11)
        ax.set_ylabel("count")
        ax.set_xlim(60, 180)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

        fs, rs = d["fe_stats"], d["ref_stats"]
        if fs and rs[0] is not None:
            ax.text(0.02, 0.97,
                    f"fe : mean={fs[0]:.4f}  std={fs[1]:.4f}  N={fs[2]:.0f}\n"
                    f"ref: mean={rs[0]:.4f}  std={rs[1]:.4f}  N={rs[2]:.0f}",
                    transform=ax.transAxes, va="top", fontsize=7, family="monospace")

        # 残差：按 bin 索引对齐，横轴借用 fe 的中心值作刻度
        axd = axes[1][col]
        diff = fe - ref
        axd.plot(d["x_fe"], diff, lw=0.9, color="crimson")
        axd.axhline(0.0, color="k", lw=0.6)
        axd.set_title(f"fe − ref (aligned by bin index)   "
                      f"max|Δ|={np.abs(diff).max():.4g}   "
                      f"Σ|Δ|/Σref={np.abs(diff).sum() / ref.sum():.4f}", fontsize=9)
        axd.set_ylabel("Δ count")
        axd.set_xlabel("angle [deg]")
        axd.set_xlim(60, 180)
        axd.grid(alpha=0.3)

    fig.suptitle(
        f"Bond angle distribution:  fe-traj ({ver})  vs  dump2analysis\n"
        f"{Path(traj).name}   |   r_cut_ab={R_CUT_AB}, r_cut_bc={R_CUT_BC}, "
        f"d_angle={D_ANGLE}   |   raw counts, NOT normalised",
        fontsize=12)
    fig.text(0.5, 0.005,
             "dump2analysis counts each A-B-C twice (O1-P-O2 and O2-P-O1); fe-traj counts it once. "
             "Counts are left unscaled on purpose.",
             ha="center", fontsize=8, style="italic")
    fig.tight_layout(rect=(0, 0.03, 1, 0.93))
    fig.savefig(png, dpi=140)
    print(f"\n图已写入 {png}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--traj", default=str(DEFAULT_TRAJ))
    p.add_argument("--outdir", default="./cmp_out")
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
    print(f"参数   : r_cut_ab={R_CUT_AB}  r_cut_bc={R_CUT_BC}  d_angle={D_ANGLE}\n")

    results = compute(traj, outdir, args.fe_traj, args.dump2analysis)

    print("\n数值摘要（原始计数，未归一化）")
    print(f"{'triplet':>10}{'Σ fe':>12}{'Σ ref':>12}{'Σfe/Σref':>11}"
          f"{'mean fe':>11}{'mean ref':>11}{'std fe':>10}{'std ref':>10}")
    for tag, d in results.items():
        fs, rs = d["fe_stats"], d["ref_stats"]
        print(f"{tag:>10}{d['y_fe'].sum():>12.0f}{d['y_ref'].sum():>12.0f}"
              f"{d['y_fe'].sum() / d['y_ref'].sum():>11.4f}"
              f"{fs[0]:>11.4f}{rs[0]:>11.4f}{fs[1]:>10.4f}{rs[1]:>10.4f}")

    plot(results, traj, outdir / "compare_angle.png")


if __name__ == "__main__":
    main()
