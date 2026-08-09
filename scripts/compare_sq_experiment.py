#!/usr/bin/env python
"""把 fe-traj 与 dump2sq 算出的 S(q) 和实验数据画在一起。

轨迹默认用 70ZnO-30P2O5（无 Al）的 NVT 轨迹：优先 tests/70Z30P00A_NVT.lammpstrj
（50 帧，8 MB，未入库），缺席时退到 tests/70Z30P00A_NVT_5.lammpstrj（5 帧，已入库）。
实验数据取 tests/Experiment.xlsx 的 Sheet1 —— 该表的四列
nq_70 / SQn_70 / xq_70 / SQx_70 正是同一组分的中子与 X 射线总结构因子。

除下面这一处已被实测证实为纯常数的定义差之外，三方曲线**不做任何缩放或平移**，
差异如实呈现。

**dump2sq 的 +1**：其 `CalcSq` 算的是 `4πρ/q · Σ_r r[g(r)−1] sin(qr) Δr`，不含
S(q) 定义里的 +1；fe-traj 与实验数据都是 q→∞ 趋于 1 的标准形式。默认给 dump2sq
整条曲线 +1 补到同一基准，`--raw` 可关掉。

**这条轨迹为什么值得单独跑**：dump2sq 的 `InitializeType` 只在第 0 帧写入
`data[i].type_new`，而 `ReadDataLammpstrj` 之后从不更新它，于是原子类型按**数组
下标**而非实际元素归类。dump 若未按 id 排序，第 1 帧起 partial 归类即失效。
examples/43Z43P15A_NPT.lammpstrj 正是这种未排序轨迹（帧间同位置同 id 仅 210/2004），
其 10 条 partial 相互相关系数高达 0.959、XRD 与中子总曲线几乎重合 —— 那份对比里
S(q) 峰位对不上，根源在此，与傅里叶变换本身无关。

本轨迹按 id 排序且帧间顺序完全稳定（5003/5003），bug 不触发，因此这是一次干净的
三方对照。脚本会自动检查排序稳定性并把结果打印出来。

用法:
    python compare_sq_experiment.py [--traj FILE] [--outdir DIR] [--raw]
"""

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import openpyxl

# ── 统一计算参数 ──────────────────────────────────────────────────────────────
R_MIN = 0.005   # fe-traj 固定值，dump2sq 向它对齐
DR = 0.01
Q_MIN = 0.1     # fe-traj 固定值（CLI 未暴露），dump2sq 向它对齐
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
    REPO / "tests" / "70Z30P00A_NVT.lammpstrj",
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


def run(cmd, expect):
    """执行外部命令并确认它真的产出了文件（这些 C 程序报错时也可能返回 0）。"""
    print("  $", " ".join(str(c) for c in cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    expect = Path(expect)
    if not expect.exists() or expect.stat().st_size == 0:
        sys.exit(
            f"命令未产出 {expect}\n"
            f"  exit={proc.returncode}\n  stdout={proc.stdout[-800:]}\n  stderr={proc.stderr[-800:]}"
        )
    return expect


def load_table(path):
    rows = [[float(x) for x in ln.split()]
            for ln in Path(path).read_text().splitlines()
            if ln.strip() and not ln.startswith("#")]
    if not rows:
        sys.exit(f"{path} 中没有数值行")
    return np.array(rows)


def fe_columns(path):
    """fe-traj 的列名行：最后一个以 '#' 开头且含 '[' 的行，制表符分隔。"""
    names = None
    for ln in Path(path).read_text().splitlines():
        if ln.startswith("#") and "[" in ln:
            names = [n.strip() for n in ln.lstrip("# ").split("\t")]
    if names is None:
        sys.exit(f"{path} 中找不到列名行")
    return names


def d2sq_columns(path):
    """dump2sq 的列名行是最后一个 '#' 行，定宽 %15s，末两列为 Sq-XRD / Sq-ND。"""
    last = None
    for ln in Path(path).read_text().splitlines():
        if ln.startswith("#"):
            last = ln
        else:
            break
    return last.lstrip("# ").split()


def fe_version(path):
    return Path(path).read_text().splitlines()[0].lstrip("# ").strip()


# ── 轨迹自检 ──────────────────────────────────────────────────────────────────

def scan_traj(traj, n_probe=3):
    """读前几帧，返回 (盒子边长, 原子数, 元素计数, 帧间同位置同 id 的比例)。

    最后一项直接决定 dump2sq 的 partial 是否可信（见模块 docstring）。
    """
    order, box, elems = [], None, None
    with open(traj) as f:
        for _ in range(n_probe):
            for _ in range(3):
                if not f.readline():
                    break
            line = f.readline()
            if not line:
                break
            n = int(line)
            f.readline()                                   # ITEM: BOX BOUNDS
            b = [[float(v) for v in f.readline().split()] for _ in range(3)]
            f.readline()                                   # ITEM: ATOMS
            rows = [f.readline().split() for _ in range(n)]
            order.append([r[0] for r in rows])
            if box is None:
                box = [hi - lo for lo, hi, *_ in b]
                elems = {}
                for r in rows:
                    elems[r[2]] = elems.get(r[2], 0) + 1
    stable = min(
        sum(a == b for a, b in zip(order[0], o)) / len(order[0])
        for o in order[1:]
    ) if len(order) > 1 else float("nan")
    return box, len(order[0]), elems, stable


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


def compute(traj, outdir, r_max, fe_bin, d2sq_bin):
    fe_out = outdir / "fe_sq.dat"
    run([fe_bin, "-m", "sq", "-i", str(traj), "-o", str(fe_out),
         "--q-max", str(Q_MAX), "--dq", str(DQ),
         "--r-max", str(r_max), "--dr", str(DR),
         "--weighting", "both"], fe_out)

    # dump2sq 把 argv 的最后一项当输入文件；输出为 <outbase>.{gr,cn,sq,coeff}
    base = outdir / "d2sq"
    run([d2sq_bin, "-o", str(base),
         f"--GRmin={R_MIN}", f"--GRmax={r_max}", f"--GRdel={DR}",
         f"--SQmin={Q_MIN}", f"--SQmax={Q_MAX}", f"--SQdel={DQ}",
         str(traj)], base.with_suffix(".sq"))

    fe = load_table(fe_out)
    fe_names = fe_columns(fe_out)
    ref = load_table(base.with_suffix(".sq"))
    ref_names = d2sq_columns(base.with_suffix(".sq"))

    dq = np.abs(fe[:, 0] - ref[:, 0]).max()
    if dq > 1e-9:
        sys.exit(f"两侧 q 网格不一致，最大偏差 {dq:.3e} —— 参数未对齐")

    return {
        "q": fe[:, 0],
        "fe": {"neutron": fe[:, fe_names.index("total_neutron")],
               "xrd": fe[:, fe_names.index("total_xrd")]},
        "ref": {"neutron": ref[:, ref_names.index("Sq-ND")],
                "xrd": ref[:, ref_names.index("Sq-XRD")]},
        # partial 之间是否相互独立 —— type_new bug 的直接体检项
        "ref_partials": ref[:, 1:ref_names.index("Sq-XRD")],
        "version": fe_version(fe_out),
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
        # fe-traj 画粗、dump2sq 细虚线叠在其上 —— 两者若重合，蓝色会从红虚线的
        # 缝隙里透出来；若只看到一条实线或一条虚线，那才是真的有差异。
        "fe": dict(color="tab:blue", lw=2.6, alpha=0.85, label="fe-traj", zorder=2),
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
        for name, y in (("fe-traj", fe), ("dump2sq", ref)):
            d = y - se_on_q
            c = styles["fe"]["color"] if name == "fe-traj" else styles["ref"]["color"]
            ax.plot(q, d, lw=0.9, color=c,
                    ls="-" if name == "fe-traj" else "--",
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
        f"S(q):  experiment  vs  fe-traj ({data['version']})  vs  dump2sq\n"
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
    p.add_argument("--fe-traj", default="fe-traj")
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
    # 两侧都必须小于最小盒长的一半：fe-traj 会 clamp 到最小面间距的一半，
    # dump2sq 则以 GRmax < lmax*0.5 决定走最小镜像分支。留一格余量。
    r_max = round((min(box) * 0.5 - DR) // DR * DR, 6)

    print(f"轨迹   : {traj}")
    print(f"实验   : {xlsx}  [{EXP_SHEET}]")
    print(f"输出   : {outdir}")
    print(f"体系   : {natoms} atoms  " + "  ".join(f"{k}:{v}" for k, v in sorted(elems.items()))
          + f"   box = {box[0]:.4f} × {box[1]:.4f} × {box[2]:.4f} Å")
    print(f"参数   : r {R_MIN}–{r_max} step {DR}   q {Q_MIN}–{Q_MAX} step {DQ}")
    print(f"对齐   : {'无（--raw）' if args.raw else f'dump2sq + {SQ_SHIFT_REF:g}'}\n")

    print(f"轨迹排序自检: 帧间同位置同 id = {stable:.4f}")
    if stable < 0.999:
        print("  ⚠ dump 未按 id 排序 —— dump2sq 的 type_new 只在第 0 帧写入，")
        print("    其后各帧的 partial 归类失效，partial 与加权总曲线均不可信。")
    else:
        print("  ✓ 顺序稳定，dump2sq 的 type_new 陈旧问题不触发。\n")

    data = compute(traj, outdir, r_max, args.fe_traj, args.dump2sq)
    exp = load_experiment(xlsx)

    # partial 相互独立性：全部雷同即为归类失效的直接证据
    m = data["q"] > 2.0
    part = data["ref_partials"][m]
    c = np.corrcoef(part.T)
    n = c.shape[0]
    print(f"dump2sq {n} 条 partial 平均非对角相关 (q>2) = {(c.sum() - n) / (n * n - n):+.3f}"
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
