"""LAMMPS dump 的轨迹自检，供 compare_*.py 共用。

核心是**帧间原子书写顺序是否稳定**这一项，它决定 `dump2sq` 的输出可不可信：

`dump2sq.c:InitializeType` 只在第 0 帧写入 `data[i].type_new`，而
`input.c:ReadDataLammpstrj` 的 sscanf 之后只更新 id/type/elem/x/y/z，**从不更新
type_new**。于是原子类型是按**数组下标**而非实际元素归类的。LAMMPS 未加
`dump_modify sort id` 时，原子书写顺序随 domain decomposition 变化，第 1 帧起
partial 归类即失效。

实测（examples/43Z43P15A_NPT.lammpstrj，帧间同位置同 id 仅 210/2004）：
10 条 partial 相互相关系数 0.959，XRD 与中子总曲线相关 0.9986 —— 加权被完全冲掉，
两条总曲线都退化成总 g(r) 的傅里叶变换（与 FT[自身 total g(r)] 相关 0.9987/0.9999）。
换一条按 id 排序的轨迹（examples/70Z30P00A_NVT.lammpstrj，5003/5003），
ferro 与 dump2sq 的总曲线 max|Δ| 降到 8e-4 ~ 2.7e-3。

`dump2analysis` 每帧按元素重新选原子，不受影响。
"""

from pathlib import Path

# 低于此比例即认为 dump 未按稳定顺序书写
ORDER_STABLE_THRESHOLD = 0.999


def scan_traj(traj, n_probe=3):
    """读前 `n_probe` 帧，返回 (盒子边长, 原子数, 元素计数, 帧间同位置同 id 比例)。"""
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
    if not order:
        raise SystemExit(f"{traj} 中读不到任何帧")
    stable = min(
        (sum(a == b for a, b in zip(order[0], o)) / len(order[0]) for o in order[1:]),
        default=float("nan"),
    )
    return box, len(order[0]), elems, stable


def report_order_stability(stable, verbose=True):
    """打印自检结论，返回 dump2sq 的输出是否可信。"""
    ok = not (stable < ORDER_STABLE_THRESHOLD)
    if not verbose:
        return ok
    print(f"轨迹排序自检: 帧间同位置同 id = {stable:.4f}")
    if ok:
        print("  ✓ 顺序稳定，dump2sq 的 type_new 陈旧问题不触发。")
    else:
        print("  ⚠ dump 未按 id 排序 —— dump2sq 的 type_new 只在第 0 帧写入，")
        print("    其后各帧的 partial 归类失效，partial 与加权总曲线均不可信。")
        print("    此轨迹上 ferro 与 dump2sq 的差异不构成对 ferro 的质疑。")
    return ok


def partial_independence(partials):
    """partial 之间的平均非对角相关系数；接近 1 即为归类失效的直接证据。"""
    import numpy as np

    c = np.corrcoef(partials.T)
    n = c.shape[0]
    if n < 2:
        return float("nan")
    return (c.sum() - n) / (n * n - n)
