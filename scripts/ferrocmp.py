"""读取 ferro 0.2.0 的分析产物、调用外部参考程序，供 compare_*.py 共用。

0.2.0 起每个分析模式只出**一份 csv**：顶部一个 `#` 注释块（版本、共享参数、
`[inputs]` 逐输入清单），下面是带表头的数据行。注释块明确是给人看的，
`pandas.read_csv(comment="#")` 会把它整块丢掉 —— 所以这里一律用 pandas 读，
不再手写逐行 `float(x) for x in line.split()`（那正是旧脚本遇到 `file` 文本列
直接 ValueError 的原因）。

四处必须知道的形状变化（相对 0.1.15 的 `.dat`）：

1. **`-o` 是文件名后缀，不是路径**。`ferro traj gr -o run1` 写出的是当前目录下的
   `gr_run1.csv`。`run_ferro()` 因此以 outdir 为工作目录调用，并把输入路径转成绝对
   路径 —— 否则产物会散落在脚本的启动目录里。

2. **gr / angle 是长表**，类型在数据列里而不在列名里。取一条曲线是
   `df[(df.center == "P") & (df.neighbor == "O")]`，见 `pick()`。
   sq 仍是宽表（主产物是两条 total，一行一个 q），两者读法不同。

3. **多输入产物带 `file` 列**。这些对拍脚本都是单输入，`one_file()` 断言其唯一后
   丢掉该列；悄悄取第一个文件是最难查的错。

4. 数值一律 `{:.6e}`，缺失值写成**空字段**（列并集下某输入没有的列）。pandas
   读回来就是 NaN，不需要额外处理 —— 但不要把 NaN 当 0 参与统计。

外部参考程序（dump2analysis / dump2sq）的输出格式没变，仍是空白分隔的数值表，
由 `load_legacy()` / `legacy_header()` 读取。
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ── 外部程序调用 ─────────────────────────────────────────────────────────────

def run(cmd, expect, cwd=None, check_returncode=False):
    """执行外部命令并确认它真的产出了文件。

    dump2analysis / dump2sq 的参数校验失败走 exit(0)，退出码不可信，只能验产物；
    ferro 的退出码可信（批内失败也会置 1），由 `check_returncode` 打开。
    """
    print("  $", " ".join(str(c) for c in cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    expect = Path(expect)
    bad = check_returncode and proc.returncode != 0
    if bad or not expect.exists() or expect.stat().st_size == 0:
        sys.exit(
            f"命令未产出 {expect}\n"
            f"  exit={proc.returncode}\n"
            f"  stdout={proc.stdout[-1200:]}\n  stderr={proc.stderr[-1200:]}"
        )
    return expect


def run_ferro(ferro_bin, argv, outdir, product):
    """跑一条 ferro 子命令，返回它写出的产物路径。

    `-o` 是后缀不是路径，故以 `outdir` 为工作目录调用；`argv` 里的输入路径必须
    已经是绝对路径。`product` 是预期文件名（`gr_<后缀>.csv` 这类），相对 outdir。
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    target = outdir / product
    if target.exists():
        target.unlink()      # 不让上一轮的残留冒充本轮产物
    return run([ferro_bin, *[str(a) for a in argv]], target,
               cwd=outdir, check_returncode=True)


# ── ferro csv 读取 ───────────────────────────────────────────────────────────

def meta_lines(path):
    """产物顶部的 `#` 注释块（已去掉前导 `# `）。参数与 `[inputs]` 清单都在里面。"""
    out = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith("#"):
            break
        out.append(line[1:].strip())
    return out


def version(path):
    """首行的 `ferro vX.Y.Z`。"""
    lines = meta_lines(path)
    return lines[0] if lines else "ferro (version unknown)"


def read(path):
    """读 ferro 产物为 DataFrame；`#` 注释块整块丢弃。"""
    df = pd.read_csv(path, comment="#")
    if df.empty:
        sys.exit(f"{path} 中没有数据行")
    return df


def one_file(df, path):
    """断言产物只来自一个输入，丢掉 `file` 列。

    这些脚本是单输入对拍。多输入时静默取第一个会让后面所有数字都对不上，
    且从图上看不出来。
    """
    if "file" not in df.columns:
        return df
    names = df["file"].unique()
    if len(names) != 1:
        sys.exit(f"{path} 含 {len(names)} 个输入 {list(names)}，对拍脚本只接受单输入")
    return df.drop(columns="file")


def pick(df, path, x, **where):
    """长表选行：`pick(df, p, "r", center="P", neighbor="O")`。

    选空即报错并列出该表实际有哪些组合 —— 拼错元素名时最需要看到的就是这个。
    返回按 `x` 升序、索引重置后的子表。
    """
    mask = np.ones(len(df), dtype=bool)
    for col, val in where.items():
        if col not in df.columns:
            sys.exit(f"{path} 中没有列 '{col}'，实际列: {list(df.columns)}")
        mask &= (df[col] == val).to_numpy()
    sub = df[mask]
    if sub.empty:
        have = df[list(where)].drop_duplicates().to_string(index=False)
        sys.exit(f"{path} 中没有 {where} 对应的行。该表实际含有：\n{have}")
    return sub.sort_values(x).reset_index(drop=True)


def column(df, path, name):
    """按列名取一列，缺失时把实际列名一并报出来（宽表用）。"""
    if name not in df.columns:
        sys.exit(f"{path} 中没有列 '{name}'，实际列: {list(df.columns)}")
    return df[name].to_numpy(dtype=float)


# ── 外部参考程序的输出（格式未变） ───────────────────────────────────────────

def load_legacy(path):
    """读空白分隔的数值表，跳过 `#` 注释与空行。dump2analysis / dump2sq 用。"""
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        rows.append([float(v) for v in line.split()])
    if not rows:
        sys.exit(f"{path} 中没有数值行")
    width = len(rows[0])
    if any(len(r) != width for r in rows):
        sys.exit(f"{path} 各行列数不一致")
    return np.array(rows)


def legacy_header(path, marker, sep=None):
    """取含 `marker` 的最后一个注释行作为列名行。"""
    names = None
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") and marker in line:
            names = line.lstrip("# ").split(sep)
    if names is None:
        sys.exit(f"{path} 中找不到含 '{marker}' 的列名行")
    return [n.strip() for n in names if n.strip()]


def same_grid(a, b, label, tol=1e-9):
    """确认两侧横轴逐点一致；不一致说明参数没对齐，继续比下去没有意义。"""
    n = min(len(a), len(b))
    d = float(np.abs(a[:n] - b[:n]).max())
    if d > tol:
        sys.exit(f"{label}: 横轴不一致，最大偏差 {d:.3e} —— 参数未对齐")
    return n
