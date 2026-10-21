"""Read ferro's analysis products and run the external reference programs; shared by
compare_*.py.

Since 0.2.0 every analysis writes **one csv**: a `#` comment block on top (version, shared
parameters, the `[inputs]` list per input), then a header and the data rows. The comment block
is for people, and `pandas.read_csv(comment="#")` drops it whole — so everything here is read
with pandas, never line by line with `float(x) for x in line.split()` (which is how the old
scripts died with a ValueError on the `file` text column).

Five changes of shape to know (relative to the `.dat` of 0.1.15):

1. **The batch tag is `-s`, and `-o` is the output directory** (since 0.3.3; in 0.2.0
   `-o` was the tag). `ferro traj gr -s run1` writes `gr_all_run1.csv` into the current
   directory. `run_ferro()` runs with outdir as the working directory and makes the
   input paths absolute, so it does not need `-o`.

2. **gr / angle are long tables**: the species are in data columns, not in column names. One
   curve is `df[(df.center == "P") & (df.neighbor == "O")]`, see `pick()`. sq is still a wide
   table (the main products are the two totals, one row per q); the two are read differently.

3. **A multi-input product has a `file` column.** These comparisons are all single-input;
   `one_file()` asserts it is unique and drops the column — silently taking the first file is
   the hardest mistake to find.

4. Numbers are `{:.6e}`; a missing value (a column another input lacks, under the column
   union) is an **empty field**. pandas reads it back as NaN with no extra handling — but do
   not count NaN as 0.

5. **Product names carry a label segment** (since 2026-08-13): `gr.csv` → `gr_P-O.csv`,
   `gr_all.csv` without a selection. Build every name with `product_name()`, never a string at
   the call site — when the label segment came, three scripts broke at once, each having
   written its own names.

The external reference programs (dump2analysis / dump2sq) kept their output format,
whitespace-separated numbers, read by `load_legacy()` / `legacy_header()`.
"""

import os
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

    批次后缀走 `-s`，产物落在当前目录，故以 `outdir` 为工作目录调用；`argv` 里的
    输入路径必须已经是绝对路径（0.3.3 起 ferro 的 `-o` 是输出目录，这里用不上它：
    cwd 就是 outdir）。`product` 是预期文件名，用 `product_name()` 拼，别手写。
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    target = outdir / product
    if target.exists():
        target.unlink()      # 不让上一轮的残留冒充本轮产物
    # cwd 换成了 outdir，`--ferro ../target/release/ferro` 这种相对路径会解析不到
    if os.sep in str(ferro_bin):
        ferro_bin = Path(ferro_bin).resolve()
    return run([str(ferro_bin), *[str(a) for a in argv]], target,
               cwd=outdir, check_returncode=True)


# ── 产物文件名 ───────────────────────────────────────────────────────────────
#
# 命名规则在 ferro-cli 的 `batch::product_name` / `batch::file_label`，这里是它的镜像。
# 手写产物名的代价已经付过一次：2026-08-13 给 traj 产物加 label 段，三个脚本同时
# 断掉。规则再变时只改这一处。

_LABEL_OK = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_+-")


def file_label(*parts):
    """有序选择的 label：`file_label("P", "O") -> "P-O"`，不给参数 -> `"all"`。

    对应 `-a`/`-b`/`-c` 这类**按写的顺序**进文件名的选择：CN(r) 有向，
    `-a P -b O` 与 `-a O -b P` 是两份不同数据，必须落到两个文件名下。
    """
    parts = [str(p) for p in parts if p is not None]
    if not parts:
        return "all"
    for p in parts:
        if not p:
            sys.exit("空的类型选择：元素或位点标签不能为空串")
        bad = [c for c in p if c not in _LABEL_OK]
        if bad:
            sys.exit(f"'{p}' 含非法字符 {bad!r}：label 会成为文件名的一段，"
                     f"ferro 只接受字母、数字、'_'、'+'、'-'")
    return "-".join(parts)


def product_name(mode, label=None, suffix=None, table=None, input=None, ext="csv"):
    """ferro 的产物文件名：`<mode>[_<table>][_<label>][_<input>][_<suffix>].<ext>`。

    * `table` 与 `mode` 同名时省略（gr / angle / sq 这类单表模式）；
      net 的六张表要显式传（`product_name("network", table="qn", ...)`）。
    * `label` 是**该模式有没有类型选择**：sq 没有，传 `None`；gr / angle 有，
      即使没筛选也得传 `file_label()` 拿到的 `"all"`。两者的区别就是差一段。
    * `input` 是逐输入产物（map 的 cube、net 导出、bader）的输入标签，单输入也带；
      traj / net 的堆叠表没有这一段。
    """
    stem = mode if table in (None, mode) else f"{mode}_{table}"
    for seg in (label, input, suffix):
        if seg:
            stem += f"_{seg}"
    return f"{stem}.{ext}"


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
