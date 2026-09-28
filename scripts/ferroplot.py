"""发表级绘图的共享层，供 plot_*.py 使用。

与 `ferrocmp.py` 的分工：那边是**对拍**（跟参考实现逐点比对，图只为看差异），
这边是**出图**（进论文的 pdf）。两者都读 ferro 的 csv，但读完之后没有共同代码，
所以是两个模块而不是一个。

## 样式

`['science', 'vibrant']` + LaTeX。`science` 是四边框 + 内向刻度，`vibrant` 是
Paul Tol 的 7 色色盲安全序列。两条实测结论：

- **不需要给字符串转义**。matplotlib 的 usetex 前导会加载 `underscore` 宏包
  （`[strings]` 选项），`43Z43P15A_NPT_5`、`P_3`、`O_b`、`50% Zn` 都原样通过。
- 仍会让**整张图渲染失败**的只有 `&`、`#`、`^` 三个字符。元素符号与位点标签不可能
  含它们（ferro 的文件名校验只放行 `[A-Za-z0-9_+-]`），只有自定义文件名里出现才会
  碰上。真碰上了就改文件名，不在这里做转义层——转义层意味着每个进图的字符串都要
  记得过一遍，漏一个就是同样的渲染失败，只是更难找。

## 颜色

颜色**恒定按 `file` 列分配**（一条轨迹一个颜色），跨子图、跨图一致。这是全部四个
脚本共同的语义：`file` 是成分/体系，是读者要追踪的那一维。

## 显示名：`labels.csv`

产物目录里一份**全部脚本共用**的 `labels.csv`（列 `key,label,show`）。成分名在
论文里全文统一，不该在 g(r) 图与 Qn 图里各改一遍。

- 第一次运行自动生成，`label` 取原名、`show` 取 1
- **行序即显示顺序**（x 轴、图例、颜色分配）；`show=0` 的不画
- 之后每次运行：已有行一字不动；本次新出现的 key 追加到末尾并打 `Note:`（打错的
  key 会以新行的形式出现在这条提示里，这是替代「报错」的那道保护）；文件里有而本次
  数据里没有的 key 静默忽略 —— 那是别的批次的成分，不是错误
- key 既可以是成分名（`file` 列），也可以是多 csv 并排时的组名（文件名后缀）

## 数据导出：`<stem>_data.csv`

每次出图都在 pdf 旁写一份长表，内容就是图上画的那些数（平滑、筛选之后的），
供按论文风格重绘。原始 `file` 与映射后的 `label` 两列都在，行序即绘制顺序。
"""

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import scienceplots  # noqa: F401  —— import 即注册样式，不直接调用

# ── 全局导出设置 ─────────────────────────────────────────────────────────────

STYLE = ["science", "vibrant"]

#: png 只为看效果，pdf 才是产物
PNG_DPI = 600

#: 42 = TrueType，保证 pdf 里的文字可被 Illustrator 选中编辑
PDF_FONTTYPE = 42


def use_style():
    """装上绘图样式。每个脚本在 main 里调一次。"""
    plt.style.use(STYLE)
    mpl.rcParams["pdf.fonttype"] = PDF_FONTTYPE
    mpl.rcParams["ps.fonttype"] = PDF_FONTTYPE


# ── 读取 ─────────────────────────────────────────────────────────────────────

def read(path):
    """读一份 ferro 产物。

    `#` 块是给人看的（版本、共享参数、`[inputs]` 清单），`comment="#"` 整块丢掉。
    缺失值是空字段 → NaN，**不要当 0 参与统计**。
    """
    path = Path(path)
    df = pd.read_csv(path, comment="#")
    if df.empty:
        raise SystemExit(f"{path} 没有数据行（只有表头？）")
    return df


def files_in(df):
    """`file` 列的取值，**保持出现顺序**。

    顺序不是字典序：它等于 `ferro -i` 的参数顺序（`expand_inputs` 保序，
    `concat_union` 按 parts 顺序拼），所以成分序列的次序是在命令行上控制的。
    分类 x 轴直接用这个顺序，不要 sort。
    """
    return list(dict.fromkeys(df["file"]))


def suffix_of(path):
    """从产物文件名取批次名，供多文件对比时的组标签用。

    `network_qn_CMD.csv` → `CMD`。取不到（`network_qn.csv` 这种没给 `-o` 的常规
    用法）就退回整个 stem —— 单文件时本就不需要为了画图去编一个 `-o`。
    """
    stem = Path(path).stem
    for mode in ("network_qn_partner", "network_qn", "network_coordination",
                 "network_composition", "network_ligand_type", "network_linkage"):
        if stem.startswith(mode + "_"):
            return stem[len(mode) + 1:]
    return stem


# ── 颜色 ─────────────────────────────────────────────────────────────────────

def color_map(names):
    """把一组名字映射到样式色序，循环取用。

    调用点一律传 `files_in(df)`，所以同一条轨迹在任何一张图里都是同一个颜色。
    """
    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    return {n: cycle[i % len(cycle)] for i, n in enumerate(names)}


def shades(base, n, lo=0.35):
    """由一个基色派生 n 级明度，用于嵌套堆积柱的内层。

    深→浅。外层色相仍然携带主结构（Qn），内层只做深浅，所以第一眼落在色相边界上
    而不是深浅上——这正是「补充而不抢焦点」的做法。
    """
    rgb = mpl.colors.to_rgb(base)
    if n <= 1:
        return [rgb]
    # 朝白色插值，最浅一级仍保留 lo 的原色，避免淡到看不见
    return [tuple(c + (1.0 - c) * (1.0 - lo) * i / (n - 1) for c in rgb) for i in range(n)]


# ── 显示名与导出 ─────────────────────────────────────────────────────────────

LABELS_FILE = "labels.csv"


class Context:
    """一张图的共享状态：显示名、可见成分及其顺序、颜色、待导出的行。

    `files` 是本次全部输入里出现过的成分（出现顺序），`groups` 是多 csv 并排时的
    组名。两者都登记进 `labels.csv`，但只有 `files` 参与颜色分配。
    """

    def __init__(self, files, outdir=None, groups=()):
        self.path = (Path(outdir) if outdir else Path.cwd()) / LABELS_FILE
        table = load_labels(self.path, [*files, *groups])
        self._label = {k: lab for k, lab, _ in table}
        shown = {k for k, _, show in table if show}
        present = set(files)
        self.names = [k for k, _, _ in table if k in present and k in shown]
        if not self.names:
            raise SystemExit(f"{self.path} 把本次全部成分都设成了 show=0，没有可画的")
        self._groups_shown = {g for g in groups if g in shown}
        self.colors = color_map(self.names)
        self.rows = []

    def shown(self, df):
        """df 里出现、且 `show=1` 的成分，按 `labels.csv` 的行序。"""
        present = set(df["file"])
        return [n for n in self.names if n in present]

    def group_shown(self, group):
        return group in self._groups_shown

    def label(self, key):
        return self._label.get(key, key)

    def record(self, **row):
        """登记一行导出数据。`label` 列由 `file` 列自动补上，紧跟在它后面。"""
        out = {}
        for k, v in row.items():
            out[k] = v
            if k == "file":
                out["label"] = self.label(v)
        self.rows.append(out)

    def record_curve(self, panel, file, quantity, x, y):
        """登记一条曲线：一点一行，列 `panel, file, label, quantity, x, y`。"""
        lab = self.label(file)
        self.rows += [{"panel": panel, "file": file, "label": lab, "quantity": quantity,
                       "x": float(a), "y": float(b)} for a, b in zip(x, y)]


def context(frames, outdir=None, groups=()):
    """frames: [(path, df)]。成分取全部产物里出现过的 `file`，保持出现顺序。"""
    files = list(dict.fromkeys(f for _, df in frames for f in files_in(df)))
    return Context(files, outdir, groups)


def load_labels(path, keys):
    """读 `labels.csv`，把本次新出现的 key 追加进去，返回 [(key, label, show)]。

    只在有新 key 时写盘；已有行原样写回（顺序、改过的 label、show 都不动）。
    """
    table = []
    if path.exists():
        with open(path, newline="", encoding="utf-8") as fh:
            for rec in csv.DictReader(fh):
                show = str(rec.get("show", "1")).strip()
                if show not in ("0", "1"):
                    raise SystemExit(f"{path}: key '{rec['key']}' 的 show 只能是 0 或 1，给的是 '{show}'")
                table.append((rec["key"], rec["label"], show == "1"))
    known = {k for k, _, _ in table}
    new = [k for k in dict.fromkeys(keys) if k not in known]
    if new:
        first = not path.exists()
        table += [(k, k, True) for k in new]
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["key", "label", "show"])
            w.writerows((k, lab, int(show)) for k, lab, show in table)
        if first:
            print(f"Note  : 已生成 {path}，改 label / show / 行序后重跑即生效")
        else:
            print(f"Note  : {path} 追加了 {len(new)} 个新 key：{', '.join(new)}")
    return table


def save(fig, stem, outdir=None, ctx=None):
    """同时写 pdf（产物）与 png（看效果），有登记的数据时再写 `<stem>_data.csv`。

    `bbox_inches="tight"` 由 `science` 样式设好，这里不重复指定，否则显式给的
    figsize 会被裁掉边距后失真。
    """
    out = Path(outdir) if outdir else Path.cwd()
    out.mkdir(parents=True, exist_ok=True)
    pdf, png = out / f"{stem}.pdf", out / f"{stem}.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=PNG_DPI)
    print(f"  -> {pdf}")
    print(f"  -> {png}")
    if ctx is not None and ctx.rows:
        data = out / f"{stem}_data.csv"
        pd.DataFrame(ctx.rows).to_csv(data, index=False)
        print(f"  -> {data}")
    return pdf


def add_outdir_arg(parser):
    """`--outdir`：出图脚本自己的输出目录。

    ferro CLI 侧这个参数已并进 `-o`（0.3.3 起 `-o` 恒为路径），这里保持
    `--outdir` 不动 —— 绘图脚本的 `-o` 另有含义（产物名后缀）。
    """
    parser.add_argument(
        "--outdir", type=Path, default=None,
        help="产物目录（不存在则创建；默认当前目录）",
    )
    return parser
