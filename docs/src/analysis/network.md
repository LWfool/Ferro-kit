# Glass Network Analysis (`ferro net`)

`ferro net` 对 MD 轨迹逐帧建立玻璃网络拓扑，把每个原子归入一个**结构化的类型**，
再在全轨迹上累加成分布。产物是五张长表，外加可选的标注轨迹。

| 量 | 表 | 含义 |
|---|---|---|
| 桥接配体数 | `bridge` | 形成子连了几个桥联配体（P 的即 Qn），可按伙伴元素分解 |
| 配体分类 | `oxy` | `O_f` / `O_n` / `O_b` / `O_t`，伙伴元素以数据列给出 |
| 配位数 | `cn` | 形成子与修饰子的总配位数分布 |
| 均值 | `mean` | 一元素一行的平均桥接数与平均配位数 |
| 桥联统计 | `linkage` | 每座桥两端的位点状态（元素 + 桥接数 + 配位数） |

> **0.2.1 起命令形状变了。** `ferro net qn` 与 `ferro net type` 已合并为叶子命令
> `ferro net`；导出退化为开关 `--export-traj`。旧标签（`P0` / `Of` / `On_P` /
> `Ob_P_P` / `X` / `Zn_f`）已全部替换，不留读取兼容层。

---

## 理论背景

### 键连关系的判定

以截断半径（cutoff）作为成键判据：两原子间距小于截断值时视为成键。距离计算使用最小
镜像约定（minimum-image convention），支持斜方和三斜晶胞：

$$d_{ij} = \bigl|\mathbf{r}_{ij} - \mathbf{M} \cdot \text{round}\!\left(\mathbf{M}^{-1}\mathbf{r}_{ij}\right)\bigr|$$

分析要求每帧均存在晶胞（PBC）。缺晶胞的帧被跳过；一帧不剩时该输入报错并跳过。

### 三种角色

参数把元素分成三类，这个划分是后面一切的前提：

| 角色 | 来源 | 参与 |
|---|---|---|
| **形成子** (former) | `--<F>-<L>=<Å>` 的左侧 | 桥接数、配位数、配体分类、桥联统计 |
| **配体** (ligand) | `--<F>-<L>=<Å>` 的右侧 | 配体分类 |
| **修饰子** (modifier) | `--modifier` 点名，截断同样用 `--<M>-<L>=<Å>` | **只**参与配位数 |

修饰子被排除在配体分类之外，是这个参数存在的全部理由：一个 Zn 挨着某个非桥氧，
不能因此把那个氧算成桥氧。点名了修饰子却不给它截断会直接报错——静默通过会让它被
当成形成子，整套氧分类随之错位。

### 配体分类

对每个配体原子 $k$，统计其形成子邻居数 $n_k$：

| $n_k$ | 标签 | 含义 |
|---|---|---|
| 0 | `O_f` | Free，游离配体 |
| 1 | `O_n` | Non-bridging，非桥配体 |
| 2 | `O_b` | Bridging，桥联配体 |
| ≥3 | `O_t` | Tricluster，三配位配体（含 Al 玻璃中的常见结构） |

**标签里不带伙伴元素。** `P–O–P` 与 `P–O–Al` 都标 `O_b`，两者的区分走 `oxy` 表的
`former_a` / `former_b` **数据列**。标签保持朴素、统计保持完整，两件事解耦：标签是
给结构文件用的，伙伴分解是给分析用的。

### 桥接配体数（P 的 Qn）

对每个形成子原子 $i$，桥接数为其邻居中被分类为 `O_b` 或 `O_t` 的配体个数：

$$n_i = \bigl|\{k \in \text{neighbors}(i) : n_k \ge 2\}\bigr|$$

三配位配体**计入**（它确实把该形成子连进了网络）。代价是
$\sum(\text{桥接数}) \neq 2 \times |\text{O\_b}|$，因为一个 `O_t` 被三边各数一次。

列名是 `n_bridge` 而不是 `qn`：这个计数对每个形成子都有定义，但 Qn 记号是四面体
形成子的惯例——**Al 没有 Qn 这一说**，描述 Al 请看 `cn` 表。对 P 而言 `n_bridge`
就是 Qn。

### 异核桥分解 $Q^n(m\mathrm{Al})$

`bridge` 表另有 `m_<X>` 列：桥接数中通向元素 X 的那部分，即文献（Brow、Eckert 等）
的 $Q^n(m\mathrm{Al})$ 记号。

- `n_bridge=2, m_Al=1, m_P=1` → $Q^2(1\mathrm{Al})$
- `groupby("n_bridge")` 即退回简单 Qn 分布

只对**恰好两个形成子**的配体成立：三配位配体没有唯一对端，计入 `n_bridge` 但不进任何
`m_<X>`，故 $\sum m \le n_\text{bridge}$，差额就是三配位桥的数目——在表里可见而非隐藏。

这个分解**无法从 `linkage` 表反推**：两座 P–O–Al 桥可能来自一个 $m_\mathrm{Al}=2$ 的 P，
也可能来自两个 $m_\mathrm{Al}=1$ 的 P。`linkage` 数的是桥，`bridge` 数的是原子。

### 桥联统计

每个连接 ≥2 个形成子的配体，把两端形成子的**位点状态**记下来：

$$\text{site} = (\text{元素},\; n_\text{bridge},\; \text{CN})$$

两端**都**带这三个字段。这一点与常见实现不同：那些实现按元素只存一个数（P 存 Qn、
Al 存 CN），于是「4 配位 Al 的桥接数是多少」问不出来。

- **规范半边**：桥联无方向，两端按 `(元素, 桥接数, CN)` 排序后小的在前，每对只存一次。
  故**行和不等于该位点的总参与度**，要算参与度得把 `_a` 与 `_b` 两列都数一遍。
- **`n_formers` 列**：普通桥氧为 2；三配位配体展开成 $C(3,2)=3$ 行并标 3。
  它们不被丢弃，因为「三配位氧连的是谁」在含 Al 体系里正是要研究的东西。

### 统计口径：`fraction` 与 `sd`

- `count`：全轨迹计数之和（整数，精确）
- `fraction`：**逐帧比例的平均**。某帧没有该 bin 时按 0 计入，不是缺失值
- `sd`：同一序列的**样本标准差**（ddof = 1），单帧时为空

> **`sd` 不是标准误。** MD 相邻帧强相关，它既不按 $1/\sqrt{N}$ 收缩，也不估计物理
> 涨落。把它读作「这个数在快照之间摆动多大」。真要误差棒应做分块平均（block
> averaging），本分析不提供。

实现用 Welford 递推而非 $\Sigma f$ / $\Sigma f^2$：后者对恒定序列会给出 ~2e-10 的
抵消噪声，而 `sd` 列里的假抖动会被读者当成物理。

---

## 快速入门

```bash
# P₂O₅ 玻璃，P-O 截断 2.4 Å
ferro net -i traj.lammpstrj --P-O=2.4

# 铝磷酸盐 + Zn 修饰子
ferro net -i traj.lammpstrj --P-O=2.4 --Al-O=2.4 --Zn-O=2.6 --modifier Zn

# 批处理：一次跑多条轨迹，结果堆成一份带 file 列的 csv
ferro net -i 'runs/*/prod.lammpstrj' --P-O=2.4 -o scan

# 只统计尾部 500 帧，8 线程
ferro net -i traj.lammpstrj --P-O=2.4 --last-n 500 --ncore 8

# 同时导出标注轨迹
ferro net -i traj.lammpstrj --P-O=2.4 --export-traj
ferro net -i traj.lammpstrj --P-O=2.4 --export-traj extxyz
```

---

## 命令行参数

### 必需

至少一个 `--<Former>-<Ligand>=<cutoff>` 形式的配对参数（元素首字母大写）：

```
--P-O=2.4        P（形成子）与 O（配体）的截断半径 2.4 Å
--Si-O=1.8       Si-O 截断 1.8 Å
--Al-O=2.4       Al-O 截断 2.4 Å
--Al-F=2.1       同一形成子可有多种配体
```

元素对写在**参数名**里，clap 建模不了，所以 `main` 在解析前先把它们从 argv 剥离。

### 可选

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-i <FILE>...` | — | 输入轨迹，多值并自展开 glob（引号括起来）。缺省时打印帮助 |
| `-o <SUFFIX>` | — | 输出**文件名后缀**（不是路径）：`network_<表>_<后缀>.csv` |
| `--last-n N` | 全部帧 | 仅使用尾部 N 帧 |
| `--ncore N` | 全部核心 | 并行线程数 |
| `--metal-units` | 关 | LAMMPS metal 单位。统计不读速度/力，**只对 `--export-traj extxyz` 有影响** |
| `--modifier E,E` | — | 只计配位数的元素，逗号分隔。须同时给出各自的截断 |
| `--export-traj [FMT]` | — | 另写标注轨迹，`lammpstrj`（默认）或 `extxyz` |

---

## 输出

五张 CSV，各带 `file` 列，`#` 注释块里是共享参数与 `[inputs]` 清单
（`pandas.read_csv(comment="#")` 会自动丢掉）。

| 文件 | 列 | 一行一个 |
|---|---|---|
| `network_bridge.csv` | `file, former, n_bridge, m_<X>…, count, fraction, sd` | (形成子, 桥接数, 伙伴分解) |
| `network_oxy.csv` | `file, type, former_a, former_b, count, fraction, sd` | (配体类型, 伙伴对) |
| `network_cn.csv` | `file, element, cn, count, fraction, sd` | (元素, 配位数) |
| `network_mean.csv` | `file, element, mean_n_bridge, mean_cn` | 元素 |
| `network_linkage.csv` | `file, elem_a, n_bridge_a, cn_a, elem_b, n_bridge_b, cn_b, n_formers, count, fraction, sd` | 桥的两端状态组合 |

三张分布表的值列语义各不相同（桥接数 / 类型标签 / 配位数），并成一张会让
`groupby` 无意义，故保持三张。`mean` 是「一元素一行」的另一种粒度，故独立成表，
而不是在每行重复一个稀疏列。

批处理时元素集不同的输入取**列并集，缺的留空（NaN），不补零、不插值**——没有 Al
的体系其 `m_Al` 列为空，而不是 0。

### 示例：`oxy` 表

```csv
file,type,former_a,former_b,count,fraction,sd
43Z43P15A,O_n,P,,2985,4.543379e-1,1.203302e-3
43Z43P15A,O_b,Al,Al,10,1.522070e-3,0.000000e0
43Z43P15A,O_b,Al,P,2693,4.098935e-1,1.154167e-3
43Z43P15A,O_b,P,P,880,1.339422e-1,0.000000e0
43Z43P15A,O_t,,,2,3.044140e-4,4.168360e-4
```

`O_t` 的伙伴列留空——它不是一对，其配对关系由 `linkage` 表的 `n_formers=3` 行给出。

### 用 pandas 分析 `linkage`

```python
import pandas as pd
d = pd.read_csv("network_linkage.csv", comment="#")

# 按元素对汇总
d.groupby(["elem_a", "elem_b"])["count"].sum()

# Al-O-Al 发生在哪种配位的 Al 之间
al = d.query("elem_a == 'Al' and elem_b == 'Al'")
al.pivot_table(index="cn_a", columns="cn_b", values="count", aggfunc="sum")

# P-O-P 的 Qn–Qn 连接矩阵
d.query("elem_a == 'P' and elem_b == 'P'").pivot_table(
    index="n_bridge_a", columns="n_bridge_b", values="count", aggfunc="sum")

# 只看真桥，排除三配位配体
d.query("n_formers == 2")
```

---

## 标注轨迹（`--export-traj`）

导出每个输入一个文件：`<输入 stem>_types[_<后缀>].<ext>`。文件名必须掺输入 stem，
否则批处理时第二个输入会覆盖第一个（与 `ferro map` 的 cube 同理）。

### 标签

| 角色 | 标签 |
|---|---|
| 形成子 | `P_0` `P_1` … `Al_2`（数字 = 桥接配体数） |
| 自由配体 | `O_f` |
| 非桥配体 | `O_n` |
| 桥联配体 | `O_b` |
| 三配位配体 | `O_t` |
| 修饰子 | `Zn`（元素符号，无角色后缀） |

全部合 `<元素>_<后缀>` 约定，按**首个下划线**拆分。修饰子不带角色后缀：旧方案
按非桥配体数 0/1/2/≥3 分档，而修饰子的实际配位数在 3–6，实测参考轨迹 97% 落进
兜底桶，没有分辨力。

### 两种格式的差别

| | `lammpstrj`（默认） | `extxyz` |
|---|---|---|
| 标签存放 | 折进 `element` 列 | 独立的 `label:S:1` 列 |
| `species` / `element` | 是标签（`P_2`） | 是纯元素（`P`） |
| 读回 | reader 按首个下划线拆成 element + label | 直接分列读取，无需猜测 |
| 下游兼容 | 不新增列，已有工具照常解析 | 新增一列 |

dump 折叠是因为它只有一列名字放得下第二个名字；extxyz 的列是自描述的，不需要将就。
**只有 `ferro net --export-traj` 会折叠**——`ferro convert` 无论标签如何都写干净的
元素符号。折叠还带守卫：标签不形如 `<元素>_…` 时不折并计数告警（CIF 的 `O1`、CP2K
的 `Fe1` 折进去后读回来就是不存在的元素 `O1`）。

### 下游选型

```bash
# 按元素选 —— 任意帧数都可以
ferro traj gr -i run_types.lammpstrj -a P -b O

# 按标签选 —— 需要单帧
ferro traj gr -i run_types.lammpstrj -x P_3 -y O_b --last-n 1
```

> **按标签选型只在单帧成立。** `g(r)` 要求逐类型的粒子数守恒，而标签是动态的——
> 一个 P 在轨迹推进中会改变自己的 Qn（实测参考轨迹 `P_3` 逐帧为
> 149 / 152 / 150 / 150 / 150）。对多帧标注轨迹按标签选型会被
> 「per-type atom counts change」错误拒绝，这是守卫而非缺陷。多帧请按元素选。

标注轨迹与统计表是互相独立的产物：想验证分类是否自洽，可对导出轨迹跑
`ferro traj gr -x P_3 -y O_b --last-n 1`，第一壳层 CN 应恰为 3——分类层说
「这个 P 有 3 个桥氧」，独立的 CN 积分路径应当数回 3。

---

## 典型参数参考

| 体系 | 形成子 | 修饰子 | 参考截断 |
|---|---|---|---|
| P₂O₅ 玻璃 | P | — | P-O: 2.3–2.4 Å |
| SiO₂ 玻璃 | Si | — | Si-O: 1.8 Å |
| Al₂O₃ | Al | — | Al-O: 2.1–2.4 Å |
| GeO₂ 玻璃 | Ge | — | Ge-O: 2.0 Å |
| ZnO–P₂O₅ 玻璃 | P | Zn | P-O: 2.4 Å，Zn-O: 2.6 Å |
| ZnO–Al₂O₃–P₂O₅ | P, Al | Zn | P-O: 2.4，Al-O: 2.4，Zn-O: 2.6 Å |

截断值应参照 g(r) 第一峰谷位置确定：

```bash
ferro traj gr -i traj.lammpstrj -a P -b O --r-max 5 --plot
```
