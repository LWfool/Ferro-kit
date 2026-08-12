# Glass Network Analysis (`ferro net`)

`ferro net` 对 MD 轨迹逐帧建立玻璃网络拓扑，把每个原子归入一个**结构化的类型**，
再在全轨迹上累加成分布。产物是六张长表，外加可选的标注轨迹。

| 量 | 表 | 含义 |
|---|---|---|
| **结构组成** | `composition` | 一物种一行：`P-Q2` `Al_4` `O_b` `Zn_4`，各占其元素的比例 |
| Qn 分布 | `qn` | Qn 形成子连了几个桥联配体 |
| 异核桥分解 | `qn_partner` | 上表再按伙伴元素拆一维，即 $Q^n(m\mathrm{Al})$ |
| 配体分类 | `ligand_type` | `O_f` / `O_n` / `O_b` / `O_t`，伙伴元素以数据列给出 |
| 配位数 | `coordination` | 形成子与修饰子的总配位数分布 |
| 桥联统计 | `linkage` | 每座桥的配体元素与两端位点状态 |

`composition` 是**其余表的摘要**，不是新的测量：想一眼看完这块玻璃的结构组成就读它，
要伙伴分解、要两端状态才去翻源表。

> **Qn 只报给 Qn 元素。** Qn 是四面体形成子的记号，默认只有 `B` / `P` / `Si`
> 出现在 `qn` 与 `qn_partner` 两张表里；**Al 之类的形成子由配位数刻画**
> （文献写作 Al[4] / Al[5] / Al[6]），它们仍然参与桥氧判定、`ligand_type` 与
> `linkage`，退出的只是这两张表的「行」。用 `--qn` 可以整体替换这个名单。

> **命令形状。** `ferro net qn` 与 `ferro net type` 已合并为叶子命令
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

### 两个必须分清的量：桥接数与配位数

对每个形成子原子 $i$：

$$n_i^\text{bridge} = \bigl|\{k \in \text{neighbors}(i) : n_k \ge 2\}\bigr|,
\qquad
\mathrm{CN}_i = \bigl|\text{neighbors}(i)\bigr|$$

**桥接数只数桥联配体，配位数数截断内的全部配体**，含非桥配体。三配位配体计入桥接数
（它确实把该形成子连进了网络），代价是
$\sum(\text{桥接数}) \neq 2 \times |\text{O\_b}|$——一个 `O_t` 被三边各数一次。

> **两者相等是体系的性质，不是定义。** 一个不带非桥配体的形成子会让
> $n^\text{bridge} = \mathrm{CN}$ 处处成立。参考轨迹里的 Al 正是如此（`ligand_type`
> 表中没有任何 `O_n, Al` 行），于是两张分布表逐档相同。换一个 Al 带非桥氧的体系，
> 两者立刻分叉。凡是读到「Al 的配位数」，认准 `coordination` 表。

### Qn 只给 Qn 形成子

$Q^n$ 是四面体形成子的记号：在配位数基本固定的位点上，一个数字就说尽了它的连接
状况。Al 不满足这个前提——它的配位数本身就是要报的量。因此：

| | 出现在 `qn` / `qn_partner` | 标签数字 | 由什么刻画 |
|---|---|---|---|
| Qn 形成子（默认 B, P, Si） | 是 | Qn（桥接数） | Qn 分布 |
| 其他形成子（Al, …） | **否** | **配位数** | `coordination` 表 |
| 修饰子（Zn, …） | 否 | 无后缀 | `coordination` 表 |

名单是**默认值**，`--qn Si,Al` 会整体替换它（不是叠加）——某个元素算不算 Qn 形成子
是体系的性质：铝硅酸盐里人们确实会报 Al 的 $Q^n(m\mathrm{Si})$。

### 异核桥分解 $Q^n(m\mathrm{Al})$

`qn_partner` 表在 `qn` 之上多一维：`m_<X>` 是桥接数中通向元素 X 的那部分，即文献
（Brow、Eckert 等）的 $Q^n(m\mathrm{Al})$ 记号。

- `qn=2, m_Al=1, m_P=1` → $Q^2(1\mathrm{Al})$
- `qn` 就是它对伙伴维度的**边际**

**伙伴分解不编码进 `label`。** 试过写成 `Q2(1Al)`，行不通：文献的 $Q^n(mX)$ 记号
表达不了「该桥的配体连了三个形成子」这种情形。参考轨迹里
`(qn=2, m_Al=1, m_P=0)` 与 `(qn=2, m_Al=1, m_P=1)` 两行都会被写成 `Q2(1Al)`——
前者有一座桥是三簇氧，$\sum m = 1 < qn$——而没有任何一列能把它们分开。那一维本就
在 `m_<X>` 列里，标签只写 `P-Q2`。

注意 Al 在这里是**伙伴而非主语**：它没有自己的行，但 `m_Al` 列照常存在。

两者是**两张表而不是一张加 `groupby`**：简单 Qn 分布是主产物，必须打开文件就能读到，
不能要求先做聚合。这与 `average` 独立成表是同一条判据——粒度不同就该分表。
`sd` 也必须各自累加而不是相加：相关项之和的方差不等于方差之和。实测参考轨迹
P 的 `qn=2` 一行，正确 `sd` 为 5.574e-3，而把 `qn_partner` 对应各行的 `sd` 相加得
9.966e-3，差近一倍。

只对**恰好两个形成子**的配体成立：三配位配体没有唯一对端，计入 `qn` 但不进任何
`m_<X>`，故 $\sum m \le n_\text{qn}$，差额就是三配位桥的数目——在表里可见而非隐藏。

这个分解**无法从 `linkage` 表反推**：两座 P–O–Al 桥可能来自一个 $m_\mathrm{Al}=2$ 的 P，
也可能来自两个 $m_\mathrm{Al}=1$ 的 P。`linkage` 数的是桥，`qn` 数的是原子。

### 桥联统计

每个连接 ≥2 个形成子的配体，记下一条

$$(\underbrace{\text{元素},\; n_\text{bridge},\; \text{CN}}_{\text{A 端}}),\;
  (\underbrace{\text{元素},\; n_\text{bridge},\; \text{CN}}_{\text{B 端}}),\;
  \text{配体元素},\; n_\text{formers}$$

两端**都**带三个字段。这一点与常见实现不同：那些实现按元素只存一个数（P 存 Qn、
Al 存 CN），于是「4 配位 Al 的桥接数是多少」问不出来。

- **`linkage` 展示列**：`Al_4-O-P_2` 这样的人可读形式，数字按各自约定（Qn 形成子取
  Qn，其余取配位数），与导出轨迹的标签同一套词汇。**筛选请用数值列**，不要正则拆
  这一列。
- **`ligand` 列**：桥中间那个原子的元素。多配体体系里 `Al-O-P` 与 `Al-F-P` 是两种桥，
  永不共用一行——合并计数会报出一个实验无法对应的数。
- **规范半边**：桥联无方向，两端按 `(元素, 桥接数, CN)` 排序后小的在前，每对只存一次。
  故**行和不等于该位点的总参与度**，要算参与度得把 `_a` 与 `_b` 两列都数一遍。
- **`n_formers` 列**：普通桥氧为 2；三配位配体展开成 $C(3,2)=3$ 行并标 3。
  它们不被丢弃，因为「三配位氧连的是谁」在含 Al 体系里正是要研究的东西。
- **`n_bridge_a/b` 保留此名而不叫 `qn_a/b`**：桥接数对每个形成子都有定义，Qn 只对
  一部分成立。这是 `linkage` 里唯一还需要给「非 Qn 形成子的桥接数」命名的地方。

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
| `--qn E,E` | `B,P,Si` | 报 Qn 的形成子，逗号分隔。**替换**默认名单而非叠加 |
| `--export-traj [FMT]` | — | 另写标注轨迹，`lammpstrj`（默认）或 `extxyz` |

---

## 输出

六张 CSV，各带 `file` 列。每个文件的 `#` 注释块里是共享参数、`[inputs]` 清单，
**以及这张表自己的逐列说明**（`pandas.read_csv(comment="#")` 会自动丢掉整块）。
所以下表只需记住哪个文件装什么，列的含义打开文件即可。

| 文件 | 一行一个 | 列 |
|---|---|---|
| `network_composition.csv` | 物种 | `label, element, count, fraction, sd` |
| `network_qn.csv` | (形成子, Qn) | `label, former, qn, count, fraction, sd` |
| `network_qn_partner.csv` | (形成子, Qn, 伙伴分解) | `label, former, qn, m_<X>…, count, fraction, sd` |
| `network_ligand_type.csv` | (配体类型, 伙伴对) | `label, type, former_a, former_b, count, fraction, sd` |
| `network_coordination.csv` | (元素, 配位数) | `element, cn, count, fraction, sd` |
| `network_linkage.csv` | (配体, 两端状态) | `linkage, ligand, elem_a, n_bridge_a, cn_a, elem_b, n_bridge_b, cn_b, n_formers, count, fraction, sd` |

`label` 列是给人读的锚点，`former` / `qn` / `cn` 等数值列是给筛选和画图用的——两者
并存而非二选一，否则「筛出 Qn ≥ 3 的」就得去切字符串。

分布表的值列语义各不相同（Qn / 类型标签 / 配位数），并成一张会让 `groupby` 无意义，
故各自成表。`qn` 与 `qn_partner`、以及 `average`，都是**粒度不同故分表**。

**没有 Qn 形成子时，前两个文件整个不出**，屏幕上打印原因。只有表头的 CSV 读起来像
「测了，结果是零」，而实情是压根没测。

批处理时元素集不同的输入取**列并集，缺的留空（NaN），不补零、不插值**——没有 Al
的体系其 `m_Al` 列为空，而不是 0。

### 示例：`network_composition.csv`

```csv
file,label,element,count,fraction,sd
43Z43P15A,P-Q0,P,25,1.344086e-2,0.000000e0
43Z43P15A,P-Q1,P,276,1.483871e-1,2.249086e-3
43Z43P15A,P-Q2,P,653,3.510753e-1,5.574312e-3
43Z43P15A,P-Q3,P,751,4.037634e-1,2.944745e-3
43Z43P15A,P-Q4,P,155,8.333333e-2,1.900825e-3
43Z43P15A,Al_4,Al,590,8.939394e-1,1.417294e-2
43Z43P15A,Al_5,Al,63,9.545455e-2,1.570943e-2
43Z43P15A,Al_6,Al,7,1.060606e-2,4.149413e-3
43Z43P15A,Zn_3,Zn,114,1.225806e-1,4.844680e-2
43Z43P15A,Zn_4,Zn,723,7.774194e-1,2.860094e-2
43Z43P15A,Zn_5,Zn,92,9.892473e-2,2.331127e-2
43Z43P15A,Zn_6,Zn,1,1.075269e-3,2.404374e-3
43Z43P15A,O_n,O,2985,4.543379e-1,1.203302e-3
43Z43P15A,O_b,O,3583,5.453577e-1,1.154167e-3
43Z43P15A,O_t,O,2,3.044140e-4,4.168360e-4
```

**分母恒为该元素的原子数**：Q2 占全部 P、`Al_4` 占全部 Al、`O_b` 占全部 O。所以
**每个元素的 `fraction` 求和为 1** —— 这是一条打开文件就能核对的恒等式。

每个元素只出现**一种刻画**：Qn 形成子出 Qn，其余形成子与修饰子出配位数，配体出类型。
P 在这里没有配位数行（它在 `network_coordination.csv` 里）。

`O_b` 那行是从 `ligand_type` 的三行（Al-Al / Al-P / P-P）聚合来的，`sd` **逐帧重新
累加**而非把三行相加——相关项之和的方差不等于方差之和。

原来的 `network_average.csv` 已删除。它的两个均值仍在每个文件的 `[inputs]` 块里
（`mean_qn P=2.40  mean_cn Al=4.12 P=4.00 Zn=3.98`），也可从本表精确复算：
$\sum_n n \cdot f_n = 2.395161$。

### 示例：`network_qn.csv`

```csv
file,label,former,qn,count,fraction,sd
43Z43P15A,P-Q0,P,0,25,1.344086e-2,0.000000e0
43Z43P15A,P-Q1,P,1,276,1.483871e-1,2.249086e-3
43Z43P15A,P-Q2,P,2,653,3.510753e-1,5.574312e-3
43Z43P15A,P-Q3,P,3,751,4.037634e-1,2.944745e-3
43Z43P15A,P-Q4,P,4,155,8.333333e-2,1.900825e-3
```

Al 不在其中——同一次运行里它出现在 `network_coordination.csv` 的 `cn` 4/5/6 三行。

### 示例：`network_ligand_type.csv`

```csv
file,label,type,former_a,former_b,count,fraction,sd
43Z43P15A,P-O_n,O_n,P,,2985,4.543379e-1,1.203302e-3
43Z43P15A,Al-O_b-Al,O_b,Al,Al,10,1.522070e-3,0.000000e0
43Z43P15A,Al-O_b-P,O_b,Al,P,2693,4.098935e-1,1.154167e-3
43Z43P15A,P-O_b-P,O_b,P,P,880,1.339422e-1,0.000000e0
43Z43P15A,O_t,O_t,,,2,3.044140e-4,4.168360e-4
```

`label` 把整行读成一句 `<former_a>-<type>-<former_b>`；缺位**不补占位符**，
所以自由氧就是 `O_f`、非桥氧是 `P-O_n`、三簇氧是 `O_t`。`former_a` / `former_b`
仍是独立列，Al-O-Al 的查询照旧是
`query("former_a=='Al' and former_b=='Al'")`，不必去切字符串。

`O_t` 的伙伴列留空——它不是一对，其配对关系由 `linkage` 表的 `n_formers=3` 行给出。
注意这里**没有 `Al-O_n` 行**：这个体系的 Al 不带非桥氧，正是它的桥接数与配位数
处处相等的原因。

> **`fraction` 的分母是该配体元素的原子数**，不是全体配体原子。这是文献 BO fraction
> 的口径；「桥氧占 O+F 的比例」不对应任何常用量。单配体体系两者相同，
> `--Al-O` + `--Al-F` 的体系才有区别。

### 示例：`network_linkage.csv`

```csv
file,linkage,ligand,elem_a,n_bridge_a,cn_a,elem_b,n_bridge_b,cn_b,n_formers,count,...
43Z43P15A,Al_4-O-Al_4,O,Al,4,4,Al,4,4,2,10,...
43Z43P15A,Al_4-O-P_3,O,Al,4,4,P,3,4,2,1277,...
43Z43P15A,Al_5-O-P_2,O,Al,5,5,P,2,4,2,83,...
```

`Al_4` 里的 4 是**配位数**，`P_3` 里的 3 是 **Qn**——这是文献自己的约定（Al[4] 对
Q³），逐文件的 `#` 头会写明。

这里用的是**原子词汇**（`P_3` 而不是 `P-Q3`）：桥联描述的是两个原子之间的连接，
而 Qn 命名的是一个含多个原子的结构单元。与导出轨迹的标签一致。

### 用 pandas 分析 `linkage`

```python
import pandas as pd
d = pd.read_csv("network_linkage.csv", comment="#")

# 按元素对汇总
d.groupby(["elem_a", "elem_b"])["count"].sum()

# Al-O-Al 发生在哪种配位的 Al 之间 —— 与 NMR 的 Al[4]/Al[5]/Al[6] 直接对照
al = d.query("elem_a == 'Al' and elem_b == 'Al' and n_formers == 2")
al.pivot_table(index="cn_a", columns="cn_b", values="count", aggfunc="sum")

# P-O-P 的 Qn–Qn 连接矩阵
d.query("elem_a == 'P' and elem_b == 'P'").pivot_table(
    index="n_bridge_a", columns="n_bridge_b", values="count", aggfunc="sum")

# 某种配位的 Al 更爱连哪种 Qn 的 P
d.query("elem_a == 'Al' and elem_b == 'P'").pivot_table(
    index="cn_a", columns="n_bridge_b", values="count", aggfunc="sum")

# 只看真桥，排除三配位配体
d.query("n_formers == 2")

# 多配体体系：O 桥与 F 桥分开看
d.groupby("ligand")["count"].sum()
```

> **交叉核对。** `linkage` 里 `Al-O-Al` 的 `n_formers=2` 计数应等于
> `ligand_type` 表 `O_b, Al, Al` 那一行；差额来自三配位氧的贡献。
> $\sum(m_\mathrm{Al} \times \text{count})$（取自 `qn_partner`）应等于 `linkage`
> 里 `Al-O-P` 的真桥计数。两条都不对时，先查截断是不是给漏了。

---

## 标注轨迹（`--export-traj`）

导出每个输入一个文件：`<输入 stem>_types[_<后缀>].<ext>`。文件名必须掺输入 stem，
否则批处理时第二个输入会覆盖第一个（与 `ferro map` 的 cube 同理）。

### 标签

| 角色 | 标签 | 数字是 |
|---|---|---|
| Qn 形成子 | `P_0` `P_1` … `Si_4` | **Qn**（桥接配体数） |
| 其他形成子 | `Al_4` `Al_5` `Al_6` | **配位数** |
| 自由配体 | `O_f` | — |
| 非桥配体 | `O_n` | — |
| 桥联配体 | `O_b` | — |
| 三配位配体 | `O_t` | — |
| 修饰子 | `Zn` | 无后缀 |

运行 `ferro net` 时这张表会按本次参数填好元素打印一次，因为哪个元素报 Qn、哪个报
配位数取决于 `--qn` 和给了哪些截断。

### 两套词汇：单元与原子

| 词汇 | 形如 | 用在 | 为什么 |
|---|---|---|---|
| **单元** | `P-Q2`、`Al_4` | `composition` / `qn` / `qn_partner` | 这三张表数的是「有多少个 Q2 单元」，Qn 本就是结构单元的记号 |
| **原子** | `P_2`、`Al_4` | `linkage`、导出轨迹 | 桥联描述的是**原子之间**的连接，而一个 Qn 单元含多个原子 |

非 Qn 形成子两套词汇恰好相同（都是 `Al_4`），因为 Al 没有单元记号可用。

> **下游选型认原子词汇。** 导出轨迹里那个 P 叫 `P_2`，不叫 `P-Q2`——轨迹标签必须能
> 被 dump reader 按首个下划线拆回 element，而 `Q2` 里的 `Q` 不是元素。所以是
> `traj gr -x P_2`，不是 `-x P-Q2`。

带元素前缀（`P-Q2` 而非 `Q2`）是因为双 Qn 形成子体系（硼磷酸盐 B+P、铝硅酸盐
Si+Al）里 `Q0`…`Q4` 会各出现两遍，而**跨元素重名没有任何相邻列能解**。

全部合 `<元素>_<后缀>` 约定，按**首个下划线**拆分。修饰子不带角色后缀：旧方案
按非桥配体数 0/1/2/≥3 分档，而修饰子的实际配位数在 3–6，实测参考轨迹 97% 落进
兜底桶，没有分辨力。

> **同样长相的数字，含义按元素而变。** 这是文献自己的读法（$Q^2$ 与 Al[4] 没人会
> 混），代价是在 Al 不带非桥氧的体系里两种口径给出同一个数字，错了也看不出来——所以
> 约定在分类时就写进类型里，而不是渲染时查表，`--qn` 换了名单，标签跟着换。

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

# 五配位 Al 周围的桥氧 —— Al-O-Al 定位研究的入口
ferro traj gr -i run_types.lammpstrj -x Al_5 -y O_b --last-n 1
```

> **按标签选型只在单帧成立。** `g(r)` 要求逐类型的粒子数守恒，而标签是动态的——
> 一个 P 在轨迹推进中会改变自己的 Qn（实测参考轨迹 `P_3` 逐帧为
> 149 / 152 / 150 / 150 / 150）。对多帧标注轨迹按标签选型会被
> 「per-type atom counts change」错误拒绝，这是守卫而非缺陷。多帧请按元素选。

标注轨迹与统计表是互相独立的产物，可以互相验证：

- `traj gr -x P_3 -y O_b --last-n 1` 的第一壳层 CN 应**恰为 3**——分类层说「这个 P
  有 3 个桥氧」，独立的 CN 积分路径应当数回 3。
- `traj gr -x Al_5 -y O_b --last-n 1` 实测参考轨迹得 4.933 而非 5.000：差的 0.067
  是一个三配位氧，它标 `O_t` 而不是 `O_b`。对得上，说明标签的数字确实是配位数。

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
