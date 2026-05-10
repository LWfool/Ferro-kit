# Glass Network Analysis (`fe-network`)

`fe-network` 对 MD 轨迹逐帧建立玻璃网络拓扑，统计以下四类量并在全轨迹上累加归一化：

| 量 | 含义 |
|---|---|
| 配位数 (CN) | 网络形成子与各配体类型的平均配位数及其分布 |
| 配体分类 | FO / NBO / BO / OBO |
| Qn 物种 | 网络形成子的桥接度分布 |
| 修饰子角色 | Free / T / B / M（可选） |

---

## 理论背景

### 键连关系的判定

以截断半径（cutoff）作为成键判据：两原子间距小于截断值时视为成键。距离计算使用最小镜像约定（minimum-image convention），支持斜方和三斜晶胞：

$$d_{ij} = \bigl|\mathbf{r}_{ij} - \mathbf{M} \cdot \text{round}\!\left(\mathbf{M}^{-1}\mathbf{r}_{ij}\right)\bigr|$$

分析要求每帧均存在晶胞（PBC）。

### 配位数 (CN)

对每个网络形成子原子 $i$（元素 $\alpha$），其与配体元素 $\beta$ 的配位数为：

$$\text{CN}^{(\alpha\beta)}_i = \bigl|\{j : \text{elem}(j)=\beta,\; d_{ij} < r_{\alpha\beta}^{\text{cut}}\}\bigr|$$

全轨迹统计：对所有帧、所有同类形成子原子的 CN 值做频率直方图，给出 CN 分布和均值。

**总 CN**：对同一形成子原子，将其对所有配体类型的 CN 求和。

### 配体分类

对每个配体原子 $k$（元素 $\beta$），统计其网络形成子邻居数 $n_k$：

| $n_k$ | 标签 | 含义 |
|---|---|---|
| 0 | `FO` | Free oxygen（游离氧） |
| 1 | `NBO(X)` | Non-bridging oxygen（非桥氧），X = 形成子元素 |
| 2 | `BO(X-Y)` | Bridging oxygen（桥氧），X、Y 字母序 |
| ≥3 | `OBO(X-Y-Z)` | Over-bridging oxygen（过桥氧） |

混合形成子体系（如 P + Si）时，标签携带具体元素信息，例如 `BO(P-Si)` 表示该桥氧连接一个 P 和一个 Si。

### Qn 物种分类

对每个网络形成子原子 $i$，$Q_n$ 表示其桥接配体数（即所连配体中被分类为 BO 或 OBO 的数量）：

$$n_i = \bigl|\{k \in \text{neighbors}(i) : \text{label}(k) \in \{\text{BO}, \text{OBO}\}\}\bigr|$$

**异核桥接标注**：若桥接配体 $k$ 的另一侧形成子与 $i$ 不同元素，则在标签中注明，例如：

- `Q3` — 3 个桥氧，全为同核（P-O-P）
- `Q3(2Al)` — 3 个桥氧，其中 2 个为 P-O-Al 异核桥
- `Q4(1Al,1Si)` — 4 个桥氧，分别含 1 个 P-O-Al 和 1 个 P-O-Si

异核元素按字母序排列，格式为 `Q{n}({count}{Elem},...)`.

### 修饰子角色分类

修饰子（modifier）阳离子（如 Zn²⁺、Na⁺）通过与非桥氧（NBO）配位锚定在网络中。对每个修饰子原子 $m$，统计其在截断半径内的 NBO 邻居数 $p_m$：

| $p_m$ | 角色标签 | 物理含义 |
|---|---|---|
| 0 | `Free` | 未与任何 NBO 配位 |
| 1 | `T` | Terminal（单端连接一个 NBO） |
| 2 | `B` | Bridging（桥接两个 NBO） |
| ≥3 | `M` | Multi（连接 3 个或更多 NBO） |

---

## 快速入门

```bash
# P₂O₅ 玻璃，P-O 截断 2.3 Å，输出 CSV
fe-network -i traj.dump --P-O=2.3

# 混合 P-Si 玻璃，多个配体对
fe-network -i traj.dump --P-O=2.3 --Si-O=1.8 -o result.csv

# 含 Zn 修饰子，分析其角色
fe-network -i traj.dump --P-O=2.3 --Zn-O=3.5 --modifier Zn

# 多修饰子
fe-network -i traj.dump --P-O=2.3 --Zn-O=3.5 --Na-O=3.2 --modifier Zn,Na

# 输出 xlsx（多 sheet）
fe-network -i traj.dump --P-O=2.3 --format xlsx -o result.xlsx

# 只统计尾部 500 帧，8 线程
fe-network -i traj.dump --P-O=2.3 --last-n 500 --ncore 8
```

---

## 命令行参数

### 必需

`fe-network` 至少需要一个 `--Former-Ligand=cutoff` 形式的配对参数（首字母大写）。

```
--P-O=2.3        P（形成子）与 O（配体）的截断半径 2.3 Å
--Si-O=1.8       Si-O 截断 1.8 Å
--Al-O=2.0       Al-O 截断 2.0 Å
--Al-F=2.1       Al-F 截断 2.1 Å（同一形成子可有多种配体）
```

### 可选参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-i <file>` | — | 输入轨迹文件（缺省时显示帮助） |
| `-o <file>` | `network` | 输出文件前缀 |
| `--format csv\|xlsx` | `csv` | 输出格式 |
| `--last-n N` | 全部帧 | 仅使用尾部 N 帧 |
| `--ncore N` | 全部核心 | 并行线程数 |
| `--metal-units` | 关 | LAMMPS metal 单位（速度 Å/ps，力 eV/Å） |
| `--modifier Elem` | — | 修饰子元素（逗号分隔，如 `Zn` 或 `Zn,Na`） |

**修饰子截断**：通过相同的 `--Elem-O=cutoff` 语法提供，系统根据 `--modifier` 列表自动路由。

---

## 输出文件

### CSV 格式（默认）

| 文件 | 内容 |
|---|---|
| `<stem>_cn.csv` | 各 (former, ligand) 配对的 CN 分布及均值，含总 CN |
| `<stem>_ligand.csv` | 各配体元素的 FO/NBO/BO/OBO 分布 |
| `<stem>_qn.csv` | 各形成子元素的 Qn 物种分布 |
| `<stem>_modifier.csv` | 各修饰子元素的角色分布（仅指定 `--modifier` 时生成） |

#### `_cn.csv` 格式

```csv
Former,Ligand,CN,Count,Fraction
P,O,3,120,0.240000
P,O,4,380,0.760000
P,O,mean,3.7600,
P,total,4,500,1.000000
P,total,mean,3.7600,
```

#### `_ligand.csv` 格式

```csv
Ligand,Class,Count,Fraction
O,FO,50,0.033333
O,NBO(P),650,0.433333
O,BO(P-P),800,0.533333
```

#### `_qn.csv` 格式

```csv
Former,Species,Count,Fraction
P,Q2,120,0.240000
P,Q3,380,0.760000
```

#### `_modifier.csv` 格式

```csv
Modifier,Role,Count,Fraction
Zn,Free,10,0.040000
Zn,T,120,0.480000
Zn,B,110,0.440000
Zn,M,10,0.040000
```

### XLSX 格式

使用 `--format xlsx` 时，所有结果写入同一文件的多个 sheet：

| Sheet | 内容 |
|---|---|
| CN | 配位数分布 |
| Ligand | 配体分类分布 |
| Qn | Qn 物种分布 |
| Modifier | 修饰子角色分布（有修饰子时才生成） |

---

## 典型参数参考

| 体系 | 形成子 | 配体 | 参考截断 |
|---|---|---|---|
| P₂O₅ 玻璃 | P | O | P-O: 2.3 Å |
| SiO₂ 玻璃 | Si | O | Si-O: 1.8 Å |
| Al₂O₃ | Al | O | Al-O: 2.1 Å |
| GeO₂ 玻璃 | Ge | O | Ge-O: 2.0 Å |
| ZnO-P₂O₅ 玻璃 | P | O | P-O: 2.3 Å，Zn（modifier）: Zn-O: 3.5 Å |

截断值应参照 g(r) 第一峰谷位置（`fe-traj -m gr`）确定。
