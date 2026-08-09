# ferro net 开发手册（原 fe-network）

## 模块位置

```
ferro-analysis/src/network/
├── mod.rs          ← 公开 API、数据结构、顶层入口、跨帧累加器
├── cn.rs           ← 邻居表构建 + CN 统计
├── ligand_class.rs ← 配体角色分类（FO/NBO/BO/OBO）
├── qn.rs           ← Qn 物种计算
└── modifier.rs     ← 修饰子角色分类（Free/T/B/M）

ferro-cli/src/bin/network.rs
                    ← CLI 入口（参数解析、路由、CSV/XLSX 输出）
```

---

## 核心数据结构

### `NetworkParams`

```rust
pub struct NetworkParams {
    // 形成子-配体截断表，键 (former_elem, ligand_elem) → 截断半径 [Å]
    pub cutoffs: CutoffTable,
    // 修饰子-配体截断表，键 (modifier_elem, ligand_elem) → 截断半径 [Å]
    pub modifier_cutoffs: CutoffTable,
}

pub type CutoffTable = BTreeMap<(String, String), f64>;
```

使用 `BTreeMap` 而非 `HashMap`，保证遍历顺序确定、输出可复现。

`NetworkParams` 的辅助方法：

| 方法 | 作用 |
|---|---|
| `formers()` | 从 `cutoffs` 提取所有形成子元素（排序） |
| `ligands()` | 从 `cutoffs` 提取所有配体元素（排序） |
| `modifiers()` | 从 `modifier_cutoffs` 提取所有修饰子元素（排序） |
| `cutoff(former, ligand)` | 查询单个形成子-配体截断 |
| `formers_for_ligand(ligand)` | 给定配体，返回所有 `(former, cutoff)` 对 |

---

### `NetworkResult`

所有帧累加后归一化的最终结果：

```rust
pub struct NetworkResult {
    // (former, ligand) → [(cn_value, count, fraction)]，按 cn_value 升序
    pub cn_dist: HashMap<(String, String), Vec<(u32, usize, f64)>>,
    // former → [(cn_value, count, fraction)]（跨所有配体的总 CN）
    pub cn_total: HashMap<String, Vec<(u32, usize, f64)>>,
    // (former, ligand) → 平均 CN
    pub mean_cn: HashMap<(String, String), f64>,
    // former → 平均总 CN
    pub mean_cn_total: HashMap<String, f64>,
    // ligand → [(class_label, count, fraction)]，按 FO < NBO < BO < OBO 排序
    pub ligand_classes: HashMap<String, Vec<(String, usize, f64)>>,
    // former → [(qn_label, count, fraction)]，按 Q0 < Q1 < ... 排序
    pub qn_species: HashMap<String, Vec<(String, usize, f64)>>,
    // modifier_elem → [(role_label, count, fraction)]，按 Free < T < B < M 排序
    pub modifier_roles: HashMap<String, Vec<(String, usize, f64)>>,
}
```

---

## 数据流（完整链路）

```
CLI 参数解析
  │  split_pair_args()  ─→  pair args (--X-O=cutoff)
  │  clap                 ─→  普通参数 (--modifier, --last-n, ...)
  │
  ├─ parse_pairs()         ─→  BTreeMap<(String,String), f64>
  ├─ modifier_elems set    ─→  路由到 former_cutoffs / modifier_cutoffs
  └─ NetworkParams::new
       │
       ▼
calc_network(&traj, &params)
  │
  ├─ 并行：traj.frames.par_iter()
  │    └─ compute_frame(frame, cell, params)
  │         ├─ cn::build_neighbor_map()        → neighbor_map
  │         ├─ ligand_class::build_ligand_nf_map() → ligand_nf_map
  │         ├─ cn::collect_cn()                → cn_by_pair
  │         ├─ cn::collect_cn_total()          → cn_total_by_former
  │         ├─ ligand_class::classify_ligands() → ligand_labels
  │         ├─ qn::calc_qn()                   → qn_labels
  │         └─ modifier::classify_modifiers()  → modifier_labels
  │
  ├─ fold()  → 每个 rayon 线程内的 Accumulator
  └─ reduce() → 合并所有 Accumulator
       └─ Accumulator::finalize() → NetworkResult

NetworkResult
  └─ write_csv() / write_xlsx()
```

---

## 子模块详解

### `cn.rs` — 邻居表与配位数

#### `build_neighbor_map`

**输入**：`Frame`，`Cell`，`NetworkParams`  
**输出**：`HashMap<(former, ligand), Vec<Vec<usize>>>`

外层 `Vec` 对应该 former 元素的原子列表（与帧内 `frame.atoms` 中该元素出现顺序一致），内层 `Vec` 为该原子的配体邻居全局索引。

算法：
1. 按元素对 `frame.atoms` 建立索引 `elem_atoms: HashMap<&str, Vec<usize>>`
2. 对每个 `(former, ligand)` 对，双重遍历 former 原子 × ligand 原子
3. 用 `cell.minimum_image()` 计算最小镜像距离，与 `cutoff²` 比较

时间复杂度：当前 O(N_former × N_ligand) 每帧。若帧内原子数极大可改为 cell-list（参考 `md/gr.rs` 的 `CellList` 实现）。

#### `collect_cn` / `collect_cn_total`

纯统计，将邻居表的 `Vec` 长度转为 CN 值，再对多配体类型求和得总 CN。

---

### `ligand_class.rs` — 配体分类

#### `build_ligand_nf_map`

**输入**：`Frame`，`Cell`，`NetworkParams`  
**输出**：`HashMap<usize, Vec<(String, usize)>>`  
（配体原子全局索引 → 其网络形成子邻居列表，每项含形成子元素名和全局索引）

逻辑与 `build_neighbor_map` 对称——从配体视角遍历所有 former。

#### `classify_ligands`

对每个配体原子调用 `make_label(nf_neighbors)`：

```
nf_neighbors 长度 → 标签
0              → "FO"
1              → "NBO(X)"              X = former 元素
2              → "BO(X-Y)"            X, Y 字母序排列
≥3             → "OBO(X-Y-Z-...)"    字母序
```

---

### `qn.rs` — Qn 物种

#### `calc_qn`

对每个 former 原子 $i$：
1. 收集 $i$ 的所有配体邻居索引（跨所有配体类型）
2. 对每个配体邻居查 `ligand_nf_map`，调用 `make_label` 判断是否为 BO/OBO
3. 统计桥接配体数 `n_bridging` 和异核 former 的出现次数 `hetero_counts`
4. 调用 `make_qn_label(n_bridging, &hetero_counts)`

**异核计数细节**：一个桥接配体对每种异元素 former 只计一次（用 `HashSet<&str> seen` 去重），避免 OBO 氧有多个同种 former 邻居时重复计数。

#### `make_qn_label`

```
n=0, hetero 空         → "Q0"
n=3, hetero 空         → "Q3"
n=3, hetero {Al→2}    → "Q3(2Al)"
n=4, hetero {Al→1, Si→1} → "Q4(1Al,1Si)"   ← 字母序
```

#### `qn_label_order`

排序键：取标签第二个字符（数字）作为 Q 值，未解析的标签返回 99。

---

### `modifier.rs` — 修饰子角色

#### `build_nbo_set`

从 `ligand_nf_map` 中提取所有恰好有 **1 个** NF 邻居的配体原子索引，组成 `HashSet<usize>`。这些原子即 NBO（Non-Bridging Oxygen/Anion）。

```rust
ligand_nf_map.iter()
    .filter(|(_, nf)| nf.len() == 1)
    .map(|(&idx, _)| idx)
    .collect()
```

#### `classify_modifiers`

对帧内每个修饰子原子 $m$：
1. 遍历 NBO 集合，计算 $m$ 到每个 NBO 的最小镜像距离
2. 小于 `modifier_cutoff` 的计入 `nbo_count`
3. 按 `nbo_count` → `role_label()`

```
0 → "Free"
1 → "T"
2 → "B"
≥3 → "M"
```

**与 `build_neighbor_map` 的区别**：修饰子分析只针对 NBO 子集，而不是全部配体原子；因此不复用 `build_neighbor_map` 的结果。

---

### `mod.rs` — 跨帧累加器 `Accumulator`

`Accumulator` 在 rayon `fold`/`reduce` 中使用，负责将逐帧的 `FrameData`（原子级标签列表）聚合为全局计数：

```
FrameData.cn_by_pair        HashMap<pair, Vec<u32>>
  └─ push → Accumulator.cn_pair    HashMap<pair, HashMap<u32, usize>>

FrameData.modifier_labels   HashMap<modifier_elem, Vec<String>>
  └─ push → Accumulator.modifier   HashMap<modifier_elem, HashMap<label, usize>>
```

`merge(other: Self)` 合并两个 `Accumulator`（rayon reduce 的归约操作），对应字段逐一合并计数。

`finalize()` 将 `HashMap<key, usize>` 转为有序的 `Vec<(label, count, fraction)>`，分别按各自的排序键排列：

| 字段 | 排序规则 |
|---|---|
| `cn_dist` / `cn_total` | CN 值升序 |
| `ligand_classes` | `ligand_class_order()`: FO(0) < NBO(1) < BO(2) < OBO(3) |
| `qn_species` | `qn_label_order()`: Q0 < Q1 < ... |
| `modifier_roles` | `modifier_role_order()`: Free(0) < T(1) < B(2) < M(3) |

---

## CLI 参数路由逻辑

```rust
// 1. 将 argv 拆分为 pair 参数 vs. 普通参数
let (pair_args, clap_args) = split_pair_args(&all_args);
// 判据: "--X..." 其中 X 首字母大写且含 "="

// 2. 解析 pair 参数 → 统一的大表
let all_pairs: BTreeMap<(String,String), f64> = parse_pairs(&pair_args)?;
// "--P-O=2.3" → ("P", "O") → 2.3

// 3. 从 --modifier Zn,Na 构建修饰子集合
let modifier_elems: HashSet<String> = ...;

// 4. 按修饰子集合路由
for ((elem, ligand), cutoff) in all_pairs {
    if modifier_elems.contains(&elem) {
        modifier_cutoffs.insert((elem, ligand), cutoff);
    } else {
        cutoffs.insert((elem, ligand), cutoff);  // 视为形成子
    }
}
```

不在 `--modifier` 中的元素自动视为网络形成子，保持向后兼容。

---

## 并行策略

```
traj.frames.par_iter()              ← rayon 帧级并行
  .filter_map(|frame| compute_frame(...))
  .fold(|| Accumulator::new(...), |mut acc, fd| { acc.push(&fd); acc })
  .reduce(|| Accumulator::new(...), |mut a, b| { a.merge(b); a })
```

- `fold`：每个 rayon 线程内本地累加，避免锁竞争
- `reduce`：线程间合并 `Accumulator`（需 `Accumulator::merge` 语义正确）
- `compute_frame` 不可变借用 `frame`，完全无锁

---

## 如何扩展

### 添加新的配体分类规则

修改 `ligand_class.rs::make_label()`，无需改动其他模块。

### 添加新的形成子统计量

1. 在 `FrameData` 中增加字段
2. `compute_frame()` 中填充
3. `Accumulator` 中增加对应计数字段（`push` + `merge`）
4. `finalize()` 中生成 `NetworkResult` 新字段

### 支持修饰子-修饰子 CN

目前修饰子只做 NBO 角色分类，不计算修饰子之间或修饰子与配体的 CN。
若需要，在 `modifier_cutoffs` 中可定义修饰子-配体对，然后在 `compute_frame` 中对修饰子调用 `build_neighbor_map` 的逻辑，将结果单独存储。

### 切换为 Cell-list 加速

`build_neighbor_map` 和 `classify_modifiers` 当前均为暴力 O(N²) 遍历。
若体系超过 ~5000 原子，可复用 `ferro-analysis/src/md/angle.rs` 中的 `CellList`（`pub(super)` 可见），将查找复杂度降至 O(N)。

---

## 错误处理

- 空轨迹 / 任意帧缺少 Cell → `calc_network` 返回 `None`，CLI 层转为 `anyhow::Error`
- 无 pair 参数 → CLI 提前 `bail!`
- `cell.minimum_image()` 返回 `Err` 的情况（奇异晶胞）→ 计算闭包中 `.expect("cell is non-singular")`（已在 `calc_network` 入口处验证 Cell 存在）
