# Bader 电荷分析 — 开发文档

参考实现：Henkelman 组 Fortran 源码（`examples/bader/`）

参考文献：
- Henkelman et al., *Comput. Mater. Sci.* **36**, 254 (2006) — on-grid / near-grid
- Sanville et al., *J. Comput. Chem.* **28**, 899 (2007) — near-grid 精化
- Tang et al., *J. Phys.: Condens. Matter* **21**, 084204 (2009) — off-grid
- Yu & Trinkle, *J. Chem. Phys.* **134**, 064111 (2011) — weight 方法

---

## 1. 理论基础

### 核心思想

Bader 分析基于电荷密度拓扑：空间中每个体积元被分配给最近的电荷密度极大值点，从而把空间划分为若干 **Bader 体积**（Bader volume）。每个体积对应一个原子，对体积内的电荷积分即得该原子的 Bader 电荷。

**Bader 体积边界**（零通量面）满足：
$$\nabla\rho(\mathbf{r}) \cdot \hat{n}(\mathbf{r}) = 0$$
即法向量方向上的电荷密度梯度为零。

### 算法基本流程

```
对每个网格点 p：
  1. 从 p 出发沿梯度方向爬升到局部极大值
  2. 记录极大值所属体积编号 volnum[p]
  
对每个 Bader 体积 v：
  3. 对体积内所有点的电荷密度求和 → volchg[v]
  
对每个 Bader 体积 v：
  4. 找到距极大值最近的原子 j → nnion[v] = j
  5. ionchg[j] += volchg[v]
```

---

## 2. 输入格式：VASP CHGCAR

### 文件结构

```
系统名称（任意字符串）
缩放因子 (scalefactor)
晶格矩阵 3×3（行：a, b, c，单位 Å × scalefactor）
离子种类名称（VASP5 格式）
各种离子数量（如 "8   4   12"）
Direct
离子分数坐标（每行一个原子）
（空行）
网格维度 N1 N2 N3
电荷密度数据（x/N1 变化最快，z/N3 最慢）
```

### 密度归一化约定（关键）

VASP 存储的值是物理密度乘以晶胞体积：
$$\rho_{\text{stored}}(\mathbf{r}) = \rho_{\text{physical}}(\mathbf{r}) \times V_{\text{cell}}$$

单位：`e`（无量纲，已乘了体积因子）。

**读取时不做除法，直接存入 `chg.rho`。** 后续各处的使用约定：

| 用途 | 公式 | 说明 |
|---|---|---|
| 真空判断 | `rho_stored / V_cell <= vacval` | 还原为 e/Å³ |
| 电荷积分 | `Σ rho_stored / nrho` | 得到体积总电子数（e） |
| 梯度上升 | 直接比较 `rho_stored` | 比例关系不变 |

推导：`total_e = Σ ρ_phys * dV = Σ (ρ_stored/V_cell) * (V_cell/nrho) = Σ ρ_stored / nrho` ✓

### CHGCAR 读取时计算的辅助量

读取完成后立即初始化以下变换矩阵（`chgcar_mod.f90` L95-128）：

```
lat2car[:,i] = dir2car[:,i] / N_i       # 格点→笛卡尔（列向量 = a_i/N_i）
car2lat = inverse(lat2car)              # 笛卡尔→格点

lat_dist[d1,d2,d3] = |lat2car @ [d1,d2,d3]|    # 26邻居的实际距离（Å）
lat_i_dist[d1,d2,d3] = 1 / lat_dist            # 距离倒数（(0,0,0)处=0）
```

---

## 3. 核心数据结构

### Fortran → Rust 映射

#### `charge_obj` → `ChargeGrid`

```rust
pub struct ChargeGrid {
    // 存储 rho_stored = rho_physical × V_cell，x 变化最快
    pub rho: Vec<f64>,           // len = nrho = N1 * N2 * N3
    pub shape: [usize; 3],       // [N1, N2, N3]
    pub nrho: usize,             // N1 * N2 * N3

    // 坐标变换（列 i = a_i/N_i，即每个格步的笛卡尔向量）
    pub lat2car: Matrix3<f64>,
    pub car2lat: Matrix3<f64>,   // = lat2car.try_inverse()

    // 26 邻居预计算距离（d1,d2,d3 ∈ {-1,0,1}，存为 [3][3][3]）
    pub lat_dist:   [[[f64; 3]; 3]; 3],
    pub lat_i_dist: [[[f64; 3]; 3]; 3],  // (0,0,0) 处为 0.0
}
```

索引约定：`rho[n1 + N1*(n2 + N2*n3)]`（x 最快，与 CHGCAR 文件顺序一致）

#### `bader_obj` → `BaderResult`

```rust
pub struct BaderResult {
    pub nvols: usize,

    // 每个格点所属体积编号（0=未分配，-1=真空临时，>0=体积编号）
    // 分析结束后 -1 → nvols+1
    pub volnum: Vec<i32>,        // len = nrho，x 最快

    // 梯度上升路径跟踪（仅分析过程中使用）
    pub known: Vec<u8>,          // 0=未知, 1=当前路径在此点, 2=已完全确定

    // 各 Bader 体积的极大值格点坐标（1-indexed）
    pub volpos_lat: Vec<[f64; 3]>,  // len = nvols
    pub volpos_dir: Vec<[f64; 3]>,  // 分数坐标
    pub volpos_car: Vec<[f64; 3]>,  // 笛卡尔坐标（Å）

    // 各体积积分电荷（e）
    pub volchg: Vec<f64>,        // len = nvols

    // 各原子结果
    pub ionchg:     Vec<f64>,    // len = nions，各原子 Bader 电荷（e）
    pub iondist:    Vec<f64>,    // len = nvols，极大值→最近原子距离（Å）
    pub nnion:      Vec<usize>,  // len = nvols，各体积对应原子编号（0-indexed）
    pub ionvol:     Vec<f64>,    // len = nions，各原子 Bader 体积（Å³）

    pub vacchg: f64,
    pub vacvol: f64,

    // 分析过程内部状态（不对外暴露）
    pub(crate) tol: f64,         // 体积电荷截断阈值（默认 1e-3）
    pub(crate) stepsize: f64,    // off-grid 步长（Å），默认 min(voxel尺寸)
}
```

#### `ions_obj` → 直接使用 `Frame`

| Fortran 字段 | Ferro 对应 |
|---|---|
| `r_car[i,:]` | `frame.atoms[i].position` |
| `r_dir[i,:]` | `frame.cell.unwrap().cartesian_to_fractional(pos)` |
| `r_lat[i,:]` | `car2lat @ r_car + org_lat`（格点坐标，用于 ion → 格点）|
| `lattice` | `frame.cell.unwrap().matrix`（行向量 a, b, c）|
| `dir2car` | `cell.matrix.transpose()` |
| `nions` | `frame.atoms.len()` |

---

## 4. 算法方法详解

所有方法的主流程相同：

```
for p in all_grid_points:
    if volnum[p] == 0:
        maximize(bdr, chg, p)          // 从 p 爬升到极大值
        path_vol = volnum[p_final]
        if path_vol == 0:              // 新极大值
            nvols += 1
            path_vol = nvols
            volpos_lat[nvols] = p_final
        for q in path:
            if volnum[q] != -1:        // 非真空
                volnum[q] = path_vol
            known[q] = 0              // 重置路径标记（known=2 保留）
```

### 4.1 On-grid 方法

**原理**：在 26 个邻居中找到"梯度校正密度"最大的方向移动。

**距离校正**：非正交/非均匀格点中邻居距离不同。直接比较密度值会偏向近邻，引入各向异性偏差。校正方式——将密度外推到单位格步距离处：

$$\rho_{\text{tmp}} = \rho_{\text{ctr}} + (\rho_{\text{nbr}} - \rho_{\text{ctr}}) \times \text{lat\_i\_dist}[d_1,d_2,d_3]$$

即 $\rho_{\text{ctr}} + \Delta\rho / d$，相当于在中心到邻居方向上单位距离处的外推密度。移动到使该值最大的邻居。

**Rust 实现** (`step_ongrid`)：

```rust
fn step_ongrid(chg: &ChargeGrid, p: &mut [usize; 3]) {
    let rho_ctr = rho_val(chg, p);
    let mut rho_max = rho_ctr;
    let mut pm = *p;
    for d1 in -1i32..=1 {
        for d2 in -1i32..=1 {
            for d3 in -1i32..=1 {
                let pt = [p[0] as i32 + d1, p[1] as i32 + d2, p[2] as i32 + d3];
                let rho_nbr = rho_val_i(chg, pt);
                let li = lat_i_dist(chg, d1, d2, d3);
                let rho_tmp = rho_ctr + (rho_nbr - rho_ctr) * li;
                if rho_tmp > rho_max {
                    rho_max = rho_tmp;
                    pm = pbc_i(chg, pt);
                }
            }
        }
    }
    *p = pm;
    // 若 p 未变，说明到达极大值
}
```

**终止**：`p` 不再变化（所有邻居校正密度 ≤ 自身）。

### 4.2 Near-grid 方法（默认）

**原理**：中心差分计算梯度方向，用浮点累积小数位移 `dr`，当累积量超过 0.5 时额外整步移动。效果：在格点网格上近似连续梯度上升，避免 on-grid 的直线偏差。

**关键细节：`dr` 是跨调用的持久状态**（Fortran `SAVE`），在每条新路径开始时（`pnum == 1`）归零，路径中持续累积。Rust 实现必须将 `dr` 存入路径状态结构体，不能作为局部变量。

**`step_neargrid` 完整逻辑**（`bader_mod.f90` L454-498）：

```
输入：p（当前格点）、dr（持久累积量，路径开始时=0）
输出：p（新格点，若不变说明到达极大值）

1. 计算梯度方向 grad = rho_grad_dir(chg, p)
   （中心差分，见 §4.5）

2. 若 max|grad| < 1e-30:
   → 调用 is_max 判断是否真极大值
      是：dr=0, return（未移动 = 极大值）
      否：fallback → step_ongrid(chg, pm)

3. 归一化：coeff = 1 / max(|grad_i|)
   gradrl = coeff * grad

4. 整数步：pm = p + round(gradrl)
5. 小数累积：dr += gradrl - round(gradrl)
6. 若 |dr_i| >= 0.5：pm += round(dr)；dr -= round(dr)

7. known[p] = 1  （标记当前点在路径上）

8. pbc(pm, npts)
9. 若 known[pm] == 1:
   → 检测到循环，fallback → step_ongrid(chg, pm)；dr = 0

10. p = pm
```

**`max_neargrid` 的两种退出条件**：

| 条件 | 原因 | 处理 |
|---|---|---|
| `p` 不变（step 后 p == path[pnum]）| 到达极大值 | EXIT，p 为新极大值 |
| `known[p] == 2`（quit_known 模式）| 已完全确定体积 | EXIT，继承该点 volnum |

注意：`known == 1` 只在 `step_neargrid` 内部触发 fallback（循环检测），不在 `max_neargrid` 层面触发退出。

**`known = 2` 的设置**：通过 `assign_surrounding_pts` → `known_volnum_ongrid` 完成：当某点的 6 个面邻居全部属于同一体积时，标记该点 `known = 2`（已完全确定）。

### 4.3 Off-grid 方法

**原理**：三线性插值实现连续空间中的梯度上升，步长为 `stepsize`（默认为最小格步边长）。

**`step_offgrid`**：
```
1. 在浮点坐标 r 处三线性插值得到 (rho, grad_car)
2. 沿 grad_car 方向步进：dr_car = grad_car * stepsize / |grad_car|
3. dr_lat = car2lat @ dr_car
4. r += dr_lat
5. pm = to_lat(r)（找最近格点，非正交格需在 2×2×2 邻居中找最近的）
6. 若 rho(pm) < rho(p)：退出（密度下降）
```

**三线性插值** `rho_grad`（`charge_mod.f90` L98-168）——在浮点格点坐标 `r = p0 + (f1,f2,f3)` 处：

```
g_i = 1 - f_i

rho00_ = rho000*g3 + rho001*f3       # 沿第3轴插值
rho01_ = rho010*g3 + rho011*f3
rho10_ = rho100*g3 + rho101*f3
rho11_ = rho110*g3 + rho111*f3

rho0__ = rho00_*g2 + rho01_*f2       # 沿第2轴
rho1__ = rho10_*g2 + rho11_*f2

rho    = rho0__*g1 + rho1__*f1       # 沿第1轴 → 插值密度

# 梯度（格点坐标系中的偏导数）
grad_lat[0] = rho1__ - rho0__
grad_lat[1] = rho_1_ - rho_0_
grad_lat[2] = rho__1 - rho__0

# 笛卡尔梯度（行向量形式）
grad_car = grad_lat @ car2lat
       等价于 car2lat^T @ grad_lat（列向量形式，物理正确）
```

### 4.4 Weight 方法（Yu-Trinkle，最高精度）

**原理**：将格点按密度降序排列，从高到低流动概率分配电荷，无格点偏差（lattice bias）。

#### Step 1：Wigner-Seitz Voronoi 分解 (`ws_voronoi`)

构建 WS 胞的邻居向量集合 `{vect}` 及对应面积权重 `{alpha}`：

```
遍历所有 3^3-1=26 个近邻方向（范围可扩展到 nRange=3，共 (2*3+1)^3-1=342 个）
筛选：R/2 在 WS 胞内部的向量保留

对每个保留方向 n，计算面积因子 alpha[n]：
  找到 WS 胞面 n 的所有顶点（三个面的交点，需在 WS 胞内）
  按绕 R[n] 方向排序顶点
  alpha[n] = (1/4 * Rmag[n]) * Σ triple_product(v_i, v_{i+1}, R[n])
```

`alpha` 的物理意义：流量 `t[n] = alpha[n] * Δρ`，正比于面 `n` 的面积投影。

#### Step 2：排序

所有格点按 `rho_stored` 降序排列（密度相等时按索引次序）。记 `chgList[n]` 为第 `n` 大的格点，`indList[pos] = n`（位置→排名反查表）。

#### Step 3：流分配（高密度→低密度遍历）

```rust
for n in 0..npts {
    // 找所有密度更高的邻居（排名 m < n）
    let mut above = vec![];
    let mut tsum = 0.0;
    for nv in 0..num_vect {
        let p = pbc(chgList[n].pos + vect[nv]);
        let m = indList[p];
        if m < n {   // m 密度更高
            let t = alpha[nv] * (chgList[m].rho - chgList[n].rho);
            above.push((m, t));
            tsum += t;
        }
    }

    if above.is_empty() {
        // 极大值 → 新 basin
        basin[n] = new_basin_id;
        volnum[chgList[n].pos] = basin[n];
        volpos_lat[basin_id] = chgList[n].pos;  // ← 用 chgList[n].pos，不是 p！
    } else {
        let tbasin = basin[above[0].0];
        let is_boundary = above.iter().any(|(m, _)| basin[*m] != tbasin || tbasin == 0);

        if is_boundary {
            // 边界点：按流比例分配给各上游 basin
            basin[n] = 0;
            for (m, t) in &above {
                prob[m][numbelow[m]] = t / tsum;
                neigh[m][numbelow[m]] = n;
                numbelow[m] += 1;
            }
            // volnum 取流量最大的上游 basin（用于可视化，不影响电荷计算）
        } else {
            // 内部点：继承上游 basin
            basin[n] = tbasin;
            volnum[chgList[n].pos] = tbasin;
        }
    }
}
```

#### Step 4：电荷积分

Fortran 原始实现是 O(N_vols × N_grid)（`weight_mod.f90` L215-239），效率较低。**Rust 推荐单遍 O(N_grid) 实现**：

```rust
// 权重数组 w[n] = 该点属于当前 basin 的权重（0~1）
// 按密度降序（n=0 是最高密度）遍历

// 初始化
let mut w = vec![0.0f64; npts];
let mut volchg = vec![0.0f64; nvols];

// 极大值点 w = 1；内部点从上游继承
for n in 0..npts {
    if basin[n] > 0 && w[n] == 0.0 {
        w[n] = if basin[n] == current_basin { 1.0 } else { 0.0 };
    }
    // 向下游传播权重
    for nb in 0..numbelow[n] {
        w[neigh[n][nb]] += prob[n][nb] * w[n];
    }
    if w[n] != 0.0 {
        let bv = basin[n];
        if bv > 0 { volchg[bv] += w[n] * chgList[n].rho; }
    }
}
// 注意：实际上需要对每个 basin 分别做一遍，或用更高效的矩阵方法
```

实际上 Fortran 的实现就是逐 basin 的，如果 N_ions 较少（通常 < 200）可以接受。Rust 可先保持相同结构，后期优化。

#### Step 5：真空点处理

```
分析结束后：volnum == -1 → nvols + 1
```

---

### 4.5 `rho_grad_dir` — Near-grid 梯度方向计算

**中心差分** + **一维极大抑制**（`charge_mod.f90` L176-218）：

```
grad_lat[i] = (rho(p + e_i) - rho(p - e_i)) / 2

// 若第 i 轴两侧都比当前点低（局部极大值的平台），清零该分量
if rho(p+e_i) < rho(p) && rho(p-e_i) < rho(p):
    grad_lat[i] = 0
```

**坐标变换**（`charge_mod.f90` L207-212）：

```fortran
rho_grad_car = MATMUL(rho_grad_lat, chg%car2lat)   ! 格点梯度 → 笛卡尔
rho_grad_dir = MATMUL(chg%car2lat, rho_grad_car)   ! 笛卡尔 → 再乘 car2lat
```

⚠️ **已知问题**：第二步对非正交晶胞理论上应用 `lat2car^T @ grad_car`，而非 `car2lat @ grad_car`。对于接近正交的晶胞（绝大多数体系）误差可忽略；对高度倾斜晶胞（如 α < 60°）建议改用 off-grid 或 weight 方法。Rust 实现可保持原始行为（复现 Fortran），在注释中标注这一近似。

---

## 5. 边缘精化（`refine_edge`）

网格方法在 Bader 体积边界处可能分配有误，需迭代精化。

**标记边缘点**（`is_vol_edge`）：检查 26 个邻居，若有邻居 volnum 与自身不同，则为边缘点。标记方式：`volnum[p] = -volnum[p]`（取负）。

**重新计算**：对所有 `volnum < 0` 的点重新执行梯度上升（`max_neargrid` 或 `max_ongrid`）。

**迭代策略**：
- `-r -1`（默认自动）：反复执行，直到 `num_reassign == 0`
- `-r N`（指定次数）：执行 N 次
- `-r -2`：仅精化第一次检测到的边缘点（不扩散）

---

## 6. 真空处理

密度低于阈值的格点标记为真空：

```
if |rho_stored[p]| / V_cell <= vacval:
    volnum[p] = -1
```

默认 `vacval = 1e-3` e/Å³。

分析结束时将 `-1` 替换为 `nvols + 1`（真空"体积"索引），电荷积分时跳过这些点。

---

## 7. 电荷积分与原子分配

### 体积电荷

```
volchg[v] = Σ_{p in vol_v} rho_stored[p] / nrho
```

其中 `nrho = N1 * N2 * N3`，结果单位为 e（总电子数贡献）。

### 原子分配（`assign_chg2atom`）

对每个 Bader 体积 `v`，找距极大值 `volpos_dir[v]` 最近的原子（考虑 PBC）：

```rust
// 用 Cell::minimum_image 处理周期性
for v in 0..nvols {
    let (ion_idx, dist) = frame.atoms.iter().enumerate()
        .map(|(j, atom)| {
            let dv = cell.minimum_image(volpos_car[v] - atom.position);
            (j, dv.norm())
        })
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap();
    nnion[v] = ion_idx;
    iondist[v] = dist;
    ionchg[ion_idx] += volchg[v];
}
```

### 原子体积

```
ionvol[j] = count(volnum == Bader_vols_of_j) * V_cell / nrho
```

---

## 9. 已知陷阱与修正

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| `ChargeGrid` 构造 | VASP 密度已乘 V_cell，不再除以体积 | 读入直接存，真空判断除以 vol，积分除以 nrho |
| `step_neargrid` | `dr` 是跨调用状态（Fortran SAVE） | 放入 `PathState.dr`，路径开始时归零 |
| Weight 极大值位置 | Fortran 存的是 `p`（最后邻居），不是极大值 | Rust 改为存 `chgList[n].pos`（极大值本身）|
| `known` 状态 | `known=1` 触发 step 内 fallback；`known=2` 触发 max 内 EXIT | 两种情况处理路径完全不同，不要混用 |
| `rho_grad_dir` | 对非正交晶胞梯度方向存在近似误差 | 高倾斜晶胞改用 weight 方法 |
| CHGCAR 读取 | 文件中密度循环顺序：n1 最快，n3 最慢 | Rust Vec 索引：`rho[n1 + N1*(n2 + N2*n3)]` |
| `lat_i_dist(0,0,0)` | 中心点距离为 0，不能取倒数 | 存为 `0.0`（Fortran 代码同样处理） |
| 体积电荷归一化 | `volchg = sum / nrho`，不是 `sum / vol` | 见 §2 归一化约定 |
| weight 边界点 volnum | 边界点 volnum 存的是"最大流"来源体积（近似） | 仅影响可视化；电荷积分仍用 w 权重，不依赖此 volnum |

---

## 10. 输出文件

| 文件 | 内容 |
|---|---|
| `ACF.dat` | 每个原子：坐标、Bader 电荷、最小表面距离、原子体积 |
| `BCF.dat` | 每个显著 Bader 体积：坐标、电荷、归属原子、距离 |
| `AVF.dat` | 每个原子对应的显著体积编号列表 |
| `BvIndex.dat` | 每个格点的 Bader 体积编号（可视化用） |
| `AtIndex.dat` | 每个格点的归属原子编号（可视化用） |

**ACF/BCF/AVF 三个 writer 不走 `Table`**（全项目唯一的例外）：它们是 Henkelman 组的
既定格式，外部工具按它解析。0.2.0 让其余 7 个 `write_*` 退役时，这三个明确保留。

已实现部分见 `progress.md`「ferro-analysis / dft」；`--outdir` 尚未接入（三个固定
文件名写在当前目录，连跑两个体系会互相覆盖），见 `plan.md` 的三项小待办。
