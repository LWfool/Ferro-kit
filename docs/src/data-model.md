# Core Data Model

## Type Hierarchy

Every ferro API accepts and returns `Trajectory`, even for single-frame files.

```
Trajectory
  └── Vec<Frame>
        ├── Vec<Atom>           — atomic species, positions, optional properties
        ├── Option<Cell>        — None = non-periodic
        ├── [bool; 3]           — pbc flags per axis
        └── Optional results    — energy, forces, stress, velocities
```

## Atom

```rust
pub struct Atom {
    pub element: String,
    pub position: Vector3<f64>,  // Å, Cartesian
    pub label: Option<String>,   // e.g. "Fe1", "Fe2"
    pub mass: Option<f64>,       // None → look up from element table
    pub magmom: Option<f64>,     // initial magnetic moment (DFT input)
    pub charge: Option<f64>,     // Bader / DDEC charge (post-processing)
}
```

The atom **index is implicit** (position in `Vec<Atom>`); no `index` field is stored to prevent inconsistency.

## Frame

```rust
pub struct Frame {
    pub atoms: Vec<Atom>,
    pub cell: Option<Cell>,              // None = non-periodic
    pub pbc: [bool; 3],
    pub charge: i32,
    pub multiplicity: u32,
    pub bonds: Option<Vec<(usize,usize)>>,
    pub energy: Option<f64>,             // eV
    pub forces: Option<Vec<Vector3<f64>>>, // eV/Å
    pub stress: Option<Matrix3<f64>>,    // eV/Å³
    pub velocities: Option<Vec<Vector3<f64>>>, // Å/fs (internal standard)
}
```

`pbc = [false,false,false]` → molecular system.  
`pbc = [true,true,false]` → surface / slab.  
`pbc = [true,true,true]` → bulk crystal or glass.

## Cell

```rust
pub struct Cell {
    pub matrix: Matrix3<f64>,  // row vectors a, b, c  [Å]
}
```

| Method | Description |
|---|---|
| `from_lengths_angles(a,b,c,α,β,γ)` | Construct from lattice parameters |
| `lengths() -> [f64; 3]` | \|a\|, \|b\|, \|c\| |
| `angles() -> [f64; 3]` | α, β, γ in degrees |
| `volume() -> f64` | Cell volume [Å³] |
| `fractional_to_cartesian(f)` | f·M |
| `cartesian_to_fractional(c)` | c·M⁻¹ |
| `wrap_position(c)` | Fold Cartesian position into [0, L) |
| `minimum_image(v)` | Apply minimum-image convention to vector v |

### Coordinate Conventions

The cell matrix stores row vectors:

$$\mathbf{M} = \begin{pmatrix} \mathbf{a} \\ \mathbf{b} \\ \mathbf{c} \end{pmatrix}$$

Fractional → Cartesian: $\mathbf{r} = \mathbf{f} \cdot \mathbf{M}$

Cartesian → Fractional: $\mathbf{f} = \mathbf{r} \cdot \mathbf{M}^{-1}$

For a triclinic cell:

$$\mathbf{M} = \begin{pmatrix}
a & 0 & 0 \\
b\cos\gamma & b\sin\gamma & 0 \\
c\cos\beta & c(\cos\alpha - \cos\beta\cos\gamma)/\sin\gamma & c\sqrt{1 - \cos^2\alpha - \cos^2\beta - \cos^2\gamma + 2\cos\alpha\cos\beta\cos\gamma}/\sin\gamma
\end{pmatrix}$$

## Trajectory

```rust
pub struct Trajectory {
    pub frames: Vec<Frame>,
    pub metadata: Option<TrajectoryMetadata>,
}
```

NPT trajectories (variable box per frame) are handled naturally: each `Frame` carries its own `Cell`.

## Element vs. Label

`Atom::element` 是**化学元素**；`Atom::label` 是可选的**位点类型**。两者分开存放，
所以一个原子既能按元素选也能按位点选，而不必二选一：

```rust
Atom { element: "P".into(), label: Some("P_3".into()), .. }
```

分析命令的两组选择参数互斥：`-a/-b/-c` 按 `element`，`-x/-y/-z` 按 `label`。

### 标签约定

位点标签一律形如 `<元素>_<后缀>`，按**首个下划线**拆分。后缀里再有下划线不影响
（`Fe_site_A` 的元素仍是 `Fe`）—— 这条对任何来源的标签都成立，不限于 `ferro net`
写的那几种。

| 标签 | 含义 |
|---|---|
| `P_0` … `P_4` | Qn 形成子，数字 = **同元素**连接数 n（P–O–P），见 `network.md` |
| `Al_4` `Al_5` `Al_6` | 非 Qn 形成子，数字 = **配位数** |
| `O_f` | 自由配体（不连任何形成子） |
| `O_n` | 非桥配体（连 1 个） |
| `O_b` | 桥联配体（连 2 个） |
| `O_t` | 三配位配体（连 ≥3 个，tricluster） |
| `Zn`、`Na`… | 修饰子，裸元素符号，无角色后缀 |

> 0.2.1 之前的旧格式（`P0` / `Of` / `On_P` / `Ob_P_P` / `X` / `Zn_f`）已全部替换，
> **不留读取兼容层**：重跑一次 `ferro net --export-traj` 即为新格式，而加格式猜测
> 反而可能把真元素误判成旧标签。

### 谁会填 `label`

| 来源 | 说明 |
|---|---|
| LAMMPS dump reader | `element` 列写成 `P_3` 时拆成 `element="P"` + `label="P_3"`，读完打印一次映射表 |
| extxyz reader | 独立的 `label:S:1` 列 |
| CIF / CP2K / QE reader | 各自的位点名（`O1`、`Fe1`）——注意这些**不合** `<元素>_<后缀>` 约定 |
| `ferro net --export-traj` | 分类结果写进 `label`，`element` 保持元素 |

### 谁会写出 `label`

**writer 从不自作主张把 `label` 折进元素列。** `ferro convert` 无论标签如何都写干净的
元素符号。只有 `ferro net --export-traj` 折叠，且只对 LAMMPS dump——那个格式只有一列
名字放得下第二个名字；extxyz 走独立的 `label:S:1` 列，两个方向都无损。

折叠带守卫：标签不形如 `<元素>_…` 时不折并计数告警。CIF 的 `O1` 折进去后再读回来，
就是一个不存在的元素 `O1`。
