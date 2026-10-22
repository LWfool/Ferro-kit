# Ferro — 项目概要

> **文档分工**：`docs/src/`（mdBook）= 怎么用，面向使用者 · `dev/` = 为什么这样定 +
> 现状 + 待办，面向开发 · `CLAUDE.md` = 给 Claude 的硬约束与导航。
> 用法问题（CLI 参数、输出列结构、数据模型字段）一律查 `docs/src/`。

## 定位

Rust 重写的计算化学后处理工具链，面向周期性体系（晶体、表面、玻璃）的 MD 轨迹分析。
目标：替代原有 Python 脚本，兼顾性能与可扩展性，未来计划集成深度学习工作流。

## 架构

```
ferro-cli / ferro-python        ← 唯一允许组合多个 crate 的入口层
    ├── ferro-core              ← 纯数据结构 + 静态参考数据
    │                             + Table（分析产物的中立载体）
    │                             + PeriodicCell（三维周期胞、最小镜像）
    ├── ferro-io       → core   ← 格式读写 + write_table（分析产物的唯一出口）
    ├── ferro-structure→ core   ← 超胞、真空层、合并、初始盒子
    └── ferro-analysis → core   ← 纯计算，不碰文件系统；结果暴露 to_tables()
                                  md/、network/（含 AtomType 分类与 NeighborList）、
                                  dft/（Bader、ChgSDF）、
                                  ml/（数据集筛选与合并）
```

中间层 crate 之间不允许互相依赖。**但共享类型可以下沉** —— 判据：一个类型该放
`ferro-core`，当且仅当两个以上中间层需要叫出它的名字。两个方向各一个范例：
`Trajectory`（io 产出、analysis 消费）与 `Table`（analysis 产出、io 消费）。分析私有的中间产物
（`GrResult`、`SqResult`…）留在 `ferro-analysis`，只有可序列化投影（`to_tables()` 出的 `Table`）跨层。

这条规则的完整论证见 `issues.md`「分析产物为什么不进 ferro-io」—— ASE 与 pymatgen
是 Python，**能**写出循环依赖却依然不写，说明这是设计选择而非语言约束。

## 命令入口（0.2.0 起；`net` 于 0.2.1 降为叶子命令）

**单二进制 `ferro`**，子命令按**产物**分组（不是按实现它的 crate）。命令树、参数与输出列
结构见 `docs/src/cli-reference.md`；`-o` / `-i` 的硬约束见 `CLAUDE.md`，不在这里再抄一份
（命令树抄在多处正是 `issues.md`「一处语义改了、四处文档跟不上」那类漂移）。

`dataset` 的产物是**目录**（DeePMD 的 system 就是目录），故它的 `-o` 是输出根目录。
0.3.3 把这一条推成了全仓规则（`-o` 恒为路径），它不再是例外。

## 编码规范

见 `CLAUDE.md`「编码约定」（注释语言、公式写法）。**不在两处各写一份** —— 这段
此前是 CLAUDE.md 那几行的副本，2026-10-17 规则改成英文时已经漂过一次。

## 版本规则

规则本身见 `CLAUDE.md`「版本规则」与 `memory.md`「版本号更新规则」。这里只记 dev 独有的：

**已知例外**：`v0.2.1`（2026-08-12，net 重构）含破坏性改动（`net qn`/`net type`
删除、标签格式全变、表名与列结构全变），按规则本应是 0.3.0，但按用户当时的
明确要求走了 patch 位。翻 git 历史时注意：**0.2.0 → 0.2.1 之间有 breaking change**。
`v0.3.2`、`v0.3.3`、`v0.3.4`、`v0.3.6` 同属此类，见下面各节。

## 锚点 tag

**`v0.1.15`** = 批处理/长表重构前的最后状态（单文件输入、`.dat` 宽表、writer 在
ferro-analysis）。此后所有分析产物的文件名、扩展名、列结构全变。

**`v0.3.2-alpha1`**（指向 `c7627ff`）= 2026-09-12 简化重构前的最后状态（按 R1–R6 的全库审计，
判据见 `issues.md`「全库简化审计发现的陷阱」）。命名带 `-alpha1` 而不是直接用 `v0.3.2`：
0.3.2 当时**尚未发版**，占掉那个名字会让将来的发版点无名可用（`v0.1.15` 当年能直接用版本号，
是因为它已经发过了）。

---

以下各版本批次**按时间倒序**。用户可见的「改了什么」（旧 → 新命令行对照、迁移清单）在手册
`changelog.md` 的对应版本节，这里不再列表；只留每批的范围、版本号为什么这样定、库 API 的变化
（不进 changelog）与「为什么」。

## 0.4.0 之后（未发版）

用户可见的变化在手册 `changelog.md` 的 `Unreleased`。

- **库 API**（不进 changelog）：2026-10-18 – 10-21 各批（审查第 3–5 步、遗留清理、第 12–14 批、瘦身 S2 / P3）
  曾把 Rust 库的变化写进 changelog，2026-10-21 按「库 API 不进 changelog」的既有分工移回这里（瘦身 D-2，用户定）。
  下面是原文（英文），按原所在的 changelog 小节标出：

- 原在「Analyses」：

  >  Library: `GrResult.params.dr` is the request, `GrResult.dr_used` the
  >   actual width (`params.dr` used to hold the actual width).

- 原在「Analyses」：

  > - Library: `ferro_analysis::{distance, angle, dihedral, radius_of_gyration, bounding_box,
  >   mean_squared_displacement, rmsd_to_first, center_of_mass_trajectory,
  >   radius_of_gyration_trajectory}` are removed. Nothing called them; they took no minimum image
  >   and did not unwrap, so `mean_squared_displacement` sat next to `calc_msd` giving a different
  >   answer. Use `calc_msd`, and `Cell::minimum_image` for distances.

- 原在「Analyses」：

  > - Library: `Table::validate` and `Table::concat_union` return `ferro_core::Result`
  >   (`ChemError::ValidationError`) instead of `Result<_, String>`.

- 原在「Analyses」：

  > - Library: `find_supercell_dims` refuses a NaN or infinite `min_length`, and `build_box` /
  >   `estimate_box_length` a NaN or infinite `density` / `min_dist`. NaN used to pass the
  >   `<= 0` checks: `min_length` was then ignored, the others panicked.

- 原在「Analyses」：

  > - Library: `build_box` takes a `seed` (last argument); the same seed gives the same box. It drew
  >   from an unseeded generator, in an element order that changed from run to run.

- 原在「Analyses」：

  > - Library: `MsdParams`, `VacfParams`, `RotCorrParams`, `BondLifeParams` and `VanHoveParams`
  >   default `dt` to NaN, and `validate` / `calc_*` refuse it as "dt must be given"; the default of
  >   1.0 fs used to compute silently. The CLI and Python already required it.

- 原在「Command line」：

  >   Library: `ferro_io::read_charge_grid` reads a `.cube` or a CHGCAR by name, the choice `bader`
  >   made on its own.

- 原在「Messages: one layout, on stderr」：

  > - Rust library: the per-format readers and `ReadReport::warnings` return
  >   `Vec<ferro_io::Notice>` (`level`, `line`, `message`) instead of `Vec<String>`.

- 原在「Masses: no 1 amu stand-in」：

  > - Library: `effective_mass()` (on `Atom`, `Site`, `AtomRef`) and `Frame::total_mass` return
  >   `Result`.

- 原在「Library: removed and moved items」：

  > ### Library: removed and moved items
  >
  > No command-line behaviour changes; Rust callers change as follows.
  >
  > - `ferro_analysis::ml` is split by task. What `ferro dataset` uses now lives under
  >   `ferro_analysis::ml::dataset` and is reached only there: `ml::filter_frames` and
  >   `ml::filter::filter_frames` are both `ml::dataset::filter_frames`, likewise `Split`, the merge
  >   functions, `DEFAULT_SEED`, `frame_table`, `filter_diagnostics` and `quantiles`. The geometry
  >   helpers (`min_pair_distance`, `count_with_coordination`, `coordination_histogram`) and the
  >   remaining diagnostics functions are no longer public.
  > - `ferro_io` readers return what they noticed alongside the frames, under the plain name: the
  >   `_with_warnings` / `_with_stats` variants are gone, and `read_lammps_dump`, `read_deepmd_npy`
  >   and `read_qe_out` return `(Trajectory, Vec<Notice>)` or `(Trajectory, AimdStats)` (take `.0`
  >   for the frames alone). `read_aimd_with_stats` is `read_aimd`.
  > - No longer reachable from `ferro_io`: `read_as`, `read_xyz`, `read_extxyz`, `read_pdb`,
  >   `read_cif`, `read_poscar`, `read_contcar`, `read_lammps_data`, `read_cp2k_md`, `read_cp2k_sp`,
  >   `read_vasp_outcar`, `read_vasprun`, `read_cp2k_inp`, `read_cp2k_restart`, `read_qe_input`,
  >   `is_extxyz`, `apply_type_map`, `write_xyz`, `index_width`, `indexed_path`. Read any format
  >   through `ferro_io::read` (its warnings and AIMD statistics are in `ReadReport`), or a program
  >   output through `read_aimd` / `read_aimd_as`.
  > - Removed, with no caller left: `ferro_analysis::VERSION`; `MsdFit::{frac_lo, frac_hi}` (the range is `MsdParams.fit_range`);
  >   `AtomType::{display_rank, is_bridging}`; `Frame::{geometric_center, center_of_mass, center,
  >   is_periodic, velocities_mut}`; `Trajectory::{frame_mut, iter_frames}`; `Atom::distance_to`;
  >   `Table::n_cols`; `FrameVerdict::is_clean` (test `flags == 0`).

- 原在「Writers: fewer silent results」：

  > - Library: `ferro_io::write` / `write_numbered` take `&WriteOptions { units, drop_partial }`
  >   in place of the `units` argument; `write_extxyz_with`, `write_lammps_dump` and
  >   `write_lammps_data` take the `drop_partial` flag.

- 原在「Network classification belongs to the analysis library」：

  > Library:
  >
  > - `Frame.charge`, `Frame.multiplicity`, `Frame.bonds`, `TrajectoryMetadata.timestep` and
  >   `Trajectory::time_at` are removed (nothing produced them).  `TrajectoryMetadata.source` is split
  >   into `title` and `program`.
  > - `AtomType`, `TypeParams`, `CutoffTable`, `classify_frame`, the network graph and `NeighborList`
  >   moved from `ferro_core` to `ferro_analysis`; `PeriodicCell` stays in core.  New
  >   `calc_network_typed` returns the per-frame types with the statistics.
  >   `ferro_structure::find_clusters` and `classify_trajectory` are removed.
  > - `ferro_analysis::ml::merge::sort_atoms` is `Trajectory::canonically_sorted`; the (Z, symbol)
  >   order is `ferro_core::data::elements::canonical_cmp`.  New `Frame::virial` and
  >   `Cell::stress_from_virial`.
  > - `ferro_analysis` and `ferro_io` are reached from the crate root only (`ferro_analysis::GrParams`,
  >   not `ferro_analysis::md::gr::GrParams`); `ml` stays a module.  The time-correlation analyses
  >   check the frame spacing themselves, so `check_frame_spacing` is no longer public.
  > - The `--ratio` split, the stitching of restart segments, the per-frame table of
  >   `collect --type inspect` and the filter's diagnostic tables are in `ferro_analysis::ml`.

- 原在「A third of the memory」：

  > Library: `Frame.atoms`, `Frame.forces` and `Frame.velocities` are no longer public fields.  Read
  > atoms with `atom(i)` / `atoms()` / `positions()`, forces with `forces()`; build frames with
  > `add_atom` / `set_atoms`, and set forces with `set_forces`, which refuses a wrong length.

- 原在「Reader warnings reach Python」：

  > Library: `calc_cube_density` and `calc_cube_jump` return `Result` and say why they failed (no
  > periodic cell, no velocities / forces, too few frames, atoms changing between frames, an empty
  > element filter) instead of `None`; `SqResult` carries its warnings.

- 原在「Python `write`: CIF, POSCAR and LAMMPS data get one file per frame」：

  >  Library callers of `write_cif` / `write_poscar` /
  > `write_lammps_data` get an error for several frames instead.

- 2026-10-22（`3baec43`）：`GrResult` 新增 `key_element`（分组键 → 元素，取自首帧原子；标签跨两种元素
  记为空串）。`calc_sq_from_gr` 改经它查散射数据：此前按标签分组时从标签字符串反推元素，`Os` 会读成锇。
  CLI 的 sq 固定按元素分组、Python 不暴露 sq，用户入口不受影响。

- 2026-10-22（`a66f0a1`、`17b7d1a`）：删除 `AngleStats` 与 `AngleResult.stats`（逐三元组 mean / std / count），
  `to_tables` 仍只出 `angle` 一张表（`a66f0a1` 曾加 `stats` 表，`17b7d1a` 撤掉）；`AngleResult` /
  `VanHoveResult::meta_lines` 不再含逐输入的量（angle 的类型清单与 `[statistics]`，vanhove 的 τ、time、
  atoms、origins、elements）。`md/util.rs::unwrapped_series` 为 crate 内共用（`76e4a71`），不是公开项。

- 2026-10-22（`5686678`、`1201faf`）：新增 `FilterParams::validate` 与 `TypeParams::validate`；`filter_frames`
  与 `classify_frame(_detailed)`（经它的 `calc_network`）入口先调。此前库调用方的 NaN / 负阈值静默关判据、
  NaN 截断零成键，现返回 `ValidationError`；0 仍合法。`write_deepmd_npy_bounds` 拒收不首尾相接的边界（`2d02d2c`）。

## `v0.4.0`（2026-10-17 打 tag，发版）

输入生成迁到 ferro-python 收尾（`decisions.md`「`ferro.jobs`」）与 dataset 审查的修复。删了
`ferro job` 子命令与 `ferro-workflow` crate，用户定升次版本位 0.4.0。用户侧变化见 changelog 的 `v0.4.0` 一节。

- **库 API**：删 `ferro_core::spin` 与 `elements::{group_number, valence_electrons, is_transition_metal}`；
  `merge::sort_atoms` 逐帧各取置换（此前取第 0 帧的置换套全部帧）

## `v0.3.6`（2026-10-15 打 tag，发版）

**一次发出两批**：从未单独发版的 `0.3.5`（下一节）与 2026-10-14 的 2026-10-11 全库审查 15 条修复。
两批都含破坏性改动，按用户选择走 patch 位（同 `v0.3.2`、`v0.3.4`）。`v0.3.5-alpha` 是中途锚点，
不是发版点。用户侧迁移清单见 changelog 的 `v0.3.6` 一节。

- **库 API**（不进 changelog）：`Cell` 新增 `origin` 字段（字面量构造须补）；`AimdStats` 新增 `n_no_cell`；
  `rust-version = 1.89`；`calc_sq_from_gr` 拒收 r_min≠0 的 `GrResult`

挂起未修的条目见 `plan/suspended.md`。

## 2026-10-03 – 10-09 的一批（版本号 0.3.5，**随 `v0.3.6` 发版**）

2026-10-09 打 tag `v0.3.5-alpha`（中途锚点，不是发版点）。2026-10-03 四路审查（对照 ASE /
MDAnalysis / dpdata / pymatgen / CP2K 源码）严重与中两级的修复，39 条，逐条见 `decisions.md`
「2026-10-03 四路审查：严重 / 中两级已修条目」。用户可见的变化在 changelog 的 `v0.3.6` 一节。

- **库 API**（不进 changelog）：`ferro_core::classify_frame[_detailed]`、`build_network_graph`、
  `ferro_structure::classify_trajectory`、`ml::geometry` 两个配位函数、`cutoff_scan`、
  `pooled_coordination` 改返回 `Result`；`check_frame_spacing` 由 CLI 下沉为
  `ferro_analysis::md::check_frame_spacing`；`diagnostics::distribution_table` 去掉 `name` 参数

## `v0.3.4`（2026-10-03 打 tag，发版）

**一次发出三批**：从未单独发版的 `0.3.2`、`0.3.3`（含 2026-09-21 ~ 09-27 挂在 0.3.3 下的
各批），以及 2026-10-03 的物理审查修复（下一节）。三批都含破坏性改动，按用户要求
走 patch 位。`v0.3.2-alpha1` / `v0.3.3-alpha1` 两个 tag 是中途锚点，不是发版点。
用户侧迁移清单见 changelog 的 `v0.3.4` 一节。

## 2026-10-03 的一批（随 `v0.3.4` 发版）

2026-10-02 物理审查的中等问题修复（审查 M-A – M-F、低 1 – 6）。最大的一条是 `--metal-units` →
`--units real|metal`、dump 含速度/力列时必填：两种体系速度差 10³、力差 23 倍，漏写开关是静默的
数量级错误，同 M1 对 atom style 的裁定 —— 不猜。中子散射长度表与 NIST 逐项核对（用户选「先全表
核对」），B 保持 ¹¹B（用户裁定：硼玻璃中子实验用富集样品）。

- **库 API**：`read_lammps_dump` / `write_lammps_dump` 收 `Option<LammpsUnits>`，`LammpsUnits` 去掉
  `Default`；`find_supercell_dims` 与 `calc_cube_radius` 改返回 `Result`
- 手册补了审查第四节的 6 条口径说明；第 1 条（Lorch 窗）待加功能，第 6、7 条等 Kabsch 修复
  （`plan/suspended.md` S-A）时一起写
- 上线即抓到 `tests/70Z30P00A_NVT_5.lammpstrj` 的步距是 12000/12000/13000/12000 —— 这条 fixture
  不能再喂给时间相关分析

## 2026-09-26 的一批（版本号**仍是 0.3.3**，未发版）

依赖精简：`traj gr/sq/msd/angle` 的 `--plot` 移除，连带删掉 plotters 与从未调用过的 serde derive，
运行时依赖 138 → 97（取舍见 `decisions.md`「依赖精简」）。09-27 在同一个 0.3.3 里追加：msd /
vacf / rotcorr 改全原点 FFT（数值全变）、msd 头部只留共享参数、新增 `traj bondlife` 与 rotcorr 的
`--vector bond` / `--legendre`。

## 2026-09-22 的第二批（版本号**仍是 0.3.3**，未发版）

起因是实测 `ferro dataset collect -i 1Al/*.out` 报「different compositions」而
两个化学式打出来一模一样。根因是 collect 用**逐原子原始顺序**比较成分、却用
**计数**渲染报错，而 CP2K 按输入文件列原子的顺序写坐标表，同一体系的两次单点
计算原子顺序常常不同。这批让 collect 逐文件规范排序、新增 `collect --format`、
25 页帮助精简为五段、手册新增 `changelog.md`。

### 为什么无条件重排、不给开关

查了本机 deepmd 环境（`~/.miniforge3/envs/deepmd`）的两份实现，两边都在做同一
件事：

- **DeePMD-kit 3.1.3 加载时就重排**。`deepmd/utils/data.py` 的 `_make_idx_map`
  默认 `sort_atoms=True`，用 `np.lexsort((idx, atom_type))` 按类型稳定排序；
  文档字符串写明非 mixed type 的 descriptor **要求**如此。磁盘顺序对训练没有
  影响
- **dpdata 1.0.2 的 `append` 自动排**。`system.py:491` 在 `atom_types` 不一致时
  调 `sort_atom_types()` 重排两边，注释写着 "allow to append a system with
  different atom_types order" —— 正是 collect 要做的事

顺带对上了一个口径：ferro 的 `type_map` 按 `(Z, 符号)` 排，重排后 `type.raw`
单调不减，于是 DeePMD 内部那个 lexsort 成为**恒等置换**。

给开关等于承认「不排」也是一种合法产物，而它唯一的效果是让同目录的两个文件因为
原子顺序不同而无法合并。`merge` 一直是无条件排的，collect 不排才是那个例外。

**mixed type 是例外**：DeePMD 在那条路上明确不排序（`real_atom_types` 逐帧各异）。
`plan/io.md` 的 mixed type 待办已记下这一条。

### 为什么嗅探没有按最初的要求删掉

用户最初提的是「删除嗅探，完全依赖 `--format`」，理由是版本多、易错、难维护。
查证后留下了：`sniff` 只有 40 行且**不含任何版本知识** —— 担心的那部分全在四个
reader 的块匹配里，删掉 `sniff` 一行都减不掉。留着它换来两样东西：与 `--format`
不一致时的一行 NOTE（否则手滑写错与有意覆盖长得一样），以及「同目录 OUTCAR +
vasprun.xml」那条守卫的比较对象（`--format` 一给，`stats.format` 就恒等于它）。

### 晶胞精度：明确选了低精度的那个源

`CELL| Vector` 的分量只打 3 位小数，而 `|a|` 与角度打 6 位。从长度+角度重建更
精确但会**强制 a 轴沿 x**，且 ≤7.1 连 `|a|` 都只有 3 位、根本无处可重建。实测
体积相对误差 1.16e-4，20 GPa 下 virial 绝对偏差 3.8e-3 eV，低于 virial 训练
RMSE 一到两个数量级——**按矢量行读，不复杂化**。

## 2026-09-22 的一批（版本号**仍是 0.3.3**，未发版）

起因是 `private/cp2kdata`（robinzyb/cp2kdata，MIT），核心需求是**单点能数据的
读取**——把一批 FP 计算的输出收成训练集。范围只取能量 / 力 / 坐标 / 盒子 /
应力五项，与 MD 端一致；cube、pdos、Mulliken / Hirshfeld 布居、振动频率、偶极、
DFT+U 占据数**全部不取**（要么 Ferro 已有更稳健的实现，要么 `Frame` 里没有存放
处）。GEO_OPT / CELL_OPT 也不取——cp2kdata 里那两条路的坐标根本不来自 `.out`，
而是去 glob `*pos*.xyz`，属于已划在范围外的 CP2K 多文件布局。

这批新增 `readers/cp2k_sp.rs`（`AimdFormat::Cp2kSp`），`cp2k_out.rs` 改名 `cp2k_md.rs`，
`sniff` 按 `GLOBAL| Run type` 二次分派，`AimdStats` 加 `version` / `version_note`，MD 路补两处老版本的块。

### 为什么单点单开一个文件而不是在 `cp2k_md.rs` 里分支

用户明确要求**不共用单位与数据抽取**（"强行兼容可能会给未来造成隐患"）。事实上
两者也确实无可共用之处：MD 的帧锚点是 `MD| Step number`、坐标与力是两个 xyz 块；
单点的锚点是 `ENERGY| Total FORCE_EVAL`、坐标与力是两张表。两个 `version_note`
函数看着像重复，但它们判定的是不同的东西（MD 需要 xyz 块与 cell 行，单点需要坐标
表与力表），将来会各自分岔——按 `CLAUDE.md` 的 R6，没有同一个变化原因就不合。

### 版本策略：提示而不是拒绝

用户最初提的是白名单 + 拒绝，查证后收敛成**照常解析 + 打提示**。查了 cp2kdata
全部 28 份 fixture，格式断代是三段：

| | ≤ 7.1 | 8.1 – 2024 | 2025 + |
|---|---|---|---|
| 能量 | `energy (a.u.):` | `energy [a.u.]:` | `energy [hartree]` |
| 力 | `ATOMIC FORCES in [a.u.]` | 同左 | `FORCES| Atomic forces` |
| 应力 | ` STRESS TENSOR [GPa]` | `STRESS| Analytical …` | 同左 |

即「2025/2026 + 可能还有 2023/2024」这个白名单里装着**两种**布局，不是一种。

**最终口径**：明确实现 **2023–2024** 与 **2025–2026** 两代。2025–2026 静默读
（2026 沿用 2025 布局，据用户确认，未实测）；**2024 及之前每文件打一行 `NOTE:`
点名版本，然后照常提取**。7.1 及更早"过老，暂不考虑"，解析代码与 fixture 保留
但不列为支持 —— `cp2k_sp_v61.out` 是 `ATOMIC FORCES` 块唯一的真实文件样例，而
2023–2024 用的正是那个块。

两代的差别只有两处（坐标表、kind 块、`CELL|` 三者两代同形，已在真实
2023.1/2023.2/2024.1 输出上逐块核对）：

| | 2023–2024 | 2025–2026 |
|---|---|---|
| 能量 | `energy [a.u.]:` | `energy [hartree]` |
| 力 | `ATOMIC FORCES in [a.u.]` | `FORCES| Atomic forces [hartree/bohr]` |

2023–2024 **没有单点 fixture**（手上那两代的样例全是 AIMD），所以它落在"打提示"
那一侧；块组合由 `the_2023_2024_generation_parses_to_the_same_numbers` 钉住。

硬拒绝被否掉的理由：它会把**没见过但格式没变**的将来版本一并挡掉，而
`cp2k_md.rs` 的模块文档里早就写着相反的判断（token 匹配而非版本表，且
`STRESS_UNIT` 是输入关键字、不是版本的函数）。

### kind 名进 `label`，`type_map` 仍按 element

CP2K 允许同元素多 kind。查 fixture 时发现两种：`Fe1`/`Fe2`（自旋初猜），以及
v7.1 那份里 kind 名叫 `Al` 而坐标表 Element 列是 `Fe 26` 的取代体系。

`element` 恒取坐标表里的真元素，kind 名只在与元素不同时进 `label`。
`type_map.raw` 仍按 element 建，于是 `Fe1`/`Fe2` 合成一个训练类型——DeePMD
模型本就是按元素参数化的。但这扔掉了用户特意做出的区分，所以 `collect` 每个
system 打印一次映射表，照 LAMMPS dump reader 的先例。

**dpdata 的 CP2K 插件默认相反**（`true_symbols=False`，拿 kind 名当 `atom_names`），
所以同一份输出两边建出的 type_map 可能不同。

## 2026-09-21 的一批（版本号**仍是 0.3.3**，未发版）

**含破坏性改动，但按用户明确要求不升版本号** —— 这一批没有自己的版本号，
它与 `v0.3.2`、`v0.3.3` 挤在同一个 `0.3.3` 里。翻 git 历史时注意：**`0.3.3` 这个版本号下
有多批改动，且都是破坏性的**。

起因是把 `private/cp2k_grep_strInfo.py` 的能力并进来。该脚本的九项能力里 Ferro
已有六项且实现更稳健，真正合进来的是 **LAMMPS 导出 + 逐帧标量序列**（`collect --type inspect`）；
CP2K 的多文件布局与 DeePMD mixed type 明确划在范围外（仍在 `plan/io.md`）。同批修了 LAMMPS dump
三斜盒子行（改 `*_bound`），`Frame` 加 `temperature` / `step`。

推翻的旧判据两条，都写进了 `issues.md`：

- **「`collect` 的 `-o` 必填」**（0.3.2 定）。原判据反对的是默认值 `.`（会把 npy
  撒进正在工作的目录），而跟着输入走的默认值撒不到别处去。结论改为：**默认值必须
  跟着输入走，不能是 `.`**
- **「`collect` 不加后缀」**（0.3.3 定，理由是产物是没筛过的原始数据、不该标成
  训练集）。改为加 `.db`（database），但它**不是**划分后缀，不进 `SPLIT_SUFFIXES`

## `v0.3.3`（2026-09-20，未发版）

**含破坏性改动但按用户明确要求走 patch 位** —— 与 `v0.2.1`、`v0.3.2` 同一类例外。
主体是 `-o` 恒为路径（`--outdir` 删除）、批次标记改 `-s/--suffix`、`--mkdir` 与创建确认、
`bader` 报告改写到输入旁边、`dataset` 产物默认带 `.train`。

**内部**：`ferro-io` 的 reader/writer 路径参数由 `&str` 统一为 `&Path`；
`BaderResult` 的三个 `write_*` 改为 `*_text()` 返回字符串，落盘移到 CLI ——
`ferro-analysis` 不再碰文件系统。新增 `tests/CHGCAR_2atoms` 与
`ferro-cli/tests/bader_reports.rs`。


## `v0.3.2`（2026-08-26，未发版）

**含破坏性改动但按用户明确要求走 patch 位** —— 与 `v0.2.1` 同一类例外。翻 git
历史时注意：**0.3.1 → 0.3.2 之间 `dataset collect` 的产物布局不兼容**（分组单位由文件改为
目录，collect 拼**同一次运行**的碎片、merge 合**不同运行**），extxyz 的 `stress=` 变号
（旧行为是错的，不是另一种约定）。

新增 `ferro doc`、VASP AIMD 读取、`dataset filter` / `merge` 的 `--type` 与 `--ratio`。两条判据：

- `collect` 改为**按文件头几行的横幅嗅探**格式：VASP 写的叫 `OUTCAR` 无扩展名，`.out` 又太通用，
  按名字判必然开一串特例。`Cp2kOutStats` 下沉改名为 `AimdStats`（第二个使用者到齐，符合下沉判据）；
  新增依赖 `quick-xml`，流式不建 DOM —— AIMD 的 vasprun 常有几百 MB
- GPUMD/NEP 的 `train.xyz` 就是 extxyz，故导出**没有新命令** —— 缺的只是「读 system 目录」，
  而 filter/merge 本来就在读。重复帧**不去重**（重启重叠段位置速度相同，能量力也相同），
  但 step 区间打出来让重叠可见

## `v0.3.1`：机器学习数据集（2026-08-26 发版）

新增 `ferro dataset` 三步（`collect` / `filter` / `merge`），配套
`ferro-io` 的 CP2K out reader 与 DeePMD npy 读写、`ferro-analysis/src/ml/`、
`ferro-core/src/array_order.rs`。**全部是新增，无破坏性改动** —— 唯一动到既有
行为的是 `HARTREE_TO_EV` 由旧值 27.211396132 改为 CODATA 2018 的
27.211386245988，而它此前全项目无使用点。

三个跨命令的约定在这批里定下，其他地方也适用：

| 约定 | 出处 |
|---|---|
| **落盘矩阵一律行优先**，取九个数走 `matrix3_row_major`，禁用 nalgebra 的 `as_slice()`（它是列优先，且对称张量下这个错误完全静默） | `array_order.rs` |
| **`Frame::stress` 的符号 = 正为压缩**，与 CP2K/VASP/QE 输出一致、与 ASE/GPUMD 的 stress 相反；`virial = stress × V` 不变号 | `frame.rs` 字段文档 |
| **单位从输出文本自读**，认不出报错不默认 —— CP2K 的 `STRESS_UNIT` 是输入关键字，单位不是版本的函数 | `readers/cp2k_md.rs` |
| **AIMD 的能量取自由能**（CP2K 的 `FORCE_EVAL`、VASP 的 `free energy TOTEN`），不取 `sigma->0` —— 力是自由能对坐标的导数，配 `sigma->0` 等于给模型两半不同的泛函 | `readers/vasp_outcar.rs` |

`v0.3.1` 是纯新增，故升 patch 位。

## `v0.3.0`（2026-08-20 发版）

`v0.2.1` → `v0.3.0` 累积了三批破坏性改动，此前按用户要求一直压着不发版，等配套的
Python 绘图脚本跟进后一并升：net 重做（六张表改名、Al 退出 Qn 表、标签数字换义，2026-08-12）、
traj 产物命名加 label 段（2026-08-13）、net 口径改为文献 $Q^n_m$（**数值全变**，2026-08-20）。
**旧产物与旧命令行都不兼容**，对照表见 changelog 的 `v0.3.0` 一节。

$Q^n_m$ 口径的依据是 arXiv 2510.13545 等文献原文：n 只数同核桥，m 数异核桥，总桥由两者相加；
铝磷酸盐文献专门指出「n 是总桥氧数」是领域内的已知误解。与实验对比时的两条口径选择见手册
`analysis/network.md`「Comparing with experiment」。

三批合起来按规则升次版本位，故 `0.2.1 → 0.3.0`（跳过 0.2.x）。
