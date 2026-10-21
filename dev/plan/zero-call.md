# 零调用清单

> 生产代码零调用、尚未决定删留的公开项。从 `dev/plan.md` 按需读取。

## 优先级中

### 零调用清单（2026-09-12 实测快照，判断待定）

这不是待删清单。**多数功能还没有进入实际使用、没有调试过**，因此没法区分
「需求确定但 CLI 层还没接」与「写代码过程中的遗留」——`cube_jump` 就是前者的
现成例子（429 行实现 + 10 个测试 + 手册页都在位，缺的只是 CLI 入口）。

留这份快照是为了**将来能做 diff**：等这批功能跑起来之后重扫一次，「用起来了」
与「始终没人碰」自动分开。现在删掉任何一条，都是拿判断力换行数。

按 2026-09-12 的实测（727 个生产函数，扫描排除 `#[cfg(test)]`）：

| 类别 | 位置 | 行 |
|---|---|---|
| ~~整文件零引用~~ | ~~`ferro-analysis/src/trajectory_analysis.rs`、`geometry.rs`~~ 2026-10-19 删（遗留清理 5.8，`6ed8046`） | — |
| ~~整文件零引用~~ | ~~`ferro-workflow/src/templates.rs`（其中 `orca_sp_template` 是 Ferro 不支持的 ORCA）~~ 2026-10-17 随 crate 删除 | — |
| ~~整文件零引用~~ | ~~`args/corr.rs`（`CorrMode`）· `args/traj.rs`（`TrajMode`）~~ 2026-10-21 删（瘦身 S2.3，用户定） | — |
| 单位转换链 | **2026-10-21 更正**：`convert_length`（cp2k / qe_out reader）、`convert_time`（LAMMPS dump 读写）都已有生产调用，移出清单。`convert_energy` 已活（cp2k_md / cp2k_sp 的能量单位，第 3 批），但 `EnergyUnit::{KcalPerMol, KJPerMol, Wavenumber}` 三个变体从不被构造（`from_name` 只出 `Hartree` / `EV`），连带 `KCAL_MOL_TO_EV` / `KJ_MOL_TO_EV` / `WAVENUMBER_TO_EV` 与 `EV_TO_KJ_MOL` / `EV_TO_WAVENUMBER` 只服务它们（`EV_TO_KCAL_MOL` 另有 LAMMPS real 单位在用）；`LengthUnit::{Nanometer, Picometer}` 同理。**2026-10-21 用户定：保留**（刚按 CODATA 2018 校正，属有意留的扩展点） | ~35 |
| 被绕过的算法 | `ferro-analysis/src/dft/bader_grid.rs:375 max_neargrid` —— 见下一条 | 20 |
| ~~取代方案只落实了一半~~ | ~~`AtomType::display_rank` / `is_bridging`~~ 2026-10-21 删（瘦身 S2.9，用户定） | — |
| ~~基础访问器~~ | ~~`Frame::geometric_center` / `is_periodic` / `center_of_mass` / `center` / `velocities_mut`、`Trajectory::frame_mut` / `iter_frames`、`Atom::distance_to`、`Table::n_cols`、`all_frames_have_cells` / `cell_of`、`FrameVerdict::is_clean`~~ 2026-10-21 删（瘦身 S2.8，用户定）。`Cell::wrap_position` 保留（msd 测试造夹具用、手册 data-model 有记载、merge.rs 有互指注释）；`ChargeGrid::lat_dist_i` 归 Bader 挂起组 | — |
| ~~未使用依赖~~ | ~~`ferro-workflow` 的 `serde`~~ —— 2026-09-26 随 serde 整体移除 | — |
| 2026-10-11 重扫新增（审查 G-2） | `Cell::wrap_position`（只有测试用；`Frame::wrap_all` 自己按 pbc 逐轴实现）；`cube_sdf.rs:141-142 ClusterSnapshot::anchor_idx`（`#[allow(dead_code)]` 压住的只写字段）；~~`pub mod corr`~~（2026-10-21 随 `CorrMode` 删）。**已失效**：`Frame::wrap_all` 自 `48bf11e`（审查 M5）起被 `cube_reference_frame` 调用，已移出基础访问器行 | — |
| 2026-10-20 扫描新增（A-17） | ~~`Frame::center_of_mass` / `center` / `velocities_mut`~~（2026-10-21 删）；`ferro-structure` 的 `find_supercell_dims` / `build_box` / `estimate_box_length` 只有测试调用（无 CLI / Python 入口，同 `cube_jump` 属「需求在、入口没接」）。第 12 批改了其中四个（A-1、A-3、A-8、A-11、A-12），行为由测试钉住。**2026-10-21 用户定：box_builder 保留**（完整功能，缺的只是入口） | — |

**方法上要记住的一条**：`cargo clippy` 零警告不说明没有死代码 —— 库 crate 里
`pub` 项不触发 `dead_code` lint，上面全部躲过了它。查零调用要按名字逐个 grep
全仓，不能靠编译器。

~~`Atom`/`Frame`/`Cell`/`Trajectory` 上的 `derive(Serialize, Deserialize)`~~ ——
**2026-09-26 已删**：全仓没有任何格式后端（serde_json / bincode 都没有），derive
从未被调用。连带摘掉三处 nalgebra 的 `serde-serialize`，`Cargo.lock` 里 serde 系
归零。将来要做检查点或 Python pickle 时，四个 derive + 两行 Cargo.toml 即可加回。
