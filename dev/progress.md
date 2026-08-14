# 当前进度

> 各命令的用法与输出列结构见 `docs/src/`；踩过的坑见 `issues.md`；
> 本文件只记**现状**：什么已完成、代码在哪、验证到什么程度。

## 测试总数：439 个（全部通过，clippy 零警告）

| Crate | 测试数 |
|---|---|
| ferro-core | 84 |
| ferro-io | 65 |
| ferro-structure | 70 |
| ferro-analysis | 163 |
| ferro-workflow | 23 |
| ferro-cli（lib 29 + 集成 5） | 34 |

版本号 **0.2.1**（workspace 统一；ferro-python 已同步并复核编译通过）。
`v0.2.1` 之后另有一批**尚未发版**的破坏性改动，清单见 `overview.md`。

## 锚点 tag

**`v0.1.15`** = 批处理/长表重构前的最后状态（单文件输入、`.dat` 宽表、writer 在
ferro-analysis）。此后所有分析产物的文件名、扩展名、列结构全变。

## scripts/

两组脚本，共享层各一份。都读 ferro 的 csv，但读完之后没有共同代码，故不合并：

| 组 | 共享层 | 脚本 | 用途 |
|---|---|---|---|
| 对拍 | `ferrocmp.py` | `compare_rdf/angle/sq/sq_experiment.py` | 跟 dump2analysis / dump2sq 逐点比 |
| 出图 | `ferroplot.py` | `plot_gr/angle/sq/net.py` | 发表级 pdf（+ png 看效果） |

对拍侧的产物名一律由 `ferrocmp.product_name()` 拼（`batch::out_path` 的镜像），
调用点不写文件名字符串；四个脚本均已按 label 段复跑验证。

出图侧样式 `['science','vibrant']` + LaTeX + 四边框；每个脚本顶部一个 `CFG` 配置块。
`plot_net.py` 的 x 轴是**成分**（`file` 列），100 % 堆积柱；`--partner` 展开成
色相 = Qn、同色系明度 = m_<X> 的嵌套条带。多 csv（`-o` 的 suffix 区分 CMD / MLMD）
画成同一刻度下并排多根柱。

---

## 各 Crate 完成状态

### ferro-core

- `Atom`、`Frame`、`Trajectory`、`Cell` 核心数据结构
- **`table.rs`**：`Table` + `enum Column { Num, Text }` —— 分析产物跨层的中立载体
  - `concat_union(label, parts)`：跨文件堆叠 + **列并集**，缺的列留空（NaN）
  - `to_comment_lines()`：渲染成对齐文本，供 `#` 元信息块内嵌
  - `Column::cell()` 定义全项目输出约定：`{:.6e}`，**NaN 渲染为空字段**
- `cell.rs`：`cartesian_to_fractional` / `wrap_position` / `minimum_image` 均返回 `Result`
  - **`interplanar_spacings()` / `minimum_image_cutoff()`**：dᵢ = 1/‖(Mᵀ)⁻¹.row(i)‖，
    最小镜像的正确几何上界（非正交晶胞 d < L，用边长会高估）；`calc_gr` 与
    `CellList::build` 共用同一份
- **`spin.rs`**：`guess_spin`（magmom 求和 → 氧化态+Hund → 电子数奇偶下限，三级回退）、
  `assign_oxidation_states`、`total_electron_count`、`parity_min_multiplicity`；过渡金属取高自旋
  - 验证：ZnP₂O₆→0（抗磁）、MnS→5（d⁵ 高自旋）、Fe₂O₃→10
- **`cluster.rs`**：`build_network_graph`（former–ligand 邻接 + 配体分类 + 个人 Qn +
  连通分量）、`connected_components`（通用并查集）。供 `ferro-structure::find_clusters`
  与 `ferro-analysis::cube_sdf` 复用
- **`network_type.rs`**：`enum AtomType`（`Former{elem,bridging,cn,bridges_to,qn}` /
  `Ligand{elem,partners}` / `Modifier{elem,cn}` / `Other`）、`TypeParams`、
  `classify_frame[_detailed]`
  - **`label()` 是全项目唯一把类型渲染成文本的地方**；`class_rank()` / `display_rank()`
    取代三个 `*_label_order` 与 `type_sort_key`；`is_bridging()` 取代前缀匹配
  - `classify_frame_detailed` 另返回 ligand→former 邻接（`FrameTypes`）供 linkage 使用
- `data/elements.rs`：含 `symbol_to_z`、`group_number`、`valence_electrons`、
  `is_transition_metal`、**`split_element_label`**（按第一个下划线拆分，返回
  `LabelSplit::{Plain, Split, Unknown}`；不用贪婪前缀匹配，避免 `Pb`→铅的静默误判）
- `data/compounds.rs`、`data/qn_elements.rs`（默认 `{B,P,Si}`）
- `CubeData`、`charge_grid.rs`、`units.rs`、`error.rs`
- `Frame::unique_elements()`（替代 5 处重复实现）

### ferro-io

- **`writers/table.rs`**：`write_table(&Table, path, TableFormat)` —— 分析产物的唯一出口，
  取代 ferro-analysis 里原来的 7 个 `write_*`。ragged 表在**创建文件之前**拒绝
- **`lammps_dump.rs`**：`element` 列写位点标签时拆成 element + label，读取结束打印一次
  映射表。writer 的整数 `type` 列**在整条轨迹上确定一次**（0.2.1）
- **`extxyz.rs`**：`label:S:1` 列（读写两侧）。extxyz 的列自描述，故标签走自己一列、
  `species` 保持纯元素 —— 与 dump「没地方放第二个名字只能折进 element 列」相反
- **`cube.rs`**：`read_cube`（可视化）+ `read_cube_as_chg`（Bader 用：Bohr→Å、索引转置、
  密度缩放 `rho_stored = ρ_cube × V_cell_Bohr`），共用 `parse_header()`
- 其余格式：XYZ、PDB、CIF、VASP、CHGCAR、lammps_data、CP2K、QE

### ferro-structure

- `supercell.rs`、`vacuum.rs`、`merge.rs`
- `box_builder.rs`：随机放置 + soft-core 弛豫，O(N) cell-list 加速
  - `parse_formula` 栈式解析，支持嵌套 `()` / `[]`：`Ca3(PO4)2`、`(NH4)2SO4`、`K4[Fe(CN)6]`
  - `resolve_component`：先查 COMPOUNDS，查不到则当化学式解析 → 可建库外化合物的盒子
  - 拒绝空式、空组、零计数、括号错配、未知元素
- `typing.rs`：`classify_trajectory`、`apply_type_labels`
- `cluster.rs`：`find_clusters`（复用 `ferro_core::connected_components`）

### ferro-analysis / md

**0.1.15 之后分析层不再碰文件系统。** 各结果类型提供
`to_tables() -> Vec<(String, Table)>` + `meta_lines() -> Vec<String>`。

| 文件 | 分析 | 要点 |
|---|---|---|
| `gr.rs` | g(r) + CN(r) | CellList 加速；**n² 个有序对**；逐帧体积归一化 |
| `sq.rs` | S(q) | 只对规范半边做 FT；Faber-Ziman 上三角权重 |
| `msd.rs` | 均方位移 | NPT 安全，atom-major 并行；`fit_diffusion` 出 D 与 R² |
| `angle.rs` | 键角分布 | `AngleParams::ends` 按用户写的 `-a`/`-c` 顺序分派 cutoff |
| `vanhove.rs` `vacf.rs` `rotcorr.rs` | 自关联 / 速度自关联 / 转动相关 | |
| `cube_density.rs` | 3D 密度/速度/力分布 | |
| `cube_sdf.rs` | 团簇 SDF（Kabsch 对齐） | 用 `ferro_core::build_network_graph`，已去 petgraph |

`gr.rs` 的字段：`rho_g = ⟨ρ_f·g_f⟩`（供 S(q) 逐帧变换重构）、`volume_std`（两遍算法）、
`rho = N·⟨1/V⟩`。粒子数守恒校验：逐帧比对总数与分组计数，不符即 `Err`。
排序为 `(elem_z, 字符串)` 二级比较。`GroupBy::{Element, Label}`。

`sq.rs` 的 `to_tables(gr)` 恒写全部规范半边；加权 partial `w_ij(q)·S_ij(q)` 逐点相加
恰等于对应 total。散射因子查表失败时告警。

### ferro-analysis / network

单文件 `mod.rs`，依赖 `ferro_core::classify_frame_detailed`。
结果 `NetworkResult` 出六张长表（`composition` / `qn` / `qn_partner` / `ligand_type` /
`coordination` / `linkage`），每张自带 `#` 头。

**表的含义、列结构与 pandas 用法见 `docs/src/analysis/network.md`（467 行，最完整的一份）。**
设计判据与踩过的坑见 `issues.md`「network 重构（0.2.1）编码陷阱」。这里只记实现事实：

- 参数 `TypeParams`（`cutoffs` + `modifier_cutoffs` + `qn_elements`），`--qn` 整体替换名单
- `Bin { count, fraction, sd }`：`fraction` 是逐帧比例的平均，`sd` 是同一序列的样本
  标准差（ddof=1，缺席帧按 0 计入），用 Welford + Chan/Golub/LeVeque 成对合并（rayon 是归并）
- `linkage` 规范半边存储，两端各带元素/桥接数/配位数，`LinkKey` 含配体元素维
- 旧的 `cn.rs`、`ligand_class.rs`、`qn.rs`、`modifier.rs` 已删除（逻辑迁移至 ferro-core）

### ferro-analysis / dft

- `bader.rs`：`BaderAnalyzer` builder、`BaderResult`、ACF/BCF/AVF 输出（Henkelman 格式，
  外部工具按它解析，故**不走 `Table`**）
- `bader_grid.rs`：on-grid / near-grid / off-grid 三种梯度上升，含边缘精化
- `bader_weight.rs`：Yu-Trinkle weight 方法（WS Voronoi + 流分配）
- `chg_sdf.rs`：电荷密度团簇 SDF（Kabsch 对齐 + pull 插值旋转子格）
- 算法规格与 Fortran 逆向结论见 `bader.md`

### ferro-workflow

- `GaussianJobBuilder`
- `Cp2kJobBuilder`：energy/force/geo-opt/cell-opt/md/freq；PBE/BLYP/PBE0/B3LYP/SCAN/
  r2SCAN/HSE06 等；DFT-D3(BJ)；对角化/OT；k 点、涂抹、cube/Molden 输出；
  CSVR/NoseHoover/Langevin/NVE + NPT
  - `auto_spin` 经 `ferro_core::guess_spin` 推断多重度 + UKS（CLI 显式 `--multiplicity` 时关闭）
  - **`cp2k_basis_db`**：从 6 个 CP2K 文件解析的基组/赝势库（**2829 条**），覆盖
    PBE/SCAN/全电子（pob、-ae）；`basis()`/`potential()` 按元素+族前缀+泛函精确匹配
- `QeJobBuilder`：pw.x（scf/nscf/bands/relax/vc-relax/md/vc-md），复用 `guess_spin`

### ferro-cli

单二进制 `ferro`，八个 `fe-*` 已删除。目录结构见 `CLAUDE.md`「代码导航」，
命令与参数见 `docs/src/cli-reference.md`。实现要点：

- 每个子命令自带参数结构体（`--dt` 在 `msd` 与 `vacf` 下语义不同，不共用字段；
  旧 `fe-traj` 是 29 个字段的大结构体）。`CommonArgs` / `SelectArgs` 经
  `#[command(flatten)]` 接入
- **`batch.rs` 对结果类型泛型**，不认识任何分析类型：`expand_inputs`（自展开 glob，
  零匹配报错）、`map_inputs<T>`（串行遍历，轨迹逐条释放；帧内并行不变）、`stack<T>`、
  `write_all`、`Output { dir, label, suffix }`、`Summary`（存**预格式化文本**）
- 三级帮助全部手写在 `help.rs`（clap 的派生格式塞不下输出列结构这类段落）
- `net` 的 `--P-O=2.3` 由 `main` 在 clap 解析前从 argv 剥离
- **`plot.rs` 面板模型**：`Panel`/`Series` + 通用 `render`，一格一个量、一条曲线一个
  文件，颜色按文件跨格一致，图例只画第一格。**500 dpi**，版式按 96 dpi 编写并统一过
  `px()` 缩放。矢量 PDF 方案已验证可用但**因依赖成本回退**（见 `issues.md`）
- **不留兼容层**：输出格式同期变更，留着 `fe-traj` 会让旧脚本「跑成功」却吐出自己
  解析不了的 csv —— 静默坏数据比命令消失难查

### ferro-python（PyO3 绑定）

独立 workspace（pyo3 0.29 extension-module），用 maturin 构建。
`lib.rs` / `types.rs`（`PyTrajectory`）/ `io.rs`（11 格式按扩展名分派）/ `structure.rs` /
`analysis.rs`（`gr_pair` / `gr_all` / `msd` → `dict[str, list[float]]`）。

本 crate 无 `cargo test`（cdylib 绑定层），经 maturin + Python 冒烟测试验证。

---

## 依赖版本（2026-08-08 全量升级至最新）

| 依赖 | 版本 | 声明位置 | 备注 |
|---|---|---|---|
| `nalgebra` | 0.35 | workspace | 0.34→0.35 零代码改动；g(r)/CN 输出逐字节一致 |
| `ndarray` | 0.17 | workspace | |
| `rayon` | 1.8 | workspace | |
| `serde` | 1.0 | workspace | derive |
| `thiserror` | 2.0 | workspace | |
| `anyhow` | 1.0 | workspace | |
| `rand` | 0.10 | workspace | 0.8→0.10 改名三处，见 `issues.md` |
| `clap` | 4.5 | ferro-cli | derive |
| `plotters` | 0.3 | ferro-cli | `default-features = false` + 必须保留 `ttf`；backend 为 `bitmap` |
| `pyo3` | 0.29 | ferro-python | 0.21→0.29 仅需 `skip_from_py_object`；**运行时未验证** |

`cargo update` 在主 workspace 与 ferro-python 均已无可更新项。

---

## 已知限制

- **`ferro-python` 能编译**（2026-08-12 复核：`cargo clean && cargo check` 干净通过）。
  真问题是它作为独立 workspace 被主 workspace 的 `cargo build/test/clippy` 全部跳过，
  断裂不会被自动发现；改公共 API 后须手动补跑。待办是 **pyo3 0.29 的运行时验证**
  （本机无 maturin），优先级中
- `box_builder`：不支持水合物点记法（`CuSO4·5H2O`）—— 水应作为独立 component 传入；
  无 CLI / Python 入口，只能作为库函数调用
- `cp2k_basis_db`：源数据为 gitignore 的 examples/ 6 文件，DB 已固化为静态表
- QE 赝势仅占位 `<El>.UPF`（UPF 与 CP2K GTH 库不通用，由用户提供 `pseudo_dir`）
- **REPL / 脚本模式未实现**。`main.rs` 是可用的子命令分发器，裸 `ferro` 打印分类总览
- 代码库未纳入 rustfmt 管理，格式为手工维护（含刻意的列对齐）；是否采用见 `plan.md`
- `cube_density.rs:188` 用参考帧体积归一化，NPT 下同类偏差；需先定义「NPT 下 3D 密度图
  指什么」再动（见 `issues.md`）
- **按标签选型只在单帧成立**：`traj gr -x P_3 -y O_b` 对多帧标注轨迹会被 `calc_gr` 的
  逐类型粒子数守恒守卫拒绝（实测 `P_3`：149/152/150/150/150）。多帧请按元素选
- `ferro map chg-sdf` 的 `--cubes` 仍是多 cube 聚合成一张 SDF，与 `map` 其余模式
  「一输入一产物」相反。**优先级高**
- `ferro-python` 仍只暴露 gr/msd，未包 net
- `spin.rs`：纯共价分子（如 O₂ 三重态）回退奇偶下限，无法给出 MO 简并导致的自旋
