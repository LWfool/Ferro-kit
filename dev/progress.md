# 当前进度

> 各命令的用法与输出列结构见 `docs/src/`；踩过的坑见 `issues.md`；
> 本文件只记**现状**：什么已完成、代码在哪、验证到什么程度。

## 测试总数：574 个（全部通过，clippy 零警告）

| Crate | 测试数 |
|---|---|
| ferro-core | 95 |
| ferro-io | 110（另有 2 个 `#[ignore]`：真实 40 MB CP2K out、296 MB OUTCAR + 19.7 MB vasprun 与 dpdata 对拍，需 `-- --ignored`） |
| ferro-structure | 72 |
| ferro-analysis | 195 |
| ferro-workflow | 23 |
| ferro-cli（lib 70 + bin 2 + 集成 5） | 77 |

版本号 **0.3.2**（workspace 统一；ferro-python 已同步）。
`v0.2.1 → v0.3.0` 的三批破坏性改动清单见 `overview.md`。

`v0.3.1` 相对 `v0.3.0` **全部是新增**（`ferro dataset` 三步、CP2K out reader、
DeePMD npy 读写、`ml` 模块、`array_order`），没有破坏性变更 —— 唯一动到既有
行为的是 `HARTREE_TO_EV` 由旧值改为 CODATA 2018，而它此前全项目无使用点。

`v0.3.2`（未发版）含破坏性改动（`collect` 的产物布局），但按用户要求走 patch 位；
清单见 `overview.md`。

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
- **`network_type.rs`**：`enum AtomType`（`Former{elem,qn,n_bo,cn,bridges_to}` /
  `Ligand{elem,partners}` / `Modifier{elem,cn}` / `Other`）、`TypeParams`、
  `classify_frame[_detailed]`
  - **`label()` 是全项目唯一把类型渲染成文本的地方**；`class_rank()` / `display_rank()`
    取代三个 `*_label_order` 与 `type_sort_key`；`is_bridging()` 取代前缀匹配
  - `classify_frame_detailed` 另返回 ligand→former 邻接（`FrameTypes`）供 linkage 使用
- `data/elements.rs`：含 `symbol_to_z`、`group_number`、`valence_electrons`、
  `is_transition_metal`、**`split_element_label`**（按第一个下划线拆分，返回
  `LabelSplit::{Plain, Split, Unknown}`；不用贪婪前缀匹配，避免 `Pb`→铅的静默误判）
- `data/compounds.rs`、`data/qn_elements.rs`（默认 `{B,P,Si}`）
- `CubeData`、`charge_grid.rs`、`units.rs`（含 `AMU_ANG3_TO_G_CM3`，由 `AVOGADRO`
  导出，供 `ferro info` 报 g/cm³）、`error.rs`
- `Frame::unique_elements()`（替代 5 处重复实现）
- **`array_order.rs`**（2026-08-25）：`matrix3_row_major` / `matrix3_from_row_major`。
  nalgebra 是列优先且无开关，`as_slice()` 给出的是转置，而 npy 的 box/virial、
  extxyz 的 `Lattice`、GPUMD 的 `lattice` 全要行优先。测试用**非对称**矩阵 ——
  对称张量（stress）被转置后数值不变，这个错误只在 cell 上暴露、在 stress 上静默
- **帧选择**（2026-08-22）：`Trajectory::select(start, end, stride)` /
  `select_indices` / `spread_indices`。区间是 **0 基闭区间**，与 `ferro info`
  打印的帧号对齐；`spread_indices` 是 `--number` 的等间隔且**恒含两端**。
  `tail()` 已改用 `select` 实现，帧选择逻辑只此一处

### ferro-io

- **`writers/table.rs`**：`write_table(&Table, path, TableFormat)` —— 分析产物的唯一出口，
  取代 ferro-analysis 里原来的 7 个 `write_*`。ragged 表在**创建文件之前**拒绝
- **`lammps_dump.rs`**：`element` 列写位点标签时拆成 element + label，读取结束打印一次
  映射表。writer 的整数 `type` 列**在整条轨迹上确定一次**（0.2.1）
- **`extxyz.rs`**：`label:S:1` 列（读写两侧）。extxyz 的列自描述，故标签走自己一列、
  `species` 保持纯元素 —— 与 dump「没地方放第二个名字只能折进 element 列」相反
  - **应力（2026-08-26 修正）**：`stress=` 随 ASE（正 = 拉伸，读写两侧变号），
    `virial=` 随 QUIP/GPUMD/DeePMD（eV，正 = 压缩，除 `|det(box)|`、不变号；
    无 Lattice 报错）。两键同在则交叉校验 `virial ≈ stress × V`，不一致报错而
    不选一个。九个数**检查对称性**（规格要求对称，藉此绕开 ASE「列优先」与
    GPUMD「行优先」的相反描述）；**6 分量 Voigt 拒收**（顺序不普适且无法从数据
    检测）。写侧只写 `stress=`，不写 `virial=`
  - `Properties` 的力列 **`force` 与 `forces` 两种拼法都收** —— GPUMD 的
    `train.xyz` 用单数，此前只认复数，读 NEP 数据集会静默丢掉全部受力
  - 符号的依据是**外部实物**：fixture 由 ASE 3.29.0 生成，期望值由 dpdata 1.0.2
    的换算式独立算出（`~/.miniforge3/envs/deepmd`）。自读自写的回路验证不了符号
  - **写侧可选 `virial=`**（`write_extxyz_with` + `StressKey`，2026-08-26）：eV，
    正 = 压缩，乘体积不变号；无 cell 报错。GPUMD 两个键都认但都给时取 virial，
    故 NEP 导出走它。恒只写一个键
- **`cube.rs`**：`read_cube`（可视化）+ `read_cube_as_chg`（Bader 用：Bohr→Å、索引转置、
  密度缩放 `rho_stored = ρ_cube × V_cell_Bohr`），共用 `parse_header()`
- **`cp2k_out.rs`**（2026-08-25）：CP2K MD 的 stdout 日志（坐标/力/应力全打到
  `__STD_OUT__` 时一个文件自足）。定位靠 `MD| Step number` 锚点 + **区间内扫描**，
  不用固定行偏移（偏移随 ensemble 变，NVT/NPT_I/NPT_F 差 10/18/20 行）；区间上界
  是下一个锚点，缺块的帧宁可丢也不借下一帧的数据。单位从文本自读
  （`[hartree]` / `[bar]`），认不出报错。丢帧三类（SCF 未收敛 / 块截断 / 组成不符）
  计数由 `Cp2kOutStats` 带出。
  **文本锚点按 token 序列匹配**（`mod tag` 一张多候选表），对列宽、对齐、缩进
  与制表符免疫；数值行不写死下标（应力靠「解得出三个浮点」筛，cell 从尾部取）。
  四种排版变形 + 多一行表头 + cell 多一列，均有测试钉住结果逐位相同
- **`writers/deepmd.rs`**（2026-08-25）：DeePMD system 目录（`type.raw` +
  `type_map.raw` + `set.NNN/*.npy`）。磁盘上一律二维 `float64`；
  `virial = stress × V` 不变号；半有半无的属性直接报错。`write_deepmd_npy_sets`
  按 set 切分，**余数摊进各 set** 而不是留尾巴（410 帧按 400 切是 205+205）
- **`readers/deepmd.rs`**（2026-08-25）：上者的逆，返回 `Trajectory`（不造数据集
  专用类型，否则 analysis 要依赖 io）。**dpdata 的 float32 与 ferro 的 float64
  都收**；`virial → stress` 除以 `|det(box)|`（system 里没有 volume.npy）；
  额外键（`atom_ener` 等）**告警而非静默丢**
- **`vasp_outcar.rs` / `vasprun.rs`**（2026-08-27）：VASP AIMD 两条路。口径三方
  实证一致（dpdata 1.0.2、用户的 `private/dp_makedataliu.py`、ASE 的 OUTCAR 解析）：
  能量取 `free energy TOTEN`（不是 `sigma->0`）；`in kB` 六个数按 VASP 自己的
  **`XX YY ZZ XY YZ ZX`**（与 extxyz 规格的 `XX YY ZZ YZ XZ XY` 不同），转 eV/Å³
  **不变号**；体积走 `|det|` 而非三个对角线相乘（后者只对正交胞成立）。
  **晶胞逐帧各读各的**，缺块的帧丢弃而不继承 —— 定胞下这个 bug 逐位不可见。
  收敛判据：OUTCAR 读 VASP 自己的 `EDIFF is reached`，vasprun 只能数
  `<scstep>` 与 `NELM` 比，两者会对同一次运行给出不同的丢帧数，故判据随 stats
  打出来。vasprun 走 `quick-xml` **流式**，不建 DOM
- **`aimd.rs`**（2026-08-27）：`AimdStats`（原 `Cp2kOutStats`）+ `AimdFormat` +
  `sniff`（读头 64 行认横幅）+ `read_aimd_with_stats`
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
| `cube_radius.rs` | 硬球占据图 | `ferro map radius` |
| `cube_jump.rs` | 跳跃距离分布 | 实现完整 + 10 个测试，但**没有 CLI 入口**，`lib.rs` 的再导出清单里也漏了它（只能走 `md::calc_cube_jump`）。手册页与 `ferro doc cube-jump` 都在位。见 `plan.md` |
| `cube_sdf.rs` | 团簇 SDF（Kabsch 对齐） | 用 `ferro_core::build_network_graph`，已去 petgraph |
| `scattering_data.rs` | X 射线 / 中子散射因子表 | 供 `sq.rs` 加权 |

`gr.rs` 的字段：`rho_g = ⟨ρ_f·g_f⟩`（供 S(q) 逐帧变换重构）、`volume_std`（两遍算法）、
`rho = N·⟨1/V⟩`。粒子数守恒校验：逐帧比对总数与分组计数，不符即 `Err`。
排序为 `(elem_z, 字符串)` 二级比较。`GroupBy::{Element, Label}`。

`sq.rs` 的 `to_tables(gr)` 恒写全部规范半边；加权 partial `w_ij(q)·S_ij(q)` 逐点相加
恰等于对应 total。散射因子查表失败时告警。

### ferro-analysis / network

单文件 `mod.rs`，依赖 `ferro_core::classify_frame_detailed`。
结果 `NetworkResult` 出六张长表（`composition` / `qn` / `qn_partner` / `ligand_type` /
`coordination` / `linkage`），每张自带 `#` 头。

**Qn 口径 = 文献的 $Q^n_m$（2026-08-20 起）**：`qn` 只数**同元素**连接（P–O–P），
异核连接是 `qn_partner` 的 `m_<X>` 列，总桥连接数 = `n + Σm`。桥氧个数是第三个量
（`n_bo`，池化均值在 `[inputs]` 的 `mean_n_bo`），三簇氧算 1 个桥氧但 2 个连接，
三者只在无三簇氧时相等。改动前 `qn` 是总桥数却渲染成 `P-Q3`，与文献差 130 倍
（`P-Q3` 40.4% vs 文献 n=3 的 0.27%）。判据与文献原文见 `issues.md`。

**表的含义、列结构与 pandas 用法见 `docs/src/analysis/network.md`（467 行，最完整的一份）。**
设计判据与踩过的坑见 `issues.md`「network 重构（0.2.1）编码陷阱」。这里只记实现事实：

- 参数 `TypeParams`（`cutoffs` + `modifier_cutoffs` + `qn_elements`），`--qn` 整体替换名单
- `n_edge_sharing`：共享 ≥2 个配体的形成子对数（共边多面体），非零才告警。传统 Qn
  假设全共角，共边下一个邻居贡献两个桥氧。参考轨迹 3589 对全共角、零共边
- `Bin { count, fraction, sd }`：`fraction` 是逐帧比例的平均，`sd` 是同一序列的样本
  标准差（ddof=1，缺席帧按 0 计入），用 Welford + Chan/Golub/LeVeque 成对合并（rayon 是归并）
- `linkage` 规范半边存储，两端各带元素/同核连接数(`qn_a`)/配位数，`LinkKey` 含配体元素维
- 旧的 `cn.rs`、`ligand_class.rs`、`qn.rs`、`modifier.rs` 已删除（逻辑迁移至 ferro-core）

### ferro-analysis / ml

`filter.rs`（2026-08-25）：数据集帧筛选，与 md/network/dft 并列，纯计算。
`filter_frames(&Trajectory, &FilterParams) -> FilterResult`，判据 `Criterion`
（力 / 应力，批 2 加 O-O 与 Al6）。

- **区间作用于存活序列**，不是原始帧号；原始索引由 `keep` 带出，保留帧可追溯
- **交叉表**（`cross_tab` / `overlaps`）：逐帧对**全部**帧算判定而非只算存活帧，
  才能报出每个判据「独占抓到」多少 —— 漏斗每步只在上一步存活帧上报数，冗余判据
  在那里看着也很能干
- 阈值 0 = 关闭；给了阈值但缺该标签 → 判之前就报错
- **`geometry.rs`**：`min_pair_distance`（周期最小镜像，超出最小镜像上界报错）、
  `count_with_coordination` / `coordination_histogram`（走 `classify_frame`，读
  `cn` **字段**不解析标签）、`first_shell_cutoff`（g(r) 第一峰后的极小 = 配位壳层
  外沿，用 0.02 Å 粗 bin，细 bin 的局部极小是采样噪声）
- **`diagnostics.rs`**：只读模式的四张表 —— min d(O-O) 分布、每帧 Al6 个数、
  Al 配位分布、**rcut 敏感性扫描**。敏感性表在 49Z49P02A 上是陡坡
  （2.15→0.9%、2.45→13.1%、2.75→41.4%），在 43Z43P15A 上是平线（2.1–2.6 全
  100%）—— 同一张表给出相反提示，这是它的价值
- 四条判据同一语义：「保留含 Al6 的帧」写成「删除不含 Al6 的帧」，交叉表不分裂
- **`merge.rs`**：`composition_key`（逐原子元素序列，分组用，**不看目录名**）、
  `sort_atoms`（规范序 (Z,符号)，coord/force 跟同一置换，box/energy 不动）、
  `group_name`（`112_Al32O64Zn16`，下标是实际计数不约分）、`shuffle_order`
  （seed 默认 666，不给也可复现）
- `to_tables()` 出 funnel / criteria / overlap 三张表，**计数走预格式化文本**
  （`Column::Num` 会把 2000 渲染成 `2.000000e3`）

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
- **`cmd/dataset.rs`**：`ferro dataset collect` —— AIMD out → DeePMD system
  目录。**一目录一 system**（2026-08-26 改）：同目录的 `.out` 是同一次运行被
  重启切开的段，合并回去；命名保留目录层级（剥掉公共祖先，其余原样嵌套，
  文件 stem 不进名字），只有一组时直接写进 `-o` 本身。`-o` 必填，
  `--overwrite` 拦覆盖。文件按首个 step 排序，重复帧不去重但区间打出来。
  同目录成分不符报错（不当坏帧丢），单文件解析失败跳过且最后再报一遍。
  `ferro dataset filter` —— 按力（eV/Å）/ 应力（CLI 收 GPa）阈值筛帧，
  `-i` 收 system 目录或其上层（递归找 `type.raw`），`-o` 按相对路径重建，
  **不给 `-o` 即只读**。七张表（三张统计 + 四张诊断）经 `write_all` 出 csv，
  平铺在 `-o` 根下；诊断恒算（实测挂钟无差别），只读模式全打屏、一字不落盘，
  给了 `-o` 则只打三张统计表。`--shuffle` 在**全部判据与
  抽帧之后**打乱写出顺序（seed 默认 666）—— 抽帧要看时间序，而打乱不可逆。
  `ferro dataset merge` —— 按成分分组合并，两种模式（`shuffle` 全局打乱后按
  `--set-size` 切 / `by-source` **一个 system 一个 set**、不打乱不重切，
  并写 `sets_source.txt`）
  - **`--type deepmd|nep|extxyz` 与 train/valid/test 划分**（2026-08-26，
    filter 与 merge 共用 `OutType` / `Split` / `write_split`）：nep 与 extxyz
    一个 system（merge 是一个成分组）一个 `.xyz`，差别只在应力键；`--set-size`
    不到达它们。`--ratio` 收 `8:1:1`（权重非分数），两段即 train:test；
    不给就不划分。产物用 dpgen 的目录名
    后缀 `.train/.valid/.test`（`SPLIT_SUFFIXES` 早已存在，merge 一直在继承它）。
    成员取自打乱序、各部分内部排回帧序 —— 同 seed 逐字节可复现；比例向上取到
    至少 1 帧。四处在读第一个文件前失败：by-source + 划分、by-source + 非
    deepmd、`--suffix` + 划分、输入已带划分后缀
- **`doc.rs`**（2026-08-26）：`ferro doc` —— `docs/src/` 的 24 页经 `include_str!`
  编译进二进制（208 KB）。topic 跟子命令树同名（`ferro doc dataset filter`），
  裸命令列出全部；原样输出 markdown（零依赖），stdout 是 tty 时经 `$PAGER`
  （默认 `less -R`），否则直接打印，pager 起不来就回落。
  **支持按小节寻址**（`Page.section`）：`convert` / `info` / `bader` 没有手册
  专页，是 `cli-reference.md` 的小节，取该 `##` 到下一个同级标题 —— 99 行而
  不是整本 863 行
- **帮助页与 clap 的防漂测试**（`main.rs` 的 `mod help_sync`）：正向断言每个
  长选项都在其富文本页出现（写短名也算），反向断言页里的每个 `--xxx` 在某个
  命令上真实存在（允许交叉引用，跳过「there is no --x」这类否定陈述）。
  页文本从 `include_str!("help.rs")` 按函数名切出来 —— 25 个 `print_*` 改成
  返回字符串比这个检查本身还大。**上线当场抓到 19 处真漏**：`--metal-units`
  12 页没写、`--tau`/`--ncore` 在 rotcorr/vacf/vanhove 共 5 处没写；
  `gr` 的 `--atom-c`/`--label-z` 是有意不写（SelectArgs 与 angle 共享），
  进 `UNDOCUMENTED` 白名单并写明理由
- **帮助页全部按同一模板**（2026-08-26）：一句话用途 + 完整参数表 + 输出布局 +
  2~3 个例子 + `Full documentation:  ferro doc <topic>`。23 页里 19 页 ≤40 行；
  超标的 4 页（顶层 61 / net 54 / convert 59 含 27 行生成的格式表 /
  cp2k 47 / filter 43）大头都是参数表本身 —— **参数表不为凑行数砍**
- 三级帮助全部手写在 `help.rs`（clap 的派生格式塞不下输出列结构这类段落）。
  **叶子命令 `convert` / `info` / `bader` 也走同一模式**（2026-08-22）：`-i` 是
  `Option`，为空即 `wants_help()` → 富文本页；`-h` 仍归 clap 的参数表。两套并存
  是有意的 —— 短表是参数速查，富文本是格式清单与告警说明。`job` 是唯一自己接管
  `-h` 的命令，未动
- **`convert` 支持选帧**（2026-08-22）：`--start` / `--end` / `--stride` /
  `--number`。产物个数**由目标格式决定不设开关** —— 装得下轨迹的写一个文件，
  只装单结构的（POSCAR / data / QE）一帧一个，序号用**原轨迹帧索引**插在扩展名
  之前（`POSCAR_0000` 仍匹配前缀，可读回）。`-o` 语义未动，仍是完整路径
- **`io_dispatch::supported_formats()` 是格式清单的唯一出处**，读/写/多帧三列。
  三件事清单里各占一列而不是一句散文：CP2K 的 `.inp`/`.restart` 只读不写；
  POSCAR / LAMMPS data / QE 只写第一帧且**不警告**。有测试钉住表与 `match`
  分支一致（表里写 `-` 的格式必须真的拒绝写入），否则两处手写的事实会漂
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
| `quick-xml` | 0.38 | ferro-io | vasprun.xml；净新增 1 个 crate，零传递依赖 |
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
- **`ferro bader` 的三个 `.dat` 写在当前目录**，名字是 `<输入stem>_ACF.dat`。VASP 的
  电荷密度一律叫 `CHGCAR`，故同一目录连跑两个体系后一次静默盖掉前一次。已在帮助页
  与手册告知，`--outdir` 待做（`plan.md` 优先级高）
- **`ferro job` 只用输入的第 0 帧**，多帧输入其余帧丢弃。2026-08-24 起**会告警**
  （`multi_frame_warning`，三行：帧数、忽略数、变通命令），不再是静默的。变通是先
  `ferro convert --number N` 抽成单帧文件再逐个跑 job；让 job 自己选帧见 `plan.md`
- **`ferro info` 的密度只报首尾两帧**，不是全轨迹统计。NPT 下要 mean ± σ 请读
  `ferro traj` 产物文件头的 `# volume = <mean> +/- <std>`。元素表里查不到的符号在
  `effective_mass()` 里回退 1 amu，会把密度拉低 —— 该情形有逐符号告警，但**告警只在
  info 里有**，其他用到质量的地方（msd 的权重、vacf）没有同类提示
- `ferro-python` 仍只暴露 gr/msd，未包 net
- **`ferro dataset collect` 读 CP2K / VASP OUTCAR / vasprun.xml**，QE 待扩；三者
  **都未**注册进 `io_dispatch`（`.out` 太通用、`OUTCAR` 根本没有扩展名），故
  `ferro convert -i OUTCAR` 仍不认识它们
- **VASP 的变胞（NPT）路径只经手工构造的文本验证**：用户现有的两份真实数据
  （2000 帧 OUTCAR、425 帧 vasprun）都是定胞 NVT，且用户明确说 VASP 侧不涉及
  NPT。逐帧读胞的代码有测试钉住「第 2 帧的胞来自第 2 帧」，但没有真实变胞数据
- **不支持 ML_FF 的 OUTCAR**（`free energy ML TOTEN` / `ML FORCE`）：无样例可验，
  且 dpdata 对两者用的行偏移不同（14 vs 4），说明差别不止 token 名
- `dataset` 三步（collect / filter / merge）已齐；几何判据只覆盖最小镜像范围，
  未做多层镜像扫描（小胞体系需要时再补）
- **merge 的规范序取 (Z, 符号)**，dpdata 取字母序；两者都靠 `type_map.raw`
  自描述，但同一份数据经 ferro 与经 dpdata 合并，`type.raw` 的数字会不同
- **`ferro-python` 的格式分派是独立实现**（`ferro-python/src/io.rs`），未跟着 CLI 的
  `io_dispatch.rs` 走。2026-08 加的 `.vasp`/`.pos` 扩展名只有 CLI 认，Python 侧仍只认
  `POSCAR`/`CONTCAR` 前缀。两处 match 分支本就是分开维护，改一处不会波及另一处，
  但也意味着**会漂**——加格式时两边都要看一眼
- `spin.rs`：纯共价分子（如 O₂ 三重态）回退奇偶下限，无法给出 MO 简并导致的自旋
