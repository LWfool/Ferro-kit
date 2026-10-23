# 当前进度

> 各命令的用法与输出列结构见 `docs/src/`；踩过的坑见 `issues.md`；
> 本文件只记**现状**：什么已完成、代码在哪、验证到什么程度。

## 测试总数：933 个（全部通过，clippy 零警告；2026-10-23 `ferro build` 落地后实数）；另 ferro-python 195 个 pytest

| Crate | 测试数 |
|---|---|
| ferro-core | 77 |
| ferro-io（lib 295 + 集成 4） | 299（另有 2 个 `#[ignore]`：296 MB OUTCAR + 19.7 MB vasprun 与 dpdata 对拍、`examples/qe/` 的 396 原子 QE 7.3.1 scf 与 GIPAW / pp.x 反例，需 `-- --ignored`；2026-10-17 复跑通过） |
| ferro-structure | 78 |
| ferro-analysis | 312 |
| ferro-render | 31 |
| ferro-cli（lib 108 + bin 5 + 集成 23） | 136 |

版本号 **0.4.0**（workspace 统一；ferro-python 已同步；**2026-10-17 发版，tag `v0.4.0`**）。
`rust-version = 1.89` 自 0.3.6 起声明（依赖树最高要求，nalgebra 0.35）。
各版本批次、破坏性改动、锚点 tag 见 `overview.md`。

## tests/ 的 fixture

| 文件 | 来源 | 覆盖 |
|---|---|---|
| `43Z43P15A_NPT_5.lammpstrj` / `70Z30P00A_NVT_5.lammpstrj` | 生产轨迹等间隔取 5 帧 | 正交胞，一 NPT 一 NVT |
| `CHGCAR_2atoms` | 8×8×8 合成网格 | Bader 端到端 |
| `vasp_OUTCAR_2frames` / `vasp_vasprun_2frames.xml` | 真实 VASP 运行裁出 | 定胞 NVT；温度 0.00 / 111.55 K |
| **`cp2k_md_3frames.out`**（2026-09-21） | `examples/cp2k/md/total.out` 裁 3 帧，165 KB | CP2K 端到端。**此前一份 CP2K 样例都没有**，那条链只有代码里构造的字符串。选这三帧是因为第 2、3 帧的瞬时温度与累计平均**不相等**，取错列会当场失败 |
| **`5Al_0003_1500K_f394.out`**（2026-09-22，用户提供） | 真实 CP2K 2025.2 单点，457 原子 | 新布局：`FORCES|`、`STRESS| … [bar]`、`CELL_TOP|` 干扰项 |
| **`cp2k_sp_v61.out`**（2026-09-22） | cp2kdata 的 `test_energy_force/v6.1/normal/output`，MIT，107 KB 整份照抄 | 老布局三个块一次到齐：`energy (a.u.):`、`ATOMIC FORCES in [a.u.]`、无前缀的 ` STRESS TENSOR [GPa]`；外加 gamma=120 的三斜胞与 `Fe1`/`Fe2` 两个 kind。断言数值全取自 cp2kdata 自己的 `answer/`，是跟独立实现对表 |
| **`triclinic_2frames.lammpstrj`**（2026-09-21） | 合成 | 三斜盒子（lx/ly/lz=10/12/14，xy/xz/yz=2/-3/1），cell 与 ASE 读出的逐位相同。三斜没有真实来源，而正交胞上三斜的 bug 完全不可见。由 `lammps_dump.rs` 的 `test_triclinic_bounds_are_not_xlo_xhi` 读（瘦身 A-8） |
| **`qe75_{sp,sp_spin,sp_noconv,md,vcmd}.out`**（2026-10-17，超算实算） | QE 7.5 pw.x，6 原子 $AlPO_4$，PBE-D3，`verbosity = 'high'`，共 410 KB；`ferro-python/scripts/qe_fixtures_hpc.sh` 生成，输入与作业脚本在 `examples/qe/pwx/` | 总力后跟七组力分量；`nspin = 2`；`electron_maxstep = 3` 不收敛；5 步 SVR md 与 5 步 vc-md。能量、力与 dpdata 一致到 1e-10，其余与 ASE 一致到 1e-7 |
| **`vasp651_sp/frame_{0011,0161,0194}/`**（2026-10-13，用户提供） | VASP 6.5.1 单点（r2SCAN-D3，302 原子）200 帧中挑 3 帧，`OUTCAR` + `vasprun.xml` 原样，1.3 MB | 单点布局。全部 200 帧当时与 dpdata 1.0.2 逐帧对拍（两条路），box / coord / energy / force 逐位相同，virial 差 1.4e-6 eV（换算常数 160.2176634 vs 160.2176621）。原数据随后由用户删除，测试钉住 3 帧的 dpdata 参考值。0011 的 xml 应力非对称最大，读侧取对称部分（见 `plan.md` 归档「vasprun 的非对称应力」） |

## scripts/

只剩对拍与空间群工具；**出图已于 2026-10-20（第 15 批）迁入 ferro-python**（`ferro.plot` +
`ferro-py plot`，见下文 ferro-python 一节与手册 `python/plot.md`）。

| 组 | 共享层 | 脚本 | 用途 |
|---|---|---|---|
| 对拍 | `ferrocmp.py` | `compare_rdf/angle/sq/sq_experiment.py`、`trajcheck.py` | 跟 dump2analysis / dump2sq 逐点比 |

另有两个与上面无关的独立脚本（不读 ferro 的 csv）：`gen_spacegroups.py` 从 spglib 生成
CIF reader 的空间群数据表 `ferro-io/src/readers/spacegroups.dat`；`crosscheck_spacegroups.py`
把 ferro 的 CIF 展开与 spglib、ASE 逐设置对拍。两者都用 deepmd 环境的 python 跑。

对拍侧的产物名一律由 `ferrocmp.product_name()` 拼（`batch::product_name` 的镜像），
调用点不写文件名字符串；四个脚本均已按 label 段复跑验证（2026-08-14）。

**2026-09-20 起未复跑**：`-o` 改为输出目录之后，四个脚本传给 ferro 的批次标记
从 `-o` 改成了 `-s`（它们以 outdir 为 cwd，用不到 `-o`）。改动是机械的，但没有
dump2analysis / dump2sq 在手，无法再跑一遍对拍 —— 下次跑之前先确认这一处。

---

## 各 Crate 完成状态

### ferro-core

- `Frame`、`Trajectory`、`Cell` 核心数据结构。`Frame` 拓扑 / 坐标分离：`Site` 表（元素、标签、
  质量、磁矩）经 `Arc` 在相邻帧间共享，坐标、电荷、力、速度是帧私有数组；`Atom` 只是建帧用的
  自有记录，`AtomRef` 是帧给出的单原子视图（B15）。重排原子走 `Frame::permuted`，按 (Z, 符号)
  规范序整条轨迹重排走 `Trajectory::canonically_sorted`（逐帧各取置换，相邻帧元素序列相同时复用；
  原 `ml::merge::sort_atoms`，2026-10-18 下沉）
- **`periodic.rs`**：`PeriodicCell`（缓存逆矩阵，`to_fractional` / `minimum_image` / `spacings` /
  `mic_bound`），各分析入口建一次、逐对复用；`Cell::check(what)` 是有限性与体积 > 1e-12 Å³ 的
  唯一判据（2026-10-20，A-1 + A-13）
- **`table.rs`**：`Table` + `enum Column { Num, Text }` —— 分析产物跨层的中立载体
  - `concat_union(label, parts)`：跨文件堆叠 + **列并集**，缺的列留空（NaN）
  - `to_comment_lines()`：渲染成对齐文本，供 `#` 元信息块内嵌
  - `Column::cell()` 定义全项目输出约定：`{:.6e}`，**NaN 渲染为空字段**
- `cell.rs`：`cartesian_to_fractional` / `wrap_position` / `minimum_image` 均返回 `Result`
  - **`Cell.origin`**（2026-10-14，审查 A-12）：盒子原点（LAMMPS 的 lo），构造器一律填 0。
    只是元数据，同 ASE `celldisp`：分数坐标、wrap、最小镜像与全部分析都不用它；只有 LAMMPS
    dump / data 读写它，`make_supercell` / `add_vacuum` / `merge_frames` 认它，其余格式不读不写
  - **`interplanar_spacings()` / `minimum_image_cutoff()`**：dᵢ = 1/‖(Mᵀ)⁻¹.row(i)‖，
    最小镜像的正确几何上界（非正交晶胞 d < L，用边长会高估）；`calc_gr` 与
    `box_builder` 的 `CellList` 共用同一份
- **网络类型不在 core**：`AtomType` / `TypeParams` / `classify_frame[_detailed]`、former–ligand 图
  （`build_network_graph`）与并查集 `connected_components` 于 2026-10-18 移到
  `ferro-analysis/src/network/`（只有分析层叫它们的名字，按「两个以上中间层需要」的判据不该在
  core），见下文 ferro-analysis / network
- `data/elements.rs`：元素表（`by_symbol` / `by_number`）、`symbol_to_z`、
  **`canonical_cmp` / `canonical_order`**（(Z, 符号) 规范序的唯一出处）。`group_number` /
  `valence_electrons` / `is_transition_metal` 已随 `ferro job` 删除；`split_element_label`
  （按第一个下划线拆分，不做贪婪前缀匹配，避免 `Pb`→铅的静默误判）是 LAMMPS dump reader 的
  私有函数（`ferro-io/src/readers/lammps_dump.rs`）
- 化合物表在 `ferro-structure/src/compounds.rs`（`box_builder` 唯一使用方）；Qn 形成子默认名单
  `{B,P,Si}` 在 `ferro-analysis/src/network/qn_elements.rs`
- `CubeData`、`charge_grid.rs`、`units.rs`（含 `AMU_ANG3_TO_G_CM3`，由 `AVOGADRO`
  导出，供 `ferro info` 报 g/cm³；`ASE_FS` = `ase.units.fs`，取 ASE 默认 CODATA 2014 下的值，
  供 extxyz 的 `momenta` 换算）、`error.rs`
- `Frame::unique_elements()`（替代 5 处重复实现）
- **`array_order.rs`**（2026-08-25）：`matrix3_row_major` / `matrix3_from_row_major`。
  nalgebra 是列优先且无开关，`as_slice()` 给出的是转置，而 npy 的 box/virial、
  extxyz 的 `Lattice`、GPUMD 的 `lattice` 全要行优先。测试用**非对称**矩阵 ——
  对称张量（stress）被转置后数值不变，这个错误只在 cell 上暴露、在 stress 上静默
- **`Frame.temperature` / `Frame.step`**（2026-09-21）：逐帧标量，跟着帧走而不是
  与帧数组并行的旁路 —— 后者在 `select` / `stride` 之后立刻错位且不报错。
  `temperature` 三条 AIMD 路都填（vasprun 的是反算的），`step` 只有 CP2K 与
  OUTCAR 有。`units.rs` 连带新增 `BOLTZMANN_EV_K`
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
  - **速度**（2026-10-14，审查 A-2）：读侧 `velocities` / `vel`（GPUMD）/ `velo`（QUIP）按 Å/fs 直收，
    ASE 的 `momenta` 除以质量（`masses` 列优先）再乘 `ASE_FS`；写侧写 `masses` + `momenta`
    两列（ASE 与 GPUMD 都不认 `velocities`），ASE 3.29 `get_velocities()` 读回差 ~1e-11
  - 符号、行列序、Voigt 拒收的依据见 `issues.md`「extxyz 的 stress/virial 编码陷阱」
  - **写侧可选 `virial=`**（`write_extxyz_with` + `StressKey`，2026-08-26）：eV，
    正 = 压缩，乘体积不变号；无 cell 报错。GPUMD 两个键都认但都给时取 virial，
    故 NEP 导出走它。恒只写一个键
- **`cube.rs`**：`read_cube`（可视化）+ `read_cube_as_chg`（Bader 用：Bohr→Å、索引转置、
  密度缩放 `rho_stored = ρ_cube × V_cell_Bohr`），共用 `parse_header()`
- **`cp2k_md.rs`**（2026-08-25，2026-09-22 由 `cp2k_out.rs` 改名）：CP2K MD 的 stdout 日志（坐标/力/应力全打到
  `__STD_OUT__` 时一个文件自足）。定位靠 `MD| Step number` 锚点 + **区间内扫描**，
  不用固定行偏移（偏移随 ensemble 变，NVT/NPT_I/NPT_F 差 10/18/20 行）；区间上界
  是下一个锚点，缺块的帧宁可丢也不借下一帧的数据。单位从文本自读
  （`[hartree]` / `[bar]`），认不出报错。丢帧三类（SCF 未收敛 / 块截断 / 组成不符）
  计数由 `Cp2kOutStats` 带出。
  文本锚点按 token 序列匹配（`mod tag` 一张多候选表），数值行不写死下标，依据见
  `issues.md`「CP2K out 解析 / DeePMD 数据集编码陷阱」。2026-09-22 补两处老版本的块：`energy (a.u.):` 的圆括号写法（此前整份 bail 在
  `unknown energy unit []`）、无 `STRESS|` 前缀的 ` STRESS TENSOR [GPa]`（此前帧
  照读、stress 恒 None、数据集少一个 virial.npy，全程无错误）
  **virial 与 cp2kdata 对拍**（2026-10-17，CP2K 2025.2 NVT，`examples/cp2k/md/` 的 `total.out` 2000 帧
  与 `01_glass` 1110 帧；01_glass 下四份 `.out` 逐字节相同，实为一份）：参考值取 cp2kdata 的
  `parse_stress_tensor_list` × `CELL| Volume` ÷ dpdata 的 GPa 常数（即它的 dpdata 插件口径，不变号），
  按能量逐帧对齐。virial 相对差 ≤ 2.1e-10，符号一致。体积须取 `CELL| Volume`：`CELL| Vector` 只印
  3 位小数，三斜胞下按它算体积差 ~1e-4。ferro 不收的力计算只有第 0 步、重启重复的那一步，与
  「MD run terminated by external request」前只算了力、未成步的最后一次（SCF 也未收敛）
- **`cp2k_sp.rs`**（2026-09-22）：CP2K 单点（`ENERGY` / `ENERGY_FORCE`）。与
  `cp2k_md.rs` **不共用任何单位或抽取函数**，两者无一处布局相同。帧锚点是
  `ENERGY| Total FORCE_EVAL`：坐标与胞向前找、力与应力向后找，两侧由相邻锚点
  夹住，于是 `cat` 起来的 N 份单点输出天然读成 N 帧。两代力块（`FORCES|` 与
  `ATOMIC FORCES`）、两代应力块、三种能量括号写法都认。kind 名进 `label`，
  `element` 恒取坐标表里的真元素。`CELL|` 精确匹配，不含 `CELL_TOP|`/`CELL_REF|`。
  缺应力不丢帧（没有 virial 而已），缺力必丢（`force.npy` 不是可选项）。
  明确实现 **2023–2024** 与 **2025–2026** 两代，后者静默、前者每文件打一行
  `NOTE:` 后照常提取
- **`qe_out.rs`**（2026-10-17）：QE pw.x 输出，`scf` / `md` / `vc-md` 一个 reader。步锚点
  `Self-consistent Calculation`；能量、力、应力在本步内，**结构印在上一步末尾**（第 0 步取
  输出头：`celldm(1)` × `crystal axes`、alat 单位坐标，精度约 1e-6 Å），最后一步之后的卡片
  丢弃。力只取 `Forces acting on atoms` 下第一组（verbosity high 还印七组同格式分量）；应力取
  Ry/bohr³ 列（kbar 列只两位小数），不变号。力必需、应力可缺、一帧不剩报错（同 cp2k_sp）；
  relax / vc-relax 与多横幅拼接报错。换算用 CODATA 2018（= QE 7.5 `constants.f90`）
- **`writers/deepmd.rs`**（2026-08-25）：DeePMD system 目录（`type.raw` +
  `type_map.raw` + `set.NNN/*.npy`）。磁盘上一律二维 `float64`；
  `virial = stress × V` 不变号；半有半无的属性直接报错。`write_deepmd_npy_sets`
  按 set 切分，**余数摊进各 set** 而不是留尾巴（410 帧按 400 切是 205+205）
- **`readers/deepmd.rs`**（2026-08-25）：上者的逆，返回 `Trajectory`（不造数据集
  专用类型，否则 analysis 要依赖 io）。**dpdata 的 float32 与 ferro 的 float64
  都收**；`virial → stress` 除以 `|det(box)|`（system 里没有 volume.npy）；
  额外键（`atom_ener` 等）**告警而非静默丢**
- **`vasp_outcar.rs` / `vasprun.rs`**（2026-08-27）：VASP AIMD 两条路（各条口径的
  依据见 `issues.md`「VASP AIMD 读取编码陷阱」）：能量取 `free energy TOTEN`（不是 `sigma->0`）；`in kB` 六个数按 VASP 自己的
  **`XX YY ZZ XY YZ ZX`**（与 extxyz 规格的 `XX YY ZZ YZ XZ XY` 不同），转 eV/Å³
  **不变号**；体积走 `|det|` 而非三个对角线相乘（后者只对正交胞成立）。
  **晶胞逐帧各读各的**，缺块的帧丢弃而不继承 —— 定胞下这个 bug 逐位不可见。
  收敛判据：OUTCAR 读 VASP 自己的 `EDIFF is reached`，vasprun 只能数
  `<scstep>` 与 `NELM` 比，两者会对同一次运行给出不同的丢帧数，故判据随 stats
  打出来。vasprun 走 `quick-xml` **流式**，不建 DOM。vasprun 的应力取 $(\sigma+\sigma^T)/2$（2026-10-13）：VASP 6.5.1 写的 3×3 不对称，
  原样保留则导出 nep 的 `virial=` 被本仓 extxyz reader 拒收；与 dpdata 在该路差 ~1e-3 eV
- **`aimd.rs`**（2026-08-27）：`AimdStats`（原 `Cp2kOutStats`）+ `AimdFormat` +
  `sniff`（读头 96 行认横幅）+ `read_aimd_as(path, fmt)` + 它的嗅探包装
  `read_aimd`（2026-10-21 前叫 `read_aimd_with_stats`）。2026-09-22 拆成两个：调用方在自己更清楚格式时
  （`collect --format`）绕开嗅探。同日 `sniff` 的 CP2K 分支改为**读不到
  `GLOBAL| Run type` 就报错**并点名 `--format` —— 此前默认按 MD 走，
  一份被裁过头的单点日志会报成「MD 找不到 Step number 锚点」，与真正的原因无关。
  2026-10-17 加 `QeOut`：认 `Program PWSCF`；GIPAW、pp.x 等其他 QE 程序点名拒收
- 其余格式：XYZ、PDB、CIF、VASP、CHGCAR、lammps_data、CP2K、QE
- **2026-10-03 审查 A 组的修复**（判据见 `decisions.md` 与 `issues.md`）：CIF 无 symop 时按
  空间群符号展开（`spacegroup.rs` + spglib 生成的 530 条 Hall 设置表）、跳过无 `_atom_site_`
  的 block、无胞帧拒写；PDB 每个 MODEL 前写本帧 CRYST1，非标准取向先转坐标；LAMMPS 无胞帧盒子
  外扩 1 Å、近正交不写三斜、dump 边界按 pbc；extxyz 注释行按 ASE 规则切词；DeePMD 见
  `real_atom_types.npy` 即拒收；POSCAR / CHGCAR 元素行取 POTCAR 标签前缀；CP2K 段名整行比较、
  `&CELL_REF` 不再覆盖 `&CELL`；QE namelist 按逗号与空白切分；LAMMPS data 空文件报错
- **2026-10-11 审查 A 组的修复**（2026-10-14，随 `v0.3.6`）：LAMMPS dump / data 读写 `Cell.origin`
  （A-12，三斜 dump 的 `*_bound` 也加原点、非规范胞绕原点旋转）；LAMMPS data 段头取 `read_data`
  文档全名单、Masses / Atoms 按段名定位不依赖顺序、Masses 行须恰为 `type mass`、重复 type 报错
  （A-10）；CP2K MD 缺 cell 行计入 `AimdStats.n_no_cell` 并在报错里点名 `&CELL`，第三个同构 xyz
  块（VELOCITIES 打到 stdout）整份报错（A-17）；CIF 判重按 ASE `symprec`：逐分量 1e-3、只在同一
  非对称位点的像之间、跨位点重合报错（A-7）。A-13（QE 一行 namelist）2026-10-16 随 `ferro.jobs` 第 1 步修复

### ferro-structure

- `supercell.rs`、`vacuum.rs`、`merge.rs`
  - `merge_frames`（2026-10-09 重写放置）：**不拉伸**，两块整体平移、胞中心在面内居中；
    gap 沿界面法向量（同 `add_vacuum`）；界面夹角 > 1° 报错。共格应变匹配（ASE `stack`、
    pymatgen `from_slabs`）有意不做。2026-10-14（审查 E-10）：原子先相对 `Cell.origin` 在周期
    方向折回各自胞内，块中心取 origin + 半体对角线；非周期方向越界报错
  - `add_vacuum`：厚度须有限且 > 0，奇异胞（c=0）报错（2026-10-14，审查 E-9）；
    它与 `make_supercell` 都保留 `Cell.origin`
- `box_builder.rs`：随机放置 + soft-core 弛豫，O(N) cell-list 加速
  - `parse_formula` 栈式解析，支持嵌套 `()` / `[]`：`Ca3(PO4)2`、`(NH4)2SO4`、`K4[Fe(CN)6]`
  - `resolve_component`：先查 COMPOUNDS，查不到则当化学式解析 → 可建库外化合物的盒子
  - 拒绝空式、空组、零计数、括号错配、未知元素
- CLI 入口 `ferro build supercell | vacuum | merge | box`（`ferro-cli/src/cmd/build.rs`，2026-10-23）；
  Python 四者都有（`supercell` / `add_vacuum_layer` / `merge` / `build_box`）
- `typing.rs`：`apply_type_labels`（`ferro net --export-traj` 用）。`classify_trajectory` 与
  `cluster.rs::find_clusters` 已于 2026-10-18 删除（零调用；changelog 记为库 API 移除）

### ferro-analysis / md

**0.1.15 之后分析层不再碰文件系统。** 各结果类型提供
`to_tables() -> Vec<(String, Table)>` + `meta_lines() -> Vec<String>`。

| 文件 | 分析 | 要点 |
|---|---|---|
| `gr.rs` | g(r) + CN(r) | `NeighborList`（`neighbor.rs`）加速；**n² 个有序对**；逐帧体积归一化 |
| `sq.rs` | S(q) | 只对规范半边做 FT；Faber-Ziman 上三角权重 |
| `msd.rs` | 均方位移 | NPT 安全，atom-major 并行；`fit_diffusion` 出 D 与 R² |
| `angle.rs` | 键角分布 | `AngleParams::ends` 按用户写的 `-a`/`-c` 顺序分派 cutoff |
| `vanhove.rs` `vacf.rs` `rotcorr.rs` | 自关联 / 速度自关联 / 转动相关 | |
| `cube_density.rs` | 3D 密度/速度/力分布 | |
| `cube_radius.rs` | 硬球占据图 | `ferro map radius` |
| `cube_jump.rs` | 跳跃距离分布 | 实现完整 + 测试，2026-10-08 补 `validate()`，已从 crate 根导出（`calc_cube_jump` 等），但**没有 CLI 入口**。手册页 2026-10-21 撤出手册，原样存于 `plan/cube-jump-manual.md`。见 `plan/suspended.md`（挂起的 `map jump` 入口组） |
| `cube_sdf.rs` | 团簇 SDF（Kabsch 对齐） | 用 crate 内的 `network::build_network_graph`，已去 petgraph |
| `util.rs` | `build_avg_frame`（cube 文件头的时间平均帧）；`check_frame_spacing`（帧等间隔、不重复，2026-10-09 由 CLI 下沉，CLI 五个时间相关分析与 Python `msd` 共用） | 2026-09-12 由三份逐字相同的实现合并而来，同时修掉三份都带的 PBC bug：解缠在**参考帧 cell** 的分数坐标下做，故没跨边界的原子逐位不变 |
| `scattering_data.rs` | X 射线 / 中子散射因子表 | 供 `sq.rs` 加权 |

`gr.rs` 的字段：`rho_g = ⟨ρ_f·g_f⟩`（供 S(q) 逐帧变换重构）、`volume_std`（两遍算法）、
`rho = N·⟨1/V⟩`。粒子数守恒校验：逐帧比对总数与分组计数，不符即 `Err`。
排序走 `canonical_cmp`（`ferro-core/src/data/elements.rs`）。`GroupBy::{Element, Label}`。

2026-10-07 – 10-08 的口径修正：`vanhove --tau` 默认 N/2（同 `--max-lag`）；`GrResult.r_max_used`
存截断后的值，头部写请求值；S(q) 的 `q-min` 必须 > 0，`[inputs]` 带 `q_trunc = 2π/r_max`，
低于它告警；`check::non_negative` / `ordered` 要求有限值；vanhove / map sdf / cube_jump 入口查胞可逆。

`sq.rs` 的 `to_tables(gr)` 恒写全部规范半边；加权 partial `w_ij(q)·S_ij(q)` 逐点相加
恰等于对应 total。散射因子查表失败时告警。**变换恒从 r=0 积分**（2026-10-14，审查 B-1）：
`calc_sq_from_gr` 拒收 `r_min ≠ 0` 的 `GrResult`，`traj sq` 没有 `--r-min`（与 gr 不再共用 `GrKnobs`）。

四种 cube（density/velocity/force、radius、jump、sdf）的 `CubeData.origin` 取半个体素
（2026-10-14，审查 B-2）：值属于体素中心，cube 规范第 i 点在 origin + i·spacing。

### ferro-analysis / network

`mod.rs`（六张表的统计）+ `classify.rs`（`AtomType` / `TypeParams` / `classify_frame[_detailed]`）+
`graph.rs`（`build_network_graph`、`connected_components`）+ `qn_elements.rs`；2026-10-18 由
`ferro-core` 的 `network_type.rs` / `cluster.rs` 下移而来。
- `enum AtomType`（`Former{elem,qn,n_bo,cn,bridges_to}` / `Ligand{elem,partners}` /
  `Modifier{elem,cn}` / `Other`）；`classify_frame` 与 `build_network_graph` 入口查奇异胞并返回
  `Result`（2026-10-09）—— c=0 的 slab 此前在最小镜像的 `expect` 处 panic
  - **`label()` 是全项目唯一把类型渲染成文本的地方**；`class_rank()` 取代三个
    `*_label_order` 与 `type_sort_key`（零调用的 `display_rank()`、`is_bridging()` 2026-10-21 删）
  - `classify_frame_detailed` 另返回 ligand→former 邻接（`FrameTypes`）供 linkage 使用
结果 `NetworkResult` 出六张长表（`composition` / `qn` / `qn_partner` / `ligand_type` /
`coordination` / `linkage`），每张自带 `#` 头。

**Qn 口径 = 文献的 $Q^n_m$（2026-08-20 起）**：`qn` 只数**同元素**连接（P–O–P），
异核连接是 `qn_partner` 的 `m_<X>` 列，总桥连接数 = `n + Σm`。桥氧个数是第三个量
（`n_bo`，池化均值在 `[inputs]` 的 `mean_n_bo`），三簇氧算 1 个桥氧但 2 个连接，
三者只在无三簇氧时相等。改动前 `qn` 是总桥数却渲染成 `P-Q3`，与文献差 130 倍
（`P-Q3` 40.4% vs 文献 n=3 的 0.27%）。判据与文献原文见 `issues.md`。

**表的含义、列结构与 pandas 用法见 `docs/src/analysis/network.md`（最完整的一份）。**
设计判据与踩过的坑见 `issues.md`「network 重构（0.2.1）编码陷阱」。这里只记实现事实：

- 参数 `TypeParams`（`cutoffs` + `modifier_cutoffs` + `qn_elements`），`--qn` 整体替换名单
- **截断可写 `auto`**（2026-10-07，`cmd/net.rs::resolve_auto`）：逐输入调
  `md::first_shell_cutoffs`（一次 g(r) 查全部 auto 对），头部写 `P-O=auto`，逐输入的
  截断（三位小数，手填回去逐字节相同）与 `g_min` 进 `[inputs]`。输入里没有这对元素时
  截断置 0 照常分析（与定值截断时一致），记 `-`；元素都在但 g 越不过 1 则该文件失败。
  `build_params` 阶段 auto 对占位 NaN，只为校验元素角色与打印标签方案，进
  `calc_network` 前必被替换
- `n_edge_sharing`：共享 ≥2 个配体的形成子对数（共边多面体），非零才告警。传统 Qn
  假设全共角，共边下一个邻居贡献两个桥氧。参考轨迹 3589 对全共角、零共边
- `Bin { count, fraction, sd }`：`fraction` 是逐帧比例的平均，`sd` 是同一序列的样本
  标准差（ddof=1，缺席帧按 0 计入），用 Welford + Chan/Golub/LeVeque 成对合并（rayon 是归并）
- `linkage` 规范半边存储，两端各带元素/同核连接数(`qn_a`)/配位数，`LinkKey` 含配体元素维

### ferro-analysis / ml

2026-10-21 起 `ml` 按任务分命名空间（用户定），现有文件都在 `ml/dataset/`，公开项一律
`ferro_analysis::ml::dataset::<项>` 一条路径，子模块 crate 内可见。

`filter.rs`（2026-08-25）：数据集帧筛选，与 md/network/dft 并列，纯计算。
`filter_frames(&Trajectory, &FilterParams) -> FilterResult`，判据 `Criterion`
（力 / 应力，批 2 加 O-O 与 Al6）。

- **区间作用于存活序列**，不是原始帧号；原始索引由 `keep` 带出，保留帧可追溯
- **交叉表**（`cross_tab` / `overlaps`）：逐帧对**全部**帧算判定而非只算存活帧，
  才能报出每个判据「独占抓到」多少 —— 漏斗每步只在上一步存活帧上报数，冗余判据
  在那里看着也很能干
- 阈值 0 = 关闭；给了阈值但缺该标签 → 判之前就报错
- 判据开着时，力 / 应力峰值非有限的帧判坏（2026-10-08，此前 `fold(0.0, f64::max)` 吞掉 NaN）
- **`geometry.rs`**：`min_pair_distance`（周期最小镜像，超出最小镜像上界报错）、
  `count_with_coordination` / `coordination_histogram`（走 `classify_frame`，读
  `cn` **字段**不解析标签）。`first_shell_cutoff` 已于 2026-10-07 移到
  `md/shell.rs`（改名 `first_shell_cutoffs`，net 的 auto 与 `--al6` 共用）：0.02 Å 粗
  bin，极小在 0.10 Å 滑动平均上找，零平台取最长段中点
- **`diagnostics.rs`**：只读模式的四张表 —— min d(O-O) 分布、每帧 Al6 个数、
  Al 配位分布、**rcut 敏感性扫描**。敏感性表在 49Z49P02A 上是陡坡
  （2.15→0.9%、2.45→13.1%、2.75→41.4%），在 43Z43P15A 上是平线（2.1–2.6 全
  100%）—— 同一张表给出相反提示，这是它的价值。逐 system 的统计（O–O 分位、实际所用
  Al–O 截断）进 `[inputs]` 的列而不进表头（2026-10-08，`quantiles`）—— 表头堆叠后只剩第一个
- 四条判据同一语义：「保留含 Al6 的帧」写成「删除不含 Al6 的帧」，交叉表不分裂
- **`merge.rs`**：`composition_key`（逐原子元素序列，分组用，**不看目录名**）；原子重排已下沉为
  `Trajectory::canonically_sorted`（见 ferro-core；坐标、力、速度跟同一置换 —— 后两者
  2026-09-22 补上，此前静默错位；现在逐帧各取置换，独立写出的文件不必先各自排完）、
  `group_name`（`112_Al32O64Zn16`，下标是实际计数不约分）、`shuffle_order`
  （seed 默认 666，不给也可复现）
- `to_tables()` 出 funnel / criteria / overlap 三张表，**计数走预格式化文本**
  （`Column::Num` 会把 2000 渲染成 `2.000000e3`）

### ferro-analysis / dft

- `bader.rs`：`BaderAnalyzer` builder、`BaderResult`、`acf_text`/`bcf_text`/`avf_text`
  （Henkelman 格式，外部工具按它解析，故**不走 `Table`**）。2026-09-20 起只**渲染文本**，
  落盘由 CLI 做 —— 分析层碰文件系统的唯一例外就此消失
- **`tests/CHGCAR_2atoms`**（8×8×8 合成网格，9 KB）：CHGCAR reader → 分析 → 报告
  这条链的端到端 fixture，此前只有内存里构造网格的单元测试
- `bader_grid.rs`：on-grid / near-grid / off-grid 三种梯度上升，含边缘精化
- `bader_weight.rs`：Yu-Trinkle weight 方法（WS Voronoi + 流分配）
- `chg_sdf.rs`：电荷密度团簇 SDF（Kabsch 对齐 + pull 插值旋转子格）
- 算法规格与 Fortran 逆向结论见 `bader.md`

### ferro-cli

单二进制 `ferro`，八个 `fe-*` 已删除。目录结构见 `CLAUDE.md`「代码导航」，
命令与参数见 `docs/src/cli-reference.md`。实现要点：

- 每个子命令自带参数结构体（`--dt` 在 `msd` 与 `vacf` 下语义不同，不共用字段；
  旧 `fe-traj` 是 29 个字段的大结构体）。`CommonArgs` / `SelectArgs` 经
  `#[command(flatten)]` 接入
- **`-o` 恒为路径**（2026-09-20 起，`--outdir` 已删）：写多个产物的命令（13 个分析
  命令、`chg-sdf`、`bader`、`dataset`）指目录，只写一个文件的（`convert`）
  指文件；批次标记走 `-s/--suffix`。目录不存在时经 **`outpath.rs`** 向 stderr 问一句
  `[y/N]`，非交互环境报错并指明 `--mkdir`，询问都在读第一个输入之前
- **`batch.rs` 对结果类型泛型**，不认识任何分析类型：`expand_inputs`（自展开 glob，
  零匹配报错）、`unique_labels`（同名输入改 `_2`，撞名按小写比 —— macOS / Windows 文件系统
  不分大小写，2026-10-14 审查 E-12）、`map_inputs<T>`（串行遍历，轨迹逐条释放；帧内并行不变）、`stack<T>`、
  `write_all`、`Output { dir, label, suffix }`、`Summary`（存**预格式化文本**）
- **`cmd/dataset/`**（2026-10-19 由 `dataset.rs` 拆为 `mod` / `collect` / `filter` / `merge` / `inspect`）：
  三步的行为、参数与产物布局见手册 `docs/src/dataset/`，这里只记实现位置与用户裁定
  - collect：一目录一 system（2026-08-26）；PJM 作业日志在分组前剔除（`drop_pjm_logs` / `is_pjm_log`，
    2026-10-05，规则写死、无关闭开关，用户裁定）；失败按输入文件计，原因随条目收集、跑完统一列一次，
    目录级条目只带前面没列过的文件（2026-10-23）；每个文件读完当场按 (Z, 符号) 重排
    （`Trajectory::canonically_sorted`，2026-09-22）。`--format` 压过嗅探，但嗅探仍跑：
    不一致时打 `NOTE:`，并作「同目录 OUTCAR + vasprun.xml」守卫的比较对象（`--format` 一给，
    `stats.format` 就恒等于它，守卫会失效）
  - filter：`-f` / `--s-max` 是 `Option`，默认 20 / 10 经 `FilterCmd::f_max()/s_max()` 取 ——
    用户没给时缺该标签的 system 只关这一条，显式给了才报错（2026-10-14，审查 D-4）；应力阈值没有
    短名（`-s` 是批次标记，审查 E2-2）；诊断恒算（实测挂钟无差别）；`--shuffle` 在全部判据与
    抽帧之后（抽帧要看时间序，打乱不可逆）
  - merge：`shuffle` / `by-source` 两种模式，by-source 自 2026-10-19 起在 system 内部打乱
  - `--type deepmd|nep|extxyz` 与 `--ratio` 划分（2026-08-26）：filter 与 merge 共用
    `OutType` / `Split` / `write_split`；划分后缀用 dpgen 的目录名（`SPLIT_SUFFIXES`），
    `filter` / `merge` 产物默认带 `.train`（2026-09-20 起）
- **`cmd/build.rs`**（2026-10-23）：`build supercell | vacuum | merge | box`，只读、命名、写，结构操作全在
  `ferro-structure`。四者共用 `OutArgs`（`-o` 目录 / `-s` / `--mkdir` / `--format`，格式默认取输入的）与
  `write_product`；supercell / vacuum 走 `per_input`（`batch::map_inputs`，逐输入失败跳过 + 退出码 1），
  merge / box 只出一个产物、失败即报错。手册 `docs/src/build.md`，判据见 `decisions.md`「`ferro build`」
- **`cmd/dataset/inspect.rs`**（2026-09-21）：`dataset collect --type inspect` 的三个诊断产物，
  布局与 `_info.csv` 的列见手册 `dataset/collect.md`
- **`doc.rs`**（2026-08-26）：`ferro doc` —— `docs/src/` 除 `python/` 外的全部页经 `include_str!`
  编译进二进制（`python/` 只由 `ferro-py doc` 收录）。topic 跟子命令树同名（`ferro doc dataset filter`），
  裸命令列出全部；pager 与 `NO_COLOR` 的行为见手册 `cli-reference.md` 的 doc 一节。
  **支持按小节寻址**（`Page.section`）：`convert` / `info` / `bader` 没有手册
  专页，是 `cli-reference.md` 的小节，取该 `##` 到下一个同级标题，不是整本
- **`ferro-render` crate**（2026-09-24；2026-10-06 从 `ferro-cli/src/doc/render.rs` 拆出，
  `markdown.rs` + `latex.rs`，CLI 改动不再重编它）：`ferro doc` 的终端渲染器，**只在 stdout 是
  tty 时**介入，重定向输出仍是源文件原样（24 个整页主题逐字节核对过）。全自写，
  净新增 0 crate（`libc` 只在 `cfg(unix)` 下声明，且早已经 rand → getrandom 在树里）。
  `split_blocks` 按行首切块；`inline` 是一次扫描的
  tokenizer；表格带完整边框、单元格按词折行；行内 `$...$` 按四条语法规则转
  Unicode；块级 `$$` 用同一转换合成一行、缩进 4 格、超宽在空格处折行，
  `\begin{}` 环境（矩阵、cases）保留原文。宽度取 ioctl → `COLUMNS` → 80，正文上限 100 列。
  Windows 下纯 ASCII 样式；`rendering_loses_no_content` 对每页两种样式校验不丢字
- **帮助页与 clap 的防漂测试**（`main.rs` 的 `mod help_sync`）：正向断言每个
  长选项都在其富文本页出现（写短名也算），反向断言页里的每个 `--xxx` 在某个
  命令上真实存在（允许交叉引用，跳过「there is no --x」这类否定陈述）。
  页文本从 `include_str!("help.rs")` 按函数名切出来 —— 25 个 `print_*` 改成
  返回字符串比这个检查本身还大。有意不写的参数进 `UNDOCUMENTED` 白名单并写明理由
  （如 `gr` 的 `--atom-c`/`--label-z`：SelectArgs 与 angle 共享）。另一条（2026-10-09）：带 `default_value(_t)` 的参数 doc 注释里不许再手写
  `[default: …]` —— clap 自己会渲染，曾在 dataset 的六个参数上显示两遍
- **帮助页全部按同一模板**（2026-09-22 重定为五段）：一句话用途 +
  `Parameters:`（完整）+ `Output:`（≤4 行）+ `Examples:` +
  `Full documentation:  ferro doc <topic>`，目标 ≤30 行。判据与口径一律进手册。
  仍超 30 行的页（`dataset filter`、顶层总览）大头是参数表与命令列表本身 ——
  **不为凑行数砍**。`convert` 的格式矩阵并进 `Parameters:`（它是 `-i`/`-o` 的值域）
- 三级帮助全部手写在 `help.rs`（clap 的派生格式塞不下输出列结构这类段落）。
  **叶子命令 `convert` / `info` / `bader` 也走同一模式**（2026-08-22）：`-i` 是
  `Option`，为空即 `wants_help()` → 富文本页；`-h` 仍归 clap 的参数表。两套并存
  是有意的 —— 短表是参数速查，富文本是格式清单与告警说明。（曾经唯一自己接管 `-h` 的
  `job` 已于 2026-10-17 删除，计算输入改由 `ferro-py job` 生成）
- **`convert` 支持选帧**（2026-08-22）；产物个数由目标格式决定、不设开关（用法见手册 §convert）
- **格式清单的唯一出处是 `ferro-io/src/format.rs` 的 `Format` 表**（`can_write` / `multi_frame`）；
  `io_dispatch::supported_formats()` 由 `Format::ALL` 生成读 / 写 / 写出帧数三列，不再手写。
  输入（QE、CP2K `.inp`/`.restart`）与程序输出只读；CIF / POSCAR / LAMMPS data 一帧一个文件
- `net` 的 `--P-O=2.3` 由 `main` 在 clap 解析前从 argv 剥离，**只在 `argv[1] == "net"` 时**
  （2026-10-08；此前对所有子命令剥离，`traj gr --P-O=2.3` 被静默吞掉）；截断须有限且 ≥ 0
  （0 = 这一对不成键），由 `TypeParams::validate` 检查（2026-10-22，C-2）
- **`-s` 只收字母、数字、`_`、`+`、`-`**（`batch::check_suffix`，2026-10-08），traj / map / net /
  chg-sdf / bader 入口查；`dataset merge --suffix` 是目录后缀，另一套语义
- **写侧格式判定是 `ferro_io::Format::for_output`**（2026-10-08 起在读输入前判定；原
  `io_dispatch::out_format` 已并入 `Format` 表），`convert` 在建目录、读输入之前调
- **`convert` 写出前一律按 (Z, 符号) 重排原子**（2026-10-12，用户要求，无开关）：所有格式
  同一体系的原子顺序一致，像 POSCAR 那样同元素成块。走 `Trajectory::canonically_sorted`，
  **逐帧**各取置换。writer 未改 —— vasp / lammps_data 按首次出现分组，排完序后自然是 Z 序
- `--version` 带 commit、构建日期、profile、target、rustc 与 features（2026-10-05）
- **无绘图**：`--plot` 与 `plot.rs` 已于 2026-09-26 移除（理由见 `decisions.md`
  「依赖精简」与手册 `plotting.md`），出图走 `ferro-py plot`（2026-10-20 前是 `scripts/plot_*.py`）
- **不留兼容层**：输出格式同期变更，留着 `fe-traj` 会让旧脚本「跑成功」却吐出自己
  解析不了的 csv —— 静默坏数据比命令消失难查

### ferro-python（PyO3 绑定 + `ferro.jobs`）

计算输入只由这里生成：Rust 的 `ferro job`、`ferro-workflow` crate 与 `ferro-core` 的 `spin.rs`
已于 2026-10-17 删除（判据见 `decisions.md`「`ferro.jobs`」）。

独立 workspace（pyo3 0.29，**abi3-py311**），maturin 的 Rust/Python 混合布局：Rust 扩展是
`ferro._ferro`，`python/ferro/__init__.py` 整体再导出，`import ferro` 用法不变。
- Rust 侧：`types.rs`（`PyTrajectory`，含 `pbc` / `subset` / `frame_indices` / `energy` /
  `forces` / `stress`）、`io.rs`（`read` / `write` 经 `Format::from_name` 分派，与 CLI 共用 `ferro-io/src/format.rs` 的
  `Format` 表；`read(type_map=)`；
  `read_output`（经 aimd 嗅探，CP2K / VASP / QE 都认）、`write_table`）、`jobs.rs`（`sort_atoms`、`poscar_string`、`qe_structure_cards`、
  `cp2k_structure_sections`、`parse_qe_input_parts`、`write_namelists`）、`structure.rs`、`analysis.rs`、
  `doc.rs`（`docs/src/python/` 七页经 `include_str!` 编进扩展，`_doc_pages` / `_render_markdown`，
  渲染用 `ferro-render`；2026-10-16）
- Python 侧（运行时只用标准库，`plot/` 除外）：`config.py`（`~/.config/ferro/config.toml`，`FERRO_CONFIG` 覆盖）、
  `cli.py`（命令 `ferro-py`：`config init` / `job` / `job template` / `conv` / `conv analyze` /
  `plot` / `doc`，帮助页五段，末行指向 `ferro-py doc <主题>`）、`jobs/`：`base.py`（JobOptions、InputSet、Generator、write_batch、多重度奇偶、
  kmesh）、`vasp.py`（Incar + VaspGenerator）、`qe.py`（QeInput + QeGenerator）、`cp2k.py`
  （Section 节树 + Cp2kGenerator）、`conv.py`（收敛测试）、`templates/<程序>/<任务>/` 9 份内置模板
  、`plot/`（2026-10-20 由 `scripts/ferroplot.py` + `plot_*.py` 迁入：`common.py` 共享层 + `gr` / `sq` /
  `angle` / `msd` / `net` 各一个 `plot()`；可选依赖 `ferro[plot]` = matplotlib + pandas + scienceplots，
  另需 LaTeX；`import ferro` 不加载 matplotlib。样式 `['science','vibrant']`，产物目录共用
  `labels.csv`（`key,label,show`，行序即显示顺序），每张图旁写 `<name>_data.csv`）

测试：`ferro-python/tests/` 的 pytest（个数见文件顶部），本地用 maturin 编 wheel、
`pip install --no-deps --target <scratch>` 后 `PYTHONPATH=<scratch> python -m pytest`。
GitHub Actions（`.github/workflows/python-wheels.yml`）编 manylinux2014 x86_64 与 macOS arm64
wheel，并在 Python 3.11 / 3.13 上用编好的 wheel 跑 pytest，2026-10-16 `903de96` 全绿。
**2026-10-17 起只在推 `v*` tag 或手动触发时跑**（用户定：平时推送构建没有意义）；同一 ref 的新一轮
取消进行中的旧一轮；artifacts 保留 7 天。
CP2K 输入另用 deepmd 环境的 cp2k-input-tools 0.9.1 离线校验过（它与 numpy 2 不兼容，运行时
补 `numpy.cumproduct = numpy.cumprod` 再调；关键字表较旧，不认 `IGNORE_CONVERGENCE_FAILURE`
与 `GAPW_ACCURATE_XCINT`）。三个程序的真实运行验收：超算冒烟 9/9 通过（2026-10-16 第二轮；VASP 6.5.1、QE 7.5、CP2K 2026.2；
第一轮的两处失败与修复见 `decisions.md`「`ferro.jobs`」）。

---

## 依赖版本

| 依赖 | 版本 | 声明位置 | 备注 |
|---|---|---|---|
| `nalgebra` | 0.35 | workspace | 0.34→0.35 零代码改动；g(r)/CN 输出逐字节一致 |
| `ndarray` | 0.17 | workspace | |
| `rayon` | 1.8 | workspace | |
| `rustfft` | 6.4 | workspace | msd 的 FFT 自相关（`md/correlate.rs`）；+4 crate，纯 Rust |
| `thiserror` | 2.0 | workspace | |
| `anyhow` | 1.0 | workspace | |
| `rand` | 0.10 | workspace | 0.8→0.10 改名三处，见 `issues.md` |
| `quick-xml` | 0.38 | ferro-io | vasprun.xml；净新增 1 个 crate，零传递依赖 |
| `ndarray-npy` | 0.10 | workspace | DeePMD npy 读写 |
| `chacha20` | 0.10 | workspace | 可复现的 RNG（`ml/dataset/merge`、`box_builder` 的 `--seed`） |
| `glob` | 0.3 | workspace | `-i` 的通配展开 |
| `libc` | 0.2 | ferro-render | 终端宽度 ioctl；经 rand → getrandom 已在编译图里，净新增 0 |
| `clap` | 4.5 | ferro-cli | derive |
| `pyo3` | 0.29 | ferro-python | 0.21→0.29 仅需 `skip_from_py_object`；2026-10-16 运行时已验证（abi3-py311 wheel + pytest） |

`cargo update` 在主 workspace 与 ferro-python 均已无可更新项。

---

## 已知限制

- **`ferro-python` 是独立 workspace**，被主 workspace 的 `cargo build/test/clippy` 全部跳过：
  改公共 API 后须在 `ferro-python/` 里补跑 `cargo clippy` 与 pytest（2026-10-16 就漏过一次
  clippy 警告）。pyo3 0.29 运行时已验证（2026-10-16，见上节）
- `box_builder`：不支持水合物点记法（`CuSO4·5H2O`）—— 水应作为独立 component 传入；
  入口是 `ferro build box` 与 `ferro.build_box`（2026-10-23）
- **REPL / 脚本模式未实现**。`main.rs` 是可用的子命令分发器，裸 `ferro` 打印分类总览
- 代码库未纳入 rustfmt 管理，格式为手工维护（含刻意的列对齐）；是否采用见 `plan/cli-python.md`
- `cube_density.rs:188` 用参考帧体积归一化，NPT 下同类偏差；需先定义「NPT 下 3D 密度图
  指什么」再动（见 `issues.md`）。**2026-09-12 修 `build_avg_frame` 的 PBC 解缠时
  刻意没碰这个口径**：解缠一律走参考帧的 cell，于是没跨边界的原子逐位不变，
  改动严格局限在那个 bug 上
- **按标签选型只在单帧成立**：`traj gr -x P_3 -y O_b` 对多帧标注轨迹会被 `calc_gr` 的
  逐类型粒子数守恒守卫拒绝（实测 `P_3`：149/152/150/150/150）。多帧请按元素选
- `ferro map chg-sdf` 的 `--cubes` 仍是多 cube 聚合成一张 SDF，与 `map` 其余模式
  「一输入一产物」相反（待办见 `plan/suspended.md`）
- **`ferro bader` 的三个 `.dat` 默认写在输入文件旁边**（2026-09-20 起），名字是
  `bader_ACF_<输入stem>[_<-s 后缀>].dat`（BCF / AVF 同形，`batch::product_name` 的规则），`-o` 可改目录。这是全仓唯一不默认当前目录的
  命令：VASP 的电荷密度一律叫 `CHGCAR`，默认 cwd 时两个体系必然互相覆盖
- `ferro info` 的密度只报首尾两帧（手册 §info 已写明）。元素表里查不到质量的符号
  **报错并点名**（2026-10-20，A-3，用户定 4；此前回退 1 amu 把密度拉低）：`ferro info` 不打
  Density 行、其余照出、退出码 1
- `ferro-python` 仍只暴露 gr/msd，未包 net
- **CP2K / VASP OUTCAR / vasprun.xml / QE pw.x 的输出**已进 `Format` 表（`Format::Aimd(..)`，
  按横幅嗅探，2026-10-17），`convert` / `info` / `traj` 与 `dataset collect` 都能读
- **`--type inspect` 的产物会落进 `-i` 的 glob 射程**：重跑时被当输入收进来，认不出横幅而
  跳过 + 退出码 1，不会静默污染数据（用户须知已写进手册 `dataset/collect.md`）
- **vasprun.xml 的温度是反算量**：该格式不打印瞬时温度（只有输入参数 `TEBEG`），
  由 `<i name="kinetic">` 按 `T = 2·E_kin/(3N·k_B)` 得来，`_info.csv` 的 `#` 头
  会注明。同一次运行读 OUTCAR 与读 vasprun 得到的温度因此**来源不同** —— 与两者
  收敛判据不同、丢帧数不同是同一族问题
- **VASP 的变胞（NPT）路径只经手工构造的文本验证**：用户现有的两份真实数据
  （2000 帧 OUTCAR、425 帧 vasprun）都是定胞 NVT，且用户明确说 VASP 侧不涉及
  NPT。逐帧读胞的代码有测试钉住「第 2 帧的胞来自第 2 帧」，但没有真实变胞数据
- **不支持 ML_FF 的 OUTCAR**（`free energy ML TOTEN` / `ML FORCE`）：无样例可验，
  且 dpdata 对两者用的行偏移不同（14 vs 4），说明差别不止 token 名
- `dataset` 三步（collect / filter / merge）已齐；几何判据只覆盖最小镜像范围（上界逐帧取最紧的一帧），
  未做多层镜像扫描（小胞体系需要时再补）
- **规范序取 (Z, 符号)**（`collect` 自 2026-09-22 起也排，不只 `merge`），
  dpdata 取字母序；两者都靠 `type_map.raw` 自描述，但同一份数据经 ferro 与经
  dpdata 处理，`type.raw` 的数字会不同。
  磁盘顺序对训练没有影响：**DeePMD-kit 加载时自己按类型重排**
  （`deepmd/utils/data.py` 的 `_make_idx_map`，默认 `sort_atoms=True`，
  非 mixed type 的 descriptor 要求如此）。ferro 排完之后 `type.raw` 单调不减，
  那个 lexsort 成为恒等置换

## 时间相关分析的验证记录（2026-10-21 自手册迁入，瘦身 D-6）

手册各页的「Verification」只留一句结论与给用户复现的 numpy 片段；逐项检查与测试名原文在此。
「by hand」各行用的是当时的示例轨迹与独立 numpy 实现。

### `msd`

What has been checked, so a doubtful result can be narrowed down:

| Check | Where | Agreement |
|---|---|---|
| FFT sum vs direct $O(N^2)$ sum, NPT triclinic box with boundary crossings | unit test `test_fft_matches_brute_force_npt_triclinic` | < $10^{-9}$ Å² per lag and axis |
| same, non-periodic | `test_fft_matches_brute_force_nonperiodic` | < $10^{-9}$ Å² |
| straight-line motion across a skewed triclinic face | `test_tor_crosses_skewed_axis` | $(v\,m)^2$ to $10^{-9}$, $v$ the speed |
| TOR vs lattice view under a breathing box | `test_tor_does_not_follow_lattice_scaling` | the table above |
| components add up to the total | `test_components_sum_to_total` | $10^{-12}$ |
| two example LAMMPS trajectories (NPT, NVT) vs an independent numpy implementation | by hand, 2026-09-27 | within the csv's 7 significant digits |

### `vacf`

| Check | Where | Agreement |
|---|---|---|
| FFT vs direct all-origin sum, irregular velocities, per axis | unit test `test_fft_matches_brute_force` | $10^{-12}$ |
| trapezoidal integral (constant $C_v$ gives $c\,t/3$, 0 at $t=0$) | `test_diffusion_is_trapezoidal` | $10^{-12}$ |
| $C_v(0)=0$ leaves `vacf_norm` empty | `test_zero_velocity_gives_nan_norm_not_zero` | — |
| synthetic AR(1) velocities, 30 atoms × 200 frames, vs an independent numpy implementation | by hand, 2026-09-27 | within the csv's 7 significant digits, all columns |

### `rotcorr`

| Check | Where | Agreement |
|---|---|---|
| FFT vs direct average, every frame valid (the `gmx rotacf` case) | unit test `test_fft_matches_brute_force_all_valid` | $10^{-10}$ |
| same with one frame in five invalid per molecule | `test_fft_matches_brute_force_with_invalid_frames` | $10^{-10}$, NaN at the same lags |
| lags with no valid pair are NaN | `test_lag_without_valid_pair_is_nan` | — |
| trapezoidal integral | `test_integral_is_trapezoidal` | $10^{-10}$ |
| rigid regular $PO_4$ rotating about $z$: `sum` reports cancellation, `bond` gives the analytic $C_2$ ($\cos = \tfrac13 + \tfrac23\cos\theta$ for every bond) | `test_sum_cancels_on_tetrahedron_bond_does_not` | $10^{-10}$ |
| `bond` + $P_1$ vs direct average | `test_p1_bond_matches_brute_force` | $10^{-10}$ |
| `bond` keeps following a bond stretched past $r_\text{cut}$ | `test_bond_is_followed_beyond_r_cut` | — |
| example glass, P–O `bond` mode, $P_1$ and $P_2$ (1488 bonds = 372 P × 4), vs independent numpy | by hand, 2026-09-27 | $5\times10^{-8}$ |
| example glass trajectory, P–O at $r_\text{cut}$ = 1.8 Å (all valid) and 1.52 Å (0.5 % invalid), vs an independent numpy implementation | by hand, 2026-09-27 | $5\times10^{-8}$ (csv precision), NaN positions identical |

### `bondlife`

| Check | Where | Agreement |
|---|---|---|
| $C_I$ and $S_C$ vs counting every origin directly, irregular bonding | unit test `test_both_functions_match_brute_force` | $10^{-12}$ |
| break-and-reform counts for $C_I$, not for $S_C$ | `test_intermittent_counts_reformed_continuous_does_not` | exact |
| hysteresis stops flicker; a free pair needs $r \le r_\text{bond}$ | `test_hysteresis_suppresses_flicker`, `test_hysteresis_needs_r_bond_to_form` | exact |
| gap filling: $k$ filled, $k+1$ not, ends never | `test_fill_gaps_boundaries` | exact |
| $n(t) = n(t-1) + \text{formed} - \text{broken}$ | `test_events_conserve_bond_count` | exact |
| example glass, Zn–O (single threshold 2.4 Å; hysteresis 2.3/2.6 Å with intermittency 2) vs an independent numpy implementation | by hand, 2026-09-27 | $5\times10^{-8}$ (csv precision); events identical frame by frame |
