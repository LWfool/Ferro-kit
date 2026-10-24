# 待办：ferro-io（读写、格式、测试 fixture）

> 从 `dev/plan.md` 按需读取；定级规则见 `plan.md` 顶部。

## 优先级高

### ferro dataset：剩余小项（2026-08-25）

**额外键搬运**：`atom_ener` / `fparam` 这类 `Frame` 装不下的项，读时告警、
写时丢失。真出现时再设计（需要一条绕过 `Trajectory` 的按帧索引搬运通道）。
（原节其余各项 —— NEP 导出、filter 读回 npy、merge 的 `type_map` 对齐 —— 均已完成。）


## 优先级中

### 带速度的测试 fixture（2026-09-12 提出，2026-09-21 仍未做）

`tests/` 两条 fixture **没有速度数据**，于是 `vacf` / `vanhove` / `rotcorr`
跑不到写文件那一步 —— 这三个命令端到端从来没验过（0.2.0 改名那轮就只验了
gr/angle/msd/sq/net/map），2026-09-12 的 golden master 基线也盖不到它们。

补一条带 `vx vy vz` 列的 LAMMPS dump（5 帧即可）就能把三个命令纳入。
注意**不能**用全零速度顶替：VACF 在全零速度下是退化情形，跑得通但看不出
实现对不对，那是假的安全网。

### VASP：ML_FF 的 OUTCAR（2026-08-27 提出）

VASP 的机器学习力场跑出来的输出用 `free  energy ML TOTEN`
与 `ML FORCE`。dpdata 用一个 `ml=True` 开关切 token，但它给两者的行偏移是
`[14, 4]` —— **块结构本身就不同**，不只是换个名字。没有样例就没法验，等有输出
再补；`vasp_outcar.rs` 的锚点已经是多候选表的形状，加 token 不必重构。
（原节的另一半「`io_dispatch` 注册」已由 10-17 第 3 步的 `Format` 表解决。）

## 优先级低

### DeePMD mixed type 数据的读写（2026-08-26 提出）（2026-10-22 由高降为低：功能类，按定级规则；瘦身 E-12 / E-20，用户定）

**现状（2026-10-07，审查 A6）**：reader 见到任一 set 有 `real_atom_types.npy` 即报错拒收。
**只拒不读**是有意的：只做读侧的话 `filter` 读进来再写回标准布局（按第一帧定 `type.raw`），
mixed 数据会静默降级成坏数据、退出码 0 —— 比拒收更糟。要支持就读写与下列判据一起做。

现在 `readers/deepmd.rs` / `writers/deepmd.rs` 只做**标准 system**：`type.raw` 一份
定型，故一个 system 里所有帧的成分必须完全相同 —— 这正是 `collect`「同目录成分不符
报错」与 `merge`「按 `composition_key` 分组」两条现有约束的来源。

mixed type 布局（**先核对 DeePMD-kit 文档与 dpdata 的 `deepmd/npy/mixed` 再动手**，
以下是待验证的理解）：

| 文件 | 标准 system | mixed type |
|---|---|---|
| `type_map.raw` | 该 system 的元素表 | 全局并集 |
| `type.raw` | 逐原子真实类型 | **全 0 占位** |
| `set.NNN/real_atom_types.npy` | 无 | `(nframes, natoms)` 的整型，逐帧逐原子给真实类型 |

要点与判据：

- **价值是装下成分不同的帧**（原子数仍须相同 —— dpdata 的 mixed 也按 natoms 分
  system）。所以这件事做完，`collect` 的「成分不符报错」和 `merge` 的分组要重新定：
  是继续分组、还是给一个 `--mixed` 让它们合成一个 system。**默认不改**，因为
  mixed type 只有 DPA 系列 / 多任务训练吃得下，普通 DeePMD 训练不认
- 读侧先做：能读回 dpdata 产出的 mixed system，`real_atom_types` 经 `type_map` 映射
  回元素符号，落进 `Trajectory` 天然装得下（帧与帧的元素本来就各存各的）
- 写侧作为开关，不改默认布局。`type_map` 的并集需要**可复现的稳定序** ——
  `merge.rs` 的 `(Z, 符号)` 规范序已满足，沿用同一条（注意与 dpdata 的字母序不同，
  这条差异 `progress.md`「已知限制」已记）
- npy 是整型：现有读写走的都是 `float64`（磁盘上一律二维 f64），
  `real_atom_types.npy` 是 int，**dtype 分支是新的**，读侧要同时收 int32/int64
- **`sort_atoms` 不能无脑套用**（2026-09-22 查证）：DeePMD-kit 在 mixed type 下
  明确**不**按类型排序（`deepmd/utils/data.py` 的 `sort_atoms` 文档写着
  "except mixed types"），因为 `real_atom_types` 是逐帧各异的。而 `collect`
  现在无条件排 —— 两者相遇时要先定 mixed system 的规范序是什么
- 与 `filter` 的关系：`filter` 现在按 system 读回再写出，mixed system 经它一趟必须
  仍是 mixed，否则静默降级成「全 0 类型」的坏数据

### CP2K 的 EXTXYZ / XMOL 输出：kind 名进 species 列、`Energy=` 是 Hartree（2026-08-26 提出，2026-10-22 查源码，降为低）

**降级理由**（2026-10-22 用户定）：用户的 CP2K 数据全从 `.out` 读（`cp2k_md` / `cp2k_sp`），不开
`MOTION/PRINT/* FORMAT XMOL|EXTXYZ`。下面两条只在主动开了这些输出时才碰得到。

**CP2K 的 EXTXYZ 是一量一文件**（CP2K master `motion_utils.F`、`input_cp2k_motion_print.F`）：

| print key | `Properties` | 数值默认单位 |
|---|---|---|
| `TRAJECTORY`（`-pos-1.xyz`） | `species:S:1:pos:R:3` | Å（随 `UNIT`） |
| `VELOCITIES`（`-vel-1.xyz`） | `species:S:1:velo:R:3` | bohr/au_t |
| `FORCES`（`-frc-1.xyz`） | `species:S:1:force:R:3` | Hartree/bohr |
| `STRESS` | 无 `FORMAT` 关键字，不出 EXTXYZ | — |

标题行恒为 `Lattice="…"（Å） Properties=… pbc="…" Step=N Time=t（fs） Energy=E`。Ferro 只读得进 `-pos`
文件（另两种无 `pos` 列，reader 报错，不会读错）。两处会静默出错：

1. **`Energy=` 是 Hartree**：`write_trajectory` 把 `etot` 原样写出（晶胞那行乘了 `angstrom`，能量没换），
   Ferro 按 extxyz 惯例当 eV 读，差 27.2 倍
2. **kind 名进 species 列**：`PRINT_ATOM_KIND` 为真、或 kind 无元素符号时，`particle_methods.F` 的
   `write_particle_coordinates` 写 `&KIND` 名原文，标题行无任何标记。XMOL（走 xyz reader）同样受影响。
   实测（2026-10-22）：species 为 `Fe1` / `O_s` 的文件读入退出码 0，`convert` 出的 POSCAR 元素行就是
   `O_s   Fe1`，gr 按它们分组；ASE 读同一文件 `KeyError: 'Fe1'`

**要做时的候选修法**（当时推荐，未经用户定）：

- 1：按原始大小写认 CP2K 标题特征（`Step=`、`Time=`、`Energy=` 同时出现；ASE 写小写 `energy=`），
  换算 Hartree → eV（`units.rs` 的 `HARTREE_TO_EV`）并告警；备选是丢能量并告警、指向 `.out`
- 2：整串是元素照旧（合规文件逐字节不变）；否则查元素表推断（开头全部字母 → 前两个 → 首字母），原串
  进 `label`，告警列 `原串 -> 元素`，推不出报错。合乎「reader 缺元素时从标签推、如实匹配」（`CA1` →
  Ca 同 `Os` → 锇）。改的是 xyz / extxyz 共用的 `readers/util.rs::species_element`，两个 reader 一并生效。
  **不能照抄** dump 的「无条件按下划线拆」：合规 extxyz 的 species 本就是纯元素（`Pb` 贪婪前缀的教训见
  `issues.md`）
- 无法识别、只能进手册的：`TRAJECTORY%UNIT bohr` 时坐标是 bohr、晶胞仍是 Å

更一般的一条同属此节：xyz / extxyz 的 species 不在元素表里时一律照收（不只 CP2K 的文件），修法即上面第 2 条。

### `convert` 写 CIF 时判定并写出空间群（2026-10-11 提出，功能，低）

现状：空间群只在**读** CIF 时生效（`readers/cif.rs` 按 symop 或空间群符号展开），**写侧**
（`writers/cif.rs:9,52-53`）一律写 `'P 1'` / IT 1、全部原子显式列出，无开关。实测 NaCl（Fm-3m，
两个位点）经 `ferro convert` 一圈变成 P1、8 个原子 —— 结构正确，对称性信息丢失。

要做的是**对称性判定**：从坐标反推空间群与非对称单元，再写 `_space_group_symop_operation_xyz`
与只含非对称位点的 `_atom_site_` 表。ASE / pymatgen 都靠 spglib；Rust 侧候选 `moyo`（spglib 的
纯 Rust 重写）或 spglib 的 C 绑定 —— **新依赖，动手前先问用户**。要先定的口径：

- 容差（`symprec`，spglib 默认 1e-5 Å，pymatgen `CifWriter` 常用 0.1 Å）及是否暴露成参数
- 默认行为：仍写 P1（逐字节不变）而由开关打开判定，还是默认判定
- 标准化设置：写输入原胞的设置还是标准设置（后者会旋转 / 平移晶胞，坐标变）
- 部分占位、带标签（`label`）的位点如何参与判定
- 回读校验：写出后用本仓 reader 展开，原子数与坐标应与输入一致（可与 `scripts/crosscheck_spacegroups.py`
  对 ASE 的对拍复用）

### 深度学习工作流 I/O

MACE/NequIP 兼容格式仍未开始。（DeePMD npy 读写与 GPUMD/NEP 的 `train.xyz` 导出已完成。）

## 2026-10-21 从各审查节归拢的条目（定级沿用原节，见各条）

### 测试缺口（遗留清理第 8 批，原属高；2026-10-20 用户定暂缓）

历次审查登记的测试缺口原文如下，其中一部分可能已随后续修复补上，开工前逐条复核。带速度的 dump fixture 与下面「带速度的测试 fixture」是同一件事。

- 10-11：（2026-10-14 删去已随修复补上的项）：Bader near-grid / off-grid 对解析密度的收敛；
带速度的 dump fixture（C 用合成数据验过 vacf，仓库里仍无）。

- 10-03：（2026-10-09 删去已随修复补上的项）：workflow 只有子串断言、无一条对照 CP2K / pw.x
规则（M4–M12 全漏，随 QC 专项补）；`help_sync` 的 PAGES 不含 job 三页（约 35 个参数无防漂）；
CLI 集成测试只有 bader 与 npt；gr 三斜数值（可用 27 镜像暴力对拍）；cube_sdf ≥8 P 的运行时上界与
组内置换语义（随 B-6）；vanhove / angle 末 bin；rotcorr `--legendre 1` 头部；bondlife 对 MDA 定义
钉死；spin 锕系 / 未知元素；weight BCF 与 ACF 一致性、ACF 用 ASE / pymatgen 回读（随 Bader 重写）；
`wrap_position` 恰返 L；三斜下 `classify_frame`；dataset merge 乱序输入对 dpdata 端到端；
`--last-n 0`；分派大小写、Python 写 CONTCAR；Python 绑定零测试（D-M1 只用临时 wheel 实测过）。

- 10-01：（与上表对应）：dump writer 无「非下三角胞 + 力/速度」；POSCAR 无「元素交错 +
速度」；QE / CP2K 只测了 `{}` 与 `.TRUE.`；lammps data 无「无注释 Atoms」；无截断 dump
fixture；cube 往返只测零原点；cube_radius 无窗口超网格的用例。

**审查未覆盖**：Bader 全套（已登记项之外）、`chg_sdf`、`cube_sdf` 聚类、`cube_jump`、
`cp2k_sp`、pdb / cif writer、doc 渲染器、
ferro-python 运行时。（`--metal-units` 全链路已随 2026-10-03 的 `--units` 改造走过一遍）

### 其他 io 待办（2026-10-21 从各审查节归拢）

- **`read_cube` 与 `read_cube_as_chg` 合并成一个 cube reader**（10-17 第 5 步用户定：`read_cube` 标记保留，以后合并）
- **`ferro dataset` 审查（2026-10-17）未做的低级项**：D-5 filter 不重排原子（collect / merge 会），可写进手册；D-7 超过 1000 个 set 时目录按字典序；仍存疑：左手系晶胞
