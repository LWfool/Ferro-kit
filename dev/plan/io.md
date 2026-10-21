# 待办：ferro-io（读写、格式、测试 fixture）

> 从 `dev/plan.md` 按需读取；定级规则见 `plan.md` 顶部。

## 优先级高

### DeePMD mixed type 数据的读写（2026-08-26 提出）

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

### CP2K 的 EXTXYZ 把 atom kind 写进 species 列（2026-08-26 提出）

CP2K 新版的 `MOTION/PRINT/TRAJECTORY` 多了 `FORMAT EXTXYZ`，而它的
`PRINT_ATOM_KIND` 会**把 subsys 里的 atom kind 写进 species 列**（文档原文：只对
XMOL 与 EXTXYZ 有效）。Ferro 的 extxyz reader 假设 species 是纯元素、位点标签走
独立的 `label:S:1` 列，撞上这种文件会把 kind 当元素收下。

与 LAMMPS dump「没地方放第二个名字只能折进 element 列」是同一族问题，但**不能照抄
那边的解法**：dump 那边是无条件按下划线拆，而 extxyz 的 species 列在合规文件里就
该是纯元素，无条件拆会误伤。要定的是判据 —— 拆还是不拆、按什么拆、
`split_element_label` 的 `Unknown` 分支怎么处理（`Pb` 那类贪婪前缀误判的教训见
`issues.md`）。

CP2K 的 EXTXYZ **只写 cell + 坐标**，不写 stress/virial（力在 `PRINT/FORCES` 另一个
文件里），故这条与应力无关，是纯粹的标签映射问题。


## 优先级低

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
