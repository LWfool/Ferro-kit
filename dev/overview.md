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
    ├── ferro-core              ← 纯数据结构 + 静态参考数据 + spin（未成对电子推断）
    │                             + Table（分析产物的中立载体）
    │                             + AtomType（网络分类结果，结构化而非字符串）
    ├── ferro-io       → core   ← 格式读写 + write_table（分析产物的唯一出口）
    ├── ferro-structure→ core   ← 超胞、真空层、合并、初始盒子
    ├── ferro-analysis → core   ← 纯计算，不碰文件系统；结果暴露 to_tables()
    │                             md/、network/、dft/（Bader、ChgSDF）、
    │                             ml/（数据集筛选与合并）
    └── ferro-workflow → core   ← QC 输入生成（Gaussian / CP2K / QE，含基组赝势库）
```

中间层 crate 之间不允许互相依赖。**但共享类型可以下沉** —— 判据：一个类型该放
`ferro-core`，当且仅当两个以上中间层需要叫出它的名字。两个方向各一个范例：
`Trajectory`（io 产出、analysis 消费）与 `Table`（analysis 产出、io 消费）。

这条规则的完整论证见 `issues.md`「分析产物为什么不进 ferro-io」—— ASE 与 pymatgen
是 Python，**能**写出循环依赖却依然不写，说明这是设计选择而非语言约束。

## 命令入口（0.2.0 起；`net` 于 0.2.1 降为叶子命令）

**单二进制 `ferro`**，子命令按**产物**分组（不是按实现它的 crate）：

```
ferro traj  gr | sq | msd | angle | vacf | rotcorr | vanhove   → 堆叠 csv + 可选 PNG
ferro map   density | velocity | force | radius | sdf | chg-sdf → 逐输入一个 .cube
ferro net                                                      → 六张堆叠 csv
ferro dataset collect | filter | merge                         → DeePMD system 目录
            collect --type inspect                             → 诊断三件套（不出数据集）
ferro bader | convert | info | job
ferro doc   <topic>                                            → 手册（编译在二进制里）
```

`dataset` 的产物是**目录**（DeePMD 的 system 就是目录），故它的 `-o` 是输出根目录。
0.3.3 把这一条推成了全仓规则，它不再是例外。

`-i` 恒为多值并自展开 glob；逐文件独立分析，结果堆叠成一份带 `file` 列的 csv。
**`-o` 恒为路径**（0.3.3 起）：写多个产物的命令指目录，`convert` / `job` 指文件；
批次标记走 `-s/--suffix`。目录缺失时先问 `[y/N]`，非交互环境须给 `--mkdir`。

参数与输出列结构见 `docs/src/cli-reference.md`。

## 编码规范

见 `CLAUDE.md`「编码约定」：手册与 doc 注释用英文、`//` 内部注释用中文、
公式一律 LaTeX `$...$`。**不在两处各写一份** —— 这段此前就是 CLAUDE.md
那几行的副本。

## 版本规则

版本号在根 `Cargo.toml` 集中管理，各 crate 继承 `workspace.package.version`；
`ferro-python` 是独立 workspace，需手动同步。

**只在用户明确要求时才动版本号**，不要每次改代码就自动 +1。破坏性改动升次版本位
（0.1.15 → 0.2.0），其余升 patch 位。发布点打 annotated tag。

**已知例外**：`v0.2.1`（2026-08-12，net 重构）含破坏性改动（`net qn`/`net type`
删除、标签格式全变、表名与列结构全变），按上面的规则本应是 0.3.0，但按用户当时的
明确要求走了 patch 位。翻 git 历史时注意：**0.2.0 → 0.2.1 之间有 breaking change**。

## `v0.3.0` 的破坏性改动（2026-08-20 发版）

`v0.2.1` → `v0.3.0` 累积了三批破坏性改动，此前按用户要求一直压着不发版，等配套的
Python 绘图脚本跟进后一并升。**旧产物与旧命令行都不兼容**，三批分别是：

**net 重做（2026-08-12）**

| 改动 | 影响 |
|---|---|
| 六张表改名 | `bridge`→`qn`、`partner`→`qn_partner`、`oxy`→`ligand_type`、`cn`→`coordination`、`mean`→**删除** |
| 新增 `composition` 表 | 一物种一行的结构组成，取代 `average` |
| Al 退出 Qn 表 | 非 Qn 形成子只在 `coordination` 出现 |
| 标签数字换义 | 非 Qn 形成子由桥接数改为配位数（`Al_4` = 四配位） |
| 两套标签词汇 | 分布表 `P-Q2`（单元），`linkage` 与导出轨迹 `P_2`（原子） |
| `linkage` 加列 | `linkage`（展示）、`ligand`（配体元素，进键） |
| `ligand_type` 改列 | `label` 合并为 `Al-O_b-P`，原标签移入 `type` |
| 配体 `fraction` 分母 | 全体配体原子 → 该配体元素 |
| 新增 `--qn` | 替换默认 Qn 名单 `{B,P,Si}` |

**traj 产物命名（2026-08-13）**

| 改动 | 影响 |
|---|---|
| 文件名加 label 段 | `gr.csv` → `gr_P-O.csv` / `gr_all.csv`；六个命令全改，**含无筛选时的 `_all`** |
| 新增 `--outdir` | 进 `CommonArgs`，覆盖 11 个命令的 csv/png/cube/导出轨迹；`chg-sdf` 单独加同名参数 |
| `traj sq` 移除 `-a/-b/-x/-y` | 旧命令行直接报 `unexpected argument`；按 label 分辨的 partial 从 CLI 消失 |

**net 口径改为文献 $Q^n_m$（2026-08-20）**

| 改动 | 影响 |
|---|---|
| `qn` 只数同元素连接 | **数值全变**。`P-Q3` 由 40.4% 变 0.27%，`mean_qn` 由 2.40 变 0.95。异核桥进 `m_<X>` 列，总桥 = `n + Σm` |
| `bridges_to` 数连接非桥氧 | 三簇配体贡献 2 而非被跳过，Σm 不再亏空 |
| `qn_partner` 去掉自身元素列 | `m_P` 恒等于 `qn`，重复列 |
| `linkage` 列改名 | `n_bridge_a/b` → `qn_a/qn_b`，值改为同核连接数 |
| `[inputs]` 新增 `mean_n_bo` | 桥氧个数（旧口径那个数）与 `mean_qn` 并列 |
| 新增共边告警 | 形成子共享 ≥2 配体时提示核对 cutoff |

依据是 arXiv 2510.13545 等文献原文：n 只数同核桥，m 数异核桥，总桥由两者相加；
铝磷酸盐文献专门指出「n 是总桥氧数」是领域内的已知误解。判据与排查清单见
`issues.md`。

三批合起来按规则升次版本位，故 `0.2.1 → 0.3.0`（跳过 0.2.x）。配套的
`scripts/plot_net.py` 已同步（`a87843d`）。

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
| **单位从输出文本自读**，认不出报错不默认 —— CP2K 的 `STRESS_UNIT` 是输入关键字，单位不是版本的函数 | `readers/cp2k_out.rs` |
| **AIMD 的能量取自由能**（CP2K 的 `FORCE_EVAL`、VASP 的 `free energy TOTEN`），不取 `sigma->0` —— 力是自由能对坐标的导数，配 `sigma->0` 等于给模型两半不同的泛函 | `readers/vasp_outcar.rs` |

`v0.3.1` 是纯新增，故升 patch 位。

## `v0.3.2`（2026-08-26，未发版）

**含破坏性改动但按用户明确要求走 patch 位** —— 与 `v0.2.1` 同一类例外。翻 git
历史时注意：**0.3.1 → 0.3.2 之间 `dataset collect` 的产物布局不兼容**。

| 改动 | 影响 |
|---|---|
| `collect` 分组单位由**文件**改为**目录** | 同目录的 `.out` 合并成一个 system（它们是同一次运行被重启切开的段）。这也是 collect 与 merge 的新分界：collect 拼**同一次运行**的碎片，merge 合**不同运行** |
| `collect` 命名保留目录层级 | 剥掉公共祖先后原样嵌套，文件 stem 不进名字。`run1_total/` → `run1/`；只有一组时祖先是整条路径，直接写进 `-o` 本身 |
| `collect` 的 `-o` 改为必填 | 默认 `.` 会把 npy 撒进正在工作的目录 |
| `collect` 新增 `--overwrite` | 与 `filter` / `merge` 一致 |
| `filter` 的报告落盘 | 七张 csv 平铺在 `-o` 根下（`filter_*.csv`）。诊断表改为恒算 |
| 批内失败提示改措辞 | 不再承诺「见 `[inputs]` 块」—— dataset 的产物是目录，没有那个文件 |
| **extxyz 的 `stress=` 变号** | 写出的九个数与 0.3.1 逐个反号。旧行为不是「另一种约定」而是错的：extxyz 的 `stress=` 是 ASE 约定（正 = 拉伸），`Frame::stress` 是正 = 压缩，两侧都没变号。读侧同时改：`virial=` 不再与 `stress=` 混作一谈（差一个体积因子），非对称张量与 6 分量 Voigt 一律拒收 |

**新增**：`ferro doc`（手册编译进二进制，`include_str!` 24 页 208 KB，topic 跟
子命令树同名，tty 下走 `$PAGER`）。

**新增**：VASP AIMD 读取 —— `readers/vasp_outcar.rs` 与 `readers/vasprun.rs`，
`collect` 改为**按文件头几行的横幅嗅探**格式（VASP 写的叫 `OUTCAR` 无扩展名，
`.out` 又太通用，按名字判必然开一串特例）。`Cp2kOutStats` 下沉改名为
`AimdStats`（第二个使用者到齐，符合 CLAUDE.md 的下沉判据）。新增依赖
`quick-xml`（净新增 1 个 crate，流式不建 DOM —— AIMD 的 vasprun 常有几百 MB）。

**新增**：`dataset filter` / `merge` 的 `--type deepmd|nep|extxyz` 与
`--ratio 8:1:1`（train:valid:test，两段即 train:test；默认不划分，旧命令行行为不变）。GPUMD/NEP 的
`train.xyz` 就是 extxyz，故导出**没有新命令** —— 缺的只是「读 system 目录」，
而 filter/merge 本来就在读。划分产物沿用 dpgen 的目录名后缀 `.train/.valid/.test`。

**新增约定**：`Cp2kOutStats.steps` 带出文件的 step 区间；重复帧**不去重**（重启
重叠段位置速度相同，能量力也相同），但区间打出来让重叠可见。

## `v0.3.3`（2026-09-20，未发版）

**含破坏性改动但按用户明确要求走 patch 位** —— 与 `v0.2.1`、`v0.3.2` 同一类例外。
`v0.3.2` 亦未单独发版，两批改动会一并进入下一个发布点。

| 改动 | 影响 |
|---|---|
| **`-o` 恒为路径，`--outdir` 删除** | 13 个分析命令加 `chg-sdf` 的 `-o` 由「文件名后缀」变为「输出目录」；旧命令行 `-o cmp` 不报错，但含义变成写进 `cmp/` 目录 |
| **批次标记改 `-s/--suffix`** | `dataset filter` 的 `-s` 仍是 `--s-max`，`job` 的 `-s` 仍是 software（两者不收批次标记，不冲突） |
| **新增 `--mkdir` 与创建确认** | 目录不存在时向 stderr 问 `[y/N]`；**非交互环境（脚本、CI）不给 `--mkdir` 会报错退出** |
| `convert` 的 `-o` 可带目录 | 以分隔符结尾则报错（格式从文件名推断）；`job` 的 `-o` 以分隔符结尾则放默认名 |
| **`bader` 报告改写到输入旁边** | 默认位置由当前目录改为 `-i` 所在目录；新增 `-o` / `-s` |
| **`dataset` 产物默认带 `.train`** | `filter` / `merge` 的输出目录名多出 `.train`（`collect` 不加）；`merge` 遇混合后缀报错；`.valid`/`.test` 不能再划分 |
| `dataset` 的 `--outdir` 改名 `--output` | 短名 `-o` 不变 |

**内部**：`ferro-io` 的 reader/writer 路径参数由 `&str` 统一为 `&Path`；
`BaderResult` 的三个 `write_*` 改为 `*_text()` 返回字符串，落盘移到 CLI ——
`ferro-analysis` 不再碰文件系统。新增 `tests/CHGCAR_2atoms` 与
`ferro-cli/tests/bader_reports.rs`。

## 2026-09-21 的一批（版本号**仍是 0.3.3**，未发版）

**含破坏性改动，但按用户明确要求不升版本号** —— 这一批没有自己的版本号，
它与 `v0.3.2`、`v0.3.3` 挤在同一个 `0.3.3` 里，一并进入下一个发布点。
翻 git 历史时注意：**`0.3.3` 这个版本号下有三批改动，且都是破坏性的**。

起因是把 `private/cp2k_grep_strInfo.py` 的能力并进来。该脚本的九项能力里 Ferro
已有六项且实现更稳健，真正合进来的是 **LAMMPS 导出 + 逐帧标量序列**；CP2K 的
多文件布局与 DeePMD mixed type 明确划在范围外（仍在 `plan.md`）。

| 改动 | 影响 |
|---|---|
| **LAMMPS dump 三斜盒子行改 `*_bound`** | 读写两侧口径都变。正交胞产物**逐字节不变**；三斜胞此前写出的盒子偏小，OVITO / ASE / LAMMPS `read_dump` 读到的都是错的 |
| **`collect` 的 `-o` 由必填改可选** | 不给时产物落在 AIMD 目录**同级** |
| **`collect` 的产物目录名恒带 `.db`** | `sets/run1/` → `sets/run1.db/`；`filter`/`merge` 会剥掉它再命名自己的产物 |
| 新增 `collect --type inspect` | 只出诊断三件套，不出数据集；`-o` 在这条路上报错 |
| `Frame` 加 `temperature` / `step` 两个字段 | 三条 AIMD 路都填（vasprun 的温度是反算的，无步号） |
| extxyz 读侧收 `Temperature=` | 写侧**不写** |

推翻的旧判据两条，都写进了 `issues.md`：

- **「`collect` 的 `-o` 必填」**（0.3.2 定）。原判据反对的是默认值 `.`（会把 npy
  撒进正在工作的目录），而跟着输入走的默认值撒不到别处去。结论改为：**默认值必须
  跟着输入走，不能是 `.`**
- **「`collect` 不加后缀」**（0.3.3 定，理由是产物是没筛过的原始数据、不该标成
  训练集）。改为加 `.db`（database），但它**不是**划分后缀，不进 `SPLIT_SUFFIXES`

## 2026-09-22 的一批（版本号**仍是 0.3.3**，未发版）

起因是 `private/cp2kdata`（robinzyb/cp2kdata，MIT），核心需求是**单点能数据的
读取**——把一批 FP 计算的输出收成训练集。范围只取能量 / 力 / 坐标 / 盒子 /
应力五项，与 MD 端一致；cube、pdos、Mulliken / Hirshfeld 布居、振动频率、偶极、
DFT+U 占据数**全部不取**（要么 Ferro 已有更稳健的实现，要么 `Frame` 里没有存放
处）。GEO_OPT / CELL_OPT 也不取——cp2kdata 里那两条路的坐标根本不来自 `.out`，
而是去 glob `*pos*.xyz`，属于已划在范围外的 CP2K 多文件布局。

| 改动 | 影响 |
|---|---|
| **`cp2k_out.rs` 改名 `cp2k_md.rs`** | `read_cp2k_out*` → `read_cp2k_md*`；`AimdFormat::Cp2kOut` → `Cp2kMd`；fixture 同步改名。纯机械，行为零变化 |
| 新增 `readers/cp2k_sp.rs` | CP2K 单点（`ENERGY` / `ENERGY_FORCE`）reader，`AimdFormat::Cp2kSp` |
| `sniff` 按 `GLOBAL| Run type` 二次分派 | 同一个 `CP2K|` 横幅下两套布局，此前一律当 MD |
| `AimdStats` 加 `version` / `version_note` | 版本判定在 reader，打印在 CLI |
| **MD 路补两处老版本的块** | `energy (a.u.):` 的圆括号写法、无 `STRESS|` 前缀的 ` STRESS TENSOR [GPa]` |
| `collect` 报告加两行 | 版本提示、CP2K kind 名映射 |

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

## 2026-09-22 的第二批（版本号**仍是 0.3.3**，未发版）

起因是实测 `ferro dataset collect -i 1Al/*.out` 报「different compositions」而
两个化学式打出来一模一样。根因是 collect 用**逐原子原始顺序**比较成分、却用
**计数**渲染报错，而 CP2K 按输入文件列原子的顺序写坐标表，同一体系的两次单点
计算原子顺序常常不同。

| 改动 | 影响 |
|---|---|
| **collect 逐文件规范排序** | 产物的 `type.raw` / `coord.npy` / `force.npy` 原子顺序与旧产物不同。同成分不同顺序的文件从此能合并 |
| `sort_atoms` 补 `velocities` / `bonds` 置换 | 此前只置换 `forces`，另两个静默错位。merge 路看不见（system 不存速度与键），collect 的诊断导出会看见 |
| 成分报错点名差异元素 | `differing: Al 2 vs 3`，不再只并排两个化学式 |
| **新增 `collect --format`** | `cp2k/md` / `cp2k/sp` / `vasp/outcar` / `vasp/xml`。嗅探仍跑，降为出 NOTE 与喂守卫 |
| `sniff` 对 CP2K 缺 run type 改报错 | 此前默认按 MD 走，一份被裁过头的单点日志会报成「MD 找不到锚点」 |
| **25 页帮助全部精简为五段** | `collect` 96 → 28 行。删掉的是判据与口径，手册里逐条核对过都已有 |
| 手册新增 `changelog.md` | `ferro doc changelog`。从 v0.3.0 起，只写「旧 → 新」命令行对照 |

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
`plan.md` 的 mixed type 待办已记下这一条。

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
