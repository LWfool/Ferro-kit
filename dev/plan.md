# 后续计划

> 归档只保留**判据**与**被实测推翻的原计划** —— 「做了什么、怎么改的」翻 git 历史
> （提交号已列出），「现在是什么样」看 `progress.md` 与 `docs/src/`。

## 优先级高

### GPUMD/NEP 的 `train.xyz` 导出（**最高优先**，2026-08-26）

前置的 extxyz 符号/体积因子修正**已完成**（见归档），导出本身仍待做。这轮查证顺带
把导出要用的约定全部落实了，下次开工不必重查：

| 事实 | 出处 |
|---|---|
| `virial=` 单位 eV，**正 = 压缩**，等于 `Frame::stress × V` 不变号 | GPUMD 手册；与 dpdata 1.0.2、extxyz 规格三方一致 |
| `stress=` 单位 eV/Å³，**正 = 拉伸**（ASE 约定），写出时变号 | 同上 |
| 两键同在时 **GPUMD 取 virial** | GPUMD 手册。注意 Ferro 的 reader 更严：不一致直接报错 |
| 力列 `force:R:3` 与 `forces:R:3` **都合法** | GPUMD 手册；读侧已两种都收 |
| `lattice="ax ay az bx by bz cx cy cz"` | 与 Ferro 现有写法一致，不必改 |

仍要定的：

- **命令归属**：`ferro dataset export --format nep`（与 collect/filter/merge 同链）
  还是 `ferro convert` 的一个目标格式。判据：产物是**一个文件**而非目录，更像
  convert；但输入是 DeePMD system 目录，而 `convert` 的 `-i` 现在不收目录。倾向前者
- **写 `stress=` 还是 `virial=`**：通用 extxyz writer 定的是只写 `stress=`（两个键
  = 两处可能矛盾的事实）。NEP 侧两个都认，故沿用 `stress=` 即可，除非实测发现
  GPUMD 对 `stress=` 的处理有别
- `config_type` / `weight` 这类 NEP 可选键：先不写，需要时再加

---

### DeePMD mixed type 数据的读写（2026-08-26 提出）

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
- 与 `filter` 的关系：`filter` 现在按 system 读回再写出，mixed system 经它一趟必须
  仍是 mixed，否则静默降级成「全 0 类型」的坏数据

---

### VASP AIMD 数据读取（2026-08-26 提出）

`readers/vasp.rs` 现在只有 `read_poscar` / `read_contcar`，AIMD 轨迹无入口；
`dataset collect` 也是直接调 `read_cp2k_out_with_stats`，**没有格式分派**。
这条要同时补 reader 和 collect 的分派。

**来源三选一**：

| 来源 | 内容 | 代价 |
|---|---|---|
| `OUTCAR` | 能量 / 力 / 应力 / 逐帧晶胞全有 | 纯文本，**零新依赖**，可照搬 `cp2k_out.rs` 的 token 锚点 + 区间扫描 |
| `vasprun.xml` | 同上，且结构化 | 要拉 XML 解析依赖 |
| `XDATCAR` | **只有坐标**，无力无能量 | 做不了训练集，只能当轨迹 |

倾向 `OUTCAR` —— 与 CP2K 那条路同构，`cp2k_out.rs` 的三条经验（按 token 序列匹配
而非固定行偏移、区间上界取下一个锚点、数值行不写死下标）可直接复用。

**待核对的坑**（动手前逐条验，别照记忆写）：

- **取哪个能量**：`free  energy   TOTEN` 与 `energy(sigma->0)` 是两个数。dpdata 取
  sigma→0 那个；要与 dpdata 对拍就得同口径。这条决定要写进 reader 的 doc 注释
- **应力符号与分量顺序**：`in kB` 行是 Voigt 六分量，顺序是 **XX YY ZZ XY YZ ZX**
  （不是常见的 YZ XZ XY），单位 kB = kBar，`units.rs` 的 `PressureUnit::Kbar` 已有。
  符号是否与 Ferro 的「正 = 压缩」一致**必须实测核对**，办法是拿一个已知受压体系
  或与 dpdata 的结果对拍 —— 对称张量下符号错了不会有任何形状异常
- **逐帧晶胞**：NPT 下每步都有 `VOLUME and BASIS-vectors` 块，取 `direct lattice
  vectors` 三行，**行优先**（`matrix3_from_row_major`）
- **元素与计数**：元素名在 `POTCAR:` / `VRHFIN =` 行（VASP 5+ 有时只在开头出现一次），
  每种个数在 `ions per type =`，顺序 = POSCAR 顺序。三处都缺就报错，不猜
- **力块**：`TOTAL-FORCE (eV/Angst)` 锚点，`-----` 夹住；坐标与力同表且坐标是笛卡尔
- **截断与未收敛**：最后一帧常被中断切掉；与 `cp2k_out` 同策略——缺块的帧宁可丢，
  丢帧计数由 stats 带出
- **ML_FF 的 OUTCAR 排版有别**（多机器学习力场块、能量行不同），先只保证纯 AIMD，
  遇到再说

**CLI 接入**：`collect` 需要格式分派，而**不能靠扩展名** —— `OUTCAR` 没有扩展名，
`.out` 又太通用（这也是它至今没进 `io_dispatch` 的原因）。判据：collect 是按目录
批量的，一批里可能 CP2K 与 VASP 混，故走**文件名 + 内容嗅探**（读头若干行找
`vasp.` 版本横幅 / `CP2K|` 横幅）比加一个全局 `--format` 开关好。
`io_dispatch` 侧可按前缀 `OUTCAR` 注册只读格式（与 `POSCAR`/`CONTCAR` 的前缀判断
同一模式），让 `ferro convert -i OUTCAR -o traj.xyz` 也能用；注册时记得
**`ferro-python/src/io.rs` 是另一处独立的分派**，加格式要两边都看。

**Voigt 顺序表（`ferro-io/src/voigt.rs`）跟着这条待办建**：extxyz 那轮本打算先建，
但 6 分量在 extxyz 侧已决定拒收，表在那里没有调用者，先建等于先造一个没人用的抽象。
两条已实证的顺序可直接写进去：

| 来源 | 6 分量顺序 | 实证 |
|---|---|---|
| extxyz 规格 / ASE | `xx yy zz` `yz xz xy` | `ase/stress.py:84` |
| VASP `in kB` / GPUMD `stress_*.out` | `xx yy zz` `xy yz zx` | `ase/io/vasp_parsers/vasp_outcar_parsers.py:93` 的 `[[0,1,2,4,5,3]]` 重排 |

LAMMPS 那条**没有实证，先不写** —— 半可信的表比没有表危险。

---

### scripts/：net 剩余四张表的画法

四个发表级绘图脚本已完成。**net 六张表里还有四张没有画法**，因为用户明确说还没想好
形式，不替他猜：

| 表 | 状态 |
|---|---|
| `network_qn` / `network_qn_partner` | 已画（`plot_net.py qn` / `--partner`） |
| `network_coordination` | 已画（`plot_net.py cn`） |
| `network_composition` | **未定形式**。它是其余表的摘要，一物种一行；堆积柱会跟 qn/cn 两张图重复 |
| `network_ligand_type` | **未定形式**。伙伴对细分后类别数不定 |
| `network_linkage` | **未定形式**。热图方案已做出原型并被否，见下节 |
| 实验数据对比 | **未定形式**。多半是单点或水平参考线，与 MD 的条带不是一种几何 |

`plot_net.py` 的 `kind` 已经是子命令位置参数，加一种图就是加一个分支 + 一个
`series_*` 函数，不必重构。

#### 口径改动已跟进（2026-08-20，提交 `a87843d`）

原以为脚本会「断在缺列上」，实测**没崩**——`series_qn` 早有 `not m_cols` 的兜底分支。
真正的问题是两处**逻辑错误**，比崩溃难查：退化警告的判据（`len(m_cols)==1`）在新
口径下把最有信息的那一维误判为退化；无 `m_` 列时报「给的是 network_qn.csv」，而
Zn–P–O 这类无异核形成子的体系其 `qn_partner` 与 `qn` 列结构本就相同、无法区分。

顺带补了 `--partner` 的图例注解：图例只画 Qn 色相，明度档携带 `m_<X>` 却从无说明。
旧口径下明度维是退化的所以无所谓，现在它是真正的信息维。

**`--partner` 是否改回默认仍未决**：当初设成开关的唯一理由（单形成子体系下
`m_P ≡ qn`，展开是假的两级结构）已经消失。这属出图形式问题，与四张表的画法一起定。

#### linkage 热图原型（2026-08-15，五种版式全部被否）

做了五张原型（三角 / 镜像 / 上下三角双量 × facet 网格 / block 摘要），跑在
`tests/` 两条轨迹上。**用户看过实物后判断效果均不理想，形式重新待定。**
产物留在 <https://claude.ai/code/artifact/667e09ec-884d-4a82-8905-fbb01eb65cfb>；
脚本在会话 scratchpad，未进 `scripts/`。

下次重开时**不必重跑**的事实：

- **「热图会是三角形」这句原判断是错的**。规范半边按 `(元素, 桥接数, 配位数)` 排序，
  `elem_a ≠ elem_b` 时次序完全由元素决定，两端各自遍历全部状态 —— **异元素 block 是
  满矩阵**，只有同元素 block 才是三角。而 Al–P 这个满矩阵恰是最有物理意义的一张
  （NMR 的 Al[4]/[5]/[6] × ³¹P 的 Qn 直接对照）
- **`n_formers == 2` 是干净切口**：过滤后 Al–O–P 计数求和 = 2693，与 `ligand_type`
  的 `Al-O_b-P` 行**逐字相等**。可直接写成自检断言。不过滤则分母含义从「桥」滑到
  「配对观测」，三簇氧抽出的一对与真桥混在一起
- **跨成分可比只能靠零模型，不能靠归一化**。行归一化只除掉一端的丰度，Al 含量变化时
  数字照样变，说不清是亲和性变了还是 P 变少了。可行的是桥端随机配对：
  `p_i = e_i/2N`（`e_i` 为该标签作桥端出现次数，对角贡献两个端），
  `E_ij = 2N·p_i·p_j`（i≠j）/ `N·p_i²`（i=j），值取 `log2(O/E)`。
  `Σ_{i≤j} E_ij = N` 自洽。**`p_i` 逐成分各算各的**，丰度被分母同步吸收；边缘取自
  linkage 表自身而非 `composition` 表，否则 V 里残留丰度信息
- **实测结论（口径可信，只是画法不好）**：`43Z43P15A` 的 Al–O–Al 压制 39 倍
  （10/388），Al_4–Al_5 / Al_4–Al_6 / Al_5–Al_5 三格观测为 0 而期望 103/14/7；
  Al–P 全线 +0.4~+1.4；P–P 高 Qn 之间 −0.5~−2.4。`70Z30P00A` 的 P–P 三格 log 比全部
  贴近 0，即二元磷酸盐的链是**随机连接**的 —— 这个零结果是三元体系择优的对照
- **block 级摘要在单形成子体系下恒等于 0**：只有 P–P 一个 block 时全部桥都在里面，
  obs ≡ exp。它只在多形成子 + 三个以上成分时才有定义
- **对角线不能参与镜像**：i–i 桥只存一份，但它**贡献两个同类桥端**。边缘统计要算两次、
  矩阵取值只能算一次，两者极易混 —— 原型第一版把 `P_3–P_3` 的 162 算成了 324

### 三项小待办（2026-08-13 提出）

三件事互相独立，都不大，凑在一起是因为都属「上一次改动划在范围外的部分」：

1. **`bader` 的 `--outdir`**。它不走 `CommonArgs`，三个 `.dat` 一律写在当前目录。
   要连 `write_acf` / `write_bcf` / `write_avf` 的签名一起动，与 Henkelman 格式无关，
   只是路径。

   **原先这条写的「三个固定名字」不准确**（2026-08-22 核对代码）：名字是
   `<输入stem>_ACF.dat`，跟着输入走。但危害不变甚至更隐蔽 —— VASP 的电荷密度
   一律叫 `CHGCAR`，所以 `run1/CHGCAR` 与 `run2/CHGCAR` 在同一个工作目录跑，
   两次都写 `CHGCAR_ACF.dat`，后一次静默盖掉前一次。已在 `ferro bader` 的帮助页
   与手册里如实告知，`--outdir` 仍待做
2. **`convert` / `job` 的 `-o` 统一为文件名**，路径走 `--outdir`。现在它们的 `-o` 是完整
   路径，与其余 11 个命令的约定相反。（2026-08-22：这条差异已先在 `ferro convert`
   的帮助页与手册里写明，行为未动）

   **2026-08-25 更新**：现在共有**三种** `-o` 语义 —— 分析命令的文件名后缀、
   `convert`/`job` 的完整路径、`dataset` 的输出根**目录**。第三种不是历史包袱
   而是必然的：DeePMD 的 system 就是目录，没有「文件名」可言。所以统一的目标
   应当是**两种**（后缀 / 目录），把 `convert`/`job` 的完整路径归入前者，而不是
   指望三种收敛成一种。
3. **`ferro-io` 的 writer 路径统一 `&str` → `&Path`**。九个 writer 全收 `&str`，而
   `batch.rs` 内部已是 `PathBuf`，只能在边界 `to_string_lossy()` 转一次
   （`Output::join_str` 就是为此存在）。改动机械但会碰到 io_dispatch、ferro-python 与
   一批测试，故没有混进命名那次提交

### ferro job 从轨迹抽帧（2026-08-22 提出，告警已于 2026-08-24 补上）

`job` 现在无条件取第 0 帧（`cmd/job.rs:186`）。多帧输入会告警（见下），但**仍然只
生成第 0 帧的输入**——喂一条 500 帧轨迹得到的还是最没平衡那个构型，只是这次用户
知道了。真正的抽帧 + 批量生成仍待做。

`convert` 的 `--start/--end/--stride/--number` 已把选帧逻辑做进
`Trajectory::select_indices` / `spread_indices`，job 复用即可，不必重写。
真正要定的是**产物命名与批量语义**，与 `convert` 不完全一样：

- job 的 `-o` 现在是完整路径且有默认值（`job.gjf` / `job.inp`），多帧要变成
  `job_0000.inp` 这类。与待办 #2「`-o` 统一为文件名」是同一件事，该合并做
- 一个构型一个输入文件是显然的（QC 输入本来就一个结构一个），所以不存在
  `convert` 那种「一个文件还是 N 个」的分支 —— 恒为 N 个
- ~~**最小可用的第一步是加警告**~~ **已做（2026-08-24）**：`multi_frame_warning`
  打三行 —— 帧数、忽略了几帧、以及可直接抄的 `ferro convert --number` 命令。
  产物形态未动。**做这一步时撞出 `ferro job` 在 debug 构建下必 panic**（clap 的
  `help` 参数重名），见 `issues.md`「构建 / 工具链陷阱」

在此之前的变通：`ferro convert -i traj.dump -o conf.vasp --number 20` 抽成单帧
文件，再逐个喂给 job。

### ferro dataset：剩余小项（2026-08-25）

`collect` / `filter` / `merge` 均已落地。剩下的都不大：

- **多层周期镜像**：几何判据现在把最小镜像当上界（超出报错），没做参考脚本的
  多层扫描（`n = floor(rcut / w + 0.5)`）。当前体系盒子远大于阈值，不构成限制
- **额外键搬运**：`atom_ener` / `fparam` 这类 `Frame` 装不下的项，读时告警、
  写时丢失。真出现时再设计（需要一条绕过 `Trajectory` 的按帧索引搬运通道）
- **GPUMD/NEP 的 `train.xyz` 导出**：已单独提为本章第一条「GPUMD/NEP 的
  `train.xyz` 导出」；它的前置（extxyz 的 stress/virial 符号与体积因子）已完成

两件**不必新写**的事已经在库里：帧区间与间隔用 `Trajectory::select_indices` /
`spread_indices`（`convert` 的 `--start/--end/--stride/--number` 就是它）；
O-O 间距与 Al6 配位用 `ferro_core::classify_frame` 出的
`AtomType::Former{cn,..}`（`Al_4`/`Al_6` 的数字就是配位数）。`filter` 在 CLI 层
组合 io 与 analysis，中间层仍不互依赖。

还需要定的：`filter` 要读回 npy（`ndarray-npy` 读侧未用过）；merge 的 `type_map`
跨 system 对齐（`collect` 已按 `(Z, 符号)` 排序，这条规则可复现是对齐的前提）。

### ferro map chg-sdf 的 --cubes 拆成单文件（2026-08-11 提为高）

现在是多个 cube 聚合成**一张** SDF（`cmd/map.rs::run_chg_sdf`），与 `ferro map` 其余
模式「一个输入 → 一个 `.cube`」的语义相反，`--cubes` 也是 `-i` 之外唯一的输入参数。

阻塞点不在遍历而在**中间产物格式**：跨文件加权平均不可交换，逐文件出产物之后若还想
聚合，产物里必须带样本计数（否则 5 帧的 SDF 与 500 帧的 SDF 被等权平均）。
故顺序是：先定义带计数的中间格式 → 再拆 `--cubes` → 再接批处理遍历。

---

## 优先级中

### ferro-cli：REPL / 脚本模式（2026-08-11 由高降中）

`main.rs` 现在是子命令分发器，裸 `ferro` 打印分类总览。REPL 落地时改为
**tty 进 REPL、管道读 stdin**（`python`/`node`/`irb` 的惯例；`isatty` 判断，CI 里
`echo ... | ferro` 不会挂住）。

- 依赖 `rustyline`；三种模式：交互 REPL、脚本文件（`ferro -f workflow.mf`）、管道输入
- **必须在同一进程内链接全部命令**：`read` 之后轨迹要留在内存里给后续命令用，
  这正是 REPL 相对 shell 循环的全部价值。这条也是 0.2.0 决定合并成单二进制的理由 ——
  子进程分发做不到状态保持，而链接进来之后再包一层 `fe-*` 前端就是重复
- 脚本语法应建在**已经定型**的子命令树上（`gr -a P -b O` 直接复用 `cmd::traj` 的
  参数结构体），不要另起一套
- 注意与「批处理输入」是两件事：这里是**命令**的批处理（一个脚本跑多条命令），
  那里是**输入文件**的批处理（一条命令跑多个轨迹）。两者可叠加但互不依赖

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

### ferro-python：pyo3 0.29 运行时验证

0.21 → 0.29 只做了类型层验证。本机无 maturin，`cargo build` 在 macOS link 阶段过不了
（`extension-module` 需 `-undefined dynamic_lookup`）：

```bash
pip install maturin
cd ferro-python && maturin build --interpreter "$(which python)"
pip install target/wheels/*.whl
```

冒烟测试要覆盖：读 xyz/cif/lammpstrj、`supercell`、`write`、
**新拆的 `gr_pair` / `gr_all`（含 `by="label"`）**、`msd`。

### 是否采用 rustfmt（待定，暂缓）

收益：把格式从代码审查面里移除。代价三点：

1. 一次性约 80 文件的大 diff，`git blame` 在这些行上全部指向那一次提交
2. 会拆掉现有的刻意列对齐（`gr.rs` 的 `WelfordStats`、`angle.rs` 的 `CellList`、
   `units.rs` 的枚举表）
3. 生成文件要显式排除：`rustfmt.toml` 的 `ignore` 仅 nightly 可用，stable 上需给
   `cp2k_basis_db.rs` 的静态表加 `#[rustfmt::skip]`

若采用：单独一次纯格式提交（`style: adopt rustfmt`）+ 固定 `rustfmt.toml`，
不要混进特性或依赖升级提交。

### 搁置项

- **`vanhove` 加 `tau` 列**：现在一次只算一个 τ、写在 `#` 头里。加列后将来支持多 τ 是
  加行而非改列结构

### ferro-structure：补充结构操作

`rotate.rs` / `orient.rs` / `substitute.rs` / `disturb.rs` / `select.rs`。

移植参考 Multiwfn（`examples/Multiwfn_2026.4.10_src_Linux`）的 `otherfunc3.f90`
`geom_operation`（行 1975–3066）：

| 待办 | Multiwfn 参考 |
|---|---|
| `rotate.rs` | 菜单 3、4（绕笛卡尔轴/键/指定向量旋转、旋转矩阵） |
| `orient.rs` | 菜单 5/6/8/11（对齐键/向量/最长轴/平面到笛卡尔轴或平面） |
| `disturb.rs` | 菜单 18 / `displace_geom`(1879)（高斯随机位移，默认 σ=0.03 Å） |
| `select.rs` | `util.f90:985 str2arr`（`"2,3,7-10"` 选择语法） |
| `substitute.rs` | 无直接对应（Multiwfn 仅 15/16 加删原子） |

其余可参考项：晶胞数学在 `PBC.f90`；菜单 20 边界分子补全、22 原子折叠入胞、
25 提取团簇、28 坐标轴互换。

### ferro-workflow：VASP（POSCAR/INCAR/KPOINTS）

---

## 优先级低

### `ferro doc` 的 markdown 渲染（2026-08-26 提出，已量过代价）

`ferro doc` 现在原样打 markdown，零依赖。三个不好看的地方：表格不对齐、
长行不按终端宽换行、LaTeX 在终端里读不了。

**依赖代价已实测**（各方案的 `Cargo.lock` 与 Ferro 现有 152 个包求差）：

| 方案 | 净新增 crate | 内容 |
|---|---|---|
| 自己写 | **0 或 1** | 只有 `unicode-width`（零依赖小包） |
| `minimad`（termimad 的解析层） | **1** | 只有它自己 |
| `termimad` | **38** | `crossterm` `mio` `signal-hook` `parking_lot` `rustix` `regex` `syn`… |

termimad 那 38 个是一整套 **TUI 事件循环栈** —— 为了把静态文字打出来拉进
信号处理、异步 IO 轮询、锁原语，不成比例。`minimad` 只 +1 但**只做解析不做
布局**，而要的恰恰是布局，省不下多少。

**工作量按 3918 行手册的实际构成估**（普通段落 47% / 代码块 22% / 表格 15% /
标题 7% / LaTeX 5.5% / 列表引用 4%）：**约 250–350 行**，全在 `doc.rs` 里，
不碰别的 crate。其中表格一件事占 ~150 行。

三条已经查证、下次不必重跑的事实：

- **`unicode-width` 几乎是必然的**：表格行里 **46%（266/573）含中文**，
  `cli-reference.md` 整页是中文。按字符数算列宽会让中文表格全部错位。
  为省一个零依赖小包去手写 Unicode 宽度区间表不划算
- **终端宽度不新增编译单元**：`libc` 已在树里（经 `plotters → font-kit →
  core-foundation`），`ioctl(TIOCGWINSZ)` 可用，但要在 `Cargo.toml` 里
  声明成直接依赖
- **LaTeX 是天花板，谁都解决不了**：216 行含 `$`，形如
  `$$P = \frac{2}{3}\frac{E_kin}{V}$$`。termimad 也不认 LaTeX。自己做
  符号替换表是另一个项目，且对 5.5% 的行做**部分**正确的转换容易比原样更糟。
  无论走哪条路，`collect` 那几页的公式都还是现在这样

倾向：自己写 + `unicode-width`。收益的大头集中在表格，而表格恰是自己写最容易
做对的部分 —— 只有 GFM 管道表，无对齐标记、无嵌套。做完看实物，不行就
`git revert` 一个提交。

### 机器学习集成

- 阶段一：`linfa`（K-Means/DBSCAN/PCA）→ `ferro-analysis/src/ml/`
- 阶段二：`candle`（ONNX 推理，加载 DeepMD-kit 势函数）
- 阶段三：`burn`（纯 Rust 训练，远期）

### 深度学习工作流 I/O

DeePMD-kit 的 npy 系统目录**写**侧已完成（`writers/deepmd.rs`，2026-08-25），
读侧与 GPUMD/NEP 的 `train.xyz` 导出待做（后者是链末的 export，不是中间格式）。
MACE/NequIP 兼容格式仍未开始。

---

## 已完成（归档）

### extxyz 的 stress/virial 修正（2026-08-26，0.3.2）

四处静默缺陷，共同根因是**在没有依据的地方替用户猜了一个约定，且猜错不报错**：
`stress` 取不到就回落取 `virial`（差一个体积因子）· 两侧都没做 ASE（正 = 拉伸）与
`Frame::stress`（正 = 压缩）之间的变号 · 九个数按行优先处理而 ASE 文档说列优先 ·
`Properties` 的力列只认复数 `forces`，读 GPUMD 的 `train.xyz` 会丢掉全部受力。

**被查证推翻的原计划**（三条，都发生在动手之前）：

1. 原计划「读到 `virial=` 就报错，因为符号约定无法核实」。实测 dpdata 1.0.2 的
   `virials = -volume * stress_ase` 双向换算 + GPUMD 手册 + extxyz 规格的
   `virial -> stress` 乘 `-1/cell_vol`，三方一致，符号完全可定 —— 报错等于明知
   怎么读却拒绝读
2. 原计划「按列优先读写以对齐 ASE」。查到规格要求该张量**对称**（"fail if not
   symmetric"），而对称下行/列优先逐位相同 —— ASE 说 Fortran order、GPUMD 手册
   拼成 `vxx vxy vxz vyx ...`，两家描述相反却从无人报 bug，正是这个原因。改为
   **检查对称性**，不站队
3. 原计划「6 分量按标准 Voigt 收下并告警」。用户判断：让用户重排与重新生成 9 分量
   的工作量几乎没差别，那就取最稳的结果、正确性交回给用户 —— 改为**拒收**

`Lattice` 全程未动：实测 ASE 写出的九个数就是三个晶格矢量依次排开，Ferro 原本正确。

### 帮助页精简：其余各页（2026-08-26，0.3.2）

`dataset` 三页的模板推到了全部命令。**实际范围与「22 页」的估计不同**：

- **只有 8 页超 40 行**，其余 16 页本来就短，只需补一行手册指针
- **`ferro net` 的帮助页不在 `help.rs`**（在 `cmd/net.rs` 的 `HELP_EXTRA`，
  紧挨着它需要的 argv 剥离逻辑），所以此前**防漂测试完全没覆盖它**，而它
  有 10 个选项。已纳入
- **`convert` / `info` / `bader` 没有手册专页** —— 原以为要新写三页，实测
  `cli-reference.md`（863 行）已按命令分节覆盖了全部参数，缺的只有
  `--help`/`--input` 这类普适项

最后一条导致了唯一的设计改动：**`ferro doc` 支持按小节寻址**。`Page` 加
`section: Option<&str>`，从该 `##` 标题取到下一个同级标题。`ferro doc convert`
于是出 99 行而不是整本 863 行，手册也不必拆成一堆按命令的小文件、让「唯一
一份完整参考」碎掉。小节标题改名而表没跟上时回落到整页，另有测试钉住每个
`section` 真实存在。

**收尾时仍有 4 页超 40 行，且都不该再砍**：

| 页 | 行 | 大头是什么 |
|---|---|---|
| 顶层总览 | 61 | 它是入口地图，不是命令页 |
| `net` | 54 | 10 个选项 + 6 张输出表 |
| `convert` | 59 | 其中 **27 行是 `supported_formats()` 生成的**格式表，手写只 32 |
| `job -s cp2k` | 47 | 23 个参数，值域枚举本身就是参数表 |
| `dataset filter` | 43 | 18 行参数表 |

**判据：参数表（含值域枚举）是唯一必须完整的一段，不为压进 40 行去砍它。**
40 是目标不是硬上限。

### collect 语义重做 + filter 报告落盘 + ferro doc + 帮助精简（2026-08-26，0.3.2）

四个提交，顺序即依赖：collect 语义 → filter 报告 → `ferro doc` → 帮助精简 +
手册对账 + 防漂测试 + 本轮记录。第四步必须最后，它照着前三步的**实际**行为写。

被实测推翻或修正的：

- **「用 `_` 串联多级目录」是用户的初始提法，商议中被他自己换掉**：改为
  **保留目录结构**（`sets/a/md` 而不是 `sets/a_md`）。理由是 `filter` 的
  `-o` 本来就按相对路径重建，`find_systems` 与 `merge` 也都是递归的，
  压平反而是这条链上唯一的例外
- **「撞名报错」整条判断消失**：新规则下撞名就是「该合并」的定义，
  原来那个 `bail!("would both write the system directory ...")` 删掉。
  命名规则与分组规则合成同一条，所以这两件事不能拆成两个提交 ——
  中间态是「撞名报错但名字已变」，没法对拍
- **去重被用户判为不必要**：重启重叠段位置速度相同，能量力也相同，
  不稀释结果；重启间隔通常 5 步以内。改为不去重但**打出 step 区间**，
  让这个前提保持可检验
- **「时间序没法区分」这句原判断是错的**：`MD| Step number` 正是解析器
  挂一切的锚点，xyz 注释行还带 `time =`。两个值当时都读到了但都没存
- **「诊断只在只读模式算」的省时理由不成立**：实测 1110 帧 / 302 原子
  挂钟 0.23 s（带）vs 0.28 s（不带），rayon 跑满六核。改为恒算，
  `--no-diagnostics` 开关不必加
- **报告平铺还是塞子目录，用户的直觉对且理由更硬**：`expand_dirs` 只收
  `is_dir()`，平铺的 csv 会被后续 `merge -i clean/*` 自动滤掉，而
  `report/` 子目录反倒会被收进去当 system 候选
- **`clap_mangen` 解不了这个问题**：它从 clap 定义生成 man page，只有
  参数表 —— 正是要精简掉的那部分，散文一句都带不出来。`cargo doc` 是
  rustdoc（API 文档），也不是一回事。故 `ferro doc` 自己做
- **防漂测试当场抓到 19 处真漏**，远超预期的 1 处（`--shuffle`）。
  三次收紧判据才落到可用：短名也算写了；允许交叉引用别的命令的参数；
  跳过「there is no --x」这类否定陈述。剩下 2 处是有意不写，进白名单


### 早期（2026-05 ~ 06）

| 日期 | 内容 | 落点 |
|---|---|---|
| 05-10 | Bader 电荷分析 | `charge_grid.rs`、`chgcar.rs`、`dft/bader*.rs`；规格见 `bader.md` |
| 05-14 | network 重构 + 类型分类迁移 | `network_type.rs`、`typing.rs`、`cluster.rs` |
| 05-15 | 电荷密度团簇 SDF | `dft/chg_sdf.rs`（Kabsch + pull 插值旋转） |
| 05-15 | CP2K 输入生成 | `workflow/cp2k.rs`，混合泛函自动生成 `&HF` 块 |
| 05-16 | 未成对电子 + 基组库 + QE | `spin.rs`、`cp2k_basis_db.rs`（2829 条）、`qe.rs` |
| 05-16 | cube_sdf 迁移 + 全项目审计 | 共享原语下沉 `core/cluster.rs`，**ferro-analysis 去掉 petgraph**；审计 9 项修 8（#8 为误报） |
| 05-16 | ferro-python PyO3 绑定 | 独立 workspace，旧占位代码全部重写 |
| 05-17 | MSD 绘图 + 自扩散拟合 | `fit_diffusion`（D=slope/6，Einstein 3D）+ R²；拟合区间是**滞后时间轴的分数** |
| 06-26 | 代码审查修复三项 | g(r) r_max clamp、`bader_ongrid` 删重复循环、`box_builder` 近邻去重 |

### g(r) / CN / S(q) 重构 + element/label 拆分（2026-08-08，0.1.10–0.1.11）

提交 `caa309b`（0.1.10）· `c5a48b6`…`fc0a3b0`（0.1.11）。测试 312 → 334。

需求来源：`-a P -b O` 与 `-a O -b P` 输出完全相同，`-a`/`-b` 顺序对 CN 不起作用。

定下的语义（后续一直沿用）：

- **`gr` 对称、`cn` 有向**：`CN(A→B) = hist/(N_A·steps)`，`-a` = 中心、`-b` = 近邻
- **未加权 total S(q) 删除**：`f_i ≡ 1` 的退化情形对应不了任何实验探针，
  参考实现 `examples/code2` 中根本不存在
- **配对模型改 n² 个有序对**（3 元素 → 9），镜像对的 `gr` 数值重复照写，
  换取「每配对一组可直接提取的列」
- **拆分规则按第一个下划线**，不用贪婪前缀匹配 —— 否则 `Pb`→铅、`Po`→钋 会盖过
  「P bridging」这类伪标签意图。只作用于 LAMMPS dump；cp2k/qe 的 `extract_element`
  与 cif 的 `element_from_label` **不动**（它们处理 `Fe1`/`O2` 原生命名，字母前缀规则
  对其正确，套用下划线规则反而让 `Fe1` 退化成元素 `Fe1`）

明确不做（当时）：`neighbors_of` 每对枚举两次的 2× 浪费（纯性能）；三处
`extract_element` 的统一；旧标签格式（`P0`/`Ob`/`On_P`）的读取兼容层 —— 重跑一次即为
新格式，加格式猜测反而可能把真元素误判成旧标签。

### 依赖包全量升级（2026-08-08，0.1.12）

`ebd6147` → `915412b`。陷阱见 `issues.md`「依赖升级陷阱」。
`nalgebra` 0.34→0.35 后 g(r)/CN 输出与升级前**逐字节一致**。

### NPT 逐帧体积归一化（2026-08-08，0.1.13）

口径由「先时间平均、后归一化/变换」改为对齐 code1/code2 的「先逐帧归一化/变换、
后时间平均」。完整推导、实测差异与「明确不做」见 `issues.md`。

要点：关键恒等式 `ρ_f·g_f` 中 V_f 自行抵消 → 逐帧变换可由 `⟨ρ_f·g_f⟩` 与 `N·⟨1/V⟩`
两个时间平均量**精确**重构，无须每帧各做一次 FT，分层与 `calc_sq_from_gr` 签名不动。

fixture 由生产轨迹**等间隔**取 5 帧（连续取会丢掉体积跨度，测不出东西）。

### box_builder 括号公式（2026-08-08，0.1.14）

表面是「`parse_formula` 不支持 `Ca3(PO4)2`」，**实际阻塞点在上游**：`build_box` 只把
数据库里的 `cd.formula` 喂给解析器，而库中 26 条全是无括号有机溶剂；用户输入的
compound 一旦不在库中，`compounds::find` 返回 `None` 就报错了，**根本到不了解析器**。
只加括号支持等于写死代码。故一并做了三件事：栈式解析、`resolve_component`
（库外化合物当化学式解析）、收紧输入校验。

### 批处理 + Table 长表 + 单二进制重构（2026-08-09，0.2.0）

锚点 tag `v0.1.15`。测试 362 → 396。

值得记住的判据：

- **`-i` 恒为 `Vec`，单一代码路径**，N=1 是 N 的特例。按文件数分派会让产物形态取决于
  glob 当天匹配到几个文件，还意味着两套命名、两套绘图、两套错误处理
- **`Table` 下沉 core 而非把 writer 搬进 io**：完整论证见 `issues.md`
  「分析产物为什么不进 ferro-io」
- **表结构跟主产物的粒度走** —— 这是 gr/angle 长表而 sq 宽表的判据。决定性理由是
  不同轨迹的**元素集可能不同**（Zn-P-O 与 Al-P-O 同批跑），配对方向做宽表就要取列
  并集补空洞
- **`to_tables()` 定为固有方法不是 trait**，留待第二个消费者出现
- **不拆「纯搬家」提交**，`Table` 迁移 + 批处理 + 长表一次走完
- **`--plot` 冻结为自检用途**，不追 matplotlib。任何「加对数轴/误差棒」的需求一律
  指向 Python（长表 + `sns.lineplot(hue="file")` 一行）

实施中比原计划多出的：`angle` 多一列 `count`（整数直方图是与 dump2analysis 逐 bin
对拍的依据，只留归一化的 `p` 会废掉这条验证路径）。

### ferro net 重构：结构化类型 + 标签重做 + 合并命令 + linkage（2026-08-12，0.2.1）

原计划的「标签体系重做」与「接入批处理 + 长表化」两项，实施时发现它们不是两件事的
先后，而是同一件事的两层：**标签之所以改不动，是因为分类结果以字符串形式流转。**

| # | 提交 | 内容 |
|---|---|---|
| 1 | `5f8ac0e` | `classify_frame` 返回结构化 `AtomType`，消除**五处**标签反解析 |
| 2 | `4503d11` | 标签改 `<元素>_<后缀>`；删修饰子角色分类 |
| 3 | `f1f929d` | `net` 降为叶子命令，接 `CommonArgs` + 批处理 + 长表 |
| 4a | `801184b` | 分布表加 `sd` 列（Welford）；氧的伙伴改数据列 |
| 4b | `304bdd9` | linkage 长表 + `Q^n(mAl)` 分解 |
| 5 | `2f9021e` | 标签存 `atom.label`；extxyz 加 `label:S:1`；dump 折叠 + type 编号跨帧固定 |

**第 1 步单独成一提交是关键**：它把「信息不再经字符串往返」与「标签换格式」分开，
前者可逐字节对拍验证（22 个产物文件一致），后者才是行为变更。合到一起就没有能对拍
的中间态。

被实测推翻或修正的原计划：

- **`X` 兜底的分裂比预期严重**：计划只说「≥3 配位氧与 ≥3 NBO 修饰子塌缩」，实测单帧
  `X = 181` 里 **180 个是 Zn、1 个是氧**
- **修饰子角色分类整体删除**（原计划是保留 `Zn_f/_t/_b`）：实测 97.1% 落进兜底桶
  （0 / 1 / 26 / 903），该分档对 Zn 无分辨力
- **`qn` → `n_bridge` 这一步后来被用户纠正**：「Al 没有 Qn」的意思是 Qn 计算直接忽略
  Al，而不是把列改名。改列名保留了 Al 的行，等于用一个更含糊的名字继续报一个对 Al
  无意义的量。正确做法是让 Al 退出 Qn 表、列名退回 `qn`
- **动机比计划写的窄**：`net type 导出 → traj gr -x P_3 -y O_b 直接串起来`
  只在**单帧**成立，多帧被粒子数守恒守卫拒绝

### ferro net 产物可读性重做（2026-08-12，未发版）

起因是三条反馈：帮助文档太长；四个文件名看不懂；表里把 label 拆成 `former` +
`n_bridge` 反而更难读。追问下去牵出上面那条更根本的「Al 没有 Qn」。

提交 `450252d` · `100265d` · `e2821c6` · `c26cde6`（帮助 120→58 行）· `9536c1b`。
第二轮反馈续做：`a0d21d6`（单元/原子两套词汇、`average`→`composition`、
`ligand_type` label 合并、配体分母改逐元素）。

设计判据（grilling 里逐条过过，坑的部分见 `issues.md`）：

- **`label` 新增而非替换数值列**。label 给人读，`former`/`qn`/`cn` 给筛选和画图。
  删掉数值列会让「筛出 Qn ≥ 3」退化成字符串切分 —— 正是第 1 轮从五处调用点消灭掉的
  反向解析，不该在输出侧重新造一个
- **展示列的数字含义按元素而变是可以接受的**（`Al_4` 是配位数、`P_2` 是 Qn），因为
  那是文献自己的读法。曾以「标签要被 `traj gr -x` 机器解析」为由否掉，后来发现站不住：
  用户要的正是 `-x Al_5` 选出五配位 Al
- **详细说明进文件头**。帮助会滚走、手册在别处，`#` 头跟着文件走
- **用户否掉了自己最初提的 `Q3(2Al)` 写法**：实测参考数据里 `(qn=2,m_Al=1,m_P=0)` 与
  `(qn=2,m_Al=1,m_P=1)` 会渲染成同一个 `Q2(1Al)`，前者的 Σm 亏空 1（那座桥是三簇氧）
- **`linkage` 保持原子词汇**，理由由用户给出且比一致性论证更强：桥联表达的是
  **原子之间**的连接，Qn 是包含多个原子的**单元**

对拍（`43Z43P15A_NPT_5`）：P 的 Qn 计数 25/276/653/751/155 不变、linkage 总观测 3589
不变、`ligand_type` 与 `coordination` 逐行相同、两张 Qn 表 count 闭合。导出轨迹读回后
`traj gr -x Al_5 -y O_b` 得 CN = 4.933，缺的 0.067 是一个三簇氧（标 `O_t` 不标 `O_b`）
—— 另一条代码路径的交叉验证。测试 398 → 420。

### traj 产物命名带 label + --outdir + sq 移除选择（2026-08-13，未发版）

起因是一句具体的抱怨：`gr.csv` 看不出算的是哪一对。范围从「gr 和 angle」扩到六个命令，
并顺带拆掉了 sq 的配对选择。提交 `3062954` · `2646232` · `5234fba`。

判据：

- **label 排在 suffix 之前**。label 说的是算了什么，suffix 是批次标记；这个顺序让
  `ls gr_P-O_*` 列出同一对在各批次的结果，反过来没有对应的用法
- **两种拼法并存且各有理由**：`gr`/`angle` 按写的顺序，`--elements` 排序去重。
  看着不一致，但统一成任何一种都会错一半（详见 `issues.md`）
- **`_all` 是有代价的选择**。它让所有不带筛选的旧命令产物改名，换来「文件名一眼看出
  有没有筛选」。用户明确选了这一边
- **sq 移除选择的代价已知并接受**：按 label 分辨的 partial 从 CLI 消失（`-x/-y` 曾是
  进入 `GroupBy::Label` 的唯一入口）。理由是位点标签对应的原子数往往不足以让 partial
  显出信号。**库层 `GroupBy::Label` 不动**

实施中发现与设计讨论不符的两处：`rotcorr` 的 `--center`/`--neighbor` 其实是必填的
（讨论里假设的 `rotcorr_all.csv` 那条路根本走不到）；删掉 `SelectArgs` 的**使用**不等于
删掉参数（详见 `issues.md`）。

**`vacf` / `vanhove` / `rotcorr` 的改名未经实测验证** —— 参考 fixture 没有速度，
`vacf` 跑不到写文件那步；这三个性质本身也还没调试到。单元测试覆盖了拼名函数，
但端到端只验了 gr / angle / msd / sq / net / map。

### scripts/：四个对拍脚本修复（2026-08-11）

0.2.0 的输出变更让四个对拍脚本全部失效。典型断点：`compare_rdf.py` 的
`[float(x) for x in line.split()]` 遇到 `file` 文本列直接 `ValueError`。

**新增 `scripts/ferrocmp.py`** 收拢共用逻辑（四个脚本各有一份几乎相同的 `run` /
`load_columns` / `fe_version` / 列名行解析，形状一变就要改四遍）。陷阱见 `issues.md`
「外部脚本调用 ferro 的陷阱」。

复跑验证（数值与 `issues.md` 记录逐项吻合，说明改的只是读法不是口径）：

| 脚本 | 结果 |
|---|---|
| `compare_rdf.py` | P-O / Al-O 峰值与峰位完全相同，max\|Δ\| 5e-5 / 5e-4（dump2analysis `%12g` 的量化台阶） |
| `compare_angle.py` | Σfe/Σref = 0.5000（×2 计数约定）；`--align-binning` 下 1800 个 bin **整数零差** |
| `compare_sq.py` | 未补偿残差 q>15 均值 1.00031 / 1.00018（+1 常数）；partial 互相关 +0.959（type_new 失效） |
| `compare_sq_experiment.py` | 50 帧 rms(fe−exp) 0.0194/0.0277、max\|fe−ref\| 0.0027/0.0008、FSDP 1.95/2.05 |

顺带修：`compare_sq_experiment.py` 的 `TRAJ_CANDIDATES` 首选指向一个**不存在**的文件，
一直在静默回退到 5 帧子集。

### scripts/：对拍脚本跟进产物改名（2026-08-14）

产物加 label 段后三个脚本断了（`compare_rdf` / `compare_angle` / `compare_sq` 的 gr 段），
根因不是那三个字符串写错，而是**四个脚本各自手写产物名**。故做法是把命名规则搬进
`ferrocmp.py`：`file_label()` / `set_label()` / `product_name()` 是 `batch::out_path`
与 `batch::file_label` 的镜像，调用点只说「哪个模式、什么 label、什么后缀」。

- **`-o` 不再重复配对**。原先 `-o P-O` 是唯一区分手段，现在配对已在 label 段里，
  再塞进后缀就成了 `gr_P-O_P-O.csv`。四个脚本统一 `-o cmp` / `-o exp`
- 顺带修：`--ferro` 传相对路径必炸（`run_ferro` 以 outdir 为 cwd），含分隔符时先 `resolve()`

复跑四个脚本，数值与 2026-08-11 那轮**逐项吻合**（`max|Δ|` 5e-5/5e-4、Σfe/Σref
0.5000、q>15 均值 1.00031/1.00018、partial 互相关 +0.959、rms 0.0194/0.0277、
FSDP 1.95/2.05），说明改的只是拼名不是口径。

### scripts/：四个发表级绘图脚本（2026-08-13）

与 `compare_*.py` 的分工是清楚的：那边**对拍**（跟参考实现逐点比，图只为看差异），
这边**出图**（进论文的 pdf）。两者都读 ferro 的 csv，读完之后没有共同代码，所以是
两个共享层而不是一个。

| 脚本 | 子图 | 曲线 / 条带 |
|---|---|---|
| `plot_gr.py` | 一个 csv 一张 | 一条轨迹一条曲线；左轴 g(r) 实线、右轴 CN(r) 虚线 |
| `plot_angle.py` | 一个 csv 一张 | 同上，单轴 P(θ) |
| `plot_sq.py` | 一个 csv 一行两格 | 左 S^N(Q)、右 S^X(Q)，两条 total |
| `plot_net.py` | 一个形成子 / 元素一张 | x = 成分（`file` 列），100 % 堆积柱 |

判据：

- **net 的 x 轴是成分不是 Qn**。同一张图上看「各 Qn 占比随成分怎么变」，这是堆积柱
  相对折线的全部理由：它把「和为 1」画成图形约束。代价是**小分量看不出趋势**
  （2 % 的条带只剩一条线），要追小分量就把那列单独拉出来画折线
- **x 轴顺序 = `-i` 的参数顺序**，不从文件名解析成分数值。`43Z43P15A` 里有三个数字，
  脚本无从知道要哪个，猜错了图是错的但看不出来
- **`--partner` 是开关不是默认**。单形成子体系下 `m_P ≡ qn`，展开只会画出一条假的两级结构
- **CMD/MLMD 的方法维走文件名**，画成同一成分刻度下并排多根柱。这也是嵌套条带
  （而非双条带）方案的理由：并排那个位置要留给方法维
- **gr/angle/sq 不为方法对比设计**（用户判断是极小概率需求）；**不画误差棒**
  （`sd` 是快照间散布不是标准误）

编码陷阱见 `issues.md`「发表级绘图脚本编码陷阱」。
