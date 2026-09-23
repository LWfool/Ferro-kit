# 后续计划

> 归档只保留**判据**与**被实测推翻的原计划** —— 「做了什么、怎么改的」翻 git 历史
> （提交号已列出），「现在是什么样」看 `progress.md` 与 `docs/src/`。

## 优先级高

### `ferro doc` 的终端渲染器（2026-08-26 提出，2026-09-23 定案 —— **方案已定，可直接开工**）

手册已全英文（见归档），渲染器面对纯 ASCII，不必写 CJK 分支。

**十条判据**（2026-09-23 逐条议定，执行时不再回头）：

| # | 决策 | 理由 |
|---|---|---|
| 1 | 范围 = 只解决「自己读着难受」 | 表格是唯一确凿坏掉的：源文本全是 `\|---\|---\|` 紧凑写法，终端里必然错位 |
| 2 | **tty-only** | 渲染出的是 ANSI 控制码，进了文件就是乱码；`page_out` 的 `is_terminal()` 岔路口本来就在，不必加 `--raw` |
| 3 | **全自写**，净新增 0 crate | 见下 |
| 4 | 先翻译后渲染 | 已完成 |
| 5 | 公式：行内做、块级不做 | 块级 `\begin{pmatrix}` 无解；后补它是加一个 `Block`，无沉没成本 |
| 6 | 表格是收益主体 | 683 行表格，占手册 14.5% |

**实测依赖增量**（临时 crate 实装后与 Ferro 的 153 包求差，2026-09-22）：

| 方案 | 净新增 | 自写代码 | 冷编译 |
|---|---|---|---|
| termimad | **37**（`no-default-features` 也是 37，crossterm 是硬依赖） | ~30 行 | 7.5 s 墙钟 / 47 s CPU |
| minimad（只解析） | 1 | ~350 行 | ~1 s |
| **全自写** | **0** | ~420 行 | 0 |

「净新增 0」的两个依据：手册全英文后不需要 `unicode-width`；终端宽度不引
`terminal_size`（实测它拖进 rustix/errno/linux-raw-sys 共 **+5**），`libc` 已在树里
（经 plotters → font-kit），`ioctl(TIOCGWINSZ)` 约 15 行 unsafe。

**为什么否掉 termimad（跑过实物，不是凭感觉）**：它做对的是列宽、CJK 宽度、单元格内
折行、长段落按词折、无 `COLUMNS` 时回落 80 —— 恰恰是自己写最容易写对的部分。留给你的
是两处硬伤：手册里 **21 处 `\|` 转义**渲染成 `ENERGY\| Total`（追到解析层确认
`minimad` 正确分列但不剥反斜杠，走 minimad 方案也要自己补这一处），以及表格**没有顶边
和底边框**、样式写死在库里。为把静态文字排版拉进 crossterm/mio/signal-hook/regex
一整套 TUI 事件循环，与本仓为 `quick-xml` 只增 1 个包的判据不成比例。

**架构**（三个扩展点互不干扰）：

```
render(md, width) -> String
├── split_blocks(md) -> Vec<Block>   只看行首特征,不理解内容
│     Block::{ Heading, Para, Code, Table, Quote, List, Raw }   Raw = 认不出的原样输出
├── render_<块类型>(block, width)     一类块一个函数
└── inline(text) -> Vec<Span>        唯一的行内解析
      Span::{ Text, Bold, Code, Link }
```

加块类型 → 改 `split_blocks` + 一个 `render_*`；加行内标记 → `inline` 加一个 `Span`；
改样式 → 只动 `render_*`。手册加新内容、出现没覆盖的语法都**零改动**（走 `Raw`，
后果是「没渲染」而不是「显示错」）。

**硬约束：行内解析必须是一次扫描的 tokenizer，不是 `replace` 链。** 标记会嵌套互斥
（代码里的 `**` 不是粗体，链接文字里可以有粗体），replace 链加到第三种标记就开始互相
破坏、只能推倒重来。这是整个渲染器**唯一一处第一版必须写对**的地方，其余都可先糙后细。

**三个提交**：

| # | 内容 | 行 |
|---|---|---|
| 1 | 骨架：`split_blocks` + `Block::Raw` + 行内 tokenizer + 终端宽度 ioctl + 接进 `page_out` 的 tty 分支 | ~250 |
| 2 | 表格渲染（列宽、`\|` 转义、单元格内折行） | ~140 |
| 3 | 行内 LaTeX（`_x` `^x` `\希腊字母` `\text{}` `\mathrm{}`） | ~50 |

**两条验证**：`ferro doc <每页> > f` 与今天**逐字节相同**（tty-only 的保证，可写成
测试）；渲染输出去掉 ANSI 与边框、归一化空白后**包含原文每一个非空白字符**。

**两条实现细节**（查过，不必重查）：

- **Unicode 下标字母不全** —— 只有 `ₐₑₒₓₕₖₗₘₙₚₛₜ`，`b c d f g q w y z` 没有。
  转不了就**保留 `_x` 原样**（`\tau_c` → `τ_c`），不写成 `τc`（制造错误读法）。
  数字下标 `₀-₉` 与上标 `⁰-⁹ ⁿ ⁱ` 齐全
- **绝不自动识别化学式**，只处理显式 `$...$`。手册里 `-o cmp3`、`43Z43P15A`、
  `set.000`、`sq2` 遍地都是，自动下标会在命令行示例里误导人照着敲

行内公式实测 285 个实例、去重 **182 种**，故必须写成语法规则（`_x` / `^x` /
`\<希腊字母>` / `\text{}` / `\mathrm{}` 四条覆盖全部），不能查表。

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
- **`sort_atoms` 不能无脑套用**（2026-09-22 查证）：DeePMD-kit 在 mixed type 下
  明确**不**按类型排序（`deepmd/utils/data.py` 的 `sort_atoms` 文档写着
  "except mixed types"），因为 `real_atom_types` 是逐帧各异的。而 `collect`
  现在无条件排 —— 两者相遇时要先定 mixed system 的规范序是什么
- 与 `filter` 的关系：`filter` 现在按 system 读回再写出，mixed system 经它一趟必须
  仍是 mixed，否则静默降级成「全 0 类型」的坏数据

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

### ferro job 从轨迹抽帧（2026-08-22 提出，告警已于 2026-08-24 补上）

`job` 现在无条件取第 0 帧（`cmd/job.rs:186`）。多帧输入会告警（见下），但**仍然只
生成第 0 帧的输入**——喂一条 500 帧轨迹得到的还是最没平衡那个构型，只是这次用户
知道了。真正的抽帧 + 批量生成仍待做。

`convert` 的 `--start/--end/--stride/--number` 已把选帧逻辑做进
`Trajectory::select_indices` / `spread_indices`，job 复用即可，不必重写。
真正要定的是**产物命名与批量语义**，与 `convert` 不完全一样：

- job 的 `-o` 现在是完整路径且有默认值（`job.gjf` / `job.inp`），多帧要变成
  `job_0000.inp` 这类。`-o` 的语义已于 2026-09-20 统一（见归档），多帧命名仍待定：
  沿用 `convert` 的 `indexed_path`（序号插在扩展名之前）即可，不必另起一套
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
- ~~**GPUMD/NEP 的 `train.xyz` 导出**~~ **已做（2026-08-26）**：走
  `filter`/`merge` 的 `--type nep`，没有新命令。见归档

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

### 零调用清单（2026-09-12 实测快照，判断待定）

这不是待删清单。**多数功能还没有进入实际使用、没有调试过**，因此没法区分
「需求确定但 CLI 层还没接」与「写代码过程中的遗留」——`cube_jump` 就是前者的
现成例子（429 行实现 + 10 个测试 + 手册页都在位，缺的只是 CLI 入口）。

留这份快照是为了**将来能做 diff**：等这批功能跑起来之后重扫一次，「用起来了」
与「始终没人碰」自动分开。现在删掉任何一条，都是拿判断力换行数。

按 2026-09-12 的实测（727 个生产函数，扫描排除 `#[cfg(test)]`）：

| 类别 | 位置 | 行 |
|---|---|---|
| 整文件零引用 | `ferro-analysis/src/trajectory_analysis.rs`（文件头自称「旧接口，保持编译兼容」；4 个 `pub fn` 里 2 个连测试都没有，2 个只有自己的测试） | 113 |
| 整文件零引用 | `ferro-analysis/src/geometry.rs`（`distance`/`angle`/`dihedral`/`radius_of_gyration`/`bounding_box`，3 个测试；`docs/src/` 零提及） | 104 |
| 整文件零引用 | `ferro-workflow/src/templates.rs`（其中 `orca_sp_template` 是 Ferro 不支持的 ORCA） | 25 |
| 整文件零引用 | `ferro-cli/src/args/corr.rs`（`CorrMode`）· `args/traj.rs:5`（`TrajMode`），0.2.0 前 `--mode` 时代遗留 | 8 + 8 |
| 单位转换链 | `ferro-core/src/units.rs` 的 `convert_length` / `convert_energy` / `convert_time` **生产零调用**（只有自己的测试），连带 6 个私有方法只服务它们。四个 `convert_*` 里只有 `convert_pressure` 活着（vasprun / vasp_outcar / cp2k_out / filter / dataset 共 5 处） | ~55 |
| 被绕过的算法 | `ferro-analysis/src/dft/bader_grid.rs:375 max_neargrid` —— 见下一条 | 20 |
| 取代方案只落实了一半 | `ferro-core/src/network_type.rs:233 display_rank` 生产零调用（`class_rank` 有 3 处）。`progress.md` 说两者一起取代了三个 `*_label_order`，实际只落实了一半 | 8 |
| 基础访问器 | `Frame::geometric_center` / `atom_mut` / `wrap_all` / `is_periodic`、`Trajectory::frame_mut` / `iter_frames` / `time_at`、`Atom::distance_to`、`Table::n_cols`、`compounds::with_density`、`ClusterResult::ids`、`ml/diagnostics.rs:177,182`、`filter::is_clean`、`qn_elements::has_qn`、`charge_grid::lat_dist_i` | ~50 |
| 未使用依赖 | `ferro-workflow/Cargo.toml:8` 的 `serde`（实测删掉后 `cargo check -p ferro-workflow` 通过） | — |

**方法上要记住的一条**：`cargo clippy` 零警告不说明没有死代码 —— 库 crate 里
`pub` 项不触发 `dead_code` lint，上面全部躲过了它。查零调用要按名字逐个 grep
全仓，不能靠编译器。

另一个待拍板的：`Atom`/`Frame`/`Cell`/`Trajectory` 上的 `derive(Serialize,
Deserialize)` 全仓从未序列化过。**但有半条是无条件该做的**：`ferro-io` 与
`ferro-analysis` 的 `nalgebra` 开了 `serde-serialize` feature 而这两个 crate
本身零 serde 用法，即使保留 derive 也该摘掉这两个 feature。

### bader weight 的真空电荷取错（2026-09-20 发现，处置未定）

`bader_weight.rs:265` 用 `volchg[nvols]` 当真空电荷，而该数组在 weight 方法里是
**1 索引**的（上方几行的注释自己写着），取到的是最后一个 Bader 体积。
`bader_grid.rs` 的三条路是 0 索引、真空在 `[nvols]`，那边是对的。

实测（`tests/CHGCAR_2atoms`）：ACF 末尾的 `Total` 报 158.98 e，网格实际只有
105.99 e。逐原子的 `ionchg` 不受影响，被污染的只有 `vacchg` 与由它加出来的
`Total`，以及 ACF 的真空行。

**改一个下标就够，但会改变 weight 方法的产物**，所以没有顺手混进 `-o` 那一轮：
bug 修复的 diff 要能一眼说清「变的都该变」（判据同 `build_avg_frame` 那条）。
修的时候连带做三件事：

- `ferro-cli/tests/bader_reports.rs` 的真空断言把 `Weight` 加回来（现在注释里
  指名跳过它）
- 真空非空的情形现有 fixture 盖不到 —— `CHGCAR_2atoms` 的密度处处高于阈值，
  三条路的 `vacchg` 都是 0，所以这个 bug 只在 `Total` 行上看得出来。要钉住修复
  得再加一个带真空区的 fixture，或给现有的调 `--vacval`
- `dev/bader.md` 未提真空桶的索引基准，修完补一句

### `max_neargrid` 被绕过（2026-09-12 发现）

三件事实，处置未定：

- `bader_neargrid:542` 直接调 `step_neargrid` 把上升循环内联了，因而绕过了
  `max_neargrid`（`bader_grid.rs:375`，20 行）
- 两个兄弟都在用：`max_ongrid` ← `refine_edge:464`，`max_offgrid` ← `bader_offgrid:744`。
  三条路里只有 near-grid 这条不对称
- `dev/bader.md` 有 **3 处**写着 `max_neargrid` 的规格（含 known=1 / known=2
  两种退出条件那条陷阱），删函数要跟着改

**不预先定处置方向**：Bader 是本仓唯一「对着 Fortran 逆向出来」的模块，删还是
让 `bader_neargrid` 用回它以恢复三条路对称，该在真要动这个算法的时候连着算法
一起判断。现在定一个方向，只会让将来的人跳过这次判断。

### 带速度的测试 fixture（2026-09-12 提出，2026-09-21 仍未做）

`tests/` 两条 fixture **没有速度数据**，于是 `vacf` / `vanhove` / `rotcorr`
跑不到写文件那一步 —— 这三个命令端到端从来没验过（0.2.0 改名那轮就只验了
gr/angle/msd/sq/net/map），2026-09-12 的 golden master 基线也盖不到它们。

补一条带 `vx vy vz` 列的 LAMMPS dump（5 帧即可）就能把三个命令纳入。
注意**不能**用全零速度顶替：VACF 在全零速度下是退化情形，跑得通但看不出
实现对不对，那是假的安全网。

### `ferro map jump`：丢失的 CLI 入口（2026-09-03 审计发现）

`ferro-analysis/src/md/cube_jump.rs` **429 行实现完整、10 个单元测试、已从
`md/mod.rs` 导出**，但没有任何 CLI 入口，`lib.rs` 的 `pub use md::{...}` 清单里
也漏了它 —— 0.2.0 把八个 `fe-*` 合并成单个 `ferro` 时，`fe-cube -m jump` 没有
跟着搬过来。

**这不是死代码，是漏登记的待办**：`docs/src/analysis/cube-jump.md` 开头就写着
「仅库函数，没有 CLI 入口……待 `map jump` 补上」，`SUMMARY.md` 也挂着这一页，
即手册已经对用户承诺了。此前 `plan.md` 里却没有这一条。

补的时候按 `map` 现有形状走：`cmd/map.rs:223 run_grid` 的 `drive()` 已经把
「一输入一 `.cube`」的遍历收好了，`calc_cube_jump` 的签名与
`calc_cube_density` / `calc_cube_radius` 同形，接一个分支 + 一个 `JumpCmd`
参数结构体即可。手册页与 `doc.rs` 的 `PAGES` **都已经在位**（`ferro doc cube-jump`
现在就能打出来），缺的只有 `help.rs` 的帮助页与 `print_map_overview` 里的一行。

### VASP：ML_FF 的 OUTCAR 与 io_dispatch 注册（2026-08-27 提出）

两件被这一轮明确划在范围外的事：

- **ML_FF 的 OUTCAR**：VASP 的机器学习力场跑出来的输出用 `free  energy ML TOTEN`
  与 `ML FORCE`。dpdata 用一个 `ml=True` 开关切 token，但它给两者的行偏移是
  `[14, 4]` —— **块结构本身就不同**，不只是换个名字。没有样例就没法验，等有输出
  再补；`vasp_outcar.rs` 的锚点已经是多候选表的形状，加 token 不必重构
- **`io_dispatch` 注册**（2026-09-21 复核：**仍未做，本轮有意绕开**）：把 LAMMPS
  导出挂在 `dataset collect --type inspect` 上而不是 `convert` 上，于是下面这条
  判据一次也没被逼着定。`ferro convert -i OUTCAR -o traj.xyz` 依旧不认识。要动
  `io_dispatch.rs` 的两处 match、`supported_formats()` 那张有测试盯着的清单，以及
  `ferro-python/src/io.rs` 那处独立分派。判据要先定：`OUTCAR` 没有扩展名，而
  `io_dispatch` 现在是按扩展名（加 `POSCAR`/`CONTCAR` 的前缀特例）分派的 ——
  是给它加内容嗅探，还是只按前缀认 `OUTCAR*`

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

### ferro-cli：REPL / 脚本模式（2026-08-11 由高降中，2026-09-15 降低）

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

### ferro-python：pyo3 0.29 运行时验证（2026-09-15 由中降低）

0.21 → 0.29 只做了类型层验证。本机无 maturin，`cargo build` 在 macOS link 阶段过不了
（`extension-module` 需 `-undefined dynamic_lookup`）：

```bash
pip install maturin
cd ferro-python && maturin build --interpreter "$(which python)"
pip install target/wheels/*.whl
```

冒烟测试要覆盖：读 xyz/cif/lammpstrj、`supercell`、`write`、
**新拆的 `gr_pair` / `gr_all`（含 `by="label"`）**、`msd`。

### 是否采用 rustfmt（待定，暂缓；2026-09-15 由中降低）

收益：把格式从代码审查面里移除。代价三点：

1. 一次性约 80 文件的大 diff，`git blame` 在这些行上全部指向那一次提交
2. 会拆掉现有的刻意列对齐（`gr.rs` 的 `WelfordStats`、`angle.rs` 的 `CellList`、
   `units.rs` 的枚举表）
3. 生成文件要显式排除：`rustfmt.toml` 的 `ignore` 仅 nightly 可用，stable 上需给
   `cp2k_basis_db.rs` 的静态表加 `#[rustfmt::skip]`

若采用：单独一次纯格式提交（`style: adopt rustfmt`）+ 固定 `rustfmt.toml`，
不要混进特性或依赖升级提交。

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

### 手册英文化 + 公式统一（2026-09-23 落地）

用户决定 `docs/src` 改英文，理由是手册的读者是外人；`//` 内部注释与测试断言**保持
中文**（它们承载判据，母语更准，读者是开发者）。规则进 `CLAUDE.md`「编码约定」，
`dev/overview.md` 那段原是它的副本、改为指针。

八个提交：产物摘出 git → 规则 → 术语表 → 公式统一 → 四批翻译。
`docs/src` 25 页 4705 行现为**零中文字符、零中文标点**（894 行中文：散文 555 +
表格 290 + 代码块 49；89% 集中在 `cli-reference.md` 516 行与 `network.md` 277 行）。

定下来的判据：

- **术语表先行**（`dev/glossary.md`，65 条，用户逐条裁决）。不先定表，两页大的会译出
  两套词汇。用户的四条裁决：口径 → statement、约定 → regulation、点名拆成
  declare/name 两处（另一处「位点名」是 n-gram 断词错误）、规范半边 = canonical half
- **公式统一 `$...$`，不手写 Unicode 上下标**，但**单位符号例外**（`Å³` `Å⁻¹`
  `g/cm³`）—— 手册原本就是「公式 LaTeX、单位 Unicode」（`sq.md:38` 即
  `` $q \approx 1$–2 Å⁻¹ ``），把单位塞进 LaTeX 只会让表格更难读。`SUMMARY.md` 的
  目录标题也不动（mdBook 侧边栏不过 MathJax）
- **`docs/book/` 摘出 git**。91 个 HTML 的时间戳停在 08-13 而源文件已到 09-22，
  跟着仓走只会让每个手册提交带一批陈旧产物

**被实测推翻的**：

- **「公式在 mdBook 里已正常渲染」是错的** —— MathJax 2.7 默认定界符是 `\(...\)`，
  手册里 255 行 `$...$` 在网页上**从来没被渲染过**，一直是字面的 `$Q_n$`。新增
  `docs/theme/head.hbs` 放开 `$` 才修好；`workflow/spin.md` 的 `\(...\)` 是全书唯一
  渲染正确的写法，已统一为 `$...$`
- **中文标点是盲点** —— `，。、（）「」` 不在 `[\u4e00-\u9fff]` 里，「零中文残留」
  的检查全部漏掉它们。最后一处是替换顺序的锅：先把全角括号换成半角，`），` 那条规则
  就失配了
- **校验脚本的判据改了三次才对**（`.claude/check_tr.py`，随 TASK 文件同去留）：
  ① 行内代码逐项相等 → 把 `gr_<配对>` → `gr_<pair>` 报成违规；② 只比不含中文的项 →
  新版翻译出的英文项成了「多出来的」；③ **旧版的英文字面量必须仍然存在**（多重集包含，
  不是相等）才对。契约保护的是英文字面量，不是占位符。顺带堵了一个假阴性：在子目录里
  跑 `git show` 拿到空串，空集减空集恒过

**顺带发现未修**：`cli-reference.md:718` 写着 bader「没有 `-o`」，而同节 744 行写着
「`-o <DIR>` 收集到别处」—— 0.3.3 给 bader 加了 `-o`/`-s`，这一行漏改。翻译时**照原文
直译**没有顺手修：翻译提交要能对拍，修 bug 该是独立一个提交。

### cp2k_grep_strInfo.py 的能力合并：LAMMPS 导出 + 逐帧标量（2026-09-21）

用户在 `private/` 放了一个 626 行的脚本要求合并。**先做的是对账**：它的九项能力
里 Ferro 已有六项且实现更稳健（VASP/CP2K 读取、DeePMD 写出与 8:1:1 划分、
lammpstrj/data 导出），`outputLammpsin` 甚至引用了一个**根本没定义**的
`LAMMPSIN` 变量。真正缺的只有三样：CP2K 多文件布局、DeePMD mixed type、
逐帧标量序列。

**用户选的范围是 LAMMPS 导出 + 逐帧标量**，另两样留在优先级高那一栏。

六个提交，顺序即依赖：三斜 bug → `Frame` 字段 → 三个 reader 填温度 →
`collect` 的 `-o` 与 `.db` → `--type inspect` → 文档。

**被实测推翻或修正的**：

- **原以为缺的是 LAMMPS writer，实测 Ferro 的两个 writer 已支持三斜**，
  `cell_to_lammps` 与脚本的 `getLammpsCell` 逐行同构。真正的断点是**入口**：
  `io_dispatch` 不认 `.out`/`OUTCAR`/`vasprun.xml`
- **却发现了一处脚本对、Ferro 错的**：dump 的三斜盒子行是 `*_bound` 不是
  `xlo/xhi`，而 Ferro 读写两侧一致地错，往返测试全绿。用户当时的指示是
  「只取能力，不记录这些问题」，但这一条涉及 Ferro 本身，提出后用户判为「修，
  且记录」。**教训是对账要做到实现层，不能停在能力清单**
- **原计划的 `time_fs` 列去掉**：设计时说「step 与 time_fs 是 CP2K 锚点上现成的」，
  这话对**文件**成立、对 `Frame` 不成立 —— 两者都没被带出来。实测三条路里
  `step` 有两条、`time_fs` 只有一条，故只加 `Frame.step`
- **脚本在 VASP 路抓的温度一度被怀疑是错的**（正则 `temperature +([\d\.]+) K`
  命中的是 `EKIN_LAT` 行）。实测**它是对的**：那个括号标的是离子温度，
  `EKIN=4.268142` / `N=296` 反算得 111.56 与打印的 111.55 吻合
- **`-o` 的形态被用户改过两轮**：先定「`--type inspect` 时 `-o` 改写 inspect 位置」
  （与 bader 同构），再改成「`-o` 完全交给数据集，inspect 位置固定、给了 `-o`
  报错」。`collect` 的 `-o` 也从必填改为可选
- **后缀从 `.train` 改成 `.db`**：用户先说可能不需要后缀，最后定为全部加 `.db`

**明确划在范围外**：CP2K 多文件布局（`-pos-1.xyz`/`-frc-1.xyz`/`-1.cell`）·
DeePMD mixed type · `io_dispatch` 注册 · 抽帧参数 · `ferro info` 的「只报首尾
两帧」。脚本里**不照搬**的四处：体积用对角线乘积、CP2K 固定行偏移（只认
NVT/NPT_I）、质量表读 `~/.emacs.d` 的 elisp、倾斜折叠的除数写错。

### 原子顺序 + --format + 帮助页精简（2026-09-22 落地）

起因是实测 `collect -i 1Al/*.out` 报「different compositions」而两个化学式打出来
一模一样。落地清单见 `overview.md`「2026-09-22 的第二批」。

**被实测推翻的两条原计划**：

- **「删掉嗅探，完全依赖 `--format`」**（用户最初的要求，理由是版本多、易错、
  难维护）。查证后留下了：`sniff` 40 行且不含任何版本知识，担心的那部分全在四个
  reader 的块匹配里。留着它换来一行 NOTE 与那条一目录一格式的守卫
- **「`--format` 必填」**（我最初的推荐）。用户改为「有 format 用 format，没有就
  嗅探，探不到才报错」，旧命令行因此一条都不用改

**动手前查了本机 deepmd 环境**，两条结论都进了 `overview.md`：DeePMD-kit 加载时
自己按类型排序（所以磁盘顺序对训练无影响），dpdata 的 `append` 在 atom_types
不一致时自动排两边（所以这是领域内的既有做法，不是 ferro 的发明）。

**划在范围外**：`filter` / `merge` 的 `--format`（它们读 npy，不需要）·
`io_dispatch` 注册 `.out`（第三轮绕开）· mixed type 下的规范序（见优先级高那栏）。

### cp2kdata 的能力合并：CP2K 单点读取（2026-09-22 落地）

核心需求是**单点能数据的读取**，用于建训练集。落地清单见 `overview.md`
「2026-09-22 的一批」。

**明确划在范围外**（本轮不做，也未计划）：
- **GEO_OPT / CELL_OPT**。看着像"顺手"，其实不是：cp2kdata 里这两条路的坐标
  **不来自 `.out`**，`parse_cell_opt` 去 glob `*pos*.xyz`，`parse_geo_opt_info`
  干脆不产出坐标。做它们等于做 CP2K 多文件布局，那是另一件已划在范围外的事
- cube / pdos / Mulliken / Hirshfeld 布居 / 振动频率 / 偶极 / DFT+U 占据数。
  前两项 Ferro 已有自己的实现，后五项 `Frame` 里没有存放处
- `io_dispatch` 注册 `.out`（`convert` 仍进不去，与上一轮同样绕开）
- 2023 / 2024 的**单点** fixture。手上那两个版本的样例全是 AIMD，单点布局是
  按同版本 MD 输出推断的，因而落在"打提示"那一侧而不是"已验证"

**仍然开着的**：
- **kind 分成独立训练类型**。现在 `type_map` 恒按 element，`Fe1`/`Fe2` 合成一个。
  dpdata 默认相反。要分的话是 `collect` 加一个开关，不是改默认
- **2026 版本未经实测**，是按用户"格式沿用 2025，未修改"的说法列进已验证的。
  拿到 2026 的输出跑一遍即可确认

### 三项小待办：&Path、bader --outdir、-o 语义统一（2026-08-13 提出，2026-09-20 落地）

三件事原本互相独立，做的时候合成了一条链：路径类型统一之后，`-o` 的语义才好改。
六个提交，顺序即依赖。

**动手前推翻的原计划**：原案是「`-o` 统一为文件名、路径走 `--outdir`」，收敛到
**两种**语义（后缀 / 目录）。用户指出 `-o`/`--outdir` 本身就乱，要求改成
pathlib 那样「带目录的路径和单个文件名都收」。于是目标从「两种语义」变成**一条
规则**：

> `-o` 恒为路径。一次运行写多个产物的命令（13 个分析命令、`chg-sdf`、`bader`、
> `dataset`）指**目录**，只写一个文件的（`convert`、`job`）指**文件**。
> 批次标记腾给 `-s/--suffix`，`--outdir` 删除。

定下来的判据：

- **目录不存在时先问一句**（stderr、`[y/N]`、默认否），非交互环境直接报错并指明
  `--mkdir`。自动创建会把打错的路径变成一次「看着成功」的运行；而在脚本里连问都
  没人能答，所以非交互不能沉默地二选一
- **「只给目录」只认用户写的尾部分隔符**，不看磁盘上有没有同名目录。否则同一条
  命令行的含义会随磁盘状态漂移（今天 `out` 不存在→产物叫 `gr_out.csv`，明天有人
  建了 `out/`→写进目录），与「不按输入文件数分派两条路径」是同一条理由
- **`convert` 拒收以分隔符结尾的 `-o`**：目标格式正是从文件名推断的，目录名里没有
  它。`job` 相反，它有默认名（`job.gjf`/`job.inp`/`pw.in`），可以放进目录
- **`bader` 的报告默认写在输入旁边**，是全仓唯一不默认当前目录的命令。VASP 的
  电荷密度一律叫 `CHGCAR`，默认当前目录时 `run1/CHGCAR` 与 `run2/CHGCAR` 必然
  互相覆盖 —— 换个默认值就消掉了这个坑，比加参数让用户每次记得写更靠得住
- **`write_acf`/`write_bcf`/`write_avf` 改为返回文本**（`acf_text` 等），文件由 CLI
  写。既然要动签名，顺手清掉「分析层不碰文件系统」的唯一例外
- **reader 也跟着 writer 一起改 `&Path`**。只改 writer 会让 `ferro-io` 一半 `&str`
  一半 `&Path`，`io_dispatch` 两侧各转一次

**dataset 的默认后缀**（同一轮，用户要求）：`filter` / `merge` 的产物默认带
`.train`，划分出的三部分是 `.train`/`.valid`/`.test`，`collect` 不加（它的产物是
没筛过的原始数据）。连带三条：名字已带后缀不叠加；`merge` 遇到混合后缀**报错**而
不是默认标成训练集；`.train` 输入允许再划分（先剥后缀），`.valid`/`.test` 拒绝。

**验证**：53 个产物的逐字节对拍贯穿六个提交，只有 dataset 那步的目录名多出
`.train`，内容零差异。bader 另造了 `tests/CHGCAR_2atoms`（8×8×8 合成网格，9 KB）
补上端到端覆盖 —— 此前那条链只有内存里构造网格的单元测试。

**顺带发现的缺陷（未修）**：`bader_weight.rs` 的真空电荷取 `volchg[nvols]`，而该
数组在 weight 方法里是 1 索引的，取到的是最后一个 Bader 体积。ACF 末尾的 `Total`
因此虚高（新 fixture 上 158.98 e vs 实际 105.99 e）。`ongrid`/`neargrid` 走另一份
代码，不受影响。见 `issues.md`。

### 近距接触结构的定向构造（2026-09-01 提出，2026-09-15 废弃）

起因：MD 里原子互相穿透，训练集在阳离子对（P-P、Al-Al）短程区没有样本，势函数
在那里纯外推。原计划在 `ferro-structure` 里挑一对原子定向压缩、造扫描结构送 DFT。

**废弃理由**：改由 dpgen 主动学习去探索这些构型，Ferro 不做定向构造。原计划里
「压近的帧会被 `filter` 的 `f_max` 删掉」这条冲突随之消失 —— dpgen 产出的结构
正常流程下不再过二次筛选。

留下的事实（重开或做别的扰动时不必重查）：

- **六个工具（ASE `rattle`、pymatgen `perturb`、dpdata `perturb`、dpgen
  `create_random_disturb.py`、hiphive `mc_rattle`、ASE GA）无一做定向压缩**，全是
  各向同性随机抖动。管最小间距的只有 hiphive 与 ASE GA，方向是**阻止**近距接触
- **dpdata 的 `perturb` 是 dpgen 那个脚本的修正版**：dpgen 的方向按立方体归一化、
  偏向 ⟨111⟩，模长 U[0, dmax) 也不是球内均匀；dpdata 两处都改对了。要统计微扰用
  dpdata，别用 dpgen 的 `pert_atom`
- **排斥壁两条路配合使用**：ZBL 混入管 d→0（NEP 内置 `zbl`，DeePMD 走 `use_srtab`
  + `sw_rmin`/`sw_rmax`），数据管约 1.5–2.7 Å
- **存疑**：当时的锚点表（P-P 2.73 Å、Al-Al 2.83 Å 以下无样本）标为「训练集见过的
  最近 r」，实测对象却是 `tests/43Z43P15A_NPT_5.lammpstrj`（5 帧）。引用前先在
  真实训练集上重测

### VASP AIMD 读取：OUTCAR + vasprun.xml（2026-08-27，0.3.2）

原计划只做 OUTCAR，用户要求**两个都做**。实测 `quick-xml` 与 `roxmltree` **各自只
净新增 1 个 crate**（零传递依赖），代价不是问题；选 `quick-xml` 的**流式**是因为
DOM 的内存约为文件的 5–10 倍，而这条链的输入本来就是几百 MB。

**被查证推翻的原计划**（两条，都在动手前）：

1. 上一轮写的「dpdata 取 `energy(sigma->0)`」**是错的**。`dpdata/formats/vasp/
   outcar.py:118` 取的是 `free  energy   TOTEN`，用户的 `private/dp_makedataliu.py`
   也是（`lines[-3].split()[-2]`）。理由本身也站得住：力是自由能的导数
2. 原计划「同目录既有 OUTCAR 又有 vasprun 时优先 vasprun、回落 OUTCAR」被用户
   否掉 —— 改为**由用户指名文件**（跟 CP2K 端一样），程序不做目录级偏好判断。
   于是分派改成读头几行认横幅，而不是一张扩展名特例表

**定下来的判据**：

- **格式按内容嗅探,不按名字**。VASP 写的叫 `OUTCAR`（无扩展名），用户手上那份叫
  `50Z50P_0.970_3000K.outcar`，`.out` 又谁都可能用 —— 按名字判必然要开一串特例，
  而横幅是唯一的
- **同目录混格式报错**。真实运行目录里两个文件并存且记同一批帧，按 collect
  「同目录 = 同一次运行的分段」拼接会让数据集**悄悄翻倍**：成分一致、两个文件
  各自也都读得通，不会有别的症状
- **晶胞逐帧各读各的**，缺块的帧丢弃而不继承上一帧。定胞下继承与不继承逐位相同，
  这个 bug 只会在第一次跑变胞时显形 —— 与 `array_order.rs` 的转置同族
- **`in kB` 是 VASP 自己的 Voigt 顺序 `XX YY ZZ XY YZ ZX`**，与 extxyz 规格的
  `XX YY ZZ YZ XZ XY` 不同。三个独立来源一致（dpdata 的下标映射、用户脚本的
  `[[0,3,5,3,1,4,5,4,2]]`、ASE 的 `[[0,1,2,4,5,3]]` 重排）
- **体积用 `|det|`**，不用三个对角线相乘（用户脚本是后者，只对正交胞成立）
- 两条 VASP 路径的**收敛判据不同**（OUTCAR 有 VASP 自己的 `EDIFF is reached`，
  vasprun 只能数 `<scstep>` 与 `NELM` 比），故判据随 stats 打进报告 —— 否则同一次
  运行换个来源、丢帧数不同会没人说得清

**实测发现的两处坑**（都不会报错，只会给出错的数）：vasprun 的 `<energy>` 每个
`<scstep>` 都有一份，收敛值是最后一个（取第一个得到 33072.96 这种数）；
`<array name="atoms">` 的 `<field>` 表头也是文本，当成元素会得到 298 个物种、
每帧判 incomplete、整个文件读成空 —— 后者是实际写代码时踩到的。

**验证**：`tests/` 两份从真实文件裁出的 fixture（OUTCAR 340 KB、vasprun 188 KB，
dpdata 仍能正常读）；`#[ignore]` 的全量对拍（2000 帧 / 425 帧）与 dpdata 逐位一致，
virial 差 8.5e-9 相对量 —— dpdata 用的 eV/Å³→GPa 常数是 160.2176621，`units.rs`
用 CODATA 2018 的 160.2176634。用户那份 OUTCAR 本身是**两段拼的**
（1..1575 接 1..425，第 1576 步被中断），正好验到重启计数与截断帧丢弃。

### GPUMD/NEP 导出：`--type` + train/valid/test 划分（2026-08-26，0.3.2）

**原计划是加一个 `ferro dataset export` 子命令，被用户否掉** —— 改为在 `filter` 与
`merge` 上加 `--type`，划分标签也加在这两处。事后看这个改法更对：NEP 的 `train.xyz`
本来就是 extxyz，ferro 的 writer 早就能写，真正缺的只是「读 DeePMD system 目录」，
而 filter/merge 本来就在读；新命令会把同一件事再写一遍。

定下来的判据：

- **一个 system 一个 `.xyz`**，不并成一份 train.xyz。NEP 要的是单份，但那是一句
  `cat`；分开保留的是「抽掉某个来源不必重跑」的能力
- **划分产物用目录名后缀** `.train/.valid/.test`，不造 `train/` 子目录 ——
  `SPLIT_SUFFIXES` 早就在仓里，merge 一直在继承它，这是 dpgen/dpdata 的约定
- **成员取自打乱序**：直接切尾巴的话 test 全是轨迹末尾一段连续状态。各部分内部
  再排回帧序，同 seed 下逐字节可复现
- **比例向上取到至少 1 帧**：给了比例却拿到 0 帧，等于静默地没有验证集
- **一个 `--ratio 8:1:1` 而不是两个 ratio 参数**（用户第二轮要求）。数字是
  **权重不是分数**，`8:1:1` 与 `80:10:10` 同义，不必凑成和为 1；两段写法 `9:1`
  按 `train:test` 解释 —— NEP 要的就是这两份文件，而位置歧义靠每次运行按名字
  打出三部分的帧数来自证，第一行输出就能看见解释对不对
- `nep` 与 `extxyz` 两个值只差一个应力键（`virial=` / `stress=`）。留两个值而不是
  一个，是因为「按下游训练框架命名」比「按文件格式命名」更接近用户提问的方式

**没做、也不该在这里做的**：`weight` / `dipole` / `pol` / `bec` 这些 NEP 可选键。
`Frame` 装不下它们，硬塞要开一条绕过 `Trajectory` 的通道，与「额外键搬运」是同一件
待办。

**已知局限（写进帮助页与手册，不是缺陷）**：同一段 MD 的相邻帧高度相关，帧级 test
仍偏乐观；诚实的估计要按整个 system 划分，做法是把来源分开筛、留一个不动。

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
