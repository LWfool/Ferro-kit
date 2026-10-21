# 挂起组：Bader / cube_sdf·chg_sdf / map jump

> 这三组等对应功能改写时一并处理，平时不修。从 `dev/plan.md` 按需读取；定级规则见 `plan.md` 顶部。

## 优先级高

### ferro map chg-sdf 的 --cubes 拆成单文件（2026-08-11 提为高）

现在是多个 cube 聚合成**一张** SDF（`cmd/map.rs::run_chg_sdf`），与 `ferro map` 其余
模式「一个输入 → 一个 `.cube`」的语义相反，`--cubes` 也是 `-i` 之外唯一的输入参数。

阻塞点不在遍历而在**中间产物格式**：跨文件加权平均不可交换，逐文件出产物之后若还想
聚合，产物里必须带样本计数（否则 5 帧的 SDF 与 500 帧的 SDF 被等权平均）。
故顺序是：先定义带计数的中间格式 → 再拆 `--cubes` → 再接批处理遍历。


## 优先级中

### 2026-10-02 物理审查的三条严重问题（用户定为中，功能后续会改写，届时一并修）

三条都已用数值实验证实（临时 crate，不在仓库），均未修。

**S-A　Kabsch 旋转求成了转置**（`md/cube_sdf.rs::kabsch_rotation`）。H = Σ m·rᵀ =
U S Vᵀ 时，使 R·m ≈ r 的解是 R = V·D·Uᵀ，代码返回 `u * diag * vt`（= Rᵀ）。
`map sdf` 与 `map chg-sdf`（`dft/chg_sdf.rs` 复用，`rotate_grid` 的 pull 插值同错）
叠加的是反向旋转的团簇，置换枚举也按错误 RMSD 选。实测同一畸变 PO₄ 刚性旋转两次，
RMSD mean 1.47 Å（应 ≈ 0）。现有 `test_kabsch_known_rotation` 点全在 xy 平面、绕 z
转，H 退化且旋转与 diag(1,1,0) 可交换，两种写法结果相同，所以没抓到。修法：
`v * diag * u.transpose()` + 三维非退化回归测试（随机旋转 + 四面体，RMSD < 1e-10）。

**S-B　Bader 三条网格路的真空电荷恒为 0，真空里冒出大量伪体积**（`dft/bader_grid.rs`
on-grid `:213-232`、near-grid `:580-597`、off-grid `:763-770`）。`volchg.push(0.0)`
后立刻把这个 0 当 `vacchg` 读出；且真空在梯度上升**之后**才标，平坦/噪声真空区的
每个局部极大都成了独立 Bader 体积再分给最近原子（与 `dev/bader.md` §4 的伪代码顺序
相反）。实测 20 Å 盒两高斯峰 + 5e-4 e/Å³ 背景：网格总电荷 17.61 e，三条路 Total 都是
13.65 e，nvols 5 万量级；`vacvol` 是对的。修法：上升前 `mark_vacuum`，路径跳过
`volnum == -1`；`vacchg = Σ_{vac} rho / nrho`；补带真空区的 fixture。

**S-C　Bader weight 不守恒电荷，权重也不是 Yu–Trinkle 的**（`dft/bader_weight.rs`
`:64`、`:138-140`、`:162`）。实测两峰无真空，网格总量 15.925 e：near-grid 严格守恒并随
网格收敛，weight 丢 0.1–0.3 e 且不收敛。机理（推断）：边界点用负数存 basin，若某点所有
上游邻居都是同一个边界点，它被判为内部点并继承负 basin，不在任何 `neigh` 列表里，电荷
整体丢失。另：通量系数应为 Voronoi 面积 / |R|，代码用 |R|；WS 筛选严格不等式把面积为 0
的棱/角向量也留下了；`ionvol` 不计边界点；真空点参与流动。修法：真空标记与 basin 号
分开存、按 Voronoi 面积算 α、先标真空再流动。与下面「weight 的真空电荷取错」一并处理。

### bader weight 的真空电荷取错（2026-09-20 发现，处置未定）

`bader_weight.rs:265` 用 `volchg[nvols]` 当真空电荷，而该数组在 weight 方法里是
**1 索引**的（上方几行的注释自己写着），取到的是最后一个 Bader 体积。
`bader_grid.rs` 的三条路是 0 索引、真空在 `[nvols]` —— 下标对，但那里是硬编码的 0（见上面 S-B）。

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

### `ferro map jump`：丢失的 CLI 入口（2026-09-03 审计发现）

`ferro-analysis/src/md/cube_jump.rs` **429 行实现完整、10 个单元测试、已从
`md/mod.rs` 导出**，但没有任何 CLI 入口，`lib.rs` 的 `pub use md::{...}` 清单里
也漏了它 —— 0.2.0 把八个 `fe-*` 合并成单个 `ferro` 时，`fe-cube -m jump` 没有
跟着搬过来。

**这不是死代码，是漏登记的待办**：原手册页 `docs/src/analysis/cube-jump.md` 开头就写着
「仅库函数，没有 CLI 入口……待 `map jump` 补上」，`SUMMARY.md` 也挂着这一页，
即手册已经对用户承诺了。此前 `plan.md` 里却没有这一条。

补的时候按 `map` 现有形状走：`cmd/map.rs:223 run_grid` 的 `drive()` 已经把
「一输入一 `.cube`」的遍历收好了，`calc_cube_jump` 的签名与
`calc_cube_density` / `calc_cube_radius` 同形，接一个分支 + 一个 `JumpCmd`
参数结构体即可。缺 `help.rs` 的帮助页、`print_map_overview` 里的一行，以及手册页。

**手册页 2026-10-21 撤出手册**（瘦身 D-4，用户定）：没有 CLI 入口的功能不占手册页。原页原样存在
`dev/plan/cube-jump-manual.md`；补入口时按 `add-analysis` 技能把它搬回 `docs/src/analysis/`，
topic 名改用 `map jump`，并改正其中过时的用法示例（`-o run1` 是 0.3.4 以前把 `-o` 当批次标记的写法，
现在是 `-s`），再登记 `SUMMARY.md` 与 `doc.rs` 的 `PAGES`。

### `ferro map sdf` 导出旋转后结构：压缩 npz（2026-09-26 提出）

用 `ndarray-npy` 的 `NpzWriter::new_compressed`（默认 feature `compressed_npz`，
依赖全为纯 Rust，Linux/Windows 无系统库要求）——**所以 ndarray-npy 不关默认
feature**。选型依据：

| 格式 | Python 端 | Rust 端代价 | 结论 |
|---|---|---|---|
| npz（deflate） | `np.load`，只需 numpy | 已在树里，新增 0 | 采用 |
| npy + zstd | 另装 `zstandard` 并自拼 header | zstd 是 C 绑定 | 否 |
| HDF5 | `h5py` | 链 libhdf5，跨平台最难 | 否 |

体积主要由精度决定而非算法：f64 坐标尾数近随机，deflate 只压 10–30%；**存 f32
直接减半**（坐标误差 ~1e-5 Å，dpdata 默认也是 float32）。建议默认 f32 + 开关留
f64。数组命名与形状在实现时定。

## 历次审查登记的挂起条目（2026-10-21 从各审查节归拢，原文）

### 归类（10-19 遗留清理）

**归类（只改登记，不改代码）**：
- **Bader 改写组**：D-8、D-9、D-10、refine_edge 的负号与真空哨兵相撞、Bader 参数校验与手册不符、`dev/bader.md`
  过时与空断言、`bader_grid` 三份重复的后处理（G-7）、布局建议 5、C-D2、C-D3、Bader 归属的 MIC（口径）、
  near-grid / off-grid 收敛测试、weight 方法 BCF 与 ACF 的一致性测试
- **cube_sdf / chg_sdf 改写组**：B-3、B-4、B-6、B-10、D-11、`cube-sdf.md` 过时、chg_sdf 插值的 `rem_euclid` 越界、
  chg_sdf 跨文件平均按体积加权、chg-sdf 的参数校验 / `-s STEM` 帮助 / 中文报错、`calc_chg_sdf` /
  `calc_cluster_sdf` 返回 `Option`、≥8 个 P 的置换测试、**`NetworkGraph.former_qn` 同名不同义**（实为桥氧数，
  与 `AtomType::Former.qn` 不是一个量；唯一使用方是 cube_sdf 的族名）
- **`map jump` 入口组**：B-9

### 2026-10-17 全盘架构审查

| # | 级 | 位置 | 问题 | 修法 |
|---|---|---|---|---|
| A-3 | 轻 | `dft/bader.rs:165,194,217` | ACF / BCF / AVF 文本在 analysis 里拼，不走 `write_table` | 随 Bader 重写：改 `to_tables()`，或承认是外部格式搬 `ferro-io/src/writers/bader.rs` |
| B7 | 中（实测） | `charge_grid.rs:54-55` | `ChargeGrid::new` 对退化胞 `expect` panic：第三晶格矢量为 0 的 CHGCAR 让 `ferro bader` 退出码 134。B-1 / D-M17 修了分析侧，漏了 core 构造器 | 返回 `Result` |
| B8 | 轻（实测） | `readers/chgcar.rs:103-107` | CHGCAR 网格维数为 0 时 bader 输出 0 e、退出码 0（cube reader 已有此检查） | 「维数 ≥ 1」放进 `ChargeGrid::new`，一处守住两个 reader |
| D-C1 | 中（实测） | `cube_sdf.rs:383,161-166,634-638,169-174`；`chg_sdf.rs:382-390` | 形成子硬编码「P」：`--former Si` 时族名仍 `P1`、cube 头里 Si 写 Z=15；`--modifier Pb` 时 Pb 写 Z=15，Q1 下 Pb 被 `p_indices` 当骨架参加 n! 置换（Pb 版 2 分 19 秒未完）；`--ligand` 非 O 时元素仍写 O | 并入 cube_sdf 改写组（下面 R-2 的结构化 `Site`） |
| D-A3 | 中 | 各 `calc_*` | 返回类型四种：`Result`、`Option`（cube_density / cube_jump / chg_sdf，丢失败原因；`cube_jump.rs:136` 把 `check_same_atoms` 的错误 `.ok()?` 吞掉）、`Result<Option>`（cube_sdf）、不可失败（Bader；`BaderParams` 无 `validate`，`bader_grid.rs:722` 无胞帧 `unwrap` panic） | 统一 `Result<T>`，结果带 warnings（第 3 步） |

D 组改写组的新结构（只给结构，条目见原组）：
- **cube_sdf / chg_sdf（R-2）**：① 逐帧并行提取团簇，沿键 BFS 展开、检测逾渗（B-4）；位点用结构化枚举
  `enum Site { Former{elem,qn}, Ligand{elem,kind}, Modifier{elem} }`（消 D-C1）② 按签名分族，族内确定性选参考
  ③ 独立 `align.rs`，Kabsch 写对（S-A），置换限组内、组大时匈牙利 + Kabsch 迭代（B-6）④ 网格以参考结构为中心、
  以 Å 定界（B-3），归一化定为「每团簇每 Å³ 数密度」（B-10）⑤ 告警进结果
- **Bader（R-3）**：① 一个驱动 + 三种上升策略，真空先于上升标记（S-B），`refine_edge` 按方法分派、看 26 邻（D-8、
  `max_neargrid`）② 体积编号用 `enum Basin { Unassigned, Vacuum, Vol(u32) }`，去掉 0 / nvols+1 / 负数三种哨兵 ——
  C-D2、S-C、refine 撞真空哨兵三条都出在这里 ③ 四种方法共用一份后处理，守恒性断言 `Σionchg + vacchg == 网格总量`
  ④ `BaderParams::validate`，`run()` 返回 `Result` ⑤ **`dev/bader.md` 写的参考实现 `examples/bader/` 本地不存在**，
  重写前须把 Henkelman 源码或二进制放回作对拍基准

- 10-17 第 3 步挂起：E2-1 与 D-A1 的 cube_sdf / chg_sdf 逐团簇告警、cube_sdf / chg_sdf / Bader 的返回类型（随改写组）
- E 组 R-6（低）：`bader -i` 改收多文件走 `map_inputs`

### 2026-10-11 全库审查

| # | 位置 | 问题 | 复现 | 修法 |
|---|---|---|---|---|
| D-8 **挂起：Bader 重写** | `dft/bader_grid.rs:474`（refine_edge 内 `max_ongrid`）、`:439-449`（边缘只看 6 个面邻居） | **CLI 默认的 near-grid（`-r -1`）与 on-grid 逐位相同**：refine_edge 对边缘点跑 on-grid 上升，边缘判定也只看 6 邻（规格与 Henkelman 2009 为 26 邻、near-grid 时调 near-grid 上升）。即默认路径实际是 on-grid | 10 Å 立方胞两高斯（中心 (3.1,3.3,5.0) σ0.9 幅 1.0；(6.3,6.9,5.2) σ1.4 幅 0.6），N=60：on-grid 与 near-grid `-r -1` 均 11.91354，near-grid `-r 0` 11.88225；N=120 默认误差 0.024 e 且随加密几乎不降，`-r 0` 1e-4 | refine_edge 按方法分派，near-grid 调 `max_neargrid`（顺带解决「`max_neargrid` 被绕过」一节）；边缘判定改 26 邻；加收敛测试 |
| D-9 **挂起：Bader 重写** | `dft/bader_grid.rs:658`（`grad_lat_1`）、`:736`（`start = [i+0.5, …]`） | off-grid 三线性插值的 $\partial\rho/\partial f_2$ 写成 x 向差分，起点偏半格 | 同 D-8 体系 `-m offgrid`：N=60 nvols **21918**（应 2）；N=30 3064；原子电荷误差 0.1–0.19 e | 改 `grad_lat_1`、起点 `[i,j,k]`；单测「线性密度 ρ=y 的梯度只有 y 分量」+ 收敛测试 |

| # | 位置 | 问题 | 修法 |
|---|---|---|---|
| B-3 **挂起：cube_sdf 改写** | `md/cube_sdf.rs:614-627 build_reference_frame` | `map sdf` 参考结构多加一次 `(min+max)/2` 平移，与密度网格错位；现被 S-A 掩盖（实测仅 0.04 Å），S-A 修后显形 | 直接用 `pos`，删 `center`/`box_size` 参数；随 S-A 修，测试「参考原子落在各自密度峰上」 |
| B-4 **挂起：cube_sdf 改写** | `md/cube_sdf.rs:401-404` | 以锚原子最小镜像取局部坐标，跨度超过半盒或逾渗的团簇被折断，无告警（推断，未构造算例） | BFS 展开 + 逾渗检测；结果计数被跳过的团簇。随 B-6 / S-A 一并做 |
| D-10 **挂起：Bader 重写** | `dft/bader_weight.rs:107-113`、`:146-148` | S-C 的补充：weight 在等密度平台上 `tsum = 0` → 0/0 = NaN，电荷成片丢（floor=0.2 时 208.6 e 剩 97.2 e） | 随 S-C / Bader 重写；回归测试 Σionchg + vacchg == 网格总量 |
| D-11 **挂起：Bader 重写** | `docs/src/analysis/chg-sdf.md:85`；`dft/chg_sdf.rs:10`、`:345-347 build_family` | 手册教「除以体素体积」，应为晶胞体积（差 $N_1N_2N_3$ 倍）；写出的 cube 不是 e/Bohr³ | `build_family` 除以 V_cell_Bohr 还原 e/Bohr³ 再写；与「chg_sdf 跨文件平均 ρ·V_cell」口径一起定 |

- **B-9**（挂到「`ferro map jump`：丢失的 CLI 入口」，补入口时一并定）：`cube_jump.rs:192-217` 滑窗重叠，τ > 1 时同一次
  跳跃计 τ 次，手册却称「raw count of jump events」。二选一：窗口不重叠（步长 τ），或改称「超阈值窗口数」
- **B-10**（挂到 B-6，cube_sdf 改写时一并定）：`cube_sdf.rs:607-609` 每张网格按全体体素均值 = 1 归一，含 padding 空体素，
  `--padding` 3 → 6 Å 数值降为 1/2–1/3，不同类型不可横比。TRAVIS / VMD volmap 按每帧每团簇数密度或相对体相密度归一

- `bader_grid.rs` 三条网格路的后处理逐段复制（`:161-`、`:508-`、`:719-`，15 个连续 6 行窗口相同；G-7），S-B 修法要改三遍 →
  随 Bader 重写抽出

- 10-11 轻级 D 中归 Bader 组的：refine_edge 取负标边缘，体积 1 取负后是 −1 与真空哨兵相撞、凭空多一个体积（`bader_grid.rs:478-487`，`ionchg` 不受影响）；Bader / chg-sdf 参数校验与手册零碎不符；`dev/bader.md` 过时处与 Bader 单测的空断言

### 2026-10-03 四路审查

| # | 位置 | 问题 | 复现 | 修法 |
|---|---|---|---|---|
| C-D2 **搁置**（2026-10-08 用户定：随 Bader 重写一并修，同 S-B / S-C）| `dft/bader.rs:205` vs `bader_weight.rs:98,145` | weight 的 `volchg` 是 1 索引，`bcf_text` 按 0 索引读，**BCF 电荷列整体错一位**。「bader weight 的真空电荷取错」一节只说了 vacchg，不完整 | `ferro bader -i tests/CHGCAR_2atoms -m weight`：BCF 体积 1 电荷 0、体积 2 为 53.000（实为体积 1 的）；`Vacuum charge` = 原子 2 电荷 52.99 | 建 `BaderResult` 前转 0 索引；或按布局建议 5 复用 grid 的函数 |

| # | 位置 | 问题 / 复现 | 修法 |
|---|---|---|---|
| B-6 **搁置**（2026-10-08 用户定：该功能后续会更新，届时一并修；候选思路：组内枚举、匈牙利 + Kabsch 迭代如 ArbAlign、距离指纹配对）| `cube_sdf.rs:502` 起 `heap_permutations` | 同签名族对**全部 P** 做 n! 枚举，每排列一次 SVD；docstring（`:414`）说「组内」，代码不分组。10 元环 1.2 s、11 元 10.7 s、13 元外推约 28 min，偏磷酸盐 Q2 长链常见 | 只在同标签组内置换并设上限，超限退 Hungarian 并告警 |

| # | 位置 | 问题 / 复现 | 修法 |
|---|---|---|---|
| C-D3 **搁置**（2026-10-08 用户定：随 Bader 重写一并修，同 C-D2）| `dft/bader.rs:174-178`；`cli-reference.md:800`、`dev/bader.md` §10 | ACF 的 `MinDist` 实为「极大值到原子距离」的最小值（原子 1 输出 0），不是 Henkelman 的「到 Bader 表面最小距离」；`bader.md` 称外部工具可按 Henkelman 格式解析，实测 ASE `attach_charges` 与 pymatgen `_parse_acf` 都失败（表头 `—` 非 ASCII、无 `----` 分隔、无 `VACUUM CHARGE:`） | 实现表面距离并恢复 Henkelman 版式；或改手册写明自有格式与列含义 |

- chg_sdf 跨文件平均的是 ρ·V_cell，NPT 快照下是体积加权（推断）
- `cluster.rs:56 NetworkGraph.former_qn` 是桥氧个数，`AtomType::Former.qn` 是同元素连接数，
- Bader 归属原子的 MIC 用分数取整、不查上界，强三斜小胞取错镜像（推断，极少触发）

5. **Bader weight 复用 grid**：`bader_weight.rs:236-262` 抄了 `assign_chg2atom` 与
   `calc_atomic_vol`，改传 `&volchg[1..=nvols]` 复用（R6），同时消掉 C-D2 那一族索引错误

### 2026-10-01 全库审查

- **cube_*（除 `calc_cube_radius`，已随审查低 6 改 `Result`）、`calc_chg_sdf`、`calc_cluster_sdf` 仍返回 `Option`**（M3 未动；`calc_network` 已随 L2 改 `Result`）：
  参数已在读文件前查过，剩下的 `None` 都是轨迹层面的（缺 cell、缺速度、没有团簇），
  但 CLI 只能猜原因（`missing cell, velocities, or forces?`）。改 `Result` 是同一个做法
