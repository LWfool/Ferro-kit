# 当前进度

## 测试总数：439 个（全部通过，clippy 零警告）

| Crate | 测试数 |
|---|---|
| ferro-core | 84 |
| ferro-io | 65 |
| ferro-structure | 70 |
| ferro-analysis | 163 |
| ferro-workflow | 23 |
| ferro-cli（lib 29 + 集成 5） | 34 |

## scripts/

两组脚本，共享层各一份：

| 组 | 共享层 | 脚本 | 用途 |
|---|---|---|---|
| 对拍 | `ferrocmp.py` | `compare_rdf/angle/sq/sq_experiment.py` | 跟 dump2analysis / dump2sq 逐点比 |
| 出图 | `ferroplot.py` | `plot_gr/angle/sq/net.py` | 发表级 pdf（+ png 看效果） |

两组不合并：都读 ferro 的 csv，但读完之后没有共同代码。

出图侧样式 `['science','vibrant']` + LaTeX + 四边框；每个脚本顶部一个 `CFG` 配置块。
`plot_net.py` 的 x 轴是**成分**（`file` 列），100 % 堆积柱，`--partner` 展开成
色相 = Qn、同色系明度 = m_<X> 的嵌套条带并给每个 Qn 组加框线。多 csv（`-o` 的 suffix
区分 CMD / MLMD）画成同一刻度下并排多根柱。

**对拍脚本当前是断的** —— `traj` 产物加了 label 段，三处写死的产物名要跟进，
见 `dev/plan.md`。

## 锚点 tag

**`v0.1.15`** = 批处理/长表重构前的最后状态（单文件输入、`.dat` 宽表、writer 在
ferro-analysis）。此后所有分析产物的文件名、扩展名、列结构全变。

---

## 各 Crate 完成状态

### ferro-core
- `Atom`、`Frame`、`Trajectory`、`Cell` 核心数据结构
- **`table.rs`**：`Table` + `enum Column { Num, Text }` —— 分析产物跨层的中立载体
  - `concat_union(label, parts)`：跨文件堆叠 + **列并集**，某输入缺的列留空（NaN），
    不补零、不插值 —— 这是「元素集不同的轨迹能同批跑」的支点
  - `to_comment_lines()`：把一张小表渲染成对齐文本，供 `#` 元信息块内嵌
  - `Column::cell()` 定义全项目输出约定：`{:.6e}`，**NaN 渲染为空字段**
  - 分层判据（已写进 CLAUDE.md）：类型下沉 core 当且仅当两个以上中间层要叫出它的名字。
    `Trajectory` 管结构方向，`Table` 管产物方向，两个方向由此对称
- `CubeData`：`data: Vec<f64>` + `shape: [usize; 3]`；`idx/get/get_mut` 访问器
- 静态数据：`elements.rs`（含 `symbol_to_z`、`group_number`、`valence_electrons`、`is_transition_metal`）、`compounds.rs`
- 单位转换：`LengthUnit`、`EnergyUnit`、`PressureUnit`、`TimeUnit`
- 错误类型：`ChemError` / `ferro_core::Result<T>`
- **`spin.rs`**：从结构推断未成对电子 / 自旋多重度
  - `guess_spin`：magmom 求和 → 氧化态+Hund 规则 → 电子数奇偶下限，三级回退
  - `assign_oxidation_states`：电负性规则 + 电荷守恒（离子型固体）
  - `total_electron_count` / `parity_min_multiplicity`；过渡金属取高自旋
  - 验证：ZnP₂O₆→0（抗磁）、MnS→5（d⁵ 高自旋）、Fe₂O₃→10
- `cell.rs`：`cartesian_to_fractional` / `wrap_position` / `minimum_image` 均返回 `Result`
  - **`interplanar_spacings()` / `minimum_image_cutoff()`**：dᵢ = 1/‖(Mᵀ)⁻¹.row(i)‖，
    最小镜像的正确几何上界（非正交晶胞 d < L，用边长会高估截断）；
    `calc_gr` 与 `CellList::build` 共用同一份
- `data/elements.rs` **`split_element_label`**：按第一个下划线拆分 `<Element>_<suffix>`，
  返回 `LabelSplit::{Plain, Split, Unknown}`；不用贪婪前缀匹配，避免 `Pb`→铅 的静默误判
- **`cluster.rs`**：网络连通性共享原语
  - `build_network_graph`：former–ligand 邻接图 + 配体分类 + 个人 Qn + 连通分量
  - `connected_components`：通用并查集（分量 ID 按首见根确定）
  - 供 `ferro-structure::find_clusters` 与 `ferro-analysis::cube_sdf` 复用，消除重复
- `Frame::unique_elements()`：按首现顺序去重元素（替代 5 处重复实现）
- **`network_type.rs`**：`AtomType`、`TypeParams`、`classify_frame[_detailed]`
  - **`enum AtomType`（0.2.1）**：`Former { elem, bridging, cn, bridges_to }` /
    `Ligand { elem, partners }` / `Modifier { elem, cn }` / `Other { elem }`
    —— 分类结果不再是字符串。此前 `classify_frame` 把算出的数字渲染成标签就丢弃,
    下游各自写解析器猜回来,全项目共五个这样的反解析点(`extract_qn`、三个
    `*_label_order`、`cmd/net.rs::type_sort_key`、`structure/cluster.rs` 的
    `starts_with("Ob_")`),改标签格式需同时改六处,漏一处即静默错数据
  - **`label()` 是全项目唯一把类型渲染成文本的地方**；`class_rank()` / `display_rank()`
    取代三个 `*_label_order` 与 `type_sort_key`；`is_bridging()` 取代前缀匹配
  - 标签(0.2.1 起,全部合 `<元素>_<后缀>` 约定)：
    形成子 `P_0`/`Al_2`（数字 = 桥接配体数）、`O_f`/`O_n`/`O_b`/`O_t`、修饰子裸元素符号
  - `classify_frame_detailed` 另返回 ligand→former 邻接（`FrameTypes`），供 linkage
    统计使用；配对统计需要图而不只是逐原子类型,重算一遍邻居搜索会把最贵的一步做两次

### ferro-io
- **`writers/table.rs`**：`write_table(&Table, path, TableFormat)` —— 分析产物的唯一出口
  - 取代 ferro-analysis 里原来的 7 个 `write_*`
  - ragged 表在**创建文件之前**拒绝，不留半个文件；字段含逗号才加引号
  - `TableFormat` 目前只有 `Csv`；加容器是加 match 分支，不是加调用路径
- 格式读写：XYZ、PDB、CIF、LAMMPS dump 等
- **`lammps_dump.rs`**：`element` 列写位点标签（`P_0`/`O_b`/`Al_2`）时拆成
  element + label，读取结束打印一次映射表；前缀非法元素时整串当元素并告警。
  普通符号原样通过、`label = None`
  - writer 的整数 `type` 列**在整条轨迹上确定一次**（0.2.1）。原先在帧循环内按
    「该帧首次出现的名字顺序」重建：元素只有几种时碰巧稳定，写位点标签后某帧
    恰好没有 `P_4`，同一个 type 编号就会在不同帧指向不同的东西
- **`extxyz.rs`**：`label:S:1` 列（0.2.1，读写两侧）。extxyz 的列是自描述的，
  故标签走自己一列、`species` 保持纯元素 —— 与 dump「没地方放第二个名字只能折进
  element 列」相反，两个方向都无损
- **`cube.rs`（重构）**：
  - `read_cube(path) -> Result<CubeData>`：可视化/密度图（原有功能保留）
  - **NEW** `read_cube_as_chg(path) -> Result<(Frame, ChargeGrid)>`：供 Bader 分析用
    - 坐标：Bohr → Å，减去非零原点偏移
    - 索引转置：cube 存储（i3 最快） → ChargeGrid（n1 最快）
    - 密度缩放：`rho_stored = ρ_cube × V_cell_Bohr`
  - 共用 `parse_header()` 中间结构体，避免重复解析逻辑

### ferro-structure
- `supercell.rs`、`vacuum.rs`、`merge.rs`
- `box_builder.rs`：随机放置 + soft-core 弛豫，O(N) cell-list 加速
  - **`parse_formula` 支持嵌套括号**（0.1.14）：栈式解析，`()` / `[]` 均可且需正确配对；
    `Ca3(PO4)2`、`(NH4)2SO4`、`K4[Fe(CN)6]`、`Mg[Al(OH)4]2`
  - **`resolve_component`**：先查 COMPOUNDS（保留库中人工录入的 molecular_mass），
    查不到则把输入当化学式解析、按原子量求和 → 可直接建库外化合物的盒子
  - 拒绝空式、空组、零计数、括号错配、未知元素（原先空式返回 Ok(空表)，把报错推迟到
    "no atoms to place"）
  - 私有 `CellList::neighbors` 加 `seen` bin 去重：修复小盒子（n_bins ≤ 2）下近邻被重复计入、弛豫斥力成倍累加的问题
- **`typing.rs`**：`classify_trajectory`、`apply_type_labels`（单帧类型标注 + 结构导出）
- **`cluster.rs`**：`find_clusters`（多形成子，复用 `ferro_core::connected_components`）

### ferro-analysis / md

**0.1.15 之后：分析层不再碰文件系统。**7 个 `write_*` 全部退役，各结果类型改为
`to_tables() -> Vec<(String, Table)>` + `meta_lines() -> Vec<String>`。
`meta_lines` 只放**批内共享**的参数；逐文件才有意义的量（组成、clamp 后的 r_max）
走 CLI 的 `[inputs]` 清单，不能混进共享参数块冒充全局事实。

| 模式 | 表结构 | 列 |
|---|---|---|
| gr | 长表 | `file, r, center, neighbor, gr, cn` |
| angle | 长表 | `file, angle, end_a, center, end_c, count, p` |
| sq | 宽表 | `file, q, total_xrd, total_neutron, <pair>_sq/_xrd/_neutron` |
| msd | 直加 file 列 | `file, time, msd, msd_a, msd_b, msd_c` |
| vacf | 直加 file 列 | `file, time, vacf, vacf_x/y/z, diffusion` |
| rotcorr | 直加 file 列 | `file, time, c2, integral` |
| vanhove | 直加 file 列 | `file, r, gs` |

判据：**表结构跟主产物的粒度走**。gr/angle 的每一对/三元组都是独立主产物，且元素集
逐文件可能不同 → 长表（类型进数据列，堆叠无需对列）；sq 的主产物是两条 total
（一行一个 q），partial 只是诊断附属 → 宽表。
angle 保留原始 `count` 与归一化 `p` 两列：整数直方图是 `compare_angle.py` 与
dump2analysis 逐 bin 对拍的依据，不能只留 `p`。

- `gr.rs`：g(r) + CN(r)，CellList 加速；**n² 个有序对**输出
  - **0.1.15**：`--r-min` 上 CLI（`GrParams.r_min` 早已存在，只是没暴露）
  - **逐帧体积归一化**（0.1.13）：`g = ⟨hist·V_f⟩`，非 `⟨hist⟩·⟨V⟩`；fold 维护
    `plain`（Σcount）+ `vol_w`（Σcount·V_f）两份直方图，与 code1/code2 的 `cr` 临时数组同构
  - 新字段 `rho_g = ⟨ρ_f·g_f⟩`（供 S(q) 重构逐帧变换）、`volume_std`（两遍算法）；
    `rho` 重定义为 `N·⟨1/V⟩`（原 `N/⟨V⟩`）
  - 粒子数守恒校验：逐帧比对总数与分组计数，不符即 `Err`
  - `gr` 对称（`A-B` 与 `B-A` 逐点相同），`cn` 有向（`CN(A→B)=Σ ni·hist/(N_A·steps)`）
  - `GroupBy::{Element, Label}`：Label 模式按 `Atom::label` 分组（未设则回退 element）
  - 排序 `(elem_z, 字符串)` 二级比较 → 同 Z 标签列顺序可复现
  - **r_max clamp** 用最小**面间距**的一半（`Cell::minimum_image_cutoff`），非最短边长
  - `write_gr(result, path, pair)` 单函数：`Some` 出三列、`None` 出宽表；`write_cn` 已删
  - 已删：`PairStats` / `pair_stats` / `WelfordStats` / `GrParams.r_cut` / total
- `sq.rs`：S(q)，只对规范半边做傅里叶变换
  - **0.1.15**：`--q-min` 上 CLI（同上）
  - **逐帧口径**（0.1.13）：被积函数改 `(rho_g − rho)`、系数 `4πρ/q → 4π/q`；
    签名与分层不变，结果严格等于 code2 的逐帧变换（对拍 <1e-10）
  - `SqResult`：`sq`（未加权 partial）/ `sq_xrd` / `sq_neutron` / `total_xrd` / `total_neutron`
  - 加权 partial `w_ij(q)·S_ij(q)` 逐点相加恰等于对应 total
  - `to_tables(gr)` **恒写全部规范半边**（0.2.1 后）：pair 参数已删,CLI 侧的
    `-a/-b/-x/-y` 一并移除 —— 只留一对会藏起 `Σ w_ij·S_ij = total` 这条闭合
  - 散射因子查表失败时告警（此前静默输出垃圾值）
  - Faber-Ziman 权重改上三角遍历，修复异种对被重复累加（实测 total_xrd 曾偏大 1.59 倍）
- `msd.rs`：MSD，NPT 安全，atom-major 并行；`fit_diffusion` 线性最小二乘拟合自扩散系数 D=slope/6（Einstein 3D）+ R²；`plot_msd` 绘图（total + a/b/c + 拟合线）
- `angle.rs`：键角分布；支持 `GroupBy`，类型排序改 `(Z, 字符串)`
  - **0.1.15**：`--r-cut-ab` / `--r-cut-bc` 改按调用方写的 `-a`/`-c` 顺序分配
    （新增 `AngleParams::ends`），此前按原子序数派发、静默忽略用户顺序；
    未点名三元组时仍用规范 (Z, 符号) 顺序，两端同类型时取 min
  - **0.1.15**：新增 `--angle-min` / `--angle-max`，直方图上下界可调；上界闭区间
    （`acos(-1)` 恰为 180.0）。加 `--angle-min 0.05 --angle-max 180.05` 可逐点复现
    dump2analysis 的分箱（补 ×2 计数后 1800 个 bin 整数零差）
- `vanhove.rs`：Van Hove 自关联
- `vacf.rs`：速度自相关
- `rotcorr.rs`：转动相关 C₂(t)
- `cube_density.rs`：3D 密度/速度/力分布
- `cube_sdf.rs`：团簇 SDF（Kabsch 对齐）；团簇查找改用 `ferro_core::build_network_graph`，已去除 petgraph 依赖
- `chg_sdf.rs`：电荷密度团簇 SDF（Kabsch 对齐 + pull 插值旋转子格）

### ferro-analysis / network（已重构）
- 单文件 `mod.rs`，依赖 `ferro_core::classify_frame_detailed`
- 参数类型：`TypeParams`（来自 ferro-core，`cutoffs` + `modifier_cutoffs`）
- 参数类型另带 `qn_elements`（默认 `{B,P,Si}`，来自 `ferro_core::data::qn_elements`），
  由 `--qn` 整体替换。约定随参数走而不是渲染时查表 —— 一个类型必须携带它是在哪套
  约定下产生的，否则同一批次两个文件里 `Al_4` 含义不同
- 结果：`NetworkResult`（`qn_dist` / `qn_partner_dist` / `mean_qn` / `cn_dist` /
  `mean_cn` / `oxy_dist` / `ligand_dist` / `linkage`），`to_tables()` 出六张长表
  （`composition` / `qn` / `qn_partner` / `ligand_type` / `coordination` / `linkage`）
- **`Bin { count, fraction, sd }`**：`fraction` 是**逐帧比例的平均**，`sd` 是同一序列的
  样本标准差（ddof=1），缺席帧按 f=0 计入。`sd` **不是标准误** —— MD 相邻帧强相关，
  既不按 1/√N 收缩也不估计物理涨落
  - 用 Welford 而非 `Σf`/`Σf²`：两遍求和对恒定序列给出 ~2e-10 的抵消噪声
    （实测 P 的 Q0 档逐帧恒为 5 个），而 sd 列里的假抖动会被当成物理。
    rayon 是归并而非顺序折叠，故用 Chan/Golub/LeVeque 成对合并公式
- **Qn 只报给 Qn 形成子**：Al 从 `qn` / `qn_partner` 的**行**里整体退出，仍保留在
  `m_Al` 伙伴列、`ligand_type`、`coordination`、`linkage`。Qn 是四面体形成子的记号，
  Al 的表征量是配位数（Al[4]/Al[5]/Al[6]），给它报 Qn 等于发明一个文献里没有的量。
  非 Qn 形成子的 **label 数字改取配位数**（`Al_4` = 四配位），`AtomType::Former.qn:
  Option<u32>` 既是标志也是取值
  - **测试盲区**：参考轨迹的 Al 不带非桥氧（`ligand_type` 无 `O_n, Al` 行），故逐
    原子 `bridging == cn`，两种口径的分叉从未被触发。已补
    `test_non_qn_former_labels_by_coordination_not_bridging` 造出分叉构型
  - 两张 Qn 表在没有 Qn 形成子时**整个不写**并打印原因：只有表头的 CSV 读起来像
    「测了，结果是零」
- **`linkage` 表**：一行一种「配体元素 × 两端位点状态」组合，两端**各带元素、桥接数、
  配位数**三个字段，外加人可读的 `linkage` 展示列（`Al_4-O-P_2`，经
  `AtomType::label` 渲染，与导出轨迹同一套词汇）。参考实现（`private/coord_analysis.py`）按元素只存一个数（P 存 Qn、
  Al 存 CN，由 `QnEle` 列表决定），故「4 配位 Al 的桥接数」问不出来
  - 规范半边：桥联无方向，两端按 `(元素, 桥接数, 配位数)` 排序后小的在前，
    每对只存一次（与 sq 只做规范半边同理）。行和不等于该位点的总参与度
  - `n_formers` 列：普通桥氧为 2，三配位氧展开成 C(3,2)=3 行并标 3
  - **`ligand` 进键**（0.2.1 后补）：此前双配体体系（`--Al-O` + `--Al-F`）会把 Al-O-P
    与 Al-F-P 静默合并成一行，报出的计数无法与实验对应
  - `n_bridge_a/b` 保留此名而非 `qn_a/b`：桥接数对每个形成子都有定义，Qn 只对一部分
    成立。这是唯一还需要给「非 Qn 形成子的桥接数」命名的地方
- **`qn` 与 `qn_partner` 是两张表**（0.2.1 后修）：前者是 Qn 的边际分布（打开文件即可
  读），后者多一维伙伴分解。0.2.1 曾合成一张，理由是「`groupby` 能退回去」——结果是
  打开文件**看不到 Qn 占比**，它散在 14 行里。粒度不同就该分表，与 `average` 独立成
  表同理；`sd` 必须各自累加，相关项之和的方差不等于方差之和（实测 P 的 `qn=2`：
  正确 5.574e-3，partner 各行相加 9.966e-3）
- **`composition` 表**（0.2.1 后加，取代 `average`）：一物种一行，分母恒为**该元素**
  的原子数（Q2 占全部 P、`Al_4` 占全部 Al、`O_b` 占全部 O），故「每元素 fraction 求和
  为 1」是可核对的恒等式。每个元素只出现一种刻画，P 不再按配位数重复一遍。配体行由
  `ligand_type` 的伙伴维聚合而来，`sd` **逐帧重算**不是把源行相加
  - 原 `average` 的 `mean_qn`/`mean_cn` 未丢：在每个文件的 `[inputs]` 块里，且可从
    本表精确复算（`Σn·f_n` 实测 2.395161 / 4.116667）
- **两套标签词汇**：分布表（`composition`/`qn`/`qn_partner`）用**单元**词汇 `P-Q2`，
  `linkage` 与导出轨迹用**原子**词汇 `P_2`/`Al_4`。桥联连的是原子而 Qn 单元含多个
  原子；轨迹标签还必须能被 dump reader 拆回 element，`Q2` 里 `Q` 不是元素。
  `AtomType::label` 是原子形式，`species_label`（ferro-analysis）是单元形式并在
  非 Qn 情形委托给前者
  - 带元素前缀是因为双 Qn 形成子体系（B+P、Si+Al）里 `Q0`…`Q4` 会重名，而**跨元素
    重名没有任何相邻列能解**
  - **伙伴分解不编码进 label**：文献的 `Q^n(mX)` 表达不了「该桥的配体连三个形成子」，
    参考数据里 `(qn=2,m_Al=1,m_P=0)` 与 `(qn=2,m_Al=1,m_P=1)` 都会写成 `Q2(1Al)`
- **`ligand_type` 的 `label` 合并**为 `Al-O_b-P` / `P-O_n` / `O_f` / `O_t`（缺位不补
  占位符），原标签移入 `type` 列，`former_a`/`former_b` 保留
- **配体 `fraction` 的分母改为「该配体元素」**（0.2.1 后修）：文献 BO fraction 的口径。
  此前取全体配体原子，双配体体系里 `O_b` 会被 F 稀释。单配体体系数值不变
- **每张表自带 `#` 头**（表意 + 逐列 + 该表的坑）。为此 `batch::write_all` 由覆盖
  `Table::meta` 改为拼接，`Table::concat_union` 同步保留各部分 meta（取第一份），
  否则说明在 `stack` 阶段就丢了。帮助文本删掉的内容落在这里
- **`label` 列与数值列并存**，不是二选一：label 给人读，`former`/`qn`/`cn` 给筛选和
  画图。删掉数值列会让「筛出 Qn ≥ 3」退化成字符串切分 —— 正是这次重构从五处调用点
  消灭掉的反向解析
- **`m_<X>` 分解**（`bridges_to`）：Q^n(mAl) 记号。只对「恰好两个形成子」的配体成立，
  故 `Σm ≤ qn`，差额即三配位桥数。**无法从 linkage 反推** —— 两座 P-O-Al 桥
  可能来自一个 m_Al=2 的 P，也可能来自两个 m_Al=1 的 P；linkage 数桥，这里数原子
- 修饰子角色分类（`Zn_f`/`Zn_t`/`Zn_b`/`X`）已删除：实测参考轨迹 97.1% 落进 ≥3
  兜底桶（0 / 1 / 26 / 903），无分辨力。改用配位数描述，Zn 进 CN 表
- 旧的 `cn.rs`、`ligand_class.rs`、`qn.rs`、`modifier.rs` 已删除（逻辑迁移至 ferro-core）

### ferro-analysis / dft
- `bader.rs`：`BaderAnalyzer` builder 模式、`BaderResult`、ACF/BCF/AVF 输出
  - `write_acf` 接受 `frame: &Frame` 参数，写入真实原子坐标（修复占位符 0,0,0）
- `bader_grid.rs`：on-grid / near-grid / off-grid 三种梯度上升方法，含边缘精化
  - `bader_ongrid` 删除被立即丢弃的第一遍梯度上升循环（遗留开发产物），运行时减半，行为不变
- `bader_weight.rs`：Yu-Trinkle weight 方法（WS Voronoi + 流分配）
  - 修复边界点电荷丢失（改用 O(N_vols × N_pts) Fortran 算法）
  - 修复 volnum 边界点指向错误（用负数编码最大概率上游体积）
- `charge_grid.rs`（ferro-core）：`ChargeGrid` 数据结构
  - 修复 `lat2car` 构造错误（列 = a_i/N_i，非正交晶胞关键）；新增非正交测试
- `chgcar.rs`（ferro-io）：CHGCAR 文件读取器
- `cube.rs`（ferro-io）：Gaussian cube 文件读取器（支持 QE pp.x 输出）
- `ferro bader`（ferro-cli）：自动识别扩展名（`.cube` → `read_cube_as_chg`，其他 → `read_chgcar`）

### ferro-workflow
- `GaussianJobBuilder`：Gaussian 输入文件生成
- `Cp2kJobBuilder`：CP2K 输入文件生成
  - 支持任务：energy、force、geo-opt、cell-opt、md、freq
  - 支持泛函：PBE/BLYP/PBE0/B3LYP/RevPBE/PBEsol/SCAN/r2SCAN/HSE06 及自定义
  - 支持色散：DFT-D3、D3(BJ)；SCF 求解器：对角化、OT；k 点、涂抹、原子电荷、cube 输出、Molden 输出
  - MD 参数：CSVR/NoseHoover/Langevin/NVE 热浴，NPT 压浴
  - **`auto_spin`**：经 `ferro_core::guess_spin` 推断多重度 + UKS（CLI 显式 `--multiplicity` 时关闭）
  - **`cp2k_basis_db`**：从 6 个 CP2K 文件解析的基组/赝势数据库（2829 条）
    - 覆盖 PBE / SCAN / 全电子（pob、-ae），全元素全等级
    - `basis()` / `potential()` 按元素 + 族前缀 + 泛函精确匹配，q 价电子数一致
    - 新增全电子基组选项 `pob-dzvp`、`pob-tzvp`；启发式 `gth_nval` 降为兜底
- `QeJobBuilder`：Quantum ESPRESSO pw.x 输入文件生成
  - 任务：scf/nscf/bands/relax/vc-relax/md/vc-md（&CONTROL/&SYSTEM/&ELECTRONS/&IONS/&CELL）
  - 泛函（input_dft）：PBE/PBEsol/RevPBE/BLYP/SCAN/r2SCAN/PBE0/HSE06 及自定义
  - 平面波截断、k 点（Gamma / Monkhorst-Pack）、展宽（gaussian/mp/mv/fd）
  - 自旋复用 `guess_spin`（→ nspin/tot_magnetization），与 CP2K 路径一致

### ferro-cli

**单二进制 `ferro`（0.2.0）。**八个 `fe-*` 已删除，子命令按**产物**分组：

```
ferro traj  gr | sq | msd | angle | vacf | rotcorr | vanhove   → 堆叠 csv + 可选 PNG
ferro map   density | velocity | force | radius | sdf | chg-sdf → 逐输入一个 .cube
ferro net                                                      → 六张堆叠 csv
                                                                 + 可选标注轨迹
ferro bader | convert | info | job
```

- `fe-traj` / `fe-corr` 合并为 `traj`：七个分析共用同一条导出管线，而原来的
  「结构 vs 动力学」界限本就不成立（`msd` 是时间关联量却在结构一侧）
- `-m` 全面消失，位置即模式。**每个子命令自带参数结构体**——`--dt` 在 `msd` 与 `vacf`
  下语义不同，不再被迫共用一个字段（旧 `fe-traj` 是 29 个字段的大结构体）
- `CommonArgs` / `SelectArgs` 经 `#[command(flatten)]` 接入，消掉四份重复声明
  （`CommonArgs` 此前定义了却从未被任何 bin 使用）
- 三级帮助：`ferro` → 分类总览；`ferro traj` → 该组命令；`ferro traj gr`（无 `-i`）→
  参数、输出列结构与示例
- 裸 `ferro` 打印总览。REPL 落地后再改为「tty 进 REPL、管道读 stdin」
- `net` 的 `--P-O=2.3` 由 `main` 在 clap 解析前从 argv 剥离（元素对写在参数名里，
  clap 建模不了）
- **不留兼容层**：输出格式同期变更，留着 `fe-traj` 会让旧脚本「跑成功」却吐出自己
  解析不了的 csv——静默坏数据比命令消失难查

- **产物命名（0.2.1 后）**：`<outdir>/<mode>[_<表>][_<label>]_<后缀>.csv`
  - `Output { dir, label, suffix }` 收口三个字段（否则 `write_all` 要吃 8 个位置参数），
    配 `out_path` / `file_label`（原样拼）/ `set_label`（排序去重）
  - `label` 由类型选择填出：`gr_P-O.csv`、`angle_O-P-O.csv`、`msd_O-P.csv`、
    `rotcorr_O-H.csv`；未点名时为 `all`。**排在 suffix 之前** —— `ls gr_P-O_*`
    列出同一对在各批次的结果
  - **两条拼法，各有理由**：`gr`/`angle` 按写的顺序（`CN` 有向，`-a P -b O` 与
    `-a O -b P` 本就是两份数据）；`--elements` 排序去重（是集合，同一份数据不能
    落到两个文件名下）
  - 选中的串会成为**路径的一段**，故拼名前校验 `[A-Za-z0-9_+-]`，在读第一个文件
    之前失败。替换非法字符被否掉：那会让 `-a P/2` 与 `-a P_2` 静默写进同一个文件
  - `--outdir` 进 `CommonArgs`（11 个命令），`Output::prepare` 缺失即创建并打印，
    同样在读第一个文件之前。覆盖 csv / png / `map` 的 cube / `net --export-traj`
    的轨迹；`chg-sdf` 不走 `CommonArgs`，单独加了同名参数
  - PNG 由数据文件**路径**派生（`png_path` 改收 `&Path`），目录与 label 自动跟随
- **`traj sq` 没有类型选择**（0.2.1 后）：`SqCmd` 去掉 `SelectArgs`，
  `SqResult::to_tables` 去掉 pair 参数。主产物是两条 total，partial 的价值全在
  `Σ w_ij·S_ij = total` 这条闭合上，只留一对恰好把它藏起来。按 label 分辨的 partial
  随之消失（位点标签对应的原子数往往不足以显出信号）；**库层 `GroupBy::Label` 不动**
- **`batch.rs`**：多输入驱动，**对结果类型泛型**，不认识任何分析类型
  - `expand_inputs`：自展开 glob（`glob` crate，不支持 `{a,b}`），参数顺序 + 参数内
    字典序，跨模式按规范化路径去重，**零匹配报错**
  - `map_inputs<T>(inputs, f)`：串行遍历（轨迹是内存大头，逐条释放；帧内并行不变），
    失败跳过并收集
  - `stack<T>`：按表名分组，各组 `Table::concat_union("file", ...)`
  - `write_all` / `out_path` / `Output` / `meta_block` / `Summary`
  - `Summary` 存**预格式化文本**：`[inputs]` 是给人看的，`{:.6e}` 会把「5 帧」写成
    `5.000000e0`
- **单一代码路径**：N=1 是 N 的特例。按文件数分派会让产物形态取决于 glob 当天匹配到
  几个文件
- 失败跳过 + `[inputs]` 留原因 + **退出码非零**（否则 shell 里 `&&` 串联会把批内失败
  当成功）
- `-o` 是**文件名后缀**：`<mode>[_<子表>]_<后缀>.csv`
- **`plot.rs` 重写为面板模型**：`Panel`/`Series` + 通用 `render`，一张 PNG 分格，
  一格一个量、一条曲线一个文件，颜色按文件跨格一致，图例只画第一格。
  gr 的双 Y 轴已删（g(r) 与 CN 各占一格）
- **绘图分辨率 500 dpi**（0.2.1）：一格 2708×2083 px，约可清晰放大到 5×。
  版式按 96 dpi 像素编写、统一过 `px()` 缩放，改 `DPI` 只整体放缩不重排——
  但**每个长度都必须过 `px()`**（含字号与 `legend_area_size`），漏掉就留在未缩放的
  默认值上。矢量 PDF 方案已验证可用但**因依赖成本回退**，见 `dev/issues.md`
- **类型选择**（`traj` 的 gr/sq/angle）：`-a/-b/-c` 按 element、`-x/-y/-z` 按 label，
  两组互斥；`-a`/`-x` 为中心、`-b`/`-y` 为近邻，顺序影响 CN
- `plot_gr` 单配对（修复原先只画低 Z 作中心那一向 CN 的 bug）
- **`ferro net`（0.2.1 起为叶子命令）**：`net qn` / `net type` 合并。两者从来不是
  两个分析——都对每帧每个原子做同一次分类，`type` 只是把结果写出去而不是汇总。
  导出退化为开关 `--export-traj [lammpstrj|extxyz]`，统计恒定执行
  - 接 `CommonArgs` + `batch::map_inputs`：`-i` 多值 + glob，五张长表带 `file` 列
  - `--format`（xlsx）与 `--frame` 删除；`-o` 由「路径基名」改为全项目统一的「文件名后缀」
  - `--modifier` 点名却没给对应截断时直接报错（静默通过会让该元素被当成形成子，
    氧的分类整体错位）

### ferro-python（PyO3 绑定，已完成）
- 独立 workspace（不进主 workspace，pyo3 0.21 extension-module），用 maturin 构建
- `lib.rs`（`#[pymodule]` + `pyerr` 错误转换）、`types.rs`（`PyTrajectory`）
- `io.rs`：`read`（11 格式按扩展名分派）/ `write`，metal_units 开关
- `structure.rs`：`supercell` / `add_vacuum_layer` / `merge`
- `analysis.rs`：`gr_pair(traj, a, b, by=...)` / `gr_all(traj, by=...)` / `msd`
  → `dict[str, list[float]]`（PyO3 自动转）；已删 `r_cut`、无 total
- 验证：`cargo check`/`clippy` 零警告；`maturin build` link 通过；Python 冒烟测试全过
  （读 xyz/cif/lammpstrj、supercell 568→2272、write、gr、msd 含元素过滤）
- 注：本 crate 无 `cargo test`（cdylib 绑定层），经 maturin + Python 验证

- 版本号：**0.2.1**（workspace 统一；ferro-python 已同步并复核编译通过）。
  `v0.2.1` 之后另有一批**尚未发版**的 net 破坏性改动，清单见 `overview.md`

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
| `rand` | 0.10 | workspace | 0.8→0.10 改名三处，见 issues.md |
| `clap` | 4.5 | ferro-cli | derive |
| `plotters` | 0.3 | ferro-cli | `default-features = false` + 必须保留 `ttf`；backend 为 `bitmap` |
| `pyo3` | 0.29 | ferro-python | 0.21→0.29 仅需 `skip_from_py_object`；**运行时未验证** |

`cargo update` 在主 workspace 与 ferro-python 均已无可更新项。

---

## 已知限制

- `ferro-python` **能编译**（2026-08-12 复核：`cargo clean && cargo check` 干净通过，
  clippy 零警告，版本已同步至 0.2.1）。此前记录里的「编译断裂」经查不成立。
  真问题是它作为独立 workspace 被主 workspace 的 `cargo build/test/clippy` 全部跳过，
  断裂不会被自动发现；改公共 API 后须手动补跑。
  待办仍是 **pyo3 0.29 的运行时验证**（本机无 maturin），**优先级中**

- `box_builder`：不支持水合物点记法（`CuSO4·5H2O`）—— 水应作为独立 component 传入
- `box_builder` 无 CLI / Python 入口，目前只能作为库函数调用
- `cp2k_basis_db`：源数据为 gitignore 的 examples/ 6 文件，DB 已固化为静态表（重生成需源文件）
- QE 赝势仅占位 `<El>.UPF`（QE 用 UPF，与 CP2K GTH 数据库不通用，由用户提供 pseudo_dir）
- **REPL / 脚本模式未实现**。`main.rs` 是可用的子命令分发器（不是占位符），
  裸 `ferro` 打印分类总览；REPL 落地后改为「tty 进 REPL、管道读 stdin」
- **pyo3 0.29 升级只做了类型层验证**（`cargo check`/`clippy` 零警告）。本机无 maturin，
  `cargo build` 在 macOS link 阶段过不了，模块初始化与 `#[pyfunction]` 签名的运行时
  行为（含新拆的 `gr_pair`/`gr_all`）待 `maturin build` + Python 冒烟测试确认
- 代码库未纳入 rustfmt 管理，格式为手工维护（含刻意的列对齐）；是否采用见 dev/plan.md
- `cube_density.rs:188` 用参考帧体积归一化，NPT 下同类偏差；需先定义「NPT 下 3D 密度图
  指什么」再动（见 `dev/issues.md`「NPT 逐帧体积归一化」的「明确不做」）
- **按标签选型只在单帧成立**：`traj gr -x P_3 -y O_b` 对多帧标注轨迹会被
  `calc_gr` 的「逐类型粒子数守恒」守卫拒绝——标签 population 逐帧在变
  （实测 P_3：149/152/150/150/150）。这是 `dev/plan.md` 原先写的动机
  「导出即可直接串起来」比实际窄的地方。多帧要么按元素选（`-a P -b O`，正常），
  要么等 g(r) 支持逐帧归一化的类型计数（**未列入计划**，需先定义口径）
- `ferro map chg-sdf` 的 `--cubes` 仍是多 cube 聚合成一张 SDF，与 `map` 其余模式
  「一输入一产物」相反；拆分需先定义带样本计数的中间产物格式。
  **优先级高**（2026-08-11）
