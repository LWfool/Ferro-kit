# 已知问题 / 编码陷阱

## 编码陷阱

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| `cell.rs` | 浮点误差约 1e-15 | 断言用 `< 1e-10`，不能用 `assert_eq!` |
| `cube_density.rs` | 原子坐标放格点边界会导致归属歧义 | 测试中用格点中心 `(n + 0.5) / N * L` |
| `cell.rs` | `wrap_position` 负数行为 | 用 `rem_euclid(1.0)`，不要用 `x - x.floor()` |
| `ferro-cli` | plotters `ttf` feature | `default-features = false` 会禁掉字体支持，必须保留 |
| `ferro-cli` | plotters 借用 | `BitMapBackend::new(&path, ...)` 借用 `path`，绘图块用 `{}` 先析构再使用 `path` |
| 全局 | clippy snake_case | 变量名不能用单大写字母（如 `L`），改为 `box_len` 等 |

| `cube_density/radius/jump/sdf.rs` | ndarray `Array3` 仍在内部使用（分析阶段） | 在 `CubeData` 构造点用 `.into_iter().collect()` 转为 `Vec<f64>`；不要对外暴露 Array3 |

---

## Bader 实现陷阱（来自 Fortran 源码逆向分析）

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| CHGCAR 读取 | VASP 存储的是 `ρ × V_cell`，不是物理密度 | 读入不做归一化；真空判断除以 `vol`；电荷积分除以 `nrho` |
| `ChargeGrid` | ndarray 不可对外暴露 | 存储用 `Vec<f64>`，内部计算可用局部 ndarray |
| `step_neargrid` | Fortran `dr` 是 `SAVE` 变量，跨调用持续累积 | 放入 `PathState.dr`，新路径（pnum==1）时归零 |
| `max_neargrid` | `known=1` 与 `known=2` 处理完全不同 | known=1：step 内 fallback to ongrid；known=2：max 层 EXIT 继承 volnum |
| weight 极大值位置 | Fortran 原始代码存的是循环末尾的邻居 `p`，不是极大值本身 | Rust 改为 `chgList[n].pos`（bug 修正） |
| `rho_grad_dir` | 对高倾斜晶胞梯度方向有近似误差（car2lat 用了两次） | 正交/弱倾斜体系无影响；强倾斜晶胞改用 weight 方法 |
| CHGCAR 密度循环顺序 | 文件中 n1(x) 最快，n3(z) 最慢 | Rust Vec 索引：`rho[n1 + N1*(n2 + N2*n3)]` |
| `lat_i_dist(0,0,0)` | 中心点距离为 0，不能取倒数 | 显式存为 `0.0`（step_ongrid 跳过 d=0 方向） |

---

## network_type / typing 编码陷阱

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| `network_type.rs` `build_nf_map` | 若形成子元素同时也是配体元素（如 Si-O 且 O-O），同一原子会同时出现在 former 和 ligand 中 | 当前跳过 `fa_idx == la_idx`；实际体系不会出现，无需处理 |
| `classify_frame` 形成子标签 | 桥接数来自桥氧计数，桥氧定义为 NF 邻居数 ≥ 2（含三配位氧 `O_t`） | 0.2.1 起由 `AtomType::Former.bridging` 直接携带；**任何从标签字符串回推数字的写法都不要写**，这正是 0.2.1 重构消除的五处缺陷 |
| `ferro-cli` 解构 `cli` | `let Some(x) = cli.field` 是部分移动，后续不能再借用 `cli` | 用 `let Cli { field1, field2, .. } = cli;` 一次性解构，再分别传参 |
| `cluster.rs` Union-Find | `find()` 需要 `&mut Vec`，不能与 `local_cluster` 构建同时进行 | 先用 `(0..n).map(|li| find(...))` 收集，再组装 `Vec` |

## spin / cp2k_basis_db / qe 编码陷阱

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| `spin.rs` `assign_oxidation_states` | 单质（distinct 元素 < 2）氧化态法不适用 | 返回 `None`，由 `guess_spin` 回退奇偶下限 |
| `spin.rs` 过渡金属 | 一律高自旋；低自旋 / 多中心磁耦合无法从结构推断 | 结果为初猜须 DFT 验证；附带 warnings |
| `cp2k.rs` `auto_spin` | builder 默认 `auto_spin=true`，会覆盖 CLI 显式多重度 | job.rs 设 `auto_spin = args.multiplicity.is_none()` |
| `cp2k_basis_db` `name_rank` | 同 q 多候选（`GTH-PBE` vs `GTH-NLCC-PBE`） | 精确等于 `GTH-PBE`/`GTH-SCAN` 优先（rank 0） |
| `cp2k_basis_db` | 文件自动生成，标注「Do not edit by hand」 | 改逻辑须改生成脚本并重跑（数据源于 6 文件） |
| `cp2k.rs` `gth_nval` | 兜底价电子数不能用 DB 内 `max(q)`（Zn 有 q2/q12/q20） | 复用 `cp2k_basis_db::basis()` 默认族解析，取惯用 q |

## cluster / 审计 编码陷阱

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| `ferro-structure::find_clusters` | 用 `AtomType::is_bridging()` 当连通性判据 | 它是「**恰好**两个形成子」，三簇配体返回 `false` → 网络在连得最紧的节点处断开（实测 3 个 P 共享一个三簇氧报成 3 个团簇，应为 1 个）。连通性判据是 `partners.len() >= 2`，与 `LigandKind::Bridging`、与模块 doc 一致。`ferro-core::build_network_graph` 那份从一开始就是对的，两份实现不一致了很久 |
| `find_clusters` 的截断 | 按形成子元素取 `cutoffs.keys().find()` 的第一个 | 按 `(形成子, **该配体元素**)` 取。`BTreeMap` 序下 `("P","F") < ("P","O")`，双配体体系会拿 P-F 的 1.4 Å 去判 P-O 的 1.6 Å 键，桥氧整体判不出来。与 `LinkKey` 少一个配体元素维是同一类洞：参考体系只有 O，所以一直没暴露 |
| `ferro-core::cluster` | ferro-analysis / ferro-structure 中间层不能互依赖 | 共享原语下沉 `ferro-core`（两者共同依赖）；ferro-analysis 已去 petgraph |
| `cif.rs` 坐标回退 | 依赖远处 `ensure!` 的 `.unwrap()` 脆弱 | 本地 `else if let (Some,Some,Some)` 守卫 + `bail!`，不变量自洽 |
| `vasp.rs`/`chgcar.rs` VASP4 检测 | 守卫 `parse::<u64>` 而 unwrap `parse::<usize>`，32 位平台大值溢出 panic | 单次按 `usize` 解析为 `Option<Vec<usize>>`：`Some`→VASP4，`None`→VASP5 |
| `extract_element` (cp2k/qe reader) | `it.next().unwrap()` | 安全：上方 2 行 `is_empty()` 局部守卫（局部不变量可接受，区别于 cif 的远处不变量）|

## ferro-python 编码陷阱

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| ferro-python 构建 | `cargo build` 在 macOS link 阶段报 undefined Python 符号 | 用 maturin（提供 `-undefined dynamic_lookup`）；`cargo check` 可类型校验 |
| `Cargo.toml` | 不进主 workspace，但 `version.workspace=true` 需成员资格 | 显式 `version`/`edition` + 空 `[workspace]` 表使其独立 |
| maturin 解释器 | 默认可能选到杂散 python3.9，wheel tag 不匹配当前 python | `maturin build --interpreter "$(which python)"` |
| `pyproject.toml` | `readme="README.md"` 缺文件 → maturin 失败；失效 `python-source` | 必须有 README.md；无 python 覆盖层时删 `python-source` |
| `pyproject.toml` 的版本号 | 手动同步 `[project] version` | **已漂移过一次**：0.3.0 发版时发现它还停在 `0.1.7`，落后四个版本，而 maturin 打 wheel 用的正是这个值——装出来的包会自称 0.1.7。已改用 `dynamic = ["version"]` 让 maturin 从 Cargo.toml 取，`[project]` 下不再写 `version`。**本机无 maturin，此改动未经构建验证**，首次 `maturin build` 时留意 |
| `maturin develop` | 无激活 venv/conda 时报错 | 用 `maturin build` + `pip install <wheel>`，或先激活环境 |

## 依赖升级陷阱（2026-08-08）

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| `rand` 0.8 → 0.10 | 除 `thread_rng()`→`rng()`、`gen_range()`→`random_range()` 之外，`random_range` 还从 `Rng` trait **移到了 `RngExt`** | `use rand::Rng` 必须改成 `use rand::RngExt`；只改方法名会得到 `no method named random_range found for struct ThreadRng` |
| `pyo3` 0.29 `#[pyclass]` + `Clone` | 隐式 `FromPyObject` 改为需显式选择，不选会有 deprecated 警告 | 本仓所有绑定函数按 `&PyTrajectory` 借用，无按值提取场景 → 选 `skip_from_py_object`。选 `from_py_object` 会让每次提取克隆整条轨迹，是不该发生的深拷贝 |
| `ferro-python` 独立 workspace | 有自己的 `Cargo.lock`，`cargo build/test/update --workspace` 全部跳过它 → lockfile 长期漂移、编译断裂不被发现（`MsdParams.fit_range` 自 0.1.8 起断裂到 0.1.11 才被发现） | 每次改 `ferro-analysis` 等公共 API 后，额外跑 `cd ferro-python && cargo check`；升级依赖时额外跑 `cd ferro-python && cargo update` |
| 0.x 版本号语义 | Cargo 把 `0.x` 的**次版本位**当破坏性位，`0.34 → 0.35` 属跨主版本，`cargo update` 不会碰 | 用 `cargo update --dry-run --verbose` 看 `Unchanged ... (available: ...)` 行，或装 `cargo-edit` 后 `cargo upgrade --incompatible` |

## 构建 / 工具链陷阱

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| `cmd/job.rs` 的 `-h` | 自己声明 `help` 字段接管 `-h/--help`，却没关掉 clap 自动生成的那个。**两个同名参数只在 debug 断言里报**：`cargo build` 通过、release 正常跑，而 **debug 版 `ferro job` 一跑就 panic**（`clap_builder/debug_asserts.rs:99`）。测试全绿也发现不了 —— 23 个 workflow 测试测的是 builder，从不经过 clap 解析 | `#[command(disable_help_flag = true)]`。2026-08-24 修。**任何自定义 `-h`/`-V` 都要同时关掉 clap 的自动版**；发现它的唯一途径是拿 debug 二进制真跑一次子命令 |
| `cargo fmt`（全仓） | 代码库未用 rustfmt 管理，全仓 `cargo fmt` 会改动 ~80 文件（含 22k 行生成的 `cp2k_basis_db.rs`），淹没特性 diff | 只对改动文件 `cargo fmt -- <file>`，或跳过格式化手工对齐周边风格；用 `cargo build`/`cargo clippy` 验证。若已误跑全仓，先提交真实改动再 `git checkout -- .` 丢弃格式化 |

## g(r) / box_builder / Bader 编码陷阱

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| `gr.rs` `calc_gr` | 默认 r_max=10.005 对小盒子超过 MIC 上界，CellList 每轴退化为 1 个 cell，最小镜像失效，g(r) 尾部被错误压向 0 | `calc_gr` 内把有效 r_max clamp 到 `Cell::minimum_image_cutoff()`（最小**面间距**的一半，非最短边长），并写回 `result.params.r_max`；`with_auto_rmax` 现由 ferro-python 使用，同样走该上界 |
| `box_builder.rs` `CellList::neighbors` | 27-bin 遍历缺去重，小盒子（n_bins ≤ 2）同一 bin 被多次映射 → 近邻重复计入、弛豫斥力成倍累加 | 加 `seen` bin 去重（与分析层 `angle.rs::CellList` 语义一致） |
| `bader_grid.rs` `bader_ongrid` | 曾有第一遍梯度上升结果被 `volnum=vec![0;..]` 立即清零重算（开发遗留），运行时翻倍 | 已删，只保留全路径追踪的单遍循环；改 on-grid 时勿再引入双重计算 |

## g(r) / S(q) 审查项处置（2026-08-08）

七项中六项已修（0.1.10–0.1.13、0.2.1），教训分别落在下方各节与
`network 重构（0.2.1）编码陷阱`。两条仍需知道的：

- `angle.rs` `neighbors_of` 每对被枚举两次后靠 `j<=i` 丢一半，**2× 计算浪费但不影响
  正确性 —— 决定不处理**
- `sq.rs` 的 Faber-Ziman 权重曾全 n² 遍历 + `factor=2.0`，异种对被加两次。实测
  `examples/70Z30P00A_NVT.lammpstrj` 大 q 处 `total_xrd` 旧/新 = **1.588**、
  `total_neutron` = 1.524，修后尾部收敛到 1.007。改权重公式后**必须验大 q 是否趋于 1**

## NPT 逐帧体积归一化（2026-08-08 已修，0.1.13）

**病根是口径不一致，不是局部笔误。** 参考实现（`code1/gr.c:139-146`、
`code2/dump2sq.c:392-406`）是**逐帧归一化/变换后时间平均**；Ferro 移植时把体积提到了
循环外（`⟨hist⟩·⟨V⟩`），只对时间平均后的 g 做一次变换。NVT 下两者严格恒等
（V 恒定，运算线性可交换）；NPT 下每处交换留下一个协方差项。

**关键恒等式：ρ·g 中体积自行抵消**

```
ρ_f · g_ij,f = (N_f/V_f) · ni·hist_f·V_f/(4πr²dr·N_A·N_B')
             = ni·N·hist_f/(4πr²dr·N_A·N_B')          ← V_f 消去

⟨ρ_f(g_f−1)⟩ = rho_g − rho,   rho_g := ⟨ρ_f·g_f⟩（纯计数），rho := N·⟨1/V⟩
```

所以逐帧变换的结果可由两个时间平均量**精确**重构，无须每帧各做一次 FT。
`GrResult` 因此有 `rho_g` 字段；fold 维护两份直方图 —— `plain`（Σcount，供 CN 与
`rho_g`）与 `vol_w`（Σcount·V_f，供 g(r)）。不是近似：两个测试与字面逐帧实现对拍到
1e-10 以内。

### 实测差异（完整轨迹，`-a P -b O`）

| 轨迹 | g(r) 最大绝对差 | S(q) 最大绝对差 |
|---|---|---|
| `70Z30P00A_NVT.lammpstrj`（V 恒定） | **0**（逐字节相同） | **0**（逐字节相同） |
| `43Z43P15A_NPT.lammpstrj`（σ_V/⟨V⟩=0.641%） | 3.06e-3（相对峰值 0.0088%） | 1.06e-2（相对峰值 0.82%） |

**S(q) 的偏差比 g(r) 大约 100 倍**，主体来自 ρ 的口径（`N/⟨V⟩` → `N·⟨1/V⟩`，本例相对差
3.9e-5，与 `⟨V⟩·⟨1/V⟩−1 = 0.003927%` 吻合）而非 g 本身。**结论：这次修复真正值钱的是
S(q) 一侧；g(r) 一侧基本是正确性洁癖。** 0.1.12 及更早跑出的 NPT `.sq` 与新版会有约
1e-2 的差异，NVT 结果完全不受影响。

### 附带修复

- 粒子数守恒校验：归一化把 N_A/N_B/N 当常数提出求和号，逐帧比对计数，不符即 `Err`
  （原先 `else { continue }` 会把首帧没见过的类型静默丢弃）
- `volume_std` 用**两遍算法**：`⟨V²⟩−⟨V⟩²` 在 V≈6e4 时抵消误差放大到 1e-3，定容轨迹也
  报不出 0；改 `Σ(v−⟨V⟩)²` 后降到 1 ulp（7e-12），测试改用相对容差判定「恒容」
- 文件头 `# volume = <mean> +/- <std> Ang^3`，用户可直接判断轨迹是否 NPT

### 明确不做

- `cube_density.rs:188` 的 `ref_cell.volume()` 是同类病，但它的 NPT 问题不止体积
  （格点在逐帧变形的盒子里含义就不同），需先定义「NPT 下 3D 密度图指什么」，单独排期
- 变粒子数（GCMC / 反应 MD）：会让各 partial 来自不同系综，`Σ w_ij·S_ij = total`
  这条 0.1.10 刚立起来的恒等式失效，需独占一次「total 怎么定义」的讨论

### 一处已修正的错误论断

曾断言「标签分组与元素分组的加权总 S(q) 逐点严格相等」。**不成立**：同种对 g_AA 的
归一化用 `N_A(N_A−1)`，而标签层求和重构出的是 `N_A²`，差一个 O(1/N) 的自排除项。
纯改名（每元素一个标签）是精确相等；位点细分则留有随 N 缩小的有限尺寸差。
测试 `test_label_rename_gives_identical_totals` /
`test_label_subdivision_differs_only_by_finite_size_term` 分别钉住这两种情形。

---

## element / label 分配规则（2026-08-08 查证）

`calc_gr` / `calc_angle` 默认读 `atom.element`，`GroupBy::Label` 时读 `atom.label`
（未设则回退 element）。除此之外 `atom.label` 的消费者只有
`ferro-io/src/writers/cif.rs:57,82`。各 reader 的分配：

| Reader | `element` 来源 | `label` |
|---|---|---|
| `xyz.rs:54` | 第 1 列原样，不校验 | `None` |
| `extxyz.rs:85` | `species` 列 | `None` |
| `pdb.rs:78-88` | 第 77–78 列；行短于 78 字符时退化为原子名首字符，取不到则 `"X"` | `None` |
| `vasp.rs:41-81` | POSCAR 第 5/6 行符号 + 计数展开 | `None` |
| `lammps_data.rs:94-100` | Masses 段注释 `# Fe`，无注释则按质量反查 | `None` |
| `lammps_dump.rs` | dump 的 `element` 列按第一个下划线拆分后的元素前缀；**无该列 → `X{type}`**（`X1`/`X2`…） | 有下划线且前缀合法时为整串（`O_b_P_P`），否则 `None` |
| `cif.rs:218` | `_atom_site_type_symbol`，缺失则从 label 取字母前缀 | `_atom_site_label` 原样 |
| `cp2k.rs:97-120` | `&KIND` 映射，或 label 字母前缀 | 坐标段第 1 列原样 |
| `qe.rs:113-129` | `ATOMIC_SPECIES` label 的字母前缀 | label 原样 |

**陷阱**：

| 位置 | 陷阱 | 说明 |
|---|---|---|
| `lammps_dump.rs` | dump 缺 `element` 列时元素名为 `X1`/`X2`，`symbol_to_z` 全部返回 255 | S(q) 的 XRD 系数与中子散射长度查表全部失败，`total_xrd`/`total_neutron` 无意义。0.1.11 起 `calc_sq_from_gr` 会打印告警（此前静默）。用 `-m sq` 前仍应确认 `ITEM: ATOMS` 头含 `element` |
| 三处元素提取 | `symbol_to_z`（精确→两字节→首字节）与 `cp2k.rs:221`/`qe.rs:193` 的 `extract_element`、`cif.rs:463` 的 `element_from_label`（字母前缀 + 首字母大写）**规则不一致** | 对 `ZnA`：前者 → `Zn`，后者 → `Zna`。新增的下划线拆分规则只作用于 LAMMPS dump，不动这三处（它们处理 `Fe1`/`O2` 原生命名，字母前缀规则对其正确） |
| `symbol_to_z` 贪婪匹配 | `Pb`→铅(82)、`Po`→钋(84)、`Os`→锇(76) 会盖过"P bridging"这类伪标签意图 | 这是选择「按第一个下划线拆分」而非贪婪前缀匹配的主要理由 |
| `typing.rs:22-23` | `apply_type_labels` 把类型标签写进 **`element`** 字段而非 `label` | 伪元素轨迹即由此产生；`elem_z` 的前缀匹配就是为它服务的 |
| `cif.rs` 位点 | `Fe1`/`Fe2` 的 `element` 均为 `Fe`，g(r) 无法分辨不等价位点 | 用 `-x/-y` 按 `label` 选（CIF reader 已填 `label`）。注意 `Fe1` 不合 `<元素>_<后缀>` 约定，导出 LAMMPS dump 时**不会**被折进 element 列（0.2.1 的守卫），走 extxyz 的 `label:S:1` 列才无损 |

---

## 外部参考工具的差异与缺陷（2026-08-09 查证）

`ferro traj` 的 gr / sq / angle 以 `private/code1`（dump2analysis）与 `private/code2`
（dump2sq）为参考实现。三处差异全部定位，其中一处是参考程序的真 bug。**下次有人问
"为什么 ferro 的 S(q) 和 dump2sq 对不上"，先看这一节。**

### dump2sq 的 partial 按数组下标归类（对方的 bug，不复现）

`dump2sq.c:InitializeType` 只在**第 0 帧**写入 `data[i].type_new`；
`input.c:ReadDataLammpstrj` 的 `sscanf` 此后只更新 id/type/elem/x/y/z，**从不更新
type_new**。于是原子类型按数组下标而非实际元素归类。LAMMPS 未加 `dump_modify sort id`
时原子书写顺序随 domain decomposition 变化，第 1 帧起归类即失效。

实测对照：

| 轨迹 | 帧间同位置同 id | partial 平均互相关 (q>2) | XRD vs ND 相关 | fe 与 d2sq 的 max\|Δ\| |
|---|---|---|---|---|
| `examples/43Z43P15A_NPT.lammpstrj`（未排序） | 210/2004 | **0.959** | 0.9986 | 0.31 / 0.14 |
| `tests/70Z30P00A_NVT.lammpstrj`（已排序） | 5003/5003 | −0.118 | 0.4280 | **0.0008 / 0.0027** |

决定性验证：`.gr` 的 `Total` 列不经过 `type_new`，是全文件唯一干净的一列。拿它自己做
傅里叶变换，与 dump2sq 输出的两条总曲线相关 0.99868（XRD）/ 0.99993（ND）——
**dump2sq 的加权总 S(q) 实际就是总 g(r) 的傅里叶变换**，Faber-Ziman 加权被完全冲掉。

处置：**不复现**。`scripts/trajcheck.py` 集中记录成因与数据，`compare_sq.py` 与
`compare_sq_experiment.py` 每次运行自检并在终端与图上警示。dump2analysis 每帧按元素
重新选原子，不受影响。

### angle 的两处约定差（保留 ferro 现有行为）

| | ferro | dump2analysis |
|---|---|---|
| 计数 | 每个几何角一次（PO₄ 给 6 个 O-P-O） | 有序端对，两端**同元素**时数两次；不同元素时一次 |
| bin k 覆盖 | `[kΔθ, (k+1)Δθ)` | `[0.05+kΔθ, 0.15+kΔθ)` |

计数因子不影响 mean / std / 峰位。半格偏移源自 `angle.c:OutputAngle` 的
`min_angle = 0.05`；那套分箱会把 <0.05° 的角写到 `a[-1]`（越界），上界也延到
180.05°，故不采作默认。0.1.15 起 `--angle-min` / `--angle-max` 可调，
`--angle-min 0.05 --angle-max 180.05 --d-angle 0.1` 可逐点复现——实测 O-P-O 与
O-Al-O 各 1800 个 bin，补上 ×2 后整数计数**零差**。

### S(q) 缺 +1（在比较脚本侧补偿）

`dump2sq.c:CalcSq` 只算 `4πρ/q · Σ_r r[g−1]sin(qr)Δr`，不含定义里的 +1。实测残差在
q>15 处为 1.00031（XRD）/ 1.00018（ND），标准差 0.014 / 0.002，确认是纯常数。比较
脚本默认给 dump2sq 加 1（`--raw` 可关），因为 S(q) 定义要求大 q 趋于 1，ferro 已是
标准形式。

### 两条曾经写错、已更正的论断

- `compare_sq.py` 早先把上面第一条描述成"dump2sq 的 g(r) 有伪峰"——伪峰是症状，
  成因是 type_new 陈旧。
- `compare_angle.py` 早先把 bin 偏移写成"横轴标注约定不同、同一索引覆盖同一区间"
  ——错的，区间真的错开半格。

---

## angle cutoff 归属（2026-08-09 已修，0.1.15）

`--r-cut-ab` / `--r-cut-bc` 原先按**原子序数**派发给两端，静默忽略用户在 `-a`/`-c`
写的顺序（`-a Zn -b O -c P --r-cut-ab 2.5` 会把 2.5 给 P，因 Z=15 < 30），与意图相反。
两个 cutoff 相等时不显现，所以一直没暴露。

修法：`AngleParams.ends: Option<(String,String)>`，CLI 点名三元组时把 `-a`/`-c` 原样
传下去；未点名时（扫全部三元组）无用户顺序可依，仍用规范 (Z, 符号) 顺序。

**两端同类型时 `ends` 不生效**，一律取 `min(r_cut_ab, r_cut_bc)` —— 谁被排成 lo 端
取决于邻居枚举顺序，只有取 min 才可复现（测试 `test_symmetric_ends_ignore_order`）。

---

## Table / 批处理编码陷阱（2026-08-09）

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| `Table` 的缺失值 | 列并集下某输入缺某列，补 0 会被当成「测到了 0」 | `Column::Num` 存 `f64::NAN`，`cell()` 渲染为**空字段**，pandas 读回 `NaN`。空 ≠ 零，这条是列并集能被接受的全部理由 |
| `Summary` 的数值 | 用 `Column::Num` 会让「5 帧」写成 `5.000000e0` | `[inputs]` 是给人看的注释块，值**预格式化成文本**（`fmt_compact`：整数保持整数，其余四位有效数字，非有限值写 `-`） |
| `meta_lines()` | 把逐文件才有意义的量（组成、clamp 后的 r_max）写进共享参数块，多文件下第一个文件的组成会被当成全局事实 | 共享参数留 `meta_lines()`，逐文件的走 `Summary::note()` 进 `[inputs]` 清单 |
| 按输入数量分派单/批两条路 | 产物形态会取决于 glob 当天匹配到几个文件；平时 5 个走批处理、某天只剩 1 个就悄悄换形态 | **单一代码路径**，N=1 是 N 的特例。两条路还意味着两套命名、两套绘图、两套错误处理 |
| `expand_inputs` 零匹配 | 静默跑零个文件是最难查的失败，尤其在脚本里 | 零匹配**报错**并回显模式 |
| 批内失败 | 只打印不改退出码，shell 里 `ferro traj gr ... && next` 会把失败当成功 | 跳过 + `[inputs]` 留原因 + **退出码 1**；参数级错误仍在读第一个文件前快速失败 |
| `glob` crate | 不支持 `{a,b}` 花括号 | 让 shell 展开（展开后的字面路径原样通过）；文档写明 |
| 多输入下的 `.cube` | 文件名不带 stem 时第二个输入直接覆盖第一个 | `ferro map` 的文件名必须掺输入 stem（`density_<stem>.cube`） |
| plotters 子图 | `root.split_evenly((rows, cols))` 后每格要各自 `ChartBuilder`，且没有 figure 级图例 | 颜色按文件分配保证跨格一致，图例只画第一格 |

## 外部脚本调用 ferro 的陷阱（2026-08-11，`scripts/ferrocmp.py`）

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| `-o` 的语义 | 它是**文件名后缀不是路径**。`subprocess.run([... "-o", str(outdir/"x")])` 会写出 `gr_<绝对路径>.csv` 这种文件名，或直接失败 | 以 outdir 为**工作目录**调用（`cwd=outdir`），输入路径转绝对；产物名自己拼 `<mode>[_<表>]_<后缀>.csv` |
| 复跑的残留 | 上一轮的产物若还在，新一轮参数写错导致 ferro 没写文件时，脚本会拿旧文件当本轮结果，数字全对但结论是错的 | 调用前 `unlink()` 预期产物，再验存在 |
| 退出码 | dump2analysis / dump2sq 参数校验失败走 `exit(0)`，退出码不可信；**ferro 的可信**（批内失败也置 1） | 两种检查分开：C 程序只验产物，ferro 同时查 `returncode` |
| 读长表 | 按列名取列的老写法在长表下取不到东西；直接 `df.iloc[:, 4]` 又会在列顺序变化时静默错位 | 按数据列选行（`center` / `neighbor` / `end_a`…），选空时把该表实际有的组合打印出来 |
| `file` 列 | 单输入脚本遇到多输入产物时静默取第一个，后面所有数字都对不上且图上看不出来 | 断言 `df.file.nunique() == 1` 后再丢列 |
| 产物名 | 在调用点手写 `f"gr_{tag}.csv"` | 用 `ferrocmp.product_name()` 拼（`batch::out_path` 的镜像）。2026-08-13 加 label 段时，三个脚本各写各的名字，同时断掉 |
| `-o` 与 label 重复 | 把配对同时塞进 `-o`（`-o P-O` → `gr_P-O_P-O.csv`） | 配对已由 label 段进文件名，`-o` 只留批次标记（`cmp`）。这正是加 label 段换来的东西 |
| `--ferro` 传相对路径 | 直接交给 `subprocess` | `run_ferro` 以 outdir 为 cwd，`../target/release/ferro` 会 `FileNotFoundError`。含分隔符时先 `resolve()` |

## 产物命名 / --outdir 编码陷阱（2026-08-13）

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| 删 clap 参数 | 只改函数体、不删结构体字段 | 两处都要改。`run_sq` 不再读 `c.select` 之后，`traj sq -a P -b O` **仍被接受**并静默产出 `sq.csv`；删掉 `SqCmd` 的 `#[command(flatten)] select` 之后 clap 才报 `unexpected argument`。「参数没人用了」不等于「参数没了」 |
| 有序选择 vs 集合选择 | 用同一条规则拼名 | 分开：`gr`/`angle` 按写的顺序（`CN` 有向，`-a P -b O` 与 `-a O -b P` 是两份不同数据，必须两个文件）；`--elements` 排序去重（是集合，`O,P` 与 `P,O` 选中同一批原子，同一份数据不能落到两个文件名下）。统一成任何一种都会错一半 |
| 用户串进文件名 | 直接拼 | 先校验 `[A-Za-z0-9_+-]`，**在读第一个文件之前**报错。选中的串会成为路径的一段，`-a 'P/2'` 不校验就写到别处去了 |
| 非法字符的处理 | 替换成下划线 | 报错。替换会让 `-a P/2` 与 `-a P_2` 静默写进同一个文件 —— 与补零冒充「测到了 0」同一类错误 |
| `--outdir` 的创建时机 | 写第一个文件时才建 | `Output::prepare()` 在读第一个输入**之前**调用，与「参数级错误快速失败」一致。否则跑完一小时分析才发现路径打不开 |
| PNG 的目录与 label | 在绘图侧再接一遍 `--outdir` | `png_path` 收数据文件的**完整路径**（`&Path` 而非 `&str`），目录与 label 自动跟随。绘图侧不该知道命名规则 |
| `rotcorr` 的「无选择」分支 | 以为 `--center`/`--neighbor` 可省 | 两者都是必填（`Option` 只为「不带 `-i` 时打帮助」，随后 `bail`）。它恒有 label，`rotcorr_all.csv` 走不到，别为它写代码 |
| `write_all` 的参数个数 | 继续加位置参数 | 加到第 8 个时收进 `Output { dir, label, suffix }`。三个字段总是一起走，分开传只会让调用点一串 `None, Some(x), None` |
| `ferro-io` 的 writer 路径 | 以为都收 `&Path` | 九个 writer 全收 `&str`，`batch` 内部是 `PathBuf`，只能在边界转（`Output::join_str`）。统一成 `&Path` 是待办，见 `dev/plan.md` |

## 发表级绘图脚本编码陷阱（2026-08-13，`scripts/plot_*.py`）

样式栈是 `['science','vibrant']` + `usetex`。三条实测结论先摆前面：

| 事 | 结论 |
|---|---|
| 字符串转义 | **不需要**。matplotlib 的 usetex 前导会加载 `underscore` 宏包（`[strings]`），`43Z43P15A_NPT_5`、`P_3`、`O_b`、`50% Zn` 全部原样通过 |
| 仍会炸的字符 | `&`、`#`、`^` —— 让**整张图渲染失败**（LaTeX 报错，不是画歪）。元素符号与位点标签不可能含它们，只有自定义文件名会碰上 |
| usetex 速度 | 一张图约 0.5 s，不是瓶颈 |

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| `science` 样式 + 双轴 | `ytick.right: True` 会在右侧画**主轴的镜像刻度**，跟 twinx 的 CN 刻度叠在一起 | 主轴显式 `ax.tick_params(right=False)`，把右边让给 CN 轴 |
| 堆积柱的「段」 | 把 csv 的**每一行**当一个段 | 一个段是一个**物种**。`qn_partner` 里同一个 `(qn, m_<X>…)` 组合在每个成分各有一行，不归并的话明度档按行号分配 —— 颜色就不再编码 m 而是编码「第几行」，第二个成分整根柱子变淡 |
| `m_<X>` 列的 NaN | 照 `== key` 匹配 | 先 `fillna(0)`。列并集下某成分缺整列 `m_Al`（那个体系没有 Al），空字段读回是 NaN，匹配不上就**整根柱子静默消失**。一般不该把 NaN 当 0，这里 NaN 确实等于 0（没有 X 就没有连向 X 的桥），是有明确语义的例外 |
| 多子图共享图例 | 每格各自 `color_map(该格的类别)` | 颜色必须**跨全图按类别值分配**。否则 Al 格的 cn=4 与 Zn 格的 cn=3 都拿到色序第一个颜色，一份共享图例就是在撒谎。段标签也不能带元素名（`Al_4`），否则 Zn 格里同一个橙色挂着 `Al` 的名字 |
| 明度档数 | 按第一组的段数算 | 取**所有格、所有组**的最大值。不同成分下伙伴分解行数不同（某成分一个 `Q2(m_Al=2)` 都没有），只按第一组数会在别的组上索引越界 |
| 图例位置 | `ax.legend(bbox_to_anchor=(0.5,-0.28))` | 挂到 `fig.legend` 上。成分名是旋转的长标签，挂在 ax 上会被它们顶穿 |
| 多组并排时的标题 | 默认 pad | 组标签（CMD/MLMD）竖排在柱顶会顶到标题，`n_grp > 1` 时给 `title_pad` |
| g(r) 的 CN 右轴 | 用默认自动范围 | r 到 8 Å 时 CN 已到 110，把一壳层平台压成一条线。这正是 `CFG["cn_ylim"]` / `xlim` 存在的理由，出图时按体系设 |

## 高 DPI 绘图陷阱（2026-08-10）

`--plot` 仍出 PNG，但分辨率由 96 dpi 提到 **500 dpi**（一格 2708×2083 px）。
版式在 96 dpi 像素下编写，统一经 `px()` 缩放，故改 `DPI` 只整体放缩、不重排。

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| 只放大画布 | 画布放大了但字号/线宽/边距留在原值，图会变成「巨大画布上一堆蚂蚁字」——比不放大更难看 | **每个长度都要过 `px()`**：`caption`、`margin`、`x/y_label_area_size`、`stroke_width`、legend 色条 |
| `configure_mesh` 字号 | 刻度与轴标题字号有**独立默认值**，不跟随 `caption`；漏掉就小到看不清 | 显式 `.axis_desc_style(("sans-serif", px(14)))` + `.label_style(("sans-serif", px(12)))` |
| plotters legend 区宽 | `legend_area_size` 默认 30px 是**固定值**，不随画布缩放；`px(18)` 的色条会伸出去压在标签文字上 | 显式 `.legend_area_size(px(30))`，与色条同比缩放 |
| 内存与体积 | 2×2 的 msd 图是 5416×4166 px，RGB 缓冲约 68 MB，PNG 约 675 KB | 目前可接受；再往上调 `DPI` 前先算 `w×h×3` |
| 验证「分辨率真的生效了」 | 像素数写错时，扩展名、格式、肉眼缩略图全都看不出来 | 测试 `test_render_writes_png_at_full_resolution` 直接从 PNG 的 IHDR 读宽高（偏移 16/20，大端）断言等于 `px(520)`/`px(400)` |

### 为什么没有采用矢量输出（2026-08-10 定案）

曾实装矢量 PDF（plotters `SVGBackend` → `svg2pdf::to_pdf`，提交 `a4c3c11`）并验证可用：
真矢量、字体内嵌、放大无损。**因依赖成本回退**——`svg2pdf` 会拉进
usvg / resvg / fontdb / tiny-skia 等约 40 个 crate，与计算无关，不值得为一个自检图付。

若将来重新考虑，已知结论存档：

- plotters **没有** PDF backend（只有 bitmap 与 svg），换矢量不是改扩展名的事
- `svg2pdf` 只需 `default-features = false, features = ["text"]`（绘图 SVG 无位图无滤镜）
- `usvg` 解析不到 `font-family` 时**静默丢弃全部文字**、不报错，产物仍是合法 PDF 却没有
  任何标注——必须显式检查 `opt.fontdb.is_empty()`（读是字段，写才是 `fontdb_mut()`）
- 矢量下 `DPI` 不再是画质而是尺度：`svg2pdf` 按 `page = px × 72/DPI` 定页面大小
- 500 dpi 位图约可清晰放大到 5×，够读峰形；需要超出这个范围的，正途是走 Python 出图

## 分析产物为什么不进 ferro-io（2026-08-09 定案）

`ferro-io` 只依赖 `ferro-core`，而 `write_gr` 的入参 `GrResult` 在 `ferro-analysis`
—— 直接把 writer 搬进 io 会让分层图出现横边，违反「中间层不得互相依赖」。

成熟工具无一例外走「中立载体」：ASE 的结果对象自带 `write()` 走通用 `jsonio`；
pymatgen 的 `MSONable` + `MontyEncoder`；OVITO 的 `DataTable` 被所有 modifier 产出、
被所有 exporter 消费。ASE 与 pymatgen 是 Python，**能**写出循环依赖却依然不写，
说明这是设计选择而非语言约束。

结论：规则不用改，缺的是反方向的共享类型。补上 `Table` 之后
`io --Trajectory--> analysis` 与 `analysis --Table--> io` 对称，`GrResult` 继续留在
analysis（io 没理由认识它），跨层的只有它的可序列化投影。

---

## 单二进制重构陷阱（2026-08-09，0.2.0）

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| `cmd/net.rs` | `--P-O=2.3` 把元素对写在**参数名**里，clap 无法建模成固定 flag 集 | `main` 在 `Cli::parse_from` 之前用 `split_pair_args` 从 argv 剥离，剩下的交给 clap |
| `cmd/job.rs` | 原 `main` 按值消费 `args` 的各字段（`if let Some(m) = args.method`），改成 `&JobCmd` 后全是 E0507 | `#[derive(Clone)]` + 函数入口 `let args = args.clone();`，函数体一行不用改 |
| `main.rs` 的 `Command` 枚举 | `JobCmd` 有 30+ 字段，clippy `large_enum_variant` 报变体大小差异 | 大的那个装箱：`Job(Box<JobCmd>)` |
| 子命令帮助 | clap 的 `--help` 派生格式塞不下「输出列结构」「长表为什么长」这类段落 | 保留 `disable_help_flag` 风格的自定义 `help.rs`；`ferro traj gr` 无 `-i` 时打印富文本 |
| 三级帮助的第一级 | clap 不支持给子命令**分组**显示 | 帮助本来就是手写的，`print_overview()` 自己排版分类 |

> 已知限制清单见 `progress.md`。一条容易混淆的留在这里：
> `ferro traj sq` 的 `--q-min` 只影响**输出 q 网格**，g(r) 的积分范围由 `--r-min` /
> `--r-max` 决定。

## 与实验对比偏差大时的排查清单（2026-08-20）

计算结果与 NMR / 衍射对不上时，**先看这两条口径选择**，它们都不是 bug，
都是有意的建模决定，但都能把数字改变一个量级。按怀疑优先级排列：

### 1. Zn 被当作修饰子还是形成子（量级最大）

实测 `43Z43P15A_NPT_5`，同一条轨迹两种口径：

| | mean_qn(总桥) | O_n 非桥氧 | O_b 桥氧 | O_t 三簇 |
|---|---|---|---|---|
| `--modifier Zn`（现在的用法） | 2.40 | 2985 | 3583 | 2 |
| Zn 当形成子（`--Zn-O=2.8` 进 cutoffs） | 3.96 | 70 | 5607 | **893** |

判为 `O_n` 的氧里 **97.7% 连着 Zn**，所以这个开关几乎重新定义了整张网络。
文献立场：ZnO 是 **intermediate oxide**，四配位时进网络、八面体时作修饰子，
而该体系 72.6% 的 Zn 是四配位。

**标准做法是 Zn 作 modifier**（31P NMR 的 Qn 按 P–O–P 连接定义，Zn–O–P 不计），
ferro 的默认用法与之一致。仅当与实验偏差异常大时才回头怀疑这一条。

### 2. n 的「同核」边界取严格同元素

`qn` 只把**同元素**连接计入 n（P–O–P）。文献里这个边界实际是「同属主骨架」，
由作者按体系定：arXiv 2510.13545 把 **P–O–Si 也算进 n**（Si 与 P 同为四面体骨架），
而硼磷酸盐文献把 B 算进 m（`Qⁿ(mB)`）。

ferro 选严格同元素，因为它覆盖 `Qⁿ(mAl)` / `Qⁿ(mB)` 这两个最常见记号，
且规则一句话说得清。**代价**：若对照的文献把某个异种形成子并进了 n
（如含 Si 的体系），ferro 的 n 会系统性偏低，差额恰好是 `m_<该元素>` 列。
不丢信息，加回去即可，但对比前必须先确认对方的 n 包含哪些元素。

---

## 数据集 merge 编码陷阱（2026-08-25）

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| 需求表述有歧义时 | 挑一种读法直接实现 | 「每个 system 对应一个 set」既可读作「恒为一个」也可读作「不跨 system 地切」，两次改动才落到后者。含「每个 X 一个 Y」这类量词的需求，先把边界情形（X 比 Y 大很多时怎么办）问清楚再动手 |
| 分组依据 | 按目录名（`*.train`、`init.011`） | 名字说明不了里面装的是什么。按**逐原子的元素序列**分组，成分相同才合并；顺带让异成分体系自动分开 |
| 原子重排 | 只换 `type.raw` 不动数据 | 逐原子数组（`coord`、`force`）必须跟同一个置换走，而与原子编号无关的量（`box`、`energy`、`virial`）**不能**跟着动。漏掉 force 会让力对错原子，数值全在、量级正常、完全看不出来 |
| 规范序的选择 | 随手选一个 | ferro 取 `(Z, 符号)`（与 `collect` 的 type_map、`traj gr` 的分组序一致），dpdata 取字母序。两者都靠 `type_map.raw` 自描述所以都对，但**同一份数据经 ferro 与经 dpdata 合并，`type.raw` 的数字不同** —— 混用两条工具链时别拿旧的 `type.raw` 配新的 `type_map.raw` |
| by-source 的语义 | 拼完再按 `--set-size` 均匀切 | 切分要**在每个 system 内部**进行：set 不横跨 system，每个 set 仍可追溯到单一条件；`--set-size` 照常生效。writer 因此需要一个显式边界的入口 —— 全局均匀切表达不了「边界落在 system 边界上」 |
| 切分的余数 | 留在最后一个 set | 500 帧按 400 切要给 250+250 而不是 400+100。一大一小对训练负载和「拿一个 set 当验证集」都更差；这条对两种模式都适用 |
| 打乱的可复现性 | 用系统随机源 | seed 有默认值（666），不给 `--seed` 也能复现。数据集是要归档的，「这一版是怎么打乱的」必须答得出来 |
| 后缀继承 | 无条件加 `.train` | 只在组内**全部**输入共享同一个划分后缀时继承；混合后缀或无后缀（dpgen collect 那种）就不加，否则等于凭空断言这是训练集 |

## 数据集 filter 编码陷阱（2026-08-25）

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| 筛选判据放哪 | 「阈值比较只有三行，放 CLI 就行」 | **代码量不是分层判据**。项目所有分析方法都在 ferro-analysis，判据放 CLI 会让 `cmd/dataset.rs` 成为唯一自带算法的子命令。更要命的是连带效应：判据在 CLI 时，承载类型看起来只被 io 产出、cli 消费，于是「不下沉 core」的结论也跟着错 —— 判据回到 analysis，`Trajectory` 立刻满足「两个以上中间层要叫出它的名字」 |
| 数据集的承载类型 | 造一个 `DeepmdSystem` 放 ferro-io | analysis 的筛选函数要收它，`analysis → io` 就是横向依赖，直接违反铁律。而且类型带格式名，等 GPUMD 的 extxyz 进来就是第二个类型。用已下沉的 `Trajectory`：extxyz reader 本来就产出它，npy reader 也产出它，中立性自动成立 |
| 区间的基准 | 作用在原始帧号上 | 作用在**存活序列**上。`--start 10` 是「通过质量判据的第 10 帧」—— 在丢掉未知多少帧之后，这是唯一还讲得通的语义。原始索引由 `keep` 带出，保留帧仍可追溯 |
| 判据贡献的度量 | 看漏斗每步剔掉多少 | 漏斗每步只在**上一步的存活帧**上报数，于是一个只会重复抓别人已抓帧的判据看着也很能干。要逐帧对**全部**帧算判定，再统计「独占抓到」多少。实测 `-f 3 -s 0.2`：force 判坏 1978 帧、独占 **0**，漏斗（2000→22→0）完全看不出它冗余 |
| 表里的计数 | 用 `Column::Num` | 会把「2000 帧」渲染成 `2.000000e3`。整数计数走 `push_text` 预格式化 —— 与 `Summary` 那条是同一个坑，只是换了个地方再踩一次 |
| set 切分的余数 | `chunks(set_size)` 朴素切 | 410 帧按 400 切给出 400+10，那个 10 帧的 set 当验证集没有统计意义。余数应摊进各 set（205+205）|
| 重写 system 目录 | 直接覆盖各 npy | 帧数变少时旧的 `set.001`、`set.002` 会留在原地，读回来帧数还是旧的。写之前先清掉全部 `set.*` |
| `expand_inputs` | 拿它展开目录 | 它只收 `is_file()`，给目录会报 "no file matches"。DeePMD 的 system 是目录，用 `expand_dirs` |
| npy 的 dtype | 只按自己写出的 float64 读 | dpdata 默认写 float32。两种都要收，别的 dtype 报错而不是猜 |

## CP2K out 解析 / DeePMD 数据集编码陷阱（2026-08-25）

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| 文本锚点 | 字面匹配（`starts_with(" MD\| Step number")`） | **按 token 序列匹配**（`["MD\|","Step","number"]`，`split_whitespace` 后逐词比）。CP2K 每个版本都会动列宽与对齐，字面匹配把解析器钉死在一个版本的空白上；而空白正是 `split_whitespace` 已经扔掉的东西。所有 tag 收进 `mod tag` 一张表，**每个 tag 是候选列表**，支持某版本改了措辞就是加一行，不是在扫描器里加分支 |
| 数值行的列 | 写死下标（`f[2..5]`、`fields[2..11]`、恰好 12 字段） | **从尾部取或按可解析性筛**。应力块靠「这一行能解出三个浮点」自动跳过表头与 `1/3 Trace`、`Determinant` 摘要行，不必知道有几行表头；cell 行末位是体积、其前九个是晶胞，从尾部取则前缀增减一列都不影响 |
| xyz 块的注释行 | 要求它以 `i =` 开头 | 注释行的措辞不是 xyz 格式的约定，版本间无保证。判据改成**结构性**的：一行原子数、一行解析**不成**原子行的（这才叫注释）、然后那么多行能解析的原子行。首尾各验一行即可，中间留给真正读取时逐行解析 |
| `nalgebra::Matrix3` 落盘 | 用 `as_slice()` 取九个数 | 它是**列优先**（无开关），给出的是转置；而 npy 的 box/virial、extxyz 的 `Lattice`、GPUMD 的 `lattice` 全要行优先。走 `matrix3_row_major()`。**这个错误有一半是静默的**：stress 对称，转置后数值不变，于是你会在 cell 上发现 bug、修掉、以为 stress 也对了。测试必须用**非对称**矩阵 |
| CP2K 的应力单位 | 按版本查表推断 | `STRESS_UNIT` 是 **`&PRINT &STRESS_TENSOR` 下的输入关键字**，同一个二进制能打 bar / GPa / atm —— 单位不是版本的函数，版本表天生解不了。从 `STRESS\| Analytical stress tensor [bar]` 的方括号自读，认不出**报错**。dpdata 硬编码 `/ GPa`（`formats/cp2k/output.py`），遇到 bar 输出静默错 5 个数量级 |
| 力的单位 | 也想从文本读 | 力是唯一**没有**单位标注的量（xyz 块的注释行只有 `i = …, time = …, E = …`），只能按 a.u.（Hartree/Bohr）兜底。这是「文本没写时才回落」的唯一实例 |
| 坐标块 vs 力块 | 靠内容区分 | 两个块**逐字同构**：同原子数、同 `i = …` 注释行、同元素列。只有**顺序**能分辨（坐标在前），所以扫描必须有上界 —— 限定在 `[锚点_i, 锚点_{i+1})` 内。没有上界时，某帧缺块会让扫描跨界抓到下一帧的数据，**数值完全正常只是整体错位一帧** |
| 块定位 | 用固定行偏移（`A + 10` 这类） | 偏移随 ensemble 变（NVT/NPT_I/NPT_F 的坐标块偏移是 10/18/20），CP2K 多打一行（`MD\| Estimated peak process memory` 已经在输出里了）整张表就错位，而错位读到的是别的块的数字、不报错。改为从锚点扫描特征行；偏移只用来**校验**（与首帧不同则告警） |
| `STRESS\|` 与 `MD\| Pressure` | 当成同一个量 | `STRESS\|` 是**势能部分**，`MD\| Pressure` 含动能项。实测（`total.out` 第 1 帧）：`(2/3)E_kin/V = 14134` bar，`(1/3)Tr σ = 2345` bar，两者相加 16479 ≈ `MD\| Pressure` 的 16480.9。训练集只能用 `STRESS\|` —— 用错等于把动能烘进势函数，别的温度下系统性错压 |
| stress 的符号 | 以为 CP2K 与 ASE 同号 | CP2K / VASP / QE **三者一致，且都与 ASE 相反**（正 = 压缩）。证据：`ase/calculators/cp2k.py:363` 的 `-1.0 * stress  # cp2k uses the opposite sign`、`ase/io/vasp.py:630` 的 `stress *= -0.1*GPa`、`ase/io/espresso.py:2098` 的 `# sign convention is opposite of ase`。ferro 存原样符号，故 `virial = +stress×V` 不变号；**导出 GPUMD 的 `stress=` 时要变号**（那个关键字是 ASE 约定），`virial=` 不用 |
| 重启拼接的 out | 按序号取第 k 个 `ENERGY\|` / `STRESS\|` | 重启时会**重打一次初始值**：`total.out` 有 2002 个 `ENERGY\|` 但只有 2000 帧，按序号取在重启点之后整体错位一位。用「从锚点向上找最近的」，天然不受影响。段的标志是 `MD_INI\| MD initialization`（`GO CP2K GO` 在新版本里根本不打印，dpdata 的段分隔符不可用） |
| 初始构型 | 以为会混进数据集 | CP2K 的初始构型打印在**第一个 `MD\| Step number` 之前**，按锚点定位天然排除；且重启**不**重打构型，故不存在重复帧。实测 2000 个锚点区间内各恰好 2 个 xyz 块，无一例外 |
| DeePMD npy 的形状 | 照文档写 `(nframes, natoms, 3)` | 磁盘上**一律二维**：dpdata 落盘前统一 `reshape([nframes, -1])`，文档里的三维是逻辑形状。`energy.npy` 是 `(nframes, 1)` 而非 `(nframes,)` |
| npy 的精度 | 跟 dpdata 默认的 float32 | `collect` 是流水线的**头**，filter 与 merge 都读它；头上有损无法在下游还原。写 float64，降精度留给喂训练框架那一步 |
| 一个 system 的组成 | 把多条轨迹并进一个目录 | `type.raw` 一个目录只写一次，故各帧的原子数与类型序列必须逐项相同。合并是 `merge` 的职责，`collect` 恒为一输入一目录 |
| 输出目录名 | 用文件 stem | CP2K 日志常常**全叫 `total.out`**，靠目录区分。stem 撞车时改 `<父目录>_<stem>`，仍撞则报错 —— 静默覆盖会让你以为收集了两条轨迹 |

## collect 目录语义 / 帮助防漂编码陷阱（2026-08-26）

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| 多输入的产物命名 | 把多级路径用分隔符压平（`a_md`） | **保留目录结构**（`a/md`）。`filter` 的 `-o` 本来就按相对路径重建，`find_systems` 与 `merge` 也都递归下钻 —— 压平会是这条链上唯一的例外。命名规则取「剥掉公共祖先后剩下的层级」：公共前缀按定义不携带区分信息，剩下的必然唯一 |
| 撞名 | 报错要求用户改名 | 分组单位选对了，撞名**就是「该合并」的定义**。`collect` 按目录分组之后，同目录的两个 out 本来就该进一个 system，原来那条 `bail!` 整条消失。**这也意味着命名与分组不能拆成两个提交** —— 中间态是「撞名报错但名字已变」，没法对拍 |
| 分组键 | 用路径原样 | `a/x.out` 与 `./a/y.out` 是同一个目录却分成两组。键取 `canonicalize`，**但显示与消息用用户写下的路径** —— 规范化后是绝对路径，刷屏且认不出 |
| 合并多文件的排序 | 逐帧按 step 全局排序 | 看着最彻底，但重启若不从 checkpoint 而是新起一个 job，step 会从 0 重新计数，全局排序会把两段真实轨迹**交错洗牌**，比不排序更糟。按文件的首 step 排、文件内保持原序：正常重启下与全局排序等价，重号下退化成「按文件拼」，最坏情况不比现状差 |
| 重叠帧 | 一定要去重 | 重启只重跑 checkpoint 以来的几步，位置速度相同则能量力也相同，不稀释结果。**不去重，但把每个源文件的 step 区间打出来** —— 这个前提只有区间可见时才检验得了 |
| 一个目录里成分不一致 | 当作坏帧丢掉（复用既有的 `composition changed` 计数） | **报错并指名两个文件与各自组成**。「把两个体系放进了一个目录」是人的错误，「这一帧的元素列坏了」是数据的问题，两者的处理动作完全不同；渲染成后者会让人去查 SCF 设置 |
| 合并语义下的单文件失败 | 整个目录失败 | 跳过坏文件、用剩下的建 system。但 **system 目录看着完全正常而帧数少了**，所以那份跳过清单要在最后的汇总里**再说一遍**，不能只在中间滚过去 |
| 产物为目录树的命令 | `-o` 给个 `.` 的默认值 | 改必填。默认 `.` 会把 npy 撒进正在工作的目录；且与「名字为空时写进 `-o` 本身」叠加后，不给 `-o` 的那条路会被 `--overwrite` 恒定拦下，等于一条走不通的路 |
| 报告文件放哪 | 塞进 `<outdir>/report/` 子目录显得整齐 | **平铺**。`expand_dirs` 只收 `is_dir()`，平铺的 csv 会被后续 `merge -i clean/*` 自动滤掉，而 `report/` 子目录反倒会被收进去当 system 候选 |
| 堆叠多 system 的行标签 | 用 `label_of`（取 stem） | 嵌套结构下 `a/md` 与 `b/md` 的叶子名相同，堆起来分不出是谁。用**相对路径**。`batch::Summary` 因此加了 `failed_one` 与 `into_table_named` |
| 「诊断只在只读模式算」这类省时优化 | 凭直觉认定它贵 | **先量**。实测 1110 帧 / 302 原子挂钟 0.23 s（带诊断）vs 0.28 s（不带）—— rayon 跑满六核，差别在噪声里。这个优化的代价是那几张表永远落不了盘 |
| 批内失败的提示语 | 承诺失败原因在某个文件里 | `[inputs]` 块只有分析命令的 csv 才有，`dataset` 的产物是目录。**指错地方比不指路更糟** —— 用户会照着去找一个不存在的东西。原因每个 SKIP 行都已说过，这句的唯一职责是解释退出码 |
| clap 定义与手写帮助页 | 靠人记得两边都改 | 加测试。`main.rs` 的 `mod help_sync` 正向查「每个长选项都在页里出现」、反向查「页里的 `--xxx` 真实存在」。**上线当场抓到 19 处真漏**，而提出它的起因只有 1 处（`--shuffle`） |
| 那个测试的判据 | 一次写对 | 收紧了三次。① 页里写短名（`-a ELEM`）也算写了 —— 要断言的是参数被交代过，不是拼法；② 允许**交叉引用**别的命令的参数（`ferro net --export-traj` 出现在 gr 页里是对的），只抓哪儿都不存在的；③ 跳过「there is no --from / --to flag」这类**否定陈述** —— 那句话正是在做对的事 |
| 共享参数组的未文档化项 | 一律要求写进帮助 | `SelectArgs` 由 `gr` 与 `angle` 共用，`gr` 因此接受 `-c`/`-z` 却对它无意义。**有意不写要进白名单并写明理由**，而不是放宽断言 —— 放宽等于把这个检查关掉 |
| 帮助页文本怎么拿到测试里 | 把 25 个 `print_*` 改成返回 `&str` | 比这个检查本身还大的改动。`include_str!("help.rs")` 按函数名切出函数体即可；函数名写错会当场 panic（找不到那一页），所以映射表自己不会烂 |

## 一处语义改了、四处文档跟不上（2026-08-26）

`collect` 改成「一目录一 system」时只改了 `docs/src/dataset/collect.md`，
事后审出**四处**还写着旧语义或缺新事实：

| 位置 | 漏了什么 | 为什么容易漏 |
|---|---|---|
| `docs/src/cli-reference.md` § collect | 还写「一个输入一个 system 目录，目录名取 stem，撞车改 `<父目录>_<stem>`」 | 手册**两处**讲同一个命令：专页讲「为什么」，参考页讲「参数与语义」。改了专页很容易以为完了 |
| `docs/src/cli-reference.md` § filter | 还写「报告三张表**只打印不落盘**」 | 同上 |
| `CLAUDE.md` 代码导航 | 缺 `doc.rs`、缺 `help_sync` | 新增文件时导航表不在改动路径上 |
| `README.md` / `cli-reference.md` / `dev/overview.md` 的命令树 | 三处都缺 `ferro doc` | **同一张命令树抄在三个文件里** |

**教训**：改命令语义时要 grep 的不只是它的专页。这几处的共同点是「同一份事实
被抄在多个地方」——命令树三份、命令语义两份。防漂测试只盯 clap 与帮助页那一对，
盯不到文档之间。

**当前没有更好的办法**：把命令树做成生成的（像 `supported_formats()` 那样）
需要 markdown 侧有注入点，mdBook 的预处理器能做但要多一个构建步骤。记在这里，
下次再改命令树时**三处一起改**。

## extxyz 的 stress/virial 编码陷阱（2026-08-26）

- **自读自写的回路验证不了符号**。读侧与写侧同时漏掉变号时，往返测试逐位相同、
  全绿。符号的依据必须来自**外部实物**：fixture 由 ASE 3.29.0 生成，期望值由
  dpdata 1.0.2 的换算式算出，两者独立落到同一个数。这与 `array_order.rs` 用
  **非对称**矩阵测行优先是同一类防护 —— 对称张量下转置静默，自洽回路下符号静默
- **`stress` 与 `virial` 是两个量，不是两个名字**。差一个体积因子（virial 是 eV，
  stress 是 eV/Å³），且相对「正 = 压缩」这个内部约定**符号相反**：extxyz 的
  `stress=` 随 ASE（正 = 拉伸，读时变号），`virial=` 随 QUIP/GPUMD/DeePMD
  （正 = 压缩，除体积、不变号）。三个独立来源一致：extxyz 规格写
  `virial -> stress` 是乘 `-1/cell_vol`；GPUMD 手册明写两个键的正负各代表什么；
  dpdata 1.0.2 的 `virials = -volume * stress_ase`
- **行优先还是列优先：不要站队，检查对称性**。ASE 文档说这九个数是 Fortran order
  （列优先），GPUMD 手册把它拼成 `vxx vxy vxz vyx ...`（行优先）—— 同一个格式被两家
  描述成相反的顺序，却从没人报过 bug，因为**规格要求这个张量对称**（"fail if not
  symmetric"），对称下两种读法逐位相同。所以 reader 检查对称性而不是赌一家。
  注意 `Lattice` 不在此列：它的九个数是三个晶格矢量依次排开，两家描述一致，Ferro
  原本就是对的，**别顺手一起改**
- **6 分量 Voigt 一律拒收**。规格允许这个形式，但组件顺序不是普适的：规格与 ASE 是
  `xx yy zz yz xz xy`，VASP 的 `in kB` 行与 GPUMD 自己的 `stress_*.out` 是
  `xx yy zz xy yz zx`（同一个程序在两种文件里用两种顺序），而文件里没有任何字段
  说明产出者。六个独立浮点数下任何排列都自洽，**从数据里检测不出来**，故拒收而非
  猜一种再告警
- **`force` 与 `forces` 两种拼法都要收**。GPUMD 的 `train.xyz` 用单数，ASE 用复数，
  两者都是合法的 `Properties` 列名。只认复数时读 NEP 数据集会**一声不吭丢掉全部
  受力**（`find("forces")` 返回 `None`，不报错）
- **报错消息里也不许用 `as_slice()`**。它是列优先，会把张量打成转置，读日志的人
  据此判断「哪一侧错了」时会被带偏。走 `matrix3_row_major`

## dataset 划分 / --type 编码陷阱（2026-08-26）

- **`with_extension` 会吃掉 `.train`**。输出基名常是 `sysA.train`，
  `with_extension("xyz")` 把它变成 `sysA.xyz` —— 划分后缀被当成扩展名换掉，
  三个部分互相覆盖。拼字符串：`format!("{}.xyz", base.display())`
- **划分不能取尾巴**。按帧序切最后 10% 的话，test 全是轨迹末尾那一段连续状态，
  与训练集分布不同却又高度自相关，两头都不像。成员取自 `shuffle_order`，
  各部分内部再排回帧序（只有成员是随机的，产物仍逐字节可复现）
- **比例向上取到至少 1 帧**。`(n * ratio).round()` 在 n 小时给 0 —— 用户要了
  验证集，命令成功返回，而验证集不存在。0 帧的部分直接不写，非 0 比例则保底 1 帧；
  train 会被挤空时报错并**指名是哪个 system**
- **`--mode by-source` 与帧级划分互斥**。by-source 的全部产物就是「set 边界落在
  system 边界上」，帧级随机划分正好打碎它。同理 by-source + `--type nep`：
  extxyz 没有 set，边界无处安放。两条都在读第一个文件前失败
- **输入已带 `.train` 再划分要拦**。`merge -i 'data/*.train' --test-ratio 0.1`
  看着很合理，实际是把别人已经划好的数据集再划一次，产出 `X.train.test`。
  只看路径就能判，故在读第一条轨迹之前拦
- **`-0.0` 不是 bug 但看着像**。变号之后写出的 `stress="-0.0000000000 ..."` 会让
  读日志的人怀疑符号逻辑。`fmt` 里把 `-0.0` 归一成 `0.0`

## 帮助页精简 / `ferro doc` 小节寻址编码陷阱（2026-08-26）

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| 精简前的范围估计 | 按「页数」估 | 按**超标页数**估。23 页里只有 8 页超 40 行，其余 16 页只需补一行手册指针。「其余 22 页」这个说法把补一行和重写一页算成了同一件事 |
| 「没有手册就得新写」 | 看 `docs/src/` 有没有同名文件 | `cli-reference.md`（863 行）**已按命令分节**覆盖了全部参数。查覆盖要把它算作回落，否则会误判成「要新写三页手册」 |
| 指针指向大页 | `ferro doc convert` → 整本 `cli-reference` | 让 `ferro doc` 支持**按小节寻址**（`Page.section`，取该 `##` 到下一个同级标题）。另一条路是把参考拆成一堆按命令的小文件，那会让「唯一一份完整参考」碎掉 |
| 小节标题写错 | 输出空白 | 找不到就**回落到整页**，另有测试钉住每个 `section` 在其页里真实存在、且解析出的内容多于 5 行 |
| 防漂测试的覆盖面 | 只扫 `help.rs` | **`ferro net` 的帮助页在 `cmd/net.rs`**（`HELP_EXTRA`，紧挨它需要的 argv 剥离逻辑），于是它是唯一没被检查的命令，而它有 10 个选项。测试要显式把那个文件也 `include_str!` 进来 |
| `--P-O=2.4` 这类参数 | 断言它在 clap 上存在 | `net` 的配对参数把元素对写在**参数名**里，clap 建模不了，`main` 在解析前就从 argv 剥离了 —— clap 永远不知道它。判据取那个 `=`（别处的参数从不这么写） |
| 40 行上限 | 当成硬约束，砍参数表凑数 | **参数表（含值域枚举）是唯一必须完整的一段**。不知道 `--functional` 收哪些值就没法敲命令。收尾时 4 页超标，大头都是参数表，不该再砍 |
| `convert` 的行数 | 按总行数判断超标 | 其中 **27 行是 `supported_formats()` 生成的**格式表 —— 那正是这页存在的理由。手写部分只有 32 行 |
| zsh 里量行数 | `for c in "traj gr"; do $F $c; done` | **zsh 默认不对未加引号的变量做分词**（与 bash 相反），`ferro "traj gr"` 成了一个参数，量出来全是 7 行的 clap 报错。要写 `${=c}` |

## network 重构（0.2.1）编码陷阱

| 位置 | 陷阱 | 正确做法 |
|---|---|---|
| `Former.qn` | 当成「总桥氧数」 | 它是**同元素**连接数（P–O–P），文献 `Qⁿ_m` 的 n。总桥数 = `n + Σm`，而桥氧个数是**第三个**量（`n_bo`）：三簇氧算 1 个桥氧但 2 个连接，三者只在无三簇氧时才相等 |
| `bridges_to` | 只统计「恰好两个形成子」的配体 | 统计全部 `nf.len()>=2`，除自己外每个形成子各算一个连接。文献数的是 *connections*（arXiv 2510.13545: "the number of aluminum and boron **connections** with phosphate"），按桥氧记会让 Σm 亏空、无法与 n 相加 |
| `qn_partner` 的列 | 把形成子自身元素也出成 `m_<自己>` 列 | 排除。新口径下 `m_P ≡ qn`，是恒等重复列；文献 `Pⁿ_mAl,xB` 也没给自己留下标位。去掉不影响闭合——`m_P` 由 n 唯一决定，不是独立分组维度 |
| `linkage` 的端点列 | 列名说桥氧、值填连接数 | 端点数字列必须与 `label()` 渲染的数字**同源**，否则同一行里 `P_2` 和数字列会对不上。故改名 `qn_a/qn_b` 并填同元素连接数 |
| `qn` 表 sd 全为 0 | 当成累加器坏了 | 可能是真的：新口径下 P–O–P 骨架若逐帧不变，sd 就是 0。参考轨迹 5 帧全是 93/207/71/1（独立 Python 复算一致），涨落全在 P–O–Al 与配位数一侧。旧口径把刚性与涨落混在一个数里，所以从前 sd 非零 |
| 共边多面体 | 假设网络全部共角 | 传统 Qn 的三条前提是「全四面体、全共角、无三键氧」。共边下一个邻居贡献两个桥氧而数字看着正常，故加 `n_edge_sharing` 检测。参考体系 3589 对全共角、零共边；实测把截断放大到 3.2 Å 会报出 1117 对——**共边非零时先怀疑截断** |
| `AtomType` 的消费 | 从 `label()` 的字符串反解析数字 | 直接读字段。0.2.1 前有五处这样的解析器，散在四个 crate，改标签格式要同步改六处，漏一处即静默错数据（`extract_qn` 对 `"P_3"` 解析失败后 `unwrap_or(0)`，所有 P 变 Q0） |
| `Moments` 的 sd | 用 `Σf² − N·mean²` 求方差 | 用 Welford。恒定序列下两遍求和给出 ~2e-10 的抵消噪声，而 `sd` 列里的假抖动会被读者当成物理。rayon 是归并不是顺序折叠，故需 Chan/Golub/LeVeque 成对合并 |
| `Moments` 的缺席帧 | 只累加出现过的帧 | 缺席帧按 `f = 0` 计入，在 `finish(n_frames)` 里按同一条合并公式补进来。否则 5 帧里只出现 1 次的 bin 会被当成「比例 1.0」 |
| `Accumulator::new` 预置 | 按参数预置每个元素的空条目 | 参数里点名但轨迹里没有的元素必须在 `finalize` 丢掉，否则均值算成 `0.0`，把「不存在」伪装成「平均值为零」（实测 70Z30P00A 的 mean 表里冒出 `Al=0.00`） |
| `write_lammps_dump` 的 type 列 | 在帧循环内按首见顺序重建 | 整条轨迹确定一次。元素只有几种时逐帧碰巧稳定，写位点标签后某帧恰好没有 `P_4`，同一编号就在不同帧指向不同的东西 |
| label 折进 element 列 | 无条件折叠 | 只折 `<element>_…` 形式的标签。CIF/CP2K/QE 的位点名（`O1`、`Fe1`）折进去后读回来就是不存在的元素 `O1`；不合约定的要跳过并计数告警 |
| 按标签跑 g(r) | 以为标注轨迹能直接多帧分析 | 只在单帧成立。`calc_gr` 要求逐类型粒子数守恒，而标签 population 逐帧在变（实测 `P_3`：149/152/150/150/150）。多帧请按元素选 |
| linkage 的规范半边 | 按行求和当作某位点的总参与度 | 每对只存一次（两端按 `(元素, 桥接数, 配位数)` 排序），行和不是总参与度；要算参与度得把 `_a` 与 `_b` 两列都数一遍 |
| `bridges_to` 与三配位配体 | 期待 `Σm == qn` | 三配位配体计入 `qn` 但没有唯一对端，不进 `bridges_to`，故 `Σm ≤ qn`，差额即三配位桥数 |
| `bridging` 与 `cn` | 当成同一个数 | 两个不同的量：`cn` 数截断内**全部**配体，`bridging` 只数连了 ≥2 个形成子的。参考轨迹的 Al 不带非桥氧（`ligand_type` 无 `O_n, Al` 行），故逐原子恰好相等 —— 这是该体系的性质，不是定义，换个体系立刻分叉 |
| 非 Qn 形成子的 label | 以为数字一律是桥接数 | Qn 形成子（默认 `B,P,Si`）取 Qn，其余取**配位数**。`Al_4` 是四配位 Al。约定在分类时写进 `AtomType::Former.qn`，不在渲染时查静态表 —— `--qn` 能在运行时覆盖它 |
| `LinkKey` 的维度 | 只存两端位点状态 | 必须含配体元素。双配体体系（`--Al-O` + `--Al-F`）里 Al-O-P 与 Al-F-P 是两种桥，合并计数无法与实验对应。参考体系只有 O，所以这个洞在 0.2.1 里一直没暴露 |
| `batch::write_all` 的 meta | `table.meta = meta.clone()` 覆盖 | 拼接：公共块 + 该表自述。分析侧用 `Table::meta_line` 写的逐列说明会被覆盖掉；`Table::concat_union` 也要保留各部分 meta（取第一份），否则说明在 `stack` 阶段就丢了 |
| `composition` 的配体行 | 把 `ligand_type` 三行的 `sd` 相加 | 逐帧重新累加。相关项之和的方差不等于方差之和 —— 与 `qn`/`qn_partner` 那次同一个坑。参考数据里恰好相等，因为三个分量有两个方差为 0（常数不改变和的方差），不能当作「相加也对」的证据 |
| 配体 `fraction` 的分母 | 全体配体原子 | 该配体元素的原子数。「桥氧占全部氧」是文献 BO fraction，「占 O+F」不对应任何常用量。单配体体系看不出区别 |
| Qn 标签写成 `Q2` | 以为够短够清楚 | 带元素前缀 `P-Q2`。双 Qn 形成子体系（B+P、Si+Al）里 `Q0`…`Q4` 各出现两遍，跨元素重名没有相邻列能解 |
| 把伙伴分解写进 label | `Q2(1Al)` 这种文献写法 | 不写。`Q^n(mX)` 表达不了「该桥的配体连了三个形成子」，参考数据里 `(qn=2,m_Al=1,m_P=0)`（Σm 亏空 1）与 `(qn=2,m_Al=1,m_P=1)` 会渲染成同一个标签 |
| 单元词汇与原子词汇混用 | 让 `linkage` 也写 `P-Q2` | 桥联连的是**原子**，Qn 命名的是含多个原子的**单元**；导出轨迹更必须用原子形式，否则 dump reader 拆不出 element |
| 空分布表 | 照写只有表头的 CSV | 跳过并打印原因。只有表头的 `network_qn.csv` 读起来像「测了，Qn 是零」，而实情是压根没测 —— 与 `Accumulator` 那条 `Al=0.00` 是同一类错误 |
