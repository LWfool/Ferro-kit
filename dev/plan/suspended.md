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

**这不是死代码，是漏登记的待办**：`docs/src/analysis/cube-jump.md` 开头就写着
「仅库函数，没有 CLI 入口……待 `map jump` 补上」，`SUMMARY.md` 也挂着这一页，
即手册已经对用户承诺了。此前 `plan.md` 里却没有这一条。

补的时候按 `map` 现有形状走：`cmd/map.rs:223 run_grid` 的 `drive()` 已经把
「一输入一 `.cube`」的遍历收好了，`calc_cube_jump` 的签名与
`calc_cube_density` / `calc_cube_radius` 同形，接一个分支 + 一个 `JumpCmd`
参数结构体即可。手册页与 `doc.rs` 的 `PAGES` **都已经在位**（`ferro doc cube-jump`
现在就能打出来），缺的只有 `help.rs` 的帮助页与 `print_map_overview` 里的一行。

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
