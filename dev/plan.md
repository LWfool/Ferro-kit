# 后续计划

## 已完成（归档）

### ferro-analysis：Bader 电荷分析（2026-05-10）
- `ferro-core/src/charge_grid.rs`、`ferro-io/src/readers/chgcar.rs`
- `ferro-analysis/src/dft/`：`bader.rs`、`bader_grid.rs`、`bader_weight.rs`
- `ferro-cli/src/bin/bader.rs`（`fe-bader`）
- 测试：6 个；待优化：梯度上升并行化

### ferro-analysis/network 重构 + 原子类型分类迁移（2026-05-14）
- `ferro-core/src/network_type.rs`：`TypeParams`、`classify_frame`、标签排序函数
  - 形成子标签：`P0`/`Al2`…，氧标签：`Of`/`On_P`/`Ob_P_P`/`X`，修饰子：`Zn_f`/`Zn_t`/`Zn_b`/`X`
- `ferro-structure/src/typing.rs`：`classify_trajectory`、`apply_type_labels`
- `ferro-structure/src/cluster.rs`：`find_clusters`（Union-Find 基于桥氧连接）
- `ferro-analysis/src/network/mod.rs`：简化为单文件，输出 Qn/氧类型/CN 三张分布表
- `ferro-cli/src/bin/network.rs`：`-m Qn`（轨迹统计）/ `-m type`（单帧标注 + 可选导出）

### ferro-analysis/dft：电荷密度团簇 SDF（ChgSDF）（2026-05-15）
- `ferro-analysis/src/dft/chg_sdf.rs`：`calc_chg_sdf` / `ChgSdfParams` / `ChgSdfResult`
  - 从 QE pp.x cube 文件提取子格，Kabsch 对齐后 pull 插值旋转，累加平均
  - 4 个单元测试（均匀密度插值、子格提取、恒等旋转、180°旋转）
- `ferro-cli/src/bin/cube.rs`：`-m chg_sdf` 模式，接受 `--cubes` 替代 `-i`

### ferro-workflow：CP2K 输入文件生成（2026-05-15）
- `ferro-workflow/src/cp2k.rs`：`Cp2kJobBuilder` + 全部枚举类型
  - 支持任务、泛函、基组、色散、SCF、PBC、k 点、截断、MD 参数
  - 混合泛函（PBE0/B3LYP/HSE06）自动生成 `&HF` 块
  - SCAN/r2SCAN 通过 LIBXC 支持
  - `gth_nval()` 80+ 元素赝势价电子数查表
  - 5 个单元测试
- `ferro-cli/src/bin/job.rs`：`-s cp2k` 分支 + 全部 CLI 参数
- `ferro-cli/src/help.rs`：`print_fe_job_overview()` / `print_job_help()` / `print_job_cp2k()`

### 未成对电子估算 + CP2K 基组库 + QE 移植（2026-05-16）
- `ferro-core/src/spin.rs`：`guess_spin`/`assign_oxidation_states`/`parity_min_multiplicity`
  - magmom → 氧化态+Hund → 奇偶下限三级回退；过渡金属高自旋
  - `elements.rs` 新增 `group_number`/`valence_electrons`/`is_transition_metal`
  - `fe-job` 新增 `--charge`/`--multiplicity`/`--auto-spin`
- `ferro-workflow/src/cp2k.rs`：`Cp2kJobBuilder.auto_spin` 经 `guess_spin` 推断多重度/UKS
- `ferro-workflow/src/cp2k_basis_db.rs`：2829 条基组/赝势库（PBE/SCAN/全电子）
  - 由 examples/ 6 文件解析固化；`basis()`/`potential()` 按元素精确匹配
  - 新增 `pob-dzvp`/`pob-tzvp` 全电子选项
- `ferro-workflow/src/qe.rs`：`QeJobBuilder`（QE pw.x，scf/relax/md/...），复用 `guess_spin`
  - `fe-job -s qe` 分支 + `print_job_qe()` 帮助
- 测试 299 全过，clippy 零警告（并修复 3 处历史警告）

### cube_sdf 迁移 + 全项目审计修复（2026-05-16）
- `ferro-core/src/cluster.rs`（新建）：`NetworkGraph`/`build_network_graph`/`connected_components`
  - `cube_sdf::process_frame` 改用之，删除手工邻接 + petgraph BFS；**ferro-analysis 去掉 petgraph 依赖**
  - `find_clusters` 复用 `connected_components`，删除私有 Union-Find
- 文档手册补全：新增 `docs/src/workflow/{spin,job-builders}.md` + Input Generation 章节；cli-reference 补 QE/charge-spin/pob
- 审计修复 9 项：
  - #1 cif.rs 笛卡尔回退本地守卫；#2 删 `GromacsTopologyBuilder` 桩；#3 删死代码 `properties.rs`
  - #4 移除死 workspace 依赖 petgraph；#5 `Frame::unique_elements()` 消除 5 处重复
  - #6 QE writer/builder 跨层重复（架构所限，加交叉引用文档）；#7 `gth_nval` 改由基组库派生
  - #8 经核实为误报（angle.rs 已有 rem_euclid）；#9 修复 vasp/chgcar VASP4 计数 u64/usize 平台溢出 panic
- 测试 304 全过，clippy 零警告

### ferro-python：PyO3 绑定（2026-05-16）
- 独立 workspace + pyo3 0.21 extension-module；旧占位代码（pyo3 0.20 / 已删 `Molecule` API）全部重写
- `lib.rs`/`types.rs`/`io.rs`/`structure.rs`/`analysis.rs` 按 CLAUDE.md 规范
- `pyproject.toml`/`README.md` 修正以适配 maturin；版本随主 workspace bump 至 0.1.7
- 验证：cargo check/clippy 零警告，maturin build link 通过，Python 冒烟测试全过

### ferro-analysis/md：MSD 绘图 + 自扩散系数拟合（2026-05-17）
- `ferro-analysis/src/md/msd.rs`：新增 `MsdParams.fit_range`、`MsdFit`、`MsdResult.fit`
  - `fit_diffusion`：MSD–t 线性最小二乘，D=slope/6（Einstein 3D，仅总 D）+ R²
  - 拟合区间为 MSD 滞后时间轴的**分数**，索引映射 `i=round(f·(n-1))`；非法区间在重负载并行前快速失败
  - `write_msd`：命中拟合时追加 D/R² 文件头块（Å²/fs、cm²/s ×0.1、m²/s ×1e-5）
- `ferro-cli/src/plot.rs`：`plot_msd`（total 主线 + a/b/c 淡线 + 黑色拟合线，图例含 D/R²）
- `ferro-cli/src/bin/traj.rs`：`--fit-range FMIN,FMAX`；给 `--fit-range` 即算 D 并打印（与 `--plot` 无关），`--plot` 时叠加拟合线
- `ferro-cli/src/help.rs` / `CLAUDE.md`：补 `--fit-range` / `--plot`（msd）文档
- 5 个单元测试（精确直线、分数→索引映射、非法区间、calc 填充 fit、write 头块）；测试 309 全过，clippy 零警告
- 设计/计划：`docs/superpowers/specs|plans/2026-05-17-msd-plot-diffusion-fit*`

### 代码审查修复：g(r) r_max / box_builder / Bader（2026-06-26）
- `ferro-analysis/src/md/gr.rs`：`calc_gr` 把有效 r_max clamp 到所有帧最短晶格矢量的一半
  （MIC 上界），防止小盒子（默认 r_max=10.005 > L/2）下 g(r) 尾部被错误压向 0；clamp 值写回 result.params.r_max
- `ferro-analysis/src/dft/bader_grid.rs`：删除 `bader_ongrid` 第一遍被立即丢弃的梯度上升循环（运行时减半，行为不变）
- `ferro-structure/src/box_builder.rs`：私有 `CellList::neighbors` 加 `seen` bin 去重，修复小盒子（n_bins ≤ 2）近邻重复计数
- 新增 3 个测试（gr clamp + 控制、box_builder 小盒子无重复）；测试 312 全过，clippy 零警告
- 遗留：废弃的 `GrParams::with_auto_rmax`（定义但从未被 CLI 调用）已被 calc_gr 内联 clamp 取代，可后续清理

### g(r) / CN / S(q) 重构 + element/label 拆分（2026-08-08 已实施）

**结果**：提交 1 = `caa309b`（0.1.10）；提交 2 = `c5a48b6`…`fc0a3b0`（0.1.11）。
测试 312 → 334 全过，clippy 零警告。逐项落地情况见 `dev/progress.md` 与 `dev/issues.md`。

需求来源：`-a P -b O` 与 `-a O -b P` 输出完全相同（`gr.dat` 均为 `O-P` 列 + `total`），
`-a`/`-b` 顺序对 CN 不起作用。审查中另外发现 4 处独立问题。

#### 语义确认

- `gr` 对称：`A-B` 与 `B-A` 数值逐点相同（g(r) 无方向）
- `cn` 有向：`CN(A→B) = hist/(N_A·steps)`，`-a` = 中心原子，`-b` = 近邻原子
- 未加权 total S(q)：`f_i ≡ 1` 的退化情形，对应不了任何实验探针，
  参考实现 `examples/code2` 中不存在 → 删除
- `.sq` 的 partial 是 g_ij(r) 的傅里叶变换（未加权），散射因子只在求和成 total 时介入

#### 公式审查结论（`calc_gr` 主体正确）

| 量 | 代码 | 判定 |
|---|---|---|
| 异种 g_AB | `hist·V / (4πr²·dr·N_A·N_B·steps)` | ✅ `hist` 为无序对计数，恰等于 `Σ_{i∈A}Σ_{j∈B}` |
| 同种 g_AA | `2·hist·V / (4πr²·dr·N_A·(N_A−1)·steps)` | ✅ 因子 2 + 有限尺寸修正正确 |
| CN(A→B) | `Σ hist/(N_A·steps)`（同种乘 2） | ✅ |

---

### 提交 1 — S(q) Faber-Ziman 权重修复（先做，独立）

| 项 | 内容 |
|---|---|
| 改动 | `sq.rs:209` 内层 → `for ib in ia..elems.len()`，对齐 `examples/code2/dump2sq.c:327` 的上三角 + `n=2.0` |
| 连带 | `canonical_pair_key` 在上三角下恒等，可简化或保留作断言 |
| 测试 | `Σ_{i≤j} w_ij ≈ 1`（XRD 多 q 点 + 中子） |
| 验证 | 拿实际轨迹跑，确认 `total_xrd` 下降约 30% |
| 版本 | 0.1.10 |

---

### 提交 2 — element/label 拆分 + gr/sq 输出重构（0.1.11）

#### A. 拆分规则（按第一个下划线）

新增 `ferro_core` 统一函数：

| 输入 | element | label |
|---|---|---|
| `O` `P` `Zn` | 原样 | `None` |
| `O_b_P_P` | `O` | `Some("O_b_P_P")` |
| `Zn_f` | `Zn` | `Some("Zn_f")` |
| `foo_bar`（前缀非合法元素） | 原样 | `None` + **警告** |

优点：普通轨迹零影响；不存在 `Pb`→铅、`Po`→钋 这类贪婪匹配的静默误判。

**不动** `cp2k.rs:221` / `qe.rs:193` 的 `extract_element` 与 `cif.rs:463` 的
`element_from_label`——它们处理 `Fe1`/`O2` 这类原生命名，字母前缀规则对其正确，
套用下划线规则反而让 `Fe1` 退化成元素 `Fe1`。下划线规则只作用于 LAMMPS dump。

#### B. `lammps_dump.rs`

读取时对 `element` 列应用拆分规则，填 `element` + `label`；**打印一次映射表**
（`O_b_P_P → O`、`Zn_f → Zn` …），其后附一行说明
`-a/-b/-c use element, -x/-y/-z use label`。前缀非法者单独列出告警。

#### C. `gr.rs`

**删除**：`PairStats`、`GrResult.pair_stats`、`WelfordStats`、`GrParams.r_cut`、
`hist[n_pairs]` 总累加器与 total 归一化块、`lib.rs:19` / `md/mod.rs:23` 的 `PairStats` 导出
（`r_cut` 在 `GrParams` 中唯一用途就是 pair_stats；RDF 不是获取键长的手段）。

**新增** `GrParams.group_by: GroupBy::{Element, Label}`。

**配对模型**：6 个无序对 → **n² 个有序对**（3 元素 → 9），每对产出 `gr` + `cn`，
镜像对的 `gr` 数值重复照写（换取"每配对一组可直接提取的列"的规整结构）。

**修复**：
- 排序改 `sort_by(|a,b| (elem_z(a), a.as_str()).cmp(&(elem_z(b), b.as_str())))`
  — 消除同 Z 标签（`O_f`/`O_b_P_P`、`P_0..P_3`）列顺序随 HashSet 迭代序漂移
- r_max clamp 改用最小**面间距**；抽 `ferro_core::Cell::interplanar_spacings()`
  供 `CellList::build`（`angle.rs:79-81`）与 `calc_gr` 共用
- `sorted_keys` 去掉无条件 `push("total")`（否则移除 total 后写出全 0 幽灵列）

**输出** `write_gr(result, path, pair: Option<(&str,&str)>)`，删 `write_cn`：
- `Some(("P","O"))` → 三列 `r[Ang] \t P-O_gr \t P-O_cn`
- `None` → 宽表，按 (Z_center, Z_neighbor) **分组排序**：
  `O-O_gr, O-O_cn, O-P_gr, O-P_cn, O-Zn_gr, O-Zn_cn, P-O_gr, ...`
- 文件头删 `r_cut` 行与 pair stats 块，新增
  `# pairs: <A>-<B>_gr is symmetric, <A>-<B>_cn is directed (center=A, neighbor=B)`

#### D. `angle.rs`

同样加 `group_by`；`:229` 的 `sort_by_key(elem_z)` 改二级排序
（本次改动会让 angle 大量接触伪标签，必须修）。

#### E. `sq.rs`

- 只对 canonical 键做 FT（跳过镜像，省 3 次变换）
- 无 `total` 列
- 每对新增两条加权列 `w_ij(q)·S_ij(q)`，逐点相加恰等于对应 total
- **label 分组时**：partial 按 label 算，加权时 `label → element` 查散射因子。
  加权求和是双线性的，label 层分解可精确重构 element 层，
  故 `total_xrd`/`total_neutron` 与 element 分组结果**完全相同**，只是分解更细
- **导出规则**：
  - 无 `-x/-y`（element 宽表）→ `q, total_xrd, total_neutron` + 6 对 × 3 列 = 21 列
  - 有 `-x/-y`（label 分组）→ 只写 `q, total_xrd, total_neutron` + 指定那对 3 列 = 6 列
- 列顺序：`q[Ang^-1], total_xrd, total_neutron`，随后每对三列组
  `O-P_sq, O-P_xrd, O-P_neutron`（按 (Z,Z) 排）
- 散射因子查表失败 → 打印警告，不报错
- 重写 `sq.rs:381/398` 断言 `sq["total"]` 的两个测试

#### F. `ferro-cli`

- `-a/-b/-c` = element，`-x/-y/-z` = label，**两组互斥**，同时出现报
  `-a/-b/-c and -x/-y/-z are mutually exclusive`（`-x/-y/-z` 短参数当前未被占用）
- `-a`/`-x` = 中心，`-b`/`-y` = 近邻，顺序有意义
- 不再生成 `gr_cn.dat`；`--r-cut` 从 `-m gr` / `-m sq` 移除
- `--plot` 仅在给了配对参数时生效，否则写完数据文件打印提示并跳过（不报错）
- `plot.rs`：`plot_gr` 简化为单对（g(r) 实线左轴 + CN 淡线右轴），删 total 灰线。
  原 `pair_keys` 取自 `gr.keys()` 导致只画低 Z 中心那一向的 bug 随之消失
- `help.rs` 同步

#### G. `ferro-python`

- `gr_pair(traj, a, b, *, by="element", r_max=None, dr=0.01, r_min=0.005)`
  → `{"r", "A-B_gr", "A-B_cn"}`
- `gr_all(traj, ...)` → `r` + 9 对 18 键（恒按 element）
- 删 `r_cut`；无 total

#### H. 测试与文档

回归测试：n² 对键名完整性 · 同 Z 标签列顺序确定性 · `CN(P→O)/CN(O→P)==N_O/N_P` ·
非正交盒子 clamp 到面间距一半 · `Σw≈1` · `Σ w_ij·S_ij == total_xrd` ·
下划线拆分（含非法前缀告警）。

文档：`dev/issues.md`（NPT 体积归一化限制）· `dev/progress.md` · `CLAUDE.md` · `help.rs`。

#### I. 命名约定（已拍板）

- `.sq` 的 raw partial 列名用 `O-P_sq`（非裸 `O-P`），与 `.dat` 的 `_gr`/`_cn` 后缀统一
- `.dat` 宽表按中心元素分组（`O-*` 全部、再 `P-*`、再 `Zn-*`），不按"先全 gr 后全 cn"分区

#### J. 明确不做

- NPT 逐帧体积归一化（后续单独提交）
- `neighbors_of` 每对枚举两次的 2× 浪费（纯性能）
- cp2k / qe / cif 三处 `extract_element` 的统一（见 A 节理由）
- 旧标签格式（`P0` / `Ob` / `On_P`）的读取兼容层——重跑一次 `fe-network -m type`
  即为新格式，加格式猜测反而可能把真元素误判成旧标签

---

### 依赖包全量升级（2026-08-08 已实施）

五次提交，`ebd6147` → `915412b`，版本 0.1.11 → 0.1.12。

| 步骤 | 依赖 | 变化 | 代码改动 | 验证 |
|---|---|---|---|---|
| 1 | ferro-python lockfile | 32 包 | 0 | `cargo check` |
| 2 | `rust_xlsxwriter` | 0.80 → 0.97 | 0 | 实跑 `fe-network -m qn --format xlsx`，三表结构完整 |
| 3 | `rand` | 0.8 → 0.10 | 3 行（`RngExt` / `rng()` / `random_range`） | box_builder 19 测试全过 |
| 4 | `nalgebra` | 0.34 → 0.35 | 0 | g(r)/CN 输出与 0.34 **逐字节一致** |
| 5 | `pyo3` | 0.21 → 0.29 | 1 处 `skip_from_py_object` | 仅 `cargo check`/`clippy`，运行时未验证 |

陷阱记录见 `dev/issues.md`「依赖升级陷阱」。测试 334 全过，clippy 零警告。

---

### NPT 逐帧体积归一化（2026-08-08 已实施，0.1.13）

`gr.rs` / `sq.rs` 的口径由「先时间平均、后归一化/变换」改为对齐 code1/code2 的
「先逐帧归一化/变换、后时间平均」。完整推导、实测差异与「明确不做」见
`dev/issues.md`「NPT 逐帧体积归一化」。

要点：
- 关键恒等式 `ρ_f·g_f` 中 V_f 自行抵消 → 逐帧变换可由 `⟨ρ_f·g_f⟩` 与 `N·⟨1/V⟩`
  两个时间平均量**精确**重构，无须每帧各做一次 FT。分层与 `calc_sq_from_gr` 签名不动
- 实测：NVT 逐字节不变；NPT 下 g(r) 差 3.06e-3（0.0088%）、S(q) 差 1.06e-2（0.82%）
- 测试 334 → 347：`ferro-analysis` 单元测试 +8，`ferro-cli/tests/npt_normalisation.rs` +5
- fixture：`tests/{43Z43P15A_NPT_5,70Z30P00A_NVT_5}.lammpstrj`，各由生产轨迹**等间隔**
  取 5 帧（连续取会丢掉体积跨度，测不出东西）

### box_builder 括号公式（2026-08-08 已实施，0.1.14）

表面是「`parse_formula` 不支持 `Ca3(PO4)2`」，实际阻塞点在上游：`build_box` 只把
**数据库里的** `cd.formula` 喂给解析器，而库中 26 条全是无括号有机溶剂；用户输入的
`compound` 一旦不在库中，`compounds::find` 直接返回 `None` 就报错了，**根本到不了解析器**。
只加括号支持等于写死代码。

故一并做了三件事：
1. `parse_formula` 改栈式解析，支持嵌套 `()` / `[]`，并校验配对
2. 新增 `resolve_component`：先查库（保留人工录入的 `molecular_mass`），
   查不到则把输入当化学式解析、按原子量求和 → 库外化合物可直接建盒
3. 收紧输入校验：空式 / 空组 / 零计数 / 括号错配 / 未知元素一律 `Err`
   （原先空式返回 `Ok(空表)`，把报错推迟到 `build_box` 的 "no atoms to place"）

测试 347 → 357。遗留：无 CLI / Python 入口；不支持水合物点记法 `CuSO4·5H2O`
（水本就该作为独立 component 传入）。

---

### 批处理 + Table 长表 + 单二进制重构（2026-08-09，**已完成，0.2.0**）

锚点 tag：**`v0.1.15`** —— 单文件输入、`.dat` 宽表、writer 位于 ferro-analysis 的
最后一个状态。此后所有输出的文件名、扩展名、列结构全变。

需求：一条命令跑多个轨迹，逐文件独立分析，结果汇总成一份可直接喂 pandas 的表。

```bash
ferro traj gr -i a.lammpstrj b.lammpstrj -a P -b O          # 多文件
ferro traj gr -i 'runs/*/prod.lammpstrj' -a P -b O -o run1  # 自己展开 glob
```

#### 架构

| # | 决定 |
|---|---|
| 1 | `-i` 恒为 `Vec<PathBuf>`（`num_args = 1..`），**单一代码路径**，不按文件数分派；N=1 是 N 的特例 |
| 2 | `ferro-cli/src/batch.rs`，对结果类型**泛型**：`map_inputs<T>(inputs, f) -> (Vec<(PathBuf,T)>, Vec<Failure>)`；只管索引/循环/跳过/汇报，不认识任何分析类型 |
| 3 | 循环体 read → calc → 结果存 `Vec`、轨迹逐条释放；循环结束后统一导出 + 绘图。（结果小、轨迹大，所以攒结果、丢轨迹） |
| 4 | `ferro_core::Table` + `enum Column { Num(Vec<f64>), Text(Vec<String>) }`；结果类型加**固有方法** `to_tables() -> Vec<(String, Table)>`（**不定 trait**，留待第二个消费者出现）；`ferro_io::write_table(&Table, path, Format)` |
| 5 | 分层判据已写进 `CLAUDE.md`：跨中间层共享的类型下沉 `ferro-core`，判据 = 两个以上中间层需要叫出它的名字。`Trajectory`（结构方向）/ `Table`（产物方向）为两个范例 |
| 6 | 不拆「纯搬家」提交，`Table` 迁移 + 批处理 + 长表一次走完 |

#### 输入

| # | 决定 |
|---|---|
| 7 | 自己展开 glob（新增 `glob` crate）；**不支持** `{a,b}`（交给 shell）；参数顺序 + 参数内字典序；跨模式按规范化路径去重；**零匹配报错**（静默跑零个文件是最难查的失败） |
| 8 | 单文件失败 → 提示 + 跳过 + 继续，末尾汇总，**退出码非零**；参数级错误在读第一个文件前快速失败 |

#### 输出

| # | 决定 |
|---|---|
| 9 | **不留逐文件产物**，只出汇总；`-o` 是**后缀**，命名 `<mode>[_<子表>]_<后缀>.csv`；`.dat` 退役 |
| 10 | 格式恒为 csv，**移除 `--format`**（`ferro net` 保留它自己那个 csv/xlsx） |
| 11 | 元信息 = 一行一个输入的表（文件名、帧数、原子数、volume 均值±标准差、生效 r_max、状态/失败原因），压在数据文件顶部的 `#` 块内，**不单独出文件**（已知 pandas `comment="#"` 会丢弃，接受——这些只是辅助参考） |
| 12 | 数值 `{:.6e}` 不变；列名小写无括号（`gr`/`cn`，保证 `df.gr` 可用） |
| 13 | **表结构跟主产物的粒度走** —— 这是 14/16 不对称的判据 |
| 14 | `gr` 长表：`file, r, center, neighbor, gr, cn` |
| 15 | `angle` 长表：`file, angle, end_a, center, end_c, p` |
| 16 | `sq` 宽表：`file, q, total_xrd, total_neutron, <pair>_sq, <pair>_xrd, <pair>_neutron`。主产物是两条 total（一行一个 q），partial 为诊断附属故横向摊开；元素集不同 → **取列并集、缺的留空**（NaN ≠ 插值） |
| 17 | `msd`/`vacf`/`rotcorr`：现有列前面加 `file` 列，其余不动 |
| 18 | `ferro map` 是唯一例外：仍走同一套遍历与失败跳过，但逐文件出 `.cube`、**文件名必须带 stem**，无汇总、无图 |

长表而非宽表的决定性理由：不同轨迹的**元素集可能不同**（Zn-P-O 与 Al-P-O 同批跑），
配对方向做宽表就要取列并集补空洞。`gr`/`angle` 因此长表化；`sq` 的主产物是 total、
partial 只是附属，才保留宽表（并接受列并集留空）。

#### 绘图

| # | 决定 |
|---|---|
| 19 | 一张 PNG 内**分格子图**：`sq` 1×2（XRD \| Neutron）、`msd` 2×2（total \| a \| b \| c）、`gr` 1×2（g(r) \| CN，**删掉双 Y 轴** `plot.rs:63,71,86`）；每格 N 条曲线，**文件是唯一变量维度**，图例取 stem |
| 20 | 必须点名配对才画图；多文件不自动打开查看器 |
| 21 | `--plot` **冻结为自检用途**，不追 matplotlib；正式作图走 Python（长表 + `sns.lineplot(hue="file")` 一行）。任何「加对数轴/误差棒」的需求一律指向 Python |

#### 实施结果

七步全部落地，测试 362 → 396，clippy 零警告：

1. `ferro-core::table` —— `Table` / `Column` / `concat_union`
2. `ferro-io::write_table` —— csv writer（`TableFormat` 枚举留 xlsx 位）
3. `ferro-analysis` —— 各结果类型 `to_tables()` + `meta_lines()`；**7** 个 `write_*` 退役
   （`ferro bader` 的 ACF/BCF/AVF 不动：Henkelman 既定格式，外部工具按它解析）
4. `ferro-cli::batch` —— glob 展开、`map_inputs`、`stack`、`write_all`、`Summary`
5. `-i` 改 `Vec`，各命令接入循环
6. `plot.rs` —— 面板模型（子图 + 多文件曲线），gr 双 Y 轴删除
7. **单二进制 `ferro`**：八个 `fe-*` 删除，子命令按产物分组
   （`traj` / `map` / `net` / `bader` / `convert` / `info` / `job`），三级帮助

`angle` 比原计划多出一列 `count`：整数直方图是 `compare_angle.py` 与 dump2analysis
逐 bin 对拍的依据，只留归一化的 `p` 会废掉这条验证路径。

元信息拆两处：`meta_lines()` 只放批内共享参数，逐文件的量（组成、clamp 后的 `r_max`）
进 `[inputs]` 清单——否则第一个文件的组成会摆在全局参数区冒充全局事实。

#### 已知会破的东西

- `scripts/` 四个 Python 脚本（`compare_rdf.py:69` 是 `[float(x) for x in line.split()]`，
  多一列文本直接 `ValueError`）—— 本次不改；2026-08-11 已提为高优先级，见下
- 所有现存的 `.dat` 产物与依赖它们的个人脚本

### ferro net 重构：结构化类型 + 标签重做 + 合并命令 + linkage（2026-08-12，已完成，0.2.1）

原计划的「net type 标签体系重做」与「net 接入批处理 + 长表化」两项，实施时发现
它们不是两件事的先后，而是同一件事的两层：标签之所以改不动，是因为**分类结果以
字符串形式流转**。六个提交：

| # | 提交 | 内容 |
|---|---|---|
| 1 | `5f8ac0e` | `classify_frame` 返回结构化 `AtomType`，消除**五处**标签反解析（比预估多一处：`ferro-structure/cluster.rs` 的 `starts_with("Ob_")`）。标签文本不动 → 22 个产物文件逐字节对拍一致 |
| 2 | `4503d11` | 标签改 `<元素>_<后缀>`；删修饰子角色分类 |
| 3 | `f1f929d` | `net` 降为叶子命令，接 `CommonArgs` + 批处理 + 长表；删 `--format`/`--frame` |
| 4a | `801184b` | 分布表加 `sd` 列（Welford）；氧的伙伴改数据列；`qn` → `n_bridge` |
| 4b | `304bdd9` | linkage 长表 + `Q^n(mAl)` 分解 |
| 5 | `2f9021e` | 标签存 `atom.label`；extxyz 加 `label:S:1`；dump 折叠 + type 编号跨帧固定 |

**第 1 步单独成一提交是关键**：它把「信息不再经字符串往返」与「标签换格式」分开，
前者可逐字节对拍验证，后者才是行为变更。合到一起就没有能对拍的中间态。

实施中由实测推翻或修正的原计划：

- **`X` 兜底的分裂**：计划只说「≥3 配位氧与 ≥3 NBO 修饰子塌缩」，实测单帧
  `X = 181` 里 **180 个是 Zn、1 个是氧** —— 比预期严重得多
- **修饰子角色分类整体删除**（原计划是保留 `Zn_f/_t/_b` 不变）：实测 97.1% 落进
  兜底桶（`Zn_f` 0 / `Zn_t` 1 / `Zn_b` 26 / `X` 903），该分档对 Zn 无分辨力
- **氧标签不带伙伴后缀**，伙伴改走数据列：标签保持朴素而统计保持完整，两者解耦
- **`qn` → `n_bridge`**：Al 没有 Qn 这一说；实测该体系里 Al 的 `n_bridge` 分布与
  `cn` 分布逐行相同（Al 没有非桥氧），旧表里那三行 Al 是顶着错名的冗余
- **动机比计划写的窄**：`net type 导出 → traj gr -x P_3 -y O_b 直接串起来`
  只在**单帧**成立，多帧被 `calc_gr` 的粒子数守恒守卫拒绝（标签 population 逐帧在变）
- 顺带修：`Accumulator` 预置空条目导致不存在的元素被写成「平均值 0」；
  `write_lammps_dump` 的 type 编号逐帧重建
- 新增 linkage 表与 `m_<X>` 分解（原计划无，源自对 `private/coord_analysis.py`
  的审查与用户对 Al-O-Al 的研究需求）

测试 398 → 420，clippy 零警告。产物交叉对账全部闭合（见提交 `304bdd9` 说明）。

### scripts/ 四个对拍脚本修复（2026-08-11 已实施）

0.2.0 的输出变更（`.dat` 退役 → csv、宽表 → 长表、新增 `file` 文本列、顶部 `#` 元信息块）
让四个对拍脚本全部失效。典型断点：`compare_rdf.py` 的
`[float(x) for x in line.split()]` 遇到 `file` 文本列直接 `ValueError`。

**新增 `scripts/ferrocmp.py`** 收拢共用逻辑 —— 四个脚本各有一份几乎相同的
`run` / `load_columns` / `fe_version` / 列名行解析，形状一变就要改四遍：

| 函数 | 作用 |
|---|---|
| `run_ferro(bin, argv, outdir, product)` | **以 outdir 为工作目录**跑 ferro（`-o` 是后缀不是路径），跑前删同名残留，查退出码 + 查产物 |
| `read` / `one_file` / `pick` / `column` | `read_csv(comment="#")`；断言单输入并丢 `file` 列；长表按数据列选行；宽表按列名取列 |
| `meta_lines` / `version` | `#` 块（mean/std 这类只在注释里的量仍需解析文本） |
| `load_legacy` / `legacy_header` / `same_grid` | 参考程序侧（格式未变）+ 横轴对齐校验 |
| `run` | C 程序专用：dump2analysis / dump2sq 报错也 exit(0)，退出码不可信，只能验产物 |

要点：

- `-o` 是**后缀不是路径**是本次最大的形状变化：不改工作目录的话产物会散落在
  脚本的启动目录里，而不是 `--outdir`
- 长表取一条曲线是 `pick(df, path, "r", center="P", neighbor="O")`；选空时把该表
  实际有哪些组合打出来（元素名拼错时最需要看到的就是这个）
- `one_file` 显式拒绝多输入：静默取第一个会让后面所有数字都对不上，且图上看不出来
- `--r-min` 自 0.1.15 已上 CLI，两侧显式传 0.005 对齐，旧脚本里那条
  「fe-traj 的 r_min 固定 0.005，CLI 未暴露」的不对称说明随之删除
- `--fe-traj` 参数改名 `--ferro`；`trajcheck.py` 里的 `fe-traj` 字样一并更新

复跑验证（数值与 `dev/issues.md` 记录逐项吻合，说明改的只是读法不是口径）：

| 脚本 | 结果 |
|---|---|
| `compare_rdf.py` | P-O / Al-O 的 g(r)、CN(r) 峰值与峰位完全相同，max\|Δ\| 5e-5 / 5e-4 —— 仍是 dump2analysis `%12g` 的量化台阶 |
| `compare_angle.py` | Σfe/Σref = 0.5000（×2 计数约定）；`--align-binning` 下 1800 个 bin **整数零差** |
| `compare_sq.py` | 未补偿残差 q>15 均值 1.00031 / 1.00018（+1 常数）；partial 互相关 +0.959（type_new 失效） |
| `compare_sq_experiment.py` | 50 帧 rms(fe−exp) 0.0194/0.0277、max\|fe−ref\| 0.0027/0.0008、FSDP 1.95/2.05；5 帧回退 0.0197/0.0282 |

顺带修：`compare_sq_experiment.py` 的 `TRAJ_CANDIDATES` 首选指向
`tests/70Z30P00A_NVT.lammpstrj`，该文件**不存在**（50 帧的那条在 `examples/` 下），
一直在静默回退到 5 帧子集。改为 `examples/70Z30P00A_NVT.lammpstrj`。

---

## 优先级高

> 2026-08-12 更新：原「net type 标签重做」与「net 接入批处理」两项**已完成**，
> 归档见下方。剩下 `chg-sdf` 一项。

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
  这正是 REPL 相对 shell 循环的全部价值。这条也是 0.2.0 决定合并成单二进制的理由——
  子进程分发做不到状态保持，而链接进来之后再包一层 `fe-*` 前端就是重复
- 脚本语法应建在**已经定型**的子命令树上（`gr -a P -b O` 直接复用 `cmd::traj` 的
  参数结构体），不要另起一套
- 注意与「批处理输入」是两件事：这里是**命令**的批处理（一个脚本跑多条命令），
  那里是**输入文件**的批处理（一条命令跑多个轨迹）。两者可叠加但互不依赖
- 降为中级的理由：子命令树刚在 0.2.0 定型，脚本语法要建在它上面；先让上面四项
  高优先级把 `net` / `scripts` / 标签体系的形状敲定，避免 REPL 复用到一半又改

### ferro-python（2026-08-11 全部归为中级）

两件事，第一件是第二件的前提。

**1. ~~修复编译断裂~~ —— 2026-08-12 复核：不存在。** 计划此前写着
「`analysis.rs` 仍调用 0.2.0 删除的 `write_gr` 路径，`cargo check` 直接失败」，
实测 `cargo clean && cargo check` **干净通过、clippy 零警告**，`gr_pair` /
`gr_all` / `msd` 三个 `#[pyfunction]` 都在，没有任何 `write_*` 调用。
该条目应是 0.2.0 重构当时的推测被写成了事实，未经复核。

留下的真问题仍是**独立 workspace 导致主 workspace 的 `cargo build/test/clippy`
全部跳过它**：编译断裂不会被发现，只能靠每次改 `ferro-core` / `ferro-analysis`
公共 API 后补跑 `cd ferro-python && cargo check`。0.2.1 的 net 重构改了这两个
crate 的公共 API，已复核通过。

**2. pyo3 0.29 运行时验证。** 0.21 → 0.29 只做了类型层验证。本机无 maturin，
`cargo build` 在 macOS link 阶段过不了（`extension-module` 需 `-undefined dynamic_lookup`）：

```bash
pip install maturin
cd ferro-python && maturin build --interpreter "$(which python)"
pip install target/wheels/*.whl
```

冒烟测试要覆盖：读 xyz/cif/lammpstrj、`supercell`、`write`、
**新拆的 `gr_pair` / `gr_all`（含 `by="label"`）**、`msd`。

### 搁置项（待实测后确定）

- **`vanhove` 加 `tau` 列**：现在一次只算一个 τ、写在 `#` 头里。加列后将来支持多 τ 是
  加行而非改列结构

### 是否采用 rustfmt（待定，暂缓，2026-08-09 由高降中）

收益：把格式从代码审查面里移除。代价三点：

1. 一次性约 80 文件的大 diff，`git blame` 在这些行上全部指向那一次提交
2. 会拆掉现有的刻意列对齐（`gr.rs` 的 `WelfordStats`、`angle.rs` 的 `CellList`、
   `units.rs` 的枚举表）
3. 生成文件要显式排除：`rustfmt.toml` 的 `ignore` 仅 nightly 可用，stable 上需给
   `cp2k_basis_db.rs` 的静态表加 `#[rustfmt::skip]`

若采用：单独一次纯格式提交（`style: adopt rustfmt`）+ 固定 `rustfmt.toml`，
不要混进特性或依赖升级提交。

### ferro-structure：补充结构操作
- `rotate.rs` / `orient.rs`：旋转和晶体学取向
- `substitute.rs`：元素替换
- `disturb.rs`：随机位移（热扰动初始构型）
- `select.rs`：原子选择器

#### Multiwfn 结构操作参考（examples/Multiwfn_2026.4.10_src_Linux）
核心模块 `otherfunc3.f90` 的 `geom_operation`（行 1975–3066，约 1090 行）为
「Geometry relevant operations」主菜单 + 分派；辅以 `displace_geom`(1879)、
`extract_molclus`(3066) 及 `PBC.f90`（晶胞数学层）、`util.f90:985 str2arr`
（原子选择串 `"2,3,7-10"` 解析）。

`geom_operation` 操作目录：
- 平移/旋转/取向：1 按向量平移；2 平移使选中部分质心到原点；
  3 绕笛卡尔轴/键旋转（子菜单：X/Y/Z/绕键/绕指定向量）；4 按旋转矩阵旋转；
  5 使键平行于向量/轴；6 使向量平行于…；7 使偶极矩平行于…；
  8 使最长轴平行于…；9 镜像反演；10 中心反演；11 使选中原子平面平行于笛卡尔平面
- 编辑：12 缩放坐标(含晶胞)；13 重排原子序；15 加原子；16 删原子；17 裁剪原子
- 周期/晶胞：18 随机位移(displace_geom)；19 复制构建超胞(nrepli1/2/3，支持负方向)；
  20 边界截断分子补全(makemolwhole)；21 缩放晶胞+坐标；22 原子折叠入胞；
  23 沿晶轴平移；24 平移使选中部分居中；25 提取团簇；26 设晶胞；27 加边界原子；
  28 坐标轴互换
- 输出：-1 xyz / -2 pdb / -3 gjf / -4 cif；-9 恢复原始；-10 返回；0 可视化
- 子菜单：中心可选 几何中心/质心/核电荷中心；displace_geom 支持选原子(str2arr)、
  选坐标分量(X/Y/Z 及组合)、标准差(默认 0.03 Å)、生成数量

移植对应：
| ferro-structure 待办 | Multiwfn 参考 |
|---|---|
| `rotate.rs` | 菜单 3、4（轴/键/矩阵旋转） |
| `orient.rs` | 菜单 5/6/8/11（对齐键/向量/最长轴/平面） |
| `disturb.rs` | 菜单 18 / `displace_geom`（高斯随机位移） |
| `select.rs` | `util.f90:str2arr`（`"2,3,7-10"` 选择语法） |
| `substitute.rs` | 无直接对应（Multiwfn 仅 15/16 加删原子，元素替换需自行扩展） |
| supercell | 菜单 19（ferro 已实现） |

### ferro-workflow：扩展 QC 软件支持
- VASP（POSCAR/INCAR/KPOINTS）

---

## 优先级低

### 机器学习集成
- 阶段一：`linfa`（K-Means/DBSCAN/PCA），集成至 `ferro-analysis/src/ml/`
- 阶段二：`candle`（ONNX 推理，加载 DeepMD-kit 势函数）
- 阶段三：`burn`（纯 Rust 训练，远期）

### 深度学习工作流 I/O
- `ferro-io` 新增 DeePMD-kit 训练集格式、MACE/NequIP 兼容格式
