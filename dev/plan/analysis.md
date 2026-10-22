# 待办：ferro-analysis（MD 分析、数据集、网络）

> 从 `dev/plan.md` 按需读取；定级规则见 `plan.md` 顶部。

## 优先级中

### vanhove：默认 1 个原点 + 仍用格点视图解包裹（2026-09-27 发现）

固定 lag 的位移直方图，FFT 无用。~~`tau` 默认 `n_steps − 1` → 只有 1 个原点~~（2026-10-08 随 B-5 改为 N/2）；解包裹
仍是 `unwrap_frac`（分数坐标解包裹 × 当帧盒子，即格点视图，NPT 下有 msd 手册里写的
问题）。改默认 tau（N/2？）并换 TOR —— `msd.rs::unwrap_tor` 届时有第二个用户，
可下沉共享。

### bondlife 的后续（2026-09-27 记）

- **氢键几何判据**（供体–H–受体角度）：水–玻璃界面的 H 键需要它；现只有距离判据
- **cell list**：候选搜索是逐帧 中心数 × 邻居数 全扫；上万原子 × 上千帧会慢
- **Luzar–Chandler 速率常数** k、k′：需另定义 n(t)（键断但仍在第二阈值内的对），
  gmx hbond-legacy 用 `-r2`
- Yamamoto–Onuki 严格的端点定义（t₀ 时 r ≤ A₁、t₀+Δt 时 r ≤ A₂，不看中间）与现行
  滞后序列只在 r_break > r_bond 时略有差别；如需逐字复现可加 FFT 互相关

2026-10-10 对照业内做法补记三条：

- **`tau_int_integral` 随 `max_lag` 不收敛**：间歇 $C_I$ 有平台（约 $\sum p^2/\sum p$，键在窗内
  反复回来），尾部按 $t^{-3/2}$ 走，积分值取决于截到哪。手册应写明，或只报 1/e 交点 /
  多指数拟合积分（Schulze & Spiekermann 2023 的做法：去掉前 ~100 fs 后积分，未衰减完就
  2–3 指数拟合后积分）—— 待定
- **事件的最短持续过滤**：截断附近来回进出的 O 产生大量假事件。配位态原型（下一节）实测
  单截断比双截断核心态的计数速率虚高 44%（合成轨迹，已知真值）
- **按近邻集合的停留分析**：「近邻集合不变」判据（集合身份，不只是 CN 不变），与下一节
  共用实现

### 配位态动力学：逐中心 CN(t)、核心态转变、集合身份自相关（2026-10-10 提出，中；未落地）

用户研究 Al5→Al6（ZnO–P₂O₅–Al₂O₃ 玻璃，DPMD）。bondlife 只看原子对，缺「逐中心原子的近邻
集合时间序列」。Python 原型在仓库外（用户计算目录 `…/DC2/al_cn_kinetics/al_cn_kinetics.py`
+ `test_synthetic.py`），落地须按 `add-analysis` 技能走，且先问用户要不要进 Ferro。

**判据**（原型已实现并经合成轨迹验证）：

- 核心态：缓冲带 $[r_c-\delta, r_c+\delta)$ 内无配体、且 $r_c-\delta$ 以内恰 n 个 → 态 n；否则
  沿用上一个核心态（milestoning / 基于转变的划分，Schütte 2011）。只有进入另一核心态才计转变
- Schulze & Spiekermann 2023（GRL）的间歇配位自相关：$s_i=1$ 当 $t_0$ 时 CN=n 且 $t_0+t$ 时
  近邻**集合**与 $t_0$ 相同；集合比较用 Zobrist 异或哈希。该文指出连续型「寿命分箱」低一个
  数量级且排序会反
- 转变概率 $P_{ij}(\tau)=P(m(t+\tau)=j \mid m(t)=i)$；可逆（同一配体进出，访问前后集合相同）
  与净交换拆分；逐中心转变次数分布（不均匀性）
- $r_c$ 取 g(r) 第一极小；谷底平坦时 argmin 位置是噪声，跨输入取中位数统一

**被实测推翻的原计划**：

- 原以「$k_{56}=N/T$ + Arrhenius」为主结果。真实数据上 $k$ 随 δ 单调降、**无平台**（Al–O
  自由能垒仅 2–4 kT）；转移矩阵 `logm` 的滞后扫描 1–150 ps 也无平台 → CN 不是马尔可夫变量，
  单一速率常数不存在。改以 $\tau\ge10$ ps 的 $P_{ij}(\tau)$ 与净交换率为主（后者对 $r_c$、δ
  半宽 ±17–33%，而温度趋势 ×8–13）
- 多指数外推在刚过 1/e 的曲线上给出 $10^6$ ps 量级：只在 $C(t_{end})<1/e$ 时报寿命，且拟合 τ
  超过 2 倍窗口不报

**已对拍**：原型 lammpstrj 读取（`x/xu/xs` 列、NPT/NVT）的 Al 配位计数与 `ferro net --Al-O=2.5`
逐项一致（examples 的 NPT 文件 + 一条 500 帧 NVT）。

### VDOS：振动态密度（2026-09-27 提出，中）

VACF 的傅里叶变换，业内常用（如 LAMMPS/MDAnalysis 生态的 power spectrum）。等 vacf
改 FFT 后顺手做：窗函数（Hann？）、归一化、频率单位（THz / cm⁻¹）需先查规范。

### 共享头部冒充全批：angle / vanhove（2026-09-27 扫描；vacf、rotcorr 已修）

与 msd 同类（`fa2727c` 已修 msd）：`meta_lines()` 里用了逐文件的量，多文件时第一个
文件的值摆在共享区。

| 分析 | 写进共享区的逐文件量 |
|---|---|
| `angle` | `elements` 列表、`[statistics]` 段（每个三元组的 mean/std/count） |
| `vanhove` | `elements`、`n_atoms`、`n_origins`、`tau_frames`、`time` |

修法同 msd：共享区只留参数，其余走 `Summary` 列或 `note()`。`angle` 的统计是逐文件
逐三元组的二维量，放 `[inputs]` 会很宽，可能要单独一张表。

### S(q) 的散射数据从标签字符串反推元素（2026-10-22 发现，中）

`md/sq.rs` 的 `elem_z_local` 用 `symbol_to_z` 从分组键取 Z，再查形状因子与中子散射长度。
g(r) 按 `GroupBy::Label` 分组时键是 site 标签，`symbol_to_z` 先整串精确匹配：标签 `Os`（表面 O）
会得到锇（76），`Ca`、`Co`、`Cs`、`Nb` 同理，S(q) 静默用错数据、不出告警。定中不定高：CLI 的 sq
固定按元素分组（`cmd/traj.rs`），Python 不暴露 sq，目前只有库函数 `calc_sq_from_gr` 可达。

修法：`GrResult` 记下每个分组键对应的元素（来自 `Site.element`），sq 按它查散射数据，不再解析字符串。
`canonical_cmp` 的排序也从标签反推 Z，只影响列顺序，可一并改成读字段。先写测试：标签 `Os` 的 O 原子，
按标签与按元素分组的 S(q) 应一致（现有 `totals_by_group` 用的标签恰好不撞元素符号，没测出来）。
判据背景见 `decisions.md`「手册迁出」gr 小节。

### D 的统计误差：GLS / 贝叶斯（2026-09-27 提出，低）

`d_err` 是 gmx msd 的半窗差，本质是线性度检查。严格误差需要 MSD 协方差，见 kinisi
（PMC11736684）。先解决上面「只有 1 个原点」，否则协方差无从谈起。

### 搁置项

- **`vanhove` 加 `tau` 列**：现在一次只算一个 τ、写在 `#` 头里。加列后将来支持多 τ 是
  加行而非改列结构


## 优先级低

### ferro dataset filter：多层周期镜像（2026-08-25 提出，2026-09-25 降低）

几何判据（`--oo-min`、`-r` Al6）现在把最小镜像当上界（超出报错），没做
`private/filter_dataset.py` 的多层扫描（`n = floor(rcut / w + 0.5)`）。当前体系
盒子约 15 Å、截断约 2.4 Å，按该公式 n = 0，两种做法逐位相同；只有小胞体系才需要。
Al6 走 `ferro-core` 的 `classify_frame`，要做多层得动 net 共用的近邻搜索，面更大。

上界检查原先只看第 0 帧，NPT 缩盒时会漏判；2026-09-25 改为逐帧取最紧的一帧
（`5833ff6`），报错点名帧号。这一半已做完，剩下的只是多层扫描本身。

### 机器学习集成

- 阶段一：`linfa`（K-Means/DBSCAN/PCA）→ `ferro-analysis/src/ml/`
- 阶段二：`candle`（ONNX 推理，加载 DeepMD-kit 势函数）
- 阶段三：`burn`（纯 Rust 训练，远期）

## 2026-10-21 从各审查节归拢的条目（定级沿用原节，见各条）

### 其他 analysis 待办（2026-10-21 从各审查节归拢）

- **S(q) 的 Lorch 窗**（10-19 用户定 1：归新功能，低；做的时候再查文献）
- **分析侧逐帧消费**（10-17 第 4 步用户定「这次不做」，记作后续项）：gr / sq / map / net 边读边算、内存与帧数无关
- **bondlife 手册口径**（10-03 物理 / 口径，未改）：$S_C$ 是「比值的和」（Luzar–Chandler / gmx 同构），MDA `autocorrelation` 是「比值的平均」，lag 1 差 3e-4。手册 `bondlife.md` 应改说「与 MDA 只差平均顺序」

### 库层参数无 `validate`（10-20 扫描 C-2，登记不修）

**登记不修**：**C-2** `ml/filter.rs` 的 `FilterParams` 与 `network/classify.rs` 的 `TypeParams` 没有 `validate`，范围检查只在
CLI（`cmd/dataset/filter.rs:92-117`、`cmd/net.rs:430`）。库调用方传 NaN 时：`f_max` / `s_max` / `oo_min` 的判据静默关闭，
`al6_rcut = NaN` 删光全部帧，`TypeParams` 截断为 NaN 则零成键。**有 Python 入口时补，走 `check.rs`**。
