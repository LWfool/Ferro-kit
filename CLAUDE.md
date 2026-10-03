# CLAUDE.md

> **动手前先读 `dev/`**：`overview.md`（定位与版本）→ `progress.md`（各 crate 现状）
> → `plan.md`（待办）→ `issues.md`（编码陷阱，**改代码前必查**）。
>
> 本文件只放硬约束与导航。三份文档分工：
> **`docs/src/`** = 怎么用（用户手册，mdBook）· **`dev/`** = 为什么这样定 + 现状 + 待办 ·
> **本文件** = 你必须遵守的规则。任何用法问题（CLI 参数、输出列、数据模型）查
> `docs/src/`，不要在这里重述。
>
> 问到进度或下一步时，读 `dev/progress.md` / `dev/plan.md`，**不要凭记忆作答**。

## Build & Test

```bash
cargo build                                   # 整个 workspace
cargo test                                    # 全部测试
cargo test --package ferro-core test_name     # 单个测试
cargo clippy                                  # 必须零警告
cd ferro-python && cargo check                # 独立 workspace，主 workspace 会跳过它
```

`cargo fmt` **不要全仓跑**：本仓未纳入 rustfmt 管理，全仓格式化会改动约 80 文件
（含 22k 行生成代码），淹没真实 diff。见 `dev/issues.md`「构建 / 工具链陷阱」。

测试 fixture：`tests/`（两条 5 帧 LAMMPS 轨迹，一 NPT 一 NVT）。

## 分层铁律

```
ferro-cli / ferro-python        ← 唯一允许组合多个 crate 的层
    ├── ferro-core                纯数据结构 + 静态数据 + Table + AtomType
    ├── ferro-io        → core    格式读写；write_table 是分析产物的唯一出口
    ├── ferro-structure → core    超胞、真空层、合并、建盒
    ├── ferro-analysis  → core    纯计算，**不碰文件系统**
    └── ferro-workflow  → core    QC 输入生成
```

**中间层 crate 不得互相依赖。** 共享类型下沉而非横向依赖：

> 一个类型该放 `ferro-core`，**当且仅当两个以上中间层需要叫出它的名字**。

两个方向各一个范例：`Trajectory`（io 产出 → analysis 消费）、`Table`（analysis 产出
→ io 消费）。分析私有的中间产物（`GrResult`、`SqResult`…）留在 `ferro-analysis`，
只有可序列化投影（`to_tables()` 出的 `Table`）跨层。

## 编码约定

- **手册（`docs/src/`）→ 英文**；`//` 内部注释与测试断言消息 → **中文** ——
  它们承载判据，母语更准，读者是开发者不是用户。`///` / `//!` doc 注释
  **中英文均可**（2026-09-28 用户裁定；现存约 550 行中文 doc 注释不再算违规）
- 手册的公式与化学式**一律 LaTeX `$...$`**，不直接写 Unicode 上下标
  （写 `$P_2O_5$`，不写 `P₂O₅`）。**三处不动**：① 单位符号（`Å³`、`Å⁻¹`、
  `g/cm³`、`cm²/s`）—— 它们是符号不是公式，进 LaTeX 只会让表格更难读，
  手册原本就是「公式 LaTeX、单位 Unicode」；② 代码块与行内代码；
  ③ `SUMMARY.md` 的目录标题（mdBook 侧边栏不过 MathJax，会字面显示）
- 库 crate 用 `ferro_core::error::ChemError` / `Result<T>`；CLI 用 `anyhow::Result`
- 顶层类型恒为 `Trajectory`，单帧文件也是（`frames: vec![frame_0]`）
- **`Molecule` 类型不存在** —— `Frame` 覆盖分子与周期体系，由 `pbc: [bool; 3]` 区分
- 原子索引是**隐式的**（在 `Vec<Atom>` 里的位置），不存 `index` 字段
- **矩阵一律行优先**：`ferro-core` 里行 = 晶格矢量 / 张量行；落盘取九个数走
  `matrix3_row_major()`，**禁止对 nalgebra 类型用 `as_slice()`**（它是列优先，
  即转置；对称张量下这个错误完全静默）

内部单位（DeePMD-kit / VASP 约定）：长度 Å · 能量 eV · 力 eV/Å · 应力 eV/Å³ ·
时间 fs · 质量 amu · 电荷 e · 温度 K。转换走 `units.rs` 的枚举，不引入 `uom`。

## 版本规则

版本号在根 `Cargo.toml` 集中管理，各 crate 继承 `workspace.package.version`；
`ferro-python` 是独立 workspace，需**手动同步**。

**只在用户明确要求时才动版本号**，不要每次改代码就自动 +1。
当前版本见 `dev/progress.md` 顶部，历次版本批次与破坏性改动（含走 patch 位的例外）见
`dev/overview.md`。何时升版本号的细则见 `dev/memory.md`「版本号更新规则」。

**`-o` 恒为路径**：写多个产物的命令指目录，`convert` / `job` 指文件；批次标记走
`-s/--suffix`。目录缺失时经 `outpath.rs` 询问，非交互环境须给 `--mkdir`。

## 扩展项目

**加文件格式** → 技能 `add-format`。硬约束：`ferro-cli/src/io_dispatch.rs` 的
`read_trajectory` 与 `write_trajectory` **两处**都要加检测，漏一处那个方向就认不出新格式。

**加分析方法** → 技能 `add-analysis`。三条违反了会静默出错的：
- `meta_lines` **只放批内共享的参数**；逐输入才有意义的量走 `[inputs]` 清单，
  否则第一个文件的组成会摆在全局参数区冒充全局事实
- 帮助页**五段，且只有五段**（用途 / `Parameters:` / `Output:` / `Examples:` /
  `Full documentation:`）；判据与「为什么」一律进手册
- 加手册页后 **`ferro-cli/src/doc.rs` 的 `PAGES` 加一条**，否则帮助页末行的
  `ferro doc <topic>` 指向空

**加结构操作**：`ferro-structure/src/` 收发 `Trajectory` → `ferro-python/src/structure.rs`。
目前无 CLI 入口（`box_builder` 也是库级）。

**加 QC 目标**：`ferro-workflow/src/job_builder.rs` + `templates.rs` → `cmd/job.rs` 分支。

**加数据集判据**：`ferro-analysis/src/ml/filter.rs` 的 `Criterion` 加变体（记得
`ALL`、`name`、`bit`、`enabled` 四处）→ `FilterParams` 加参数 → `cmd/dataset.rs`
接 CLI。判据一律表达成「符合条件则**删**」，四条同一语义交叉表才不用分裂。

## 代码导航

| 位置 | 内容 |
|---|---|
| `ferro-core/src/` | `atom/cell/frame/trajectory`、`table.rs`、`network_type.rs`（`AtomType`）、`cluster.rs`、`spin.rs`、`data/`（元素、化合物、Qn 名单） |
| `ferro-analysis/src/md/` | `gr` `sq` `msd` `angle` `vacf` `rotcorr` `vanhove` `cube_density` `cube_radius` `cube_jump`（**无 CLI 入口**，见 `dev/plan.md`）`cube_sdf` `scattering_data` |
| `ferro-analysis/src/network/` | 单文件 `mod.rs`，六张表的统计 |
| `ferro-analysis/src/dft/` | `bader*`、`chg_sdf`（Bader 算法规格见 `dev/bader.md`） |
| `ferro-analysis/src/ml/` | `filter`（帧筛选 + 交叉表）、`geometry`（最小间距、配位、RDF 壳层）、`diagnostics`（只读四表）、`merge`（分组、规范序、打乱） |
| `ferro-cli/src/` | `main.rs` 子命令树 + `mod help_sync`（帮助/clap 防漂测试）、`batch.rs` 多输入驱动（对结果类型泛型）、`outpath.rs`（`-o` 的目录创建与确认）、`cmd/`、`help.rs`、`doc.rs`（`ferro doc`，手册经 `include_str!` 编译进二进制）、`plot.rs` |

**几条容易违反的**（完整清单在 `dev/issues.md`）：

- 绝不从 `AtomType::label()` 的字符串反解析数字 —— 直接读字段。0.2.1 重构消灭的
  正是散在四个 crate 的五处这种解析
- 列并集下缺失值填 `f64::NAN`（渲染为空字段），**绝不补零** —— 零是「测到了 0」
- 批内单文件失败：跳过 + 记入 `[inputs]` + **退出码 1**；参数错误在读第一个文件前快速失败
- **`ferro net` 的帮助页在 `cmd/net.rs` 的 `HELP_EXTRA`**，不在 `help.rs` —— 它紧挨
  着自己需要的 argv 剥离逻辑。改 net 的参数时别只看 `help.rs`
- 帮助页与 clap 是**两处手写同一份事实**，靠 `help_sync` 测试盯着。别放宽它的断言
  来让测试变绿：有意不写的进白名单并写明理由
- 不按输入文件数分派单/批两条路径，N=1 是 N 的特例

## 四类信息放在哪

项目的规则、记忆、技能**全部放在本项目内**，不写进用户级的 `~/.claude/` 目录。

| 类别 | 位置 | 放什么 | 什么时候读 |
|---|---|---|---|
| 规则 | 本文件 | 任何时候都要守、违反就出错的约束 | 每次会话自动加载 |
| 状态 / 计划 | `dev/progress.md`、`dev/overview.md` / `dev/plan.md` | 已经完成了什么、版本批次；接下来按什么优先级做 | 问到进度时 |
| 记忆 | `dev/memory.md` | 被纠正过的协作做法、容易答错的工具语义 | 精简文档、动版本号、提交之前 |
| 细节 | `dev/issues.md`、`dev/bader.md`、`dev/glossary.md` | 编码陷阱、算法规格、译名 | 改代码前 |
| 技能 | `.claude/skills/*/SKILL.md` | 反复做、步骤固定的操作 | 做对应操作时 |

新信息按上表归位：**进度变了改 `dev/progress.md`，不改本文件**；被用户纠正过一次的
做法，记进 `dev/memory.md`。用法类内容仍进 `docs/src/`。

**技能清单：**

| 技能 | 用途 |
|---|---|
| `add-analysis` | 新增分析方法：实现 → CLI 分支 → 帮助页 → 手册页 → `PAGES` 登记 |
| `add-format` | 新增文件格式：reader/writer → 导出 → `io_dispatch` 两处 → Python 包装 |
| `update-dev-records` | 「更新开发记录」或「发版」：定版本号 → 更新 dev/ → 同步 ferro-python → 打 tag |
