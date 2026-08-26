# CLAUDE.md

> **动手前先读 `dev/`**：`overview.md`（定位与版本）→ `progress.md`（各 crate 现状）
> → `plan.md`（待办）→ `issues.md`（编码陷阱，**改代码前必查**）。
>
> 本文件只放硬约束与导航。三份文档分工：
> **`docs/src/`** = 怎么用（用户手册，mdBook）· **`dev/`** = 为什么这样定 + 现状 + 待办 ·
> **本文件** = 你必须遵守的规则。任何用法问题（CLI 参数、输出列、数据模型）查
> `docs/src/`，不要在这里重述。

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

- `///` / `//!` doc 注释 → **英文**；`//` 内部注释 → **中文**
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
当前 `0.3.2` **未发版**，其中含破坏性改动（`dataset collect` 的产物布局）却按
用户要求走了 patch 位 —— 与 `v0.2.1` 同类例外。清单见 `dev/overview.md`。

## 扩展项目

**加文件格式**
1. `ferro-io/src/readers/<fmt>.rs` 返回 `Result<Trajectory>` + `writers/<fmt>.rs`
2. 从 `readers/mod.rs`、`writers/mod.rs` 导出
3. `ferro-cli/src/io_dispatch.rs` 加格式检测（`read_trajectory` 与 `write_trajectory` **两处**）
4. `ferro-python/src/io.rs` 加包装

**加分析方法**
1. 在 `ferro-analysis/src/<domain>/` 实现 —— **纯计算，无文件 I/O**
2. 给结果类型 `to_tables() -> Vec<(String, Table)>` 与 `meta_lines() -> Vec<String>`。
   `meta_lines` **只放批内共享的参数**；逐输入才有意义的量走 `[inputs]` 清单，
   否则第一个文件的组成会摆在全局参数区冒充全局事实
3. `ferro-cli/src/cmd/<group>.rs` 加分支：构造参数（**在读第一个文件前**校验）→
   `batch::map_inputs` → `batch::stack` → `batch::write_all`
4. `ferro-cli/src/help.rs` 加帮助并在 `print_overview` 列出。帮助页照同一模板：
   一句话用途 + **完整参数表** + 输出布局 + 2~3 个例子 +
   `Full documentation:  ferro doc <topic>`。目标 ≤40 行，**但参数表（含值域
   枚举）不为凑行数砍** —— 不知道收哪些值就没法敲命令
5. `docs/src/analysis/<name>.md` 加手册页 + `SUMMARY.md` 挂上 +
   **`ferro-cli/src/doc.rs` 的 `PAGES` 加一条**（否则第 4 步那行指针指向空）
6. `main.rs` 的 `mod help_sync` 会自动校验帮助页与 clap 一致，不必手动核对；
   有意不写的参数进 `UNDOCUMENTED` 并写明理由

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
| `ferro-analysis/src/md/` | `gr` `sq` `msd` `angle` `vacf` `rotcorr` `vanhove` `cube_density` `cube_sdf` |
| `ferro-analysis/src/network/` | 单文件 `mod.rs`，六张表的统计 |
| `ferro-analysis/src/dft/` | `bader*`、`chg_sdf`（Bader 算法规格见 `dev/bader.md`） |
| `ferro-analysis/src/ml/` | `filter`（帧筛选 + 交叉表）、`geometry`（最小间距、配位、RDF 壳层）、`diagnostics`（只读四表）、`merge`（分组、规范序、打乱） |
| `ferro-cli/src/` | `main.rs` 子命令树 + `mod help_sync`（帮助/clap 防漂测试）、`batch.rs` 多输入驱动（对结果类型泛型）、`cmd/`、`help.rs`、`doc.rs`（`ferro doc`，手册经 `include_str!` 编译进二进制）、`plot.rs` |

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
