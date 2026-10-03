---
name: add-analysis
description: 给 Ferro 新增一种分析方法（ferro traj / map / net / dataset 下的新子命令或新统计量）时使用，覆盖从 ferro-analysis 实现到 CLI 分支、帮助页、手册页、ferro doc 登记的全流程。
---

# 加分析方法

> 2026-10-03 从 `CLAUDE.md`「扩展项目」原文迁入。分层与批处理的约束见
> `CLAUDE.md`「分层铁律」「几条容易违反的」。

1. 在 `ferro-analysis/src/<domain>/` 实现 —— **纯计算，无文件 I/O**
2. 给结果类型 `to_tables() -> Vec<(String, Table)>` 与 `meta_lines() -> Vec<String>`。
   `meta_lines` **只放批内共享的参数**；逐输入才有意义的量走 `[inputs]` 清单，
   否则第一个文件的组成会摆在全局参数区冒充全局事实
3. `ferro-cli/src/cmd/<group>.rs` 加分支：构造参数（**在读第一个文件前**校验）→
   `batch::map_inputs` → `batch::stack` → `batch::write_all`
4. `ferro-cli/src/help.rs` 加帮助并在 `print_overview` 列出。帮助页照同一模板，
   **五段，且只有五段**：一句话用途 + `Parameters:`（完整，含值域枚举与默认值）
   + `Output:`（≤4 行：产物叫什么、落在哪）+ `Examples:`（2~3 条）+
   `Full documentation:  ferro doc <topic>`。目标 ≤30 行，**参数表与命令列表
   不为凑行数砍** —— 不知道收哪些值就没法敲命令。

   判据、口径、为什么这么设计**一律进手册**：帮助页答「怎么敲」，手册答
   「为什么」。2026-09-22 按这条把 25 页推平了一遍（collect 96 → 28 行）
5. `docs/src/analysis/<name>.md` 加手册页 + `SUMMARY.md` 挂上 +
   **`ferro-cli/src/doc.rs` 的 `PAGES` 加一条**（否则第 4 步那行指针指向空）
6. `main.rs` 的 `mod help_sync` 会自动校验帮助页与 clap 一致，不必手动核对；
   有意不写的参数进 `UNDOCUMENTED` 并写明理由

**收尾**：`cargo test`（含 `help_sync`）+ `cargo clippy` 零警告；改了公共 API 时
`cd ferro-python && cargo check`。
