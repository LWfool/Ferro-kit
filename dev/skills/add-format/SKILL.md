---
name: add-format
description: 给 Ferro 新增一种轨迹/结构文件格式的读写时使用（新 reader / writer、让 ferro convert 认识新扩展名、给 Python 绑定加包装）。
---

# 加文件格式

> 2026-10-03 从 `CLAUDE.md`「扩展项目」原文迁入。动手前先查 `dev/issues.md`
> 里该格式或相近格式的编码陷阱。

1. `ferro-io/src/readers/<fmt>.rs` 返回 `Result<Trajectory>` + `writers/<fmt>.rs`
2. 从 `readers/mod.rs`、`writers/mod.rs` 导出
3. `ferro-cli/src/io_dispatch.rs` 加格式检测（`read_trajectory` 与写侧的 `out_format` **两处**；`out_format` 加了 `OutFormat` 变体后，`write_trajectory` 的 match 由编译器逼着补）
4. `ferro-python/src/io.rs` 加包装

**收尾**：`cargo test` + `cargo clippy` 零警告；`cd ferro-python && cargo check`
（独立 workspace，主 workspace 会跳过它）。
