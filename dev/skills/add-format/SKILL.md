---
name: add-format
description: 给 Ferro 新增一种轨迹/结构文件格式的读写时使用（新 reader / writer、让 ferro convert 认识新扩展名、给 Python 绑定加包装）。
---

# 加文件格式

> 2026-10-03 从 `CLAUDE.md`「扩展项目」原文迁入；2026-11-01 按审查 2026-10-30 A7-M4 重写第 1–3 步。
> 动手前先查 `dev/issues.md` 里该格式或相近格式的编码陷阱。

1. `ferro-io/src/readers/<fmt>.rs` 的 reader 直接返回 `Result<(Trajectory, Vec<Notice>)>`
   （要告诉用户的事放进 `Vec<Notice>`，**不打印**；例：`readers/xyz.rs` 的 `read_xyz`）
   + `writers/<fmt>.rs`。读文本走 `readers/util.rs` 的 `read_text` / `text_lines`；末帧截断一律调
   `util::drop_truncated`（丢帧并告警，一帧不剩报错）；reader 入口给错误套一层 `parsing {path}`
2. 在 `readers/mod.rs`、`writers/mod.rs` 声明模块即可，**单个 reader 不在 `ferro-io/src/lib.rs` 根上导出**：
   外部一律经 `ferro_io::read` / `read_as`（`tests/warnings_export.rs` 也走 `read` 查告警）
3. `ferro-io/src/format.rs`：`Format` 加变体，补 `ALL`、**`LABELS`**（`--format-in` 与 Python
   `read(format=)` 的名字，`program/kind` 或无主格式的裸名；漏了就叫不出新格式，
   `a_label_reads_a_file_whatever_its_name` 会拦）、`name`、`detected_by`、`can_write`、
   `multi_frame`、`from_name`（或内容嗅探）、`read_as`、`write_one`。CLI 的 `convert` 格式表由
   `Format::ALL` 生成，CLI 与 Python 的读写都走 `ferro_io::read` / `write`，**不用再改别处**
4. 手册：`docs/src/cli-reference.md` 的 convert 格式表（含 `--format-in` 列）与 `docs/src/python/api.md` 的 `read` 表

**收尾**：`cargo test` + `cargo clippy` 零警告；`cd ferro-python && cargo check`
（独立 workspace，主 workspace 会跳过它）。
