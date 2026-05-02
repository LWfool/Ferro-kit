---
name: 记忆文件双写同步
description: 每次写入或更新 ~/.claude/projects/.../memory/ 时，同步写入项目内 .claude-memory/
type: feedback
---

每次创建或更新记忆文件时，必须同时在以下两处保持一致：

1. `~/.claude/projects/-Users-liuwei-Documents-script-Program-rust-Ferro/memory/`（系统自动加载路径）
2. `/Users/liuwei/Documents/script/Program/rust/Ferro/.claude-memory/`（项目内备份，随 git 提交）

**Why:** 用户换电脑后需从项目目录恢复记忆，两处不同步会导致记忆丢失或冲突。

**How to apply:** 写完系统路径的文件后，立即用相同内容写入项目路径。MEMORY.md 索引也需同步。敏感信息规则优先：涉及敏感内容时仍需先询问用户再写入。
