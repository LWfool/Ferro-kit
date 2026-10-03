---
name: update-dev-records
description: 用户说「更新开发记录」「更新开发进度」「发版」「打 tag」时使用。先判断这次升不升版本号，再更新 dev/ 的进度、计划、陷阱和版本批次；发版时再同步 ferro-python 版本、写 changelog、打 annotated tag。
---

# 更新开发记录 / 发版

> 2026-10-03 新建（改写）。步骤综合自 `dev/memory.md`「版本号更新规则」、
> `dev/overview.md`「版本规则」、`dev/issues.md`「ferro-python 编码陷阱」，
> 以及提交 `6e1508f`（v0.3.4 发版）的实际改动面。原理去那几处看，这里只放步骤。

## 1. 先定版本号动不动

- 用户说「更新开发记录/进度」→ 默认 patch +1；**同一句话里说了「暂不改版本」就不动**
- 用户说「发版」→ 问清版本号。破坏性改动按规则升次版本位，但用户可以要求走 patch 位
  （`v0.2.1`、`v0.3.4` 都是这类例外，要在 overview 里写明）
- 用户没提版本 → 不动。拿不准就问，不要自作主张 +1

## 2. 更新 dev/（每次都做）

1. `dev/progress.md`：顶部测试总数表（跑 `cargo test` 取实数，不要估）、版本号行、
   对应 crate 小节的现状
2. `dev/plan.md`：已落地的条目从待办挪到「已完成（归档）」。归档只留**判据**与
   **被实测推翻的原计划**，实施步骤删掉、写提交号
3. `dev/issues.md`：本轮踩到的坑按「位置 | 陷阱 | 正确做法」补进对应小节
4. `dev/overview.md`：有破坏性改动时加一节批次说明，标明「版本号**仍是 x.y.z**，未发版」
   或「随 `vx.y.z` 发版」
5. 用户可见的变化写进 `docs/src/changelog.md`（英文，写「你要改什么」，不写理由）

## 3. 升版本号（第 1 步判定要升时）

1. 根 `Cargo.toml` 的 `workspace.package.version`
2. **`ferro-python/Cargo.toml` 手动同步**（独立 workspace，不会跟着变）；
   `pyproject.toml` 是 `dynamic = ["version"]`，不用改
3. `cargo build` 与 `cd ferro-python && cargo check`，让两份 `Cargo.lock` 一起更新
4. `dev/progress.md` 的版本号行

## 4. 发版（仅当用户要求发版）

1. `dev/overview.md` 加「`vx.y.z`（日期 打 tag，发版）」一节，把之前标「未发版」的
   批次标题改成「随 `vx.y.z` 发版」
2. `docs/src/changelog.md` 的 `Unreleased` 改成 ``## `vx.y.z` — 日期``
3. 提交后打 annotated tag：`git tag -a vx.y.z -m "vx.y.z — <一句话>"`。
   中途锚点用 `-alphaN` 后缀，不要占用未发版的版本名
4. 推不推送 tag 问用户

## 5. 提交

直接在 main 上提交，只 stage 本次改动的文件，**不加任何署名行**。
