# 协作记忆

> 被用户纠正过的协作做法、容易答错的工具语义。格式：规则 → **Why** → **How to apply**。
> 新的纠正记在这里，**不写进** `~/.claude/projects/.../memory/`。
>
> 2026-10-03 由 `/organize-project` 从用户级记忆目录迁入（原文搬移，改动处已注明）。
> **没迁的两条**：`dev_records.md`（读 dev/、三份文档分工）—— 已写进 `CLAUDE.md`
> 文件头与「四类信息放在哪」，是规则不是记忆；`feedback_paid_mode.md` —— 跨项目的
> 协作习惯，与本项目无关，留在用户级。

## 文档精简的方向

要求「精简开发记录、节省上下文」时，**不要删 dev/ 里的内容**。正确做法是把重复
内容留在 `dev/`，转而精简 `CLAUDE.md`。

**Why:** 用户的开发习惯是每次会话第一件事让 Claude 读开发记录，所以 dev/ 是他
真正阅读的文件；CLAUDE.md 是给 Claude 的。删 dev/ 等于为了省 Claude 的上下文而
损害用户自己的阅读体验。2026-08-14 我第一次动手时删的是 dev/overview.md，被用户
当场拦下并给出反方向。

**How to apply:** 精简前先查全仓有几份文档（本项目有三份：`CLAUDE.md`、`dev/`、
`docs/src/` 的 mdBook 手册，当时共 7871 行）。按读者划分职责后再去重，不要凭
「哪个文件大」下手。CLAUDE.md 只留「即使没读 dev/ 也必须遵守」的硬约束。
删除时保留三类内容：**设计判据、被实测推翻的原计划、编码陷阱** —— git 历史与
用户手册都回答不了它们；可以删的是实施步骤清单和已实现的计划。
与手册**逐条重复**的「改了什么」也可以改成指针：2026-10-22 用户同意把 overview 各版本的
破坏性改动表改为 changelog 指针（瘦身 O-3），前提是先核对 changelog 对应版本节确实都有。

## 版本号更新规则

只有当用户明确说"更新开发记录"或"更新开发进度"时，才同时更新版本号（Cargo.toml patch +1）和 dev/progress.md。**但用户可显式抑制版本号变更**：若用户同时说"暂不修改版本"/"不改版本号"等，则只更新 dev 记录内容，**保持版本号不变**。

**Why:** 用户不希望每次代码修改后自动递增版本号，版本号更新应与开发记录同步进行，由用户主动触发。有时多轮工作合并到一个版本，用户会要求更新记录但暂缓 bump（2026-05-16 明确这样要求过一次）。

**How to apply:** 代码修改完成后，不要自动递增版本号。用户要求"更新开发记录/进度"时：默认一并 bump 版本 + 更新 dev/{progress,plan,issues}.md；但若用户在同一指令中明确表示暂不改版本，则只改 dev 记录内容，Cargo.toml 版本与 progress.md 版本号行均不动。

## 中文不要经 Bash 的 heredoc / 命令行参数传

本机 Bash 工具把 heredoc 或命令行参数里的中文偶尔弄成乱码：`git commit -m "中文"` 写出过
乱码提交信息，`python - <<EOF` 里含中文的脚本报 `Non-UTF-8 code`，或把乱码写进文件。

**Why:** 2026-10-16 `ferro.jobs` 一轮里出现了五六次，一次乱码提交信息已 amend 掉；其余在
提交前发现。纯 ASCII 的 heredoc 不受影响。2026-10-17 又犯两次：`cat >> test.py <<EOF` 追加含中文
注释的测试，乱码随提交进了库（`3cb0eba`，`1cad05b` 修）；`git commit -m` 中文乱码（已 amend）。
2026-10-18 又犯两次，而且**不只是中文**：heredoc 里的 Python 只要含 `→`、`ᵀ` 这类非 ASCII
字符同样报 `Non-UTF-8 code`，脚本整个没执行、文件原样未改 —— 输出被折叠时很容易当成改成功了。
2026-10-29 又犯三次（review 中级修复一轮）：`git commit -am` 中文乱码（已 amend）；heredoc
Python 改测试、拆提交两次报 `Non-UTF-8 code` 未执行，其中一次让本该分开的 R9-1 / R9-2 并进同一提交。

**How to apply:** 提交信息先用 Write 写到 scratchpad 文件再 `git commit -F`；改含中文的源码、
文档用 Edit / Write，不用 heredoc 里的 Python 替换；Bash 里 grep 中文也会报
`illegal byte sequence`，改用 Python 读文件查。改完用
`python -c "open(p, encoding='utf-8').read()"` 验一遍。

## 提交直接在 main

用户要求本项目的代码提交**直接在 main 分支进行**，不要先开特性分支再合并/PR。

**Why:** 单人项目，特性分支/PR 流程是多余开销；用户明确拒绝了 `git switch -c` 并指示「直接合并 main 分支」。

**How to apply:** 被要求「提交修改」时，直接在 main 上 `git add <具体文件> && git commit`，不要先 `git switch -c`。只在用户另行要求时才开分支/PR。注意只 stage 预期文件（避免误纳入 stray 副本如 `xxx 2.rs`）。提交信息**不加任何署名行**（迁入时改写：原文要求保留 `Co-Authored-By`，与用户全局规则冲突，以全局规则为准）。相关：「版本号更新规则」。

## CLI 的 cargo 包名是 `ferro`，不是 `ferro-cli`

`cargo test -p ferro-cli` 报 `package ID specification ... did not match any packages`，一个测试都不跑；
输出再经 `grep FAILED` 之类过滤，错误行被滤掉，看起来像「CLI 测试全过」。

**Why:** 2026-10-29 review 低级修复一轮，R5-3、R9-8、R7-7、R7-8 四步都以为跑过了 CLI 测试，实际没跑；
到 R7-9 才看见这行错误，补跑全量 `cargo test` 才确认无失败。

**How to apply:** 单测 CLI 用 `cargo test -p ferro`；拿不准就跑全量 `cargo test`。过滤输出时把
`error` 也放进 grep 的模式，别只盯 `FAILED` / `panicked`。
