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
2026-10-31 又犯一次（10-30 审查 H6）：heredoc Python 往 `issues.md` 追加含中文的表格行。这次碰巧没乱码，
事后用 `python -c` 核对过 —— 侥幸不是理由，规则照旧。
2026-10-31 又犯一次（10-30 审查 CR-5）：heredoc Python 的替换串里有 `Å`，报 `Non-UTF-8 code` 未执行，
随后的 cargo test 照常 999 过（测的是旧代码），看结果文件里第一行的 SyntaxError 才发现。纯 ASCII 的
heredoc 本轮用了十几次都没事 —— **动手前看一眼替换串里有没有 `Å`、`→`、`σ` 这类字符**，有就用 Edit。
2026-11-01 又犯两次（10-30 审查 A7-M2 的 `—`、A7-M6 的中文替换串），都没乱码，事后核对 diff 无 U+FFFD。
两次都是「同一文件改几处，想一次脚本做完」—— 几处就发几个 Edit，同一回复里并发即可。
2026-11-01 低级项那一轮又犯四次，形式变了：heredoc Python 里写「（已修）」、`sed -i` 的替换串含中文、
`git commit -m` 带中文标题、改提交信息里的数字时 `sed` 带中文。都没乱码，但规则照旧。管用的做法：
把含中文的替换写成 scratchpad 里的 `.py` **文件**（Write 写入，不经 shell），再 `python 文件`；
另备一个只收 ASCII 的替换脚本（遇非 ASCII 直接拒），它拒掉时**看清它拒了**——一次被拒没察觉，
只提交了半个修改（A4-L3，`d763063` 后由 `0f501c7` 补上）。

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

## 提交信息里的哈希一律从 `git log` 取

写 R5 基线、「见某提交」时，哈希用 `git log -1 --format=%h <ref>` 取出来再填，不凭记忆写。
amend 会换哈希：amend 过的提交，后面再引用要重新取。

**Why:** 2026-11-07 审查批 8–14 里三次写出不存在的哈希（`ec04c8b`、`0d1b5f0-era`），一次引用了
amend 前的旧哈希（`0aab0c0`，留在批 11 的提交信息里未改）。

**How to apply:** 提交信息里出现 7 位哈希的地方，先跑一次 `git log --oneline -N` 对照。

## macOS 的 `iconv` / `grep` 判 UTF-8 不可信

`iconv -f UTF-8 -t UTF-8` 会把本仓里合法的 UTF-8 文件（如 `decisions.md`、`gr.rs`）判成坏的；`grep`
在默认 locale 下报 `illegal byte sequence` 也不说明文件坏了。要验编码用 Python：
`open(f, 'rb').read().decode('utf-8')`，新增行里查 `\ufffd`。

**Why:** 2026-11-07 用 iconv 扫出 7 个「坏文件」，Python 解码全部通过，rustc 也照常编译。

