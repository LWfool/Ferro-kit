# 后续计划（大纲）

> 本文件只放大纲：定级规则、进行中的事、待办索引。**详细内容在 `dev/plan/` 的子文件里，按需读取**，
> 不要一次全读。已完成的事：「改了什么」看手册 `docs/src/changelog.md`，「为什么这样定」看
> `dev/decisions.md`，「怎么改的」翻 git 历史；「现在是什么样」看 `progress.md` 与 `docs/src/`。
>
> **定级规则**（2026-10-11 用户定）：
> - **高**：重大 bug，影响数据准确性（静默错数据、错物理结论）
> - **中**：缺陷、隐藏的不完善问题（报错不当、边角情况、口径与手册不符等）
> - **低**：功能性条目（新功能、新开关、新格式支持），**除非用户明确要求提级**
>
> 新条目写进对应子文件，并在下面的索引补一行；完成后从子文件删去，判据移进 `decisions.md`，
> 用户可见的变化写进 changelog。

## 进行中

**2026-10-21 瘦身**（报告与逐条状态：`plan/slim-2026-10-21.md`）。S1 文档改正、S2 删死代码已完成；
P1 报告入库、P2 本文件拆分（含 changelog 补缺、全仓引用改指、CLAUDE.md 与技能同步）、
P3 A-2（`_with_*` 合并回原名、收口 `ferro-io` 根上的再导出）、P4（A-5 – A-12、B-10、C-5 – C-8）、
D 组（手册压缩，7947 → 6803 行）已完成。接下来另交计划：C-1 / C-2（两份 README）、E 组其余（`dev/` 压缩）

## 待办索引

| 子文件 | 优先级 | 条目 |
|---|---|---|
| `plan/io.md` | 高 | DeePMD mixed type 读写；dataset 额外键搬运；测试缺口（第 8 批，暂缓） |
| | 中 | 带速度的测试 fixture；VASP ML_FF 的 OUTCAR |
| | 低 | `convert` 写 CIF 的空间群；MACE / NequIP 格式；`read_cube` 合并；CP2K EXTXYZ / XMOL（kind 名、Hartree 能量）；dataset 审查余项 |
| `plan/analysis.md` | 中 | bondlife 后续；配位态动力学；VDOS；D 的统计误差；bondlife 手册口径 |
| | 低 | filter 多层周期镜像；机器学习集成；S(q) Lorch 窗；分析侧逐帧消费；库层参数无 `validate`（C-2） |
| `plan/cli-python.md` | 中 | ferro-structure 补充结构操作 |
| | 低 | **结构操作的 CLI 入口（supercell / vacuum / box_builder）**；code1 / code2 引用扫描；REPL；rustfmt；ferro-py job 的遗留与存疑 |
| `plan/suspended.md` | 挂起 | Bader 改写组、cube_sdf / chg_sdf 改写组、`map jump` 入口组（随功能改写一并处理，平时不修） |
| `plan/zero-call.md` | 中 | 零调用清单 |
