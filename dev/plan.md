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

无。（2026-10-21 瘦身已于 2026-10-22 全部完成，逐条状态仍在 `plan/slim-2026-10-21.md`）

## 待办索引

| 子文件 | 优先级 | 条目 |
|---|---|---|
| `plan/io.md` | 高 | dataset 额外键搬运；测试缺口（第 8 批，暂缓） |
| | 中 | 带速度的测试 fixture；VASP ML_FF 的 OUTCAR |
| | 低 | DeePMD mixed type 读写（10-22 由高降）；`convert` 写 CIF 的空间群；MACE / NequIP 格式；`read_cube` 合并；CP2K EXTXYZ / XMOL（kind 名、Hartree 能量）；dataset 审查余项 |
| `plan/analysis.md` | 中 | bondlife 手册口径；vanhove `tau` 列（搁置） |
| | 低 | bondlife 后续、配位态动力学、VDOS（10-22 由中降）；D 的统计误差；filter 多层周期镜像；机器学习集成；S(q) Lorch 窗；分析侧逐帧消费 |
| `plan/cli-python.md` | 低 | ferro-structure 补充结构操作（10-22 由中降）；code1 / code2 引用扫描；REPL；rustfmt；ferro-py job 的遗留与存疑 |
| `plan/suspended.md` | 挂起 | Bader 改写组、cube_sdf / chg_sdf 改写组、`map jump` 入口组（随功能改写一并处理，平时不修） |
| `plan/zero-call.md` | 中 | 零调用清单 |
| `plan/review-2026-10-24.md` | 高 / 中 / 低 | 全库 code review 结果（高 4 已修 · 中 42 · 低 32）。中级已处理完（R4-7 随挂起组；R5-2 定不改；R8-3 误报）；低级已处理完（R2-9、R9-6 定不修，R3-6 / R3-8 / R4-5 随挂起组）。只剩挂起组的条目 |
| `plan/review-2026-10-30.md` | 高 / 中 / 低 | 7 代理全库审查 + `/code-review high` 复扫（高 9 · 中约 40 · 低约 90），高级 9 项与 5 项用户拍板已全部完成（2026-10-31）；中级已处理完（2026-11-01；A3-M-1 / M-2、A4-M-3、A5-M7 用户定随 Bader 挂起，A4-M-5 误报，CR-7 定不改）；低级已处理完（2026-11-01；Bader / sdf / `map jump` 挂起组随改写，几条定口径或定不改，理由在报告各条与 `decisions.md`）。遗留待办：A6-L2 的参数名由 CLI / Python 各自补上 |
| `plan/review-2026-11-06.md` | 高 / 中 / 低 | 8 路全库精简审查（过度防御 / 啰嗦 / 多重标准）：高 1 · 中 22 · 低约 150，15 项待用户拍板（D1–D15），7 个修复批次，均未动手 |
