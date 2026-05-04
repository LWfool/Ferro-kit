---
name: ferro 开发进度
description: 待办事项与编码陷阱（架构见 CLAUDE.md）
type: project
originSessionId: 0e8d82ba-c0c2-4836-af51-2fc7e824f750
---
## 待办

- `ferro-python` — PyO3 绑定（Cargo.toml 已注释掉）
- REPL / 批处理模式（ferro-cli/src/main.rs 目前是占位符）

## ferro-structure（2026-05-04 全部完成，55 测试）

- 4 个模块：supercell / vacuum / merge / box_builder
- API 设计：单帧操作（&Frame → Result<Frame>），非 Trajectory
- box_builder 新增依赖 `rand = "0.8"`
- 设计文档：`docs/superpowers/specs/2026-05-04-ferro-structure-design.md`

## 编码陷阱

- 浮点断言：`cartesian_to_fractional` 误差 ~1e-15，用 `< 1e-10`，不能 `assert_eq!`
- voxel 测试：原子坐标放格点中心 `(n+0.5)/N * L`，不要放 voxel 边界
- plotters 周期：`wrap_position` 用 `rem_euclid(1.0)` 而非 `x - x.floor()`
- plotters 字体：依赖必须开 `ttf` feature，`default-features = false` 会禁掉它
- plotters 借用：`BitMapBackend::new(&path, ...)` 借用 `path`；绘图块包进 `{}` 让 root 先析构
- snake_case：clippy 报 `non_snake_case`，变量名 `L` → `box_len`
