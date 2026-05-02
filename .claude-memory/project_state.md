---
name: ferro 开发进度
description: 待办事项与非显而易见的编码陷阱（架构见 CLAUDE.md）
type: project
originSessionId: 0e8d82ba-c0c2-4836-af51-2fc7e824f750
---
## 待办

- `ferro-structure` — **supercell 已完成**；vacuum / merge / box 估算待做
- `ferro-python` — PyO3 绑定（Cargo.toml 已注释掉）
- REPL / 批处理模式（ferro-cli/src/main.rs 目前是占位符）

## ferro-structure 实现说明

- crate 已创建并加入 workspace（2026-05-02）
- `src/supercell.rs` 提供两个公开函数：
  - `make_supercell(frame: &Frame, nx, ny, nz) -> Result<Frame>`：操作单帧；bonds 按副本偏移重映射；energy/forces/stress/velocities 不复制；charge 按副本数线性缩放；multiplicity 重置为 1
  - `find_supercell_dims(cell, n_atoms, min_length, min_atoms) -> [usize; 3]`：先用 ceil(min_length/L) 满足长度约束，再递增最短轴满足原子数约束（对应 private/make_supercell.py 的算法）
- API 设计：单帧而非 Trajectory（用户决策）
- 12 个单元测试全绿，clippy 零 warning

## 编码陷阱

- 浮点断言：`cartesian_to_fractional` 误差 ~1e-15，用 `< 1e-10`，不能 `assert_eq!`
- voxel 测试：原子坐标放格点中心 `(n+0.5)/N * L`，不要放 voxel 边界，否则落在哪个格子里取决于浮点误差
- plotters 位置周期：`wrap_position` 用 `rem_euclid(1.0)` 而非 `x - x.floor()`；`floor` 在 ~1e-15 误差下会在 x=0 边界跳变
- plotters 字体：依赖必须开 `ttf` feature，`default-features = false` 会把它禁掉，macOS 上文字渲染 panic
- plotters 借用：`BitMapBackend::new(&path, ...)` 借用 `path`；绘图块包进 `{}`，让 root 在 `Ok(path)` 前析构
