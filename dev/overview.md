# Ferro — 项目概要

## 定位

Rust 重写的计算化学后处理工具链，面向周期性体系（晶体、表面、玻璃）的 MD 轨迹分析。目标：替代原有 Python 脚本，兼顾性能与可扩展性，未来计划集成深度学习工作流。

## 架构（Cargo workspace，严格分层）

```
ferro-cli / ferro-python        ← 唯一允许组合多个 crate 的入口层
    ├── ferro-core              ← 纯数据结构 + 静态参考数据 + spin（未成对电子推断）
    │                             + Table（分析产物的中立载体）
    │                             + AtomType（网络分类结果，结构化而非字符串）
    ├── ferro-io       → core   ← 格式读写 + write_table（分析产物的唯一出口）
    ├── ferro-structure→ core   ← 超胞、真空层、合并、初始盒子
    ├── ferro-analysis → core   ← 纯计算，不碰文件系统；结果暴露 to_tables()
    │                             md/（gr、sq、msd、angle、vacf、rotcorr、vanhove、cube）
    │                             network/（桥接数、配体分类、CN、桥联统计）、dft/（Bader、ChgSDF）
    └── ferro-workflow → core   ← QC 输入生成（Gaussian / CP2K / QE，含基组赝势库）
```

中间层 crate 之间不允许互相依赖。**但共享类型可以下沉**——判据：一个类型该放
`ferro-core`，当且仅当两个以上中间层需要叫出它的名字。两个方向各一个范例：
`Trajectory`（io 产出、analysis 消费）与 `Table`（analysis 产出、io 消费）。

## 命令入口（0.2.0 起；`net` 于 0.2.1 降为叶子命令）

**单二进制 `ferro`**，子命令按**产物**分组（不是按实现它的 crate）：

```
ferro traj  gr | sq | msd | angle | vacf | rotcorr | vanhove   → 堆叠 csv + 可选 PNG
ferro map   density | velocity | force | radius | sdf | chg-sdf → 逐输入一个 .cube
ferro net                                                      → 六张堆叠 csv
                                                                 + 可选标注轨迹
ferro bader | convert | info | job
```

`-i` 恒为多值并自展开 glob；逐文件独立分析，结果堆叠成一份带 `file` 列的 csv。
`-o` 是**文件名后缀**，不是路径。

## 核心数据模型

- 顶层类型：`Trajectory { frames: Vec<Frame> }`，单帧文件也用此类型
- `Frame`：`atoms: Vec<Atom>`、`cell: Option<Cell>`、`pbc: [bool; 3]`、能量/力/应力/速度
- `Atom`：`element`、`position`（Å，Cartesian）、`label`、`mass`、`magmom`、`charge`
- `Cell`：`matrix: Matrix3<f64>`（行向量 a, b, c，单位 Å）

## 内部单位

| 量 | 单位 |
|---|---|
| 长度 | Å |
| 能量 | eV |
| 力 | eV/Å |
| 时间 | fs |
| 质量 | amu |
| 电荷 | e |
| 温度 | K |

## 编码规范

- `///` / `//!` doc 注释 → 英文
- `//` 内部注释 → 中文

## 版本规则

版本号在 `Cargo.toml`（workspace 根）集中管理，各 crate 继承
`workspace.package.version`；`ferro-python` 是独立 workspace，需手动同步。

**只在用户明确要求时才动版本号**，不要每次改代码就自动 +1。破坏性改动升次版本位
（0.1.15 → 0.2.0），其余升 patch 位。发布点打 annotated tag（`v0.1.15`、`v0.2.0`）。

**已知例外**：`v0.2.1`（2026-08-12，net 重构）含破坏性改动（`net qn`/`net type`
删除、标签格式全变、表名与列结构全变），按上面的规则本应是 0.3.0，但按用户当时的
明确要求走了 patch 位。翻 git 历史时注意：**0.2.0 → 0.2.1 之间有 breaking change**。
