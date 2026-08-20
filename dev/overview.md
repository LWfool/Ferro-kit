# Ferro — 项目概要

> **文档分工**：`docs/src/`（mdBook）= 怎么用，面向使用者 · `dev/` = 为什么这样定 +
> 现状 + 待办，面向开发 · `CLAUDE.md` = 给 Claude 的硬约束与导航。
> 用法问题（CLI 参数、输出列结构、数据模型字段）一律查 `docs/src/`。

## 定位

Rust 重写的计算化学后处理工具链，面向周期性体系（晶体、表面、玻璃）的 MD 轨迹分析。
目标：替代原有 Python 脚本，兼顾性能与可扩展性，未来计划集成深度学习工作流。

## 架构

```
ferro-cli / ferro-python        ← 唯一允许组合多个 crate 的入口层
    ├── ferro-core              ← 纯数据结构 + 静态参考数据 + spin（未成对电子推断）
    │                             + Table（分析产物的中立载体）
    │                             + AtomType（网络分类结果，结构化而非字符串）
    ├── ferro-io       → core   ← 格式读写 + write_table（分析产物的唯一出口）
    ├── ferro-structure→ core   ← 超胞、真空层、合并、初始盒子
    ├── ferro-analysis → core   ← 纯计算，不碰文件系统；结果暴露 to_tables()
    │                             md/、network/、dft/（Bader、ChgSDF）、ml/（未来）
    └── ferro-workflow → core   ← QC 输入生成（Gaussian / CP2K / QE，含基组赝势库）
```

中间层 crate 之间不允许互相依赖。**但共享类型可以下沉** —— 判据：一个类型该放
`ferro-core`，当且仅当两个以上中间层需要叫出它的名字。两个方向各一个范例：
`Trajectory`（io 产出、analysis 消费）与 `Table`（analysis 产出、io 消费）。

这条规则的完整论证见 `issues.md`「分析产物为什么不进 ferro-io」—— ASE 与 pymatgen
是 Python，**能**写出循环依赖却依然不写，说明这是设计选择而非语言约束。

## 命令入口（0.2.0 起；`net` 于 0.2.1 降为叶子命令）

**单二进制 `ferro`**，子命令按**产物**分组（不是按实现它的 crate）：

```
ferro traj  gr | sq | msd | angle | vacf | rotcorr | vanhove   → 堆叠 csv + 可选 PNG
ferro map   density | velocity | force | radius | sdf | chg-sdf → 逐输入一个 .cube
ferro net                                                      → 六张堆叠 csv
ferro bader | convert | info | job
```

`-i` 恒为多值并自展开 glob；逐文件独立分析，结果堆叠成一份带 `file` 列的 csv。
`-o` 是**文件名后缀**，不是路径；路径走 `--outdir`。

参数与输出列结构见 `docs/src/cli-reference.md`。

## 编码规范

- `///` / `//!` doc 注释 → 英文
- `//` 内部注释 → 中文

## 版本规则

版本号在根 `Cargo.toml` 集中管理，各 crate 继承 `workspace.package.version`；
`ferro-python` 是独立 workspace，需手动同步。

**只在用户明确要求时才动版本号**，不要每次改代码就自动 +1。破坏性改动升次版本位
（0.1.15 → 0.2.0），其余升 patch 位。发布点打 annotated tag。

**已知例外**：`v0.2.1`（2026-08-12，net 重构）含破坏性改动（`net qn`/`net type`
删除、标签格式全变、表名与列结构全变），按上面的规则本应是 0.3.0，但按用户当时的
明确要求走了 patch 位。翻 git 历史时注意：**0.2.0 → 0.2.1 之间有 breaking change**。

## 未发版的破坏性改动（`v0.2.1` 之后）

按用户要求暂不动版本号，等配套的 Python 绘图脚本完成后一并升。届时应计入：

**net 重做（2026-08-12）**

| 改动 | 影响 |
|---|---|
| 六张表改名 | `bridge`→`qn`、`partner`→`qn_partner`、`oxy`→`ligand_type`、`cn`→`coordination`、`mean`→**删除** |
| 新增 `composition` 表 | 一物种一行的结构组成，取代 `average` |
| Al 退出 Qn 表 | 非 Qn 形成子只在 `coordination` 出现 |
| 标签数字换义 | 非 Qn 形成子由桥接数改为配位数（`Al_4` = 四配位） |
| 两套标签词汇 | 分布表 `P-Q2`（单元），`linkage` 与导出轨迹 `P_2`（原子） |
| `linkage` 加列 | `linkage`（展示）、`ligand`（配体元素，进键） |
| `ligand_type` 改列 | `label` 合并为 `Al-O_b-P`，原标签移入 `type` |
| 配体 `fraction` 分母 | 全体配体原子 → 该配体元素 |
| 新增 `--qn` | 替换默认 Qn 名单 `{B,P,Si}` |

**traj 产物命名（2026-08-13）**

| 改动 | 影响 |
|---|---|
| 文件名加 label 段 | `gr.csv` → `gr_P-O.csv` / `gr_all.csv`；六个命令全改，**含无筛选时的 `_all`** |
| 新增 `--outdir` | 进 `CommonArgs`，覆盖 11 个命令的 csv/png/cube/导出轨迹；`chg-sdf` 单独加同名参数 |
| `traj sq` 移除 `-a/-b/-x/-y` | 旧命令行直接报 `unexpected argument`；按 label 分辨的 partial 从 CLI 消失 |

**net 口径改为文献 $Q^n_m$（2026-08-20）**

| 改动 | 影响 |
|---|---|
| `qn` 只数同元素连接 | **数值全变**。`P-Q3` 由 40.4% 变 0.27%，`mean_qn` 由 2.40 变 0.95。异核桥进 `m_<X>` 列，总桥 = `n + Σm` |
| `bridges_to` 数连接非桥氧 | 三簇配体贡献 2 而非被跳过，Σm 不再亏空 |
| `qn_partner` 去掉自身元素列 | `m_P` 恒等于 `qn`，重复列 |
| `linkage` 列改名 | `n_bridge_a/b` → `qn_a/qn_b`，值改为同核连接数 |
| `[inputs]` 新增 `mean_n_bo` | 桥氧个数（旧口径那个数）与 `mean_qn` 并列 |
| 新增共边告警 | 形成子共享 ≥2 配体时提示核对 cutoff |

依据是 arXiv 2510.13545 等文献原文：n 只数同核桥，m 数异核桥，总桥由两者相加；
铝磷酸盐文献专门指出「n 是总桥氧数」是领域内的已知误解。判据与排查清单见
`issues.md`。

按仓库自己的规则这三批合起来应升次版本位（0.3.0）。
