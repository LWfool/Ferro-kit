# 待办：ferro-cli / ferro-structure / ferro-python / 工具链

> 从 `dev/plan.md` 按需读取；定级规则见 `plan.md` 顶部。

## 优先级中

### ferro-structure：补充结构操作

`rotate.rs` / `orient.rs` / `substitute.rs` / `disturb.rs` / `select.rs`。

移植参考 Multiwfn（`examples/Multiwfn_2026.4.10_src_Linux`）的 `otherfunc3.f90`
`geom_operation`（行 1975–3066）：

| 待办 | Multiwfn 参考 |
|---|---|
| `rotate.rs` | 菜单 3、4（绕笛卡尔轴/键/指定向量旋转、旋转矩阵） |
| `orient.rs` | 菜单 5/6/8/11（对齐键/向量/最长轴/平面到笛卡尔轴或平面） |
| `disturb.rs` | 菜单 18 / `displace_geom`(1879)（高斯随机位移，默认 σ=0.03 Å） |
| `select.rs` | `util.f90:985 str2arr`（`"2,3,7-10"` 选择语法） |
| `substitute.rs` | 无直接对应（Multiwfn 仅 15/16 加删原子） |

其余可参考项：晶胞数学在 `PBC.f90`；菜单 20 边界分子补全、22 原子折叠入胞、
25 提取团簇、28 坐标轴互换。


## 优先级低

### 结构操作的 CLI 入口：supercell / vacuum / box_builder（2026-10-21 用户定登记）

三者都只有库函数，CLI 没有入口：`make_supercell` / `find_supercell_dims`、`add_vacuum`、
`build_box` / `estimate_box_length`（`ferro-structure/src/lib.rs`）。前两者另有 Python 绑定
（`ferro-python/src/structure.rs` 的 `supercell` / `add_vacuum_layer`），`build_box` 连 Python 也没有。
命令形态待定（例如 `ferro build supercell | vacuum | box`）；按定级规则属功能，低。
做的时候 `merge_frames`（Python 有 `merge`）一并考虑。

### 代码注释中的 code1 / code2 引用统一扫描（2026-10-09 提出）

`private/code1`（dump2analysis）与 `code2`（dump2sq）是私有参考代码，**不是行业标准做法** ——
vanhove 的归一化就不除 dr，与 GROMACS / pymatgen 相差 1/dr（手册已于本日改为对照这两者）。
注释里写「follows code1/xxx.c」会让读者以为 ferro 的口径以它为准。2026-10-09 实数：
`angle.rs` 5、`gr.rs` 6、`sq.rs` 5、`scattering_data.rs` 4、`cube_radius.rs` 2、`vanhove.rs` 2
（模块文档与 `calc_vanhove` 文档，均写「follows code1/vanhove.c」）、`tests/npt_normalisation.rs` 1。
逐条判：口径已偏离 code1 的改写为描述 ferro 自己的做法（必要时对照行业实现）；只是记录
「对拍时如何与它对齐」的保留。手册里 gr/msd/vacf/rotcorr 四页的 code1 提及是历史与对拍说明，已看过、不动。

**2026-10-11 审查 G-9 的逐条判断**（25 处与上面实数吻合；改 11、删 2、留 12）：
- **改**：`vanhove.rs:8,138`（优先：写「照 code1」却给出除以 dr 的公式，code1 不除；改为对照 GROMACS / pymatgen）；
  `cube_radius.rs:7`（路径错，应为 `private/code1`，且 M4 后搜索窗已不是 code1 的全遍历）；`angle.rs:20,170`；
  **`angle.rs:489`（用户可见：进 angle 产物 csv 头部，把私有代码名写进用户数据文件；改中性表述）**；
  `gr.rs:5,363,304`（改引 Allen & Tildesley）；`sq.rs:6,107`（改引 Faber & Ziman 1965 / Keen 2001）、`sq.rs:315`
  （(2−δ_ij) 对称求和是 Faber–Ziman 加权的定义）
- **删**：`gr.rs:4`、`gr.rs:314`（纯历史，无判据）
- **留**：`angle.rs:16,268`（计数差一半的对拍说明）；`cube_radius.rs:354`（测试的 brute-force 对照）；
  `gr.rs:575`、`npt_normalisation.rs:34`、`sq.rs:153,705`（对拍基准）；`scattering_data.rs:3,10,18,320`（数据出处）

### ferro-cli：REPL / 脚本模式（2026-08-11 由高降中，2026-09-15 降低）

`main.rs` 现在是子命令分发器，裸 `ferro` 打印分类总览。REPL 落地时改为
**tty 进 REPL、管道读 stdin**（`python`/`node`/`irb` 的惯例；`isatty` 判断，CI 里
`echo ... | ferro` 不会挂住）。

- 依赖 `rustyline`；三种模式：交互 REPL、脚本文件（`ferro -f workflow.mf`）、管道输入
- **必须在同一进程内链接全部命令**：`read` 之后轨迹要留在内存里给后续命令用，
  这正是 REPL 相对 shell 循环的全部价值。这条也是 0.2.0 决定合并成单二进制的理由 ——
  子进程分发做不到状态保持，而链接进来之后再包一层 `fe-*` 前端就是重复
- 脚本语法应建在**已经定型**的子命令树上（`gr -a P -b O` 直接复用 `cmd::traj` 的
  参数结构体），不要另起一套
- 注意与「批处理输入」是两件事：这里是**命令**的批处理（一个脚本跑多条命令），
  那里是**输入文件**的批处理（一条命令跑多个轨迹）。两者可叠加但互不依赖

### 是否采用 rustfmt（待定，暂缓；2026-09-15 由中降低）

收益：把格式从代码审查面里移除。代价三点：

1. 一次性约 80 文件的大 diff，`git blame` 在这些行上全部指向那一次提交
2. 会拆掉现有的刻意列对齐（`gr.rs` 的 `WelfordStats`、`angle.rs` 的 `CellList`、
   `units.rs` 的枚举表）
3. 生成文件要显式排除：`rustfmt.toml` 的 `ignore` 仅 nightly 可用，stable 上需给
   `cp2k_basis_db.rs` 的静态表加 `#[rustfmt::skip]`

若采用：单独一次纯格式提交（`style: adopt rustfmt`）+ 固定 `rustfmt.toml`，
不要混进特性或依赖升级提交。
