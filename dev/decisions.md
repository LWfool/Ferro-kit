# 设计判据与用户定

> 已执行完的审查与计划留下的**判据、用户定、被实测推翻的原计划**，原文搬自 `dev/plan.md`（2026-10-21）。
> 「改了什么」看手册 `docs/src/changelog.md`，「怎么改的」翻 git 历史。按需读取，不必通读。

## 手册迁出（2026-10-21 瘦身 D 组）

手册只写「是什么 / 怎么用」。下面是从 `docs/src/` 撤下、仍有留档价值的段落，**原文照搬**（英文），
按来源页分小节。从代码能直接读出的实现细节（rayon 的 `par_iter` / `fold`、源文件路径、耗时）
撤下时直接删了，没有搬。

### `analysis/angle.md`「Implementation Notes」（D-5）

> **Linked-Cell List.** Brute-force neighbour search is $O(N^2)$ per frame.  ferro uses a linked-cell list to reduce this to O(N·k) where $k$ is the average number of atoms in the search volume.
>
> For each central atom B, only the 27 cells within ±1 cell in each direction are searched. The cell size is chosen so that a sphere of radius $r_\text{cut,max}$ fits within the search volume:
>
> $$n_k = \left\lfloor \frac{1}{r_\text{cut,max} \cdot \|({\mathbf{M}^T})^{-1}\|_k} \right\rfloor$$
>
> This is exact for triclinic cells.

### `analysis/gr.md`「Implementation Notes」（D-5）

> - Pseudo-element labels (e.g. `"P0"`, `"Ob"`) are supported: `elem_z` resolves them by stripping numeric/alphabetic suffixes.
> - Column ordering uses atomic number Z as the primary sort key, with string labels as tiebreaker for deterministic ordering of pseudo-elements.

**存疑，未处理**：`elem_z`（`md/gr.rs`）与 `elem_z_local`（`md/sq.rs`）都是从标签字符串剥后缀反推元素，
与 CLAUDE.md「绝不从标签字符串反解析」的精神相冲（那条针对 `AtomType::label()` 的数字，这里是元素）。
`Site` 本来就带 `element` 字段，排序可以直接读它。另外 g(r) 产物早已是长表，「Column ordering」的说法
需要按现状核对。

### `dataset/collect.md`「Checked against dpdata」与三处开发口吻（D-8）

手册里只留结论表与两条用户须知。下面是撤下前的整节原文（含当时的四条论证与跨格式的 virial 符号核对），
以及同页删去的三句（2024 以前版本提示的由来、kind 拆分「可以做成开关」、`-o` 旧默认 `.` 被否决的理由）。
注意：cp2kdata 现在用的是用户自己 fork 并更新的版本（见下文「瘦身 P4 的判据」），这里的 0.7.4 是当时的对照版本。

> ## Checked against dpdata
>
> Every quantity below was compared, number by number, against
> [cp2kdata](https://github.com/robinzyb/cp2kdata) 0.7.4 reading the same file —
> the reference implementation of CP2K support for
> [dpdata](https://github.com/deepmodeling/dpdata). CP2K 2025.2 throughout.
>
> ### Single point — 457 atoms, via dpdata
>
> `dpdata.LabeledSystem(out, fmt="cp2kdata/e_f")` against the `deepmd/npy` system
> `collect` writes. Atoms were paired **explicitly by position** rather than
> sorted, so a permutation could not hide in the comparison; the pairing came out
> as the identity and every element matched.
>
> | | max\|diff\| | relative |
> |---|---|---|
> | energies | 3.2e-10 | 1.4e-15 |
> | cells | 0 | 0 |
> | coords | 0 | 0 |
> | forces | 4.4e-11 | 5.7e-12 |
> | virials | 2.8e-14 | 2.1e-16 |
>
> ### AIMD — 2000 frames, 112 atoms
>
> | | result |
> |---|---|
> | energy, all 2000 frames | max\|diff\| **0** |
> | stress, all 2000 frames | relative 1.2e-15 |
> | coords, 2000 × 112 × 3 | **0** |
> | forces, 2000 × 112 × 3 | **0** |
> | duplicate frames at the restart seam | none |
>
> ### Four things the comparison makes explicit
>
> **The force residual is a constant, not the parser.** `BOHR_TO_ANG` is written
> to ten digits (`0.529_177_210_9`) where CODATA 2018 has
> `0.529177210903`, a relative difference of 5.669e-12 — which is the force
> residual, to the digit. Redo the same comparison using ferro's own constants on
> both sides and the forces come out bit-identical. The energy residual, 1.4e-15,
> is ordinary floating-point noise: `HARTREE_TO_EV` *is* the CODATA 2018 value and
> dpdata simply carries two more digits from its own derivation.
>
> **The initial configuration of an MD segment is not collected.** Frames are
> anchored on `MD| Step number`, and CP2K dumps the starting structure *before*
> the first such line. For the reference run that is 1 frame out of 2001. The
> restart point is skipped too, but it has no trajectory block of its own, so
> nothing complete is lost there — which is also why the seam produces no
> duplicate frames.
>
> **Coordinates and forces are not verifiable through cp2kdata.** Its
> `parse_pos_xyz` / `parse_frc_xyz` read `*-pos-*.xyz` and `*frc*.xyz` *sibling
> files*; given a self-contained log it raises `No atomic coordinates found in
> cp2k output`. The figures above therefore come from an independent extractor
> written against the xyz format itself, cross-checked with the `i = N, time = …,
> E = …` comment that CP2K writes into every block — which is what pins each frame
> to its energy without either parser having a say. If you want that leg covered
> by a third-party implementation instead, have CP2K also write `-pos-1.xyz` and
> `-frc-1.xyz`; dpdata reads those.
>
> **The cell is read at full precision, cp2kdata's at three decimals.** ferro
> takes the 12-field line after the force block; cp2kdata takes `CELL| Vector`,
> which CP2K prints rounded:
>
> ```
> cp2kdata  [11.415      11.415       8.029    ]
> ferro     [11.41476108 11.41495281  8.02853232]
> ```
>
> Rounding ferro's to three decimals reproduces cp2kdata's exactly. Under NVT
> cp2kdata repeats that rounded cell for every frame, so its volumes — and hence
> its virials — carry a ~1e-4 relative offset that ferro does not.
>
> ### The virial sign, checked across formats
>
> Agreeing with cp2kdata alone cannot rule out both being wrong the same way, so
> the same comparison was run on a VASP `OUTCAR` through dpdata's own reader,
> which shares no code with cp2kdata. Energies, coordinates, forces and cells come
> out identical and the virial trace agrees in sign and magnitude
> (+1170.9134 eV both ways).
>
> The 8e-9 residual there is neither tool's parser either: VASP prints the virial
> twice, once as `Total` (eV, directly) and once as `in kB`, and **both** ferro
> and dpdata reconstruct it from `in kB` × volume. That path inherits the coarser
> printed precision — 385.09371 eV against the 385.09354 eV on the `Total` line.
>
> It exists because ferro has no *single-point* output from those releases to
> test against; the block shapes were read off real 2023.1 / 2023.2 / 2024.1 AIMD
> logs instead.
>
> If you need the kinds split
> > into separate types, say so and it can become a flag.
>
> An earlier default of `.` was rejected for scattering `.npy` files through
> whatever directory you happened to be in; a default that follows the input
> cannot do that.

## 瘦身 P4 的判据（2026-10-21）

- **dpdata 不当 CP2K 的独立对照**（用户定，A-9）：dpdata 读 CP2K 输出靠 cp2kdata 插件，原作者已
  不再维护，现用的是用户自己 fork 并更新的版本（用户测试时也与 deepmd 环境中的 cp2kdata 对比过）。
  因此 `cp2k_md.rs` 那条只打印的探针直接删除，没有去补 dpdata 参考值；CP2K 的端到端覆盖靠
  `tests/` 里的三份 CP2K fixture。VASP 与 QE 的 dpdata 对拍不受影响。
- **最小镜像三份、MIC 上界两份不合**（A-11，R6）：契约不同（可失败与否、是否缓存逆矩阵、是否只支持
  立方盒）；`box_builder` 那份在 relax 内循环里，换成 `PeriodicCell` 会改变舍入，同一 seed 出的盒子
  就变了。各处已加互指注释。`gr` 改调 `PeriodicCell::mic_bound` 后可删 `Cell::minimum_image_cutoff`，
  但要动 analysis，未做。
- **`.gitignore` 保留 `.env`、`config/secrets.yml`、`*.key`**（C-6）：同为模板残留、本仓也没有这些文件，
  但它们是防误提交密钥的兜底，删掉没有收益；删的只是 `CLAUDE.md` 与会误伤 `lib/`、`var/` 的 Python 模板条目。
- **出图只合并面板网格**（C-5）：各模块 4 行的 `plot()` 入口不合（收益为负，net 已分岔）；
  gr / angle 的 `build_figure` 合并网格后只剩三行循环与图例，不再共享。

## `ferro-io` 的公开面：告警随返回，按格式的 reader 不上根（2026-10-21，瘦身 P3，用户定）

- **告警与统计是读取结果的一部分**，reader 只有一个名字，返回 `(Trajectory, Vec<Notice>)` 或
  `(Trajectory, AimdStats)`。原先「原名丢弃告警 + `_with_warnings` / `_with_stats` 带告警」两套并存，
  丢弃版本只在测试里有人调，库调用方却最容易顺手拿到它。
- **根上只导出仓外真有调用的项**。按格式的 reader（xyz、pdb、poscar…）一律经 `read` / `read_as`
  走，告警在 `ReadReport` 里；仓外在用的 `read_lammps_dump`、`read_deepmd_npy`、`read_qe_out`、
  `read_aimd` 保留。`ReadReport` 出现在 `read` 的返回类型里，保留导出。
- **撤回 B-13**：当时为让库调用方拿到 xyz / extxyz / pdb / lammps_data 的告警而导出了四个
  `_with_warnings`；`read` 同样带出这些告警，两条路并存只会让 dead_code 失效。用户在 A（撤回）/
  B（保留四个导出）之间选 A，`ferro-io/tests/warnings_export.rs` 改为经 `read` 断言。

## 2026-10-19 四次审查遗留清理（10-11 / 10-03 / 10-01 / 09-27）

范围是下面四节里尚未修复的条目。涉及挂起组的条目只归类、不修。每批开工前先把逐条计划交用户审核；
每条修复单独提交，并带回归测试。不改行为的提交，golden master 须逐字节一致。

- **作废**：依附已删的 `ferro-workflow` / `cmd/job.rs` / `spin.rs` / `writers/qe.rs` 的条目（10-03 D 轻级的 job /
  CP2K / Gaussian 项、PBE0 待核实、QC 测试缺口，D-1）；英文 `//` 注释（10-17 规则已改）；09-27 第 2 条（第 6 步
  D3 已做）；已随 10-17 八步修掉的条目（chg-sdf `--qn` 默认值、traj 总览漏 bondlife、`is_lammps_data` /
  `Format`、`classify_trajectory`、`find_clusters` 的 MIC、近邻 O(N²) 与 `minimum_image` 求逆、`CellList`）。
  10-03 口径里 CP2K / QE 的生成口径（TTF slab、force/md 不请求应力、PBE0 用 COULOMB、bands 用 gamma、MD 步长）
  已对照 `ferro.jobs` 核过，**无一重犯**：非 3D 周期报错（`base.py:194`）、三任务模板都写 `STRESS_TENSOR`、`--timestep`
  对 QE 换算成 `dt`、模板无 HFX、无 bands / nscf 任务
- **这次不做**：R1–R6 逐函数审计（09-27）、code1 / code2 引用扫描（留在低优先级）

**用户定（2026-10-19）**：
1. S(q) 默认 `q_min` 不改（已有告警），Lorch 窗归新功能（低优先级），做的时候再查文献
2. map velocity / force 的空体素仍写 0，手册写明
3. `EV_TO_KCAL_MOL` 改用 CODATA 2018 的值
4. `merge --mode by-source` **每个体系内部也按 `--seed` 打乱**，再切 set（破坏性改动，改手册 `merge.md`）
5. 删 `trajectory_analysis.rs` 和 `geometry.rs`（零调用；MSD 不解包裹，距离 / 角度不取最小镜像）
6. dump：`element` 列优先，没有时用 `typelabel`；`ix iy iz` 不接
7. QE 的 `starting_magnetization` 不再存进 `magmom`（它是比例，不是 μB）
8. **不再支持 VASP4**：POSCAR / CHGCAR 没有元素行时报错；xyz 首列是原子序数时转成元素符号
9. **CP2K `MULTIPLE_UNIT_CELL` 按 CP2K 的规则支持**（`cell_methods.F:556`、`topology_multiple_unit_cell.F`）：
   `&CELL` 只放大晶胞，`&TOPOLOGY` 只复制坐标；后者不是 1 1 1 时两者必须一致，否则报错；
   复制顺序 k 外 j 中 i 内，块内保持原子原序。`SCALED` 坐标相对哪个胞，实现前读源码确认
10. **分箱统一按 MDAnalysis 的做法**（gr / vanhove / angle）：`n = round(range / dr)`，实际宽度取 `range / n`，末 bin
    止于 `r_max`，文件头写实际宽度。gr 默认 `r_max` 来自盒子，默认产物的 r 网格会微变（破坏性改动）。依据：MDAnalysis
    `InterRDF` 与 ASE `get_rdf` 用 bin 个数 + 范围，天然整除；GROMACS 的 bin 数取四舍五入，末 bin 可能越过 `r_max`；
    ferro 原先 gr 向下取整、vanhove / angle 向上取整
11. CN(r) 仍标在 bin 中心（与 MDAnalysis `get_cdf` 同口径；GROMACS `-cn` 标在上沿，二者差 dr/2），手册写明
12. xyz / extxyz 末帧截断时对齐 dump：丢掉末帧并告警，告警点名帧号

### 各批执行中的用户定与被推翻的原计划（摘自原分批表，提交号见 git log）

- **第 1 批（io 读侧）**：顺带修了一条同类的静默错误：PDB 用 `END` 分隔、不写 MODEL 的轨迹（CP2K / VMD 的写法）
  原先会并成一帧，现在按 `END` 分帧（同 ASE）。`MULTIPLE_UNIT_CELL` 对 5495 份 CP2K regtest / benchmark 输入实测无
  panic；写成 `${NREP}` 的预处理变量报错（不展开 `@SET`）
- **第 2 批（io 写侧）**：用户定：extxyz 的 `{:.10}` **不改**（ASE 写 8 位小数），手册写明精度；部分原子带电荷 /
  磁矩时**报错**，`convert --drop-partial` / `ferro.write(drop_partial=True)` 强制通过时整列不写（LAMMPS data 的
  电荷列必须有，全部写 0）。实现上与计划的出入：同一规则也用到 dump 与 data writer（三处同一变化原因，共用
  `writers::partial_column`）；`ferro_io::write` 改收 `&WriteOptions`，与 `ReadOptions` 对称
- **第 3 批（io 布局）**：用户定：**单位名 → 单位的对照表下沉 core**（`LengthUnit` / `EnergyUnit` /
  `PressureUnit::from_name`，固定不变的事实），提取规则（方括号、数值前一个 token、QE 的 `alat` / `crystal`）仍留在
  各 reader，cp2k_md / sp 补互指注释；`TableFormat` 单变体**保留**（公开 API、扩展点），注释记下。实现上与计划的
  出入：LAMMPS dump 的 reader / writer 自带 `0.04336410` 的 kcal 常数（布局建议 2 点名的那处，计划时被截断的 grep
  漏看，误报为已改），并入 3.6 改调 `units.rs`
- **第 4 批（analysis 数值与口径）**：分箱规则（bin 数向上取整、宽度收窄止于上界、整除时原样用步长）写进
  CLAUDE.md；sq 保持 g(r) 从 0 积分、默认 10.005（用户定，不为整格改默认值）。4.2（用户选改代码）：平滑窗由实际
  7 格 = 0.14 Å 改为记录所定的 5 格 = 0.10 Å；测试玻璃上 Zn-O / Al-O 的自动截断移动一两格（0.02–0.04 Å），P-O 与
  干净空隙不动
- **第 5 批（analysis 布局）**：用户定：map density / velocity / force / radius 选不到原子改为报错（原写全 0 格点、
  「0 atoms」）；**用户定：暂不换 TOR**（NPT 下 vanhove / map jump 仍是格点视图，待用户了解后再定）；
  `unwrapped_frac_sum` 边读边求和、互指注明不合；sq 的两个包装函数与 `elem_z_local` 各有 2–3 处调用，按 R1 保留
  （与审查建议不同）
- **第 6 批（CLI / render / structure）**：6.1 `\frac` 后紧跟因子（字母数字、`(` `[` `{` `_` `^`、命令）时整个分式加
  括号，后面是空格、右括号、运算符时不加。6.3 **python/ 子页不收**：doc.rs 既有说明写明只由 `ferro-py doc` 收录
  （两处收录会走样），与计划时的推荐相反，按既有设计。6.5 **与计划不同**：帮助页不改走 stderr（要改 25 个 print
  函数），只打一行错误指向帮助页。6.7 `build_box` 加 `seed`（ChaCha12，同 ml/merge）。6.8 **用户定**：by-source
  每个体系内按 `--seed` 打乱再切 set。6.11 采用推荐 (a)：`collect -o sys` 单组写 `sys/<目录名>.db/`，与多组一致
- **第 7 批（文档与脚本）**：7.9 **用户定**：installation.md 的命令表删除，改为指向 `ferro` 总览。7.14 **用户定**：
  ferro-cli crate 级 `allow(rustdoc::broken_intra_doc_links, rustdoc::invalid_html_tags)` 并写明原因（clap doc 即
  `--help`）
- **第 9 批（10-20 扫描，io 读侧）**：**用户定**：截断后一帧不剩报错（四个 reader 共用 `util::drop_truncated`）；END
  分帧的末块缺 END 但原子数一致时保留并告警。实现上与计划的出入：B-2 的判据收窄为「最后一行列数不足或数字读不出」
  才算截断，整行完整但内容错（原子序数 200）恒报错 —— 计划时没写，首版把它也当截断，测试当场抓到
- **第 10 批（io 写侧）**：**用户定：报错帧号统一 0 起**（同 extxyz / dump / deepmd / `partial_column`，胞的两条原为
  1 起）。被拒的 DeePMD 重写保住原系统（set 边界、逐项全有全无、力个数、virial 全部算完才动目录）。顺带注意到
  未修：`write_deepmd_npy_bounds` 的边界检查只比总帧数，重叠或越界的 bounds（如 `[(0,3),(5,7)]`，共 5 帧）能过
  检查后在切片处 panic；只有 CLI 内部生成 bounds，暂无外部入口
- **第 11 批（输出口径）**：**用户定**：与文件无关的消息以程序名（`ferro` / `ferro-py`）作来源；顶层致命错误也统一
  成一行 `ferro: FerroError: {e:#}`（anyhow 的 `Caused by:` 链改用 `: ` 连接）；多行消息合并成一行；产物表头里的
  `warning:` meta 行（`ml/diagnostics.rs`）是文件内容不是通道，不动。ferro-io 的 `Notice { level, line, message }`
  只有 io 产出、调用方都在顶层，按分层铁律留在 io。ferro-py 的 `_show_warning` 与 `report.rs` 是同一口径的两份
  实现，互指不合（R6）
- **第 12 批（core / structure）**：12.2 **用户定**：缺质量时 `ferro info` 不打 Density 行、其余照出、FerroError
  点名并退出码 1；inspect 的 density 列写 NaN、`#` 块写 `error:`。12.10 A-11 **用户定报错**：判据留 1e-3 相对容差；
  实测水与 $Ca_3(PO_4)_2$ 盒子（15–1500 原子、三个种子）500 步内全部收敛、100 步不够，文档建议步数 100 → 500。
  12.7 `EV_TO_WAVENUMBER` 旧值 8065.5440 对不上任何一版 CODATA（出处不明，注释如实写）。12.5 A-16 判为不合，三处
  互指。**用户定（2026-10-20）**：非 dump 内容的 `.lammpstrj` 读成 0 帧**不是问题、不修**（不会有人把别的格式写成
  `.lammpstrj` 扩展名）
- **第 13 批（analysis）**：用户评：C-1 / C-3 严格说不算问题（曲线固定，间隔够细不改结果），同意修。13.1 实现上与
  计划的出入：`dr` 列走 `Summary::note` 写全精度，不走数值列 —— `fmt_compact` 只留四位，会把 0.0019994 显示成
  0.0020。13.7 `PeriodicCell::to_fractional` 与 `cartesian_to_fractional` 同为 `matrix.transpose()` 的逆，计划时担心
  的 1 ulp 换格不存在。13.10 C-12 **用户定只改注释**。13.11 C-13 **用户定混用即报错**：第一个缺速度 / 力的帧点名
  报错（同 vacf）
- **第 14 批（cli / render）**：**用户定**：D-1 filter 照抄 convert 的三条检查、互指（两处不到 Rule of Three）；D-10
  `--method` 改 `ValueEnum`，拼错由 clap 报错、退出码 2。`ferro doc` 管道输出不渲染，CLI golden master 覆盖不到渲染，
  要用伪终端（`os.openpty` + `PAGER=cat`）比对。新增 `ferro_io::read_charge_grid`（网格格式不进 `Format` 表）
- **第 15 批（python / scripts）**：15.2 E-2 **用户定**：QE vc-relax 要求模板写 `ecutwfc`（`ecutrho` 照旧兜底）；CP2K
  的 `--vary-cell` 与 NPT_I / NPT_F **必须用自己的模板**且写 `&MGRID CUTOFF` —— 内置模板写死了 `CUTOFF 400`，只查
  「有没有」形同虚设，故原样的内置模板（含 `job template` 导出未改）也拒。**15B 用户定**：只迁「出图」进 `ferro.plot`，
  加 `ferro-py plot <gr|sq|angle|msd|net>`；对拍（`ferrocmp` / `compare_*` / `trajcheck`）与空间群工具留在 `scripts/`。
  形态：`-o` 恒为目录、原脚本的 `-o STEM` 改 `-s`；可选依赖 `ferro[plot]`，`import ferro` 不加载 matplotlib；
  样式改用 context manager，库调用不改调用方的 rcParams；`suffix_of` 与 `composition_of` 并排并互指
  `batch::product_name`（Python 侧无绑定，不合）
- **第 16 批（文档）**：**16.15 用户定（F-20）：Python `gr_pair` / `gr_all` 默认改为 CLI 的 `dr=0.002`、`r_min=0.001`**

## 2026-10-20 全库扫描

**用户定（2026-10-20）**：
1. E-2：QE / CP2K 的变胞任务也要求模板写截断能（Pulay 应力），同 VASP
2. E-1：VASP 的 PREC **只收 `Normal` 与 `Accurate`**（High 等已淘汰），其余写法报错
3. 输出口径统一为 Python 风格 `<来源>[:<行号>]: <类别>: <消息>`，一律走 stderr，三个类别合并管理：
   `FerroNote`（说明）、`FerroWarning`（可能有问题）、`FerroError`（错误）。例：
   `traj.pdb:57: FerroWarning: last MODEL (frame 2) has no ENDMDL; kept`、
   `sets/a.db: FerroNote: no virial.npy; stress criterion switched off`
4. A-3：查不到质量的元素名**报错**并点名，由用户改输入（不回退 1 amu、不按前缀猜）
5. C-2：**只登记**，有 Python 入口时再补（见下）
6. 全部等级逐条修；每批开工前逐条计划交用户审，每条单独提交并带回归测试

## 2026-10-17 全盘架构审查

六个子代理按领域审：A 分层与依赖 · B ferro-core 数据模型 · C ferro-io · D ferro-analysis · E ferro-cli ·
F ferro-python。已登记的条目不重报，延伸处注明原编号。复现 fixture 是一次性的（会话 scratch），未进仓。
**主会话亲自复现 / 读码确认**：A-7（pbc，见下）、C-6、C-3、A-11。

**总体结论**：分层健全，不推倒。中间层无横向依赖（Cargo.toml 与 `use` 都核对过）；analysis / structure /
render 不碰文件系统，**没有功能模块绕过 core 直接调 io**；CLI 内无格式解析、几乎无物理计算；Python 与主体库
分工基本守住，绑定层是薄封装。问题集中在：core 收了不该收的算法、同一事实散在入口层、告警与格式识别
各有多套、内存模型。**不重写**：`batch.rs` 的泛型形状、`md/correlate.rs`、gr 的 NPT 逐帧归一、`ml/merge`、
vasprun / OUTCAR / cif + spacegroup / extxyz 应力、`ferro.jobs` 的 `base.py` 三层。

**用户定（2026-10-18）**：
- 全部登记；按下面的执行顺序逐步做，**每一步的执行计划先交用户审核**，批准后才动手
- **`Atom` 存储方式一并改**（不止记方向）
- **电荷保留**。来源只有：LAMMPS data 的 `full` / `charge` 类型、lammpstrj 的 `q` 列（逐原子 `Atom.charge`），
  以及 Bader 分析产出的电荷。supercell / merge 没有 CLI 入口（只经 Python 的 structure 绑定）。
  `Frame.charge`（整数总电荷）的具体形态在第 5 步的执行计划里定
- **删 `Frame.multiplicity`、`Frame.bonds`、`TrajectoryMetadata.timestep`**（无生产者）
- 降级一条：E2-2「`filter -s 5` 静默当 5 GPa」—— 帮助页本来写的是 `-s, --s-max <GPA>`，属与全仓
  「`-s` 是批次后缀」的约定撞名，不是静默出错，归第 6 步的命名统一

第 1 步的两条用户定：**除 lammpstrj / xyz / extxyz / pdb 这类轨迹格式外都是单帧**（C-3：库 writer 遇多帧
报错，Python `ferro.write` 与 `ferro convert` 一样一帧一个文件，命名函数 `indexed_path` 下沉 ferro-io 共用）；
LAMMPS data 同一元素的原子质量不一致时 writer 报错（一个 type 一个质量）。

第 2 步的用户定：**非三维周期体系一律报错**，不按轴支持、不另定 gr / sq 归一化口径（「理论上不会有
非三维周期体系」）；pbc 只在 `PeriodicCell::new` 一处判，FFF 与 TTF 同等对待。structure 的
`box_builder` 不并入 `NeighborList`（R6，两处注释互指）。随之 `NeighborList` 的使用方只剩 core 的
网络分类与 analysis：第 5 步 A-1 把网络分类搬进 analysis 时，`NeighborList` 跟着搬，`PeriodicCell`
留 core（io / structure 以后量距离也要它）。挂起组没动内部：cube_sdf / chg_sdf 只加入口检查（截断上界），
Bader 仍用 `Cell::minimum_image`。

第 8 步的用户定：C-11 sp 向 md 对齐（锚点上方最近的 SCF 行说了算，没有就保留）；D-2 **不合**
（cp2k_md / sp 不共用代码，三行函数不抵跨文件依赖），四处互指；D-3 **合**（CHGCAR 头部就是 POSCAR，推翻
布局建议 10 的「不合」），`vasp.rs::parse_poscar_header`；HBAR2_2M 先查再改：VASP 的 HSQDTM =
RYTOEV·AUTOA² = 3.8100199（RYTOEV 经 VASP 论坛版主确认，AUTOA / HSQDTM 的定义取自复述 constant.inc 的
VaspBandUnfolding，未对照 VASP 源码本身），随机 2 万个立方胞里 3 个网格变化。实现上与计划的出入：CHGCAR
负缩放因子的报错文案随 POSCAR 统一（计划写「各自保持」，为此多传一个参数不值）；clippy.toml 的路径写
`nalgebra::Matrix::as_slice`（变异验证生效，数组 / Vec 的 as_slice 不误报）。仍开着的只有挂起组：
Bader（A-3、B7、B8）与 cube_sdf / chg_sdf 两个改写组。

第 7 步的用户定：D1 认不出元素的 UPF **跳过**（缺赝势时报错里点名被跳过的文件）；D2 `HBAR2_2M`、
`RY_TIME_FS` 是程序侧的量，**留在 Python** 并注明理由，`RY_TIME_FS` 改用 QE `constants.f90` 的精确值；
D3 E2-4 下沉的 dataset 函数**这次不做** Python 绑定（无调用方）；D4 `qe_structure_cards` 返回
`[(卡片名, 原文)]`，破坏性变化接受。实现上与计划的出入：`kmesh` 不改签名，新绑定是模块函数
`ferro.interplanar_spacings(cell)`（收 3×3 胞，测试照旧喂矩阵），不是 `Trajectory` 的方法；gr 的
`py.detach` 放进共用的 `run_gr` 里经 `Python::attach` 取 token（给 `gr_pair` 加 `py` 参数超出 clippy
的 7 参数上限）。**留给第 8 步**：`readers/cp2k.rs` 自带 `BOHR` 常数与 `units.rs` 重复；`ferro-py
--start -1` 等负数选帧参数抛 `OverflowError` traceback（argparse 收 int、绑定收 usize，改前即如此）；
`HBAR2_2M = 3.80998` 是 CODATA 的 ħ²/2mₑ，按记忆 VASP 自己的 HSQDTM = RYTOEV·AUTOA² ≈ 3.81002，
预测 NGX 若要与 VASP 逐格一致应取后者 —— 未查源码核实，待用户定是否查。

第 6 步的用户定：D1 bader 跟命名规则（`bader_ACF_<input>`）；D2 merge 的 `--suffix` 改名 `--dir-suffix`；
D3 map 单输入也带输入名；D4 帮助页全部统一列宽；D5 sdf 段序 `sdf[_famN]_<type>_<input>[_<s>]`。
实现上与计划的出入：子步骤顺序改为先命名、后合并驱动（命名统一后 map 不再分单 / 多输入，驱动才能
直接合并）；net 导出按统计表的 mode 名落为 `network_types_<input>`，不是计划里写的 `net_types`；
公共参数块用 `const` + `{COMMON_OPTS}` 内联捕获，不用 `concat!`；net 帮助页搬进 help.rs；
合并驱动后 net 的 `Labels:` 块打在 `Inputs:` 之前。新断言上线即抓到 chg-sdf `--qn` 默认值写错、
traj 总览漏 bondlife，并发现 filter 页例子仍写 `-s 8`（6.4 漏改），一并修。

第 5 步的用户定：**删 `Frame.charge`**（连同 multiplicity / bonds / timestep；逐原子电荷不动，Python 的
`charge(frame)` / `multiplicity(frame)` 一并删）；`metadata.source` 拆成 **`title` + `program`** 两个自由文本，
writer 只回写 title（AIMD 输出转 xyz / POSCAR 注释行因此不再是程序名）；`calc_network` **遇无胞帧报错**，
net 只分类一遍（`calc_network_typed`）；E2-4 只搬纯计算（`Split`、`stitch_segments`、`frame_table`、
`filter_diagnostics`），Python 绑定留第 7 步；A-12 收口后暴露的零调用项逐条定：cube_jump **保留为公开 API**、
`read_cube` **标记保留、以后与 `read_cube_as_chg` 合并成一个 cube reader**（待办）、`VERSION` 移到
`ferro_analysis::VERSION`、`l_f_adj` 字段 / `has_qn` / compounds 的 density·cas·`with_density` 删。
实现上与计划的出入：B11 的 `check_same_periodicity()` 没加（「有胞 / 无胞一致」的双向检查只有 deepmd writer
一处，R1 不成立）；5.4 拆成两次提交；4a 的提交说明写「Python 的 vacf 等以前不查」不准确 —— Python 只暴露
msd，而 msd 本来就查，Python 侧行为没有变。

第 4 步的用户定：原子存储取**方案 C**（每帧 `Arc` 共享位点表 + 自己的坐标 / 电荷 / 力 / 速度数组，
不在 `Trajectory` 上立「各帧原子相同」的不变量）；逐原子数组**私有**、只经访问器；**分析侧逐帧消费**
（gr / sq / map / net 边读边算、内存与帧数无关）**这次不做**，记作后续项。实现上：「先迁访问器、再换表示」
走不通（AoS 下 `positions()` 给不出切片），4.3 一次提交；`Frame::permuted` 顺带落地了第 5 步 A-2 的置换本体
（`sort_atoms` 已改用它），A-2 余下的是让 CLI / Python 的调用点不再经 `ml::merge`。

第 3 步的用户定：`convert` / `info` / `traj` 新开放的读入格式为 **AIMD 输出与 DeePMD 目录**（cube、CHGCAR
不开）；**非 UTF-8 严格报错**（报行号列号，不替换坏字节）；cube_sdf / chg_sdf 与 Bader 都要重构，告警与返回
类型**随改写组挂起**。实现上与计划的一处出入：拆帧没有全并进 `write` —— convert 按输入帧号命名、Python 按
位置命名，于是 io 给 `write`（按位置）与 `write_numbered`（给定编号）两个入口，共用一份循环。

Bader 与 cube_sdf 两个改写组仍各自挂起，D 组给出的新结构记在下面「D 组」R-2 / R-3。
第 2 步与第 4 步都要改 `Frame` 的几何读取路径：第 2 步的原语接口按「坐标切片 + 胞 + pbc」设计，
第 4 步改 `Atom` 时不必再改一遍调用点。

B 组重设计：R-1 周期几何原语（= 第 2 步）；R-2 拓扑与坐标分离
`Trajectory { topology: Vec<Site{element,label,mass,charge,magmom}>, frames: Vec<FrameData{positions,cell,pbc,forces,velocities,energy,stress,…}> }`，
允许逐帧覆盖 topology（GCMC / DeePMD mixed type 天然对应），`check_same_atoms` 变成构造期不变量（= 第 4 步）；
低成本缓解是 dump 流式读 + 元素 interned。**不建议**：`Frame.stress` 强类型化；合并 `CubeData` 与 `ChargeGrid`
（两者分别对齐 cube 与 CHGCAR 文件序，合并会搞乱 Bader 下标；建议 doc 互指）。

E 组重设计：R-2 产物命名统一 —— 现有四种排列（traj `<mode>_<label>_<suffix>.csv`；map `<mode>_<suffix>_<input>.cube`，
实测 `density_x_prod.cube`；net 导出 / bader `<input>_…_<suffix>`；sdf 单输入时 `-s` 直接替换前缀，实测 `x_Ob.cube`），
建议定一条 `<mode>[_<表>][_<选择>][_<输入>]_<suffix>.<ext>`，map 恒带输入名（延伸合规扫描 2、10-03 D 轻）。
破坏性，`scripts/ferrocmp.py product_name()` 与手册同步。R-6（低）：`bader -i` 改收多文件走 `map_inputs`。

## `ferro.jobs`：输入生成迁到 ferro-python（2026-10-16 定案）

生成输入不吃性能、改动最勤，Python 改完不用重编、用户也能直接改。2026-10-15–16 grilling
定案（对照 ASE `calculators/cp2k.py`、cp2k-input-tools 0.9.1、pymatgen `io.core` 的设计与依赖）：

| 方面 | 决定 |
|---|---|
| 定位 | 简易产生通用计算输入；复杂任务交给专业工具。留功能扩展点（固定原子、Gaussian、初始速度暂不做） |
| 范围 | VASP / QE / CP2K × SP / opt / AIMD；只做周期体系，无胞报错 |
| 位置 | ferro-python 的 maturin 混合布局：`import ferro`，`ferro.jobs`，命令 `ferro-py`（子命令 `job`，以后 `plot`） |
| 依赖 | 运行时**只用标准库**，Python ≥ 3.11，abi3 wheel；pytest 只作开发依赖。ASE / pymatgen / cp2k-input-tools 均有 scipy、matplotlib、pydantic 等硬依赖，超算离线难装，不依赖 |
| 分工 | 结构块（POSCAR、QE 三张卡片、CP2K `&CELL`/`&COORD`）由 Rust 写成字符串，格式规范只一份；参数块（INCAR、namelist、CP2K 节树）在 Python。QE namelist 解析由 Rust 暴露，全仓一份（同一语法、同一变化原因） |
| 模板 | 模板即原生输入文件；内置 9 份，用户模板目录同名覆盖；未知键原样保留（不提示：三个程序都没有可靠的全量键表）。2026-10-16 用户定的用法：`ferro-py job template -s … -t … -o DIR` 导出通用模板（头注释写明程序改写哪些键、哪些须自己加），用户改完以 `--template` 读回再合并结构；命令行显式给的量覆盖模板，没给的沿用模板 |
| 默认值（2026-10-16 用户定） | 内置模板泛函 PBE-D3(BJ)；AIMD 1 fs × 2000 步；k 点 kspacing 0.5 Å⁻¹，VASP 写成 INCAR 的 `KSPACING`、不生成 KPOINTS，网格 1×1×1 时提示可用 vasp_gam；SCF 不收敛照跑（CP2K `IGNORE_CONVERGENCE_FAILURE`）只在 aimd；并行参数内置模板不写，模板头注释与生成时提示用户自己加。其余默认值在模板导出后由用户对着模板审 |
| 系综 | VASP、QE 只做 NVT（QE `vc-md` 不支持热浴；npt 报错）；CP2K 支持 NVT / NPT_I / NPT_F，`--ensemble nvt\|npt_i\|npt_f`，NPT_F 的 `VIRIAL` 不开放（默认全部方向，要单方向改模板），NPT 写 `STRESS_TENSOR ANALYTICAL`（D-M4） |
| VASP meta-GGA（2026-10-16 用户定） | 模板写了 METAGGA 时 ALGO = All + ISEARCH = 1，并显式写 NGXF/NGYF/NGZF = `--fine-grid` × 粗网格（默认 2；含 Zn 可能要 3，只提示）；MAGMOM 等结构相关数组键不检查（氧化物玻璃用不到）。`--fine-grid` 一律照做只提示：PREC 非 Accurate 时粗网格按 3/2·G_cut 算，GGA 也照写 NGXF，CP2K 打开 USE_FINER_GRID，QE 忽略 |
| QE 赝势（2026-10-16 用户定） | 不用 SSSP json，索引 `[qe] pseudo_dir` 里的全部 UPF，元素 / 价电子 / 建议截断取 PP_HEADER；一元素多份时 `[qe.pseudos]` 指定 |
| CP2K 方法（2026-10-16 用户定） | GPW（GTH）与 GAPW 全电子都支持，由模板的 `&QS METHOD` 决定；全电子基组 pob-TZVP-rev2 系列（不涉及四五周期以后的元素）、`POTENTIAL ALL`；用户发现 GTH-SCAN 部分元素有问题（已报开发者），全电子更安全；`GAPW_ACCURATE_XCINT` 写进模板、默认注释掉。用户用 CP2K 2025.2 与 2026.2 |
| 收敛测试（2026-10-16 新增） | `ferro-py conv` 生成一组单点（每档一个子目录），`ferro-py conv analyze DIR` 读回输出，以最高档（或 `--ref`）为参考，报总能量偏差（meV/atom）、压力偏差、最大力偏差、逐原子力偏差的 RMS；支持多帧（选帧同 job），每帧每档一行 CSV，终端表取各帧最大。阈值只作提示，取值查惯例。VASP 先测细网格 NGXF/NGYF/NGZF 再测 ENCUT（r2SCAN 的经验，PBE 未定）；CP2K 测 CUTOFF、REL_CUTOFF 与 `USE_FINER_GRID` 开 / 关（不扫 FINE_GRID_FACTOR；开细网格后 CUTOFF 可能不必那么高，但不一定所有元素都要细网格），GAPW 时加 `GAPW_ACCURATE_XCINT` 开 / 关；QE 的单点能量 / 力 / 应力提取与截断能测试放下一批。读输出与写表在 Rust（复用 OUTCAR / CP2K reader，暴露 `write_table`），偏差算术在 Python |
| 程序改动的量 | 结构、电荷、多重度（奇偶校验）、k 点（按 kspacing，模板写死的不动）、元素列表；AIMD 加温度/步数/时间步/系综；opt 加是否变胞 |
| 自旋 | 默认按奇偶取 1 / 2；显式值与奇偶矛盾报错。未成对电子的方法未定，`guess_spin` 暂不暴露 |
| 配置 | `~/.config/ferro/config.toml`（`FERRO_CONFIG` 覆盖）：POTCAR、CP2K data、QE 赝势、模板目录；运行时读 CP2K 自带数据文件，不迁 `cp2k_basis_db.rs` |
| 批量 | 选帧沿用 `convert`（默认全部帧）；`-o` 指目录，每帧一个子目录；同组成共用一份 POTCAR + 相对符号链接，组成不一致报错 |
| LAMMPS | dump reader 加 `--type-map`；剩 `X*` 占位元素时生成任务报错 |
| convert | QE 只读不写（破坏性） |
| 默认参数 | 以用户脚本（`private/` 下 `dp_makeAIMDin.py`、`make_NMR.py`、`make_vasp_label.py`，覆盖 VASP / CP2K 的 AIMD 与 SP）为准；opt 与 QE 无参考脚本，查惯例补齐，出入处商量着定 |
| 验收 | CP2K 本机用 cp2k-input-tools 离线校验；9 个算例在超算冒烟（CentOS 7，glibc 2.17 → manylinux2014 x86_64，GitHub Actions 编） |
| 收尾 | 追平后删 ferro-workflow 与 `ferro job`（✅ 2026-10-17） |

- **CP2K opt**：2026.2 删了 `MOTION/CELL_OPT/TYPE`，写了即在读输入时 ABORT。`&CELL_OPT` 节总在
  模板里，所以不变胞的 GEO_OPT 也中止。旧版默认值即 `DIRECT_CELL_OPT`（cp2k-input-tools 自带的
  CP2K 9.0 关键字表核实），删掉对 2025.2 不变。**以后写 CP2K 模板不要写这个键。**
- **QE aimd**：第 2 步 `checkallsym` 中止 —— SVR 热浴的随机速度破坏初始对称性；aimd 模板加
  `nosym = .true.`。
- 编译器冲突：VASP / QE 要 intel，CP2K 要 gcc + openmpi，同一作业脚本里装不下 → 拆两个作业，
  用 `SKIP` 分程序、各用独立 WORKDIR（脚本开头会清空 `smoke_report.txt`）。脚本不用改。
- VASP OUTCAR 有一条「KPOINTS 文件找不到，回退自动生成」警告 —— 设计如此（k 点走 INCAR 的
  `KSPACING`），不是问题。收敛测试只生成未跑（未设 `RUN_CONV=1`）。
- 已确认（2026-10-17，用户查 2025.2 手册）：CP2K **2025.2 没有** `GAPW_ACCURATE_XCINT`，只 2026.x 有。
  模板与 `job-cp2k.md` 本已注明 2026.x；`conv.md` 的 `accurate-xcint` 一行补注。

用户已审过并认可 VASP / CP2K / QE 的 9 份内置模板；push 每次都要先问。

**2026-10-17 QE 这一批**（用户定：`ferro dataset collect` 也收 QE 的单点与 AIMD 输出；QE 的 AIMD / SP
主力仍是 VASP、CP2K，本地只有 GIPAW 流程的产物，无带力的 pw.x 输出）：

- ✅ `ferro-py conv -s qe --scan ecutwfc=…`；ecutrho 来自 UPF 建议且低于 4 × ecutwfc 时提到 4 倍并提示，
  模板写的照用只提示，≤ ecutwfc 报错（判据查 QE 7.3.1 `PW/src/input.f90` `set_cutoff`：dual ≤ 1 才
  errore，< 4 只 infomsg）。`conv analyze` 读 QE 等 reader
- ✅ **超算算 fixture**（2026-10-17，QE 7.5，intel 2025.1.3 + impi，16 进程）：
  `ferro-python/scripts/qe_fixtures_hpc.sh`，结构 `exp_603719.cif`（6 原子 AlPO₄），五个算例全部
  `JOB DONE`，放在 `examples/qe/pwx/{sp,sp_spin,sp_noconv,md,vcmd}/`（pw.in + pw.out；vcmd 另有
  pw.x 自己写的逐步文件 `pw.avec` 等，可核对逐步晶胞），作业脚本与日志在 `job/`。首跑全部
  GENERATION FAILED：作业脚本的 `conda activate` 前缺 `conda init`（用户修好重跑）
- ✅ pw.x 输出 reader `readers/qe_out.rs`（`04e95a0`）+ `aimd.rs` 嗅探；collect `--format qe/pw`
  （`49befab`）；`read_output` / `conv analyze` 读 QE、`pyerr` 带出原因链（`9d9c8bc`）。核对：能量、力
  与 dpdata 一致到 1e-10，其余与 ASE 一致到 1e-7（ASE 能量的 1.2e-6 相对偏差是它写死 CODATA 2006）；
  collect 产物经 dpdata 读回逐帧对上。变异验证三项（读第二组力、应力翻号、MD 帧错一位）都被抓到。
  **与原计划的两处出入**：① 原计划「没开 tprnfor 的帧照收」，实测 collect 会写出没有 `force.npy`
  的 system，改为力必需（同 cp2k_sp）；② 不读 md 每步的温度 —— pw.x 打印的温度属于哪一帧
  有歧义，宁缺勿错。没做 relax / vc-relax（报错），也没有它们的样例

**2026-10-16 ferro-python 审查**（Opus 子代理只读审查，19 条 + 6 条存疑；复现脚本未入库）。

用户裁定（以后照此，不再重提）：
- PREC 写全称不缩写 → 缩写报错；不同组成不得输出到同一 `-o`；**不支持 `GAPW_XC`**（配赝势时
  力有错，用户有限差分验证过）；`--vary-cell` 不替用户定 ENCUT（1.3 × ENMAX 不够），模板未写即报错
- VASP 流程：先固定 ENCUT 测细网格，再带定好的 `--fine-grid` 测 ENCUT（粗网格随之变）
- r2SCAN 下 VASP 默认开 LASPH、LMIXTAU，不用写；CP2K 2025.2 不必写 `GAMMA_CENTERED`，`MONKHORST-PACK`
  只在小胞时打开（现行为：1×1×1 不写 `&KPOINTS`），偶数网格不以 Γ 为中心不算问题
- CP2K 基组与赝势须同 q 配对；GTH-PBE 可不写电荷，GTH-SCAN 必须写明（生成一律写 `-qN`）

已修（各一提交，pytest 125 → 133）：P-1 `ffc5c9f`、P-2 `40535de`、P-3 `e23ba40`、P-5 `cec54a1`、
P-6 `80f9646`、P-7 `511820b`（强制 `ibrav = 0`、删胞参数键）、P-9 `6d90a19`、P-10 `706cc64`、
P-11 `a282757`。P-4（ecutrho）、P-8（CP2K 偶数网格偏离 Γ）按上面裁定不改。

`RECOMMENDED_POTCAR` 个别条目（Fe、Er…）未对照 vasp.at 当前表：**用户定不核对**（2026-10-17）。

**2026-10-17 `ferro dataset` 审查**（子代理只读审查：原子重排后力、坐标、盒子、能量、应力是否逐帧
逐原子对应；对照 dpdata 1.0.2 环境，构造可反查身份的数据跑 collect / filter / merge / convert）。
结论：**无逐帧、逐原子错配**。编号只在本节内有效：

| 编号 | 级别 | 问题 | 处置 |
|---|---|---|---|
| D-1 | 高 | dpdata 的 nopbc system（有 `box.npy`）经 filter / merge 丢 `nopbc`，团簇变周期 | ✅ `cfb228d` |
| D-2 | 中 | collect 用第 0 帧置换套全部帧 | ✅ `7f5fa81`（`sort_atoms` 逐帧；实际格式触发不了） |
| D-3 | 低 | filter 保留帧号、`--shuffle` 置换不落盘，产物无法追溯回原始帧 | **用户定不做** |
| D-4 | 低 | `fparam` / `aparam` / `spin` 等额外 npy 只告警即丢 | **用户定不做** |
| D-5 | 低 | filter 不重排原子（collect / merge 会） | 未做，可写进手册 |
| D-6 | 低 | extxyz 无 `masses` 列时 momenta → 速度用 ferro 质量表（P 相对误差 7.7e-6） | 可忽略 |
| D-7 | 低 | 超过 1000 个 set 时目录按字典序 | 未做 |

测试缺口已补（`673d732`）：非对称 virial 的盘上布局、merge 打乱与 filter 筛帧的逐帧逐原子对应。
仍存疑：左手系晶胞。CP2K MD 的 virial 已核对（2026-10-17，见 `progress.md`「`cp2k_md.rs`」）。

## 2026-10-11 全库审查

六个子代理按领域只读审查：A `ferro-io` · B 结构分析（gr/sq/angle/cube_*）· C 动力学（msd/vacf/
rotcorr/vanhove/bondlife）· D network/dft/ml/core · E cli 驱动/render/structure · G 横切（R1–R6、
死代码、手册一致性、scripts/、dev/ 漂移）。`ferro-python`、`ferro-workflow`、`cmd/job.rs` 按用户
要求**不在范围内**。对照 deepmd 环境的 ase 3.29 / dpdata 1.0.2 / MDAnalysis / pymatgen 源码与
`examples/external/cp2k-2026.2/src`，可联网查文献。复现 fixture 是一次性的（会话 scratch），未进仓；复现
写到可重做的程度。A-18 与 G-6 是同一条，合并。四条登记勘误（G-2、G-3、G-4、G-15）已直接改在
`plan.md` / `progress.md` / `issues.md` 的原条目里，G-3 并入零调用清单的判断。

**本轮核对一致、不必再查**：Waasmaier–Kirfel 形状因子（对 ASE；ASE 自己 Cu 常数项、O 一个系数
有误）与 pymatgen（Sn 是 pymatgen 错）；中子散射长度 99 项对 pymatgen（只 B 不同，有意 ¹¹B）；
Faber–Ziman 权重；msd / vacf / rotcorr / vanhove / bondlife 对合成轨迹解析解与 numpy 独立实现
（csv 精度 5e-8 内，含 metal 单位带速度 dump 的 vacf C_v(0) 与 Green–Kubo 平台）；changelog
`Unreleased` 覆盖 0.3.5 全部破坏性改动；手册 LaTeX 约定零违规；PAGES 覆盖全部手册页；
traj / net 产物文件名与列名与手册一致；绘图脚本对当前产物可跑通；分层、`as_slice`、`Molecule`、
`label()` 反解析零违规；collect → filter → dpdata 读回一致。

- **cp2k_md / cp2k_sp 的五个重复函数**（A-18 / G-6）：`version_note`、`unit_in_brackets`、`energy_unit`、`energy_to_ev`、
  `stress_unit`（`cp2k_md.rs:122-215` 与 `cp2k_sp.rs:81-153`），报错文案已开始漂移。**注意 2026-09-22 用户明确要求两者
  「不共用单位与数据抽取」**（`overview.md`），故不合；按 R6 补两处互指注释，并把已漂开的文案修平

## 2026-10-03 四路审查

四个 Fable 5.1 子代理按领域只读审查（A 格式读写 · B 轨迹分析 · C network/dft/ml/core ·
D cli/structure/workflow/跨 crate），对照库在 `~/.miniforge3/envs/deepmd`（ase 3.29、
dpdata 1.0.2、MDAnalysis 2.10、pymatgen 2026.7.16），CP2K 口径对照 `examples/external/cp2k-2026.2/src`。
复现 fixture 是一次性的（会话 scratch），未进仓；复现写到可重做的程度。
**「严重」级 21 条主会话已用 debug 版 `ferro` 逐条复现**（B-2 只复现了 sq 一路）；
中 / 轻级与布局建议是子代理实测或读码的结论，标「推断」的未实测。
已在本文件或 `issues.md` 登记的条目不重复。

**已对拍一致、不必再查**：g(r)/CN（含三斜 27 镜像暴力、NPT 逐帧归一）、S(q) 各 partial 与
XRD/中子加权、msd（对 MDA EinsteinMSD fft）、vacf、rotcorr、vanhove、bondlife、angle、
`map density` 积分；OUTCAR（2000 帧）、vasprun（425 帧）、dump（NPT + 三斜）、CP2K 单点
（对 dpdata）、CIF 对称展开（三个空间群）、extxyz 应力往返、CP2K `.inp` reader 跑 5502 份
regtest 零 panic；超胞（对 ASE / pymatgen）；元素质量（86 种）；Cell 运算；`net` 六张表
（对 numpy）；dataset collect/merge/filter 产物（对 dpdata 读回）。

**QC 输入生成组**（**2026-10-17 `ferro-workflow` 与 `cmd/job.rs` 已删，本组不再修**；保留作 ferro-py job 的不重犯清单，第 15–17 步已对照）：D-M4、D-M6–D-M12（CP2K / QE / Gaussian 输入生成）**挂起作专项**（2026-10-08 用户定：大工程，先记录）。动手前须逐条对照 `examples/external/cp2k-2026.2/src` 与 QE、Gaussian 官方文档核实——D-M10 等几条只是推断；修完应能让生成的输入真跑过一次各软件的输入解析。

| # | 位置 | 问题 / 依据 | 修法 |
|---|---|---|---|
| D-M4 | `cp2k.rs write_force_eval` | NPT 不写 `STRESS_TENSOR`，CP2K 启动报错（`md_run.F:331-344`） | NPT 时写 `ANALYTICAL` |
| D-M6 | `cp2k.rs --smear` | 不写 `ADDED_MOS`，CP2K 报错（`qs_environment.F:2247`）；`--scf ot` 时 smear 静默丢。手册 `spin.md:76` 示例正中 | 写 `ADDED_MOS`；OT + smear 报错 |
| D-M7 | `cp2k.rs` PBC | 分子 / `--pbc z` 不写 `POISSON_SOLVER`，默认 PERIODIC 报错；分子 `&CELL` 无 ABC。`test_energy_molecular` 断言的正是这份不可运行的输入 | 非 3D 写 MT / WAVELET + ABC，或报错 |
| D-M8 | `cp2k.rs` | `--kpoints` + `--scf ot` CP2K 报错（`qs_scf_initialization.F:886-887`），builder 不拦 | 提前拒绝 |
| D-M9 | `ferro-workflow/src/qe.rs build` | 无胞仍写 `ibrav = 0` 不写 `CELL_PARAMETERS`，pw.x 拒收；`test_scf_basic` 断言的正是它 | 报错或要求给盒子 |
| D-M10 | `qe.rs` vc-md（推断，本机无 QE 源码） | 写 `ion_dynamics='verlet'`，INPUT_PW 中 vc-md 只支持 `'beeman'` | 改 `beeman` |
| D-M11 | `ferro-workflow/src/job_builder.rs` Gaussian | 周期结构丢胞、不写 `TV`，静默变团簇计算（ASE 写三行 TV） | 写 TV 或报错 |
| D-M12 | `cmd/job.rs` | 与所选软件无关的参数静默忽略（QE 不读 `--cutoff`/`--md-timestep`/`--thermostat`/`--pbc`）；`--task` 等枚举值读完输入、打印 auto-spin 之后才校验 | clap `ValueEnum`；不适用参数报错 |

- ~~`merge` 侧向是插真空条带而非 ASE `stack` 的应变匹配~~（2026-10-09 随 D-M2 定：保留不拉伸，
  共格界面不是它的用途）；`box_builder` 按原子随机放、软约束
  （SiO₂ ×200 实测 241 对略小于 `min_dist`）

10. **决定不合、需两处互指注释**（R6）：`cell_to_lammps`/`lammps_cell_matrix`/`bounding_box`
    在 lammps_data 与 lammps_dump writer 逐字两份（**互指注释已随 `a34b488` 加上**）；POSCAR 头部
    解析在 `vasp.rs` 与 `chgcar.rs` 两份且已漂过一次（`chgcar.rs:71` 仍是「`d` 开头才算 Direct」
    旧规则，与 A7 一起修）；msd 与 bondlife 的逐帧 (Mᵀ, Mᵀ⁻¹) 预计算；`filter_split`/`merge_split`；
    `job.rs:204` 与 `convert.rs:155` 的 `-o` 尾分隔判断；`vacuum.rs`/`merge.rs` 轴名解析；
    rayon `build_global` 两处；io 与 workflow 两份 QE writer（分层所限不能合，先修平 nspin 漂移）。
    ferro-structure 4 处手写完整 `Frame { … }` 字面量**有意不合**：新增字段时编译器逼每个构造点表态

## 2026-10-01 全库审查

子代理（Fable 5.1）只读审查，主会话用 debug 版 `ferro` 逐条复现。复现 fixture
是一次性的，未进仓；每条的「复现」写到可以重做的程度。陷阱的一般化教训见
`issues.md`「2026-10-01 审查：读写两侧的静默错误」。**每条修复单独提交并带回归测试**。

**严重**（静默产出错数据，或 panic 打断整个批次）

| # | 位置 | 问题 | 复现 | 修法 |
|---|---|---|---|---|
| ~~S1~~ **已修**（`f1c7908` reader 丢残末帧 + 告警、中间帧不完整报错；`58f9c90` `Trajectory::check_same_atoms` 挂到六个时间相关分析）| `readers/lammps_dump.rs:58,108,124-127` | 截断的 dump（MD 中断，常见）不报错：残帧照收；`lines[i]` 无越界检查 | `head -n 9000 tests/43Z43P15A_NPT_5.lammpstrj`：第 4 帧 939 原子无告警；`traj msd` panic `msd.rs:133`（退出码 101，批内不再是「跳过 + 退出码 1」）；`net` 退出码 0、`[inputs]` 写 `atoms 2004 ok`，mean_cn 4.00 → 3.97；截断在 BOX BOUNDS 前则读取本身 panic `:58` | 原子行不足 n 时报错（或丢末帧并告警）；`lines[i]` 改 `get` + `bail`。另在 core 加「各帧原子数与元素序列一致」守卫，msd/vacf/bondlife/rotcorr/angle/net/map 共用 —— 现在只有 gr 有 |
| ~~S2~~ **已修**（`c4f7508`，期望值对 ASE `Prism.vector_to_lammps`）| `writers/lammps_dump.rs:82-115` | 坐标转进 LAMMPS 规范胞（下三角），速度与力**未同步旋转** | extxyz `Lattice="0 4 0 -4 0 0 0 0 4"`、原子 (1,2,3)、v=(0.5,0,0)、f=(1,0,0) → dump 坐标 (2,−1,3) 已转，vx=0.5、fx=23.06 仍在 x（应在 −y）。影响 `convert -o *.lammpstrj` 与 `collect --type inspect` | 求 R = L_lmpᵀ·(Mᵀ)⁻¹，v、f 同样左乘 R；测试用非下三角胞 + 力 |
| ~~S3~~ **已修**（`5882fd4`；复核时另发现 reader 不认速度块的模式行、ASE 写的速度被静默丢弃，`7f333aa` 一并修）| `writers/vasp.rs:38-51` | POSCAR 坐标按元素分组写，速度按原序写 | O(v=.1) Si(v=.2) O(v=.3) → 坐标序 O O Si，速度序 .1 .2 .3，后两个原子速度互换 | 速度用与坐标相同的分组循环 |
| ~~S4~~ **已修**（`29ff9d6`，口径照 pw.x `cell_base.f90`；alat 坐标、不写选项一并支持；L6 的 qe.rs 部分同修）| `readers/qe.rs:174-179`（及 `:69-72`、`:105-108`） | 卡片单位只认 `{...}` 写法，`CELL_PARAMETERS bohr`、`ATOMIC_POSITIONS crystal`、`(crystal)` 一律回落 Å；`alat` 的拒绝也只对 `{alat}` 生效；QE 无选项时的默认并不是 Å | Si 原胞 `CELL_PARAMETERS bohr` + `ATOMIC_POSITIONS crystal` → 胞读成 5.13 Å，第二个 Si 在 (0.25,0.25,0.25) Å，无告警 | 首行去掉括号后不区分大小写匹配关键字；无选项按 QE 默认或报错 |
| ~~S5~~ **已修**（`a2c013a`，照 CP2K `parser_get_logical`；另修 &COORD UNIT 非 bohr 当 Å、&CELL 认 `[bohr]` 而非不存在的 UNIT 关键字）| `readers/cp2k.rs:84` | `&COORD` 的 `SCALED` 只认 `.TRUE.`；单写 `SCALED`（CP2K lone keyword 即真）、`SCALED T` 判为假 | 同一输入：`SCALED` / `SCALED T` → 坐标 0.25（当 Å）；`SCALED .TRUE.` → 1.3575 | 按 CP2K 逻辑值规则：缺省值、`T` `TRUE` `.TRUE.` `YES` `ON` 均真 |
| ~~S6~~ **已修**（`3d3f6ed`，`ferro_io::is_extxyz` 两处分派共用；写 `.xyz` 仍是纯 XYZ、丢晶胞，未动）| `cli/io_dispatch.rs:18`、`ferro-python/src/io.rs:52` | `.xyz` 一律走纯 XYZ reader，extxyz 的 Lattice/能量/力/应力/速度全丢。**ferro 自己** `dataset --type nep|extxyz` 写的就是 `.xyz`（`cmd/dataset.rs:250`） | 同一 extxyz 文件改名 `.xyz` → `info` 报 `Cell: none (non-periodic)` | 读 `.xyz` 时看第 2 行，含 `Lattice=` 或 `Properties=` 即转 extxyz reader。Python 侧同步（两处分派会漂，见 `progress.md`） |
| ~~S7~~ **已修**（`2abe03a`，按用户裁定：同名依次加 `_2`、`_3` 并打 Note，不拼父目录、不改产物落点。用户原以为 map/net 的产物跟着输入目录走，实际是当前目录或平铺进 `-o`，与手册一致）| `cli/batch.rs:29 label_of`、`cmd/map.rs:183-190 stem_for`、`cmd/net.rs:259` | 多输入时只取 stem，不同目录的同名文件撞名。**与已登记的「`map` 产物名随输入数变」不是一条**：那条是 N=1 与 N>1 两套命名，这条是 N>1 下静默覆盖 | `-i 'runs/*/prod.lammpstrj'`（两个不同体系）：`map density` 两次都写 `density_prod.cube`，后者覆盖前者，退出码 0；`traj gr` 的 csv 10004 行 `file` 列全是 `prod` | 复用 `collect` 的规则（`issues.md:448`）：撞名改 `<父目录>_<stem>`，仍撞报错，在读第一个文件前判。与合规扫描第 2 条（`map` 命名）一起定 |

**中**

| # | 位置 | 问题 | 复现 | 修法 |
|---|---|---|---|---|
| ~~M1~~ **已修**（`a4c2b6d`，按用户裁定：**不看注释、不按列数猜**，style 由 `--atom-style` / Python `atom_style` 显式给出，缺了在读第一个文件前报错；列数只收基本列或 +3，原子数核对 `N atoms`；image flag 展开同 LAMMPS / ASE。修法栏的「按列数推断」作废）| `readers/lammps_data.rs:112,158` | `Atoms` 段无 style 注释时按 full 解析（注释在 LAMMPS 里是可选的） | 5 列 atomic → `Atoms: 0` 无报错；8 列（带 image flag）→ 第 2 原子 (1,1,1) 读成 (1,0,0) | 无注释按列数推断（5/8 atomic、6/9 charge、7/10 full），歧义报错；列不足报错而非跳过 |
| ~~M2~~ **已修**（`6dff5b1`，用户裁定缩放坐标按 LAMMPS 定义加 lo，与 ASE 不同；顺带读 `mass`、坏字段报错。`typelabel`、`ix iy iz` 未接，见下）| `readers/lammps_dump.rs:182` | 无可识别坐标列时坐标静默全 0 | `ITEM: ATOMS id type element xsu ysu zsu` → 两原子都在 (0,0,0) | 支持 `xsu/ysu/zsu`；一列都找不到报错 |
| ~~M3~~ **已修**（`aa82fb2`，12 个 `XxxParams::validate` 共用 `ferro-analysis/src/check.rs`，CLI 在建目录、读文件前调用；angle / vanhove / sq 改 `Result`。cube_*、network、chg_sdf 仍返回 `Option`，见下）| `md/gr.rs:222`、`md/sq.rs:118`、`md/cube_density.rs:84`、`cmd/traj.rs:521` | 参数未在读文件前校验，违反 `CLAUDE.md` | `--dr 0` / `--dq 0` / `--nx 0` 均 panic（101）；`--r-min 5 --r-max 3` 读完文件才逐个 skipped；`angle --d-angle 0` 报误导的「empty trajectory?」（`calc_angle` 返回 `Option`，三种失败揉成一种） | CLI 构造参数时统一查 `> 0`、`r_min < r_max`、`n ≥ 1`；`calc_angle` 改 `Result` |
| ~~M4~~ **已修**（`c35567d`，每轴下标列表，窗口 ≥ 网格时取全部；对照 C 全遍历与 numpy 独立实现。默认参数 50³/0.7 Å 不受影响，测试轨迹 `--radius 15` 修复前多计约 0.9%）| `md/cube_radius.rs:78-80,118-125` | 搜索窗 `2s+1` 超过网格点数时 `rem_euclid` 把多个偏移折到同一体素，重复计数。**不只粗网格**：40³ 网格、`--radius 4.9`（10 Å 胞）也中 | 单原子单帧，10 Å 立方胞，`--radius 4.9`：4³ 网格 max=8；40³ 网格 max=2、sum=31852（解析值 ≈ 31540）。单帧单原子任一体素应 ≤ 1 | 对折叠后的 `(ix,iy,iz)` 去重（同 `box_builder` CellList 那次），或把窗口截到 n |
| ~~M5~~ **已修**（`694f4e7` 按轴 pbc 基础设施；`48bf11e` cube 参考结构折回，**只在 cube 里做的特例**。读入时平移原点 / 折回的方案实测会扰动 NPT MSD，已否决，见 issues 与手册 msd 页）| `md/cube_density.rs:202`（origin 恒 0）、`md/util.rs:40`（平均帧不折回） | LAMMPS 盒子 `lo≠0` 时原点被丢，密度按分数坐标折进 [0,L)，原子按原坐标写出。**与已登记的 NPT 体积口径不是一条** | `tests/70Z30P00A_NVT_5`（lo=1.82）：cube 原点 0、边长 39.22 Å，5003 原子里 659 个在网格外。周期意义上密度与原子一致，**可视化效果未验** | 平均帧折回 [0,1) 后写，或把盒子原点带进 `CubeData.origin`（reader 也要保留 lo） |

## 库 crate 的错误类型：`ferro-io` 用 `anyhow`（2026-09-28）

违反 `CLAUDE.md`「库 crate 用 `ChemError`」。2026-09-28 摸底后的处置：

- **`ferro-io` 继续用 `anyhow`（用户已同意）**。30 个文件约 173 处（`context` /
  `with_context` 103、`bail!` / `ensure!` 67），另靠 `?` 自动吞 `ParseFloatError`、
  quick-xml、npy 的错误；15 处测试断言错误文本。理由：
  - 调用方只有 CLI（`{e:#}` 打整条链）与 ferro-python（只取 Display），全仓零
    downcast，**没有人按变体分支** —— 类型化错误的价值在这里是零
  - `ChemError` 的变体装的全是 `String`，迁过去仍是一段话，只多一个
    `Parse error: ` 前缀；反而丢掉 `.context()` 的「文件 → 帧 → 原因」链，那是读
    AIMD 输出时最需要的诊断
  - 该改的时机：出现真要按类别分支的调用方（ferro-python 映射 `FileNotFoundError`、
    `collect` 区分「格式认不出」与「文件损坏」）。届时设计**真正带类型**的错误，
    不是套现在的 `ChemError`

## 自旋多重度：CP2K 源码核对与行业做法（2026-09-30 / 10-08，留作以后推断未成对电子的参考）

**CP2K 自己不推断多重度**（`examples/external/cp2k-2026.2`）。`qs_environment.F:2013` 在
`MULTIPLICITY` 未给（=0）时只按价电子数奇偶取 1 / 2；`:2032` 奇数电子不开 UKS
直接 abort（有 smear 时例外）；`:2090` 按 $N_\alpha = (N+M-1)/2$、
$N_\beta = (N-M+1)/2$ 切分，$N+M-1$ 为奇数报「try a different multiplicity」。
旁路只有 `RELAX_MULTIPLICITY`（按 Aufbau 原理占据，需 `ADDED_MOS`）、`SMEAR` 下的
`FIXED_MAGNETIC_MOMENT`、`&KIND/&BS`（只改局域初猜，不改总 M）。

**已核对一致的**（不必再查）：

- 奇偶口径：ferro 用 $\sum Z - q$，CP2K 用 $\sum q_\text{val} - q$。`cp2k_basis_db.rs`
  全部赝势的芯电子数 $Z - q_\text{val}$ 无一为奇数，两边奇偶恒同
- `guess_spin` 经 `reconcile_parity` 出的 M 恒满足 $N+M-1$ 为偶数；`mult > 1` 即写
  `UKS`，故自动路径碰不到 CP2K 的两条 abort

**行业做法**（2026-10-08 调研，修 C-D1 前用户要求查）：没有软件单凭结构可靠给出未成对数 ——
它取决于配体场（高 / 低自旋）与磁性中心间排布，结构看不出。分层：① 奇偶下限（CP2K 缺省、
pymatgen `Molecule` 缺省）；② 用户 / 磁矩给定（ASE 写 Gaussian 取 Σ 初始磁矩 + 1；Gaussian、
ORCA、xtb `--uhf` 都要用户给）；③ 氧化态 + 电子计数（pymatgen `oxi_state_guesses` 按 ICSD
概率排序、`BVAnalyzer` 按键价和；`get_crystal_field_spin` 只做 d 区，高 / 低自旋与八面体 /
四面体由**用户指定**）；④ 配位环境判高低自旋（cell2mol 2024，仅 3d 单核，d⁴–d⁸ 才有歧义，
97–98%）；⑤ 周期体系**不定总自旋**：MP 给偏大初始 MAGMOM（Ce 5、Eu 10）、VASP 缺省
1 μB/原子且不设 `NUPDOWN`，由 SCF 弛豫；磁序靠枚举 FM / AFM / 亚铁磁再比能量（Horton 2019
npj Comput Mater，atomate）。MP 的 POTCAR 对 Pr–Lu 多数用 `Ln_3`（4f 在核内）。
ferro 做的是 ①–③，C-D1 后仍缺：多解取「总和最大」而非按统计概率（`assign_oxidation_states`
文档说「调用方给歧义警告」，实际没有调用方给）；共价分子只给奇偶（$O_2$ 报单重态，OpenBabel
从 xyz 读同样报 1）。第 1、4 条的裁决应参照 ⑤。

**实测线索**：用户曾用网页版 Claude（约 Opus 4.6/4.7）生成过 MnS 夹杂 MnO 的 CP2K
输入，多重度 500 多，看着异常但 CP2K 跑得很顺、SCF 收敛特别顺利，细节已不可考。
与第 4 条吻合：Mn²⁺ 是 d⁵，逐个取 5 个未成对电子，约 100 个 Mn 就是 M ≈ 501。铁磁态
全部自旋平行、不存在阻挫，SCF 好收敛是可以预期的；但 MnO、MnS 实际都是**反铁磁**，
收敛顺利不等于找到了基态。定第 4 条时可拿这个体系对照铁磁与反铁磁（`&BS` 初猜）的能量。

## 已完成项的判据（0.3.0 以后；0.3.0 以前的归档已删，见 git 历史）

### `dataset collect`：输入是当前目录下的裸文件名时，产物位置错（2026-10-17 发现并修复，`7735386`）

`Path::new("a.out").parent()` 是 `Some("")` 而非 `None`，`canonicalize("")` 失败，分组键成了空路径。
后果：带 `-o` 且另有子目录输入时，子目录组的 `rel` 成了绝对路径，产物甩开 `-o` 写到输入旁；不带 `-o`
时写出 `./system.db`。修法：空父路径当 `.`（陷阱见 `issues.md`）。

**判据（用户定，方案 A）**：多组时若某组目录恰为公共祖先，祖先上退一级，该组以目录名命名、全部
嵌套在 `-o` 下（`md/x.out md/sub/y.out -o out` → `out/md.db`、`out/md/sub.db`）。否决的方案 B 是报错
让用户分两次跑。只有一组时仍写进 `-o` 本身，不变。

### 2026-10-11 全库审查：15 条已修（2026-10-14，随 `v0.3.6`）

范围由用户 2026-10-14 圈定（另 8 条挂起，见优先级高的同名一节）。每条一个提交，均先写复现测试再修。

| # | 判据 | 提交 |
|---|---|---|
| B-1 | 方案二选一：补解析项 $-\rho\int_0^a r\sin(qr)\,dr$ 要假设 $[0,r_{min}]$ 无原子对，r_min 越过第一峰时静默错；**定为 sq 恒从 0 积分、删 `--r-min`**（用户选，与 TRAVIS / ISAACS 等从 0 起的直方图及教科书下限一致）。库层同守，拒收 r_min≠0 的 `GrResult` | `93b0ad4` |
| A-12 | 原点放 `Cell` 不放 `Frame`（用户选）：放 Frame 则每个产生新 Frame 的地方都得记得带，漏一处就静默丢。只作元数据（同 ASE `celldisp`），分析一律不用 —— 读入时平移已在 M5 实测扰动 NPT MSD 被否决。范围只限 dump / data（用户选），cube 的盒子原点未接 | `3c775ca` |
| B-2 | 值属于体素中心，cube 规范第 i 点在 origin + i·spacing，故 origin 取半体素；四处同一变化原因但一行，按 R6 不合 | `0dc1387` |
| A-10 | 段头名单取 LAMMPS `read_data` 文档全集；Masses / Atoms 按段名定位（段序任意）；Masses 行恰两项否则报错 | `dc1cf1a` |
| E-10 | 放置按胞面算，原子须先相对 `Cell.origin` 折回；非周期方向越界无合理折回，报错（自定） | `57ead5d` |
| E-12 | 撞名键小写、显示原样；`net --export-traj` 经同一函数一并覆盖 | `7685d40` |
| A-2 | 写侧列名定为 ASE 的 `momenta`（用户选）并带 `masses`：ASE 只认 momenta，`velocities`/`vel` 在 ASE 里是普通数组。换算常数取 ASE 默认 CODATA 2014 的 `units.fs`，不用本仓 2018 常数自算。顺带了结中级「extxyz 把 momenta 当速度读」 | `66e9340` |
| E-9 | 奇异胞报错，不做 ASE 式「按原子跨度补 c」（那是新功能，自定） | `03c7fa7` |
| D-4 | 区分「用户给了」与「默认值」只能靠 `Option`（derive 拿不到 value_source）；默认值缺标签时关掉该条并 Note，显式给了才报错 | `2cb7e29` |
| A-17 | 缺 cell 行单独计数（`n_no_cell`），全丢时报错点名 `&CELL`；第三个同构 xyz 块整份报错 —— VELOCITIES 若排在力块前会被当力读 | `f350185` |
| C-9 / D-2 / E-1 / E-3 | 只改手册 / 帮助页。E-3 顺带改了 6 个分析手册页同样的 8 处 `-o` 旧写法；sdf / chg-sdf 的 `-o` 本是目录不动 | `055596b` `7bfcd87` `a3154fe` `0807f17` |
| E-15 | 下限取依赖树（两个 workspace 各自 `cargo metadata`）最高 `rust-version` = 1.89；clippy `incompatible_msrv` 零告警，**未装 1.89 工具链实编** | `1b67d80` |
| A-7 | 判重同 ASE `equivalent_sites`：逐分量 1e-3、只在同一位点的像之间；跨位点重合报错（部分占位本已拒收，不多拒任何此前能读的结构）。`crosscheck_spacegroups.py` 复跑与 10-03 基线相同 | `dd76f07` |

### vasprun 的非对称应力：导出 nep 后 ferro 读不回自己的文件（2026-10-13 发现并修复）

VASP 6.5.1 写进 vasprun.xml 的 `<varray name="stress">` 不对称（`tests/vasp651_sp/frame_0011`：
xy = −0.95709、yx = −0.95684 kB；200 帧中相对非对称最大 4.4e-5，virial 绝对差最大 2.0e-3 eV），OUTCAR
的 `in kB` 只打 6 个数故看不见。原样保留时 `dataset filter --type nep` 写出的 `virial=` 被本仓
`readers/extxyz.rs::parse_tensor9` 的对称检查拒收。

**用户选的修法**：vasprun 读侧取 $(\sigma+\sigma^T)/2$。与 dpdata（原样保留 9 个数）在 vasprun 路上因此差
最多 ~1e-3 eV，OUTCAR 路不受影响。否掉的两条：extxyz 读侧放宽容差（规范要求对称，放宽等于替上游的
不对称背书）、写侧对称化（不对称的张量仍留在 `Frame` 里，别的出口照样会漏）。

**OUTCAR 与对称化后的 vasprun 仍差半个非对称量**（3 帧 fixture 实测 2.0e-3 → 1.0e-3 eV）：OUTCAR 的 `in kB`
不取平均，而是各挑一个分量（frame_0011 的 XY 取 yx、ZX 取 xz）。VASP 输出本身的属性，不是读取错误。

**遗留**：`parse_tensor9` 的容差是 `TENSOR_TOL · max(|v|, 1)`，stress 量级 < 1 eV/Å³ 时退化成绝对 1e-6 ——
同一张量写成 `stress=` 能过、写成 `virial=` 被拒。本次修复后无实例触发，口径未动。

### 2026-10-03 四路审查：严重 / 中两级已修条目（2026-10-03 – 10-09，0.3.5）

每条都先写在旧代码上失败的测试。未修的（搁置、QC 专项、轻级、布局建议）留在「优先级高」原节。

| # | 问题 | 提交与判据 |
|---|---|---|
| A1 | 无 symop 循环一律当 P1，H-M / IT 号写着非 P1 也不报；CIF2 点号标签 `_space_group_symop.operation_xyz` 同路 | `acb282d` 数据表、`7638dc0` 查表、`1a76d49` 接入；判据见 `issues.md`「CIF 空间群符号展开」 |
| A2 | CRYST1 只按第 0 帧写一次，reader 套到所有帧 | `c8ee8c8`；顺带修了同函数里的取向错位：CRYST1 只存六参数，非标准取向的胞原样写坐标，读回错位 —— 现经分数坐标转到 a 沿 x、b 在 xy 平面，同 ASE `standard_form`；对拍 ASE 写出逐位相同。**遗留**：ASE 把 `ENDMDL` 后的 `END` 读成第 N+1 个空帧，ASE 自己不写 `END`，未改 |
| A3 | 无胞帧取包围盒尺寸写成 `0..L`、坐标不平移、边界硬写 `pp` | `a34b488`；无胞帧盒子取包围盒每侧外扩 1 Å，坐标原样、lo≠0；未在真 LAMMPS 里读过 data —— 本机无 LAMMPS |
| A4 | 未闭合引号 `&inner[end+1..]` 越界 panic | `9ebb707` 逐字符解析，规则同 ASE `key_val_str_to_dict`，一并支持 `[]` `{}` 括号值与 `\` 转义；有意与 ASE 不同两处：未闭合报错、`a="" b=1` 读作两键，见 `parse_comment` 文档 |
| A5 | 不带 `=` 的裸键（ASE 读作 True）与下一个键名粘连，Lattice/Properties 被静默丢弃 | `9ebb707`；同 A4 |
| A6 | mixed type system（`type.raw` 全 0 + `real_atom_types.npy`）读入即全成 `type_map[0]`；本文件「DeePMD mixed type」一节只写了「未实现」 | `e46a5ab` 拒收；完整支持仍是下面「DeePMD mixed type」一节 |
| B-1 | 奇异胞 `expect("cell is non-singular")` panic。ASE 二维材料约定（c=0、pbc TTF）即触发；gr/msd 等同一文件是「singular 跳过 + 退出码 1」 | `bb25ba7`；`check::invertible_cells` 在 vanhove / cluster SDF / cube_jump 入口逐帧查，`calc_cluster_sdf` 改返回 `Result<Option>`。实测 c=0 slab 过全部 13 个 CLI 分析，只有这两个 panic |
| B-2 | 只有 `positive` 查 `is_finite`，两者放行 `+inf` | `a5d1be8`；报错措辞改为「must be a finite number …」 |
| B-3 | q=0 捷径返回 1.0，而公式极限是 $1+4\pi\rho\int r^2(g-1)\,dr$ | `98b5e97`；**用户否了「取 q→0 极限」**：有限盒子的 S(0) 无物理意义，改为 q-min 必须 > 0；并加 `[inputs]` 的 `q_trunc = 2π/r_max` 与低于它时的告警。查证的行业做法与出处见手册 `sq.md`「The low-q limit」 |
| C-D1 | 镧系 `group_number` 为 None 落进主族分支，算出垃圾未成对数且无告警；CP2K/QE 默认 auto-spin | `4a34108`；用户要求先查行业做法，结论记在下面「自旋多重度」一节。范围扩到同一函数的 4d/5d：一律低自旋、d⁸ 取平面四方 0；氧化态凑不平单独告警。镧系假定 4f 在价层，加告警 |
| D-S1 | `PRESSURE 1.01325E+05 # bar`，CP2K 该关键字单位就是 bar → 约 10 GPa | `a92e4b9`；用户定：`--pressure P`（bar）取代 `--barostat`，给出即 NPT_F，不给为 NVT；只拦非有限值，负压允许 |
| D-S2 | `&HF` 缺 `&INTERACTION_POTENTIAL POTENTIAL_TYPE SHORTRANGE / OMEGA 0.11`，按 COULOMB 全程 25%，得到不存在的泛函且能跑 | `bd21bc2`；用户注：杂化泛函一般用不到，作功能完善 |
| D-S3 | `--thermostat none`（帮助写 NVE）仍写 `ENSEMBLE NVT`，CP2K 默认恒温器 NOSE | `bd21bc2`；系综按恒温器 × 压强查表：csvr/nose → NVT/NPT_F，none → NVE/NPE_F，langevin → `ENSEMBLE LANGEVIN` + `&MD/&LANGEVIN`，与 `--pressure` 同给报错。**顺带修**：原 langevin 写进 `THERMOSTAT%TYPE`，CP2K 只认 NOSE/CSVR/GLE/AD_LANGEVIN（`input_cp2k_thermostats.F:237`） |
| D-S4 | `--export-traj lammpstrj` 恒以 real 单位写出，无视 `--units`；`inspect.rs` 恒写 metal，两出口口径相反 | `e78e609`；沿用 `common.read.units`。**未改**：`inspect.rs` 是 CP2K 产物导出、无 `--units`，恒写 metal，属有意的固定口径 |
| D-S5 | 写 `TIMECON_CSVR`，CP2K `&CSVR` 只有 `TIMECON`（`input_cp2k_thermostats.F:580`，无别名）；默认 `--task md` 即走此路，推断 CP2K 解析即中止（本机无 CP2K，未实跑） | `bd21bc2`；第 8 轮对照 CP2K 源码时新发现 |
| A7 | VASP 6 元素行带 POTCAR 哈希 `Na_pv/6a2f546d` 原样当元素名（ASE 读 Na） | `d7876a1`；VASP 6.4.2 HDF5 编译才写 `标签/哈希`，6.4.3 起官方称已解决；5.4.4 纯符号、VASP 4 无元素行均有测试钉住。顺带：`chgcar.rs` 坐标类型只认 d 开头，改为同 POSCAR 的 C/K 规则 |
| A8 | 任一 data 块无胞（常见 `data_global`）整份报错；ASE / pymatgen 跳过 | `714cac6` |
| A9 | 无胞帧写成只有 Cartn 的 CIF，ferro 自己读不回（ASE 读成 pbc=False） | `714cac6`；2026-10-08 用户定：无胞结构只写 xyz —— 改为 **writer 报错**，同 POSCAR writer，读侧不改 |
| A10 | 77–78 列为空时元素为 `""`；坐标坏或行短于 54 时原子静默丢弃；多字节字符切片 panic（推断） | `a6740d5`；用户定按 PDB 列对齐从原子名推元素，不照 ASE 先试两字母。顺带：77–78 列大写 `FE` 原样当元素名，改为规范成 `Fe` |
| A11 | 段名前缀匹配：`&COLVAR` 里的 `&COORDINATION` 先命中 `&coord`，合法输入报错（regtest `QS/regtest-gpw-2-2/H2O-meta.inp` 等 ≥9 份） | `276b465`；用户定整行比较，`&KIND` 带 kind 名故比首词；剥注释后补 trim_end。新发现：`&CELL_REF` 的 ABC 覆盖 `&CELL`，`parse_cell_section` 改为只读本层。regtest `H2O-meta.inp` 实读通过 |
| A12 | namelist 用空白分隔（`ibrav=0 nat=3 ntyp=1`）时整串成了一个值：`nat` 校验静默跳过、非零 ibrav 不被拒、`starting_magnetization` 丢失 | `307cd4d` |
| A13 | 空文件 `lines[i]` 越界 panic | `d734593` |
| A14 | 三斜判定 `!= 0.0`，cos 90° 残留 1e-16 使正交盒写成三斜（`convert tests/70Z30P00A_NVT_5.lammpstrj -o x.lammpstrj` 即中） | `a34b488`；`cell_to_lammps` 把相对最大边长 1e-10 以下的倾斜置 0，两份同改 |
| A15 | 边界恒写 `pp pp pp`，不看 `frame.pbc`；reader 却尊重标志（TTF slab 往返成 TTT） | `a34b488`；只改 dump，data 文件无边界标志 |
| B-4 | `[inputs]` 的 `frames` 填 `r.r.len()`（bin 数）；angle 的 `atoms` 填 `r.elements.len()`（实 2004 写 4） | `08eaa0a`；NPT fixture 实跑 5 帧 / 2004 原子，与 `ferro info` 一致 |
| B-5 | vanhove `--tau` 帮助写「default: half trajectory」，实为 `n_frames−1`（与手册一致）。`help_sync` 不比默认值 | `5d20746`；用户选默认改 N/2：同 msd/vacf/rotcorr/bondlife 的 `--max-lag` 与 gmx `-acflen`；pymatgen `VanHoveAnalysis` 是另一套设计 —— 扫一串 lag、每个 lag 固定 50 个原点。破坏性，进 changelog |
| B-7 | 头部写截断值却标「requested」，批内第一个文件的截断值冒充全批请求值；测试 `test_meta_lines_report_clamped_rmax_and_composition` 钉住的是错误语义 | `da11aed`；`GrResult.r_max_used`，`params.r_max` 留请求值；未截断的 gr / sq 产物逐字节不变，三斜 fixture 只有头部一行由 4.8136 变回请求值 10.005 |
| cube_jump 校验 | 无 `validate()`（M3 那批漏了）：`nx=0` 在 `:95` 除零；`threshold ≤ 0` 全算跳跃；`:180`「与 msd.rs 一致」已过时（msd 已改 TOR） | `13044de`；`CubeJumpParams::validate`，`calc_cube_jump` 入口也调，接进 `check.rs` 测试宏 |
| C-D4 | 力 / 应力含 NaN 的帧通过筛选进训练集（`fold(0.0, f64::max)` 吞 NaN，`m > f_max` 对 NaN 为假） | `c8970dd`；判据开着时非有限峰值判坏，峰值保留 NaN 供报告；连带修了下面「轻」C 里的 `-f nan` / `--oo-min nan` / `--al6 nan` 与诊断下溢 |
| C-D5 | `--al6` 自动截断取峰后 3 Å 窗口内第一个严格最小 bin，g 恒 0 的平台上取到第一个零 bin。`collect tests/cp2k_md_3frames.out` → `filter --al6`：cutoff 2.05，cn5=2 / cn6=46；numpy 截断 2.2–2.6 得 cn6=48 | `deac6ae` 取最长零平台中点、移入 `md/shell.rs`；`e9a88a3` 改在 0.10 Å 滑动平均上找极小，实测 2.05 → 2.57 Å。**遗留**：手册 `dataset/filter.md` 的「2.45 Å」与三对极小表是旧算法在参考体系上测的，数据不在仓内，未重测 |
| C-D6 | 多 system 报告 `concat_union` 只留第一份 meta，逐 system 统计（O–O 分位、自动 rcut）冒充全批 | `162be7d`；O–O 分位与实际 Al–O 截断改作 `[inputs]` 列，表头不再写；`diagnostics::quantiles` 拆出，只读模式顺带补打分位行，原先只在 csv 表头里 |
| D-M1 | 仍默认 `dt=1.0`、不做 `check_frame_spacing`，与 CLI 已改必填的口径漂开 | `ca24322`；`check_frame_spacing` 下沉为 `ferro_analysis::md::check_frame_spacing`，CLI 五处与 Python `msd` 共用；`dt` 必填；wheel 实测：NVT fixture 步距 12000/13000 不均被拦，CLI 同句 |
| D-M2 | ① 只给 B 居中，A 窄时不居中（与文档不符）；② B 沿自身单位矢量平移、胞用 A 的，倾角不同时剪切错位；③ gap 沿矢量量，与 `add_vacuum` 的垂直间隙口径不一 | `58b47f6`；用户定保留不拉伸的填充式拼接，不做 ASE `stack` / pymatgen `from_slabs` 的共格拉伸；两块整体平移、胞中心面内居中，gap 按面间距，界面不平行（>1°）报错。原 13 个测试中正交胞的结果不变，三斜测试改按法向高度断言 |
| D-M17 | 奇异胞（c=0 slab）在 `expect("cell must be non-singular")` panic。CLI `net` 在前面已拦，Python 绑定直调仍会 panic；`dft/chg_sdf` 经 `process_frame` 同路（cube 文件的胞，推断难触发） | `99d574c`；核实：Python 并未暴露这几个函数，原描述「Python 直调会 panic」不成立，暴露面只在库公开 API。`classify_frame(_detailed)`、`build_network_graph` 入口查奇异胞并返回 `Result`；`find_clusters` 返回 `None`，`classify_trajectory`、ml geometry 两个配位函数、`cutoff_scan`、`pooled_coordination` 向上透传；network、filter、cube_sdf 入口已查，内部 `expect` 写明由谁查过；`chg_sdf` 入口遇奇异胞返回 `None` |
| D-M3 | `--P-O=nan` 通过校验：P 的 cn=1314、全 Q0，退出码 0 | `b71b71c`；nan / inf / ≤0 一律报错 |
| D-M5 | `&THERMOSTAT TYPE LANGEVIN` 非法（合法值 NOSE/CSVR/GLE/AD_LANGEVIN）；CP2K 是 `ENSEMBLE LANGEVIN` + `MD/&LANGEVIN` | `bd21bc2`，随 D-S3 一并改 |
| D-M13 | `-s` 不校验字符：`-s 'a/../../escaped'` 写到 `-o` 之外 | `7f36b12`；`batch::check_suffix`，在 `CommonArgs::inputs`、chg-sdf、bader 三处入口调用。**未覆盖**：`dataset merge --suffix` 是目录后缀，合法值含 `.`，另一套语义，未动 |
| D-M15 | `--start/--end` 校验前就建 `-o` 目录；写侧格式读完整条输入后才查 | `d66793f`；参数检查移到建目录前；写侧格式判定抽成 `io_dispatch::out_format`，`write_trajectory` 改走它，convert 读输入前先调 |
| D-M16 | 对所有子命令剥离 `--Xx-Yy=v`，`traj gr --P-O=2.3` 静默吞掉 | `b71b71c`；判据放进 `split_pair_args` 本身，有测试 |

### scripts/：net 剩余表的画法 + labels.csv + 数据导出（2026-09-28 落地）

原待办「net 剩余四张表的画法」经一轮 grilling 定案，结论与判据：

| 表 / 项 | 结论 |
|---|---|
| `network_composition` | **不单独画**。它 = qn ∪ coordination ∪ 配体四类，唯一独有的配体四类进了 `bridge` 左格 |
| `network_ligand_type` | 被 `bridge` 右格取代。数据源改用 `linkage`：`ligand_type` 的 `O_t` 只有一行、不分伙伴，而用户选了「三簇配体按 C(k,2) 对连接照计」 |
| `network_linkage` | `linkmap` 热图（用户指定形式：Al4/5/6 × Q0…Q3），2 个形成子 3 张、3 个 6 张 |
| 实验数据对比 | **结案，不做**。用户用导出的 `<stem>_data.csv` 自行与实验数据合图 |
| `--partner` | **删除**。两级堆积柱在实物上难以读懂；`qn_partner` 表不动 |

用户定下的口径（逐条问过，别再重开）：X-O-X 分母只含形成子间连接；三簇配体照计；
Qn 轴取 n（P–O–P 个数）；每张热图自归一、格里写一位小数、色标按列共用；同种元素只画
上三角；多成分排成一张大图（行 = 成分）；轴取全部成分并集；全部在 Python 做，不改 Rust；
比例由 count 求和后相除，不画误差棒（合并行后逐帧均值无从恢复，sd 不能相加）。

`labels.csv` 的设计由用户提出：首跑按实际输入生成模板、之后读它覆盖默认。于是
「文件里有、数据里没有 → 报错」不再成立（文件跨批次、跨脚本共用），改为静默忽略；
打错的 key 以「追加新 key」的 `Note:` 暴露。覆盖全部 `plot_*.py`，用户不要
`--labels` 手动指定（只出英文版）。

验证：默认 labels 下 gr/sq/angle/msd/qn/cn 六张 png 与改动前逐字节相同；
`bridge` 的 X-O-X 总数 3589 = 3583 个 O_b + 2 个 O_t × 3；热图抽查 Al4–Q0 =
863/2697、Al4–Al4 = 10/12。fixture 无真三形成子体系，6 列排版是把 Zn 当形成子演示的。

**被否的前一版热图（2026-08-15，五种版式）留下的、仍然成立的事实**——将来若要做
「跨成分的择优连接」，从这里起步，不必重跑：

- 行归一化只除掉一端的丰度，跨成分比较要靠**零模型**：桥端随机配对
  `p_i = e_i/2N`、`E_ij = 2N·p_i·p_j`（i≠j）/ `N·p_i²`（i=j），取 `log2(O/E)`，
  `p_i` 逐成分各算、边缘取自 linkage 自身
- 实测：`43Z43P15A` 的 Al–O–Al 压制 39 倍；`70Z30P00A` 的 P–P 三格 log 比贴近 0，
  即二元磷酸盐的链是随机连接的
- 对角线贡献两个同类桥端：边缘统计算两次、矩阵取值只算一次（原型第一版在此翻倍）
- 现行 `linkmap` 是**每图自归一的占比**，不是 O/E；两者回答的问题不同

### `ferro traj bondlife`：键寿命与成键断键事件（2026-09-27 落地）

- 用户要研究水与玻璃之间化学键的形成与断裂。查证后同时给两种相关函数：
  间歇 C_I（Luzar–Chandler / gmx hbond -ac，FFT）与连续 S_C（Rapaport / MDAnalysis，
  按成段长度精确计算），外加逐帧 n_bonds / formed / broken 事件表
- 阈值：`--r-bond` 必填；`--r-break` 形成滞后区（Yamamoto–Onuki 的 A₁/A₂，读过原文
  式 3.9–3.14：3D 取 A₁ = A₂ = 1.5σ，结果对阈值不敏感只要落在两峰之间）
- `--intermittency k`：只补两段成键之间 ≤ k 帧的缺口，通向轨迹两端的不补；
  只作用于 S_C 与事件，C_I 本就容忍断开
- 寿命两种：梯形积分末值（下界）与 1/e 首次下穿（降不到给空，不外推）
- 验证：单测逐原点暴力；examples Zn–O 两组参数 vs 独立 numpy，事件逐帧相等。
  该玻璃 Zn–O 在轨迹内 C_I 仍 ≈ 0.97，1/e 寿命为空 —— 常温玻璃的正常情形
- 测试 fixture 一度把多个 O 放在同一个 Si 周围不同高度，只有一个真在阈值内，
  暴力参考却按设计键长算 —— 改为每组 Si–O 独立、组间相距 6 Å

### rotcorr：`--vector bond` 与 `--legendre 1|2`（2026-09-27 落地）

- 起因：用户研究玻璃；sum 模式对 PO₄/SiO₄ 这类四面体键矢量求和几乎抵消，量到的是
  畸变不是转动。bond 模式 = `gmx rotacf -d`：第 0 帧定对、按身份跟踪、不判断断键
  （用户选 GROMACS 做法，理由：常温玻璃不断键，有出处）；默认仍为 sum
- P₁ 同一流程三分量；τ₁/τ₂ = 3（转动扩散）vs → 1（大角跳跃）是机理判据
- **测试逼出一个真 bug**：正四面体合向量剩 ~1e-16 舍入残差，绝对阈值 1e-30 放行，
  C₂ ≈ 0.99 纯噪声。sum 模式改相对阈值 |Σd| < 1e-6·Σ|d| 判抵消
- 曾以为用户提供的 Mamedov 2005（JCP 123, 124515）是 NMR 转动弛豫对比，通读后是
  ³¹P 化学位移 + Raman 的**静态结构**研究（环/链比例），与 rotcorr 无关，用户决定忽略。
  若将来要与之对比，缺的是 P–O–P 网络的环统计
- 验证：刚性旋转正四面体解析值；examples P–O bond 模式 P₁/P₂ vs 独立 numpy

### vacf / rotcorr 全原点 FFT + 梯形积分（2026-09-27 落地）

- 同 msd：`--tau` 默认整条轨迹 → 1 个原点。改 FFT，`--max-lag` 默认 N/2（= gmx `-acflen`）
- **积分**：左矩形 → 梯形（GROMACS `print_and_integrate` 源码注释「Use trapezoidal
  rule」；MDAnalysis transport-analysis `scipy.integrate.trapezoid`）。D 旧值偏大 C(0)dt/6
- **rotcorr**：6 分量 FFT 与 gmx `autocorr.cpp` 的 P2 求和逐项一致（对角 1.5、非对角 3、
  −0.5(N−m)）；ferro 另有「无邻居帧」，加 0/1 指示 χ 的自相关当分母 —— 只平均两端都
  有效的配对（MDAnalysis waterdynamics 取交集同义，但它先分子后原点两层平均，不能 FFT）
- vacf 加 `vacf_norm`（gmx velacc 默认归一化）
- 手册写明 rotcorr「键向量求和」对 PO₄ 这类对称中心几乎抵消，量到的是畸变而非转动
- 验证：单测 vs 暴力（rotcorr 含无效帧）；vacf 合成 AR(1) extxyz、rotcorr examples P–O
  两种 r_cut（全有效 / 0.5% 无效）vs 独立 numpy，csv 精度内一致
- 事故记录：一次 Python 就地改写把 `traj.rs` 截成 0 字节（`open(p,'w').write(open(p).read())`
  先截断后读），从 git 恢复后重做。该文件当时无未提交改动，未丢工作

### MSD 全原点 FFT + TOR 解包裹 + 笛卡尔分量（2026-09-27 落地）

- **起因**：CLI 无 `--tau` → 窗口 = 整条轨迹 → 只有 1 个原点，`--shift` 无效
- **算法**：Calandrini 2011（MDAnalysis `fft=True` 同源），`MSD = S1 − 2·S2`；序列先
  减均值防相消。公共部分 `md/correlate.rs`（rustfft，+4 crate）
- **解包裹**：TOR（Bullerjahn 2023 式 2，**后一帧**盒子、`⌊x+½⌋`），读论文原文核对过；
  替换 code1 的平均盒矩阵（格点视图，且无法套进 FFT）。NVT 下两者恒等
- **分量**：`msd_x/y/z`（笛卡尔，三斜下之和也等于总量），废 `msd_a/b/c`
- **参数**：删 `--shift`，加 `--max-lag`（默认 N/2，与 vacf/rotcorr 帮助页原本的说法一致）；
  `--fit-range` 的比例改为相对 lag 轴
- **验证**：单测 FFT vs 暴力 < 1e-9（NPT 三斜、非周期）；examples 两条轨迹 vs 独立
  numpy（ASE 读入 + TOR + 暴力平均）在 csv 7 位有效数字内一致。NPT 例子 R² 0.31 → 0.97
- 用户**尚未实际使用 msd**，手册 `analysis/msd.md` 按用户要求把 max-lag / TOR /
  分量三处写到可以独立排查的程度（含可复现的 numpy 参考代码与排障表）

### scripts/plot_msd.py + MSD 头部修正（2026-09-27 落地）

- **头部修正**：`meta_lines` 只留 shift/dt/elements/拟合比例；原子数、origins 与全部
  拟合结果进 `[inputs]`，新增 `t_lo t_hi points slope intercept d_err`
- **`d_err`**：gmx msd 的半窗差（两半共用中点，半窗 < 2 点给 NaN）。查过的做法：
  gmx msd（10–90% 窗、半窗差）、MDAnalysis（linregress + 双对数斜率 1 诊断）、
  pymatgen（带截距最小二乘）、kinisi（OLS 低估误差，主张 GLS/贝叶斯）
- **脚本**：拟合线由 `[inputs]` 的 slope/intercept 画，不重做拟合；图例挂在子图下方
  （带误差的标签比面板宽）；`--loglog` 的参考线锚在拟合窗口起点、不越过数据末端
- 顺带修了手册：Å²/fs → cm²/s 写成了 $10^{-16}$（应为 0.1）；Rust 示例调用的
  `write_msd` 早已不存在；上一批漏改的 3 处「optional PNG」

### 依赖精简：serde 与 `--plot` 移除（2026-09-26 落地）

逐个依赖量了「用了多少 / 带进多少传递依赖」后做了两件事：

- **serde 全删**：四个 derive 全仓无格式后端调用；三处 nalgebra `serde-serialize`
  一并摘，`Cargo.lock` 里 serde 系归零。连带删 `ferro-workflow` 的 serde、
  `ferro-analysis` 的 thiserror（声明未用）
- **`--plot` 与 plotters 移除**：`ttf` 在 Linux 要系统 fontconfig/freetype（集群常
  缺，macOS 走 CoreText 所以一直没暴露）；无 PDF backend，矢量化要再加约 40 个
  crate；出图早已由 `scripts/plot_*.py` 承担。`batch::write_all` 的返回路径只有绘图
  在用，改为 `Result<()>`。手册新增 `plotting.md` 说明出图途径与理由

运行时依赖（去重）macOS 138 → 97，Linux 目标 → 99。gr/sq/msd/angle 在两条 tests
轨迹上的 csv 与终端输出改前改后逐字节一致。

**评估后保留的**：nalgebra / clap / rayon / quick-xml / glob / anyhow / libc（重写
不划算或零传递依赖）；ndarray + ndarray-npy（手写 npy 的边角情况留维护风险，且
npz 将来要用，见上方待办）；rand（手写会改变同 seed 的洗牌顺序）。


### `ferro doc` 的终端渲染器（2026-08-26 提出，2026-09-23 定案，2026-09-24 落地）

三个提交按原计划：骨架 `fdf3446` → 表格 `92193eb` → 行内 LaTeX `8e5d80b`。
定案时的判据全部照做：tty-only、全自写净新增 0 crate、公式行内做块级不做、
行内解析是一次扫描的 tokenizer、只处理显式 `$...$` 绝不自动识别化学式。
依赖对比（termimad +37 / minimad +1 / 自写 0）与否掉 termimad 的理由见 git 历史
里本节的旧版（`abf51cd`）。

**动手时补的约束**（用户提出）：

- **Windows 要能用**，效果差可以接受。于是加了 `Style { ansi, unicode }`，Windows
  两者皆关：老 conhost 既不认 SGR，代码页也未必显示框线字符。ioctl 放 `cfg(unix)`，
  Windows 宽度退回 `COLUMNS` / 80。公式转换挂在 `unicode` 上 —— 希腊字母在那里
  同样是问号
- **CentOS 7 要能跑**。`libc` crate 只是 C 声明，版本号不决定 glibc 下限，
  `ioctl` 是 glibc 自古就有的符号；实际下限由 Rust 工具链定（目前 2.17 = CentOS 7）。
  又因 Cargo 按 semver 统一，写老版本号也会被锁到 rand → getrandom 已用的 0.2.189，故声明
  写最宽的 `"0.2"`

**被实测推翻或修正的**：

- **「代码片段不断行」**（第 1 步定的）在表格里站不住：net 的列名清单、changelog
  的整条命令是单个代码片段，把表格撑到 135 列。改为放得下就不断、比整行还长才在
  内部空格处断，朴素样式下反引号只留首尾各一个
- **「渲染后包含原文每个字符且顺序不变」这条验证对折行的表格不成立**：屏幕逐行读
  时相邻列的片段交错。改为两层 —— 宽度给足（表格不折行）时比顺序，80 列时比字符
  多重集。故意改丢一个字符的变异测试确认它能报出来
- **`\text{}` 的下标一律不转**。按「每个字母都有下标才转」的原规则，`N_\text{frames}
  \cdot N_\text{atoms}` 渲染成 `N_{frames} · Nₐₜₒₘₛ` —— 两个同类的词因为字母凑巧
  齐不齐而写法不同。词标签不是变量，保留 `_{atoms}`；变量字母照转（`wᵢⱼ` `Qⁿₘ`）
- 下标字母表在 plan 列的 `ₐₑₒₓₕₖₗₘₙₚₛₜ` 之外补了 **`ᵢ ⱼ ᵣ ᵤ ᵥ`**（Unicode 有，
  `r_i` `w_{ij}` 在手册里常见）
- 希腊字母后的空格只在紧跟字母数字时吞掉（`\Delta r` → `Δr`），否则
  `\tau \leq N` 会成 `τ≤ N`

**明确没做**：块级 `$$` 公式（判据 5）；`\overline` / `\overrightarrow` 共 6 处留
TeX 原文（组合字符在表格里会让列宽差一格）；手册里 `.md` 链接不映射成
`ferro doc <topic>`（链接目标原样显示在括号里）。

**未验证**：Windows 目标的编译（本机没装该 target，非 unix 分支只有一个返回 `None`
的函数）；Linux 产物的 glibc 符号版本（本机是 macOS，只能按上面的推理）。

### 手册英文化 + 公式统一（2026-09-23 落地）

用户决定 `docs/src` 改英文，理由是手册的读者是外人；`//` 内部注释与测试断言**保持
中文**（它们承载判据，母语更准，读者是开发者）。规则进 `CLAUDE.md`「编码约定」，
`dev/overview.md` 那段原是它的副本、改为指针。

八个提交：产物摘出 git → 规则 → 术语表 → 公式统一 → 四批翻译。
`docs/src` 25 页 4705 行现为**零中文字符、零中文标点**（894 行中文：散文 555 +
表格 290 + 代码块 49；89% 集中在 `cli-reference.md` 516 行与 `network.md` 277 行）。

定下来的判据：

- **术语表先行**（`dev/glossary.md`，65 条，用户逐条裁决）。不先定表，两页大的会译出
  两套词汇。用户的四条裁决：口径 → statement、约定 → regulation、点名拆成
  declare/name 两处（另一处「位点名」是 n-gram 断词错误）、规范半边 = canonical half
- **公式统一 `$...$`，不手写 Unicode 上下标**，但**单位符号例外**（`Å³` `Å⁻¹`
  `g/cm³`）—— 手册原本就是「公式 LaTeX、单位 Unicode」（`sq.md:38` 即
  `` $q \approx 1$–2 Å⁻¹ ``），把单位塞进 LaTeX 只会让表格更难读。`SUMMARY.md` 的
  目录标题也不动（mdBook 侧边栏不过 MathJax）
- **`docs/book/` 摘出 git**。91 个 HTML 的时间戳停在 08-13 而源文件已到 09-22，
  跟着仓走只会让每个手册提交带一批陈旧产物

**被实测推翻的**：

- **「公式在 mdBook 里已正常渲染」是错的** —— MathJax 2.7 默认定界符是 `\(...\)`，
  手册里 255 行 `$...$` 在网页上**从来没被渲染过**，一直是字面的 `$Q_n$`。新增
  `docs/theme/head.hbs` 放开 `$` 才修好；`workflow/spin.md` 的 `\(...\)` 是全书唯一
  渲染正确的写法，已统一为 `$...$`
- **中文标点是盲点** —— `，。、（）「」` 不在 `[\u4e00-\u9fff]` 里，「零中文残留」
  的检查全部漏掉它们。最后一处是替换顺序的锅：先把全角括号换成半角，`），` 那条规则
  就失配了
- **校验脚本的判据改了三次才对**（`.claude/check_tr.py`，随 TASK 文件同去留）：
  ① 行内代码逐项相等 → 把 `gr_<配对>` → `gr_<pair>` 报成违规；② 只比不含中文的项 →
  新版翻译出的英文项成了「多出来的」；③ **旧版的英文字面量必须仍然存在**（多重集包含，
  不是相等）才对。契约保护的是英文字面量，不是占位符。顺带堵了一个假阴性：在子目录里
  跑 `git show` 拿到空串，空集减空集恒过

**顺带发现未修**：`cli-reference.md:718` 写着 bader「没有 `-o`」，而同节 744 行写着
「`-o <DIR>` 收集到别处」—— 0.3.3 给 bader 加了 `-o`/`-s`，这一行漏改。翻译时**照原文
直译**没有顺手修：翻译提交要能对拍，修 bug 该是独立一个提交。

### cp2k_grep_strInfo.py 的能力合并：LAMMPS 导出 + 逐帧标量（2026-09-21）

用户在 `private/` 放了一个 626 行的脚本要求合并。**先做的是对账**：它的九项能力
里 Ferro 已有六项且实现更稳健（VASP/CP2K 读取、DeePMD 写出与 8:1:1 划分、
lammpstrj/data 导出），`outputLammpsin` 甚至引用了一个**根本没定义**的
`LAMMPSIN` 变量。真正缺的只有三样：CP2K 多文件布局、DeePMD mixed type、
逐帧标量序列。

**用户选的范围是 LAMMPS 导出 + 逐帧标量**，另两样留在优先级高那一栏。

六个提交，顺序即依赖：三斜 bug → `Frame` 字段 → 三个 reader 填温度 →
`collect` 的 `-o` 与 `.db` → `--type inspect` → 文档。

**被实测推翻或修正的**：

- **原以为缺的是 LAMMPS writer，实测 Ferro 的两个 writer 已支持三斜**，
  `cell_to_lammps` 与脚本的 `getLammpsCell` 逐行同构。真正的断点是**入口**：
  `io_dispatch` 不认 `.out`/`OUTCAR`/`vasprun.xml`
- **却发现了一处脚本对、Ferro 错的**：dump 的三斜盒子行是 `*_bound` 不是
  `xlo/xhi`，而 Ferro 读写两侧一致地错，往返测试全绿。用户当时的指示是
  「只取能力，不记录这些问题」，但这一条涉及 Ferro 本身，提出后用户判为「修，
  且记录」。**教训是对账要做到实现层，不能停在能力清单**
- **原计划的 `time_fs` 列去掉**：设计时说「step 与 time_fs 是 CP2K 锚点上现成的」，
  这话对**文件**成立、对 `Frame` 不成立 —— 两者都没被带出来。实测三条路里
  `step` 有两条、`time_fs` 只有一条，故只加 `Frame.step`
- **脚本在 VASP 路抓的温度一度被怀疑是错的**（正则 `temperature +([\d\.]+) K`
  命中的是 `EKIN_LAT` 行）。实测**它是对的**：那个括号标的是离子温度，
  `EKIN=4.268142` / `N=296` 反算得 111.56 与打印的 111.55 吻合
- **`-o` 的形态被用户改过两轮**：先定「`--type inspect` 时 `-o` 改写 inspect 位置」
  （与 bader 同构），再改成「`-o` 完全交给数据集，inspect 位置固定、给了 `-o`
  报错」。`collect` 的 `-o` 也从必填改为可选
- **后缀从 `.train` 改成 `.db`**：用户先说可能不需要后缀，最后定为全部加 `.db`

**明确划在范围外**：CP2K 多文件布局（`-pos-1.xyz`/`-frc-1.xyz`/`-1.cell`）·
DeePMD mixed type · `io_dispatch` 注册 · 抽帧参数 · `ferro info` 的「只报首尾
两帧」。脚本里**不照搬**的四处：体积用对角线乘积、CP2K 固定行偏移（只认
NVT/NPT_I）、质量表读 `~/.emacs.d` 的 elisp、倾斜折叠的除数写错。

### 原子顺序 + --format + 帮助页精简（2026-09-22 落地）

起因是实测 `collect -i 1Al/*.out` 报「different compositions」而两个化学式打出来
一模一样。落地清单见 `overview.md`「2026-09-22 的第二批」。

**被实测推翻的两条原计划**：

- **「删掉嗅探，完全依赖 `--format`」**（用户最初的要求，理由是版本多、易错、
  难维护）。查证后留下了：`sniff` 40 行且不含任何版本知识，担心的那部分全在四个
  reader 的块匹配里。留着它换来一行 NOTE 与那条一目录一格式的守卫
- **「`--format` 必填」**（我最初的推荐）。用户改为「有 format 用 format，没有就
  嗅探，探不到才报错」，旧命令行因此一条都不用改

**动手前查了本机 deepmd 环境**，两条结论都进了 `overview.md`：DeePMD-kit 加载时
自己按类型排序（所以磁盘顺序对训练无影响），dpdata 的 `append` 在 atom_types
不一致时自动排两边（所以这是领域内的既有做法，不是 ferro 的发明）。

**划在范围外**：`filter` / `merge` 的 `--format`（它们读 npy，不需要）·
`io_dispatch` 注册 `.out`（第三轮绕开）· mixed type 下的规范序（见优先级高那栏）。

### cp2kdata 的能力合并：CP2K 单点读取（2026-09-22 落地）

核心需求是**单点能数据的读取**，用于建训练集。落地清单见 `overview.md`
「2026-09-22 的一批」。

**明确划在范围外**（本轮不做，也未计划）：
- **GEO_OPT / CELL_OPT**。看着像"顺手"，其实不是：cp2kdata 里这两条路的坐标
  **不来自 `.out`**，`parse_cell_opt` 去 glob `*pos*.xyz`，`parse_geo_opt_info`
  干脆不产出坐标。做它们等于做 CP2K 多文件布局，那是另一件已划在范围外的事
- cube / pdos / Mulliken / Hirshfeld 布居 / 振动频率 / 偶极 / DFT+U 占据数。
  前两项 Ferro 已有自己的实现，后五项 `Frame` 里没有存放处
- `io_dispatch` 注册 `.out`（`convert` 仍进不去，与上一轮同样绕开）
- 2023 / 2024 的**单点** fixture。手上那两个版本的样例全是 AIMD，单点布局是
  按同版本 MD 输出推断的，因而落在"打提示"那一侧而不是"已验证"

**仍然开着的**：
- **kind 分成独立训练类型**。现在 `type_map` 恒按 element，`Fe1`/`Fe2` 合成一个。
  dpdata 默认相反。要分的话是 `collect` 加一个开关，不是改默认
- **2026 版本未经实测**，是按用户"格式沿用 2025，未修改"的说法列进已验证的。
  拿到 2026 的输出跑一遍即可确认

### 三项小待办：&Path、bader --outdir、-o 语义统一（2026-08-13 提出，2026-09-20 落地）

三件事原本互相独立，做的时候合成了一条链：路径类型统一之后，`-o` 的语义才好改。
六个提交，顺序即依赖。

**动手前推翻的原计划**：原案是「`-o` 统一为文件名、路径走 `--outdir`」，收敛到
**两种**语义（后缀 / 目录）。用户指出 `-o`/`--outdir` 本身就乱，要求改成
pathlib 那样「带目录的路径和单个文件名都收」。于是目标从「两种语义」变成**一条
规则**：

> `-o` 恒为路径。一次运行写多个产物的命令（13 个分析命令、`chg-sdf`、`bader`、
> `dataset`）指**目录**，只写一个文件的（`convert`、`job`）指**文件**。
> 批次标记腾给 `-s/--suffix`，`--outdir` 删除。

定下来的判据：

- **目录不存在时先问一句**（stderr、`[y/N]`、默认否），非交互环境直接报错并指明
  `--mkdir`。自动创建会把打错的路径变成一次「看着成功」的运行；而在脚本里连问都
  没人能答，所以非交互不能沉默地二选一
- **「只给目录」只认用户写的尾部分隔符**，不看磁盘上有没有同名目录。否则同一条
  命令行的含义会随磁盘状态漂移（今天 `out` 不存在→产物叫 `gr_out.csv`，明天有人
  建了 `out/`→写进目录），与「不按输入文件数分派两条路径」是同一条理由
- **`convert` 拒收以分隔符结尾的 `-o`**：目标格式正是从文件名推断的，目录名里没有
  它。`job` 相反，它有默认名（`job.gjf`/`job.inp`/`pw.in`），可以放进目录
- **`bader` 的报告默认写在输入旁边**，是全仓唯一不默认当前目录的命令。VASP 的
  电荷密度一律叫 `CHGCAR`，默认当前目录时 `run1/CHGCAR` 与 `run2/CHGCAR` 必然
  互相覆盖 —— 换个默认值就消掉了这个坑，比加参数让用户每次记得写更靠得住
- **`write_acf`/`write_bcf`/`write_avf` 改为返回文本**（`acf_text` 等），文件由 CLI
  写。既然要动签名，顺手清掉「分析层不碰文件系统」的唯一例外
- **reader 也跟着 writer 一起改 `&Path`**。只改 writer 会让 `ferro-io` 一半 `&str`
  一半 `&Path`，`io_dispatch` 两侧各转一次

**dataset 的默认后缀**（同一轮，用户要求）：`filter` / `merge` 的产物默认带
`.train`，划分出的三部分是 `.train`/`.valid`/`.test`，`collect` 不加（它的产物是
没筛过的原始数据）。连带三条：名字已带后缀不叠加；`merge` 遇到混合后缀**报错**而
不是默认标成训练集；`.train` 输入允许再划分（先剥后缀），`.valid`/`.test` 拒绝。

**验证**：53 个产物的逐字节对拍贯穿六个提交，只有 dataset 那步的目录名多出
`.train`，内容零差异。bader 另造了 `tests/CHGCAR_2atoms`（8×8×8 合成网格，9 KB）
补上端到端覆盖 —— 此前那条链只有内存里构造网格的单元测试。

**顺带发现的缺陷（未修）**：`bader_weight.rs` 的真空电荷取 `volchg[nvols]`，而该
数组在 weight 方法里是 1 索引的，取到的是最后一个 Bader 体积。ACF 末尾的 `Total`
因此虚高（新 fixture 上 158.98 e vs 实际 105.99 e）。`ongrid`/`neargrid` 走另一份
代码，不受影响。见 `issues.md`。

### 近距接触结构的定向构造（2026-09-01 提出，2026-09-15 废弃）

起因：MD 里原子互相穿透，训练集在阳离子对（P-P、Al-Al）短程区没有样本，势函数
在那里纯外推。原计划在 `ferro-structure` 里挑一对原子定向压缩、造扫描结构送 DFT。

**废弃理由**：改由 dpgen 主动学习去探索这些构型，Ferro 不做定向构造。原计划里
「压近的帧会被 `filter` 的 `f_max` 删掉」这条冲突随之消失 —— dpgen 产出的结构
正常流程下不再过二次筛选。

留下的事实（重开或做别的扰动时不必重查）：

- **六个工具（ASE `rattle`、pymatgen `perturb`、dpdata `perturb`、dpgen
  `create_random_disturb.py`、hiphive `mc_rattle`、ASE GA）无一做定向压缩**，全是
  各向同性随机抖动。管最小间距的只有 hiphive 与 ASE GA，方向是**阻止**近距接触
- **dpdata 的 `perturb` 是 dpgen 那个脚本的修正版**：dpgen 的方向按立方体归一化、
  偏向 ⟨111⟩，模长 U[0, dmax) 也不是球内均匀；dpdata 两处都改对了。要统计微扰用
  dpdata，别用 dpgen 的 `pert_atom`
- **排斥壁两条路配合使用**：ZBL 混入管 d→0（NEP 内置 `zbl`，DeePMD 走 `use_srtab`
  + `sw_rmin`/`sw_rmax`），数据管约 1.5–2.7 Å
- **存疑**：当时的锚点表（P-P 2.73 Å、Al-Al 2.83 Å 以下无样本）标为「训练集见过的
  最近 r」，实测对象却是 `tests/43Z43P15A_NPT_5.lammpstrj`（5 帧）。引用前先在
  真实训练集上重测

### VASP AIMD 读取：OUTCAR + vasprun.xml（2026-08-27，0.3.2）

原计划只做 OUTCAR，用户要求**两个都做**。实测 `quick-xml` 与 `roxmltree` **各自只
净新增 1 个 crate**（零传递依赖），代价不是问题；选 `quick-xml` 的**流式**是因为
DOM 的内存约为文件的 5–10 倍，而这条链的输入本来就是几百 MB。

**被查证推翻的原计划**（两条，都在动手前）：

1. 上一轮写的「dpdata 取 `energy(sigma->0)`」**是错的**。`dpdata/formats/vasp/
   outcar.py:118` 取的是 `free  energy   TOTEN`，用户的 `private/dp_makedataliu.py`
   也是（`lines[-3].split()[-2]`）。理由本身也站得住：力是自由能的导数
2. 原计划「同目录既有 OUTCAR 又有 vasprun 时优先 vasprun、回落 OUTCAR」被用户
   否掉 —— 改为**由用户指名文件**（跟 CP2K 端一样），程序不做目录级偏好判断。
   于是分派改成读头几行认横幅，而不是一张扩展名特例表

**定下来的判据**：

- **格式按内容嗅探,不按名字**。VASP 写的叫 `OUTCAR`（无扩展名），用户手上那份叫
  `50Z50P_0.970_3000K.outcar`，`.out` 又谁都可能用 —— 按名字判必然要开一串特例，
  而横幅是唯一的
- **同目录混格式报错**。真实运行目录里两个文件并存且记同一批帧，按 collect
  「同目录 = 同一次运行的分段」拼接会让数据集**悄悄翻倍**：成分一致、两个文件
  各自也都读得通，不会有别的症状
- **晶胞逐帧各读各的**，缺块的帧丢弃而不继承上一帧。定胞下继承与不继承逐位相同，
  这个 bug 只会在第一次跑变胞时显形 —— 与 `array_order.rs` 的转置同族
- **`in kB` 是 VASP 自己的 Voigt 顺序 `XX YY ZZ XY YZ ZX`**，与 extxyz 规格的
  `XX YY ZZ YZ XZ XY` 不同。三个独立来源一致（dpdata 的下标映射、用户脚本的
  `[[0,3,5,3,1,4,5,4,2]]`、ASE 的 `[[0,1,2,4,5,3]]` 重排）
- **体积用 `|det|`**，不用三个对角线相乘（用户脚本是后者，只对正交胞成立）
- 两条 VASP 路径的**收敛判据不同**（OUTCAR 有 VASP 自己的 `EDIFF is reached`，
  vasprun 只能数 `<scstep>` 与 `NELM` 比），故判据随 stats 打进报告 —— 否则同一次
  运行换个来源、丢帧数不同会没人说得清

**实测发现的两处坑**（都不会报错，只会给出错的数）：vasprun 的 `<energy>` 每个
`<scstep>` 都有一份，收敛值是最后一个（取第一个得到 33072.96 这种数）；
`<array name="atoms">` 的 `<field>` 表头也是文本，当成元素会得到 298 个物种、
每帧判 incomplete、整个文件读成空 —— 后者是实际写代码时踩到的。

**验证**：`tests/` 两份从真实文件裁出的 fixture（OUTCAR 340 KB、vasprun 188 KB，
dpdata 仍能正常读）；`#[ignore]` 的全量对拍（2000 帧 / 425 帧）与 dpdata 逐位一致，
virial 差 8.5e-9 相对量 —— dpdata 用的 eV/Å³→GPa 常数是 160.2176621，`units.rs`
用 CODATA 2018 的 160.2176634。用户那份 OUTCAR 本身是**两段拼的**
（1..1575 接 1..425，第 1576 步被中断），正好验到重启计数与截断帧丢弃。

### GPUMD/NEP 导出：`--type` + train/valid/test 划分（2026-08-26，0.3.2）

**原计划是加一个 `ferro dataset export` 子命令，被用户否掉** —— 改为在 `filter` 与
`merge` 上加 `--type`，划分标签也加在这两处。事后看这个改法更对：NEP 的 `train.xyz`
本来就是 extxyz，ferro 的 writer 早就能写，真正缺的只是「读 DeePMD system 目录」，
而 filter/merge 本来就在读；新命令会把同一件事再写一遍。

定下来的判据：

- **一个 system 一个 `.xyz`**，不并成一份 train.xyz。NEP 要的是单份，但那是一句
  `cat`；分开保留的是「抽掉某个来源不必重跑」的能力
- **划分产物用目录名后缀** `.train/.valid/.test`，不造 `train/` 子目录 ——
  `SPLIT_SUFFIXES` 早就在仓里，merge 一直在继承它，这是 dpgen/dpdata 的约定
- **成员取自打乱序**：直接切尾巴的话 test 全是轨迹末尾一段连续状态。各部分内部
  再排回帧序，同 seed 下逐字节可复现
- **比例向上取到至少 1 帧**：给了比例却拿到 0 帧，等于静默地没有验证集
- **一个 `--ratio 8:1:1` 而不是两个 ratio 参数**（用户第二轮要求）。数字是
  **权重不是分数**，`8:1:1` 与 `80:10:10` 同义，不必凑成和为 1；两段写法 `9:1`
  按 `train:test` 解释 —— NEP 要的就是这两份文件，而位置歧义靠每次运行按名字
  打出三部分的帧数来自证，第一行输出就能看见解释对不对
- `nep` 与 `extxyz` 两个值只差一个应力键（`virial=` / `stress=`）。留两个值而不是
  一个，是因为「按下游训练框架命名」比「按文件格式命名」更接近用户提问的方式

**没做、也不该在这里做的**：`weight` / `dipole` / `pol` / `bec` 这些 NEP 可选键。
`Frame` 装不下它们，硬塞要开一条绕过 `Trajectory` 的通道，与「额外键搬运」是同一件
待办。

**已知局限（写进帮助页与手册，不是缺陷）**：同一段 MD 的相邻帧高度相关，帧级 test
仍偏乐观；诚实的估计要按整个 system 划分，做法是把来源分开筛、留一个不动。

### extxyz 的 stress/virial 修正（2026-08-26，0.3.2）

四处静默缺陷，共同根因是**在没有依据的地方替用户猜了一个约定，且猜错不报错**：
`stress` 取不到就回落取 `virial`（差一个体积因子）· 两侧都没做 ASE（正 = 拉伸）与
`Frame::stress`（正 = 压缩）之间的变号 · 九个数按行优先处理而 ASE 文档说列优先 ·
`Properties` 的力列只认复数 `forces`，读 GPUMD 的 `train.xyz` 会丢掉全部受力。

**被查证推翻的原计划**（三条，都发生在动手之前）：

1. 原计划「读到 `virial=` 就报错，因为符号约定无法核实」。实测 dpdata 1.0.2 的
   `virials = -volume * stress_ase` 双向换算 + GPUMD 手册 + extxyz 规格的
   `virial -> stress` 乘 `-1/cell_vol`，三方一致，符号完全可定 —— 报错等于明知
   怎么读却拒绝读
2. 原计划「按列优先读写以对齐 ASE」。查到规格要求该张量**对称**（"fail if not
   symmetric"），而对称下行/列优先逐位相同 —— ASE 说 Fortran order、GPUMD 手册
   拼成 `vxx vxy vxz vyx ...`，两家描述相反却从无人报 bug，正是这个原因。改为
   **检查对称性**，不站队
3. 原计划「6 分量按标准 Voigt 收下并告警」。用户判断：让用户重排与重新生成 9 分量
   的工作量几乎没差别，那就取最稳的结果、正确性交回给用户 —— 改为**拒收**

`Lattice` 全程未动：实测 ASE 写出的九个数就是三个晶格矢量依次排开，Ferro 原本正确。

### 帮助页精简：其余各页（2026-08-26，0.3.2）

`dataset` 三页的模板推到了全部命令。**实际范围与「22 页」的估计不同**：

- **只有 8 页超 40 行**，其余 16 页本来就短，只需补一行手册指针
- **`ferro net` 的帮助页不在 `help.rs`**（在 `cmd/net.rs` 的 `HELP_EXTRA`，
  紧挨着它需要的 argv 剥离逻辑），所以此前**防漂测试完全没覆盖它**，而它
  有 10 个选项。已纳入
- **`convert` / `info` / `bader` 没有手册专页** —— 原以为要新写三页，实测
  `cli-reference.md`（863 行）已按命令分节覆盖了全部参数，缺的只有
  `--help`/`--input` 这类普适项

最后一条导致了唯一的设计改动：**`ferro doc` 支持按小节寻址**。`Page` 加
`section: Option<&str>`，从该 `##` 标题取到下一个同级标题。`ferro doc convert`
于是出 99 行而不是整本 863 行，手册也不必拆成一堆按命令的小文件、让「唯一
一份完整参考」碎掉。小节标题改名而表没跟上时回落到整页，另有测试钉住每个
`section` 真实存在。

**收尾时仍有 4 页超 40 行，且都不该再砍**：

| 页 | 行 | 大头是什么 |
|---|---|---|
| 顶层总览 | 61 | 它是入口地图，不是命令页 |
| `net` | 54 | 10 个选项 + 6 张输出表 |
| `convert` | 59 | 其中 **27 行是 `supported_formats()` 生成的**格式表，手写只 32 |
| `job -s cp2k` | 47 | 23 个参数，值域枚举本身就是参数表 |
| `dataset filter` | 43 | 18 行参数表 |

**判据：参数表（含值域枚举）是唯一必须完整的一段，不为压进 40 行去砍它。**
40 是目标不是硬上限。

### collect 语义重做 + filter 报告落盘 + ferro doc + 帮助精简（2026-08-26，0.3.2）

四个提交，顺序即依赖：collect 语义 → filter 报告 → `ferro doc` → 帮助精简 +
手册对账 + 防漂测试 + 本轮记录。第四步必须最后，它照着前三步的**实际**行为写。

被实测推翻或修正的：

- **「用 `_` 串联多级目录」是用户的初始提法，商议中被他自己换掉**：改为
  **保留目录结构**（`sets/a/md` 而不是 `sets/a_md`）。理由是 `filter` 的
  `-o` 本来就按相对路径重建，`find_systems` 与 `merge` 也都是递归的，
  压平反而是这条链上唯一的例外
- **「撞名报错」整条判断消失**：新规则下撞名就是「该合并」的定义，
  原来那个 `bail!("would both write the system directory ...")` 删掉。
  命名规则与分组规则合成同一条，所以这两件事不能拆成两个提交 ——
  中间态是「撞名报错但名字已变」，没法对拍
- **去重被用户判为不必要**：重启重叠段位置速度相同，能量力也相同，
  不稀释结果；重启间隔通常 5 步以内。改为不去重但**打出 step 区间**，
  让这个前提保持可检验
- **「时间序没法区分」这句原判断是错的**：`MD| Step number` 正是解析器
  挂一切的锚点，xyz 注释行还带 `time =`。两个值当时都读到了但都没存
- **「诊断只在只读模式算」的省时理由不成立**：实测 1110 帧 / 302 原子
  挂钟 0.23 s（带）vs 0.28 s（不带），rayon 跑满六核。改为恒算，
  `--no-diagnostics` 开关不必加
- **报告平铺还是塞子目录，用户的直觉对且理由更硬**：`expand_dirs` 只收
  `is_dir()`，平铺的 csv 会被后续 `merge -i clean/*` 自动滤掉，而
  `report/` 子目录反倒会被收进去当 system 候选
- **`clap_mangen` 解不了这个问题**：它从 clap 定义生成 man page，只有
  参数表 —— 正是要精简掉的那部分，散文一句都带不出来。`cargo doc` 是
  rustdoc（API 文档），也不是一回事。故 `ferro doc` 自己做
- **防漂测试当场抓到 19 处真漏**，远超预期的 1 处（`--shuffle`）。
  三次收紧判据才落到可用：短名也算写了；允许交叉引用别的命令的参数；
  跳过「there is no --x」这类否定陈述。剩下 2 处是有意不写，进白名单
