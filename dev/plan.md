# 后续计划

> 归档只保留**判据**与**被实测推翻的原计划** —— 「做了什么、怎么改的」翻 git 历史
> （提交号已列出），「现在是什么样」看 `progress.md` 与 `docs/src/`。

## 优先级高

### 2026-10-03 四路审查发现（对照 ASE / MDAnalysis / dpdata / pymatgen / CP2K 源码；已修 A1–A15、B-1–B-3、C-D1、D-S1–D-S5，C-D2 搁置，其余未修）

四个 Fable 5.1 子代理按领域只读审查（A 格式读写 · B 轨迹分析 · C network/dft/ml/core ·
D cli/structure/workflow/跨 crate），对照库在 `~/.miniforge3/envs/deepmd`（ase 3.29、
dpdata 1.0.2、MDAnalysis 2.10、pymatgen 2026.7.16），CP2K 口径对照 `examples/cp2k-2026.2/src`。
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

**严重**（静默错数据 / 错物理，或 panic 打断批次）

| # | 位置 | 问题 | 复现 | 修法 |
|---|---|---|---|---|
| ~~A1~~ **已修**（`acb282d` 数据表、`7638dc0` 查表、`1a76d49` 接入；判据见 `issues.md`「CIF 空间群符号展开」）| `readers/cif.rs:334-335 collect_symops` | 无 symop 循环一律当 P1，H-M / IT 号写着非 P1 也不报；CIF2 点号标签 `_space_group_symop.operation_xyz` 同路 | NaCl CIF 只写 `_symmetry_space_group_name_H-M 'F m -3 m'` + Na1、Cl1：ferro 2 原子、密度 0.54 g/cm³；ASE `Cl4Na4` | 无 symop 而 H-M / IT 号非 P1 时 bail；点号标签加进 `TAGS` |
| ~~A2~~ **已修**（`c8ee8c8`；顺带修了同函数里的取向错位：CRYST1 只存六参数，非标准取向的胞原样写坐标，读回错位 —— 现经分数坐标转到 a 沿 x、b 在 xy 平面，同 ASE `standard_form`；对拍 ASE 写出逐位相同。**遗留**：ASE 把 `ENDMDL` 后的 `END` 读成第 N+1 个空帧，ASE 自己不写 `END`，未改）| `writers/pdb.rs:19-24` | CRYST1 只按第 0 帧写一次，reader 套到所有帧 | `convert -i tests/43Z43P15A_NPT_5.lammpstrj -o npt.pdb`：1 个 CRYST1（ASE 写 5 个）；第 4 帧体积 27886.6 → 29687.3 Å³ | 每个 MODEL 前写本帧 CRYST1（reader 已认逐帧） |
| ~~A3~~ **已修**（`a34b488`；无胞帧盒子取包围盒每侧外扩 1 Å，坐标原样、lo≠0；未在真 LAMMPS 里读过 data —— 本机无 LAMMPS）| `writers/lammps_dump.rs:48-54,69`、`writers/lammps_data.rs:37-50` | 无胞帧取包围盒尺寸写成 `0..L`、坐标不平移、边界硬写 `pp` | 水分子 xyz → `.lammpstrj` → `info`：Volume 0、β=γ=NaN、PBC 全 true，原子在盒外 | 无胞时坐标减 min、盒子加余量、写 `ff`；或报错要求先给胞 |
| ~~A4~~ **已修**（`9ebb707` 逐字符解析，规则同 ASE `key_val_str_to_dict`，一并支持 `[]` `{}` 括号值与 `\` 转义；有意与 ASE 不同两处：未闭合报错、`a="" b=1` 读作两键，见 `parse_comment` 文档） | `readers/extxyz.rs:203-204 read_value` | 未闭合引号 `&inner[end+1..]` 越界 panic | 注释行 `Lattice="5 0 0 0 5 0 0 0 5 Properties=…`（缺右引号）：退出码 101 | 无闭合引号 bail 并点名帧号 |
| ~~A5~~ **已修**（同 A4）| `readers/extxyz.rs:180-196 parse_comment` | 不带 `=` 的裸键（ASE 读作 True）与下一个键名粘连，Lattice/Properties 被静默丢弃 | `energy=-1.5 is_relaxed Lattice="5 0 0 0 5 0 0 0 5" Properties=…`：ferro `Cell: none`；ASE 5 Å 立方 | 按空白切词，无 `=` 的词记作标志 |
| ~~A6~~ **已修**（`e46a5ab` 拒收；完整支持仍是下面「DeePMD mixed type」一节）| `readers/deepmd.rs:62,107` | mixed type system（`type.raw` 全 0 + `real_atom_types.npy`）读入即全成 `type_map[0]`；本文件「DeePMD mixed type」一节只写了「未实现」 | dpdata `to('deepmd/npy/mixed')`（源 `tests/vasp_OUTCAR_2frames`）→ `dataset filter --type extxyz`：594 行全是 O，退出码 0 | 见到 `real_atom_types.npy` 就 bail，直到真正实现 |
| ~~B-1~~ **已修**（`bb25ba7`；`check::invertible_cells` 在 vanhove / cluster SDF / cube_jump 入口逐帧查，`calc_cluster_sdf` 改返回 `Result<Option>`。实测 c=0 slab 过全部 13 个 CLI 分析，只有这两个 panic）| `md/vanhove.rs:213`；`map sdf` 经 `cube_sdf.rs:354,400` → `ferro-core/src/cluster.rs:106`；`cube_jump.rs:149` | 奇异胞 `expect("cell is non-singular")` panic。ASE 二维材料约定（c=0、pbc TTF）即触发；gr/msd 等同一文件是「singular 跳过 + 退出码 1」 | c=0 slab extxyz 进 `traj vanhove --dt 1`：退出码 101 | 入口逐帧查可逆，或 `expect` 改 `?` |
| ~~B-2~~ **已修**（`a5d1be8`；报错措辞改为「must be a finite number …」）| `ferro-analysis/src/check.rs:22 non_negative`、`:28 ordered` | 只有 `positive` 查 `is_finite`，两者放行 `+inf` | `traj sq --q-max inf`：`sq.rs:130` panic（debug；release 回卷成 0 → 只有表头、退出码 0）。子代理另报 `vanhove --r-max inf`、`map sdf --sigma inf` 退出码 134、`--padding inf` 写出 −inf 原点的 cube（主会话未复现） | 一律要求有限 |
| ~~B-3~~ **已修**（`98b5e97`；**用户否了「取 q→0 极限」**：有限盒子的 S(0) 无物理意义，改为 q-min 必须 > 0；并加 `[inputs]` 的 `q_trunc = 2π/r_max` 与低于它时的告警。查证的行业做法与出处见手册 `sq.md`「The low-q limit」）| `md/sq.rs:167` | q=0 捷径返回 1.0，而公式极限是 $1+4\pi\rho\int r^2(g-1)\,dr$ | `traj sq -i tests/70Z30P00A_NVT_5.lammpstrj --q-min 0 --q-max 0.03 --dq 0.01`：q=0 各列 1.000，q=0.01 处 total_xrd −0.45、O-O 0.600 | q=0 用 $\sin(qr)/q \to r$ |
| ~~C-D1~~ **已修**（`4a34108`；用户要求先查行业做法，结论记在下面「自旋多重度」一节。范围扩到同一函数的 4d/5d：一律低自旋、d⁸ 取平面四方 0；氧化态凑不平单独告警。镧系假定 4f 在价层，加告警）| `ferro-core/src/spin.rs:200-213` | 镧系 `group_number` 为 None 落进主族分支，算出垃圾未成对数且无告警；CP2K/QE 默认 auto-spin | 萤石 CeO₂ 进 `job -s cp2k`：`MULTIPLICITY 9 / UKS`（应 1）。Ce₂O₃ 6（实 2）、Gd₂O₃ 6（实 14）、EuO 1（实 7） | f 区分支 `n_f = Z−54−ox` 再 `hund(n_f, 7)`；至少回退奇偶下限并告警 |
| C-D2 **搁置**（2026-10-08 用户定：随 Bader 重写一并修，同 S-B / S-C）| `dft/bader.rs:205` vs `bader_weight.rs:98,145` | weight 的 `volchg` 是 1 索引，`bcf_text` 按 0 索引读，**BCF 电荷列整体错一位**。「bader weight 的真空电荷取错」一节只说了 vacchg，不完整 | `ferro bader -i tests/CHGCAR_2atoms -m weight`：BCF 体积 1 电荷 0、体积 2 为 53.000（实为体积 1 的）；`Vacuum charge` = 原子 2 电荷 52.99 | 建 `BaderResult` 前转 0 索引；或按布局建议 5 复用 grid 的函数 |
| ~~D-S1~~ **已修**（`a92e4b9`；用户定：`--pressure P`（bar）取代 `--barostat`，给出即 NPT_F，不给为 NVT；只拦非有限值，负压允许）| `ferro-workflow/src/cp2k.rs` `write_motion` 压浴 | `PRESSURE 1.01325E+05 # bar`，CP2K 该关键字单位就是 bar → 约 10 GPa | `job -s cp2k -i examples/30Z70P.cif --task md --barostat` | 1 atm = `1.01325`；或加 `--pressure` |
| ~~D-S2~~ **已修**（`bd21bc2`；用户注：杂化泛函一般用不到，作功能完善）| `cp2k.rs` HSE06 分支 | `&HF` 缺 `&INTERACTION_POTENTIAL POTENTIAL_TYPE SHORTRANGE / OMEGA 0.11`，按 COULOMB 全程 25%，得到不存在的泛函且能跑 | `--functional hse06`：无 `INTERACTION_POTENTIAL`；对照 CP2K `tests/QS/regtest-hybrid-3/CH3-hybrid-HSE06-lsd.inp:36-39` | 补段 + 断言测试 |
| ~~D-S3~~ **已修**（`bd21bc2`；系综按恒温器 × 压强查表：csvr/nose → NVT/NPT_F，none → NVE/NPE_F，langevin → `ENSEMBLE LANGEVIN` + `&MD/&LANGEVIN`，与 `--pressure` 同给报错。**顺带修**：原 langevin 写进 `THERMOSTAT%TYPE`，CP2K 只认 NOSE/CSVR/GLE/AD_LANGEVIN（`input_cp2k_thermostats.F:237`））| `cp2k.rs` MD 分支 | `--thermostat none`（帮助写 NVE）仍写 `ENSEMBLE NVT`，CP2K 默认恒温器 NOSE | `--task md --thermostat none` → `ENSEMBLE NVT`（`input_cp2k_thermostats.F:234-238`） | 写 `ENSEMBLE NVE`；与 barostat 同开时报错或 NPE |
| ~~D-S4~~ **已修**（`e78e609`；沿用 `common.read.units`。**未改**：`inspect.rs` 是 CP2K 产物导出、无 `--units`，恒写 metal，属有意的固定口径）| `ferro-cli/src/cmd/net.rs export_labelled` | `--export-traj lammpstrj` 恒以 real 单位写出，无视 `--units`；`inspect.rs` 恒写 metal，两出口口径相反 | 3 原子 dump vx=1.0，`net --units metal --P-O=2.0 --export-traj`：产物 vx=0.001 | 沿用 `common.read.units`，同 `convert` |
| ~~D-S5~~ **已修**（`bd21bc2`；第 8 轮对照 CP2K 源码时新发现）| `cp2k.rs` CSVR 分支 | 写 `TIMECON_CSVR`，CP2K `&CSVR` 只有 `TIMECON`（`input_cp2k_thermostats.F:580`，无别名）；默认 `--task md` 即走此路，推断 CP2K 解析即中止（本机无 CP2K，未实跑） | `job -s cp2k --task md` | 改 `TIMECON`；一次性脚本把 21 份生成输入的全部小节 / 关键字名对照 CP2K 源码定义，仅此一处不存在 |

**中**

格式读写（A）：

| # | 位置 | 问题 / 复现 | 修法 |
|---|---|---|---|
| ~~A7~~ **已修**（`d7876a1`；VASP 6.4.2 HDF5 编译才写 `标签/哈希`，6.4.3 起官方称已解决；5.4.4 纯符号、VASP 4 无元素行均有测试钉住。顺带：`chgcar.rs` 坐标类型只认 d 开头，改为同 POSCAR 的 C/K 规则）| `readers/vasp.rs:57`、`readers/chgcar.rs:57` | VASP 6 元素行带 POTCAR 哈希 `Na_pv/6a2f546d` 原样当元素名（ASE 读 Na） | 照 `ase/io/vasp.py:249-250` 取 `split('/')[0].split('_')[0]`，两处同改 |
| ~~A8~~ **已修**（`714cac6`）| `readers/cif.rs:480` | 任一 data 块无胞（常见 `data_global`）整份报错；ASE / pymatgen 跳过 | 跳过无 `_atom_site_` 循环的块，全跳过再报错 |
| ~~A9~~ **已修**（`714cac6`；2026-10-08 用户定：无胞结构只写 xyz —— 改为 **writer 报错**，同 POSCAR writer，读侧不改）| `writers/cif.rs:73-99` vs reader | 无胞帧写成只有 Cartn 的 CIF，ferro 自己读不回（ASE 读成 pbc=False） | 无胞参数但有 Cartn 列时按非周期读 |
| ~~A10~~ **已修**（`a6740d5`；用户定按 PDB 列对齐从原子名推元素，不照 ASE 先试两字母。顺带：77–78 列大写 `FE` 原样当元素名，改为规范成 `Fe`）| `readers/pdb.rs:45,82,92-94` | 77–78 列为空时元素为 `""`；坐标坏或行短于 54 时原子静默丢弃；多字节字符切片 panic（推断） | 空元素回落原子名；坏坐标 bail 点名行号；按字节切 |
| ~~A11~~ **已修**（`276b465`；用户定整行比较，`&KIND` 带 kind 名故比首词；剥注释后补 trim_end。新发现：`&CELL_REF` 的 ABC 覆盖 `&CELL`，`parse_cell_section` 改为只读本层。regtest `H2O-meta.inp` 实读通过）| `readers/cp2k.rs:66,74` | 段名前缀匹配：`&COLVAR` 里的 `&COORDINATION` 先命中 `&coord`，合法输入报错（regtest `QS/regtest-gpw-2-2/H2O-meta.inp` 等 ≥9 份） | 段名按词精确比较；`&cell` 同改 |
| ~~A12~~ **已修**（`307cd4d`）| `readers/qe.rs:175,29,152` | namelist 用空白分隔（`ibrav=0 nat=3 ntyp=1`）时整串成了一个值：`nat` 校验静默跳过、非零 ibrav 不被拒、`starting_magnetization` 丢失 | 按 `key=value` token 解析；ibrav / nat 解析失败报错 |
| ~~A13~~ **已修**（`d734593`）| `readers/lammps_data.rs:60-61` | 空文件 `lines[i]` 越界 panic | `get` + bail |
| ~~A14~~ **已修**（`a34b488`；`cell_to_lammps` 把相对最大边长 1e-10 以下的倾斜置 0，两份同改）| `writers/lammps_dump.rs:56`、`writers/lammps_data.rs:46` | 三斜判定 `!= 0.0`，cos 90° 残留 1e-16 使正交盒写成三斜（`convert tests/70Z30P00A_NVT_5.lammpstrj -o x.lammpstrj` 即中） | 倾斜 / 边长 < 1e-10 置 0 |
| ~~A15~~ **已修**（`a34b488`；只改 dump，data 文件无边界标志）| `writers/lammps_dump.rs:64,69` | 边界恒写 `pp pp pp`，不看 `frame.pbc`；reader 却尊重标志（TTF slab 往返成 TTT） | 按 pbc 写 `pp`/`ff` |

轨迹分析（B）：

**B 组进度**：B-4、B-5、B-7 已修，B-6 搁置（2026-10-08）；B 组表内全部处理完。下一步 D-M17；D-M4、D-M6–12 挂起作专项。照例逐条先给方案。

| # | 位置 | 问题 / 复现 | 修法 |
|---|---|---|---|
| ~~B-4~~ **已修**（NPT fixture 实跑 5 帧 / 2004 原子，与 `ferro info` 一致）| `cmd/traj.rs:740`（vanhove）、`run_angle` | `[inputs]` 的 `frames` 填 `r.r.len()`（bin 数）；angle 的 `atoms` 填 `r.elements.len()`（实 2004 写 4） | 结果结构体加 `n_frames`；填原子数 |
| ~~B-5~~ **已修**（用户选默认改 N/2：同 msd/vacf/rotcorr/bondlife 的 `--max-lag` 与 gmx `-acflen`；pymatgen `VanHoveAnalysis` 是另一套设计 —— 扫一串 lag、每个 lag 固定 50 个原点。破坏性，进 changelog）| `traj.rs:171`、`help.rs:644` | vanhove `--tau` 帮助写「default: half trajectory」，实为 `n_frames−1`（与手册一致）。`help_sync` 不比默认值 | 帮助页先对齐实现，或随「vanhove 默认 1 个原点」一节改默认 |
| B-6 **搁置**（2026-10-08 用户定：该功能后续会更新，届时一并修；候选思路：组内枚举、匈牙利 + Kabsch 迭代如 ArbAlign、距离指纹配对）| `cube_sdf.rs:445-478 heap_permutations` | 同签名族对**全部 P** 做 n! 枚举，每排列一次 SVD；docstring（`:411`）说「组内」，代码不分组。10 元环 1.2 s、11 元 10.7 s、13 元外推约 28 min，偏磷酸盐 Q2 长链常见 | 只在同标签组内置换并设上限，超限退 Hungarian 并告警 |
| ~~B-7~~ **已修**（`GrResult.r_max_used`，`params.r_max` 留请求值；未截断的 gr / sq 产物逐字节不变，三斜 fixture 只有头部一行由 4.8136 变回请求值 10.005）| `gr.rs:495`（`params.r_max` 在 `:417` 被覆盖） | 头部写截断值却标「requested」，批内第一个文件的截断值冒充全批请求值；测试 `test_meta_lines_report_clamped_rmax_and_composition` 钉住的是错误语义 | 请求值留在 params，截断值另设字段 |
| ~~—~~ **已修**（`CubeJumpParams::validate`，`calc_cube_jump` 入口也调，接进 `check.rs` 测试宏）| `cube_jump.rs` | 无 `validate()`（M3 那批漏了）：`nx=0` 在 `:95` 除零；`threshold ≤ 0` 全算跳跃；`:180`「与 msd.rs 一致」已过时（msd 已改 TOR） | 补 validate 接进 `check.rs` 测试宏 |

network / dft / ml / core（C）：

| # | 位置 | 问题 / 复现 | 修法 |
|---|---|---|---|
| C-D3 **搁置**（2026-10-08 用户定：随 Bader 重写一并修，同 C-D2）| `dft/bader.rs:174-178`；`cli-reference.md:800`、`dev/bader.md` §10 | ACF 的 `MinDist` 实为「极大值到原子距离」的最小值（原子 1 输出 0），不是 Henkelman 的「到 Bader 表面最小距离」；`bader.md` 称外部工具可按 Henkelman 格式解析，实测 ASE `attach_charges` 与 pymatgen `_parse_acf` 都失败（表头 `—` 非 ASCII、无 `----` 分隔、无 `VACUUM CHARGE:`） | 实现表面距离并恢复 Henkelman 版式；或改手册写明自有格式与列含义 |
| ~~C-D4~~ **已修**（判据开着时非有限峰值判坏，峰值保留 NaN 供报告；连带修了下面「轻」C 里的 `-f nan` / `--oo-min nan` / `--al6 nan` 与诊断下溢）| `ml/filter.rs:219-239` | 力 / 应力含 NaN 的帧通过筛选进训练集（`fold(0.0, f64::max)` 吞 NaN，`m > f_max` 对 NaN 为假） | 非有限一律判坏 |
| ~~C-D5~~ **已修**（`deac6ae` 取最长零平台中点、移入 `md/shell.rs`；`e9a88a3` 改在 0.10 Å 滑动平均上找极小，实测 2.05 → 2.57 Å。**遗留**：手册 `dataset/filter.md` 的「2.45 Å」与三对极小表是旧算法在参考体系上测的，数据不在仓内，未重测）| `ml/geometry.rs:224-236` | `--al6` 自动截断取峰后 3 Å 窗口内第一个严格最小 bin，g 恒 0 的平台上取到第一个零 bin。`collect tests/cp2k_md_3frames.out` → `filter --al6`：cutoff 2.05，cn5=2 / cn6=46；numpy 截断 2.2–2.6 得 cn6=48 | 平台取中点或末端；补零平台测试 |
| ~~C-D6~~ **已修**（O–O 分位与实际 Al–O 截断改作 `[inputs]` 列，表头不再写；`diagnostics::quantiles` 拆出，只读模式顺带补打分位行，原先只在 csv 表头里）| `dataset.rs:1063-1068,1311`、`diagnostics.rs:108`、`table.rs:166` | 多 system 报告 `concat_union` 只留第一份 meta，逐 system 统计（O–O 分位、自动 rcut）冒充全批 | 统计改数据列，或走 `Summary::note` |

入口 / 结构 / 工作流（D）：

**QC 输入生成组**：D-M4、D-M6–D-M12（CP2K / QE / Gaussian 输入生成）**挂起作专项**（2026-10-08 用户定：大工程，先记录）。动手前须逐条对照 `examples/cp2k-2026.2/src` 与 QE、Gaussian 官方文档核实——D-M10 等几条只是推断；修完应能让生成的输入真跑过一次各软件的输入解析。

| # | 位置 | 问题 / 依据 | 修法 |
|---|---|---|---|
| ~~D-M1~~ **已修**（`check_frame_spacing` 下沉为 `ferro_analysis::md::check_frame_spacing`，CLI 五处与 Python `msd` 共用；`dt` 必填；wheel 实测：NVT fixture 步距 12000/13000 不均被拦，CLI 同句）| `ferro-python/src/analysis.rs msd` | 仍默认 `dt=1.0`、不做 `check_frame_spacing`，与 CLI 已改必填的口径漂开 | `dt` 必填；间隔检查下沉共用 |
| ~~D-M2~~ **已修**（用户定保留不拉伸的填充式拼接，不做 ASE `stack` / pymatgen `from_slabs` 的共格拉伸；两块整体平移、胞中心面内居中，gap 按面间距，界面不平行（>1°）报错。原 13 个测试中正交胞的结果不变，三斜测试改按法向高度断言）| `ferro-structure/src/merge.rs` | ① 只给 B 居中，A 窄时不居中（与文档不符）；② B 沿自身单位矢量平移、胞用 A 的，倾角不同时剪切错位；③ gap 沿矢量量，与 `add_vacuum` 的垂直间隙口径不一 | 两块都按新胞分数坐标放；gap 按面间距 |
| D-M17 | `ferro-structure::find_clusters` → `classify_frame` → `network_type.rs:353,394`（2026-10-07 修 B-1 时发现） | 奇异胞（c=0 slab）在 `expect("cell must be non-singular")` panic。CLI `net` 在前面已拦，Python 绑定直调仍会 panic；`dft/chg_sdf` 经 `process_frame` 同路（cube 文件的胞，推断难触发） | 入口查可逆（同 `check::invertible_cells`；它是 ferro-analysis 私有，structure 那边要么自查要么下沉 core） |
| ~~D-M3~~ **已修**（nan / inf / ≤0 一律报错）| `cmd/net.rs parse_pairs` | `--P-O=nan` 通过校验：P 的 cn=1314、全 Q0，退出码 0 | `is_finite() && > 0` |
| D-M4 | `cp2k.rs write_force_eval` | NPT 不写 `STRESS_TENSOR`，CP2K 启动报错（`md_run.F:331-344`） | NPT 时写 `ANALYTICAL` |
| ~~D-M5~~ **已修**（`bd21bc2`，随 D-S3 一并改）| `cp2k.rs` Langevin | `&THERMOSTAT TYPE LANGEVIN` 非法（合法值 NOSE/CSVR/GLE/AD_LANGEVIN）；CP2K 是 `ENSEMBLE LANGEVIN` + `MD/&LANGEVIN` | 按 CP2K 写法 |
| D-M6 | `cp2k.rs --smear` | 不写 `ADDED_MOS`，CP2K 报错（`qs_environment.F:2247`）；`--scf ot` 时 smear 静默丢。手册 `spin.md:76` 示例正中 | 写 `ADDED_MOS`；OT + smear 报错 |
| D-M7 | `cp2k.rs` PBC | 分子 / `--pbc z` 不写 `POISSON_SOLVER`，默认 PERIODIC 报错；分子 `&CELL` 无 ABC。`test_energy_molecular` 断言的正是这份不可运行的输入 | 非 3D 写 MT / WAVELET + ABC，或报错 |
| D-M8 | `cp2k.rs` | `--kpoints` + `--scf ot` CP2K 报错（`qs_scf_initialization.F:886-887`），builder 不拦 | 提前拒绝 |
| D-M9 | `ferro-workflow/src/qe.rs build` | 无胞仍写 `ibrav = 0` 不写 `CELL_PARAMETERS`，pw.x 拒收；`test_scf_basic` 断言的正是它 | 报错或要求给盒子 |
| D-M10 | `qe.rs` vc-md（推断，本机无 QE 源码） | 写 `ion_dynamics='verlet'`，INPUT_PW 中 vc-md 只支持 `'beeman'` | 改 `beeman` |
| D-M11 | `ferro-workflow/src/job_builder.rs` Gaussian | 周期结构丢胞、不写 `TV`，静默变团簇计算（ASE 写三行 TV） | 写 TV 或报错 |
| D-M12 | `cmd/job.rs` | 与所选软件无关的参数静默忽略（QE 不读 `--cutoff`/`--md-timestep`/`--thermostat`/`--pbc`）；`--task` 等枚举值读完输入、打印 auto-spin 之后才校验 | clap `ValueEnum`；不适用参数报错 |
| ~~D-M13~~ **已修**（`batch::check_suffix`，在 `CommonArgs::inputs`、chg-sdf、bader 三处入口调用。**未覆盖**：`dataset merge --suffix` 是目录后缀，合法值含 `.`，另一套语义，未动）| `batch.rs out_path`（及 bader、chg-sdf 的 `-s`） | `-s` 不校验字符：`-s 'a/../../escaped'` 写到 `-o` 之外 | 复用 `label_char_ok`，读文件前校验 |
| ~~D-M15~~ **已修**（参数检查移到建目录前；写侧格式判定抽成 `io_dispatch::out_format`，`write_trajectory` 改走它，convert 读输入前先调）| `cmd/convert.rs run` | `--start/--end` 校验前就建 `-o` 目录；写侧格式读完整条输入后才查 | 先查参数与格式再建目录 |
| ~~D-M16~~ **已修**（判据放进 `split_pair_args` 本身，有测试）| `main.rs split_pair_args` | 对所有子命令剥离 `--Xx-Yy=v`，`traj gr --P-O=2.3` 静默吞掉 | 只在 `argv[1]=="net"` 时剥离 |

**轻**

- A：`lammps_data.rs:135` 质量解析失败 `unwrap_or(1.0)` → 当成 H；`:208` 电荷 0 转 `None`；
  `lammps_dump.rs:131` 新版 general triclinic `BOX BOUNDS abc origin` 会按 lo/hi/tilt 误读（推断）；
  vasprun 空 `<i>` 不清 `reading_energy_name`（推断，最坏丢帧）；空 / 垃圾文件读成 0 帧返回 Ok，
  报错不点明「不是该格式」；extxyz writer 不写 `masses`、电荷 / 磁矩缺失补 0；`writers/qe.rs`
  无胞不写 `CELL_PARAMETERS`、写 `tot_magnetization` 不写 `nspin=2`（推断）；PDB 坐标 / 序号超宽错列；
  `extxyz.rs parse_properties` 的 `count` 解析失败 `unwrap_or(1)`
- B：`vanhove.rs:185`、`angle.rs:289` bin 数用 `ceil`，宽度不整除时末 bin 按整格归一、中心越界
  （`vanhove --r-max 1.0 --dr 0.07` 末 bin p_r 低估约 3.5 倍）；rotcorr `--legendre 1` 头部仍写
  `P2` / 「of c2」（`traj.rs:661`、`rotcorr.rs:329`）；`sq.rs:346 to_tables` 参数 `gr` 未用、`Result`
  永不 Err；rotcorr Sum 模式串行暴力（`:171-189`，推断）；msd lag 0 输出 −1.7e-18
- C：`cell.rs:147-160` 浮点 `rem_euclid` 可返回 L（`(-1e-17).rem_euclid(8.0)=8`），chg_sdf 三线性插值
  由此越界（推断）→ 下标再 `% n`；`spin.rs:62-69` 未知元素 Z=255 计入电子数（元素表只到 Rn，
  UO₂ 给双重态）；~~`dataset.rs:982,1003` `-f nan`、`--oo-min nan` 静默关判据、`--al6 nan` 删光
  后才失败；`diagnostics.rs:101-104` 全非有限时 `len()-1` 下溢~~（随 C-D4 修）；`diagnostics.rs``:41` 实际可采
  2·max−1 帧；三处注释错位 / 过时（`dataset.rs:281-293` 两段 `///` 挂错、`group_by_directory`
  两条注释矛盾、`network_type.rs:21` 仍写「digit = bridging ligands」）
- D：`--last-n 0` 不提前拦；`job` 不带 `-s` 退出码 0 不出文件；CP2K `r2scan + d3` 无 D3 参数
  （`qs_dispersion_utils.F`）；`--cp2k-basis` 拼错落进 `Custom` 写不存在的基组名；chg-sdf 帮助页
  `-s STEM` 与实际产物名不符，`map sdf` 的 `-s` 替换 stem 而 `density` 追加；chg-sdf 两条用户报错
  是中文；`box_builder` HashMap 序 + 无种子 RNG 不可复现；Gaussian 恒写 `%chk=job.chk`

**物理 / 口径（可能有意，待定）**

- bondlife 的 $S_C$ 是「比值的和」（Luzar–Chandler / gmx 同构），MDA `autocorrelation` 是「比值的
  平均」，lag 1 差 3e-4。手册 `bondlife.md:51-55,208` 应改说「与 MDA 只差平均顺序」
- S(q) 无窗函数，默认 `q_min=0.1` 低于 $2\pi/r_{max}\approx0.63$ Å⁻¹，NVT 上 total_xrd(0.1)=−0.42。
  建议 $q<2\pi/r_{max}$ 告警或提供 Lorch 调制。**2026-10-07 随 B-3 已加告警**（`[inputs]` 的 `q_trunc`）；剩两项待定：
  默认 `q_min` 改为逐输入 auto（$2\pi/r_{max}$，破坏性：默认产物行数变）、Lorch 窗选项（改峰形，新功能）
- `map velocity/force` 空体素写 0.0（`cube_density.rs:207`），与「测到 0」不可分；手册注明或另出计数网格
- chg_sdf 跨文件平均的是 ρ·V_cell，NPT 快照下是体积加权（推断）
- `cluster.rs:56 NetworkGraph.former_qn` 是桥氧个数，`AtomType::Former.qn` 是同元素连接数，
  core 里同名不同义 → 改名 `n_bo`
- Bader 归属原子的 MIC 用分数取整、不查上界，强三斜小胞取错镜像（推断，极少触发）
- merge `--mode by-source` 静默忽略 `--seed`，filter 遇无用 `--seed` 报错，口径不一；filter
  几何判据不读 `Frame.pbc`，手册 dataset 页未写
- `EV_TO_KCAL_MOL=23.0605419`（CODATA 2018 为 23.0605478，差 2.6e-7）
- VASP4 POSCAR 元素给 `X1`（ASE 会推）；xyz 首列是原子序数时 ferro 存 `"8"`（ASE 转 O）；
  QE `starting_magnetization` 存进 `magmom` 再经 extxyz 导出会被当 μB
- CP2K：TTF slab 映射成 `PERIODIC XYZ`；`--task force/md` 不请求 `STRESS_TENSOR`，自家 AIMD
  进 `collect` 无 virial；PBE0 / B3LYP 周期体系用 COULOMB 而非 TRUNCATED。QE bands/nscf 用
  `K_POINTS gamma`；QE MD 时间步固定 20 a.u. 不接 `--md-timestep`
- `merge` 侧向是插真空条带而非 ASE `stack` 的应变匹配；`box_builder` 按原子随机放、软约束
  （SiO₂ ×200 实测 241 对略小于 `min_dist`）

**代码布局建议**（判据注明；做与不做待定）

1. **格式分派下沉 ferro-io**：CLI `io_dispatch` 读 / 写 / `holds_multiple_frames` /
   `supported_formats` + Python 读 / 写共 6 张手写表，已在大小写（CLI 敏感、Python 转小写）、
   CONTCAR、`.lammps`、`.vasp` 上漂开。建 `enum Format` + `Format::detect(path)`（R6 第三处、
   同一变化原因；io 自己的职责，不违反分层）。A、D 两路独立提出
2. **单位常数改调 `units.rs`**：`cp2k.rs:18,357`、`qe.rs:7` 的 `BOHR=0.52917721`，
   `lammps_dump.rs` 的 `KCAL_TO_EV` / writer 的 `EV_TO_KCAL`，`box_builder::estimate_box_length`
   的 `N_A`（R6 复用已有；`CLAUDE.md`「转换走 units.rs」）。A、C、D 三路都提到
3. **md/util 合并**（R6，均 ≥3 处同一变化原因）：`--elements` 选原子 8 处（msd/vacf/vanhove/
   cube_jump 取下标、cube_density/cube_radius 计数与内联）→ `select_atoms`；格点视图解包裹三份
   （`vanhove.rs:105`、`cube_jump.rs:84`、`util.rs:63`）→ 下沉 `msd::unwrap_tor`；cube spacing
   矩阵三份（`cube_density.rs:222`、`cube_radius.rs:217`、`cube_jump.rs:208`），`voxel_idx` 顺带
4. **搬移**：`cube_reference_frame`（在 cube_density，radius/jump 共用）与 `CellList`（在 angle，
   gr 共用）下沉 `md/util.rs`；`spin.rs` 搬 ferro-workflow、`build_network_graph`/`NetworkGraph`/
   `LigandKind` 搬 ferro-analysis（下沉判据反向：只有一个中间层叫它的名字）
5. **Bader weight 复用 grid**：`bader_weight.rs:236-262` 抄了 `assign_chg2atom` 与
   `calc_atomic_vol`，改传 `&volchg[1..=nvols]` 复用（R6），同时消掉 C-D2 那一族索引错误
6. **同文件已有函数直接复用**（R6）：`network_type.rs` 的 `label()`/`class_rank()` 改调
   `site_digit()`；`Trajectory::select` 改调 `select_indices`；`find_clusters` 改用
   `classify_frame_detailed` 的 `ligand_formers`（实测同为 569 团簇、56→19 ms；第 3 步注释
   「恰好两个」与代码 ≥2 相反）；sq 的 `build_xrd_weights_at_q`/`build_neutron_weights` 分母
   移进 `weights_from_factors` 后内联，`elem_z_local` 内联（R3）
7. **内联 / 删死代码**：`ferro-structure/src/typing.rs classify_trajectory` 内联进 `net.rs`
   （R1 单调用 + R3；跳过无胞帧造成下标错位）；`args/cube.rs CubeCliMode`（从未被 clap 解析，
   三个变体从未构造）；`batch::Output::join_str` 零调用
8. **拆分 / 改名**：`cmd/dataset.rs`（1921 行）拆 `cmd/dataset/{mod,collect,filter,merge,split}.rs`，
   `cmd/inspect.rs` 搬入（R2a）；`ferro-workflow/src/job_builder.rs` → `gaussian.rs`，`lib.rs`
   改显式导出；`job.rs` 约 10 处字符串 match 改 clap `ValueEnum`（同 D-M12）；`ml/geometry.rs`
   的生产代码夹在两个测试模块之间，移到前面
9. **性能**：`Cell::minimum_image` 每次调用 `try_inverse`，用在 network_type 三处、`cluster.rs`、
   Bader 归属的 O(N²) 循环里 → core 加预存逆的 MIC 辅助（与 `issues.md:658` 否决的「并入现有
   方法」不是一件事）；`classify_formers` 每个形成子重算 `params.ligands()`（`:385`）
10. **决定不合、需两处互指注释**（R6）：`cell_to_lammps`/`lammps_cell_matrix`/`bounding_box`
    在 lammps_data 与 lammps_dump writer 逐字两份（**互指注释已随 `a34b488` 加上**）；POSCAR 头部
    解析在 `vasp.rs` 与 `chgcar.rs` 两份且已漂过一次（`chgcar.rs:71` 仍是「`d` 开头才算 Direct」
    旧规则，与 A7 一起修）；msd 与 bondlife 的逐帧 (Mᵀ, Mᵀ⁻¹) 预计算；`filter_split`/`merge_split`；
    `job.rs:204` 与 `convert.rs:155` 的 `-o` 尾分隔判断；`vacuum.rs`/`merge.rs` 轴名解析；
    rayon `build_global` 两处；io 与 workflow 两份 QE writer（分层所限不能合，先修平 nspin 漂移）。
    ferro-structure 4 处手写完整 `Frame { … }` 字面量**有意不合**：新增字段时编译器逼每个构造点表态
11. **库 crate 错误类型**：`Table::validate`、`concat_union` 返回 `Result<_, String>`（轻）
12. **注释语言**：`cp2k.rs` 12 处、`qe.rs` 9 处、`box_builder.rs` 6 处英文 `//`；断言消息
    `batch.rs:512`、`supercell.rs:322`、`box_builder.rs:531`、`dataset.rs:1854` 为英文
13. **过时注释**：`map.rs drive` 文档仍提 `--outdir`；`traj.rs` 模块文档「seven analyses … PNG」；
    两份 `qe.rs` 写 `fe-job`/`fe-convert`；`net.rs fmt_means` 的文档挂到 `warn_edge_sharing`；
    `cube_density.rs:12-16,41-44` 写 `--grid`/`--mode`；`util.rs:1` 自称「三处以上共用」而
    `build_avg_frame` 只有 1 个调用者（有测试，R2c 保留）；`dataset.rs` 模块文档「filter/merge
    not implemented yet」

**待核实**（第 8 轮顺带记下，未查）：PBE0 / B3LYP 在周期体系下用全程 COULOMB 的 HFX，CP2K
通常要求 `TRUNCATED` + `CUTOFF_RADIUS`，ferro 不写。

**测试缺口**：workflow 只有子串断言、无一条对照 CP2K / pw.x 规则（S1–S3、M4–M10 全漏）；
`help_sync` 的 PAGES 不含 job 三页（约 35 个参数无防漂）；CLI 集成测试只有 bader 与 npt；
`check.rs` 无 inf；sq q=0；奇异胞进 vanhove / map sdf；gr 三斜数值（可用 27 镜像暴力对拍）；
cube_sdf ≥8 P 的运行时上界与组内置换语义；vanhove / angle 末 bin；rotcorr `--legendre 1` 头部；
bondlife 对 MDA 定义钉死；spin 镧系 / 锕系 / 未知元素；weight BCF 与 ACF 一致性、ACF 用 ASE /
pymatgen 回读；filter 的 NaN 帧、多 system 报告头、零平台曲线；`wrap_position` 恰返 L；三斜下
`classify_frame`；merge 乱序输入对 dpdata 端到端；extxyz 未闭合引号 / 裸键；CIF 仅 H-M、
`data_global`、非周期自写自读；PDB 逐帧 CRYST1；LAMMPS 无胞写出、正交不写三斜、slab pbc 往返、
空 data；VASP 6 哈希元素行、CHGCAR `Fractional`；mixed type 应被拒；CP2K `&COORDINATION` 在
`&COORD` 前（regtest H2O-meta.inp）；QE 空白分隔 namelist；分派大小写、Python 写 CONTCAR；
merge 窄 A / 倾角不同 / 垂直间隙；`-s` 非法字符、`--last-n 0`、`--P-O=nan`；Python 绑定零测试。

### 2026-10-01 全库审查发现（19 条，全部复核属实；表内 19 条均已修，表后「待定口径」仍开着）

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

M1 期间另见、未修（待定口径）：

- **`.lammps` 扩展名两处口径相反**：CLI `io_dispatch` 当 LAMMPS data，Python
  `ferro.read` 当 LAMMPS dump（`ferro-python/src/io.rs` 的 `"lammpstrj" | "dump" | "lammps"`）。
  同一个文件两个入口读出两种东西
- **dump 的 `typelabel` 与 `ix iy iz` 未接**（M2 时按 LAMMPS 手册逐列过了一遍）：`typelabel`
  要先定它与 `element` 谁优先；`ix iy iz` 对 msd 无影响（msd 自己按最小像展开）。其余列
  （`mol`、`mu*`、`c_*`/`f_*`/`v_*` 等）Ferro 数据模型无对应字段，忽略是对的
- **cube_*（除 `calc_cube_radius`，已随审查低 6 改 `Result`）、`calc_chg_sdf`、`calc_cluster_sdf` 仍返回 `Option`**（M3 未动；`calc_network` 已随 L2 改 `Result`）：
  参数已在读文件前查过，剩下的 `None` 都是轨迹层面的（缺 cell、缺速度、没有团簇），
  但 CLI 只能猜原因（`missing cell, velocities, or forces?`）。改 `Result` 是同一个做法
- **`Atoms` 段写在 `Masses` 之前**时 Masses 不会被读，元素退化为 `X<type>`
  （`lammps_data.rs` 按 Masses → Atoms 的顺序单遍扫描）。LAMMPS 不规定段序
- **`ferro-structure::find_clusters` 未查最小镜像上界**（L2 时见到）：它在 structure 层、
  返回 `Option`、以 `Frame` 为单位，`check::within_minimum_image` 在 analysis 层够不着。
  目前只有库 API 与测试调用，无 CLI 入口
- **extxyz 的 `momenta` 列被当成速度**（L6 时见到）：`readers/extxyz.rs` 找不到
  `velocities` 就读 `momenta` 存进 `velocities`。ASE 的 momenta 是 $m\,v$（amu·Å/fs 量纲
  按 ASE 内部单位），直接当速度会让 vacf 等差一个质量因子和单位换算。要么按质量除回来，
  要么不认这一列

**轻**

| # | 位置 | 问题 | 备注 |
|---|---|---|---|
| ~~L1~~ **已修**（`read_cube` 保持绝对坐标；减原点只留在 Bader 用的 `read_cube_as_chg`；加非零原点两次往返测试）| `readers/cube.rs:87-95` vs `writers/cube.rs:44-50` | reader 把原子坐标减去 origin，writer 当绝对坐标写（origin 照写）→ 非零原点的 cube 每往返一次原子平移 −O | 读码确认。CLI/Python 均不调 `read_cube`，只影响库 API；现有往返测试是零原点 |
| ~~L2~~ **已修**（`check::within_minimum_image` 逐帧取最紧上界，超出报错并点名帧；angle / bondlife / rotcorr / network 入口调用，filter 原有的同一段改为调它。`calc_network` 顺带改返回 `Result`。gr 照旧截断 `r_max`）| angle（`CellList` 每轴 `.max(1)`）、`network_type.rs:352,393,452`、`bondlife`、`rotcorr` | 只有 gr 与 dataset filter 检查最小镜像上界；其余截断超过 MIC 上界时只取单一镜像，静默漏邻居 | 逻辑确认；只在小胞 + 大截断时触发 |
| ~~L3~~ **已修**（`chacha20 0.10.1` 设为直接依赖，用户批准；`shuffle_order(20, 666/42)` 换前后逐位相同，结果硬编码进测试）| `ml/merge.rs:129` | `StdRng` 不可移植：rand 0.10.2 源码 `rngs/std.rs` 明写「any future library version may replace the algorithm」，与「seed 默认 666、可归档复现」冲突 | 改用 `chacha20::ChaCha12Rng`（当前 StdRng 的实现，已在树里）固定版本，同 seed 结果不变 |
| ~~L4~~ **已修**（续行改为单个 `\`，测试断言消息无反斜杠与连续空格）| `cmd/dataset.rs:191` | 普通字符串里写 `\\` + 换行，错误消息里多一个字面反斜杠和一段缩进 | 同 `issues.md`「用 Python heredoc 改 Rust 源码时」那一族 |
| ~~L5~~ **已修**（判断移到 `FilterParams::enabled()`，`filter_frames` 与 `FilterResult::enabled()` 都调它）| `ml/filter.rs:270` | `filter_frames` 内联了一份与 `FilterResult::enabled()`（`:306`）相同的 match | match 是穷尽的，加变体会编译报错，**漏改不会静默**；风险是两处条件改得不一致。按 R6「同文件已有函数直接复用」 |
| ~~L6~~ **已修**（extxyz 的力、速度、坐标与 `charges`/`masses`/`magmoms`、注释行的 `energy`/`temperature` 解析失败一律报错，点名帧号、行号、字段 —— 后两组按用户同意顺带修；cube 体数据坏值报错并点名行号）| ~~`cp2k.rs:104-106`~~（随 S5 修）、~~`qe.rs:120-122`~~（随 S4 修）、`extxyz.rs:121-129`（力、速度）、`readers/cube.rs:105`（体数据） | `parse().unwrap_or(0.0)`：解析失败冒充「测到了 0」。QE 的 Fortran 写法 `0.25d0` 即中 | 读码确认；改为报错并点名行号 |
| ~~L7~~ **已修**（`a2c013a`，删分支）| `readers/cp2k.rs:124-137` | 把 `&COORD` 第 5–7 列当速度（注释写「restart 里」）。CP2K 的 `&COORD` 第 5 列是分子名，速度在独立的 `&VELOCITY` 段（bohr/au_time，不是 bohr/fs），reader 不读该段 | 读码确认；实际几乎不触发（第 5 列是字符串则解析失败跳过），触发则单位错且速度数组长度可能与原子数不一致。删分支，或正式读 `&VELOCITY` |

**测试缺口**（与上表对应）：dump writer 无「非下三角胞 + 力/速度」；POSCAR 无「元素交错 +
速度」；QE / CP2K 只测了 `{}` 与 `.TRUE.`；lammps data 无「无注释 Atoms」；无截断 dump
fixture；cube 往返只测零原点；cube_radius 无窗口超网格的用例。

**审查未覆盖**：Bader 全套（已登记项之外）、`chg_sdf`、`cube_sdf` 聚类、`cube_jump`、
`spin.rs`、`cp2k_basis_db`、workflow 模板、`cp2k_sp`、pdb / cif writer、doc 渲染器、
ferro-python 运行时。（`--metal-units` 全链路已随 2026-10-03 的 `--units` 改造走过一遍）

### 全仓规则合规扫描的遗留违规（2026-09-27 扫描）

对照 `CLAUDE.md` 逐条扫的结果。**已查无违规**：分层依赖、`clippy` 零警告、
`as_slice()`、`Molecule`、`label()` 反解析、列并集补 NaN、手册页三处登记、
`ferro-python` 版本同步。以下按严重程度排，每条独立可做：

1. ~~`ferro-io` 与 `ferro-workflow` 用 `anyhow`~~ —— 2026-09-28 降为中优先级，
   见「优先级中」的「库 crate 的错误类型」
2. **`ferro map` 的产物名随输入数变**：`cmd/map.rs` 的 `stem_for` 与
   `let multi = inputs.len() > 1` —— 单输入 `density.cube`，多输入
   `density_<stem>.cube`。正是归档里「`-i` 恒为 `Vec`，单一代码路径」点名的
   「两套命名」。修法待定（恒带 stem / 只看 `-s`），属破坏性改动，要同步手册与
   changelog
3. ~~帮助页偏离五段模板~~ —— **2026-09-28 已修**（`d4fc999`）。net 的 Output
   11 → 3 行、Qn 口径说明删去（`network.md` 已有）；job 三页的分组降为
   `Parameters:` 下的缩进小标题；`convert` 的格式矩阵是 `-i`/`-o` 的值域，
   并入 `Parameters:` 而不单列一段（用户认可）；info 补 `Output:`；bader 的
   「为什么写在输入旁」移交手册。仍超 30 行的（`job -s cp2k` 46、`convert` 55、
   gr sq angle rotcorr、`map sdf`、`dataset filter`/`merge`）全是参数表，按规则不砍
4. ~~中文 `///` / `//!` 约 550 行~~ —— **2026-09-28 用户裁定允许中文 doc 注释**，
   `CLAUDE.md` 规则已同步改写，不再算违规，别再重扫
5. ~~手册化学式用 Unicode 下标~~ —— **2026-09-28 已修**（`9145e7a`）：rotcorr.md
   三处、changelog.md 一处，外加 changelog 的裸文本 `Q^n_m`。`docs/src` 剩下的
   唯一一处 `α₁` 在行内代码里，按规定不动
6. **`ferro-analysis/src/trajectory_analysis.rs`**（提示）：自称「旧接口，保持编译
   兼容」，经 `lib.rs` `pub use *` 全量导出，仓内零调用。删或留待定，可并入中优先级
   的「零调用清单」一起判

未做：全局规则「何时该有一个独立函数」（R1~R6）的逐函数审计。

### DeePMD mixed type 数据的读写（2026-08-26 提出）

**现状（2026-10-07，审查 A6）**：reader 见到任一 set 有 `real_atom_types.npy` 即报错拒收。
**只拒不读**是有意的：只做读侧的话 `filter` 读进来再写回标准布局（按第一帧定 `type.raw`），
mixed 数据会静默降级成坏数据、退出码 0 —— 比拒收更糟。要支持就读写与下列判据一起做。

现在 `readers/deepmd.rs` / `writers/deepmd.rs` 只做**标准 system**：`type.raw` 一份
定型，故一个 system 里所有帧的成分必须完全相同 —— 这正是 `collect`「同目录成分不符
报错」与 `merge`「按 `composition_key` 分组」两条现有约束的来源。

mixed type 布局（**先核对 DeePMD-kit 文档与 dpdata 的 `deepmd/npy/mixed` 再动手**，
以下是待验证的理解）：

| 文件 | 标准 system | mixed type |
|---|---|---|
| `type_map.raw` | 该 system 的元素表 | 全局并集 |
| `type.raw` | 逐原子真实类型 | **全 0 占位** |
| `set.NNN/real_atom_types.npy` | 无 | `(nframes, natoms)` 的整型，逐帧逐原子给真实类型 |

要点与判据：

- **价值是装下成分不同的帧**（原子数仍须相同 —— dpdata 的 mixed 也按 natoms 分
  system）。所以这件事做完，`collect` 的「成分不符报错」和 `merge` 的分组要重新定：
  是继续分组、还是给一个 `--mixed` 让它们合成一个 system。**默认不改**，因为
  mixed type 只有 DPA 系列 / 多任务训练吃得下，普通 DeePMD 训练不认
- 读侧先做：能读回 dpdata 产出的 mixed system，`real_atom_types` 经 `type_map` 映射
  回元素符号，落进 `Trajectory` 天然装得下（帧与帧的元素本来就各存各的）
- 写侧作为开关，不改默认布局。`type_map` 的并集需要**可复现的稳定序** ——
  `merge.rs` 的 `(Z, 符号)` 规范序已满足，沿用同一条（注意与 dpdata 的字母序不同，
  这条差异 `progress.md`「已知限制」已记）
- npy 是整型：现有读写走的都是 `float64`（磁盘上一律二维 f64），
  `real_atom_types.npy` 是 int，**dtype 分支是新的**，读侧要同时收 int32/int64
- **`sort_atoms` 不能无脑套用**（2026-09-22 查证）：DeePMD-kit 在 mixed type 下
  明确**不**按类型排序（`deepmd/utils/data.py` 的 `sort_atoms` 文档写着
  "except mixed types"），因为 `real_atom_types` 是逐帧各异的。而 `collect`
  现在无条件排 —— 两者相遇时要先定 mixed system 的规范序是什么
- 与 `filter` 的关系：`filter` 现在按 system 读回再写出，mixed system 经它一趟必须
  仍是 mixed，否则静默降级成「全 0 类型」的坏数据

---

### ferro job 从轨迹抽帧（2026-08-22 提出，告警已于 2026-08-24 补上）

`job` 现在无条件取第 0 帧（`cmd/job.rs:186`）。多帧输入会告警（见下），但**仍然只
生成第 0 帧的输入**——喂一条 500 帧轨迹得到的还是最没平衡那个构型，只是这次用户
知道了。真正的抽帧 + 批量生成仍待做。

`convert` 的 `--start/--end/--stride/--number` 已把选帧逻辑做进
`Trajectory::select_indices` / `spread_indices`，job 复用即可，不必重写。
真正要定的是**产物命名与批量语义**，与 `convert` 不完全一样：

- job 的 `-o` 现在是完整路径且有默认值（`job.gjf` / `job.inp`），多帧要变成
  `job_0000.inp` 这类。`-o` 的语义已于 2026-09-20 统一（见归档），多帧命名仍待定：
  沿用 `convert` 的 `indexed_path`（序号插在扩展名之前）即可，不必另起一套
- 一个构型一个输入文件是显然的（QC 输入本来就一个结构一个），所以不存在
  `convert` 那种「一个文件还是 N 个」的分支 —— 恒为 N 个
- ~~**最小可用的第一步是加警告**~~ **已做（2026-08-24）**：`multi_frame_warning`
  打三行 —— 帧数、忽略了几帧、以及可直接抄的 `ferro convert --number` 命令。
  产物形态未动。**做这一步时撞出 `ferro job` 在 debug 构建下必 panic**（clap 的
  `help` 参数重名），见 `issues.md`「构建 / 工具链陷阱」

在此之前的变通：`ferro convert -i traj.dump -o conf.vasp --number 20` 抽成单帧
文件，再逐个喂给 job。

### ferro dataset：剩余小项（2026-08-25）

`collect` / `filter` / `merge` 均已落地。剩下的都不大：

- **额外键搬运**：`atom_ener` / `fparam` 这类 `Frame` 装不下的项，读时告警、
  写时丢失。真出现时再设计（需要一条绕过 `Trajectory` 的按帧索引搬运通道）
- ~~**GPUMD/NEP 的 `train.xyz` 导出**~~ **已做（2026-08-26）**：走
  `filter`/`merge` 的 `--type nep`，没有新命令。见归档

两件**不必新写**的事已经在库里：帧区间与间隔用 `Trajectory::select_indices` /
`spread_indices`（`convert` 的 `--start/--end/--stride/--number` 就是它）；
O-O 间距与 Al6 配位用 `ferro_core::classify_frame` 出的
`AtomType::Former{cn,..}`（`Al_4`/`Al_6` 的数字就是配位数）。`filter` 在 CLI 层
组合 io 与 analysis，中间层仍不互依赖。

还需要定的：`filter` 要读回 npy（`ndarray-npy` 读侧未用过）；merge 的 `type_map`
跨 system 对齐（`collect` 已按 `(Z, 符号)` 排序，这条规则可复现是对齐的前提）。

### ferro map chg-sdf 的 --cubes 拆成单文件（2026-08-11 提为高）

现在是多个 cube 聚合成**一张** SDF（`cmd/map.rs::run_chg_sdf`），与 `ferro map` 其余
模式「一个输入 → 一个 `.cube`」的语义相反，`--cubes` 也是 `-i` 之外唯一的输入参数。

阻塞点不在遍历而在**中间产物格式**：跨文件加权平均不可交换，逐文件出产物之后若还想
聚合，产物里必须带样本计数（否则 5 帧的 SDF 与 500 帧的 SDF 被等权平均）。
故顺序是：先定义带计数的中间格式 → 再拆 `--cubes` → 再接批处理遍历。

---

## 优先级中

### 2026-10-02 物理审查的三条严重问题（用户定为中，功能后续会改写，届时一并修）

三条都已用数值实验证实（临时 crate，不在仓库），均未修。

**S-A　Kabsch 旋转求成了转置**（`md/cube_sdf.rs::kabsch_rotation`）。H = Σ m·rᵀ =
U S Vᵀ 时，使 R·m ≈ r 的解是 R = V·D·Uᵀ，代码返回 `u * diag * vt`（= Rᵀ）。
`map sdf` 与 `map chg-sdf`（`dft/chg_sdf.rs` 复用，`rotate_grid` 的 pull 插值同错）
叠加的是反向旋转的团簇，置换枚举也按错误 RMSD 选。实测同一畸变 PO₄ 刚性旋转两次，
RMSD mean 1.47 Å（应 ≈ 0）。现有 `test_kabsch_known_rotation` 点全在 xy 平面、绕 z
转，H 退化且旋转与 diag(1,1,0) 可交换，两种写法结果相同，所以没抓到。修法：
`v * diag * u.transpose()` + 三维非退化回归测试（随机旋转 + 四面体，RMSD < 1e-10）。

**S-B　Bader 三条网格路的真空电荷恒为 0，真空里冒出大量伪体积**（`dft/bader_grid.rs`
on-grid `:213-232`、near-grid `:580-597`、off-grid `:763-770`）。`volchg.push(0.0)`
后立刻把这个 0 当 `vacchg` 读出；且真空在梯度上升**之后**才标，平坦/噪声真空区的
每个局部极大都成了独立 Bader 体积再分给最近原子（与 `dev/bader.md` §4 的伪代码顺序
相反）。实测 20 Å 盒两高斯峰 + 5e-4 e/Å³ 背景：网格总电荷 17.61 e，三条路 Total 都是
13.65 e，nvols 5 万量级；`vacvol` 是对的。修法：上升前 `mark_vacuum`，路径跳过
`volnum == -1`；`vacchg = Σ_{vac} rho / nrho`；补带真空区的 fixture。

**S-C　Bader weight 不守恒电荷，权重也不是 Yu–Trinkle 的**（`dft/bader_weight.rs`
`:64`、`:138-140`、`:162`）。实测两峰无真空，网格总量 15.925 e：near-grid 严格守恒并随
网格收敛，weight 丢 0.1–0.3 e 且不收敛。机理（推断）：边界点用负数存 basin，若某点所有
上游邻居都是同一个边界点，它被判为内部点并继承负 basin，不在任何 `neigh` 列表里，电荷
整体丢失。另：通量系数应为 Voronoi 面积 / |R|，代码用 |R|；WS 筛选严格不等式把面积为 0
的棱/角向量也留下了；`ionvol` 不计边界点；真空点参与流动。修法：真空标记与 basin 号
分开存、按 Voronoi 面积算 α、先标真空再流动。与下面「weight 的真空电荷取错」一并处理。

### 库 crate 的错误类型：`ferro-io` / `ferro-workflow` 用 `anyhow`（2026-09-27 扫描，2026-09-28 由高降中）

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
- **`ferro-workflow` 搁置**：三个文件只有 `use anyhow::Result`，唯一错误源是
  `writeln!` 写进 `String` 的 `fmt::Error`（实际不会发生），迁不迁行为都不变。
  若要迁：`ChemError` 加 `From<std::fmt::Error>`、删依赖，`cmd/job.rs` 不用改
- **未做**：`CLAUDE.md` 那条规则仍是「一律 `ChemError`」，尚未写明 io 例外

### `ferro job` 的自旋多重度：对照 CP2K 源码后的四处待定（2026-09-30）

**CP2K 自己不推断多重度**（`examples/cp2k-2026.2`）。`qs_environment.F:2013` 在
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

**待定的四条**（2026-09-30 grilling 提出，用户尚未裁决）：

1. **CP2K / QE 默认开 auto-spin**（`job.rs:292`，手册 `cli-reference.md:301`），与
   CP2K 自己「只按奇偶给 1 或 2」相反。选项：保留 / 改保守默认（破坏性）/ 保留但 M > 2
   时醒目提示
2. **显式 `--multiplicity` 不校验奇偶**（`job.rs:254`）：矛盾值照写，CP2K 启动时才报错。
   倾向于写文件前报错，放三个 builder 共用处（Gaussian / QE 同样要求这个奇偶条件）
3. **`--smear` 下 M 只是初值**：`FIXED_MAGNETIC_MOMENT` 默认 −100，
   `qs_mo_occupation.F:288` 走 `set_mo_occupation_3`，α/β 共用一个费米能级，总磁矩在
   SCF 里自由变化。ferro 的 `&SMEAR` 不写该关键字，手册 `spin.md:76` 的
   `--auto-spin --smear` 示例因此名不副实。选项：只改手册 / 自动锁 M−1 / 加开关
4. **多磁性中心按铁磁叠加**：`ion_unpaired` 逐离子取高自旋后直接求和，是**铁磁上限**，
   现有 warning 只说了「单离子高自旋」没说离子间排布

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

### 零调用清单（2026-09-12 实测快照，判断待定）

这不是待删清单。**多数功能还没有进入实际使用、没有调试过**，因此没法区分
「需求确定但 CLI 层还没接」与「写代码过程中的遗留」——`cube_jump` 就是前者的
现成例子（429 行实现 + 10 个测试 + 手册页都在位，缺的只是 CLI 入口）。

留这份快照是为了**将来能做 diff**：等这批功能跑起来之后重扫一次，「用起来了」
与「始终没人碰」自动分开。现在删掉任何一条，都是拿判断力换行数。

按 2026-09-12 的实测（727 个生产函数，扫描排除 `#[cfg(test)]`）：

| 类别 | 位置 | 行 |
|---|---|---|
| 整文件零引用 | `ferro-analysis/src/trajectory_analysis.rs`（文件头自称「旧接口，保持编译兼容」；4 个 `pub fn` 里 2 个连测试都没有，2 个只有自己的测试） | 113 |
| 整文件零引用 | `ferro-analysis/src/geometry.rs`（`distance`/`angle`/`dihedral`/`radius_of_gyration`/`bounding_box`，3 个测试；`docs/src/` 零提及） | 104 |
| 整文件零引用 | `ferro-workflow/src/templates.rs`（其中 `orca_sp_template` 是 Ferro 不支持的 ORCA） | 25 |
| 整文件零引用 | `ferro-cli/src/args/corr.rs`（`CorrMode`）· `args/traj.rs:5`（`TrajMode`），0.2.0 前 `--mode` 时代遗留 | 8 + 8 |
| 单位转换链 | `ferro-core/src/units.rs` 的 `convert_length` / `convert_energy` / `convert_time` **生产零调用**（只有自己的测试），连带 6 个私有方法只服务它们。四个 `convert_*` 里只有 `convert_pressure` 活着（vasprun / vasp_outcar / cp2k_out / filter / dataset 共 5 处） | ~55 |
| 被绕过的算法 | `ferro-analysis/src/dft/bader_grid.rs:375 max_neargrid` —— 见下一条 | 20 |
| 取代方案只落实了一半 | `ferro-core/src/network_type.rs:233 display_rank` 生产零调用（`class_rank` 有 3 处）。`progress.md` 说两者一起取代了三个 `*_label_order`，实际只落实了一半 | 8 |
| 基础访问器 | `Frame::geometric_center` / `atom_mut` / `wrap_all` / `is_periodic`、`Trajectory::frame_mut` / `iter_frames` / `time_at`、`Atom::distance_to`、`Table::n_cols`、`compounds::with_density`、`ClusterResult::ids`、`ml/diagnostics.rs:177,182`、`filter::is_clean`、`qn_elements::has_qn`、`charge_grid::lat_dist_i` | ~50 |
| ~~未使用依赖~~ | ~~`ferro-workflow` 的 `serde`~~ —— 2026-09-26 随 serde 整体移除 | — |

**方法上要记住的一条**：`cargo clippy` 零警告不说明没有死代码 —— 库 crate 里
`pub` 项不触发 `dead_code` lint，上面全部躲过了它。查零调用要按名字逐个 grep
全仓，不能靠编译器。

~~`Atom`/`Frame`/`Cell`/`Trajectory` 上的 `derive(Serialize, Deserialize)`~~ ——
**2026-09-26 已删**：全仓没有任何格式后端（serde_json / bincode 都没有），derive
从未被调用。连带摘掉三处 nalgebra 的 `serde-serialize`，`Cargo.lock` 里 serde 系
归零。将来要做检查点或 Python pickle 时，四个 derive + 两行 Cargo.toml 即可加回。

### bader weight 的真空电荷取错（2026-09-20 发现，处置未定）

`bader_weight.rs:265` 用 `volchg[nvols]` 当真空电荷，而该数组在 weight 方法里是
**1 索引**的（上方几行的注释自己写着），取到的是最后一个 Bader 体积。
`bader_grid.rs` 的三条路是 0 索引、真空在 `[nvols]` —— 下标对，但那里是硬编码的 0（见上面 S-B）。

实测（`tests/CHGCAR_2atoms`）：ACF 末尾的 `Total` 报 158.98 e，网格实际只有
105.99 e。逐原子的 `ionchg` 不受影响，被污染的只有 `vacchg` 与由它加出来的
`Total`，以及 ACF 的真空行。

**改一个下标就够，但会改变 weight 方法的产物**，所以没有顺手混进 `-o` 那一轮：
bug 修复的 diff 要能一眼说清「变的都该变」（判据同 `build_avg_frame` 那条）。
修的时候连带做三件事：

- `ferro-cli/tests/bader_reports.rs` 的真空断言把 `Weight` 加回来（现在注释里
  指名跳过它）
- 真空非空的情形现有 fixture 盖不到 —— `CHGCAR_2atoms` 的密度处处高于阈值，
  三条路的 `vacchg` 都是 0，所以这个 bug 只在 `Total` 行上看得出来。要钉住修复
  得再加一个带真空区的 fixture，或给现有的调 `--vacval`
- `dev/bader.md` 未提真空桶的索引基准，修完补一句

### `max_neargrid` 被绕过（2026-09-12 发现）

三件事实，处置未定：

- `bader_neargrid:542` 直接调 `step_neargrid` 把上升循环内联了，因而绕过了
  `max_neargrid`（`bader_grid.rs:375`，20 行）
- 两个兄弟都在用：`max_ongrid` ← `refine_edge:464`，`max_offgrid` ← `bader_offgrid:744`。
  三条路里只有 near-grid 这条不对称
- `dev/bader.md` 有 **3 处**写着 `max_neargrid` 的规格（含 known=1 / known=2
  两种退出条件那条陷阱），删函数要跟着改

**不预先定处置方向**：Bader 是本仓唯一「对着 Fortran 逆向出来」的模块，删还是
让 `bader_neargrid` 用回它以恢复三条路对称，该在真要动这个算法的时候连着算法
一起判断。现在定一个方向，只会让将来的人跳过这次判断。

### 带速度的测试 fixture（2026-09-12 提出，2026-09-21 仍未做）

`tests/` 两条 fixture **没有速度数据**，于是 `vacf` / `vanhove` / `rotcorr`
跑不到写文件那一步 —— 这三个命令端到端从来没验过（0.2.0 改名那轮就只验了
gr/angle/msd/sq/net/map），2026-09-12 的 golden master 基线也盖不到它们。

补一条带 `vx vy vz` 列的 LAMMPS dump（5 帧即可）就能把三个命令纳入。
注意**不能**用全零速度顶替：VACF 在全零速度下是退化情形，跑得通但看不出
实现对不对，那是假的安全网。

### `ferro map jump`：丢失的 CLI 入口（2026-09-03 审计发现）

`ferro-analysis/src/md/cube_jump.rs` **429 行实现完整、10 个单元测试、已从
`md/mod.rs` 导出**，但没有任何 CLI 入口，`lib.rs` 的 `pub use md::{...}` 清单里
也漏了它 —— 0.2.0 把八个 `fe-*` 合并成单个 `ferro` 时，`fe-cube -m jump` 没有
跟着搬过来。

**这不是死代码，是漏登记的待办**：`docs/src/analysis/cube-jump.md` 开头就写着
「仅库函数，没有 CLI 入口……待 `map jump` 补上」，`SUMMARY.md` 也挂着这一页，
即手册已经对用户承诺了。此前 `plan.md` 里却没有这一条。

补的时候按 `map` 现有形状走：`cmd/map.rs:223 run_grid` 的 `drive()` 已经把
「一输入一 `.cube`」的遍历收好了，`calc_cube_jump` 的签名与
`calc_cube_density` / `calc_cube_radius` 同形，接一个分支 + 一个 `JumpCmd`
参数结构体即可。手册页与 `doc.rs` 的 `PAGES` **都已经在位**（`ferro doc cube-jump`
现在就能打出来），缺的只有 `help.rs` 的帮助页与 `print_map_overview` 里的一行。

### VASP：ML_FF 的 OUTCAR 与 io_dispatch 注册（2026-08-27 提出）

两件被这一轮明确划在范围外的事：

- **ML_FF 的 OUTCAR**：VASP 的机器学习力场跑出来的输出用 `free  energy ML TOTEN`
  与 `ML FORCE`。dpdata 用一个 `ml=True` 开关切 token，但它给两者的行偏移是
  `[14, 4]` —— **块结构本身就不同**，不只是换个名字。没有样例就没法验，等有输出
  再补；`vasp_outcar.rs` 的锚点已经是多候选表的形状，加 token 不必重构
- **`io_dispatch` 注册**（2026-09-21 复核：**仍未做，本轮有意绕开**）：把 LAMMPS
  导出挂在 `dataset collect --type inspect` 上而不是 `convert` 上，于是下面这条
  判据一次也没被逼着定。`ferro convert -i OUTCAR -o traj.xyz` 依旧不认识。要动
  `io_dispatch.rs` 的两处 match、`supported_formats()` 那张有测试盯着的清单，以及
  `ferro-python/src/io.rs` 那处独立分派。判据要先定：`OUTCAR` 没有扩展名，而
  `io_dispatch` 现在是按扩展名（加 `POSCAR`/`CONTCAR` 的前缀特例）分派的 ——
  是给它加内容嗅探，还是只按前缀认 `OUTCAR*`

### CP2K 的 EXTXYZ 把 atom kind 写进 species 列（2026-08-26 提出）

CP2K 新版的 `MOTION/PRINT/TRAJECTORY` 多了 `FORMAT EXTXYZ`，而它的
`PRINT_ATOM_KIND` 会**把 subsys 里的 atom kind 写进 species 列**（文档原文：只对
XMOL 与 EXTXYZ 有效）。Ferro 的 extxyz reader 假设 species 是纯元素、位点标签走
独立的 `label:S:1` 列，撞上这种文件会把 kind 当元素收下。

与 LAMMPS dump「没地方放第二个名字只能折进 element 列」是同一族问题，但**不能照抄
那边的解法**：dump 那边是无条件按下划线拆，而 extxyz 的 species 列在合规文件里就
该是纯元素，无条件拆会误伤。要定的是判据 —— 拆还是不拆、按什么拆、
`split_element_label` 的 `Unknown` 分支怎么处理（`Pb` 那类贪婪前缀误判的教训见
`issues.md`）。

CP2K 的 EXTXYZ **只写 cell + 坐标**，不写 stress/virial（力在 `PRINT/FORCES` 另一个
文件里），故这条与应力无关，是纯粹的标签映射问题。

### vanhove：默认 1 个原点 + 仍用格点视图解包裹（2026-09-27 发现）

固定 lag 的位移直方图，FFT 无用。~~`tau` 默认 `n_steps − 1` → 只有 1 个原点~~（2026-10-08 随 B-5 改为 N/2）；解包裹
仍是 `unwrap_frac`（分数坐标解包裹 × 当帧盒子，即格点视图，NPT 下有 msd 手册里写的
问题）。改默认 tau（N/2？）并换 TOR —— `msd.rs::unwrap_tor` 届时有第二个用户，
可下沉共享。

### extxyz 把 `momenta` 当速度读（2026-09-27 发现）

`readers/extxyz.rs` 找不到 `velocities` 列时退到 `momenta`，但直接当速度存，没有除以
质量（ASE 的 `momenta` 是 m·v，且时间单位是 ASE 自己的）。ASE 写 extxyz 时带的正是
`momenta`，所以 ASE 产出的轨迹进 `traj vacf` 会得到错误量纲。修法：除以
`effective_mass()` 并换算 ASE 时间单位，或拒绝 `momenta` 并报错。

### bondlife 的后续（2026-09-27 记）

- **氢键几何判据**（供体–H–受体角度）：水–玻璃界面的 H 键需要它；现只有距离判据
- **cell list**：候选搜索是逐帧 中心数 × 邻居数 全扫；上万原子 × 上千帧会慢
- **Luzar–Chandler 速率常数** k、k′：需另定义 n(t)（键断但仍在第二阈值内的对），
  gmx hbond-legacy 用 `-r2`
- Yamamoto–Onuki 严格的端点定义（t₀ 时 r ≤ A₁、t₀+Δt 时 r ≤ A₂，不看中间）与现行
  滞后序列只在 r_break > r_bond 时略有差别；如需逐字复现可加 FFT 互相关

### VDOS：振动态密度（2026-09-27 提出，中）

VACF 的傅里叶变换，业内常用（如 LAMMPS/MDAnalysis 生态的 power spectrum）。等 vacf
改 FFT 后顺手做：窗函数（Hann？）、归一化、频率单位（THz / cm⁻¹）需先查规范。

### 共享头部冒充全批：angle / vanhove（2026-09-27 扫描；vacf、rotcorr 已修）

与 msd 同类（`fa2727c` 已修 msd）：`meta_lines()` 里用了逐文件的量，多文件时第一个
文件的值摆在共享区。

| 分析 | 写进共享区的逐文件量 |
|---|---|
| `angle` | `elements` 列表、`[statistics]` 段（每个三元组的 mean/std/count） |
| `vanhove` | `elements`、`n_atoms`、`n_origins`、`tau_frames`、`time` |

修法同 msd：共享区只留参数，其余走 `Summary` 列或 `note()`。`angle` 的统计是逐文件
逐三元组的二维量，放 `[inputs]` 会很宽，可能要单独一张表。

### D 的统计误差：GLS / 贝叶斯（2026-09-27 提出，低）

`d_err` 是 gmx msd 的半窗差，本质是线性度检查。严格误差需要 MSD 协方差，见 kinisi
（PMC11736684）。先解决上面「只有 1 个原点」，否则协方差无从谈起。

### `ferro map sdf` 导出旋转后结构：压缩 npz（2026-09-26 提出）

用 `ndarray-npy` 的 `NpzWriter::new_compressed`（默认 feature `compressed_npz`，
依赖全为纯 Rust，Linux/Windows 无系统库要求）——**所以 ndarray-npy 不关默认
feature**。选型依据：

| 格式 | Python 端 | Rust 端代价 | 结论 |
|---|---|---|---|
| npz（deflate） | `np.load`，只需 numpy | 已在树里，新增 0 | 采用 |
| npy + zstd | 另装 `zstandard` 并自拼 header | zstd 是 C 绑定 | 否 |
| HDF5 | `h5py` | 链 libhdf5，跨平台最难 | 否 |

体积主要由精度决定而非算法：f64 坐标尾数近随机，deflate 只压 10–30%；**存 f32
直接减半**（坐标误差 ~1e-5 Å，dpdata 默认也是 float32）。建议默认 f32 + 开关留
f64。数组命名与形状在实现时定。

### 搁置项

- **`vanhove` 加 `tau` 列**：现在一次只算一个 τ、写在 `#` 头里。加列后将来支持多 τ 是
  加行而非改列结构

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

### ferro-workflow：VASP（POSCAR/INCAR/KPOINTS）

---

## 优先级低

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

### ferro-python：pyo3 0.29 运行时验证（2026-09-15 由中降低）

0.21 → 0.29 只做了类型层验证。本机无 maturin，`cargo build` 在 macOS link 阶段过不了
（`extension-module` 需 `-undefined dynamic_lookup`）：

```bash
pip install maturin
cd ferro-python && maturin build --interpreter "$(which python)"
pip install target/wheels/*.whl
```

冒烟测试要覆盖：读 xyz/cif/lammpstrj、`supercell`、`write`、
**新拆的 `gr_pair` / `gr_all`（含 `by="label"`）**、`msd`。

### 是否采用 rustfmt（待定，暂缓；2026-09-15 由中降低）

收益：把格式从代码审查面里移除。代价三点：

1. 一次性约 80 文件的大 diff，`git blame` 在这些行上全部指向那一次提交
2. 会拆掉现有的刻意列对齐（`gr.rs` 的 `WelfordStats`、`angle.rs` 的 `CellList`、
   `units.rs` 的枚举表）
3. 生成文件要显式排除：`rustfmt.toml` 的 `ignore` 仅 nightly 可用，stable 上需给
   `cp2k_basis_db.rs` 的静态表加 `#[rustfmt::skip]`

若采用：单独一次纯格式提交（`style: adopt rustfmt`）+ 固定 `rustfmt.toml`，
不要混进特性或依赖升级提交。

### ferro dataset filter：多层周期镜像（2026-08-25 提出，2026-09-25 降低）

几何判据（`--oo-min`、`-r` Al6）现在把最小镜像当上界（超出报错），没做
`private/filter_dataset.py` 的多层扫描（`n = floor(rcut / w + 0.5)`）。当前体系
盒子约 15 Å、截断约 2.4 Å，按该公式 n = 0，两种做法逐位相同；只有小胞体系才需要。
Al6 走 `ferro-core` 的 `classify_frame`，要做多层得动 net 共用的近邻搜索，面更大。

上界检查原先只看第 0 帧，NPT 缩盒时会漏判；2026-09-25 改为逐帧取最紧的一帧
（`5833ff6`），报错点名帧号。这一半已做完，剩下的只是多层扫描本身。

### 机器学习集成

- 阶段一：`linfa`（K-Means/DBSCAN/PCA）→ `ferro-analysis/src/ml/`
- 阶段二：`candle`（ONNX 推理，加载 DeepMD-kit 势函数）
- 阶段三：`burn`（纯 Rust 训练，远期）

### 深度学习工作流 I/O

DeePMD-kit 的 npy 系统目录**写**侧已完成（`writers/deepmd.rs`，2026-08-25），
读侧与 GPUMD/NEP 的 `train.xyz` 导出待做（后者是链末的 export，不是中间格式）。
MACE/NequIP 兼容格式仍未开始。

---

## 已完成（归档）

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


### 早期（2026-05 ~ 06）

| 日期 | 内容 | 落点 |
|---|---|---|
| 05-10 | Bader 电荷分析 | `charge_grid.rs`、`chgcar.rs`、`dft/bader*.rs`；规格见 `bader.md` |
| 05-14 | network 重构 + 类型分类迁移 | `network_type.rs`、`typing.rs`、`cluster.rs` |
| 05-15 | 电荷密度团簇 SDF | `dft/chg_sdf.rs`（Kabsch + pull 插值旋转） |
| 05-15 | CP2K 输入生成 | `workflow/cp2k.rs`，混合泛函自动生成 `&HF` 块 |
| 05-16 | 未成对电子 + 基组库 + QE | `spin.rs`、`cp2k_basis_db.rs`（2829 条）、`qe.rs` |
| 05-16 | cube_sdf 迁移 + 全项目审计 | 共享原语下沉 `core/cluster.rs`，**ferro-analysis 去掉 petgraph**；审计 9 项修 8（#8 为误报） |
| 05-16 | ferro-python PyO3 绑定 | 独立 workspace，旧占位代码全部重写 |
| 05-17 | MSD 绘图 + 自扩散拟合 | `fit_diffusion`（D=slope/6，Einstein 3D）+ R²；拟合区间是**滞后时间轴的分数** |
| 06-26 | 代码审查修复三项 | g(r) r_max clamp、`bader_ongrid` 删重复循环、`box_builder` 近邻去重 |

### g(r) / CN / S(q) 重构 + element/label 拆分（2026-08-08，0.1.10–0.1.11）

提交 `caa309b`（0.1.10）· `c5a48b6`…`fc0a3b0`（0.1.11）。测试 312 → 334。

需求来源：`-a P -b O` 与 `-a O -b P` 输出完全相同，`-a`/`-b` 顺序对 CN 不起作用。

定下的语义（后续一直沿用）：

- **`gr` 对称、`cn` 有向**：`CN(A→B) = hist/(N_A·steps)`，`-a` = 中心、`-b` = 近邻
- **未加权 total S(q) 删除**：`f_i ≡ 1` 的退化情形对应不了任何实验探针，
  参考实现 `examples/code2` 中根本不存在
- **配对模型改 n² 个有序对**（3 元素 → 9），镜像对的 `gr` 数值重复照写，
  换取「每配对一组可直接提取的列」
- **拆分规则按第一个下划线**，不用贪婪前缀匹配 —— 否则 `Pb`→铅、`Po`→钋 会盖过
  「P bridging」这类伪标签意图。只作用于 LAMMPS dump；cp2k/qe 的 `extract_element`
  与 cif 的 `element_from_label` **不动**（它们处理 `Fe1`/`O2` 原生命名，字母前缀规则
  对其正确，套用下划线规则反而让 `Fe1` 退化成元素 `Fe1`）

明确不做（当时）：`neighbors_of` 每对枚举两次的 2× 浪费（纯性能）；三处
`extract_element` 的统一；旧标签格式（`P0`/`Ob`/`On_P`）的读取兼容层 —— 重跑一次即为
新格式，加格式猜测反而可能把真元素误判成旧标签。

### 依赖包全量升级（2026-08-08，0.1.12）

`ebd6147` → `915412b`。陷阱见 `issues.md`「依赖升级陷阱」。
`nalgebra` 0.34→0.35 后 g(r)/CN 输出与升级前**逐字节一致**。

### NPT 逐帧体积归一化（2026-08-08，0.1.13）

口径由「先时间平均、后归一化/变换」改为对齐 code1/code2 的「先逐帧归一化/变换、
后时间平均」。完整推导、实测差异与「明确不做」见 `issues.md`。

要点：关键恒等式 `ρ_f·g_f` 中 V_f 自行抵消 → 逐帧变换可由 `⟨ρ_f·g_f⟩` 与 `N·⟨1/V⟩`
两个时间平均量**精确**重构，无须每帧各做一次 FT，分层与 `calc_sq_from_gr` 签名不动。

fixture 由生产轨迹**等间隔**取 5 帧（连续取会丢掉体积跨度，测不出东西）。

### box_builder 括号公式（2026-08-08，0.1.14）

表面是「`parse_formula` 不支持 `Ca3(PO4)2`」，**实际阻塞点在上游**：`build_box` 只把
数据库里的 `cd.formula` 喂给解析器，而库中 26 条全是无括号有机溶剂；用户输入的
compound 一旦不在库中，`compounds::find` 返回 `None` 就报错了，**根本到不了解析器**。
只加括号支持等于写死代码。故一并做了三件事：栈式解析、`resolve_component`
（库外化合物当化学式解析）、收紧输入校验。

### 批处理 + Table 长表 + 单二进制重构（2026-08-09，0.2.0）

锚点 tag `v0.1.15`。测试 362 → 396。

值得记住的判据：

- **`-i` 恒为 `Vec`，单一代码路径**，N=1 是 N 的特例。按文件数分派会让产物形态取决于
  glob 当天匹配到几个文件，还意味着两套命名、两套绘图、两套错误处理
- **`Table` 下沉 core 而非把 writer 搬进 io**：完整论证见 `issues.md`
  「分析产物为什么不进 ferro-io」
- **表结构跟主产物的粒度走** —— 这是 gr/angle 长表而 sq 宽表的判据。决定性理由是
  不同轨迹的**元素集可能不同**（Zn-P-O 与 Al-P-O 同批跑），配对方向做宽表就要取列
  并集补空洞
- **`to_tables()` 定为固有方法不是 trait**，留待第二个消费者出现
- **不拆「纯搬家」提交**，`Table` 迁移 + 批处理 + 长表一次走完
- **`--plot` 冻结为自检用途**，不追 matplotlib。任何「加对数轴/误差棒」的需求一律
  指向 Python（长表 + `sns.lineplot(hue="file")` 一行）。**2026-09-26 进一步整体移除**，
  见归档「依赖精简」

实施中比原计划多出的：`angle` 多一列 `count`（整数直方图是与 dump2analysis 逐 bin
对拍的依据，只留归一化的 `p` 会废掉这条验证路径）。

### ferro net 重构：结构化类型 + 标签重做 + 合并命令 + linkage（2026-08-12，0.2.1）

原计划的「标签体系重做」与「接入批处理 + 长表化」两项，实施时发现它们不是两件事的
先后，而是同一件事的两层：**标签之所以改不动，是因为分类结果以字符串形式流转。**

| # | 提交 | 内容 |
|---|---|---|
| 1 | `5f8ac0e` | `classify_frame` 返回结构化 `AtomType`，消除**五处**标签反解析 |
| 2 | `4503d11` | 标签改 `<元素>_<后缀>`；删修饰子角色分类 |
| 3 | `f1f929d` | `net` 降为叶子命令，接 `CommonArgs` + 批处理 + 长表 |
| 4a | `801184b` | 分布表加 `sd` 列（Welford）；氧的伙伴改数据列 |
| 4b | `304bdd9` | linkage 长表 + `Q^n(mAl)` 分解 |
| 5 | `2f9021e` | 标签存 `atom.label`；extxyz 加 `label:S:1`；dump 折叠 + type 编号跨帧固定 |

**第 1 步单独成一提交是关键**：它把「信息不再经字符串往返」与「标签换格式」分开，
前者可逐字节对拍验证（22 个产物文件一致），后者才是行为变更。合到一起就没有能对拍
的中间态。

被实测推翻或修正的原计划：

- **`X` 兜底的分裂比预期严重**：计划只说「≥3 配位氧与 ≥3 NBO 修饰子塌缩」，实测单帧
  `X = 181` 里 **180 个是 Zn、1 个是氧**
- **修饰子角色分类整体删除**（原计划是保留 `Zn_f/_t/_b`）：实测 97.1% 落进兜底桶
  （0 / 1 / 26 / 903），该分档对 Zn 无分辨力
- **`qn` → `n_bridge` 这一步后来被用户纠正**：「Al 没有 Qn」的意思是 Qn 计算直接忽略
  Al，而不是把列改名。改列名保留了 Al 的行，等于用一个更含糊的名字继续报一个对 Al
  无意义的量。正确做法是让 Al 退出 Qn 表、列名退回 `qn`
- **动机比计划写的窄**：`net type 导出 → traj gr -x P_3 -y O_b 直接串起来`
  只在**单帧**成立，多帧被粒子数守恒守卫拒绝

### ferro net 产物可读性重做（2026-08-12，未发版）

起因是三条反馈：帮助文档太长；四个文件名看不懂；表里把 label 拆成 `former` +
`n_bridge` 反而更难读。追问下去牵出上面那条更根本的「Al 没有 Qn」。

提交 `450252d` · `100265d` · `e2821c6` · `c26cde6`（帮助 120→58 行）· `9536c1b`。
第二轮反馈续做：`a0d21d6`（单元/原子两套词汇、`average`→`composition`、
`ligand_type` label 合并、配体分母改逐元素）。

设计判据（grilling 里逐条过过，坑的部分见 `issues.md`）：

- **`label` 新增而非替换数值列**。label 给人读，`former`/`qn`/`cn` 给筛选和画图。
  删掉数值列会让「筛出 Qn ≥ 3」退化成字符串切分 —— 正是第 1 轮从五处调用点消灭掉的
  反向解析，不该在输出侧重新造一个
- **展示列的数字含义按元素而变是可以接受的**（`Al_4` 是配位数、`P_2` 是 Qn），因为
  那是文献自己的读法。曾以「标签要被 `traj gr -x` 机器解析」为由否掉，后来发现站不住：
  用户要的正是 `-x Al_5` 选出五配位 Al
- **详细说明进文件头**。帮助会滚走、手册在别处，`#` 头跟着文件走
- **用户否掉了自己最初提的 `Q3(2Al)` 写法**：实测参考数据里 `(qn=2,m_Al=1,m_P=0)` 与
  `(qn=2,m_Al=1,m_P=1)` 会渲染成同一个 `Q2(1Al)`，前者的 Σm 亏空 1（那座桥是三簇氧）
- **`linkage` 保持原子词汇**，理由由用户给出且比一致性论证更强：桥联表达的是
  **原子之间**的连接，Qn 是包含多个原子的**单元**

对拍（`43Z43P15A_NPT_5`）：P 的 Qn 计数 25/276/653/751/155 不变、linkage 总观测 3589
不变、`ligand_type` 与 `coordination` 逐行相同、两张 Qn 表 count 闭合。导出轨迹读回后
`traj gr -x Al_5 -y O_b` 得 CN = 4.933，缺的 0.067 是一个三簇氧（标 `O_t` 不标 `O_b`）
—— 另一条代码路径的交叉验证。测试 398 → 420。

### traj 产物命名带 label + --outdir + sq 移除选择（2026-08-13，未发版）

起因是一句具体的抱怨：`gr.csv` 看不出算的是哪一对。范围从「gr 和 angle」扩到六个命令，
并顺带拆掉了 sq 的配对选择。提交 `3062954` · `2646232` · `5234fba`。

判据：

- **label 排在 suffix 之前**。label 说的是算了什么，suffix 是批次标记；这个顺序让
  `ls gr_P-O_*` 列出同一对在各批次的结果，反过来没有对应的用法
- **两种拼法并存且各有理由**：`gr`/`angle` 按写的顺序，`--elements` 排序去重。
  看着不一致，但统一成任何一种都会错一半（详见 `issues.md`）
- **`_all` 是有代价的选择**。它让所有不带筛选的旧命令产物改名，换来「文件名一眼看出
  有没有筛选」。用户明确选了这一边
- **sq 移除选择的代价已知并接受**：按 label 分辨的 partial 从 CLI 消失（`-x/-y` 曾是
  进入 `GroupBy::Label` 的唯一入口）。理由是位点标签对应的原子数往往不足以让 partial
  显出信号。**库层 `GroupBy::Label` 不动**

实施中发现与设计讨论不符的两处：`rotcorr` 的 `--center`/`--neighbor` 其实是必填的
（讨论里假设的 `rotcorr_all.csv` 那条路根本走不到）；删掉 `SelectArgs` 的**使用**不等于
删掉参数（详见 `issues.md`）。

**`vacf` / `vanhove` / `rotcorr` 的改名未经实测验证** —— 参考 fixture 没有速度，
`vacf` 跑不到写文件那步；这三个性质本身也还没调试到。单元测试覆盖了拼名函数，
但端到端只验了 gr / angle / msd / sq / net / map。

### scripts/：四个对拍脚本修复（2026-08-11）

0.2.0 的输出变更让四个对拍脚本全部失效。典型断点：`compare_rdf.py` 的
`[float(x) for x in line.split()]` 遇到 `file` 文本列直接 `ValueError`。

**新增 `scripts/ferrocmp.py`** 收拢共用逻辑（四个脚本各有一份几乎相同的 `run` /
`load_columns` / `fe_version` / 列名行解析，形状一变就要改四遍）。陷阱见 `issues.md`
「外部脚本调用 ferro 的陷阱」。

复跑验证（数值与 `issues.md` 记录逐项吻合，说明改的只是读法不是口径）：

| 脚本 | 结果 |
|---|---|
| `compare_rdf.py` | P-O / Al-O 峰值与峰位完全相同，max\|Δ\| 5e-5 / 5e-4（dump2analysis `%12g` 的量化台阶） |
| `compare_angle.py` | Σfe/Σref = 0.5000（×2 计数约定）；`--align-binning` 下 1800 个 bin **整数零差** |
| `compare_sq.py` | 未补偿残差 q>15 均值 1.00031 / 1.00018（+1 常数）；partial 互相关 +0.959（type_new 失效） |
| `compare_sq_experiment.py` | 50 帧 rms(fe−exp) 0.0194/0.0277、max\|fe−ref\| 0.0027/0.0008、FSDP 1.95/2.05 |

顺带修：`compare_sq_experiment.py` 的 `TRAJ_CANDIDATES` 首选指向一个**不存在**的文件，
一直在静默回退到 5 帧子集。

### scripts/：对拍脚本跟进产物改名（2026-08-14）

产物加 label 段后三个脚本断了（`compare_rdf` / `compare_angle` / `compare_sq` 的 gr 段），
根因不是那三个字符串写错，而是**四个脚本各自手写产物名**。故做法是把命名规则搬进
`ferrocmp.py`：`file_label()` / `set_label()` / `product_name()` 是 `batch::out_path`
与 `batch::file_label` 的镜像，调用点只说「哪个模式、什么 label、什么后缀」。

- **`-o` 不再重复配对**。原先 `-o P-O` 是唯一区分手段，现在配对已在 label 段里，
  再塞进后缀就成了 `gr_P-O_P-O.csv`。四个脚本统一 `-o cmp` / `-o exp`
- 顺带修：`--ferro` 传相对路径必炸（`run_ferro` 以 outdir 为 cwd），含分隔符时先 `resolve()`

复跑四个脚本，数值与 2026-08-11 那轮**逐项吻合**（`max|Δ|` 5e-5/5e-4、Σfe/Σref
0.5000、q>15 均值 1.00031/1.00018、partial 互相关 +0.959、rms 0.0194/0.0277、
FSDP 1.95/2.05），说明改的只是拼名不是口径。

### scripts/：四个发表级绘图脚本（2026-08-13）

与 `compare_*.py` 的分工是清楚的：那边**对拍**（跟参考实现逐点比，图只为看差异），
这边**出图**（进论文的 pdf）。两者都读 ferro 的 csv，读完之后没有共同代码，所以是
两个共享层而不是一个。

| 脚本 | 子图 | 曲线 / 条带 |
|---|---|---|
| `plot_gr.py` | 一个 csv 一张 | 一条轨迹一条曲线；左轴 g(r) 实线、右轴 CN(r) 虚线 |
| `plot_angle.py` | 一个 csv 一张 | 同上，单轴 P(θ) |
| `plot_sq.py` | 一个 csv 一行两格 | 左 S^N(Q)、右 S^X(Q)，两条 total |
| `plot_net.py` | 一个形成子 / 元素一张 | x = 成分（`file` 列），100 % 堆积柱 |

判据：

- **net 的 x 轴是成分不是 Qn**。同一张图上看「各 Qn 占比随成分怎么变」，这是堆积柱
  相对折线的全部理由：它把「和为 1」画成图形约束。代价是**小分量看不出趋势**
  （2 % 的条带只剩一条线），要追小分量就把那列单独拉出来画折线
- **x 轴顺序 = `-i` 的参数顺序**，不从文件名解析成分数值。`43Z43P15A` 里有三个数字，
  脚本无从知道要哪个，猜错了图是错的但看不出来
- **`--partner` 是开关不是默认**。单形成子体系下 `m_P ≡ qn`，展开只会画出一条假的两级结构
- **CMD/MLMD 的方法维走文件名**，画成同一成分刻度下并排多根柱。这也是嵌套条带
  （而非双条带）方案的理由：并排那个位置要留给方法维
- **gr/angle/sq 不为方法对比设计**（用户判断是极小概率需求）；**不画误差棒**
  （`sd` 是快照间散布不是标准误）

编码陷阱见 `issues.md`「发表级绘图脚本编码陷阱」。
