# CLI Reference

**一个二进制 `ferro`**，子命令按**产物**分组（不是按实现它的 crate）。
0.2.0 起原来的八个 `fe-*` 二进制已全部删除，不留兼容层——输出格式同期变更，留着
`fe-traj` 会让旧脚本「跑成功」却吐出自己解析不了的 csv，静默坏数据比命令消失难查。

```
ferro traj  gr | sq | msd | angle | vacf | rotcorr | vanhove   → 堆叠 csv + 可选 PNG
ferro map   density | velocity | force | radius | sdf | chg-sdf → 逐输入一个 .cube
ferro net                                                      → 五张堆叠 csv
                                                                 + 可选标注轨迹
ferro bader | convert | info | job
```

**帮助分三级**：`ferro` 列出分组；`ferro traj` 列出该组命令；`ferro traj gr`（不给
`-i`）打印该命令的参数、输出列结构与示例。

| 旧命令（0.2.0 前） | 新命令 |
|---|---|
| `fe-traj -m gr` | `ferro traj gr` |
| `fe-corr -m vacf` | `ferro traj vacf` |
| `fe-cube -m density` | `ferro map density` |
| `fe-cube -m chg_sdf` | `ferro map chg-sdf` |
| `fe-network -m Qn` | `ferro net` |
| `fe-network -m type` | `ferro net --export-traj` |
| `fe-bader` / `fe-convert` / `fe-info` / `fe-job` | `ferro bader` / `convert` / `info` / `job` |

---

## Common Flags

`traj` / `map` / `net` 的每个命令都 flatten 了同一组 `CommonArgs`：

| Flag | Description |
|---|---|
| `-i <FILE>...` | 输入文件，**多值**并自展开 glob 模式（引号括起来，让 shell 别动它） |
| `-o <SUFFIX>` | 输出**文件名后缀**，不是路径：产物落在 `<命令>[_<表>]_<后缀>.csv` |
| `--last-n N` | 只用尾部 N 帧（跳过平衡段） |
| `--ncore N` | 并行线程数（默认全部核心） |
| `--metal-units` | LAMMPS metal 单位（速度 Å/ps，力 eV/Å）。**只影响速度与力**，坐标与晶胞两种单位下都是 Å，故对 gr/sq/msd/angle/rotcorr/vanhove/net 无影响 |

### 批处理

`-i` 恒为多值。**只有一条代码路径**：单输入是 N=1 的特例，不是特殊模式——按文件数
分派会让产物形态取决于 glob 当天匹配到几个文件。

```bash
ferro traj gr -i 'runs/*/prod.lammpstrj' -a P -b O -o scan
```

每个输入独立分析，结果堆叠成**一份**带 `file` 列的 csv。元素集不同的输入取列并集，
**缺的留空（NaN），不补零、不插值**。失败的输入被跳过、在输出的 `[inputs]` 块里留下
原因、并使**退出码为 1**（否则 shell 里 `&&` 串联会把批内失败当成功）。
`{a,b}` 花括号不支持，交给 shell 展开。

产物是逐输入的命令（`ferro map` 的 cube、`ferro net --export-traj` 的轨迹）例外：
文件名必须掺输入 stem，否则第二个输入会覆盖第一个。

### 类型选择（gr / sq / angle）

两组互斥：

| Flag | 含义 |
|---|---|
| `-a` / `-b` / `-c` | 按 `Atom::element` 选（元素） |
| `-x` / `-y` / `-z` | 按 `Atom::label` 选（位点标签） |

第一个槽是中心（pair）或端原子 A（triplet）。**顺序有意义**：`g(r)` 对称但 `CN` 有向。

---

## `ferro convert`

格式转换。

```bash
ferro convert -i input.xyz -o output.pdb
ferro convert -i input.cif -o POSCAR
```

支持格式：`.xyz`、`.extxyz`、`.pdb`、`.cif`、`POSCAR`/`CONTCAR`、LAMMPS dump、
LAMMPS data、CP2K、QE。

**元素列始终写干净的元素符号**，无论 `Atom::label` 是什么。只有
`ferro net --export-traj` 会把标签折进 LAMMPS dump 的元素列。

---

## `ferro info`

打印结构 / 轨迹摘要（晶胞参数、原子数、元素组成）。

```bash
ferro info -i input.xyz
ferro info -i traj.lammpstrj --last-n 1
```

读取带位点标签的 LAMMPS dump 时会打印一次 element/label 的拆分映射表。

---

## `ferro job`

为 **Gaussian**、**CP2K**、**Quantum ESPRESSO** 生成输入文件。不给 `-s` 打印总览；
给了 `-s <software>` 但不给 `-i` 打印该软件的专属帮助。导览见
[Job Builders](workflow/job-builders.md)。

```bash
ferro job                                    # 总览
ferro job -s cp2k                            # CP2K 专属帮助
ferro job -i input.xyz -s gaussian -m B3LYP -b 6-31G* -o job.gjf
ferro job -i input.xyz -s cp2k --task geo-opt --functional pbe --dispersion d3bj
ferro job -i Fe2O3.cif -s qe --auto-spin --kpoints 4 4 4 -o pw.in
```

### Charge / Spin（三种目标共用）

| Flag | Default | Description |
|---|---|---|
| `--charge` | (from file) | 覆盖体系总电荷（在自旋推断之前生效） |
| `--multiplicity` | (from file) | 强制多重度 2S+1；优先级最高，会关掉 auto-spin |
| `--auto-spin` | off（cp2k/qe 默认开） | 从结构推断多重度 |

推断链（magmom → 氧化态 + Hund 规则 → 电子数奇偶下限）见
[Spin Estimation](workflow/spin.md)。

### Gaussian

| Flag | Default | Description |
|---|---|---|
| `-m <method>` | (required) | DFT 方法，如 `B3LYP`、`PBE0` |
| `-b <basis>` | (required) | 基组，如 `6-31G*`、`def2-TZVP` |
| `-o <file>` | `job.gjf` | 输出文件 |

### CP2K

#### 任务与电子结构

| Flag | Default | Candidates |
|---|---|---|
| `--task` | `energy` | `energy`, `force`, `geo-opt`, `cell-opt`, `md`, `freq` |
| `--functional` | `pbe` | `pbe`, `blyp`, `pbe0`, `b3lyp`, `revpbe`, `pbesol`, `scan`, `r2scan`, `hse06` |
| `--cp2k-basis` | `dzvp-molopt-sr` | `dzvp-molopt-sr`, `tzvp-molopt`, `tzv2p-molopt`, `dzvp-gth`, `tzvp-gth`, `pob-dzvp`, `pob-tzvp`（全电子），或任意自定义字符串 |
| `--dispersion` | `none` | `none`, `d3`, `d3bj` |
| `--scf` | `diag` | `diag`（金属/大体系）, `ot`（绝缘体） |
| `--pbc` | (auto) | `xyz`, `z`, `none`；省略时从晶胞自动判断 |
| `--kpoints` | (none) | 三个整数，如 `--kpoints 2 2 2` |
| `--cutoff` | `400` | 平面波截断 [Ry] |
| `--rel-cutoff` | `50` | 相对截断 [Ry] |
| `--smear` | off | 启用 Fermi–Dirac 展宽 |

#### 输出

| Flag | Default | Candidates |
|---|---|---|
| `--atom-charge` | `none` | `none`, `mulliken`, `hirshfeld`, `hirshfeld-i` |
| `--cube` | `none` | `none`, `density`, `elf`, `hartree` |
| `--molden` | off | 导出 Molden 轨道文件 |
| `--project` | `ferro` | CP2K project 名 |

#### MD（仅 `--task md`）

| Flag | Default | Description |
|---|---|---|
| `--md-steps` | `10000` | MD 步数 |
| `--md-timestep` | `1.0` | 步长 [fs] |
| `--temperature` | `298.15` | 温度 [K] |
| `--thermostat` | `csvr` | `csvr`, `nose`, `langevin`, `none` |
| `--traj-freq` | `100` | 轨迹写出频率 [步] |
| `--barostat` | off | 启用 NPT 压浴 |

> 基组与赝势名**逐元素**从 2829 条数据库解析（PBE / SCAN / 全电子，价电子数 `q`
> 一致）。`--cp2k-basis` 选族，元素专属名自动填入。见
> [Job Builders](workflow/job-builders.md#precise-basis--pseudopotential-matching)。

### Quantum ESPRESSO

```bash
ferro job -i crystal.cif -s qe
ferro job -i metal.cif -s qe --smearing mp --kpoints 8 8 8
ferro job -i slab.xyz -s qe --qe-task relax --qe-functional scan -o pw.in
```

| Flag | Default | Candidates / Description |
|---|---|---|
| `--qe-task` | `scf` | `scf`, `nscf`, `bands`, `relax`, `vc-relax`, `md`, `vc-md` |
| `--qe-functional` | `pbe` | `pbe`, `pbesol`, `revpbe`, `blyp`, `scan`, `r2scan`, `pbe0`, `hse06` |
| `--ecutwfc` | `50` | 平面波截断 [Ry] |
| `--smearing` | `none` | `none`, `gaussian`, `mp`, `mv`, `fd`（金属用 mp/mv） |
| `--kpoints` | (Gamma) | 三个整数 → Monkhorst-Pack 网格 |
| `--pseudo-dir` | `./pseudo` | 赝势目录（`<El>.UPF`） |
| `--md-steps` | `10000` | MD 步数（`--qe-task md`/`vc-md`） |
| `--temperature` | `298.15` | MD 目标温度 [K] |
| `-o <file>` | `pw.in` | 输出文件 |

`ibrav = 0`；晶胞按 `CELL_PARAMETERS angstrom` 从结构写出。自旋走共用推断器 →
`nspin` / `tot_magnetization`。

---

## `ferro traj`

七个轨迹分析，共用同一条导出管线：一份长表或宽表 csv + 可选 PNG。

```bash
ferro traj <command> -i traj.lammpstrj [flags] -o <suffix>
```

### `gr` — 径向分布函数

$g(r)$ 与配位数 $\text{CN}(r)$。

```bash
ferro traj gr -i traj.lammpstrj -a P -b O --r-max 10.0 --dr 0.002 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--r-min` | 0.001 | 最小半径 [Å] |
| `--r-max` | 10.005 | 最大半径 [Å]；clamp 到最小**面间距**的一半（非最短边长） |
| `--dr` | 0.002 | 分箱宽度 [Å] |
| `--plot` | off | 另写 PNG（两格：g(r) \| CN(r)），需要指定配对 |

**长表**：`file, r, center, neighbor, gr, cn`。类型进数据列，故元素集不同的轨迹可直接
堆叠；不给 `-a/-b` 是加行而不是加列。`gr` 对称（`A-B` == `B-A`），`cn` 有向
（`CN(A→B)`），这个区别写进了 `center`/`neighbor` 两列而不是文档注脚。

### `sq` — 结构因子

对 $g(r)$ 做傅里叶变换得到 $S(q)$。

```bash
ferro traj sq -i traj.lammpstrj --q-max 25.0 --dq 0.02 --weighting both -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--q-min` | 0.1 | 最小 $q$ [Å⁻¹] |
| `--q-max` | 25.0 | 最大 $q$ [Å⁻¹] |
| `--dq` | 0.02 | $q$ 分箱宽度 [Å⁻¹] |
| `--weighting` | `both` | `none`, `xrd`, `neutron`, `both` |
| `--plot` | off | 另写 PNG（两格：XRD \| Neutron） |

`gr` 的 `--r-min` / `--r-max` / `--dr` 同样生效——它们决定被变换的那条 $g(r)$ 的范围。

**宽表**：`file, q, total_xrd, total_neutron`，其后每个配对三列
（`_sq` / `_xrd` / `_neutron`，只出规范半边）。主产物是两条 total（一行一个 $q$），
加权 partial 是能求和还原 total 的诊断分解。

### `msd` — 均方位移

```bash
ferro traj msd -i traj.lammpstrj --dt 2.0 --shift 10 --elements Li --fit-range 0.3,0.8 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--dt` | 1.0 | 步长 [fs] |
| `--shift` | 1 | 时间原点间隔 [帧] |
| `--elements` | (全部) | 逗号分隔的元素过滤 |
| `--fit-range` | (无) | `FMIN,FMAX` 线性拟合窗口（轨迹分数）→ 自扩散系数 D |
| `--plot` | off | 另写 PNG（2×2：total \| a \| b \| c） |

给了 `--fit-range` 即计算并打印 $D = \text{slope}/6$ 与 $R^2$（与 `--plot` 无关）。

### `angle` — 键角分布

```bash
ferro traj angle -i traj.lammpstrj -a O -b P -c O --r-cut-ab 2.4 --r-cut-bc 2.4 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--r-cut-ab` | 2.3 | 端 A 到中心 B 的截断 [Å] —— A 是 `-a`/`-x` 给的那个 |
| `--r-cut-bc` | 2.3 | 端 C 到中心 B 的截断 [Å] —— C 是 `-c`/`-z` 给的那个 |
| `--angle-min` | 0.0 | 直方图下界 [°] |
| `--angle-max` | 180.0 | 直方图上界 [°]（闭区间） |
| `--d-angle` | 0.1 | 分箱宽度 [°] |
| `--plot` | off | 另写 PNG，图例含 mean ± std |

不指定三元组时两个截断回落到规范 (Z, 符号) 顺序；两端同类型时两者都取
`min(--r-cut-ab, --r-cut-bc)`。区间外的角被**丢弃**，不只是不显示。

**长表**：`file, angle, end_a, center, end_c, count, p`。保留整数 `count` 与归一化
`p` 两列——整数直方图是与 `dump2analysis` 逐 bin 对拍的依据。详见
[Bond Angle Distribution](analysis/angle.md)。

### `vacf` — 速度自相关

```bash
ferro traj vacf -i traj.lammpstrj --dt 2.0 --elements Li --metal-units -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--dt` | 1.0 | 步长 [fs] |
| `--shift` | 1 | 时间原点间隔 [帧] |
| `--tau` | (全部) | 滞后窗口 [帧] |
| `--elements` | (全部) | 元素过滤 |

列：`file, time, vacf, vacf_x, vacf_y, vacf_z, diffusion`

### `rotcorr` — 转动相关

分子取向矢量的 $C_2(t)$。

```bash
ferro traj rotcorr -i traj.lammpstrj --center P --neighbor O --r-cut 2.4 --dt 2.0 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--center` | (required) | 中心原子元素 |
| `--neighbor` | (required) | 近邻原子元素 |
| `--r-cut` | 1.2 | 成键搜索截断 [Å] |
| `--dt` | 1.0 | 步长 [fs] |
| `--shift` | 1 | 时间原点间隔 [帧] |
| `--tau` | (全部) | 滞后窗口 [帧] |

列：`file, time, c2, integral`

### `vanhove` — Van Hove 自关联

```bash
ferro traj vanhove -i traj.lammpstrj --tau 500 --dt 2.0 --r-max 8.0 --dr 0.02 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--tau` | (末帧) | 滞后 [帧] |
| `--dt` | 1.0 | 步长 [fs] |
| `--shift` | 1 | 时间原点间隔 [帧] |
| `--r-max` | 10.0 | 最大位移 [Å] |
| `--dr` | 0.01 | 分箱宽度 [Å] |
| `--elements` | (全部) | 元素过滤 |

列：`file, r, gs`

### 绘图

`--plot` 出一张 PNG 分格，**一格一个量、一条曲线一个输入文件**；颜色按文件跨格一致，
图例只画第一格。500 dpi，一格 2708×2083 px。

`--plot` **冻结在自查质量**，不会去追 matplotlib：数据是长表 csv，一行 seaborn 就是
一张正经图（`sns.lineplot(data=df, x="r", y="gr", hue="file")`），对数轴、误差棒、
主题这些属于 Python。

---

## `ferro map`

3-D 空间分布图（Gaussian cube 格式）。**逐输入一个 `.cube`**，没有可堆叠的表也没有图，
故文件名必带输入 stem（`density_<stem>.cube`）。

```bash
ferro map <command> -i traj.lammpstrj [flags] -o <suffix>
```

### `density` — 原子数密度

```bash
ferro map density -i traj.lammpstrj --nx 80 --ny 80 --nz 80 --elements Li -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--nx/ny/nz` | 50 | 网格维度 |
| `--elements` | (全部) | 元素过滤 |

### `velocity` — 逐体素平均速率

需要轨迹带速度（LAMMPS metal dump 请加 `--metal-units`）。

```bash
ferro map velocity -i traj.lammpstrj --metal-units -o run1
```

### `force` — 逐体素平均力大小

需要轨迹带力。

```bash
ferro map force -i traj.lammpstrj -o run1
```

### `radius` — 硬球占据

```bash
ferro map radius -i traj.lammpstrj --elements Li --radius 0.7 --nx 100 --ny 100 --nz 100 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--radius` | 0.7 | 硬球半径 [Å] |
| `--nx/ny/nz` | 50 | 网格维度 |
| `--elements` | (全部) | 元素过滤 |

### `sdf` — 团簇 SDF

```bash
ferro map sdf -i traj.lammpstrj --qn 3 --former P --ligand O --cutoff-fl 2.4 \
    --modifier Zn --cutoff-ml 2.8 --grid-res 0.1 --sigma 1.5 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--qn` | 3 | 目标 $Q_n$ 级别（0–3） |
| `--former` | `P` | 网络形成子元素 |
| `--ligand` | `O` | 桥联配体元素 |
| `--cutoff-fl` | 2.4 | 形成子–配体截断 [Å] |
| `--modifier` | (无) | 修饰子阳离子元素 |
| `--cutoff-ml` | 2.8 | 修饰子–配体截断 [Å] |
| `--grid-res` | 0.1 | 体素尺寸 [Å] |
| `--sigma` | 1.5 | 高斯展宽 [体素] |
| `--padding` | 3.0 | 网格边界余量 [Å] |
| `--rmsd-warn` | 0.5 | RMSD 告警阈值 [Å] |

产物：逐原子类型一个 `<stem>_<label>.cube`（多族时 `<stem>_fam<N>_<label>.cube`）。

### `chg-sdf` — 电荷密度 SDF

从一组 QE `pp.x` 电荷密度 cube 计算 Qn 团簇周围取向平均的电子密度。**不用 `-i`**，
改接 `--cubes`。

```bash
ferro map chg-sdf --cubes frame_000.cube frame_001.cube frame_002.cube \
    --qn 2 --former P --ligand O --cutoff-fl 2.4 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--cubes <files…>` | (required) | QE pp.x cube 文件，一帧一个 |
| `--qn` | `3` | 目标 Qn 级别（0–3） |
| `--former` | `P` | 网络形成子元素 |
| `--ligand` | `O` | 桥联配体元素 |
| `--cutoff-fl` | `2.4` | 形成子–配体截断 [Å] |
| `--modifier` | (无) | 修饰子阳离子元素 |
| `--cutoff-ml` | `2.8` | 修饰子–配体截断 [Å] |
| `--chg-padding` | `6.0` | 子网格边界余量 [Å] |
| `--rmsd-warn` | `0.5` | 对齐 RMSD 告警阈值 [Å] |

产物：`<stem>_Q<n>.cube`（一签名族一个文件）。算法见
[Averaged Charge-Density SDF](analysis/chg-sdf.md)。

> `--cubes` 是 `ferro map` 里唯一「多输入聚合成**一张** SDF」的模式，与该组其余
> 「一输入一产物」的语义相反。拆分需先定义带样本计数的中间产物格式（跨文件加权平均
> 不可交换），见 `dev/plan.md`。

---

## `ferro net`

玻璃网络拓扑：桥接配体数（P 的 Qn）、配体分类、配位数、桥联统计。

```bash
ferro net -i traj.lammpstrj --P-O=2.4
ferro net -i traj.lammpstrj --P-O=2.4 --Al-O=2.4 --Zn-O=2.6 --modifier Zn
ferro net -i 'runs/*/prod.lammpstrj' --P-O=2.4 -o scan
ferro net -i traj.lammpstrj --P-O=2.4 --last-n 500 --export-traj
```

### 配对参数（必需，至少一个）

截断用 `--<Former>-<Ligand>=<cutoff>` 格式（元素首字母大写）。元素对写在**参数名**
里，clap 建模不了，故 `main` 在解析前先把它们从 argv 剥离。

```
--P-O=2.4     P-O 截断 2.4 Å（P 为形成子，O 为配体）
--Al-O=2.4    同一体系可有多个形成子
--Al-F=2.1    同一形成子可有多种配体
--Zn-O=2.6    配合 --modifier Zn 时视为修饰子–配体截断
```

### 普通参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-i <FILE>...` | — | 输入轨迹（缺省显示帮助） |
| `-o <SUFFIX>` | — | 输出文件名后缀 |
| `--last-n N` | 全部 | 仅用尾部 N 帧 |
| `--ncore N` | 全部核心 | 线程数 |
| `--metal-units` | 关 | 统计不读速度/力；只对 `--export-traj extxyz` 有影响 |
| `--modifier E,E` | — | **只计配位数**的元素，不参与桥接计数与配体分类。须同时给出各自的截断，否则报错 |
| `--export-traj [FMT]` | — | 另写标注轨迹：`lammpstrj`（默认）或 `extxyz` |

### 输出

五张 csv，各带 `file` 列：

| 文件 | 列 |
|---|---|
| `network_bridge.csv` | `file, former, n_bridge, m_<X>…, count, fraction, sd` |
| `network_oxy.csv` | `file, type, former_a, former_b, count, fraction, sd` |
| `network_cn.csv` | `file, element, cn, count, fraction, sd` |
| `network_mean.csv` | `file, element, mean_n_bridge, mean_cn` |
| `network_linkage.csv` | `file, elem_a, n_bridge_a, cn_a, elem_b, n_bridge_b, cn_b, n_formers, count, fraction, sd` |

标签：`P_0`…`P_4` / `Al_2`（数字 = 桥接配体数）、`O_f` / `O_n` / `O_b` / `O_t`、
修饰子为裸元素符号。全部合 `<元素>_<后缀>` 约定。

`--export-traj` 逐输入写 `<输入 stem>_types[_<后缀>].<ext>`。

详细说明（口径、`sd` 的含义、pandas 分析范例、两种导出格式的差别、按标签选型的
单帧限制）见 [Glass Network Analysis](analysis/network.md)。

---

## `ferro bader`

从 DFT 电荷密度做 Bader 电荷分解。支持 VASP CHGCAR 与 Gaussian/QE cube。

```bash
ferro bader -i CHGCAR                      # VASP CHGCAR
ferro bader -i charge.cube                 # Gaussian/QE cube
ferro bader -i CHGCAR -o bader             # 自定义输出 stem
ferro bader -i CHGCAR --method weight      # Yu-Trinkle weight 方法
```

| Flag | Default | Description |
|---|---|---|
| `-i <file>` | (required) | 输入（`.cube` → cube reader；其余 → CHGCAR reader） |
| `-o <stem>` | `bader` | 输出文件 stem |
| `--method` | `ongrid` | `ongrid`, `neargrid`, `offgrid`, `weight` |

### 输出文件

| 文件 | 内容 |
|---|---|
| `<stem>_ACF.dat` | Atomic Charges File —— 逐原子 Bader 电荷、体积、到表面的最小距离 |
| `<stem>_BCF.dat` | Bader Charge File —— 逐 Bader 体积的电荷、体积、坐标 |
| `<stem>_AVF.dat` | Atomic Volume File —— 原子 → Bader 体积索引映射 |

这三个是 Henkelman 组的 bader 格式，外部工具在解析，故不随其余产物迁到 csv。

---

## 输出格式约定

除 `ferro map`（cube）与 `ferro bader`（ACF/BCF/AVF）外，所有产物都是**一份 csv**，
数据上方有一个 `#` 注释块，放共享参数与 `[inputs]` 清单（逐输入的帧数、原子数、
体积、状态）。

那个块是给人看的——`pandas.read_csv(comment="#")` 会丢掉它，所以**脚本必须解析的
东西一律是列**，不会藏在注释里。

数值格式统一 `{:.6e}`，**NaN 渲染为空字段**（列并集下某输入缺的列就是空，不是 0）。
