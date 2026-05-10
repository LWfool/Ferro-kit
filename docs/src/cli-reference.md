# CLI Reference

ferro provides six command-line binaries.  All trajectory binaries share common flags for frame selection and parallelism.  Run any binary without `-i` to print mode-specific help.

## Common Flags

| Flag | Description |
|---|---|
| `-i <file>` | Input file |
| `-o <file>` | Output file (default: mode-specific name) |
| `--last-n N` | Use only the last N frames |
| `--ncore N` | Parallel threads (default: all cores) |
| `--metal-units` | LAMMPS metal units (velocities Å/ps, forces eV/Å) |

---

## `fe-convert`

Format conversion.

```bash
fe-convert -i input.xyz -o output.pdb
fe-convert -i input.cif -o POSCAR
```

Supported input formats: `.xyz`, `.pdb`, `.cif`, LAMMPS dump  
Supported output formats: `.xyz`, `.pdb`, `POSCAR`, LAMMPS dump

---

## `fe-info`

Print structure summary (cell parameters, atom counts, element list).

```bash
fe-info -i input.xyz
fe-info -i traj.dump --last-n 1
```

---

## `fe-job`

Generate QC software input files.

```bash
fe-job -i input.xyz -s gaussian -m B3LYP -b 6-31G* -o job.gjf
fe-job -i input.xyz -s gromacs -o topology.top
```

| Flag | Description |
|---|---|
| `-s <software>` | Target software: `gaussian`, `gromacs` |
| `-m <method>` | DFT method, e.g. `B3LYP`, `PBE` |
| `-b <basis>` | Basis set, e.g. `6-31G*`, `def2-TZVP` |

---

## `fe-traj`

Structural analysis of MD trajectories.

```bash
fe-traj -m <mode> -i traj.dump [flags] -o output
```

### Modes

#### `gr` — Radial Distribution Function

Computes partial and total $g(r)$ plus coordination numbers $\text{CN}(r)$.

```bash
fe-traj -m gr -i traj.dump --r-max 10.0 --dr 0.01 --r-cut 2.3 -o gr.dat
```

| Flag | Default | Description |
|---|---|---|
| `--r-max` | 10.005 | Maximum radius [Å] |
| `--dr` | 0.01 | Bin width [Å] |
| `--r-cut` | 2.3 | CN integration cutoff [Å] |

Output: `<stem>.dat` (g(r)) and `<stem>_cn.dat` (coordination numbers)

#### `sq` — Structure Factor

Computes $S(q)$ via Fourier transform of $g(r)$.

```bash
fe-traj -m sq -i traj.dump --q-max 25.0 --dq 0.05 --weighting xrd -o sq.dat
```

| Flag | Default | Description |
|---|---|---|
| `--q-max` | 25.0 | Maximum $q$ [Å⁻¹] |
| `--dq` | 0.05 | $q$ bin width [Å⁻¹] |
| `--weighting` | `both` | `xrd`, `neutron`, or `both` |

#### `msd` — Mean Squared Displacement

Computes MSD and directional components; supports NPT and non-periodic trajectories.

```bash
fe-traj -m msd -i traj.dump --dt 2.0 --shift 10 --elements Li -o msd.dat
```

| Flag | Default | Description |
|---|---|---|
| `--dt` | 1.0 | Timestep [fs] |
| `--shift` | 1 | Origin spacing [frames] |
| `--elements` | (all) | Comma-separated element filter |

#### `angle` — Bond Angle Distribution

Computes $P(\theta)$ for all A–B–C triplets.

```bash
fe-traj -m angle -i traj.dump --r-cut-ab 2.3 --r-cut-bc 2.3 --d-angle 0.1 -o angle.dat
```

| Flag | Default | Description |
|---|---|---|
| `--r-cut-ab` | 2.3 | A–B bond cutoff [Å] |
| `--r-cut-bc` | 2.3 | C–B bond cutoff [Å] |
| `--d-angle` | 0.1 | Histogram bin width [°] |

---

## `fe-corr`

Correlation function analysis.

```bash
fe-corr -m <mode> -i traj.dump [flags] -o output
```

### Modes

#### `vacf` — Velocity Autocorrelation Function

Computes VACF and the running Green-Kubo diffusion integral.

```bash
fe-corr -m vacf -i traj.dump --dt 2.0 --elements Li --metal-units -o vacf.dat
```

| Flag | Default | Description |
|---|---|---|
| `--dt` | 1.0 | Timestep [fs] |
| `--shift` | 1 | Origin spacing [frames] |
| `--tau` | (all) | Lag window [frames] |
| `--elements` | (all) | Element filter |

Output columns: `time[fs]`, `vacf[v²]`, `vacf_x`, `vacf_y`, `vacf_z`, `diffusion[v²·fs]`

#### `rotcorr` — Rotational Autocorrelation

Computes $C_2(t)$ for molecular orientation vectors.

```bash
fe-corr -m rotcorr -i traj.dump --center P --neighbor O --r-cut 2.4 --dt 2.0 -o rotcorr.dat
```

| Flag | Default | Description |
|---|---|---|
| `--center` | (required) | Central atom element |
| `--neighbor` | (required) | Neighbour atom element |
| `--r-cut` | 1.2 | Bond search cutoff [Å] |
| `--dt` | 1.0 | Timestep [fs] |
| `--shift` | 1 | Origin spacing [frames] |
| `--tau` | (all) | Lag window [frames] |

Output columns: `time[fs]`, `C(t)`, `integral[fs]`

#### `vanhove` — Van Hove Self-Correlation

Computes $G_s(r, \tau)$ displacement histogram.

```bash
fe-corr -m vanhove -i traj.dump --tau 500 --dt 2.0 --r-max 8.0 --dr 0.02 -o vanhove.dat
```

| Flag | Default | Description |
|---|---|---|
| `--tau` | (last frame) | Lag [frames] |
| `--dt` | 1.0 | Timestep [fs] |
| `--shift` | 1 | Origin spacing [frames] |
| `--r-max` | 10.0 | Max displacement [Å] |
| `--dr` | 0.01 | Bin width [Å] |
| `--elements` | (all) | Element filter |

Output columns: `r[Å]`, `Gs(r,tau)`

---

## `fe-cube`

3-D spatial distribution maps (Gaussian cube format).

```bash
fe-cube -m <mode> -i traj.dump [flags] -o output
```

### Modes

#### `density` — Atomic Number Density

```bash
fe-cube -m density -i traj.dump --nx 80 --ny 80 --nz 80 --elements Li -o li.cube
```

| Flag | Default | Description |
|---|---|---|
| `--nx/ny/nz` | 50 | Grid dimensions |
| `--elements` | (all) | Element filter |

#### `velocity` — Mean Speed per Voxel

Requires trajectory with velocities (use `--metal-units` for LAMMPS metal dumps).

```bash
fe-cube -m velocity -i traj.dump --metal-units -o velocity.cube
```

#### `force` — Mean Force Magnitude per Voxel

Requires trajectory with forces.

```bash
fe-cube -m force -i traj.dump -o force.cube
```

#### `radius` — Hard-Sphere Occupancy

```bash
fe-cube -m radius -i traj.dump --elements Li --radius 0.7 --nx 100 --ny 100 --nz 100 -o li_radius.cube
```

| Flag | Default | Description |
|---|---|---|
| `--radius` | 0.7 | Hard-sphere radius [Å] |
| `--nx/ny/nz` | 50 | Grid dimensions |
| `--elements` | (all) | Element filter |

#### `sdf` — Cluster SDF

```bash
fe-cube -m sdf -i traj.dump --qn 3 --former P --ligand O --cutoff-fl 2.4 \
         --modifier Zn --cutoff-ml 2.8 --grid-res 0.1 --sigma 1.5 -o sdf
```

| Flag | Default | Description |
|---|---|---|
| `--qn` | 3 | Target $Q_n$ level (0–3) |
| `--former` | `P` | Network-former element |
| `--ligand` | `O` | Bridging-ligand element |
| `--cutoff-fl` | 2.4 | Former–ligand cutoff [Å] |
| `--modifier` | (none) | Modifier cation element |
| `--cutoff-ml` | 2.8 | Modifier–ligand cutoff [Å] |
| `--grid-res` | 0.1 | Voxel size [Å] |
| `--sigma` | 1.5 | Gaussian broadening [voxels] |
| `--padding` | 3.0 | Grid margin [Å] |
| `--rmsd-warn` | 0.5 | RMSD warning threshold [Å] |

Output: `<stem>_<label>.cube` per atom type (multiple families: `<stem>_fam<N>_<label>.cube`).

---

## `fe-network`

Glass network analysis: CN distribution, ligand classification (FO/NBO/BO/OBO), Qn speciation, modifier cation roles.

```bash
# 基础用法
fe-network -i traj.dump --P-O=2.3

# 混合形成子 + xlsx 输出
fe-network -i traj.dump --P-O=2.3 --Si-O=1.8 --format xlsx -o result.xlsx

# 含修饰子角色分析
fe-network -i traj.dump --P-O=2.3 --Zn-O=3.5 --modifier Zn

# 多修饰子
fe-network -i traj.dump --P-O=2.3 --Zn-O=3.5 --Na-O=3.2 --modifier Zn,Na
```

### 配对参数（Pair Arguments）

截断参数使用 `--Former-Ligand=cutoff` 格式（首字母大写）：

```
--P-O=2.3     P-O 截断 2.3 Å（P 为形成子，O 为配体）
--Si-O=1.8    Si-O 截断 1.8 Å
--Zn-O=3.5    当配合 --modifier Zn 时，视为修饰子-配体截断
```

### 普通参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-i <file>` | — | 输入轨迹（缺省显示帮助） |
| `-o <file>` | `network` | 输出前缀 |
| `--format csv\|xlsx` | `csv` | 输出格式 |
| `--last-n N` | 全部 | 仅用尾部 N 帧 |
| `--ncore N` | 全部核心 | 线程数 |
| `--metal-units` | 关 | LAMMPS metal 单位 |
| `--modifier Elem` | — | 修饰子元素（逗号分隔） |

### 输出文件

| 文件 | 内容 |
|---|---|
| `<stem>_cn.csv` | CN 分布及均值 |
| `<stem>_ligand.csv` | FO / NBO / BO / OBO 分布 |
| `<stem>_qn.csv` | Qn 物种分布 |
| `<stem>_modifier.csv` | Free / T / B / M 分布（有 `--modifier` 时） |

XLSX 格式将上述内容写入同一文件的多个 sheet（CN、Ligand、Qn、Modifier）。

详细说明见 [Glass Network Analysis](analysis/network.md)。
