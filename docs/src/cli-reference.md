# CLI Reference

nexflux provides six command-line binaries.  All trajectory binaries share common flags for frame selection and parallelism.  Run any binary without `-i` to print mode-specific help.

## Common Flags

| Flag | Description |
|---|---|
| `-i <file>` | Input file |
| `-o <file>` | Output file (default: mode-specific name) |
| `--last-n N` | Use only the last N frames |
| `--ncore N` | Parallel threads (default: all cores) |
| `--metal-units` | LAMMPS metal units (velocities Å/ps, forces eV/Å) |

---

## `nex-convert`

Format conversion.

```bash
nex-convert -i input.xyz -o output.pdb
nex-convert -i input.cif -o POSCAR
```

Supported input formats: `.xyz`, `.pdb`, `.cif`, LAMMPS dump  
Supported output formats: `.xyz`, `.pdb`, `POSCAR`, LAMMPS dump

---

## `nex-info`

Print structure summary (cell parameters, atom counts, element list).

```bash
nex-info -i input.xyz
nex-info -i traj.dump --last-n 1
```

---

## `nex-job`

Generate QC software input files.

```bash
nex-job -i input.xyz -s gaussian -m B3LYP -b 6-31G* -o job.gjf
nex-job -i input.xyz -s gromacs -o topology.top
```

| Flag | Description |
|---|---|
| `-s <software>` | Target software: `gaussian`, `gromacs` |
| `-m <method>` | DFT method, e.g. `B3LYP`, `PBE` |
| `-b <basis>` | Basis set, e.g. `6-31G*`, `def2-TZVP` |

---

## `nex-traj`

Structural analysis of MD trajectories.

```bash
nex-traj -m <mode> -i traj.dump [flags] -o output
```

### Modes

#### `gr` — Radial Distribution Function

Computes partial and total $g(r)$ plus coordination numbers $\text{CN}(r)$.

```bash
nex-traj -m gr -i traj.dump --r-max 10.0 --dr 0.01 --r-cut 2.3 -o gr.dat
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
nex-traj -m sq -i traj.dump --q-max 25.0 --dq 0.05 --weighting xrd -o sq.dat
```

| Flag | Default | Description |
|---|---|---|
| `--q-max` | 25.0 | Maximum $q$ [Å⁻¹] |
| `--dq` | 0.05 | $q$ bin width [Å⁻¹] |
| `--weighting` | `both` | `xrd`, `neutron`, or `both` |

#### `msd` — Mean Squared Displacement

Computes MSD and directional components; supports NPT and non-periodic trajectories.

```bash
nex-traj -m msd -i traj.dump --dt 2.0 --shift 10 --elements Li -o msd.dat
```

| Flag | Default | Description |
|---|---|---|
| `--dt` | 1.0 | Timestep [fs] |
| `--shift` | 1 | Origin spacing [frames] |
| `--elements` | (all) | Comma-separated element filter |

#### `angle` — Bond Angle Distribution

Computes $P(\theta)$ for all A–B–C triplets.

```bash
nex-traj -m angle -i traj.dump --r-cut-ab 2.3 --r-cut-bc 2.3 --d-angle 0.1 -o angle.dat
```

| Flag | Default | Description |
|---|---|---|
| `--r-cut-ab` | 2.3 | A–B bond cutoff [Å] |
| `--r-cut-bc` | 2.3 | C–B bond cutoff [Å] |
| `--d-angle` | 0.1 | Histogram bin width [°] |

---

## `nex-corr`

Correlation function analysis.

```bash
nex-corr -m <mode> -i traj.dump [flags] -o output
```

### Modes

#### `vacf` — Velocity Autocorrelation Function

Computes VACF and the running Green-Kubo diffusion integral.

```bash
nex-corr -m vacf -i traj.dump --dt 2.0 --elements Li --metal-units -o vacf.dat
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
nex-corr -m rotcorr -i traj.dump --center P --neighbor O --r-cut 2.4 --dt 2.0 -o rotcorr.dat
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
nex-corr -m vanhove -i traj.dump --tau 500 --dt 2.0 --r-max 8.0 --dr 0.02 -o vanhove.dat
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

## `nex-cube`

3-D spatial distribution maps (Gaussian cube format).

```bash
nex-cube -m <mode> -i traj.dump [flags] -o output
```

### Modes

#### `density` — Atomic Number Density

```bash
nex-cube -m density -i traj.dump --nx 80 --ny 80 --nz 80 --elements Li -o li.cube
```

| Flag | Default | Description |
|---|---|---|
| `--nx/ny/nz` | 50 | Grid dimensions |
| `--elements` | (all) | Element filter |

#### `velocity` — Mean Speed per Voxel

Requires trajectory with velocities (use `--metal-units` for LAMMPS metal dumps).

```bash
nex-cube -m velocity -i traj.dump --metal-units -o velocity.cube
```

#### `force` — Mean Force Magnitude per Voxel

Requires trajectory with forces.

```bash
nex-cube -m force -i traj.dump -o force.cube
```

#### `radius` — Hard-Sphere Occupancy

```bash
nex-cube -m radius -i traj.dump --elements Li --radius 0.7 --nx 100 --ny 100 --nz 100 -o li_radius.cube
```

| Flag | Default | Description |
|---|---|---|
| `--radius` | 0.7 | Hard-sphere radius [Å] |
| `--nx/ny/nz` | 50 | Grid dimensions |
| `--elements` | (all) | Element filter |

#### `sdf` — Cluster SDF

```bash
nex-cube -m sdf -i traj.dump --qn 3 --former P --ligand O --cutoff-fl 2.4 \
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

## `nex-network`

Glass network analysis: Qn distribution, connectivity statistics.

```bash
nex-network -i traj.dump --P-O 2.3 --format csv -o network.csv
nex-network -i traj.dump --P-O 2.3 --format xlsx -o network.xlsx
```

| Flag | Description |
|---|---|
| `--P-O <r>` | P–O bond cutoff [Å] |
| `--format csv\|xlsx` | Output format |
