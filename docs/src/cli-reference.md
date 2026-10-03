# CLI Reference

**One binary, `ferro`**, with the subcommands grouped by **output** rather than by the crate that
implements them.  Since 0.2.0 the eight original `fe-*` binaries have all been removed, with no
compatibility layer — the output formats changed at the same time, and keeping `fe-traj` would let an old script "succeed" while emitting a csv it cannot parse itself.  Silently bad data is harder to track down than a command that is gone.

```
ferro traj  gr | sq | msd | angle | vacf | rotcorr | vanhove | bondlife   → stacked csv
ferro map   density | velocity | force | radius | sdf | chg-sdf → one .cube per input
ferro net                                                      → six stacked csv
                                                                 + optional labelled trajectory
ferro bader | convert | info | job
ferro doc   <topic>                                            → the manual (compiled into the binary)
```

**The help has three levels**: `ferro` lists the groups; `ferro traj` lists the commands in that group;
`ferro traj gr` (without `-i`) prints that command's parameters, output column structure and examples.
The leaf commands `convert` / `info` / `bader` work the same way: without `-i` they print their own page.

Two kinds of help coexist and each has its use: **the bare command** gives the rich page (format lists,
output structure, warning notes), **`-h`** gives clap's short parameter table.  `job` is the only exception — it takes over `-h` itself.

| Old command (before 0.2.0) | New command |
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

Every command under `traj` / `map` / `net` flattens the same `CommonArgs` (`convert` / `info` /
`job` / `bader` / `dataset` do **not** take it; each has its own parameters):

| Flag | Description |
|---|---|
| `-i <FILE>...` | input files; takes **several values** and expands glob patterns itself (quote them so the shell leaves them alone) |
| `-o <DIR>` | the output goes into this directory, the current directory by default.  When it does not exist ferro **asks first** (`[y/N]`, on stderr); `--mkdir` skips the question |
| `-s <SUFFIX>` | batch marker, appended to the output name: `<command>[_<table>][_<label>]_<suffix>.csv` |
| `--mkdir` | create the `-o` directory without asking.  **Required in non-interactive environments (scripts, CI)**, which otherwise error out |
| `--last-n N` | use only the last N frames (skipping the equilibration stage) |
| `--ncore N` | number of parallel threads (all cores by default) |
| `--units <UNITS>` | `real` (velocity Å/fs, force kcal/mol/Å) or `metal` (velocity Å/ps, force eV/Å): the LAMMPS unit system of a dump's velocity/force columns.  **Required whenever a dump has `vx`/`fx`-type columns** — a dump does not record its `units`, and the two systems differ by $10^3$ in velocity and ~23 in force, so the value is never guessed (DeePMD runs in `metal`).  Coordinates and cells are Å under either system; a dump with positions only reads without it |
| `--atom-style <STYLE>` | `atomic`, `charge` or `full`: the column layout of a LAMMPS data file's `Atoms` section.  **Required whenever an input is `.data` / `.lmp` / `.lammps`**, checked before the first file is read; see the format notes below |

### Batch processing

`-i` always takes several values.  There is **only one code path**: a single input is the N=1 case, not a
special mode — dispatching on the file count would make the shape of the output depend on how many files the glob happened to match that day.

```bash
ferro traj gr -i 'runs/*/prod.lammpstrj' -a P -b O -s scan
```

Each input is analysed independently and the results are stacked into **one** csv with a `file` column.
Inputs with different element sets take the union of the columns; **what is missing stays empty (NaN),
never padded with zeros, never interpolated**.  A failing input is skipped, leaves its reason in the `[inputs]` block of the output, and makes the **exit code 1** (otherwise an `&&` chain in the shell would take a failure inside the batch for success).
A parameter out of range (`--dr 0`, `--r-min 5 --r-max 3`, `--nx 0`, a step or cutoff of `nan` or `inf`,
`--angle-max 200`, …) is not a per-input failure: it stops the run **before the first file is read and
before `-o` is created**, naming the parameter and its allowed range.
Brace expansion `{a,b}` is not supported; leave that to the shell.

Commands whose output is per-input (the cubes of `ferro map`, the trajectory of `ferro net --export-traj`)
are the exception: the input stem has to be part of the file name, otherwise the second input would overwrite the first.  `-o` applies to these two as well.

The `file` column and those file names use the input's **stem**.  When two inputs share a stem
(`runs/700K/prod.lammpstrj`, `runs/900K/prod.lammpstrj`), the second and later become `prod_2`, `prod_3`, … in
input order, and ferro prints a `Note` naming each renamed file.

### Output naming

```
<-o directory>/<command>[_<table>][_<label>]_<suffix>.csv
```

`label` says **what was computed** and is filled in by the type selection; `-s` is the batch marker.  The
label comes before the suffix, so `ls gr_P-O_*` lists the results for one pair across all batches.

| Command | Where the label comes from | Example |
|---|---|---|
| `traj gr` | `-a/-b` or `-x/-y` | `gr_P-O.csv`, `gr_P_3-O_b.csv`, `gr_all.csv` with no selection |
| `traj angle` | `-a/-b/-c` or `-x/-y/-z` | `angle_O-P-O.csv`, `angle_all.csv` with no selection |
| `traj msd` / `vacf` / `vanhove` | `--elements`, **sorted and deduplicated** | `msd_O-P.csv`, `msd_all.csv` with no selection |
| `traj rotcorr` | `--center`-`--neighbor` | `rotcorr_O-H.csv` (both are required, so `all` is never reached) |
| `traj bondlife` | `--center`-`--neighbor` | `bondlife_Si-O.csv` and `bondlife_events_Si-O.csv` |
| `traj sq` | none (see below) | `sq.csv` |
| `ferro net`, `ferro map` | none | `network_qn.csv`, `density.cube` |

Two rules that are easily confused:

- **`gr` / `angle` join the parts in the order you wrote them**, so `-a P -b O` and `-a O -b P` land in two
  files.  That is correct: `g(r)` is symmetric but `CN` is directed, so the two are genuinely different data.
- **`--elements` is sorted before joining**, because it is a set: `O,P` and `P,O` select the same atoms, and
  one set of data must not end up under two file names.

The selected elements/labels become **part of a path**, so the character set `[A-Za-z0-9_+-]` is validated
before the name is assembled, and a violation errors out **before the first file is read**.  Substituting
underscores was rejected — it would let `-a P/2` and `-a P_2` write silently into the same file.

### Type selection (gr / angle)

The two groups are mutually exclusive:

| Flag | Meaning |
|---|---|
| `-a` / `-b` / `-c` | select by `Atom::element` (the element) |
| `-x` / `-y` / `-z` | select by `Atom::label` (the site label) |

The first slot is the centre (for a pair) or end atom A (for a triplet).  **The order matters**: `g(r)` is symmetric but `CN` is directed.

**`traj sq` has no type selection**: `-a/-b` and `-x/-y` have been removed.  The primary output of $S(q)$ is
the two totals, and the partials are a diagnostic decomposition that adds back to the total
($\sum w_{ij}S_{ij} = \mathrm{total}$); keeping only one pair would hide exactly that closure, and selecting
columns in pandas is enough to look at one pair.  Partials resolved by label went with them — the atom
count behind a single site label is usually too small for its partial to show any signal.  `GroupBy::Label` in the library is untouched.

---

## `ferro convert`

Format conversion.  The format on both sides is decided by the **file name**; there is no `--from` / `--to`.

```bash
ferro convert                              # without -i: prints the table below
ferro convert -i input.xyz -o output.pdb
ferro convert -i input.cif -o POSCAR
ferro convert -i traj.lammpstrj -o traj.extxyz --units metal
```

| Format | Recognised by | Read | Write | Frames written |
|---|---|:-:|:-:|---|
| XYZ | `.xyz` | y | y | all |
| extended XYZ | `.extxyz`, or a `.xyz` whose comment line declares `Lattice=` / `Properties=` | y | y | all |
| PDB | `.pdb` | y | y | all (MODEL records) |
| CIF | `.cif` | y | y | all (several data blocks) |
| LAMMPS dump | `.dump` `.lammpstrj` | y | y | all |
| VASP | `.vasp` `.pos`, or a `POSCAR*` / `CONTCAR*` prefix | y | y | **first frame only** |
| LAMMPS data | `.lammps` `.data` `.lmp` | y | y | **first frame only** |
| QE (pw.x) | `.in` `.qe` | y | y | **first frame only** |
| CP2K input | `.inp` | y | — | — |
| CP2K restart | `.restart` | y | — | — |

Things that trip people up:

- **The two CP2K input formats are read-only.**  To generate a CP2K input use `ferro job -s cp2k`,
  which writes a complete calculation setup rather than bare coordinates.
- **"First frame only" is silent**: writing a 500-frame trajectory as a POSCAR gives frame 0, with no error.
- Writing to the name `CONTCAR` produces content in **POSCAR format**.
- **A `.xyz` can be extended XYZ.**  ASE, GPUMD's `train.xyz`, CP2K's EXTXYZ output and ferro's own
  `dataset --type nep|extxyz` all write extended XYZ under `.xyz`, so a `.xyz` whose comment line (line 2)
  declares `Lattice=` or `Properties=` is read with the extended XYZ reader; anything else stays plain XYZ.
  Writing to `.xyz` still produces plain XYZ (no cell); write `.extxyz` to keep the cell.
- **A CIF with partial occupancy is refused.**  A structure holds whole atoms, so a mixed site
  (`Mg 0.5` / `Fe 0.5`) would become two overlapping atoms and a vacancy site (`O 0.9`) a full
  atom.  The error names every site with `_atom_site_occupancy` below 1; order the structure
  first (one species per site, a supercell for vacancies).
  An unknown (`?`) or malformed coordinate is an error naming the site, not an atom at the origin.
- **A CIF without symmetry operations is expanded from its space-group symbol.**  Explicit operations
  (`_space_group_symop_operation_xyz`, its CIF2 form `_space_group_symop.operation_xyz`, or
  `_symmetry_equiv_pos_as_xyz`) are used as given; otherwise the Hall symbol, the H-M symbol or the
  IT number is looked up, and only a file with none of them is read as P1.  When a symbol fits more
  than one setting — an origin choice (`F d -3 m`), or `C m m e`, whose two settings differ by an
  origin shift — the read **fails** rather than guess, because the wrong origin keeps the atom count
  plausible and every position wrong.  Append the setting (`F d -3 m :2`), give
  `_space_group_name_Hall`, or list the operations.  The monoclinic unique axis and the hexagonal vs.
  rhombohedral cell of an R group are taken from the cell parameters.  A symbol that contradicts
  the IT number, or a malformed operation, is an error.
- VASP files often have no extension, so **both the prefix and the extension are recognised**: `POSCAR`,
  `CONTCAR`, `conf.vasp` and `conf.pos` all go through the same reader/writer pair.
- **A LAMMPS data file needs `--atom-style atomic|charge|full`.**  The layout of the `Atoms` section cannot
  be read off the file: the `Atoms # full` comment is optional, and the column count is ambiguous
  (`charge` and `molecular` both have 6 columns).  Ferro guesses neither, and ignores the comment even
  when it is present.  A line whose column count does not match the style (with or without the three
  image flags) is an error, as is an atom count different from the header's `N atoms`.  Image flags
  are applied, `x + ix·a + iy·b + iz·c`, as LAMMPS `read_data` and ASE do.  ferro's own data writer
  uses `full`, so reading its output back takes `--atom-style full`.

**LAMMPS dump with a triclinic cell**: the box lines carry LAMMPS's `*_bound` (the bounding box after
tilting), not `xlo/xhi`.  Before 2026-09-21 ferro treated them as `xlo/xhi` on both read and write, which
was self-consistent but handed OVITO / ASE / LAMMPS `read_dump` a box that was too small.
**Orthorhombic output is unaffected** (with all three tilt factors 0 the two forms are bit-for-bit identical).

**Whether velocities and forces survive depends on both sides supporting them**: `.dump` to `.xyz` silently
drops the velocities, because plain XYZ has nowhere to put them.  Convert to `.extxyz` to keep them.

### Frame selection

```bash
ferro convert -i traj.dump -o sub.extxyz --start 100            # skip the relaxation stage
ferro convert -i traj.dump -o sub.extxyz --start 100 --end 199
ferro convert -i traj.dump -o POSCAR --stride 50                # one every 50 frames
ferro convert -i traj.dump -o conf.lmp --number 20              # 20 frames at even intervals
```

| Flag | Default | Description |
|---|---|---|
| `--start N` | `0` | first frame, **0-based, inclusive** |
| `--end N` | last frame | last frame, **0-based, inclusive** |
| `--stride N` | `1` | take one frame every N within `[start, end]` |
| `--number N` | — | take **N frames at even intervals** within `[start, end]`, both ends included; mutually exclusive with `--stride` |

Three points of semantics to remember:

- **Closed interval, 0-based**, matching the frame numbers `ferro info` prints: when `info` shows the last
  frame as `Frame 4`, `--end 4` reaches it.  (A half-open interval would need `--end 5`, which does not line up with `info`, so it was not adopted.)
- **`--stride` and `--number` cannot be given together**; clap rejects that while parsing.  One is an
  interval and the other a total, so the combination states two intentions at once; silently ignoring one of them would be the worse choice.
- **`--number` always includes both ends.**  The last frame is often the best-equilibrated configuration, and
  a fixed-stride walk would systematically miss it.  When more frames are asked for than exist, each frame is given once; nothing is duplicated to make up the count.

### One file or N files

**The target format decides; there is no switch**:

| Target format | Output |
|---|---|
| holds a trajectory (`.xyz` `.extxyz` `.pdb` `.cif` `.dump`) | **one** multi-frame file |
| holds a single structure (`POSCAR` `.vasp`/`.pos` `.lmp`/`.data` `.in`/`.qe`) | **one file per frame** |

Writing 20 frames as POSCAR can only ever be 20 files, so there is no switch for the user to remember.

The index is inserted **before the extension** and is the **frame index in the original trajectory** (not
"the n-th one selected"), so the output maps straight back onto the trajectory:

```
-o POSCAR    --stride 2   →  POSCAR_0000      POSCAR_0002      POSCAR_0004
-o conf.vasp --number 3   →  conf_0000.vasp   conf_0002.vasp   conf_0004.vasp
```

The index is zero-padded to at least 4 digits so that `ls` sorts by frame order.  Names recognised by a
**prefix**, such as `POSCAR`, can still be read back with an index attached (`POSCAR_0002` still matches
`POSCAR*`).  **When only one frame is selected, one file is written with no index**, whatever the format.

`-i` currently accepts a **single file** only.  Several inputs plus frame selection would make the output
names of different trajectories collide; that needs the input stem in the file name as well, which is a separate piece of work.

**The element column always carries clean element symbols**, whatever `Atom::label` holds.  Only
`ferro net --export-traj` folds the label into the element column of a LAMMPS dump.

**Reading a LAMMPS dump, columns are taken by the names on the `ITEM: ATOMS` line**, in any order.
Coordinates come from the first complete set of `x y z`, `xs ys zs`, `xu yu zu`, `xsu ysu zsu` (ASE's
priority); a dump with none of them is an error rather than atoms at the origin.  Scaled columns are
converted by LAMMPS's own definition $x = x_{lo} + s \cdot L$, so an atom reads the same whichever kind
was dumped (ASE leaves out $x_{lo}$ for scaled columns and so disagrees with its own `x` reading when the
box does not start at 0).  Also read: `id` (atom order), `type`, `element`, `mass`, `q`, `vx vy vz`,
`fx fy fz`; every other column (`mol`, `ix iy iz`, `c_*`, `f_*`, `v_*`, …) is ignored.  A value that
does not parse is an error naming the atom line and column, except on the very last line of the file,
where it is a truncated final frame and is dropped with a warning.

| Flag | Default | Description |
|---|---|---|
| `-i <file>` | (required) | input file (a single one); omit it to print the format table |
| `-o <file>` | (required) | output file, may include directories (`-o out/run1/x.extxyz`).  A missing parent directory is asked about first; `--mkdir` skips the question.  **A trailing `/` is an error** — the target format is inferred from the file name.  When several frames are written the index goes into the file-name part and the path is unchanged |
| `--mkdir` | off | create the parent directory of `-o` without asking |
| `--start` / `--end` / `--stride` / `--number` | see above | frame selection |
| `--units <UNITS>` | — | `real` or `metal`: LAMMPS units of dump velocities/forces, on read and on write.  Required when they are present |
| `--atom-style <STYLE>` | (none) | `atomic` / `charge` / `full`; required when `-i` is a LAMMPS data file |

---

## `ferro info`

Prints a summary of a structure or trajectory: frame count, element composition, cell parameters, volume
and **mass density**.  The readable formats are exactly those of `ferro convert`.

```bash
ferro info                          # without -i: prints this page
ferro info -i input.xyz
ferro info -i traj.lammpstrj
```

Per-frame report — **the first and the last frame only**, not every frame:

| Row | Content |
|---|---|
| `Atoms` | total count plus the per-element composition |
| `Cell` | a b c (Å) and α β γ (°); `none (non-periodic)` for a non-periodic system |
| `Volume` | Å³ |
| `Density` | **g/cm³** = Σ(atomic masses) / cell volume.  Masses are taken from the file when given explicitly, otherwise from the element table |
| `PBC` | the periodicity flag per axis |
| `Energy` / `Forces` / `Velocities` | whether the frame carries them |

Two edge cases for the density:

- **Without a cell the row is not printed at all**, and no `n/a` placeholder is written — no volume means
  no density, and a placeholder reads like a measurement.
- **Unknown elements pull the density down.**  A symbol that is not in the element table (a stray site
  label, or the `X` that a too-short PDB line degenerates into) falls back to 1 amu in `effective_mass()`,
  with no symptom other than a smaller number.  The density row is therefore followed by a warning naming how many atoms and which symbols triggered the fallback:

  ```
  Density: 0.0765 g/cm³
           WARNING: 2 atom(s) not in the element table (Xx×2) counted as 1 amu — the density is too low
  ```

  **When the warning is there, do not use that number.**

A different volume on the first and the last frame means an NPT trajectory, and the density drifts with it.
For the mean ± σ over the whole trajectory, read `# volume = <mean> +/- <std>` from the header of any `ferro traj` output.

Reading a LAMMPS dump that carries site labels prints the element/label split mapping once.

| Flag | Default | Description |
|---|---|---|
| `-i <file>` | (required) | input file; omit it to print this page |
| `--units <UNITS>` | — | `real` or `metal`: LAMMPS units of dump velocities/forces.  Required when they are present |
| `--atom-style <STYLE>` | (none) | `atomic` / `charge` / `full`; required when `-i` is a LAMMPS data file |

---

## `ferro job`

Generates input files for **Gaussian**, **CP2K** and **Quantum ESPRESSO**.  Without `-s` it prints an
overview; with `-s <software>` but no `-i` it prints that software's own help.  For a guided tour see
[Job Builders](workflow/job-builders.md).

```bash
ferro job                                    # overview
ferro job -s cp2k                            # CP2K-specific help
ferro job -i input.xyz -s gaussian -m B3LYP -b 6-31G* -o job.gjf
ferro job -i input.xyz -s cp2k --task geo-opt --functional pbe --dispersion d3bj
ferro job -i Fe2O3.cif -s qe --auto-spin --kpoints 4 4 4 -o pw.in
```

**Only frame 0 of the input is used** (one structure per input file).  Handing it a multi-frame trajectory
prints three `[warn]` lines (the frame count, how many frames were ignored, and a frame-selection command
to copy), but still generates the input for frame 0 only — and frame 0 is often the least relaxed configuration.  To pick a particular frame, or to generate a batch, extract them with `ferro convert` first:

```bash
ferro convert -i traj.dump -o conf.vasp --number 20       # extract 20 configurations
for f in conf_*.vasp; do
  ferro job -i "$f" -s cp2k --task energy -o "${f%.vasp}.inp"
done
```

### Charge / Spin (shared by all three targets)

| Flag | Default | Description |
|---|---|---|
| `--charge` | (from file) | override the total charge of the system (applied before the spin is inferred) |
| `--multiplicity` | (from file) | force the multiplicity 2S+1; highest priority, and it turns auto-spin off |
| `--auto-spin` | off (on by default for cp2k/qe) | infer the multiplicity from the structure |

The inference chain (magmom → oxidation states + Hund's rule → parity bound on the electron count) is
described in [Spin Estimation](workflow/spin.md).

### Gaussian

| Flag | Default | Description |
|---|---|---|
| `-m <method>` | (required) | DFT method, e.g. `B3LYP`, `PBE0` |
| `-b <basis>` | (required) | basis set, e.g. `6-31G*`, `def2-TZVP` |
| `-o <file>` | `job.gjf` | output file |

### CP2K

#### Task and electronic structure

| Flag | Default | Candidates |
|---|---|---|
| `--task` | `energy` | `energy`, `force`, `geo-opt`, `cell-opt`, `md`, `freq` |
| `--functional` | `pbe` | `pbe`, `blyp`, `pbe0`, `b3lyp`, `revpbe`, `pbesol`, `scan`, `r2scan`, `hse06` |
| `--cp2k-basis` | `dzvp-molopt-sr` | `dzvp-molopt-sr`, `tzvp-molopt`, `tzv2p-molopt`, `dzvp-gth`, `tzvp-gth`, `pob-dzvp`, `pob-tzvp` (all-electron), or any custom string |
| `--dispersion` | `none` | `none`, `d3`, `d3bj` |
| `--scf` | `diag` | `diag` (metals / large systems), `ot` (insulators) |
| `--pbc` | (auto) | `xyz`, `z`, `none`; inferred from the cell when omitted |
| `--kpoints` | (none) | three integers, e.g. `--kpoints 2 2 2` |
| `--cutoff` | `400` | plane-wave cutoff [Ry] |
| `--rel-cutoff` | `50` | relative cutoff [Ry] |
| `--smear` | off | enable Fermi–Dirac smearing |

#### Output

| Flag | Default | Candidates |
|---|---|---|
| `--atom-charge` | `none` | `none`, `mulliken`, `hirshfeld`, `hirshfeld-i` |
| `--cube` | `none` | `none`, `density`, `elf`, `hartree` |
| `--molden` | off | export a Molden orbital file |
| `--project` | `ferro` | CP2K project name |

#### MD (`--task md` only)

| Flag | Default | Description |
|---|---|---|
| `--md-steps` | `10000` | number of MD steps |
| `--md-timestep` | `1.0` | time step [fs] |
| `--temperature` | `298.15` | temperature [K] |
| `--thermostat` | `csvr` | `csvr`, `nose`, `langevin`, `none` |
| `--traj-freq` | `100` | trajectory output frequency [steps] |
| `--barostat` | off | enable an NPT barostat |

> Basis-set and pseudopotential names are resolved **per element** from a 2829-entry database (PBE / SCAN /
> all-electron, with a consistent valence-electron count `q`).  `--cp2k-basis` picks the family and the
> element-specific names are filled in automatically.  See [Job Builders](workflow/job-builders.md#precise-basis--pseudopotential-matching).

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
| `--ecutwfc` | `50` | plane-wave cutoff [Ry] |
| `--smearing` | `none` | `none`, `gaussian`, `mp`, `mv`, `fd` (mp/mv for metals) |
| `--kpoints` | (Gamma) | three integers → a Monkhorst-Pack grid |
| `--pseudo-dir` | `./pseudo` | pseudopotential directory (`<El>.UPF`) |
| `--md-steps` | `10000` | number of MD steps (`--qe-task md`/`vc-md`) |
| `--temperature` | `298.15` | target MD temperature [K] |
| `-o <file>` | `pw.in` | output file |

`ibrav = 0`; the cell is written from the structure as `CELL_PARAMETERS angstrom`.  The spin goes through
the shared inference chain → `nspin` / `tot_magnetization`.

---

## `ferro traj`

Seven trajectory analyses sharing one export pipeline: a single long or wide csv.

```bash
ferro traj <command> -i traj.lammpstrj [flags] -o <suffix>
```

### `gr` — radial distribution function

$g(r)$ and the coordination number $\text{CN}(r)$.

```bash
ferro traj gr -i traj.lammpstrj -a P -b O --r-max 10.0 --dr 0.002 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--r-min` | 0.001 | minimum radius [Å] |
| `--r-max` | 10.005 | maximum radius [Å]; clamped to half the smallest **interplanar spacing** (not the shortest edge length) |
| `--dr` | 0.002 | bin width [Å] |

**Long table**: `file, r, center, neighbor, gr, cn`.  The types go into data columns, so trajectories with
different element sets stack directly; omitting `-a/-b` adds rows rather than columns.  `gr` is symmetric
(`A-B` == `B-A`) and `cn` is directed (`CN(A→B)`), a distinction written into the `center`/`neighbor` columns rather than into a footnote.

### `sq` — structure factor

$S(q)$ from the Fourier transform of $g(r)$.

```bash
ferro traj sq -i traj.lammpstrj --q-max 25.0 --dq 0.02 --weighting both -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--q-min` | 0.1 | minimum $q$ [Å⁻¹] |
| `--q-max` | 25.0 | maximum $q$ [Å⁻¹] |
| `--dq` | 0.02 | $q$ bin width [Å⁻¹] |
| `--weighting` | `both` | `none`, `xrd`, `neutron`, `both` |

The `--r-min` / `--r-max` / `--dr` of `gr` apply here too — they set the range of the $g(r)$ being transformed.

**Wide table**: `file, q, total_xrd, total_neutron`, then three columns per pair (`_sq` / `_xrd` /
`_neutron`, canonical half only).  The primary output is the two totals (one $q$ per row); the weighted
partials are a diagnostic decomposition that sums back to the total.

### `msd` — mean squared displacement

```bash
ferro traj msd -i traj.lammpstrj --dt 2.0 --elements Li --fit-range 0.3,0.8 -o run1
```

Every lag is averaged over all time origins (FFT); periodic inputs are unwrapped with the TOR scheme;
columns `file, time, msd, msd_x, msd_y, msd_z` (Cartesian).  The conventions are spelled out in
[MSD](analysis/msd.md).

| Flag | Default | Description |
|---|---|---|
| `--dt` | **required** | time between stored frames [fs] = MD time step × dump interval.  When the trajectory carries step numbers (LAMMPS dump, CP2K, OUTCAR), duplicated or unevenly spaced frames are refused |
| `--max-lag` | N/2 | longest lag [frames], `1..N-1` |
| `--elements` | (all) | comma-separated element filter |
| `--fit-range` | (none) | `FMIN,FMAX` linear-fit window as fractions of the lag axis (`0..max-lag`) → the self-diffusion coefficient D |

Giving `--fit-range` computes $D = \text{slope}/6$, its error `d_err` and $R^2$, prints them, and adds
`t_lo, t_hi, points, slope, intercept, d_ang2_per_fs, d_err, r2` to each input's row in the `[inputs]` list.
Plot with `scripts/plot_msd.py`, which draws the fitted line from those columns.

### `angle` — bond angle distribution

```bash
ferro traj angle -i traj.lammpstrj -a O -b P -c O --r-cut-ab 2.4 --r-cut-bc 2.4 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--r-cut-ab` | 2.3 | cutoff from end A to centre B [Å] — A is the one given by `-a`/`-x` |
| `--r-cut-bc` | 2.3 | cutoff from end C to centre B [Å] — C is the one given by `-c`/`-z` |
| `--angle-min` | 0.0 | lower bound of the histogram [°] |
| `--angle-max` | 180.0 | upper bound of the histogram [°] (inclusive) |
| `--d-angle` | 0.1 | bin width [°] |

Without a triplet the two cutoffs fall back to the canonical (Z, symbol) order; when both ends are the
same type they both take `min(--r-cut-ab, --r-cut-bc)`.  Angles outside the range are **discarded**, not merely hidden.

**Long table**: `file, angle, end_a, center, end_c, count, p`.  Both the integer `count` and the normalised
`p` are kept — the integer histogram is what the bin-by-bin cross-check against `dump2analysis` rests on.
See [Bond Angle Distribution](analysis/angle.md) for details.

### `vacf` — velocity autocorrelation

```bash
ferro traj vacf -i traj.lammpstrj --dt 2.0 --elements Li --units metal -o run1
```

Every lag averages all time origins (FFT); `diffusion` is the trapezoidal Green–Kubo running integral.
Details in [VACF](analysis/vacf.md).

| Flag | Default | Description |
|---|---|---|
| `--dt` | **required** | time between stored frames [fs] = MD time step × dump interval.  When the trajectory carries step numbers (LAMMPS dump, CP2K, OUTCAR), duplicated or unevenly spaced frames are refused |
| `--max-lag` | N/2 | longest lag [frames], `1..N-1` |
| `--elements` | (all) | element filter |

Columns: `file, time, vacf, vacf_norm, vacf_x, vacf_y, vacf_z, diffusion`

### `rotcorr` — rotational correlation

$C_2(t)$ of a molecular orientation vector.

```bash
ferro traj rotcorr -i water.lammpstrj --center O --neighbor H --r-cut 1.2 --dt 2.0 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--center` | (required) | element of the centre atom |
| `--neighbor` | (required) | element of the neighbour atom |
| `--r-cut` | 1.2 | cutoff for the bond search [Å] |
| `--vector` | `sum` | `sum` = the centre's bonds summed, per frame; `bond` = each bond of frame 0 followed by atom identity (`gmx rotacf -d`) — use for tetrahedra |
| `--legendre` | 2 | order of $P_\ell$: 1 or 2 |
| `--dt` | **required** | time between stored frames [fs] = MD time step × dump interval.  When the trajectory carries step numbers (LAMMPS dump, CP2K, OUTCAR), duplicated or unevenly spaced frames are refused |
| `--max-lag` | N/2 | longest lag [frames], `1..N-1` |

Every lag averages all (molecule, origin) pairs valid at both ends (FFT); `integral` is trapezoidal.
Details in [Rotational correlation](analysis/rotcorr.md).

Columns: `file, time, c2, integral` (`c1` with `--legendre 1`)

### `bondlife` — bond lifetimes and bond events

Intermittent $C_I(t)$ and continuous $S_C(t)$ bond correlation functions over all time origins, and
bonds formed / broken per frame.  Details in [Bond lifetimes](analysis/bondlife.md).

```bash
ferro traj bondlife -i glass_water.lammpstrj --center Si --neighbor O --r-bond 2.2 --dt 1000 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--center` | (required) | centre element |
| `--neighbor` | (required) | neighbour element, may equal `--center` |
| `--r-bond` | (required) | a free pair bonds at $r \le$ this [Å] |
| `--r-break` | `--r-bond` | a bond survives while $r \le$ this [Å] |
| `--intermittency` | 0 | fill breaks of at most this many frames ($S_C$, events) |
| `--dt` | **required** | time between stored frames [fs] = MD time step × dump interval.  When the trajectory carries step numbers (LAMMPS dump, CP2K, OUTCAR), duplicated or unevenly spaced frames are refused |
| `--max-lag` | N/2 | longest lag [frames] |

Columns: `bondlife_<C>-<N>.csv` → `file, time, c_int, s_cont`; `bondlife_events_<C>-<N>.csv` →
`file, time, n_bonds, formed, broken`

### `vanhove` — Van Hove self-correlation

```bash
ferro traj vanhove -i traj.lammpstrj --tau 500 --dt 2.0 --r-max 8.0 --dr 0.02 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--tau` | (last frame) | lag [frames] |
| `--dt` | **required** | time between stored frames [fs] = MD time step × dump interval.  When the trajectory carries step numbers (LAMMPS dump, CP2K, OUTCAR), duplicated or unevenly spaced frames are refused |
| `--shift` | 1 | spacing between time origins [frames] |
| `--r-max` | 10.0 | maximum displacement [Å] |
| `--dr` | 0.01 | bin width [Å] |
| `--elements` | (all) | element filter |

Columns: `file, r, gs`

### Plotting

ferro writes csv only; figures are made in Python. See [Plotting](plotting.md).

---

## `ferro map`

3-D spatial distributions in Gaussian cube format.  **One `.cube` per input**, with no stackable table, so the file name always carries the input stem (`density_<stem>.cube`).

```bash
ferro map <command> -i traj.lammpstrj [flags] -o <suffix>
```

### `density` — atom number density

```bash
ferro map density -i traj.lammpstrj --nx 80 --ny 80 --nz 80 --elements Li -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--nx/ny/nz` | 50 | grid dimensions |
| `--elements` | (all) | element filter |

### `velocity` — mean speed per voxel

Requires a trajectory carrying velocities (a LAMMPS dump needs `--units real|metal`).

```bash
ferro map velocity -i traj.lammpstrj --units metal -o run1
```

### `force` — mean force magnitude per voxel

Requires a trajectory carrying forces.

```bash
ferro map force -i traj.lammpstrj -o run1
```

### `radius` — hard-sphere occupancy

```bash
ferro map radius -i traj.lammpstrj --elements Li --radius 0.7 --nx 100 --ny 100 --nz 100 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--radius` | 0.7 | hard-sphere radius [Å] |
| `--nx/ny/nz` | 50 | grid dimensions |
| `--elements` | (all) | element filter |

### `sdf` — cluster SDF

```bash
ferro map sdf -i traj.lammpstrj --qn 3 --former P --ligand O --cutoff-fl 2.4 \
    --modifier Zn --cutoff-ml 2.8 --grid-res 0.1 --sigma 1.5 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--qn` | 3 | target $Q_n$ level (0–3) |
| `--former` | `P` | network former element |
| `--ligand` | `O` | bridging ligand element |
| `--cutoff-fl` | 2.4 | network former–ligand cutoff [Å] |
| `--modifier` | (none) | modifier cation element |
| `--cutoff-ml` | 2.8 | modifier–ligand cutoff [Å] |
| `--grid-res` | 0.1 | voxel size [Å] |
| `--sigma` | 1.5 | Gaussian broadening [voxels] |
| `--padding` | 3.0 | grid boundary margin [Å] |
| `--rmsd-warn` | 0.5 | RMSD warning threshold [Å] |

Output: one `<stem>_<label>.cube` per atom type (`<stem>_fam<N>_<label>.cube` when there are several families).

### `chg-sdf` — charge-density SDF

Computes the orientationally averaged electron density around Qn clusters from a set of QE `pp.x`
charge-density cubes.  It does **not** use `-i`; it takes `--cubes` instead.

```bash
ferro map chg-sdf --cubes frame_000.cube frame_001.cube frame_002.cube \
    --qn 2 --former P --ligand O --cutoff-fl 2.4 -o run1
```

| Flag | Default | Description |
|---|---|---|
| `--cubes <files…>` | (required) | QE pp.x cube files, one per frame |
| `--qn` | `3` | target Qn level (0–3) |
| `--former` | `P` | network former element |
| `--ligand` | `O` | bridging ligand element |
| `--cutoff-fl` | `2.4` | network former–ligand cutoff [Å] |
| `--modifier` | (none) | modifier cation element |
| `--cutoff-ml` | `2.8` | modifier–ligand cutoff [Å] |
| `--chg-padding` | `6.0` | sub-grid boundary margin [Å] |
| `--rmsd-warn` | `0.5` | alignment RMSD warning threshold [Å] |

Output: `<stem>_Q<n>.cube` (one file per signature family).  For the algorithm see
[Averaged Charge-Density SDF](analysis/chg-sdf.md).

> `--cubes` is the only mode in `ferro map` that aggregates several inputs into **one** SDF, the opposite
> of the "one input, one output" semantics of the rest of the group.  Splitting it up requires first
> defining an intermediate format that carries the sample count (a weighted average across files does not commute); see `dev/plan.md`.

---

## `ferro net`

Glass network topology: bridging-ligand count (the Qn of P), ligand classification, coordination number and linkage statistics.

```bash
ferro net -i traj.lammpstrj --P-O=2.4
ferro net -i traj.lammpstrj --P-O=2.4 --Al-O=2.4 --Zn-O=2.6 --modifier Zn
ferro net -i 'runs/*/prod.lammpstrj' --P-O=2.4 -o scan
ferro net -i traj.lammpstrj --P-O=2.4 --last-n 500 --export-traj
ferro net -i traj.lammpstrj --Al-O=2.4 --Si-O=2.0 --qn Si,Al
```

### Pair parameters (required, at least one)

Cutoffs use the form `--<Former>-<Ligand>=<cutoff>` (element symbols capitalised).  The element pair lives
in the **parameter name**, which clap cannot model, so `main` strips them out of argv before parsing.

```
--P-O=2.4     P-O cutoff 2.4 Å (P is the network former, O the ligand)
--Al-O=2.4    one system may have several network formers
--Al-F=2.1    one network former may have several ligand species
--Zn-O=2.6    read as a modifier-ligand cutoff when --modifier Zn is given
```

### Ordinary parameters

| Parameter | Default | Notes |
|---|---|---|
| `-i <FILE>...` | — | input trajectories (shows the help when omitted) |
| `-o <DIR>` | current directory | output directory; ferro asks first when it does not exist, `--mkdir` skips the question |
| `-s <SUFFIX>` | — | batch marker, appended to the output name |
| `--mkdir` | off | create `-o` without asking (required in non-interactive environments) |
| `--last-n N` | all | use only the last N frames |
| `--ncore N` | all cores | number of threads |
| `--units <UNITS>` | — | `real` or `metal`.  The statistics read neither velocities nor forces, but a dump that has them is refused without it |
| `--modifier E,E` | — | elements that count towards **coordination number only**, taking no part in the bridge count or the ligand classification.  Their cutoffs must be given as well, otherwise ferro errors out |
| `--qn E,E` | `B,P,Si` | network formers to report Qn for.  **Replaces** the default list rather than adding to it; naming a non-former, or an element already taken by `--modifier`, errors out |
| `--export-traj [FMT]` | — | also write a labelled trajectory: `lammpstrj` (default) or `extxyz` |

### Output

Six csv files, each with a `file` column.  **The `#` header of every file carries its own column-by-column
description**, which `pandas.read_csv(comment="#")` drops as a whole.

| File | What it holds |
|---|---|
| `network_composition.csv` | **an overview of the structural composition**: `P-Q2` `Al_4` `O_b` `Zn_4`, each as a fraction of its own element (summing to 1 per element) |
| `network_qn.csv` | the Qn distribution, readable as soon as the file is opened |
| `network_qn_partner.csv` | the same, split by partner element, i.e. $Q^n(m\mathrm{Al})$ |
| `network_ligand_type.csv` | the ligand classification; `label` reads as `Al-O_b-P` |
| `network_coordination.csv` | the coordination-number distribution (network formers + modifiers) |
| `network_linkage.csv` | how the bridges connect: the ligand element plus the site state at both ends |

**Qn is reported only for Qn network formers** (`B,P,Si` by default).  Network formers such as Al are
characterised by coordination number and do not appear as rows in the first two files, but they are still
in the `m_Al` column, in `ligand_type` and in `linkage`.  When there is no Qn network former the first two files are not written at all, and the reason is printed to the screen.

The labels come in **two vocabularies**: the distribution tables (`composition` / `qn` / `qn_partner`) use
the **unit** vocabulary `P-Q2`, because they count structural units; `linkage` and the exported trajectory
use the **atom** vocabulary `P_2` / `Al_4`, because a linkage joins atoms and a trajectory label must be
splittable back into an element.  For non-Qn network formers the two coincide (`Al_4`, where the number is
the **coordination number**).  Ligands are `O_f` / `O_n` / `O_b` / `O_t`, modifiers are bare element symbols.  The table is printed once per run for the parameters actually given.  **Downstream `-x/-y` use the atom vocabulary.**

`--export-traj` writes `<input stem>_types[_<suffix>].<ext>` per input.

For the details (the statements behind the numbers, what `sd` means, pandas examples, how the two export
formats differ, and the single-frame restriction on selecting by label) see [Glass Network Analysis](analysis/network.md).

---

## `ferro bader`

Bader charge decomposition from a DFT charge density.  VASP CHGCAR and Gaussian/QE cube are supported.
A cube must be in Bohr (positive voxel counts) and hold one value per grid point; an Å-unit cube
(negative voxel count), a multi-orbital cube, or one whose data count differs from $N_1 N_2 N_3$ is
refused rather than read as a shifted or empty grid.  An orbital cube (negative atom count) with a
single orbital is read, skipping its orbital list.

```bash
ferro bader                                # without -i: prints the methods and the output description
ferro bader -i CHGCAR                      # VASP CHGCAR
ferro bader -i charge.cube                 # Gaussian/QE cube
ferro bader -i CHGCAR --method weight      # the Yu-Trinkle weight method
ferro bader -i CHGCAR --refine 3 --vacval 1e-4
```

| Flag | Default | Description |
|---|---|---|
| `-i <file>` | (required) | input (`.cube` → the cube reader; anything else → the CHGCAR reader) |
| `-m, --method` | `neargrid` | `ongrid` \| `neargrid` \| `offgrid` \| `weight` |
| `-r, --refine` | `-1` | edge refinement: `-1` automatic, `-2` a single pass, `N` for N passes |
| `-v, --vacval` | `1e-3` | vacuum density threshold [e/Å³] |

**There is no `-o`**: the output file names are decided by the stem of the input file (see below).

### Choosing among the four methods

| Method | Notes |
|---|---|
| `neargrid` | the default.  Gradient ascent with an accumulated off-grid correction plus edge refinement; accurate for ordinary cells |
| `ongrid` | the cheapest, ascending steepest between grid points only.  The basin surfaces come out stepped and the charges are systematically a little off |
| `offgrid` | an interpolated gradient: slower, but with no grid bias |
| `weight` | Yu-Trinkle: a grid point is **split by weight** across several basins according to the flux, rather than assigned whole.  **Use it for strongly tilted (non-orthogonal) cells** — the gradient directions that on/near grid rely on have a known approximation error there |

### Output files

Three `.dat` files in Henkelman format, **named after the stem of the input file**:

| File | Content |
|---|---|
| `<input stem>_ACF.dat` | Atomic Charges File — the Bader charge, volume and minimum distance to the surface, per atom |
| `<input stem>_BCF.dat` | Bader Charge File — the charge, volume and coordinates of each Bader volume |
| `<input stem>_AVF.dat` | Atomic Volume File — the atom → Bader volume index mapping |

These three follow the bader format of the Henkelman group and are parsed by external tools, which is why they did not move to csv with the rest of the outputs.

The three reports are **written next to the input file by default** (`-i run1/CHGCAR` →
`run1/CHGCAR_ACF.dat`), which makes this the only command in the repository that does not default to the
current directory.  VASP charge densities are all called `CHGCAR`, and only landing next to the input, in
their own run directories, do they avoid overwriting each other.  `-o <DIR>` collects them elsewhere, and `-s <SUFFIX>` tells two runs with different parameters on the same input apart (`<stem>_ACF_<suffix>.dat`).

---

## `ferro dataset`

A three-step pipeline for machine-learning training sets.  What sets it apart from the other commands:
the output is a **directory** (a DeePMD system is a directory), so `-o` is the output **root directory**;
and it does not take `CommonArgs`.  The output of `filter` and `merge` carries a `.train` suffix by
default, and the three parts `--ratio` splits out are `.train` / `.valid` / `.test`; a name already ending
in one of them does not get another.  The output of `collect` carries `.db` (raw collected data), which is
**not** a split suffix: `filter` / `merge` strip it before naming their own output — `md.db` filters into `md.train`, not `md.db.train`.

```
ferro dataset collect   AIMD output      → DeePMD system directory (--type inspect for diagnostics)
ferro dataset filter    system directory → a filtered system directory
ferro dataset merge     several systems  → merged by composition
```

All three steps read and write the same directory format, and **none of them modifies its own input**.

### `collect` — AIMD output into a dataset

| Flag | Description |
|---|---|
| `-i <FILE>...` | CP2K MD log / CP2K single-point output / VASP OUTCAR / vasprun.xml; globs are supported |
| `-o <DIR>` | output root directory, **optional**; without it the output lands next to each input directory |
| `--type <WHAT>` | `deepmd` (a DeePMD system) \| `inspect` (diagnostics only, no dataset) [deepmd] |
| `--mkdir` | create `-o` without asking (required when there is no terminal) |
| `--overwrite` | allow writing into an existing non-empty directory |

**One system per input directory.**  The `.out` files in one directory are segments of a single run cut
apart by restarts, and are joined back together — this is also the dividing line between collect and
merge: collect joins the fragments of **one run**, merge combines **different runs**.  The same applies to
a batch of single points: one directory holds several `ENERGY_FORCE` outputs of the same composition, each file contributing one frame.

**Which reader is used is decided by the content of the file**, not by its name.  CP2K writes two
completely different layouts under the same `CP2K|` banner, told apart by `GLOBAL| Run type`:
`ENERGY` / `ENERGY_FORCE` go to the single-point reader, everything else to the MD reader.

The directory name is **what is left after the common ancestor is stripped**, nested as it was rather than flattened, and the file stem does not go into the name:

| `-i` | Output |
|---|---|
| `/data/md/*.out` (no `-o`) | `/data/md.db/` |
| `run*/*.out -o sets` | `sets/run1.db/`, `sets/run2.db/` |
| `/s/a/md/x.out /s/b/md/x.out -o sets` | `sets/a/md.db/`, `sets/b/md.db/` |
| `*.out -o sys` (a single directory) | `sys.db/` itself |

A common prefix carries no distinguishing information by definition, so what remains after stripping it is
necessarily unique and a name collision is no longer an error — a collision is the definition of "these belong together".

The output directory name **always carries `.db`** (database), marking it as raw collected data.  It is
**not** a split suffix: `.train`/`.valid`/`.test` say which part of a split something is, `.db` says where
the data came from.  Keeping it out of that table was a deliberate trade-off — the guard that forbids
splitting a `.valid` again reads exactly that table, and `.db` would be rejected along with them.
`filter` / `merge` strip it before naming their own output, so `md.db` filters into `md.train` rather than
`md.db.train`; two stacked suffixes would also leave `merge`'s shared-suffix check (which looks at the last segment only) without a unique answer.

Without `-o` the output is written **next to** the AIMD directory.  The earlier rule that "`-o` is
mandatory" was aimed at a default of `.`, which would scatter npy files into the directory being worked in; a default that follows the input cannot scatter them anywhere else.

**`--type inspect`**: writes diagnostics only and no dataset.  The three files land in
`<AIMD directory>/ferro_inspect/`, and `-o` **errors out** on this path (there is no dataset to place):

| File | Content |
|---|---|
| `<directory name>.lammpstrj` | all frames, for looking at |
| `<directory name>.data` | the **last frame**, for continuing a run (frame 0 is the initial configuration you fed in) |
| `<directory name>_info.csv` | per frame: `frame` `step` `temperature` `energy` `volume` `density` `source`; the run summary is in the `#` header |

The temperature is read directly from CP2K and from OUTCAR; vasprun.xml does not print it, so it is
back-calculated from the ionic kinetic energy as `T = 2·E_kin/(3N·k_B)`, which the `#` header notes.  Missing values are rendered as **empty fields**, never padded with zeros.

Files are ordered by their first `MD| Step number`, and the order within a file is kept.  **Overlapping
frames are not deduplicated** (a restart only reruns the few steps since the checkpoint, and identical
positions and velocities give identical energies and forces), but the step range of every source file is
printed so that this premise stays checkable.  Single points have no step number and fall back to **ordering by file name**; that column then reads `N single point(s)` instead of `steps ?`.

Inconsistent compositions in one directory **error out immediately**, naming both files, rather than being
dropped as bad frames — that is human error, not a problem with the data.  A single out file that fails to
parse is skipped and the system is built from the rest; the skip list is reported again at the end and the exit code is set to 1.

**MD** requires CP2K to print coordinates, forces and stress all to `__STD_OUT__`, which makes one out
file self-contained.  **Single points** need no `&MOTION`: the coordinates and forces come from the two
tables that `PRINT_LEVEL MEDIUM` prints anyway.  The frame anchor is `ENERGY| Total FORCE_EVAL`, so N
single-point outputs concatenated with `cat` read as N frames.

Units are read from the text itself (`[hartree]` / `[bar]`), and an unrecognised one **errors out** rather
than defaulting — `STRESS_UNIT` is a CP2K input keyword, and one version can emit bar, GPa or atm.
Forces are the only quantity with no unit annotation and fall back to a.u.

**CP2K versions**: two generations are implemented — **2023–2024** (`energy [a.u.]:` plus the
`ATOMIC FORCES` table) and **2025–2026** (`energy [hartree]` plus the `FORCES|` block); the coordinate
table, the kind block and `CELL|` have the same shape in both.  2025–2026 is read silently; **2024 and
earlier get one extra `NOTE:` line per file naming the version, and are then extracted as usual** — the
note asks for a spot-check of one frame, it is not a refusal.  The parenthesised energy and the unprefixed stress block of pre-8.1 versions are still recognised, but those versions are out of scope.

**CP2K kind names**: with several kinds per element (`Fe1`/`Fe2`), `element` takes the real element and
the kind name goes into `label`, while `type_map.raw` is still built from elements, so the two merge into
one training type.  The mapping is printed once per system.  The CP2K plugin of dpdata defaults to the opposite (taking the kind name as the element).

**Cross-check against dpdata**: both single points and AIMD were compared item by item against cp2kdata
0.7.4 (CP2K 2025.2).  Across all 2000 AIMD frames the energies, coordinates and forces are bit-for-bit
identical and the stress differs by 1.2e-15 relative; all five single-point quantities agree.  For the full
numbers and the four accompanying notes (where the force residual comes from, the first MD frame, why coordinates and forces cannot be verified through cp2kdata, and cell precision) see "Checked against dpdata" in `ferro doc dataset collect`.

Frames are dropped for three reasons, and they are **always counted**: SCF not converged / a truncated
block (including a force count that does not match the atom count) / a mismatched composition.  A single
point without stress is **not** a dropped frame — a run without `STRESS_TENSOR` is still good data, it merely has no virial; a missing force is always fatal, since `force.npy` is not optional in DeePMD.

Output:

```
<outdir>/<name>.db/
  type.raw          the type index per atom, 0-based
  type_map.raw      element symbols, ordered by (Z, symbol)
  set.000/coord.npy (nframes, natoms*3)  Å
          box.npy   (nframes, 9)         Å, row-major
          energy.npy(nframes, 1)         eV
          force.npy (nframes, natoms*3)  eV/Å
          virial.npy(nframes, 9)         eV = stress × V
```

Everything on disk is **two-dimensional float64**.  dpdata defaults to float32 and ferro does not follow —
this is the head of the pipeline and everything downstream reads it, so precision lost here cannot be recovered.

### `filter` — select frames by quality

| Flag | Default | Description |
|---|---|---|
| `-i <DIR>...` | | system directories, or a directory containing them (`type.raw` is searched for recursively) |
| `-o <DIR>` | | output root directory, rebuilt along the path relative to `-i`; **omitting it means read-only** |
| `-f, --f-max <EV_PER_A>` | 20.0 | drop the frame when the largest force **vector magnitude** in it exceeds this; 0 disables |
| `-s, --s-max <GPA>` | 10.0 | drop the frame when the largest absolute value among its 9 stress components exceeds this; 0 disables |
| `--oo-min [<DMIN>]` | off / 2.0 when bare | drop the frame when its smallest O–O distance falls below this |
| `--al6 [<RCUT>]` | off / automatic when bare | keep only frames containing a 6-coordinate Al; given bare, the cutoff is the outer edge of the first Al–O RDF shell |
| `--start <N>` | 0 | start of the range (the index among the **surviving frames**, 0-based and inclusive) |
| `--end <N>` | last frame | end of the range (0-based, **inclusive**) |
| `--stride <N>` | 1 | take one out of every N surviving frames |
| `-N, --number <N>` | | take this many frames at even intervals, both ends included; mutually exclusive with `--stride` |
| `--shuffle` | off | shuffle before writing, **after every criterion and after the frame sampling** |
| `--seed <N>` | 666 | the seed for `--shuffle`; giving it without `--shuffle` errors out |
| `--set-size <N>` | 400 | frames per output set; 0 means no splitting |
| `--overwrite` | | allow writing into an existing non-empty directory |

The funnel (narrowing step by step):

```
all frames → |F|max → |σ|max → min d(O-O) → Al6 → [range / sampling] → shuffle
```

**The range and the sampling act on the index among the surviving frames**, not on the original frame
number — after an unknown number of frames has been dropped, that is the only semantics that still makes sense.  The report always gives the original frame numbers.

A threshold of 0 disables that criterion: an explicit zero states "do not judge", which no small positive number can express.

The report has three tables: `[funnel]` for what survives each step, `[criteria]` for how many frames each
criterion rejects and how many it rejects **exclusively**, and `[overlap]` for the pairwise overlap.  **The
exclusive count is the evidence that a criterion is worth having** — each funnel step only counts against the survivors of the previous one, so a criterion that merely re-catches what others already caught still looks busy there.

There are four diagnostic tables as well: the distribution of min d(O–O), the number of Al6 per frame, the
Al coordination distribution, and an **rcut sensitivity scan**.  The last one matters most — on one system
it can climb steeply from 0.9% to 41.4%, while on another it is a flat 100% line.  All four are computed
unconditionally: measured on 1110 frames of 302 atoms the wall-clock time is the same as without them, and computing them only in read-only mode would mean they never reach disk.

Printing and writing are separate:

| | Screen | Written to disk |
|---|---|---|
| without `-o` (read-only) | all seven printed | **nothing at all** |
| with `-o` | only the three statistics tables | all seven written |

The seven csv files sit **flat** in the `-o` root (`filter_funnel.csv`, `filter_rcut_scan.csv`, …), go
through the same writer as every other output, and carry their own `#` header and `[inputs]` list.  Several
systems are stacked into one file, with the row label being the **path relative to `-i`** in the `system`
column — in a nested layout `a/md` and `b/md` have the same leaf name and could not be told apart once stacked.

Keeping them flat rather than in a `report/` subdirectory is deliberate: `expand_dirs` only picks up
directories, so a later `merge -i clean/*` filters the flat csv files out by itself, whereas a `report/` would be picked up as a system candidate.

### `merge` — combine by composition

| Flag | Default | Description |
|---|---|---|
| `-i <DIR>...` | | the system directories to combine; globs are supported |
| `-o <DIR>` | | output root directory, one subdirectory per composition |
| `--mode <MODE>` | shuffle | `shuffle` \| `by-source` |
| `--seed <N>` | 666 | the seed for `shuffle`; unused by `by-source` |
| `--set-size <N>` | 400 | frames per output set; 0 means no splitting |
| `--suffix <EXT>` | see right | force the output directory suffix.  By default: inherit the suffix when the group shares one, use `.train` when none has one; **a mixture of different parts errors out** rather than labelling test data as a training set |
| `--overwrite` | | allow writing into an existing non-empty directory |

**Grouping does not look at directory names** — `init.011` says nothing about what is inside.  Grouping is
by the per-atom element sequence, and only identical compositions are combined.  The output directory is
named `<atom count>_<formula>` (`112_Al32O64Zn16`), where the subscripts are the actual counts, not reduced.

The systems in a group may order their atoms differently: on merging they are brought to the canonical
order `(Z, symbol)`, with **the per-atom arrays (coord, force) following the same permutation**, while
quantities independent of atom numbering (box, energy, virial) are carried over unchanged.  DP is invariant under a permutation of atom numbering; this changes the notation, not the physics.

| Mode | Behaviour |
|---|---|
| `shuffle` | concatenate everything of the same composition → shuffle by seed → split by `--set-size` |
| `by-source` | no mixing and no shuffling; each system is split **within itself**, no set spans two systems, and the correspondence is written to `sets_source.txt` |

Both modes spread the remainder evenly: 500 frames split by 400 gives 250+250, not 400+100.

`filter --shuffle` and `merge --mode shuffle` are **one or the other, not a sequence**: shuffle in merge to
get sets that mix several sources, shuffle in filter when the dataset goes straight to the trainer.

---

## `ferro doc`

Every page of this manual is compiled into the binary via `include_str!` (24 pages, 208 KB), so a ferro
installed with `cargo install` carries it too.

```bash
ferro doc                          # list every topic
ferro doc dataset filter           # read one page
ferro doc net > net.md             # redirected output is not paged and is a clean file
```

**Topics carry the same names as the subcommand tree** (`dataset filter`, `traj gr`, `net`), so the
`Full documentation:` line at the end of every help page is the next command to type rather than a path to
go looking for.  Pages that do not correspond to a single command use flat names (`data-model`,
`installation`, `python`, `cli-reference`).  Short forms such as `gr`, `filter` and `network` have aliases.

`convert` / `info` / `bader` have no page of their own — they are **sections of this page**, and `ferro doc`
addresses sections, taking that `##` heading up to the next one of the same level, so it prints a few dozen lines rather than the whole book.

**Redirected or piped, the markdown is printed exactly as written** — `ferro doc net > net.md` is the
source file, byte for byte.  **On a terminal it is rendered first** and then goes through `$PAGER`
(`less -R` by default), the behaviour of `git` and `man`.  A missing pager, or one that will not start,
falls back to printing rather than failing.

Rendering reflows prose to the terminal width (at most 100 columns), draws tables with borders and wraps
their cells to fit, and turns inline formulas such as `$\alpha_1$` into `α₁`.  A sub- or superscript
that has no Unicode form stays as written (`τ_c`, `r_{min}`), and so do display formulas (`$$...$$`) and
any command the renderer does not know.  The width comes from the terminal, then `$COLUMNS`, then 80.

On Windows, and wherever `NO_COLOR` is set, the output is plainer: Windows gets ASCII only (no colour,
`+--+` table borders, formulas left as TeX), because an older console may show neither escape codes nor
box-drawing characters; `NO_COLOR` turns off colour only.

---

## Output format conventions

Apart from `ferro map` (cube) and `ferro bader` (ACF/BCF/AVF), every output is **one csv** with a `#`
comment block above the data holding the shared parameters and the `[inputs]` list (frame count, atom
count, volume and status per input).

That block is for people to read — `pandas.read_csv(comment="#")` drops it, so **anything a script has to
parse is always a column** and is never hidden in a comment.

Numbers are formatted uniformly as `{:.6e}`, and **NaN renders as an empty field** (under a column union, a column an input lacks is empty, not 0).
