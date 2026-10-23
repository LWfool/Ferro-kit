# CLI Reference

**One binary, `ferro`**, with the subcommands grouped by **output** rather than by the crate that
implements them.

```
ferro traj  gr | sq | msd | angle | vacf | rotcorr | vanhove | bondlife   → stacked csv
ferro map   density | velocity | force | radius | sdf | chg-sdf → one .cube per input
ferro net                                                      → six stacked csv
                                                                 + optional labelled trajectory
ferro bader                                                    → ACF / BCF / AVF tables (Henkelman's layout)
ferro convert | info                                           → structure files / a summary
ferro build supercell | vacuum | merge | box                   → one structure file per input
ferro dataset collect | filter | merge                         → DeePMD system directories
ferro doc   <topic>                                            → the manual (compiled into the binary)
```

**The help has three levels**: `ferro` lists the groups; `ferro traj` lists the commands in that group;
`ferro traj gr` (without `-i`) prints that command's parameters, output column structure and examples.
The leaf commands `convert` / `info` / `bader` work the same way: without `-i` they print their own page.

Two kinds of help coexist and each has its use: **the bare command** gives the rich page (format lists,
output structure, warning notes), **`-h`** gives clap's short parameter table. Options given without
an input (`ferro traj gr --dr 0.1`) are an error, not a request for the page: the command exits non-zero.

---

## Common Flags

Every command under `traj` / `map` / `net` takes these flags (`convert` / `info` / `bader` / `dataset`
do **not**; each has its own parameters):

| Flag | Description |
|---|---|
| `-i <FILE>...` | input files; takes **several values** and expands glob patterns itself (quote them so the shell leaves them alone) |
| `-o <DIR>` | the output goes into this directory, the current directory by default.  When it does not exist ferro **asks first** (`[y/N]`, on stderr); `--mkdir` skips the question |
| `-s <SUFFIX>` | batch marker, the last segment of every product name: `<command>[_<table>][_<selection>][_<input>]_<suffix>.<ext>` |
| `--mkdir` | create the `-o` directory without asking.  **Required in non-interactive environments (scripts, CI)**, which otherwise error out |
| `--last-n N` | use only the last N frames (skipping the equilibration stage) |
| `--ncore N` | number of parallel threads (all cores by default) |
| `--units <UNITS>` | `real` (velocity Å/fs, force kcal/mol/Å) or `metal` (velocity Å/ps, force eV/Å): the LAMMPS unit system of a dump's velocity/force columns.  **Required whenever a dump has `vx`/`fx`-type columns** — a dump does not record its `units`, and the two systems differ by $10^3$ in velocity and ~23 in force, so the value is never guessed (DeePMD runs in `metal`).  Coordinates and cells are Å under either system; a dump with positions only reads without it |
| `--type-map <MAP>` | LAMMPS dump without an `element` or `typelabel` column: type number → element, e.g. `1=Zr,2=P,3=O`.  Without it such atoms read as `X1`, `X2`, …; a dump that names its elements in either column refuses it |
| `--atom-style <STYLE>` | `atomic`, `charge` or `full`: the column layout of a LAMMPS data file's `Atoms` section.  **Required whenever an input is `.data` / `.lmp` / `.lammps`**, checked before the first file is read; see the format notes below |

### Batch processing

`-i` always takes several values; a single input is simply N=1 and gives output of the same shape.

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
input order, and ferro prints a `FerroNote` naming each renamed file.  Stems are compared ignoring case
(`Prod` and `prod` clash), because the macOS and Windows file systems do.

### Output naming

```
<-o directory>/<command>[_<table>][_<selection>][_<input>][_<suffix>].<ext>
```

`selection` says **what was computed** and is filled in by the type selection; `input` names the input
file and appears only in per-input products (`ferro map`, `net --export-traj`, `bader`, `build`); `-s` is the batch
marker.  The selection comes before the suffix, so `ls gr_P-O_*` lists the results for one pair across all
batches.

| Command | Where the selection comes from | Example |
|---|---|---|
| `traj gr` | `-a/-b` or `-x/-y` | `gr_P-O.csv`, `gr_P_3-O_b.csv`, `gr_all.csv` with no selection |
| `traj angle` | `-a/-b/-c` or `-x/-y/-z` | `angle_O-P-O.csv`, `angle_all.csv` with no selection |
| `traj msd` / `vacf` / `vanhove` | `--elements`, **sorted and deduplicated** | `msd_O-P.csv`, `msd_all.csv` with no selection |
| `traj rotcorr` | `--center`-`--neighbor` | `rotcorr_O-H.csv` (both are required, so `all` is never reached) |
| `traj bondlife` | `--center`-`--neighbor` | `bondlife_Si-O.csv` and `bondlife_events_Si-O.csv` |
| `traj sq` | none (see below) | `sq.csv` |
| `ferro net`, `ferro map` | none | `network_qn.csv`, `density_<input>.cube` |
| `ferro build` | none | `supercell_<input>.vasp`, `merge_<A>_<B>.cif`, `box.vasp` |

Two rules that are easily confused:

- **`gr` / `angle` join the parts in the order you wrote them**, so `-a P -b O` and `-a O -b P` land in two
  files.  That is correct: `g(r)` is symmetric but `CN` is directed, so the two are genuinely different data.
- **`--elements` is sorted before joining**, because it is a set: `O,P` and `P,O` select the same atoms, and
  one set of data must not end up under two file names.

The selected elements/labels become **part of a path**, so the character set `[A-Za-z0-9_+-]` is validated
before the name is assembled, and a violation errors out **before the first file is read**.

### Type selection (gr / angle)

The two groups are mutually exclusive:

| Flag | Meaning |
|---|---|
| `-a` / `-b` / `-c` | select by `Atom::element` (the element) |
| `-x` / `-y` / `-z` | select by `Atom::label` (the site label) |

The first slot is the centre (for a pair) or end atom A (for a triplet).  **The order matters**: `g(r)` is symmetric but `CN` is directed.

**`traj sq` has no type selection**: every pair is always written; see [S(q)](analysis/sq.md#no-type-selection).

---

## Messages

Everything ferro says outside its products is **one line on stderr**, in the layout of Python's warnings:

```
<source>[:<line>]: <category>: <message>
```

```
traj.pdb:57: FerroWarning: last MODEL (frame 2) has no ENDMDL; kept
sets/a.db: FerroNote: no virial; the default stress criterion is off for it
runs/b/prod.lammpstrj: FerroError: trajectory is empty; skipped
ferro: FerroError: 1 input(s) failed; see the messages above
```

| Part | Meaning |
|---|---|
| source | the input path as you typed it (a product's path when the message is about the file being written), or `ferro` for a message about no file in particular |
| line | the line of that file the message is about, when there is one |
| `FerroNote` | worth knowing; nothing looks wrong |
| `FerroWarning` | something may be wrong with the input or the result; read it before using the numbers |
| `FerroError` | failed: one input skipped in a batch, or the whole run stopped (exit code 1) |

stdout carries only the results and progress (`[2/5] runs/b/prod.lammpstrj`, the tables of `info` and
`dataset`), so `> log.txt` keeps the run's report and `2> messages.txt` its messages, and
`grep FerroWarning` finds every warning of a batch.  A message is never split over several lines.
`ferro-py` prints the same layout; from Python the notes and warnings arrive as `ferro.FerroNote` /
`ferro.FerroWarning` (see the Python API page).

---

## `ferro convert`

Format conversion.  `-o` is a **directory** (created if missing) and `--format EXT` names the format to
write; both are required.  Products are named after the input: `traj.lammpstrj` with `--format extxyz`
gives `DIR/traj.extxyz`, and an input without an extension (`POSCAR`) keeps its whole name.
The input format is decided by the **file name**; there is no `--from`.
The one exception: a name that says nothing about its format — an `OUTCAR`, a
`vasprun.xml`, a CP2K or pw.x `.out` log — is recognised by the program banner inside, and a
directory holding `type.raw` is a DeePMD system.

```bash
ferro convert                              # without -i: prints the table below
ferro convert -i input.xyz -o out --format pdb                     # out/input.pdb
ferro convert -i input.cif -o . --format vasp                      # ./input.vasp
ferro convert -i traj.lammpstrj -o out --format extxyz --units metal
ferro convert -i traj.lammpstrj -o cifs --format cif --units metal   # cifs/traj_0000.cif ...
```

| Format | Recognised by | Read | Write | Frames written |
|---|---|:-:|:-:|---|
| XYZ | `.xyz` | y | y | all |
| extended XYZ | `.extxyz`, or a `.xyz` whose comment line declares `Lattice=` / `Properties=` | y | y | all |
| PDB | `.pdb` | y | y | all (MODEL records) |
| LAMMPS dump | `.dump` `.lammpstrj` | y | y | all |
| CIF | `.cif` | y | y | **one per file** |
| VASP | `.vasp` `.poscar` `.pos`, or a `POSCAR*` / `CONTCAR*` prefix | y | y | **one per file** |
| LAMMPS data | `.lammps` `.data` `.lmp` | y | y | **one per file** |
| QE (pw.x) | `.in` `.qe` | y | — | — |
| CP2K input | `.inp` | y | — | — |
| CP2K restart | `.restart` | y | — | — |
| VASP OUTCAR | banner (an `OUTCAR` has no extension) | y | — | — |
| vasprun.xml | banner (`.xml`) | y | — | — |
| CP2K output | banner (`.out`), MD or single point by its run type | y | — | — |
| QE output | banner (`.out`), pw.x only | y | — | — |
| DeePMD | a directory holding `type.raw` | y | — | — |

Extensions are compared **case-insensitively**: `CONF.XYZ` is XYZ, `OUT.CIF` is CIF.

Things that trip people up:

- **Atoms are written in atomic-number order, in every format.**  A file listing Zn, Al, P, O is
  written as all O, then Al, P, Zn, the way a POSCAR groups its species, so one system converted
  from any program comes out in the same order.  Atoms of one element keep their relative order,
  and forces, velocities and charges move with their atoms.  The LAMMPS atom ids and type numbers
  follow the new order (`O` becomes type 1).  There is no switch to keep the input order.
- **The two CP2K input formats and QE input are read-only.**  To generate either input use
  `ferro-py job -s cp2k|qe`, which writes a complete calculation setup rather than bare coordinates.
- **A program output is read as `ferro dataset collect` reads it**: frames whose SCF did not
  converge, and an incomplete last step, are dropped, and the counts are printed as a warning.
  A file with no usable frame is an error.  Program outputs and DeePMD systems are read-only too.
- **Bytes that are not UTF-8 text are an error** naming the file, line, column and byte (a
  binary file, or another program's output mixed into a log through `2>&1`).  Nothing is replaced
  and read past.
- **"One per file" splits the trajectory**: 500 frames of `traj.dump` with `--format cif` give
  `traj_0000.cif` … `traj_0499.cif`, the index being the original frame number.  CIF is split too, although the format
  allows several `data_` blocks in one file: VESTA and Materials Studio show only the first block.
  Thin the frames out with `--stride` / `--number`.
- **A `.xyz` can be extended XYZ.**  ASE, GPUMD's `train.xyz`, CP2K's EXTXYZ output and ferro's own
  `dataset --type nep|extxyz` all write extended XYZ under `.xyz`, so a `.xyz` whose comment line (line 2)
  declares `Lattice=` or `Properties=` is read with the extended XYZ reader; anything else stays plain XYZ.
  Writing to `.xyz` still produces plain XYZ (no cell); write `.extxyz` to keep the cell.
- **CIF.**  Explicit symmetry operations (`_space_group_symop_operation_xyz`, its CIF2 form, or
  `_symmetry_equiv_pos_as_xyz`) are used as given; otherwise the Hall symbol, the H-M symbol or the IT
  number is looked up, and a file with none of them reads as P1.  A symbol that fits more than one setting
  (`F d -3 m`, `C m m e`) is an error rather than a guess: append the setting (`F d -3 m :2`), give
  `_space_group_name_Hall`, or list the operations.  The monoclinic unique axis and the R-group cell come
  from the cell parameters; images of one site within $10^{-3}$ (fractional) are one atom, two sites on
  one position an error.  Partial occupancy (`_atom_site_occupancy` < 1) is refused, naming every such
  site — order the structure first; an unknown (`?`) or malformed coordinate is an error naming the site.
  Blocks without an `_atom_site_` loop (a `data_global` block) are skipped.  Writing a CIF needs a cell;
  write a non-periodic structure as `.xyz` or `.extxyz`.
- VASP files often have no extension, so **both the prefix and the extension are recognised**: `POSCAR`,
  `CONTCAR`, `conf.vasp`, `conf.poscar` and `conf.pos` all go through the same reader/writer pair.
- **A LAMMPS data file needs `--atom-style atomic|charge|full`.**  The layout of the `Atoms` section cannot
  be read off the file: the `Atoms # full` comment is optional, and the column count is ambiguous
  (`charge` and `molecular` both have 6 columns).  Ferro guesses neither, and ignores the comment even
  when it is present.  A line whose column count does not match the style (with or without the three
  image flags) is an error, as is an atom count different from the header's `N atoms`.  Image flags
  are applied, `x + ix·a + iy·b + iz·c`, as LAMMPS `read_data` and ASE do.  ferro's own data writer
  uses `full`, so reading its output back takes `--atom-style full`.

**LAMMPS dump with a triclinic cell**: the box lines carry LAMMPS's `*_bound` (the bounding box after
tilting), not `xlo/xhi`; with all three tilt factors 0 the two forms are identical.

**LAMMPS dump with a general triclinic cell** (`dump_modify triclinic/general yes`, header
`ITEM: BOX BOUNDS abc origin`): the three box lines are the edge vectors $\mathbf{A}$, $\mathbf{B}$,
$\mathbf{C}$, each followed by one component of the origin, and they are read as such; positions,
velocities and forces are taken as written, already in that frame.  ferro writes the restricted form
(`*_bound` and tilts), rotating positions, velocities and forces into LAMMPS's frame together.

**Element names in a LAMMPS dump**: the `element` column is used; without it, the `typelabel` column
(LAMMPS type labels, set with `labelmap`); with both, `element` wins.  Either is split the same way
when it carries site labels (`P_3` → element `P`, label `P_3`).  With neither, atoms read as `X1`,
`X2`, … by type number unless `--type-map` names them.

**LAMMPS output of a structure without a cell** (a molecule from `.xyz`): the box is the atoms' bounding
box widened by 1 Å on every side, and the coordinates are written unchanged, so `xlo` is generally not 0.
The 1 Å has no physical meaning; it only keeps the box non-degenerate (a planar molecule has zero thickness)
and keeps atoms off the box faces.  A dump writes the boundary flags from the frame's periodicity, `pp` or
`ff` per direction, so a cell-less frame gets `ff ff ff` and a slab keeps `pp pp ff`.  A data file has no
boundary flags; set `boundary` in the LAMMPS input.

**Whether velocities and forces survive depends on both sides supporting them**: `.dump` to `.xyz` silently
drops the velocities, because plain XYZ has nowhere to put them.  Convert to `.extxyz` to keep them: they are
written as ASE's `momenta` column with a `masses` column next to it, so ASE's `get_velocities()` returns them.

**Extended XYZ numbers are written fixed-point with 10 decimals** (ASE writes 8): positions, cell,
forces, energies and stresses keep every digit above $10^{-10}$, far below DFT precision.  Values smaller
than that, and the last bits of a float64, do not survive a round trip.

**A charge or magnetic-moment column that only some atoms carry is refused** (extended XYZ, LAMMPS dump
and data): a structure merged from a charged and an uncharged part would otherwise have the missing values
written as 0, which reads back as "measured 0".  `--drop-partial` leaves the column out instead; LAMMPS data,
whose charge column every atom must have, then writes every charge as 0.  A structure in which no atom has a
charge is unaffected: extended XYZ and dumps write no charge column, LAMMPS data writes 0.

### Frame selection

```bash
ferro convert -i traj.dump -o sub --format extxyz --start 100   # skip the relaxation stage
ferro convert -i traj.dump -o sub --format extxyz --start 100 --end 199
ferro convert -i traj.dump -o vasp --format vasp --stride 50    # one every 50 frames
ferro convert -i traj.dump -o lmp --format lmp --number 20      # 20 frames at even intervals
```

| Flag | Default | Description |
|---|---|---|
| `--start N` | `0` | first frame, **0-based, inclusive** |
| `--end N` | last frame | last frame, **0-based, inclusive** |
| `--stride N` | `1` | take one frame every N within `[start, end]` |
| `-N, --number N` | — | take **N frames at even intervals** within `[start, end]`, both ends included; mutually exclusive with `--stride` (the same short name as `dataset filter`) |

Three points of semantics to remember:

- **Closed interval, 0-based**, matching the frame numbers `ferro info` prints: when `info` shows the last
  frame as `Frame 4`, `--end 4` reaches it.
- **`--stride` and `--number` cannot be given together**; clap rejects that while parsing.  One is an
  interval and the other a total, so the combination states two intentions at once; silently ignoring one of them would be the worse choice.
- **`--number` always includes both ends.**  The last frame is often the best-equilibrated configuration, and
  a fixed-stride walk would systematically miss it.  When more frames are asked for than exist, each frame is given once; nothing is duplicated to make up the count.

### One file or N files

**The target format decides; there is no switch**:

| Target format | Output |
|---|---|
| holds a trajectory (`.xyz` `.extxyz` `.pdb` `.dump`/`.lammpstrj`) | **one** multi-frame file |
| holds a single structure (`.cif` `POSCAR` `.vasp`/`.poscar`/`.pos` `.lmp`/`.data`/`.lammps`) | **one file per frame** |

The inputs (QE `.in`/`.qe`, CP2K `.inp`/`.restart`), the program outputs and DeePMD systems are read-only.

Writing 20 frames as POSCAR can only ever be 20 files, so there is no switch for the user to remember.

The index is inserted **before the extension** and is the **frame index in the original trajectory** (not
"the n-th one selected"), so the output maps straight back onto the trajectory:

```
-i traj.dump --format vasp --stride 2   →  traj_0000.vasp   traj_0002.vasp   traj_0004.vasp
-i traj.dump --format cif  --number 3   →  traj_0000.cif    traj_0002.cif    traj_0004.cif
```

The index is zero-padded to at least 4 digits so that `ls` sorts by frame order.  **When only one frame is selected, one file is written with no index**, whatever the format.

`-i` currently accepts a **single file** only.  Several inputs plus frame selection would make the output
names of different trajectories collide.

**The element column carries clean element symbols**, with one exception: the LAMMPS dump writer folds
a label of the form `<element>_<suffix>` (`P_3`, `O_b`) into its element column, and the dump reader
splits it back.  Labels such as `Fe1` are not folded.

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
| `-o <dir>` | (required) | output **directory**; missing directories are created without asking.  Products are named after the input: `DIR/<stem>.<EXT>`, or `DIR/<stem>_0000.<EXT>` … when several frames go to a one-structure format |
| `--format <EXT>` | (required) | format to write, as an extension (`cif`, `vasp`, `extxyz`, …) |
| `--start` / `--end` / `--stride` / `--number` | see above | frame selection |
| `--units <UNITS>` | — | `real` or `metal`: LAMMPS units of dump velocities/forces, on read and on write.  Required when they are present |
| `--type-map <MAP>` | — | LAMMPS dump without an `element` or `typelabel` column: type number → element, e.g. `1=Zr,2=P,3=O`.  Without it such atoms read as `X1`, `X2`, …; a dump that names its elements in either column refuses it |
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
- **An atom with no mass leaves the row out.**  A symbol that is not in the element table (a stray site
  label, or a placeholder such as `X` written by another program) has no mass unless the file gives one.
  Ferro does not guess one (no 1 amu stand-in, no reading `Fe1` as Fe): the density row is left out, the
  rest of the report is printed, an error names how many atoms and which symbols, and the exit code is 1:

  ```
  box.xyz: FerroError: frame 0: no mass for 2 atom(s) (Xx×2): not in the element table and no mass given; density not computed
  ```

  Correct the element names in the input. The same atoms are an error wherever a mass is needed:
  extended XYZ `momenta` on reading, the `masses` / `momenta` columns on writing, and the `Masses`
  section of a LAMMPS data file.

A different volume on the first and the last frame means an NPT trajectory, and the density drifts with it.
For the mean ± σ over the whole trajectory, read `# volume = <mean> +/- <std>` from the header of any `ferro traj` output.

Reading a LAMMPS dump that carries site labels prints the element/label split mapping once.

| Flag | Default | Description |
|---|---|---|
| `-i <file>` | (required) | input file; omit it to print this page |
| `--units <UNITS>` | — | `real` or `metal`: LAMMPS units of dump velocities/forces.  Required when they are present |
| `--type-map <MAP>` | — | LAMMPS dump without an `element` or `typelabel` column: type number → element, e.g. `1=Zr,2=P,3=O`.  Without it such atoms read as `X1`, `X2`, …; a dump that names its elements in either column refuses it |
| `--atom-style <STYLE>` | (none) | `atomic` / `charge` / `full`; required when `-i` is a LAMMPS data file |

---

## `ferro traj`

Trajectory analyses.  Each writes **one** csv for all its inputs (see
[Output format conventions](#output-format-conventions)); figures are drawn by `ferro-py plot`, see
[Plotting](plotting.md).  Parameters, output columns and conventions are on each command's page.

| Command | Computes | Manual |
|---|---|---|
| `gr` | radial distribution function $g(r)$ and coordination number $CN(r)$ | [g(r)](analysis/gr.md) |
| `sq` | structure factor $S(q)$, X-ray and neutron weighted | [S(q)](analysis/sq.md) |
| `msd` | mean squared displacement and the self-diffusion coefficient | [MSD](analysis/msd.md) |
| `angle` | bond angle distribution | [Bond angles](analysis/angle.md) |
| `vacf` | velocity autocorrelation and its Green–Kubo integral | [VACF](analysis/vacf.md) |
| `rotcorr` | rotational correlation $C_1(t)$ / $C_2(t)$ of an orientation vector | [Rotational correlation](analysis/rotcorr.md) |
| `bondlife` | intermittent and continuous bond correlations, bonds formed / broken per frame | [Bond lifetimes](analysis/bondlife.md) |
| `vanhove` | van Hove self-correlation $G_s(r, \tau)$ | [Van Hove](analysis/vanhove.md) |

---

## `ferro map`

3-D spatial distributions in Gaussian cube format.  **One `.cube` per input**, with no stackable table.
Every name carries its input's name, a single input included (`density_<input>.cube`), so that none
overwrites another and the name does not depend on how many files a glob matched.  `-s TAG` comes last:
`density_<input>_TAG.cube`.  All products follow one rule,
`<mode>[_<table>][_<selection>][_<input>][_<suffix>].<ext>`.

| Command | Computes | Manual |
|---|---|---|
| `density` | atom number density per voxel | [3-D density maps](analysis/cube-density.md) |
| `velocity` / `force` | mean speed / mean force magnitude per voxel | [3-D density maps](analysis/cube-density.md) |
| `radius` | hard-sphere occupancy | [Hard-sphere occupancy](analysis/cube-radius.md) |
| `sdf` | cluster spatial distribution function around $Q_n$ clusters | [Cluster SDF](analysis/cube-sdf.md) |
| `chg-sdf` | orientation-averaged charge density around $Q_n$ clusters; takes `--cubes` (QE pp.x cubes) instead of `-i` and writes **one** SDF from all of them | [Charge-density SDF](analysis/chg-sdf.md) |

---

## `ferro net`

Glass network topology: $Q^n$ of the network formers, ligand types, coordination numbers and linkages,
written as six stacked csv; cutoffs are given per pair as `--P-O=2.4` (or `auto`).  See
[Glass network analysis](analysis/network.md).

```bash
ferro net -i traj.lammpstrj --P-O=2.4 --Al-O=2.4 --Zn-O=2.6 --modifier Zn
```

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
| `-o, --output <DIR>` | next to the input | directory for the three reports |
| `-s, --suffix <SUFFIX>` | — | tag appended to the report names: `bader_ACF_<stem>_<suffix>.dat` |
| `--mkdir` | — | create the output directory without asking (required when there is no terminal) |

### Choosing among the four methods

| Method | Notes |
|---|---|
| `neargrid` | the default.  Gradient ascent with an accumulated off-grid correction plus edge refinement; accurate for ordinary cells |
| `ongrid` | the cheapest, ascending steepest between grid points only.  The basin surfaces come out stepped and the charges are systematically a little off |
| `offgrid` | an interpolated gradient: slower, but with no grid bias |
| `weight` | Yu-Trinkle: a grid point is **split by weight** across several basins according to the flux, rather than assigned whole.  **Use it for strongly tilted (non-orthogonal) cells** — the gradient directions that on/near grid rely on have a known approximation error there |

### Output files

Three `.dat` files in Henkelman format, `bader_<kind>_<input stem>[_<suffix>].dat`:

| File | Content |
|---|---|
| `bader_ACF_<input stem>.dat` | Atomic Charges File — the Bader charge, volume and minimum distance to the surface, per atom |
| `bader_BCF_<input stem>.dat` | Bader Charge File — the charge, volume and coordinates of each Bader volume |
| `bader_AVF_<input stem>.dat` | Atomic Volume File — the atom → Bader volume index mapping |

These three follow the bader format of the Henkelman group and are parsed by external tools.

The three reports are **written next to the input file by default** (`-i run1/CHGCAR` →
`run1/bader_ACF_CHGCAR.dat`), which makes this the only command in the repository that does not default to the
current directory.  VASP charge densities are all called `CHGCAR`, and only landing next to the input, in
their own run directories, do they avoid overwriting each other.  `-o <DIR>` collects them elsewhere, and `-s <SUFFIX>` tells two runs with different parameters on the same input apart (`bader_ACF_<stem>_<suffix>.dat`).

---

## `ferro build`

Makes structures: one product per input, `<-o dir>/<op>_<input>[_<suffix>].<ext>`, in the
input's format unless `--format` names another. It does not take the common flags `--last-n` /
`--ncore`.

| Command | Does |
|---|---|
| `supercell` | replicate a cell, given (`--dims`) or sized by `--min-atoms` / `--min-length` |
| `vacuum` | add a vacuum layer along one lattice vector |
| `merge` | stack two structures with a gap, `merge_<A>_<B>` |
| `box` | pack molecules at random into a cubic box of a given density, `box.vasp` |

Parameters and the reasoning behind them: [Building structures](build.md).

## `ferro dataset`

A three-step pipeline for machine-learning training sets.  What sets it apart from the other commands:
the output is a **directory** (a DeePMD system is a directory), so `-o` is the output **root directory**;
and it does not take the common flags.  The output of `filter` and `merge` carries a `.train` suffix by
default, and the three parts `--ratio` splits out are `.train` / `.valid` / `.test`; a name already ending
in one of them does not get another.  The output of `collect` carries `.db` (raw collected data), which is
**not** a split suffix: `filter` / `merge` strip it before naming their own output — `md.db` filters into `md.train`, not `md.db.train`.

All three steps read and write the same directory format, and **none of them modifies its own input**.

| Command | Does | Manual |
|---|---|---|
| `collect` | AIMD output → DeePMD system directory (`--type inspect` for diagnostics only) | [Collecting AIMD output](dataset/collect.md) |
| `filter` | system directory → a filtered system directory | [Filtering a dataset](dataset/filter.md) |
| `merge` | several systems → merged by composition | [Merging datasets](dataset/merge.md) |

---

## `ferro doc`

Every page of this manual is compiled into the binary via `include_str!`, so a ferro installed with
`cargo install` carries it too. The `ferro-py` pages (`python/`) are compiled into ferro-py instead and
read with `ferro-py doc`.

```bash
ferro doc                          # list every topic
ferro doc dataset filter           # read one page
ferro doc net > net.md             # redirected output is not paged and is a clean file
```

**Topics carry the same names as the subcommand tree** (`dataset filter`, `traj gr`, `net`), so the
`Full documentation:` line at the end of every help page is the next command to type rather than a path to
go looking for.  Pages that do not correspond to a single command use flat names (`data-model`,
`installation`, …); `ferro doc` with no argument lists every one.  Short forms such as `gr`, `filter` and `network` have aliases.

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
count, volume and status per input).  Several inputs stack into that one table, told apart by a `file`
column.

That block is for people to read — `pandas.read_csv(comment="#")` drops it, so **anything a script has to
parse is always a column** and is never hidden in a comment.

Numbers are formatted uniformly as `{:.6e}`, and **NaN renders as an empty field** (under a column union, a column an input lacks is empty, not 0).
