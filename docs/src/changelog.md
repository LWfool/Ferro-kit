# Changelog

What changed for *you* — the command lines you have to edit, the products you
have to regenerate. The reasoning behind each decision lives in the development notes,
not here.

Versions from `v0.3.0` onwards are listed. Anything produced by `0.2.x` or
earlier is incompatible with today's ferro in file names, column structure and
command line alike; regenerate it rather than migrating it.

---

## Unreleased

### Readers: fewer silent results

- A file in which no frame is found is an error (`no frame found — is this really a … file?`).
  An empty `.xyz` / `.extxyz` / `.pdb` / `.lammpstrj`, or a file of another format under one of
  those names, used to read as 0 frames with exit status 0.
- XYZ and extended XYZ cut off in their last frame (an interrupted MD run) keep the complete
  frames and drop the last one with a warning naming it, as LAMMPS dumps already did; the whole
  file used to be rejected.
- An atomic number in the first column of a plain XYZ file is read as that element (`8` → O),
  as ASE does; it used to become an element named `8`.
- Extended XYZ: per-atom charges are also read from the `charge` and `initial_charges`
  columns ASE writes (only `charges` was). A `Properties` count that is not a positive integer
  is an error; it used to be taken as 1, shifting every later column.
- PDB: blocks closed by `END` with no `MODEL` records (CP2K and VMD trajectories) are read as
  one frame each; they used to merge into a single frame holding every block's atoms. A last
  `MODEL` with no `ENDMDL` is kept with a warning when it has as many atoms as the frame
  before it, and dropped with a warning otherwise; it used to vanish silently.
  `CRYST1 1 1 1 90 90 90`, the PDB placeholder for a structure with no unit cell, is read as
  non-periodic with a warning; it used to become a periodic 1 Å box.
- CIF: a symmetry operation with a term ferro does not recognise (`2x`, `y*1`) is an error
  naming the term; the term used to be dropped, misplacing the atoms it generated.

### Python: `ferro.FerroError`, and `ferro-py` checks before it reads

- What ferro's Rust libraries report is raised as `ferro.FerroError`.  It is a subclass of
  `RuntimeError`, so `except RuntimeError` keeps working; `except ferro.FerroError` no longer
  swallows Python's own errors.  `read`, `write`, `read_output`, `gr_pair`, `gr_all` and `msd`
  release the GIL while they work.
- `qe_structure_cards` returns `[(card name, text)]` instead of one text holding three cards.
- New: `ferro.interplanar_spacings(cell)` and `ferro.EV_ANG3_TO_KBAR` (CODATA 2018; `conv
  analyze` used the 2014 value, a difference in the ninth digit).
- `ferro-py job` / `conv` report a bad frame selection, `--units` or `--atom-style` before they
  read the first file, and make `-o` only once every input has been read and generated: with
  `--mkdir`, a failed run no longer leaves an empty directory behind.
- `ferro-py job -s qe --timestep` converts with QE's own Rydberg time unit (0.0483776865 fs,
  was 0.048378): `dt` for 2 fs is now 41.3414, was 41.3411.
- A UPF in `pseudo_dir` whose element neither its header nor its name gives is skipped, not a
  traceback; it is named when an element then has no pseudopotential.  A template's `nraise`
  that is not a whole number is named in the error.
- A negative `--start`, `--end`, `--stride` or `--number` to `ferro-py job` / `conv` is a usage
  error (exit 2), not an `OverflowError` traceback.
- The VASP fine grid (`NGXF`…) is predicted with VASP's own ħ²/2mₑ, 3.81002 eV·Å² (`HSQDTM`),
  instead of the CODATA 3.80998; only a cell that sits right on a grid-size boundary changes.

### Readers

- A CP2K single-point output with no `SCF run` line at all (FIST and other runs without SCF) is
  read, as the MD reader already did; before, every frame was dropped as unconverged.
- CP2K input coordinates and cells given in bohr convert with CODATA 2018 (0.5291772109 Å), as
  everywhere else in ferro; they moved in the ninth significant digit.

### One naming rule for every product

Every product is now named `<command>[_<table>][_<selection>][_<input>][_<suffix>].<ext>`.  The
stacked csv files of `traj` and `net` keep their names; the per-input products change, and a
single input now carries its input name too, so the name no longer depends on how many files a
glob matched:

| Command | Before | Now |
|---|---|---|
| `map density` (also velocity, force, radius) | `density.cube`, `density_<s>.cube`, `density_<s>_<input>.cube` | `density_<input>[_<s>].cube` |
| `map sdf` | `sdf_<type>.cube`, `<s>_<type>.cube`, `<s>_<input>_fam<N>_<type>.cube` | `sdf[_fam<N>]_<type>_<input>[_<s>].cube` |
| `map chg-sdf` | `chg_sdf[_<s>][_fam<N>]_Q<n>.cube` | `chg_sdf[_fam<N>]_Q<n>[_<s>].cube` |
| `net --export-traj` | `<input>_types[_<s>].<ext>` | `network_types_<input>[_<s>].<ext>` |
| `bader` | `<input>_ACF[_<s>].dat` | `bader_ACF_<input>[_<s>].dat` (also BCF, AVF) |

`-s` is the batch tag everywhere, so two `dataset` options that used the name differently are
renamed:

- `dataset filter -s` is gone; write `--s-max`.  `filter -s 5` used to mean 5 GPa.
- `dataset merge --suffix` is now `--dir-suffix`: it names the output directory's suffix (`.train`).

Also:

- `map` with every input failing is an error (`every input failed`), as `traj` and `net` already were.
- The `[inputs]` block lists a failed input in its place, not at the end, so it matches the
  `[k/N]` lines on screen.  The `skipped:` line of a failed input goes to stderr.
- `ferro net` prints its label scheme before the `Inputs:` line.
- `dataset collect --type inspect` warns, on screen and in the `#` block of the frame table, when
  an element missing from the element table makes the density too low (it counts as 1 amu), as
  `ferro info` does.  A plain file where a system directory should go is refused with the usual
  `--overwrite` message.
- The help pages share one layout and one block for the common options; `map chg-sdf` showed
  `--qn` defaulting to 2, it is 3.

### Network classification belongs to the analysis library

Nothing changes on the command line except four corrections:

- `ferro net` refuses a trajectory in which a frame has no cell.  The statistics used to leave such
  frames out silently while `--export-traj` refused them; each frame is now classified once and
  feeds both.
- Converting a program output (CP2K, pw.x, OUTCAR, vasprun.xml) to xyz, POSCAR or LAMMPS data no
  longer writes the program name as the structure's title or comment line.
- `dataset collect --type inspect` on a vasprun.xml run now carries the note that its temperature
  column is derived from the kinetic energy.  It was meant to and never did.
- Writing a DeePMD system whose frames carry a stress in a zero-volume cell is an error; it used to
  write a zero virial.  The `--ratio` errors of `dataset filter` / `merge` now begin with
  `Validation error:`, like the other parameter errors.

Python: `Trajectory.charge(frame)` and `multiplicity(frame)` are gone; they always returned 0 and 1.

Library:

- `Frame.charge`, `Frame.multiplicity`, `Frame.bonds`, `TrajectoryMetadata.timestep` and
  `Trajectory::time_at` are removed (nothing produced them).  `TrajectoryMetadata.source` is split
  into `title` and `program`.
- `AtomType`, `TypeParams`, `CutoffTable`, `classify_frame`, the network graph and `NeighborList`
  moved from `ferro_core` to `ferro_analysis`; `PeriodicCell` stays in core.  New
  `calc_network_typed` returns the per-frame types with the statistics.
  `ferro_structure::find_clusters` and `classify_trajectory` are removed.
- `ferro_analysis::ml::merge::sort_atoms` is `Trajectory::canonically_sorted`; the (Z, symbol)
  order is `ferro_core::data::elements::canonical_cmp`.  New `Frame::virial` and
  `Cell::stress_from_virial`.
- `ferro_analysis` and `ferro_io` are reached from the crate root only (`ferro_analysis::GrParams`,
  not `ferro_analysis::md::gr::GrParams`); `ml` stays a module.  The time-correlation analyses
  check the frame spacing themselves, so `check_frame_spacing` is no longer public.
- The `--ratio` split, the stitching of restart segments, the per-frame table of
  `collect --type inspect` and the filter's diagnostic tables are in `ferro_analysis::ml`.

### A third of the memory

A trajectory now stores what an atom is (element, label, mass) once for a run of frames, and only
positions, charges, forces and velocities per frame; the large text readers (LAMMPS dump, extxyz,
CP2K and pw.x logs) no longer hold the whole file as text while parsing, and `map density` /
`map radius` / `traj gr` sum frames in batches instead of keeping one grid per frame.  Reading a
200-frame, 2004-atom dump took 80 MB and takes 19 MB; `map density` on it 468 MB → 71 MB, `traj gr`
227 MB → 57 MB.  Every result is unchanged, byte for byte.

Library: `Frame.atoms`, `Frame.forces` and `Frame.velocities` are no longer public fields.  Read
atoms with `atom(i)` / `atoms()` / `positions()`, forces with `forces()`; build frames with
`add_atom` / `set_atoms`, and set forces with `set_forces`, which refuses a wrong length.

### `convert`, `info` and `traj` read program outputs and DeePMD systems

A VASP `OUTCAR` / `vasprun.xml`, a CP2K or pw.x `.out` log and a DeePMD system directory used to
be readable only by `dataset collect` / `dataset filter`: `ferro info -i vasprun.xml` answered
`Unsupported .xml`.  Every command that reads a structure now takes them, the outputs recognised by
their program banner and read the way `collect` reads them (unconverged frames dropped and
counted).  Extensions are compared case-insensitively (`.XYZ`, `.CIF`), and the Python `read` /
`write` use the same format table as `ferro convert`, so the two can no longer disagree.

### Text that is not UTF-8 is an error at its line

One byte that is not UTF-8 used to reject a whole file with `cannot open FILE` — although it had
opened fine — or `stream did not contain valid UTF-8`, without a position.  The error now names the
file, line, column and byte.

### Reader warnings reach Python

The readers no longer print.  The CLI shows the same messages as before (the LAMMPS dump site-label
mapping no longer ends with a line about `-a/-b/-c`); Python's `read` raises them as `UserWarning`,
where they used to go straight to stderr.  `read_output` also returns `version_note`,
`step_range`, `layout_drift` and `restarts`, which it used to drop — a CP2K version ferro has no
fixture for was flagged on the command line and silent in Python.

Library: `calc_cube_density` and `calc_cube_jump` return `Result` and say why they failed (no
periodic cell, no velocities / forces, too few frames, atoms changing between frames, an empty
element filter) instead of `None`; `SqResult` carries its warnings.

### A box that is not periodic in all three directions is an error

A frame with a cell but `pbc` other than `T T T` — extxyz `pbc="F F F"` or `"T T F"`, a LAMMPS
dump with an `f` or `s` boundary, CP2K `PERIODIC XY` / `NONE`, a DeePMD system with `nopbc` —
used to be measured as a bulk by every analysis: a molecule in a 3 Å box found its neighbours
through images that do not exist, and `traj gr` reported them with exit code 0. Such an input
now fails with `only 3-D periodic cells are supported`, in `traj`, `map`, `net` and the
geometric criteria of `dataset filter`. Results for fully periodic input are unchanged, byte for
byte.

`dataset filter --al6` also checks the cutoff scan of its diagnostics against the
minimum-image bound now, and `map sdf` checks `--cutoff-fl` / `--cutoff-ml`; both used to
measure past it without a word.

### Python `write`: CIF, POSCAR and LAMMPS data get one file per frame

`ferro.write` of a trajectory with several frames to a CIF, POSCAR or LAMMPS data
path used to write only the first frame (POSCAR, data) or all frames into one CIF.
It now writes one file per frame, `conf_0000.cif`, `conf_0001.cif`, …, as
`ferro convert` does. Library callers of `write_cif` / `write_poscar` /
`write_lammps_data` get an error for several frames instead.

### extxyz without a `pos` or a species column is an error

A frame whose `Properties` has no `pos` column used to read every atom at the
origin, and one with no `species` column named every atom `X`. The first is now an
error; for the second, ASE's atomic-number column `Z` gives the element, and a frame
with neither is an error.

### `dataset merge --suffix` must look like `.train`

The suffix is appended to the output directory names, and one containing `/` or `..`
wrote outside `-o`. It must now be `.` followed by letters, digits, `_`, `+` or `-`;
anything else is refused before any input is read.

### `dataset collect`: a VASP output with no usable frame is skipped, not joined

An OUTCAR or vasprun.xml whose every step was dropped (an unfinished first step,
say) used to be read as zero frames and joined to the other files of its directory;
the composition check then failed and the good files were lost with it. It is now
skipped with `no usable frame`, as CP2K and QE outputs already were.

### Site labels survive a LAMMPS dump round trip

Writing a dump now folds `<element>_<suffix>` site labels (`P_3`, `O_b`) into the
element column, as `ferro net --export-traj` always did. A `convert` from a labelled
dump, or from an extxyz with a `label` column, to a dump used to drop every label;
`ferro.write` to a dump keeps them the same way.

### LAMMPS data keeps the masses of its `Masses` section

The mass given for each atom type is now the atoms' mass, not only a hint for the
element: isotopes and coarse-grained masses survive (`ferro info` densities, a
`convert` back to data). Writing a data file whose atoms of one element carry
different masses is an error, since a type has one mass.

### `ferro dataset collect` reads Quantum ESPRESSO pw.x output

`scf`, `md` and `vc-md` outputs are recognised by their `Program PWSCF` banner, or
named with `--format qe/pw` (`ferro doc collect`). An scf needs `tprnfor = .true.`:
frames without forces are dropped, as for CP2K.

### `ferro-py conv -s qe --scan ecutwfc=…`

QE convergence tests scan `ecutwfc` (`ferro-py doc conv`); `conv analyze` reads the
pw.x outputs.

### `ferro-py job -s qe`: `ecutrho` is at least 4 × `ecutwfc`

When `ecutrho` comes from the UPF suggestions and that is below 4 × `ecutwfc`, the job now
writes 4 × `ecutwfc`. Regenerate QE inputs whose template sets a high `ecutwfc` but no
`ecutrho`. An `ecutrho` at or below `ecutwfc` in the template is now an error.

## `v0.4.0` — 2026-10-17

### `dataset filter` / `merge` kept non-periodic systems periodic

A non-periodic DeePMD system (a `nopbc` marker) that also carries `box.npy` — dpdata writes one
for every cluster — lost its `nopbc` marker when `filter` or `merge` wrote it back, so DeePMD-kit
trained it as a periodic system in that box. The marker is now written whenever no frame is
periodic. Regenerate filtered or merged cluster datasets.

### `ferro job` is removed; use `ferro-py job`

`ferro job` and the `ferro-workflow` crate behind it are gone. Calculation inputs come from
`ferro-py job` (VASP, Quantum ESPRESSO and CP2K; single point, optimisation, AIMD), see
`ferro-py doc job`. What `ferro job` did and `ferro-py job` does not: Gaussian input, CP2K
`--task freq`, QE `nscf` / `bands`, CP2K `--atom-charge` / `--cube` / `--molden`, and the
multiplicity guess (`--auto-spin`). `ferro doc job` and `ferro doc spin` are gone with them.

### `convert` no longer writes QE input

`ferro convert -o x.in` (and `.qe`, and `ferro.write(traj, "x.in")` in Python) wrote a pw.x
input with placeholder `&CONTROL` / `&SYSTEM` / `&ELECTRONS`, a `prefix` copied from the
source file name and `K_POINTS gamma` — a file that looked runnable but was not a calculation
setup. QE input is now read-only, like CP2K input: `.in` / `.qe` is refused on write. Generate a
pw.x input with `ferro-py job -s qe`.

### Python: file names mean what they mean on the command line

`ferro.read("x.lammps")` read the file as a dump, while the command line reads `.lammps` as a
LAMMPS data file. Python now follows the command line: `.lammps` is read as data (and needs
`atom_style`) and can be written as data. Likewise `.vasp` / `.pos` are now read and written as POSCAR, and
a `CONTCAR*` name can be written (POSCAR content), as on the command line.

### `.poscar` is a VASP structure

Both the command line and Python now read and write `x.poscar` as POSCAR, alongside `.vasp`,
`.pos` and the `POSCAR*` / `CONTCAR*` names.

### `--type-map` for LAMMPS dumps without an element column

A dump with only `id type x y z` columns read its atoms as `X1`, `X2`, … with no way to name
them. `--type-map 1=Zr,2=P,3=O` (on every command that reads a structure) maps type numbers to
elements. A type missing from the map, an unknown symbol, or a dump that already has an element
column is an error.

### QE input read: one-line namelists, undeclared labels

A namelist written on one line (`&SYSTEM ibrav=0, celldm(1)=10.26, nat=2, ntyp=1 /`, valid
Fortran that pw.x accepts) was ignored as a whole: `celldm(1)` was lost and a `CELL_PARAMETERS`
without an option fell back to Bohr, giving a cell over 10× too small with exit code 0. It is now
read. A missing `&SYSTEM` or `ibrav`, and an `ATOMIC_POSITIONS` label not declared in
`ATOMIC_SPECIES`, are errors.

---

## `v0.3.6` — 2026-10-15

### CIF symmetry expansion: duplicates judged as in ASE

Symmetry images were merged when closer than 0.01 in fractional distance, against every atom of
the same element. Two distinct sites closer than that lost one of them without notice, and images
0.001–0.01 apart were taken for one atom. The rule is now ASE's: images of the **same** site merge
when every fractional coordinate agrees within $10^{-3}$, and two sites landing on one position is an
error naming both.

### `dataset collect` (CP2K MD): the manual's input template prints the cell; missing cells are named

The CP2K MD template in the manual left out `&MOTION/&PRINT/&CELL`, so a run set up from it had
every frame dropped as "incomplete" without saying what was missing. The template now prints the
cell, and frames without a cell line are counted separately ("no cell line") with a pointer to
the fix. A log that also prints `&VELOCITIES` to `__STD_OUT__` is refused: its third block of the
same shape could be taken for the forces.

### `dataset filter`: a system without virial no longer fails at the default settings

`-s` defaults to 10 GPa, so a dataset without `virial.npy` failed under the default
settings, with an error saying a stress threshold had been given. At its default the
criterion is now switched off for such a system, with a `Note`; `-f` likewise for a
system without forces. An explicit `-s` / `-f` still makes a missing label an error.

### Python `add_vacuum_layer`: NaN thickness and singular cells are errors

A `NaN` thickness, or a cell with a zero lattice vector (ASE's $c = 0$ convention for 2D
materials), returned a frame whose cell was all `NaN`. Both are now errors; give a 2D
structure a full 3-D cell first.

### Extended XYZ velocities: GPUMD and ASE files read correctly; written for ASE (breaking)

Reading dropped the velocities of a GPUMD file (`vel`) or a QUIP file (`velo`) without notice,
and took ASE's `momenta` — mass × velocity in ASE's time unit — as a velocity in Å/fs, off by
the mass and a factor of about 10. `vel` and `velo` are now velocities in Å/fs, and `momenta`
is divided by the `masses` column (or the element mass) and converted. Writing used a
`velocities` column that neither ASE nor GPUMD recognises; velocities are now written as
`momenta` with a `masses` column, so ASE's `get_velocities()` returns them, and ferro reads
them back unchanged.

### Inputs whose names differ only in case no longer overwrite each other

Two inputs named `Prod.lammpstrj` and `prod.lammpstrj` (in different directories) kept their
stems, so on macOS and Windows their `map` cubes and `net --export-traj` files overwrote each
other. Stems are now compared ignoring case; the second becomes `prod_2`, with a `Note`.

### Python `merge`: atoms outside their cell no longer land in the other block

Each block is placed by its cell, but atoms were moved as they were: one outside its cell
(a LAMMPS box with `xlo ≠ 0`, or unwrapped coordinates) followed the block into the other
one — with a 2 Å gap two atoms ended up 0.18 Å apart. Atoms are now folded into their own
cell first, the box origin is respected, and an atom outside its cell along a non-periodic
axis is an error.

### LAMMPS data: every section is recognised; bad `Masses` lines are errors

The reader knew ten section names. Any other section after `Masses` (`PairIJ Coeffs`,
`Dihedral Coeffs`, …) was read as more masses: its lines became new atom types, and on a
water file both O and H turned into `X` with no error. All section names of LAMMPS
`read_data` are now recognised, sections are found in any order (a `Masses` written after
`Atoms` used to be skipped), and a `Masses` line that is not exactly `type mass`, or a type
given twice, is an error naming the line.

### `map` cube files: the grid is no longer shifted by half a voxel

Each cube value is a voxel average (or, for `map radius`, a test at the voxel centre), but the
header origin was written as the voxel corner, so a viewer drew the whole grid half a voxel off —
about 0.4 Å per axis at the default 50³ grid on a 40 Å box. The origin is now half a voxel in,
which puts every point at its voxel centre. Applies to `map density/velocity/force/radius/sdf`
and the library's jump map; the grid values themselves are unchanged.

### LAMMPS dump / data keep the box origin

A box whose `xlo ylo zlo` is not 0 used to lose its origin on reading: the writers put the box at
0 while the coordinates kept their values, so after `ferro convert` hundreds of atoms sat outside
the box (225 of 5003 on a sample trajectory with `lo = 1.82` Å). `Cell` now carries an `origin`,
which the LAMMPS readers fill and the LAMMPS writers write back; supercells and added vacuum keep
it. Analyses do not use it, so their results are unchanged.

### `traj sq`: `--r-min` removed; the transform starts at $r = 0$ (breaking)

The Fourier transform integrated only over the $g(r)$ bins, so with `--r-min` above 0 it
dropped the $-\rho$ contribution of $[0, r_{min}]$ and shifted the whole $S(q)$ without any
warning (`--r-min 1.2` moved the neutron-weighted $S(0.7)$ of a ZnO–$P_2O_5$ glass from 0.09
to 0.64). `traj sq` no longer takes `--r-min`, and the library's `calc_sq_from_gr` rejects a
$g(r)$ whose `r_min` is not 0. At the old default of 0.001 Å the change is below
$10^{-3}$. `traj gr` keeps `--r-min`.

### `convert`: CIF is written one frame per file; `-o` can be a directory (breaking)

- Several frames written to `.cif` now give `conf_0000.cif`, `conf_0001.cif`, …
  like POSCAR, instead of one file with several `data_` blocks, which VESTA and
  Materials Studio read as its first frame only.
- `--format EXT` makes `-o` a directory: `-o cifs/ --format cif` writes
  `cifs/frame_0000.cif` …
- `--mkdir` is removed from `convert`; a missing output directory is created
  without asking. The other commands keep `--mkdir`.
- Every format is now written with its atoms sorted by atomic number (Zn, Al,
  P, O → O, Al, P, Zn). Positions are unchanged; only the order is. LAMMPS atom
  ids and type numbers follow the new order, so scripts that pick atoms by id
  or type need checking.

### Smaller changes

- `dataset collect` skips PJM job logs (`x<N>-<job>.<id>.out`) before grouping,
  with one `NOTE` line. They used to count as failed inputs and set exit code 1.
- `ferro --version` also prints the commit, build date, profile, target, rustc
  version and enabled features.
- `traj msd` / `vacf` / `rotcorr` / `vanhove` / `bondlife` without `-i` show the
  help page again, instead of a "missing `--dt`" error.
- `ferro doc` renders display formulas (`$$ … $$`) in the terminal.
- `dataset` help pages no longer print some defaults twice.
- `dataset collect` symmetrises the stress read from vasprun.xml. VASP 6.5.1
  writes it slightly asymmetric, and `filter` / `merge --type nep` then wrote
  a `virial=` that ferro itself refused to read back. Off-diagonal virials
  move by up to about 1 meV.

### Python `merge`: placement fixed for unequal and tilted cells (breaking)

- The narrower block is now centred in the interface plane whichever it is;
  before, only the second one was.
- The second block used to be shifted along its own lattice vectors inside a
  cell built from the first one's. When the two tilts differed, its bottom face
  missed the first block's top face. It now sits exactly `gap` above it.
- `gap` is measured along the interface normal, as in `add_vacuum`. Along a
  tilted join axis it used to be measured along the lattice vector, so the real
  gap came out smaller.
- Interfaces that are not parallel (more than 1° apart) are now an error, and so
  is a non-finite `gap`.

Orthogonal cells of equal in-plane size give the same result as before.

### Python `msd`: `dt` is required (breaking)

`ferro.msd(traj)` used to assume 1 fs between frames and skip the frame-spacing
check that `traj msd` makes. `dt` now has no default, and a trajectory with
repeated or unevenly spaced steps is rejected, as on the command line.

### `dataset filter`: per-system statistics move to `[inputs]` (breaking)

With several systems, the stacked report kept only the first system's header, so
its minimum O–O quantiles and its Al–O cutoff read as if they held for the whole
batch. `filter_min_oo.csv` and `filter_al_cn.csv` no longer carry those lines.
Instead `[inputs]` has one value per system: `min_oo_min`, `min_oo_p1`,
`min_oo_p50`, `min_oo_p99`, `min_oo_max` and `min_oo_mean` with `--oo-min`, and
`al6_rcut` (the cutoff actually used) with `--al6`. The quantile line is now also
printed in read-only mode.

### `convert`: bad arguments no longer leave an empty directory

`--stride 0`, `--number 0`, `--start` past `--end` and an unknown output extension
are now reported before `-o`'s directory is created and before the input is read.
The unknown extension used to be reported only after the whole input was read.

### `-s` accepts only letters, digits, `_`, `+` and `-`

The batch tag goes into the output file name, and `-s 'a/../../x'` used to write
outside `-o`. Other characters are now an error before any input is read, as they
already were for `-a` / `-b` / `--elements`. Applies to `traj`, `map`, `net`,
`map chg-sdf` and `bader`.

### Element-pair cutoffs: `nan` rejected, and only `ferro net` takes them

`--P-O=nan` was accepted and every pair counted as bonded; it is now an error, as
`inf` and non-positive values are. A `--X-Y=value` given to any other subcommand
(`ferro traj gr --P-O=2.3`) was dropped silently; clap now reports it as an
unexpected argument.

### `dataset filter`: frames with NaN forces or stresses

A frame whose forces or stress contain NaN passed the `-f` / `-s` criteria and went
into the training set. With the criterion on, such a frame is now dropped (inf
already was). `-f`, `-s`, `--oo-min` and `--al6` given `nan` or `inf` are now an
error before any input is read; `nan` used to switch the criterion off silently.

### `traj vanhove`: `--tau` defaults to half the trajectory

Without `--tau`, the lag was the last frame, so frame 0 was the only time origin. It is
now $N/2$ frames ($N/2$ origins with the default `--shift 1`), as `--max-lag` in msd /
vacf / rotcorr / bondlife and `gmx velacc -acflen`. Runs without `--tau` give different
numbers; to reproduce an old product pass `--tau` = frames − 1. Runs with `--tau` are
unchanged.

### `traj vanhove` / `traj angle`: `[inputs]` frames and atoms

`[inputs]` listed the number of `r` bins as `frames` (vanhove) and the number of
element types as `atoms` (angle). They are now the frame and atom counts that
`ferro info` reports. Data columns are unchanged.

### `traj gr` / `traj sq`: the header's requested `r_max`

When `--r-max` was clamped to half the smallest interplanar spacing, the header line
`r_max = … (requested; …)` showed the clamped value of the first input instead of the
request. It now shows the request; the clamped value per input stays in `[inputs]`.
Data are unchanged.

### Reading VASP 6.4.2 element lines, CIF, PDB, CP2K and QE input

- POSCAR / CONTCAR / CHGCAR element lines written by VASP 6.4.2 built with HDF5
  (`Na_pv/6a2f546d`, `Cs/`) were taken verbatim as element names. The element is
  now the part before `/` and `_`, as ASE reads it; plain symbols (VASP 5.4.4)
  and count-only lines (VASP 4) read as before.
- A CHGCAR whose coordinate line is not `Direct` (e.g. `Fractional`) was read as
  Cartesian; it now follows VASP's rule (only `C`/`K` mean Cartesian), as POSCAR
  already did.
- A CIF with a block that has no atom sites (`data_global`) failed as a whole;
  such blocks are now skipped.
- PDB: an element in columns 77–78 written in capitals (`FE`, as RCSB files do)
  was taken as the element name `FE`; it is now `Fe`. Blank columns 77–78 gave an
  empty element; the element now comes from the atom name, following the PDB
  column alignment (` CA ` is carbon, `CA  ` is calcium). An `ATOM`/`HETATM` line
  with a bad or missing coordinate was dropped silently; it is now an error naming
  the line. A non-ASCII title no longer crashes the reader.
- CP2K input: section names were matched by prefix, so a `&COLVAR` holding
  `&COORDINATION` before `&COORD` made the read fail (CP2K's own
  `H2O-meta.inp`). An `ABC` inside `&CELL_REF` replaced the `&CELL` one; it is
  now ignored.
- QE input: a namelist written with blanks between assignments
  (`ibrav=0 nat=3 ntyp=1`, legal Fortran) was read as one value, so the `nat`
  check was skipped, a non-zero `ibrav` got through and `starting_magnetization`
  was lost. Blanks and commas now both separate; an `ibrav` or `nat` that does
  not parse is an error.
- **Writing a frame without a cell to `.cif` is now an error** (write `.xyz` or
  `.extxyz`). It used to produce a block ferro itself could not read back.

### `ferro net --export-traj lammpstrj` honours `--units`

The exported dump always wrote velocities and forces in LAMMPS `real` units,
whatever `--units` said: a `metal` dump with `vx = 1.0` came back as `0.001`. The
export now uses the `--units` it was read with, as `convert` does; a trajectory
with velocities or forces and no `--units` is an error naming the flag.
Coordinate-only exports are unchanged.

### `ferro job -s cp2k`: MD ensembles, the CSVR keyword, HSE06

- `--task md` with the default CSVR thermostat wrote `TIMECON_CSVR`, a keyword
  CP2K does not have (it is `TIMECON`), so the input stopped at parsing.
- `--thermostat none` wrote `ENSEMBLE NVT`, which CP2K runs with its default
  Nosé–Hoover thermostat. It now writes `NVE` (`NPE_F` with `--pressure`).
- `--thermostat langevin` put `LANGEVIN` into `&THERMOSTAT TYPE`, which CP2K
  rejects. It now writes `ENSEMBLE LANGEVIN` with `&LANGEVIN GAMMA`. CP2K's
  Langevin ensemble has no barostat, so combining it with `--pressure` is an error.
- `--functional hse06` used the full-range Coulomb operator for its 25% exact
  exchange, a functional that is not HSE06. It now writes
  `&INTERACTION_POTENTIAL` with `POTENTIAL_TYPE SHORTRANGE` and `OMEGA 0.11`.

Regenerate MD and HSE06 inputs made by earlier versions.

### `ferro job -s cp2k`: `--barostat` is replaced by `--pressure`

`--barostat` wrote `PRESSURE 1.01325E+05`, but CP2K reads that keyword in bar:
the run was set to about 10 GPa, not 1 atm. The flag is gone; `--pressure P` gives
the pressure in bar and switches MD to NPT with a flexible cell (`--pressure 1.01325`
is 1 atm). Without it MD stays NVT. Regenerate any NPT input made with
`--barostat`.

### `ferro job` auto-spin: lanthanides and 4d/5d metals

The estimated multiplicity changes for compounds of La–Lu and of the 4d/5d
transition metals. Lanthanides were treated as main-group elements (fluorite
$CeO_2$ came out as multiplicity 9 with `UKS`); they now fill 4f by Hund's rule
($CeO_2$ 1, $Ce_2O_3$ 3, $Gd_2O_3$ 15, $EuO$ 8). 4d/5d ions were estimated
high-spin like 3d ions; they are now low-spin, with $d^8$ square-planar ($K_2PtCl_4$
3 → 1, $RuO_2$ 5 → 3). CP2K and QE inputs generated without `--multiplicity` for
such systems should be regenerated. A new warning names the 4f-in-core case, and
a compound whose oxidation states cannot be balanced (mixed valence such as
$Fe_3O_4$) now says so instead of reporting it as non-ionic.

### `ferro net` cutoffs can be `auto`

`--P-O=auto` takes the cutoff from each input's own g(r) (first minimum behind the
first peak). New; command lines with numbers are unaffected. With `auto` the
`[inputs]` block gains `cutoff` and `g_min` columns and the shared header shows
`P-O=auto`.

### `dataset filter --al6` (bare): the derived cutoff changes

The first-minimum search now works on a 0.10 Å moving average and, inside an empty
gap, takes the gap's midpoint instead of its first bin. The old search could land on
an isolated empty bin in the tail of the peak (2.05 Å instead of 2.57 Å on the
CP2K test trajectory). Datasets filtered with bare `--al6` may keep or drop
different frames; rerun the filter if that matters.

### A CIF without symmetry operations is expanded from its space group

A CIF that named its space group (`_symmetry_space_group_name_H-M 'F m -3 m'`,
an IT number or a Hall symbol) but listed no operations was read as P1: the
asymmetric unit only, e.g. 2 atoms of NaCl instead of 8, with exit code 0. The
symbol is now looked up and expanded. An ambiguous symbol (an origin choice such
as `F d -3 m`) is an error asking for `:1`/`:2`, the Hall symbol, or the
operations. A malformed symmetry operation used to be dropped silently; it is
now an error naming it.

### `traj sq`: no $q = 0$, and a warning below $2\pi/r_{max}$

`--q-min 0` used to give $S(0) = 1$ by fiat, a value disconnected from the
neighbouring points. The transform of a finite box's $g(r)$ has no physical meaning
at $q = 0$, so `--q-min` must now be greater than 0. The `[inputs]` block gains a
`q_trunc` column, $2\pi/r_{max}$, below which the truncation of $g(r)$ dominates
$S(q)$, and a warning is printed when `--q-min` lies below it. With the default
`--r-max 10` and `--q-min 0.1` that warning appears on every run; the data rows are
unchanged.

### `inf` is refused as a range bound; singular cells no longer crash `vanhove` / `map sdf`

`--q-max inf`, `--r-max inf`, `--sigma inf`, `--padding inf` and the like slipped
past the range checks: `traj sq` crashed or wrote a file with only a header and exit
code 0. They are now parameter errors, reported before any file is read, as the
reference already said. A frame whose cell is singular (ASE's convention for 2-D
materials, $c = 0$) crashed `traj vanhove` and `map sdf` with exit code 101 and
aborted the whole batch; it is now an error that skips that input, like every other
command already did.

### Mixed-type DeePMD systems are refused

A system in dpdata's mixed-type layout (`real_atom_types.npy` in its sets) was read
with its placeholder `type.raw`, so every atom became the first element of
`type_map.raw` and `dataset filter` / `merge` wrote that out with exit code 0. Such
a system is now an error naming the set. Convert it to the standard `deepmd/npy`
layout first.

### LAMMPS output: boxes for cell-less structures, orthogonal boxes, boundary flags

Three fixes to the dump and data writers. A structure without a cell used to get a
box of the bounding box's size placed at the origin while the coordinates stayed
put, so atoms sat outside it (and a planar molecule got a box of zero thickness).
The box is now the bounding box widened by 1 Å on every side, around the unchanged
coordinates. An orthogonal cell built from angles was written as triclinic because
$\cos 90^\circ$ is not exactly 0 in floating point; tilts below $10^{-10}$ of the
longest edge are now written as 0, which drops the tilt line from such boxes and
leaves atom lines unchanged. A dump's boundary flags were always `pp pp pp`; they
now follow the frame's periodicity, so a slab reads back as a slab.

### extxyz comment lines are split the way ASE splits them

A bare key without `=` (ASE reads it as `True`) used to glue itself to the key after
it, so a comment line such as `energy=-1.5 is_relaxed Lattice="…" Properties=…`
silently lost its `Lattice` and `Properties`. Bare keys are now read as `T`. Values
in `[…]` or `{…}` and backslash-escaped quotes are understood as well. An
unclosed quote or bracket used to crash ferro and is now an error naming the
frame and line.

### PDB output: one `CRYST1` per model, coordinates in the cell's standard orientation

A multi-frame `.pdb` carried only the first frame's cell, so every frame of an NPT
trajectory read back with frame 0's box. Each `MODEL` is now preceded by its own
`CRYST1`. Because `CRYST1` holds only $a, b, c, \alpha, \beta, \gamma$, a reader
rebuilds the cell with $a$ along $x$ and $b$ in the $xy$ plane; coordinates of a cell
in any other orientation (a POSCAR, for instance) used to be written unrotated and
no longer matched the cell. They are now rotated into that orientation, as ASE does.
LAMMPS cells are already in it, so their coordinates are unchanged. Regenerate
`.pdb` files written from NPT trajectories or non-LAMMPS cells.

## `v0.3.4` — 2026-10-03

This release carries everything since `v0.3.1`: the `0.3.2` and `0.3.3` batches,
which were never released on their own, and the fixes of 2026-10-03. All three
contain breaking changes; read the whole section before upgrading scripts.

### `--metal-units` → `--units real|metal`, required for dumps with velocities/forces

A LAMMPS dump does not record its `units`, and real and metal differ by $10^3$ in
velocity and ~23 in force. The old default (real) silently mis-scaled every
metal-unit dump read without `--metal-units` — which is every DeePMD run.
`--metal-units` is gone; a dump with `vx`/`fx`-type columns now needs
`--units real` or `--units metal`, on every command, or the read fails. Dumps with
positions only are unaffected. Python: `metal_units=True` → `units="metal"`.

| before | now |
|---|---|
| `ferro traj vacf -i t.dump --metal-units` | `ferro traj vacf -i t.dump --units metal` |
| `ferro traj gr -i t.dump …` (dump has `vx vy vz`) | add `--units real` or `--units metal` |

### `--dt` is required for the time-correlation commands; uneven frames are refused

`traj msd`, `vacf`, `rotcorr`, `bondlife` and `vanhove` no longer default `--dt`
to 1.0 fs. It is the time between **stored** frames (MD time step × dump
interval); a forgotten `--dt` used to rescale the time axis and $D$ silently by
the dump interval. LAMMPS dumps now keep their `TIMESTEP`, and when step numbers
are present (dump, CP2K, OUTCAR) a file with a duplicated frame (a restart
written twice) or a changed dump interval is skipped with a message naming the
frames. `convert` to a dump writes the original step numbers.

### A CIF with partial occupancy is an error

`_atom_site_occupancy` was ignored: a mixed site became two atoms at the same
position, silently. Such a file is now refused, naming every partial site.
An unknown (`?`) or malformed coordinate used to put the atom at the origin;
it is now an error naming the site.

### Malformed cube headers are errors

A negative voxel count (Å units) used to give an empty grid and a Bader total of
0 e with exit code 0; an orbital cube's orbital list shifted the whole grid;
surplus data was cut off. All three now fail with a message, except a
single-orbital cube, which is read correctly.

### `traj vanhove`: the `gs` column is now `p_r`, a density in Å⁻¹

The old `gs` was the probability per bin, so its values scaled with `--dr`, and
displacements beyond `--r-max` were dropped while the total still claimed to be
1. The column is now `p_r` $= 4\pi r^2 G_s = $ old `gs` / `dr`, directly comparable
with the Gaussian reference, and `[inputs]` gains `outside_fraction` (the
share beyond `--r-max`, warned above 1 %).

### Neutron scattering lengths: Te, Eu, Hf corrected

Checked entry by entry against the NIST table: Te 5.68 → 5.80 fm, Eu 5.30 →
7.22 fm, Hf 7.77 → 7.7 fm. `total_neutron` changes for systems containing them.
B stays at the $^{11}B$ value (6.65 fm) on purpose, now stated in the manual and
the csv header.

### Single-atom species: self-pair $g(r)$ and $S(q)$ are empty, not 0

With one atom of a species, `gr` for its self-pair was 0 everywhere and $S_{AA}(q)$
oscillated wildly. Both are now empty (NaN); the weighted totals treat that pair
as $S_{AA} = 1$.

### CHGCAR: an unreadable density value is an error

A malformed number inside the density grid used to become 0. It now fails with
its position. Fortran's `0.1234-100` (exponent without `E`) is read correctly.

### `add_vacuum_layer`: thickness is the perpendicular gap

On a tilted axis the vacuum used to be added to the vector length, leaving a real
gap of thickness·cos θ. The gap is now exactly `thickness`, as in ASE.

### `map radius`: radius beyond the minimum-image bound is refused

A `--radius` larger than half the cell's smallest interplanar spacing used to
under-count (only the nearest image); it is now an error naming the frame. The
search window follows each frame's cell, so a shrinking NPT box is fully
covered.

### `-o` is always a path

`--outdir` is gone. `-o` names a **directory** for every command whose run
writes several products (the 13 analysis commands, `map chg-sdf`, `bader`,
`dataset`), and the **output file** for `convert` and `job`, which write
exactly one. A batch tag goes to `-s` / `--suffix`.

| before | now |
|---|---|
| `ferro traj gr -i t.dump -a P -b O -o cmp` | `ferro traj gr -i t.dump -a P -b O -s cmp` |
| `ferro traj gr … --outdir scan` | `ferro traj gr … -o scan` |
| `ferro dataset filter -i raw --outdir clean` | `ferro dataset filter -i raw -o clean` |

`-o cmp` does not fail — it now means *write into the directory `cmp/`*. Check
every script that passed a tag to `-o`.

A missing directory is offered for creation on stderr; **a non-interactive run
(script, CI) must pass `--mkdir`** or it exits with an error.

### `ferro bader` writes beside its input

The three reports default to the **input file's own directory**, not the
current one, and are named after the input stem: `run1/CHGCAR` gives
`run1/CHGCAR_ACF.dat`. VASP calls every charge density `CHGCAR`, so the old
default made two runs overwrite each other. `-o` collects them elsewhere, `-s`
tags them.

### `dataset` product names

`collect` products carry `.db`, `filter` and `merge` products carry `.train`
unless `--ratio` splits them into `.train` / `.valid` / `.test`.

| before | now |
|---|---|
| `sets/run1/` | `sets/run1.db/` |
| `clean/md/` | `clean/md.train/` |

`merge` refuses a mix of split suffixes, and `.valid` / `.test` cannot be split
again. `collect`'s `--outdir` became `--output` (`-o` unchanged), and its `-o`
is optional: without it, products land **beside** each AIMD directory.

### One system per directory

`collect` groups by **directory**, not by file: the `.out` files of one
directory are the restart segments of one run and are reassembled into one
system. Two compositions in one directory is an error, not a dropped frame.

Atom order is **no longer** one of those differences. Every file is sorted into
canonical `(Z, symbol)` order as it is read, so two single points of one
material that merely list their atoms differently now collect into one system.
Products therefore have a different atom order than `0.3.1` produced — training
is unaffected (`type_map.raw` is self-describing, and DeePMD-kit re-sorts by
type on load anyway), but the files are not byte-identical.

### `collect --format`

```bash
ferro dataset collect -i 'sp/*.log' --format cp2k/sp
```

`cp2k/md` | `cp2k/sp` | `vasp/outcar` | `vasp/xml`. The banner is still read;
`--format` overrides it and a mismatch prints one `NOTE:` line. A CP2K file
with no `GLOBAL| Run type` line is now an error naming this flag, where it used
to be read as MD.

### LAMMPS triclinic dump boxes

The three box lines of a **triclinic** dump are written as the bounding box
(`xlo_bound`/`xhi_bound`), which is what the format specifies. Orthogonal cells
are byte-identical to before; triclinic ones were written too small, so OVITO,
ASE and LAMMPS `read_dump` all read the wrong cell. **Regenerate triclinic
dumps written by an earlier ferro.**

### extxyz stress sign

`stress=` follows the ASE convention (positive = tension). The nine numbers are
the negative of what `0.3.1` wrote. `virial=` is unchanged (eV, positive =
compression). A file carrying both is cross-checked and rejected if they
disagree; a 6-component Voigt vector is refused.

### `--plot` is gone

`traj gr`, `sq`, `msd` and `angle` no longer take `--plot`; passing it is now an
`unexpected argument` error. Drop the flag and plot the csv with the scripts in
`scripts/` or your own Python — see [Plotting](plotting.md). The csv products are
unchanged byte for byte.

| before | now |
|---|---|
| `ferro traj msd -i t.dump --dt 1.0 --plot` | `ferro traj msd -i t.dump --dt 1.0` |

### `scripts/plot_net.py --partner` is gone

The two-level stacked bars (Qn hue, $m_X$ shading) were hard to read and have been
removed; `plot_net.py qn` now draws the plain Qn distribution only. `ferro net`
still writes `network_qn_partner.csv`, so filter $Q^n(m\mathrm{Al})$ in pandas.

| before | now |
|---|---|
| `python plot_net.py qn network_qn_partner.csv --partner` | `python plot_net.py qn network_qn.csv` |

### `traj vacf` / `rotcorr` average every origin; trapezoidal integrals

Same cause as `msd` below: `--tau` defaulted to the whole trajectory, so only one time origin was used.

| before | now |
|---|---|
| `--shift N`, `--tau N` | `--shift` removed; `--max-lag N` (default half the trajectory) |
| `diffusion` / `integral`: rectangular sum, too large by $C(0)\,dt/6$ (vacf) and $dt/2$ (rotcorr) | trapezoidal, as GROMACS and MDAnalysis |
| rotcorr: frames without a neighbour still counted in the denominator | only pairs valid at both ends are averaged; lags without any are empty |
| — | new column `vacf_norm` = $C_v/C_v(0)$ |
| `[inputs]` `origins` | `max_lag`, `min_origins`; vacf adds `diffusion_end`, `species`; rotcorr adds `valid_fraction` |

`vanhove` keeps `--tau` and `--shift` for now.

### `traj rotcorr --vector bond`, `--legendre 1|2`

New, nothing to migrate: `--vector bond` follows each centre–neighbour bond of the first frame
(`gmx rotacf -d`) — the way to see tetrahedra ($PO_4$, $SiO_4$) rotate, whose summed bonds cancel.
`--legendre 1` gives $C_1$ (column `c1`).  In the default `sum` mode, bonds that cancel exactly are now
treated as no orientation instead of a vector made of rounding noise.

### `traj msd` averages every origin; `--shift` → `--max-lag`; `msd_x/y/z`

MSD numbers change for every input.  Earlier versions averaged a **single** time origin (the CLI had
no `--tau`, so the window was the whole trajectory) and unwrapped NPT boxes in the lattice view.  Now
every lag averages all origins (FFT) and periodic inputs are unwrapped with the TOR scheme; see
[MSD](analysis/msd.md) for every detail.

| before | now |
|---|---|
| `--shift N` | removed — every frame is an origin |
| lag axis = the whole trajectory | `--max-lag`, default half the trajectory |
| `--fit-range` fractions of the whole trajectory | fractions of the lag axis — the same numbers cover half the time span by default |
| columns `msd_a, msd_b, msd_c` (crystal axes, did not sum to the total in a triclinic box) | `msd_x, msd_y, msd_z` (Cartesian, sum to the total) |
| `[inputs]` `origins` | `max_lag`, `min_origins`; `frames` is now the real frame count |
| Python `ferro.msd(t, shift=, tau=)` | `ferro.msd(t, max_lag=)` |

### MSD fit results moved to `[inputs]`

With several inputs the `msd` header used to show the **first** file's atom count, origins, slope, $D$
and $R^2$ as if they held for the whole batch. The header now keeps only shared parameters; per-input
values are columns of the `[inputs]` list, which also gains `t_lo`, `t_hi`, `points`, `slope`,
`intercept` and the new error `d_err`. Data columns are unchanged. Scripts that grepped `D (total)`
from the header must read the `d_ang2_per_fs` column of `[inputs]` instead.

### New, nothing to migrate

`ferro doc <topic>` (the manual, compiled into the binary; rendered with tables and
formulas on a terminal, the plain source when redirected) · VASP `OUTCAR` and
`vasprun.xml` reading · CP2K single-point reading · `dataset collect --type
inspect` (diagnostics, no dataset) · `--type nep|extxyz` and `--ratio 8:1:1`
on `filter` and `merge` · `--qn` on `net` · `scripts/plot_msd.py` · `traj bondlife` (bond lifetimes and bond formation / breaking events) · `plot_net.py bridge` and `linkmap` · `labels.csv` display names and the `<stem>_data.csv` export in every plotting script.

### Help pages are shorter

Every page is now five sections: what it does, the full parameter table, the
output layout, examples, and a pointer to the manual. Everything else moved
here into the manual. No command line changed.

---

## `v0.3.1` — machine-learning datasets (2026-08-26)

`ferro dataset collect | filter | merge`, DeePMD npy reading and writing, and
the CP2K MD reader. **Purely additive** — no existing behaviour changed.

---

## `v0.3.0` — glass network rework (2026-08-20)

Three rounds of breaking changes released together.

### `ferro net`: six tables, renamed and restructured

| before | now |
|---|---|
| `bridge` | `qn` |
| `partner` | `qn_partner` |
| `oxy` | `ligand_type` |
| `cn` | `coordination` |
| `mean` | removed |
| — | `composition` (new) |

`net qn` and `net type` are gone; `net` is a single leaf command. The number in
a non-Qn former's label is now its **coordination** number (`Al_4` = four
coordinated), not its bridge count.

### `ferro net`: `n` counts homopolar bridges only

The Qn convention follows the literature's $Q^n_m$: **`n` counts P–O–P bridges
only**, heteropolar ones go to the `m_<X>` columns of `qn_partner`, and the
total is `n + sum(m)`. **Every number in the `qn` table changed** — one
reference run moved from 40.4 % `P-Q3` to 0.27 %, and `mean_qn` from 2.40 to
0.95. Plots and comparisons built on the old numbers have to be redone.

### `ferro traj`: product names carry the selection

`gr.csv` became `gr_P-O.csv`, and `gr_all.csv` when nothing is selected. All
six commands changed, including the unselected case.

`traj sq` lost `-a` / `-b` / `-x` / `-y`: every pair is always written, because
the partials sum back to the weighted totals and keeping one pair would hide
exactly that. Filter columns in pandas instead.
