# Changelog

What changed for *you* — the command lines you have to edit, the products you
have to regenerate. The reasoning behind each decision lives in the development notes,
not here.

Versions from `v0.3.0` onwards are listed. Anything produced by `0.2.x` or
earlier is incompatible with today's ferro in file names, column structure and
command line alike; regenerate it rather than migrating it.

---

## Unreleased

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
