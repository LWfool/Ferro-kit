# ferro-py job: Calculation Inputs

`ferro-py job` turns structures into ready-to-run inputs for VASP, Quantum
ESPRESSO (pw.x) and CP2K — a single point (`sp`), a geometry optimisation
(`opt`) or an NVT / NpT molecular dynamics run (`aimd`).

```bash
ferro-py job -s vasp -t sp -i POSCAR -o sp/
ferro-py job -s cp2k -t aimd -i traj.lammpstrj --number 20 --temperature 1500 -o md/
ferro-py job -s qe -t opt -i 'cifs/*.cif' --vary-cell --template my_relax.in -o relax/
```

It is meant for the common calculation, done reliably: a periodic structure,
one of three tasks, sensible defaults. Anything more specialised — fixed atoms,
initial velocities, NEB, phonons, molecules — is a job for the program's own
tools or a hand-edited input.

This page covers what is common to the three programs. What each one writes,
and its limits, is on its own page:

| Program | Files written | Page |
|---|---|---|
| VASP | `INCAR`, `POSCAR`, `POTCAR` | `ferro-py doc job vasp` |
| Quantum ESPRESSO | `pw.in` | `ferro-py doc job qe` |
| CP2K | `cp2k.inp` | `ferro-py doc job cp2k` |

## The template is the input

Every job starts from a **template, which is a native input file of the
program** — an `INCAR`, a `pw.in`, a `cp2k.inp`. ferro-py rewrites a small,
fixed set of quantities in it and leaves every other line as written,
comments included. A setting you know from your own runs therefore goes into
the template exactly as you would write it for the program; ferro-py does not
need to understand it, and keys it has never heard of pass through untouched.

Which template a job uses:

1. `--template FILE` on the command line;
2. otherwise `<dir>/<program>/<task>/<file>` under the configured
   `[templates] dir` (`ferro-py doc config`), when that file exists;
3. otherwise the built-in template.

The nine built-in templates (three programs × three tasks) are generic
PBE-D3(BJ) inputs. To start your own, export one:

```bash
ferro-py job template -s vasp -t aimd -o my_aimd/     # writes my_aimd/INCAR
# edit my_aimd/INCAR, then
ferro-py job -s vasp -t aimd -i POSCAR --template my_aimd/INCAR -o md/
```

The header comment of each template lists what ferro-py rewrites in it and what
is deliberately missing (the parallel settings — they belong to your machine).

## What ferro-py rewrites

| Quantity | Source | Notes |
|---|---|---|
| Structure | the input file | cell and positions; atoms sorted by (Z, symbol), as `ferro convert` does |
| Element list | the structure | POTCARs / UPF files / `&KIND`s, in the structure's element order |
| Charge | `--charge N`  [0] | |
| Multiplicity | `--multiplicity M`  [by parity] | below |
| k-point mesh | `--kspacing K`  [template, else 0.5 Å⁻¹] | below |
| Temperature, steps, time step | `--temperature`, `--steps`, `--timestep` (fs) | `aimd` only; otherwise the template's |
| Ensemble | `--ensemble nvt\|npt_i\|npt_f` | `aimd` only; NpT is CP2K only and, like `--vary-cell`, needs the cutoff in your template |
| Variable cell | `--vary-cell` | `opt` only; the template must set the cutoff (Pulay stress), see the program pages |
| Fine grid | `--fine-grid F` | VASP and CP2K; see their pages |

A quantity given on the command line overrides the template; one not given
keeps the template's value.

### Charge and multiplicity

The electron count is the sum of the **valence electrons of the
pseudopotentials actually chosen** (ZVAL of the POTCARs, `z_valence` of the UPF
files, the GTH electron count of CP2K), minus `--charge`. Not the atomic numbers:
a lanthanide POTCAR with the 4f shell frozen in the core (`Gd_3`) has an odd
core, and counting $Z$ would get the parity wrong.

- Without `--multiplicity`: 1 for an even count, 2 for an odd one.
- With it: the value is checked against the count's parity, and an impossible
  one (multiplicity 1 with an odd count) is an error.

A multiplicity above 1 switches the job to spin-polarised; each program page
says with which keys. ferro-py does not guess magnetic ground states: for a
high-spin system give `--multiplicity` yourself.

### k-points

`--kspacing K` is the largest allowed spacing between k-points, in Å⁻¹ with the
$2\pi$ included — the definition of VASP's `KSPACING`. Along each reciprocal
vector $\mathbf b_i$ the mesh has

$$N_i = \max\left(1, \left\lceil \frac{|\mathbf b_i|}{K} \right\rceil\right), \qquad \mathbf b_1 = 2\pi\,\frac{\mathbf a_2 \times \mathbf a_3}{\mathbf a_1 \cdot (\mathbf a_2 \times \mathbf a_3)}$$

points. Without `--kspacing`, a mesh written in the template is kept; a template
without one gets 0.5 Å⁻¹. A large cell (a glass of a few hundred atoms) comes
out as a single Γ point. Every job prints the mesh it got.

## Inputs and frames

`-i` takes one or more files in any format `ferro` reads; a quoted glob
(`'cifs/*.cif'`) is expanded by ferro-py. Only periodic structures are accepted:
a file without a cell, or not periodic along all three axes, is an error.

From a trajectory **every frame** becomes a job unless you select some, with the
same options as `ferro convert`:

| Option | Meaning |
|---|---|
| `--start N` / `--end N` | first / last frame, 0-based, both inclusive |
| `--stride N` | every Nth frame within the range |
| `--number N` | N frames spread evenly, both ends kept |

LAMMPS inputs need what `ferro convert` needs: `--atom-style` for a data file,
`--units` for a dump with velocities or forces, and `--type-map 1=Zr,2=O,…` for a
dump without an element column. A dump whose atoms are still placeholders
(`X1`, `X2`) is an error until `--type-map` names them.

## Output layout

`-o DIR` receives one sub-directory per job:

```
DIR/POSCAR/                 single-structure file: named after the file
DIR/traj_0000/              trajectory: <file stem>_<frame>, frame 0-based, 4 digits
DIR/traj_0004/
DIR/POTCAR                  VASP, several jobs: one shared copy …
DIR/traj_0000/POTCAR -> ../POTCAR      … linked from each job
```

- A missing `DIR` is created after asking; without a terminal, pass `--mkdir`.
- Nothing is written until every job has been generated: one failure (an
  element without a POTCAR, a frame without a cell) leaves `DIR` untouched.
- Jobs in one `DIR` must share their composition and element order, so that one
  POTCAR serves them all. Different compositions go to different `-o`
  directories; mixing them in one call is an error, and so is generating into a
  `DIR` that already holds the POTCAR of another composition.
- The links are relative (`../POTCAR`), so the whole directory can be moved or
  copied with `rsync -a` to the cluster.

Under each job ferro-py prints what it decided, for example

```
Wrote: md/traj_0000/  INCAR POSCAR POTCAR->../POTCAR
  POTCAR: O P Zn
  1110 electrons, multiplicity 1
  PREC Accurate, KSPACING 0.5 -> mesh 1 x 1 x 1 (Gamma only: vasp_gam can run it)
  no NCORE / NPAR / KPAR in the INCAR: add your machine's parallel settings
```

Read these lines: they are where a default you did not expect shows up.

## Errors are raised before files are read

Option errors — a task-specific option on the wrong task (`--temperature` on an
`sp`), an unknown ensemble, NpT for VASP or QE, a bad `--type-map` — are reported
before the first structure is read, so a typo does not cost the time of reading a
long trajectory.
