# ferro-py job: CP2K

```bash
ferro-py job -s cp2k -t sp   -i POSCAR -o sp/
ferro-py job -s cp2k -t opt  -i POSCAR --vary-cell --template my_opt.inp -o relax/
ferro-py job -s cp2k -t aimd -i md.lammpstrj --number 20 --ensemble npt_i --template my_md.inp -o md/
```

Each job directory holds one `cp2k.inp`. What is common to the three programs —
templates, charge and multiplicity, frames, the output layout — is in
`ferro-py doc job`. The built-in templates were run with CP2K 2026.2; they avoid
keywords that 2025.2 does not know, except where noted.

## Method: GPW or GAPW

The template's `&QS METHOD` decides how atoms are described:

| `METHOD` | Electrons | Basis family | Potential |
|---|---|---|---|
| `GPW` (default) | valence, GTH pseudopotentials | `[cp2k] basis_family`  [`DZVP-MOLOPT-SR-GTH`] | `[cp2k] potential_family`  [`GTH-PBE`] |
| `GAPW` | all electrons | `[cp2k] ae_basis_family`  [`pob-TZVP-rev2`] | `ALL` |
| `GAPW_XC` | — | — | **not supported**: an error |

`GAPW_XC` is refused because its forces fail a finite-difference check, and CP2K
aborts outright when it is combined with all-electron potentials.

For an all-electron run, change the template's `METHOD` to `GAPW`,
`BASIS_SET_FILE_NAME` to `BASIS_pob` and `POTENTIAL_FILE_NAME` to
`ALL_POTENTIALS`, and consider uncommenting `GAPW_ACCURATE_XCINT` (a CP2K 2026.x
keyword; it only acts with `GAPW`). The pob-TZVP-rev2 entries from Rb on carry an
ECP rather than all electrons; such an element is an error under `GAPW`.

## Basis sets and potentials

The data files are the ones the template names in `BASIS_SET_FILE_NAME` and
`POTENTIAL_FILE_NAME` (several of each are allowed), looked up in `[cp2k]
data_dir` (`ferro-py doc config`). For each element:

- **`GPW`**: the basis is `<basis_family>-qN` and the potential
  `<potential_family>-qN` **with the same N** — basis set and pseudopotential
  must describe the same number of valence electrons. When the family has
  several `-qN` entries for an element (Rh: `-q9` and `-q17`), only those with a
  matching potential count; if more than one is left, or none, the job is an
  error that asks you to name them in `[cp2k.kinds]`. The potential is always
  written with its `-qN`: GTH-PBE would also accept the bare family name, but
  GTH-SCAN potentials must be named with their charge.
- **`GAPW`**: the all-electron entry of `<ae_basis_family>` and `POTENTIAL ALL`.
- **`[cp2k.kinds]`** names the basis and / or the potential of an element
  directly and overrides the above. A basis without `-qN` needs its potential
  named there too.

The electron count used for the multiplicity is the sum of the occupations on
the line below the potential's header in the data file.

## What is written

| Where | Keyword / section | When | Value |
|---|---|---|---|
| `&FORCE_EVAL &SUBSYS` | `&CELL`, `&COORD`, one `&KIND` per element | always | replaces the template's `&CELL`, `&COORD`, `&KIND`, `&TOPOLOGY` |
| `&DFT` | `CHARGE`, `MULTIPLICITY` | always | `--charge`, the multiplicity |
| `&DFT` | `UKS` | multiplicity above 1 | `T` |
| `&DFT &KPOINTS` | `SCHEME MONKHORST-PACK N1 N2 N3` | a mesh larger than 1×1×1 | below |
| `&DFT &XC &XC_GRID` | `USE_FINER_GRID` | `--fine-grid` | `T` |
| `&GLOBAL` | `RUN_TYPE` | `opt --vary-cell` | `CELL_OPT` (the template's `&CELL_OPT` applies) |
| `&FORCE_EVAL` | `STRESS_TENSOR` | `--vary-cell`, or an NpT ensemble | `ANALYTICAL` |
| `&MOTION &MD` | `TEMPERATURE`, `STEPS`, `TIMESTEP` | `aimd` options | K, steps, fs |
| `&MOTION &MD` | `ENSEMBLE` | `aimd --ensemble` | `NVT`, `NPT_I`, `NPT_F` |
| `&MOTION &MD &BAROSTAT` | `PRESSURE`, `TIMECON` | NpT and no `&BAROSTAT` in the template | 1.01325 bar, 1000 fs |

## k-points

Without `--kspacing`, a `&KPOINTS` section written in the template is kept.
Otherwise the mesh from the spacing (`ferro-py doc job`; default 0.5 Å⁻¹)
decides:

- **1×1×1** — the usual case for a cell of a few hundred atoms: no `&KPOINTS`,
  the Γ point only, and the template's OT solver is kept.
- **Larger** — a small cell: `&KPOINTS` with `SCHEME MONKHORST-PACK N1 N2 N3`.
  OT works at Γ only, so the template's `&OT` and `&OUTER_SCF` are replaced by
  standard diagonalisation with Broyden mixing (`ALPHA 0.4`, `NBROYDEN 10`) and
  `MAX_SCF 128`.

`MONKHORST-PACK` is the standard Monkhorst–Pack mesh: along a direction with an
even $N$ it does not contain Γ (for $N = 2$ the points are $\pm 1/4$), while the
VASP and QE meshes of the same `--kspacing` do. Once the mesh is converged this
makes no difference to the results; it does mean that an unconverged mesh is not
the same sampling in the three programs.

A template that smears occupations (`&SMEAR`) must use diagonalisation, not OT,
and set `ADDED_MOS`; either mistake is an error.

## Fine grid

`--fine-grid` switches `USE_FINER_GRID` on; its numeric factor is not used,
because CP2K has only the switch (2026.x adds `FINE_GRID_FACTOR`, default 4,
which the template can set). For a functional that is not a meta-GGA a note says
the finer grid is probably unnecessary.

## Optimisation

`opt` runs `GEO_OPT`; with `--vary-cell`, `CELL_OPT` with the analytic stress.

**A varying cell — `--vary-cell` or an NpT ensemble — needs your own template
with `&MGRID CUTOFF`**; the built-in template (`CUTOFF 400`, a generic guess) and
a template without `CUTOFF` are errors. The cell changes under a fixed grid, so
too low a cutoff gives Pulay stress and a cell that comes out too small. Find the
value with `ferro-py conv -s cp2k --scan cutoff=…` (`ferro-py doc conv`).
The template's `&CELL_OPT` has **no `TYPE` keyword**: CP2K 2026.2 removed it
(direct cell optimisation, the old default, is the only mode left) and aborts
when it is present.

## MD

`aimd` offers NVT, NPT_I (isotropic cell) and NPT_F (fully flexible cell).
For NpT the cutoff must come from your template (see Optimisation), the analytic
stress is switched on and, when the template has no `&BAROSTAT`, one is added at
1.01325 bar with a 1000 fs time constant. NPT_F
couples all cell directions; to restrict them, set `VIRIAL` in the template's
`&BAROSTAT`.

The aimd template sets `IGNORE_CONVERGENCE_FAILURE` in `&SCF`, so a step whose
SCF does not converge does not stop the trajectory (CP2K 2024.1 or later). The
sp and opt templates do not: there an unconverged SCF should stop the run.

Initial velocities are not written; CP2K draws them at `TEMPERATURE`.

## The built-in templates

All three: PBE with D3(BJ) (`&VDW_POTENTIAL`, `dftd3.dat`), `METHOD GPW`,
`BASIS_MOLOPT` + `GTH_POTENTIALS`, `CUTOFF 400` Ry, `REL_CUTOFF 60` Ry,
`EPS_DEFAULT 5.0E-11`, OT (`PRECONDITIONER FULL_ALL`, `MINIMIZER DIIS`,
`ALGORITHM IRAC`) with `EPS_SCF 5.0E-06`, `MAX_SCF 25` and an `&OUTER_SCF` of 5,
`STRESS_TENSOR ANALYTICAL`, Mulliken and Hirshfeld output off.
`GAPW_ACCURATE_XCINT` and `&XC_GRID USE_FINER_GRID` are present but commented
out.

| Task | Specific settings |
|---|---|
| `sp` | `RUN_TYPE ENERGY_FORCE`; forces and stress printed |
| `opt` | `RUN_TYPE GEO_OPT`; `&GEO_OPT` and `&CELL_OPT` with BFGS, `MAX_ITER 200`, `MAX_FORCE 4.5E-4` Hartree/bohr; `&CELL_OPT` `PRESSURE_TOLERANCE 100` bar, `EXTERNAL_PRESSURE 1.01325` bar |
| `aimd` | `RUN_TYPE MD`; NVT, 2000 steps of 1 fs at 300 K, Nosé–Hoover thermostat with `TIMECON 150` fs; `EXTRAPOLATION ASPC`, order 3; `IGNORE_CONVERGENCE_FAILURE`; positions, forces, cell and stress printed every step, restart file every 5 steps |

No parallel settings are written: CP2K takes them from the MPI launcher and
`OMP_NUM_THREADS` of your job script.

The cutoffs of the template are a starting point; test them with
`ferro-py conv -s cp2k` (`ferro-py doc conv`).
