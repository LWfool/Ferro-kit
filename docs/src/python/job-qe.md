# ferro-py job: Quantum ESPRESSO

```bash
ferro-py job -s qe -t sp   -i POSCAR -o sp/
ferro-py job -s qe -t opt  -i 'cifs/*.cif' --vary-cell -o relax/
ferro-py job -s qe -t aimd -i POSCAR --temperature 1000 --timestep 1 -o md/
```

Each job directory holds one `pw.in` for `pw.x`. What is common to the three
programs — templates, charge and multiplicity, frames, the output layout — is in
`ferro-py doc job`.

The header comment of the template is not copied into `pw.in` (it is written for
whoever edits the template), and neither are comments inside the namelists;
comments between the cards are kept.

## Pseudopotentials

ferro-py reads every `*.upf` file in `[qe] pseudo_dir` (`ferro-py doc config`)
and takes from each file's `PP_HEADER` (UPF v1 or v2) the element, the number of
valence electrons and the suggested cutoffs. Then, per element:

- the file named in `[qe.pseudos]`, if any;
- otherwise the one file of that element in the directory;
- several files and no `[qe.pseudos]` entry is an error listing them; none is an
  error too.

The functional is the one of the UPF files: the template does not set
`input_dft`. `pseudo_dir` is written into `&CONTROL` as an absolute path.

## What is written

| Where | Key / card | When | Value |
|---|---|---|---|
| `&CONTROL` | `pseudo_dir` | always | `[qe] pseudo_dir` |
| `&SYSTEM` | `nat`, `ntyp` | always | from the structure |
| `&SYSTEM` | `ibrav` | always | 0; `celldm(1..6)`, `A`, `B`, `C`, `cosAB`, `cosAC`, `cosBC` are removed |
| `&SYSTEM` | `ecutwfc`, `ecutrho` | the template does not set it | the largest suggestion among the UPF files, each key on its own |
| `&SYSTEM` | `tot_charge` | always | `--charge` |
| `&SYSTEM` | `nspin`, `tot_magnetization` | multiplicity above 1 | 2, multiplicity − 1 |
| cards | `ATOMIC_SPECIES` | always | element, mass, UPF file |
| cards | `CELL_PARAMETERS {angstrom}` | always | the cell |
| cards | `ATOMIC_POSITIONS {angstrom}` | always | Cartesian positions |
| cards | `K_POINTS automatic` | `--kspacing`, or the template has no `K_POINTS` | `N1 N2 N3  0 0 0` (Γ-centred) |
| `&CONTROL` | `calculation` | `opt --vary-cell` | `'vc-relax'` |
| `&CELL` | `cell_dynamics`, `press_conv_thr` | `opt --vary-cell` and the template has no `&CELL` | `'bfgs'`, 0.5 kbar |
| `&IONS` | `tempw` | `aimd --temperature` | K |
| `&CONTROL` | `nstep` | `aimd --steps` | steps |
| `&CONTROL` | `dt` | `aimd --timestep` | fs ÷ 0.048378 (Rydberg atomic units: 1 fs = 20.67) |

The cell always comes from `CELL_PARAMETERS`, which is why `ibrav` is forced to
0: with any other `ibrav`, pw.x builds the cell from `celldm` and ignores the
card, so a template with `ibrav = 2` would silently compute a different cell.

A `K_POINTS` card written in the template (`gamma`, an explicit list, a shifted
mesh) is kept unless `--kspacing` is given.

With multiplicity 1, `nspin` stays as the template has it.

## Cutoffs: test `ecutwfc`

A template without `ecutwfc` / `ecutrho` gets the largest value suggested in the
UPF headers, and the job prints which. The suggestion is a starting point, not a
converged value:

- **`ecutwfc` must be tested** for your system;
- `ecutrho` converges easily; the largest UPF suggestion is normally enough.

Write the tested values into your template. A convergence test for QE
(`ferro-py conv -s qe`) is not offered yet; `ferro-py doc conv` covers VASP and
CP2K.

## MD

`aimd` is **NVT only**: `vc-md` has no thermostat, so `--ensemble npt_i|npt_f`
is an error (it is a CP2K option). The built-in template uses `verlet` with the
stochastic velocity-rescaling thermostat (`ion_temperature = 'svr'`,
`nraise = 150`) and `nosym = .true.`.

`nosym` is required for MD: the thermostat's random velocities break the
symmetry of the initial structure, and pw.x stops at the second step
(`checkallsym: some of the original symmetry operations not satisfied`). Keep it
in your own aimd template.

`nraise` counts **steps**, not fs: 150 steps are 150 fs at the template's 1 fs
time step, and become 300 fs with `--timestep 2`. Adjust `nraise` in the template
when you change the time step and want the same thermostat time constant.

## `--fine-grid`

Ignored, with a note: pw.x has no separate fine grid — the dense grid follows
`ecutrho`.

## The built-in templates

All three: `ibrav = 0`, D3(BJ) (`vdw_corr = 'dft-d3'`, `dftd3_version = 4`),
`tstress` and `tprnfor` on, `outdir = './tmp/'`, `disk_io = 'low'`,
`diagonalization = 'david'`, `electron_maxstep = 200`, no cutoffs (taken from the
UPF files, see above).

| Task | Specific settings |
|---|---|
| `sp` | `calculation = 'scf'`, `conv_thr = 1e-10` |
| `opt` | `calculation = 'relax'` (`'vc-relax'` with `--vary-cell`), `nstep = 200`, `etot_conv_thr = 1e-4`, `forc_conv_thr = 1e-3`, `ion_dynamics = 'bfgs'`, `conv_thr = 1e-10` |
| `aimd` | `calculation = 'md'`, `nstep = 2000`, `dt = 20.67` (1 fs), `ion_dynamics = 'verlet'`, `ion_temperature = 'svr'`, `tempw = 300`, `nraise = 150`, `nosym = .true.`, `conv_thr = 1e-8` |

No parallel settings are written: they belong on the `pw.x` command line
(`-nk`, `-nd`, …) of your job script.
