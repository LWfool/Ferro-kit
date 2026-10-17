# ferro-py conv: Convergence Tests

`ferro-py conv` writes the same single points for every value of one parameter;
you run them; `ferro-py conv analyze` reads the outputs back and compares each
value with the most accurate one.

```bash
ferro-py conv -s vasp --scan fine-grid=2,2.5,3,4 -i POSCAR --template r2scan.INCAR -o ng/
ferro-py conv -s vasp --scan encut=500:900:100 --fine-grid 3 -i md.lammpstrj -o encut/
ferro-py conv -s cp2k --scan cutoff=300:700:100 -i POSCAR -o cutoff/
ferro-py conv -s qe --scan ecutwfc=40:100:10 -i POSCAR -o ecut/
# … run every job …
ferro-py conv analyze encut/
```

## What can be scanned

| Program | `--scan` | Values | What changes in the sp template |
|---|---|---|---|
| VASP | `encut=…` | eV | `ENCUT` |
| VASP | `fine-grid=…` | × NG | `NGXF/NGYF/NGZF` $= F \times$ the coarse grid (`ferro-py doc job vasp`) |
| CP2K | `cutoff=…` | Ry | `&MGRID CUTOFF` |
| CP2K | `rel-cutoff=…` | Ry | `&MGRID REL_CUTOFF` |
| CP2K | `finer-grid` | off, on | `&XC_GRID USE_FINER_GRID` |
| CP2K | `accurate-xcint` | off, on | `&QS GAPW_ACCURATE_XCINT` (`METHOD GAPW` only; CP2K 2026.x, absent in 2025.2) |
| QE | `ecutwfc=…` | Ry | `&SYSTEM ecutwfc` |

Values are `START:STOP:STEP` (both ends included) or a list `a,b,c`; at least two
distinct values, used in increasing order. The two CP2K switches take no values.

The jobs are single points built from the **sp template** (yours with
`--template`, else the configured or built-in one) exactly as `ferro-py job`
builds them, with the scanned value written in. Every other option of `job`
applies: `--charge`, `--multiplicity`, `--kspacing`, frame selection. A test on
several frames of a trajectory (`--number 3`) is safer than on one structure.

`--fine-grid` fixes the fine grid at every value — VASP's factor, CP2K's switch.
It cannot be combined with `--scan fine-grid` or `--scan finer-grid`, which set
the same thing.

For QE only `ecutwfc` is scanned. `ecutrho` stays as the template writes it, or
is the largest UPF suggestion; at a value where the suggestion is below
4 × `ecutwfc`, that value gets 4 × `ecutwfc` instead (`ferro-py doc job qe`).

## In which order

**VASP with a meta-GGA (r2SCAN).** The fine grid and `ENCUT` both matter, and
the fine grid is a multiple of the coarse grid that `ENCUT` sets:

1. Write a reasonable `ENCUT` into the sp template and scan the fine grid:
   `--scan fine-grid=2,2.5,3,4`. Pick the smallest factor that is converged.
2. Scan `ENCUT` with that factor fixed: `--scan encut=… --fine-grid F`. Each
   `ENCUT` recomputes the coarse grid, and the fine grid follows it.

Write both results into your templates. A variable-cell optimisation needs the
tested `ENCUT` in its template (`ferro-py doc job vasp`).

**VASP with a GGA.** The fine grid rarely matters; scan `ENCUT`.

**CP2K.** Scan `cutoff`, then `rel-cutoff` at the chosen cutoff, then
`finer-grid` (for a meta-GGA). With the finer grid switched on, a lower `CUTOFF`
may be enough, though not necessarily for every element. For `METHOD GAPW`,
scan `accurate-xcint` as well.

**QE.** Scan `ecutwfc`. `ecutrho` converges easily; the largest UPF suggestion
(or 4 × `ecutwfc`, whichever is larger) is normally enough.

## Output layout

```
DIR/ferro-conv.json        parameter, values and every job (read by analyze)
DIR/POTCAR                 VASP: one copy, linked from every job
DIR/encut_500/traj_0000/
DIR/encut_500/traj_0004/
DIR/encut_600/traj_0000/
…
```

`DIR/<parameter>_<value>/<structure>/`, each job laid out as `ferro-py job`
writes it. Run every job in its own directory; the file names `analyze` looks
for are those of the built-in templates (below).

## conv analyze

```bash
ferro-py conv analyze encut/
ferro-py conv analyze ng/ --ref 4 --df-rms 2
```

`analyze` reads `DIR/ferro-conv.json` — never the directory names — and the
output of every job:

| Program | Output files tried, in order |
|---|---|
| VASP | `OUTCAR`, `vasprun.xml` |
| CP2K | `cp2k.out`, then any `*.out`, `*.log` |
| QE | `pw.out`, then any `*.out` |

The **reference** is the highest value (`on` for a switch), or `--ref VALUE`.
For each structure and value, against the same structure at the reference:

| Column | Definition | Unit |
|---|---|---|
| `dE_meV_atom` | $\Delta E = (E - E_{\mathrm{ref}}) / N$ | meV/atom |
| `dP_kbar` | $\Delta P = P - P_{\mathrm{ref}}$, $P = \mathrm{tr}\,\sigma / 3$ | kbar |
| `dF_max_meV_A` | $\max_i \lvert \mathbf F_i - \mathbf F_i^{\mathrm{ref}} \rvert$ | meV/Å |
| `dF_rms_meV_A` | $\sqrt{\frac{1}{N} \sum_i \lvert \mathbf F_i - \mathbf F_i^{\mathrm{ref}} \rvert^2}$ | meV/Å |
| `within` | `yes` when every difference is within its threshold | |

The last complete step of each output is used. A job that cannot be compared
says why in `status`: `no output`, `SCF not converged`, `no complete step`,
`unreadable (…)`, or `reference …` when the reference job is the one missing.

### Thresholds

| Option | Default | Basis |
|---|---|---|
| `--de` | 1 meV/atom | the common criterion for total energies (SSSP uses 2 meV/atom for cohesive energies) |
| `--dp` | 1 kbar | SSSP's tolerance for the pressure, which matters for `vc-relax` |
| `--df-rms` | 5 meV/Å | about a tenth of the force error of typical machine-learned potentials |
| `--df-max` | 10 meV/Å | a looser bound for the single worst atom |

The thresholds only mark rows; they decide nothing. Set them to what your
application needs: training data for a potential needs tighter forces than a
relaxation does.

### What is printed and written

The terminal shows, per value, the **worst frame** of each difference and whether
all frames are within the thresholds, then the first (lowest) value for which
they are. The CSV, `DIR/conv_<parameter>.csv` unless `-o FILE`, has one row per
value and structure:

```
value, structure, frame, status, n_atoms, energy_eV, dE_meV_atom, pressure_kbar,
dP_kbar, dF_max_meV_A, dF_rms_meV_A, within
```

"The first value within every threshold" is the lowest value that passes, not a
proof of convergence. When the differences do not fall steadily — a value that
passes followed by a higher one that does not — look at the table and choose the
value above which every row passes.
