# ferro-py job: VASP

```bash
ferro-py job -s vasp -t sp   -i POSCAR -o sp/
ferro-py job -s vasp -t opt  -i POSCAR --vary-cell --template my_opt/INCAR -o relax/
ferro-py job -s vasp -t aimd -i traj.lammpstrj --number 10 --temperature 1500 -o md/
```

Each job directory holds `INCAR`, `POSCAR` and `POTCAR`. There is **no
`KPOINTS` file**: the mesh is the `KSPACING` tag of the INCAR (below). What is
common to the three programs — templates, charge and multiplicity, frames, the
output layout — is in `ferro-py doc job`.

## POSCAR and POTCAR

`POSCAR` is written in direct (fractional) coordinates, with the element line;
atoms are sorted by (Z, symbol), so each element appears once and in the same
order as in the POTCAR.

`POTCAR` is the concatenation of one POTCAR per element, in that order, read
from `[vasp] potcar_dir` (`ferro-py doc config`) as `<potcar_dir>/<name>/POTCAR`.
The `<name>` of an element is, in this order:

1. its entry in `[vasp.potcar]`;
2. the PBE POTCAR that VASP recommends for it (potpaw.64, the vasp.at list of
   available pseudopotentials) — table below;
3. the element symbol itself.

When `<potcar_dir>/<name>/POTCAR` does not exist, the error lists the variants
of that element that do exist (`Zn`, `Zn_pv`, `Zn_sv_GW`, …), to pick one in
`[vasp.potcar]`.

| Elements | Recommended POTCAR |
|---|---|
| Li K Ca Sc Ti V Rb Sr Y Zr Nb Mo Cs Ba W Fr Ra | `<El>_sv` |
| Na Cr Mn Fe Co Ni Cu Tc Ru Rh Hf Ta | `<El>_pv` |
| Ga Ge In Sn Tl Pb Bi Po | `<El>_d` |
| Pr Nd Pm Sm Gd Tb Dy Ho Tm Lu | `<El>_3` |
| Eu Er Yb | `<El>_2` |
| every other element | `<El>` |

## What is written into the INCAR

| Tag | When | Value |
|---|---|---|
| `NELECT` | `--charge` is not 0 | $\sum$ ZVAL − charge |
| `ISPIN`, `NUPDOWN` | multiplicity above 1 | `ISPIN = 2`, `NUPDOWN` = multiplicity − 1 |
| `NUPDOWN` | `--multiplicity 1` given and the template has `ISPIN = 2` | 0 |
| `KSPACING` | always | `--kspacing`, else the template's, else 0.5 |
| `ISIF` | `opt --vary-cell` | 3 |
| `TEBEG` (and `TEEND` if present) | `aimd --temperature` | the temperature |
| `NSW`, `POTIM` | `aimd --steps`, `--timestep` | steps, fs |
| `ALGO`, `ISEARCH` | the template sets `METAGGA` | `All`, 1 |
| `NGXF`, `NGYF`, `NGZF` | meta-GGA or `--fine-grid` | below |

Everything else is the template's. The tag lines ferro-py changes keep their
comments; tags it adds go at the end.

### k-points

`KSPACING` with VASP's default `KGAMMA = .TRUE.` gives a Γ-centred mesh; VASP
computes the mesh itself from the tag, and ferro-py prints the same mesh under
each job together with the `PREC` in use:

```
PREC Accurate, KSPACING 0.5 -> mesh 1 x 1 x 1 (Gamma only: vasp_gam can run it)
```

A Γ-only mesh can run with `vasp_gam`, which is faster than `vasp_std`.

## Precision

`PREC` must be written **in full** — `Low`, `Medium`, `High`, `Normal`,
`Single`, `SingleN` or `Accurate` (any case). An abbreviation such as `PREC = Acc`
is an error: ferro-py computes the fine grid from `PREC` (below) and does not
guess what a short form means. Without a `PREC` tag VASP's default `Normal` is
assumed. The built-in templates use `Accurate`.

## Variable-cell optimisation needs ENCUT

`opt --vary-cell` sets `ISIF = 3`, so the cell relaxes with a fixed set of
plane waves. Too low a cutoff then gives Pulay stress and a cell that comes out
too small. **The template must therefore set `ENCUT`**; without it the job is an
error. Find the value with a convergence test (`ferro-py doc conv`) and write it
into your opt template. For a fixed cell, a template without `ENCUT` runs at
VASP's default, the largest ENMAX of the POTCARs.

## Meta-GGA (r2SCAN) and the fine grid

To run r2SCAN, replace `GGA = PE` in the template by `METAGGA = R2SCAN` and, for
D3(BJ), add the r2SCAN damping parameters (`VDW_S8`, `VDW_A1`, `VDW_A2`) next to
`IVDW = 12`. `LASPH` and `LMIXTAU` need not be written: VASP switches them on for
r2SCAN itself. ferro-py then:

- sets `ALGO = All` and `ISEARCH = 1` — the Davidson default converges poorly for
  meta-GGAs;
- writes an explicit fine grid, because the meta-GGA energy is sensitive to the
  grid on which $\tau$ is integrated (egg-box noise in the forces).

The fine grid is a factor $F$ times VASP's coarse FFT grid:

$$N_i^{\mathrm{fine}} = \mathrm{good}\left(\lceil F \cdot N_i \rceil\right), \qquad N_i = \mathrm{good}\left(\mathrm{round}\left(\frac{k\, G_{\mathrm{cut}} |\mathbf a_i|}{2\pi}\right)\right), \qquad G_{\mathrm{cut}} = \sqrt{\frac{E_{\mathrm{cut}}}{\hbar^2/2m_e}}$$

with $k = 4$ for `PREC = Accurate` or `High` and $k = 3$ otherwise,
$\hbar^2/2m_e = 3.80998$ eV·Å², $E_{\mathrm{cut}}$ the template's `ENCUT` (else
the largest ENMAX), and $\mathrm{good}(n)$ the smallest even number $\geq n$
whose only prime factors are 2, 3, 5, 7 — the sizes VASP itself picks.

| Situation | Fine grid |
|---|---|
| `METAGGA` set, no `NGXF` in the template, no `--fine-grid` | written with $F = 2$ |
| `METAGGA` set, `NGXF` in the template, no `--fine-grid` | the template's |
| `--fine-grid F` | written with that $F$, overriding the template |
| GGA with `--fine-grid F` | written as asked, with a note that a GGA rarely needs it |

With `PREC = Accurate`, $F = 2$ is what VASP would use by default anyway; the
point of writing it is that the convergence test can vary it. Systems with Zn
may need $F = 3$ (a note says so). Which $F$ is enough is a question for
`ferro-py conv -s vasp --scan fine-grid=…` (`ferro-py doc conv`): test the fine
grid first at a fixed `ENCUT`, then `ENCUT` at the chosen grid.

Array tags tied to the atom order (`MAGMOM`, `LDAUL`, …) are not rewritten; a
template that has them must match the sorted element order.

## MD

`aimd` is **NVT only** (`--ensemble npt_i|npt_f` is a CP2K option and an error
here). The built-in template uses the Nosé–Hoover thermostat (`MDALGO = 2`,
`SMASS = 1`), 2000 steps of 1 fs at 300 K, with `ISIF = 2` so the stress is still
written at every step.

## The built-in templates

All three: `PREC = Accurate`, `LREAL = Auto`, `ISMEAR = 0` with `SIGMA = 0.05`,
`ISYM = 0`, PBE (`GGA = PE`) with D3(BJ) (`IVDW = 12`), `KSPACING = 0.5`,
`ENCUT` unset.

| Task | Specific settings |
|---|---|
| `sp` | `IBRION = -1`, `NSW = 0`, `ISIF = 2` (stress computed), `EDIFF = 1E-05` |
| `opt` | `IBRION = 2`, `NSW = 200`, `ISIF = 2` (3 with `--vary-cell`), `EDIFF = 1E-05`, `EDIFFG = -0.02` |
| `aimd` | `IBRION = 0`, `MDALGO = 2`, `SMASS = 1`, `NSW = 2000`, `POTIM = 1.0`, `TEBEG = 300`, `ISIF = 2`, `EDIFF = 1E-04`, `LWAVE` / `LCHARG` on for restarts |

No parallel settings (`NCORE`, `NPAR`, `KPAR`) are written: they depend on the
machine. A job whose INCAR has none of them prints a reminder.
