# Job Builders

`fe-job` generates ready-to-run input files for quantum-chemistry / DFT codes from any structure `ferro` can read (XYZ, CIF, PDB, POSCAR, LAMMPS dump, …).

```bash
fe-job                                 # overview of supported software
fe-job -s cp2k                         # CP2K-specific help
fe-job -i struct.cif -s qe -o pw.in    # generate a QE input
```

| Target | `-s` value | Output | Notes |
|---|---|---|---|
| Gaussian 16/09 | `gaussian` | `.gjf` | molecular DFT, method + basis |
| CP2K (Quickstep) | `cp2k` | `.inp` | periodic DFT/MD; GPW; precise basis DB |
| Quantum ESPRESSO | `qe` | `.in` | plane-wave `pw.x`; scf/relax/md/bands |

All three share charge/spin handling (see [Spin Estimation](spin.md)).

---

## Gaussian

```bash
fe-job -i mol.xyz -s gaussian -m B3LYP -b 6-31G* -o job.gjf
fe-job -i radical.xyz -s gaussian --charge 0 --multiplicity 2
```

`-m` method and `-b` basis set are written verbatim into the route section.

## CP2K

A full GPW/Quickstep input: `&GLOBAL`, `&FORCE_EVAL` (`&SUBSYS`, `&DFT`), and a task-specific `&MOTION`/`&VIBRATIONAL_ANALYSIS` block.

```bash
# Geometry optimisation, PBE-D3(BJ)
fe-job -i glass.xyz -s cp2k --task geo-opt --functional pbe --dispersion d3bj

# AIMD, 1500 K, CSVR thermostat
fe-job -i glass.xyz -s cp2k --task md --temperature 1500 --md-steps 50000

# Hybrid PBE0 with OT, all-electron pob basis
fe-job -i crystal.cif -s cp2k --functional pbe0 --scf ot --cp2k-basis pob-tzvp
```

Tasks: `energy`, `force`, `geo-opt`, `cell-opt`, `md`, `freq`.
Functionals: PBE/BLYP/revPBE/PBEsol (GGA), PBE0/B3LYP/HSE06 (hybrid, auto `&HF` block), SCAN/r²SCAN (meta-GGA via LIBXC).

### Precise basis & pseudopotential matching

CP2K basis-set and pseudopotential names are **element- and valence-specific** (e.g. `DZVP-MOLOPT-SR-GTH-q16` for Fe vs `-q6` for O).  `ferro` ships a database of **2829 entries** parsed from the official CP2K files (`BASIS_MOLOPT`, `BASIS_MOLOPT_UCL`, `BASIS_MOLOPT_UZH`, `BASIS_pob`, `POTENTIAL`, `POTENTIAL_UZH`), covering:

- **PBE / GGA** MOLOPT basis sets, all elements, all levels
- **SCAN** MOLOPT basis sets
- **All-electron** `pob` basis sets (`pob-dzvp`, `pob-tzvp`)

For each `&KIND` the builder looks up the exact basis name and a pseudopotential whose valence count `q` matches, preferring canonical names (`GTH-PBE` over `GTH-NLCC-PBE`).  Unknown/exotic elements fall back to a heuristic `q`-guess.  `--cp2k-basis` selects the family; the per-element exact name is resolved automatically.

## Quantum ESPRESSO

A complete `pw.x` input with `ibrav = 0`: `&CONTROL`, `&SYSTEM`, `&ELECTRONS`, plus `&IONS`/`&CELL` for relaxation/MD, then `ATOMIC_SPECIES`, `CELL_PARAMETERS angstrom`, `ATOMIC_POSITIONS angstrom`, `K_POINTS`.

```bash
# SCF, Gamma point
fe-job -i crystal.cif -s qe

# Metal: Methfessel-Paxton smearing + k-mesh
fe-job -i metal.cif -s qe --smearing mp --kpoints 8 8 8

# SCAN relaxation
fe-job -i slab.xyz -s qe --qe-task relax --qe-functional scan -o pw.in
```

Tasks: `scf`, `nscf`, `bands`, `relax`, `vc-relax`, `md`, `vc-md`.
Functionals map to `input_dft` (PBE uses the pseudopotential's own functional).
Spin reuses the same [`guess_spin`](spin.md) estimator → `nspin` / `tot_magnetization`.
Pseudopotentials are referenced as `<Element>.UPF` in `--pseudo-dir` (UPF files are user-supplied; the CP2K GTH database does not apply to QE).

---

## Charge & spin (all targets)

| Flag | Effect |
|---|---|
| `--charge INT` | Override total system charge (applied before spin estimation) |
| `--multiplicity INT` | Force 2S+1 — highest priority, disables auto-spin |
| `--auto-spin` | Estimate multiplicity from structure ([details](spin.md)) |

See the [CLI Reference](../cli-reference.md#fe-job) for the complete flag tables.
