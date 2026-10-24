# ferro

A modular computational-chemistry post-processing toolkit written in Rust, aimed at
**periodic systems** — crystals, surfaces and glasses — and at the MD trajectories
they produce.

MIT licensed. The current version and test count are in [`dev/progress.md`](dev/progress.md).

> **Documentation**
> [`docs/src/`](docs/src) is the user manual (mdBook): every command, every output
> column, the data model. [`dev/`](dev) records the design rationale, the current
> state and the open work. [`CLAUDE.md`](CLAUDE.md) holds the hard rules for
> contributors (and for Claude).

## Workspace layout

```
ferro/
├── ferro-core/       # Data structures, static data, Table, PeriodicCell
├── ferro-io/         # Format readers/writers; write_table is the one exit for analysis output
├── ferro-structure/  # Supercells, vacuum layers, merging, box building, type labelling
├── ferro-analysis/   # Pure computation — never touches the filesystem
├── ferro-render/     # Terminal rendering for `ferro doc` (markdown + LaTeX); depends on no ferro crate
├── ferro-cli/        # The single `ferro` binary
└── ferro-python/     # PyO3 bindings (its own workspace, built with maturin)
```

Dependencies are strictly one-directional, and **middle-layer crates must not depend on
each other**:

```
ferro-cli / ferro-python        ← the only layer allowed to combine crates
    ├── ferro-core
    ├── ferro-io        → core
    ├── ferro-structure → core
    ├── ferro-analysis  → core
    └── ferro-render              (depends on no ferro crate)
```

Shared types sink into `ferro-core` rather than crossing sideways: a type belongs there
exactly when two or more middle layers need to name it. `Trajectory` (produced by io,
consumed by analysis) and `Table` (produced by analysis, consumed by io) are the two
directions of that rule.

## Build

```bash
cargo build --release        # → target/release/ferro
cargo test                   # whole workspace
cargo clippy                 # must be warning-free
```

`ferro-python` is a separate workspace and is skipped by the commands above; check it
with `cd ferro-python && cargo check` after changing any public API.

Do **not** run `cargo fmt` over the whole repository — this codebase is not under
rustfmt management and a full pass rewrites ~115 files (see `dev/issues.md`).

## Command line

One binary, `ferro`, with subcommands grouped by *product* rather than by the crate
that implements them. Running `ferro`, `ferro <group>` or `ferro <group> <command>`
with no input prints progressively more detailed help.

```
ferro traj  gr | sq | msd | angle | vacf | rotcorr | vanhove | bondlife → stacked CSV
ferro map   density | velocity | force | radius | sdf | chg-sdf → one .cube per input
ferro net                                                       → six stacked CSVs
ferro dataset collect | filter | merge                         → DeePMD system directories
ferro bader | convert | info
ferro doc   <topic>                                            → the manual, compiled into the binary
```

```bash
# g(r) + CN(r) for one pair, over every production run
ferro traj gr -i 'runs/*/prod.lammpstrj' -a P -b O --r-max 8.0 -o scan

# X-ray and neutron weighted structure factor
ferro traj sq -i traj.lammpstrj --weighting both

# self-diffusion from the MSD
ferro traj msd -i traj.lammpstrj --dt 2.0 --elements Li

# glass network topology; Zn treated as a modifier
ferro net -i traj.lammpstrj --P-O=2.4 --Al-O=2.4 --Zn-O=2.6 --modifier Zn

# Bader charge partitioning
ferro bader -i CHGCAR -m neargrid

# format conversion and quick inspection
ferro convert -i traj.lammpstrj -o out --format extxyz
ferro info -i NaCl.cif

# calculation inputs (VASP / QE / CP2K) come from the Python package
ferro-py job -s cp2k -t opt -i slab.cif -o opt/
```

### Batch input and output naming

`-i` always takes several files and expands glob patterns itself (quote them, so the
shell does not). Each input is analysed independently and the results stack into **one**
CSV with a `file` column — `N = 1` is just a special case of `N`, never a separate code
path. A failing input is skipped, listed in the `[inputs]` block of the output, and sets
exit code 1; argument errors fail before the first file is read.

```
<-o dir>/<command>[_<table>][_<label>][_<input>][_<suffix>].csv
```

`-o` is where products go — a path, for every command: the output **directory** for the
runs that write several products, the output **file** for `convert`. A missing
directory is created after a `[y/N]` prompt, or straight away with `--mkdir` (which
scripts and CI must pass, having no terminal to answer). `-s` tags the batch. The label segment
says what was analysed, so `traj gr -a P -b O` writes `gr_P-O.csv` and an unfiltered run
writes `gr_all.csv`. Missing values under a column union are written as empty fields
(`NaN` on read-back) — never as zero.

ferro writes csv only and draws no figures; publication figures are drawn by
`ferro-py plot` (see `docs/src/python/plot.md`).

## What it computes

| Group | Methods |
|---|---|
| Trajectory | g(r) + CN(r), S(q) (Faber-Ziman, XRD/neutron weighted), MSD + self-diffusion D, bond-angle distribution, VACF + Green-Kubo, rotational correlation C₂(t), Van Hove Gs(r,τ), bond lifetimes |
| Spatial maps | number density, velocity field, force field, hard-sphere occupancy, cluster SDF (Kabsch-aligned), charge-density cluster SDF |
| Network topology | composition, Qⁿ speciation, partner decomposition Qⁿ_m, ligand types, coordination numbers, bridge linkage |
| Charges | Bader partitioning — on-grid, near-grid, off-grid, Yu-Trinkle weight (ACF/BCF/AVF output) |
| ML datasets | DeePMD system collection, frame filtering, merging and train/valid/test splitting; GPUMD/NEP export |
| Structure | supercells, vacuum layers, merging, formula-driven box building, type labelling |
| Input generation | VASP, CP2K, Quantum ESPRESSO pw.x — single point, optimisation, AIMD, convergence tests; in the Python package (`ferro-py job`, `ferro-py conv`) |

Everything is normalised per frame before time-averaging, which matters under NPT; the
derivation is in `dev/issues.md`.

### A note on Qⁿ

Since 0.3.0 the network analysis follows the literature's $Q^n_m$ convention: **n counts
homopolar connections only** (P–O–P), heteropolar ones appear as `m_<X>` columns, and the
total number of bridging connections is `n + Σm`. The number of bridging *oxygens* is a
third quantity again — a tricluster is one oxygen but two connections. Results from
0.2.x are not comparable; see `dev/overview.md`.

## File formats

Structures and trajectories: XYZ / extended XYZ, PDB, CIF, POSCAR / CONTCAR, LAMMPS data and dump,
Gaussian cube, VASP CHGCAR, CP2K and QE inputs (read-only), and the AIMD outputs of VASP
(OUTCAR, vasprun.xml), CP2K and QE, plus DeePMD npy systems. The authoritative table —
what is recognised by which name, read or written — is printed by `ferro convert` with no
input and lives in [`docs/src/cli-reference.md`](docs/src/cli-reference.md).

Analysis products are CSV only, all written through the single `write_table` exit.

## Internal units

DeePMD-kit / VASP convention. Conversions go through the enums in `units.rs`.

| Quantity | Unit | | Quantity | Unit |
|---|---|---|---|---|
| Length | Å | | Time | fs |
| Energy | eV | | Mass | amu |
| Force | eV/Å | | Charge | e |
| Stress | eV/Å³ | | Temperature | K |

## Core data model

See [`docs/src/data-model.md`](docs/src/data-model.md).

## Python bindings

Built with [maturin](https://www.maturin.rs/) — `cargo build` from the workspace root
does not apply, because `ferro-python` is its own workspace.

```bash
cd ferro-python
maturin build --interpreter "$(which python)"
pip install target/wheels/*.whl
```

```python
import ferro

t  = ferro.read("traj.lammpstrj", units="metal")
sc = ferro.supercell(t, 2, 2, 1)
ferro.write(sc, "POSCAR", format="vasp")

g = ferro.gr_pair(t, "P", "O", by="element", r_max=10.0, dr=0.02)
d = ferro.msd(t, dt=2.0, elements=["Li"])
```

The full API is in [`docs/src/python/api.md`](docs/src/python/api.md); the `ferro-py` commands
are read with `ferro-py doc`.

## Python scripts

`scripts/` holds the cross-checks against the reference C implementations
(`compare_rdf` / `compare_angle` / `compare_sq` / `compare_sq_experiment`, sharing
`ferrocmp.py`) and the space-group table tools (`gen_spacegroups` / `crosscheck_spacegroups`).
The publication figures moved into ferro-python as `ferro-py plot`.

`trajcheck.py` records why one reference implementation's partial S(q) cannot be
reproduced (its atom typing goes stale after frame 0) — read it before concluding that
ferro disagrees.

## Test fixtures and known limitations

The fixture list is in [`dev/progress.md`](dev/progress.md); open work and known
limitations are indexed in [`dev/plan.md`](dev/plan.md).
