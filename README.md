# ferro

A modular computational-chemistry post-processing toolkit written in Rust, aimed at
**periodic systems** — crystals, surfaces and glasses — and at the MD trajectories
they produce.

**Version 0.3.0** · 444 tests · MIT

> **Documentation**
> [`docs/src/`](docs/src) is the user manual (mdBook): every command, every output
> column, the data model. [`dev/`](dev) records the design rationale, the current
> state and the open work. [`CLAUDE.md`](CLAUDE.md) holds the hard rules for
> contributors (and for Claude).

## Workspace layout

```
ferro/
├── ferro-core/       # Data structures, static data, Table, AtomType, spin
├── ferro-io/         # Format readers/writers; write_table is the one exit for analysis output
├── ferro-structure/  # Supercells, vacuum layers, merging, box building, typing, clustering
├── ferro-analysis/   # Pure computation — never touches the filesystem
├── ferro-workflow/   # QC input generation (Gaussian / CP2K / Quantum ESPRESSO)
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
    └── ferro-workflow  → core
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
rustfmt management and a full pass rewrites ~80 files (see `dev/issues.md`).

## Command line

One binary, `ferro`, with subcommands grouped by *product* rather than by the crate
that implements them. Running `ferro`, `ferro <group>` or `ferro <group> <command>`
with no input prints progressively more detailed help.

```
ferro traj  gr | sq | msd | angle | vacf | rotcorr | vanhove   → stacked CSV (+ optional PNG)
ferro map   density | velocity | force | radius | sdf | chg-sdf → one .cube per input
ferro net                                                       → six stacked CSVs
ferro bader | convert | info | job
```

```bash
# g(r) + CN(r) for one pair, over every production run
ferro traj gr -i 'runs/*/prod.lammpstrj' -a P -b O --r-max 8.0 -o scan

# X-ray and neutron weighted structure factor
ferro traj sq -i traj.lammpstrj --weighting both --plot

# self-diffusion from the MSD
ferro traj msd -i traj.lammpstrj --dt 2.0 --elements Li

# glass network topology; Zn treated as a modifier
ferro net -i traj.lammpstrj --P-O=2.4 --Al-O=2.4 --Zn-O=2.6 --modifier Zn

# Bader charge partitioning
ferro bader -i CHGCAR -m neargrid

# format conversion and quick inspection
ferro convert -i traj.lammpstrj -o traj.extxyz
ferro info -i NaCl.cif

# quantum-chemistry input files
ferro job -s cp2k -i slab.cif --task geo-opt --functional pbe
```

### Batch input and output naming

`-i` always takes several files and expands glob patterns itself (quote them, so the
shell does not). Each input is analysed independently and the results stack into **one**
CSV with a `file` column — `N = 1` is just a special case of `N`, never a separate code
path. A failing input is skipped, listed in the `[inputs]` block of the output, and sets
exit code 1; argument errors fail before the first file is read.

```
<outdir>/<command>[_<table>][_<label>]_<suffix>.csv
```

`--outdir` is where products go; `-o` is a **suffix, not a path**. The label segment
says what was analysed, so `traj gr -a P -b O` writes `gr_P-O.csv` and an unfiltered run
writes `gr_all.csv`. Missing values under a column union are written as empty fields
(`NaN` on read-back) — never as zero.

`--plot` writes a 500 dpi PNG for a quick look. It is deliberately frozen at
self-inspection quality; publication figures go through the Python scripts below.

## What it computes

| Group | Methods |
|---|---|
| Trajectory | g(r) + CN(r), S(q) (Faber-Ziman, XRD/neutron weighted), MSD + self-diffusion D, bond-angle distribution, VACF + Green-Kubo, rotational correlation C₂(t), Van Hove Gs(r,τ) |
| Spatial maps | number density, velocity field, force field, hard-sphere occupancy, cluster SDF (Kabsch-aligned), charge-density cluster SDF |
| Network topology | composition, Qⁿ speciation, partner decomposition Qⁿ_m, ligand types, coordination numbers, bridge linkage |
| Charges | Bader partitioning — on-grid, near-grid, off-grid, Yu-Trinkle weight (ACF/BCF/AVF output) |
| Structure | supercells, vacuum layers, merging, formula-driven box building, type labelling, cluster detection |
| Input generation | Gaussian, CP2K (with a 2829-entry basis/pseudopotential database), Quantum ESPRESSO pw.x |

Everything is normalised per frame before time-averaging, which matters under NPT; the
derivation is in `dev/issues.md`.

### A note on Qⁿ

Since 0.3.0 the network analysis follows the literature's $Q^n_m$ convention: **n counts
homopolar connections only** (P–O–P), heteropolar ones appear as `m_<X>` columns, and the
total number of bridging connections is `n + Σm`. The number of bridging *oxygens* is a
third quantity again — a tricluster is one oxygen but two connections. Results from
0.2.x are not comparable; see `dev/overview.md`.

## File formats

| Format | Read | Write |
|---|:---:|:---:|
| XYZ | ✓ | ✓ |
| Extended XYZ (with `label:S:1`) | ✓ | ✓ |
| PDB | ✓ | ✓ |
| CIF | ✓ | ✓ |
| POSCAR / CONTCAR | ✓ | ✓ (POSCAR) |
| LAMMPS data | ✓ | ✓ |
| LAMMPS dump | ✓ | ✓ |
| CP2K `.inp` / `.restart` | ✓ | — |
| Quantum ESPRESSO `.in` | ✓ | ✓ |
| Gaussian cube | ✓ | ✓ |
| VASP CHGCAR | ✓ | — |

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

`Trajectory` is always the top-level type, even for a single-frame file. There is no
`Molecule` type — `Frame` covers molecules and periodic systems alike, distinguished by
`pbc`. Atom indices are implicit: an atom's identity is its position in `Vec<Atom>`.

```rust
pub struct Trajectory {
    pub frames: Vec<Frame>,
    pub metadata: TrajectoryMetadata,
}

pub struct Frame {
    pub atoms: Vec<Atom>,
    pub cell: Option<Cell>,                 // None = non-periodic
    pub pbc: [bool; 3],
    pub charge: i32,
    pub multiplicity: u32,
    pub bonds: Option<Vec<(usize, usize)>>,
    pub energy: Option<f64>,                // eV
    pub forces: Option<Vec<Vector3<f64>>>,  // eV/Å
    pub stress: Option<Matrix3<f64>>,       // eV/Å³
    pub velocities: Option<Vec<Vector3<f64>>>, // Å/fs
}

pub struct Atom {
    pub element: String,
    pub position: Vector3<f64>,             // Å, Cartesian
    pub label: Option<String>,              // site label, e.g. "P_3", "Fe1"
    pub mass: Option<f64>,                  // amu; None = element table value
    pub magmom: Option<f64>,                // μ_B
    pub charge: Option<f64>,                // e
}
```

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

t  = ferro.read("traj.lammpstrj", metal_units=True)
sc = ferro.supercell(t, 2, 2, 1)
ferro.write(sc, "POSCAR")

g = ferro.gr_pair(t, "P", "O", by="element", r_max=10.0, dr=0.02)
d = ferro.msd(t, dt=2.0, elements=["Li"])
```

Currently exposed: `read` / `write` / `supercell` / `add_vacuum_layer` / `merge` /
`gr_pair` / `gr_all` / `msd`, plus the `Trajectory` methods. Network analysis is not
wrapped yet. See `docs/src/python.md`.

## Python scripts

`scripts/` holds two independent groups, each with its own shared layer:

| Group | Shared layer | Scripts | Purpose |
|---|---|---|---|
| Cross-checking | `ferrocmp.py` | `compare_rdf` / `compare_angle` / `compare_sq` / `compare_sq_experiment` | Point-by-point comparison against the reference C implementations |
| Plotting | `ferroplot.py` | `plot_gr` / `plot_angle` / `plot_sq` / `plot_net` | Publication-grade PDF (SciencePlots + LaTeX) |

`trajcheck.py` records why one reference implementation's partial S(q) cannot be
reproduced (its atom typing goes stale after frame 0) — read it before concluding that
ferro disagrees.

## Test fixtures

`tests/` carries two five-frame LAMMPS trajectories, one NPT and one NVT. The NPT one
exists so that per-frame volume normalisation is actually exercised; its frames are
sampled at even intervals, because consecutive frames do not span enough volume to test
anything.

## Known limitations

- `ferro map chg-sdf` still aggregates several cubes into one SDF, unlike every other
  `map` mode; splitting it needs an intermediate format carrying sample counts first
- `bader` writes `ACF.dat` / `BCF.dat` / `AVF.dat` into the current directory under fixed
  names — running two systems in a row overwrites them
- `convert` and `job` still take a full path in `-o`, opposite to the other 11 commands
- Selecting by site label (`-x` / `-y`) only works on a single frame: g(r) requires a
  fixed particle count per type, and labels change as the run evolves
- The pyo3 0.29 bindings compile but have not been verified at runtime
- No REPL / script mode yet

The full list, with the reasoning behind each, is in `dev/progress.md` and `dev/plan.md`.
