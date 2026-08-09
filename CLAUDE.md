# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

```bash
cargo build                                                    # build entire workspace
cargo build --release
cargo build --package ferro-core                             # single crate
cargo test                                                     # all tests
cargo test --package ferro-io                                # single crate
cargo test --package ferro-core test_basic_molecule          # single test
cargo fmt
cargo clippy
cargo check

# CLI (dev)
cargo run --bin fe-convert -- -i input.xyz -o output.pdb
cargo run --bin fe-info    -- -i input.xyz
cargo run --bin fe-job     -- -i input.xyz -s gaussian -m B3LYP -o job.gjf
cargo run --bin fe-traj    -- -m gr  -i traj.dump
cargo run --bin fe-traj    -- -m gr  -i 'runs/*/prod.dump' -a P -b O -o scan   # batch
cargo run --bin fe-cube    -- -m sdf -i traj.dump --qn 3 --modifier Zn

# Python bindings (requires maturin; ferro-python not yet in workspace)
cd ferro-python && maturin develop
```

Test fixtures: `tests/` (water.xyz, water.pdb, water.cif)

---

## Architecture

Cargo workspace with a strict layered dependency graph. Middle-layer crates must NOT depend on each other — only the top-layer entry points combine them.

```
ferro-cli / ferro-python        ← only layer that combines multiple crates
    ├── ferro-core                ← pure data structures + static reference data
    ├── ferro-io        → core    ← format readers/writers
    ├── ferro-structure → core    ← supercell, vacuum, merge, box estimation
    ├── ferro-analysis  → core    ← md/, dft/ (future), ml/ (future)
    └── ferro-workflow  → core    ← QC software input file builders
```

**Shared types go down, not sideways.** The rule above forbids `ferro-io → ferro-analysis`,
but it never forbade the two from *sharing* types — the sharing point is `ferro-core`.

> A type belongs in `ferro-core` **iff two or more middle layers need to name it.**

Two worked examples, one per direction of data flow:

| Direction | Shared type | Producer | Consumer |
|---|---|---|---|
| structures in | `Trajectory` | `ferro-io` | `ferro-analysis` |
| results out | `Table` | `ferro-analysis` | `ferro-io` |

Analysis-private intermediates (`GrResult`, `SqResult`, …) stay in `ferro-analysis`;
only their serialisable projection (`Table`, via `to_tables()`) crosses a layer boundary.
This is the same pattern as pymatgen's `MSONable` and OVITO's `DataTable` — neither
lets its io module import its analysis module.

### Crate responsibilities

| Crate | Role |
|---|---|
| `ferro-core` | `Atom`, `Frame`, `Trajectory`, `Cell`, `Table`; static element/compound data; error types; unit conversion |
| `ferro-io` | Format readers (`read_xyz`, `read_pdb`, ...) and writers; returns/accepts `Trajectory`. `write_table` is the single exit point for analysis results |
| `ferro-structure` | Supercell, vacuum layer, merge, initial box estimation from compound data |
| `ferro-analysis` | Pure computation — no filesystem. Results expose `to_tables() -> Vec<(String, Table)>`. Sub-modules: `md/` (g(r), S(q), MSD, angle, VACF, rotcorr, VanHove, cube density, cluster SDF), `dft/` (future), `ml/` (future) |
| `ferro-workflow` | QC input file builders: `GaussianJobBuilder`, `GromacsTopologyBuilder`, etc. |
| `ferro-cli` | CLI + REPL + batch mode (shared interpreter); `batch.rs` drives multi-input runs, generic over the result type |
| `ferro-python` | PyO3 wrappers only; pure Rust libs have zero Python awareness |

---

## Core Data Model

**Primary use case is periodic systems** (crystals, surfaces). Non-periodic (molecular) systems are secondary.

**Always use `Trajectory` as the top-level type**, even for single-frame files.
- Single-frame file → `Trajectory { frames: vec![frame_0] }`
- MD trajectory → `Trajectory { frames: vec![frame_0, frame_1, ...] }`

All module APIs accept/return `Trajectory`. All core types derive `Clone`.

### Atom
```rust
pub struct Atom {
    pub element: String,
    pub position: Vector3<f64>,    // Å, Cartesian
    pub label: Option<String>,     // human-readable tag e.g. "Fe1", "Fe2"
    pub mass: Option<f64>,         // None = look up from elements table
    pub magmom: Option<f64>,       // initial magnetic moment (DFT input)
    pub charge: Option<f64>,       // Bader/DDEC charge (post-processing result)
}
```

Atom index is **implicit** (position in `Vec<Atom>`). No `index` field stored to avoid inconsistency.

### Frame (≈ ASE Atoms)
```rust
pub struct Frame {
    pub atoms: Vec<Atom>,
    pub cell: Option<Cell>,                      // None = non-periodic
    pub pbc: [bool; 3],                          // periodic in x/y/z
    pub charge: i32,                             // total system charge
    pub multiplicity: u32,                       // 2S+1; unpaired electrons = multiplicity-1
    pub bonds: Option<Vec<(usize, usize)>>,      // optional bond list (i,j)
    // post-processing results (written back after calculation):
    pub energy: Option<f64>,
    pub forces: Option<Vec<Vector3<f64>>>,
    pub stress: Option<Matrix3<f64>>,
    pub velocities: Option<Vec<Vector3<f64>>>,
}
```

`pbc` controls periodicity: `[false,false,false]` = molecular, `[true,true,true]` = crystal, `[true,true,false]` = surface/slab.

**`Molecule` struct does not exist** — `Frame` covers all cases.

### ferro-core file structure
```
ferro-core/src/
├── lib.rs
├── atom.rs
├── cell.rs
├── frame.rs
├── trajectory.rs
├── data/
│   ├── mod.rs
│   ├── elements.rs  # static: symbol, atomic number, mass, oxidation states, electron config
│   └── compounds.rs # static: name, formula, molecular_mass, density — for box estimation
├── units.rs
└── error.rs
```

### Cell
```rust
pub struct Cell {
    pub matrix: Matrix3<f64>,  // row vectors a, b, c in Å
}
```

| Method | Purpose |
|---|---|
| `from_matrix` / `from_lengths_angles` | constructors |
| `lengths() -> [f64; 3]` | a, b, c |
| `angles() -> [f64; 3]` | α, β, γ in degrees |
| `volume() -> f64` | |
| `fractional_to_cartesian` / `cartesian_to_fractional` | coordinate transforms |
| `wrap_position` | fold Cartesian position back into box |
| `minimum_image` | minimum image convention for PBC distance/vector |

NPT trajectories (varying box per frame) are naturally handled: `cell: Option<Cell>` lives inside `Frame`, so each frame carries its own cell.

### Units
Internal standard follows **DeePMD-kit / VASP convention**:

| Quantity | Internal unit |
|---|---|
| Length | Å |
| Energy | eV |
| Force | eV/Å |
| Stress | eV/Å³ |
| Time | fs |
| Mass | amu |
| Charge | e |
| Temperature | K |

Self-implemented (no `uom` or external unit crates). Enum-based conversion:
```rust
pub enum LengthUnit   { Angstrom, Bohr, Nanometer }
pub enum EnergyUnit   { EV, Hartree, KcalPerMol, KJPerMol, Wavenumber }
pub enum PressureUnit { EVPerAng3, GPa, Kbar }
pub enum TimeUnit     { Femtosecond, Picosecond }

pub fn convert_length(value: f64, from: LengthUnit, to: LengthUnit) -> f64
pub fn convert_energy(value: f64, from: EnergyUnit, to: EnergyUnit) -> f64
pub fn convert_pressure(value: f64, from: PressureUnit, to: PressureUnit) -> f64
```

### Error handling
- Library crates: `ferro_core::error::ChemError` / `ferro_core::Result<T>`
- CLI: `anyhow::Result`

---

## Static Reference Data (ferro-core/data/)

### elements.rs
Per-element: symbol, atomic number, atomic mass, common oxidation states, electron configuration, electronegativity.

### compounds.rs
Used by `ferro-structure` to estimate initial MD simulation box size:
```rust
pub struct CompoundData {
    pub name: &'static str,
    pub formula: &'static str,
    pub molecular_mass: f64,       // g/mol
    pub density: Option<f64>,      // g/cm³ at standard conditions; None for gases
    pub cas: Option<&'static str>,
}
```
Box estimation logic (V = Σ n_i·M_i / ρ_mix·Nₐ) lives in `ferro-structure`, not here.

---

## Execution Modes

### 1. One-shot CLI
```bash
ferro convert -i a.xyz -o b.pdb
```

### 2. Interactive REPL (rustyline)
```bash
ferro
ferro> read water.xyz
ferro> supercell 2 2 1
ferro> write POSCAR
```

### 3. Batch / script mode
Same interpreter as REPL, non-interactive. For shell script integration.
```bash
ferro -f workflow.mf
echo -e "read water.xyz\nsupercell 2 2 1\nwrite POSCAR" | ferro
```

### 4. Python library
`import ferro` via PyO3 in `ferro-python`.

### ferro-cli internal structure
```
ferro-cli/src/
├── main.rs              # REPL / batch mode detection (placeholder)
├── lib.rs               # re-exports args, help, io_dispatch, plot
├── args/
│   ├── common.rs
│   ├── traj.rs          # TrajMode, SqWeightingCli
│   ├── corr.rs
│   └── cube.rs
├── help.rs              # print_fe_traj_overview / print_traj_help / print_corr_help / print_cube_help
├── io_dispatch.rs       # read_trajectory (format detection)
├── plot.rs              # plot_gr / plot_angle / plot_sq / open_plot  (plotters)
└── bin/
    ├── convert.rs
    ├── info.rs
    ├── job.rs
    ├── traj.rs
    ├── corr.rs
    ├── cube.rs
    └── network.rs
```

---

## Python Bindings (ferro-python)

All PyO3 glue lives here. Library crates have zero Python awareness.

```
ferro-python/src/
├── lib.rs        # #[pymodule] entry
├── types.rs      # PyTrajectory (#[pyclass] wrapping inner: Trajectory)
├── io.rs
├── structure.rs
└── analysis.rs
```

Return types: `Vec<f64>`, `HashMap<String, Vec<f64>>` — PyO3 converts automatically to Python list/dict. No numpy or polars Rust crates.

```toml
# ferro-python/Cargo.toml — minimal deps
[dependencies]
pyo3 = { version = "0.21", features = ["extension-module"] }
ferro-core      = { path = "../ferro-core" }
ferro-io        = { path = "../ferro-io" }
ferro-structure = { path = "../ferro-structure" }
ferro-analysis  = { path = "../ferro-analysis" }
```

---

## CLI Binaries Reference

| Binary | Mode / flag | Purpose |
|---|---|---|
| `fe-convert` | `-i <in> -o <out>` | Format conversion |
| `fe-info` | `-i <file>` | Print structure info |
| `fe-job` | `-i <file> -s gaussian [-m B3LYP] [-b 6-31G*]` | Generate QC input file |
| `fe-traj` | `-m gr\|sq\|msd\|angle` | Structural analysis |
| `fe-corr` | `-m vacf\|rotcorr\|vanhove` | Correlation functions |
| `fe-cube` | `-m density\|velocity\|force\|sdf` | Spatial distribution maps |
| `fe-network` | `--P-O=2.3 [--format csv\|xlsx]` | Glass network analysis |

Common flags shared by all trajectory binaries: `--last-n N`, `--ncore N`, `--metal-units`.

**Batch input.** `-i` takes one or more files and expands glob patterns itself
(`-i 'runs/*/prod.dump'` — quote it so the shell leaves it alone; `{a,b}` braces are
*not* supported, let the shell do those). There is **one code path**: each input is
analysed independently and the results are stacked into a single CSV with a `file`
column. A single input is the N=1 case, not a special mode — dispatching on the number
of inputs would make the output shape depend on what a glob happened to match that day.

A failing input is reported, skipped, listed in the output's `[inputs]` block with its
reason, and makes the **exit code 1**. Parameter errors still fail fast, before the
first file is read.

**fe-traj help system** (two-level):
- `fe-traj` / `fe-traj -h` → overview: available modes + common flags
- `fe-traj -m <MODE>` (no `-i`) → mode-specific parameters + examples

### `fe-traj` — detailed flags

```
# type selection — two mutually exclusive groups
-a ELEM   [gr/sq] centre element; [angle] end atom A
-b ELEM   [gr/sq] neighbour element; [angle] centre atom B   (requires -a)
-c ELEM   [angle] end atom C   (requires -a -b)
-x LABEL  [gr/sq] centre site label; [angle] end atom A
-y LABEL  [gr/sq] neighbour site label; [angle] centre atom B  (requires -x)
-z LABEL  [angle] end atom C   (requires -x -y)

# gr-specific
--r-max FLOAT   max r [Å]                  default 10.005
                (clamped to half the smallest interplanar spacing)
--dr    FLOAT   bin width [Å]              default 0.01

# sq-specific
--q-max     FLOAT   max q [Å⁻¹]           default 25.0
--dq        FLOAT   q bin width [Å⁻¹]     default 0.05
--weighting ENUM    none|xrd|neutron|both  default both

# msd-specific
--dt    FLOAT   timestep [fs]              default 1.0
--shift INT     time-origin stride         default 1
--elements El,El,...  filter elements
--fit-range FMIN,FMAX  linear-fit window (traj. fractions) → self-diffusion D

# angle-specific
--r-cut-ab FLOAT   A→B cutoff [Å]         default 2.3
--r-cut-bc FLOAT   B→C cutoff [Å]         default 2.3
--d-angle  FLOAT   bin width [°]           default 0.1

# plot (gr / sq / angle / msd)
--plot   write PNG next to the data file (opened in a viewer only for a single input)
```

**`-o` is a name suffix, not a path**: results land in `<mode>[_<table>]_<suffix>.csv`
(`gr_run1.csv`, `network_qn_run1.csv`). Per-input files are not produced — the batch
summary is the product.

**element vs. label.** `Atom::element` is the chemical element; `Atom::label` is an
optional site type. The LAMMPS dump reader splits an `element` column entry written as
`<Element>_<suffix>` (`P_0`, `O_b_P_P`, `Zn_f`) at the **first underscore** — element
prefix into `element`, whole string into `label` — and prints the mapping once on read.
Plain symbols pass through untouched with `label = None`. `-a/-b/-c` select by element,
`-x/-y/-z` by label; only `ferro-io`'s CIF / CP2K / QE readers also populate `label`,
from their native site names.

### Output format

Every table-producing mode writes **one CSV**, with a `#` comment block above the data
holding the shared parameters and an `[inputs]` table (one row per input: frames, atoms,
volume, the clamped `r_max`, composition, status). That block is for humans —
`pandas.read_csv(comment="#")` drops it, so anything a script must parse is a column.

**The table's shape follows its primary product's granularity.** That is why `gr` and
`sq` differ:

**`gr` / `angle` — long format.** `file, r, center, neighbor, gr, cn` (and
`file, angle, end_a, center, end_c, count, p`). Types live in *data columns*, so runs
with different element sets stack with no column alignment; omitting `-a/-b` adds rows,
never columns. `g(r)` is symmetric (`A-B` == `B-A`); `CN(r)` is directed, which the
`center`/`neighbor` split writes into the table structure rather than a doc comment.
`angle` keeps a raw `count` alongside the normalised `p` because the integer histogram
is what `scripts/compare_angle.py` checks against dump2analysis bin for bin.

**`sq` — wide format.** `file, q, total_xrd, total_neutron`, then `_sq` / `_xrd` /
`_neutron` per pair (canonical half only — S(q) has no directed counterpart). The
primary product is the pair of totals, which are one value per `q`; the weighted
partials `w_ij(q)·S_ij(q)` are a diagnostic decomposition that sums back to the totals.
Inputs with different element sets contribute different pair columns: the union is
taken and **gaps are left empty (NaN), never interpolated or zero-filled**.

Label-resolved totals match element-resolved ones up to an O(1/N) term: same-type
partials normalise by `N_A(N_A−1)` whereas summing label partials reconstructs `N_A²`.

**Plot behaviour** — one PNG split into panels, **one panel per quantity, one curve per
input file**; colours are per file and consistent across panels, and only the first
panel carries a legend:
- `gr`: two panels (g(r) | CN(r)); needs a named pair, otherwise skipped with a note
- `sq`: two panels (XRD | Neutron) — different quantities, so no shared axis
- `angle`: one panel for the named triplet, mean ± std in the legend
- `msd`: 2×2 (total | a | b | c); `--fit-range` prints D = slope/6 and R² per input

`--plot` is **frozen at self-check quality** and will not chase matplotlib. The data is
long-format CSV, so a real figure is one line of seaborn
(`sns.lineplot(data=df, x="r", y="gr", hue="file")`); log axes, error bars and themes
belong in Python.

**`fe-cube` is the one exception** to the batch output rules: its product is a 3-D grid
file, so there is nothing to stack and nothing to plot. It still shares the traversal and
failure handling (error handling should not exist in two versions), but writes one
`.cube` per input, and the file name **must** carry the input stem
(`density_<stem>.cube`) or a second input would overwrite the first.

`fe-bader`'s `write_acf` / `write_bcf` / `write_avf` are likewise untouched: ACF.dat is
the Henkelman-group bader format that external tools parse.

### `fe-cube -m sdf` — Cluster SDF

Identifies Qn-type clusters (connected components of network-former atoms sharing bridging ligands), aligns each to a reference via Kabsch + permutation enumeration, and outputs per-atom-type 3D probability density as Gaussian cube files.

```
ferro-analysis/src/md/cube_sdf.rs
```

Key types:
```rust
pub struct ClusterSdfParams {
    pub former: String,              // network former element, e.g. "P"
    pub ligand: String,              // bridging ligand element, e.g. "O"
    pub target_qn: u8,               // 0/1/2/3 — determined by max individual Qn in component
    pub former_ligand_cutoff: f64,
    pub modifier: Option<String>,    // modifier cation, e.g. Some("Zn")
    pub modifier_cutoff: f64,
    pub grid_res: f64,               // Å/voxel
    pub sigma: f64,                  // Gaussian broadening in voxels
    pub padding: f64,                // grid boundary margin in Å
    pub rmsd_warn_threshold: f64,
}
```

Atom-type labels: `P0/P1/P2/P3` (individual Qn), `Of/On/Ob` (O connectivity), modifier element symbol.
Output: `<stem>_<label>.cube` per atom type (multi-family: `<stem>_fam<N>_<label>.cube`).

### `ferro-analysis/src/md/` file structure

| File | Analysis |
|---|---|
| `gr.rs` | Radial distribution function g(r) + coordination number CN(r) |
| `sq.rs` | Structure factor S(q) via Fourier transform of g(r) |
| `msd.rs` | Mean square displacement (time-shift average, NPT-safe) |
| `angle.rs` | Bond angle distribution P(θ) |
| `vanhove.rs` | Van Hove self-correlation Gs(r, τ) |
| `vacf.rs` | Velocity autocorrelation function |
| `rotcorr.rs` | Rotational correlation C₂(t) for molecular bond vectors |
| `cube_density.rs` | 3D spatial density / velocity / force maps |
| `cube_sdf.rs` | Cluster SDF — Kabsch alignment + per-type density accumulation |

---

## Extending the Project

### Add a file format
1. `ferro-io/src/readers/<fmt>.rs` — return `Result<Trajectory>`
2. `ferro-io/src/writers/<fmt>.rs`
3. Export from `readers/mod.rs`, `writers/mod.rs`
4. Add format detection in `ferro-cli/src/commands/io.rs`
5. Add wrapper in `ferro-python/src/io.rs`

### Add a structure operation
1. Implement in `ferro-structure/src/` — takes/returns `Trajectory`
2. `ferro-cli/src/commands/structure.rs`
3. `ferro-python/src/structure.rs`

### Add an analysis method
1. Implement in `ferro-analysis/src/<domain>/`
2. `ferro-cli/src/commands/analysis.rs`
3. `ferro-python/src/analysis.rs`

### Add a QC software target
1. Builder in `ferro-workflow/src/job_builder.rs`
2. Templates in `ferro-workflow/src/templates.rs`
3. CLI branch in `ferro-cli/src/commands/`

---

## Key Dependencies

| Crate | Purpose |
|---|---|
| `nalgebra` | 3D vectors/matrices; coordinates are `nalgebra::Vector3<f64>`, cell is `Matrix3<f64>` |
| `ndarray` | Multi-dimensional arrays for bulk trajectory data |
| `rayon` | Data parallelism |
| `thiserror` | Derive macros for `ChemError` |
| `anyhow` | Error propagation in CLI |
| `clap` | CLI argument parsing |
| `rustyline` | Readline-style input for REPL (to be added to workspace deps) |
| `pyo3` | Python bindings (ferro-python only) |
