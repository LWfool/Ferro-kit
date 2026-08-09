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

# CLI (dev) — one binary, subcommands grouped by what they produce
cargo run --bin ferro -- traj gr -i traj.dump -a P -b O
cargo run --bin ferro -- traj gr -i 'runs/*/prod.dump' -a P -b O -o scan   # batch
cargo run --bin ferro -- map density -i traj.dump
cargo run --bin ferro -- net qn -i traj.dump --P-O=2.3
cargo run --bin ferro -- convert -i input.xyz -o output.pdb
cargo run --bin ferro -- job -i input.xyz -s gaussian -m B3LYP -o job.gjf

# Python bindings (separate workspace, requires maturin; currently does not compile)
cd ferro-python && maturin develop
```

Test fixtures: `tests/` (two 5-frame LAMMPS trajectories — one NPT, one NVT — plus `Experiment.xlsx`)

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
├── table.rs         # Table / Column — analysis results on their way to ferro-io
├── charge_grid.rs   # ChargeGrid (Bader)
├── cube_data.rs     # CubeData (3-D grids)
├── cluster.rs       # build_network_graph / connected_components (shared primitive)
├── network_type.rs  # TypeParams / classify_frame
├── spin.rs          # guess_spin / oxidation states
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

**Implemented:** one-shot CLI only.

```bash
ferro traj gr -i traj.dump -a P -b O          # one input
ferro traj gr -i 'runs/*/prod.dump' -o scan   # many, stacked into one csv
```

**Not implemented** — `main.rs` is a subcommand dispatcher; bare `ferro` prints the
overview. The REPL / script interpreter (`rustyline`, `ferro -f workflow.mf`, piped
stdin, bare `ferro` entering an interactive session on a tty) is the next large feature;
see `dev/plan.md`. It must link every command into one process — a stateful
`read` → `gr` → `sq` session cannot shell out per command without re-reading the
trajectory — which is also why 0.2.0 collapsed the `fe-*` binaries into one.

**Python library:** `import ferro` via PyO3 in `ferro-python` — **currently does not
compile**, see the Python Bindings section below.

### ferro-cli internal structure
```
ferro-cli/src/
├── main.rs              # the ferro binary: subcommand tree + dispatch
├── lib.rs               # re-exports args, batch, cmd, help, io_dispatch, plot
├── args/
│   ├── common.rs        # CommonArgs (-i/-o/--last-n/--ncore/--metal-units) + SelectArgs (-a..-z)
│   ├── traj.rs          # SqWeightingCli
│   ├── corr.rs
│   └── cube.rs          # CubeCliMode (still the analysis-facing enum)
├── batch.rs             # multi-input driving, generic over the result type
├── cmd/
│   ├── traj.rs          # gr/sq/msd/angle/vacf/rotcorr/vanhove — one Args struct each
│   ├── map.rs           # density/velocity/force/radius/sdf/chg-sdf
│   ├── net.rs           # qn/type + the --P-O=2.3 argv pre-pass
│   ├── bader.rs, convert.rs, info.rs, job.rs
│   └── mod.rs
├── help.rs              # print_overview / print_*_overview / per-command help text
├── io_dispatch.rs       # read_trajectory / read_trajectory_tail (format detection)
└── plot.rs              # Panel/Series + render; plot_gr / plot_sq / plot_msd / plot_angle
```

`batch.rs` is the multi-input driver and is **generic over the result type** —
`map_inputs<T>(inputs, f) -> (Vec<(PathBuf, T)>, Vec<Failure>)`. It owns indexing,
looping, skipping and reporting; it never names `GrResult` or any other analysis type,
and the analysis crates never learn that batching exists. Deciding *which* scalars are
worth summarising is analysis semantics, so that stays in `cmd/`.


---

## Python Bindings (ferro-python)

> **Currently broken.** `analysis.rs` still calls the `write_gr` path that 0.2.0 removed,
> so `cd ferro-python && cargo check` fails. `ferro-python` is a **separate workspace**
> with its own lockfile, so the main workspace's `cargo build`/`test`/`clippy` skip it
> entirely — that is exactly how it drifted. Run its `cargo check` after touching any
> public API in `ferro-analysis` or `ferro-core`.

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
pyo3 = { version = "0.29", features = ["extension-module"] }
ferro-core      = { path = "../ferro-core" }
ferro-io        = { path = "../ferro-io" }
ferro-structure = { path = "../ferro-structure" }
ferro-analysis  = { path = "../ferro-analysis" }
```

---

## CLI Reference

**One binary, `ferro`.** Subcommands are grouped by **what they produce**, not by which
crate implements them — the `fe-*` binaries are gone.

| Command | Purpose | Product |
|---|---|---|
| `ferro traj gr\|sq\|msd\|angle\|vacf\|rotcorr\|vanhove` | Trajectory analysis | one stacked CSV + optional PNG |
| `ferro map density\|velocity\|force\|radius\|sdf\|chg-sdf` | Spatial distribution maps | one `.cube` grid **per input** |
| `ferro net qn\|type` | Glass network topology (`--P-O=2.3`) | distribution tables / labelled structure |
| `ferro bader` | Bader charge partitioning | ACF/BCF/AVF (Henkelman format) |
| `ferro convert` / `info` / `job` | Structure I/O, info, QC input files | files / stdout |

`traj` merges the old `ferro traj` + `ferro traj`: all seven share one output pipeline, and the
"structural vs dynamic" line never held anyway (`msd` is a time correlation that sat on
the structural side). `-m` is gone everywhere — position is the mode. Each subcommand
owns its own argument struct, so `--dt` can mean different things for `msd` and `vacf`
without sharing a field; `CommonArgs` / `SelectArgs` are `#[command(flatten)]`-ed in.

**Help is three levels**: `ferro` → the groups; `ferro traj` → that group's commands;
`ferro traj gr` (no `-i`) → that command's parameters, output columns and examples.

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

### `ferro traj` — detailed flags

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

**`ferro map` is the one exception** to the batch output rules: its product is a 3-D grid
file, so there is nothing to stack and nothing to plot. It still shares the traversal and
failure handling (error handling should not exist in two versions), but writes one
`.cube` per input, and the file name **must** carry the input stem
(`density_<stem>.cube`) or a second input would overwrite the first.

`ferro bader`'s `write_acf` / `write_bcf` / `write_avf` are likewise untouched: ACF.dat is
the Henkelman-group bader format that external tools parse.

### `ferro map sdf` — Cluster SDF

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
4. Add format detection in `ferro-cli/src/io_dispatch.rs` (both `read_trajectory` and `write_trajectory`)
5. Add wrapper in `ferro-python/src/io.rs`

### Add a structure operation
1. Implement in `ferro-structure/src/` — takes/returns `Trajectory`
2. `ferro-python/src/structure.rs`
3. No CLI surface yet (`box_builder` is library-only too)

### Add an analysis method
1. Implement in `ferro-analysis/src/<domain>/` — **pure computation, no file I/O**
2. Give the result type `to_tables() -> Vec<(String, Table)>` and `meta_lines() -> Vec<String>`.
   Put only batch-wide parameters in `meta_lines`; anything that varies per input belongs
   in the `[inputs]` summary, or it will masquerade as a global fact
3. Add a subcommand arm in `ferro-cli/src/cmd/<group>.rs`: build params (validating before
   any file is read) → `batch::map_inputs` → `batch::stack` → `batch::write_all`
4. Add the help text in `ferro-cli/src/help.rs` and list it in `print_overview`
5. `ferro-python/src/analysis.rs`

### Add a QC software target
1. Builder in `ferro-workflow/src/job_builder.rs`
2. Templates in `ferro-workflow/src/templates.rs`
3. CLI branch in `ferro-cli/src/cmd/job.rs`

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
| `glob` | Input pattern expansion for batch runs (`ferro-cli`) |
| `rustyline` | REPL input — **not a dependency yet**; add when the REPL lands |
| `pyo3` | Python bindings (ferro-python only) |
