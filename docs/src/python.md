# Python Bindings

`ferro` ships first-class Python bindings (PyO3) covering the core workflow:
**read → transform → analyze**.  The pure-Rust crates have zero Python
awareness; the bindings live in the separate `ferro-python` crate, which is its
own Cargo workspace and is built with [maturin](https://www.maturin.rs/) rather
than `cargo build`.

## Install

```bash
pip install maturin
cd ferro-python
maturin develop                 # build + install into the active venv/conda env
# or, without an active environment:
maturin build --release --interpreter "$(which python)"
pip install target/wheels/ferro-*.whl
```

> `cargo build` from the workspace root does **not** build this crate (and a
> standalone `cargo build` fails at the macOS link step on undefined Python
> symbols — that is expected; maturin supplies the correct link flags).
>
> Because it is a separate workspace with its own lockfile, the root
> `cargo build` / `test` / `clippy` skip it entirely — a break in the bindings is
> not caught by the main test suite.  After changing any public API in
> `ferro-core` or `ferro-analysis`, run `cd ferro-python && cargo check`.
>
> **Status (0.3.0):** type-level checks pass cleanly (`cargo clean && cargo check`,
> clippy zero warnings).  The pyo3 0.21 → 0.29 upgrade has **not** been verified at
> runtime — this machine has no maturin, so module initialisation and the
> `#[pyfunction]` signatures are unconfirmed in practice.

## Quick start

```python
import ferro

t = ferro.read("traj.lammpstrj", metal_units=True)   # -> Trajectory
print(len(t), t.n_atoms(), t.elements())

sc = ferro.supercell(t, 2, 2, 1)
ferro.write(sc, "POSCAR")

g = ferro.gr_pair(t, "P", "O", r_max=10.0, dr=0.02)   # dict[str, list[float]]
d = ferro.msd(t, dt=2.0, elements=["Li"])
```

## `Trajectory`

Returned by `ferro.read`; wraps the Rust `Trajectory` (single-frame files use
this type too).

| Method | Returns | Description |
|---|---|---|
| `n_frames()` | `int` | Number of frames |
| `n_atoms()` | `int \| None` | Atom count (None if empty) |
| `elements()` | `list[str]` | First-frame elements, first-appearance order |
| `positions(frame)` | `list[(float,float,float)]` | Cartesian coords [Å] |
| `symbols(frame)` | `list[str]` | Element symbols (matches `positions`) |
| `cell(frame)` | `list[list[float]] \| None` | 3×3 cell rows [Å]; None if non-periodic |
| `charge(frame)` | `int` | Total charge |
| `multiplicity(frame)` | `int` | Spin multiplicity 2S+1 |
| `tail(n)` | `Trajectory` | New trajectory with the last `n` frames |
| `len(t)` | `int` | = `n_frames()` |

Frame indices are 0-based; out-of-range raises `IndexError`.

## I/O

### `read(path, metal_units=False) -> Trajectory`

Format auto-detected from the file name / extension:

| Extension / name | Format |
|---|---|
| `.xyz` / `.extxyz` | XYZ / extended XYZ |
| `.pdb` / `.cif` | PDB / CIF |
| `POSCAR*` / `CONTCAR*` | VASP |
| `.in` / `.qe` | Quantum ESPRESSO input |
| `.inp` / `.restart` | CP2K input / restart |
| `.lammpstrj` / `.dump` / `.lammps` | LAMMPS dump (`metal_units` switches real↔metal) |
| `.data` / `.lmp` | LAMMPS data |

### `write(traj, path, metal_units=False)`

Writes by extension: `xyz`, `extxyz`, `pdb`, `cif`, `POSCAR`, `in`/`qe`,
`data`/`lmp`, `lammpstrj`/`dump`.

## Structure operations

| Function | Description |
|---|---|
| `supercell(traj, nx, ny, nz)` | Per-frame supercell → new `Trajectory` |
| `add_vacuum_layer(traj, axis, thickness)` | Add vacuum along `"x"`/`"y"`/`"z"` |
| `merge(a, b, axis, gap)` | Merge first frames of two trajectories |

## Analysis

Analysis functions return `dict[str, list[float]]` (PyO3 converts the Rust
`HashMap<String, Vec<f64>>` automatically — no NumPy dependency).

### `gr_pair(traj, a, b, by="element", r_max=None, dr=0.01, r_min=0.005)`
### `gr_all(traj, by="element", r_max=None, dr=0.01, r_min=0.005)`

Radial distribution function and coordination number.  `gr_pair` computes one
named pair, `gr_all` every ordered pair in a single pass.  `r_max=None`
auto-selects the minimum-image bound (half the smallest **interplanar spacing**,
not the shortest cell vector — the two differ for non-orthogonal cells).

`by="element"` groups on `Atom::element`, `by="label"` on `Atom::label`.

Returned keys:

- `"r"` — bin-centre radii [Å]
- `"<A>-<B>_gr"` — partial g(r); symmetric, `A-B` equals `B-A`
- `"<A>-<B>_cn"` — directed cumulative CN(A→B)

`gr_all` gives n² ordered pairs for n types.  There is **no total**: an unweighted
total is the degenerate `f_i ≡ 1` case, corresponds to no experimental probe, and
was removed in 0.1.11.

```python
g = ferro.gr_pair(t, "P", "O", r_max=8.0, dr=0.05)
import matplotlib.pyplot as plt
plt.plot(g["r"], g["P-O_gr"])

# 按位点标签分组（标注轨迹）
g = ferro.gr_pair(t, "P_3", "O_b", by="label", r_max=5.0)
```

### `msd(traj, dt=1.0, shift=1, tau=None, elements=None)`

Mean squared displacement (time-origin averaged, NPT-safe).

Returned keys: `"time"` [fs], `"msd"` (total), `"msd_a"`, `"msd_b"`, `"msd_c"`
(crystal axes for periodic systems, x/y/z otherwise).

```python
d = ferro.msd(t, dt=2.0, elements=["Li"])
```

## Errors

Rust errors (I/O failures, missing cells, empty trajectories) surface as Python
exceptions (`RuntimeError` / `IndexError` / `ValueError`) with the underlying
message preserved.
