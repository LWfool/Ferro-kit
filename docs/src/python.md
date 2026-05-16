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

## Quick start

```python
import ferro

t = ferro.read("traj.lammpstrj", metal_units=True)   # -> Trajectory
print(len(t), t.n_atoms(), t.elements())

sc = ferro.supercell(t, 2, 2, 1)
ferro.write(sc, "POSCAR")

g = ferro.gr(t, r_max=10.0, dr=0.02)        # dict[str, list[float]]
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

### `gr(traj, r_max=None, dr=0.01, r_cut=2.3, r_min=0.005)`

Radial distribution function and coordination number.  `r_max=None`
auto-selects half the shortest cell vector of the first frame.

Returned keys:

- `"r"` — bin-centre radii [Å]
- `"gr:<El1-El2>"`, `"gr:total"` — partial / total g(r)
- `"cn:<center-neighbor>"` — directed cumulative CN(r)

```python
g = ferro.gr(t, r_max=8.0, dr=0.05)
import matplotlib.pyplot as plt
plt.plot(g["r"], g["gr:total"])
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
