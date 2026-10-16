# Python Bindings

`ferro` ships first-class Python bindings (PyO3) covering the core workflow:
**read → transform → analyze**.  The pure-Rust crates have zero Python
awareness; the bindings live in the separate `ferro-python` crate, which is its
own Cargo workspace and is built with [maturin](https://www.maturin.rs/) rather
than `cargo build`.

## Install

The package is one wheel per platform (abi3: any CPython ≥ 3.11) with the Rust
code compiled in and nothing else to install, so it runs on an offline cluster.
GitHub Actions builds `manylinux2014` x86_64 (CentOS 7 and later) and macOS arm64
wheels; download the artifact and

```bash
pip install --no-index ferro-*.whl      # `import ferro` and the ferro-py command
```

To build one yourself:

```bash
pip install maturin
cd ferro-python
maturin build --release                 # -> target/wheels/ferro-*.whl
```

> `cargo build` from the workspace root does **not** build this crate: it is a
> separate workspace with its own lockfile, so the root `cargo build` / `test` /
> `clippy` skip it, and a break in the bindings is not caught by the main test
> suite. After changing a public API of `ferro-core`, `ferro-io` or
> `ferro-analysis`, run `cd ferro-python && cargo check`.

The command-line side — `ferro-py job`, `conv`, `config` — has its own pages:
run `ferro-py doc` for the list.

## Quick start

```python
import ferro

t = ferro.read("traj.lammpstrj", units="metal")   # -> Trajectory
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
| `pbc(frame)` | `(bool,bool,bool)` | Periodicity along a, b, c |
| `energy(frame)` | `float \| None` | Total energy [eV]; None for structure files |
| `forces(frame)` | `list[(float,float,float)] \| None` | Forces [eV/Å], same order as `positions` |
| `stress(frame)` | `list[list[float]] \| None` | Stress tensor [eV/Å³], row-major, **positive = compression** (the sign CP2K, VASP and QE print) |
| `subset(indices)` | `Trajectory` | New trajectory with the given frames, in the given order |
| `frame_indices(start=None, end=None, stride=None, number=None)` | `list[int]` | The frame choice of `ferro convert`: 0-based, `end` inclusive, all frames by default |
| `len(t)` | `int` | = `n_frames()` |

Frame indices are 0-based; out-of-range raises `IndexError`.

## I/O

### `read(path, units=None, atom_style=None, type_map=None) -> Trajectory`

Format auto-detected from the file name / extension:

| Extension / name | Format |
|---|---|
| `.xyz` / `.extxyz` | XYZ / extended XYZ |
| `.pdb` / `.cif` | PDB / CIF |
| `POSCAR*` / `CONTCAR*` / `.vasp` / `.poscar` / `.pos` | VASP |
| `.in` / `.qe` | Quantum ESPRESSO input |
| `.inp` / `.restart` | CP2K input / restart |
| `.lammpstrj` / `.dump` | LAMMPS dump; **`units` (`"real"` or `"metal"`) is required when it has velocity/force columns**, since a dump does not record its unit system |
| `.data` / `.lmp` / `.lammps` | LAMMPS data; **`atom_style` (`"atomic"`, `"charge"` or `"full"`) is required**, never guessed from the `Atoms` comment or the column count |

A dump without an element column reads its atoms as `X1`, `X2`, …; `type_map={1: "Zr", 2: "O"}`
names them, exactly as `--type-map` on the command line.  A type missing from the map, an
unknown symbol, or `type_map` on anything but a dump is an error.

### `write(traj, path, units=None)`

Writes by extension: `xyz`, `extxyz`, `pdb`, `cif`, `POSCAR*`/`CONTCAR*`/`vasp`/`poscar`/`pos` (POSCAR format),
`data`/`lmp`/`lammps`, `lammpstrj`/`dump`.  QE input (`in`/`qe`) is read-only.  Writing velocities or forces to a dump requires `units`.

### `read_output(path) -> (Trajectory, dict)`

Reads a calculation output — VASP `OUTCAR` / `vasprun.xml`, CP2K single-point or MD output —
recognised by its content, as `ferro dataset collect` does.  Frames carry energy, forces and
stress.  Frames that were dropped (SCF not converged, incomplete) are not in the trajectory; the
dict counts them (`scf_failed`, `incomplete`, …) next to `format`, `steps`, `kept` and `version`.

### `write_table(path, columns, meta=[])`

Writes a CSV in the layout of every other ferro product: `# `-prefixed `meta` lines, a header,
the data; NaN becomes an empty field.  `columns` is a list of `(name, values)`, each column all
numbers or all strings, all of the same length.

### Building blocks for input generation

These return text for a calculation input; atoms are written in the order they are in, so call
`sort_atoms` first to get the (Z, symbol) order that `ferro convert` writes.

| Function | Description |
|---|---|
| `sort_atoms(traj)` | Per frame, atoms in (Z, symbol) order → new `Trajectory` |
| `poscar_string(traj, frame_index=0, comment="ferro")` | POSCAR text of one frame |
| `qe_structure_cards(traj, pseudos, frame_index=0)` | pw.x `ATOMIC_SPECIES` / `CELL_PARAMETERS` / `ATOMIC_POSITIONS`; `pseudos` maps every element to a pseudopotential file name.  A frame without a cell is an error |
| `cp2k_structure_sections(traj, frame_index=0)` | CP2K `&CELL` + `&COORD`; a frame without a cell is an error |
| `parse_qe_input_parts(text)` | pw.x input → `([(namelist, [(key, raw value)])], cards text)`; the same parser `read` uses |
| `write_namelists(namelists)` | The inverse: namelists back to text |

## Structure operations

| Function | Description |
|---|---|
| `supercell(traj, nx, ny, nz)` | Per-frame supercell → new `Trajectory` |
| `add_vacuum_layer(traj, axis, thickness)` | Add vacuum along the `"x"`/`"y"`/`"z"` lattice vector; `thickness` is the **perpendicular** gap (a tilted axis is lengthened by thickness / cos θ, as in ASE).  Atoms are not moved or centred.  A singular cell (c = 0 for a 2D material) or a non-finite thickness is an error |
| `merge(a, b, axis, gap)` | Stack the first frames of two trajectories along `axis`: neither is strained, both are centred in the interface plane, `gap` is measured along the interface normal (as in `add_vacuum`); the interfaces must be parallel. Atoms are folded into their own cell first; one outside it along a non-periodic axis is an error |

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

`gr_all` gives $n^2$ ordered pairs for n types.  There is **no total**: an unweighted
total is the degenerate `f_i ≡ 1` case, corresponds to no experimental probe, and
was removed in 0.1.11.

```python
g = ferro.gr_pair(t, "P", "O", r_max=8.0, dr=0.05)
import matplotlib.pyplot as plt
plt.plot(g["r"], g["P-O_gr"])

# group by site label (labelled trajectory)
g = ferro.gr_pair(t, "P_3", "O_b", by="label", r_max=5.0)
```

### `msd(traj, dt, max_lag=None, elements=None)`

Mean squared displacement: every lag averaged over all time origins (FFT), periodic systems unwrapped
with the TOR scheme.  `dt` is the time between stored frames in fs (MD timestep × dump interval) and
is required, as `--dt` is on the command line: a guessed default rescales the whole time axis
silently.  A trajectory that carries step numbers must be evenly spaced with no repeated frames,
otherwise a `RuntimeError` names the offending pair.  `max_lag` in frames, default half the trajectory.  Details in
[MSD](analysis/msd.md).

Returned keys: `"time"` [fs], `"msd"` (total), `"msd_x"`, `"msd_y"`, `"msd_z"` (Cartesian components,
which add up to the total).

```python
d = ferro.msd(t, dt=2.0, elements=["Li"])
```

## Errors

Rust errors (I/O failures, missing cells, empty trajectories) surface as Python
exceptions (`RuntimeError` / `IndexError` / `ValueError`) with the underlying
message preserved.
