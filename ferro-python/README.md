# ferro (Python bindings)

PyO3 bindings for the **ferro** computational-chemistry toolkit (MD trajectory
post-processing for periodic systems). The pure-Rust crates have zero Python
awareness; this package is the only PyO3 glue layer.

## Build / install

Requires [maturin](https://www.maturin.rs/) (`pip install maturin`).
`cargo build` from the workspace root does **not** apply — this crate is its
own workspace and must be built by maturin:

```bash
cd ferro-python
maturin develop            # build + install into the active venv
# or
maturin build --release    # produce a wheel in target/wheels/
```

## Usage

```python
import ferro

t = ferro.read("traj.lammpstrj", metal_units=True)   # -> Trajectory
print(len(t), t.n_atoms(), t.elements())

sc = ferro.supercell(t, 2, 2, 1)
ferro.write(sc, "POSCAR")

g = ferro.gr(t, r_max=10.0, dr=0.02)        # dict[str, list[float]]:
#   "r", "gr:<El1-El2>", "gr:total", "cn:<center-neighbor>"

d = ferro.msd(t, dt=2.0, elements=["Li"])   # "time","msd","msd_a/b/c"
```

### API surface

| Function | Description |
|---|---|
| `read(path, metal_units=False)` | Auto-detect format → `Trajectory` |
| `write(traj, path, metal_units=False)` | Write by extension |
| `supercell(traj, nx, ny, nz)` | Per-frame supercell |
| `add_vacuum_layer(traj, axis, thickness)` | Add vacuum along `x`/`y`/`z` |
| `merge(a, b, axis, gap)` | Merge first frames |
| `gr(traj, r_max=None, dr=0.01, r_cut=2.3, r_min=0.005)` | g(r) + CN(r) |
| `msd(traj, dt=1.0, shift=1, tau=None, elements=None)` | MSD(t) |

`Trajectory` methods: `n_frames()`, `n_atoms()`, `elements()`,
`positions(frame)`, `symbols(frame)`, `cell(frame)`, `charge(frame)`,
`multiplicity(frame)`, `tail(n)`, `len()`.

Analysis functions return `dict[str, list[float]]` (PyO3 auto-converts).
