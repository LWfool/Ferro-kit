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
maturin build --interpreter "$(which python)"   # wheel → target/wheels/
pip install target/wheels/*.whl
# or, with a venv/conda env already active:
maturin develop
```

Pin the interpreter explicitly: without `--interpreter`, maturin may pick a stray
`python3.9` and tag the wheel for the wrong version.

## Usage

```python
import ferro

t = ferro.read("traj.lammpstrj", metal_units=True)   # -> Trajectory
print(len(t), t.n_atoms(), t.elements())

sc = ferro.supercell(t, 2, 2, 1)
ferro.write(sc, "POSCAR")

# one pair -> {"r", "P-O_gr", "P-O_cn"}
g = ferro.gr_pair(t, "P", "O", r_max=10.0, dr=0.02)

# every ordered pair in one pass -> {"r", "<A>-<B>_gr", "<A>-<B>_cn", ...}
ga = ferro.gr_all(t, by="label")

d = ferro.msd(t, dt=2.0, elements=["Li"])            # "time","msd","msd_a/b/c"
```

### API surface

| Function | Description |
|---|---|
| `read(path, metal_units=False)` | Auto-detect format → `Trajectory` |
| `write(traj, path, metal_units=False)` | Write by extension |
| `supercell(traj, nx, ny, nz)` | Per-frame supercell |
| `add_vacuum_layer(traj, axis, thickness)` | Add vacuum along `x`/`y`/`z` |
| `merge(a, b, axis, gap)` | Merge the first frames of two systems |
| `gr_pair(traj, a, b, by="element", r_max=None, dr=0.01, r_min=0.005)` | g(r) + CN(r) for one pair |
| `gr_all(traj, by="element", r_max=None, dr=0.01, r_min=0.005)` | g(r) + CN(r) for all n² ordered pairs |
| `msd(traj, dt=1.0, shift=1, tau=None, elements=None)` | MSD(t), time-origin averaged |

`Trajectory` methods: `n_frames()`, `n_atoms()`, `elements()`,
`positions(frame)`, `symbols(frame)`, `cell(frame)`, `charge(frame)`,
`multiplicity(frame)`, `tail(n)`, `len()`, `repr()`.

Analysis functions return `dict[str, list[float]]` (PyO3 auto-converts).
`by="label"` groups by site label instead of element — valid on a single frame only,
since g(r) requires a fixed particle count per type.

Network analysis (`ferro net`) is not wrapped yet.

## Notes

- `__version__` is taken from `Cargo.toml`; keep it in step with the root workspace
  version by hand — this crate is not a member of it
- `gr(...)` was split into `gr_pair` / `gr_all` in 0.1.10; the old name is gone
