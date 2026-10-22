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

The wheel is **abi3** (CPython ≥ 3.11): one wheel per platform serves every
Python version from 3.11 on. The package mixes Rust and Python: the Rust
extension is `ferro._ferro`, pure-Python modules live in `python/ferro/`, and
`import ferro` re-exports the extension, so the usage below is unchanged.

## Usage

```python
import ferro

t = ferro.read("traj.lammpstrj", units="metal")   # -> Trajectory
print(len(t), t.n_atoms(), t.elements())

sc = ferro.supercell(t, 2, 2, 1)
ferro.write(sc, "POSCAR")

# one pair -> {"r", "P-O_gr", "P-O_cn"}
g = ferro.gr_pair(t, "P", "O", r_max=10.0, dr=0.02)

# every ordered pair in one pass -> {"r", "<A>-<B>_gr", "<A>-<B>_cn", ...}
ga = ferro.gr_all(t, by="label")

d = ferro.msd(t, dt=2.0, elements=["Li"])            # "time","msd","msd_x/y/z"
```

The full API — every function, the `Trajectory` methods, the warning classes — is in
[`docs/src/python/api.md`](../docs/src/python/api.md), also readable as `ferro-py doc api`;
the `ferro-py job` / `conv` / `config` / `plot` commands have their own pages there.

Tests: `pytest ferro-python/tests` with the built wheel on the path.

## Notes

- `__version__` is taken from `Cargo.toml`; keep it in step with the root workspace
  version by hand — this crate is not a member of it
