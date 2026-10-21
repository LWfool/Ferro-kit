# ferro-python Overview

`ferro-python` is the Python side of ferro: a package named `ferro` that gives
`import ferro` in Python and the `ferro-py` command on the shell. This page says
what it is made of and why it is split the way it is. How to use it is in its own
manual, compiled into the package:

```bash
ferro-py doc              # list of pages: api, config, job, job vasp|qe|cp2k, conv, plot
```

The pages are deliberately not repeated here — `ferro doc` and `ferro-py doc`
each carry their own manual, so neither can drift from a copy in the other.

## Basic facts

| | |
|---|---|
| Package | `ferro` (wheel), version taken from `ferro-python/Cargo.toml` |
| Python | CPython ≥ 3.11, one **abi3** wheel per platform for every version |
| Wheels | `manylinux2014` x86_64 (glibc 2.17: CentOS 7 and later) and macOS arm64, built by GitHub Actions and tested on Python 3.11 and 3.13 |
| Install | `pip install --no-index ferro-*.whl` — nothing else, so an offline cluster works |
| Gives | `import ferro`; the `ferro-py` command; the manual, `ferro-py doc` |
| Does not need | the Rust `ferro` binary, the source tree, any network access |
| Source | `ferro-python/`: a Cargo workspace of its own, built with maturin, not by the root `cargo build` |

## Dependencies

| Side | Dependency | Why |
|---|---|---|
| Rust | `ferro-core`, `ferro-io`, `ferro-structure`, `ferro-analysis` | data model, file formats, structure operations, analysis — the same code as the `ferro` binary |
| Rust | `ferro-render` | renders the manual for `ferro-py doc`, as it does for `ferro doc` |
| Rust | `pyo3` 0.29 (`abi3-py311`) | the Python binding |
| Python, at run time | the standard library only (`tomllib`, `argparse`, `subprocess`, …) | clusters are offline; ASE, pymatgen and cp2k-input-tools all pull in SciPy, Matplotlib or pydantic |
| Python, optional | `ferro[plot]`: matplotlib, pandas, scienceplots (plus a LaTeX installation) | only for `ferro-py plot` / `ferro.plot`; `import ferro` does not load them |
| Build | `maturin` ≥ 1.0, < 2.0 | builds the mixed Rust / Python wheel |
| Tests | `pytest` | development only, not a dependency of the wheel |

Input generation lives only in Python (below): the Rust `ferro job` and its
`ferro-workflow` crate were removed in 0.4.0.

## Responsibilities

The split follows one rule: **a file format is specified once, in Rust; what
changes often and what users want to edit is in Python.**

| Rust extension `ferro._ferro` | Python package `ferro` |
|---|---|
| reading and writing structures and trajectories (`read`, `write`, the `Trajectory` type) | `ferro.jobs`: input generation for VASP, QE and CP2K (`sp`, `opt`, `aimd`) |
| the structure blocks of input files: POSCAR, QE's three cards, CP2K's `&CELL` / `&COORD` | the parameter blocks: INCAR, the pw.x namelists, the CP2K section tree |
| parsing and writing QE namelists (one grammar for the reader and the generator) | the templates: nine native input files, users' own templates by name |
| reading program outputs: energy, forces, stress (`read_output`) and writing ferro's tables (`write_table`) | convergence tests and their analysis arithmetic (`ferro-py conv`) |
| structure operations (`supercell`, `add_vacuum_layer`, `merge`) and analysis (`gr_pair`, `gr_all`, `msd`) | the configuration file, `~/.config/ferro/config.toml` |
| the manual pages and their terminal rendering (`ferro-py doc`) | the `ferro-py` command line and its help pages |
| | `ferro.plot`: figures from ferro's csv tables (`ferro-py plot`) |

Input generation moved to Python because it is not performance-bound, changes
more often than anything else, and a user can read and edit it without a Rust
toolchain. What it cannot afford to get wrong twice — how a POSCAR or a QE card is
laid out, how an OUTCAR is read — stays in the one Rust implementation that the
`ferro` binary also uses.

## Architecture

```
  shell: ferro-py job | conv | config | plot | doc Python: import ferro
              │                                            │
              ▼                                            │
  ferro/cli.py ─── ferro/config.py                         │
              │                                            │
              ▼                                            │
  ferro/jobs/                                              │
    base.py     Generator → InputSet; options,             │
                multiplicity, k-mesh, batch writing        │
    vasp.py  qe.py  cp2k.py   one Generator per program    │
    conv.py     convergence tests                          │
    templates/  <program>/<task>/<input file>              │
  ferro/plot/   gr sq angle msd net (optional matplotlib)  │
              │                                            │
              ▼                                            ▼
  ferro._ferro  (Rust, PyO3; ferro/__init__.py re-exports it)
    types.rs  io.rs  jobs.rs  structure.rs  analysis.rs  doc.rs
              │
              ▼
  ferro-core · ferro-io · ferro-structure · ferro-analysis · ferro-render
```

- `ferro/__init__.py` re-exports the extension, so `import ferro` exposes the Rust
  functions directly; `ferro.jobs` is pure Python on top of them.
- A generator turns a template, one frame of a trajectory and a few options into
  an `InputSet` (the files of one job directory). Adding a program means adding a
  generator and its templates; adding a feature means adding an option that each
  generator applies.
- `ferro-python` is, with `ferro-cli`, the only layer that combines several ferro
  crates; the crates below it know nothing of Python.

## Where to read next

| Question | Page |
|---|---|
| Python API: `read`, `Trajectory`, analysis functions | `ferro-py doc api` (also `docs/src/python/api.md`) |
| Configuration: POTCARs, UPF files, CP2K data | `ferro-py doc config` |
| Generating inputs | `ferro-py doc job`, then `job vasp`, `job qe`, `job cp2k` |
| Convergence tests | `ferro-py doc conv` |
| Figures from ferro's tables | `ferro-py doc plot` |
