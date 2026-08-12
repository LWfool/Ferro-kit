# Installation

## Requirements

- Rust 1.75+ (`rustup` recommended)
- Cargo (bundled with Rust)

## Build from Source

```bash
git clone https://github.com/ferro/ferro
cd ferro
cargo build --release
```

The build produces **one** binary, `target/release/ferro`. (Before 0.2.0 there were
eight `fe-*` binaries; they are gone, with no compatibility shims — see
[CLI Reference](cli-reference.md) for the migration table.)

| Command group | Purpose | Product |
|---|---|---|
| `ferro traj gr\|sq\|msd\|angle\|vacf\|rotcorr\|vanhove` | Trajectory analysis | one stacked CSV + optional PNG |
| `ferro map density\|velocity\|force\|radius\|sdf\|chg-sdf` | 3-D spatial maps | one `.cube` per input |
| `ferro net` | Glass network topology | five stacked CSVs + optional labelled trajectory |
| `ferro bader` | Bader charge partitioning | ACF/BCF/AVF |
| `ferro convert` / `info` / `job` | Structure I/O, info, QC input files | files / stdout |

Run `ferro` with no arguments for the overview, `ferro <group>` for that group's
commands, and `ferro <group> <command>` without `-i` for its parameters.

## Python Bindings

`ferro-python` is its own workspace and is built with maturin (not the root
`cargo build`).  `maturin develop` needs an active virtualenv / conda env:

```bash
pip install maturin
cd ferro-python
maturin develop                                   # into the active env
# or, without an active env:
maturin build --release --interpreter "$(which python)"
pip install target/wheels/ferro-*.whl
```

See [Python Bindings](python.md) for the full API reference.

## Generating This Documentation

```bash
cd docs
mdbook build          # HTML output → docs/book/
mdbook-pandoc         # PDF output → docs/book/ferro.pdf
```

Requires [mdbook](https://rust-lang.github.io/mdBook/), [mdbook-pandoc](https://github.com/max-heller/mdbook-pandoc), pandoc, and a LaTeX distribution (XeLaTeX recommended).
