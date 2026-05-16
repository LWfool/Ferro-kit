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

The compiled CLI binaries are placed in `target/release/`:

| Binary | Purpose |
|---|---|
| `fe-convert` | Format conversion |
| `fe-info` | Print structure information |
| `fe-job` | Generate QC input files |
| `fe-traj` | Structural analysis (g(r), S(q), MSD, angle) |
| `fe-corr` | Correlation functions (VACF, rotcorr, van Hove) |
| `fe-cube` | 3-D spatial maps (density, jump, radius, SDF) |
| `fe-network` | Glass network analysis (Qn, CN, ligand classification) |

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
