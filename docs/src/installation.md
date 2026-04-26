# Installation

## Requirements

- Rust 1.75+ (`rustup` recommended)
- Cargo (bundled with Rust)

## Build from Source

```bash
git clone https://github.com/nexflux/nexflux
cd nexflux
cargo build --release
```

The compiled CLI binaries are placed in `target/release/`:

| Binary | Purpose |
|---|---|
| `nex-convert` | Format conversion |
| `nex-info` | Print structure information |
| `nex-job` | Generate QC input files |
| `nex-traj` | Structural analysis (g(r), S(q), MSD, angle) |
| `nex-corr` | Correlation functions (VACF, rotcorr, van Hove) |
| `nex-cube` | 3-D spatial maps (density, jump, radius, SDF) |
| `nex-network` | Glass network analysis (Qn, CN, ligand classification) |

## Python Bindings

```bash
cd nexflux-python
pip install maturin
maturin develop
```

## Generating This Documentation

```bash
cd docs
mdbook build          # HTML output → docs/book/
mdbook-pandoc         # PDF output → docs/book/nexflux.pdf
```

Requires [mdbook](https://rust-lang.github.io/mdBook/), [mdbook-pandoc](https://github.com/max-heller/mdbook-pandoc), pandoc, and a LaTeX distribution (XeLaTeX recommended).
