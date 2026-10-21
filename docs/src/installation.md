# Installation

## Requirements

- Rust 1.89+ (`rustup` recommended; nalgebra 0.35 needs it)
- Cargo (bundled with Rust)

## Build from Source

```bash
git clone https://github.com/ferro/ferro
cd ferro
cargo build --release
```

The build produces **one** binary, `target/release/ferro`.

Run `ferro` with no arguments for the overview of every command, `ferro <group>` for
that group's commands, and `ferro <group> <command>` without `-i` for its parameters.

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

See [Python Bindings](python/api.md) for the full API reference.

## Generating This Documentation

```bash
cd docs
mdbook build          # HTML output → docs/book/
mdbook-pandoc         # PDF output → docs/book/ferro.pdf
```

Requires [mdbook](https://rust-lang.github.io/mdBook/), [mdbook-pandoc](https://github.com/max-heller/mdbook-pandoc), pandoc, and a LaTeX distribution (XeLaTeX recommended).
