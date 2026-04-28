# Introduction

**ferro** is a Rust library and command-line toolkit for post-processing molecular dynamics (MD) trajectories and computational chemistry data.  It is designed for periodic systems (crystals, surfaces, glasses) and covers:

- **Structural analysis** — g(r), S(q), bond angle distributions, Qn speciation
- **Dynamical analysis** — MSD, VACF, van Hove correlation, rotational correlation
- **Spatial maps** — 3-D density/velocity/force grids, jump-distance maps, hard-sphere occupancy, cluster SDF
- **I/O** — readers and writers for LAMMPS dump, VASP POSCAR/OUTCAR, XYZ, PDB, CIF, QE, CP2K, Gaussian cube

## Design Goals

| Goal | Implementation |
|---|---|
| Correct periodic-boundary handling | Minimum-image convention; fractional-coordinate unwrapping for NPT |
| Multi-component support | Partial g(r), directed CN, per-element filters throughout |
| Performance | Rayon data parallelism; linked-cell neighbour lists for O(N) searches |
| Extensibility | Strict layered crate architecture; no cross-dependencies between middle layers |

## Architecture

```
ferro-cli / ferro-python        ← only layer combining multiple crates
    ├── ferro-core                ← Atom, Frame, Trajectory, Cell; static data; units; errors
    ├── ferro-io        → core    ← format readers / writers
    ├── ferro-structure → core    ← supercell, vacuum, merge, box estimation
    ├── ferro-analysis  → core    ← md/, dft/ (future), ml/ (future)
    └── ferro-workflow  → core    ← QC software input builders
```

## Internal Units

| Quantity | Unit |
|---|---|
| Length | Å |
| Energy | eV |
| Force | eV/Å |
| Stress | eV/Å³ |
| Time | fs |
| Mass | amu |
| Charge | e |
| Temperature | K |
