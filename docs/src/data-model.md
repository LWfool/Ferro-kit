# Core Data Model

## Type Hierarchy

Every nexflux API accepts and returns `Trajectory`, even for single-frame files.

```
Trajectory
  └── Vec<Frame>
        ├── Vec<Atom>           — atomic species, positions, optional properties
        ├── Option<Cell>        — None = non-periodic
        ├── [bool; 3]           — pbc flags per axis
        └── Optional results    — energy, forces, stress, velocities
```

## Atom

```rust
pub struct Atom {
    pub element: String,
    pub position: Vector3<f64>,  // Å, Cartesian
    pub label: Option<String>,   // e.g. "Fe1", "Fe2"
    pub mass: Option<f64>,       // None → look up from element table
    pub magmom: Option<f64>,     // initial magnetic moment (DFT input)
    pub charge: Option<f64>,     // Bader / DDEC charge (post-processing)
}
```

The atom **index is implicit** (position in `Vec<Atom>`); no `index` field is stored to prevent inconsistency.

## Frame

```rust
pub struct Frame {
    pub atoms: Vec<Atom>,
    pub cell: Option<Cell>,              // None = non-periodic
    pub pbc: [bool; 3],
    pub charge: i32,
    pub multiplicity: u32,
    pub bonds: Option<Vec<(usize,usize)>>,
    pub energy: Option<f64>,             // eV
    pub forces: Option<Vec<Vector3<f64>>>, // eV/Å
    pub stress: Option<Matrix3<f64>>,    // eV/Å³
    pub velocities: Option<Vec<Vector3<f64>>>, // Å/fs (internal standard)
}
```

`pbc = [false,false,false]` → molecular system.  
`pbc = [true,true,false]` → surface / slab.  
`pbc = [true,true,true]` → bulk crystal or glass.

## Cell

```rust
pub struct Cell {
    pub matrix: Matrix3<f64>,  // row vectors a, b, c  [Å]
}
```

| Method | Description |
|---|---|
| `from_lengths_angles(a,b,c,α,β,γ)` | Construct from lattice parameters |
| `lengths() -> [f64; 3]` | \|a\|, \|b\|, \|c\| |
| `angles() -> [f64; 3]` | α, β, γ in degrees |
| `volume() -> f64` | Cell volume [Å³] |
| `fractional_to_cartesian(f)` | f·M |
| `cartesian_to_fractional(c)` | c·M⁻¹ |
| `wrap_position(c)` | Fold Cartesian position into [0, L) |
| `minimum_image(v)` | Apply minimum-image convention to vector v |

### Coordinate Conventions

The cell matrix stores row vectors:

$$\mathbf{M} = \begin{pmatrix} \mathbf{a} \\ \mathbf{b} \\ \mathbf{c} \end{pmatrix}$$

Fractional → Cartesian: $\mathbf{r} = \mathbf{f} \cdot \mathbf{M}$

Cartesian → Fractional: $\mathbf{f} = \mathbf{r} \cdot \mathbf{M}^{-1}$

For a triclinic cell:

$$\mathbf{M} = \begin{pmatrix}
a & 0 & 0 \\
b\cos\gamma & b\sin\gamma & 0 \\
c\cos\beta & c(\cos\alpha - \cos\beta\cos\gamma)/\sin\gamma & c\sqrt{1 - \cos^2\alpha - \cos^2\beta - \cos^2\gamma + 2\cos\alpha\cos\beta\cos\gamma}/\sin\gamma
\end{pmatrix}$$

## Trajectory

```rust
pub struct Trajectory {
    pub frames: Vec<Frame>,
    pub metadata: Option<TrajectoryMetadata>,
}
```

NPT trajectories (variable box per frame) are handled naturally: each `Frame` carries its own `Cell`.

## Pseudo-Element Labels

nexflux supports pseudo-element labels for sub-classified atoms (e.g. after network analysis):

| Label | Meaning |
|---|---|
| `P0/P1/P2/P3` | Phosphorus with Qn connectivity |
| `Of` | Free oxygen (0 P neighbours) |
| `On` | Non-bridging oxygen (1 P neighbour) |
| `Ob` | Bridging oxygen (≥2 P neighbours) |
| `Zn`, `Na`, … | Modifier cation symbols |

Atomic-number lookup (`elem_z`) resolves pseudo-labels by stripping suffixes: `"P3"` → P (Z=15), `"Ob"` → O (Z=8).
