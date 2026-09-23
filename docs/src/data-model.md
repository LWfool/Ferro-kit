# Core Data Model

## Type Hierarchy

Every ferro API accepts and returns `Trajectory`, even for single-frame files.

```
Trajectory
  └── Vec<Frame>
        ├── Vec<Atom>           — atomic species, positions, optional properties
        ├── Option<Cell>        — None = non-periodic
        ├── [bool; 3]           — pbc flags per axis
        └── Optional results    — energy, forces, stress, velocities,
                                  temperature, step
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
    pub temperature: Option<f64>,        // K, instantaneous ionic temperature
    pub step: Option<i64>,               // MD step number, when the engine prints one
}
```

`temperature` and `step` are filled by the AIMD readers and are `None` for
everything else (a POSCAR, a CIF, a relaxation). They live on the frame rather
than in arrays beside the trajectory because they are per-frame quantities:
`select` / `--stride` must carry them along, and a parallel array would shift out
of step with the frames without any error.

What fills them:

| reader | `temperature` | `step` |
|---|---|---|
| CP2K out | `MD\| Temperature [K]`, the instantaneous column | `MD\| Step number` |
| VASP OUTCAR | the value in `(temperature X K)` — despite the label it is the *ionic* temperature | ionic iteration |
| vasprun.xml | **derived**: `T = 2·E_kin/(3N·k_B)` from `<i name="kinetic">`, because the format does not print it | — |
| extxyz | `Temperature=` when present (read only; the writer never emits it) | — |

`pbc = [false,false,false]` → molecular system.  
`pbc = [true,true,false]` → surface / slab.

### Sign of the stress

`stress` is always **positive = compression** (the sign CP2K / VASP / QE print), stored row-major.
The virial follows from it directly and **without a sign flip**: `virial = stress × V` (eV),
which is exactly what DeePMD, QUIP and GPUMD call `virial`.

**ASE is the opposite** (positive = tension), so ferro flips the sign at the boundary when reading or writing ASE-family formats:

| Format and key | What is in the file | What ferro does |
|---|---|---|
| extxyz `stress=` | eV/Å³, positive = tension | **sign flipped** on both read and write |
| extxyz `virial=` | eV, positive = compression | divided by `\|det(box)\|`, **no sign flip**; errors out when there is no `Lattice` |
| DeePMD `virial.npy` | eV, positive = compression | `stress × V`, no sign flip |
| VASP `in kB` line | kBar, positive = compression, Voigt order `XX YY ZZ XY YZ ZX` | unit conversion only, **no sign flip** |

When an extxyz frame carries both `stress=` and `virial=`, the two must agree (`virial ≈ stress × V`),
otherwise ferro errors out — it will not pick one for you.  The tensor must also be **symmetric**
(the extxyz specification requires it), and the **6-component Voigt form is rejected**: its component
order is not consistent across programs (the specification and ASE use `xx yy zz yz xz xy`, VASP's
`in kB` and GPUMD's `stress_*.out` use `xx yy zz xy yz zx`), and nothing in the file says which one wrote it.  Rewrite it as 9 numbers.
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
| `cartesian_to_fractional(c)` | c·M^{-1} |
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

## Element vs. Label

`Atom::element` is the **chemical element**; `Atom::label` is an optional **site type**.  They are stored
separately, so an atom can be selected by element or by site without having to choose one or the other:

```rust
Atom { element: "P".into(), label: Some("P_3".into()), .. }
```

The two groups of selection flags on the analysis commands are mutually exclusive: `-a/-b/-c` select by `element`, `-x/-y/-z` by `label`.

### Label convention

A site label always has the form `<element>_<suffix>` and is split at the **first underscore**.  Further
underscores in the suffix do not matter (`Fe_site_A` still has element `Fe`) — this holds for labels from
any source, not only the ones `ferro net` writes.

| Label | Meaning |
|---|---|
| `P_0` … `P_4` | Qn network former; the number is the count of **same-element** connections n (P–O–P), see `network.md` |
| `Al_4` `Al_5` `Al_6` | non-Qn network former; the number is the **coordination number** |
| `O_f` | free ligand (bonded to no network former) |
| `O_n` | non-bridging ligand (bonded to one) |
| `O_b` | bridging ligand (bonded to two) |
| `O_t` | three-coordinate ligand (bonded to ≥3, tricluster) |
| `Zn`, `Na`… | modifier; bare element symbol, no role suffix |

> The pre-0.2.1 formats (`P0` / `Of` / `On_P` / `Ob_P_P` / `X` / `Zn_f`) have all been replaced,
> **with no read-compatibility layer**: rerunning `ferro net --export-traj` once produces the new
> format, whereas guessing the format could misread a real element as an old label.

### What fills `label`

| Source | Notes |
|---|---|
| LAMMPS dump reader | an `element` column written as `P_3` is split into `element="P"` + `label="P_3"`; the mapping is printed once after reading |
| extxyz reader | its own `label:S:1` column |
| CIF / CP2K inp / QE reader | their own site names (`O1`, `Fe1`) — note these do **not** follow the `<element>_<suffix>` convention |
| CP2K single-point out reader | the `&KIND` name (from the `ATOMIC KIND INFORMATION` block), filled **only when it differs from the element**.  Several kinds per element is common practice (`Fe1`/`Fe2` carrying two magnetic-moment guesses), and substituted systems have been seen where the kind name is simply another element symbol — so `element` always takes the real element from the coordinate table, and the kind name goes into `label` only |
| `ferro net --export-traj` | the classification goes into `label`, `element` stays the element |

### What writes `label` out

**No writer folds `label` into the element column on its own.**  `ferro convert` writes clean element
symbols whatever the labels are.  Only `ferro net --export-traj` folds, and only for LAMMPS dump — that
format has just one name column to put a second name in; extxyz uses its own `label:S:1` column, lossless in both directions.

The fold has a guard: a label not of the form `<element>_…` is not folded, and those occurrences are
counted and warned about.  A CIF `O1` folded in and read back would be a non-existent element `O1`.
