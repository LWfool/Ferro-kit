# 3-D Spatial Density Maps

## Theory

The spatial density map divides the simulation box into an $n_x \times n_y \times n_z$ voxel grid and accumulates a time-averaged scalar quantity for each voxel.  Three modes are supported:

| Mode | Accumulated quantity | Unit |
|---|---|---|
| `Density` | atomic number density | atoms/Å³ |
| `Velocity` | mean atomic speed $|\mathbf{v}|$ | Å/fs |
| `Force` | mean atomic force magnitude $|\mathbf{f}|$ | eV/Å |

The result is a Gaussian cube file that can be visualised directly in VESTA, VMD, or similar tools.

### Density Mode

Each atom is mapped to a voxel by converting its Cartesian position to fractional coordinates:

$$\mathbf{f}_j = \mathbf{M}^{-1} \mathbf{r}_j, \quad f_{j,k} = f_{j,k} \bmod 1$$

The voxel index along axis $k$ is:

$$i_k = \lfloor f_{j,k} \cdot n_k \rfloor \bmod n_k$$

After accumulating counts over all atoms and frames, the number density in voxel $(i, j, k)$ is:

$$\rho_{ijk} = \frac{\text{count}_{ijk}}{N_\text{frames} \cdot V_\text{voxel}}$$

where the voxel volume is:

$$V_\text{voxel} = \frac{V_\text{cell}}{n_x \cdot n_y \cdot n_z}$$

### Velocity / Force Mode

For `Velocity` and `Force` modes the average magnitude is computed per voxel:

$$\langle |\mathbf{q}| \rangle_{ijk} = \frac{\sum_{\text{visits}} |\mathbf{q}|}{\text{count}_{ijk}}$$

where $\mathbf{q}$ is the velocity or force vector of the visiting atom.  Voxels with no visits are set to zero.

### Voxel Spacing Matrix

The output cube file encodes the voxel step vectors.  For a general triclinic cell with matrix $\mathbf{M}$ (rows = $\mathbf{a}, \mathbf{b}, \mathbf{c}$), the spacing matrix is:

$$\mathbf{S} = \begin{pmatrix} \mathbf{a}/n_x \\ \mathbf{b}/n_y \\ \mathbf{c}/n_z \end{pmatrix}$$

This ensures correct spatial mapping for both orthogonal and triclinic simulation cells.

### Time-Averaged Structure

The cube file header includes a time-averaged atomic structure (mean position over all frames), which serves as a reference geometry for visualisation.

The grid always starts at the origin and spans the cell, $[0, L)$ along each axis, and every atom is binned by
its position folded back into the cell.  The reference structure is folded the same way: an atom whose mean
position lies outside the cell along a periodic axis is moved back by whole lattice vectors, so that the atoms
and the grid line up in a viewer.  A LAMMPS box with `xlo` $\neq 0$ is the usual case: its atoms sit in
$[x_{lo}, x_{lo} + L)$, and without the fold the part beyond $L$ would be drawn outside the grid.

**This fold is specific to the cube outputs** (`map density`, `map velocity`, `map force`, `map radius`).
Everywhere else ferro keeps the coordinates exactly as read, neither shifted nor folded, because folding
during reading would perturb time-correlation analyses of NPT trajectories (see the MSD page).

## Parameters

```rust
pub struct CubeDensityParams {
    pub nx: usize,                        // grid divisions along a axis; default: 50
    pub ny: usize,                        // grid divisions along b axis; default: 50
    pub nz: usize,                        // grid divisions along c axis; default: 50
    pub elements: Option<Vec<String>>,    // None = all atoms
    pub mode: CubeMode,                   // Density | Velocity | Force; default: Density
}
```

## Output

A Gaussian cube file (`.cube`):
- Header: time-averaged atomic positions + unit cell
- Data: $n_x \times n_y \times n_z$ scalar field

The cube format is directly accepted by VESTA, VMD, Ovito, and most electronic-structure visualisation packages.

### Units in the file — read before integrating

- **Geometry is in Bohr** (the cube convention): origin, voxel vectors and atom positions.
- **Data is in the table's unit**, atoms/Å³ for `density`, Å/fs or eV/Å for `velocity` / `force` — not
  converted to Bohr.  Viewers only colour the values, so isosurfaces are fine; but a tool that
  *integrates* the grid (VESTA's or Bader-style integration, which assume e/Bohr³) multiplies by a voxel
  volume in Bohr³, and the total comes out too large by $(1\ \text{Å}/1\ \text{Bohr})^{3} \approx 6.748$.
  To integrate, multiply the values by the voxel volume in Å³ yourself (the sum over the grid then gives
  the mean number of atoms in the box).
- In `velocity` and `force` modes a voxel **no atom ever visited is written as 0**: the cube format has
  no missing-value marker.  0 there means "no data", not "measured at rest" — do not average those
  voxels in.

## Usage

```bash
# Number density of all atoms
ferro map density -i traj.dump --nx 80 --ny 80 --nz 80 -o density.cube

# Li-only density
ferro map density -i traj.dump --elements Li -o li_density.cube

# Time-averaged velocity magnitude
ferro map velocity -i traj.dump --units metal -o velocity.cube

# Time-averaged force magnitude
ferro map force -i traj.dump -o force.cube
```

```rust
use ferro_analysis::md::{CubeDensityParams, CubeMode, calc_cube_density};
use ferro_io::write_cube;

let params = CubeDensityParams {
    nx: 80, ny: 80, nz: 80,
    elements: Some(vec!["Li".into()]),
    mode: CubeMode::Density,
};
let result = calc_cube_density(&traj, &params).unwrap();
write_cube("li_density.cube", &result.cube).unwrap();
```

## Implementation Notes

- Parallelism: per-frame `par_iter`.  Each frame independently produces `(count, value_sum)` arrays; results are reduced by element-wise addition.
- Fractional coordinates are folded with `rem_euclid(1.0)` to handle atoms outside the nominal box (e.g. from NPT fluctuations or wrapped dump files).
- `Velocity` mode requires `frame.velocities` to be populated; `Force` mode requires `frame.forces`.  Frames missing the required data are silently skipped.
- For NPT trajectories the spacing matrix is taken from the first periodic frame; the density is therefore referenced to that cell geometry.
