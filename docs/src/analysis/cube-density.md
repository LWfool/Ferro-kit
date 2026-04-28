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

## Usage

```bash
# Number density of all atoms
fe-cube -m density -i traj.dump --nx 80 --ny 80 --nz 80 -o density.cube

# Li-only density
fe-cube -m density -i traj.dump --elements Li -o li_density.cube

# Time-averaged velocity magnitude
fe-cube -m velocity -i traj.dump --metal-units -o velocity.cube

# Time-averaged force magnitude
fe-cube -m force -i traj.dump -o force.cube
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
