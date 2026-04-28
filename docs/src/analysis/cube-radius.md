# Hard-Sphere Occupancy

## Theory

The hard-sphere occupancy map treats each atom as a sphere of radius $r$.  For each (frame, atom) pair, every voxel whose centre lies within $r$ of the atom is incremented by 1.  The resulting 3-D grid represents the cumulative time-averaged occupancy: high values indicate regions of space that are frequently occupied by the selected atom type.

Unlike the density mode (which bins atoms by position), this method smears each atom over a finite volume, producing smoother maps suitable for visualising preferred coordination sites or channel geometry in disordered materials.

### Voxel Inclusion Criterion

For an atom at fractional coordinate $\mathbf{f}_\text{atom}$ and a voxel with centre at fractional coordinate $\mathbf{v}$:

1. Compute the fractional displacement with minimum-image convention:
   $$\Delta\mathbf{f} = \mathbf{f}_\text{atom} - \mathbf{v}, \quad \Delta f_k \leftarrow \Delta f_k - \text{round}(\Delta f_k)$$

2. Convert to Cartesian using the cell matrix transpose:
   $$\Delta\mathbf{r} = \mathbf{M}^T \Delta\mathbf{f}$$

3. Increment the voxel count if $|\Delta\mathbf{r}| < r$.

The voxel centre fractional coordinates are:

$$v_k = \frac{i_k + 0.5}{n_k}, \quad i_k \in \{0, \ldots, n_k - 1\}$$

### Bounding-Box Optimisation

Brute-force comparison of every atom against every voxel is $O(N \cdot n_x n_y n_z)$.  Instead, for each atom, only the voxels within a bounding box of half-width

$$s_k = \left\lceil r \cdot \left\|(\mathbf{M}^T)^{-1}_k\right\| \cdot n_k \right\rceil + 1$$

around the atom's home voxel are checked.  The row norm of $(\mathbf{M}^T)^{-1}$ gives the maximum fractional offset in direction $k$ corresponding to a Cartesian distance of $r$, making the bound exact for triclinic cells.

## Parameters

```rust
pub struct CubeRadiusParams {
    pub nx: usize,                        // grid divisions along a axis; default: 100
    pub ny: usize,                        // grid divisions along b axis; default: 100
    pub nz: usize,                        // grid divisions along c axis; default: 100
    pub radius: f64,                      // hard-sphere radius [Å]; default: 0.7
    pub elements: Option<Vec<String>>,    // None = all atoms
}
```

## Output

A Gaussian cube file with raw occupancy counts per voxel.  Values are not normalised.  To convert to fractional occupancy, divide by $N_\text{frames} \cdot N_\text{atoms}$.

## Usage

```bash
# Li occupancy map, 100³ grid, radius 0.7 Å
fe-cube -m radius -i traj.dump --elements Li --radius 0.7 --nx 100 --ny 100 --nz 100 -o li_radius.cube

# All-atom occupancy
fe-cube -m radius -i traj.dump --radius 1.0 -o all_radius.cube
```

```rust
use ferro_analysis::md::{CubeRadiusParams, calc_cube_radius};
use ferro_io::write_cube;

let params = CubeRadiusParams {
    nx: 100, ny: 100, nz: 100,
    radius: 0.7,
    elements: Some(vec!["Li".into()]),
};
let result = calc_cube_radius(&traj, &params).unwrap();
write_cube("li_radius.cube", &result.cube).unwrap();
```

## Comparison with Density Mode

| Feature | Density mode | Radius mode |
|---|---|---|
| Resolution | Depends on bin width $1/n_k$ | Depends on $r$ and $n_k$ |
| Smoothing | None (point-like atoms) | Hard-sphere smearing |
| Normalisation | atoms/Å³ | Raw counts |
| Typical use | Structural analysis | Site visualisation |

## Implementation Notes

- Parallelism: per-frame `par_iter`.  Each frame independently produces a count array; results are reduced by element-wise addition.
- The bounding-box optimisation reduces the inner loop from $O(n_x n_y n_z)$ to $O(s_x s_y s_z)$ per atom, giving roughly $10\times$ speedup at default settings.
- Fractional coordinates are folded into $[0, 1)$ with `rem_euclid` before the bounding-box centre is computed, ensuring correct periodic wrapping.
