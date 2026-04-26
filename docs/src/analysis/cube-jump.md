# Jump Distance Distribution

## Theory

The jump distribution map identifies spatial hotspots of large atomic displacements (jump events).  For each atom and each time window $[t,\, t+\tau]$, the true Cartesian displacement is computed from unwrapped fractional coordinates.  If the displacement exceeds a threshold $d_\text{min}$, the event is recorded at a voxel position on the 3-D grid.

This analysis is particularly useful for characterising ion hopping in glassy or disordered materials: the resulting cube shows *where* jumps originate, terminate, or pass through, enabling the identification of preferred migration corridors.

### Fractional-Coordinate Unwrapping

To obtain the minimum-image displacement correctly for periodic (and NPT) trajectories, fractional coordinates are unwrapped in-place before computing differences.  For a single atom over all frames:

$$f_{j,k}^{(t)} \leftarrow f_{j,k}^{(t)} - \text{round}\!\left(f_{j,k}^{(t)} - f_{j,k}^{(t-1)}\right), \quad k \in \{a, b, c\}$$

This step-by-step correction removes periodic-boundary crossings from the fractional trajectory.

### Displacement Calculation

For each time window $[t,\, t+\tau]$, the unwrapped fractional displacement is:

$$\Delta\mathbf{f} = \mathbf{f}^{(t+\tau)} - \mathbf{f}^{(t)}$$

The true Cartesian displacement uses the average cell matrix of the two endpoint frames (consistent with the MSD NPT treatment):

$$\Delta\mathbf{r} = \overline{\mathbf{M}}^T\, \Delta\mathbf{f}, \quad \overline{\mathbf{M}} = \frac{\mathbf{M}^{(t)} + \mathbf{M}^{(t+\tau)}}{2}$$

A jump event is recorded when $|\Delta\mathbf{r}| > d_\text{min}$.

### Recording Position

Each jump event is recorded at the atom's wrapped fractional position.  Three choices are available:

| Option | Recorded position | Typical use |
|---|---|---|
| `Start` | $\mathbf{f}^{(t)}$ wrapped | Identify jump origins |
| `End` | $\mathbf{f}^{(t+\tau)}$ wrapped | Identify jump destinations |
| `Midpoint` | $(\mathbf{f}^{(t)} + \mathbf{f}^{(t+\tau)})/2$ wrapped | Highlight migration pathways |

The voxel index is computed by folding the fractional coordinate into $[0, 1)$ and multiplying by the grid size.

### Grid Value

Each voxel stores the raw count of jump events recorded at that location.  No normalisation is applied; to convert to a probability density, divide by the total event count.

## Parameters

```rust
pub struct CubeJumpParams {
    pub nx: usize,                        // grid divisions along a axis; default: 50
    pub ny: usize,                        // grid divisions along b axis; default: 50
    pub nz: usize,                        // grid divisions along c axis; default: 50
    pub tau: usize,                       // displacement lag [frames]; default: 1
    pub threshold: f64,                   // minimum displacement [Å]; default: 1.0
    pub elements: Option<Vec<String>>,    // None = all atoms
    pub record_at: JumpPosition,          // Start | End | Midpoint; default: Start
}
```

## Output

A Gaussian cube file with raw jump-event counts per voxel.  Typical post-processing:

1. Load in VESTA with isosurface rendering to visualise migration channels.
2. Normalise by total jump count to obtain the spatial probability distribution.

## Usage

```bash
# Jump origins for Li, lag = 1 frame, threshold = 1.0 Å
nex-cube -m jump -i traj.dump --elements Li --tau 1 --threshold 1.0 -o jump.cube

# Jump midpoints (migration pathway visualisation)
nex-cube -m jump -i traj.dump --elements Li --record-at midpoint -o jump_mid.cube
```

```rust
use nexflux_analysis::md::{CubeJumpParams, JumpPosition, calc_cube_jump};
use nexflux_io::write_cube;

let params = CubeJumpParams {
    nx: 50, ny: 50, nz: 50,
    tau: 1,
    threshold: 1.0,
    elements: Some(vec!["Li".into()]),
    record_at: JumpPosition::Start,
};
let result = calc_cube_jump(&traj, &params).unwrap();
write_cube("jump.cube", &result.cube).unwrap();
println!("{} jump events recorded", result.n_jumps);
```

## Implementation Notes

- Parallelism: per-atom `par_iter`.  Each atom independently unwraps its fractional trajectory and accumulates jump events into a local grid; results are reduced by element-wise addition.
- Fractional coordinates for all frames are pre-extracted into a `Vec<Vec<[f64;3]>>` before parallelisation (shared read-only, zero copying per thread).
- For NPT trajectories, cell matrices are also cached per frame so that the average $\overline{\mathbf{M}}$ can be computed without locking.
- Atoms absent from some frames (length-mismatch guard) are silently skipped.
