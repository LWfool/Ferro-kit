# Mean Squared Displacement (MSD)

## Theory

The mean squared displacement (MSD) quantifies how far atoms diffuse over time.  In the diffusive regime it grows linearly, and the self-diffusion coefficient $D$ is extracted via the Einstein relation.

### Einstein Relation

$$\text{MSD}(t) = \langle |\mathbf{r}(t_0 + t) - \mathbf{r}(t_0)|^2 \rangle$$

$$D = \lim_{t \to \infty} \frac{\text{MSD}(t)}{6t}$$

(factor of 6 for 3-D isotropic diffusion; use 2 for 1-D or 4 for 2-D).

### Time-Shift Averaging

To reduce statistical noise, the MSD is averaged over all possible time origins $t_0$ separated by `shift` frames:

$$\text{MSD}(\tau) = \frac{1}{N_\text{origins}} \sum_{p} \frac{1}{N_\text{atoms}} \sum_j \left|\mathbf{r}_j(p+\tau) - \mathbf{r}_j(p)\right|^2$$

where $p$ runs over $\{0, \text{shift}, 2\cdot\text{shift}, \ldots\}$ subject to $p + \tau \leq N_\text{frames}$.

### Periodic Boundary Handling

Atoms crossing periodic boundaries must be **unwrapped** to obtain true displacements.

**Algorithm** (matching code1/msd.c `EstimateMSD`):

1. Convert Cartesian coordinates to fractional: $\mathbf{f}_j^{(t)} = \mathbf{r}_j^{(t)} \cdot \mathbf{M}^{-1}$
2. Unwrap: for each step $t > 0$ and each fractional component $k$:
   $$f_{j,k}^{(t)} \leftarrow f_{j,k}^{(t)} - \text{round}\!\left(f_{j,k}^{(t)} - f_{j,k}^{(t-1)}\right)$$
3. Compute unwrapped Cartesian displacement $\Delta\mathbf{r}_j = \mathbf{\Delta f}_j \cdot \overline{\mathbf{M}}$

### NPT Trajectories

For variable-box (NPT) simulations, the box matrix changes each frame.  The total MSD uses the **average** of the origin- and endpoint-frame matrices:

$$\overline{\mathbf{M}} = \frac{\mathbf{M}^{(p)} + \mathbf{M}^{(p+\tau)}}{2}$$

The directional MSD along axis $k$ uses the endpoint cell parameter $L_k^{(p+\tau)}$ (simplified approximation matching code1):

$$\text{MSD}_k(\tau) = \langle (\Delta f_k)^2 \rangle \cdot \left(L_k^{(p+\tau)}\right)^2$$

### Directional MSD

For periodic systems, `msd_a`, `msd_b`, `msd_c` give displacements along the three crystal axes.  For non-periodic systems they correspond to Cartesian $x$, $y$, $z$.

## Parameters

```rust
pub struct MsdParams {
    pub tau: Option<usize>,           // window size [frames]; None = all frames
    pub shift: usize,                 // origin spacing [frames]; default: 1
    pub dt: f64,                      // time step [fs]; default: 1.0
    pub elements: Option<Vec<String>>,// None = all atoms
}
```

## Output

Columns: `time[fs]`, `msd[Å²]`, `msd_a[Å²]`, `msd_b[Å²]`, `msd_c[Å²]`

## Usage

```bash
nex-traj -m msd -i traj.dump --dt 2.0 --shift 10 -o output.msd
```

```rust
use nexflux_analysis::md::{MsdParams, calc_msd, write_msd};

let params = MsdParams { tau: Some(1000), shift: 10, dt: 2.0,
    elements: Some(vec!["Li".into()]) };
let result = calc_msd(&traj, &params).unwrap();
write_msd(&result, "output.msd").unwrap();
```

## Extracting the Diffusion Coefficient

Fit the linear region (avoiding the ballistic regime at short $t$ and the noise-dominated long-$t$ tail):

$$D = \frac{1}{6} \cdot \frac{d\,\text{MSD}}{d\,t}$$

To convert from Å²/fs to cm²/s: multiply by $10^{-16}$.

## Implementation Notes

- Parallelism: per-origin `par_iter`; `unwrap_frac` runs serially before parallelisation (sequential dependency).
- Non-periodic path: Cartesian coordinates are used directly without unwrapping.
