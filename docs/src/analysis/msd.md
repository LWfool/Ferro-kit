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
    pub fit_range: Option<(f64, f64)>,// linear-fit window as fractions; None = no fit
}
```

## Output

`msd_<element>[_<suffix>].csv` (elements sorted and deduplicated; `msd_all…` without `--elements`): `file, time, msd, msd_a, msd_b, msd_c`
(time in fs, squared displacement in Å²).

The `#` header holds only what the whole batch shares: `shift`, `dt`, the element selection and the fit
window as fractions.  Everything that differs per input is a column of the `[inputs]` list: `frames`,
`atoms`, `origins`, `species`, and — when `--fit-range FMIN,FMAX` is given — `t_lo`, `t_hi` (the window in fs),
`points`, `slope` (Å²/fs), `intercept` (Å²), `d_ang2_per_fs`, `d_err` and `r2`.  $D$ is also printed to stdout
in Å²/fs, cm²/s and m²/s.

All output is **one** csv; multiple inputs are stacked into a single table with a `file` column.  The `#` comment block holds the shared parameters and the `[inputs]` list
(`pandas.read_csv(comment="#")` drops it).  `-o` takes the **output directory**; the batch suffix goes to `-s`.


Figures are made from this csv outside ferro; see [Plotting](../plotting.md).

## Usage

```bash
ferro traj msd -i traj.dump --dt 2.0 --shift 10 -o run1
```

```rust
use ferro_analysis::md::{MsdParams, calc_msd};

let params = MsdParams { tau: Some(1000), shift: 10, dt: 2.0,
    elements: Some(vec!["Li".into()]), fit_range: Some((0.3, 0.8)) };
let result = calc_msd(&traj, &params).unwrap();
if let Some(fit) = &result.fit {
    println!("D = {:.3e} ± {:.1e} Å²/fs", fit.d_ang2_per_fs, fit.d_err);
}
```

## Extracting the Diffusion Coefficient

Fit the linear region (avoiding the ballistic regime at short $t$ and the noise-dominated long-$t$ tail):

$$D = \frac{1}{6} \cdot \frac{d\,\text{MSD}}{d\,t}$$

The fit is ordinary least squares $\text{MSD} = \text{slope}\cdot t + \text{intercept}$ over the points
$i \in [\operatorname{round}(f_\min (n-1)),\ \operatorname{round}(f_\max (n-1))]$, the same model as `gmx msd`.
To convert from Å²/fs to cm²/s multiply by $0.1$ (1 Å² = $10^{-16}$ cm², 1 fs = $10^{-15}$ s).

Check the window on a log-log plot (`scripts/plot_msd.py --loglog`): in the diffusive regime the MSD has
slope 1, and the window should sit inside that stretch — past the ballistic start and before the poorly
averaged tail.

### Uncertainty

`d_err` follows `gmx msd`: the window is split at its midpoint (shared by both halves), each half is fitted
separately, and $d_{err} = |D_1 - D_2|$.  It is empty when a half holds fewer than 2 points.

It is a **linearity check more than a confidence interval**.  A large `d_err` relative to $D$ means the MSD
bends inside the window — the window is not in the diffusive regime, or the tail is noise.  It is not a
statistical error bar: MSD points at different lags are strongly correlated, so least-squares errors
of any kind understate the true uncertainty.  For a rigorous one use generalized least squares or
Bayesian regression with the MSD covariance (the `kinisi` package), or the spread of $D$ over independent
runs.

## Implementation Notes

- Parallelism: per-origin `par_iter`; `unwrap_frac` runs serially before parallelisation (sequential dependency).
- Non-periodic path: Cartesian coordinates are used directly without unwrapping.
