# Rotational Autocorrelation Function C₂(t)

## Theory

The second-order rotational autocorrelation function $C_2(t)$ measures how quickly molecular orientations decorrelate over time.  It is defined using the second-order Legendre polynomial $P_2$:

$$C_2(t) = \langle P_2(\cos\theta(t)) \rangle = \left\langle \frac{3\cos^2\theta(t) - 1}{2} \right\rangle$$

where $\theta(t)$ is the angle between the orientation vector at time $t_0$ and time $t_0 + t$.

### Orientation Vector

For each centre atom $c$, the orientation vector $\mathbf{u}_c(t)$ is the vectorial sum of all bond vectors to neighbour atoms within the cutoff radius $r_\text{cut}$:

$$\mathbf{u}_c(t) = \sum_{n \in \text{neighbours}(c,\, r_\text{cut})} (\mathbf{r}_n - \mathbf{r}_c)_\text{min-image}$$

This generalises to arbitrary A-centre-B structures (e.g. water: O as centre, H as neighbour; phosphate tetrahedra: P as centre, O as neighbour).

### P₂ Calculation

For each centre atom $j$ and time origin $p$:

$$P_2\!\left[\cos\theta_j(p, \tau)\right] = \frac{3(\mathbf{u}_j(p) \cdot \mathbf{u}_j(p+\tau))^2}{2|\mathbf{u}_j(p)|^2 |\mathbf{u}_j(p+\tau)|^2} - \frac{1}{2}$$

Atoms with zero orientation vector (no neighbours) are excluded.

### Time-Shift Averaging

$$C_2(\tau) = \frac{1}{N_\text{origins}} \sum_p \frac{1}{N_\text{mol}} \sum_j P_2\!\left[\cos\theta_j(p, \tau)\right]$$

### Physical Properties

| Property | Value |
|---|---|
| $C_2(0)$ | Always = 1 (angle with itself is 0°) |
| $C_2(\infty)$ | → 0 for freely rotating molecules |
| Range | $[-0.5,\ 1]$ |

$C_2(t)$ for a perfect rotor decays exponentially: $C_2(t) = e^{-t/\tau_c}$, where $\tau_c$ is the rotational correlation time.

### Rotational Correlation Time

The rotational correlation time $\tau_c$ is obtained from the running integral:

$$\tau_c(t) = \int_0^t C_2(t') \, dt'$$

which converges to $\tau_c$ as $t \to \infty$.

## Parameters

```rust
pub struct RotCorrParams {
    pub center: String,      // centre element, e.g. "O"
    pub neighbor: String,    // neighbour element, e.g. "H"
    pub r_cut: f64,          // cutoff radius [Å]; default: 1.2
    pub tau: Option<usize>,  // window [frames]; None = all frames
    pub shift: usize,        // origin spacing; default: 1
    pub dt: f64,             // time step [fs]; default: 1.0
}
```

## Output

Columns: `time[fs]`, `C(t)`, `integral[fs]`

## Usage

```bash
fe-corr -m rotcorr -i traj.dump --center P --neighbor O --rcut 2.4 --dt 2.0 -o output.rotcorr
```

```rust
use ferro_analysis::md::{RotCorrParams, calc_rotcorr, write_rotcorr};

let params = RotCorrParams {
    center: "P".into(), neighbor: "O".into(),
    r_cut: 2.4, tau: Some(500), shift: 5, dt: 2.0,
};
let result = calc_rotcorr(&traj, &params).unwrap();
write_rotcorr(&result, "output.rotcorr").unwrap();
```

## Implementation Notes

- Orientation vectors are precomputed for all frames before the parallel loop.
- For periodic systems, the minimum-image convention is applied to bond vectors.
- Parallelism: per-origin `par_iter`.
