# Rotational Autocorrelation Function $C_2(t)$

## Theory

The second-order rotational autocorrelation function $C_2(t)$ measures how quickly molecular orientations
decorrelate over time.  It is defined with the second-order Legendre polynomial $P_2(x) = (3x^2-1)/2$:

$$C_2(t) = \big\langle P_2\big(\hat{\mathbf{u}}(t_0)\cdot\hat{\mathbf{u}}(t_0+t)\big) \big\rangle$$

where $\hat{\mathbf{u}}$ is the unit orientation vector of a molecule.

| Property | Value |
|---|---|
| $C_2(0)$ | 1 |
| $C_2(\infty)$ | → 0 for isotropic tumbling |
| Range | $[-0.5,\ 1]$ |

For rotational diffusion $C_2(t) = e^{-t/\tau_c}$, where $\tau_c$ is the rotational correlation time.

## Orientation vector

For each centre atom $c$ the orientation vector is the sum of its bond vectors to every
neighbour-element atom within $r_\text{cut}$ (minimum image in periodic cells):

$$\mathbf{u}_c(t) = \sum_{n \in \text{neighbours}(c,\, r_\text{cut},\, t)} (\mathbf{r}_n - \mathbf{r}_c)_\text{min-image}, \qquad \hat{\mathbf{u}}_c = \mathbf{u}_c / |\mathbf{u}_c|$$

- The neighbour set is **found afresh in every frame**, so an exchange of neighbours (a proton hop, a
  bond breaking) changes the vector.
- A centre with **no neighbour** within $r_\text{cut}$ in a frame has no orientation there: that
  (molecule, frame) is **invalid**.
- Centres are chosen by element in the first frame.

**Difference from GROMACS.**  `gmx rotacf` defines the vector from fixed atoms — a pair $i\!-\!j$, or the
normal $\mathbf{ij}\times\mathbf{jk}$ of a triplet — so every molecule has a vector in every frame.
ferro's summed-bond vector suits a single bond (O–H in a hydroxide) or a lone pair direction (the H–O–H
bisector of water, which is the sum of the two O–H bonds).

**Symmetric centres cancel.**  For a regular tetrahedron (P with four O, Si with four O) the four bond
vectors sum to nearly zero; $\mathbf{u}$ is then set by the tetrahedron's *distortion*, not its
orientation, and $C_2$ measures how fast the distortion changes.  That is a legitimate quantity, but not
the tumbling of the unit.  To follow the tumbling itself, pick a single bond whose partner is unique
within $r_\text{cut}$.

## Averaging: valid pairs only

At lag $m$, a (molecule, origin) pair counts only when the vector is valid at **both** ends, and every
such pair has the same weight:

$$C_2(m) = \frac{\sum_{j}\sum_{t:\ \text{valid}(j,t)\,\wedge\,\text{valid}(j,t+m)} P_2\big(\hat{\mathbf{u}}_j(t)\cdot\hat{\mathbf{u}}_j(t+m)\big)}{\#\{(j,t) \text{ valid at both ends}\}}$$

MDAnalysis waterdynamics (`WaterOrientationalRelaxation`) keeps the same pairs — the molecules present at
both $t_0$ and $t_0+t$ — but averages over molecules first and over origins second.  ferro pools all
valid pairs instead, which is what lets the sum go through an FFT.  When every frame is valid the two
coincide.

A lag with **no** valid pair is left empty (NaN), not written as 0.  `valid_fraction` in `[inputs]` is
the share of (molecule, frame) with a vector; well below 1 means $r_\text{cut}$ does not hold the bond.

### How it is computed: FFT

$(\hat{\mathbf{u}}(t)\cdot\hat{\mathbf{u}}(t+m))^2 = \sum_{ik} q_{ik}(t)\,q_{ik}(t+m)$ with
$q_{ik} = \hat u_i \hat u_k$, a sum of products at the two times.  Setting $q = 0$ on invalid frames and
letting $\chi(t) \in \{0, 1\}$ mark valid frames,

$$C_2(m) = \frac{\tfrac32 \sum_{i \le k} w_{ik}\,\mathrm{AC}[q_{ik}](m) \;-\; \tfrac12\,\mathrm{AC}[\chi](m)}{\mathrm{AC}[\chi](m)}, \qquad w_{ii} = 1,\ w_{i<k} = 2$$

where $\mathrm{AC}[x](m) = \sum_j \sum_t x_j(t)\,x_j(t+m)$ is the all-origin autocorrelation sum
(per molecule by FFT, then summed).  $\mathrm{AC}[\chi](m)$ is exactly the number of valid pairs; it is
rounded back to an integer before dividing.

With every frame valid, $\mathrm{AC}[\chi](m) = N_\text{mol}(N-m)$ and this is term for term the
`gmx rotacf -P 2` sum in GROMACS `autocorr.cpp` (weight 1.5 on the three diagonal products, 3 on the
three off-diagonal ones, $-0.5\,(N-m)$).

## Time origins and the lag axis

As in [MSD](msd.md#time-origins-and-the-lag-axis): every lag uses every origin; the longest lag is
`--max-lag`, default $\lfloor N/2 \rfloor$ (the default of `gmx rotacf -acflen`), allowed $1\ldots N-1$;
`min_origins` $= N - m_\max$.

## Rotational correlation time

`integral` is the running integral $\int_0^t C_2\,dt'$ by the **trapezoidal rule**
($I(0)=0$, $I(m) = I(m-1) + \tfrac{dt}{2}(C_2(m-1)+C_2(m))$) — the rule behind GROMACS' "Correlation time
(integral over corrfn)".  It levels off at $\tau_c$ once $C_2$ has decayed.  If $C_2$ has not reached zero
by the longest lag, the integral underestimates $\tau_c$; fitting an exponential to $C_2$ is the
alternative (`gmx rotacf -fitfn exp`).

If any lag is empty (NaN), the integral is empty from there on.

## Parameters

```rust
pub struct RotCorrParams {
    pub center: String,           // centre element, e.g. "O"
    pub neighbor: String,         // neighbour element, e.g. "H"
    pub r_cut: f64,               // bond search cutoff [Å]; default: 1.2
    pub max_lag: Option<usize>,   // longest lag [frames], 1..=N-1; None = N/2
    pub dt: f64,                  // time between stored frames [fs]; default: 1.0
}
```

| CLI flag | Default | Meaning |
|---|---|---|
| `--center` | required | centre element |
| `--neighbor` | required | neighbour element |
| `--r-cut` | 1.2 | bond search cutoff [Å] |
| `--dt` | 1.0 | time between stored frames [fs] |
| `--max-lag` | $N/2$ | longest lag in frames |

## Output

`rotcorr_<center>-<neighbor>[_<suffix>].csv` (both are required, so `all` is never reached):
`file, time, c2, integral` (time in fs, integral in fs).

`[inputs]`: `frames`, `atoms` (the number of centres), `max_lag`, `min_origins`, `valid_fraction`.

All output is **one** csv; multiple inputs are stacked into a single table with a `file` column.  The `#`
comment block holds the shared parameters and the `[inputs]` list (`pandas.read_csv(comment="#")` drops
it).  `-o` takes the **output directory**; the batch suffix goes to `-s`.

## Usage

```bash
ferro traj rotcorr -i water.dump --center O --neighbor H --r-cut 1.2 --dt 2.0 -o run1
```

```rust
use ferro_analysis::md::{RotCorrParams, calc_rotcorr};

let params = RotCorrParams {
    center: "O".into(), neighbor: "H".into(), r_cut: 1.2, max_lag: Some(500), dt: 2.0,
};
let result = calc_rotcorr(&traj, &params)?;
```

## Verification

| Check | Where | Agreement |
|---|---|---|
| FFT vs direct average, every frame valid (the `gmx rotacf` case) | unit test `test_fft_matches_brute_force_all_valid` | $10^{-10}$ |
| same with one frame in five invalid per molecule | `test_fft_matches_brute_force_with_invalid_frames` | $10^{-10}$, NaN at the same lags |
| lags with no valid pair are NaN | `test_lag_without_valid_pair_is_nan` | — |
| trapezoidal integral | `test_integral_is_trapezoidal` | $10^{-10}$ |
| example glass trajectory, P–O at $r_\text{cut}$ = 1.8 Å (all valid) and 1.52 Å (0.5 % invalid), vs an independent numpy implementation | by hand, 2026-09-27 | $5\times10^{-8}$ (csv precision), NaN positions identical |

The numpy reference, given orientation vectors `U` of shape (frames, molecules, 3) with zero rows where
there is no neighbour:

```python
import numpy as np
n2 = (U ** 2).sum(axis=2); ok = n2 >= 1e-30
u = np.where(ok[..., None], U / np.sqrt(np.where(ok, n2, 1))[..., None], 0)
N = len(U)
c2 = []
for m in range(N // 2 + 1):
    both = ok[m:] & ok[:N - m]
    c = (u[m:] * u[:N - m]).sum(axis=2)
    c2.append((1.5 * c ** 2 - 0.5)[both].mean() if both.any() else np.nan)
```

## Troubleshooting

| Symptom | First thing to check |
|---|---|
| `no X atom has a Y neighbor within r_cut` | `--r-cut` shorter than the bond, or the element names |
| `valid_fraction` well below 1 | `--r-cut` at the edge of the bond-length distribution; raise it past the first g(r) peak |
| $C_2$ decays implausibly fast for a symmetric unit | summed bonds nearly cancel — see [Symmetric centres cancel](#orientation-vector) |
| empty `c2` at long lags | no molecule valid at both ends that far apart |
| values differ from ferro 0.3.2 or earlier | expected — see below |

## Differences from earlier versions

Up to 0.3.2 (following `code1/rotcorr.c`):

- **One origin.**  `--tau` defaulted to the whole trajectory, so only frame 0 was an origin; `--shift`
  had no effect.  Now every lag averages all origins; `--shift` is gone and `--tau` became `--max-lag`.
- **Invalid frames pulled $C_2$ down.**  Pairs without a vector were skipped, but the average was still
  divided by the full number of origins, and an origin with no valid pair counted as 0.  Now only valid
  pairs enter the numerator and the denominator.
- **Rectangular integral.**  The running integral was too large by $C_2(0)\,dt/2 = dt/2$.  Now
  trapezoidal.

## Implementation Notes

- Orientation vectors are computed for all frames first (per frame, a scan over all atoms for each
  centre), then correlated per molecule in parallel.
- FFT: `ferro-analysis/src/md/correlate.rs` (crate `rustfft`), seven autocorrelations per molecule.
