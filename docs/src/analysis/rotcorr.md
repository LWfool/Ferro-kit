# Rotational Autocorrelation Function $C_\ell(t)$

## Theory

The rotational autocorrelation functions measure how quickly the orientation of a structural unit
decorrelates.  With $\hat{\mathbf{u}}$ a unit vector fixed to the unit and $P_\ell$ the Legendre
polynomials,

$$C_\ell(t) = \big\langle P_\ell\big(\hat{\mathbf{u}}(t_0)\cdot\hat{\mathbf{u}}(t_0+t)\big) \big\rangle, \qquad P_1(x) = x, \quad P_2(x) = \tfrac12(3x^2-1)$$

`--legendre` picks $\ell$ (default 2); the output column is `c1` or `c2` accordingly.

| Property | $C_1$ | $C_2$ |
|---|---|---|
| $C_\ell(0)$ | 1 | 1 |
| $C_\ell(\infty)$, isotropic | 0 | 0 |
| Range | $[-1,\ 1]$ | $[-0.5,\ 1]$ |
| rotational diffusion, coefficient $D_r$ | $e^{-2D_r t}$ | $e^{-6D_r t}$ |

**Why both.**  For small-step rotational diffusion $\tau_1/\tau_2 = 3$, from the $\ell(\ell+1)$ in the
exponents.  Large-angle jumps decorrelate both orders at once and push the ratio towards 1.  The ratio
of the two correlation times is therefore a standard diagnostic of the reorientation mechanism.  $C_2$ is
also the order probed by second-rank experiments (NMR relaxation, depolarised light scattering), $C_1$
the order of dielectric relaxation for a dipole along $\hat{\mathbf{u}}$.

`gmx rotacf -P 1` / `-P 2` compute the same two functions.

## Orientation vector: `--vector sum` or `bond`

| | `sum` (default) | `bond` |
|---|---|---|
| unit | one centre atom | one centre–neighbour bond |
| vector | sum of the centre's bonds within $r_\text{cut}$ | that bond |
| neighbours | searched afresh in every frame | fixed in the **first** frame |
| a bond stretching past $r_\text{cut}$ later | dropped (frame may become invalid) | still followed |
| use for | a single bond (O–H of a hydroxide), the H–O–H bisector of water | tetrahedra ($PO_4$, $SiO_4$, $AlO_4$), any unit whose bonds cancel |
| GROMACS | — | `gmx rotacf -d` with the pairs as index |

### `sum`

For each centre atom $c$ the orientation vector is the sum of its bond vectors to every
neighbour-element atom within $r_\text{cut}$ (minimum image in periodic cells):

$$\mathbf{u}_c(t) = \sum_{n \in \text{neighbours}(c,\, r_\text{cut},\, t)} (\mathbf{r}_n - \mathbf{r}_c)_\text{min-image}, \qquad \hat{\mathbf{u}}_c = \mathbf{u}_c / |\mathbf{u}_c|$$

- The neighbour set is **found afresh in every frame**, so an exchange of neighbours (a proton hop, a
  bond breaking) changes the vector.
- A centre with **no neighbour** within $r_\text{cut}$ in a frame has no orientation there: that
  (molecule, frame) is **invalid**.
- Centres are chosen by element in the first frame.

- **Cancellation is invalid.**  When the bonds cancel to within $10^{-6}$ of the summed bond lengths,
  $|\mathbf{u}| < 10^{-6}\sum|\mathbf{d}|$, what is left is floating-point residue whose direction is
  noise; that frame is treated as invalid too.  (An absolute threshold does not catch this: a perfect
  tetrahedron leaves a residue of $\sim 10^{-16}$ Å whose square still passes $10^{-30}$, and gave a
  meaningless $C_2 \approx 0.99$ in testing.)

**Symmetric centres.**  For a tetrahedron (P with four O, Si with four O) the four bonds nearly cancel;
in a real glass the residue is set by the tetrahedron's *distortion*, not its orientation, so `sum`
measures how fast the distortion changes — not the tumbling of the unit.  Use `bond` for that.  If every
frame cancels, `sum` stops with an error that says so.

### `bond`

Every pair (centre $c$, neighbour $n$) with $|\mathbf{r}_n - \mathbf{r}_c|_\text{min-image} < r_\text{cut}$ **in
the first frame** becomes one unit.  From then on the two atoms are followed by identity:

$$\mathbf{u}_{cn}(t) = (\mathbf{r}_n(t) - \mathbf{r}_c(t))_\text{min-image in frame } t$$

- This is `gmx rotacf -d`: GROMACS takes atom pairs from an index file and follows them whatever their
  distance; here the pairs are read off the first frame instead of an index file.  Every frame is valid,
  so `valid_fraction` is 1 and the averages are exactly the GROMACS ones.
- A $PO_4$ tetrahedron contributes its four P–O bonds as four units; the average over them is the
  tetrahedron's reorientation, and because the four bonds point in different directions, a rotation about
  any axis shows up.
- `units` in `[inputs]` is the number of bonds (1488 for 372 P with four O each); `atoms` is still the
  number of centres.
- **A bond that breaks is still followed** (the chosen GROMACS convention).  In a glass at room
  temperature bonds do not break on MD time scales, so this does not arise.  In a melt, after an exchange
  the vector points to the former partner; how often that happens is what a bond-breaking correlation
  function measures (planned, see `dev/plan.md`).
- **Minimum image per frame.**  If the two atoms ever drift more than half a box length apart (only
  possible after the bond has broken), the minimum image flips the vector to the other side.  Again a
  melt-only concern.
- The pair search is done once, so `bond` is also much faster than `sum`, which scans all atoms for every
  centre in every frame.

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

and for $\ell = 1$ simply

$$C_1(m) = \frac{\sum_i \mathrm{AC}[\hat u_i](m)}{\mathrm{AC}[\chi](m)}$$

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

`integral` is the running integral $\int_0^t C_\ell\,dt'$ by the **trapezoidal rule**
($I(0)=0$, $I(m) = I(m-1) + \tfrac{dt}{2}(C_2(m-1)+C_2(m))$) — the rule behind GROMACS' "Correlation time
(integral over corrfn)".  It levels off at the correlation time $\tau_\ell$ once $C_\ell$ has decayed.  If it has not reached zero
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
    pub vector: RotVector,        // Sum (default) | Bond
    pub legendre: Legendre,       // P1 | P2 (default)
}
```

| CLI flag | Default | Meaning |
|---|---|---|
| `--center` | required | centre element |
| `--neighbor` | required | neighbour element |
| `--r-cut` | 1.2 | bond search cutoff [Å]; in `bond` mode applied to the first frame only |
| `--vector` | `sum` | `sum` or `bond`, see [Orientation vector](#orientation-vector---vector-sum-or-bond) |
| `--legendre` | 2 | order $\ell$: 1 or 2 |
| `--dt` | 1.0 | time between stored frames [fs] |
| `--max-lag` | $N/2$ | longest lag in frames |

## Output

`rotcorr_<center>-<neighbor>[_<suffix>].csv` (both are required, so `all` is never reached):
`file, time, c2, integral` — `c1` in place of `c2` with `--legendre 1` (time in fs, integral in fs).

`[inputs]`: `frames`, `atoms` (centres), `units` (centres for `sum`, bonds for `bond`), `max_lag`,
`min_origins`, `valid_fraction`.

**The file name does not carry `--vector` or `--legendre`.**  Two runs on the same input into the same
directory overwrite each other; tell them apart with `-s` (e.g. `-s bondP1`).

All output is **one** csv; multiple inputs are stacked into a single table with a `file` column.  The `#`
comment block holds the shared parameters and the `[inputs]` list (`pandas.read_csv(comment="#")` drops
it).  `-o` takes the **output directory**; the batch suffix goes to `-s`.

## Usage

```bash
ferro traj rotcorr -i water.dump --center O --neighbor H --r-cut 1.2 --dt 2.0 -o run1

# PO4 tetrahedra in a glass: each P–O bond (first g(r) minimum as r-cut), both orders
ferro traj rotcorr -i glass.dump --center P --neighbor O --r-cut 1.8 --vector bond --legendre 1 -s P1
ferro traj rotcorr -i glass.dump --center P --neighbor O --r-cut 1.8 --vector bond --legendre 2 -s P2
```

```rust
use ferro_analysis::md::{Legendre, RotCorrParams, RotVector, calc_rotcorr};

let params = RotCorrParams {
    center: "P".into(), neighbor: "O".into(), r_cut: 1.8, max_lag: Some(500), dt: 2.0,
    vector: RotVector::Bond, legendre: Legendre::P2,
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
| rigid regular $PO_4$ rotating about $z$: `sum` reports cancellation, `bond` gives the analytic $C_2$ ($\cos = \tfrac13 + \tfrac23\cos\theta$ for every bond) | `test_sum_cancels_on_tetrahedron_bond_does_not` | $10^{-10}$ |
| `bond` + $P_1$ vs direct average | `test_p1_bond_matches_brute_force` | $10^{-10}$ |
| `bond` keeps following a bond stretched past $r_\text{cut}$ | `test_bond_is_followed_beyond_r_cut` | — |
| example glass, P–O `bond` mode, $P_1$ and $P_2$ (1488 bonds = 372 P × 4), vs independent numpy | by hand, 2026-09-27 | $5\times10^{-8}$ |
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
| `no X atom has an orientation in any frame` | `--r-cut` shorter than the bond, the element names, or (with `sum`) bonds that cancel — use `--vector bond` |
| `no X-Y bond within r_cut in the first frame` | `--r-cut` or element names; `bond` mode searches frame 0 only |
| `valid_fraction` well below 1 | `--r-cut` at the edge of the bond-length distribution; raise it past the first g(r) peak |
| $C_2$ decays implausibly fast for a symmetric unit | with `sum`, the bonds nearly cancel and it measures distortion — use `--vector bond` |
| a second run replaced the first | same file name for any `--vector` / `--legendre`; add `-s` |
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

- Orientation vectors are computed for all frames first — `sum`: per frame a scan over all atoms for each
  centre; `bond`: one pair search in frame 0, then one vector per pair and frame — then correlated per
  unit in parallel.
- FFT: `ferro-analysis/src/md/correlate.rs` (crate `rustfft`), seven autocorrelations per molecule.
