# Mean Squared Displacement (MSD)

## Theory

The mean squared displacement (MSD) quantifies how far atoms diffuse over time.  In the diffusive regime
it grows linearly, and the self-diffusion coefficient $D$ is extracted via the Einstein relation.

### Einstein Relation

$$\text{MSD}(t) = \langle |\mathbf{r}(t_0 + t) - \mathbf{r}(t_0)|^2 \rangle$$

$$D = \lim_{t \to \infty} \frac{\text{MSD}(t)}{6t}$$

(factor of 6 for 3-D isotropic diffusion; use 2 for 1-D or 4 for 2-D).

The average $\langle\cdot\rangle$ runs over atoms **and** over time origins $t_0$.  How ferro takes each
of those, and how it turns wrapped positions into the $\mathbf{r}(t)$ above, is spelled out below —
these are the places to look first when a number is in doubt.

## Time origins and the lag axis

### Every lag uses every origin

With $N$ frames and lag $m$ (in frames), there are $N - m$ pairs of frames $m$ apart, and ferro averages
over **all** of them:

$$\text{MSD}(m) = \frac{1}{N_\text{atoms}} \sum_j \frac{1}{N-m} \sum_{t=0}^{N-1-m} \left|\mathbf{r}_j(t+m) - \mathbf{r}_j(t)\right|^2$$

So short lags are averaged over almost $N$ origins and the longest lag over the fewest.  This is the
"windowed" MSD of MDAnalysis and `gmx msd`.  There is no origin stride: every frame is an origin.

### Longest lag: `--max-lag`, default $N/2$

The output lag axis is $m = 0, 1, \ldots, m_\max$ with $m_\max$ = `--max-lag`, or $\lfloor N/2 \rfloor$
when it is not given.  The time column is $t = m \cdot dt$.

| Quantity | Value | Where to read it |
|---|---|---|
| rows per input | $m_\max + 1$ | the csv |
| origins at lag $m$ | $N - m$ | — |
| origins at the longest lag | $N - m_\max$ | `min_origins` in `[inputs]` |
| allowed `--max-lag` | $1 \ldots N-1$ | a value outside fails that input, not the batch |

Why stop at $N/2$: past the middle of the trajectory the number of origins drops below half, and the
last lags rest on a handful of frame pairs — the curve turns noisy there and says little.  Stopping at
$N/2$ keeps at least $N/2$ origins at every lag.  `--max-lag N-1` gives the full axis when you want to
look at the tail anyway.

Because $m_\max$ follows each input's frame count, a batch of trajectories of different lengths gets
different lag axes; each input's $m_\max$ is the `max_lag` column of `[inputs]`.  `--last-n` cuts the
trajectory **before** any of this, so $N$ is the frame count after the cut.

**`--fit-range` is a fraction of this lag axis**, not of the trajectory.  With $N = 101$ frames the default
axis is $m = 0 \ldots 50$, and `--fit-range 0.3,0.8` fits $m = 15 \ldots 40$ (see
[the fit](#extracting-the-diffusion-coefficient) for the rounding).  The same fractions on
`--max-lag 100` would fit $m = 30 \ldots 80$.

### How it is computed: FFT

Summing all origins directly costs $O(N^2)$ per atom.  ferro uses the FFT algorithm of Calandrini et al.
(2011, the one behind MDAnalysis' `fft=True`), $O(N \log N)$, applied separately to each atom and each
Cartesian axis $x(t)$:

$$\text{MSD}_x(m) = S_1(m) - 2\,S_2(m)$$

$$S_2(m) = \frac{1}{N-m} \sum_{t=0}^{N-1-m} x(t)\,x(t+m)$$

$S_2$ is the autocorrelation, computed by zero-padding $x$ to a power of two $\geq 2N$ (so the circular
correlation of the FFT equals the linear one), transforming, taking $|X|^2$ and transforming back.
$S_1$ comes from a recursion over the squared positions $D(t) = x(t)^2$:

$$Q(0) = 2\sum_t D(t), \quad Q(m) = Q(m-1) - D(m-1) - D(N-m), \quad S_1(m) = \frac{Q(m)}{N-m}$$

Before this, each series has its mean subtracted.  The MSD does not change under a constant shift, but
$S_1$ and $2S_2$ are both of order $|x|^2$ and their difference is the small MSD; centring keeps them
near the size of the atom's excursion instead of its distance from the origin, so no significant digits
are lost to the subtraction.

The result is the exact all-origin average, not an approximation: the tests compare it against the direct
$O(N^2)$ sum to $10^{-9}$ Å² (see [Verification](#verification)).

## Periodic boundaries: TOR unwrapping

Dump files hold **wrapped** positions $\mathbf{w}(t)$, folded back into the box.  Displacements need
**unwrapped** positions $\mathbf{u}(t)$, which keep counting when an atom leaves through one face and
enters through the opposite one.  ferro builds them with the toroidal-view-preserving (TOR) scheme of
von Bülow, Bullerjahn and Hummer (2020), recommended by Bullerjahn et al. (2023) for diffusion
coefficients at constant pressure:

$$\mathbf{u}(0) = \mathbf{w}(0), \qquad \mathbf{u}(i+1) = \mathbf{u}(i) + \operatorname{mic}_{i+1}\!\big(\mathbf{w}(i+1) - \mathbf{w}(i)\big)$$

In words: each step adds the **shortest** displacement between the two wrapped positions, measured in the
box of the **later** frame $i+1$.  In one dimension this is eq 2 of Bullerjahn et al. (2023):

$$u_{i+1} = u_i + (w_{i+1} - w_i) - \left\lfloor \frac{w_{i+1} - w_i}{L_{i+1}} + \frac{1}{2} \right\rfloor L_{i+1}$$

For a general (triclinic) cell with lattice vectors as the rows of $\mathbf{M}_{i+1}$, the minimum image is
taken in fractional coordinates:

$$\mathbf{f} = \mathbf{M}_{i+1}^{-\mathsf{T}}\,\Delta\mathbf{w}, \qquad \mathbf{f} \leftarrow \mathbf{f} - \lfloor \mathbf{f} + \tfrac{1}{2} \rfloor, \qquad \operatorname{mic}_{i+1}(\Delta\mathbf{w}) = \mathbf{M}_{i+1}^{\mathsf{T}}\,\mathbf{f}$$

The floor-plus-half is written as in the paper; it differs from rounding only at exactly $\pm 1/2$.

### Why TOR and not the lattice view

The obvious alternative — unwrap fractional coordinates, then multiply by the current box — follows the
periodic *lattice*.  Under NPT the barostat rescales the whole lattice, so an atom that has crossed $n$
boundaries is moved by $n$ times the change in box length every time the box breathes, although it has
not moved relative to its neighbours.  The error grows with how far the atom has wandered, and inflates
$D$ in long NPT runs (von Bülow et al. 2020).

A concrete case, which is also a unit test: a cubic box of 10 Å, an atom at fractional $x = 0.95$ crosses
the face (next frame $x = 0.05$), then stays at $x = 0.05$ while the box goes to 11 Å and back to 10 Å.

| step | box | TOR $\Delta x$ | lattice-view $\Delta x$ |
|---|---|---|---|
| 0 → 1 | 10 → 10 | +1.00 | +1.00 |
| 1 → 2 | 10 → 11 | +0.05 | +1.05 |
| 2 → 3 | 11 → 10 | −0.05 | −1.05 |

TOR sees the atom sit still apart from the barostat's small rescaling of its wrapped position; the lattice
view sees it jump by a whole box-length change each way.

**Under NVT the two are identical**: with a constant box, the minimum image of each step is exactly the
fractional unwrap.  Only NPT (or any varying box) results depend on the choice.

### The coordinates are used as read

ferro neither shifts nor folds the coordinates it reads (the cube outputs fold their reference structure,
and that fold stays inside the cube).  For TOR this matters under NPT:

- **Moving the box origin to 0 frame by frame is wrong.**  LAMMPS rescales an NPT box about its centre, so
  `xlo` changes every frame, $x_{lo}(t) = c - L(t)/2$.  Subtracting it moves every atom by the same
  $\Delta L/2$ per step, a displacement that never happened, which adds up to $(L(t) - L(0))/2$.  It is
  bounded, so $D$ from a long fit survives, but the short-lag MSD is biased; on the 5-frame NPT test
  trajectory the lag-1 MSD rose by 57 %.  MDAnalysis's LAMMPS reader does subtract `xlo` per frame
  (since 2.4.0), and its `NoJump` transformation is the lattice-view scheme above, so its NPT MSD is not
  comparable with ferro's.
- **Folding during reading also perturbs it, slightly.**  Folding into $[0, L)$ instead of the box LAMMPS
  wrapped into adds boundary crossings, and under a varying box each one leaves a residual of order
  $\Delta L$; on the same trajectory the MSD moved by 0.03-0.2 %.  Under NVT neither operation changes the
  MSD.

### Conditions and limits

- **The dump interval must be short enough** that no atom moves half a box length (in any cell direction)
  between two stored frames.  Otherwise the shortest image is the wrong one and the MSD silently loses
  that step.  This holds for every unwrapping scheme, not only TOR; it is the first thing to check when
  an MSD looks too small.
- **Triclinic cells**: the fractional-coordinate minimum image is exact as long as the step is small
  compared with the cell, which the previous point already requires.
- **Atoms, not molecules.**  TOR is applied to each atom on its own.  Its unwrapped coordinates do not
  preserve bonds across the boundary (a molecule straddling a face can look stretched), which does not
  affect an *atomic* MSD.  A molecular (centre-of-mass) MSD would need whole molecules first; ferro does
  not compute one.
- **Already unwrapped input** (e.g. LAMMPS `xu yu zu`) is fine: each step is already the short one, so TOR
  reproduces it.
- **Every frame needs a cell** once the first frame has one; a frame without a cell, or with a singular
  cell, fails that input.
- **Non-periodic input** (no cell in the first frame) is used as it is, with no unwrapping.
- **The per-axis `pbc` flags are not read**: a frame with a cell is unwrapped in all three directions.
  For a slab this is harmless along the vacuum (atoms do not cross it), but report $D$ from the in-plane
  components (see below), not from the total.
- **No centre-of-mass drift removal.**  If the whole system drifts (a thermostat that does not conserve
  momentum, or a poorly zeroed initial velocity), the drift adds $v_\text{COM}^2 t^2$ to every atom's MSD
  and inflates $D$.  `gmx msd -rmcomm` and LAMMPS `compute msd com yes` subtract it; ferro does not.
  Check that the COM displacement is negligible against the MSD, or remove the drift before analysis.

## Components: Cartesian $x$, $y$, $z$

`msd_x`, `msd_y`, `msd_z` are the MSD of the Cartesian components of the displacement, each averaged the
same way as the total:

$$\text{MSD}_x(m) = \frac{1}{N_\text{atoms}} \sum_j \frac{1}{N-m} \sum_t \big(u_{j,x}(t+m) - u_{j,x}(t)\big)^2$$

and **$\text{msd} = \text{msd}_x + \text{msd}_y + \text{msd}_z$ exactly, for any cell shape**.

- For an orthorhombic box the Cartesian axes coincide with $a$, $b$, $c$, so the components are the
  crystal-axis components.
- For a triclinic box they are **not** along the crystal axes.  To study diffusion along a skewed lattice
  direction, project the displacement onto that direction yourself; a projection onto non-orthogonal axes
  would not add up to the total, which is why ferro does not report one.
- ferro fits only the total.  A 1-D coefficient from a component is $D_x = \text{slope}_x / 2$.

## Parameters

```rust
pub struct MsdParams {
    pub max_lag: Option<usize>,       // longest lag [frames], 1..=N-1; None = N/2
    pub dt: f64,                      // time step [fs]; default: 1.0
    pub elements: Option<Vec<String>>,// None = all atoms
    pub fit_range: Option<(f64, f64)>,// linear-fit window as fractions of the lag axis; None = no fit
}
```

| CLI flag | Default | Meaning |
|---|---|---|
| `--dt` | **required** | time between stored frames [fs] = MD time step × dump interval.  When the trajectory carries step numbers (LAMMPS dump, CP2K, OUTCAR), duplicated or unevenly spaced frames are refused |
| `--max-lag` | $N/2$ | longest lag in frames |
| `--elements` | all | which atoms enter the average, chosen by element in the first frame |
| `--fit-range` | off | `FMIN,FMAX` window as fractions of the lag axis → $D$ |
| `--last-n` | all | keep only the last N frames before anything else |

## Output

`msd_<element>[_<suffix>].csv` (elements sorted and deduplicated; `msd_all…` without `--elements`):
`file, time, msd, msd_x, msd_y, msd_z` (time in fs, squared displacement in Å²).

The `#` header holds only what the whole batch shares: `max lag` when given, the origin, unwrapping and
axis conventions, `dt`, the element selection, and the fit window as fractions.  Everything that differs
per input is a column of the `[inputs]` list:

| Column | Meaning |
|---|---|
| `frames` | $N$, after `--last-n` |
| `atoms` | atoms averaged |
| `max_lag` | $m_\max$ of this input |
| `min_origins` | $N - m_\max$, origins at the longest lag |
| `species` | elements present among those atoms |
| `t_lo`, `t_hi` | fit window in fs (with `--fit-range`) |
| `points` | lags in the fit window |
| `slope` | Å²/fs |
| `intercept` | Å² |
| `d_ang2_per_fs` | $D = \text{slope}/6$ |
| `d_err` | see [Uncertainty](#uncertainty) |
| `r2` | $R^2$ of the fit |

$D$ is also printed to stdout in Å²/fs, cm²/s and m²/s.

All output is **one** csv; multiple inputs are stacked into a single table with a `file` column.  The `#`
comment block holds the shared parameters and the `[inputs]` list (`pandas.read_csv(comment="#")` drops
it).  `-o` takes the **output directory**; the batch suffix goes to `-s`.

Figures are made from this csv outside ferro; see [Plotting](../plotting.md).  `scripts/plot_msd.py`
draws the fitted line from the `[inputs]` columns.

## Usage

```bash
ferro traj msd -i traj.dump --dt 2.0 --fit-range 0.3,0.8 -o run1
ferro traj msd -i 'runs/*.lammpstrj' --elements Li --dt 1000 --max-lag 500
```

```rust
use ferro_analysis::md::{MsdParams, calc_msd};

let params = MsdParams { max_lag: Some(1000), dt: 2.0,
    elements: Some(vec!["Li".into()]), fit_range: Some((0.3, 0.8)) };
let result = calc_msd(&traj, &params).unwrap();
if let Some(fit) = &result.fit {
    println!("D = {:.3e} ± {:.1e} Å²/fs", fit.d_ang2_per_fs, fit.d_err);
}
```

## Extracting the Diffusion Coefficient

Fit the linear region (avoiding the ballistic regime at short $t$ and the noise-dominated long-$t$ tail):

$$D = \frac{1}{6} \cdot \frac{d\,\text{MSD}}{d\,t}$$

The fit is ordinary least squares $\text{MSD} = \text{slope}\cdot t + \text{intercept}$ over the lags
$m \in [\operatorname{round}(f_\min\, m_\max),\ \operatorname{round}(f_\max\, m_\max)]$, the same model as
`gmx msd`.  To convert from Å²/fs to cm²/s multiply by $0.1$ (1 Å² = $10^{-16}$ cm², 1 fs = $10^{-15}$ s).

Check the window on a log-log plot (`scripts/plot_msd.py --loglog`): in the diffusive regime the MSD has
slope 1, and the window should sit inside that stretch — past the ballistic start and before the poorly
averaged tail.

### Isotropy and finite size

- **$D = \text{slope}/6$ assumes isotropic diffusion.**  For a slab, a channel structure or any anisotropic
  system use the components: $D_\alpha = \text{slope}_\alpha/2$ from `msd_x`, `msd_y`, `msd_z`
  (for in-plane diffusion in a slab, $D_\parallel = \text{slope}(\text{msd}_x + \text{msd}_y)/4$).
- **$D$ from a periodic box is systematically too small** because of the hydrodynamic self-interaction
  with the periodic images.  For a cubic box of side $L$ the Yeh–Hummer correction is
  $D_\infty = D_\text{PBC} + \xi k_B T / (6\pi\eta L)$ with $\xi \approx 2.837$ and $\eta$ the shear viscosity.
  It matters most for liquids in small boxes; ferro does not apply it, so compare $D$ between
  runs of the same box size, or correct it yourself.

### Uncertainty

`d_err` follows `gmx msd`: the window is split at its midpoint (shared by both halves), each half is fitted
separately, and $d_{err} = |D_1 - D_2|$.  It is empty when a half holds fewer than 2 points.

It is a **linearity check more than a confidence interval**.  A large `d_err` relative to $D$ means the MSD
bends inside the window — the window is not in the diffusive regime, or the tail is noise.  It is not a
statistical error bar: MSD points at different lags are strongly correlated, so least-squares errors
of any kind understate the true uncertainty.  For a rigorous one use generalized least squares or
Bayesian regression with the MSD covariance (the `kinisi` package), or the spread of $D$ over independent
runs.

## Verification

What has been checked, so a doubtful result can be narrowed down:

| Check | Where | Agreement |
|---|---|---|
| FFT sum vs direct $O(N^2)$ sum, NPT triclinic box with boundary crossings | unit test `test_fft_matches_brute_force_npt_triclinic` | < $10^{-9}$ Å² per lag and axis |
| same, non-periodic | `test_fft_matches_brute_force_nonperiodic` | < $10^{-9}$ Å² |
| straight-line motion across a skewed triclinic face | `test_tor_crosses_skewed_axis` | $(v\,m)^2$ to $10^{-9}$, $v$ the speed |
| TOR vs lattice view under a breathing box | `test_tor_does_not_follow_lattice_scaling` | the table above |
| components add up to the total | `test_components_sum_to_total` | $10^{-12}$ |
| two example LAMMPS trajectories (NPT, NVT) vs an independent numpy implementation | by hand, 2026-09-27 | within the csv's 7 significant digits |

The independent check reads the dump with ASE, unwraps with TOR, and averages every origin directly —
no FFT.  To repeat it on your own trajectory (`python` with numpy and ASE):

```python
import numpy as np
from ase.io import read

frames = read("traj.lammpstrj", index=":", format="lammps-dump-text")
w = np.array([f.get_positions() for f in frames])      # (N, atoms, 3), wrapped
cells = np.array([f.cell.array for f in frames])        # rows = lattice vectors
u = np.empty_like(w); u[0] = w[0]
for i in range(1, len(w)):                              # TOR, box of the later frame
    frac = (w[i] - w[i - 1]) @ np.linalg.inv(cells[i])
    frac -= np.floor(frac + 0.5)
    u[i] = u[i - 1] + frac @ cells[i]
N, max_lag = len(u), len(u) // 2
msd = [((u[m:] - u[:N - m]) ** 2).sum(axis=2).mean() for m in range(max_lag + 1)]
```

`msd` should match the `msd` column (with `--dt` only rescaling the time axis).  Restrict the atoms with
`frames[0].get_chemical_symbols()` to compare an `--elements` run.

## Troubleshooting

| Symptom | First thing to check |
|---|---|
| MSD far smaller than expected, or flat after a jump | dump interval too long — atoms move half a box between frames (see [Conditions](#conditions-and-limits)) |
| $D$ differs from ferro 0.3.2 or earlier | expected: earlier versions averaged a single time origin and used a lattice-view unwrap; see below |
| `--fit-range` now covers a different time span | the fractions are of the lag axis, which now ends at $N/2$ |
| `d_err` comparable to $D$ | window not in the linear regime, or glassy plateau with $D \approx 0$; look at `--loglog` |
| `msd_x`+`msd_y`+`msd_z` ≠ `msd` | should not happen; report it |
| an input fails with `max-lag must be within …` | `--max-lag` ≥ that input's frame count |

## Differences from earlier versions

Up to 0.3.2 ferro followed `code1/msd.c`, and results are not comparable:

- **Origins.**  The CLI had no `--tau`, so the window was the whole trajectory and only **one** origin
  (frame 0) fitted in it; every lag was a single displacement sample.  `--shift` had no effect.
- **Unwrapping.**  Fractional coordinates were unwrapped, and each displacement converted with the mean of
  the origin- and end-frame box matrices — a lattice-view scheme with the NPT problem described above.
- **Components.**  `msd_a/b/c` were fractional displacements times the end-frame axis length; they did not
  add up to the total for a triclinic box.  They are replaced by `msd_x/y/z`.

## Implementation Notes

- Parallelism: per atom (`par_iter`).  Each worker thread plans nothing — the FFT is planned once for the
  trajectory length and shared; each thread owns only its scratch buffers.
- Unwrapping is also per atom: a step depends on the previous step of the same atom only.
- Cost: $O(N_\text{atoms}\,N\log N)$ for the FFT plus $O(N_\text{atoms}\,N)$ for the unwrap.  2004 atoms ×
  101 frames runs in about 0.05 s.
- The FFT autocorrelation lives in `ferro-analysis/src/md/correlate.rs` (crate `rustfft`), shared with
  the other time-correlation functions as they move over.
