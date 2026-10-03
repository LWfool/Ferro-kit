# Velocity Autocorrelation Function (VACF)

## Theory

The velocity autocorrelation function (VACF) $C_v(t)$ measures how quickly atomic velocities decorrelate
over time.  Its time integral gives the self-diffusion coefficient (Green–Kubo), and its Fourier
transform gives the vibrational density of states (VDOS).

### Definition

$$C_v(t) = \langle \mathbf{v}(t_0) \cdot \mathbf{v}(t_0 + t) \rangle = \frac{1}{N_\text{atoms}} \sum_j \langle \mathbf{v}_j(t_0) \cdot \mathbf{v}_j(t_0 + t) \rangle_{t_0}$$

Expanding into Cartesian components:

$$C_v(t) = C_x(t) + C_y(t) + C_z(t), \qquad C_x(t) = \frac{1}{N_\text{atoms}} \sum_j \langle v_{x,j}(t_0)\, v_{x,j}(t_0+t) \rangle_{t_0}$$

At $t = 0$: $C_v(0) = \langle v^2 \rangle = 3 k_B T / m$ (equipartition), so $C_v(0)$ is a quick check on
the temperature and the velocity unit.

Velocities have no periodic jumps, so **no unwrapping** is involved.

## Time origins and the lag axis

This follows the same conventions as [MSD](msd.md#time-origins-and-the-lag-axis); the essentials:

- **Every lag uses every origin.**  With $N$ frames, lag $m$ averages the $N-m$ pairs of frames $m$ apart:

$$C_v(m) = \frac{1}{N_\text{atoms}} \sum_j \frac{1}{N-m} \sum_{t=0}^{N-1-m} \mathbf{v}_j(t)\cdot\mathbf{v}_j(t+m)$$

- **Longest lag** `--max-lag`, default $\lfloor N/2 \rfloor$ — the default of `gmx velacc -acflen`.  The
  output has $m_\max + 1$ rows; the longest lag rests on $N - m_\max$ origins (`min_origins` in
  `[inputs]`).  Allowed values are $1 \ldots N-1$; a value outside fails that input only.
- **FFT.**  The sum over origins is an autocorrelation, computed per atom and axis by zero-padding to a
  power of two $\geq 2N$, transforming, taking $|V|^2$ and transforming back — $O(N \log N)$, and exactly
  equal to the direct sum (tests compare the two to $10^{-12}$).  This is what `gmx velacc` and
  MDAnalysis transport-analysis (`fft=True`) do.

## Normalised VACF

`vacf_norm` $= C_v(t) / C_v(0)$ — the form `gmx velacc` writes by default (`-normalize`).  It starts at 1
and makes runs at different temperatures or of different elements comparable in shape.  When
$C_v(0) = 0$ (every selected atom at rest) it is left empty rather than written as 0.

## Self-diffusion: Green–Kubo

$$D = \frac{1}{3} \int_0^\infty C_v(t) \, dt$$

`diffusion` is the **running** integral $D(t) = \tfrac13 \int_0^t C_v\,d\tau$, evaluated with the
trapezoidal rule:

$$I(0) = 0, \qquad I(m) = I(m-1) + \frac{dt}{2}\big(C_v(m-1) + C_v(m)\big), \qquad D(m) = \frac{I(m)}{3}$$

This is the rule GROMACS uses for correlation integrals (`print_and_integrate`, "Use trapezoidal rule")
and MDAnalysis transport-analysis uses for `self_diffusivity_gk` (`scipy.integrate.trapezoid`).

**Reading $D$ off the column.**  $D(t)$ rises, may overshoot, and levels off once $C_v$ has decayed; $D$
is the plateau.  At long $t$ the integral accumulates noise from the poorly averaged tail and wanders —
take the plateau, not the last value.  `diffusion_end` in `[inputs]` is the last value only for a quick
comparison across inputs.

Units: $D$ in Å²/fs; multiply by 0.1 for cm²/s.

## Units

Velocities are used in the internal unit **Å/fs**:

- LAMMPS dumps: `--units real|metal` is required, because a dump does not record its unit system.
  `real` velocities are already Å/fs; `metal` (Å/ps) are converted.  Giving the wrong one makes $C_v$
  off by $10^6$ — check $C_v(0)$ against $3k_BT/m$.
- Extended XYZ: the `velocities` column is taken as it is, assumed Å/fs.

Then $C_v$ is in Å²/fs² and `diffusion` in Å²/fs.

## Parameters

```rust
pub struct VacfParams {
    pub max_lag: Option<usize>,        // longest lag [frames], 1..=N-1; None = N/2
    pub dt: f64,                       // time between stored frames [fs]; default: 1.0
    pub elements: Option<Vec<String>>, // None = all atoms
}
```

| CLI flag | Default | Meaning |
|---|---|---|
| `--dt` | **required** | time between stored frames [fs] = MD time step × dump interval.  When the trajectory carries step numbers (LAMMPS dump, CP2K, OUTCAR), duplicated or unevenly spaced frames are refused |
| `--max-lag` | $N/2$ | longest lag in frames |
| `--elements` | all | atoms averaged, chosen by element in the first frame |
| `--last-n` | all | keep only the last N frames before anything else |
| `--units <UNITS>` | — | `real` or `metal`: LAMMPS units of the dump velocities (required for a dump) |

## Output

`vacf_<element>[_<suffix>].csv` (elements sorted and deduplicated; `vacf_all…` without `--elements`):
`file, time, vacf, vacf_norm, vacf_x, vacf_y, vacf_z, diffusion` (time in fs).

The `#` header holds only what the whole batch shares (max lag when given, origin convention, `dt`,
element selection, units, integration rule).  Per input, `[inputs]` lists:

| Column | Meaning |
|---|---|
| `frames` | $N$, after `--last-n` |
| `atoms` | atoms averaged |
| `max_lag` | $m_\max$ of this input |
| `min_origins` | $N - m_\max$ |
| `diffusion_end` | $D$ at the longest lag (Å²/fs) — not the plateau |
| `species` | elements among those atoms |

All output is **one** csv; multiple inputs are stacked into a single table with a `file` column.  The `#`
comment block holds the shared parameters and the `[inputs]` list (`pandas.read_csv(comment="#")` drops
it).  `-o` takes the **output directory**; the batch suffix goes to `-s`.

## Usage

```bash
ferro traj vacf -i traj.dump --dt 2.0 --elements Li -o run1
ferro traj vacf -i traj.lammpstrj --dt 5.0 --units metal --max-lag 400
```

```rust
use ferro_analysis::md::{VacfParams, calc_vacf};

let params = VacfParams { max_lag: Some(500), dt: 2.0, elements: Some(vec!["Li".into()]) };
let result = calc_vacf(&traj, &params)?;
```

## Verification

| Check | Where | Agreement |
|---|---|---|
| FFT vs direct all-origin sum, irregular velocities, per axis | unit test `test_fft_matches_brute_force` | $10^{-12}$ |
| trapezoidal integral (constant $C_v$ gives $c\,t/3$, 0 at $t=0$) | `test_diffusion_is_trapezoidal` | $10^{-12}$ |
| $C_v(0)=0$ leaves `vacf_norm` empty | `test_zero_velocity_gives_nan_norm_not_zero` | — |
| synthetic AR(1) velocities, 30 atoms × 200 frames, vs an independent numpy implementation | by hand, 2026-09-27 | within the csv's 7 significant digits, all columns |

The numpy reference, for velocities `v` of shape (frames, atoms, 3) in Å/fs:

```python
import numpy as np
N, dt = len(v), 2.0
vacf = np.array([(v[m:] * v[:N - m]).sum(axis=2).mean() for m in range(N // 2 + 1)])
D = np.concatenate([[0], np.cumsum(0.5 * dt * (vacf[1:] + vacf[:-1]))]) / 3
```

## Troubleshooting

| Symptom | First thing to check |
|---|---|
| $C_v(0)$ far from $3k_BT/m$ | velocity unit: metal-unit dump read with `--units real` (factor $10^6$) |
| `frame k has no velocities` | the dump lacks `vx vy vz`, or only some frames carry them |
| $D$ differs from ferro 0.3.2 or earlier | expected — see below |
| `diffusion` never levels off | trajectory too short for $C_v$ to decay, or `--max-lag` too small |
| `diffusion` drifts at long $t$ | tail noise; read the plateau, shorten `--max-lag` |

## Differences from earlier versions

Up to 0.3.2 (following `code1/velcorr.c`):

- **One origin.**  `--tau` defaulted to the whole trajectory (the help page said half), so only frame 0
  was an origin; `--shift` had no effect.  Now every lag averages all origins; `--shift` is gone and
  `--tau` became `--max-lag`.
- **Rectangular integral.**  $D(t) = \tfrac13\sum_{i\le n} C_v(i)\,dt$ counted the first point in full,
  so $D$ was too large by $C_v(0)\,dt/6$ at every $t$.  Now trapezoidal.
- New column `vacf_norm`.

## Not yet: VDOS

$$g(\nu) \propto \int_0^\infty C_v(t) \cos(2\pi \nu t) \, dt$$

ferro does not compute the vibrational density of states yet (`dev/plan.md`).

## Implementation Notes

- Parallelism: per atom; the FFT is planned once per trajectory length and shared, each thread owns its
  buffers (`ferro-analysis/src/md/correlate.rs`, crate `rustfft`).
