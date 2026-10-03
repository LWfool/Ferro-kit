# Bond Lifetimes and Bond Events (`ferro traj bondlife`)

How long bonds between two elements live, and when they form and break.  Built for questions such as
"how long does an Si–O bond to a water oxygen last" or "when do Al–O bonds form at the glass–water
interface".  Bonds are defined by distance only (a contact criterion); there is no angle criterion,
so this is not a hydrogen-bond analysis.

## Which pairs are followed

A **candidate** is a (centre, neighbour) pair — centre of element `--center`, neighbour of element
`--neighbor` — that comes within `--r-bond` in **at least one frame**.  Distances use the minimum image
in each frame's own cell.  With `--center` equal to `--neighbor` (e.g. O–O), each pair is counted once.
The minimum image is exact only for distances up to half the smallest interplanar spacing of the cell; a `--r-break` (or `--r-bond` when `--r-break` is absent) beyond that bound in **any** frame (NPT boxes shrink) is rejected with the frame named, rather than silently missing neighbours in further images.

Pairs that never come within `--r-bond` never bond and do not matter to any quantity below; leaving them
out only saves work.

## When a pair counts as bonded: the existence function

Each candidate gets an existence function $h(t) \in \{0, 1\}$ per frame, with **hysteresis**:

$$h(t) = \begin{cases} [\,r(t) \le r_\text{bond}\,] & h(t-1) = 0 \text{ (free: must come close to bond)} \\ [\,r(t) \le r_\text{break}\,] & h(t-1) = 1 \text{ (bonded: stays bonded until it moves beyond } r_\text{break}) \end{cases}$$

and $h(0) = [\,r(0) \le r_\text{bond}\,]$.

- `--r-break` defaults to `--r-bond`, which makes it a single threshold — the contact criterion of
  `gmx hbond -contact`.
- With `--r-break` > `--r-bond`, a bond whose length vibrates around one threshold is not counted as
  breaking and re-forming on every vibration.  This is the two-threshold idea of Yamamoto and Onuki
  (1998): a pair is bonded at $r \le A_1\sigma$ and broken at $r > A_2\sigma$.
- **Choosing the thresholds.**  `--r-bond` at the first minimum of the centre–neighbour g(r); `--r-break`
  a little further, but below the second peak.  Yamamoto and Onuki found the result insensitive to both
  as long as they lie between the first and second peaks — check that yours is too by moving them.

## Two correlation functions

Both average over **all** pairs and **all** time origins, lag $m = 0 \ldots m_\max$
(`--max-lag`, default $N/2$; see [MSD](msd.md#time-origins-and-the-lag-axis) for the lag conventions).

### Intermittent: $C_I(t)$ — `c_int`

$$C_I(m) = \frac{\sum_\text{pairs} \sum_{t=0}^{N-1-m} h(t)\,h(t+m)}{\sum_\text{pairs} \sum_{t=0}^{N-1-m} h(t)}$$

The fraction of bonds present at $t$ that are present again at $t+m$ — **whether or not they broke in
between**.  This is the hydrogen-bond correlation of Luzar and Chandler (1996) and the `gmx hbond -ac`
function ("autocorrelations of the existence functions").  It decays when partners separate for good,
so it measures how long a pair stays a pair.  Computed by FFT autocorrelation of $h$.

### Continuous: $S_C(t)$ — `s_cont`

$$S_C(m) = \frac{\#\{(\text{pair}, t):\ h = 1 \text{ at every frame of } [t,\ t+m]\}}{\sum_\text{pairs} \sum_{t=0}^{N-1-m} h(t)}$$

The fraction of bonds present at $t$ that stay present **without a break** through $t+m$ — Rapaport's
(1983) continuous lifetime, and MDAnalysis' survival-probability `autocorrelation`.  It measures how long
one uninterrupted bonding episode lasts, and is always $\le C_I$.

Computed exactly from the lengths of the bonded stretches: a stretch of $L$ frames supplies $L - m$
origins to lag $m$ (none when $L \le m$).

**The continuous function depends on the frame interval.**  A break shorter than the dump interval is not
seen, one just longer ends the stretch.  That is what `--intermittency` is for.

### `--intermittency k`

Breaks of at most $k$ frames **between two bonded stretches** are filled before $S_C$ and the events are
computed — the `intermittency` parameter of MDAnalysis (Gowers and Carbone 2015), and Rapaport's
tolerance time.  A break that runs into the start or end of the trajectory is never filled: whether it
would close is unknown.  $C_I$ is computed on the unfilled $h$, since it tolerates breaks by definition.

Hysteresis handles flicker in distance, intermittency handles flicker in time; they can be combined.

### Which one to read

| | $C_I$ | $S_C$ |
|---|---|---|
| a bond breaks and re-forms with the same partner | still alive | ended |
| typical time scale | long: partner exchange | short: one episode |
| sensitive to dump interval | little | yes — use `--intermittency` |
| water H-bond literature reports | both, usually side by side | |

For bond formation and breaking between water and a glass surface, $C_I$ tells how long a given
Si–O(water) pair persists; $S_C$ how long each bonding episode lasts; the gap between them is the amount
of breaking and re-forming with the same partner.

## Lifetimes

Each input's `[inputs]` row carries two estimates for each function:

| Column | Definition |
|---|---|
| `tau_int_integral`, `tau_cont_integral` | trapezoidal integral of the function up to the longest lag — the integral lifetime $\langle\tau\rangle = \int_0^\infty C\,dt$ used for H-bond lifetimes |
| `tau_int_1e`, `tau_cont_1e` | first lag at which the function falls below $1/e$, linearly interpolated — the $\tau_b$ of Yamamoto and Onuki, $N_\text{bond}(\tau_b) = N_\text{bond}(0)/e$ |

Both are in fs (via `--dt`).  **Read them with care:**

- If the function has not decayed by the longest lag, the integral is only a lower bound, and the $1/e$
  lifetime is **empty** — not extrapolated.  In a glass at room temperature this is the normal case:
  bonds outlive the trajectory.
- Decay is often not exponential (Yamamoto and Onuki fit a stretched exponential
  $\exp[-(t/\tau_b)^{a'}]$ at low temperature); then the integral and the $1/e$ time differ, and neither
  is "the" lifetime.  Fit the `c_int` / `s_cont` columns yourself if you need a model.

## Bond events

`bondlife_events_<C>-<N>.csv` has one row per **frame** (not per lag):

| Column | Meaning |
|---|---|
| `time` | frame time [fs] |
| `n_bonds` | bonds present |
| `formed` | bonds present now but not in the previous frame |
| `broken` | bonds present in the previous frame but not now |

They use the same existence function, with gap filling, so $n(t) = n(t-1) + \text{formed}(t) -
\text{broken}(t)$ exactly.  `formed` and `broken` are **empty** at the first frame, which has no previous
frame.  A burst of `formed` marks when bonding happens; `mean_bonds` in `[inputs]` is the average of
`n_bonds`.

## Parameters

```rust
pub struct BondLifeParams {
    pub center: String,          // e.g. "Si"
    pub neighbor: String,        // e.g. "O"; may equal center
    pub r_bond: f64,             // a free pair bonds at r <= r_bond [Å]
    pub r_break: Option<f64>,    // a bond survives while r <= r_break [Å]; None = r_bond
    pub intermittency: usize,    // fill breaks of <= k frames (S_C and events); 0 = off
    pub max_lag: Option<usize>,  // longest lag [frames], 1..=N-1; None = N/2
    pub dt: f64,                 // time between stored frames [fs]
}
```

| CLI flag | Default | Meaning |
|---|---|---|
| `--center` | required | centre element |
| `--neighbor` | required | neighbour element (may equal `--center`) |
| `--r-bond` | required | bond-forming distance [Å] |
| `--r-break` | `--r-bond` | bond-breaking distance [Å], must be ≥ `--r-bond` |
| `--intermittency` | 0 | breaks of at most this many frames are filled |
| `--dt` | **required** | time between stored frames [fs] = MD time step × dump interval.  When the trajectory carries step numbers (LAMMPS dump, CP2K, OUTCAR), duplicated or unevenly spaced frames are refused |
| `--max-lag` | $N/2$ | longest lag [frames] |

## Output

- `bondlife_<C>-<N>[_<suffix>].csv`: `file, time, c_int, s_cont`
- `bondlife_events_<C>-<N>[_<suffix>].csv`: `file, time, n_bonds, formed, broken`
- `[inputs]`: `frames`, `atoms` (centres), `candidates`, `mean_bonds`, `max_lag`, `min_origins`, and the
  four lifetimes above.

Empty fields (NaN) mean "not defined", never zero: a lag with no bonded origin, a lifetime that did not
occur within the lag axis, the events of the first frame.

## Usage

```bash
# Si–O bonds in a glass–water system; first g(r) minimum at 2.2 Å
ferro traj bondlife -i glass_water.dump --center Si --neighbor O --r-bond 2.2 --dt 1000

# O–H with a hysteresis band and 2-frame tolerance
ferro traj bondlife -i run.dump --center O --neighbor H --r-bond 1.25 --r-break 1.4 --intermittency 2
```

## Verification

| Check | Where | Agreement |
|---|---|---|
| $C_I$ and $S_C$ vs counting every origin directly, irregular bonding | unit test `test_both_functions_match_brute_force` | $10^{-12}$ |
| break-and-reform counts for $C_I$, not for $S_C$ | `test_intermittent_counts_reformed_continuous_does_not` | exact |
| hysteresis stops flicker; a free pair needs $r \le r_\text{bond}$ | `test_hysteresis_suppresses_flicker`, `test_hysteresis_needs_r_bond_to_form` | exact |
| gap filling: $k$ filled, $k+1$ not, ends never | `test_fill_gaps_boundaries` | exact |
| $n(t) = n(t-1) + \text{formed} - \text{broken}$ | `test_events_conserve_bond_count` | exact |
| example glass, Zn–O (single threshold 2.4 Å; hysteresis 2.3/2.6 Å with intermittency 2) vs an independent numpy implementation | by hand, 2026-09-27 | $5\times10^{-8}$ (csv precision); events identical frame by frame |

The numpy reference, for distances `D` of shape (frames, candidate pairs):

```python
import numpy as np
N, P = D.shape
h = np.zeros_like(D, bool); on = np.zeros(P, bool)
for t in range(N):                                  # hysteresis
    on = np.where(on, D[t] <= r_break, D[t] <= r_bond); h[t] = on
c_int = [(h[:N-m] & h[m:]).sum() / h[:N-m].sum() for m in range(N // 2 + 1)]
hf = h.copy()                                       # fill gaps here for --intermittency > 0
s_cont = []
for m in range(N // 2 + 1):
    cont = np.ones((N - m, P), bool)
    for u in range(m + 1): cont &= hf[u:N - m + u]
    s_cont.append(cont.sum() / hf[:N - m].sum())
```

## Troubleshooting

| Symptom | First thing to check |
|---|---|
| `no X-Y pair within r_bond in any frame` | `--r-bond` below the bond length, or the element names |
| `c_int` ≈ 1 throughout, lifetimes empty | bonds outlive the trajectory (normal for a glass); the integral is a lower bound |
| `s_cont` falls much faster than `c_int` | frequent breaking and re-forming — or flicker at the threshold: raise `--r-break` or `--intermittency` and see whether it persists |
| results move with `--r-bond` | threshold not in the g(r) minimum between the first two peaks |
| very many `candidates` | `--r-bond` reaches into the second shell |
| slow on a large system | every centre–neighbour pair is scanned in every frame (no cell list yet) |

## References

- D. C. Rapaport, *Mol. Phys.* **50**, 1151 (1983) — continuous vs intermittent lifetimes.
- A. Luzar and D. Chandler, *Nature* **379**, 55 (1996) — $\langle h(0)h(t)\rangle/\langle h\rangle$.
- H. Yamamoto and A. Onuki, *Phys. Rev. E* **58**, 3515 (1998) — bond breakage with two thresholds,
  $\tau_b$ at $1/e$ (eqs 3.9–3.14).
- GROMACS `gmx hbond -ac -contact`; MDAnalysis `lib.correlations.autocorrelation` with `intermittency`.

## Implementation Notes

- Two passes: every frame is scanned for pairs within `--r-bond` (parallel over frames), then each
  candidate's distance series is built (parallel over candidates).
- $C_I$: FFT autocorrelation of $h$ (`correlate.rs`), rounded back to integer counts; $S_C$: stretch
  lengths, $O(\text{stretches} \times m_\max)$.
- Not implemented: angle (hydrogen-bond) criteria, cell lists, the Luzar–Chandler rate constants
  $k$, $k'$.
