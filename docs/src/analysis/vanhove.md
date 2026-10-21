# Van Hove Self-Correlation Function

## Theory

The van Hove self-correlation function $G_s(r, \tau)$ is the probability distribution of single-particle displacements over a time interval $\tau$.  It provides a more complete picture of atomic motion than the MSD alone: it can reveal non-Gaussian dynamics, heterogeneous mobility, and jump diffusion.

### Definition

$$G_s(r, \tau) = \frac{1}{N} \sum_j \langle \delta\!\left(r - |\mathbf{r}_j(t_0 + \tau) - \mathbf{r}_j(t_0)|\right) \rangle_{t_0}$$

ferro writes the **radial probability density** $P(r, \tau) = 4\pi r^2 G_s(r, \tau)$, in Å⁻¹: the histogram count of bin $i$ divided by $N_\text{origins} N_\text{atoms} \Delta r$.  It does not depend on the bin width, and

$$\int_{r_\min}^{r_\max} P(r, \tau)\,dr = 1 - f_\text{out},$$

where $f_\text{out}$ is the fraction of displacements outside $[r_\min, r_\max)$.  Those displacements stay in the denominator — they are part of the distribution, just not drawn — so a cut tail shows up as a missing integral rather than being spread over the visible range.  $f_\text{out}$ is listed per input as `outside_fraction` in `[inputs]`, with a warning above 1 %: at long $\tau$ in a liquid, raise `--r-max`.

### Conventions in other codes

`p_r` is the same quantity as the `-or` output of GROMACS `gmx vanhove` (count divided by origins, atoms and bin width, in nm⁻¹; its first bin is centred on $r = 0$ and therefore half as wide) and as `gsrt` of pymatgen's `VanHoveAnalysis` (labelled $4\pi r^2 G_s$; each displacement is smeared with a Gaussian instead of binned).  Codes that divide the histogram only by origins and atoms output the probability mass per bin, $P(r_i)\,\Delta r$: it sums to 1, has no unit, and is smaller than ferro's `p_r` by the factor $\Delta r$ — 100 times at the default $\Delta r = 0.01$ Å.  Multiply `p_r` by `--dr` to compare with such output.

### Gaussian Reference

For a purely diffusive Gaussian process:

$$P^\text{Gaussian}(r, \tau) = 4\pi r^2 \left(\frac{1}{4\pi D\tau}\right)^{3/2} \exp\!\left(-\frac{r^2}{4D\tau}\right)$$

which is directly comparable with the `p_r` column (same units, Å⁻¹).

Deviations from this Gaussian form — e.g. a secondary peak at large $r$ — indicate heterogeneous dynamics or discrete jump events.

### Non-Gaussian Parameter

The non-Gaussian parameter $\alpha_2(\tau)$ quantifies the deviation from Gaussian behaviour:

$$\alpha_2(\tau) = \frac{3\langle r^4(\tau)\rangle}{5\langle r^2(\tau)\rangle^2} - 1$$

$\alpha_2 = 0$ for a Gaussian distribution; $\alpha_2 > 0$ indicates fat tails (fast-moving particles).  ferro does not compute $\alpha_2$; with $f_\text{out} \approx 0$ it can be estimated from the output as $\langle r^n \rangle \approx \sum_i r_i^n P(r_i)\,\Delta r$.

### Algorithm

1. Convert Cartesian positions to fractional coordinates.
2. Unwrap the fractional coordinates frame by frame: a jump of more than half a cell between neighbouring frames is taken as a boundary crossing and removed.  This is the **lattice view**, not the TOR scheme that `msd` uses — see [Periodic boundaries](#periodic-boundaries) below.
3. Convert back to absolute Cartesian positions with each frame's own cell.
4. For each time origin $p$ (every `--shift` frames, while $p + \tau$ is still inside the trajectory) and each selected atom $j$:
   $$r = |\mathbf{r}_j(p + \tau) - \mathbf{r}_j(p)|$$
   Accumulate into histogram bin $\lfloor (r - r_\min) / \Delta r \rfloor$ if $r_\min \le r < r_\max$.
5. Normalise: $P(r_i) = \text{count}(r_i) / (N_\text{origins} \cdot N_\text{atoms} \cdot \Delta r)$.

$\Delta r$ is `--dr` when $[r_\min, r_\max]$ holds a whole number of bins, otherwise slightly
narrower, $(r_\max - r_\min) / \lceil (r_\max - r_\min) / \texttt{dr} \rceil$, so the last bin ends
exactly at $r_\max$; the `dr` header line reports the width used.

## Parameters

```rust
pub struct VanHoveParams {
    pub tau: Option<usize>,            // lag [frames], 1..=n_frames-1; None = n_frames / 2
    pub shift: usize,                  // origin spacing [frames]; default: 1
    pub dt: f64,                       // time between stored frames [fs]; required (Default is NaN)
    pub r_min: f64,                    // default: 0.0 Å
    pub r_max: f64,                    // default: 10.0 Å
    pub dr: f64,                       // bin width [Å]; default: 0.01
    pub elements: Option<Vec<String>>, // None = all atoms
}
```

These are the library defaults.  On the command line `--dt` is required, and `r_min` is not exposed: `ferro traj vanhove` always uses $r_\min = 0$.  A `--tau` of the frame count or more is an error, as `--max-lag` is for the other time correlations.

| Flag | Default | Description |
|---|---|---|
| `--tau` | N/2 | lag [frames] |
| `--dt` | **required** | time between stored frames [fs] = MD time step × dump interval.  When the trajectory carries step numbers (LAMMPS dump, CP2K, OUTCAR), duplicated or unevenly spaced frames are refused |
| `--shift` | 1 | spacing between time origins [frames] |
| `--r-max` | 10.0 | maximum displacement [Å] |
| `--dr` | 0.01 | bin width [Å] |
| `--elements` | (all) | element filter |

## Output

`vanhove_<element>[_<suffix>].csv` (elements sorted and deduplicated; `vanhove_all…` without `--elements`): `file, r, p_r` — $P(r) = 4\pi r^2 G_s$ in Å⁻¹, $\int P\,dr = 1 - f_\text{out}$.

The lag time $\tau$ is recorded only in the `#` header block (once in frames, once in fs).  One $\tau$ per run;
when several $\tau$ are supported later a `tau` column will be added — that adds rows rather than changing the column structure.

One csv for all inputs, stacked with a `file` column; see [Output format conventions](../cli-reference.md#output-format-conventions).

## Usage

```bash
ferro traj vanhove -i traj.dump --tau 500 --dt 2.0 -s run1
```

```rust
use ferro_analysis::{VanHoveParams, calc_vanhove};

let params = VanHoveParams {
    tau: Some(500), shift: 1, dt: 2.0,
    r_min: 0.0, r_max: 8.0, dr: 0.02,
    elements: Some(vec!["Li".into()]),
};
let result = calc_vanhove(&traj, &params).unwrap();
println!("cut tail: {:.1} %", result.outside_fraction * 100.0);
let tables = result.to_tables();   // [("vanhove", Table{r, p_r})]
```

## Interpreting Results

| Feature | Interpretation |
|---|---|
| Single peak near $r = 0$ | Localised, non-diffusing atoms |
| Broad peak with Gaussian shape | Normal diffusion |
| Bimodal distribution | Coexistence of slow and fast populations |
| Secondary peak at $r \approx$ jump length | Discrete jump mechanism |

## Scope: self part only

ferro computes the **self part** $G_s$ (each atom against its own earlier position).  The **distinct
part** $G_d(r,\tau)$ — an atom against *other* atoms' earlier positions, which at $\tau = 0$ reduces to
$\rho\,g(r)$ and describes how the neighbour shell relaxes — is not computed, so the full
$G(r,\tau) = G_s + G_d$ is not available.  For $\tau = 0$ use `traj gr`.

## Periodic boundaries

Periodic trajectories are unwrapped in the lattice view (step 2 above): fractional coordinates are unwrapped, then multiplied by the current cell.  **Under NVT this is identical to the TOR scheme of `msd`.**  Under NPT it is not: the barostat's rescaling moves an atom that has crossed $n$ boundaries by $n$ times the change in box length, which widens $G_s$ at long $\tau$.  The mechanism and a worked example are in the [MSD page](msd.md#why-tor-and-not-the-lattice-view).  Switching vanhove to TOR is planned.  For non-periodic systems, raw Cartesian distances are used directly.
