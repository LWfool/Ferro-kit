# Structure Factor S(q)

## Theory

The static structure factor $S(q)$ is computed from the radial distribution function $g(r)$ via a Fourier sine transform (Faber-Ziman formalism).

### Partial S(q)

$$S_{\alpha\beta}(q) = 1 + \frac{4\pi\rho}{q} \int_0^\infty r\left[g_{\alpha\beta}(r) - 1\right] \sin(qr) \, dr$$

In discretised form:

$$S_{\alpha\beta}(q) = 1 + \frac{4\pi\rho}{q} \sum_{r_i} r_i \left[g_{\alpha\beta}(r_i) - 1\right] \sin(q r_i) \, \Delta r$$

where $\rho = N/V$ is the total number density.

### Weighted Total S(q)

For experimental comparison, the total S(q) can be weighted by scattering factors.

**Faber-Ziman weights:**

$$S_\text{weighted}(q) = \sum_{\alpha \leq \beta} w_{\alpha\beta}(q) \cdot S_{\alpha\beta}(q)$$

$$w_{\alpha\beta}(q) = \frac{(2 - \delta_{\alpha\beta}) \, c_\alpha c_\beta f_\alpha(q) f_\beta(q)}{\left[\sum_\gamma c_\gamma f_\gamma(q)\right]^2}$$

where $c_\alpha = N_\alpha / N$ is the mole fraction of species $\alpha$, and $f_\alpha(q)$ is the scattering factor.

- **X-ray** (`Xrd`): $f_\alpha(q)$ are the $q$-dependent atomic form factors.
- **Neutron** (`Neutron`): $f_\alpha$ is replaced by the $q$-independent coherent scattering length $b_\alpha^\text{coh}$.

### Physical Interpretation

$S(q)$ characterises structural correlations at length scale $2\pi/q$. Key features:

| Feature | Interpretation |
|---|---|
| First sharp diffraction peak (FSDP) at $q \approx 1$–2 Å⁻¹ | Medium-range order (network periodicity ~3–5 Å) |
| Principal peak at $q \approx 2$–3 Å⁻¹ | Nearest-neighbour distance |
| $S(q) \to 1$ as $q \to \infty$ | Loss of structural correlations at short wavelengths |

## Parameters

```rust
pub struct SqParams {
    pub q_min: f64,              // default: 0.1 Å⁻¹
    pub q_max: f64,              // default: 25.0 Å⁻¹
    pub dq: f64,                 // default: 0.02 Å⁻¹
    pub weighting: SqWeighting,  // default: None
}

pub enum SqWeighting {
    None,     // equal weights (Faber-Ziman, no form factors)
    Xrd,      // X-ray atomic form factors
    Neutron,  // neutron coherent scattering lengths
    Both,     // compute both simultaneously
}
```

| CLI flag | Field | Default |
|---|---|---|
| `--q-min` | `q_min` | 0.1 |
| `--q-max` | `q_max` | 25.0 |
| `--dq` | `dq` | 0.02 |
| `--weighting` | `weighting` | `both` |

The $g(r)$ that feeds the transform is controlled by the `gr` flags — `--r-min`, `--r-max`,
`--dr` — which therefore set the integration range and hence the truncation ripple.

## Comparison with `code2/dump2sq`

ferro follows `dump2sq`'s formula, with one deliberate difference and one caveat about its
input.

**The `+1`.** `dump2sq`'s `CalcSq` computes only
$\frac{4\pi\rho}{q}\sum_r r[g(r)-1]\sin(qr)\,\Delta r$ and omits the leading $1$ of the
definition, so its curves sit one unit low and do not approach 1 at large $q$. ferro writes
the standard form. Measured on a full trajectory, the residual is a pure constant:
1.00031 (XRD) / 1.00018 (neutron) for $q > 15$, with standard deviations 0.014 / 0.002.

**Atom ordering in the input dump.** `dump2sq.c:InitializeType` writes `data[i].type_new` on
frame 0 only; `ReadDataLammpstrj` never refreshes it, so atoms are classified by **array
index** rather than by element. If the dump was written without `dump_modify sort id`, its
partials — and both weighted totals — are meaningless from frame 1 onwards. Symptoms: all
partials collapse onto one another (mean pairwise correlation 0.959 on a 2004-atom test
case) and the XRD and neutron totals become nearly identical (correlation 0.9986), because
the Faber-Ziman weighting is washed out and both totals degenerate into the Fourier
transform of the *total* $g(r)$.

On an id-sorted dump the two programs agree closely. Measured on a 5003-atom,
50-frame NVT trajectory with identical parameters:

| | max&#124;Δ&#124; vs `dump2sq` | rms vs experiment |
|---|---|---|
| Neutron | 0.0027 | 0.0194 (both) |
| X-ray | 0.0008 | 0.0277 (both) |

`scripts/compare_sq.py` and `scripts/compare_sq_experiment.py` run this check automatically
and refuse to let an unsorted dump be read as evidence against ferro; the shared logic lives
in `scripts/trajcheck.py`.

## Output

```bash
fe-traj -m sq -i traj.dump -o output
# writes output.gr, output.sq
```

The `.sq` file header records both the g(r) parameters (used as input) and the S(q) parameters.  
Column ordering matches the `.gr` file. Additional columns `total_xrd` and/or `total_neutron` are appended when weighting is requested.

## Usage

```rust
use ferro_analysis::md::{GrParams, SqParams, SqWeighting, calc_gr, calc_sq_from_gr, write_sq};

let gr = calc_gr(&traj, &GrParams::with_auto_rmax(&traj)).unwrap();
let sq = calc_sq_from_gr(&gr, &SqParams {
    q_min: 0.5, q_max: 20.0, dq: 0.02,
    weighting: SqWeighting::Neutron,
});
write_sq(&gr, &sq, "output.sq").unwrap();
```

## Implementation Notes

- The Fourier transform is parallelised over $q$ values with `rayon::par_iter`.
- Scattering data (form factors and neutron lengths) are stored as static tables in `ferro_analysis::md::scattering_data`.
- For a single-element system, `total_xrd ≈ total` because all Faber-Ziman weights reduce to 1.
