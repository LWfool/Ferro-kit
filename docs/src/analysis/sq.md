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
    pub dq: f64,                 // default: 0.05 Å⁻¹
    pub weighting: SqWeighting,  // default: None
}

pub enum SqWeighting {
    None,     // equal weights (Faber-Ziman, no form factors)
    Xrd,      // X-ray atomic form factors
    Neutron,  // neutron coherent scattering lengths
    Both,     // compute both simultaneously
}
```

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
