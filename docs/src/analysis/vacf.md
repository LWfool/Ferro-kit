# Velocity Autocorrelation Function (VACF)

## Theory

The velocity autocorrelation function (VACF) $C_v(t)$ measures how quickly atomic velocities decorrelate over time.  Its Fourier transform gives the vibrational density of states (VDOS), which encodes information about lattice dynamics and diffusion mechanisms.

### Definition

$$C_v(t) = \langle \mathbf{v}(t_0) \cdot \mathbf{v}(t_0 + t) \rangle = \frac{1}{N} \sum_j \langle \mathbf{v}_j(t_0) \cdot \mathbf{v}_j(t_0 + t) \rangle_{t_0}$$

Expanding into Cartesian components:

$$C_v(t) = \frac{1}{N} \sum_j \left[ v_{x,j}(0) v_{x,j}(t) + v_{y,j}(0) v_{y,j}(t) + v_{z,j}(0) v_{z,j}(t) \right]$$

At $t = 0$: $C_v(0) = \langle v^2 \rangle = 3 k_B T / m$ (equipartition theorem).

### Self-Diffusion Coefficient

The Green-Kubo relation connects the VACF to the self-diffusion coefficient $D$:

$$D = \frac{1}{3} \int_0^\infty C_v(t) \, dt$$

In practice, the running integral is computed as:

$$D(t) = \frac{1}{3} \sum_{i=0}^{n} C_v(i \cdot \Delta t) \cdot \Delta t$$

which converges to $D$ as $t \to \infty$ (rectangular approximation, same as code1).

### Vibrational Density of States

The VDOS $g(\nu)$ is obtained via Fourier transform of the VACF:

$$g(\nu) \propto \int_0^\infty C_v(t) \cos(2\pi \nu t) \, dt$$

Peaks in $g(\nu)$ correspond to characteristic vibrational modes.

### Time-Shift Averaging

$$C_v(\tau) = \frac{1}{N_\text{origins}} \sum_p \frac{1}{N_\text{atoms}} \sum_j \mathbf{v}_j(p) \cdot \mathbf{v}_j(p + \tau)$$

Unlike position-based analyses, **no coordinate unwrapping is needed** — velocities do not have periodic boundary jumps.

## Parameters

```rust
pub struct VacfParams {
    pub tau: Option<usize>,            // window [frames]; None = all frames
    pub shift: usize,                  // origin spacing [frames]; default: 1
    pub dt: f64,                       // time step [fs]; default: 1.0
    pub elements: Option<Vec<String>>, // None = all atoms
}
```

## Output

Columns: `time[fs]`, `vacf[v²]`, `vacf_x[v²]`, `vacf_y[v²]`, `vacf_z[v²]`, `diffusion[v²·fs]`

### Unit Note

Velocities are stored in whatever unit the source trajectory uses.  LAMMPS metal-unit dump files store velocities in Å/ps.  To convert to internal standard units (Å/fs), multiply by $10^{-3}$.  As a consequence:
- $C_v$ in metal units: (Å/ps)² → multiply by $10^{-6}$ to get (Å/fs)²
- $D$ in metal units: Å²/ps → multiply by $10^{-3}$ to get Å²/fs → multiply by $10^{-16}$ to get cm²/s

This conversion will be applied automatically in a future IO unit normalisation update.

## Usage

```bash
nex-corr -m vacf -i traj.dump --dt 2.0 --elements Li -o output.vacf
```

```rust
use nexflux_analysis::md::{VacfParams, calc_vacf, write_vacf};

let params = VacfParams { tau: Some(500), dt: 2.0, shift: 1,
    elements: Some(vec!["Li".into()]) };
let result = calc_vacf(&traj, &params).unwrap();
write_vacf(&result, "output.vacf").unwrap();
```

## Implementation Notes

- Requires `frame.velocities` to be populated (returns `None` if any frame lacks velocities).
- Parallelism: per-origin `par_iter`.
