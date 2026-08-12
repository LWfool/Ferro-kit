# Van Hove Self-Correlation Function

## Theory

The van Hove self-correlation function $G_s(r, \tau)$ is the probability distribution of single-particle displacements over a time interval $\tau$.  It provides a more complete picture of atomic motion than the MSD alone: it can reveal non-Gaussian dynamics, heterogeneous mobility, and jump diffusion.

### Definition

$$G_s(r, \tau) = \frac{1}{N} \sum_j \langle \delta\!\left(r - |\mathbf{r}_j(t_0 + \tau) - \mathbf{r}_j(t_0)|\right) \rangle_{t_0}$$

In practice this is discretised as a histogram $P(r_i, \tau)$ normalised so that $\sum_i P(r_i, \tau) = 1$ (discrete probability mass function).

### Gaussian Reference

For a purely diffusive Gaussian process:

$$G_s^\text{Gaussian}(r, \tau) = \left(\frac{1}{4\pi D\tau}\right)^{3/2} \exp\!\left(-\frac{r^2}{4D\tau}\right) \cdot 4\pi r^2$$

Deviations from this Gaussian form — e.g. a secondary peak at large $r$ — indicate heterogeneous dynamics or discrete jump events.

### Non-Gaussian Parameter

The non-Gaussian parameter $\alpha_2(\tau)$ quantifies the deviation from Gaussian behaviour:

$$\alpha_2(\tau) = \frac{3\langle r^4(\tau)\rangle}{5\langle r^2(\tau)\rangle^2} - 1$$

$\alpha_2 = 0$ for a Gaussian distribution; $\alpha_2 > 0$ indicates fat tails (fast-moving particles).

### Algorithm

Follows code1/vanhove.c (`EstimateVanHove`):

1. Convert Cartesian positions to fractional coordinates.
2. Unwrap fractional coordinates (same procedure as MSD).
3. Convert back to absolute Cartesian positions.
4. For each time origin $p$ and each selected atom $j$:
   $$r = |\mathbf{r}_j(p + \tau) - \mathbf{r}_j(p)|$$
   Accumulate into histogram bin $\lfloor r / \Delta r \rfloor$.
5. Normalise: $P(r_i) = \text{count}(r_i) / (N_\text{origins} \cdot N_\text{atoms})$.

## Parameters

```rust
pub struct VanHoveParams {
    pub tau: Option<usize>,            // lag [frames]; None = n_frames - 1
    pub shift: usize,                  // origin spacing [frames]; default: 1
    pub dt: f64,                       // time step [fs]; default: 1.0
    pub r_min: f64,                    // default: 0.0 Å
    pub r_max: f64,                    // default: 10.0 Å
    pub dr: f64,                       // bin width [Å]; default: 0.01
    pub elements: Option<Vec<String>>, // None = all atoms
}
```

## Output

`vanhove[_<suffix>].csv`：`file, r, gs`（归一化为 $\sum g_s = 1$）。

滞后时间 $\tau$ 只在 `#` 头块里记录（帧数与 fs 各一份）。一次只算一个 $\tau$；
将来支持多 $\tau$ 时会加 `tau` 列——那是加行而不是改列结构。

所有产物是**一份** csv，多输入时堆叠成一张表并加 `file` 列；`#` 注释块里是共享参数与
`[inputs]` 清单（`pandas.read_csv(comment="#")` 会丢掉）。`-o` 给的是**文件名后缀**。


## Usage

```bash
ferro traj vanhove -i traj.dump --tau 500 --dt 2.0 -o run1
```

```rust
use ferro_analysis::md::{VanHoveParams, calc_vanhove, write_vanhove};

let params = VanHoveParams {
    tau: Some(500), shift: 1, dt: 2.0,
    r_min: 0.0, r_max: 8.0, dr: 0.02,
    elements: Some(vec!["Li".into()]),
};
let result = calc_vanhove(&traj, &params).unwrap();
write_vanhove(&result, "output.vanhove").unwrap();
```

## Interpreting Results

| Feature | Interpretation |
|---|---|
| Single peak near $r = 0$ | Localised, non-diffusing atoms |
| Broad peak with Gaussian shape | Normal diffusion |
| Bimodal distribution | Coexistence of slow and fast populations |
| Secondary peak at $r \approx$ jump length | Discrete jump mechanism |

## Implementation Notes

- Parallelism: per-origin `par_iter`.
- For periodic systems, fractional-coordinate unwrapping (shared with MSD) is applied before computing displacements.
- For non-periodic systems, raw Cartesian distances are used directly.
