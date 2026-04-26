# Radial Distribution Function g(r)

## Theory

The radial distribution function (RDF) $g_{\alpha\beta}(r)$ describes the probability of finding a particle of type $\beta$ at distance $r$ from a particle of type $\alpha$, relative to an ideal gas at the same number density.

### Partial g(r)

For a system with $N_\alpha$ atoms of type $\alpha$ and $N_\beta$ atoms of type $\beta$ in a volume $V$:

$$g_{\alpha\beta}(r) = \frac{V}{N_\alpha N_\beta} \frac{\langle n_{\alpha\beta}(r, r+\Delta r)\rangle}{4\pi r^2 \Delta r}$$

where $n_{\alpha\beta}(r, r+\Delta r)$ is the number of $\beta$ atoms in the shell $[r, r+\Delta r)$ around any $\alpha$ atom.

In practice, the histogram count $C_{\alpha\beta}$ is accumulated and normalised as follows:

**Same species** ($\alpha = \beta$):
$$g_{\alpha\alpha}(r_i) = \frac{2 \, C_{\alpha\alpha}(r_i)}{4\pi r_i^2 \Delta r \cdot \frac{N_\alpha - 1}{V} \cdot N_\alpha \cdot N_\text{frames}}$$

**Different species** ($\alpha \neq \beta$):
$$g_{\alpha\beta}(r_i) = \frac{C_{\alpha\beta}(r_i)}{4\pi r_i^2 \Delta r \cdot \frac{N_\beta}{V} \cdot N_\alpha \cdot N_\text{frames}}$$

The bin centre is $r_i = r_\text{min} + (i + 0.5) \Delta r$.

### Total g(r)

All atoms are treated as a single species:

$$g_\text{total}(r_i) = \frac{2 \, C_\text{total}(r_i)}{4\pi r_i^2 \Delta r \cdot \frac{N-1}{V} \cdot N \cdot N_\text{frames}}$$

### Coordination Number CN(r)

The cumulative coordination number gives the average number of neighbours within distance $r$.

**Directed CN** (`"A-B"` = average B atoms within $r$ around each A):

$$\text{CN}_{A \to B}(r) = \sum_{r_i \leq r} \frac{m \cdot C_{AB}(r_i)}{N_A \cdot N_\text{frames}}$$

where $m = 2$ for same-species pairs and $m = 1$ for cross-species pairs.  
For $A \neq B$: the reverse direction `"B-A"` is also computed with $N_B$ as the denominator.

### Minimum-Image Convention

All pair distances use the minimum-image convention:

$$\mathbf{r}_{ij}^\text{min} = \mathbf{r}_{ij} - \mathbf{M} \cdot \text{round}\!\left(\mathbf{M}^{-1} \mathbf{r}_{ij}\right)$$

This is correct for orthorhombic and triclinic cells as long as $r_\text{max} < L_\text{min}/2$.

## Parameters

```rust
pub struct GrParams {
    pub r_min: f64,   // default: 0.005 Å
    pub r_max: f64,   // default: 10.005 Å  (use with_auto_rmax for safety)
    pub dr: f64,      // default: 0.01 Å
    pub r_cut: f64,   // default: 2.3 Å  (for short-range bond statistics)
}
```

`r_max` must satisfy $r_\text{max} < L_\text{min}/2$.  Use `GrParams::with_auto_rmax(&traj)` to set it automatically from the first frame.

## Output

| File | Content |
|---|---|
| `.gr` | Columns: `r[Å]`, then symmetric partial g(r) in periodic-table order, `total` last |
| `.cn` | Columns: `r[Å]`, then directed CN pairs in (Z_center, Z_neighbor) order |

Header lines record all parameters, atom counts, average volume, and number density.

### Bond Statistics

Within `r_cut`, the mean and standard deviation of pair distances are computed and printed in the `.cn` header:

$$\bar{d} = \frac{1}{N_\text{bonds}} \sum_k d_k, \quad \sigma = \sqrt{\frac{1}{N_\text{bonds}} \sum_k (d_k - \bar{d})^2}$$

## Usage

```bash
nex-traj -m gr -i traj.dump -o output
# writes output.gr, output.cn
```

```rust
use nexflux_analysis::md::{GrParams, calc_gr, write_gr, write_cn};

let params = GrParams::with_auto_rmax(&traj);
let result = calc_gr(&traj, &params).unwrap();
write_gr(&result, "output.gr").unwrap();
write_cn(&result, "output.cn").unwrap();
```

## Implementation Notes

- Parallelism: per-frame `par_iter` with `fold`/`reduce` accumulation.
- Pseudo-element labels (e.g. `"P0"`, `"Ob"`) are supported: `elem_z` resolves them by stripping numeric/alphabetic suffixes.
- Column ordering uses atomic number Z as the primary sort key, with string labels as tiebreaker for deterministic ordering of pseudo-elements.
