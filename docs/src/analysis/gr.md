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
    pub r_min: f64,       // default: 0.001 Å
    pub r_max: f64,       // default: 10.005 Å  (use with_auto_rmax for safety)
    pub dr: f64,          // default: 0.002 Å
    pub group_by: GroupBy, // resolve partials over elements or site labels
}
```

| CLI flag | Field | Default |
|---|---|---|
| `--r-min` | `r_min` | 0.001 |
| `--r-max` | `r_max` | 10.005 |
| `--dr` | `dr` | 0.002 |
| `-a`/`-b` (element) or `-x`/`-y` (label) | `group_by` | all pairs |

`r_max` must satisfy $r_\text{max} < L_\text{min}/2$; it is clamped internally to half the
smallest interplanar spacing seen across all frames.  Use `GrParams::with_auto_rmax(&traj)`
to set it from the first frame.

Bin $i$ covers $[r_\text{min} + i\,\Delta r,\ r_\text{min} + (i{+}1)\Delta r)$ and is labelled
at its centre, $r_i = r_\text{min} + (i + \tfrac{1}{2})\Delta r$.

To land on the same grid as `code1/gr.c` and `code2/dump2sq.c`, pass the same values on both
sides — e.g. `--r-min 0.005 --dr 0.01` here against `--rmin 0.005 --dr 0.01` there puts bin
centres at 0.01, 0.02, … for both.  The defaults are deliberately finer than that.

## Output

`gr[_<suffix>].csv`，**长表**：

| 列 | 含义 |
|---|---|
| `file` | 输入文件 stem |
| `r` | 半径 [Å] |
| `center` / `neighbor` | 该行属于哪一个有序对 |
| `gr` | $g(r)$ |
| `cn` | $CN(r)$ |

类型进**数据列**而不是列名：元素集不同的轨迹可直接堆叠而不需要对列，不给 `-a/-b`
是加行而不是加列。`gr` 对称（`A-B` 与 `B-A` 逐点相同），`cn` 有向
（`CN(A→B) = Σ hist/(N_A·steps)`）——这个区别写进了表结构而不是文档注脚。

所有产物是**一份** csv，多输入时堆叠成一张表并加 `file` 列；`#` 注释块里是共享参数与
`[inputs]` 清单（`pandas.read_csv(comment="#")` 会丢掉）。`-o` 给的是**文件名后缀**。


$g(r)$ is **symmetric** — `A-B` and `B-A` are pointwise identical.  $CN(r)$ is **directed** —
`A-B` is the average number of B around each A, so the order of `-a` / `-b` matters.

Header lines record all parameters, atom counts, average volume, and number density.

## Usage

```bash
ferro traj gr -i traj.dump -o run1

# one pair, narrowed range
ferro traj gr -a P -b O -i traj.dump --r-min 0.001 --r-max 15.0 --dr 0.002 -o po
```

```rust
use ferro_analysis::md::{GrParams, calc_gr, write_gr};

let params = GrParams::with_auto_rmax(&traj);
let result = calc_gr(&traj, &params).unwrap();
write_gr(&result, "output.gr", None).unwrap();
```

## Implementation Notes

- Parallelism: per-frame `par_iter` with `fold`/`reduce` accumulation.
- Pseudo-element labels (e.g. `"P0"`, `"Ob"`) are supported: `elem_z` resolves them by stripping numeric/alphabetic suffixes.
- Column ordering uses atomic number Z as the primary sort key, with string labels as tiebreaker for deterministic ordering of pseudo-elements.
