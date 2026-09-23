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

`gr_<pair>[_<suffix>].csv` (`gr_all…` when no pair is named), **long table**:

| Column | Meaning |
|---|---|
| `file` | Input file stem |
| `r` | Radius [Å] |
| `center` / `neighbor` | Which ordered pair the row belongs to |
| `gr` | $g(r)$ |
| `cn` | $CN(r)$ |

Types go into **data columns** rather than column names: trajectories with different element sets stack
directly without aligning columns, and omitting `-a/-b` adds rows rather than columns.  `gr` is symmetric
(`A-B` and `B-A` are identical point by point), `cn` is directed (`CN(A→B) = Σ hist/(N_A·steps)`) — that distinction is written into the table structure, not into a footnote.

All output is **one** csv; multiple inputs are stacked into a single table with a `file` column.  The `#` comment block holds the shared parameters and the `[inputs]` list
(`pandas.read_csv(comment="#")` drops it).  `-o` takes the **output directory**; the batch suffix goes to `-s`.


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

## Selecting by site label holds over a single frame only

`-x`/`-y` resolve partials over `Atom::label` instead of the element. That works
on one frame and is rejected on a labelled trajectory: `calc_gr` guards the
per-type particle count frame by frame, and labels shift as the run evolves — on
the reference trajectory `P_3` goes 149 / 152 / 150 / 150 / 150 over five frames.

The guard is not a limitation of the implementation but of the quantity: a
partial g(r) normalises by $N_A N_B$, and a count that changes between frames
has no single value to normalise by. Use `--last-n 1` for one frame, or select
by element with `-a`/`-b`.

## Implementation Notes

- Parallelism: per-frame `par_iter` with `fold`/`reduce` accumulation.
- Pseudo-element labels (e.g. `"P0"`, `"Ob"`) are supported: `elem_z` resolves them by stripping numeric/alphabetic suffixes.
- Column ordering uses atomic number Z as the primary sort key, with string labels as tiebreaker for deterministic ordering of pseudo-elements.
