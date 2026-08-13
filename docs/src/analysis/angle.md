# Bond Angle Distribution

## Theory

The bond angle distribution $P(\theta)$ describes the probability of observing a three-body angle A–B–C, where B is the central atom.  It is a key structural descriptor for network-forming glasses (e.g. P–O–P, O–P–O angles in phosphate glasses).

### Angle Calculation

For a triplet (A, B, C) with B as centre:

$$\theta = \arccos\!\left(\frac{\overrightarrow{BA} \cdot \overrightarrow{BC}}{|\overrightarrow{BA}||\overrightarrow{BC}|}\right)$$

The vectors $\overrightarrow{BA}$ and $\overrightarrow{BC}$ use the minimum-image convention for periodic cells.

### Selection Criteria

An atom pair (A, B) forms a bond if $|\overrightarrow{BA}| < r_\text{cut,AB}$.  
Similarly for (B, C): $|\overrightarrow{BC}| < r_\text{cut,BC}$.

**Which physical end is "A"** is decided by the order *you* wrote it. When the triplet is
named on the command line, `--r-cut-ab` goes to the atom given by `-a` / `-x` and
`--r-cut-bc` to the one given by `-c` / `-z`, regardless of atomic number:

```bash
# Zn gets 2.5, P gets 2.0 — exactly as written
ferro traj angle -a Zn -b O -c P -i traj.dump --r-cut-ab 2.5 --r-cut-bc 2.0
```

> **Changed in v0.1.15.** Earlier versions handed the two cutoffs out by atomic number,
> so the command above silently gave 2.5 to P (Z = 15 < Z(Zn) = 30) — the opposite of
> what was written. Results are unaffected whenever the two cutoffs are equal.

Without a named triplet (scanning every triplet at once) there is no user order to follow,
so the canonical (Z, symbol) order is used: `r_cut_ab` → lower-Z end, `r_cut_bc` → higher-Z end.

When both end atoms are the same type, which one gets sorted first depends on the neighbour
enumeration order, so $r_\text{cut} = \min(r_\text{cut,AB}, r_\text{cut,BC})$ is applied
to both — the only assignment that gives a reproducible result.

### Canonical Key Convention

For each triplet, the key `"ElemA-ElemCenter-ElemC"` is assigned with $Z(\text{ElemA}) \leq Z(\text{ElemC})$.  When both end atoms have the same Z, lexicographic order on the label string is used.  This ensures each angle type appears exactly once.

**No double-counting**: the enumeration requires $\text{idx}(A) < \text{idx}(C)$ to avoid counting (A,B,C) and (C,B,A) as separate events. Each *geometric* angle is therefore counted once — a PO₄ tetrahedron contributes 6 O–P–O angles, not 12.

### Histogram Range

Counts land in bins of width `d_angle` starting at `angle_min`; bin $i$ is labelled at its
centre, $\theta_i = \text{angle\_min} + (i + \tfrac{1}{2})\,\Delta\theta$. The upper edge is
**closed**: a perfectly collinear triplet gives $\arccos(-1) = 180.0$ exactly, which a
half-open interval would discard entirely.

Angles outside $[\text{angle\_min}, \text{angle\_max}]$ are dropped, so narrowing the window
also narrows what is counted — it is not merely a plotting range.

### Comparison with `dump2analysis -m angle`

Two convention differences, both documented here so the numbers can be read side by side:

| | ferro | `dump2analysis` |
|---|---|---|
| Counting | each geometric angle once | ordered end pairs — **twice** when both ends are the same element (once when they differ) |
| Bin $k$ covers | $[k\Delta\theta,\ (k{+}1)\Delta\theta)$ | $[0.05 + k\Delta\theta,\ 0.15 + k\Delta\theta)$ |

The factor of two cancels out of the mean, the standard deviation and the peak positions.
The bin offset is half a bin; `dump2analysis`'s grid misses $[0, 0.05)$ (angles there index
`a[-1]`, out of bounds) and extends past 180°, so ferro does not adopt it as a default. It
can be reproduced exactly when needed:

```bash
ferro traj angle -a O -b P -c O -i traj.dump \
        --angle-min 0.05 --angle-max 180.05 --d-angle 0.1
```

Verified on `examples/43Z43P15A_NPT.lammpstrj`: after accounting for the ×2 counting
convention, all 1800 bins of both O–P–O and O–Al–O match `dump2analysis` as **exact
integers** (zero differing bins). See `scripts/compare_angle.py --align-binning`.

### Statistics

From the raw count histogram $h(\theta_i)$:

$$\bar{\theta} = \frac{\sum_i \theta_i \, h_i}{\sum_i h_i}, \quad \sigma = \sqrt{\frac{\sum_i (\theta_i - \bar{\theta})^2 h_i}{\sum_i h_i}}$$

## Parameters

```rust
pub struct AngleParams {
    pub r_cut_ab: f64,   // cutoff from end A to centre [Å]; default: 2.3
    pub r_cut_bc: f64,   // cutoff from end C to centre [Å]; default: 2.3
    pub angle_min: f64,  // histogram lower edge [degrees]; default: 0.0
    pub angle_max: f64,  // histogram upper edge [degrees]; default: 180.0
    pub d_angle: f64,    // histogram bin width [degrees]; default: 0.1
    pub group_by: GroupBy,       // resolve triplets over elements or site labels
    pub ends: Option<(String, String)>,  // (A, C) as the caller wrote them; see above
}
```

| CLI flag | Field | Default |
|---|---|---|
| `--r-cut-ab` | `r_cut_ab` | 2.3 |
| `--r-cut-bc` | `r_cut_bc` | 2.3 |
| `--angle-min` | `angle_min` | 0.0 |
| `--angle-max` | `angle_max` | 180.0 |
| `--d-angle` | `d_angle` | 0.1 |
| `-a`/`-b`/`-c` (element) or `-x`/`-y`/`-z` (label) | `group_by`, `ends` | all triplets |

## Output

`angle_<三元组>[_<suffix>].csv`（未点名三元组时为 `angle_all…`），**长表**：`file, angle, end_a, center, end_c, count, p`。

三元组进**数据列**，故元素集不同的轨迹可直接堆叠。同时保留整数 `count` 与归一化 `p`
两列：整数直方图是与 `dump2analysis` 逐 bin 对拍的依据，只留 `p` 就对不了。

所有产物是**一份** csv，多输入时堆叠成一张表并加 `file` 列；`#` 注释块里是共享参数与
`[inputs]` 清单（`pandas.read_csv(comment="#")` 会丢掉）。`-o` 给的是**文件名后缀**。


## Usage

```bash
ferro traj angle -i traj.dump --r-cut-ab 2.3 --r-cut-bc 2.3 -o run1

# one triplet, narrowed to the tetrahedral region
ferro traj angle -a O -b P -c O -i traj.dump \
        --angle-min 90 --angle-max 130 --d-angle 0.2 -o opo
```

```rust
use ferro_analysis::md::{AngleParams, calc_angle, write_angle};

let params = AngleParams {
    r_cut_ab: 2.4,
    r_cut_bc: 2.4,
    d_angle: 0.5,
    ..AngleParams::default()
};
let result = calc_angle(&traj, &params).unwrap();
write_angle(&result, "output.angle").unwrap();
```

## Implementation Notes

### Linked-Cell List

Brute-force neighbour search is O(N²) per frame.  ferro uses a linked-cell list to reduce this to O(N·k) where $k$ is the average number of atoms in the search volume.

For each central atom B, only the 27 cells within ±1 cell in each direction are searched. The cell size is chosen so that a sphere of radius $r_\text{cut,max}$ fits within the search volume:

$$n_k = \left\lfloor \frac{1}{r_\text{cut,max} \cdot \|({\mathbf{M}^T})^{-1}\|_k} \right\rfloor$$

This is exact for triclinic cells.

### Parallelism

Computation is parallelised per frame with `rayon::par_iter().fold().reduce()`.
