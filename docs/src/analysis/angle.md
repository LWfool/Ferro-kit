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

When both end atoms are the same element: $r_\text{cut} = \min(r_\text{cut,AB}, r_\text{cut,BC})$.

### Canonical Key Convention

For each triplet, the key `"ElemA-ElemCenter-ElemC"` is assigned with $Z(\text{ElemA}) \leq Z(\text{ElemC})$.  When both end atoms have the same Z, lexicographic order on the label string is used.  This ensures each angle type appears exactly once.

**No double-counting**: the enumeration requires $\text{idx}(A) < \text{idx}(C)$ to avoid counting (A,B,C) and (C,B,A) as separate events.

### Statistics

From the raw count histogram $h(\theta_i)$:

$$\bar{\theta} = \frac{\sum_i \theta_i \, h_i}{\sum_i h_i}, \quad \sigma = \sqrt{\frac{\sum_i (\theta_i - \bar{\theta})^2 h_i}{\sum_i h_i}}$$

## Parameters

```rust
pub struct AngleParams {
    pub r_cut_ab: f64,   // cutoff for lower-Z end to centre [Å]; default: 2.3
    pub r_cut_bc: f64,   // cutoff for higher-Z end to centre [Å]; default: 2.3
    pub d_angle: f64,    // histogram bin width [degrees]; default: 0.1
}
```

## Output

Columns: `angle[deg]`, then one column per triplet key in (Z_center, Z_left, Z_right) order.

The header lists mean, standard deviation, and total count per triplet.

## Usage

```bash
nex-traj -m angle -i traj.dump --rcut-ab 2.3 --rcut-bc 2.3 -o output.angle
```

```rust
use nexflux_analysis::md::{AngleParams, calc_angle, write_angle};

let params = AngleParams { r_cut_ab: 2.4, r_cut_bc: 2.4, d_angle: 0.5 };
let result = calc_angle(&traj, &params).unwrap();
write_angle(&result, "output.angle").unwrap();
```

## Implementation Notes

### Linked-Cell List

Brute-force neighbour search is O(N²) per frame.  nexflux uses a linked-cell list to reduce this to O(N·k) where $k$ is the average number of atoms in the search volume.

For each central atom B, only the 27 cells within ±1 cell in each direction are searched. The cell size is chosen so that a sphere of radius $r_\text{cut,max}$ fits within the search volume:

$$n_k = \left\lfloor \frac{1}{r_\text{cut,max} \cdot \|({\mathbf{M}^T})^{-1}\|_k} \right\rfloor$$

This is exact for triclinic cells.

### Parallelism

Computation is parallelised per frame with `rayon::par_iter().fold().reduce()`.
