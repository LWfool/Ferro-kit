# ferro-structure: vacuum / merge / box builder

Date: 2026-05-04

## Overview

Three new modules for `ferro-structure`, following the existing `supercell.rs` pattern (single-frame API, return `Result<Frame>`, `energy/forces/stress/velocities` not copied).

## 1. `vacuum.rs` — Add vacuum layer

### Signature

```rust
pub fn add_vacuum(frame: &Frame, axis: &str, thickness: f64) -> Result<Frame>
```

### Behavior

- `axis`: "x" / "y" / "z" → mapped to cell matrix row index (0/1/2)
- `thickness`: must be > 0, in Å
- The specified axis's cell row vector is scaled by `(L + thickness) / L`
- Atomic positions unchanged
- `pbc` unchanged
- `bonds` preserved (indices unchanged)
- `energy`, `forces`, `stress`, `velocities` set to `None` (copied from source, not meaningful for modified structure)
- `charge`, `multiplicity` preserved

### Errors

- `ValidationError` if frame has no cell
- `ValidationError` if axis is not "x"/"y"/"z"
- `ValidationError` if thickness <= 0

## 2. `merge.rs` — Merge two frames along an axis

### Signature

```rust
pub fn merge_frames(frame_a: &Frame, frame_b: &Frame, axis: &str, gap: f64) -> Result<Frame>
```

### Behavior

- Both frames must have cells
- `axis`: "x" / "y" / "z"
- `gap`: vacuum gap at interface (Å), must be >= 0
- Merged cell along join axis: `L_a + gap + L_b`
- Other two axes: `max(L_a, L_b)`; the shorter structure is centered on that axis
- `frame_b` atoms translated along join axis by `L_a + gap`
- `frame_b` atoms centered on the other two axes (offset = `(max_L - L_b) / 2`)
- Bonds not merged (indices shift)
- `energy/forces/stress/velocities` → `None`
- `charge` = sum of both frames' charges
- `multiplicity` = 1

### Errors

- `ValidationError` if either frame has no cell
- `ValidationError` if axis is not "x"/"y"/"z"
- `ValidationError` if gap < 0

## 3. `box_builder.rs` — Mixed-compound box construction

### Data types

```rust
pub struct Component {
    pub compound: String,    // name or formula, looked up in COMPOUNDS database
    pub n_molecules: usize,
}
```

### Functions

```rust
/// Estimate cubic box edge length (Å) from compound list and target density.
pub fn estimate_box_length(components: &[Component], density: f64) -> Result<f64>

/// Build a mixed-compound box: estimate size → random placement → soft-core relaxation.
pub fn build_box(components: &[Component], density: f64, min_dist: f64) -> Result<Frame>
```

### `estimate_box_length`

- Look up each compound in `ferro_core::data::compounds::COMPOUNDS`
- `total_mass = Σ(n_i × M_i)` (g/mol)
- `V = total_mass / (density × N_A)` (cm³) → convert to Å³ (× 10^24)
- `L = V^(1/3)` for cubic box

### `build_box`

1. Parse each compound's formula → element stoichiometry
   - e.g., P₂O₅ with n=100 → {P: 200, O: 500}
   - Formula parser: handle single-letter (C, N, O) and two-letter (Zn, Fe, Cl) elements, optional count
2. Aggregate element counts across all components
3. Call `estimate_box_length` for cubic L
4. Create cell: `Cell::from_lengths_angles(L, L, L, 90, 90, 90)`
5. Place atoms randomly: each atom gets `position = random::<Vector3<f64>>() * L` (uniform in [0, L))
6. Soft-core relaxation (steepest descent):
   - For each step (default 100):
     - For each atom i, compute net repulsive force from atoms j within `min_dist`:
       - `r_ij = |pos_j - pos_i|` (minimum image convention)
       - If `r_ij < min_dist`: `force_i += (1 - r_ij/min_dist) * direction_ij`
     - Update: `pos_i += step_size * force_i`
     - `step_size` starts at 0.1 Å, can be constant
   - Wrap positions into cell after each step
   - Early exit if max displacement < threshold
7. Set `pbc = [true; 3]`, `charge = 0`, `multiplicity = 1`
8. Return Frame

### Performance

- Cell list for neighbor search: divide box into bins of size `min_dist`
- O(N × steps) with cell list, O(N² × steps) without
- 5000–10000 atoms: estimated < 5 seconds with cell list

### Errors

- `ValidationError` if compound not found in COMPOUNDS
- `ValidationError` if density <= 0
- `ValidationError` if any component has n_molecules = 0
- `ValidationError` if min_dist <= 0

## Testing strategy

### vacuum.rs
- Add vacuum to cubic cell, verify cell length increased correctly
- Add vacuum along each axis (x/y/z)
- Verify atomic positions unchanged
- Verify pbc unchanged
- Error: no cell, invalid axis, negative thickness

### merge.rs
- Merge two cubic cells, verify merged cell dimensions
- Verify frame_b atoms translated correctly
- Verify shorter structure centered on non-join axes
- Error: no cell on either frame, invalid axis, negative gap

### box_builder.rs
- `estimate_box_length`: known compound (water), verify volume calculation
- Formula parser: test various formulas (H2O, P2O5, ZnO, C6H12O6)
- `build_box`: verify atom count matches stoichiometry
- `build_box`: verify no atoms closer than `min_dist` after relaxation
- `build_box`: verify all atoms within cell after relaxation
- Integration test: build water box, verify reasonable density

## File structure

```
ferro-structure/src/
├── lib.rs              # add mod vacuum, merge, box_builder; re-export
├── supercell.rs        # existing
├── vacuum.rs           # new
├── merge.rs            # new
└── box_builder.rs      # new
```

## Dependencies

No new external dependencies needed. Uses existing `ferro-core` types and `nalgebra`.

`rayon` may be added later for parallel relaxation if performance needs it, but single-threaded should suffice for the initial implementation.
