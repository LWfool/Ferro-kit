# Findings

## Context

- ferro-structure currently has only `supercell.rs` (done, 12 tests, clippy clean)
- `lib.rs` re-exports supercell functions only
- CLI and Python not yet integrated with ferro-structure
- No Python reference code for vacuum/merge/box_builder
- `CompoundData` in `ferro-core/src/data/compounds.rs` has: name, formula, molecular_mass, density, CAS
- COMPOUNDS formulas are all simple (no parentheses): H2O, P2O5, ZnO, CH3OH, C2H5OH, etc.
- `elements.rs`: `by_symbol(symbol) -> Option<&ElementData>` with `atomic_mass: f64`
- `compounds.rs`: `find(query) -> Option<&CompoundData>` by name or formula (case-insensitive)
- Cell matrix: row vectors a, b, c; row index 0=x, 1=y, 2=z

## Key patterns from supercell.rs

- Functions take `&Frame`, return `Result<Frame>`
- Validate cell exists: `frame.cell.as_ref().ok_or_else(|| ChemError::ValidationError(...))`
- Validate dimensions: check for zero
- `energy/forces/stress/velocities` set to `None` on output
- `bonds` preserved with index remapping where needed
- `charge` scaled or combined, `multiplicity` reset to 1 where appropriate
