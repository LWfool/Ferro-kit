# Task Plan: ferro-structure — vacuum / merge / box builder

## Status: DONE

## Phases

### Phase 1: vacuum.rs
- [x] Create `vacuum.rs` with `add_vacuum` function
- [x] Input validation (cell exists, axis valid, thickness > 0)
- [x] Cell matrix scaling along specified axis
- [x] Update `lib.rs` to export vacuum module
- [x] Unit tests (each axis, error cases, pbc unchanged) — 12 tests, all green
- [x] `cargo clippy` clean

### Phase 2: merge.rs
- [x] Create `merge.rs` with `merge_frames` function
- [x] Input validation (both cells exist, axis valid, gap >= 0)
- [x] Cell dimension calculation (join axis = sum + gap, others = max)
- [x] Atom translation and centering
- [x] Update `lib.rs` to export merge module
- [x] Unit tests (cubic cells, different sizes, triclinic, error cases) — 13 tests, all green
- [x] `cargo clippy` clean

### Phase 3: box_builder.rs
- [x] Create `box_builder.rs` with `Component` struct
- [x] Formula parser (simple formulas: C, H2O, P2O5, ZnO, CH3OH)
- [x] `estimate_box_length` function
- [x] Random atom placement
- [x] Soft-core relaxation with cell list
- [x] `build_box` function (full pipeline)
- [x] Update `lib.rs` to export box_builder module
- [x] Unit tests (formula parsing, box estimation, relaxation convergence) — 18 tests, all green
- [x] `cargo clippy` clean

### Phase 4: Integration
- [ ] `cargo test --package ferro-structure` all green
- [ ] Update CLAUDE.md if needed
- [ ] Commit

## Decisions

| Decision | Choice | Reason |
|----------|--------|--------|
| vacuum axis type | `&str` ("x"/"y"/"z") | User preference |
| vacuum pbc | Don't change | User preference |
| vacuum centering | Don't center atoms | User preference |
| merge gap | User-specified Å | User preference |
| merge cell other axes | max of two structures | Prevents overlap |
| box shape | Cubic | User preference |
| min distance | Unified (not per-element pair) | User preference |
| placement strategy | Random + soft-core relaxation | User preference |
| relaxation steps | ~100, step_size 0.1 Å | Fast enough for 10k atoms |

## Dependencies

- No new external crates needed
- Uses: `ferro_core::frame::Frame`, `ferro_core::cell::Cell`, `ferro_core::atom::Atom`
- Uses: `ferro_core::data::compounds` for box_builder
- Uses: `nalgebra::Vector3`
