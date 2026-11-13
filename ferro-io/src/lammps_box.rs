//! LAMMPS's box convention, shared by the data and dump readers and writers.
//!
//! LAMMPS keeps a restricted triclinic box: a along x, b in the xy plane, the cell matrix
//! lower triangular, `[[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]]` in rows. The four files used
//! to build that matrix themselves, two of them with a verbatim copy of `cell_to_lammps` and
//! `bounding_box` (review 2026-10-30 A1-L13: one convention, so one place).

use ferro_core::{Cell, Frame};
use nalgebra::Matrix3;

/// The cell of LAMMPS box parameters, rows `[lx 0 0] [xy ly 0] [xz yz lz]`.
pub(crate) fn lammps_cell(lx: f64, ly: f64, lz: f64, xy: f64, xz: f64, yz: f64) -> Cell {
    Cell::from_matrix(Matrix3::new(lx, 0.0, 0.0, xy, ly, 0.0, xz, yz, lz))
}

/// LAMMPS box parameters `(lx, ly, lz, xy, xz, yz)`.
pub(crate) type BoxParams = (f64, f64, f64, f64, f64, f64);

/// A cell's LAMMPS box parameters, from its lengths and angles.
pub(crate) fn cell_to_lammps(cell: &Cell) -> BoxParams {
    let [a, b, c] = cell.lengths();
    let [alpha, beta, gamma] = cell.angles();
    let (al, be, ga) = (alpha.to_radians(), beta.to_radians(), gamma.to_radians());
    let lx = a;
    let xy = b * ga.cos();
    let xz = c * be.cos();
    // Both writers pass the cell through Cell::check first (volume > 1e-6 a·b·c), so
    // ly >= 1e-6 b and lz >= 1e-6 c: no square root below is of a negative number and
    // ly is no divisor near 0 (review 2026-11-10 B3-3)
    let ly = (b * b - xy * xy).sqrt();
    let yz = (b * c * al.cos() - xy * xz) / ly;
    let lz = (c * c - xz * xz - yz * yz).sqrt();
    // The cos of 90 degrees leaves ~1e-16 L, which would make an orthogonal box triclinic:
    // tilts below 1e-10 of the longest edge are set to 0
    let tol = 1e-10 * lx.max(ly).max(lz);
    let snap = |t: f64| if t.abs() < tol { 0.0 } else { t };
    let (xy, xz, yz) = (snap(xy), snap(xz), snap(yz));
    (lx, ly, lz, xy, xz, yz)
}

/// How far the box of a frame without a cell reaches past the atoms on each side (Å); the
/// coordinates are written as they are (dump and data both allow lo != 0). The margin has
/// no physical meaning: it keeps lo < hi (a flat molecule is 0 thick along one axis) and
/// keeps atoms off the hi face, where an `f` boundary would drop them.
const CELLLESS_PAD: f64 = 1.0;

/// The LAMMPS box the data and dump writers write for `frame`: its origin and
/// `(lx, ly, lz, xy, xz, yz)`. A frame without a cell gets the atoms' bounding box,
/// [`CELLLESS_PAD`] wider on each side (review 2026-11-06 D14: one copy, not one per writer).
pub(crate) fn frame_box(frame: &Frame) -> ([f64; 3], BoxParams) {
    match &frame.cell {
        Some(cell) => (cell.origin.into(), cell_to_lammps(cell)),
        None => {
            let mut mn = [f64::MAX; 3];
            let mut mx = [f64::MIN; 3];
            for a in frame.atoms() {
                for k in 0..3 {
                    mn[k] = mn[k].min(a.position[k]);
                    mx[k] = mx[k].max(a.position[k]);
                }
            }
            let l = |k: usize| mx[k] - mn[k] + 2.0 * CELLLESS_PAD;
            (mn.map(|m| m - CELLLESS_PAD), (l(0), l(1), l(2), 0.0, 0.0, 0.0))
        }
    }
}

/// The LAMMPS tilts `(xy, xz, yz)` of the same lattice on a reduced basis, each at most
/// half the edge it leans along — LAMMPS `read_data` refuses a box skewed further
/// (review R2-4). Each step adds a lattice vector to another: `c -= k·b` moves `yz` and
/// `xz`, then `c -= k·a` and `b -= k·a` move `xz` and `xy`. Only the box is reduced; the
/// atoms keep their Cartesian positions, which the same lattice still holds.
pub(crate) fn reduce_tilts(lx: f64, ly: f64, xy: f64, xz: f64, yz: f64) -> (f64, f64, f64) {
    // ly > 0: a cell passed `check` before either writer gets here, so b is not along a and
    // ly = b |sin gamma| > 0; a frame without one has ly >= 2 CELLLESS_PAD (the `else 0`
    // that stood here could not be reached, review 2026-11-06 A3-8)
    let k = (yz / ly).round();
    let (xz, yz) = (xz - k * xy, yz - k * ly);
    let xz = xz - (xz / lx).round() * lx;
    let xy = xy - (xy / lx).round() * lx;
    (xy, xz, yz)
}
