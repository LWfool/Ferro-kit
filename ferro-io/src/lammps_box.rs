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

/// A cell's LAMMPS box parameters `(lx, ly, lz, xy, xz, yz)`, from its lengths and angles.
pub(crate) fn cell_to_lammps(cell: &Cell) -> (f64, f64, f64, f64, f64, f64) {
    let [a, b, c] = cell.lengths();
    let [alpha, beta, gamma] = cell.angles();
    let (al, be, ga) = (alpha.to_radians(), beta.to_radians(), gamma.to_radians());
    let lx = a;
    let xy = b * ga.cos();
    let xz = c * be.cos();
    let ly = (b * b - xy * xy).max(0.0).sqrt();
    let yz = if ly > 1e-10 { (b * c * al.cos() - xy * xz) / ly } else { 0.0 };
    let lz = (c * c - xz * xz - yz * yz).max(0.0).sqrt();
    // The cos of 90 degrees leaves ~1e-16 L, which would make an orthogonal box triclinic:
    // tilts below 1e-10 of the longest edge are set to 0
    let tol = 1e-10 * lx.max(ly).max(lz);
    let snap = |t: f64| if t.abs() < tol { 0.0 } else { t };
    let (xy, xz, yz) = (snap(xy), snap(xz), snap(yz));
    (lx, ly, lz, xy, xz, yz)
}

/// The atoms' bounding box `(xlo, xhi, ylo, yhi, zlo, zhi)`, for a frame without a cell.
pub(crate) fn bounding_box(frame: &Frame) -> (f64, f64, f64, f64, f64, f64) {
    let mut mn = [f64::MAX; 3];
    let mut mx = [f64::MIN; 3];
    for a in frame.atoms() {
        for k in 0..3 {
            mn[k] = mn[k].min(a.position[k]);
            mx[k] = mx[k].max(a.position[k]);
        }
    }
    (mn[0], mx[0], mn[1], mx[1], mn[2], mx[2])
}
