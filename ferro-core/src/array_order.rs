//! Row-major (C order) conversion at the memory boundary.
//!
//! Every array that leaves Rust memory — a `.npy` file, an extxyz `Lattice=`
//! string, a DeePMD `box.npy` row — is written in **row-major order**, and
//! every matrix in `ferro-core` means row = lattice vector (or row = tensor
//! row). That convention is settled; what is NOT settled by the type system is
//! the physical layout: `nalgebra` stores `Matrix3` **column-major** and offers
//! no switch, so `Matrix3::as_slice()` yields `XX YX ZX XY YY ZY XZ YZ ZZ` —
//! the transpose of what every file format expects.
//!
//! Indexing (`m[(i, j)]`) hides the layout; only `as_slice()` leaks it. So the
//! rule is: never call `as_slice()` on a nalgebra type for I/O — call
//! [`matrix3_row_major`] instead, and keep this the single place where the
//! transposition is decided.
//!
//! Half of that mistake is silent: a symmetric tensor (stress) is unchanged by
//! transposition, so a wrong `as_slice()` shows up on the cell matrix and never
//! on the stress — one gets caught, the other ships. The test below therefore
//! uses an asymmetric matrix on purpose.

use nalgebra::Matrix3;

/// Flattens a 3×3 matrix in row-major order: `XX XY XZ YX YY YZ ZX ZY ZZ`.
///
/// This is the order used by DeePMD's `box.npy` / `virial.npy`, by the extxyz
/// `Lattice=` keyword and by GPUMD's `lattice=` — i.e. by every format ferro
/// reads or writes.
pub fn matrix3_row_major(m: &Matrix3<f64>) -> [f64; 9] {
    [
        m[(0, 0)], m[(0, 1)], m[(0, 2)],
        m[(1, 0)], m[(1, 1)], m[(1, 2)],
        m[(2, 0)], m[(2, 1)], m[(2, 2)],
    ]
}

/// Rebuilds a 3×3 matrix from nine row-major values.
///
/// The inverse of [`matrix3_row_major`]; `Matrix3::from_row_slice` already does
/// this, but going through a named function keeps both directions of the
/// convention visible in one module.
pub fn matrix3_from_row_major(v: &[f64; 9]) -> Matrix3<f64> {
    Matrix3::from_row_slice(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The matrix is asymmetric on purpose: a symmetric one cannot detect a
    /// transposition, which is exactly why this bug survives on stress tensors.
    fn asymmetric() -> Matrix3<f64> {
        Matrix3::from_row_slice(&[
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
        ])
    }

    #[test]
    fn row_major_follows_xx_xy_xz_order() {
        let m = asymmetric();
        assert_eq!(
            matrix3_row_major(&m),
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
        );
    }

    /// Pins the difference from `as_slice()` so nobody "simplifies" the helper
    /// away: nalgebra's own slice is the transpose of what files want.
    #[test]
    fn nalgebra_as_slice_is_the_transpose() {
        let m = asymmetric();
        assert_eq!(
            m.as_slice(),
            &[1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 9.0]
        );
        assert_ne!(matrix3_row_major(&m).as_slice(), m.as_slice());
    }

    #[test]
    fn round_trip_is_identity() {
        let m = asymmetric();
        assert_eq!(matrix3_from_row_major(&matrix3_row_major(&m)), m);
    }
}
