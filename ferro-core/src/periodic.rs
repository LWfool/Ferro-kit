//! Periodic geometry of a 3-D periodic frame: the minimum image.
//!
//! Every analysis that measures distances in a box goes through [`PeriodicCell`].
//! Building one is where the box is validated — singular, non-finite and not
//! fully periodic cells are refused there — so the hot loops that use it cannot
//! fail. Only 3-D periodic systems are supported: a frame that carries a box but
//! whose `pbc` is not `[true; 3]` (an LAMMPS `f`/`s` boundary, CP2K `PERIODIC XY`,
//! a DeePMD `nopbc` system) is an error, never quietly treated as a bulk.

use nalgebra::{Matrix3, Vector3};

use crate::cell::Cell;
use crate::error::{ChemError, Result};
use crate::frame::Frame;

/// A validated 3-D periodic cell with its inverse computed once.
#[derive(Debug, Clone)]
pub struct PeriodicCell {
    /// `Mᵀ`: fractional → Cartesian (`cart = Mᵀ·frac`, rows of `M` are lattice vectors)
    mt: Matrix3<f64>,
    /// `(Mᵀ)⁻¹`: Cartesian → fractional; its rows are the reciprocal lattice vectors
    inv: Matrix3<f64>,
}

impl PeriodicCell {
    /// Validates `cell` and `pbc`: the cell must pass [`Cell::check`] and all three
    /// directions be periodic.
    pub fn new(cell: &Cell, pbc: [bool; 3]) -> Result<Self> {
        if pbc != [true; 3] {
            let flag = |p: bool| if p { 'T' } else { 'F' };
            return Err(ChemError::ValidationError(format!(
                "only 3-D periodic cells are supported, got pbc = \"{} {} {}\"",
                flag(pbc[0]), flag(pbc[1]), flag(pbc[2])
            )));
        }
        cell.check("cell matrix")?;
        let mt = cell.matrix.transpose();
        let inv = mt
            .try_inverse()
            .ok_or_else(|| ChemError::ValidationError("cell matrix is singular".into()))?;
        Ok(Self { mt, inv })
    }

    /// The frame's cell; `Err` when it has none or the cell is refused by [`Self::new`].
    pub fn of(frame: &Frame) -> Result<Self> {
        let cell = frame.cell.as_ref().ok_or_else(|| {
            ChemError::ValidationError("frame has no periodic cell".into())
        })?;
        Self::new(cell, frame.pbc)
    }

    pub fn to_fractional(&self, cart: Vector3<f64>) -> Vector3<f64> {
        self.inv * cart
    }

    pub fn to_cartesian(&self, frac: Vector3<f64>) -> Vector3<f64> {
        self.mt * frac
    }

    /// Minimum-image form of a displacement: each fractional component in [-0.5, 0.5].
    /// Exact for distances below [`Self::mic_bound`].
    ///
    /// Kept apart from `Cell::minimum_image` (fallible, inverts per call) and
    /// `box_builder::minimum_image_ortho` (cubic only) on purpose: the contracts differ
    /// (slim scan A-11).
    pub fn minimum_image(&self, diff: Vector3<f64>) -> Vector3<f64> {
        let f = self.inv * diff;
        self.mt * f.map(|c| c - c.round())
    }

    /// `(Mᵀ)⁻¹`; row `i` is reciprocal lattice vector `i` (without 2π), of norm `1/dᵢ`.
    pub fn reciprocal(&self) -> &Matrix3<f64> {
        &self.inv
    }

    /// Interplanar spacings d₁, d₂, d₃ \[Å\]; see [`Cell::interplanar_spacings`].
    pub fn spacings(&self) -> [f64; 3] {
        [0, 1, 2].map(|i| 1.0 / self.inv.row(i).norm())
    }

    /// Largest distance the minimum image measures correctly: half the smallest
    /// interplanar spacing. Same formula as `Cell::minimum_image_cutoff`, whose one
    /// production caller (`gr`) could move here (slim scan A-11).
    pub fn mic_bound(&self) -> f64 {
        let d = self.spacings();
        d[0].min(d[1]).min(d[2]) * 0.5
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cubic(a: f64) -> Cell {
        Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap()
    }

    #[test]
    fn only_finite_invertible_3d_periodic_cells_are_accepted() {
        assert!(PeriodicCell::new(&cubic(10.0), [true; 3]).is_ok());
        let slab = PeriodicCell::new(&cubic(10.0), [true, true, false]).unwrap_err();
        assert!(slab.to_string().contains("pbc = \"T T F\""), "{slab}");
        assert!(PeriodicCell::new(&cubic(10.0), [false; 3]).is_err());

        let mut flat = cubic(10.0);
        flat.matrix[(2, 2)] = 0.0;
        assert!(PeriodicCell::new(&flat, [true; 3]).unwrap_err().to_string().contains("singular"));
        let mut nan = cubic(10.0);
        nan.matrix[(0, 1)] = f64::NAN;
        assert!(PeriodicCell::new(&nan, [true; 3]).unwrap_err().to_string().contains("not finite"));

        let mut frame = Frame::new();
        assert!(PeriodicCell::of(&frame).is_err());
        frame.cell = Some(cubic(10.0));
        frame.pbc = [true; 3];
        assert!(PeriodicCell::of(&frame).is_ok());
    }

    #[test]
    fn minimum_image_matches_cell() {
        let cell = Cell::from_lengths_angles(7.0, 9.0, 11.0, 75.0, 85.0, 95.0).unwrap();
        let pc = PeriodicCell::new(&cell, [true; 3]).unwrap();
        let d = Vector3::new(6.1, -8.2, 4.4);
        assert_eq!(pc.minimum_image(d), cell.minimum_image(d).unwrap());
        assert_eq!(pc.mic_bound(), cell.minimum_image_cutoff().unwrap());
        let f = Vector3::new(0.3, 0.5, 0.7);
        assert!((pc.to_fractional(pc.to_cartesian(f)) - f).norm() < 1e-12);
    }
}
