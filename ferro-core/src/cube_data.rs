//! `CubeData` — complete data carrier for a Gaussian cube file.
//!
//! A cube file contains two inseparable parts: the atomic structure (single frame) and
//! 3-D volumetric data (electron density, wavefunction, etc.).
//! A dedicated struct is used to avoid polluting the general-purpose [`Frame`].

use nalgebra::{Matrix3, Vector3};
use crate::Frame;

/// Container for a Gaussian cube file: atomic structure plus volumetric grid data.
///
/// All length quantities are in **Ångström** (internal standard units).
///
/// The volumetric data is stored as a flat `Vec<f64>` in C-order (row-major),
/// axis order X → Y → Z (outer loop X, middle loop Y, inner loop Z),
/// matching the cube format convention.
/// Use [`idx`](CubeData::idx) to convert (ix, iy, iz) to a flat index.
#[derive(Debug, Clone)]
pub struct CubeData {
    /// Atomic structure: atoms, cell, pbc (single frame)
    pub frame: Frame,
    /// Volumetric data, flat C-order (X outer, Z inner), length = nx × ny × nz
    pub data: Vec<f64>,
    /// Grid dimensions `[nx, ny, nz]`
    pub shape: [usize; 3],
    /// Grid origin \[Å\]
    pub origin: Vector3<f64>,
    /// Voxel vectors \[Å\]: row `i` is the step vector along grid axis `i`
    pub spacing: Matrix3<f64>,
}

impl CubeData {
    /// Grid dimensions as `(nx, ny, nz)`.
    pub fn shape(&self) -> (usize, usize, usize) {
        (self.shape[0], self.shape[1], self.shape[2])
    }

    /// Flat index for voxel `(ix, iy, iz)` in C-order layout.
    #[inline]
    pub fn idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix * self.shape[1] * self.shape[2] + iy * self.shape[2] + iz
    }

    /// Read the value at voxel `(ix, iy, iz)`.
    #[inline]
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.data[self.idx(ix, iy, iz)]
    }

    /// Mutable reference to the value at voxel `(ix, iy, iz)`.
    #[inline]
    pub fn get_mut(&mut self, ix: usize, iy: usize, iz: usize) -> &mut f64 {
        let i = self.idx(ix, iy, iz);
        &mut self.data[i]
    }
}
