//! `ChargeGrid` — charge density grid for Bader analysis (VASP CHGCAR convention).
//!
//! Stores `rho_stored = rho_physical × V_cell` with x-axis (N1) fastest,
//! matching CHGCAR file order: `rho[n1 + N1*(n2 + N2*n3)]`.

use nalgebra::{Matrix3, Vector3};
use crate::Cell;

/// Charge density grid with precomputed coordinate transforms for Bader analysis.
///
/// Density values follow the VASP convention: stored as `ρ_phys × V_cell` (units: e),
/// not normalized.  Key usage rules:
///   - Vacuum test: `rho_stored / V_cell <= vacval`
///   - Charge integral: `Σ rho_stored / nrho`
///   - Gradient ascent: compare `rho_stored` directly (proportional to physical density)
#[derive(Debug, Clone)]
pub struct ChargeGrid {
    /// Charge density values, `rho_stored = rho_physical × V_cell`.
    /// Indexed as `rho[n1 + N1*(n2 + N2*n3)]` (x fastest, z slowest).
    pub rho: Vec<f64>,
    /// Grid dimensions `[N1, N2, N3]`.
    pub shape: [usize; 3],
    /// Total number of grid points: `N1 * N2 * N3`.
    pub nrho: usize,
    /// Lattice-step → Cartesian transform.  Column `i` = `a_i / N_i` (Å per grid step).
    pub lat2car: Matrix3<f64>,
    /// Cartesian → lattice-step transform.  `= lat2car.try_inverse()`.
    pub car2lat: Matrix3<f64>,
    /// Precomputed distances for the 26+1 neighbors: `lat_dist[d1+1][d2+1][d3+1]` in Å.
    /// `d_i ∈ {-1, 0, 1}`; `(0,0,0)` entry = 0.0.
    pub lat_dist: [[[f64; 3]; 3]; 3],
    /// Inverse distances: `1 / lat_dist`.  `(0,0,0)` entry = 0.0 (not infinity).
    pub lat_i_dist: [[[f64; 3]; 3]; 3],
}

impl ChargeGrid {
    /// Construct a `ChargeGrid` from density data, grid dimensions, and the unit cell.
    ///
    /// `rho` must have length `N1 * N2 * N3`, stored with x (N1) fastest.
    /// The cell matrix rows are the lattice vectors a, b, c in Å.
    pub fn new(rho: Vec<f64>, shape: [usize; 3], cell: &Cell) -> Self {
        let [n1, n2, n3] = shape;
        let nrho = n1 * n2 * n3;
        assert_eq!(rho.len(), nrho, "rho length {} != N1*N2*N3 = {}", rho.len(), nrho);

        // lat2car: column i = a_i / N_i  (each lattice step in Cartesian)
        let m = cell.matrix;
        let lat2car = Matrix3::new(
            m[(0, 0)] / n1 as f64, m[(0, 1)] / n1 as f64, m[(0, 2)] / n1 as f64,
            m[(1, 0)] / n2 as f64, m[(1, 1)] / n2 as f64, m[(1, 2)] / n2 as f64,
            m[(2, 0)] / n3 as f64, m[(2, 1)] / n3 as f64, m[(2, 2)] / n3 as f64,
        );
        let car2lat = lat2car.try_inverse()
            .expect("lat2car is singular — cell matrix or grid dimensions are degenerate");

        // Precompute neighbor distances
        let mut lat_dist = [[[0.0_f64; 3]; 3]; 3];
        let mut lat_i_dist = [[[0.0_f64; 3]; 3]; 3];
        for di in -1i32..=1 {
            for dj in -1i32..=1 {
                for dk in -1i32..=1 {
                    let idx = [(di + 1) as usize, (dj + 1) as usize, (dk + 1) as usize];
                    let vec = lat2car * Vector3::new(di as f64, dj as f64, dk as f64);
                    let d = vec.norm();
                    lat_dist[idx[0]][idx[1]][idx[2]] = d;
                    lat_i_dist[idx[0]][idx[1]][idx[2]] = if d > 0.0 { 1.0 / d } else { 0.0 };
                }
            }
        }

        Self { rho, shape, nrho, lat2car, car2lat, lat_dist, lat_i_dist }
    }

    /// Periodic density query: folds `p` via PBC and returns `rho[p]`.
    #[inline]
    pub fn rho_val(&self, p: [i32; 3]) -> f64 {
        let [n1, n2, n3] = self.pbc_i(p);
        self.rho[n1 + self.shape[0] * (n2 + self.shape[1] * n3)]
    }

    /// Fold lattice coordinates into `[0, N_i)` via periodic boundary conditions.
    #[inline]
    pub fn pbc_i(&self, p: [i32; 3]) -> [usize; 3] {
        std::array::from_fn(|i| p[i].rem_euclid(self.shape[i] as i32) as usize)
    }

    /// Distance between center and neighbor at offset `(d1, d2, d3)` in Å.
    /// `d_i ∈ {-1, 0, 1}`.
    #[inline]
    pub fn lat_dist_i(&self, d1: i32, d2: i32, d3: i32) -> f64 {
        self.lat_dist[(d1 + 1) as usize][(d2 + 1) as usize][(d3 + 1) as usize]
    }

    /// Inverse distance at offset `(d1, d2, d3)`.  Returns 0.0 for `(0,0,0)`.
    #[inline]
    pub fn lat_i_dist_i(&self, d1: i32, d2: i32, d3: i32) -> f64 {
        self.lat_i_dist[(d1 + 1) as usize][(d2 + 1) as usize][(d3 + 1) as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Cell;

    fn cubic_cell(a: f64) -> Cell {
        Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap()
    }

    #[test]
    fn test_charge_grid_construction() {
        let cell = cubic_cell(10.0);
        let shape = [5, 5, 5];
        let rho = vec![1.0; 125];
        let chg = ChargeGrid::new(rho, shape, &cell);
        assert_eq!(chg.nrho, 125);
        assert_eq!(chg.shape, [5, 5, 5]);
    }

    #[test]
    fn test_rho_val_basic() {
        let cell = cubic_cell(10.0);
        let shape = [3, 3, 3];
        let mut rho = vec![0.0; 27];
        rho[0] = 42.0; // (0,0,0)
        rho[1] = 99.0; // (1,0,0)
        let chg = ChargeGrid::new(rho, shape, &cell);
        assert!((chg.rho_val([0, 0, 0]) - 42.0).abs() < 1e-15);
        assert!((chg.rho_val([1, 0, 0]) - 99.0).abs() < 1e-15);
    }

    #[test]
    fn test_rho_val_pbc() {
        let cell = cubic_cell(10.0);
        let shape = [3, 3, 3];
        let mut rho = vec![0.0; 27];
        rho[0] = 42.0; // (0,0,0)
        let chg = ChargeGrid::new(rho, shape, &cell);
        // PBC wrap: (-1, 0, 0) → (2, 0, 0)
        assert!((chg.rho_val([-1, 0, 0]) - chg.rho_val([2, 0, 0])).abs() < 1e-15);
        // PBC wrap: (3, 0, 0) → (0, 0, 0)
        assert!((chg.rho_val([3, 0, 0]) - 42.0).abs() < 1e-15);
    }

    #[test]
    fn test_pbc_i() {
        let cell = cubic_cell(10.0);
        let shape = [10, 10, 10];
        let rho = vec![0.0; 1000];
        let chg = ChargeGrid::new(rho, shape, &cell);
        assert_eq!(chg.pbc_i([0, 0, 0]), [0, 0, 0]);
        assert_eq!(chg.pbc_i([-1, -1, -1]), [9, 9, 9]);
        assert_eq!(chg.pbc_i([10, 10, 10]), [0, 0, 0]);
        assert_eq!(chg.pbc_i([5, -3, 12]), [5, 7, 2]);
    }

    #[test]
    fn test_lat_dist_cubic() {
        // Cubic 10 Å, 5×5×5 grid → step = 2 Å per axis
        let cell = cubic_cell(10.0);
        let shape = [5, 5, 5];
        let rho = vec![0.0; 125];
        let chg = ChargeGrid::new(rho, shape, &cell);
        // Face neighbor distance = 2.0 Å
        assert!((chg.lat_dist_i(1, 0, 0) - 2.0).abs() < 1e-10);
        assert!((chg.lat_dist_i(0, 1, 0) - 2.0).abs() < 1e-10);
        assert!((chg.lat_dist_i(0, 0, 1) - 2.0).abs() < 1e-10);
        // (0,0,0) distance = 0
        assert!((chg.lat_dist_i(0, 0, 0)).abs() < 1e-15);
        // Body diagonal: sqrt(3) * 2
        let expected = (3.0_f64).sqrt() * 2.0;
        assert!((chg.lat_dist_i(1, 1, 1) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_lat_i_dist_cubic() {
        let cell = cubic_cell(10.0);
        let shape = [5, 5, 5];
        let rho = vec![0.0; 125];
        let chg = ChargeGrid::new(rho, shape, &cell);
        // (0,0,0) → 0.0 (not infinity)
        assert!((chg.lat_i_dist_i(0, 0, 0)).abs() < 1e-15);
        // face neighbor → 1/2 = 0.5
        assert!((chg.lat_i_dist_i(1, 0, 0) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_index_consistency() {
        // Verify rho_val matches direct indexing
        let cell = cubic_cell(6.0);
        let shape = [3, 3, 3];
        let rho: Vec<f64> = (0..27).map(|i| i as f64).collect();
        let chg = ChargeGrid::new(rho, shape, &cell);
        for n3 in 0..3 {
            for n2 in 0..3 {
                for n1 in 0..3 {
                    let expected = (n1 + 3 * (n2 + 3 * n3)) as f64;
                    assert!((chg.rho_val([n1, n2, n3]) - expected).abs() < 1e-15);
                }
            }
        }
    }
}
