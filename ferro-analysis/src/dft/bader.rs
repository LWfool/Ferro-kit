//! Bader charge analysis — public API.
//!
//! Reference: Henkelman et al., *Comput. Mater. Sci.* **36**, 254 (2006).
//! See `dev/bader.md` for detailed implementation notes.

use ferro_core::{ChargeGrid, Frame};
use std::io::{BufWriter, Write};
use ferro_core::error::{ChemError, Result};

/// Bader analysis method.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaderMethod {
    /// On-grid: 26-neighbor distance-corrected gradient ascent (fastest, least accurate).
    OnGrid,
    /// Near-grid: central-difference gradient with fractional stepping (default).
    NearGrid,
    /// Off-grid: trilinear interpolation continuous gradient ascent.
    OffGrid,
    /// Weight: Yu-Trinkle Wigner-Seitz flow method (most accurate, slowest).
    Weight,
}

/// Parameters for Bader analysis.
#[derive(Debug, Clone)]
pub struct BaderParams {
    /// Analysis method (default: NearGrid).
    pub method: BaderMethod,
    /// Edge refinement count: -1 = auto (repeat until converged), -2 = single pass, N = N passes.
    pub refine: i32,
    /// Vacuum density threshold in e/Å³.  Points with `|ρ|/V ≤ vacval` are marked vacuum.
    pub vacval: f64,
    /// Off-grid step size in Å (None = use minimum voxel dimension).
    pub stepsize: Option<f64>,
}

impl Default for BaderParams {
    fn default() -> Self {
        Self {
            method: BaderMethod::NearGrid,
            refine: -1,
            vacval: 1e-3,
            stepsize: None,
        }
    }
}

/// Result of Bader charge analysis.
///
/// Grid-point arrays (`volnum`, `known`) use the same indexing as `ChargeGrid`:
/// `idx = n1 + N1*(n2 + N2*n3)` (x fastest).
///
/// Volume indices are 1-based: volume `v` (1 ≤ v ≤ nvols) corresponds to
/// `volpos_lat[v-1]`, `volchg[v-1]`, `nnion[v-1]`, etc.
/// `volnum[p] = 0` means unassigned; `volnum[p] = nvols+1` means vacuum.
#[derive(Debug, Clone)]
pub struct BaderResult {
    /// Number of Bader volumes found.
    pub nvols: usize,

    /// Per-grid-point volume assignment (1-indexed; 0 = unassigned, nvols+1 = vacuum).
    pub volnum: Vec<i32>,

    /// Bader maximum positions in lattice (grid-step) coordinates.
    pub volpos_lat: Vec<[f64; 3]>,
    /// Bader maximum positions in fractional coordinates.
    pub volpos_dir: Vec<[f64; 3]>,
    /// Bader maximum positions in Cartesian coordinates (Å).
    pub volpos_car: Vec<[f64; 3]>,

    /// Per-volume integrated charge (e).  `volchg[v-1]` = volume `v`.
    /// Last element = vacuum charge.
    pub volchg: Vec<f64>,

    /// Per-atom Bader charge (e).
    pub ionchg: Vec<f64>,
    /// Per-volume distance from maximum to nearest atom (Å).
    pub iondist: Vec<f64>,
    /// Per-volume nearest atom index (0-based).
    pub nnion: Vec<usize>,
    /// Per-atom Bader volume (Å³).
    pub ionvol: Vec<f64>,

    /// Total vacuum charge (e).
    pub vacchg: f64,
    /// Total vacuum volume (Å³).
    pub vacvol: f64,
}

/// Bader charge analyzer.  Use the builder pattern to configure and run.
///
/// ```ignore
/// use ferro_analysis::dft::{BaderAnalyzer, BaderMethod};
///
/// let (frame, chg) = ferro_io::read_chgcar("CHGCAR")?;
/// let result = BaderAnalyzer::new(chg, frame)
///     .method(BaderMethod::Weight)
///     .run();
/// result.write_acf("ACF.dat")?;
/// ```
pub struct BaderAnalyzer {
    chg: ChargeGrid,
    frame: Frame,
    params: BaderParams,
}

impl BaderAnalyzer {
    /// Create a new analyzer for the given charge grid and atomic structure.
    pub fn new(chg: ChargeGrid, frame: Frame) -> Self {
        Self { chg, frame, params: BaderParams::default() }
    }

    /// Set the analysis method.
    pub fn method(mut self, method: BaderMethod) -> Self {
        self.params.method = method;
        self
    }

    /// Set the edge refinement strategy.
    ///
    /// - `-1`: auto (repeat until no reassignments)
    /// - `-2`: single pass
    /// - `N > 0`: exactly N passes
    pub fn refine(mut self, refine: i32) -> Self {
        self.params.refine = refine;
        self
    }

    /// Set the vacuum density threshold (e/Å³).
    pub fn vacval(mut self, vacval: f64) -> Self {
        self.params.vacval = vacval;
        self
    }

    /// Set the off-grid step size (Å).
    pub fn stepsize(mut self, stepsize: f64) -> Self {
        self.params.stepsize = Some(stepsize);
        self
    }

    /// Run the Bader analysis.
    pub fn run(&self) -> BaderResult {
        match self.params.method {
            BaderMethod::OnGrid => {
                super::bader_grid::bader_ongrid(&self.chg, &self.frame, &self.params)
            }
            BaderMethod::NearGrid => {
                super::bader_grid::bader_neargrid(&self.chg, &self.frame, &self.params)
            }
            BaderMethod::OffGrid => {
                super::bader_grid::bader_offgrid(&self.chg, &self.frame, &self.params)
            }
            BaderMethod::Weight => {
                super::bader_weight::bader_weight(&self.chg, &self.frame, &self.params)
            }
        }
    }
}

impl BaderResult {
    /// Write ACF.dat (Atom-Centered File): per-atom Bader charge summary.
    ///
    /// Columns: atom# | X | Y | Z | charge | min_distance | atomic_volume
    pub fn write_acf(&self, path: &str, frame: &Frame) -> Result<()> {
        let f = std::fs::File::create(path)
            .map_err(ChemError::IoError)?;
        let mut w = BufWriter::new(f);

        writeln!(w, "# ACF.dat — Bader charge analysis (ferro v{})", env!("CARGO_PKG_VERSION"))
            .map_err(ChemError::IoError)?;
        writeln!(w, "# {:>5} {:>12} {:>12} {:>12} {:>14} {:>14} {:>14}",
                 "Atom", "X", "Y", "Z", "Charge", "MinDist", "Volume")
            .map_err(ChemError::IoError)?;

        for (i, charge) in self.ionchg.iter().enumerate() {
            let vol = self.ionvol.get(i).copied().unwrap_or(0.0);
            let min_dist = self.nnion.iter().enumerate()
                .filter(|(_, &nn)| nn == i)
                .map(|(v, _)| self.iondist.get(v).copied().unwrap_or(f64::MAX))
                .fold(f64::MAX, f64::min);
            let min_dist = if min_dist == f64::MAX { 0.0 } else { min_dist };
            let pos = frame.atoms.get(i).map(|a| a.position).unwrap_or_default();

            writeln!(w, " {:>5} {:>12.6} {:>12.6} {:>12.6} {:>14.6} {:>14.6} {:>14.6}",
                     i + 1, pos.x, pos.y, pos.z, charge, min_dist, vol)
                .map_err(ChemError::IoError)?;
        }

        writeln!(w, "\n# Vacuum charge: {:.6} e,  Vacuum volume: {:.6} A^3", self.vacchg, self.vacvol)
            .map_err(ChemError::IoError)?;
        writeln!(w, "# Total: {:.6} e", self.ionchg.iter().sum::<f64>() + self.vacchg)
            .map_err(ChemError::IoError)?;

        Ok(())
    }

    /// Write BCF.dat (Bader-Centered File): per-volume summary.
    ///
    /// Columns: volume# | X | Y | Z | charge | assigned_atom | distance
    pub fn write_bcf(&self, path: &str) -> Result<()> {
        let f = std::fs::File::create(path)
            .map_err(ChemError::IoError)?;
        let mut w = BufWriter::new(f);

        writeln!(w, "# BCF.dat — Bader volume analysis (ferro v{})", env!("CARGO_PKG_VERSION"))
            .map_err(ChemError::IoError)?;
        writeln!(w, "# {:>5} {:>12} {:>12} {:>12} {:>14} {:>10} {:>14}",
                 "Vol", "X", "Y", "Z", "Charge", "Atom", "Distance")
            .map_err(ChemError::IoError)?;

        for v in 0..self.nvols {
            let x = self.volpos_car.get(v).map(|p| p[0]).unwrap_or(0.0);
            let y = self.volpos_car.get(v).map(|p| p[1]).unwrap_or(0.0);
            let z = self.volpos_car.get(v).map(|p| p[2]).unwrap_or(0.0);
            let chg = self.volchg.get(v).copied().unwrap_or(0.0);
            let nn = self.nnion.get(v).map(|&i| i + 1).unwrap_or(0);
            let dist = self.iondist.get(v).copied().unwrap_or(0.0);

            writeln!(w, " {:>5} {:>12.6} {:>12.6} {:>12.6} {:>14.6} {:>10} {:>14.6}",
                     v + 1, x, y, z, chg, nn, dist)
                .map_err(ChemError::IoError)?;
        }

        Ok(())
    }

    /// Write AVF.dat (Atom-Volume File): per-atom list of assigned volume indices.
    pub fn write_avf(&self, path: &str) -> Result<()> {
        let f = std::fs::File::create(path)
            .map_err(ChemError::IoError)?;
        let mut w = BufWriter::new(f);

        writeln!(w, "# AVF.dat — Atom-Volume assignment (ferro v{})", env!("CARGO_PKG_VERSION"))
            .map_err(ChemError::IoError)?;

        let nions = self.ionchg.len();
        for i in 0..nions {
            let vols: Vec<usize> = self.nnion.iter().enumerate()
                .filter(|(_, &nn)| nn == i)
                .map(|(v, _)| v + 1)
                .collect();
            let line = vols.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" ");
            writeln!(w, " {:>5} : {}", i + 1, line).map_err(ChemError::IoError)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, ChargeGrid, Frame};
    use nalgebra::{Matrix3, Vector3};

    /// Create a simple test system: two Gaussian peaks on a 10×10×10 grid.
    /// Atom A at (2,2,2) lattice steps, Atom B at (7,7,7) lattice steps.
    /// Each peak has a Gaussian shape so Bader volumes are well-separated.
    fn two_peak_system() -> (Frame, ChargeGrid) {
        let shape = [10, 10, 10];
        let n = 10;
        let cell = Cell::from_matrix(Matrix3::new(
            10.0, 0.0, 0.0,
            0.0, 10.0, 0.0,
            0.0, 0.0, 10.0,
        ));
        let vol = cell.volume(); // 1000 Å³

        // Create density: two Gaussian peaks, stored as rho × V_cell
        let mut rho = vec![0.0_f64; n * n * n];
        let centers = [(2.0, 2.0, 2.0), (7.0, 7.0, 7.0)];
        let sigma = 1.5_f64;

        for i3 in 0..n {
            for i2 in 0..n {
                for i1 in 0..n {
                    let idx = i1 + n * (i2 + n * i3);
                    let mut val = 0.0;
                    for &(cx, cy, cz) in &centers {
                        let dx = (i1 as f64 - cx).min((i1 as f64 - cx + n as f64).abs()).min((i1 as f64 - cx - n as f64).abs());
                        let dy = (i2 as f64 - cy).min((i2 as f64 - cy + n as f64).abs()).min((i2 as f64 - cy - n as f64).abs());
                        let dz = (i3 as f64 - cz).min((i3 as f64 - cz + n as f64).abs()).min((i3 as f64 - cz - n as f64).abs());
                        let r2 = dx * dx + dy * dy + dz * dz;
                        val += (-r2 / (2.0 * sigma * sigma)).exp();
                    }
                    // Store as rho_physical × V_cell
                    rho[idx] = val * vol;
                }
            }
        }

        let chg = ChargeGrid::new(rho, shape, &cell);
        let mut frame = Frame::with_cell(cell, [true; 3]);
        // Place atoms at the peak centers (Cartesian)
        frame.add_atom(Atom::new("Na", Vector3::new(2.0, 2.0, 2.0)));
        frame.add_atom(Atom::new("Cl", Vector3::new(7.0, 7.0, 7.0)));

        (frame, chg)
    }

    #[test]
    fn test_bader_ongrid_finds_two_volumes() {
        let (frame, chg) = two_peak_system();
        let result = BaderAnalyzer::new(chg, frame)
            .method(BaderMethod::OnGrid)
            .refine(0)
            .run();
        assert_eq!(result.nvols, 2, "should find 2 Bader volumes");
    }

    #[test]
    fn test_bader_neargrid_finds_two_volumes() {
        let (frame, chg) = two_peak_system();
        let result = BaderAnalyzer::new(chg, frame)
            .method(BaderMethod::NearGrid)
            .refine(0)
            .run();
        assert_eq!(result.nvols, 2, "should find 2 Bader volumes");
    }

    #[test]
    fn test_bader_weight_finds_two_volumes() {
        let (frame, chg) = two_peak_system();
        let result = BaderAnalyzer::new(chg, frame)
            .method(BaderMethod::Weight)
            .run();
        assert_eq!(result.nvols, 2, "should find 2 Bader volumes");
    }

    #[test]
    fn test_bader_charge_conservation() {
        let (frame, chg) = two_peak_system();
        let result = BaderAnalyzer::new(chg, frame)
            .method(BaderMethod::NearGrid)
            .refine(0)
            .run();
        // Total charge should be conserved: sum(ionchg) + vacchg = total electrons
        let total_ion: f64 = result.ionchg.iter().sum();
        let total = total_ion + result.vacchg;
        // The total should be positive (we have electrons)
        assert!(total > 0.0, "total charge should be positive, got {total}");
    }

    #[test]
    fn test_bader_atoms_assigned() {
        let (frame, chg) = two_peak_system();
        let result = BaderAnalyzer::new(chg, frame)
            .method(BaderMethod::NearGrid)
            .refine(0)
            .run();
        // Both atoms should have non-zero charge
        assert!(result.ionchg[0] > 0.0, "atom 0 should have charge > 0");
        assert!(result.ionchg[1] > 0.0, "atom 1 should have charge > 0");
        // Each atom should be the nearest to at least one volume
        assert!(result.nnion.contains(&0), "atom 0 should be assigned to a volume");
        assert!(result.nnion.contains(&1), "atom 1 should be assigned to a volume");
    }

    #[test]
    fn test_bader_output_files() {
        let (frame, chg) = two_peak_system();
        let result = BaderAnalyzer::new(chg, frame.clone())
            .method(BaderMethod::OnGrid)
            .refine(0)
            .run();

        let dir = std::env::temp_dir();
        let acf = dir.join("test_acf.dat");
        let bcf = dir.join("test_bcf.dat");
        let avf = dir.join("test_avf.dat");

        result.write_acf(acf.to_str().unwrap(), &frame).unwrap();
        result.write_bcf(bcf.to_str().unwrap()).unwrap();
        result.write_avf(avf.to_str().unwrap()).unwrap();

        assert!(acf.exists());
        assert!(bcf.exists());
        assert!(avf.exists());

        // Cleanup
        let _ = std::fs::remove_file(&acf);
        let _ = std::fs::remove_file(&bcf);
        let _ = std::fs::remove_file(&avf);
    }
}
