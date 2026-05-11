//! CHGCAR format reader (VASP charge density file).
//!
//! Returns `(Frame, ChargeGrid)` — the structural data and the volumetric charge density.
//! Density values are stored as-is (`ρ × V_cell`), not normalized, per VASP convention.

use ferro_core::{Atom, Cell, ChargeGrid, Frame};
use nalgebra::{Matrix3, Vector3};
use anyhow::{ensure, Context, Result};

/// Read a VASP CHGCAR file, returning the structural frame and charge density grid.
pub fn read_chgcar(path: &str) -> Result<(Frame, ChargeGrid)> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path}"))?;
    parse_chgcar(&content).with_context(|| format!("parsing {path}"))
}

fn parse_chgcar(content: &str) -> Result<(Frame, ChargeGrid)> {
    let mut lines = content.lines();
    let mut next = |what: &str| -> Result<&str> {
        lines.next().with_context(|| format!("unexpected EOF before {what}"))
    };

    // ── POSCAR-style header ────────────────────────────────────────────────
    let _comment = next("comment")?.trim().to_string();

    let scale: f64 = next("scaling factor")?
        .trim().parse().context("invalid scaling factor")?;
    ensure!(scale > 0.0, "negative scaling factor is not supported");

    let mut m = [[0.0_f64; 3]; 3];
    for (i, row) in m.iter_mut().enumerate() {
        let v = floats(next(&format!("lattice vector {i}"))?, 3)
            .with_context(|| format!("invalid lattice vector {i}"))?;
        *row = [v[0] * scale, v[1] * scale, v[2] * scale];
    }
    let cell = Cell::from_matrix(Matrix3::new(
        m[0][0], m[0][1], m[0][2],
        m[1][0], m[1][1], m[1][2],
        m[2][0], m[2][1], m[2][2],
    ));

    // VASP4 vs VASP5 detection
    let line5 = next("element/count line")?.trim();
    let (elements, counts): (Vec<String>, Vec<usize>) =
        if line5.split_whitespace().all(|t| t.parse::<u64>().is_ok()) {
            let c: Vec<usize> = line5.split_whitespace().map(|s| s.parse().unwrap()).collect();
            let e = (1..=c.len()).map(|i| format!("X{i}")).collect();
            (e, c)
        } else {
            let e: Vec<String> = line5.split_whitespace().map(|s| s.to_string()).collect();
            let c: Vec<usize> = next("atom counts")?
                .split_whitespace()
                .map(|s| s.parse().context("invalid count"))
                .collect::<Result<_>>()?;
            ensure!(e.len() == c.len(), "element/count length mismatch");
            (e, c)
        };

    // Optional "Selective dynamics"
    let mut coord_type = next("coordinate type")?;
    if coord_type.trim().to_lowercase().starts_with('s') {
        coord_type = next("coordinate type after Selective dynamics")?;
    }
    let is_direct = coord_type.trim().to_lowercase().starts_with('d');

    // Atom positions
    let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
    let _total: usize = counts.iter().sum();
    for (elem, &count) in elements.iter().zip(counts.iter()) {
        for _ in 0..count {
            let v = floats(next("coordinate")?, 3).context("invalid coordinate")?;
            let pos = if is_direct {
                cell.fractional_to_cartesian(Vector3::new(v[0], v[1], v[2]))
            } else {
                Vector3::new(v[0] * scale, v[1] * scale, v[2] * scale)
            };
            frame.add_atom(Atom::new(elem.as_str(), pos));
        }
    }

    // ── Charge density grid ────────────────────────────────────────────────

    // Skip blank lines between atom coordinates and grid dimensions
    let grid_line = loop {
        let line = lines.next()
            .context("unexpected EOF before grid dimensions")?;
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            break trimmed;
        }
    };

    let dims: Vec<usize> = grid_line.split_whitespace()
        .map(|s| s.parse().context("invalid grid dimension"))
        .collect::<Result<_>>()?;
    ensure!(dims.len() >= 3, "grid dimension line needs at least 3 integers");
    let shape = [dims[0], dims[1], dims[2]];
    let nrho = shape[0] * shape[1] * shape[2];

    // Read all remaining density values (whitespace-separated, x fastest)
    let rho: Vec<f64> = lines
        .flat_map(|l| l.split_whitespace().map(|s| s.parse::<f64>().unwrap_or(0.0)))
        .collect();
    ensure!(
        rho.len() >= nrho,
        "charge density data too short: got {}, expected {}",
        rho.len(), nrho
    );
    let rho = rho[..nrho].to_vec();

    let chg = ChargeGrid::new(rho, shape, &cell);
    Ok((frame, chg))
}

fn floats(line: &str, min: usize) -> Result<Vec<f64>> {
    let v: Vec<f64> = line.split_whitespace()
        .map_while(|s| s.parse::<f64>().ok())
        .collect();
    ensure!(v.len() >= min, "expected ≥{min} floats, got {}", v.len());
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Minimal CHGCAR: 2×2×2 grid, 1 atom (NaCl-like, single Na at origin)
    const MINIMAL_CHGCAR: &str = "\
NaCl minimal
  1.00000000
     5.64000000   0.00000000   0.00000000
     0.00000000   5.64000000   0.00000000
     0.00000000   0.00000000   5.64000000
  Na
  1
Direct
  0.000000  0.000000  0.000000

  2  2  2
 1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0
";

    // CHGCAR with VASP5 element line and selective dynamics
    const VASP5_CHGCAR: &str = "\
VASP5 test
  1.00000000
     4.00000000   0.00000000   0.00000000
     0.00000000   4.00000000   0.00000000
     0.00000000   0.00000000   4.00000000
  Fe  O
  1  1
Selective dynamics
Direct
  0.000000  0.000000  0.000000   T  T  T
  0.500000  0.500000  0.500000   T  T  T

  2  2  2
 10.0  20.0  30.0  40.0  50.0  60.0  70.0  80.0
";

    fn tmp(name: &str, content: &str) -> String {
        let p = std::env::temp_dir().join(name);
        std::fs::write(&p, content).unwrap();
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn test_read_frame_atoms() {
        let (frame, _) = parse_chgcar(MINIMAL_CHGCAR).unwrap();
        assert_eq!(frame.n_atoms(), 1);
        assert_eq!(frame.atom(0).element, "Na");
    }

    #[test]
    fn test_read_cell() {
        let (frame, _) = parse_chgcar(MINIMAL_CHGCAR).unwrap();
        let cell = frame.cell.as_ref().unwrap();
        let [a, b, c] = cell.lengths();
        assert!((a - 5.64).abs() < 1e-6);
        assert!((b - 5.64).abs() < 1e-6);
        assert!((c - 5.64).abs() < 1e-6);
    }

    #[test]
    fn test_read_grid_shape() {
        let (_, chg) = parse_chgcar(MINIMAL_CHGCAR).unwrap();
        assert_eq!(chg.shape, [2, 2, 2]);
        assert_eq!(chg.nrho, 8);
    }

    #[test]
    fn test_read_density_values() {
        let (_, chg) = parse_chgcar(MINIMAL_CHGCAR).unwrap();
        // x fastest: rho[0]=1, rho[1]=2, rho[2]=3, rho[3]=4, ...
        assert!((chg.rho[0] - 1.0).abs() < 1e-10);
        assert!((chg.rho[1] - 2.0).abs() < 1e-10);
        assert!((chg.rho[4] - 5.0).abs() < 1e-10);
        assert!((chg.rho[7] - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_read_vasp5_elements() {
        let (frame, _) = parse_chgcar(VASP5_CHGCAR).unwrap();
        assert_eq!(frame.n_atoms(), 2);
        assert_eq!(frame.atom(0).element, "Fe");
        assert_eq!(frame.atom(1).element, "O");
    }

    #[test]
    fn test_read_vasp5_positions() {
        let (frame, chg) = parse_chgcar(VASP5_CHGCAR).unwrap();
        // Fe at origin
        let p0 = frame.atom(0).position;
        assert!(p0.norm() < 1e-6);
        // O at (0.5, 0.5, 0.5) → (2.0, 2.0, 2.0) Å
        let p1 = frame.atom(1).position;
        assert!((p1.x - 2.0).abs() < 1e-6);
        assert!((p1.y - 2.0).abs() < 1e-6);
        assert!((p1.z - 2.0).abs() < 1e-6);
        // Grid
        assert_eq!(chg.shape, [2, 2, 2]);
    }

    #[test]
    fn test_read_from_file() {
        let path = tmp("test_chgcar.CHGCAR", MINIMAL_CHGCAR);
        let (frame, chg) = read_chgcar(&path).unwrap();
        assert_eq!(frame.n_atoms(), 1);
        assert_eq!(chg.shape, [2, 2, 2]);
    }
}
