//! CHGCAR format reader (VASP charge density file).
//!
//! Returns `(Frame, ChargeGrid)` — the structural data and the volumetric charge density.
//! Density values are stored as-is (`ρ × V_cell`), not normalized, per VASP convention.

use std::path::Path;
use ferro_core::{Atom, Cell, ChargeGrid, Frame};
use nalgebra::{Matrix3, Vector3};
use anyhow::{ensure, Context, Result};
use super::util::floats;
use super::vasp::potcar_symbol;

/// Read a VASP CHGCAR file, returning the structural frame and charge density grid.
pub fn read_chgcar(path: &Path) -> Result<(Frame, ChargeGrid)> {
    let path_ = path.display();
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path_}"))?;
    parse_chgcar(&content).with_context(|| format!("parsing {path_}"))
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
    // Parse line5 once as usize counts: Some ⇒ VASP4, None ⇒ VASP5.
    // (Parsing directly as usize avoids the prior u64-check / usize-unwrap
    // mismatch that could panic on 32-bit targets for large values.)
    let line5_counts: Option<Vec<usize>> =
        line5.split_whitespace().map(|t| t.parse::<usize>().ok()).collect();
    let (elements, counts): (Vec<String>, Vec<usize>) =
        if let Some(c) = line5_counts {
            let e = (1..=c.len()).map(|i| format!("X{i}")).collect();
            (e, c)
        } else {
            // 元素行规则与 POSCAR 相同（含 6.4.2 的 `标签/哈希`），见 vasp.rs
            let e: Vec<String> = line5.split_whitespace().map(potcar_symbol).collect();
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
    // VASP 的规则：首字符 C/c/K/k 为 Cartesian，其余一律 Direct（同 vasp.rs）
    let is_direct = !coord_type.trim_start().starts_with(['C', 'c', 'K', 'k']);

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

    // 只取前 nrho 个数（x 最快）。之后是 augmentation 段（含文字）与自旋密度块，不属于
    // 这张网格；但前 nrho 个里的坏值要报错 —— 以前 `unwrap_or(0.0)` 让它冒充「这里密度为 0」
    let mut rho: Vec<f64> = Vec::with_capacity(nrho);
    'grid: for (k, l) in lines.enumerate() {
        for tok in l.split_whitespace() {
            if rho.len() == nrho { break 'grid; }
            let v = parse_fortran_float(tok).with_context(|| format!(
                "invalid charge density value {tok:?} on line {} after the grid dimensions \
                 (value {} of {nrho})", k + 1, rho.len() + 1
            ))?;
            rho.push(v);
        }
    }
    ensure!(
        rho.len() == nrho,
        "charge density data too short: got {}, expected {}",
        rho.len(), nrho
    );

    let chg = ChargeGrid::new(rho, shape, &cell);
    Ok((frame, chg))
}


/// 解析一个数，兼容 Fortran 在三位指数时省掉 `E` 的写法（`0.1234-100` = 0.1234E-100）。
/// VASP 的 CHGCAR 在真空区的极小密度上会写出这种数。
fn parse_fortran_float(tok: &str) -> Result<f64> {
    if let Ok(v) = tok.parse::<f64>() {
        return Ok(v);
    }
    // 第 0 位之后最后一个正负号前补 E；前面不能已经有 E/e（那是真坏值）
    let split = tok.rfind(['+', '-']).filter(|&i| i > 0 && !tok[..i].contains(['E', 'e']));
    split
        .and_then(|i| format!("{}E{}", &tok[..i], &tok[i..]).parse::<f64>().ok())
        .ok_or_else(|| anyhow::anyhow!("not a number"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bad_density_value_is_an_error_but_trailing_sections_are_not() {
        // 网格内的坏值：以前静默变 0
        let bad = MINIMAL_CHGCAR.replace(" 5.0  6.0", " 5.0  6.O");
        let msg = format!("{:#}", parse_chgcar(&bad).unwrap_err());
        assert!(msg.contains("\"6.O\"") && msg.contains("value 6 of 8"), "应点名坏值与位置，实际：{msg}");
        // Fortran 三位指数省掉 E 的写法照常读
        let fortran = MINIMAL_CHGCAR.replace(" 5.0  6.0", " 0.5-100  6.0");
        let (_, chg) = parse_chgcar(&fortran).unwrap();
        assert!(chg.rho.iter().any(|&v| v > 0.0 && v < 1e-99), "0.5-100 应读成 5e-101");
        assert_eq!(parse_fortran_float("-0.25-101").unwrap(), -0.25e-101);
        assert!(parse_fortran_float("1.0E-3-2").is_err());
        // 网格之后的 augmentation 段是文字，不该报错
        let tail = format!("{MINIMAL_CHGCAR}augmentation occupancies   1  4\n 0.1 0.2 0.3 0.4\n");
        let (_, chg) = parse_chgcar(&tail).unwrap();
        assert_eq!(chg.shape, [2, 2, 2]);
    }

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

    use crate::testutil::write_tmp as tmp;

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
    fn test_element_line_and_coordinate_mode() {
        // 6.4.2 的 `标签/哈希`；`Fractional` 按 VASP 规则是 Direct（以前只认 d 开头）
        let text = VASP5_CHGCAR
            .replace("  Fe  O\n", "  Fe_pv/6a2f546d  O/\n")
            .replace("Direct\n", "Fractional\n");
        let (frame, _) = parse_chgcar(&text).unwrap();
        assert_eq!(frame.atom(0).element, "Fe");
        assert_eq!(frame.atom(1).element, "O");
        assert!((frame.atom(1).position.x - 2.0).abs() < 1e-6, "Fractional 应按 Direct 读");
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
