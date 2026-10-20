//! Vacuum layer addition.

use ferro_core::cell::Cell;
use ferro_core::error::{ChemError, Result};
use ferro_core::frame::Frame;

/// Add a vacuum layer along the specified axis.
///
/// - `axis`: "x", "y", or "z" (maps to cell matrix row 0, 1, 2).
/// - `thickness`: vacuum thickness in Å (must be > 0), measured **perpendicular** to the
///   two other lattice vectors — the real gap between periodic images, also for a tilted axis.
///   Atoms are not moved or centred: a slab crossing the cell boundary is split by the gap.
///
/// Returns a new `Frame` with the cell expanded along the given axis.
/// Atomic positions, `pbc`, `bonds`, `charge`, and `multiplicity` are preserved.
/// `energy`, `forces`, `stress`, `velocities` are set to `None`.
///
/// # Errors
/// - `ValidationError` if the frame has no cell.
/// - `ValidationError` if `axis` is not "x"/"y"/"z".
/// - `ValidationError` if `thickness` is not a finite number > 0.
/// - `ValidationError` if the cell is singular (e.g. c = 0, the ASE convention for 2D materials)
///   or not finite (`Cell::check`).
pub fn add_vacuum(frame: &Frame, axis: &str, thickness: f64) -> Result<Frame> {
    let axis_idx = match axis {
        "x" => 0,
        "y" => 1,
        "z" => 2,
        _ => {
            return Err(ChemError::ValidationError(format!(
                "invalid axis '{}', expected \"x\", \"y\", or \"z\"",
                axis
            )));
        }
    };

    // 写成「合法则放行」，NaN 落到报错分支（审查 E-9：以前 NaN 厚度静默得到 NaN 晶胞）
    if !(thickness > 0.0 && thickness.is_finite()) {
        return Err(ChemError::ValidationError(format!(
            "vacuum thickness must be a finite number > 0, got {}",
            thickness
        )));
    }

    let cell = frame.cell.as_ref().ok_or_else(|| {
        ChemError::ValidationError("add_vacuum requires a periodic frame with a cell".into())
    })?;
    // A singular cell (ASE's c = 0 for a 2D material) has spacing 0/0; filling c in
    // from the atoms' extent is a feature of its own
    cell.check("add_vacuum: the cell (a 2D material needs a full 3D cell first)")?;

    // 沿该晶格矢量拉长，使它对面的面间距（两个相邻晶面的垂直距离）增加 thickness。
    // 以前是让矢量长度增加 thickness，c 轴倾斜时实际真空间隙只有 thickness·cosθ。
    // 面间距 d_i = V / |a_j × a_k| 与 |a_i| 成正比，所以按 d_i 求缩放比即可（同 ASE 的
    // vacuum / cosθ）
    let m = &cell.matrix;
    let (j, k) = ((axis_idx + 1) % 3, (axis_idx + 2) % 3);
    let face = m.row(j).cross(&m.row(k)).norm();
    let spacing = m.determinant().abs() / face;
    let scale = (spacing + thickness) / spacing;

    let mut new_m = *m;
    for j in 0..3 {
        new_m[(axis_idx, j)] *= scale;
    }

    // Same atoms and bonds; the results of the old box do not carry over
    let mut out = frame.clone();
    // Only one lattice vector stretches; the box origin stays (review A-12)
    out.cell = Some(Cell { matrix: new_m, origin: cell.origin });
    out.set_forces(None)?;
    out.set_velocities(None)?;
    out.energy = None;
    out.stress = None;
    out.temperature = None;
    out.step = None;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::atom::Atom;
    use nalgebra::Vector3;

    fn cubic_frame(a: f64) -> Frame {
        let cell = Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(1.0, 2.0, 3.0)));
        frame
    }

    #[test]
    fn test_add_vacuum_z() {
        let frame = cubic_frame(10.0);
        let out = add_vacuum(&frame, "z", 5.0).unwrap();
        let [la, lb, lc] = out.cell.as_ref().unwrap().lengths();
        assert!((la - 10.0).abs() < 1e-10);
        assert!((lb - 10.0).abs() < 1e-10);
        assert!((lc - 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_add_vacuum_x() {
        let frame = cubic_frame(10.0);
        let out = add_vacuum(&frame, "x", 8.0).unwrap();
        let [la, lb, lc] = out.cell.as_ref().unwrap().lengths();
        assert!((la - 18.0).abs() < 1e-10);
        assert!((lb - 10.0).abs() < 1e-10);
        assert!((lc - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_add_vacuum_y() {
        let frame = cubic_frame(10.0);
        let out = add_vacuum(&frame, "y", 3.0).unwrap();
        let [la, lb, lc] = out.cell.as_ref().unwrap().lengths();
        assert!((la - 10.0).abs() < 1e-10);
        assert!((lb - 13.0).abs() < 1e-10);
        assert!((lc - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_tilted_axis_gets_full_perpendicular_gap() {
        // β = 60°：c 与 ab 面夹 60°，面间距 = 10·sin60°。加 5 Å 后面间距应恰好 +5
        let cell = Cell::from_lengths_angles(8.0, 8.0, 10.0, 90.0, 60.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("O", Vector3::zeros()));
        let spacing = |f: &Frame| {
            let m = f.cell.as_ref().unwrap().matrix;
            m.determinant().abs() / m.row(0).cross(&m.row(1)).norm()
        };
        let out = add_vacuum(&frame, "z", 5.0).unwrap();
        assert!((spacing(&out) - spacing(&frame) - 5.0).abs() < 1e-10,
            "垂直间隙应增加 5 Å，实际 {}", spacing(&out) - spacing(&frame));
        let [_, _, lc] = out.cell.as_ref().unwrap().lengths();
        assert!((lc - (10.0 + 5.0 / 60f64.to_radians().sin())).abs() < 1e-10, "c 应拉长 5/sin60°");
    }

    // 审查 A-12：LAMMPS 读进来的 lo 不能被加真空清零
    #[test]
    fn test_box_origin_preserved() {
        let mut frame = cubic_frame(10.0);
        let o = Vector3::new(1.5, -2.0, 3.0);
        frame.cell.as_mut().unwrap().origin = o;
        let out = add_vacuum(&frame, "z", 5.0).unwrap();
        assert_eq!(out.cell.as_ref().unwrap().origin, o, "加真空应保留盒子原点");
    }

    #[test]
    fn test_positions_unchanged() {
        let frame = cubic_frame(10.0);
        let out = add_vacuum(&frame, "z", 5.0).unwrap();
        let pos = out.positions()[0];
        assert!((pos.x - 1.0).abs() < 1e-10);
        assert!((pos.y - 2.0).abs() < 1e-10);
        assert!((pos.z - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_pbc_unchanged() {
        let frame = cubic_frame(10.0);
        let out = add_vacuum(&frame, "z", 5.0).unwrap();
        assert_eq!(out.pbc, [true; 3]);
    }

    #[test]
    fn test_atom_charges_preserved() {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("O", Vector3::zeros()));
        frame.set_charge(0, Some(-0.8));
        let out = add_vacuum(&frame, "z", 5.0).unwrap();
        assert_eq!(out.atom(0).charge, Some(-0.8));
    }

    #[test]
    fn test_results_cleared() {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("O", Vector3::zeros()));
        frame.energy = Some(-100.0);
        frame.set_forces(Some(vec![Vector3::zeros()])).unwrap();
        let out = add_vacuum(&frame, "z", 5.0).unwrap();
        assert!(out.energy.is_none());
        assert!(out.forces().is_none());
    }

    #[test]
    fn test_triclinic_cell() {
        let cell = Cell::from_lengths_angles(5.0, 6.0, 7.0, 80.0, 85.0, 95.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Si", Vector3::new(1.0, 1.0, 1.0)));
        let out = add_vacuum(&frame, "z", 10.0).unwrap();
        let lengths = out.cell.as_ref().unwrap().lengths();
        // c 拉长到垂直间隙 +10 Å：17.2099…，期望值取自 ASE 3.26 `add_vacuum`
        //（旧实现让 |c| 本身 +10 = 17，倾斜时实际间隙不足 10）
        assert!((lengths[2] - 17.20991466012953).abs() < 1e-9, "c = {}", lengths[2]);
        // x, y 不变
        assert!((lengths[0] - 5.0).abs() < 1e-10);
        assert!((lengths[1] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_error_no_cell() {
        let mut frame = Frame::new();
        frame.add_atom(Atom::new("Fe", Vector3::zeros()));
        assert!(add_vacuum(&frame, "z", 5.0).is_err());
    }

    #[test]
    fn test_error_invalid_axis() {
        let frame = cubic_frame(10.0);
        assert!(add_vacuum(&frame, "w", 5.0).is_err());
    }

    #[test]
    fn test_error_zero_thickness() {
        let frame = cubic_frame(10.0);
        assert!(add_vacuum(&frame, "z", 0.0).is_err());
    }

    // 审查 E-9：NaN / inf 厚度与奇异胞（ASE 二维材料的 c=0 约定）以前静默返回 NaN 晶胞
    #[test]
    fn test_error_non_finite_thickness() {
        let frame = cubic_frame(10.0);
        for t in [f64::NAN, f64::INFINITY] {
            assert!(add_vacuum(&frame, "z", t).is_err(), "厚度 {t} 应报错");
        }
    }

    #[test]
    fn test_error_singular_cell() {
        let mut frame = cubic_frame(10.0);
        frame.cell.as_mut().unwrap().matrix[(2, 2)] = 0.0;
        assert!(add_vacuum(&frame, "z", 5.0).is_err());
    }

    #[test]
    fn test_error_negative_thickness() {
        let frame = cubic_frame(10.0);
        assert!(add_vacuum(&frame, "z", -1.0).is_err());
    }

    // Review A-13: a NaN cell passed the determinant test (NaN < 1e-12 is false)
    #[test]
    fn test_error_non_finite_cell() {
        let mut frame = cubic_frame(10.0);
        frame.cell.as_mut().unwrap().matrix[(0, 1)] = f64::NAN;
        let e = add_vacuum(&frame, "z", 5.0).unwrap_err().to_string();
        assert!(e.contains("not finite"), "{e}");
    }
}
