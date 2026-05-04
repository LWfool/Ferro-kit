//! Merge two frames along an axis with a vacuum gap at the interface.

use nalgebra::Vector3;

use ferro_core::cell::Cell;
use ferro_core::error::{ChemError, Result};
use ferro_core::frame::Frame;

/// Merge two frames along the specified axis with a vacuum gap at the interface.
///
/// - `axis`: "x", "y", or "z" (the direction along which structures are stacked).
/// - `gap`: vacuum gap thickness in Å at the interface (must be >= 0).
///
/// The merged cell along the join axis has length `L_a + gap + L_b`.
/// On the other two axes, the length is `max(L_a, L_b)`; the shorter structure
/// is centered on those axes.
///
/// `bonds` are not merged (indices shift). `energy`, `forces`, `stress`,
/// `velocities` are set to `None`. `charge` = sum of both. `multiplicity` = 1.
///
/// # Errors
/// - `ValidationError` if either frame has no cell.
/// - `ValidationError` if `axis` is not "x"/"y"/"z".
/// - `ValidationError` if `gap < 0`.
pub fn merge_frames(frame_a: &Frame, frame_b: &Frame, axis: &str, gap: f64) -> Result<Frame> {
    let join = match axis {
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

    if gap < 0.0 {
        return Err(ChemError::ValidationError(format!(
            "gap must be >= 0, got {}",
            gap
        )));
    }

    let cell_a = frame_a.cell.as_ref().ok_or_else(|| {
        ChemError::ValidationError("merge_frames: frame_a has no cell".into())
    })?;
    let cell_b = frame_b.cell.as_ref().ok_or_else(|| {
        ChemError::ValidationError("merge_frames: frame_b has no cell".into())
    })?;

    let lens_a = cell_a.lengths();
    let lens_b = cell_b.lengths();

    // ── 合并后各轴长度 ─────────────────────────────────────────────────────────
    let mut new_lens = [0.0; 3];
    for i in 0..3 {
        if i == join {
            new_lens[i] = lens_a[i] + gap + lens_b[i];
        } else {
            new_lens[i] = lens_a[i].max(lens_b[i]);
        }
    }

    // ── 构建新 cell 矩阵（缩放原 a 的行向量）──────────────────────────────────
    let m_a = &cell_a.matrix;
    let mut new_m = *m_a;
    for i in 0..3 {
        let scale = new_lens[i] / lens_a[i];
        for j in 0..3 {
            new_m[(i, j)] *= scale;
        }
    }

    // ── frame_b 平移 ──────────────────────────────────────────────────────────
    // 拼接轴方向：平移 L_a + gap
    // 其他轴方向：居中偏移 (max_L - L_b) / 2
    let row_unit = |row: usize| -> Vector3<f64> {
        let r = cell_b.matrix.row(row);
        let v = Vector3::new(r[0], r[1], r[2]);
        v / v.norm()
    };

    let mut offset = Vector3::zeros();
    for i in 0..3 {
        if i == join {
            offset += row_unit(i) * (lens_a[i] + gap);
        } else {
            offset += row_unit(i) * ((new_lens[i] - lens_b[i]) / 2.0);
        }
    }

    let atoms_b: Vec<_> = frame_b
        .atoms
        .iter()
        .map(|atom| {
            let mut a = atom.clone();
            a.position += offset;
            a
        })
        .collect();

    // ── 组装 ──────────────────────────────────────────────────────────────────
    let mut atoms = frame_a.atoms.clone();
    atoms.extend(atoms_b);

    Ok(Frame {
        atoms,
        cell: Some(Cell::from_matrix(new_m)),
        pbc: frame_a.pbc,
        charge: frame_a.charge + frame_b.charge,
        multiplicity: 1,
        bonds: None,
        energy: None,
        forces: None,
        stress: None,
        velocities: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::atom::Atom;

    fn cubic_frame(a: f64, atoms: Vec<(&str, f64, f64, f64)>) -> Frame {
        let cell = Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        for (elem, x, y, z) in atoms {
            frame.add_atom(Atom::new(elem, Vector3::new(x, y, z)));
        }
        frame
    }

    #[test]
    fn test_merge_z_basic() {
        let a = cubic_frame(10.0, vec![("Fe", 1.0, 1.0, 1.0)]);
        let b = cubic_frame(10.0, vec![("O", 5.0, 5.0, 5.0)]);
        let out = merge_frames(&a, &b, "z", 2.0).unwrap();
        let [la, lb, lc] = out.cell.as_ref().unwrap().lengths();
        assert!((la - 10.0).abs() < 1e-10);
        assert!((lb - 10.0).abs() < 1e-10);
        // z = 10 + 2 + 10 = 22
        assert!((lc - 22.0).abs() < 1e-10);
    }

    #[test]
    fn test_merge_z_atom_positions() {
        let a = cubic_frame(10.0, vec![("Fe", 0.0, 0.0, 0.0)]);
        let b = cubic_frame(10.0, vec![("O", 0.0, 0.0, 0.0)]);
        let out = merge_frames(&a, &b, "z", 2.0).unwrap();
        // frame_a 原子位置不变
        let pos_a = out.atoms[0].position;
        assert!((pos_a.z - 0.0).abs() < 1e-10);
        // frame_b 原子 z 平移 10 + 2 = 12
        let pos_b = out.atoms[1].position;
        assert!((pos_b.z - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_merge_different_sizes_centering() {
        // a: 10x10x10, b: 6x6x6 → 非拼接轴 max=10, b 居中偏移 (10-6)/2=2
        let a = cubic_frame(10.0, vec![("Fe", 0.0, 0.0, 0.0)]);
        let b = cubic_frame(6.0, vec![("O", 0.0, 0.0, 0.0)]);
        let out = merge_frames(&a, &b, "z", 0.0).unwrap();
        let pos_b = out.atoms[1].position;
        // x, y 居中偏移 2
        assert!((pos_b.x - 2.0).abs() < 1e-10);
        assert!((pos_b.y - 2.0).abs() < 1e-10);
        // z 平移 10
        assert!((pos_b.z - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_merge_x_axis() {
        let a = cubic_frame(10.0, vec![("Fe", 0.0, 0.0, 0.0)]);
        let b = cubic_frame(10.0, vec![("O", 0.0, 0.0, 0.0)]);
        let out = merge_frames(&a, &b, "x", 3.0).unwrap();
        let [la, lb, lc] = out.cell.as_ref().unwrap().lengths();
        assert!((la - 23.0).abs() < 1e-10);
        assert!((lb - 10.0).abs() < 1e-10);
        assert!((lc - 10.0).abs() < 1e-10);
        // frame_b 原子 x 平移 10 + 3 = 13
        let pos_b = out.atoms[1].position;
        assert!((pos_b.x - 13.0).abs() < 1e-10);
    }

    #[test]
    fn test_merge_atom_count() {
        let a = cubic_frame(10.0, vec![("Fe", 0.0, 0.0, 0.0), ("Fe", 1.0, 0.0, 0.0)]);
        let b = cubic_frame(10.0, vec![("O", 0.0, 0.0, 0.0)]);
        let out = merge_frames(&a, &b, "z", 1.0).unwrap();
        assert_eq!(out.atoms.len(), 3);
    }

    #[test]
    fn test_merge_charge_sum() {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut a = Frame::with_cell(cell.clone(), [true; 3]);
        a.charge = 2;
        a.add_atom(Atom::new("Fe", Vector3::zeros()));
        let mut b = Frame::with_cell(cell, [true; 3]);
        b.charge = -1;
        b.add_atom(Atom::new("O", Vector3::zeros()));
        let out = merge_frames(&a, &b, "z", 1.0).unwrap();
        assert_eq!(out.charge, 1);
        assert_eq!(out.multiplicity, 1);
    }

    #[test]
    fn test_merge_results_cleared() {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut a = Frame::with_cell(cell.clone(), [true; 3]);
        a.add_atom(Atom::new("Fe", Vector3::zeros()));
        a.energy = Some(-100.0);
        let mut b = Frame::with_cell(cell, [true; 3]);
        b.add_atom(Atom::new("O", Vector3::zeros()));
        let out = merge_frames(&a, &b, "z", 1.0).unwrap();
        assert!(out.energy.is_none());
        assert!(out.forces.is_none());
        assert!(out.stress.is_none());
        assert!(out.velocities.is_none());
    }

    #[test]
    fn test_merge_zero_gap() {
        let a = cubic_frame(10.0, vec![("Fe", 0.0, 0.0, 0.0)]);
        let b = cubic_frame(10.0, vec![("O", 0.0, 0.0, 0.0)]);
        let out = merge_frames(&a, &b, "z", 0.0).unwrap();
        let lc = out.cell.as_ref().unwrap().lengths()[2];
        assert!((lc - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_merge_triclinic() {
        let cell_a = Cell::from_lengths_angles(5.0, 6.0, 7.0, 80.0, 85.0, 95.0).unwrap();
        let mut a = Frame::with_cell(cell_a, [true; 3]);
        a.add_atom(Atom::new("Si", Vector3::zeros()));
        let cell_b = Cell::from_lengths_angles(4.0, 5.0, 6.0, 80.0, 85.0, 95.0).unwrap();
        let mut b = Frame::with_cell(cell_b, [true; 3]);
        b.add_atom(Atom::new("O", Vector3::zeros()));
        let out = merge_frames(&a, &b, "z", 1.0).unwrap();
        let lengths = out.cell.as_ref().unwrap().lengths();
        // z: 7 + 1 + 6 = 14
        assert!((lengths[2] - 14.0).abs() < 1e-10);
        // x: max(5, 4) = 5
        assert!((lengths[0] - 5.0).abs() < 1e-10);
        // y: max(6, 5) = 6
        assert!((lengths[1] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_error_no_cell_a() {
        let mut a = Frame::new();
        a.add_atom(Atom::new("Fe", Vector3::zeros()));
        let b = cubic_frame(10.0, vec![("O", 0.0, 0.0, 0.0)]);
        assert!(merge_frames(&a, &b, "z", 1.0).is_err());
    }

    #[test]
    fn test_error_no_cell_b() {
        let a = cubic_frame(10.0, vec![("Fe", 0.0, 0.0, 0.0)]);
        let mut b = Frame::new();
        b.add_atom(Atom::new("O", Vector3::zeros()));
        assert!(merge_frames(&a, &b, "z", 1.0).is_err());
    }

    #[test]
    fn test_error_invalid_axis() {
        let a = cubic_frame(10.0, vec![("Fe", 0.0, 0.0, 0.0)]);
        let b = cubic_frame(10.0, vec![("O", 0.0, 0.0, 0.0)]);
        assert!(merge_frames(&a, &b, "w", 1.0).is_err());
    }

    #[test]
    fn test_error_negative_gap() {
        let a = cubic_frame(10.0, vec![("Fe", 0.0, 0.0, 0.0)]);
        let b = cubic_frame(10.0, vec![("O", 0.0, 0.0, 0.0)]);
        assert!(merge_frames(&a, &b, "z", -1.0).is_err());
    }
}
