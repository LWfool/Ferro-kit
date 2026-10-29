//! Merge two frames along an axis with a vacuum gap at the interface.

use nalgebra::{Matrix3, Vector3};

use ferro_core::atom::Atom;
use ferro_core::cell::Cell;
use ferro_core::error::{ChemError, Result};
use ferro_core::frame::Frame;

/// Merge two frames along the specified axis with a vacuum gap at the interface.
///
/// - `axis`: "x", "y", or "z" (the direction along which structures are stacked).
/// - `gap`: vacuum gap thickness in Å at the interface (must be >= 0), measured
///   **perpendicular** to the interface — the same convention as `add_vacuum`.
///
/// Neither block is strained: both are moved rigidly. This is a padded stack for
/// building a starting box, not a coherent interface (ASE `stack` and pymatgen
/// `Interface.from_slabs` strain the second block onto the first one's lattice
/// instead, which is impossible when the in-plane sizes differ a lot).
///
/// The merged cell keeps the directions of `frame_a`'s lattice vectors:
/// - each in-plane vector is as long as the longer of the two blocks';
/// - the join vector is stretched until the cell's height along the interface
///   normal is `h_a + gap + h_b`, with `h` each block's interplanar spacing.
///
/// `frame_a`'s bottom face sits at height 0 and `frame_b`'s at `h_a + gap`. In the
/// interface plane both blocks are centred: each cell's centre is moved onto the
/// merged cell's centre line.
///
/// Each block is placed by its cell, so its atoms are first folded into that cell
/// (relative to `Cell::origin`) along its periodic axes; the box may start anywhere,
/// as a LAMMPS box with `xlo ≠ 0` does.
///
/// The atoms of `frame_a` come first, then those of `frame_b`, each with its label,
/// mass, magnetic moment and charge. The result takes `frame_a`'s `pbc` and a box origin
/// of 0; `energy`, `forces`, `stress`, `velocities`, `temperature` and `step` are `None`.
///
/// # Errors
/// - `ValidationError` if either frame has no cell, or one refused by `Cell::check`
///   (singular or not finite).
/// - `ValidationError` if `axis` is not "x"/"y"/"z".
/// - `ValidationError` if `gap` is negative or not finite.
/// - `ValidationError` if the two interfaces are not parallel (normals more than
///   1° apart): a gap along the normal has no meaning then.
/// - `ValidationError` if an atom lies outside its cell along a non-periodic axis.
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

    if !(gap >= 0.0 && gap.is_finite()) {
        return Err(ChemError::ValidationError(format!(
            "gap must be a finite number >= 0, got {}",
            gap
        )));
    }

    let cell_a = frame_a.cell.as_ref().ok_or_else(|| {
        ChemError::ValidationError("merge_frames: frame_a has no cell".into())
    })?;
    let cell_b = frame_b.cell.as_ref().ok_or_else(|| {
        ChemError::ValidationError("merge_frames: frame_b has no cell".into())
    })?;

    cell_a.check("merge_frames: frame_a's cell")?;
    cell_b.check("merge_frames: frame_b's cell")?;

    let (p, q) = ((join + 1) % 3, (join + 2) % 3);
    let rows = |c: &Cell| -> [Vector3<f64>; 3] { [0, 1, 2].map(|i| c.matrix.row(i).transpose()) };
    let (ra, rb) = (rows(cell_a), rows(cell_b));

    // 界面法向朝拼接矢量一侧；h = 拼接矢量在法向上的投影 = 面间距
    // Both cells passed `Cell::check`, so h = V / |face| is finite and > 0
    let normal = |r: &[Vector3<f64>; 3]| -> (Vector3<f64>, f64) {
        let n = r[p].cross(&r[q]);
        let h = n.dot(&r[join]) / n.norm();
        (n.normalize() * h.signum(), h.abs())
    };
    let (n_a, h_a) = normal(&ra);
    let (n_b, h_b) = normal(&rb);
    if n_a.dot(&n_b) < 1f64.to_radians().cos() {
        return Err(ChemError::ValidationError(format!(
            "merge_frames: the two interfaces are not parallel ({:.2}° apart along axis {axis}); \
             a gap measured along the normal is undefined",
            n_a.dot(&n_b).clamp(-1.0, 1.0).acos().to_degrees()
        )));
    }

    // ── 新胞：沿 A 的方向；面内取两块中较长者，拼接矢量拉到法向高度 H ──────────
    let height = h_a + gap + h_b;
    let mut new_rows = ra;
    for i in [p, q] {
        let (la, lb) = (ra[i].norm(), rb[i].norm());
        new_rows[i] = ra[i] * (la.max(lb) / la);
    }
    new_rows[join] = ra[join] * (height / h_a);
    let new_m = Matrix3::from_rows(&new_rows.map(|r| r.transpose()));

    // ── 整体平移：块的胞中心挪到新胞面内中心线上、该块自身的半高处 ────────────
    // 中心线上高度 z 的点 = 面内两矢量各一半 + 拼接矢量的 z/H
    // 块的胞中心 = 盒子原点 + 半条体对角线（LAMMPS 的 lo ≠ 0 时原点不在 0，审查 E-10）
    let centre_line = (new_rows[p] + new_rows[q]) * 0.5;
    let shift = |c: &Cell, r: &[Vector3<f64>; 3], bottom: f64, h: f64| -> Vector3<f64> {
        let target = centre_line + new_rows[join] * ((bottom + h / 2.0) / height);
        target - (c.origin + (r[0] + r[1] + r[2]) * 0.5)
    };
    let (shift_a, shift_b) = (shift(cell_a, &ra, 0.0, h_a), shift(cell_b, &rb, h_a + gap, h_b));

    // Placement goes by the cell faces, so atoms are first folded into the cell about
    // the box origin; otherwise an atom outside moves with its block into the other one
    // (review E-10: gap 2 Å, nearest distance 0.18 Å measured). An atom outside along a
    // non-periodic axis has no sensible fold: an error rather than two blocks overlapping.
    // Not `Frame::wrap_all` (folds about 0, skips non-periodic axes silently) nor
    // `Cell::wrap_position` (ignores pbc and origin); see the latter (review A-16)
    const WRAP_TOL: f64 = 1e-10;
    let moved = |frame: &Frame, cell: &Cell, t: Vector3<f64>, name: &str| -> Result<Vec<Atom>> {
        let mut atoms = Vec::with_capacity(frame.n_atoms());
        let mut stranded = 0usize;
        // inverted once, not per atom: the same product as cartesian_to_fractional (review R6-7)
        let inv = cell.matrix.transpose().try_inverse()
            .ok_or_else(|| ChemError::ValidationError("cell matrix is singular".into()))?;
        for atom in frame.atoms() {
            let f = inv * (atom.position - cell.origin);
            let mut a = atom.to_atom();
            for i in 0..3 {
                if f[i] >= -WRAP_TOL && f[i] < 1.0 + WRAP_TOL { continue; }
                if frame.pbc[i] {
                    a.position -= cell.matrix.row(i).transpose() * f[i].floor();
                } else {
                    stranded += 1;
                    break;
                }
            }
            a.position += t;
            atoms.push(a);
        }
        if stranded > 0 {
            return Err(ChemError::ValidationError(format!(
                "merge_frames: {stranded} atom(s) of {name} lie outside its cell along a \
                 non-periodic axis, where they cannot be folded back"
            )));
        }
        Ok(atoms)
    };

    // ── 组装 ──────────────────────────────────────────────────────────────────
    let mut atoms = moved(frame_a, cell_a, shift_a, "frame_a")?;
    atoms.extend(moved(frame_b, cell_b, shift_b, "frame_b")?);

    let mut frame = Frame::with_cell(Cell::from_matrix(new_m), frame_a.pbc);
    frame.set_atoms(atoms);
    Ok(frame)
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let pos_a = out.positions()[0];
        assert!((pos_a.z - 0.0).abs() < 1e-10);
        // frame_b 原子 z 平移 10 + 2 = 12
        let pos_b = out.positions()[1];
        assert!((pos_b.z - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_merge_different_sizes_centering() {
        // a: 10x10x10, b: 6x6x6 → 非拼接轴 max=10, b 居中偏移 (10-6)/2=2
        let a = cubic_frame(10.0, vec![("Fe", 0.0, 0.0, 0.0)]);
        let b = cubic_frame(6.0, vec![("O", 0.0, 0.0, 0.0)]);
        let out = merge_frames(&a, &b, "z", 0.0).unwrap();
        let pos_b = out.positions()[1];
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
        let pos_b = out.positions()[1];
        assert!((pos_b.x - 13.0).abs() < 1e-10);
    }

    // 审查 E-10：块按胞面放置，原子不在胞内时两块静默重叠
    #[test]
    fn test_atoms_outside_cell_are_folded_first() {
        let a = cubic_frame(10.0, vec![("Fe", 5.0, 5.0, 9.0)]);
        // 分数 z = -0.9，折回后 0.1：应落在 B 块底部 12 + 1，而不是 12 − 9 = 3（A 块内部）
        let b = cubic_frame(10.0, vec![("O", 5.0, 5.0, -9.0)]);
        let out = merge_frames(&a, &b, "z", 2.0).unwrap();
        assert!((out.positions()[1].z - 13.0).abs() < 1e-10, "实得 {}", out.positions()[1].z);
    }

    // LAMMPS 的 lo ≠ 0：块在 [lo, lo+L)，按盒子的真实位置放
    #[test]
    fn test_box_origin_is_respected() {
        let mut a = cubic_frame(10.0, vec![("Fe", 5.0, 5.0, 2.0)]);
        let mut b = cubic_frame(10.0, vec![("O", 5.0, 5.0, 1.92)]);
        a.cell.as_mut().unwrap().origin = Vector3::new(0.0, 0.0, 1.82);
        b.cell.as_mut().unwrap().origin = Vector3::new(0.0, 0.0, 1.82);
        let out = merge_frames(&a, &b, "z", 2.0).unwrap();
        assert!((out.positions()[0].z - 0.18).abs() < 1e-10, "A：实得 {}", out.positions()[0].z);
        assert!((out.positions()[1].z - 12.1).abs() < 1e-10, "B：实得 {}", out.positions()[1].z);
    }

    // 非周期方向越出盒子的原子没有合理的折回方式，报错而不是让两块重叠
    #[test]
    fn test_error_atom_outside_non_periodic_axis() {
        let a = cubic_frame(10.0, vec![("Fe", 5.0, 5.0, 5.0)]);
        let mut b = cubic_frame(10.0, vec![("O", 5.0, 5.0, -3.0)]);
        b.pbc = [true, true, false];
        assert!(merge_frames(&a, &b, "z", 2.0).is_err());
    }

    #[test]
    fn test_merge_atom_count() {
        let a = cubic_frame(10.0, vec![("Fe", 0.0, 0.0, 0.0), ("Fe", 1.0, 0.0, 0.0)]);
        let b = cubic_frame(10.0, vec![("O", 0.0, 0.0, 0.0)]);
        let out = merge_frames(&a, &b, "z", 1.0).unwrap();
        assert_eq!(out.n_atoms(), 3);
    }

    #[test]
    fn test_merge_keeps_atom_charges() {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut a = Frame::with_cell(cell.clone(), [true; 3]);
        a.add_atom(Atom::new("Fe", Vector3::zeros()));
        a.set_charge(0, Some(2.0));
        let mut b = Frame::with_cell(cell, [true; 3]);
        b.add_atom(Atom::new("O", Vector3::zeros()));
        b.set_charge(0, Some(-1.0));
        let out = merge_frames(&a, &b, "z", 1.0).unwrap();
        assert_eq!(out.atom(0).charge, Some(2.0));
        assert_eq!(out.atom(1).charge, Some(-1.0));
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
        assert!(out.forces().is_none());
        assert!(out.stress.is_none());
        assert!(out.velocities().is_none());
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
        let (h_a, h_b) = (
            a.cell.as_ref().unwrap().interplanar_spacings().unwrap()[2],
            b.cell.as_ref().unwrap().interplanar_spacings().unwrap()[2],
        );
        let out = merge_frames(&a, &b, "z", 1.0).unwrap();
        let cell = out.cell.as_ref().unwrap();
        let lengths = cell.lengths();
        // z：法向高度 = h_a + 1 + h_b（gap 沿法向，同 add_vacuum）
        let h = cell.interplanar_spacings().unwrap()[2];
        assert!((h - (h_a + 1.0 + h_b)).abs() < 1e-10, "实际 {h}");
        // x: max(5, 4) = 5
        assert!((lengths[0] - 5.0).abs() < 1e-10);
        // y: max(6, 5) = 6
        assert!((lengths[1] - 6.0).abs() < 1e-10);
        // 角度沿用 A 的
        let [al, be, ga] = cell.angles();
        assert!((al - 80.0).abs() < 1e-8 && (be - 85.0).abs() < 1e-8 && (ga - 95.0).abs() < 1e-8);
    }

    fn frame_rows(rows: [f64; 9], elem: &str) -> Frame {
        let cell = Cell::from_matrix(nalgebra::Matrix3::from_row_slice(&rows));
        let mut f = Frame::with_cell(cell, [true; 3]);
        f.add_atom(Atom::new(elem, Vector3::zeros()));
        f
    }

    /// 审查 D-M2 ①：文档说较小的一块居中，A 较窄时却没动
    #[test]
    fn test_narrower_a_is_centred_too() {
        let a = cubic_frame(6.0, vec![("Fe", 0.0, 0.0, 0.0)]);
        let b = cubic_frame(10.0, vec![("O", 0.0, 0.0, 0.0)]);
        let out = merge_frames(&a, &b, "z", 0.0).unwrap();
        let pa = out.positions()[0];
        assert!((pa - Vector3::new(2.0, 2.0, 0.0)).norm() < 1e-10, "A 应在面内居中，实际 {pa:?}");
        let pb = out.positions()[1];
        assert!((pb - Vector3::new(0.0, 0.0, 6.0)).norm() < 1e-10, "B 应紧贴 A 顶面，实际 {pb:?}");
    }

    /// 审查 D-M2 ②：B 曾沿自身晶格矢量平移，倾角与 A 不同时底面落不到 A 顶面
    #[test]
    fn test_tilted_b_sits_on_top_of_a() {
        let a = frame_rows([10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0], "Fe");
        let b = frame_rows([10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 3.0, 0.0, 10.0], "O");
        let out = merge_frames(&a, &b, "z", 2.0).unwrap();
        // B 原点在它的底面上：底面应在 A 顶面（z=10）之上 2 Å
        let pb = out.positions()[1];
        assert!((pb.z - 12.0).abs() < 1e-10, "B 底面应在 z = 12，实际 {pb:?}");
        // 面内按胞中心居中：B 中心 x = (10+3)/2 = 6.5，挪到 5
        assert!((pb.x + 1.5).abs() < 1e-10 && pb.y.abs() < 1e-10, "B 应按胞中心居中，实际 {pb:?}");
    }

    /// 审查 D-M2 ③：gap 曾沿晶格矢量量，倾斜拼接轴下真实间隙变小；
    /// 改为与 add_vacuum 相同的面间距口径
    #[test]
    fn test_gap_is_measured_along_the_normal() {
        let rows = [10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 3.0, 0.0, 10.0];
        let out = merge_frames(&frame_rows(rows, "Fe"), &frame_rows(rows, "O"), "z", 2.0).unwrap();
        let h = out.cell.as_ref().unwrap().interplanar_spacings().unwrap()[2];
        assert!((h - 22.0).abs() < 1e-10, "法向高度应为 10 + 2 + 10，实际 {h}");
        // 拼接矢量沿 A 的方向
        let c = out.cell.as_ref().unwrap().matrix.row(2).transpose();
        assert!((c - Vector3::new(6.6, 0.0, 22.0)).norm() < 1e-10, "实际 {c:?}");
        let pb = out.positions()[1];
        assert!((pb.z - 12.0).abs() < 1e-10, "实际 {pb:?}");
    }

    /// 两块的界面不平行时「沿法向的 gap」无定义
    #[test]
    fn test_error_non_parallel_interfaces() {
        let a = frame_rows([10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0], "Fe");
        let b = frame_rows([10.0, 0.0, 0.0, 0.0, 10.0, 3.0, 0.0, 0.0, 10.0], "O");
        let err = merge_frames(&a, &b, "z", 1.0).unwrap_err().to_string();
        assert!(err.contains("parallel"), "实际 {err}");
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
        // Review A-14: the check refused NaN and infinity, the message only spoke of >= 0
        for bad in [f64::NAN, f64::INFINITY] {
            let e = merge_frames(&a, &b, "z", bad).unwrap_err().to_string();
            assert!(e.contains("finite number >= 0"), "{e}");
        }
    }

    // Review A-13: the merge's own singularity test now goes through `Cell::check`,
    // which names the frame and also refuses a non-finite cell
    #[test]
    fn test_error_singular_or_non_finite_cell() {
        let a = cubic_frame(10.0, vec![("Fe", 0.0, 0.0, 0.0)]);
        let mut flat = cubic_frame(10.0, vec![("O", 0.0, 0.0, 0.0)]);
        flat.cell.as_mut().unwrap().matrix[(0, 0)] = 0.0;
        let e = merge_frames(&a, &flat, "z", 1.0).unwrap_err().to_string();
        assert!(e.contains("frame_b's cell is singular"), "{e}");
        let mut nan = cubic_frame(10.0, vec![("O", 0.0, 0.0, 0.0)]);
        nan.cell.as_mut().unwrap().matrix[(0, 1)] = f64::NAN;
        let e = merge_frames(&nan, &a, "z", 1.0).unwrap_err().to_string();
        assert!(e.contains("frame_a's cell is not finite"), "{e}");
    }
}
