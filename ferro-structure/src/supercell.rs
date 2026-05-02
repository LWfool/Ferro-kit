//! Supercell construction.

use nalgebra::{Matrix3, Vector3};

use ferro_core::cell::Cell;
use ferro_core::error::{ChemError, Result};
use ferro_core::frame::Frame;

// ── 自动搜索 ──────────────────────────────────────────────────────────────────

/// Find the smallest `[nx, ny, nz]` satisfying the given box-size and
/// atom-count constraints.
///
/// - `min_length`: minimum required supercell length along each axis (Å).
///   Pass `0.0` to skip the length constraint.
/// - `min_atoms`: minimum required total atom count.
///   Pass `0` to skip the atom-count constraint.
///
/// When both constraints are `0` / `0.0`, returns `[1, 1, 1]`.
///
/// # Algorithm
/// 1. Set `ni = ceil(min_length / Li)` for each axis (length constraint).
/// 2. Iteratively increment the shortest axis until
///    `nx * ny * nz * n_atoms >= min_atoms` (atom-count constraint).
///    Incrementing the shortest axis keeps the box as isotropic as possible.
pub fn find_supercell_dims(
    cell: &Cell,
    n_atoms: usize,
    min_length: f64,
    min_atoms: usize,
) -> [usize; 3] {
    let lengths = cell.lengths();

    // 步骤 1：用长度约束确定每轴最小倍数
    let ceil_dim = |l: f64| -> usize {
        if min_length > 0.0 {
            ((min_length / l).ceil() as usize).max(1)
        } else {
            1
        }
    };
    let mut n = [ceil_dim(lengths[0]), ceil_dim(lengths[1]), ceil_dim(lengths[2])];

    // 步骤 2：递增最短轴直到原子数满足要求
    while n[0] * n[1] * n[2] * n_atoms.max(1) < min_atoms {
        let super_lengths = [
            n[0] as f64 * lengths[0],
            n[1] as f64 * lengths[1],
            n[2] as f64 * lengths[2],
        ];
        // 找最短轴（f64 partial_cmp 不会 panic，因为 lengths 来自 norm()）
        let shortest = super_lengths
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        n[shortest] += 1;
    }

    n
}

// ── 扩胞 ──────────────────────────────────────────────────────────────────────

/// Build a supercell by replicating `frame` by `(nx, ny, nz)` along the
/// three lattice vectors.
///
/// # Errors
/// - `ValidationError` if the frame has no cell (`cell` is `None`).
/// - `ValidationError` if any of `nx`, `ny`, `nz` is zero.
///
/// # Notes
/// - Per-atom scalar fields (`label`, `mass`, `magmom`, `charge`) are
///   copied verbatim to every replica.
/// - Intra-cell bonds (if present) are remapped; cross-cell bonds are not
///   added (periodic images are handled by the cell geometry itself).
/// - `energy`, `forces`, `stress`, `velocities` are **not** copied: they
///   belong to a specific calculation result, not to the structure.
/// - Frame `charge` scales linearly with the number of replicas.
/// - Frame `multiplicity` is reset to 1 (supercells are typically used as
///   MD starting structures where spin is re-specified separately).
pub fn make_supercell(frame: &Frame, nx: usize, ny: usize, nz: usize) -> Result<Frame> {
    // ── 输入校验 ──────────────────────────────────────────────────────────────
    if nx == 0 || ny == 0 || nz == 0 {
        return Err(ChemError::ValidationError(
            format!("supercell dimensions must be ≥ 1, got ({nx}, {ny}, {nz})"),
        ));
    }
    let cell = frame.cell.as_ref().ok_or_else(|| {
        ChemError::ValidationError(
            "make_supercell requires a periodic frame with a cell".into(),
        )
    })?;

    // ── 新格子矩阵（各行独立缩放）────────────────────────────────────────────
    let m = &cell.matrix;
    #[rustfmt::skip]
    let new_matrix = Matrix3::new(
        m[(0,0)] * nx as f64, m[(0,1)] * nx as f64, m[(0,2)] * nx as f64,
        m[(1,0)] * ny as f64, m[(1,1)] * ny as f64, m[(1,2)] * ny as f64,
        m[(2,0)] * nz as f64, m[(2,1)] * nz as f64, m[(2,2)] * nz as f64,
    );
    let new_cell = Cell::from_matrix(new_matrix);

    // ── 晶格行向量（用于 Cartesian 偏移计算）──────────────────────────────────
    let a_vec = Vector3::new(m[(0, 0)], m[(0, 1)], m[(0, 2)]);
    let b_vec = Vector3::new(m[(1, 0)], m[(1, 1)], m[(1, 2)]);
    let c_vec = Vector3::new(m[(2, 0)], m[(2, 1)], m[(2, 2)]);

    // ── 扩展原子列表 ──────────────────────────────────────────────────────────
    let n_orig = frame.atoms.len();
    let n_replicas = nx * ny * nz;
    let mut new_atoms = Vec::with_capacity(n_orig * n_replicas);

    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let offset = ix as f64 * a_vec + iy as f64 * b_vec + iz as f64 * c_vec;
                for atom in &frame.atoms {
                    let mut new_atom = atom.clone();
                    new_atom.position += offset;
                    new_atoms.push(new_atom);
                }
            }
        }
    }

    // ── 重映射键（仅胞内键，按副本偏移下标）──────────────────────────────────
    let new_bonds = frame.bonds.as_ref().map(|bonds| {
        (0..n_replicas)
            .flat_map(|replica| {
                let offset = replica * n_orig;
                bonds.iter().map(move |&(i, j)| (i + offset, j + offset))
            })
            .collect::<Vec<_>>()
    });

    // ── 组装新帧 ──────────────────────────────────────────────────────────────
    Ok(Frame {
        atoms: new_atoms,
        cell: Some(new_cell),
        pbc: frame.pbc,
        charge: frame.charge * n_replicas as i32,
        multiplicity: 1,
        bonds: new_bonds,
        energy: None,
        forces: None,
        stress: None,
        velocities: None,
    })
}

// ── 测试 ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::atom::Atom;

    /// 构造一个正交晶胞，含两个原子（分数坐标 0,0,0 和 0.5,0.5,0.5）
    fn bcc_frame(a: f64) -> Frame {
        let cell = Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(a * 0.5, a * 0.5, a * 0.5)));
        frame
    }

    #[test]
    fn test_identity_supercell() {
        let frame = bcc_frame(2.87);
        let sc = make_supercell(&frame, 1, 1, 1).unwrap();
        assert_eq!(sc.atoms.len(), frame.atoms.len());
        let [la, lb, lc] = sc.cell.as_ref().unwrap().lengths();
        assert!((la - 2.87).abs() < 1e-10);
        assert!((lb - 2.87).abs() < 1e-10);
        assert!((lc - 2.87).abs() < 1e-10);
    }

    #[test]
    fn test_2x2x2_atom_count_and_cell() {
        let a = 2.87;
        let frame = bcc_frame(a);
        let sc = make_supercell(&frame, 2, 2, 2).unwrap();

        // 原子数：2 * 2*2*2 = 16
        assert_eq!(sc.atoms.len(), 16);

        // 格子边长：各方向乘以倍数
        let [la, lb, lc] = sc.cell.as_ref().unwrap().lengths();
        assert!((la - 2.0 * a).abs() < 1e-10);
        assert!((lb - 2.0 * a).abs() < 1e-10);
        assert!((lc - 2.0 * a).abs() < 1e-10);
    }

    #[test]
    fn test_2x3x1_mixed() {
        let a = 3.0;
        let frame = bcc_frame(a);
        let sc = make_supercell(&frame, 2, 3, 1).unwrap();
        assert_eq!(sc.atoms.len(), 2 * 2 * 3 * 1);
        let [la, lb, lc] = sc.cell.as_ref().unwrap().lengths();
        assert!((la - 2.0 * a).abs() < 1e-10);
        assert!((lb - 3.0 * a).abs() < 1e-10);
        assert!((lc - a).abs() < 1e-10);
    }

    #[test]
    fn test_atom_positions_correct() {
        // 单原子正交晶胞，边长 4 Å，原子在原点
        let cell = Cell::from_lengths_angles(4.0, 4.0, 4.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));

        let sc = make_supercell(&frame, 2, 2, 1).unwrap();
        // 4 个副本的原子坐标应为 (0,0,0), (4,0,0), (0,4,0), (4,4,0)
        let mut xs: Vec<f64> = sc.atoms.iter().map(|a| a.position.x).collect();
        let mut ys: Vec<f64> = sc.atoms.iter().map(|a| a.position.y).collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((xs[0] - 0.0).abs() < 1e-10);
        assert!((xs[1] - 0.0).abs() < 1e-10);
        assert!((xs[2] - 4.0).abs() < 1e-10);
        assert!((xs[3] - 4.0).abs() < 1e-10);
        assert!((ys[0] - 0.0).abs() < 1e-10);
        assert!((ys[1] - 0.0).abs() < 1e-10);
        assert!((ys[2] - 4.0).abs() < 1e-10);
        assert!((ys[3] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_charge_scales() {
        let cell = Cell::from_lengths_angles(4.0, 4.0, 4.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.charge = 2;
        frame.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));
        let sc = make_supercell(&frame, 3, 1, 1).unwrap();
        assert_eq!(sc.charge, 6);
    }

    #[test]
    fn test_bonds_remapped() {
        let cell = Cell::from_lengths_angles(4.0, 4.0, 4.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("H", Vector3::new(1.0, 0.0, 0.0)));
        frame.bonds = Some(vec![(0, 1)]);

        let sc = make_supercell(&frame, 2, 1, 1).unwrap();
        let bonds = sc.bonds.as_ref().unwrap();
        // 2 副本 * 1 键 = 2 键
        assert_eq!(bonds.len(), 2);
        assert!(bonds.contains(&(0, 1)));
        assert!(bonds.contains(&(2, 3)));
    }

    #[test]
    fn test_error_no_cell() {
        let mut frame = Frame::new();
        frame.add_atom(Atom::new("Fe", Vector3::zeros()));
        assert!(make_supercell(&frame, 2, 2, 2).is_err());
    }

    #[test]
    fn test_error_zero_dim() {
        let frame = bcc_frame(2.87);
        assert!(make_supercell(&frame, 0, 1, 1).is_err());
    }

    #[test]
    fn test_find_dims_min_length_only() {
        // a=b=c=5 Å，min_length=12 Å → 需要 ceil(12/5)=3 → [3,3,3]
        let cell = Cell::from_lengths_angles(5.0, 5.0, 5.0, 90.0, 90.0, 90.0).unwrap();
        let dims = find_supercell_dims(&cell, 1, 12.0, 0);
        assert_eq!(dims, [3, 3, 3]);
    }

    #[test]
    fn test_find_dims_min_atoms_only() {
        // a=5,b=6,c=7；2 atoms；min_atoms=20 → 1x1x1 给 2，不够
        // 递增最短轴（5 Å，即 x）：2x1x1=4, 3x1x1=6, 4x1x1=8, 5x1x1=10
        // 然后 y 也只有 6<10，比较后继续增 x：6x1x1=12, ...
        // 实际上：min_length=0 时先 [1,1,1]，然后每次增最短的超晶胞轴
        // 最短 super_length: 1*5=5 < 1*6=6 < 1*7=7 → 增 x
        // [2,1,1]: 2*5=10>6>7 → 最短 y=6 → 增 y
        // [2,2,1]: 10,12,7 → 最短 z → 增 z
        // [2,2,2]: 10,12,14 → 2*2*2*2=16 < 20 → 最短 x=10 → 增 x
        // [3,2,2]: 15,12,14 → 最短 y=12 → 增 y
        // [3,3,2]: 15,18,14 → 最短 z → 增 z
        // [3,3,3]: 15,18,21 → 3*3*3*2=54 ≥ 20 ✓
        let cell = Cell::from_lengths_angles(5.0, 6.0, 7.0, 90.0, 90.0, 90.0).unwrap();
        let dims = find_supercell_dims(&cell, 2, 0.0, 20);
        assert!(dims[0] * dims[1] * dims[2] * 2 >= 20);
        // 同时验证各轴 super_length 尽量平衡（最大/最小比 < 2）
        let cell_lengths = cell.lengths();
        let super_lengths: Vec<f64> = dims.iter().zip(cell_lengths.iter()).map(|(&n, &l)| n as f64 * l).collect();
        let max_l = super_lengths.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_l = super_lengths.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(max_l / min_l < 2.5, "box too anisotropic: {:?}", super_lengths);
    }

    #[test]
    fn test_find_dims_both_constraints() {
        // a=b=c=3 Å，4 atoms，min_length=8, min_atoms=50
        // 长度约束：ceil(8/3)=3 → [3,3,3] → 3*3*3*4=108 ≥ 50 ✓
        let cell = Cell::from_lengths_angles(3.0, 3.0, 3.0, 90.0, 90.0, 90.0).unwrap();
        let dims = find_supercell_dims(&cell, 4, 8.0, 50);
        assert!(dims[0] >= 3 && dims[1] >= 3 && dims[2] >= 3);
        assert!(dims[0] * dims[1] * dims[2] * 4 >= 50);
    }

    #[test]
    fn test_find_dims_no_constraint() {
        let cell = Cell::from_lengths_angles(5.0, 5.0, 5.0, 90.0, 90.0, 90.0).unwrap();
        assert_eq!(find_supercell_dims(&cell, 4, 0.0, 0), [1, 1, 1]);
    }
}
