//! Supercell construction.

use nalgebra::{Matrix3, Vector3};

use ferro_core::cell::Cell;
use ferro_core::error::{ChemError, Result};
use ferro_core::frame::Frame;

/// Largest supercell built, in atoms: past it a typo in `min_length` or the multiples
/// used to overflow `usize` or abort on the allocation instead of reporting. Also the
/// cap on a whole supercelled trajectory, which is held in memory at once.
pub const MAX_SUPERCELL_ATOMS: usize = 100_000_000;

fn too_many(what: &str) -> ChemError {
    ChemError::ValidationError(format!("{what} would exceed {MAX_SUPERCELL_ATOMS} atoms"))
}

// ── 自动搜索 ──────────────────────────────────────────────────────────────────

/// Find the smallest `[nx, ny, nz]` satisfying the given box-size and
/// atom-count constraints.
///
/// - `min_length`: minimum required **interplanar spacing** of the supercell along each
///   axis (Å) — the perpendicular thickness, which is what bounds a minimum-image cutoff
///   (`r_cut ≤ min_length / 2`). For a non-orthogonal cell it is smaller than the vector
///   length. Pass `0.0` to skip the length constraint.
/// - `min_atoms`: minimum required total atom count.
///   Pass `0` to skip the atom-count constraint.
///
/// When both constraints are `0` / `0.0`, returns `[1, 1, 1]`.
///
/// # Algorithm
/// 1. Set `ni = ceil(min_length / di)` for each axis, `di` the interplanar spacing.
/// 2. Iteratively increment the thinnest axis until
///    `nx * ny * nz * n_atoms >= min_atoms` (atom-count constraint).
///    Incrementing the thinnest axis keeps the box as isotropic as possible.
///
/// # Errors
/// - `ValidationError` if the cell is singular or not finite (`Cell::check`).
/// - `ValidationError` if `min_length` is negative, NaN or infinite.
/// - `ValidationError` if the result would hold more than 1e8 atoms.
pub fn find_supercell_dims(
    cell: &Cell,
    n_atoms: usize,
    min_length: f64,
    min_atoms: usize,
) -> Result<[usize; 3]> {
    // NaN would fall through to "no constraint", infinity to a usize::MAX multiple, and a
    // negative length to "no constraint" as well.
    if !(min_length.is_finite() && min_length >= 0.0) {
        return Err(ChemError::ValidationError(format!("min_length must be finite and >= 0, got {min_length}")));
    }
    if min_atoms > MAX_SUPERCELL_ATOMS {
        return Err(too_many(&format!("min_atoms {min_atoms}")));
    }
    // Interplanar spacings, not vector lengths: a minimum-image cutoff needs the
    // perpendicular thickness, and in a triclinic cell d < |a|, so lengths would give
    // too small a supercell. Repeating n times makes the spacing exactly n·d
    let lengths = cell.interplanar_spacings()?;
    // An empty cell counts as one atom in every cap below, so the multiples stay bounded
    // (and the f64 -> usize cast cannot saturate) even with no atoms to count
    let per_cell = n_atoms.max(1);
    // Names the supercell that broke the cap, not the constraint that led to it: the input
    // alone can be over it, or the last slab of the min_atoms search can carry it over
    let oversized = |d: [f64; 3]| ChemError::ValidationError(format!(
        "a {:.0}x{:.0}x{:.0} supercell of {n_atoms} atoms would hold {:.0} atoms, over {MAX_SUPERCELL_ATOMS}",
        d[0], d[1], d[2], d[0] * d[1] * d[2] * per_cell as f64));

    // Step 1: the smallest multiple per axis that meets min_length, bounded in f64
    // before the cast to usize can saturate
    let ceil_dim = |l: f64| (min_length / l).ceil().max(1.0);
    let nf = [ceil_dim(lengths[0]), ceil_dim(lengths[1]), ceil_dim(lengths[2])];
    if nf[0] * nf[1] * nf[2] * per_cell as f64 > MAX_SUPERCELL_ATOMS as f64 {
        return Err(oversized(nf));
    }
    let mut n = nf.map(|v| v as usize);

    // Step 2: grow the thinnest axis until there are enough atoms
    while n[0] * n[1] * n[2] * per_cell < min_atoms {
        let super_lengths = [
            n[0] as f64 * lengths[0],
            n[1] as f64 * lengths[1],
            n[2] as f64 * lengths[2],
        ];
        // The thinnest axis (partial_cmp cannot panic: the cell passed `Cell::check`,
        // so every spacing is finite)
        let shortest = super_lengths
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        n[shortest] += 1;
    }
    // The last step adds a whole slab, which can carry the count past min_atoms
    if n[0] * n[1] * n[2] * per_cell > MAX_SUPERCELL_ATOMS {
        return Err(oversized(n.map(|v| v as f64)));
    }

    Ok(n)
}

// ── 扩胞 ──────────────────────────────────────────────────────────────────────

/// Checks that supercelling every frame by `dims` stays within [`MAX_SUPERCELL_ATOMS`]
/// in total. `make_supercell` caps one frame; callers that build all frames at once
/// call this first, or a long trajectory of frames each under the cap still exhausts
/// memory (review fb7d8c3-1).
///
/// `atoms_per_frame` gives each input frame's atom count.
pub fn check_supercell_total(atoms_per_frame: impl IntoIterator<Item = usize>, dims: [usize; 3]) -> Result<()> {
    let (mut frames, mut atoms) = (0usize, Some(0usize));
    for n in atoms_per_frame {
        frames += 1;
        atoms = atoms.and_then(|a| a.checked_add(n));
    }
    let [nx, ny, nz] = dims;
    let total = atoms.and_then(|a| a.checked_mul(nx)).and_then(|v| v.checked_mul(ny)).and_then(|v| v.checked_mul(nz));
    if total.is_none_or(|t| t > MAX_SUPERCELL_ATOMS) {
        let shown = total.map_or("more than usize can hold".into(), |t| t.to_string());
        return Err(ChemError::ValidationError(format!(
            "a {nx}x{ny}x{nz} supercell of {frames} frame(s) would hold {shown} atoms in total, \
             over {MAX_SUPERCELL_ATOMS}; supercell a single frame instead"
        )));
    }
    Ok(())
}

/// Build a supercell by replicating `frame` by `(nx, ny, nz)` along the
/// three lattice vectors.
///
/// # Errors
/// - `ValidationError` if the frame has no cell (`cell` is `None`), or one refused by
///   `Cell::check` (singular or not finite).
/// - `ValidationError` if any of `nx`, `ny`, `nz` is zero.
/// - `ValidationError` if the supercell would hold more than 1e8 atoms.
///
/// # Notes
/// - Per-atom scalar fields (`label`, `mass`, `magmom`, `charge`) are
///   copied verbatim to every replica.
/// - `energy`, `forces`, `stress`, `velocities`, `temperature` and `step` are
///   **not** copied: they belong to a specific calculation result, not to the
///   structure.
/// - `pbc` and the box origin are kept; replicas extend along +a / +b / +c.
pub fn make_supercell(frame: &Frame, nx: usize, ny: usize, nz: usize) -> Result<Frame> {
    // ── Input checks ─────────────────────────────────────────────────────────
    if nx == 0 || ny == 0 || nz == 0 {
        return Err(ChemError::ValidationError(
            format!("supercell dimensions must be ≥ 1, got ({nx}, {ny}, {nz})"),
        ));
    }
    let total = nx.checked_mul(ny).and_then(|v| v.checked_mul(nz))
        .and_then(|v| v.checked_mul(frame.n_atoms()));
    if total.is_none_or(|t| t > MAX_SUPERCELL_ATOMS) {
        return Err(too_many(&format!("a {nx}x{ny}x{nz} supercell")));
    }
    let cell = frame.cell.as_ref().ok_or_else(|| {
        ChemError::ValidationError(
            "make_supercell requires a periodic frame with a cell".into(),
        )
    })?;
    cell.check("make_supercell: the cell")?;

    // ── 新格子矩阵（各行独立缩放）────────────────────────────────────────────
    let m = &cell.matrix;
    #[rustfmt::skip]
    let new_matrix = Matrix3::new(
        m[(0,0)] * nx as f64, m[(0,1)] * nx as f64, m[(0,2)] * nx as f64,
        m[(1,0)] * ny as f64, m[(1,1)] * ny as f64, m[(1,2)] * ny as f64,
        m[(2,0)] * nz as f64, m[(2,1)] * nz as f64, m[(2,2)] * nz as f64,
    );
    // 副本沿 +a/+b/+c 展开，盒子原点不动（审查 A-12：from_matrix 会把 LAMMPS 的 lo 清零）
    let new_cell = Cell { matrix: new_matrix, origin: cell.origin };

    // ── 晶格行向量（用于 Cartesian 偏移计算）──────────────────────────────────
    let a_vec = Vector3::new(m[(0, 0)], m[(0, 1)], m[(0, 2)]);
    let b_vec = Vector3::new(m[(1, 0)], m[(1, 1)], m[(1, 2)]);
    let c_vec = Vector3::new(m[(2, 0)], m[(2, 1)], m[(2, 2)]);

    // ── 扩展原子列表 ──────────────────────────────────────────────────────────
    let n_orig = frame.n_atoms();
    let n_replicas = nx * ny * nz;
    let mut new_atoms = Vec::with_capacity(n_orig * n_replicas);

    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                let offset = ix as f64 * a_vec + iy as f64 * b_vec + iz as f64 * c_vec;
                for atom in frame.atoms() {
                    let mut new_atom = atom.to_atom();
                    new_atom.position += offset;
                    new_atoms.push(new_atom);
                }
            }
        }
    }

    // ── 组装新帧 ──────────────────────────────────────────────────────────────
    let mut out = Frame::with_cell(new_cell, frame.pbc);
    out.set_atoms(new_atoms);
    Ok(out)
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
    fn a_negative_min_length_is_rejected() {
        let cell = bcc_frame(2.87).cell.unwrap();
        assert!(find_supercell_dims(&cell, 2, -1.0, 0).is_err());
    }

    /// Requests past MAX_SUPERCELL_ATOMS (1e8) are an error, not an overflow or a huge allocation.
    #[test]
    fn an_oversized_supercell_is_rejected() {
        let frame = bcc_frame(2.87);
        let cell = frame.cell.clone().unwrap();
        assert!(find_supercell_dims(&cell, 2, 1e9, 0).is_err());
        assert!(find_supercell_dims(&cell, 2, 0.0, usize::MAX).is_err());
        assert!(make_supercell(&frame, 1 << 22, 1 << 22, 1 << 22).is_err());
        assert!(make_supercell(&frame, 1000, 1000, 1000).is_err());
    }

    /// Each frame is under the cap, the trajectory is not (review fb7d8c3-1).
    #[test]
    fn a_trajectory_over_the_total_cap_is_rejected() {
        let err = check_supercell_total([1000; 2000], [10, 10, 10]).unwrap_err().to_string();
        assert!(err.contains("2000 frame(s)") && err.contains("2000000000"), "{err}");
        assert!(check_supercell_total([usize::MAX, 1], [1, 1, 1]).is_err());
        assert!(check_supercell_total([1, 1], [usize::MAX, 2, 1]).is_err());
    }

    #[test]
    fn a_trajectory_at_the_total_cap_passes() {
        assert!(check_supercell_total([1000; 100], [10, 10, 10]).is_ok());
        assert!(check_supercell_total([1000; 101], [10, 10, 10]).is_err());
    }

    /// The error names the supercell that is too big, not min_length / min_atoms that
    /// were within bounds (review fb7d8c3-4).
    #[test]
    fn an_oversized_supercell_error_names_the_supercell() {
        let cell = bcc_frame(2.87).cell.unwrap();
        // The input alone is over the cap; min_length 0 played no part
        let err = find_supercell_dims(&cell, 150_000_000, 0.0, 0).unwrap_err().to_string();
        assert!(err.contains("1x1x1 supercell of 150000000 atoms") && !err.contains("min_length"), "{err}");
        // min_atoms is within the cap, the last slab carries the count over it
        let err = find_supercell_dims(&cell, 60_000_000, 0.0, 100_000_000).unwrap_err().to_string();
        assert!(err.contains("2x1x1 supercell of 60000000 atoms would hold 120000000"), "{err}");
    }

    #[test]
    fn test_identity_supercell() {
        let frame = bcc_frame(2.87);
        let sc = make_supercell(&frame, 1, 1, 1).unwrap();
        assert_eq!(sc.n_atoms(), frame.n_atoms());
        let [la, lb, lc] = sc.cell.as_ref().unwrap().lengths();
        assert!((la - 2.87).abs() < 1e-10);
        assert!((lb - 2.87).abs() < 1e-10);
        assert!((lc - 2.87).abs() < 1e-10);
    }

    // 审查 A-12：LAMMPS 读进来的 lo 不能被超胞清零
    #[test]
    fn test_box_origin_preserved() {
        let mut frame = bcc_frame(2.87);
        let o = Vector3::new(1.5, -2.0, 3.0);
        frame.cell.as_mut().unwrap().origin = o;
        let sc = make_supercell(&frame, 2, 1, 3).unwrap();
        assert_eq!(sc.cell.as_ref().unwrap().origin, o, "超胞应保留盒子原点");
    }

    #[test]
    fn test_2x2x2_atom_count_and_cell() {
        let a = 2.87;
        let frame = bcc_frame(a);
        let sc = make_supercell(&frame, 2, 2, 2).unwrap();

        // 原子数：2 * 2*2*2 = 16
        assert_eq!(sc.n_atoms(), 16);

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
        assert_eq!(sc.n_atoms(), 2 * 2 * 3);
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
        let mut xs: Vec<f64> = sc.atoms().map(|a| a.position.x).collect();
        let mut ys: Vec<f64> = sc.atoms().map(|a| a.position.y).collect();
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
    fn test_every_replica_keeps_the_atom_charges() {
        let cell = Cell::from_lengths_angles(4.0, 4.0, 4.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));
        frame.set_charge(0, Some(-0.8));
        let sc = make_supercell(&frame, 3, 1, 1).unwrap();
        assert!(sc.atoms().all(|a| a.charge == Some(-0.8)));
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
    fn test_min_length_uses_interplanar_spacing() {
        // γ = 60° 的六方型胞：|a| = |b| = 5，但 a、b 方向的面间距只有 5·sin60° ≈ 4.33。
        // 要求 14 Å：按长度给 3（3·5 = 15），实际厚度 3·4.33 = 12.99 不够；按面间距给 4
        let cell = Cell::from_lengths_angles(5.0, 5.0, 20.0, 90.0, 90.0, 60.0).unwrap();
        let dims = find_supercell_dims(&cell, 1, 14.0, 0).unwrap();
        assert_eq!(dims, [4, 4, 1], "最小镜像要的是垂直厚度");
        let d = cell.interplanar_spacings().unwrap();
        for i in 0..3 {
            assert!(dims[i] as f64 * d[i] >= 14.0, "轴 {i} 的超胞厚度不足");
        }
    }

    #[test]
    fn test_find_dims_min_length_only() {
        // a=b=c=5 Å，min_length=12 Å → 需要 ceil(12/5)=3 → [3,3,3]
        let cell = Cell::from_lengths_angles(5.0, 5.0, 5.0, 90.0, 90.0, 90.0).unwrap();
        let dims = find_supercell_dims(&cell, 1, 12.0, 0).unwrap();
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
        let dims = find_supercell_dims(&cell, 2, 0.0, 20).unwrap();
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
        let dims = find_supercell_dims(&cell, 4, 8.0, 50).unwrap();
        assert!(dims[0] >= 3 && dims[1] >= 3 && dims[2] >= 3);
        assert!(dims[0] * dims[1] * dims[2] * 4 >= 50);
    }

    #[test]
    fn test_find_dims_no_constraint() {
        let cell = Cell::from_lengths_angles(5.0, 5.0, 5.0, 90.0, 90.0, 90.0).unwrap();
        assert_eq!(find_supercell_dims(&cell, 4, 0.0, 0).unwrap(), [1, 1, 1]);
    }

    #[test]
    fn a_min_length_that_is_not_finite_is_an_error() {
        let cell = Cell::from_lengths_angles(5.0, 5.0, 5.0, 90.0, 90.0, 90.0).unwrap();
        for bad in [f64::NAN, f64::INFINITY] {
            assert!(find_supercell_dims(&cell, 4, bad, 0).is_err(), "{bad}");
        }
    }

    // Review A-1: a NaN cell gave [1, 1, 1] or [50, 1, 1] instead of an error
    #[test]
    fn a_cell_that_is_not_finite_or_singular_is_an_error() {
        let ok = Cell::from_lengths_angles(5.0, 5.0, 5.0, 90.0, 90.0, 90.0).unwrap();
        let mut nan = ok.clone();
        nan.matrix[(0, 1)] = f64::NAN;
        let mut flat = ok.clone();
        flat.matrix[(2, 2)] = 0.0;
        for cell in [nan, flat] {
            assert!(find_supercell_dims(&cell, 4, 12.0, 0).is_err());
            assert!(find_supercell_dims(&cell, 4, 0.0, 50).is_err());
            let mut frame = Frame::with_cell(cell, [true; 3]);
            frame.add_atom(Atom::new("Fe", Vector3::zeros()));
            let e = make_supercell(&frame, 2, 1, 1).unwrap_err().to_string();
            assert!(e.contains("make_supercell: the cell is"), "{e}");
        }
    }
}
