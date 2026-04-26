//! Hard-sphere spatial occupancy map.
//!
//! For each voxel, counts how many (frame, atom) pairs have the atom
//! within `radius` Å of the voxel centre.  Applies the minimum-image
//! convention for periodic cells.
//!
//! Algorithm mirrors `examples/code1/cube_radius.c`:
//! 1. voxel centre in fractional: `(ix+0.5)/nx`, …
//! 2. fractional displacement atom → voxel, minimum-image applied
//! 3. convert to Cartesian: `cart = matrix.transpose() * frac_disp`
//! 4. if `|cart| < radius`: `rho[ix][iy][iz] += 1`
//!
//! Rust 侧改进：
//! - rayon 帧级并行
//! - 对每个原子计算 bounding-box，只遍历可能在 radius 内的 voxel

use nexflux_core::{CubeData, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use ndarray::Array3;
use rayon::prelude::*;

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Parameters for the hard-sphere spatial occupancy calculation.
#[derive(Debug, Clone)]
pub struct CubeRadiusParams {
    /// Grid divisions along a axis (CLI: `--nx`)
    pub nx: usize,
    /// Grid divisions along b axis (CLI: `--ny`)
    pub ny: usize,
    /// Grid divisions along c axis (CLI: `--nz`)
    pub nz: usize,
    /// Hard-sphere cutoff radius \[Å\] (CLI: `--radius`)
    pub radius: f64,
    /// Elements to include (`None` = all atoms; CLI: `--elements Li,Na`)
    pub elements: Option<Vec<String>>,
}

impl Default for CubeRadiusParams {
    fn default() -> Self {
        Self { nx: 100, ny: 100, nz: 100, radius: 0.7, elements: None }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Result of a hard-sphere occupancy cube calculation.
pub struct CubeRadiusResult {
    /// 3-D voxel data (raw counts) + time-averaged structure
    pub cube: CubeData,
    /// Number of frames that contributed to the grid
    pub n_frames: usize,
    /// Number of selected atoms (counted from first valid frame)
    pub n_atoms: usize,
    pub params: CubeRadiusParams,
}

// ─── 内部辅助 ────────────────────────────────────────────────────────────────

/// Compute the voxel search extent along each crystal axis.
///
/// From `frac = (M^T)^{-1} * cart`, the norm of row i of `(M^T)^{-1}` gives the maximum
/// fractional-coordinate offset in direction i when the Cartesian distance equals `radius`.
/// Multiply by the voxel count along that axis, ceil, and add 1 as a safety margin.
fn search_range(
    mat: &Matrix3<f64>,
    radius: f64,
    nx: usize,
    ny: usize,
    nz: usize,
) -> (i64, i64, i64) {
    // (M^T)^{-1} : Cartesian → 分数坐标的变换矩阵
    let m_t_inv = mat
        .transpose()
        .try_inverse()
        .unwrap_or_else(|| *mat); // 奇异矩阵时退化（不应发生）
    let sx = (radius * m_t_inv.row(0).norm() * nx as f64).ceil() as i64 + 1;
    let sy = (radius * m_t_inv.row(1).norm() * ny as f64).ceil() as i64 + 1;
    let sz = (radius * m_t_inv.row(2).norm() * nz as f64).ceil() as i64 + 1;
    (sx, sy, sz)
}

/// 处理单帧，返回该帧的计数贡献。
///
/// 返回 `None` 表示帧没有周期性 cell（跳过该帧）。
fn process_frame(
    frame: &Frame,
    params: &CubeRadiusParams,
    search: (i64, i64, i64),
) -> Option<Array3<f64>> {
    let cell = frame.cell.as_ref()?;
    let (nx, ny, nz) = (params.nx, params.ny, params.nz);
    let (sx, sy, sz) = search;
    let radius2 = params.radius * params.radius;
    // cart_disp = mat_t * frac_disp（与 C 代码 x = dx*A[0][0] + … 等价）
    let mat_t = cell.matrix.transpose();

    let mut rho = Array3::<f64>::zeros((nx, ny, nz));

    for atom in &frame.atoms {
        if let Some(elems) = &params.elements {
            if !elems.contains(&atom.element) { continue; }
        }

        // 原子分数坐标，归一化到 [0, 1)
        let frac = cell.cartesian_to_fractional(atom.position);
        let fx = frac.x.rem_euclid(1.0);
        let fy = frac.y.rem_euclid(1.0);
        let fz = frac.z.rem_euclid(1.0);

        // 原子所在 voxel 的整数索引
        let cx = (fx * nx as f64).floor() as i64;
        let cy = (fy * ny as f64).floor() as i64;
        let cz = (fz * nz as f64).floor() as i64;

        // 只遍历 bounding-box 内的 voxel
        for dix in -sx..=sx {
            for diy in -sy..=sy {
                for diz in -sz..=sz {
                    // 周期性折叠 voxel 索引
                    let ix = ((cx + dix).rem_euclid(nx as i64)) as usize;
                    let iy = ((cy + diy).rem_euclid(ny as i64)) as usize;
                    let iz = ((cz + diz).rem_euclid(nz as i64)) as usize;

                    // voxel 中心的分数坐标
                    let vx = (ix as f64 + 0.5) / nx as f64;
                    let vy = (iy as f64 + 0.5) / ny as f64;
                    let vz = (iz as f64 + 0.5) / nz as f64;

                    // 分数位移，应用最小像约定（对应 C 代码 dx = dx - rint(dx)）
                    let d = Vector3::new(fx - vx, fy - vy, fz - vz);
                    let dfrac = d.map(|v| v - v.round());

                    // 转为 Cartesian 并判断距离（对应 C 代码坐标转换 + if r < RADIUS）
                    if (mat_t * dfrac).norm_squared() < radius2 {
                        rho[[ix, iy, iz]] += 1.0;
                    }
                }
            }
        }
    }
    Some(rho)
}

/// Build a time-averaged representative frame for the cube file header.
fn build_avg_frame(traj: &Trajectory) -> Frame {
    let ref_f = traj.frames.first().unwrap();
    let n = ref_f.atoms.len();
    let mut pos_sum = vec![Vector3::<f64>::zeros(); n];
    let mut valid = 0usize;
    for frame in &traj.frames {
        if frame.atoms.len() != n { continue; }
        for (s, a) in pos_sum.iter_mut().zip(frame.atoms.iter()) {
            *s += a.position;
        }
        valid += 1;
    }
    let mut out = ref_f.clone();
    if valid > 0 {
        for (a, s) in out.atoms.iter_mut().zip(pos_sum.iter()) {
            a.position = s / valid as f64;
        }
    }
    out
}

// ─── 主函数 ──────────────────────────────────────────────────────────────────

/// Calculate a hard-sphere spatial occupancy map from a trajectory.
///
/// Each voxel accumulates the total count of (frame, atom) pairs where
/// the atom lies within `params.radius` Å of the voxel centre.
///
/// Returns `None` if no frame with a periodic cell is found.
pub fn calc_cube_radius(
    traj: &Trajectory,
    params: &CubeRadiusParams,
) -> Option<CubeRadiusResult> {
    let (nx, ny, nz) = (params.nx, params.ny, params.nz);

    let ref_frame = traj.frames.iter().find(|f| f.cell.is_some())?;
    let ref_cell = ref_frame.cell.as_ref().unwrap();

    let search = search_range(&ref_cell.matrix, params.radius, nx, ny, nz);

    // 帧级并行
    let results: Vec<Array3<f64>> = traj
        .frames
        .par_iter()
        .filter_map(|f| process_frame(f, params, search))
        .collect();

    if results.is_empty() { return None; }

    let n_frames = results.len();

    // 累加所有帧的贡献
    let data = results.into_iter().fold(
        Array3::<f64>::zeros((nx, ny, nz)),
        |mut acc, rho| { acc += &rho; acc },
    );

    let n_atoms = ref_frame
        .atoms
        .iter()
        .filter(|a| match &params.elements {
            Some(elems) => elems.contains(&a.element),
            None => true,
        })
        .count();

    // spacing: 每个 voxel 对应的晶格步长矩阵（行 i = 第 i 轴单步向量）
    let m = ref_cell.matrix;
    let spacing = Matrix3::new(
        m[(0, 0)] / nx as f64, m[(0, 1)] / nx as f64, m[(0, 2)] / nx as f64,
        m[(1, 0)] / ny as f64, m[(1, 1)] / ny as f64, m[(1, 2)] / ny as f64,
        m[(2, 0)] / nz as f64, m[(2, 1)] / nz as f64, m[(2, 2)] / nz as f64,
    );

    let cube = CubeData {
        frame: build_avg_frame(traj),
        data,
        origin: Vector3::zeros(),
        spacing,
    };

    Some(CubeRadiusResult { cube, n_frames, n_atoms, params: params.clone() })
}

// ─── 测试 ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use nexflux_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    fn make_traj(positions: Vec<(f64, f64, f64)>, elements: Vec<&str>) -> Trajectory {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        for (pos, elem) in positions.iter().zip(elements.iter()) {
            frame.add_atom(Atom::new(*elem, Vector3::new(pos.0, pos.1, pos.2)));
        }
        let mut traj = Trajectory::new();
        traj.frames.push(frame);
        traj
    }

    #[test]
    fn test_no_cell_returns_none() {
        let mut frame = Frame::default();
        frame.add_atom(Atom::new("Li", Vector3::zeros()));
        let mut traj = Trajectory::new();
        traj.frames.push(frame);
        let params = CubeRadiusParams { nx: 10, ny: 10, nz: 10, ..Default::default() };
        assert!(calc_cube_radius(&traj, &params).is_none());
    }

    #[test]
    fn test_atom_marks_nearby_voxels() {
        // 10Å 立方盒子，10×10×10 网格 → voxel 边长 1Å
        // 原子在 (5,5,5)，radius=1.2Å → 应标记以 (5,5,5) 为中心的若干 voxel
        let traj = make_traj(vec![(5.0, 5.0, 5.0)], vec!["Li"]);
        let params = CubeRadiusParams {
            nx: 10, ny: 10, nz: 10, radius: 1.2, elements: None,
        };
        let res = calc_cube_radius(&traj, &params).unwrap();
        // 原子所在 voxel (5,5,5)（中心在 5.5Å，原子在 5.0Å，距离 0.5Å < 1.2Å）→ 应被标记
        assert!(res.cube.data[[5, 5, 5]] > 0.0);
        // 远处 voxel (0,0,0)（中心在 0.5Å）距离 > 1.2Å → 应为 0
        assert_eq!(res.cube.data[[0, 0, 0]], 0.0);
    }

    #[test]
    fn test_multi_frame_accumulates() {
        // 同一位置的两帧 → 计数翻倍
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for _ in 0..2 {
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("Li", Vector3::new(5.0, 5.0, 5.0)));
            traj.frames.push(frame);
        }
        let params = CubeRadiusParams {
            nx: 10, ny: 10, nz: 10, radius: 1.2, elements: None,
        };
        let res1 = calc_cube_radius(&{
            let mut t = Trajectory::new();
            t.frames.push(traj.frames[0].clone());
            t
        }, &params).unwrap();
        let res2 = calc_cube_radius(&traj, &params).unwrap();
        // 两帧计数应恰好是单帧的两倍
        for ix in 0..10 { for iy in 0..10 { for iz in 0..10 {
            assert!((res2.cube.data[[ix,iy,iz]] - 2.0 * res1.cube.data[[ix,iy,iz]]).abs() < 1e-10);
        }}}
        assert_eq!(res2.n_frames, 2);
    }

    #[test]
    fn test_element_filter() {
        // 两个原子：Li 在 (5,5,5)，O 在 (1,1,1)；只统计 Li
        let traj = make_traj(vec![(5.0, 5.0, 5.0), (1.0, 1.0, 1.0)], vec!["Li", "O"]);
        let params = CubeRadiusParams {
            nx: 10, ny: 10, nz: 10,
            radius: 1.2,
            elements: Some(vec!["Li".to_string()]),
        };
        let res = calc_cube_radius(&traj, &params).unwrap();
        // O 所在区域（voxel (1,1,1) 中心 1.5Å）距离 O (1,1,1) 约 0.87Å < 1.2Å
        // 但 Li 与之距离 ≈ 6.9Å > 1.2Å，所以 voxel (1,1,1) 应为 0
        assert_eq!(res.cube.data[[1, 1, 1]], 0.0);
        // Li 附近 voxel (5,5,5) 应有计数
        assert!(res.cube.data[[5, 5, 5]] > 0.0);
        assert_eq!(res.n_atoms, 1);
    }

    #[test]
    fn test_periodic_wrap() {
        // 原子在 (10.5, 5.0, 5.0) → 分数 (1.05, 0.5, 0.5) → 折叠到 (0.05, 0.5, 0.5)
        // 等同于原子在 (0.5, 5.0, 5.0)，voxel (0,5,5) 中心在 0.5Å，距离 0Å < radius
        let traj = make_traj(vec![(10.5, 5.0, 5.0)], vec!["Li"]);
        let params = CubeRadiusParams {
            nx: 10, ny: 10, nz: 10, radius: 1.2, elements: None,
        };
        let res = calc_cube_radius(&traj, &params).unwrap();
        // voxel (0,5,5) 的中心在 (0.5, 5.5, 5.5)Å，原子在 (0.5, 5.0, 5.0)Å
        // 距离 = sqrt(0 + 0.25 + 0.25) ≈ 0.71Å < 1.2Å → 应被标记
        assert!(res.cube.data[[0, 5, 5]] > 0.0);
    }

    #[test]
    fn test_search_range_cubic() {
        // 10Å 立方盒子，10 voxel，radius=1.5Å → 每轴搜索 ceil(1.5/1.0)+1 = 3
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let (sx, sy, sz) = search_range(&cell.matrix, 1.5, 10, 10, 10);
        assert!(sx >= 2 && sy >= 2 && sz >= 2);
    }

    #[test]
    fn test_grid_shape() {
        let traj = make_traj(vec![(5.0, 5.0, 5.0)], vec!["Li"]);
        let params = CubeRadiusParams { nx: 4, ny: 5, nz: 6, radius: 0.7, elements: None };
        let res = calc_cube_radius(&traj, &params).unwrap();
        assert_eq!(res.cube.shape(), (4, 5, 6));
    }
}
