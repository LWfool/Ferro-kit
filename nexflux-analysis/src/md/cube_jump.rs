//! Jump distance spatial distribution.
//!
//! For each atom at each time step, computes the displacement over a lag
//! `tau` frames using unwrapped fractional coordinates.  Events where
//! |Δr| exceeds `threshold` are recorded as "jumps" and accumulated into a
//! 3D voxel grid at the atom's wrapped position (start, end, or midpoint).
//!
//! Typical use case: identify spatial hotspots of large-displacement (jump) events in ionic glasses.
//! Algorithm:
//!   1. For each atom, extract fractional coordinates across the full trajectory and unwrap (remove PBC jumps).
//!   2. For each time window [t, t+tau], compute the true Cartesian displacement from the unwrapped Δfrac.
//!   3. If |Δr| > threshold, record a jump event by incrementing the voxel at the wrapped position.
//!
//! Parallelism: per-atom par_iter; each atom independently unwraps and accumulates, then results are reduced.

use nexflux_core::{CubeData, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use ndarray::Array3;
use rayon::prelude::*;

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Position at which a jump event is recorded on the voxel grid.
#[derive(Debug, Clone, PartialEq)]
pub enum JumpPosition {
    /// Record at the atom's position at the beginning of the jump window
    Start,
    /// Record at the atom's position at the end of the jump window
    End,
    /// Record at the fractional midpoint of the unwrapped trajectory segment
    Midpoint,
}

/// Parameters for the jump distance spatial distribution.
#[derive(Debug, Clone)]
pub struct CubeJumpParams {
    /// Grid divisions along a axis
    pub nx: usize,
    /// Grid divisions along b axis
    pub ny: usize,
    /// Grid divisions along c axis
    pub nz: usize,
    /// Displacement lag \[frames\]
    pub tau: usize,
    /// Minimum displacement to qualify as a jump \[Å\]
    pub threshold: f64,
    /// Elements to include (`None` = all atoms)
    pub elements: Option<Vec<String>>,
    /// Where on the trajectory to record the jump event
    pub record_at: JumpPosition,
}

impl Default for CubeJumpParams {
    fn default() -> Self {
        Self {
            nx: 50, ny: 50, nz: 50,
            tau: 1,
            threshold: 1.0,
            elements: None,
            record_at: JumpPosition::Start,
        }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Result of a jump spatial distribution calculation.
pub struct CubeJumpResult {
    /// 3-D voxel grid (raw jump-event counts) + time-averaged structure
    pub cube: CubeData,
    /// Total number of frames in the trajectory
    pub n_frames: usize,
    /// Number of selected atoms
    pub n_atoms: usize,
    /// Total number of recorded jump events
    pub n_jumps: usize,
    pub params: CubeJumpParams,
}

// ─── 内部辅助 ────────────────────────────────────────────────────────────────

/// 原地 unwrap 单原子的分数坐标序列，消除 PBC 穿越引起的跳变。
fn unwrap_single(frac: &mut [[f64; 3]]) {
    for i in 1..frac.len() {
        for k in 0..3 {
            let delta = frac[i][k] - frac[i - 1][k];
            frac[i][k] -= delta.round();
        }
    }
}

/// 分数坐标折叠到 [0, 1) 并映射到 voxel 索引。
fn voxel_idx(fx: f64, fy: f64, fz: f64, nx: usize, ny: usize, nz: usize) -> (usize, usize, usize) {
    let ix = (fx.rem_euclid(1.0) * nx as f64).floor() as usize % nx;
    let iy = (fy.rem_euclid(1.0) * ny as f64).floor() as usize % ny;
    let iz = (fz.rem_euclid(1.0) * nz as f64).floor() as usize % nz;
    (ix, iy, iz)
}

/// 构造时间平均帧（供 cube 文件头使用）。
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

/// Calculate a spatial jump-distance distribution from a trajectory.
///
/// Each voxel accumulates the number of jump events (displacements exceeding
/// `params.threshold` over `params.tau` frames) recorded at that location.
/// The displacement uses the minimum-image convention via fractional-coordinate
/// unwrapping, making the result correct for periodic and NPT trajectories.
///
/// Returns `None` if:
/// - no frame with a periodic cell is found,
/// - the trajectory has fewer than `tau + 1` frames, or
/// - no atoms match the element filter.
pub fn calc_cube_jump(
    traj: &Trajectory,
    params: &CubeJumpParams,
) -> Option<CubeJumpResult> {
    let n_frames = traj.frames.len();
    if n_frames < params.tau + 1 { return None; }

    let ref_frame = traj.frames.iter().find(|f| f.cell.is_some())?;
    let ref_cell = ref_frame.cell.as_ref().unwrap();

    let (nx, ny, nz) = (params.nx, params.ny, params.nz);
    let threshold2 = params.threshold * params.threshold;

    // 确定要处理的原子索引（以 ref_frame 为准）
    let selected_atoms: Vec<usize> = ref_frame
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| match &params.elements {
            Some(elems) => elems.contains(&a.element),
            None => true,
        })
        .map(|(i, _)| i)
        .collect();

    if selected_atoms.is_empty() { return None; }
    let n_atoms = selected_atoms.len();

    // 预提取所有帧的分数坐标和盒子矩阵转置（NPT：各帧用自身 cell）
    let frac_all: Vec<Vec<[f64; 3]>> = traj.frames.iter().map(|frame| {
        let cell = frame.cell.as_ref().unwrap_or(ref_cell);
        frame.atoms.iter().map(|a| {
            let f = cell.cartesian_to_fractional(a.position);
            [f.x, f.y, f.z]
        }).collect()
    }).collect();

    let mat_t_all: Vec<Matrix3<f64>> = traj.frames.iter().map(|f| {
        f.cell.as_ref().unwrap_or(ref_cell).matrix.transpose()
    }).collect();

    // 以原子为粒度并行，各原子独立 unwrap + 统计
    let (data, n_jumps) = selected_atoms
        .par_iter()
        .map(|&atom_idx| {
            // 变长轨迹安全检查（正常 MD 轨迹不会触发）
            if frac_all.iter().any(|f| atom_idx >= f.len()) {
                return (Array3::<f64>::zeros((nx, ny, nz)), 0usize);
            }

            let mut frac_atom: Vec<[f64; 3]> = frac_all.iter()
                .map(|f| f[atom_idx])
                .collect();

            unwrap_single(&mut frac_atom);

            let mut local = Array3::<f64>::zeros((nx, ny, nz));
            let mut local_jumps = 0usize;

            for t in 0..(n_frames - params.tau) {
                let f0 = frac_atom[t];
                let f1 = frac_atom[t + params.tau];

                // unwrapped Δfrac → Cartesian（NPT：用两帧矩阵的平均，与 msd.rs 一致）
                let dfrac = Vector3::new(f1[0] - f0[0], f1[1] - f0[1], f1[2] - f0[2]);
                let avg_mat_t = (mat_t_all[t] + mat_t_all[t + params.tau]) * 0.5;

                if (avg_mat_t * dfrac).norm_squared() < threshold2 { continue; }

                // 记录位置使用 wrapped 分数坐标
                let (fx, fy, fz) = match params.record_at {
                    JumpPosition::Start    => (f0[0], f0[1], f0[2]),
                    JumpPosition::End      => (f1[0], f1[1], f1[2]),
                    JumpPosition::Midpoint => (
                        (f0[0] + f1[0]) * 0.5,
                        (f0[1] + f1[1]) * 0.5,
                        (f0[2] + f1[2]) * 0.5,
                    ),
                };

                let (ix, iy, iz) = voxel_idx(fx, fy, fz, nx, ny, nz);
                local[[ix, iy, iz]] += 1.0;
                local_jumps += 1;
            }
            (local, local_jumps)
        })
        .reduce(
            || (Array3::<f64>::zeros((nx, ny, nz)), 0),
            |(mut a, ja), (b, jb)| { a += &b; (a, ja + jb) },
        );

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

    Some(CubeJumpResult { cube, n_frames, n_atoms, n_jumps, params: params.clone() })
}

// ─── 测试 ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use nexflux_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    fn make_cell() -> Cell {
        Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap()
    }

    /// 构造含两帧、单原子的轨迹（10Å 立方盒子）。
    fn two_frame_traj(pos0: (f64,f64,f64), pos1: (f64,f64,f64), elem: &str) -> Trajectory {
        let cell = make_cell();
        let mut traj = Trajectory::new();
        for pos in [pos0, pos1] {
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            f.add_atom(Atom::new(elem, Vector3::new(pos.0, pos.1, pos.2)));
            traj.frames.push(f);
        }
        traj
    }

    #[test]
    fn test_no_cell_returns_none() {
        let mut frame = Frame::default();
        frame.add_atom(Atom::new("Li", Vector3::zeros()));
        let mut traj = Trajectory::new();
        traj.frames.push(frame.clone());
        traj.frames.push(frame);
        let params = CubeJumpParams { tau: 1, ..Default::default() };
        assert!(calc_cube_jump(&traj, &params).is_none());
    }

    #[test]
    fn test_too_few_frames_returns_none() {
        let cell = make_cell();
        let mut f = Frame::with_cell(cell, [true; 3]);
        f.add_atom(Atom::new("Li", Vector3::zeros()));
        let mut traj = Trajectory::new();
        traj.frames.push(f);
        let params = CubeJumpParams { tau: 1, ..Default::default() };
        assert!(calc_cube_jump(&traj, &params).is_none());
    }

    #[test]
    fn test_no_jump_below_threshold() {
        // 原子移动 0.5Å，threshold=1.0 → 无跳跃事件
        let traj = two_frame_traj((5.0, 5.0, 5.0), (5.5, 5.0, 5.0), "Li");
        let params = CubeJumpParams {
            nx: 10, ny: 10, nz: 10,
            tau: 1, threshold: 1.0,
            ..Default::default()
        };
        let res = calc_cube_jump(&traj, &params).unwrap();
        assert_eq!(res.n_jumps, 0);
        assert_eq!(res.cube.data.sum(), 0.0);
    }

    #[test]
    fn test_jump_recorded_at_start() {
        // 原子 (2.5,5.5,5.5) → (5.5,5.5,5.5)，距离 3Å > 1.0
        // 起始 frac (0.25,0.55,0.55) → voxel (2,5,5)（中心位置，远离 voxel 边界）
        let traj = two_frame_traj((2.5, 5.5, 5.5), (5.5, 5.5, 5.5), "Li");
        let params = CubeJumpParams {
            nx: 10, ny: 10, nz: 10,
            tau: 1, threshold: 1.0,
            record_at: JumpPosition::Start,
            ..Default::default()
        };
        let res = calc_cube_jump(&traj, &params).unwrap();
        assert_eq!(res.n_jumps, 1);
        assert_eq!(res.cube.data[[2, 5, 5]], 1.0);
    }

    #[test]
    fn test_jump_recorded_at_end() {
        // JumpPosition::End → 终点 frac (0.55,0.55,0.55) → voxel (5,5,5)
        let traj = two_frame_traj((2.5, 5.5, 5.5), (5.5, 5.5, 5.5), "Li");
        let params = CubeJumpParams {
            nx: 10, ny: 10, nz: 10,
            tau: 1, threshold: 1.0,
            record_at: JumpPosition::End,
            ..Default::default()
        };
        let res = calc_cube_jump(&traj, &params).unwrap();
        assert_eq!(res.n_jumps, 1);
        assert_eq!(res.cube.data[[5, 5, 5]], 1.0);
    }

    #[test]
    fn test_jump_recorded_at_midpoint() {
        // (3.5,5.5,5.5) → (7.5,5.5,5.5)，位移 4Å（< 半盒 5Å，不触发 unwrap）
        // 中点 unwrapped frac ((0.35+0.75)/2=0.55, 0.55, 0.55) → voxel (5,5,5)
        let traj = two_frame_traj((3.5, 5.5, 5.5), (7.5, 5.5, 5.5), "Li");
        let params = CubeJumpParams {
            nx: 10, ny: 10, nz: 10,
            tau: 1, threshold: 1.0,
            record_at: JumpPosition::Midpoint,
            ..Default::default()
        };
        let res = calc_cube_jump(&traj, &params).unwrap();
        assert_eq!(res.n_jumps, 1);
        assert_eq!(res.cube.data[[5, 5, 5]], 1.0);
    }

    #[test]
    fn test_pbc_unwrap_gives_min_image_displacement() {
        // 原子 (0.5,5,5) → (9.5,5,5)：直接距离 9Å，最小像距离 1Å（跨越 -x 边界）
        // threshold=2.0：若未正确 unwrap 则 9Å 被记录；正确 unwrap 后 1Å < 2Å 不记录
        let traj = two_frame_traj((0.5, 5.0, 5.0), (9.5, 5.0, 5.0), "Li");
        let params = CubeJumpParams {
            nx: 10, ny: 10, nz: 10,
            tau: 1, threshold: 2.0,
            ..Default::default()
        };
        let res = calc_cube_jump(&traj, &params).unwrap();
        assert_eq!(res.n_jumps, 0, "PBC unwrap 应给出 1Å 位移，而非直接距离 9Å");
    }

    #[test]
    fn test_element_filter() {
        // Li 原子大幅跳跃，O 原子不动；分别用不同 filter 验证
        let cell = make_cell();
        let mut traj = Trajectory::new();
        let mut f0 = Frame::with_cell(cell.clone(), [true; 3]);
        f0.add_atom(Atom::new("Li", Vector3::new(1.0, 5.0, 5.0)));
        f0.add_atom(Atom::new("O",  Vector3::new(5.0, 5.0, 5.0)));
        let mut f1 = Frame::with_cell(cell, [true; 3]);
        f1.add_atom(Atom::new("Li", Vector3::new(8.0, 5.0, 5.0))); // 7Å
        f1.add_atom(Atom::new("O",  Vector3::new(5.0, 5.0, 5.0))); // 不动
        traj.frames.push(f0);
        traj.frames.push(f1);

        let base = CubeJumpParams {
            nx: 10, ny: 10, nz: 10, tau: 1, threshold: 1.0, ..Default::default()
        };

        // 只统计 O → 无跳跃
        let res_o = calc_cube_jump(&traj, &CubeJumpParams {
            elements: Some(vec!["O".to_string()]), ..base.clone()
        }).unwrap();
        assert_eq!(res_o.n_jumps, 0);
        assert_eq!(res_o.n_atoms, 1);

        // 只统计 Li → 1 次跳跃
        let res_li = calc_cube_jump(&traj, &CubeJumpParams {
            elements: Some(vec!["Li".to_string()]), ..base
        }).unwrap();
        assert_eq!(res_li.n_jumps, 1);
    }

    #[test]
    fn test_multi_window_accumulates() {
        // 3 帧：Li 在 (2.5,5.5,5.5), (5.5,5.5,5.5), (8.5,5.5,5.5)
        // tau=1 → 2 个窗口均超过阈值（3Å > 1.0Å）
        // 中心位置确保 voxel 映射远离边界（frac 0.25/0.55/0.85 → voxel 2/5/8）
        let cell = make_cell();
        let mut traj = Trajectory::new();
        for x in [2.5f64, 5.5, 8.5] {
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            f.add_atom(Atom::new("Li", Vector3::new(x, 5.5, 5.5)));
            traj.frames.push(f);
        }
        let params = CubeJumpParams {
            nx: 10, ny: 10, nz: 10, tau: 1, threshold: 1.0, ..Default::default()
        };
        let res = calc_cube_jump(&traj, &params).unwrap();
        assert_eq!(res.n_jumps, 2);
        assert_eq!(res.cube.data[[2, 5, 5]], 1.0); // 第一个窗口的起点 voxel
        assert_eq!(res.cube.data[[5, 5, 5]], 1.0); // 第二个窗口的起点 voxel
        assert_eq!(res.n_frames, 3);
    }

    #[test]
    fn test_grid_shape() {
        let traj = two_frame_traj((2.0, 5.0, 5.0), (8.0, 5.0, 5.0), "Li");
        let params = CubeJumpParams {
            nx: 4, ny: 5, nz: 6, tau: 1, threshold: 1.0, ..Default::default()
        };
        let res = calc_cube_jump(&traj, &params).unwrap();
        assert_eq!(res.cube.shape(), (4, 5, 6));
    }
}
