//! Spatial density distribution (Gaussian cube format).
//!
//! Divides the simulation box into an nx×ny×nz grid and accumulates the time-averaged
//! distribution of selected atoms in each voxel.  Three modes are supported:
//!   - `Density`  — atomic number density \[atoms/Å³\]
//!   - `Velocity` — average speed |v| per voxel \[Å/fs\] (requires `frame.velocities`)
//!   - `Force`    — average force magnitude |f| per voxel \[eV/Å\] (requires `frame.forces`)
//!
//! Output is a [`CubeData`] that can be written directly via `ferro_io::write_cube`
//! to produce a Gaussian cube file for visualisation in VESTA, VMD, etc.
//!
//! **CLI parameters** (for future ferro-cli integration):
//!   `--elements Li,Na`              — include only specified elements (default: all)
//!   `--grid 50 50 50`               — grid dimensions (default: 50×50×50)
//!   `--mode density|velocity|force` — accumulation mode (default: density)
//!
//! Parallelism: per-frame `par_iter`; each frame independently produces (count, value_sum) arrays, then reduced.

use crate::check;
use ferro_core::{CubeData, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use ndarray::Array3;
use rayon::prelude::*;
use super::util::build_avg_frame;

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Spatial distribution mode.
#[derive(Debug, Clone, PartialEq)]
pub enum CubeMode {
    /// Time-averaged number density \[atoms/Å³\]
    Density,
    /// Time-averaged speed |v| per voxel \[Å/fs\] (requires `frame.velocities`)
    Velocity,
    /// Time-averaged force magnitude |f| per voxel \[eV/Å\] (requires `frame.forces`)
    Force,
}

/// Parameters for spatial cube density calculation.
///
/// CLI mapping:
/// - `nx/ny/nz` ← `--grid nx ny nz`
/// - `elements` ← `--elements Li,Na`
/// - `mode`     ← `--mode density|velocity|force`
#[derive(Debug, Clone)]
pub struct CubeDensityParams {
    /// Grid divisions along a axis (CLI: `--grid nx ...`)
    pub nx: usize,
    /// Grid divisions along b axis
    pub ny: usize,
    /// Grid divisions along c axis
    pub nz: usize,
    /// Elements to include (`None` = all atoms; CLI: `--elements Li,Na`)
    pub elements: Option<Vec<String>>,
    /// Quantity to accumulate on the grid
    pub mode: CubeMode,
}

impl CubeDensityParams {
    /// Value ranges that do not depend on the trajectory; the CLI calls this
    /// before reading the first file, and the `calc_*` entry calls it again.
    pub fn validate(&self) -> ferro_core::Result<()> {
        check::at_least_one("nx", self.nx)?;
        check::at_least_one("ny", self.ny)?;
        check::at_least_one("nz", self.nz)
    }
}

impl Default for CubeDensityParams {
    fn default() -> Self {
        Self { nx: 50, ny: 50, nz: 50, elements: None, mode: CubeMode::Density }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Result of a spatial cube density calculation.
///
/// `cube` contains the 3-D grid data and the time-averaged atomic structure.
/// Pass it directly to `ferro_io::write_cube` to produce a Gaussian cube file.
pub struct CubeDensityResult {
    /// 3-D voxel data + time-averaged structure
    pub cube: CubeData,
    /// Number of frames that contributed to the grid
    pub n_frames: usize,
    /// Number of selected atoms per frame (counted from first valid frame)
    pub n_atoms: usize,
    pub params: CubeDensityParams,
}

// ─── cube 参考结构（特例）─────────────────────────────────────────────────────

/// cube 文件里的参考结构：时间平均帧，再把周期方向上越出晶胞的原子折回 `[0, 1)`。
///
/// **特例：只有 cube 产物这样做，与 Ferro 其余各处不同**（审查 M5，用户裁定
/// 2026-10-02）。别处一律保留 reader 读到的原始坐标、不平移也不折回 —— 读入时折回
/// 会给 NPT 的 MSD 带来扰动（TOR 展开下额外折回留下 ΔL 量级残差），按帧减 lo 更会
/// 引入虚假整体漂移，见 `dev/issues.md`。cube 不一样：网格恒铺在 `[0, L)`（原点 0），
/// 体素归属也按折回后的分数坐标算，参考结构若不折回，LAMMPS 盒子 `lo ≠ 0` 时会有原子
/// 画在网格外（`tests/70Z30P00A_NVT_5` 曾有 659 个）。密度数值不受影响。
///
/// `cube_density`、`cube_radius`、`cube_jump` 三处共用。
pub(super) fn cube_reference_frame(traj: &Trajectory) -> Frame {
    let mut frame = build_avg_frame(traj);
    frame.wrap_all();
    frame
}

// ─── 内部辅助 ────────────────────────────────────────────────────────────────

/// Convert fractional coordinate to voxel indices, wrapping periodically.
fn voxel_idx(frac: Vector3<f64>, nx: usize, ny: usize, nz: usize) -> (usize, usize, usize) {
    let ix = (frac.x.rem_euclid(1.0) * nx as f64).floor() as usize % nx;
    let iy = (frac.y.rem_euclid(1.0) * ny as f64).floor() as usize % ny;
    let iz = (frac.z.rem_euclid(1.0) * nz as f64).floor() as usize % nz;
    (ix, iy, iz)
}

/// Process one frame into (count, value_sum) accumulator arrays.
///
/// Returns `None` if the frame has no cell or lacks the required data for the chosen mode.
fn process_frame(
    frame: &Frame,
    params: &CubeDensityParams,
    nx: usize,
    ny: usize,
    nz: usize,
) -> Option<(Array3<f64>, Array3<f64>)> {
    let cell = frame.cell.as_ref()?;
    match params.mode {
        CubeMode::Velocity if frame.velocities.is_none() => return None,
        CubeMode::Force if frame.forces.is_none() => return None,
        _ => {}
    }

    let mut count = Array3::<f64>::zeros((nx, ny, nz));
    let mut value_sum = Array3::<f64>::zeros((nx, ny, nz));

    for (i, atom) in frame.atoms.iter().enumerate() {
        if let Some(elems) = &params.elements {
            if !elems.contains(&atom.element) { continue; }
        }
        let frac = cell.cartesian_to_fractional(atom.position).ok()?;
        let (ix, iy, iz) = voxel_idx(frac, nx, ny, nz);
        match params.mode {
            CubeMode::Density => {
                count[[ix, iy, iz]] += 1.0;
            }
            CubeMode::Velocity => {
                let v = frame.velocities.as_ref().unwrap()[i];
                value_sum[[ix, iy, iz]] += v.norm();
                count[[ix, iy, iz]] += 1.0;
            }
            CubeMode::Force => {
                let f = frame.forces.as_ref().unwrap()[i];
                value_sum[[ix, iy, iz]] += f.norm();
                count[[ix, iy, iz]] += 1.0;
            }
        }
    }
    Some((count, value_sum))
}


// ─── 主函数 ──────────────────────────────────────────────────────────────────

/// Calculate spatial distribution on a voxel grid from a trajectory.
///
/// Returns `None` if:
/// - no frame with a periodic cell is found, or
/// - no frame has the required velocity/force data for the selected mode.
pub fn calc_cube_density(
    traj: &Trajectory,
    params: &CubeDensityParams,
) -> Option<CubeDensityResult> {
    // 返回 Option，报不出原因；CLI 已在读文件前用 validate 报过错，这里只防 panic
    params.validate().ok()?;
    let (nx, ny, nz) = (params.nx, params.ny, params.nz);

    let ref_frame = traj.frames.iter().find(|f| f.cell.is_some())?;
    let ref_cell = ref_frame.cell.as_ref().unwrap();

    // Parallel over frames
    let results: Vec<(Array3<f64>, Array3<f64>)> = traj
        .frames
        .par_iter()
        .filter_map(|f| process_frame(f, params, nx, ny, nz))
        .collect();

    if results.is_empty() { return None; }

    let n_frames = results.len();

    let (count, value_sum) = results.into_iter().fold(
        (Array3::<f64>::zeros((nx, ny, nz)), Array3::<f64>::zeros((nx, ny, nz))),
        |(mut c, mut v), (dc, dv)| { c += &dc; v += &dv; (c, v) },
    );

    let voxel_vol = ref_cell.volume() / (nx * ny * nz) as f64;

    let data = match params.mode {
        CubeMode::Density => count.mapv(|c| c / (n_frames as f64 * voxel_vol)),
        CubeMode::Velocity | CubeMode::Force => Array3::from_shape_fn(
            (nx, ny, nz),
            |(ix, iy, iz)| {
                let c = count[[ix, iy, iz]];
                if c > 0.0 { value_sum[[ix, iy, iz]] / c } else { 0.0 }
            },
        ),
    };

    let n_atoms = ref_frame
        .atoms
        .iter()
        .filter(|a| match &params.elements {
            Some(elems) => elems.contains(&a.element),
            None => true,
        })
        .count();

    // Spacing matrix: row i = cell.row(i) / ni  [Å]
    let m = ref_cell.matrix;
    let spacing = Matrix3::new(
        m[(0, 0)] / nx as f64, m[(0, 1)] / nx as f64, m[(0, 2)] / nx as f64,
        m[(1, 0)] / ny as f64, m[(1, 1)] / ny as f64, m[(1, 2)] / ny as f64,
        m[(2, 0)] / nz as f64, m[(2, 1)] / nz as f64, m[(2, 2)] / nz as f64,
    );

    let cube = CubeData {
        // 特例：参考结构折回盒内，见 cube_reference_frame
        frame: cube_reference_frame(traj),
        data: data.into_iter().collect(),
        shape: [nx, ny, nz],
        origin: Vector3::zeros(),
        spacing,
    };

    Some(CubeDensityResult { cube, n_frames, n_atoms, params: params.clone() })
}

// ─── 测试 ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    /// 10×10×10 Å 立方盒子，原子在 (1,1,1)
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
        frame.add_atom(Atom::new("H", Vector3::zeros()));
        let mut traj = Trajectory::new();
        traj.frames.push(frame);
        let params = CubeDensityParams { nx: 2, ny: 2, nz: 2, ..Default::default() };
        assert!(calc_cube_density(&traj, &params).is_none());
    }

    #[test]
    fn test_density_single_atom() {
        // 10³ box, 2³ grid → voxel_vol = 125 Å³
        // 格点中心位置：voxel (0,0,0) 中心 = frac (0.25,0.25,0.25) = cart (2.5,2.5,2.5)
        // density[0,0,0] = 1 / (1 frame × 125 Å³) = 0.008
        let traj = make_traj(vec![(2.5, 2.5, 2.5)], vec!["Li"]);
        let params = CubeDensityParams { nx: 2, ny: 2, nz: 2, ..Default::default() };
        let res = calc_cube_density(&traj, &params).unwrap();
        assert!((res.cube.get(0, 0, 0) - 0.008).abs() < 1e-10);
        assert_eq!(res.cube.get(1, 0, 0), 0.0);
        assert_eq!(res.n_frames, 1);
        assert_eq!(res.n_atoms, 1);
    }

    #[test]
    fn test_density_multi_frame_averages() {
        // 2 identical frames → same density as 1 frame (2 counts / 2 frames / vox_vol)
        // 格点中心：voxel (0,0,0) → (2.5,2.5,2.5)
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for _ in 0..2 {
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("Li", Vector3::new(2.5, 2.5, 2.5)));
            traj.frames.push(frame);
        }
        let params = CubeDensityParams { nx: 2, ny: 2, nz: 2, ..Default::default() };
        let res = calc_cube_density(&traj, &params).unwrap();
        assert!((res.cube.get(0, 0, 0) - 0.008).abs() < 1e-10);
        assert_eq!(res.n_frames, 2);
    }

    #[test]
    fn test_density_element_filter() {
        // 格点中心位置：H → voxel (0,0,0) 中心 (2.5,2.5,2.5)；O → voxel (1,1,1) 中心 (7.5,7.5,7.5)
        // 仅统计 H：voxel (1,1,1) 应保持为 0
        let traj = make_traj(vec![(2.5, 2.5, 2.5), (7.5, 7.5, 7.5)], vec!["H", "O"]);
        let params = CubeDensityParams {
            nx: 2, ny: 2, nz: 2,
            elements: Some(vec!["H".to_string()]),
            ..Default::default()
        };
        let res = calc_cube_density(&traj, &params).unwrap();
        assert!((res.cube.get(0, 0, 0) - 0.008).abs() < 1e-10);
        assert_eq!(res.cube.get(1, 1, 1), 0.0);
        assert_eq!(res.n_atoms, 1);
    }

    #[test]
    fn test_periodic_wrap() {
        // 格点中心 (2.5,2.5,2.5) 超出盒子一个周期：frac (1.25,0.25,0.25) → rem_euclid → (0.25,0.25,0.25) → voxel (0,0,0)
        let traj = make_traj(vec![(12.5, 2.5, 2.5)], vec!["Li"]);
        let params = CubeDensityParams { nx: 2, ny: 2, nz: 2, ..Default::default() };
        let res = calc_cube_density(&traj, &params).unwrap();
        assert!((res.cube.get(0, 0, 0) - 0.008).abs() < 1e-10);
    }

    #[test]
    fn test_velocity_mode() {
        // 格点中心：voxel (0,0,0) → (2.5,2.5,2.5)；velocity (2,0,0) → |v|=2
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Li", Vector3::new(2.5, 2.5, 2.5)));
        frame.velocities = Some(vec![Vector3::new(2.0, 0.0, 0.0)]);
        let mut traj = Trajectory::new();
        traj.frames.push(frame);
        let params = CubeDensityParams {
            nx: 2, ny: 2, nz: 2,
            mode: CubeMode::Velocity,
            ..Default::default()
        };
        let res = calc_cube_density(&traj, &params).unwrap();
        assert!((res.cube.get(0, 0, 0) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_velocity_mode_no_velocities_returns_none() {
        let traj = make_traj(vec![(2.5, 2.5, 2.5)], vec!["Li"]);
        let params = CubeDensityParams {
            nx: 2, ny: 2, nz: 2,
            mode: CubeMode::Velocity,
            ..Default::default()
        };
        assert!(calc_cube_density(&traj, &params).is_none());
    }

    #[test]
    fn test_force_mode() {
        // 格点中心：voxel (0,0,0) → (2.5,2.5,2.5)；force (0,3,4) → |f|=5
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("O", Vector3::new(2.5, 2.5, 2.5)));
        frame.forces = Some(vec![Vector3::new(0.0, 3.0, 4.0)]);
        let mut traj = Trajectory::new();
        traj.frames.push(frame);
        let params = CubeDensityParams {
            nx: 2, ny: 2, nz: 2,
            mode: CubeMode::Force,
            ..Default::default()
        };
        let res = calc_cube_density(&traj, &params).unwrap();
        assert!((res.cube.get(0, 0, 0) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_grid_shape() {
        let traj = make_traj(vec![(1.0, 1.0, 1.0)], vec!["Li"]);
        let params = CubeDensityParams { nx: 2, ny: 3, nz: 4, ..Default::default() };
        let res = calc_cube_density(&traj, &params).unwrap();
        assert_eq!(res.cube.shape(), (2, 3, 4));
    }

    #[test]
    fn test_spacing_matrix() {
        // 10Å cubic box, 2×2×2 grid → each voxel step = 5Å along each axis
        let traj = make_traj(vec![(1.0, 1.0, 1.0)], vec!["Li"]);
        let params = CubeDensityParams { nx: 2, ny: 2, nz: 2, ..Default::default() };
        let res = calc_cube_density(&traj, &params).unwrap();
        let s = res.cube.spacing;
        assert!((s[(0, 0)] - 5.0).abs() < 1e-10); // a/2 along x
        assert!((s[(1, 1)] - 5.0).abs() < 1e-10); // b/2 along y
        assert!((s[(2, 2)] - 5.0).abs() < 1e-10); // c/2 along z
    }

    #[test]
    fn test_reference_frame_is_wrapped_but_the_density_is_not_changed() {
        // 模拟 LAMMPS 盒子 lo = 1.8：原子绝对坐标在 [1.8, 11.8)，有的 x > 10
        let outside = make_traj(vec![(10.5, 3.0, 3.0), (5.0, 11.2, 4.0), (2.0, 2.0, 2.0)], vec!["O", "O", "H"]);
        let folded = make_traj(vec![(0.5, 3.0, 3.0), (5.0, 1.2, 4.0), (2.0, 2.0, 2.0)], vec!["O", "O", "H"]);
        let params = CubeDensityParams { nx: 5, ny: 5, nz: 5, ..Default::default() };
        let a = calc_cube_density(&outside, &params).unwrap();
        let b = calc_cube_density(&folded, &params).unwrap();
        // 体素归属本来就按折回后的分数坐标：密度逐位相同
        assert_eq!(a.cube.data, b.cube.data);
        // 参考结构折回网格 [0, 10)，盒内原子不动
        for (x, y) in a.cube.frame.atoms.iter().zip(&b.cube.frame.atoms) {
            assert!((x.position - y.position).norm() < 1e-12, "{:?} vs {:?}", x.position, y.position);
        }
        assert_eq!(a.cube.frame.atoms[2].position, Vector3::new(2.0, 2.0, 2.0));
    }
}
