//! Spatial density distribution (Gaussian cube format).
//!
//! Divides the simulation box into an nx×ny×nz grid and accumulates the time-averaged
//! distribution of selected atoms in each voxel.  Three modes are supported:
//!   - `Density`  — atomic number density \[atoms/Å³\]
//!   - `Velocity` — average speed |v| per voxel \[Å/fs\] (requires `frame.velocities()`)
//!   - `Force`    — average force magnitude |f| per voxel \[eV/Å\] (requires `frame.forces()`)
//!
//! Output is a [`CubeData`] that can be written directly via `ferro_io::write_cube`
//! to produce a Gaussian cube file for visualisation in VESTA, VMD, etc.
//!
//! **CLI**: `ferro map density|velocity|force` picks the mode; `--nx --ny --nz` set the
//! grid (default 50 each) and `--elements Li,Na` the atoms (default: all).
//!
//! Parallelism: per-frame `par_iter`; each frame independently produces (count, value_sum) arrays, then reduced.

use crate::check;
use ferro_core::error::ChemError;
use ferro_core::{CubeData, Frame, Trajectory};
use nalgebra::Vector3;
use ndarray::Array3;
use super::util::{cube_from_grid, par_map_in_order, select_atoms};

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Spatial distribution mode.
#[derive(Debug, Clone, PartialEq)]
pub enum CubeMode {
    /// Time-averaged number density \[atoms/Å³\]
    Density,
    /// Time-averaged speed |v| per voxel \[Å/fs\] (requires `frame.velocities()`)
    Velocity,
    /// Time-averaged force magnitude |f| per voxel \[eV/Å\] (requires `frame.forces()`)
    Force,
}

/// Parameters for spatial cube density calculation.
///
/// CLI mapping:
/// - `nx/ny/nz` ← `--nx --ny --nz`
/// - `elements` ← `--elements Li,Na`
/// - `mode`     ← the subcommand, `map density|velocity|force`
#[derive(Debug, Clone)]
pub struct CubeDensityParams {
    /// Grid divisions along a axis (CLI: `--nx`)
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
        CubeMode::Velocity if frame.velocities().is_none() => return None,
        CubeMode::Force if frame.forces().is_none() => return None,
        _ => {}
    }

    let mut count = Array3::<f64>::zeros((nx, ny, nz));
    let mut value_sum = Array3::<f64>::zeros((nx, ny, nz));

    for (i, atom) in frame.atoms().enumerate() {
        if let Some(elems) = &params.elements {
            if !elems.iter().any(|e| e == atom.element) { continue; }
        }
        let frac = cell.cartesian_to_fractional(atom.position).ok()?;
        let (ix, iy, iz) = voxel_idx(frac, nx, ny, nz);
        match params.mode {
            CubeMode::Density => {
                count[[ix, iy, iz]] += 1.0;
            }
            CubeMode::Velocity => {
                let v = frame.velocities().unwrap()[i];
                value_sum[[ix, iy, iz]] += v.norm();
                count[[ix, iy, iz]] += 1.0;
            }
            CubeMode::Force => {
                let f = frame.forces().unwrap()[i];
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
/// Errors if a parameter is out of range, a frame's cell is not 3-D periodic, no frame
/// has a periodic cell, or no frame carries the velocities / forces the mode needs.
pub fn calc_cube_density(
    traj: &Trajectory,
    params: &CubeDensityParams,
) -> ferro_core::Result<CubeDensityResult> {
    params.validate()?;
    check::periodic_cells(traj, None)?;
    let (nx, ny, nz) = (params.nx, params.ny, params.nz);

    let ref_frame = traj.frames.iter().find(|f| f.cell.is_some())
        .ok_or_else(|| ChemError::ValidationError("no frame has a periodic cell".into()))?;
    let ref_cell = ref_frame.cell.as_ref().unwrap();
    let n_atoms = select_atoms(ref_frame, params.elements.as_deref())?.len();

    // Parallel over frames, summed in frame order
    let mut count = Array3::<f64>::zeros((nx, ny, nz));
    let mut value_sum = Array3::<f64>::zeros((nx, ny, nz));
    let mut n_frames = 0usize;
    par_map_in_order(&traj.frames, |f| process_frame(f, params, nx, ny, nz), |r| {
        if let Some((dc, dv)) = r {
            count += &dc;
            value_sum += &dv;
            n_frames += 1;
        }
    });

    if n_frames == 0 {
        let what = if matches!(params.mode, CubeMode::Velocity) { "velocities" } else { "forces" };
        return Err(ChemError::ValidationError(format!(
            "no frame with a periodic cell carries {what}, which mode {:?} needs", params.mode
        )));
    }

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

    let cube = cube_from_grid(traj, ref_cell, data);

    Ok(CubeDensityResult { cube, n_frames, n_atoms, params: params.clone() })
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
    fn test_no_cell_is_an_error() {
        let mut frame = Frame::default();
        frame.add_atom(Atom::new("H", Vector3::zeros()));
        let mut traj = Trajectory::new();
        traj.frames.push(frame);
        let params = CubeDensityParams { nx: 2, ny: 2, nz: 2, ..Default::default() };
        let err = calc_cube_density(&traj, &params).err().unwrap().to_string();
        assert!(err.contains("no frame has a periodic cell"), "{err}");
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

    // 审查 B-2：cube 规范里第 i 点在 origin + i·spacing。体素平均值属于体素中心，
    // origin 写 0 时整张网格偏半个体素（这里 1.25 Å/轴）
    #[test]
    fn test_grid_point_sits_on_voxel_centre() {
        let atom = Vector3::new(3.75, 6.25, 1.25); // 4³ 网格、spacing 2.5 时体素 (1,2,0) 的中心
        let traj = make_traj(vec![(atom.x, atom.y, atom.z)], vec!["Li"]);
        let params = CubeDensityParams { nx: 4, ny: 4, nz: 4, ..Default::default() };
        let cube = calc_cube_density(&traj, &params).unwrap().cube;
        assert!(cube.get(1, 2, 0) > 0.0, "原子应落在体素 (1,2,0)");
        let at = cube.origin + cube.spacing.transpose() * Vector3::new(1.0, 2.0, 0.0);
        assert!((at - atom).norm() < 1e-10, "峰值格点应在原子上，实得 {at:?}");
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
        frame.set_velocities(Some(vec![Vector3::new(2.0, 0.0, 0.0)])).unwrap();
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
    fn test_velocity_mode_no_velocities_is_an_error() {
        let traj = make_traj(vec![(2.5, 2.5, 2.5)], vec!["Li"]);
        let params = CubeDensityParams {
            nx: 2, ny: 2, nz: 2,
            mode: CubeMode::Velocity,
            ..Default::default()
        };
        let err = calc_cube_density(&traj, &params).err().unwrap().to_string();
        assert!(err.contains("carries velocities"), "{err}");
    }

    #[test]
    fn test_force_mode() {
        // 格点中心：voxel (0,0,0) → (2.5,2.5,2.5)；force (0,3,4) → |f|=5
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("O", Vector3::new(2.5, 2.5, 2.5)));
        frame.set_forces(Some(vec![Vector3::new(0.0, 3.0, 4.0)])).unwrap();
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
        for (x, y) in a.cube.frame.atoms().zip(b.cube.frame.atoms()) {
            assert!((x.position - y.position).norm() < 1e-12, "{:?} vs {:?}", x.position, y.position);
        }
        assert_eq!(a.cube.frame.positions()[2], Vector3::new(2.0, 2.0, 2.0));
    }

    /// An element filter that matches nothing is an error, not a grid of zeros
    #[test]
    fn an_empty_selection_is_an_error() {
        let traj = make_traj(vec![(5.0, 5.0, 5.0)], vec!["Li"]);
        let params = CubeDensityParams { nx: 4, ny: 4, nz: 4, elements: Some(vec!["Na".into()]), ..Default::default() };
        let err = calc_cube_density(&traj, &params).err().unwrap().to_string();
        assert!(err.contains("no atom matches the element filter"), "{err}");
    }
}
