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

use ferro_core::{CubeData, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use ndarray::Array3;
use rayon::prelude::*;

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
        let frac = cell.cartesian_to_fractional(atom.position);
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

/// Calculate spatial distribution on a voxel grid from a trajectory.
///
/// Returns `None` if:
/// - no frame with a periodic cell is found, or
/// - no frame has the required velocity/force data for the selected mode.
pub fn calc_cube_density(
    traj: &Trajectory,
    params: &CubeDensityParams,
) -> Option<CubeDensityResult> {
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
        frame: build_avg_frame(traj),
        data,
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
        // 1 atom at (1,1,1) → frac (0.1,0.1,0.1) → voxel (0,0,0)
        // density[0,0,0] = 1 / (1 frame × 125 Å³) = 0.008
        let traj = make_traj(vec![(1.0, 1.0, 1.0)], vec!["Li"]);
        let params = CubeDensityParams { nx: 2, ny: 2, nz: 2, ..Default::default() };
        let res = calc_cube_density(&traj, &params).unwrap();
        let d = &res.cube.data;
        assert!((d[[0, 0, 0]] - 0.008).abs() < 1e-10);
        assert_eq!(d[[1, 0, 0]], 0.0);
        assert_eq!(res.n_frames, 1);
        assert_eq!(res.n_atoms, 1);
    }

    #[test]
    fn test_density_multi_frame_averages() {
        // 2 identical frames → same density as 1 frame (2 counts / 2 frames / vox_vol)
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for _ in 0..2 {
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("Li", Vector3::new(1.0, 1.0, 1.0)));
            traj.frames.push(frame);
        }
        let params = CubeDensityParams { nx: 2, ny: 2, nz: 2, ..Default::default() };
        let res = calc_cube_density(&traj, &params).unwrap();
        assert!((res.cube.data[[0, 0, 0]] - 0.008).abs() < 1e-10);
        assert_eq!(res.n_frames, 2);
    }

    #[test]
    fn test_density_element_filter() {
        // H at (1,1,1) → voxel (0,0,0), O at (6,6,6) → voxel (1,1,1)
        // Filter to H only: O voxel should remain 0
        let traj = make_traj(vec![(1.0, 1.0, 1.0), (6.0, 6.0, 6.0)], vec!["H", "O"]);
        let params = CubeDensityParams {
            nx: 2, ny: 2, nz: 2,
            elements: Some(vec!["H".to_string()]),
            ..Default::default()
        };
        let res = calc_cube_density(&traj, &params).unwrap();
        let d = &res.cube.data;
        assert!((d[[0, 0, 0]] - 0.008).abs() < 1e-10);
        assert_eq!(d[[1, 1, 1]], 0.0);
        assert_eq!(res.n_atoms, 1);
    }

    #[test]
    fn test_periodic_wrap() {
        // Atom at (11,1,1): frac = (1.1,0.1,0.1) → rem_euclid → (0.1,0.1,0.1) → voxel (0,0,0)
        // Same result as atom at (1,1,1)
        let traj = make_traj(vec![(11.0, 1.0, 1.0)], vec!["Li"]);
        let params = CubeDensityParams { nx: 2, ny: 2, nz: 2, ..Default::default() };
        let res = calc_cube_density(&traj, &params).unwrap();
        assert!((res.cube.data[[0, 0, 0]] - 0.008).abs() < 1e-10);
    }

    #[test]
    fn test_velocity_mode() {
        // Atom at (1,1,1), velocity (2,0,0) → |v|=2 → grid[0,0,0] = 2.0
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Li", Vector3::new(1.0, 1.0, 1.0)));
        frame.velocities = Some(vec![Vector3::new(2.0, 0.0, 0.0)]);
        let mut traj = Trajectory::new();
        traj.frames.push(frame);
        let params = CubeDensityParams {
            nx: 2, ny: 2, nz: 2,
            mode: CubeMode::Velocity,
            ..Default::default()
        };
        let res = calc_cube_density(&traj, &params).unwrap();
        assert!((res.cube.data[[0, 0, 0]] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_velocity_mode_no_velocities_returns_none() {
        // Frame has no velocities → None for Velocity mode
        let traj = make_traj(vec![(1.0, 1.0, 1.0)], vec!["Li"]);
        let params = CubeDensityParams {
            nx: 2, ny: 2, nz: 2,
            mode: CubeMode::Velocity,
            ..Default::default()
        };
        assert!(calc_cube_density(&traj, &params).is_none());
    }

    #[test]
    fn test_force_mode() {
        // Force (0,3,4) → |f| = 5 → grid[0,0,0] = 5.0
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("O", Vector3::new(1.0, 1.0, 1.0)));
        frame.forces = Some(vec![Vector3::new(0.0, 3.0, 4.0)]);
        let mut traj = Trajectory::new();
        traj.frames.push(frame);
        let params = CubeDensityParams {
            nx: 2, ny: 2, nz: 2,
            mode: CubeMode::Force,
            ..Default::default()
        };
        let res = calc_cube_density(&traj, &params).unwrap();
        assert!((res.cube.data[[0, 0, 0]] - 5.0).abs() < 1e-10);
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
}
