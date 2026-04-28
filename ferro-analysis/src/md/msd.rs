//! Mean squared displacement (MSD) calculation and output.
//!
//! Workflow: `calc_msd` → `write_msd`.
//! Algorithm follows code1/msd.c (`EstimateMSD`):
//!   1. Convert each frame's Cartesian coordinates to fractional coordinates.
//!   2. Unwrap: detect fractional-coordinate jumps (|Δ| > 0.5) and correct cross-boundary displacements.
//!   3. Time-shift averaging: step = shift, window = tau.
//!   4. NPT support: total MSD uses the average of the origin- and endpoint-frame cell matrices.
//!
//! Parallelism: per time-origin par_iter; each origin computed independently then reduced.

use rayon::prelude::*;
use std::collections::BTreeSet;
use std::io::{BufWriter, Write};
use ferro_core::Trajectory;
use super::gr::VERSION;

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Parameters for MSD calculation.
#[derive(Debug, Clone)]
pub struct MsdParams {
    /// Lag window size in frames (`None` = use all frames)
    pub tau: Option<usize>,
    /// Time shift between origins in frames (default: 1)
    pub shift: usize,
    /// Time step per frame \[fs\] (default: 1.0)
    pub dt: f64,
    /// Elements to include (`None` = all atoms)
    pub elements: Option<Vec<String>>,
}

impl Default for MsdParams {
    fn default() -> Self {
        MsdParams { tau: None, shift: 1, dt: 1.0, elements: None }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Result of an MSD calculation.
///
/// For periodic trajectories the directional columns `msd_a/b/c` represent
/// displacement along the three crystal axes.  For non-periodic trajectories
/// they correspond to the Cartesian x/y/z axes.
#[derive(Debug, Clone)]
pub struct MsdResult {
    /// Time axis \[fs\]; `time[i] = i * dt`
    pub time: Vec<f64>,
    /// Total MSD \[Å²\]
    pub msd: Vec<f64>,
    /// Directional MSD along the a-axis (or x) \[Å²\]
    pub msd_a: Vec<f64>,
    /// Directional MSD along the b-axis (or y) \[Å²\]
    pub msd_b: Vec<f64>,
    /// Directional MSD along the c-axis (or z) \[Å²\]
    pub msd_c: Vec<f64>,
    /// Number of atoms included in the calculation
    pub n_atoms: usize,
    /// Number of time origins averaged
    pub n_origins: usize,
    pub params: MsdParams,
    /// Elements included (sorted alphabetically)
    pub elements: Vec<String>,
}

// ─── 内部辅助 ─────────────────────────────────────────────────────────────────

/// Unwrap fractional coordinates in-place to remove periodic-boundary jumps.
///
/// Checks the fractional-coordinate difference between adjacent frames:
/// |Δ| > 0.5 indicates a boundary crossing; corrected by subtracting round(Δ).
pub(super) fn unwrap_frac(frac: &mut [Vec<[f64; 3]>]) {
    let n_steps = frac.len();
    if n_steps < 2 { return; }
    let n_atoms = frac[0].len();
    for j in 0..n_atoms {
        for i in 1..n_steps {
            for k in 0..3 {
                let delta = frac[i][j][k] - frac[i - 1][j][k];
                frac[i][j][k] -= delta.round();
            }
        }
    }
}

// ─── 计算 ────────────────────────────────────────────────────────────────────

/// Compute mean squared displacement with time-shift averaging.
///
/// The algorithm matches code1/msd.c `EstimateMSD`:
/// 1. Convert atom Cartesian coordinates to fractional (periodic case only).
/// 2. Unwrap fractional coordinates across periodic boundaries.
/// 3. Average over all time origins spaced `params.shift` frames apart,
///    computed in parallel (one task per origin).
///
/// Returns `None` if:
/// - The trajectory has fewer than 2 frames
/// - No atoms match the element filter
/// - `tau` exceeds the trajectory length
/// - Any frame is missing a cell (periodic path only)
///
/// # NPT handling
/// Total MSD uses the average of the origin- and endpoint-frame cell matrices
/// to convert fractional displacements to Cartesian. Directional MSD uses the
/// endpoint cell parameter (same simplified approximation as code1/msd.c).
pub fn calc_msd(traj: &Trajectory, params: &MsdParams) -> Option<MsdResult> {
    let n_steps = traj.n_frames();
    if n_steps < 2 { return None; }

    // 确定参与计算的原子下标（按第一帧筛选元素）
    let ref_frame = traj.first()?;
    let atom_indices: Vec<usize> = ref_frame.atoms.iter().enumerate()
        .filter(|(_, a)| match &params.elements {
            Some(elems) => elems.contains(&a.element),
            None => true,
        })
        .map(|(i, _)| i)
        .collect();
    if atom_indices.is_empty() { return None; }
    let n_atoms = atom_indices.len();

    let tau = params.tau.unwrap_or(n_steps).min(n_steps).max(1);
    let shift = params.shift.max(1);

    // 收集元素列表（用于输出文件头）
    let elements: Vec<String> = {
        let mut set = BTreeSet::new();
        for &i in &atom_indices { set.insert(ref_frame.atoms[i].element.clone()); }
        set.into_iter().collect()
    };

    if ref_frame.cell.is_some() {
        calc_msd_periodic(traj, &atom_indices, n_atoms, tau, shift, params, elements)
    } else {
        calc_msd_nonperiodic(traj, &atom_indices, n_atoms, tau, shift, params, elements)
    }
}

/// MSD for periodic boundary conditions (fractional-coordinate unwrapping path, parallelised per time origin).
fn calc_msd_periodic(
    traj: &Trajectory,
    atom_indices: &[usize],
    n_atoms: usize,
    tau: usize,
    shift: usize,
    params: &MsdParams,
    elements: Vec<String>,
) -> Option<MsdResult> {
    let n_steps = traj.n_frames();

    // 先校验所有帧都有 cell
    if traj.frames.iter().any(|f| f.cell.is_none()) { return None; }

    // 构建 frac[step][atom_local] = [fa, fb, fc]（串行，顺序依赖无法并行）
    let mut frac: Vec<Vec<[f64; 3]>> = traj.frames.iter().map(|frame| {
        let cell = frame.cell.as_ref().unwrap();
        atom_indices.iter().map(|&i| {
            let f = cell.cartesian_to_fractional(frame.atoms[i].position);
            [f.x, f.y, f.z]
        }).collect()
    }).collect();

    // Unwrap 分数坐标（顺序依赖，串行）
    unwrap_frac(&mut frac);

    // 枚举所有 origin 的起始帧
    let p_values: Vec<usize> = (0..)
        .map(|k: usize| k * shift)
        .take_while(|&p| p + tau <= n_steps)
        .collect();
    let n_origins = p_values.len();
    if n_origins == 0 { return None; }

    // 并行计算各 origin 的局部累积，每个 origin 产生 Vec<[f64;4]>(tau)
    let accum: Vec<[f64; 4]> = p_values.par_iter()
        .map(|&p| {
            let cell_orig = traj.frames[p].cell.as_ref().unwrap();
            let mut local = vec![[0.0f64; 4]; tau];
            for i in 0..tau {
                let cell_end = traj.frames[p + i].cell.as_ref().unwrap();
                // NPT：取两端盒子矩阵的平均（同 code1 的做法）
                let avg_mat = (cell_end.matrix + cell_orig.matrix) * 0.5;
                let [a_len, b_len, c_len] = cell_end.lengths();
                let mut sum = [0.0f64; 4];
                for j in 0..n_atoms {
                    let dx = frac[p + i][j][0] - frac[p][j][0];
                    let dy = frac[p + i][j][1] - frac[p][j][1];
                    let dz = frac[p + i][j][2] - frac[p][j][2];
                    // 分数位移 → Cartesian（avg_mat 行向量 = a,b,c）
                    let cx = dx*avg_mat[(0,0)] + dy*avg_mat[(1,0)] + dz*avg_mat[(2,0)];
                    let cy = dx*avg_mat[(0,1)] + dy*avg_mat[(1,1)] + dz*avg_mat[(2,1)];
                    let cz = dx*avg_mat[(0,2)] + dy*avg_mat[(1,2)] + dz*avg_mat[(2,2)];
                    sum[0] += cx*cx + cy*cy + cz*cz;
                    // 各轴分量：分数位移 × endpoint 轴长（简化近似，同 code1）
                    sum[1] += dx*dx * a_len*a_len;
                    sum[2] += dy*dy * b_len*b_len;
                    sum[3] += dz*dz * c_len*c_len;
                }
                local[i] = [
                    sum[0] / n_atoms as f64,
                    sum[1] / n_atoms as f64,
                    sum[2] / n_atoms as f64,
                    sum[3] / n_atoms as f64,
                ];
            }
            local
        })
        .reduce(
            || vec![[0.0f64; 4]; tau],
            |mut a, b| {
                for i in 0..tau { for k in 0..4 { a[i][k] += b[i][k]; } }
                a
            },
        );

    Some(build_result(accum, tau, n_origins, n_atoms, elements, params))
}

/// MSD for non-periodic (molecular) systems (Cartesian coordinates directly, parallelised per time origin).
fn calc_msd_nonperiodic(
    traj: &Trajectory,
    atom_indices: &[usize],
    n_atoms: usize,
    tau: usize,
    shift: usize,
    params: &MsdParams,
    elements: Vec<String>,
) -> Option<MsdResult> {
    let n_steps = traj.n_frames();

    // 收集各帧 Cartesian 坐标（非周期不需要 unwrap）
    let cart: Vec<Vec<[f64; 3]>> = traj.frames.iter().map(|frame| {
        atom_indices.iter().map(|&i| {
            let p = &frame.atoms[i].position;
            [p.x, p.y, p.z]
        }).collect()
    }).collect();

    let p_values: Vec<usize> = (0..)
        .map(|k: usize| k * shift)
        .take_while(|&p| p + tau <= n_steps)
        .collect();
    let n_origins = p_values.len();
    if n_origins == 0 { return None; }

    let accum: Vec<[f64; 4]> = p_values.par_iter()
        .map(|&p| {
            let mut local = vec![[0.0f64; 4]; tau];
            for i in 0..tau {
                let mut sum = [0.0f64; 4];
                for j in 0..n_atoms {
                    let dx = cart[p + i][j][0] - cart[p][j][0];
                    let dy = cart[p + i][j][1] - cart[p][j][1];
                    let dz = cart[p + i][j][2] - cart[p][j][2];
                    sum[0] += dx*dx + dy*dy + dz*dz;
                    sum[1] += dx*dx;
                    sum[2] += dy*dy;
                    sum[3] += dz*dz;
                }
                local[i] = [
                    sum[0] / n_atoms as f64,
                    sum[1] / n_atoms as f64,
                    sum[2] / n_atoms as f64,
                    sum[3] / n_atoms as f64,
                ];
            }
            local
        })
        .reduce(
            || vec![[0.0f64; 4]; tau],
            |mut a, b| {
                for i in 0..tau { for k in 0..4 { a[i][k] += b[i][k]; } }
                a
            },
        );

    Some(build_result(accum, tau, n_origins, n_atoms, elements, params))
}

/// Build an `MsdResult` from the parallel-reduction accumulation array.
fn build_result(
    accum: Vec<[f64; 4]>,
    tau: usize,
    n_origins: usize,
    n_atoms: usize,
    elements: Vec<String>,
    params: &MsdParams,
) -> MsdResult {
    let inv = 1.0 / n_origins as f64;
    let time:  Vec<f64> = (0..tau).map(|i| i as f64 * params.dt).collect();
    let msd:   Vec<f64> = (0..tau).map(|i| accum[i][0] * inv).collect();
    let msd_a: Vec<f64> = (0..tau).map(|i| accum[i][1] * inv).collect();
    let msd_b: Vec<f64> = (0..tau).map(|i| accum[i][2] * inv).collect();
    let msd_c: Vec<f64> = (0..tau).map(|i| accum[i][3] * inv).collect();
    MsdResult { time, msd, msd_a, msd_b, msd_c, n_atoms, n_origins, params: params.clone(), elements }
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

/// Write MSD data to a tab-separated text file (`.msd`).
///
/// Columns: `time[fs]`, `msd[Ang^2]`, `msd_a[Ang^2]`, `msd_b[Ang^2]`, `msd_c[Ang^2]`
///
/// For periodic trajectories a/b/c are the crystal-axis directions.
/// For non-periodic trajectories a/b/c correspond to Cartesian x/y/z.
pub fn write_msd(result: &MsdResult, path: &str) -> std::io::Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(path)?);

    writeln!(w, "# ferro v{}", VERSION)?;
    writeln!(w, "# Mean Squared Displacement (MSD)")?;
    writeln!(w, "# {}", "-".repeat(60))?;
    writeln!(w, "# tau     = {} frames", result.time.len())?;
    writeln!(w, "# shift   = {} frames", result.params.shift)?;
    writeln!(w, "# dt      = {} fs", result.params.dt)?;
    writeln!(w, "# atoms   = {}", result.n_atoms)?;
    writeln!(w, "# origins = {}", result.n_origins)?;
    write!(w,   "# elements:")?;
    for e in &result.elements { write!(w, " {}", e)?; }
    writeln!(w)?;
    writeln!(w, "# {}", "-".repeat(60))?;
    writeln!(w, "# time[fs]\tmsd[Ang^2]\tmsd_a[Ang^2]\tmsd_b[Ang^2]\tmsd_c[Ang^2]")?;

    for i in 0..result.time.len() {
        writeln!(w, "{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}",
            result.time[i], result.msd[i],
            result.msd_a[i], result.msd_b[i], result.msd_c[i])?;
    }
    Ok(())
}

// ─── 测试 ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    /// 构建 n 帧的 NVT 轨迹：单个 Fe 原子每帧沿 x 移动 v Å
    fn make_traj_linear(a: f64, v: f64, n: usize) -> Trajectory {
        let cell = Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for i in 0..n {
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("Fe", Vector3::new(i as f64 * v, 0.0, 0.0)));
            traj.add_frame(frame);
        }
        traj
    }

    /// 构建静态轨迹：n 帧全相同（3×3×3 Fe sc）
    fn make_traj_static(n: usize) -> Trajectory {
        let a = 2.87_f64;
        let side = 3.0 * a;
        let cell = Cell::from_lengths_angles(side, side, side, 90.0, 90.0, 90.0).unwrap();
        let mut ref_frame = Frame::with_cell(cell, [true; 3]);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    ref_frame.add_atom(Atom::new(
                        "Fe",
                        Vector3::new(i as f64 * a, j as f64 * a, k as f64 * a),
                    ));
                }
            }
        }
        let mut traj = Trajectory::new();
        for _ in 0..n { traj.add_frame(ref_frame.clone()); }
        traj
    }

    #[test]
    fn test_msd_static_is_zero() {
        // 所有帧完全相同 → 所有 lag 的 MSD 应为 0
        let traj = make_traj_static(8);
        let result = calc_msd(&traj, &MsdParams::default()).unwrap();
        for &v in &result.msd {
            assert!(v.abs() < 1e-10, "static MSD should be 0, got {}", v);
        }
    }

    #[test]
    fn test_msd_linear_motion() {
        // 单原子以 v=0.3 Å/step 沿 x 匀速运动，MSD[lag] = (v * lag)²
        let a = 20.0;
        let v = 0.3;
        let n = 6;
        let traj = make_traj_linear(a, v, n);
        let result = calc_msd(&traj, &MsdParams {
            tau: Some(n), shift: 1, dt: 1.0, elements: None,
        }).unwrap();

        for (lag, &msd_val) in result.msd.iter().enumerate() {
            let expected = (v * lag as f64).powi(2);
            assert!(
                (msd_val - expected).abs() < 1e-8,
                "lag {}: expected {:.6e}, got {:.6e}", lag, expected, msd_val,
            );
        }
        // a 方向应等于 total（运动沿 x=a 轴）
        for (lag, &msd_a) in result.msd_a.iter().enumerate() {
            let expected = (v * lag as f64).powi(2);
            assert!((msd_a - expected).abs() < 1e-8, "msd_a lag {}: {}", lag, msd_a);
        }
    }

    #[test]
    fn test_msd_unwrap_across_boundary() {
        // 原子从接近边界处出发，以 0.5 Å/step 运动，会跨越周期边界
        // unwrap 后 MSD 应等于 (v * lag)²
        let a = 5.0;
        let v = 0.5;
        let x0 = 4.8_f64;
        let n = 5;
        let cell = Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for i in 0..n {
            let x_raw = x0 + i as f64 * v;
            let x_wrapped = x_raw - (x_raw / a).floor() * a;
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("Fe", Vector3::new(x_wrapped, 0.0, 0.0)));
            traj.add_frame(frame);
        }
        let result = calc_msd(&traj, &MsdParams {
            tau: Some(n), shift: 1, dt: 1.0, elements: None,
        }).unwrap();

        for (lag, &msd_val) in result.msd.iter().enumerate() {
            let expected = (v * lag as f64).powi(2);
            assert!(
                (msd_val - expected).abs() < 1e-8,
                "unwrap lag {}: expected {:.6e}, got {:.6e}", lag, expected, msd_val,
            );
        }
    }

    #[test]
    fn test_msd_element_filter() {
        // 轨迹含 Fe 和 Li，只计算 Li 的 MSD
        let a = 10.0;
        let v_li = 0.4;
        let cell = Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for i in 0..5 {
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("Fe", Vector3::new(1.0, 1.0, 1.0)));
            frame.add_atom(Atom::new("Li", Vector3::new(i as f64 * v_li, 0.0, 0.0)));
            traj.add_frame(frame);
        }

        let result = calc_msd(&traj, &MsdParams {
            tau: Some(5), shift: 1, dt: 1.0,
            elements: Some(vec!["Li".to_string()]),
        }).unwrap();

        assert_eq!(result.n_atoms, 1);
        assert!(result.elements.contains(&"Li".to_string()));
        assert!(!result.elements.contains(&"Fe".to_string()));

        for (lag, &msd_val) in result.msd.iter().enumerate() {
            let expected = (v_li * lag as f64).powi(2);
            assert!((msd_val - expected).abs() < 1e-8, "Li MSD lag {}: {}", lag, msd_val);
        }
    }

    #[test]
    fn test_msd_time_shift_averaging() {
        // shift=1, tau=3, 10 帧 → origins: p=0..7 → 8 origins
        let traj = make_traj_static(10);
        let result = calc_msd(&traj, &MsdParams {
            tau: Some(3), shift: 1, dt: 2.0, elements: None,
        }).unwrap();
        assert_eq!(result.n_origins, 8);
        assert_eq!(result.time.len(), 3);
        assert!((result.time[0] - 0.0).abs() < 1e-10);
        assert!((result.time[1] - 2.0).abs() < 1e-10);
        assert!((result.time[2] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_write_msd() {
        use std::io::Read;
        let traj = make_traj_linear(10.0, 0.2, 5);
        let result = calc_msd(&traj, &MsdParams::default()).unwrap();
        let path = "/tmp/test_ferro.msd";
        write_msd(&result, path).expect("write_msd failed");

        let mut content = String::new();
        std::fs::File::open(path).unwrap().read_to_string(&mut content).unwrap();
        assert!(content.starts_with("# ferro v"), "missing version header");
        assert!(content.contains("# time[fs]"), "missing column header");
        assert!(content.contains("MSD"), "missing MSD title");
    }
}
