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
use ferro_core::{Table, Trajectory};
use ferro_core::error::ChemError;

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
    /// Linear-fit window as fractions of the MSD lag-time axis
    /// (`(fmin, fmax)`, `0 <= fmin < fmax <= 1`). `None` = no fit.
    pub fit_range: Option<(f64, f64)>,
}

impl Default for MsdParams {
    fn default() -> Self {
        MsdParams { tau: None, shift: 1, dt: 1.0, elements: None, fit_range: None }
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
    /// Linear-fit / self-diffusion result (`None` unless `fit_range` was set)
    pub fit: Option<MsdFit>,
}

/// Linear-fit result for self-diffusion coefficient extraction.
///
/// `D = slope / 6` (Einstein relation, 3-D isotropic). Unit conversions:
/// `D[cm²/s] = d_ang2_per_fs · 0.1`, `D[m²/s] = d_ang2_per_fs · 1e-5`.
#[derive(Debug, Clone)]
pub struct MsdFit {
    /// Lower fraction of the lag-time axis used for the fit
    pub frac_lo: f64,
    /// Upper fraction of the lag-time axis used for the fit
    pub frac_hi: f64,
    /// First time point of the fit window \[fs\]
    pub t_lo: f64,
    /// Last time point of the fit window \[fs\]
    pub t_hi: f64,
    /// Fitted slope of MSD vs time \[Å²/fs\]
    pub slope: f64,
    /// Fitted intercept \[Å²\]
    pub intercept: f64,
    /// Self-diffusion coefficient `slope / 6` \[Å²/fs\]
    pub d_ang2_per_fs: f64,
    /// Uncertainty of `d_ang2_per_fs` \[Å²/fs\]: `|D(first half) − D(second half)|` of
    /// the fit window, as `gmx msd` reports it. `NaN` when a half holds fewer than
    /// 2 points. Only meaningful if the MSD is linear over the whole window.
    pub d_err: f64,
    /// Coefficient of determination of the linear fit
    pub r2: f64,
    /// Number of points used in the fit
    pub n_points: usize,
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

    // 转置为 atom-major：各原子时间序列互相独立，可并行处理
    let mut by_atom: Vec<Vec<[f64; 3]>> = (0..n_atoms)
        .map(|j| frac.iter().map(|step| step[j]).collect())
        .collect();

    by_atom.par_iter_mut().for_each(|coords| {
        for i in 1..n_steps {
            let (prev, curr_and_later) = coords.split_at_mut(i);
            for (c, p) in curr_and_later[0].iter_mut().zip(prev[i - 1].iter()) {
                *c -= (*c - *p).round();
            }
        }
    });

    // 转置回 step-major
    for (i, step) in frac.iter_mut().enumerate() {
        for j in 0..n_atoms {
            step[j] = by_atom[j][i];
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
/// Returns `Err` if:
/// - The trajectory has fewer than 2 frames
/// - No atoms match the element filter
/// - `tau` exceeds the trajectory length
/// - Any frame is missing a cell (periodic path only)
///
/// # NPT handling
/// Total MSD uses the average of the origin- and endpoint-frame cell matrices
/// to convert fractional displacements to Cartesian. Directional MSD uses the
/// endpoint cell parameter (same simplified approximation as code1/msd.c).
pub fn calc_msd(traj: &Trajectory, params: &MsdParams) -> ferro_core::Result<MsdResult> {
    let n_steps = traj.n_frames();
    if n_steps < 2 {
        return Err(ChemError::ValidationError("trajectory requires at least 2 frames".into()));
    }

    // Fail fast on an obviously bad fit-range before the heavy parallel loop.
    if let Some((fmin, fmax)) = params.fit_range {
        if !(0.0..=1.0).contains(&fmin) || !(0.0..=1.0).contains(&fmax) || fmin >= fmax {
            return Err(ChemError::ValidationError(format!(
                "fit-range must satisfy 0 <= fmin < fmax <= 1, got [{fmin}, {fmax}]"
            )));
        }
    }

    // 确定参与计算的原子下标（按第一帧筛选元素）
    let ref_frame = traj.first()
        .ok_or_else(|| ChemError::ValidationError("empty trajectory".into()))?;
    let atom_indices: Vec<usize> = ref_frame.atoms.iter().enumerate()
        .filter(|(_, a)| match &params.elements {
            Some(elems) => elems.contains(&a.element),
            None => true,
        })
        .map(|(i, _)| i)
        .collect();
    if atom_indices.is_empty() {
        return Err(ChemError::ValidationError("no atoms match the element filter".into()));
    }
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
) -> ferro_core::Result<MsdResult> {
    let n_steps = traj.n_frames();

    // 先校验所有帧都有 cell
    if traj.frames.iter().any(|f| f.cell.is_none()) {
        return Err(ChemError::ValidationError("all frames must have a periodic cell for MSD".into()));
    }

    // 构建 frac[step][atom_local] = [fa, fb, fc]（串行，顺序依赖无法并行）
    let mut frac: Vec<Vec<[f64; 3]>> = traj.frames.iter().map(|frame| {
        let cell = frame.cell.as_ref().unwrap();
        atom_indices.iter().map(|&i| {
            let f = cell.cartesian_to_fractional(frame.atoms[i].position)
                .expect("cell is non-singular");
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
    if n_origins == 0 {
        return Err(ChemError::ValidationError("tau exceeds trajectory length".into()));
    }

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
                for (f_end, f_orig) in frac[p + i].iter().zip(frac[p].iter()) {
                    let dx = f_end[0] - f_orig[0];
                    let dy = f_end[1] - f_orig[1];
                    let dz = f_end[2] - f_orig[2];
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

    build_result(accum, tau, n_origins, n_atoms, elements, params)
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
) -> ferro_core::Result<MsdResult> {
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
    if n_origins == 0 {
        return Err(ChemError::ValidationError("tau exceeds trajectory length".into()));
    }

    let accum: Vec<[f64; 4]> = p_values.par_iter()
        .map(|&p| {
            let mut local = vec![[0.0f64; 4]; tau];
            for i in 0..tau {
                let mut sum = [0.0f64; 4];
                for (c_end, c_orig) in cart[p + i].iter().zip(cart[p].iter()) {
                    let dx = c_end[0] - c_orig[0];
                    let dy = c_end[1] - c_orig[1];
                    let dz = c_end[2] - c_orig[2];
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

    build_result(accum, tau, n_origins, n_atoms, elements, params)
}

/// Ordinary least-squares fit of total MSD vs time over a fractional window
/// of the lag-time axis. `frac = (fmin, fmax)` with `0 <= fmin < fmax <= 1`
/// mapped to indices `i_lo = round(fmin·(n-1))`, `i_hi = round(fmax·(n-1))`.
///
/// Returns slope/intercept, `D = slope / 6` (Einstein, 3-D isotropic), the
/// fit `R²` and `d_err` from refitting each half of the window. Errors on
/// invalid range, length mismatch, or a window with fewer than 2 points /
/// zero x-variance.
pub fn fit_diffusion(
    time: &[f64],
    msd: &[f64],
    frac: (f64, f64),
) -> ferro_core::Result<MsdFit> {
    let (fmin, fmax) = frac;
    if !(0.0..=1.0).contains(&fmin) || !(0.0..=1.0).contains(&fmax) || fmin >= fmax {
        return Err(ChemError::ValidationError(format!(
            "fit-range must satisfy 0 <= fmin < fmax <= 1, got [{fmin}, {fmax}]"
        )));
    }
    let n = time.len();
    if n != msd.len() {
        return Err(ChemError::ValidationError(
            "fit_diffusion: time/msd length mismatch".into(),
        ));
    }
    if n < 2 {
        return Err(ChemError::ValidationError(
            "MSD curve has fewer than 2 points; cannot fit".into(),
        ));
    }
    let last = (n - 1) as f64;
    let i_lo = (fmin * last).round() as usize;
    let i_hi = ((fmax * last).round() as usize).min(n - 1);
    if i_hi <= i_lo || (i_hi - i_lo + 1) < 2 {
        return Err(ChemError::ValidationError(format!(
            "fit-range [{fmin}, {fmax}] selects fewer than 2 points (i_lo={i_lo}, i_hi={i_hi})"
        )));
    }

    let xs = &time[i_lo..=i_hi];
    let ys = &msd[i_lo..=i_hi];
    let (slope, intercept) = ols(xs, ys).ok_or_else(|| {
        ChemError::ValidationError("degenerate fit window (zero x-variance)".into())
    })?;

    // gmx msd 的误差：窗口两半各拟合一次，D 之差。两半共用中点，
    // 每半至少 2 点，否则给 NaN 而不是报错 —— D 本身仍然有效
    let mid = (i_lo + i_hi) / 2;
    let half_d = |a: usize, b: usize| {
        (b > a).then(|| ols(&time[a..=b], &msd[a..=b])).flatten().map(|(k, _)| k / 6.0)
    };
    let d_err = match (half_d(i_lo, mid), half_d(mid, i_hi)) {
        (Some(d1), Some(d2)) => (d1 - d2).abs(),
        _ => f64::NAN,
    };

    let m = xs.len() as f64;
    let mean_y = ys.iter().sum::<f64>() / m;
    let ss_tot: f64 = ys.iter().map(|y| (y - mean_y).powi(2)).sum();
    let ss_res: f64 = xs
        .iter()
        .zip(ys)
        .map(|(x, y)| (y - (slope * x + intercept)).powi(2))
        .sum();
    let r2 = if ss_tot.abs() < f64::EPSILON {
        1.0
    } else {
        1.0 - ss_res / ss_tot
    };

    Ok(MsdFit {
        frac_lo: fmin,
        frac_hi: fmax,
        t_lo: xs[0],
        t_hi: xs[xs.len() - 1],
        slope,
        intercept,
        d_ang2_per_fs: slope / 6.0,
        d_err,
        r2,
        n_points: xs.len(),
    })
}

/// Ordinary least squares `y = slope·x + intercept`. `None` when x has no variance.
fn ols(xs: &[f64], ys: &[f64]) -> Option<(f64, f64)> {
    let m = xs.len() as f64;
    let sx: f64 = xs.iter().sum();
    let sy: f64 = ys.iter().sum();
    let sxx: f64 = xs.iter().map(|x| x * x).sum();
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| x * y).sum();
    let denom = m * sxx - sx * sx;
    if denom.abs() < f64::EPSILON {
        return None;
    }
    let slope = (m * sxy - sx * sy) / denom;
    Some((slope, (sy - slope * sx) / m))
}

/// Build an `MsdResult` from the parallel-reduction accumulation array,
/// computing the diffusion fit when `params.fit_range` is set.
fn build_result(
    accum: Vec<[f64; 4]>,
    tau: usize,
    n_origins: usize,
    n_atoms: usize,
    elements: Vec<String>,
    params: &MsdParams,
) -> ferro_core::Result<MsdResult> {
    let inv = 1.0 / n_origins as f64;
    let time:  Vec<f64> = (0..tau).map(|i| i as f64 * params.dt).collect();
    let msd:   Vec<f64> = (0..tau).map(|i| accum[i][0] * inv).collect();
    let msd_a: Vec<f64> = (0..tau).map(|i| accum[i][1] * inv).collect();
    let msd_b: Vec<f64> = (0..tau).map(|i| accum[i][2] * inv).collect();
    let msd_c: Vec<f64> = (0..tau).map(|i| accum[i][3] * inv).collect();
    let fit = match params.fit_range {
        Some(fr) => Some(fit_diffusion(&time, &msd, fr)?),
        None => None,
    };
    Ok(MsdResult {
        time, msd, msd_a, msd_b, msd_c,
        n_atoms, n_origins, params: params.clone(), elements, fit,
    })
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

/// Write MSD data to a tab-separated text file (`.msd`).
///
/// Columns: `time[fs]`, `msd[Ang^2]`, `msd_a[Ang^2]`, `msd_b[Ang^2]`, `msd_c[Ang^2]`
///
/// For periodic trajectories a/b/c are the crystal-axis directions.
impl MsdResult {
    /// Projects the result into the table the writers consume.
    ///
    /// Mean squared displacement: `time, msd, msd_a, msd_b, msd_c`.
    /// The `file` column is added by the caller when stacking several inputs
    /// (see `ferro_core::Table::concat_union`).
    pub fn to_tables(&self) -> Vec<(String, Table)> {
        let mut t = Table::new();
        t.push_num("time", self.time.clone())
            .push_num("msd", self.msd.clone())
            .push_num("msd_a", self.msd_a.clone())
            .push_num("msd_b", self.msd_b.clone())
            .push_num("msd_c", self.msd_c.clone());
        vec![("msd".to_string(), t)]
    }

    /// Parameter block for the comment header above the data.
    ///
    /// Only what the whole batch shares. Atom count, time origins, the fit's time
    /// window and everything fitted differ per input and go to the `[inputs]` list.
    pub fn meta_lines(&self) -> Vec<String> {
        let p = &self.params;
        let mut v = Vec::new();
        if let Some(tau) = p.tau {
            v.push(format!("tau      = {tau} frames"));
        }
        v.push(format!("shift    = {} frames", p.shift));
        v.push(format!("dt       = {} fs", p.dt));
        v.push(match &p.elements {
            Some(els) => format!("elements = {}", els.join(" ")),
            None => "elements = all".to_string(),
        });
        match p.fit_range {
            Some((lo, hi)) => {
                v.push(format!("fit      = [{lo:.2}, {hi:.2}] of each input's lag-time axis"));
                v.push("D        = slope / 6 [Ang^2/fs]; x0.1 -> cm^2/s".to_string());
                v.push("d_err    = |D(first half) - D(second half)| of the window (gmx msd)".to_string());
            }
            None => v.push("fit      = off".to_string()),
        }
        v
    }
}

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
            tau: Some(n), shift: 1, dt: 1.0, elements: None, fit_range: None,
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
            tau: Some(n), shift: 1, dt: 1.0, elements: None, fit_range: None,
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
            elements: Some(vec!["Li".to_string()]), fit_range: None,
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
            tau: Some(3), shift: 1, dt: 2.0, elements: None, fit_range: None,
        }).unwrap();
        assert_eq!(result.n_origins, 8);
        assert_eq!(result.time.len(), 3);
        assert!((result.time[0] - 0.0).abs() < 1e-10);
        assert!((result.time[1] - 2.0).abs() < 1e-10);
        assert!((result.time[2] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_to_tables_columns() {
        let traj = make_traj_linear(10.0, 0.2, 5);
        let result = calc_msd(&traj, &MsdParams::default()).unwrap();
        let (name, t) = result.to_tables().remove(0);
        assert_eq!(name, "msd");
        assert_eq!(t.names(), vec!["time", "msd", "msd_a", "msd_b", "msd_c"]);
        assert_eq!(t.n_rows(), result.time.len());
        assert!(t.validate().is_ok());
    }

    #[test]
    fn test_fit_diffusion_exact_line() {
        // msd = 6*D*t + c with known D → recovered slope/6 == D, R² == 1
        let d_true = 1.5e-4_f64; // Å²/fs
        let c = 0.7_f64;
        let dt = 2.0_f64;
        let n = 500;
        let time: Vec<f64> = (0..n).map(|i| i as f64 * dt).collect();
        let msd: Vec<f64> = time.iter().map(|&t| 6.0 * d_true * t + c).collect();

        let fit = fit_diffusion(&time, &msd, (0.2, 0.9)).unwrap();
        assert!((fit.d_ang2_per_fs - d_true).abs() < 1e-12,
            "D: expected {d_true}, got {}", fit.d_ang2_per_fs);
        assert!((fit.slope - 6.0 * d_true).abs() < 1e-12);
        assert!((fit.intercept - c).abs() < 1e-9);
        assert!((fit.r2 - 1.0).abs() < 1e-12, "R² = {}", fit.r2);
    }

    #[test]
    fn test_fit_range_index_mapping() {
        // n=11, dt=10 → time 0..100; frac (0.3,0.8) → indices 3..=8
        let dt = 10.0;
        let n = 11;
        let time: Vec<f64> = (0..n).map(|i| i as f64 * dt).collect();
        let msd: Vec<f64> = time.clone(); // slope 1
        let fit = fit_diffusion(&time, &msd, (0.3, 0.8)).unwrap();
        assert_eq!(fit.n_points, 6);
        assert!((fit.t_lo - 30.0).abs() < 1e-12);
        assert!((fit.t_hi - 80.0).abs() < 1e-12);
        assert!((fit.slope - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_d_err_zero_on_exact_line() {
        let time: Vec<f64> = (0..101).map(|i| i as f64).collect();
        let msd: Vec<f64> = time.iter().map(|&t| 6.0 * 2e-3 * t + 0.4).collect();
        let fit = fit_diffusion(&time, &msd, (0.1, 0.9)).unwrap();
        assert!(fit.d_err.abs() < 1e-12, "直线上两半的 D 应相同，d_err = {}", fit.d_err);
    }

    #[test]
    fn test_d_err_is_half_window_difference() {
        // 折线：t<=50 斜率 6，t>=50 斜率 12（在 t=50 连续）。窗口 [0,100] 的中点正好是
        // 拐点，两半各自是精确直线：D1 = 1、D2 = 2，d_err 必须恰为 1
        let time: Vec<f64> = (0..101).map(|i| i as f64).collect();
        let msd: Vec<f64> = time
            .iter()
            .map(|&t| if t <= 50.0 { 6.0 * t } else { 300.0 + 12.0 * (t - 50.0) })
            .collect();
        let fit = fit_diffusion(&time, &msd, (0.0, 1.0)).unwrap();
        assert!((fit.d_err - 1.0).abs() < 1e-9, "d_err 应为 |1 - 2| = 1，实为 {}", fit.d_err);
    }

    #[test]
    fn test_d_err_nan_when_half_too_short() {
        // 窗口只有 2 点：D 可算，但中点等于一端，半窗只剩 1 点
        let time = vec![0.0, 1.0, 2.0, 3.0];
        let msd = vec![0.0, 6.0, 12.0, 18.0];
        let fit = fit_diffusion(&time, &msd, (0.0, 0.34)).unwrap();
        assert_eq!(fit.n_points, 2);
        assert!((fit.d_ang2_per_fs - 1.0).abs() < 1e-12);
        assert!(fit.d_err.is_nan(), "半窗不足 2 点时 d_err 应为 NaN");
    }

    #[test]
    fn test_fit_range_invalid() {
        let time: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let msd = time.clone();
        assert!(fit_diffusion(&time, &msd, (0.8, 0.3)).is_err()); // fmin>=fmax
        assert!(fit_diffusion(&time, &msd, (-0.1, 0.5)).is_err()); // out of range
        assert!(fit_diffusion(&time, &msd, (0.2, 1.5)).is_err()); // out of range
        let t2 = vec![0.0, 1.0];
        let m2 = vec![0.0, 1.0];
        assert!(fit_diffusion(&t2, &m2, (0.0, 0.001)).is_err()); // <2 points
    }

    #[test]
    fn test_calc_msd_populates_fit() {
        let traj = make_traj_static(10);

        let none = calc_msd(&traj, &MsdParams::default()).unwrap();
        assert!(none.fit.is_none());

        let with = calc_msd(&traj, &MsdParams {
            fit_range: Some((0.0, 1.0)),
            ..MsdParams::default()
        }).unwrap();
        assert!(with.fit.is_some());
        let f = with.fit.unwrap();
        assert!((f.d_ang2_per_fs).abs() < 1e-12); // static traj → D = 0
    }

    #[test]
    fn test_meta_lines_hold_only_batch_shared_values() {
        let traj = make_traj_linear(10.0, 0.2, 6);
        let result = calc_msd(&traj, &MsdParams {
            fit_range: Some((0.0, 1.0)),
            ..MsdParams::default()
        }).unwrap();
        let meta = result.meta_lines().join("\n");
        assert!(meta.contains("fit      = [0.00, 1.00]"), "缺拟合窗口比例:\n{meta}");
        // 拟合出的数与原子数逐文件不同，写进共享区就会让第一个文件冒充全批
        let d = format!("{:.6e}", result.fit.as_ref().unwrap().d_ang2_per_fs);
        assert!(!meta.contains(&d), "D 的数值不该出现在共享区:\n{meta}");
        assert!(!meta.contains("atoms"), "原子数不该出现在共享区:\n{meta}");

        let plain = calc_msd(&traj, &MsdParams::default()).unwrap();
        assert!(plain.meta_lines().join("\n").contains("fit      = off"));
    }
}
