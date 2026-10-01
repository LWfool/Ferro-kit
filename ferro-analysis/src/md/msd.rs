//! Mean squared displacement (MSD) and self-diffusion fit.
//!
//! 1. Periodic trajectories are unwrapped with the TOR scheme (von Bülow et al.
//!    2020; Bullerjahn et al. 2023): each step adds the minimum-image displacement
//!    in the **later** frame's box, so NPT box fluctuations do not inflate the MSD.
//!    Non-periodic trajectories are used as they are.
//! 2. MSD at every lag `m` averages over **all** `N − m` time origins, computed per
//!    atom and Cartesian axis in `O(N log N)` with the FFT algorithm of
//!    Calandrini et al. (2011): `MSD(m) = S1(m) − 2·S2(m)`.
//!
//! Parallelism: per atom (each atom's three series are independent).

use rayon::prelude::*;
use std::collections::BTreeSet;
use ferro_core::{Table, Trajectory};
use ferro_core::error::ChemError;
use nalgebra::Vector3;

use super::correlate::{resolve_max_lag, AutocorrPlan};

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Parameters for MSD calculation.
#[derive(Debug, Clone)]
pub struct MsdParams {
    /// Longest lag in frames, `1 ..= n_frames − 1` (`None` = `n_frames / 2`)
    pub max_lag: Option<usize>,
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
        MsdParams { max_lag: None, dt: 1.0, elements: None, fit_range: None }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Result of an MSD calculation.
///
/// `msd_x/y/z` are Cartesian components and add up to `msd` for any cell shape.
#[derive(Debug, Clone)]
pub struct MsdResult {
    /// Lag-time axis \[fs\]; `time[m] = m · dt`, `m = 0 ..= max_lag`
    pub time: Vec<f64>,
    /// Total MSD \[Å²\]
    pub msd: Vec<f64>,
    /// Cartesian x component \[Å²\]
    pub msd_x: Vec<f64>,
    /// Cartesian y component \[Å²\]
    pub msd_y: Vec<f64>,
    /// Cartesian z component \[Å²\]
    pub msd_z: Vec<f64>,
    /// Frames in the trajectory
    pub n_frames: usize,
    /// Number of atoms included in the calculation
    pub n_atoms: usize,
    /// Time origins at the longest lag, `n_frames − max_lag` — the fewest of any
    /// lag (lag `m` averages `n_frames − m`)
    pub min_origins: usize,
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

/// Unwrapped Cartesian series, `[atom][axis][frame]`.
type Series = Vec<[Vec<f64>; 3]>;

/// TOR unwrapping (Bullerjahn et al. 2023, eq 2), in 3-D:
/// `u[i+1] = u[i] + mic_{i+1}(w[i+1] − w[i])`, `u[0] = w[0]`, where `mic_{i+1}`
/// is the minimum image in frame `i+1`'s cell, taken in fractional coordinates as
/// `f − ⌊f + 1/2⌋`. A step is far below half a cell, so this is exact for any
/// cell shape.
fn unwrap_tor(traj: &Trajectory, atom_indices: &[usize]) -> ferro_core::Result<Series> {
    // 每帧的 (Mᵀ, (Mᵀ)⁻¹)：行优先 matrix 的行是晶格矢量，cart = Mᵀ·frac
    let boxes = traj.frames.iter().map(|f| {
        let m = f.cell.as_ref()
            .ok_or_else(|| ChemError::ValidationError(
                "all frames must have a periodic cell for MSD".into()))?
            .matrix.transpose();
        let inv = m.try_inverse()
            .ok_or_else(|| ChemError::ValidationError("cell matrix is singular".into()))?;
        Ok((m, inv))
    }).collect::<ferro_core::Result<Vec<_>>>()?;

    Ok(atom_indices.par_iter().map(|&i| {
        let n = traj.n_frames();
        let mut s = [Vec::with_capacity(n), Vec::with_capacity(n), Vec::with_capacity(n)];
        let mut w_prev = traj.frames[0].atoms[i].position;
        let mut u = w_prev;
        for (k, frame) in traj.frames.iter().enumerate() {
            if k > 0 {
                let w = frame.atoms[i].position;
                let (m, inv) = &boxes[k];
                let f = inv * (w - w_prev);
                // ⌊x + 1/2⌋ 而非 round()：与论文式 2 逐字一致（两者只在恰好 ±0.5 处不同）
                let f_mic = f.map(|c| c - (c + 0.5).floor());
                u += m * f_mic;
                w_prev = w;
            }
            for (axis, v) in s.iter_mut().zip([u.x, u.y, u.z]) {
                axis.push(v);
            }
        }
        s
    }).collect())
}

/// MSD summed over atoms, per axis: `out[axis][m]`, `m = 0 ..= max_lag`, each lag
/// averaged over its `N − m` origins.
///
/// Calandrini et al. (2011): `MSD(m) = S1(m) − 2·S2(m)` with
/// `S2(m) = Σ_t x[t]·x[t+m] / (N − m)` (FFT) and `S1` from the recursion
/// `Q(m) = Q(m−1) − x[m−1]² − x[N−m]²`, `Q(0) = 2Σx²`, `S1(m) = Q(m) / (N − m)`.
fn msd_sums(series: &Series, max_lag: usize) -> [Vec<f64>; 3] {
    let n = series[0][0].len();
    let plan = AutocorrPlan::new(n);
    let zero = || [vec![0.0; max_lag + 1], vec![0.0; max_lag + 1], vec![0.0; max_lag + 1]];
    series.par_iter()
        .map_init(|| (plan.worker(), vec![0.0; n], vec![0.0; max_lag + 1]), |(ac, x, s2), atom| {
            let mut out = zero();
            for (axis, raw) in atom.iter().enumerate() {
                // MSD 与平移无关；减去均值让 S1 与 2·S2 的量级接近原子的活动范围，
                // 否则两个 ~|r|² 的大数相减会吃掉有效位
                let mean = raw.iter().sum::<f64>() / n as f64;
                for (xi, &r) in x.iter_mut().zip(raw) {
                    *xi = r - mean;
                }
                ac.sums(x, s2);
                let mut q = 2.0 * x.iter().map(|v| v * v).sum::<f64>();
                for m in 0..=max_lag {
                    if m > 0 {
                        q -= x[m - 1] * x[m - 1] + x[n - m] * x[n - m];
                    }
                    let origins = (n - m) as f64;
                    out[axis][m] = (q - 2.0 * s2[m]) / origins;
                }
            }
            out
        })
        .reduce(zero, |mut a, b| {
            for (aa, bb) in a.iter_mut().zip(&b) {
                for (x, y) in aa.iter_mut().zip(bb) {
                    *x += y;
                }
            }
            a
        })
}

// ─── 计算 ────────────────────────────────────────────────────────────────────

/// Computes the MSD, averaging each lag over all time origins.
///
/// Returns `Err` if:
/// - The trajectory has fewer than 2 frames
/// - `max_lag` is outside `1 ..= n_frames − 1`
/// - No atoms match the element filter
/// - A periodic trajectory has a frame without a cell, or a singular cell
/// - `fit_range` is invalid or selects fewer than 2 points
pub fn calc_msd(traj: &Trajectory, params: &MsdParams) -> ferro_core::Result<MsdResult> {
    let n_frames = traj.n_frames();
    traj.check_same_atoms()?;
    let max_lag = resolve_max_lag(n_frames, params.max_lag)?;
    // 在重计算之前先挡掉明显错误的拟合窗口
    if let Some((fmin, fmax)) = params.fit_range {
        if !(0.0..=1.0).contains(&fmin) || !(0.0..=1.0).contains(&fmax) || fmin >= fmax {
            return Err(ChemError::ValidationError(format!(
                "fit-range must satisfy 0 <= fmin < fmax <= 1, got [{fmin}, {fmax}]"
            )));
        }
    }
    // 按第一帧筛选元素
    let ref_frame = &traj.frames[0];
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
    let elements: Vec<String> = atom_indices.iter()
        .map(|&i| ref_frame.atoms[i].element.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();

    let series: Series = if ref_frame.cell.is_some() {
        unwrap_tor(traj, &atom_indices)?
    } else {
        atom_indices.iter().map(|&i| {
            let p: Vec<Vector3<f64>> = traj.frames.iter().map(|f| f.atoms[i].position).collect();
            [p.iter().map(|v| v.x).collect(), p.iter().map(|v| v.y).collect(),
             p.iter().map(|v| v.z).collect()]
        }).collect()
    };

    let [sx, sy, sz] = msd_sums(&series, max_lag);
    let inv = 1.0 / n_atoms as f64;
    let msd_x: Vec<f64> = sx.iter().map(|v| v * inv).collect();
    let msd_y: Vec<f64> = sy.iter().map(|v| v * inv).collect();
    let msd_z: Vec<f64> = sz.iter().map(|v| v * inv).collect();
    let msd: Vec<f64> = (0..=max_lag).map(|m| msd_x[m] + msd_y[m] + msd_z[m]).collect();
    let time: Vec<f64> = (0..=max_lag).map(|m| m as f64 * params.dt).collect();

    let fit = match params.fit_range {
        Some(fr) => Some(fit_diffusion(&time, &msd, fr)?),
        None => None,
    };
    Ok(MsdResult {
        time, msd, msd_x, msd_y, msd_z,
        n_frames, n_atoms, min_origins: n_frames - max_lag,
        params: params.clone(), elements, fit,
    })
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

// ─── 输出函数 ────────────────────────────────────────────────────────────────

impl MsdResult {
    /// Projects the result into the table the writers consume.
    ///
    /// Mean squared displacement: `time, msd, msd_x, msd_y, msd_z`.
    /// The `file` column is added by the caller when stacking several inputs
    /// (see `ferro_core::Table::concat_union`).
    pub fn to_tables(&self) -> Vec<(String, Table)> {
        let mut t = Table::new();
        t.push_num("time", self.time.clone())
            .push_num("msd", self.msd.clone())
            .push_num("msd_x", self.msd_x.clone())
            .push_num("msd_y", self.msd_y.clone())
            .push_num("msd_z", self.msd_z.clone());
        vec![("msd".to_string(), t)]
    }

    /// Parameter block for the comment header above the data.
    ///
    /// Only what the whole batch shares. Atom count, time origins, the fit's time
    /// window and everything fitted differ per input and go to the `[inputs]` list.
    pub fn meta_lines(&self) -> Vec<String> {
        let p = &self.params;
        let mut v = Vec::new();
        v.push(match p.max_lag {
            Some(m) => format!("max lag  = {m} frames"),
            None => "max lag  = half of each input's frames (see [inputs])".to_string(),
        });
        v.push("origins  = all: lag m averages frames - m origins (FFT, Calandrini 2011)".to_string());
        v.push("unwrap   = TOR for periodic inputs: minimum image in the later frame's cell (Bullerjahn 2023)".to_string());
        v.push("axes     = Cartesian x/y/z; msd = msd_x + msd_y + msd_z".to_string());
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

    #[test]
    fn test_msd_rejects_frames_with_different_atoms() {
        // 末帧少原子（截断轨迹）此前在按第 0 帧原子下标取坐标时越界 panic
        let mut traj = make_traj_static(5);
        traj.frames[4].atoms.truncate(10);
        let err = calc_msd(&traj, &MsdParams::default()).expect_err("原子数不一致应报错而不是 panic");
        assert!(err.to_string().contains("frame 4 has 10 atoms"), "{err}");
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

    /// 暴力参考：对给定的逐帧位置（已解包裹），每个 lag 对全部原点直接求平均
    fn brute_msd(pos: &[Vec<Vector3<f64>>], max_lag: usize) -> Vec<[f64; 3]> {
        let n = pos.len();
        let na = pos[0].len() as f64;
        (0..=max_lag).map(|m| {
            let mut acc = [0.0; 3];
            for t in 0..n - m {
                for (a, b) in pos[t + m].iter().zip(&pos[t]) {
                    let d = a - b;
                    acc[0] += d.x * d.x;
                    acc[1] += d.y * d.y;
                    acc[2] += d.z * d.z;
                }
            }
            acc.map(|v| v / ((n - m) as f64 * na))
        }).collect()
    }

    /// 独立实现的 TOR：逐步加后一帧盒子下的最小镜像位移（走 Cell::minimum_image）
    fn tor_reference(traj: &Trajectory) -> Vec<Vec<Vector3<f64>>> {
        let mut u: Vec<Vector3<f64>> = traj.frames[0].atoms.iter().map(|a| a.position).collect();
        let mut out = vec![u.clone()];
        for k in 1..traj.n_frames() {
            let cell = traj.frames[k].cell.as_ref().unwrap();
            for (j, uj) in u.iter_mut().enumerate() {
                let d = traj.frames[k].atoms[j].position - traj.frames[k - 1].atoms[j].position;
                *uj += cell.minimum_image(d).unwrap();
            }
            out.push(u.clone());
        }
        out
    }

    /// 伪随机但可复现的数（不引 rand：测试只需要「不规则」）
    fn jitter(k: usize) -> f64 {
        ((k as f64 * 12.9898).sin() * 43758.5453).fract() - 0.5
    }

    /// NPT 三斜轨迹：盒子逐帧伸缩，原子做随机游走并被包裹回盒内
    fn make_traj_npt_triclinic(n: usize, n_atoms: usize) -> Trajectory {
        let mut traj = Trajectory::new();
        let mut cart: Vec<Vector3<f64>> = (0..n_atoms)
            .map(|j| Vector3::new(jitter(j) * 8.0 + 4.0, jitter(j + 99) * 8.0 + 4.0, jitter(j + 7) * 8.0 + 4.0))
            .collect();
        for i in 0..n {
            let s = 1.0 + 0.03 * jitter(1000 + i);
            let cell = Cell::from_lengths_angles(8.0 * s, 8.5 * s, 9.0 * s, 80.0, 95.0, 105.0).unwrap();
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            for (j, c) in cart.iter_mut().enumerate() {
                *c += Vector3::new(jitter(5 * i + j), jitter(7 * i + 3 * j + 1), jitter(11 * i + j + 2)) * 1.2;
                frame.add_atom(Atom::new("Fe", cell.wrap_position(*c).unwrap()));
            }
            traj.add_frame(frame);
        }
        traj
    }

    #[test]
    fn test_fft_matches_brute_force_npt_triclinic() {
        // 盒子伸缩 + 三斜 + 跨边界：FFT 全原点平均必须与暴力求和逐点一致
        let traj = make_traj_npt_triclinic(40, 5);
        let r = calc_msd(&traj, &MsdParams { max_lag: Some(39), ..MsdParams::default() }).unwrap();
        let reference = brute_msd(&tor_reference(&traj), 39);
        for (m, want) in reference.iter().enumerate() {
            for (axis, got) in [&r.msd_x, &r.msd_y, &r.msd_z].iter().enumerate() {
                assert!((got[m] - want[axis]).abs() < 1e-9,
                    "lag {m} 轴 {axis}: fft {} vs brute {}", got[m], want[axis]);
            }
        }
    }

    #[test]
    fn test_fft_matches_brute_force_nonperiodic() {
        let mut traj = Trajectory::new();
        for i in 0..25 {
            let mut f = Frame::new();
            for j in 0..3 {
                let k = i * 3 + j;
                f.add_atom(Atom::new("O", Vector3::new(jitter(k) * 5.0 + i as f64 * 0.1, jitter(k + 50), jitter(k + 90) * 2.0)));
            }
            traj.add_frame(f);
        }
        let r = calc_msd(&traj, &MsdParams { max_lag: Some(20), ..MsdParams::default() }).unwrap();
        let pos: Vec<Vec<Vector3<f64>>> = traj.frames.iter()
            .map(|f| f.atoms.iter().map(|a| a.position).collect()).collect();
        for (m, want) in brute_msd(&pos, 20).iter().enumerate() {
            let total: f64 = want.iter().sum();
            assert!((r.msd[m] - total).abs() < 1e-9, "lag {m}: {} vs {total}", r.msd[m]);
        }
    }

    #[test]
    fn test_components_sum_to_total() {
        let traj = make_traj_npt_triclinic(20, 4);
        let r = calc_msd(&traj, &MsdParams::default()).unwrap();
        for m in 0..r.msd.len() {
            let s = r.msd_x[m] + r.msd_y[m] + r.msd_z[m];
            assert!((r.msd[m] - s).abs() < 1e-12, "lag {m}: 三斜盒子下分量之和也必须等于总量");
        }
    }

    #[test]
    fn test_tor_crosses_skewed_axis() {
        // 固定三斜盒子里匀速直线运动、每步都可能跨斜边界：解包裹后 MSD = (|v|·m)²
        let cell = Cell::from_lengths_angles(6.0, 6.0, 6.0, 70.0, 75.0, 60.0).unwrap();
        let v = Vector3::new(0.9, -0.7, 1.1);
        let mut traj = Trajectory::new();
        for i in 0..12 {
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            let p = Vector3::new(5.5, 0.2, 5.8) + v * i as f64;
            f.add_atom(Atom::new("Li", cell.wrap_position(p).unwrap()));
            traj.add_frame(f);
        }
        let r = calc_msd(&traj, &MsdParams { max_lag: Some(11), ..MsdParams::default() }).unwrap();
        for (m, &got) in r.msd.iter().enumerate() {
            let want = (v.norm() * m as f64).powi(2);
            assert!((got - want).abs() < 1e-9, "lag {m}: {got} vs {want}");
        }
    }

    #[test]
    fn test_tor_does_not_follow_lattice_scaling() {
        // 原子先跨一次边界（x: 9.5 → 10.5，包裹后 0.5），之后分数坐标不动而盒子 10→11→10。
        // 格点视图（分数坐标解包裹 × 当帧盒长）：1.05 × 11 = 11.55，盒子一胀就凭空多走
        // 1.05 Å，且随镜像编号放大；TOR 只加最小镜像位移：+0.05、−0.05。
        let frames = [(10.0, 0.95), (10.0, 0.05), (11.0, 0.05), (10.0, 0.05)];
        let mut traj = Trajectory::new();
        for (len, fx) in frames {
            let cell = Cell::from_lengths_angles(len, len, len, 90.0, 90.0, 90.0).unwrap();
            let mut fr = Frame::with_cell(cell.clone(), [true; 3]);
            fr.add_atom(Atom::new("Na", cell.fractional_to_cartesian(Vector3::new(fx, 0.5, 0.5))));
            traj.add_frame(fr);
        }
        let r = calc_msd(&traj, &MsdParams { max_lag: Some(1), ..MsdParams::default() }).unwrap();
        // lag 1 的三个原点，x 位移 TOR: 1.0, 0.05, -0.05
        let tor = (1.0 + 0.05_f64.powi(2) * 2.0) / 3.0;
        let lattice = (1.0 + 1.05_f64.powi(2) * 2.0) / 3.0;
        assert!((r.msd_x[1] - tor).abs() < 1e-9,
            "msd_x lag 1 应为 TOR 的 {tor}，实为 {}（格点视图会得 {lattice}）", r.msd_x[1]);
    }

    #[test]
    fn test_msd_linear_motion() {
        // 单原子沿 x 匀速 v：每个原点的位移都是 v·m，全原点平均仍为 (v·m)²
        let v = 0.3;
        let traj = make_traj_linear(20.0, v, 6);
        let r = calc_msd(&traj, &MsdParams { max_lag: Some(5), ..MsdParams::default() }).unwrap();
        for (m, (&tot, &x)) in r.msd.iter().zip(&r.msd_x).enumerate() {
            let want = (v * m as f64).powi(2);
            assert!((tot - want).abs() < 1e-9, "lag {m}: {tot}");
            assert!((x - want).abs() < 1e-9, "msd_x lag {m}: {x}");
        }
    }

    #[test]
    fn test_msd_unwrap_across_boundary() {
        let (a, v, n) = (5.0, 0.5, 5);
        let cell = Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for i in 0..n {
            let x_raw = 4.8 + i as f64 * v;
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("Fe", Vector3::new(x_raw - (x_raw / a).floor() * a, 0.0, 0.0)));
            traj.add_frame(frame);
        }
        let r = calc_msd(&traj, &MsdParams { max_lag: Some(4), ..MsdParams::default() }).unwrap();
        for (m, &got) in r.msd.iter().enumerate() {
            assert!((got - (v * m as f64).powi(2)).abs() < 1e-9, "lag {m}: {got}");
        }
    }

    #[test]
    fn test_msd_element_filter() {
        let (a, v_li) = (10.0, 0.4);
        let cell = Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for i in 0..5 {
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("Fe", Vector3::new(1.0, 1.0, 1.0)));
            frame.add_atom(Atom::new("Li", Vector3::new(i as f64 * v_li, 0.0, 0.0)));
            traj.add_frame(frame);
        }
        let r = calc_msd(&traj, &MsdParams {
            max_lag: Some(4), elements: Some(vec!["Li".to_string()]), ..MsdParams::default()
        }).unwrap();
        assert_eq!(r.n_atoms, 1);
        assert_eq!(r.elements, vec!["Li".to_string()]);
        for (m, &got) in r.msd.iter().enumerate() {
            assert!((got - (v_li * m as f64).powi(2)).abs() < 1e-9, "Li lag {m}: {got}");
        }
    }

    #[test]
    fn test_default_max_lag_is_half_and_origins() {
        // 11 帧 → 默认 max_lag = 5，时间轴 0..=5，最长 lag 的原点数 11 − 5 = 6
        let traj = make_traj_static(11);
        let r = calc_msd(&traj, &MsdParams { dt: 2.0, ..MsdParams::default() }).unwrap();
        assert_eq!(r.time.len(), 6);
        assert!((r.time[5] - 10.0).abs() < 1e-12);
        assert_eq!(r.min_origins, 6);
        assert_eq!(r.n_frames, 11);
    }

    #[test]
    fn test_max_lag_out_of_range() {
        let traj = make_traj_static(6);
        for bad in [0, 6, 99] {
            assert!(calc_msd(&traj, &MsdParams { max_lag: Some(bad), ..MsdParams::default() }).is_err(),
                "max_lag = {bad} 应被拒绝");
        }
        assert!(calc_msd(&traj, &MsdParams { max_lag: Some(5), ..MsdParams::default() }).is_ok());
    }

    #[test]
    fn test_to_tables_columns() {
        let traj = make_traj_linear(10.0, 0.2, 5);
        let result = calc_msd(&traj, &MsdParams::default()).unwrap();
        let (name, t) = result.to_tables().remove(0);
        assert_eq!(name, "msd");
        assert_eq!(t.names(), vec!["time", "msd", "msd_x", "msd_y", "msd_z"]);
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
