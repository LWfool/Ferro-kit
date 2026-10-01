//! Velocity autocorrelation function (VACF) and Green–Kubo self-diffusion.
//!
//! `Cv(m) = (1/N_atoms) Σⱼ (1/(N−m)) Σₜ vⱼ(t)·vⱼ(t+m)` — every lag averaged over all
//! `N − m` time origins, computed per atom and Cartesian axis with the FFT
//! autocorrelation in [`super::correlate`] (as `gmx velacc` and MDAnalysis
//! transport-analysis do).
//!
//! `diffusion(t) = (1/3)∫₀ᵗ Cv(τ)dτ` by the running trapezoidal rule; it levels off
//! at the self-diffusion coefficient D.
//!
//! Velocities are used as stored in `frame.velocities`, in the internal unit Å/fs
//! (the LAMMPS reader converts metal-unit Å/ps with `--metal-units`). `Cv` is then in
//! Å²/fs² and `diffusion` in Å²/fs.
//!
//! Parallelism: per atom.

use rayon::prelude::*;
use std::collections::BTreeSet;
use ferro_core::{Table, Trajectory};
use ferro_core::error::ChemError;

use super::correlate::{cumulative_trapezoid, resolve_max_lag, AutocorrPlan};

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Parameters for velocity autocorrelation function calculation.
#[derive(Debug, Clone)]
pub struct VacfParams {
    /// Longest lag in frames, `1 ..= n_frames − 1` (`None` = `n_frames / 2`)
    pub max_lag: Option<usize>,
    /// Time step per frame \[fs\] (default: 1.0)
    pub dt: f64,
    /// Elements to include (`None` = all atoms)
    pub elements: Option<Vec<String>>,
}

impl Default for VacfParams {
    fn default() -> Self {
        VacfParams { max_lag: None, dt: 1.0, elements: None }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Result of a velocity autocorrelation function calculation.
///
/// `vacf[m] = vacf_x[m] + vacf_y[m] + vacf_z[m]`.
#[derive(Debug, Clone)]
pub struct VacfResult {
    /// Lag-time axis \[fs\]; `time[m] = m · dt`, `m = 0 ..= max_lag`
    pub time: Vec<f64>,
    /// Total VACF Cv(t) \[Å²/fs²\]
    pub vacf: Vec<f64>,
    /// `vacf / vacf[0]` (NaN when `vacf[0] = 0`) — the normalised form `gmx velacc`
    /// writes by default
    pub vacf_norm: Vec<f64>,
    /// x-component of VACF \[Å²/fs²\]
    pub vacf_x: Vec<f64>,
    /// y-component of VACF \[Å²/fs²\]
    pub vacf_y: Vec<f64>,
    /// z-component of VACF \[Å²/fs²\]
    pub vacf_z: Vec<f64>,
    /// Running Green–Kubo integral `(1/3)∫₀ᵗ Cv dτ` (trapezoidal) \[Å²/fs\]
    pub diffusion: Vec<f64>,
    /// Frames in the trajectory
    pub n_frames: usize,
    /// Number of atoms included
    pub n_atoms: usize,
    /// Time origins at the longest lag, `n_frames − max_lag`
    pub min_origins: usize,
    pub params: VacfParams,
    /// Element types included, sorted alphabetically
    pub elements: Vec<String>,
}

// ─── 计算 ────────────────────────────────────────────────────────────────────

/// Computes the VACF, averaging each lag over all time origins.
///
/// Returns `Err` if the trajectory has fewer than 2 frames, `max_lag` is outside
/// `1 ..= n_frames − 1`, any frame lacks velocities, or no atom matches the
/// element filter.
pub fn calc_vacf(traj: &Trajectory, params: &VacfParams) -> ferro_core::Result<VacfResult> {
    let n_frames = traj.n_frames();
    traj.check_same_atoms()?;
    let max_lag = resolve_max_lag(n_frames, params.max_lag)?;
    if let Some(k) = traj.frames.iter().position(|f| f.velocities.is_none()) {
        return Err(ChemError::ValidationError(format!(
            "frame {k} has no velocities; VACF needs velocities in every frame"
        )));
    }

    // 按第一帧筛选原子下标
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

    // 逐原子逐轴：Σₜ v(t)v(t+m) 的全原点和，再在原子间累加
    let plan = AutocorrPlan::new(n_frames);
    let zero = || [vec![0.0; max_lag + 1], vec![0.0; max_lag + 1], vec![0.0; max_lag + 1]];
    let sums = atom_indices.par_iter()
        .map_init(|| (plan.worker(), vec![0.0; n_frames], vec![0.0; max_lag + 1]), |(ac, x, s), &i| {
            let mut out = zero();
            for (axis, o) in out.iter_mut().enumerate() {
                for (xk, f) in x.iter_mut().zip(&traj.frames) {
                    *xk = f.velocities.as_ref().unwrap()[i][axis];
                }
                ac.sums(x, s);
                o.copy_from_slice(s);
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
        });

    // 第 m 个 lag 有 N−m 个原点
    let per_lag = |v: &[f64]| -> Vec<f64> {
        v.iter().enumerate().map(|(m, s)| s / ((n_frames - m) as f64 * n_atoms as f64)).collect()
    };
    let [vacf_x, vacf_y, vacf_z] = [per_lag(&sums[0]), per_lag(&sums[1]), per_lag(&sums[2])];
    let vacf: Vec<f64> = (0..=max_lag).map(|m| vacf_x[m] + vacf_y[m] + vacf_z[m]).collect();
    // C(0) = 0 只在全体静止时出现，归一化无意义 —— 给 NaN（空字段），不给 0
    let vacf_norm: Vec<f64> = vacf.iter()
        .map(|v| if vacf[0] != 0.0 { v / vacf[0] } else { f64::NAN })
        .collect();
    let diffusion: Vec<f64> = cumulative_trapezoid(&vacf, params.dt).iter().map(|v| v / 3.0).collect();
    let time: Vec<f64> = (0..=max_lag).map(|m| m as f64 * params.dt).collect();

    Ok(VacfResult {
        time, vacf, vacf_norm, vacf_x, vacf_y, vacf_z, diffusion,
        n_frames, n_atoms, min_origins: n_frames - max_lag,
        params: params.clone(), elements,
    })
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

impl VacfResult {
    /// Projects the result into the table the writers consume.
    ///
    /// `time, vacf, vacf_norm, vacf_x, vacf_y, vacf_z, diffusion`.
    /// The `file` column is added by the caller when stacking several inputs
    /// (see `ferro_core::Table::concat_union`).
    pub fn to_tables(&self) -> Vec<(String, Table)> {
        let mut t = Table::new();
        t.push_num("time", self.time.clone())
            .push_num("vacf", self.vacf.clone())
            .push_num("vacf_norm", self.vacf_norm.clone())
            .push_num("vacf_x", self.vacf_x.clone())
            .push_num("vacf_y", self.vacf_y.clone())
            .push_num("vacf_z", self.vacf_z.clone())
            .push_num("diffusion", self.diffusion.clone());
        vec![("vacf".to_string(), t)]
    }

    /// Parameter block for the comment header above the data.
    ///
    /// Only what the whole batch shares; frames, atoms, the lag range and origins
    /// differ per input and go to the `[inputs]` list.
    pub fn meta_lines(&self) -> Vec<String> {
        let p = &self.params;
        vec![
            match p.max_lag {
                Some(m) => format!("max lag   = {m} frames"),
                None => "max lag   = half of each input's frames (see [inputs])".to_string(),
            },
            "origins   = all: lag m averages frames - m origins (FFT)".to_string(),
            format!("dt        = {} fs", p.dt),
            match &p.elements {
                Some(els) => format!("elements  = {}", els.join(" ")),
                None => "elements  = all".to_string(),
            },
            "units     = velocities Ang/fs; vacf Ang^2/fs^2; diffusion Ang^2/fs (x0.1 -> cm^2/s)".to_string(),
            "vacf_norm = vacf / vacf(0)".to_string(),
            "diffusion = (1/3) * running trapezoidal integral of vacf (Green-Kubo)".to_string(),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Frame, Trajectory};
    use nalgebra::Vector3;

    /// 伪随机但可复现（不引 rand：测试只需要「不规则」）
    fn jitter(k: usize) -> f64 {
        ((k as f64 * 12.9898).sin() * 43758.5453).fract() - 0.5
    }

    /// 每帧每个原子给定速度的轨迹
    fn traj_from(vel: impl Fn(usize, usize) -> Vector3<f64>, elems: &[&str], n: usize) -> Trajectory {
        let mut traj = Trajectory::new();
        for t in 0..n {
            let mut frame = Frame::new();
            for e in elems {
                frame.add_atom(Atom::new(*e, Vector3::zeros()));
            }
            frame.velocities = Some((0..elems.len()).map(|j| vel(t, j)).collect());
            traj.add_frame(frame);
        }
        traj
    }

    #[test]
    fn test_fft_matches_brute_force() {
        // 不规则速度，每个 lag 对全部原点暴力求平均，逐轴比
        let (n, na) = (33, 4);
        let traj = traj_from(|t, j| Vector3::new(jitter(3 * t + j), jitter(5 * t + 7 * j + 1), jitter(t * 11 + j + 3)),
            &["O"; 4], n);
        let r = calc_vacf(&traj, &VacfParams { max_lag: Some(n - 1), ..Default::default() }).unwrap();
        let v = |t: usize, j: usize| traj.frames[t].velocities.as_ref().unwrap()[j];
        for m in 0..n {
            let mut want = [0.0; 3];
            for t in 0..n - m {
                for j in 0..na {
                    let (a, b) = (v(t, j), v(t + m, j));
                    want[0] += a.x * b.x;
                    want[1] += a.y * b.y;
                    want[2] += a.z * b.z;
                }
            }
            let want = want.map(|w| w / ((n - m) * na) as f64);
            for (axis, got) in [&r.vacf_x, &r.vacf_y, &r.vacf_z].iter().enumerate() {
                assert!((got[m] - want[axis]).abs() < 1e-12, "lag {m} 轴 {axis}: {} vs {}", got[m], want[axis]);
            }
        }
    }

    #[test]
    fn test_constant_velocity_is_flat_and_norm_is_one() {
        let traj = traj_from(|_, _| Vector3::new(2.0, 1.0, 0.5), &["Li", "Li"], 10);
        let r = calc_vacf(&traj, &VacfParams::default()).unwrap();
        for (m, (&v, &nv)) in r.vacf.iter().zip(&r.vacf_norm).enumerate() {
            assert!((v - 5.25).abs() < 1e-10, "lag {m}: {v}");
            assert!((nv - 1.0).abs() < 1e-12, "lag {m}: vacf_norm {nv}");
        }
    }

    #[test]
    fn test_zero_velocity_gives_nan_norm_not_zero() {
        // C(0)=0 时归一化没有定义：给 NaN（空字段），不能冒充「测到了 0」
        let traj = traj_from(|_, _| Vector3::zeros(), &["Li"], 6);
        let r = calc_vacf(&traj, &VacfParams::default()).unwrap();
        assert!(r.vacf.iter().all(|v| v.abs() < 1e-15));
        assert!(r.vacf_norm.iter().all(|v| v.is_nan()));
        assert!(r.diffusion.iter().all(|v| v.abs() < 1e-15));
    }

    #[test]
    fn test_diffusion_is_trapezoidal() {
        // 恒定速度 → C 恒为 c，梯形积分 = c·t，diffusion = c·t/3；
        // 左矩形法会多出 c·dt/3，这里 m=0 必须恰为 0
        let traj = traj_from(|_, _| Vector3::new(1.0, 0.0, 0.0), &["Li"], 8);
        let dt = 2.0;
        let r = calc_vacf(&traj, &VacfParams { dt, ..Default::default() }).unwrap();
        for (m, &d) in r.diffusion.iter().enumerate() {
            let want = m as f64 * dt / 3.0;
            assert!((d - want).abs() < 1e-12, "lag {m}: {d} vs {want}");
        }
    }

    #[test]
    fn test_element_filter() {
        let traj = traj_from(|t, j| if j == 0 { Vector3::new(1.0, 0.0, 0.0) }
                                    else { Vector3::new((t as f64).cos(), (t as f64).sin(), 0.0) },
            &["Fe", "Li"], 10);
        let r = calc_vacf(&traj, &VacfParams { elements: Some(vec!["Fe".into()]), ..Default::default() }).unwrap();
        assert_eq!(r.n_atoms, 1);
        assert_eq!(r.elements, vec!["Fe".to_string()]);
        assert!(r.vacf.iter().all(|v| (v - 1.0).abs() < 1e-10));
    }

    #[test]
    fn test_default_max_lag_and_origins() {
        let traj = traj_from(|_, _| Vector3::new(1.0, 0.0, 0.0), &["Li"], 11);
        let r = calc_vacf(&traj, &VacfParams::default()).unwrap();
        assert_eq!(r.time.len(), 6);
        assert_eq!(r.min_origins, 6);
        assert_eq!(r.n_frames, 11);
    }

    #[test]
    fn test_missing_velocities_names_the_frame() {
        let mut traj = traj_from(|_, _| Vector3::new(1.0, 0.0, 0.0), &["Li"], 4);
        traj.frames[2].velocities = None;
        let err = calc_vacf(&traj, &VacfParams::default()).unwrap_err().to_string();
        assert!(err.contains("frame 2"), "错误信息应指出是哪一帧：{err}");
    }

    #[test]
    fn test_to_tables_columns() {
        let traj = traj_from(|_, _| Vector3::new(1.0, 0.0, 0.0), &["Li"], 6);
        let r = calc_vacf(&traj, &VacfParams::default()).unwrap();
        let (name, t) = r.to_tables().remove(0);
        assert_eq!(name, "vacf");
        assert_eq!(t.names(), vec!["time", "vacf", "vacf_norm", "vacf_x", "vacf_y", "vacf_z", "diffusion"]);
        assert_eq!(t.n_rows(), r.time.len());
        assert!(!r.meta_lines().join("\n").contains("atoms"), "原子数是逐文件的量，不进共享区");
    }
}
