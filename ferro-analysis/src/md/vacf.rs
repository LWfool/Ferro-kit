//! Velocity autocorrelation function (VACF) calculation and output.
//!
//! Cv(t) = <v(0)·v(t)> = (1/N) Σⱼ [vx_j(0)vx_j(t) + vy_j(0)vy_j(t) + vz_j(0)vz_j(t)]
//!
//! Workflow: `calc_vacf` → `write_vacf`.
//! Algorithm follows code1/velcorr.c (`EstimateVelocityCorr`):
//!   1. Read velocity vectors from all frames (no coordinate unwrapping needed — velocities have no PBC issues).
//!   2. For each time origin p (spaced by shift), compute velocity dot products for all lags i ∈ [0, tau).
//!   3. Average over all origins.
//!   4. Running integral D(t) = Σᵢ Cv(i)·dt / 3 → converges to self-diffusion coefficient D as t → ∞.
//!
//! Parallelism: per time-origin par_iter; each origin computed independently then reduced.
//!
//! **Unit note**: velocities use the raw values stored in `frame.velocities`.
//! The IO layer has not yet converted LAMMPS metal-unit velocities (Å/ps) to the internal standard (Å/fs).
//! For metal-unit trajectories, Cv is in (Å/ps)² and D is in Å²/ps;
//! multiply by 1e-6 (Cv) or 1e-3 (D) to obtain standard units.
//! This will be resolved when IO unit normalisation is implemented.

use rayon::prelude::*;
use std::collections::BTreeSet;
use std::io::{BufWriter, Write};
use ferro_core::Trajectory;
use super::gr::VERSION;

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Parameters for velocity autocorrelation function calculation.
#[derive(Debug, Clone)]
pub struct VacfParams {
    /// Correlation window size in frames (`None` = all frames)
    pub tau: Option<usize>,
    /// Time shift between origins in frames (default: 1)
    pub shift: usize,
    /// Time step per frame \[fs\] (default: 1.0)
    pub dt: f64,
    /// Elements to include (`None` = all atoms)
    pub elements: Option<Vec<String>>,
}

impl Default for VacfParams {
    fn default() -> Self {
        VacfParams { tau: None, shift: 1, dt: 1.0, elements: None }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Result of a velocity autocorrelation function calculation.
///
/// `vacf[i]` = `vacf_x[i] + vacf_y[i] + vacf_z[i]` = Cv(i·dt) in velocity² units.
/// `diffusion[i]` is the running integral `Σₖ₌₀ⁱ vacf[k]·dt / 3`, which converges
/// to the self-diffusion coefficient D as i → ∞.
#[derive(Debug, Clone)]
pub struct VacfResult {
    /// Time axis \[fs\]; `time[i] = i * dt`
    pub time: Vec<f64>,
    /// Total VACF Cv(t) = Cv_x + Cv_y + Cv_z \[vel²\]
    pub vacf: Vec<f64>,
    /// x-component of VACF \[vel²\]
    pub vacf_x: Vec<f64>,
    /// y-component of VACF \[vel²\]
    pub vacf_y: Vec<f64>,
    /// z-component of VACF \[vel²\]
    pub vacf_z: Vec<f64>,
    /// Running integral ∫₀ᵗ Cv(τ)dτ / 3 \[vel²·fs\]
    pub diffusion: Vec<f64>,
    /// Number of atoms included
    pub n_atoms: usize,
    /// Number of time origins averaged
    pub n_origins: usize,
    pub params: VacfParams,
    /// Element types included, sorted alphabetically
    pub elements: Vec<String>,
}

// ─── 计算 ────────────────────────────────────────────────────────────────────

/// Compute the velocity autocorrelation function with time-shift averaging.
///
/// The algorithm matches code1/velcorr.c `EstimateVelocityCorr`.
/// Unlike position-based analyses, no coordinate unwrapping is needed.
///
/// Returns `None` if:
/// - The trajectory has fewer than 2 frames
/// - Any frame is missing `velocities`
/// - No atoms match the element filter
/// - `tau` exceeds the trajectory length
pub fn calc_vacf(traj: &Trajectory, params: &VacfParams) -> Option<VacfResult> {
    let n_steps = traj.n_frames();
    if n_steps < 2 { return None; }

    // 校验所有帧都有速度数据
    if traj.frames.iter().any(|f| f.velocities.is_none()) { return None; }

    // 按第一帧筛选原子下标
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

    // 构建速度数组 vel[step][atom_local] = [vx, vy, vz]
    let vel: Vec<Vec<[f64; 3]>> = traj.frames.iter().map(|frame| {
        let vels = frame.velocities.as_ref().unwrap();
        atom_indices.iter().map(|&i| {
            let v = &vels[i];
            [v.x, v.y, v.z]
        }).collect()
    }).collect();

    // 枚举所有 origin：p + tau <= n_steps（内层循环访问 vel[p+i]，i ∈ [0,tau)）
    let p_values: Vec<usize> = (0..)
        .map(|k: usize| k * shift)
        .take_while(|&p| p + tau <= n_steps)
        .collect();
    let n_origins = p_values.len();
    if n_origins == 0 { return None; }

    // 并行计算各 origin 的局部累积 [tau][3]（x/y/z 三分量）
    let accum: Vec<[f64; 3]> = p_values.par_iter()
        .map(|&p| {
            let mut local = vec![[0.0f64; 3]; tau];
            for i in 0..tau {
                let mut sum = [0.0f64; 3];
                for j in 0..n_atoms {
                    sum[0] += vel[p][j][0] * vel[p + i][j][0];
                    sum[1] += vel[p][j][1] * vel[p + i][j][1];
                    sum[2] += vel[p][j][2] * vel[p + i][j][2];
                }
                local[i] = [
                    sum[0] / n_atoms as f64,
                    sum[1] / n_atoms as f64,
                    sum[2] / n_atoms as f64,
                ];
            }
            local
        })
        .reduce(
            || vec![[0.0f64; 3]; tau],
            |mut a, b| {
                for i in 0..tau {
                    a[i][0] += b[i][0];
                    a[i][1] += b[i][1];
                    a[i][2] += b[i][2];
                }
                a
            },
        );

    // 归一化：除以 n_origins
    let inv = 1.0 / n_origins as f64;
    let time:   Vec<f64> = (0..tau).map(|i| i as f64 * params.dt).collect();
    let vacf_x: Vec<f64> = (0..tau).map(|i| accum[i][0] * inv).collect();
    let vacf_y: Vec<f64> = (0..tau).map(|i| accum[i][1] * inv).collect();
    let vacf_z: Vec<f64> = (0..tau).map(|i| accum[i][2] * inv).collect();
    let vacf:   Vec<f64> = (0..tau).map(|i| vacf_x[i] + vacf_y[i] + vacf_z[i]).collect();

    // Running integral D(t) = Σᵢ Cv(i)·dt / 3 (rectangular approximation matching code1)
    let mut diffusion = vec![0.0f64; tau];
    let mut running = 0.0;
    for i in 0..tau {
        running += vacf[i] * params.dt / 3.0;
        diffusion[i] = running;
    }

    Some(VacfResult {
        time, vacf, vacf_x, vacf_y, vacf_z, diffusion,
        n_atoms, n_origins, params: params.clone(), elements,
    })
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

/// Write VACF data to a tab-separated text file (`.vacf`).
///
/// Columns: `time[fs]`, `vacf[v^2]`, `vacf_x`, `vacf_y`, `vacf_z`, `diffusion[v^2*fs]`
///
/// The velocity unit depends on the source trajectory.  For trajectories read
/// from LAMMPS metal-unit dump files, velocities are in Å/ps (not yet converted
/// to the internal standard Å/fs); this will be fixed when IO unit normalisation
/// is implemented.
pub fn write_vacf(result: &VacfResult, path: &str) -> std::io::Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(path)?);

    writeln!(w, "# ferro v{}", VERSION)?;
    writeln!(w, "# Velocity Autocorrelation Function (VACF)")?;
    writeln!(w, "# {}", "-".repeat(60))?;
    writeln!(w, "# tau     = {} frames", result.time.len())?;
    writeln!(w, "# shift   = {} frames", result.params.shift)?;
    writeln!(w, "# dt      = {} fs", result.params.dt)?;
    writeln!(w, "# atoms   = {}", result.n_atoms)?;
    writeln!(w, "# origins = {}", result.n_origins)?;
    write!(w,   "# elements:")?;
    for e in &result.elements { write!(w, " {}", e)?; }
    writeln!(w)?;
    writeln!(w, "# NOTE: velocity unit matches source file (see IO unit normalisation TODO)")?;
    writeln!(w, "# {}", "-".repeat(60))?;
    writeln!(w, "# time[fs]\tvacf[v^2]\tvacf_x[v^2]\tvacf_y[v^2]\tvacf_z[v^2]\tdiffusion[v^2*fs]")?;

    for i in 0..result.time.len() {
        writeln!(w, "{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}",
            result.time[i],
            result.vacf[i],
            result.vacf_x[i],
            result.vacf_y[i],
            result.vacf_z[i],
            result.diffusion[i])?;
    }
    Ok(())
}

// ─── 测试 ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Frame, Trajectory};
    use nalgebra::Vector3;

    /// 构建速度恒定的轨迹：每帧每个原子速度相同
    fn make_const_vel_traj(vx: f64, vy: f64, vz: f64, n: usize) -> Trajectory {
        let mut traj = Trajectory::new();
        for _ in 0..n {
            let mut frame = Frame::new();
            frame.add_atom(Atom::new("Li", Vector3::zeros()));
            frame.add_atom(Atom::new("Li", Vector3::zeros()));
            frame.velocities = Some(vec![
                Vector3::new(vx, vy, vz),
                Vector3::new(vx, vy, vz),
            ]);
            traj.add_frame(frame);
        }
        traj
    }

    /// 构建零速度轨迹
    fn make_zero_vel_traj(n: usize) -> Trajectory {
        make_const_vel_traj(0.0, 0.0, 0.0, n)
    }

    /// 构建含 Fe 和 Li 的混合速度轨迹
    fn make_mixed_elem_traj(n: usize) -> Trajectory {
        let mut traj = Trajectory::new();
        for i in 0..n {
            let mut frame = Frame::new();
            frame.add_atom(Atom::new("Fe", Vector3::zeros()));
            frame.add_atom(Atom::new("Li", Vector3::zeros()));
            let t = i as f64;
            frame.velocities = Some(vec![
                Vector3::new(1.0, 0.0, 0.0),         // Fe 速度恒定
                Vector3::new(t.cos(), t.sin(), 0.0),  // Li 速度随时间变化
            ]);
            traj.add_frame(frame);
        }
        traj
    }

    #[test]
    fn test_vacf_constant_velocity() {
        // 恒定速度 → Cv(t) = v² 对所有 t 均相等（VACF 平坦）
        let vx = 2.0_f64;
        let vy = 1.0_f64;
        let vz = 0.5_f64;
        let n = 10;
        let traj = make_const_vel_traj(vx, vy, vz, n);
        let params = VacfParams { tau: Some(5), ..Default::default() };
        let res = calc_vacf(&traj, &params).unwrap();

        let cv0 = vx*vx + vy*vy + vz*vz; // = v²，每个原子相同
        for (i, &v) in res.vacf.iter().enumerate() {
            assert!((v - cv0).abs() < 1e-10,
                "lag {}: expected Cv={:.4}, got {:.4}", i, cv0, v);
        }
    }

    #[test]
    fn test_vacf_zero_velocity() {
        // 零速度 → Cv(t) = 0 for all t
        let traj = make_zero_vel_traj(8);
        let res = calc_vacf(&traj, &VacfParams::default()).unwrap();
        for &v in &res.vacf {
            assert!(v.abs() < 1e-10, "zero-vel VACF should be 0, got {}", v);
        }
        for &d in &res.diffusion {
            assert!(d.abs() < 1e-10, "zero-vel diffusion should be 0, got {}", d);
        }
    }

    #[test]
    fn test_vacf_cv0_equals_mean_square_velocity() {
        // Cv(0) = (1/N) Σⱼ |vⱼ|² = 均方速度
        let vx = 3.0_f64;
        let vy = 4.0_f64;
        let vz = 0.0_f64;
        let traj = make_const_vel_traj(vx, vy, vz, 5);
        let res = calc_vacf(&traj, &VacfParams { tau: Some(3), ..Default::default() }).unwrap();
        let expected_cv0 = vx*vx + vy*vy + vz*vz; // = 25
        assert!((res.vacf[0] - expected_cv0).abs() < 1e-10,
            "Cv(0) = {:.4}, expected {:.4}", res.vacf[0], expected_cv0);
    }

    #[test]
    fn test_vacf_directional_components() {
        // vx=2, vy=0, vz=0 → vacf_x = 4, vacf_y = vacf_z = 0, total = 4
        let traj = make_const_vel_traj(2.0, 0.0, 0.0, 6);
        let res = calc_vacf(&traj, &VacfParams { tau: Some(3), ..Default::default() }).unwrap();
        for i in 0..3 {
            assert!((res.vacf_x[i] - 4.0).abs() < 1e-10, "vacf_x[{}] = {}", i, res.vacf_x[i]);
            assert!(res.vacf_y[i].abs() < 1e-10,          "vacf_y[{}] = {}", i, res.vacf_y[i]);
            assert!(res.vacf_z[i].abs() < 1e-10,          "vacf_z[{}] = {}", i, res.vacf_z[i]);
            assert!((res.vacf[i] - 4.0).abs() < 1e-10,   "vacf[{}] = {}", i, res.vacf[i]);
        }
    }

    #[test]
    fn test_vacf_element_filter() {
        // 只计算 Fe：Fe 速度恒定 (1,0,0) → Cv(t) = 1.0 for all t
        let traj = make_mixed_elem_traj(10);
        let res = calc_vacf(&traj, &VacfParams {
            tau: Some(5), shift: 1, dt: 1.0,
            elements: Some(vec!["Fe".to_string()]),
        }).unwrap();
        assert_eq!(res.n_atoms, 1);
        assert!(res.elements.contains(&"Fe".to_string()));
        assert!(!res.elements.contains(&"Li".to_string()));
        for &v in &res.vacf {
            assert!((v - 1.0).abs() < 1e-10, "Fe VACF should be 1.0, got {}", v);
        }
    }

    #[test]
    fn test_vacf_n_origins() {
        // n=10, tau=3, shift=1 → p + 3 <= 10 → p ∈ {0..7} → 8 origins
        let traj = make_zero_vel_traj(10);
        let res = calc_vacf(&traj, &VacfParams { tau: Some(3), shift: 1, ..Default::default() }).unwrap();
        assert_eq!(res.n_origins, 8);
    }

    #[test]
    fn test_vacf_diffusion_accumulates() {
        // 恒定速度 v²=1，dt=2.0 → D(t) = Σᵢ v²·dt/3 = i·(1·2/3)
        let traj = make_const_vel_traj(1.0, 0.0, 0.0, 8);
        let res = calc_vacf(&traj, &VacfParams { tau: Some(4), dt: 2.0, ..Default::default() }).unwrap();
        for (i, &d) in res.diffusion.iter().enumerate() {
            let expected = (i as f64 + 1.0) * 1.0 * 2.0 / 3.0;
            assert!((d - expected).abs() < 1e-10,
                "diffusion[{}] = {:.6}, expected {:.6}", i, d, expected);
        }
    }

    #[test]
    fn test_vacf_missing_velocities_returns_none() {
        // frame 缺少 velocities → 应返回 None
        let mut traj = Trajectory::new();
        let mut f1 = Frame::new();
        f1.add_atom(Atom::new("Li", Vector3::zeros()));
        f1.velocities = Some(vec![Vector3::new(1.0, 0.0, 0.0)]);
        let mut f2 = Frame::new();
        f2.add_atom(Atom::new("Li", Vector3::zeros()));
        // f2.velocities 故意留 None
        traj.add_frame(f1);
        traj.add_frame(f2);

        let res = calc_vacf(&traj, &VacfParams::default());
        assert!(res.is_none(), "should return None when velocities missing");
    }

    #[test]
    fn test_write_vacf() {
        use std::io::Read;
        let traj = make_const_vel_traj(1.0, 0.0, 0.0, 8);
        let res = calc_vacf(&traj, &VacfParams { tau: Some(4), ..Default::default() }).unwrap();
        let path = "/tmp/test_ferro.vacf";
        write_vacf(&res, path).expect("write_vacf failed");

        let mut content = String::new();
        std::fs::File::open(path).unwrap().read_to_string(&mut content).unwrap();
        assert!(content.starts_with("# ferro v"), "missing version header");
        assert!(content.contains("VACF"), "missing VACF title");
        assert!(content.contains("# time[fs]"), "missing column header");
        assert!(content.contains("diffusion"), "missing diffusion column");
    }
}
