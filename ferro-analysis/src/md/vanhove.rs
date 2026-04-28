//! Van Hove self-correlation function Gs(r, τ) calculation and output.
//!
//! Gs(r, τ) gives the probability distribution of atomic displacement distances over time interval τ
//! (discrete PMF, sum = 1).
//!
//! Workflow: `calc_vanhove` → `write_vanhove`.
//! Algorithm follows code1/vanhove.c (`EstimateVanHove`):
//!   1. Fractional coordinate conversion + unwrapping (identical to msd.rs).
//!   2. Convert back to absolute Cartesian coordinates.
//!   3. For each time origin p, compute |r(p+τ) − r(p)| and accumulate into a histogram.
//!   4. Normalise: gs[i] /= n_origins × n_atoms.
//!
//! Parallelism: per time-origin par_iter; each origin computed independently then reduced.

use rayon::prelude::*;
use std::collections::BTreeSet;
use std::io::{BufWriter, Write};
use ferro_core::Trajectory;
use super::gr::VERSION;
use super::msd::unwrap_frac;

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Parameters for van Hove self-correlation function calculation.
#[derive(Debug, Clone)]
pub struct VanHoveParams {
    /// Time lag in frames (`None` = `n_frames − 1`)
    pub tau: Option<usize>,
    /// Time shift between origins in frames (default: 1)
    pub shift: usize,
    /// Time step per frame \[fs\] (default: 1.0)
    pub dt: f64,
    /// Minimum displacement for histogram \[Å\] (default: 0.0)
    pub r_min: f64,
    /// Maximum displacement for histogram \[Å\] (default: 10.0)
    pub r_max: f64,
    /// Bin width \[Å\] (default: 0.01)
    pub dr: f64,
    /// Elements to include (`None` = all atoms)
    pub elements: Option<Vec<String>>,
}

impl Default for VanHoveParams {
    fn default() -> Self {
        VanHoveParams {
            tau: None, shift: 1, dt: 1.0,
            r_min: 0.0, r_max: 10.0, dr: 0.01,
            elements: None,
        }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Result of a van Hove self-correlation calculation.
///
/// `gs` is a discrete probability distribution of atomic displacements:
/// `gs[i]` is the fraction of (atom, origin) pairs whose displacement fell
/// in bin `i`.  The sum of all `gs` values equals 1.0.
#[derive(Debug, Clone)]
pub struct VanHoveResult {
    /// Bin-centre displacement values \[Å\]
    pub r: Vec<f64>,
    /// Gs(r, τ): normalized displacement histogram (sum = 1)
    pub gs: Vec<f64>,
    /// Actual τ used \[frames\]
    pub tau_frames: usize,
    /// Actual τ in physical time \[fs\]
    pub time: f64,
    /// Number of atoms included
    pub n_atoms: usize,
    /// Number of time origins averaged
    pub n_origins: usize,
    pub params: VanHoveParams,
    /// Element types included, sorted by atomic number
    pub elements: Vec<String>,
}

// ─── 计算 ────────────────────────────────────────────────────────────────────

/// Compute the van Hove self-correlation function Gs(r, τ) with time-shift averaging.
///
/// The algorithm matches code1/vanhove.c `EstimateVanHove`:
/// 1. Convert atom Cartesian coordinates to fractional (periodic) or keep Cartesian (non-periodic).
/// 2. Unwrap fractional coordinates across periodic boundaries.
/// 3. Convert back to absolute Cartesian positions.
/// 4. For each time origin p (spaced by `shift`), compute the scalar displacement
///    `|r(p+τ) − r(p)|` for every selected atom and accumulate into a histogram.
/// 5. Normalise: `gs[i] = count[i] / (n_origins × n_atoms)`.
///
/// Returns `None` if:
/// - The trajectory has fewer than 2 frames
/// - No atoms match the element filter
/// - `τ` exceeds `n_frames − 1`
/// - Any frame is missing a cell (periodic path only)
pub fn calc_vanhove(traj: &Trajectory, params: &VanHoveParams) -> Option<VanHoveResult> {
    let n_steps = traj.n_frames();
    if n_steps < 2 { return None; }

    // 按第一帧筛选参与计算的原子下标
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

    // tau 不超过 n_steps-1；None 时使用全轨迹长度
    let tau = params.tau.unwrap_or(n_steps - 1).min(n_steps - 1).max(1);
    let shift = params.shift.max(1);

    // 收集元素列表（用于输出文件头）
    let elements: Vec<String> = {
        let mut set = BTreeSet::new();
        for &i in &atom_indices { set.insert(ref_frame.atoms[i].element.clone()); }
        set.into_iter().collect()
    };

    let n_bins = ((params.r_max - params.r_min) / params.dr).ceil() as usize;
    if n_bins == 0 { return None; }

    // origin 起始帧：访问 cart[p] 和 cart[p+tau]，均需 < n_steps
    let p_values: Vec<usize> = (0..)
        .map(|k: usize| k * shift)
        .take_while(|&p| p + tau < n_steps)
        .collect();
    let n_origins = p_values.len();
    if n_origins == 0 { return None; }

    // 构建绝对 Cartesian 坐标数组（周期系统经过 unwrap）
    let cart: Vec<Vec<[f64; 3]>> = if ref_frame.cell.is_some() {
        // 周期系统：分数坐标 unwrap 后转回绝对坐标
        if traj.frames.iter().any(|f| f.cell.is_none()) { return None; }

        let mut frac: Vec<Vec<[f64; 3]>> = traj.frames.iter().map(|frame| {
            let cell = frame.cell.as_ref().unwrap();
            atom_indices.iter().map(|&i| {
                let f = cell.cartesian_to_fractional(frame.atoms[i].position);
                [f.x, f.y, f.z]
            }).collect()
        }).collect();
        unwrap_frac(&mut frac);

        // 用各帧自身的 cell 矩阵将分数坐标转回绝对坐标（NPT 自然支持）
        frac.iter().zip(traj.frames.iter())
            .map(|(step_frac, frame)| {
                let m = &frame.cell.as_ref().unwrap().matrix;
                step_frac.iter().map(|&[fa, fb, fc]| [
                    fa * m[(0,0)] + fb * m[(1,0)] + fc * m[(2,0)],
                    fa * m[(0,1)] + fb * m[(1,1)] + fc * m[(2,1)],
                    fa * m[(0,2)] + fb * m[(1,2)] + fc * m[(2,2)],
                ]).collect()
            }).collect()
    } else {
        // 非周期：直接使用 Cartesian 坐标，无需 unwrap
        traj.frames.iter().map(|frame| {
            atom_indices.iter().map(|&i| {
                let p = &frame.atoms[i].position;
                [p.x, p.y, p.z]
            }).collect()
        }).collect()
    };

    // 并行计算各 origin 的局部直方图，reduce 合并
    let hist: Vec<u64> = p_values.par_iter()
        .map(|&p| {
            let mut local = vec![0u64; n_bins];
            for j in 0..n_atoms {
                let dx = cart[p + tau][j][0] - cart[p][j][0];
                let dy = cart[p + tau][j][1] - cart[p][j][1];
                let dz = cart[p + tau][j][2] - cart[p][j][2];
                let r = (dx*dx + dy*dy + dz*dz).sqrt();
                if r >= params.r_min && r < params.r_max {
                    let bin = ((r - params.r_min) / params.dr) as usize;
                    if bin < n_bins { local[bin] += 1; }
                }
            }
            local
        })
        .reduce(
            || vec![0u64; n_bins],
            |mut a, b| { for i in 0..n_bins { a[i] += b[i]; } a },
        );

    // 归一化：bin 计数 / (n_origins × n_atoms)
    let norm = 1.0 / (n_origins as f64 * n_atoms as f64);
    let r: Vec<f64> = (0..n_bins)
        .map(|i| params.r_min + (i as f64 + 0.5) * params.dr)
        .collect();
    let gs: Vec<f64> = hist.iter().map(|&c| c as f64 * norm).collect();

    Some(VanHoveResult {
        r, gs,
        tau_frames: tau,
        time: tau as f64 * params.dt,
        n_atoms, n_origins,
        params: params.clone(),
        elements,
    })
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

/// Write van Hove self-correlation function to a tab-separated text file (`.vanhove`).
///
/// Two columns: `r[Å]` and `Gs(r,tau)`.  Values are normalised so that
/// the discrete sum over all bins equals 1.0.
pub fn write_vanhove(result: &VanHoveResult, path: &str) -> std::io::Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(path)?);

    writeln!(w, "# ferro v{}", VERSION)?;
    writeln!(w, "# van Hove Self-Correlation Function Gs(r, tau)")?;
    writeln!(w, "# {}", "-".repeat(60))?;
    writeln!(w, "# tau     = {} frames", result.tau_frames)?;
    writeln!(w, "# time    = {} fs", result.time)?;
    writeln!(w, "# shift   = {} frames", result.params.shift)?;
    writeln!(w, "# dr      = {} Ang", result.params.dr)?;
    writeln!(w, "# r_min   = {} Ang", result.params.r_min)?;
    writeln!(w, "# r_max   = {} Ang", result.params.r_max)?;
    writeln!(w, "# atoms   = {}", result.n_atoms)?;
    writeln!(w, "# origins = {}", result.n_origins)?;
    write!(w,   "# elements:")?;
    for e in &result.elements { write!(w, " {}", e)?; }
    writeln!(w)?;
    writeln!(w, "# {}", "-".repeat(60))?;
    writeln!(w, "# r[Ang]\tGs(r,tau)")?;

    for (r, gs) in result.r.iter().zip(result.gs.iter()) {
        writeln!(w, "{:.6e}\t{:.6e}", r, gs)?;
    }
    Ok(())
}

// ─── 测试 ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    /// 静态轨迹：n 帧完全相同，位移恒为 0
    fn make_static_traj(n: usize) -> Trajectory {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut ref_frame = Frame::with_cell(cell, [true; 3]);
        ref_frame.add_atom(Atom::new("Fe", Vector3::new(1.0, 1.0, 1.0)));
        ref_frame.add_atom(Atom::new("Fe", Vector3::new(3.0, 3.0, 3.0)));
        let mut traj = Trajectory::new();
        for _ in 0..n { traj.add_frame(ref_frame.clone()); }
        traj
    }

    /// 匀速轨迹：单个 Fe 原子每帧沿 x 移动 v Å（大盒子，不跨界）
    fn make_linear_traj(a: f64, v: f64, n: usize) -> Trajectory {
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
    fn test_vanhove_static_zero_displacement() {
        // 静态轨迹，所有位移 = 0 → 全部计数落在 bin 0（r ∈ [0, dr)）
        let traj = make_static_traj(10);
        let params = VanHoveParams { tau: Some(3), shift: 1, dr: 0.01, ..Default::default() };
        let res = calc_vanhove(&traj, &params).unwrap();

        // bin 0 应包含所有计数（归一化后 = 1.0）
        assert!((res.gs[0] - 1.0).abs() < 1e-10,
            "static: gs[0] should be 1.0, got {}", res.gs[0]);
        // 其余 bin 应为 0
        let rest: f64 = res.gs[1..].iter().sum();
        assert!(rest.abs() < 1e-10, "static: non-zero bins beyond 0: {}", rest);
        // 归一化检查：sum(gs) == 1
        let total: f64 = res.gs.iter().sum();
        assert!((total - 1.0).abs() < 1e-10, "normalization: sum={}", total);
    }

    #[test]
    fn test_vanhove_known_displacement() {
        // 单原子以 v=1.0 Å/step 匀速运动，tau=3 → 所有原点的位移均为 3.0 Å
        // 使用 dr=0.5（精确可表示的 2 的幂次分量）避免浮点边界问题
        let a = 50.0;
        let v = 1.0_f64;
        let tau = 3usize;
        let dr = 0.5;
        let traj = make_linear_traj(a, v, 20);
        let params = VanHoveParams {
            tau: Some(tau), shift: 1, dt: 1.0,
            r_min: 0.0, r_max: 10.0, dr,
            elements: None,
        };
        let res = calc_vanhove(&traj, &params).unwrap();

        // 归一化检查
        let sum: f64 = res.gs.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10, "normalization: sum = {}", sum);

        // 所有概率质量应集中在 r = 3.0 Å 附近（单一 bin）
        let expected_r = tau as f64 * v;
        let mean_r: f64 = res.r.iter().zip(res.gs.iter()).map(|(&r, &g)| r * g).sum();
        assert!((mean_r - expected_r).abs() < dr,
            "mean displacement {:.4} ≈ expected {:.1}", mean_r, expected_r);
    }

    #[test]
    fn test_vanhove_unwrap_across_boundary() {
        // 原子从接近边界处出发，每步 0.6 Å 跨越周期边界
        // unwrap 后位移应等于 tau * v
        let a = 5.0;
        let v = 0.6_f64;
        let x0 = 4.8_f64;
        let tau = 2usize;
        let dr = 0.01;
        let n = 10;
        let cell = Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for i in 0..n {
            let x_raw = x0 + i as f64 * v;
            let x_wrapped = x_raw - (x_raw / a).floor() * a;
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("Li", Vector3::new(x_wrapped, 0.0, 0.0)));
            traj.add_frame(frame);
        }

        let params = VanHoveParams {
            tau: Some(tau), shift: 1, dt: 1.0,
            r_min: 0.0, r_max: 5.0, dr,
            elements: None,
        };
        let res = calc_vanhove(&traj, &params).unwrap();

        // 期望位移 = tau * v = 1.2 Å
        // 用加权均值验证（避免浮点 bin 边界问题）
        let expected_r = tau as f64 * v;
        let sum: f64 = res.gs.iter().sum();
        assert!((sum - 1.0).abs() < 1e-10, "normalization: sum = {}", sum);
        let mean_r: f64 = res.r.iter().zip(res.gs.iter()).map(|(&r, &g)| r * g).sum();
        assert!((mean_r - expected_r).abs() < dr,
            "unwrap mean displacement {:.4} ≈ {:.4}", mean_r, expected_r);
    }

    #[test]
    fn test_vanhove_element_filter() {
        // 含 Fe 和 Li 的轨迹；只计算 Li，Fe 不动，Li 以 0.5 Å/step 运动
        let a = 20.0;
        let v = 0.5_f64;
        let tau = 2usize;
        let dr = 0.01;
        let cell = Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for i in 0..10 {
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("Fe", Vector3::new(1.0, 1.0, 1.0)));
            frame.add_atom(Atom::new("Li", Vector3::new(i as f64 * v, 0.0, 0.0)));
            traj.add_frame(frame);
        }

        let params = VanHoveParams {
            tau: Some(tau), shift: 1, dt: 1.0,
            r_min: 0.0, r_max: 5.0, dr,
            elements: Some(vec!["Li".to_string()]),
        };
        let res = calc_vanhove(&traj, &params).unwrap();

        assert_eq!(res.n_atoms, 1, "should include only Li");
        assert!(res.elements.contains(&"Li".to_string()));
        assert!(!res.elements.contains(&"Fe".to_string()));

        // Li 位移 = tau * v = 1.0 Å；用加权均值验证
        let expected_r = tau as f64 * v;
        let mean_r: f64 = res.r.iter().zip(res.gs.iter()).map(|(&r, &g)| r * g).sum();
        assert!((mean_r - expected_r).abs() < dr,
            "Li mean displacement {:.4} ≈ {:.4}", mean_r, expected_r);
    }

    #[test]
    fn test_vanhove_normalization() {
        // 任意轨迹的 gs 之和应等于 1（PMF 归一化）
        let traj = make_static_traj(8);
        let params = VanHoveParams { tau: Some(2), ..Default::default() };
        let res = calc_vanhove(&traj, &params).unwrap();
        let total: f64 = res.gs.iter().sum();
        assert!((total - 1.0).abs() < 1e-10, "PMF sum = {}", total);
    }

    #[test]
    fn test_vanhove_n_origins() {
        // n=10, tau=3, shift=1 → p ∈ {0,1,...,6} → 7 origins
        let traj = make_static_traj(10);
        let params = VanHoveParams { tau: Some(3), shift: 1, ..Default::default() };
        let res = calc_vanhove(&traj, &params).unwrap();
        // take_while(p+3 < 10) → p ∈ {0,1,2,3,4,5,6} → 7 origins
        assert_eq!(res.n_origins, 7);
    }

    #[test]
    fn test_vanhove_time_value() {
        let traj = make_static_traj(10);
        let params = VanHoveParams { tau: Some(5), dt: 2.0, ..Default::default() };
        let res = calc_vanhove(&traj, &params).unwrap();
        assert_eq!(res.tau_frames, 5);
        assert!((res.time - 10.0).abs() < 1e-10, "time = tau * dt = 10 fs");
    }

    #[test]
    fn test_write_vanhove() {
        use std::io::Read;
        let traj = make_linear_traj(50.0, 0.1, 10);
        let params = VanHoveParams { tau: Some(5), ..Default::default() };
        let res = calc_vanhove(&traj, &params).unwrap();
        let path = "/tmp/test_ferro.vanhove";
        write_vanhove(&res, path).expect("write_vanhove failed");

        let mut content = String::new();
        std::fs::File::open(path).unwrap().read_to_string(&mut content).unwrap();
        assert!(content.starts_with("# ferro v"), "missing version header");
        assert!(content.contains("van Hove"), "missing title");
        assert!(content.contains("# r[Ang]"), "missing column header");
    }
}
