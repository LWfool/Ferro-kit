//! Van Hove self-correlation function Gs(r, τ) calculation and output.
//!
//! Gs(r, τ) gives the probability distribution of atomic displacement distances over time interval τ.
//! The output is the radial probability density P(r) = 4πr²·Gs(r, τ) [Å⁻¹], so that
//! ∫P dr = 1 − (fraction of displacements outside [r_min, r_max)).
//!
//! Workflow: `calc_vanhove` → `VanHoveResult::to_tables`.
//! Algorithm follows code1/vanhove.c (`EstimateVanHove`), except for the unwrapping:
//!   1. TOR unwrapping (`util::unwrapped_series`, as `msd`); code1 used the lattice
//!      view (unwrapped fraction × current box), which turns NPT box breathing into
//!      displacement.
//!   2. For each time origin p, compute |r(p+τ) − r(p)| and accumulate into a histogram.
//!   3. Normalise: p_r[i] = count[i] / (n_origins × n_atoms × dr).
//!
//! Parallelism: per time-origin par_iter; each origin computed independently then reduced.

use crate::check;
use super::correlate::resolve_max_lag;
use super::util::{check_frame_spacing, select_atoms, unwrapped_series};
use rayon::prelude::*;
use std::collections::BTreeSet;
use ferro_core::error::ChemError;
use ferro_core::{Table, Trajectory};

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Parameters for van Hove self-correlation function calculation.
#[derive(Debug, Clone)]
pub struct VanHoveParams {
    /// Time lag in frames (`None` = `n_frames / 2`)
    pub tau: Option<usize>,
    /// Time shift between origins in frames (default: 1)
    pub shift: usize,
    /// Time step per frame \[fs\]. No default: `Default` leaves it NaN and the entry rejects that
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

impl VanHoveParams {
    /// Value ranges that do not depend on the trajectory; the CLI calls this
    /// before reading the first file, and the `calc_*` entry calls it again.
    pub fn validate(&self) -> ferro_core::Result<()> {
        check::frame_interval(self.dt)?;
        check::at_least_one_if_given("tau", self.tau)?;
        check::at_least_one("shift", self.shift)?;
        check::non_negative("r-min", self.r_min)?;
        check::positive("dr", self.dr)?;
        check::ordered("r-min", self.r_min, "r-max", self.r_max)?;
        check::holds_a_bin("r range", self.r_min, self.r_max, "dr", self.dr)
    }
}

impl Default for VanHoveParams {
    fn default() -> Self {
        VanHoveParams {
            tau: None, shift: 1, dt: f64::NAN,
            r_min: 0.0, r_max: 10.0, dr: 0.01,
            elements: None,
        }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Result of a van Hove self-correlation calculation.
///
/// `p_r` is the radial probability density of atomic displacements,
/// P(r) = 4πr²·Gs(r, τ) in Å⁻¹: `p_r[i] · dr` is the fraction of (atom, origin)
/// pairs whose displacement fell in bin `i`.  It is independent of `dr`, so it can be
/// compared against the Gaussian reference directly, and Σ p_r·dr = 1 − `outside_fraction`.
#[derive(Debug, Clone)]
pub struct VanHoveResult {
    /// Bin-centre displacement values \[Å\]
    pub r: Vec<f64>,
    /// P(r) = 4πr²·Gs(r, τ) \[Å⁻¹\]
    pub p_r: Vec<f64>,
    /// Fraction of (atom, origin) displacements outside \[r_min, r_max): missing from `p_r`.
    /// A long τ in a liquid pushes the tail past `r_max`; this says how much was cut.
    pub outside_fraction: f64,
    /// Actual τ used \[frames\]
    pub tau_frames: usize,
    /// Actual τ in physical time \[fs\]
    pub time: f64,
    /// Frames in the trajectory
    pub n_frames: usize,
    /// Number of atoms included
    pub n_atoms: usize,
    /// Number of time origins averaged
    pub n_origins: usize,
    pub params: VanHoveParams,
    /// Element types included, sorted alphabetically
    pub elements: Vec<String>,
}

// ─── 计算 ────────────────────────────────────────────────────────────────────

/// Compute the van Hove self-correlation function Gs(r, τ) with time-shift averaging.
///
/// The algorithm matches code1/vanhove.c `EstimateVanHove` except for step 1:
/// 1. Unwrap periodic trajectories with TOR (`util::unwrapped_series`); keep
///    Cartesian positions as they are for non-periodic ones.
/// 2. For each time origin p (spaced by `shift`), compute the scalar displacement
///    `|r(p+τ) − r(p)|` for every selected atom and accumulate into a histogram.
/// 3. Normalise: `p_r[i] = count[i] / (n_origins × n_atoms × dr)`; displacements outside
///    `[r_min, r_max)` stay in the denominator and are reported as `outside_fraction`.
///
/// Returns `Err`, each with its own message, if:
/// - a parameter is out of range ([`VanHoveParams::validate`])
/// - the trajectory has fewer than 2 frames, or its frames hold different atoms
/// - no atoms match the element filter
/// - some frames have a cell and others do not
///
/// - `τ` is not within `1 ..= n_frames − 1` (as `--max-lag` is checked by the other
///   time-correlation analyses; it used to be clamped to `n_frames − 1` without a word)
pub fn calc_vanhove(traj: &Trajectory, params: &VanHoveParams) -> ferro_core::Result<VanHoveResult> {
    params.validate()?;
    check_frame_spacing(traj)?;
    let n_steps = traj.n_frames();
    if n_steps < 2 {
        return Err(ChemError::ValidationError(format!("trajectory requires at least 2 frames, got {n_steps}")));
    }
    check::same_atoms(traj)?;
    let cells = check::periodic_cells(traj, None)?;

    // 按第一帧筛选参与计算的原子下标
    let ref_frame = &traj.frames[0];
    let atom_indices = select_atoms(ref_frame, params.elements.as_deref())?;
    let n_atoms = atom_indices.len();

    // The same rule as max-lag: None takes half the trajectory (n_steps - 1 would leave a
    // single time origin), and τ must lie in 1 ..= n_steps − 1
    let tau = resolve_max_lag("tau", n_steps, params.tau)?;
    let shift = params.shift.max(1);

    // 收集元素列表（用于输出文件头）
    let elements: Vec<String> = {
        let mut set = BTreeSet::new();
        for &i in &atom_indices { set.insert(ref_frame.element(i).to_string()); }
        set.into_iter().collect()
    };

    let (n_bins, dr) = super::util::bins(params.r_min, params.r_max, params.dr);
    if n_bins == 0 {
        // validate 已要求 dr 有限且 > 0、r_max > r_min，这里只是防御
        return Err(ChemError::ValidationError(format!("dr = {} leaves no bin", params.dr)));
    }
    // Everything downstream (histogram, normalisation, header) uses the width actually used
    let params = &VanHoveParams { dr, ..params.clone() };

    // origin 起始帧：访问 cart[p] 和 cart[p+tau]，均需 < n_steps
    let p_values: Vec<usize> = (0..)
        .map(|k: usize| k * shift)
        .take_while(|&p| p + tau < n_steps)
        .collect();
    let n_origins = p_values.len();
    if n_origins == 0 {
        // tau <= n_steps - 1 and shift >= 1 make p = 0 always fit; kept against later changes
        return Err(ChemError::ValidationError(format!("no time origin fits tau = {tau} in {n_steps} frames")));
    }

    // Continuous Cartesian series (TOR-unwrapped when periodic), [atom][axis][frame]
    let series = unwrapped_series(traj, &cells, &atom_indices);

    // 并行计算各 origin 的局部直方图，reduce 合并
    let hist: Vec<u64> = p_values.par_iter()
        .map(|&p| {
            let mut local = vec![0u64; n_bins];
            for [x, y, z] in &series {
                let dx = x[p + tau] - x[p];
                let dy = y[p + tau] - y[p];
                let dz = z[p + tau] - z[p];
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

    // 归一化成概率密度：计数 / (样本数 × dr)。分母含落在范围外的样本 —— 它们是分布的
    // 一部分，只是没画出来；除以落在范围内的样本数会把被截的尾巴摊进可见部分。
    // 以前不除 dr，输出是每个 bin 的概率质量，数值随 dr 变，没法和 4πr²Gs 的高斯参考式比
    let n_samples = n_origins as f64 * n_atoms as f64;
    let norm = 1.0 / (n_samples * params.dr);
    let r: Vec<f64> = (0..n_bins)
        .map(|i| params.r_min + (i as f64 + 0.5) * params.dr)
        .collect();
    let p_r: Vec<f64> = hist.iter().map(|&c| c as f64 * norm).collect();
    let inside: u64 = hist.iter().sum();
    let outside_fraction = 1.0 - inside as f64 / n_samples;

    Ok(VanHoveResult {
        r, p_r, outside_fraction,
        tau_frames: tau,
        time: tau as f64 * params.dt,
        n_frames: n_steps, n_atoms, n_origins,
        params: params.clone(),
        elements,
    })
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

impl VanHoveResult {
    /// Projects the result into the table the writers consume.
    ///
    /// Van Hove self-correlation at one tau: `r, p_r` (P(r) = 4πr²Gs, Å⁻¹). The tau itself is metadata for now.
    /// The `file` column is added by the caller when stacking several inputs
    /// (see `ferro_core::Table::concat_union`).
    pub fn to_tables(&self) -> Vec<(String, Table)> {
        let mut t = Table::new();
        t.push_num("r", self.r.clone())
            .push_num("p_r", self.p_r.clone());
        vec![("vanhove".to_string(), t)]
    }

    /// Parameter block for the comment header above the data.
    pub fn meta_lines(&self) -> Vec<String> {
        vec![
            format!("tau     = {} frames", self.tau_frames),
            format!("time    = {} fs", self.time),
            format!("shift   = {} frames", self.params.shift),
            format!("dr      = {} Ang", self.params.dr),
            format!("r_min   = {} Ang", self.params.r_min),
            format!("r_max   = {} Ang", self.params.r_max),
            format!("atoms   = {}", self.n_atoms),
            format!("origins = {}", self.n_origins),
            format!("elements: {}", self.elements.join(" ")),
        ]
    }
}

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

    /// A range narrower than one step is refused up front, as gr does (review C-11);
    /// `bins` would otherwise make one bin narrower than asked
    #[test]
    fn test_range_narrower_than_one_bin_is_an_error() {
        let p = VanHoveParams { r_min: 0.0, r_max: 0.05, dr: 0.1, dt: 1.0, ..Default::default() };
        let err = p.validate().unwrap_err().to_string();
        assert!(err.contains("r range [0, 0.05] is narrower than one dr = 0.1 bin"), "{err}");
    }

    #[test]
    fn test_vanhove_static_zero_displacement() {
        // 静态轨迹，所有位移 = 0 → 全部计数落在 bin 0（r ∈ [0, dr)）
        let traj = make_static_traj(10);
        let params = VanHoveParams { tau: Some(3), shift: 1, dr: 0.01, dt: 1.0, ..Default::default() };
        let res = calc_vanhove(&traj, &params).unwrap();

        // bin 0 应包含所有计数：密度 = 1/dr
        assert!((res.p_r[0] * 0.01 - 1.0).abs() < 1e-10,
            "static: p_r[0]·dr should be 1.0, got {}", res.p_r[0] * 0.01);
        // 其余 bin 应为 0
        let rest: f64 = res.p_r[1..].iter().sum();
        assert!(rest.abs() < 1e-10, "static: non-zero bins beyond 0: {}", rest);
        assert_eq!(res.outside_fraction, 0.0);
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
        let sum: f64 = res.p_r.iter().sum::<f64>() * dr;
        assert!((sum - 1.0).abs() < 1e-10, "normalization: sum = {}", sum);

        // 所有概率质量应集中在 r = 3.0 Å 附近（单一 bin）
        let expected_r = tau as f64 * v;
        let mean_r: f64 = res.r.iter().zip(res.p_r.iter()).map(|(&r, &p)| r * p * dr).sum();
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
        let sum: f64 = res.p_r.iter().sum::<f64>() * dr;
        assert!((sum - 1.0).abs() < 1e-10, "normalization: sum = {}", sum);
        let mean_r: f64 = res.r.iter().zip(res.p_r.iter()).map(|(&r, &p)| r * p * dr).sum();
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
        let mean_r: f64 = res.r.iter().zip(res.p_r.iter()).map(|(&r, &p)| r * p * dr).sum();
        assert!((mean_r - expected_r).abs() < dr,
            "Li mean displacement {:.4} ≈ {:.4}", mean_r, expected_r);
    }

    #[test]
    fn test_vanhove_density_is_dr_independent_and_reports_cut_tail() {
        // 两个原子：一个不动，一个每帧走 1.5 Å；tau=2 → 位移 0 与 3 Å 各占一半
        let cell = Cell::from_lengths_angles(50.0, 50.0, 50.0, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for i in 0..6 {
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("Li", Vector3::new(1.0, 1.0, 1.0)));
            frame.add_atom(Atom::new("Li", Vector3::new(1.0 + 1.5 * i as f64, 5.0, 5.0)));
            traj.add_frame(frame);
        }
        let run = |r_max: f64, dr: f64| calc_vanhove(&traj, &VanHoveParams {
            tau: Some(2), r_max, dr, dt: 1.0, ..Default::default()
        }).unwrap();

        // 范围盖住全部位移：积分 = 1，无溢出；换 dr 后 r=0 处的密度按 1/dr 变，积分不变
        for dr in [0.5, 0.25] {
            let res = run(10.0, dr);
            let integral: f64 = res.p_r.iter().sum::<f64>() * dr;
            assert!((integral - 1.0).abs() < 1e-10, "dr={dr}: ∫P dr 应为 1，实际 {integral}");
            assert_eq!(res.outside_fraction, 0.0);
        }
        // r_max=2 截掉 3 Å 那一半：积分 0.5，溢出比例 0.5（以前静默报 Σ=0.5 冒充归一）
        let res = run(2.0, 0.5);
        let integral: f64 = res.p_r.iter().sum::<f64>() * 0.5;
        assert!((integral - 0.5).abs() < 1e-10, "截尾后 ∫P dr 应为 0.5，实际 {integral}");
        assert!((res.outside_fraction - 0.5).abs() < 1e-10, "溢出比例应为 0.5，实际 {}", res.outside_fraction);
    }

    #[test]
    fn test_vanhove_n_origins() {
        // n=10, tau=3, shift=1 → p ∈ {0,1,...,6} → 7 origins
        let traj = make_static_traj(10);
        let params = VanHoveParams { tau: Some(3), shift: 1, dt: 1.0, ..Default::default() };
        let res = calc_vanhove(&traj, &params).unwrap();
        // take_while(p+3 < 10) → p ∈ {0,1,2,3,4,5,6} → 7 origins
        assert_eq!(res.n_origins, 7);
        // [inputs] 的 frames 列读它；旧代码填的是 bin 数，见审查 B-4
        assert_eq!(res.n_frames, 10);
    }

    #[test]
    fn test_vanhove_default_tau_is_half_trajectory() {
        // 不给 tau 时取 n/2（同 correlate.rs 的 max_lag 默认与 gmx 的 -acflen），
        // 而不是 n-1 —— 后者只剩 1 个原点
        let traj = make_static_traj(10);
        let res = calc_vanhove(&traj, &VanHoveParams { dt: 1.0, ..Default::default() }).unwrap();
        assert_eq!(res.tau_frames, 5, "默认 tau 应为 10/2");
        assert_eq!(res.n_origins, 5, "p ∈ {{0..4}} 满足 p+5 < 10");
        // 两帧时 n/2 = 1，仍有 1 个原点
        let res = calc_vanhove(&make_static_traj(2), &VanHoveParams { dt: 1.0, ..Default::default() }).unwrap();
        assert_eq!((res.tau_frames, res.n_origins), (1, 1));
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
    fn test_to_tables_columns() {
        let traj = make_static_traj(10);
        let res = calc_vanhove(&traj, &VanHoveParams { dt: 1.0, ..Default::default() }).unwrap();
        let (name, t) = res.to_tables().remove(0);
        assert_eq!(name, "vanhove");
        assert_eq!(t.names(), vec!["r", "p_r"]);
        assert_eq!(t.n_rows(), res.r.len());
        assert!(t.validate().is_ok());
        // tau is metadata, not a column, for now (dev/plan/analysis.md, shelved items)
        assert!(res.meta_lines().join("\n").contains("tau"));
    }

    /// `--r-max 1.0 --dr 0.07` used to give 15 bins of 0.07, the last one [0.98, 1.05)
    /// collecting only up to 1.0 but normalised as a full bin: 3.5 times too low.
    #[test]
    fn test_uneven_range_has_no_partial_last_bin() {
        // A displacement of 0.99 Å falls in the last bin
        let traj = make_linear_traj(50.0, 0.33, 20);
        let params = VanHoveParams {
            tau: Some(3), shift: 1, dt: 1.0, r_min: 0.0, r_max: 1.0, dr: 0.07, elements: None,
        };
        let res = calc_vanhove(&traj, &params).unwrap();
        assert_eq!(res.p_r.len(), 15);
        let w = res.params.dr;
        assert!((w * 15.0 - 1.0).abs() < 1e-12, "width {w}");
        assert!(*res.r.last().unwrap() < 1.0, "last centre inside r_max");
        // Every sample is in the last bin: density × width = 1
        assert!((res.p_r[14] * w - 1.0).abs() < 1e-12, "{}", res.p_r[14] * w);
        assert!(res.meta_lines().join("\n").contains(&format!("dr      = {w} Ang")));
    }

    /// NPT: the box breathing must not show up as displacement. The atom crosses a face
    /// (x: 9.5 → 10.5, wrapped to 0.5), then sits still in fractional coordinates while
    /// the box goes 10 → 11 → 10. The lattice view (unwrapped fraction × current box)
    /// reads 1.05 × 11 = 11.55 and moves the atom 1.05 Å per breath; TOR adds only the
    /// minimum-image steps, 1.0, +0.05, −0.05. Same frames as msd's
    /// `test_tor_does_not_follow_lattice_scaling`.
    #[test]
    fn test_tor_does_not_follow_lattice_scaling() {
        let frames = [(10.0, 0.95), (10.0, 0.05), (11.0, 0.05), (10.0, 0.05)];
        let mut traj = Trajectory::new();
        for (len, fx) in frames {
            let cell = Cell::from_lengths_angles(len, len, len, 90.0, 90.0, 90.0).unwrap();
            let mut fr = Frame::with_cell(cell.clone(), [true; 3]);
            // y = z = 0 so the breathing box does not move the atom along them
            fr.add_atom(Atom::new("Na", cell.fractional_to_cartesian(Vector3::new(fx, 0.0, 0.0))));
            traj.add_frame(fr);
        }
        let params = VanHoveParams { tau: Some(1), dt: 1.0, r_max: 2.0, dr: 0.01, ..Default::default() };
        let res = calc_vanhove(&traj, &params).unwrap();
        // Probability mass with bin centre in [lo, hi): the 0.05 Å steps sit on a bin edge
        let mass = |lo: f64, hi: f64| -> f64 {
            res.r.iter().zip(&res.p_r).filter(|(r, _)| (lo..hi).contains(*r)).map(|(_, p)| p * res.params.dr).sum()
        };
        // Three origins: 1.0 once, 0.05 twice; nothing at the lattice view's 1.05
        assert!((mass(0.03, 0.07) - 2.0 / 3.0).abs() < 1e-9, "mass near 0.05 Å: {}", mass(0.03, 0.07));
        assert!((mass(0.98, 1.02) - 1.0 / 3.0).abs() < 1e-9, "mass near 1.0 Å: {}", mass(0.98, 1.02));
        assert_eq!(mass(1.03, 1.07), 0.0, "lattice-view displacement 1.05 Å must be empty");
    }

    /// `--tau` past the last frame used to be clamped silently; the other time
    /// correlations refuse `--max-lag` past it.
    #[test]
    fn test_tau_beyond_the_trajectory_is_an_error() {
        let traj = make_linear_traj(50.0, 1.0, 10);
        let p = |tau| VanHoveParams { tau: Some(tau), shift: 1, dt: 1.0, ..Default::default() };
        let err = calc_vanhove(&traj, &p(10)).unwrap_err().to_string();
        assert!(err.contains("tau must be within 1..=9 for 10 frames, got 10"), "{err}");
        assert_eq!(calc_vanhove(&traj, &p(9)).unwrap().tau_frames, 9);
    }
}
