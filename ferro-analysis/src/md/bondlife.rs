//! Bond lifetimes: intermittent and continuous bond correlation functions, and
//! per-frame bond formation / breaking events.
//!
//! A *candidate* is a center–neighbor pair that comes within `r_bond` in at least one
//! frame. Its existence function h(t) ∈ {0, 1} follows a hysteresis rule: a free pair
//! bonds when `r ≤ r_bond`; a bonded pair stays bonded while `r ≤ r_break`
//! (`r_break ≥ r_bond`, equal by default — then it is the single-threshold contact
//! criterion of `gmx hbond -contact`).
//!
//! - Intermittent (Luzar–Chandler; `gmx hbond -ac`):
//!   `C_I(m) = Σ h(t)h(t+m) / Σ h(t)`, sums over all pairs and all origins `t ≤ N−1−m`.
//! - Continuous (Rapaport; MDAnalysis `autocorrelation`):
//!   `S_C(m)` = the same ratio with the numerator counting only origins whose bond
//!   exists at **every** frame of `[t, t+m]`. Gaps of at most `intermittency` frames
//!   inside a bonded stretch are filled first (MDAnalysis `intermittency`).
//! - Events: bonds present, formed and broken per frame, on the gap-filled h.
//!
//! Parallelism: per frame for distances, per candidate for the correlations.

use rayon::prelude::*;
use nalgebra::Matrix3;
use ferro_core::{Table, Trajectory};
use ferro_core::error::ChemError;

use super::correlate::{cumulative_trapezoid, resolve_max_lag, AutocorrPlan};

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Parameters for the bond-lifetime analysis.
#[derive(Debug, Clone)]
pub struct BondLifeParams {
    /// Element of the center atoms (e.g. `"Si"`)
    pub center: String,
    /// Element of the neighbor atoms (e.g. `"O"`); may equal `center`
    pub neighbor: String,
    /// A free pair becomes bonded at `r ≤ r_bond` \[Å\]
    pub r_bond: f64,
    /// A bonded pair stays bonded while `r ≤ r_break` \[Å\] (`None` = `r_bond`)
    pub r_break: Option<f64>,
    /// Breaks of at most this many frames inside a bonded stretch are filled
    /// (continuous function and events only; default 0)
    pub intermittency: usize,
    /// Longest lag in frames, `1 ..= n_frames − 1` (`None` = `n_frames / 2`)
    pub max_lag: Option<usize>,
    /// Time between stored frames \[fs\]
    pub dt: f64,
}

impl Default for BondLifeParams {
    fn default() -> Self {
        BondLifeParams {
            center: "Si".into(),
            neighbor: "O".into(),
            r_bond: 2.0,
            r_break: None,
            intermittency: 0,
            max_lag: None,
            dt: 1.0,
        }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Result of the bond-lifetime analysis.
#[derive(Debug, Clone)]
pub struct BondLifeResult {
    /// Lag-time axis \[fs\], `m = 0 ..= max_lag`
    pub time: Vec<f64>,
    /// Intermittent bond correlation `C_I` (NaN at a lag with no bonded origin)
    pub c_int: Vec<f64>,
    /// Continuous bond correlation `S_C`
    pub s_cont: Vec<f64>,
    /// Frame times \[fs\], all frames
    pub frame_time: Vec<f64>,
    /// Bonds present in each frame (gap-filled)
    pub n_bonds: Vec<f64>,
    /// Bonds formed since the previous frame (NaN at frame 0)
    pub formed: Vec<f64>,
    /// Bonds broken since the previous frame (NaN at frame 0)
    pub broken: Vec<f64>,
    /// Frames in the trajectory
    pub n_frames: usize,
    /// Atoms of the center element
    pub n_centers: usize,
    /// Pairs within `r_bond` in at least one frame
    pub n_candidates: usize,
    /// Time origins at the longest lag, `n_frames − max_lag`
    pub min_origins: usize,
    pub params: BondLifeParams,
}

// ─── 内部辅助 ─────────────────────────────────────────────────────────────────

/// Fills breaks of at most `k` frames that lie **between** two bonded stretches.
/// A break running into either end of the trajectory is left alone — it cannot be
/// known to close.
fn fill_gaps(h: &mut [bool], k: usize) {
    if k == 0 {
        return;
    }
    let mut t = 0;
    while t < h.len() {
        if h[t] {
            t += 1;
            continue;
        }
        let start = t;
        while t < h.len() && !h[t] {
            t += 1;
        }
        if start > 0 && t < h.len() && t - start <= k {
            h[start..t].iter_mut().for_each(|v| *v = true);
        }
    }
}

/// Existence function with hysteresis: bond at `d ≤ r_bond`, keep while `d ≤ r_break`.
fn existence(d: &[f64], r_bond: f64, r_break: f64) -> Vec<bool> {
    let mut on = false;
    d.iter().map(|&r| {
        on = if on { r <= r_break } else { r <= r_bond };
        on
    }).collect()
}

/// Correlation time as the running trapezoidal integral's last value, and as the
/// lag where the function first drops below 1/e (linear interpolation; NaN if it
/// never does within the lag axis).
pub fn lifetimes(c: &[f64], dt: f64) -> (f64, f64) {
    let integral = *cumulative_trapezoid(c, dt).last().unwrap_or(&f64::NAN);
    let target = (-1.0f64).exp();
    let cross = c.windows(2).enumerate().find(|(_, w)| w[0] >= target && w[1] < target)
        .map(|(m, w)| (m as f64 + (w[0] - target) / (w[0] - w[1])) * dt)
        .unwrap_or(f64::NAN);
    (integral, cross)
}

// ─── 计算 ────────────────────────────────────────────────────────────────────

/// Runs the bond-lifetime analysis.
///
/// Returns `Err` if the trajectory has fewer than 2 frames, `max_lag` is outside
/// `1 ..= n_frames − 1`, the thresholds are not `0 < r_bond ≤ r_break`, an element
/// is absent, or no pair comes within `r_bond` in any frame.
pub fn calc_bondlife(traj: &Trajectory, params: &BondLifeParams) -> ferro_core::Result<BondLifeResult> {
    let n_frames = traj.n_frames();
    let max_lag = resolve_max_lag(n_frames, params.max_lag)?;
    let r_bond = params.r_bond;
    let r_break = params.r_break.unwrap_or(r_bond);
    if !(r_bond > 0.0 && r_break >= r_bond) {
        return Err(ChemError::ValidationError(format!(
            "need 0 < r_bond <= r_break, got r_bond = {r_bond}, r_break = {r_break}")));
    }

    let ref_frame = &traj.frames[0];
    let of = |el: &str| -> Vec<usize> {
        ref_frame.atoms.iter().enumerate().filter(|(_, a)| a.element == el).map(|(i, _)| i).collect()
    };
    let (centers, neighbors) = (of(&params.center), of(&params.neighbor));
    for (list, el) in [(&centers, &params.center), (&neighbors, &params.neighbor)] {
        if list.is_empty() {
            return Err(ChemError::ValidationError(format!("no {el} atom in the first frame")));
        }
    }
    let same = params.center == params.neighbor;

    // 每帧 (Mᵀ, (Mᵀ)⁻¹)：分数坐标取最小镜像；无盒子的帧直接用笛卡尔差
    let boxes: Vec<Option<(Matrix3<f64>, Matrix3<f64>)>> = traj.frames.iter().map(|f| {
        f.cell.as_ref().map(|c| {
            let m = c.matrix.transpose();
            m.try_inverse().map(|inv| (m, inv))
                .ok_or_else(|| ChemError::ValidationError("cell matrix is singular".into()))
        }).transpose()
    }).collect::<ferro_core::Result<_>>()?;
    let dist = |k: usize, i: usize, j: usize| -> f64 {
        let f = &traj.frames[k];
        let d = f.atoms[j].position - f.atoms[i].position;
        match &boxes[k] {
            Some((m, inv)) => {
                let fr = inv * d;
                (m * fr.map(|c| c - (c + 0.5).floor())).norm()
            }
            None => d.norm(),
        }
    };

    // 第一遍：任一帧进入 r_bond 的原子对都是候选（同元素时只取 i < j）
    let mut candidates: Vec<(usize, usize)> = (0..n_frames).into_par_iter().flat_map_iter(|k| {
        let mut hit = Vec::new();
        for &i in &centers {
            for &j in &neighbors {
                if i == j || (same && j < i) { continue; }
                if dist(k, i, j) <= r_bond { hit.push((i, j)); }
            }
        }
        hit
    }).collect();
    candidates.sort_unstable();
    candidates.dedup();
    if candidates.is_empty() {
        return Err(ChemError::ValidationError(format!(
            "no {}-{} pair within r_bond = {r_bond} Ang in any frame", params.center, params.neighbor)));
    }

    // 第二遍：每个候选的存在序列 h（滞后规则）与补缺口后的 h_fill
    let series: Vec<(Vec<bool>, Vec<bool>)> = candidates.par_iter().map(|&(i, j)| {
        let d: Vec<f64> = (0..n_frames).map(|k| dist(k, i, j)).collect();
        let h = existence(&d, r_bond, r_break);
        let mut hf = h.clone();
        fill_gaps(&mut hf, params.intermittency);
        (h, hf)
    }).collect();

    // 间歇型：FFT 自相关；连续型：成段长度 L 对 lag m 贡献 max(0, L − m)
    let plan = AutocorrPlan::new(n_frames);
    let zero = || [vec![0.0; max_lag + 1], vec![0.0; max_lag + 1], vec![0.0; max_lag + 1], vec![0.0; max_lag + 1]];
    let [num_i, den_i, num_c, den_c] = series.par_iter()
        .map_init(|| (plan.worker(), vec![0.0; n_frames], vec![0.0; max_lag + 1]), |(ac, x, s), (h, hf)| {
            let mut out = zero();
            for (xt, &b) in x.iter_mut().zip(h) { *xt = if b { 1.0 } else { 0.0 }; }
            ac.sums(x, s);
            out[0].copy_from_slice(s);
            // 分母：lag m 可用的原点 t ≤ N−1−m 中成键的个数
            for (series, den) in [(h, 1usize), (hf, 3)] {
                let mut acc: f64 = series.iter().filter(|&&b| b).count() as f64;
                for m in 0..=max_lag {
                    if m > 0 && series[n_frames - m] { acc -= 1.0; }
                    out[den][m] = acc;
                }
            }
            let mut t = 0;
            while t < n_frames {
                if !hf[t] { t += 1; continue; }
                let a = t;
                while t < n_frames && hf[t] { t += 1; }
                let len = t - a;
                for (m, v) in out[2].iter_mut().take(len.min(max_lag + 1)).enumerate() {
                    *v += (len - m) as f64;
                }
            }
            out
        })
        .reduce(zero, |mut a, b| {
            for (aa, bb) in a.iter_mut().zip(&b) {
                for (x, y) in aa.iter_mut().zip(bb) { *x += y; }
            }
            a
        });
    // FFT 的和带舍入：0/1 序列的相关应为整数
    let ratio = |n: &[f64], d: &[f64], round: bool| -> Vec<f64> {
        n.iter().zip(d).map(|(&a, &b)| {
            let a = if round { a.round() } else { a };
            if b < 0.5 { f64::NAN } else { a / b }
        }).collect()
    };
    let c_int = ratio(&num_i, &den_i, true);
    let s_cont = ratio(&num_c, &den_c, false);

    // 事件：补缺口后的 h 上逐帧计数
    let mut n_bonds = vec![0.0; n_frames];
    let mut formed = vec![0.0; n_frames];
    let mut broken = vec![0.0; n_frames];
    for (_, hf) in &series {
        for k in 0..n_frames {
            if hf[k] { n_bonds[k] += 1.0; }
            if k > 0 && hf[k] && !hf[k - 1] { formed[k] += 1.0; }
            if k > 0 && !hf[k] && hf[k - 1] { broken[k] += 1.0; }
        }
    }
    // 第 0 帧没有「上一帧」：形成/断裂无定义，给 NaN 而不是 0
    formed[0] = f64::NAN;
    broken[0] = f64::NAN;

    Ok(BondLifeResult {
        time: (0..=max_lag).map(|m| m as f64 * params.dt).collect(),
        c_int, s_cont,
        frame_time: (0..n_frames).map(|k| k as f64 * params.dt).collect(),
        n_bonds, formed, broken,
        n_frames,
        n_centers: centers.len(),
        n_candidates: candidates.len(),
        min_origins: n_frames - max_lag,
        params: params.clone(),
    })
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

impl BondLifeResult {
    /// Two tables: `bondlife` (`time, c_int, s_cont`) and `events`
    /// (`time, n_bonds, formed, broken`, one row per frame).
    pub fn to_tables(&self) -> Vec<(String, Table)> {
        let mut c = Table::new();
        c.push_num("time", self.time.clone())
            .push_num("c_int", self.c_int.clone())
            .push_num("s_cont", self.s_cont.clone());
        let mut e = Table::new();
        e.push_num("time", self.frame_time.clone())
            .push_num("n_bonds", self.n_bonds.clone())
            .push_num("formed", self.formed.clone())
            .push_num("broken", self.broken.clone());
        vec![("bondlife".to_string(), c), ("events".to_string(), e)]
    }

    /// `(tau_int_integral, tau_int_1e, tau_cont_integral, tau_cont_1e)` \[fs\] —
    /// see [`lifetimes`].
    pub fn taus(&self) -> (f64, f64, f64, f64) {
        let (ii, i1) = lifetimes(&self.c_int, self.params.dt);
        let (ci, c1) = lifetimes(&self.s_cont, self.params.dt);
        (ii, i1, ci, c1)
    }

    /// Mean number of bonds per frame.
    pub fn mean_bonds(&self) -> f64 {
        self.n_bonds.iter().sum::<f64>() / self.n_frames as f64
    }

    /// Parameter block for the comment header above the data (batch-shared values only).
    pub fn meta_lines(&self) -> Vec<String> {
        let p = &self.params;
        vec![
            format!("center       = {}", p.center),
            format!("neighbor     = {}", p.neighbor),
            format!("r_bond       = {} Ang  (a free pair bonds at r <= r_bond)", p.r_bond),
            format!("r_break      = {} Ang  (a bonded pair stays bonded while r <= r_break)",
                p.r_break.unwrap_or(p.r_bond)),
            format!("intermittency = {} frames  (breaks this short are filled; s_cont and events only)",
                p.intermittency),
            match p.max_lag {
                Some(m) => format!("max lag      = {m} frames"),
                None => "max lag      = half of each input's frames (see [inputs])".to_string(),
            },
            format!("dt           = {} fs", p.dt),
            "c_int        = sum h(t)h(t+m) / sum h(t), all origins (Luzar-Chandler, intermittent)".to_string(),
            "s_cont       = bonded at every frame of [t, t+m] / bonded at t (Rapaport, continuous)".to_string(),
            "tau_*_int    = trapezoidal integral to max lag; tau_*_1e = first crossing of 1/e".to_string(),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame};
    use nalgebra::Vector3;

    /// n_o 组 Si–O，第 k 组的键长在第 t 帧为 `r(t, k)`（沿 x）。各组沿 y 相隔 6 Å，
    /// 跨组距离恒 > 6 Å，不会成为候选
    fn traj_from(n: usize, n_o: usize, r: impl Fn(usize, usize) -> f64) -> Trajectory {
        let cell = Cell::from_lengths_angles(40.0, 40.0, 40.0, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for t in 0..n {
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            for k in 0..n_o {
                let y = 5.0 + 6.0 * k as f64;
                f.add_atom(Atom::new("Si", Vector3::new(20.0, y, 20.0)));
                f.add_atom(Atom::new("O", Vector3::new(20.0 + r(t, k), y, 20.0)));
            }
            traj.add_frame(f);
        }
        traj
    }

    fn params() -> BondLifeParams {
        BondLifeParams { r_bond: 2.0, ..Default::default() }
    }

    /// 暴力参考：逐原点逐 lag 直接数
    fn brute(h: &[Vec<bool>], max_lag: usize, continuous: bool) -> Vec<f64> {
        let n = h[0].len();
        (0..=max_lag).map(|m| {
            let (mut num, mut den) = (0.0, 0.0);
            for s in h {
                for t in 0..n - m {
                    if !s[t] { continue; }
                    den += 1.0;
                    let ok = if continuous { (t..=t + m).all(|u| s[u]) } else { s[t + m] };
                    if ok { num += 1.0; }
                }
            }
            if den == 0.0 { f64::NAN } else { num / den }
        }).collect()
    }

    /// 伪随机但可复现
    fn jitter(k: usize) -> f64 {
        ((k as f64 * 12.9898).sin() * 43758.5453).fract() - 0.5
    }

    #[test]
    fn test_both_functions_match_brute_force() {
        // 3 个 O 在 2.0 Å 阈值附近乱跳：间歇、连续两种都与逐原点直接计数一致
        let (n, no) = (40, 3);
        let traj = traj_from(n, no, |t, k| 2.0 + 0.6 * jitter(7 * t + 13 * k));
        let r = calc_bondlife(&traj, &BondLifeParams { max_lag: Some(39), ..params() }).unwrap();
        let h: Vec<Vec<bool>> = (0..no)
            .map(|k| (0..n).map(|t| 2.0 + 0.6 * jitter(7 * t + 13 * k) <= 2.0).collect()).collect();
        for (m, (&ci, &sc)) in r.c_int.iter().zip(&r.s_cont).enumerate() {
            let (wi, wc) = (brute(&h, 39, false)[m], brute(&h, 39, true)[m]);
            assert!((ci - wi).abs() < 1e-12 || (ci.is_nan() && wi.is_nan()), "C_I lag {m}: {ci} vs {wi}");
            assert!((sc - wc).abs() < 1e-12 || (sc.is_nan() && wc.is_nan()), "S_C lag {m}: {sc} vs {wc}");
        }
    }

    #[test]
    fn test_permanent_bond_is_one() {
        let traj = traj_from(10, 1, |_, _| 1.6);
        let r = calc_bondlife(&traj, &params()).unwrap();
        assert!(r.c_int.iter().chain(&r.s_cont).all(|&v| (v - 1.0).abs() < 1e-12));
        // 永不断的键在 lag 轴内降不到 1/e：交点留空，不外推
        assert!(lifetimes(&r.c_int, 1.0).1.is_nan());
    }

    #[test]
    fn test_intermittent_counts_reformed_continuous_does_not() {
        // 键：成、成、断、成、成 —— lag 3 从 t=0 到 t=3 两端都成键：间歇算存活，连续不算
        let d = [1.8, 1.8, 2.5, 1.8, 1.8, 1.8, 1.8, 1.8];
        let traj = traj_from(8, 1, |t, _| d[t]);
        let r = calc_bondlife(&traj, &BondLifeParams { max_lag: Some(3), ..params() }).unwrap();
        // 间歇 lag 3：原点 t∈{0,1,3,4}（t≤4 且成键），两端成键的是 0→3、1→4、3→6、4→7 → 4/4
        assert!((r.c_int[3] - 1.0).abs() < 1e-12, "C_I(3) = {}", r.c_int[3]);
        // 连续 lag 3：只有 t=3、4 的 [t, t+3] 全程成键 → 2/4
        assert!((r.s_cont[3] - 0.5).abs() < 1e-12, "S_C(3) = {}", r.s_cont[3]);
    }

    #[test]
    fn test_hysteresis_suppresses_flicker() {
        // 在 1.9 ↔ 2.1 之间抖动：单阈值 2.0 会反复断连；r_break = 2.3 时一直成键
        let traj = traj_from(10, 1, |t, _| if t % 2 == 0 { 1.9 } else { 2.1 });
        let single = calc_bondlife(&traj, &params()).unwrap();
        let hyst = calc_bondlife(&traj, &BondLifeParams { r_break: Some(2.3), ..params() }).unwrap();
        assert!(single.formed.iter().skip(1).sum::<f64>() > 0.0);
        assert_eq!(hyst.formed.iter().skip(1).sum::<f64>(), 0.0);
        assert!(hyst.s_cont.iter().all(|&v| (v - 1.0).abs() < 1e-12));
    }

    #[test]
    fn test_hysteresis_needs_r_bond_to_form() {
        // 起始在 2.1（介于 r_bond 与 r_break 之间）：未成键的对要到 r_bond 以内才成键
        let d = [2.1, 2.1, 1.9, 2.1, 2.4];
        let h = existence(&d, 2.0, 2.3);
        assert_eq!(h, vec![false, false, true, true, false]);
    }

    #[test]
    fn test_fill_gaps_boundaries() {
        let mut h = vec![true, false, false, true, false, false, false, true, false];
        fill_gaps(&mut h, 2);
        // 2 帧缺口补上；3 帧不补；末尾缺口通向轨迹尽头，不补
        assert_eq!(h, vec![true, true, true, true, false, false, false, true, false]);
        let mut lead = vec![false, true];
        fill_gaps(&mut lead, 5);
        assert_eq!(lead, vec![false, true], "开头的缺口不知道从何时开始，不补");
    }

    #[test]
    fn test_intermittency_changes_only_continuous() {
        let d = [1.8, 1.8, 2.5, 1.8, 1.8, 1.8, 1.8, 1.8];
        let traj = traj_from(8, 1, |t, _| d[t]);
        let r = calc_bondlife(&traj, &BondLifeParams { max_lag: Some(3), intermittency: 1, ..params() }).unwrap();
        assert!((r.s_cont[3] - 1.0).abs() < 1e-12, "补上 1 帧缺口后连续型应全程存活");
        assert!((r.c_int[3] - 1.0).abs() < 1e-12);
        assert!(r.formed[1..].iter().all(|&v| v == 0.0), "补缺口后没有断后重连事件");
    }

    #[test]
    fn test_events_conserve_bond_count() {
        let traj = traj_from(30, 4, |t, k| 2.0 + 0.8 * jitter(3 * t + 5 * k + 1));
        let r = calc_bondlife(&traj, &params()).unwrap();
        assert!(r.formed[0].is_nan() && r.broken[0].is_nan());
        for k in 1..30 {
            assert_eq!(r.n_bonds[k], r.n_bonds[k - 1] + r.formed[k] - r.broken[k], "frame {k}");
        }
    }

    #[test]
    fn test_same_element_counts_each_pair_once() {
        let cell = Cell::from_lengths_angles(20.0, 20.0, 20.0, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for _ in 0..4 {
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            f.add_atom(Atom::new("O", Vector3::new(5.0, 5.0, 5.0)));
            f.add_atom(Atom::new("O", Vector3::new(7.0, 5.0, 5.0)));
            traj.add_frame(f);
        }
        let r = calc_bondlife(&traj, &BondLifeParams {
            center: "O".into(), neighbor: "O".into(), r_bond: 3.0, ..Default::default()
        }).unwrap();
        assert_eq!(r.n_candidates, 1);
        assert_eq!(r.n_bonds[0], 1.0);
    }

    #[test]
    fn test_minimum_image_across_the_box() {
        // Si 在 0.5，O 在 39.0：盒长 40 下相距 1.5 Å
        let cell = Cell::from_lengths_angles(40.0, 40.0, 40.0, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for _ in 0..3 {
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            f.add_atom(Atom::new("Si", Vector3::new(0.5, 1.0, 1.0)));
            f.add_atom(Atom::new("O", Vector3::new(39.0, 1.0, 1.0)));
            traj.add_frame(f);
        }
        assert_eq!(calc_bondlife(&traj, &params()).unwrap().n_candidates, 1);
    }

    #[test]
    fn test_lifetimes() {
        // C = e^{-m/4}，dt = 1：1/e 交点在 m = 4（线性插值有小偏差）
        let c: Vec<f64> = (0..20).map(|m| (-(m as f64) / 4.0).exp()).collect();
        let (_, t1e) = lifetimes(&c, 1.0);
        assert!((t1e - 4.0).abs() < 0.1, "{t1e}");
    }

    #[test]
    fn test_bad_thresholds_and_no_pair() {
        let traj = traj_from(4, 1, |_, _| 3.0);
        assert!(calc_bondlife(&traj, &params()).unwrap_err().to_string().contains("no Si-O pair"));
        assert!(calc_bondlife(&traj, &BondLifeParams { r_break: Some(1.0), ..params() }).is_err());
    }

    #[test]
    fn test_tables() {
        let traj = traj_from(6, 1, |_, _| 1.6);
        let r = calc_bondlife(&traj, &params()).unwrap();
        let t = r.to_tables();
        assert_eq!(t[0].1.names(), vec!["time", "c_int", "s_cont"]);
        assert_eq!(t[1].0, "events");
        assert_eq!(t[1].1.names(), vec!["time", "n_bonds", "formed", "broken"]);
        assert_eq!(t[1].1.n_rows(), 6);
    }
}
