//! Rotational autocorrelation function C₂(t).
//!
//!   C₂(m) = ⟨P₂(û(t)·û(t+m))⟩,  P₂(x) = (3x² − 1)/2
//!
//! Two orientation vectors ([`RotVector`]):
//! - `Sum`: u_c(t) = Σ_{n ∈ neighbors(c, r_cut, t)} minimum_image(r_n − r_c) — the bonds
//!   from each center to every neighbor-element atom within `r_cut`, found afresh in
//!   every frame. A center with no neighbor in a frame has no orientation there (an
//!   *invalid* frame).
//! - `Bond`: every center–neighbor pair within `r_cut` **in the first frame** is one
//!   unit, followed by atom identity for the whole run whatever its length later —
//!   `gmx rotacf -d`. Every frame is valid.
//!
//! Order ([`Legendre`]): `C₁ = ⟨û(t)·û(t+m)⟩` or `C₂ = ⟨P₂(û(t)·û(t+m))⟩`.
//!
//! Every lag averages over all valid (molecule, origin) pairs — pairs whose vectors
//! are valid at both ends. With `q_ij = û_i û_j` (zero on invalid frames) and the
//! validity indicator χ, `(û(t)·û(t+m))² = Σ_ij q_ij(t) q_ij(t+m)`, so
//!
//!   C₂(m) = [ (3/2) Σ_ij w_ij AC[q_ij](m) − (1/2) AC[χ](m) ] / AC[χ](m)
//!   C₁(m) = Σ_i AC[û_i](m) / AC[χ](m)
//!
//! with `w = 1` on the diagonal, 2 off it, and AC the all-origin autocorrelation
//! sum (FFT, [`super::correlate`]). With every frame valid this is exactly the
//! `gmx rotacf -P 2` sum (`autocorr.cpp`: 1.5 on the diagonal, 3 off it, −0.5·(N−m)).
//!
//! Parallelism: per molecule.

use rayon::prelude::*;
use ferro_core::{Table, Trajectory};
use ferro_core::error::ChemError;

use super::correlate::{cumulative_trapezoid, resolve_max_lag, AutocorrPlan};

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Which vector stands for a unit's orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RotVector {
    /// Sum of the center's bonds within `r_cut`, neighbors found per frame
    #[default]
    Sum,
    /// One center–neighbor bond per unit, pairs fixed in the first frame (`gmx rotacf -d`)
    Bond,
}

/// Order of the Legendre polynomial in `C_ℓ(t) = ⟨P_ℓ(û(0)·û(t))⟩`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Legendre {
    /// `P₁(x) = x`
    P1,
    /// `P₂(x) = (3x² − 1)/2`
    #[default]
    P2,
}

/// Parameters for rotational autocorrelation function calculation.
#[derive(Debug, Clone)]
pub struct RotCorrParams {
    /// Element symbol of the center atom (e.g. `"O"` for water)
    pub center: String,
    /// Element symbol of the neighbor atoms (e.g. `"H"` for water)
    pub neighbor: String,
    /// Cutoff radius for neighbor search \[Å\] (default: 1.2)
    pub r_cut: f64,
    /// Longest lag in frames, `1 ..= n_frames − 1` (`None` = `n_frames / 2`)
    pub max_lag: Option<usize>,
    /// Time step per frame \[fs\] (default: 1.0)
    pub dt: f64,
    /// Orientation vector (default: `Sum`)
    pub vector: RotVector,
    /// Legendre order (default: `P2`)
    pub legendre: Legendre,
}

impl Default for RotCorrParams {
    fn default() -> Self {
        RotCorrParams {
            center: "O".to_string(),
            neighbor: "H".to_string(),
            r_cut: 1.2,
            max_lag: None,
            dt: 1.0,
            vector: RotVector::Sum,
            legendre: Legendre::P2,
        }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Result of a rotational autocorrelation function calculation.
///
/// `rotcorr[m] = C_ℓ(m·dt)` — in `[−0.5, 1]` for P₂, `[−1, 1]` for P₁ — and `NaN`
/// at a lag with no valid pair.
#[derive(Debug, Clone)]
pub struct RotCorrResult {
    /// Lag-time axis \[fs\]; `time[m] = m · dt`, `m = 0 ..= max_lag`
    pub time: Vec<f64>,
    /// Rotational correlation function C_ℓ(t)
    pub rotcorr: Vec<f64>,
    /// Running trapezoidal integral ∫₀ᵗ C₂ dτ \[fs\]
    pub integral: Vec<f64>,
    /// Frames in the trajectory
    pub n_frames: usize,
    /// Center atoms of the center element
    pub n_centers: usize,
    /// Units correlated: center atoms (`Sum`) or bonds (`Bond`)
    pub n_units: usize,
    /// Time origins at the longest lag, `n_frames − max_lag`
    pub min_origins: usize,
    /// Fraction of (unit, frame) with a valid orientation vector (1 for `Bond`)
    pub valid_fraction: f64,
    pub params: RotCorrParams,
}

// ─── 计算 ────────────────────────────────────────────────────────────────────

/// Computes C₂(t), averaging each lag over all valid (molecule, origin) pairs.
///
/// Returns `Err` if the trajectory has fewer than 2 frames, `max_lag` is outside
/// `1 ..= n_frames − 1`, no atom is of the center element, or no center has a
/// neighbor within `r_cut` in any frame.
pub fn calc_rotcorr(traj: &Trajectory, params: &RotCorrParams) -> ferro_core::Result<RotCorrResult> {
    let n_frames = traj.n_frames();
    let max_lag = resolve_max_lag(n_frames, params.max_lag)?;

    // 确定参与计算的 center 原子下标（按第一帧筛选）
    let ref_frame = &traj.frames[0];
    let center_indices: Vec<usize> = ref_frame.atoms.iter().enumerate()
        .filter(|(_, a)| a.element == params.center)
        .map(|(i, _)| i)
        .collect();
    if center_indices.is_empty() {
        return Err(ChemError::ValidationError(format!(
            "no {} atom in the first frame", params.center)));
    }

    let has_cell = ref_frame.cell.is_some();
    let r_cut2 = params.r_cut * params.r_cut;

    // 当帧的键矢量：周期体系取最小镜像
    let bond = |frame: &ferro_core::Frame, from: usize, to: usize| {
        let d = frame.atoms[to].position - frame.atoms[from].position;
        match (&frame.cell, has_cell) {
            (Some(cell), true) => cell.minimum_image(d).expect("cell is non-singular"),
            _ => d,
        }
    };
    let is_neighbor = |i: usize, ci: usize| -> bool {
        i != ci && ref_frame.atoms[i].element == params.neighbor
    };

    // orient[帧][单元] = 取向向量（未归一化；零向量 = 该帧无效）
    let orient: Vec<Vec<[f64; 3]>> = match params.vector {
        RotVector::Sum => traj.frames.iter().map(|frame| {
            center_indices.iter().map(|&ci| {
                let mut u = [0.0_f64; 3];
                let mut len_sum = 0.0;
                for ni in (0..frame.atoms.len()).filter(|&ni| is_neighbor(ni, ci)) {
                    let d = bond(frame, ci, ni);
                    if d.norm_squared() < r_cut2 {
                        u[0] += d.x;
                        u[1] += d.y;
                        u[2] += d.z;
                        len_sum += d.norm();
                    }
                }
                // 对称单元的键相互抵消时，合向量只剩 ~1e-16 的舍入残差，方向是噪声。
                // 相对键长之和判抵消（绝对阈值 1e-30 挡不住它），按无效帧处理
                let n = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
                if n < 1e-6 * len_sum { [0.0; 3] } else { u }
            }).collect()
        }).collect(),
        RotVector::Bond => {
            // 第一帧定键，之后按原子身份一直跟踪（gmx rotacf -d）：不再看距离
            let pairs: Vec<(usize, usize)> = center_indices.iter()
                .flat_map(|&ci| (0..ref_frame.atoms.len())
                    .filter(move |&ni| is_neighbor(ni, ci))
                    .map(move |ni| (ci, ni)))
                .filter(|&(ci, ni)| bond(ref_frame, ci, ni).norm_squared() < r_cut2)
                .collect();
            if pairs.is_empty() {
                return Err(ChemError::ValidationError(format!(
                    "no {}-{} bond within r_cut = {} Ang in the first frame",
                    params.center, params.neighbor, params.r_cut)));
            }
            traj.frames.par_iter().map(|frame| {
                pairs.iter().map(|&(ci, ni)| {
                    let d = bond(frame, ci, ni);
                    [d.x, d.y, d.z]
                }).collect()
            }).collect()
        }
    };

    let n_units = orient[0].len();

    // 逐分子：单位向量的 6 个分量 q_ij（无效帧为 0）与有效性指示 χ，各做全原点自相关
    let plan = AutocorrPlan::new(n_frames);
    let zero = || (vec![0.0; max_lag + 1], vec![0.0; max_lag + 1]);
    // 分量 (i, j, 权重)：P₂ 取 û_iû_j（非对角 ×2），P₁ 取 û_i（j = None）
    let comps: &[(usize, Option<usize>, f64)] = match params.legendre {
        Legendre::P2 => &[(0, Some(0), 1.0), (1, Some(1), 1.0), (2, Some(2), 1.0),
                          (0, Some(1), 2.0), (1, Some(2), 2.0), (2, Some(0), 2.0)],
        Legendre::P1 => &[(0, None, 1.0), (1, None, 1.0), (2, None, 1.0)],
    };
    let (sq, cnt) = (0..n_units).into_par_iter()
        .map_init(|| (plan.worker(), vec![0.0; n_frames], vec![0.0; max_lag + 1]), |(ac, x, s), k| {
            let unit: Vec<Option<[f64; 3]>> = orient.iter().map(|step| {
                let [a, b, c] = step[k];
                let n2 = a * a + b * b + c * c;
                // 零向量：该帧没有邻居，无取向可言
                (n2 >= 1e-30).then(|| { let n = n2.sqrt(); [a / n, b / n, c / n] })
            }).collect();
            let (mut sq, mut cnt) = zero();
            for &(i, j, w) in comps {
                for (xt, u) in x.iter_mut().zip(&unit) {
                    *xt = u.map_or(0.0, |u| u[i] * j.map_or(1.0, |j| u[j]));
                }
                ac.sums(x, s);
                for (o, v) in sq.iter_mut().zip(s.iter()) {
                    *o += w * v;
                }
            }
            for (xt, u) in x.iter_mut().zip(&unit) {
                *xt = if u.is_some() { 1.0 } else { 0.0 };
            }
            ac.sums(x, s);
            cnt.copy_from_slice(s);
            (sq, cnt)
        })
        .reduce(zero, |(mut a, mut b), (c, d)| {
            for (x, y) in a.iter_mut().zip(&c) { *x += y; }
            for (x, y) in b.iter_mut().zip(&d) { *x += y; }
            (a, b)
        });

    let n_valid = orient.iter()
        .map(|step| step.iter().filter(|u| u.iter().map(|c| c * c).sum::<f64>() >= 1e-30).count())
        .sum::<usize>();
    if n_valid == 0 {
        return Err(ChemError::ValidationError(format!(
            "no {} atom has an orientation in any frame: either no {} neighbor within \
             r_cut = {} Ang, or its bonds cancel (a symmetric unit such as a tetrahedron; \
             use --vector bond)",
            params.center, params.neighbor, params.r_cut)));
    }

    // 有效配对数是 0/1 序列的自相关，FFT 下带舍入：四舍五入回整数；为 0 的 lag 给 NaN
    let rotcorr: Vec<f64> = sq.iter().zip(&cnt).map(|(&q, &k)| {
        let k = k.round();
        match (k < 1.0, params.legendre) {
            (true, _) => f64::NAN,
            (false, Legendre::P2) => (1.5 * q - 0.5 * k) / k,
            (false, Legendre::P1) => q / k,
        }
    }).collect();
    let time: Vec<f64> = (0..=max_lag).map(|m| m as f64 * params.dt).collect();
    let integral = cumulative_trapezoid(&rotcorr, params.dt);

    Ok(RotCorrResult {
        time, rotcorr, integral,
        n_frames,
        n_centers: center_indices.len(),
        n_units,
        min_origins: n_frames - max_lag,
        valid_fraction: n_valid as f64 / (n_units * n_frames) as f64,
        params: params.clone(),
    })
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

impl RotCorrResult {
    /// Projects the result into the table the writers consume.
    ///
    /// Rotational correlation: `time, c2, integral` (`c1` in place of `c2` for P₁).
    /// The `file` column is added by the caller when stacking several inputs
    /// (see `ferro_core::Table::concat_union`).
    pub fn to_tables(&self) -> Vec<(String, Table)> {
        let mut t = Table::new();
        t.push_num("time", self.time.clone())
            .push_num(match self.params.legendre { Legendre::P1 => "c1", Legendre::P2 => "c2" },
                self.rotcorr.clone())
            .push_num("integral", self.integral.clone());
        vec![("rotcorr".to_string(), t)]
    }

    /// Parameter block for the comment header above the data.
    ///
    /// Only what the whole batch shares; molecule counts, origins and the valid
    /// fraction differ per input and go to the `[inputs]` list.
    pub fn meta_lines(&self) -> Vec<String> {
        let p = &self.params;
        vec![
            format!("center   = {}", p.center),
            format!("neighbor = {}", p.neighbor),
            format!("r_cut    = {} Ang", p.r_cut),
            match p.vector {
                RotVector::Sum => "vector   = sum: sum of center->neighbor bonds within r_cut, found per frame",
                RotVector::Bond => "vector   = bond: each center-neighbor pair within r_cut in frame 0, followed by atom identity (gmx rotacf -d)",
            }.to_string(),
            match p.legendre {
                Legendre::P1 => "order    = P1: c = <u(t).u(t+m)>",
                Legendre::P2 => "order    = P2: c = <(3 (u(t).u(t+m))^2 - 1) / 2>",
            }.to_string(),
            match p.max_lag {
                Some(m) => format!("max lag  = {m} frames"),
                None => "max lag  = half of each input's frames (see [inputs])".to_string(),
            },
            "origins  = all: each lag averages every (molecule, origin) pair valid at both ends (FFT)".to_string(),
            format!("dt       = {} fs", p.dt),
            "integral = running trapezoidal integral of c2 [fs]".to_string(),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    /// 构建单水分子轨迹，所有帧取向相同（O 在原点，两个 H 固定）
    fn make_fixed_orientation(n: usize) -> Trajectory {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for _ in 0..n {
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));
            frame.add_atom(Atom::new("H", Vector3::new(0.96, 0.0,  0.0)));
            frame.add_atom(Atom::new("H", Vector3::new(-0.24, 0.93, 0.0)));
            traj.add_frame(frame);
        }
        traj
    }

    /// 构建两帧轨迹：t=0 取向沿 x，t=1 取向沿 y（cosθ=0 → P₂=-0.5）
    fn make_perpendicular_orientations() -> Trajectory {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();

        // 帧0：H 沿 +x 方向
        let mut f0 = Frame::with_cell(cell.clone(), [true; 3]);
        f0.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));
        f0.add_atom(Atom::new("H", Vector3::new(1.0, 0.0, 0.0)));
        f0.add_atom(Atom::new("H", Vector3::new(1.0, 0.0, 0.0)));
        traj.add_frame(f0);

        // 帧1：H 沿 +y 方向（与帧0 垂直）
        let mut f1 = Frame::with_cell(cell.clone(), [true; 3]);
        f1.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));
        f1.add_atom(Atom::new("H", Vector3::new(0.0, 1.0, 0.0)));
        f1.add_atom(Atom::new("H", Vector3::new(0.0, 1.0, 0.0)));
        traj.add_frame(f1);

        traj
    }

    /// 伪随机但可复现
    fn jitter(k: usize) -> f64 {
        ((k as f64 * 12.9898).sin() * 43758.5453).fract() - 0.5
    }

    /// 3 个 O 各带一个 H，H 方向逐帧乱转；`invalid(t, k)` 为真的帧把 H 挪出截断
    fn tumbling(n: usize, invalid: impl Fn(usize, usize) -> bool) -> Trajectory {
        let cell = Cell::from_lengths_angles(30.0, 30.0, 30.0, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for t in 0..n {
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            for k in 0..3 {
                let o = Vector3::new(5.0 + 8.0 * k as f64, 10.0, 10.0);
                let d = Vector3::new(jitter(7 * t + k), jitter(11 * t + 3 * k + 1), jitter(13 * t + k + 2))
                    .normalize();
                let h = if invalid(t, k) { o + d * 4.0 } else { o + d * 0.96 };
                f.add_atom(Atom::new("O", o));
                f.add_atom(Atom::new("H", h));
            }
            traj.add_frame(f);
        }
        traj
    }

    /// 暴力参考：逐分子取单位向量（无邻居 → None），对每个 lag 在两端都有效的全部
    /// (分子, 原点) 配对上直接平均 P₂
    fn brute_c2(traj: &Trajectory, r_cut: f64, max_lag: usize) -> Vec<f64> {
        let n = traj.n_frames();
        let unit = |t: usize, k: usize| -> Option<Vector3<f64>> {
            let f = &traj.frames[t];
            let o = f.atoms[2 * k].position;
            let mut u = Vector3::zeros();
            for a in f.atoms.iter().filter(|a| a.element == "H") {
                let d = f.cell.as_ref().unwrap().minimum_image(a.position - o).unwrap();
                if d.norm() < r_cut { u += d; }
            }
            (u.norm_squared() >= 1e-30).then(|| u.normalize())
        };
        (0..=max_lag).map(|m| {
            let (mut sum, mut cnt) = (0.0, 0usize);
            for t in 0..n - m {
                for k in 0..3 {
                    if let (Some(a), Some(b)) = (unit(t, k), unit(t + m, k)) {
                        let c = a.dot(&b);
                        sum += 0.5 * (3.0 * c * c - 1.0);
                        cnt += 1;
                    }
                }
            }
            if cnt == 0 { f64::NAN } else { sum / cnt as f64 }
        }).collect()
    }

    #[test]
    fn test_fft_matches_brute_force_all_valid() {
        // 全部帧有效：即 gmx rotacf -P 2 的情形
        let traj = tumbling(30, |_, _| false);
        let r = calc_rotcorr(&traj, &RotCorrParams { max_lag: Some(29), ..Default::default() }).unwrap();
        for (m, want) in brute_c2(&traj, 1.2, 29).iter().enumerate() {
            assert!((r.rotcorr[m] - want).abs() < 1e-10, "lag {m}: {} vs {want}", r.rotcorr[m]);
        }
        assert!((r.valid_fraction - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_fft_matches_brute_force_with_invalid_frames() {
        // 部分帧没有邻居：只在两端都有效的配对上平均，分母是有效配对数而不是原点数
        let traj = tumbling(30, |t, k| (t + 2 * k) % 5 == 0);
        let r = calc_rotcorr(&traj, &RotCorrParams { max_lag: Some(29), ..Default::default() }).unwrap();
        for (m, want) in brute_c2(&traj, 1.2, 29).iter().enumerate() {
            if want.is_nan() {
                assert!(r.rotcorr[m].is_nan(), "lag {m}: 无有效配对应为 NaN");
            } else {
                assert!((r.rotcorr[m] - want).abs() < 1e-10, "lag {m}: {} vs {want}", r.rotcorr[m]);
            }
        }
        assert!(r.valid_fraction < 1.0);
    }

    #[test]
    fn test_lag_without_valid_pair_is_nan() {
        // 分子只在第 0、1 帧有效：lag ≥ 2 没有两端都有效的配对 → NaN，而不是 0
        let traj = tumbling(6, |t, _| t >= 2);
        let r = calc_rotcorr(&traj, &RotCorrParams { max_lag: Some(5), ..Default::default() }).unwrap();
        assert!((r.rotcorr[0] - 1.0).abs() < 1e-10);
        assert!(r.rotcorr[2..].iter().all(|c| c.is_nan()), "{:?}", r.rotcorr);
    }

    #[test]
    fn test_c0_always_one() {
        let traj = make_fixed_orientation(8);
        let res = calc_rotcorr(&traj, &RotCorrParams::default()).unwrap();
        assert!((res.rotcorr[0] - 1.0).abs() < 1e-10, "C(0) = {}", res.rotcorr[0]);
    }

    #[test]
    fn test_fixed_orientation_flat() {
        let traj = make_fixed_orientation(10);
        let res = calc_rotcorr(&traj, &RotCorrParams::default()).unwrap();
        for (i, &c) in res.rotcorr.iter().enumerate() {
            assert!((c - 1.0).abs() < 1e-10, "C({i}) = {c}");
        }
    }

    #[test]
    fn test_perpendicular_gives_minus_half() {
        let traj = make_perpendicular_orientations();
        let res = calc_rotcorr(&traj, &RotCorrParams { r_cut: 1.5, max_lag: Some(1), ..Default::default() }).unwrap();
        assert!((res.rotcorr[0] - 1.0).abs() < 1e-10, "C(0)={}", res.rotcorr[0]);
        assert!((res.rotcorr[1] + 0.5).abs() < 1e-10, "C(1) = {}", res.rotcorr[1]);
    }

    #[test]
    fn test_no_neighbor_anywhere_is_an_error() {
        let traj = tumbling(4, |_, _| true);
        let err = calc_rotcorr(&traj, &RotCorrParams::default()).unwrap_err().to_string();
        assert!(err.contains("r_cut"), "应说明是截断内找不到邻居：{err}");
    }

    #[test]
    fn test_integral_is_trapezoidal() {
        // C 恒为 1 → 梯形积分 = t；左矩形法会在 m=0 就给出 dt
        let traj = make_fixed_orientation(8);
        let res = calc_rotcorr(&traj, &RotCorrParams { dt: 2.0, ..Default::default() }).unwrap();
        for (m, &ig) in res.integral.iter().enumerate() {
            assert!((ig - 2.0 * m as f64).abs() < 1e-10, "integral[{m}] = {ig}");
        }
    }

    #[test]
    fn test_default_max_lag_and_origins() {
        let res = calc_rotcorr(&make_fixed_orientation(11), &RotCorrParams::default()).unwrap();
        assert_eq!(res.time.len(), 6);
        assert_eq!(res.min_origins, 6);
        assert_eq!(res.n_frames, 11);
    }

    /// 正四面体 PO₄（P 在原点），每帧绕 z 轴转 `w` 弧度
    fn rotating_tetrahedron(n: usize, w: f64) -> Trajectory {
        let cell = Cell::from_lengths_angles(30.0, 30.0, 30.0, 90.0, 90.0, 90.0).unwrap();
        let s = 1.55 / 3f64.sqrt();
        let base = [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]];
        let c = Vector3::new(15.0, 15.0, 15.0);
        let mut traj = Trajectory::new();
        for t in 0..n {
            let (sn, cs) = (w * t as f64).sin_cos();
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            f.add_atom(Atom::new("P", c));
            for b in base {
                let (x, y) = (b[0] * s, b[1] * s);
                f.add_atom(Atom::new("O", c + Vector3::new(cs * x - sn * y, sn * x + cs * y, b[2] * s)));
            }
            traj.add_frame(f);
        }
        traj
    }

    #[test]
    fn test_sum_cancels_on_tetrahedron_bond_does_not() {
        let traj = rotating_tetrahedron(12, 0.3);
        let sum = RotCorrParams { center: "P".into(), neighbor: "O".into(), r_cut: 1.8, ..Default::default() };
        let err = calc_rotcorr(&traj, &sum).unwrap_err().to_string();
        assert!(err.contains("cancel") && err.contains("--vector bond"), "应提示键向量抵消：{err}");

        // 绕 z 转 θ：与 z 轴夹角 β 的键，cos = cos²β + sin²β·cosθ；正四面体四个键 cos²β 都是 1/3
        let bond = RotCorrParams { vector: RotVector::Bond, max_lag: Some(11), ..sum };
        let r = calc_rotcorr(&traj, &bond).unwrap();
        assert_eq!(r.n_units, 4);
        assert_eq!(r.n_centers, 1);
        for (m, &got) in r.rotcorr.iter().enumerate() {
            let x = 1.0 / 3.0 + 2.0 / 3.0 * (0.3 * m as f64).cos();
            let want = 0.5 * (3.0 * x * x - 1.0);
            assert!((got - want).abs() < 1e-10, "lag {m}: {got} vs {want}");
        }
    }

    #[test]
    fn test_p1_bond_matches_brute_force() {
        // 3 个 O 各带一个 H：bond 与 sum 的单元相同；P1 = ⟨û(t)·û(t+m)⟩ 暴力平均
        let traj = tumbling(24, |_, _| false);
        let r = calc_rotcorr(&traj, &RotCorrParams {
            vector: RotVector::Bond, legendre: Legendre::P1, max_lag: Some(23), ..Default::default()
        }).unwrap();
        let unit = |t: usize, k: usize| {
            let f = &traj.frames[t];
            f.cell.as_ref().unwrap()
                .minimum_image(f.atoms[2 * k + 1].position - f.atoms[2 * k].position).unwrap().normalize()
        };
        for m in 0..=23 {
            let mut acc = 0.0;
            for t in 0..24 - m { for k in 0..3 { acc += unit(t, k).dot(&unit(t + m, k)); } }
            let want = acc / (3 * (24 - m)) as f64;
            assert!((r.rotcorr[m] - want).abs() < 1e-10, "lag {m}: {} vs {want}", r.rotcorr[m]);
        }
        assert_eq!(r.to_tables()[0].1.names(), vec!["time", "c1", "integral"]);
    }

    #[test]
    fn test_p1_antiparallel_is_minus_one() {
        let traj = make_perpendicular_orientations();
        let r = calc_rotcorr(&traj, &RotCorrParams {
            r_cut: 1.5, legendre: Legendre::P1, max_lag: Some(1), ..Default::default()
        }).unwrap();
        assert!(r.rotcorr[1].abs() < 1e-10, "垂直时 P1 = 0，实为 {}", r.rotcorr[1]);
    }

    #[test]
    fn test_bond_is_followed_beyond_r_cut() {
        // 第 0 帧键长 0.96 定键；之后 H 退到 4 Å：sum 模式判为无效，bond 模式照 GROMACS 继续跟踪
        let traj = tumbling(6, |t, _| t >= 2);
        let r = calc_rotcorr(&traj, &RotCorrParams {
            vector: RotVector::Bond, max_lag: Some(5), ..Default::default()
        }).unwrap();
        assert!((r.valid_fraction - 1.0).abs() < 1e-12);
        assert!(r.rotcorr.iter().all(|c| c.is_finite()), "bond 模式每个 lag 都有值：{:?}", r.rotcorr);
        let sum = calc_rotcorr(&traj, &RotCorrParams { max_lag: Some(5), ..Default::default() }).unwrap();
        assert!(sum.valid_fraction < 1.0);
    }

    #[test]
    fn test_bond_needs_a_bond_in_the_first_frame() {
        let traj = tumbling(4, |t, _| t == 0);
        let err = calc_rotcorr(&traj, &RotCorrParams { vector: RotVector::Bond, ..Default::default() })
            .unwrap_err().to_string();
        assert!(err.contains("first frame"), "{err}");
    }

    #[test]
    fn test_to_tables_columns() {
        let traj = make_fixed_orientation(8);
        let res = calc_rotcorr(&traj, &RotCorrParams::default()).unwrap();
        let (name, t) = res.to_tables().remove(0);
        assert_eq!(name, "rotcorr");
        assert_eq!(t.names(), vec!["time", "c2", "integral"]);
        assert_eq!(t.n_rows(), res.time.len());
        assert!(t.validate().is_ok());
        let meta = res.meta_lines().join("\n");
        assert!(meta.contains("r_cut"));
        assert!(!meta.contains("molecules"), "分子数是逐文件的量，不进共享区");
    }
}
