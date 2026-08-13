//! Structure factor S(q) calculation and output.
//!
//! Workflow: compute g(r) first (`calc_gr`), then call `calc_sq_from_gr` to obtain S(q).
//! Output: two files — `.gr` (written by the gr module) and `.sq` (written here).
//!
//! Formula (Faber-Ziman, from code2/dump2sq.c CalcSq):
//!   S_ij(q) = 1 + (4πρ/q) Σ_r r[g_ij(r)−1] sin(qr) Δr

use super::gr::{sorted_keys, GrResult};
use ferro_core::error::ChemError;
use ferro_core::Table;
use super::scattering_data::{form_factor_xrd, neutron_bcoh};
use rayon::prelude::*;
use std::collections::BTreeMap;

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Scattering weighting scheme for the total S(q).
///
/// `None`    — equal weights (original Faber-Ziman without form factors)
/// `Xrd`     — X-ray: weighted by q-dependent atomic form factors f(q)
/// `Neutron` — neutron: weighted by q-independent coherent scattering lengths bcoh
/// `Both`    — compute both XRD and neutron weighted totals simultaneously
#[derive(Debug, Clone, PartialEq, Default)]
pub enum SqWeighting {
    #[default]
    None,
    Xrd,
    Neutron,
    Both,
}

/// Parameters for S(q) calculation.
#[derive(Debug, Clone)]
pub struct SqParams {
    /// Minimum q value \[Å⁻¹\] (default: 0.1)
    pub q_min: f64,
    /// Maximum q value \[Å⁻¹\] (default: 25.0)
    pub q_max: f64,
    /// q step size \[Å⁻¹\] (default: 0.02)
    pub dq: f64,
    /// Scattering factor weighting for weighted total S(q) (default: None)
    pub weighting: SqWeighting,
}

impl Default for SqParams {
    fn default() -> Self {
        SqParams { q_min: 0.1, q_max: 25.0, dq: 0.02, weighting: SqWeighting::None }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Result of an S(q) calculation.
///
/// S(q) is a symmetric quantity with no directed counterpart, so — unlike `GrResult` —
/// only the canonical half of the pairs is kept: key `"A-B"` with A before B in
/// `GrResult.elements` order.
///
/// The weighted partials are the decomposition of the corresponding total:
/// `Σ_pairs sq_xrd[pair][qi] == total_xrd[qi]` exactly.
#[derive(Debug, Clone)]
pub struct SqResult {
    /// q values \[Å⁻¹\]
    pub q: Vec<f64>,
    /// Unweighted partial S_ij(q) — the Fourier sine transform of g_ij(r)
    pub sq: BTreeMap<String, Vec<f64>>,
    /// XRD-weighted partial w_ij(q)·S_ij(q); empty unless XRD weighting was requested
    pub sq_xrd: BTreeMap<String, Vec<f64>>,
    /// Neutron-weighted partial w_ij·S_ij(q); empty unless neutron weighting was requested
    pub sq_neutron: BTreeMap<String, Vec<f64>>,
    /// Σ over pairs of `sq_xrd` — comparable with an X-ray diffraction pattern
    pub total_xrd: Option<Vec<f64>>,
    /// Σ over pairs of `sq_neutron` — comparable with a neutron diffraction pattern
    pub total_neutron: Option<Vec<f64>>,
    pub params: SqParams,
    /// Number density used \[Å⁻³\]
    pub rho: f64,
}

// ─── 计算 ────────────────────────────────────────────────────────────────────

/// Resolve element symbol → atomic number Z, with pseudo-element fallback.
///
/// Tries exact match, then 2-char (uppercase+lowercase), then 1-char.
/// Returns 0 for unrecognised symbols.
fn elem_z_local(symbol: &str) -> usize {
    let z = ferro_core::data::elements::symbol_to_z(symbol);
    if z == 255 { 0 } else { z as usize }
}

/// Compute S(q) from a pre-computed g(r) result via Fourier sine transform.
///
/// Formula (code2 `CalcSq`):
///   S_ij(q) = 1 + (4πρ/q) Σ_r  r [g_ij(r) − 1] sin(qr) Δr
///
/// Only the canonical half of `gr.gr` is transformed — S(q) is symmetric, so the mirror
/// keys `"B-A"` would be exact duplicates of `"A-B"`.
///
/// When `params.weighting` is `Xrd`, `Neutron`, or `Both`, each pair additionally gets
/// its weighted contribution `w_ij(q)·S_ij(q)`, and their sum is the total:
///
/// S_weighted(q) = Σ_{i≤j} w_ij(q) · S_ij(q)
///
/// XRD weights:  w_ij = (2−δᵢⱼ)·cᵢcⱼfᵢ(q)fⱼ(q) / [Σₖ cₖfₖ(q)]²
/// Neutron weights: same formula with fᵢ → bcoh_i (q-independent)
///
/// Scattering factors are looked up from the atomic number, which is resolved through
/// the element prefix — so site labels (`O_f`, `P_0`) map to the same factor as their
/// element.
///
/// A label-resolved decomposition therefore reproduces the element-resolved total, but
/// only up to an O(1/N) term: same-type partials normalise by `N_A(N_A−1)` whereas summing
/// the label-level partials reconstructs `N_A²`. Pure relabelling (one label per element)
/// is exact; subdividing a site into several labels leaves a finite-size difference that
/// vanishes as the system grows.
pub fn calc_sq_from_gr(gr: &GrResult, params: &SqParams) -> SqResult {
    let n_q = ((params.q_max - params.q_min) / params.dq).floor() as usize + 1;
    let q_vals: Vec<f64> = (0..n_q)
        .map(|i| params.q_min + i as f64 * params.dq)
        .collect();

    let rho = gr.rho;
    let dr = gr.params.dr;
    let pi4 = 4.0 * std::f64::consts::PI;

    // 规范半边：gr.elements 已排序，取 i ≤ j
    let canonical: Vec<String> = (0..gr.elements.len())
        .flat_map(|i| {
            (i..gr.elements.len())
                .map(move |j| (i, j))
        })
        .map(|(i, j)| format!("{}-{}", gr.elements[i], gr.elements[j]))
        .filter(|k| gr.gr.contains_key(k))
        .collect();

    // ── 计算各 partial S_ij(q) ────────────────────────────────────────────────
    //
    // 参考实现（code2/dump2sq.c）在帧循环内部用**该帧的** ρ_f 与 g_f 各做一次变换，
    // 最后才时间平均（`CalcSq` 逐帧调用 + `TimeAverage`）。直接对时间平均后的 g 变换
    // 并不等价：NPT 下会丢掉 Cov(ρ_f, FT[g_f])。
    //
    // 但变换对 g 是线性的，而被积函数里的 ρ_f·g_f 中体积恰好抵消，于是
    //     ⟨ρ_f·(g_f − 1)⟩ = ⟨ρ_f·g_f⟩ − ⟨ρ_f⟩ = gr.rho_g − gr.rho
    // 逐帧变换的结果可以由两个时间平均量精确重构，只做一次变换即可 —— 与逐帧实现
    // 严格相等（测试 `test_folded_matches_literal_per_frame_transform` 钉住这一点），
    // 而非近似。NVT 下 `rho_g ≡ rho·gr`，整个表达式退化回 `4πρ/q·∫r(g−1)sin(qr)dr`。
    let mut sq_map: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for label in &canonical {
        let rho_g_vals = &gr.rho_g[label];
        let sq_vals: Vec<f64> = q_vals.par_iter().map(|&qi| {
            if qi.abs() < 1e-10 { return 1.0; }
            let prefactor = pi4 / qi;
            let integral: f64 = gr.r.iter().zip(rho_g_vals.iter())
                .map(|(&ri, &rg)| ri * (rg - rho) * (qi * ri).sin() * dr)
                .sum();
            1.0 + prefactor * integral
        }).collect();
        sq_map.insert(label.clone(), sq_vals);
    }

    // ── 散射因子加权 ─────────────────────────────────────────────────────────
    let want_xrd = matches!(params.weighting, SqWeighting::Xrd | SqWeighting::Both);
    let want_neu = matches!(params.weighting, SqWeighting::Neutron | SqWeighting::Both);

    let mut sq_xrd: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut sq_neutron: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut total_xrd = None;
    let mut total_neutron = None;

    let n_total: usize = gr.element_counts.values().sum();
    if (want_xrd || want_neu) && n_total > 0 {
        // 浓度 cᵢ = Nᵢ / N_total
        let elems: Vec<(&str, f64, usize)> = gr.elements.iter()
            .map(|e| {
                let c = gr.element_counts.get(e).copied().unwrap_or(0) as f64
                    / n_total as f64;
                (e.as_str(), c, elem_z_local(e))
            })
            .collect();
        warn_missing_scattering_data(&elems);

        if want_neu {
            let w = build_neutron_weights(&elems, &gr.gr);
            let (parts, total) = apply_weights(&sq_map, n_q, |key, _qi| {
                w.get(key).copied().unwrap_or(0.0)
            });
            sq_neutron = parts;
            total_neutron = Some(total);
        }
        if want_xrd {
            // XRD 权重随 q 变化，逐 q 预先构造
            let w_per_q: Vec<BTreeMap<String, f64>> = q_vals.par_iter()
                .map(|&q| build_xrd_weights_at_q(&elems, q, &gr.gr))
                .collect();
            let (parts, total) = apply_weights(&sq_map, n_q, |key, qi| {
                w_per_q[qi].get(key).copied().unwrap_or(0.0)
            });
            sq_xrd = parts;
            total_xrd = Some(total);
        }
    }

    SqResult {
        q: q_vals,
        sq: sq_map,
        sq_xrd,
        sq_neutron,
        total_xrd,
        total_neutron,
        params: params.clone(),
        rho,
    }
}

/// Multiply each partial by its weight and accumulate the total, so that
/// `Σ_pairs parts[pair][qi] == total[qi]` holds by construction.
fn apply_weights<F>(
    sq_map: &BTreeMap<String, Vec<f64>>,
    n_q: usize,
    weight: F,
) -> (BTreeMap<String, Vec<f64>>, Vec<f64>)
where
    F: Fn(&str, usize) -> f64,
{
    let mut parts: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut total = vec![0.0f64; n_q];
    for (key, vals) in sq_map {
        let mut contrib = vec![0.0f64; n_q];
        for qi in 0..n_q {
            contrib[qi] = weight(key, qi) * vals[qi];
            total[qi] += contrib[qi];
        }
        parts.insert(key.clone(), contrib);
    }
    (parts, total)
}

/// Warn when a type has no entry in the scattering-factor tables.
///
/// Reachable whenever the element column holds something the periodic table does not
/// recognise — e.g. a LAMMPS dump without an `element` column, where types degrade to
/// `X1`/`X2`. The weighted totals would otherwise be silently meaningless.
fn warn_missing_scattering_data(elems: &[(&str, f64, usize)]) {
    let missing: Vec<&str> = elems.iter()
        .filter(|(_, _, z)| *z == 0)
        .map(|(e, _, _)| *e)
        .collect();
    if !missing.is_empty() {
        eprintln!(
            "[ferro] warning: no scattering data for {} — weighted total S(q) is meaningless. \
             Check that the trajectory carries real element symbols.",
            missing.join(", ")
        );
    }
}

/// Build Faber-Ziman XRD weights w_ij at a single q value.
///
/// w_ij = (2−δᵢⱼ)·cᵢcⱼfᵢ(q)fⱼ(q) / <f(q)>²
/// Keys are the symmetric pair labels present in `gr_keys`.
fn build_xrd_weights_at_q(
    elems: &[(&str, f64, usize)],
    q: f64,
    gr_gr: &BTreeMap<String, Vec<f64>>,
) -> BTreeMap<String, f64> {
    let f: Vec<f64> = elems.iter().map(|(_, _, z)| form_factor_xrd(*z, q)).collect();
    let f_avg: f64 = elems.iter().zip(f.iter()).map(|((_, c, _), fi)| c * fi).sum();
    let denom = f_avg * f_avg;
    weights_from_factors(elems, &f, denom, gr_gr)
}

/// Build q-independent neutron weights w_ij.
fn build_neutron_weights(
    elems: &[(&str, f64, usize)],
    gr_gr: &BTreeMap<String, Vec<f64>>,
) -> BTreeMap<String, f64> {
    let b: Vec<f64> = elems.iter().map(|(_, _, z)| neutron_bcoh(*z)).collect();
    let b_avg: f64 = elems.iter().zip(b.iter()).map(|((_, c, _), bi)| c * bi).sum();
    let denom = b_avg * b_avg;
    weights_from_factors(elems, &b, denom, gr_gr)
}

/// Shared weight builder given per-element factor values and denominator.
fn weights_from_factors(
    elems: &[(&str, f64, usize)],
    factors: &[f64],
    denom: f64,
    gr_gr: &BTreeMap<String, Vec<f64>>,
) -> BTreeMap<String, f64> {
    let mut w: BTreeMap<String, f64> = BTreeMap::new();
    if denom.abs() < 1e-30 { return w; }
    // 上三角遍历（ib 从 ia 起），每个对称对只访问一次，异种对乘 (2−δᵢⱼ)。
    // 与参考实现 code2/dump2sq.c:326-338 的 `for k=j` 循环一致。
    // 全 n² 遍历会让异种对被累加两次（w_ij 翻倍、Σw > 1、加权总 S(q) 整体偏大）。
    for (ia, (ea, ca, _)) in elems.iter().enumerate() {
        for (ib, (eb, cb, _)) in elems.iter().enumerate().skip(ia) {
            // `elems` 与 `GrResult.elements` 同序，故 ia ≤ ib 时 "ea-eb" 即 g(r) 的对称对 key
            let key = format!("{ea}-{eb}");
            if !gr_gr.contains_key(&key) { continue; }
            let factor = if ia == ib { 1.0 } else { 2.0 };
            *w.entry(key).or_insert(0.0) +=
                factor * ca * cb * factors[ia] * factors[ib] / denom;
        }
    }
    w
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

impl SqResult {
    /// Projects the result into the wide table the writers consume.
    ///
    /// Columns: `q, total_xrd, total_neutron`, then `<pair>_sq / _xrd / _neutron` per pair.
    ///
    /// Wide rather than long — the mirror image of `GrResult::to_tables` — because the
    /// **primary product here is the pair of totals**, which are one value per `q`.
    /// The partials are a diagnostic decomposition (`Σ_pairs w_ij·S_ij == total`), so they
    /// spread sideways rather than forcing the totals to repeat once per pair. Trajectories
    /// with different element sets end up with different pair columns; stacking them takes
    /// the column union and leaves the gaps empty (`ferro_core::Table::concat_union`).
    ///
    /// **Every** canonical pair is written, always. There is no pair filter: keeping
    /// one pair next to the totals hides the `Σ w_ij·S_ij == total` closure that is the
    /// only reason the partials are here at all.
    pub fn to_tables(&self, gr: &GrResult) -> Result<Vec<(String, Table)>, ChemError> {
        let _ = gr;
        let keys: Vec<String> = sorted_keys(&self.sq);

        let mut t = Table::new();
        t.push_num("q", self.q.clone());
        if let Some(v) = &self.total_xrd { t.push_num("total_xrd", v.clone()); }
        if let Some(v) = &self.total_neutron { t.push_num("total_neutron", v.clone()); }

        let n = self.q.len();
        let take = |m: &BTreeMap<String, Vec<f64>>, k: &str| -> Vec<f64> {
            m.get(k).cloned().unwrap_or_else(|| vec![f64::NAN; n])
        };
        for k in &keys {
            t.push_num(format!("{k}_sq"), take(&self.sq, k));
            if self.total_xrd.is_some() {
                t.push_num(format!("{k}_xrd"), take(&self.sq_xrd, k));
            }
            if self.total_neutron.is_some() {
                t.push_num(format!("{k}_neutron"), take(&self.sq_neutron, k));
            }
        }
        Ok(vec![("sq".to_string(), t)])
    }

    /// Parameter block for the comment header, including the g(r) it was transformed from.
    pub fn meta_lines(&self, gr: &GrResult) -> Vec<String> {
        let mut v = vec!["[g(r) parameters]".to_string()];
        v.extend(gr.meta_lines());
        v.push("[S(q) parameters]".to_string());
        v.push(format!("q_min   = {} Ang^-1", self.params.q_min));
        v.push(format!("q_max   = {} Ang^-1", self.params.q_max));
        v.push(format!("dq      = {} Ang^-1", self.params.dq));
        v.push("partials: <A>-<B>_sq unweighted; _xrd / _neutron are w_ij(q)*S_ij(q),".to_string());
        v.push("          which sum over pairs to total_xrd / total_neutron".to_string());
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::gr::{GrParams, GroupBy, calc_gr};
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    fn make_sc_fe(n: usize) -> Frame {
        let a = 2.87_f64;
        let side = n as f64 * a;
        let cell = Cell::from_lengths_angles(side, side, side, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    frame.add_atom(Atom::new(
                        "Fe",
                        Vector3::new(i as f64 * a, j as f64 * a, k as f64 * a),
                    ));
                }
            }
        }
        frame
    }

    #[test]
    fn test_sq_q_axis_length() {
        let traj = Trajectory::from_frame(make_sc_fe(3));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 3.9, dr: 0.01, ..Default::default()
        }).unwrap();
        let sq_res = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 0.5, q_max: 10.0, dq: 0.1, ..Default::default()
        });
        let expected = ((10.0_f64 - 0.5) / 0.1).floor() as usize + 1;
        assert_eq!(sq_res.q.len(), expected);
    }

    #[test]
    fn test_sq_keeps_only_canonical_half_of_gr_pairs() {
        // S(q) 对称、无有向对应物，故只保留规范半边，不复制镜像键
        let traj = Trajectory::from_frame(make_multi_crystal(4, &["O", "Si"]));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 5.9, dr: 0.01, ..Default::default()
        }).unwrap();
        let sq_res = calc_sq_from_gr(&gr_res, &SqParams::default());

        // g(r) 有 n² = 4 个有序对，S(q) 只有 n(n+1)/2 = 3 个
        assert_eq!(gr_res.gr.len(), 4);
        assert_eq!(sq_res.sq.len(), 3);
        for k in ["O-O", "O-Si", "Si-Si"] {
            assert!(sq_res.sq.contains_key(k), "missing {k}");
        }
        assert!(!sq_res.sq.contains_key("Si-O"), "mirror key should be skipped");
        assert!(!sq_res.sq.contains_key("total"), "unweighted total must be gone");
    }

    #[test]
    fn test_sq_large_q_tail_near_one() {
        let traj = Trajectory::from_frame(make_sc_fe(4));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 5.5, dr: 0.01, ..Default::default()
        }).unwrap();
        let sq_res = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 1.0, q_max: 30.0, dq: 0.1, ..Default::default()
        });
        let vals = &sq_res.sq["Fe-Fe"];
        let n = vals.len();
        let tail_mean: f64 = vals[n / 2..].iter().sum::<f64>() / (n / 2) as f64;
        assert!(
            tail_mean > 0.5 && tail_mean < 1.5,
            "tail mean S(q) = {:.3}, expected ~1",
            tail_mean
        );
    }

    #[test]
    fn test_xrd_total_equals_partial_for_single_element() {
        let traj = Trajectory::from_frame(make_sc_fe(3));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 3.9, dr: 0.01, ..Default::default()
        }).unwrap();
        let sq_res = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 1.0, q_max: 10.0, dq: 0.5,
            weighting: SqWeighting::Xrd,
        });
        let xrd = sq_res.total_xrd.as_ref().expect("total_xrd should be present");
        // 单元素体系权重恒为 1，加权总应等于该唯一 partial
        let part = &sq_res.sq["Fe-Fe"];
        let max_diff = xrd.iter().zip(part.iter()).map(|(a,b)| (a-b).abs()).fold(0.0_f64, f64::max);
        assert!(max_diff < 1e-9, "single-element XRD total should equal the partial, diff={max_diff:.4}");
        assert!(sq_res.total_neutron.is_none(), "neutron not requested");
    }

    #[test]
    fn test_neutron_total_equals_partial_for_single_element() {
        let traj = Trajectory::from_frame(make_sc_fe(3));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 3.9, dr: 0.01, ..Default::default()
        }).unwrap();
        let sq_res = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 1.0, q_max: 10.0, dq: 0.5,
            weighting: SqWeighting::Neutron,
        });
        let neu = sq_res.total_neutron.as_ref().expect("total_neutron should be present");
        let part = &sq_res.sq["Fe-Fe"];
        let max_diff = neu.iter().zip(part.iter()).map(|(a,b)| (a-b).abs()).fold(0.0_f64, f64::max);
        assert!(max_diff < 1e-9, "single-element neutron total should equal the partial, diff={max_diff:.4}");
        assert!(sq_res.total_xrd.is_none(), "xrd not requested");
    }

    /// 构造多组分晶体（元素按 (i+j+k) % n_elems 轮换，a=3.0 Å）
    fn make_multi_crystal(n: usize, elems: &[&str]) -> Frame {
        let a = 3.0_f64;
        let side = n as f64 * a;
        let cell = Cell::from_lengths_angles(side, side, side, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let pos = Vector3::new(i as f64 * a, j as f64 * a, k as f64 * a);
                    frame.add_atom(Atom::new(elems[(i + j + k) % elems.len()], pos));
                }
            }
        }
        frame
    }

    /// 复现 `calc_sq_from_gr` 内部的 `elems` 构造（浓度 + 原子序数）
    fn build_elems(gr: &ferro_core::Result<crate::GrResult>) -> Vec<(&str, f64, usize)> {
        let gr = gr.as_ref().unwrap();
        let n_total: usize = gr.element_counts.values().sum();
        gr.elements.iter()
            .map(|e| {
                let c = gr.element_counts.get(e).copied().unwrap_or(0) as f64 / n_total as f64;
                (e.as_str(), c, elem_z_local(e))
            })
            .collect()
    }

    #[test]
    fn test_faber_ziman_weights_sum_to_one() {
        // Faber-Ziman 归一化：Σ_{i≤j} w_ij = (Σ cᵢfᵢ)² / <f>² = 1。
        // 回归保护：曾因全 n² 遍历 + factor=2.0 让异种对被累加两次，Σw > 1，
        // 导致 total_xrd / total_neutron 整体偏大（二元体系约 1.4 倍）。
        for elems_list in [&["O", "Si"][..], &["O", "P", "Zn"][..]] {
            let traj = Trajectory::from_frame(make_multi_crystal(4, elems_list));
            let gr_res = calc_gr(&traj, &GrParams {
                r_min: 0.1, r_max: 5.9, dr: 0.01, ..Default::default()
            });
            let elems = build_elems(&gr_res);
            let gr_gr = &gr_res.as_ref().unwrap().gr;

            let neu = build_neutron_weights(&elems, gr_gr);
            let neu_sum: f64 = neu.values().sum();
            assert!(
                (neu_sum - 1.0).abs() < 1e-10,
                "{elems_list:?}: Σ neutron w_ij = {neu_sum:.6}, expected 1"
            );

            for q in [0.5_f64, 2.0, 5.0, 12.0, 25.0] {
                let xrd = build_xrd_weights_at_q(&elems, q, gr_gr);
                let xrd_sum: f64 = xrd.values().sum();
                assert!(
                    (xrd_sum - 1.0).abs() < 1e-10,
                    "{elems_list:?}: Σ XRD w_ij(q={q}) = {xrd_sum:.6}, expected 1"
                );
            }
        }
    }

    #[test]
    fn test_hetero_pair_weight_not_double_counted() {
        // 二元体系异种对权重应为 2·c_A·c_B·f_A·f_B/<f>²，而非其两倍。
        let traj = Trajectory::from_frame(make_multi_crystal(4, &["O", "Si"]));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 5.9, dr: 0.01, ..Default::default()
        });
        let elems = build_elems(&gr_res);
        let w = build_neutron_weights(&elems, &gr_res.as_ref().unwrap().gr);

        let (_, c_o, z_o) = elems[0];
        let (_, c_si, z_si) = elems[1];
        let b_o = neutron_bcoh(z_o);
        let b_si = neutron_bcoh(z_si);
        let b_avg = c_o * b_o + c_si * b_si;
        let expected = 2.0 * c_o * c_si * b_o * b_si / (b_avg * b_avg);

        let got = w["O-Si"];
        assert!(
            (got - expected).abs() < 1e-12,
            "w(O-Si) = {got:.6}, expected {expected:.6} (2× 说明异种对被重复累加)"
        );
    }

    #[test]
    fn test_weighted_partials_sum_to_total() {
        // 加权 partial 是 total 的精确分解 —— 这正是导出这些列的意义
        let traj = Trajectory::from_frame(make_multi_crystal(4, &["O", "P", "Zn"]));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 5.9, dr: 0.01, ..Default::default()
        }).unwrap();
        let sq = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 0.5, q_max: 20.0, dq: 0.25, weighting: SqWeighting::Both,
        });

        let tx = sq.total_xrd.as_ref().unwrap();
        let tn = sq.total_neutron.as_ref().unwrap();
        for qi in 0..sq.q.len() {
            let sx: f64 = sq.sq_xrd.values().map(|v| v[qi]).sum();
            let sn: f64 = sq.sq_neutron.values().map(|v| v[qi]).sum();
            assert!((sx - tx[qi]).abs() < 1e-12, "q={}: Σ xrd partials {sx} != total {}", sq.q[qi], tx[qi]);
            assert!((sn - tn[qi]).abs() < 1e-12, "q={}: Σ neutron partials {sn} != total {}", sq.q[qi], tn[qi]);
        }
    }

    /// 按元素轮换构造晶体，并给指定元素的原子打上位点标签。
    ///
    /// `subdivide = false` 时每个元素只对应一个标签（纯改名）；
    /// `true` 时把 O 拆成两个位点。
    fn make_labelled_crystal(n: usize, subdivide: bool) -> Frame {
        let a = 3.0_f64;
        let side = n as f64 * a;
        let cell = Cell::from_lengths_angles(side, side, side, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        let mut n_o = 0usize;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let pos = Vector3::new(i as f64 * a, j as f64 * a, k as f64 * a);
                    let elem = ["O", "P", "Zn"][(i + j + k) % 3];
                    let mut atom = Atom::new(elem, pos);
                    atom.label = Some(match elem {
                        "O" if subdivide => {
                            n_o += 1;
                            if n_o.is_multiple_of(2) { "O_b_P_P" } else { "O_f" }.to_string()
                        }
                        "O" => "O_f".to_string(),
                        "P" => "P_0".to_string(),
                        _   => "Zn_f".to_string(),
                    });
                    frame.add_atom(atom);
                }
            }
        }
        frame
    }

    fn totals_by_group(frame: Frame) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let traj = Trajectory::from_frame(frame);
        let sq_params = SqParams {
            q_min: 0.5, q_max: 20.0, dq: 0.25, weighting: SqWeighting::Both,
        };
        let run = |by| {
            let gr = calc_gr(&traj, &GrParams {
                r_min: 0.1, r_max: 5.9, dr: 0.01, group_by: by,
            }).unwrap();
            let sq = calc_sq_from_gr(&gr, &sq_params);
            (sq.total_xrd.unwrap(), sq.total_neutron.unwrap())
        };
        let (ex, en) = run(GroupBy::Element);
        let (lx, ln) = run(GroupBy::Label);
        (ex, en, lx, ln)
    }

    fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()).fold(0.0_f64, f64::max)
    }

    /// 合成 NPT 轨迹：同一晶体逐帧按 `scales` 缩放，体积随之变化。
    fn make_npt_traj(n: usize, elems: &[&str], scales: &[f64]) -> Trajectory {
        let base = make_multi_crystal(n, elems);
        let base_cell = base.cell.as_ref().unwrap().matrix;
        let mut traj = Trajectory::new();
        for &s in scales {
            let mut f = Frame::with_cell(Cell::from_matrix(base_cell * s), [true; 3]);
            for a in &base.atoms {
                f.add_atom(Atom::new(&a.element, a.position * s));
            }
            traj.add_frame(f);
        }
        traj
    }

    #[test]
    fn test_folded_matches_literal_per_frame_transform() {
        // 参考实现 code2/dump2sq.c 在帧循环内逐帧做傅里叶变换（`CalcSq` 每帧一次，
        // 末尾 `TimeAverage`）。本实现只做一次变换，靠 ⟨ρ_f·g_f⟩ 重构 —— 这里证明
        // 两者严格相等，而非近似。
        //
        // 单帧调用 calc_gr + calc_sq_from_gr 就是「该帧的 g_f 与 ρ_f 各变换一次」，
        // 对全部帧取算术平均即得字面逐帧结果。
        let traj = make_npt_traj(4, &["O", "Si"], &[1.0, 1.09, 0.93, 1.13, 0.90]);
        let gp = GrParams { r_min: 0.1, r_max: 5.0, dr: 0.02, ..Default::default() };
        let sp = SqParams { q_min: 0.5, q_max: 12.0, dq: 0.1, ..Default::default() };

        let folded = calc_sq_from_gr(&calc_gr(&traj, &gp).unwrap(), &sp);

        let mut literal: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        for frame in &traj.frames {
            let single = Trajectory::from_frame(frame.clone());
            let s = calc_sq_from_gr(&calc_gr(&single, &gp).unwrap(), &sp);
            for (k, v) in &s.sq {
                let e = literal.entry(k.clone()).or_insert_with(|| vec![0.0; v.len()]);
                for (a, b) in e.iter_mut().zip(v.iter()) { *a += b; }
            }
        }
        let nf = traj.frames.len() as f64;
        for v in literal.values_mut() { for x in v.iter_mut() { *x /= nf; } }

        assert_eq!(folded.sq.len(), literal.len());
        for (k, v) in &folded.sq {
            let d = max_abs_diff(v, &literal[k]);
            assert!(d < 1e-12, "pair {k}: folded vs literal per-frame S(q) differs by {d:.3e}");
        }
    }

    #[test]
    fn test_npt_transform_differs_from_average_volume_shortcut() {
        // 反向断言：对时间平均后的 g 直接变换（旧口径 4πρ/q·∫r(g−1)sin(qr)dr，
        // ρ=N/⟨V⟩）与逐帧口径**不同**。没有这条，上面的对拍在两种口径恰好等价时
        // 也会通过。
        let traj = make_npt_traj(4, &["O", "Si"], &[1.0, 1.09, 0.93, 1.13, 0.90]);
        let gp = GrParams { r_min: 0.1, r_max: 5.0, dr: 0.02, ..Default::default() };
        let sp = SqParams { q_min: 0.5, q_max: 12.0, dq: 0.1, ..Default::default() };
        let gr = calc_gr(&traj, &gp).unwrap();
        let folded = calc_sq_from_gr(&gr, &sp);

        let pi4 = 4.0 * std::f64::consts::PI;
        let rho_old = gr.element_counts.values().sum::<usize>() as f64 / gr.avg_volume;
        let mut worst = 0.0_f64;
        for (k, gr_vals) in &gr.gr {
            if !folded.sq.contains_key(k) { continue; }
            let old: Vec<f64> = folded.q.iter().map(|&qi| {
                if qi.abs() < 1e-10 { return 1.0; }
                let integral: f64 = gr.r.iter().zip(gr_vals.iter())
                    .map(|(&ri, &g)| ri * (g - 1.0) * (qi * ri).sin() * gp.dr)
                    .sum();
                1.0 + pi4 * rho_old / qi * integral
            }).collect();
            worst = worst.max(max_abs_diff(&folded.sq[k], &old));
        }
        assert!(worst > 1e-6,
            "per-frame and average-volume S(q) must differ under NPT (max diff {worst:.3e})");
    }

    #[test]
    fn test_label_rename_gives_identical_totals() {
        // 标签只是改名（每个元素一个标签）时，加权总必须逐点相同 ——
        // 这检验散射因子是经元素前缀查表的，O_f / P_0 / Zn_f 拿到的是 O / P / Zn 的因子。
        let (ex, en, lx, ln) = totals_by_group(make_labelled_crystal(4, false));
        assert!(max_abs_diff(&ex, &lx) < 1e-12, "total_xrd differs under pure relabelling");
        assert!(max_abs_diff(&en, &ln) < 1e-12, "total_neutron differs under pure relabelling");
    }

    #[test]
    fn test_label_subdivision_differs_only_by_finite_size_term() {
        // 位点细分时两者并非逐点严格相等：同种对 g_AA 的归一化用 N_A(N_A−1)，
        // 而标签层求和重构出的是 N_A²，差一个 O(1/N_A) 的自排除项。
        // 该差随体系增大而缩小 —— 这里断言"确实随 N 变小"，把它钉成有限尺寸效应而非 bug。
        let (ex4, en4, lx4, ln4) = totals_by_group(make_labelled_crystal(4, true));
        let (ex6, en6, lx6, ln6) = totals_by_group(make_labelled_crystal(6, true));

        let dx4 = max_abs_diff(&ex4, &lx4);
        let dx6 = max_abs_diff(&ex6, &lx6);
        let dn4 = max_abs_diff(&en4, &ln4);
        let dn6 = max_abs_diff(&en6, &ln6);

        assert!(dx4 > 0.0 && dn4 > 0.0, "sanity: subdivision should not be exactly equal");
        assert!(dx6 < dx4, "xrd deviation should shrink with N: {dx6:.3e} !< {dx4:.3e}");
        assert!(dn6 < dn4, "neutron deviation should shrink with N: {dn6:.3e} !< {dn4:.3e}");
        assert!(dx4 < 0.05 && dn4 < 0.05, "deviation should stay small: xrd {dx4:.4}, neutron {dn4:.4}");
    }

    #[test]
    fn test_to_tables_always_writes_every_pair() {
        let traj = Trajectory::from_frame(make_multi_crystal(4, &["O", "Si"]));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 5.9, dr: 0.05, ..Default::default()
        }).unwrap();
        let sq = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 1.0, q_max: 10.0, dq: 0.5, weighting: SqWeighting::Both,
        });

        // 没有配对过滤:q + 2 条 total + 3 对 × 3 列 = 12 列,恒定如此。
        // 只留一对会藏起 Σ w_ij·S_ij == total 的闭合,而那是 partial 存在的唯一理由
        let (_, t) = sq.to_tables(&gr_res).unwrap().remove(0);
        assert_eq!(t.n_cols(), 12);
        assert_eq!(t.n_rows(), sq.q.len(), "宽表 → 一行一个 q,配对增多只加列");
        assert!(t.validate().is_ok());
        assert!(!t.names().contains(&"total"), "未加权 total 不该复活");
        for n in ["O-Si_sq", "O-Si_xrd", "O-Si_neutron", "O-O_sq", "Si-Si_sq"] {
            assert!(t.names().contains(&n), "{n} 应在表中");
        }
    }

    #[test]
    fn test_partials_sum_to_total_in_the_table() {
        // 0.1.10 立起来的恒等式必须在导出层依然成立
        let traj = Trajectory::from_frame(make_multi_crystal(4, &["O", "Si"]));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 5.9, dr: 0.05, ..Default::default()
        }).unwrap();
        let sq = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 1.0, q_max: 10.0, dq: 0.5, weighting: SqWeighting::Both,
        });
        let (_, t) = sq.to_tables(&gr_res).unwrap().remove(0);

        let col = |name: &str| match t.column(name).unwrap() {
            ferro_core::Column::Num(v) => v.clone(),
            _ => panic!("{name} must be numeric"),
        };
        let total = col("total_xrd");
        let parts: Vec<Vec<f64>> = t
            .names()
            .iter()
            .filter(|n| n.ends_with("_xrd") && **n != "total_xrd")
            .map(|n| col(n))
            .collect();
        for i in 0..total.len() {
            let sum: f64 = parts.iter().map(|p| p[i]).sum();
            assert!((sum - total[i]).abs() < 1e-9, "q index {i}: {sum} != {}", total[i]);
        }
    }

    #[test]
    fn test_meta_lines_carry_both_parameter_blocks() {
        let traj = Trajectory::from_frame(make_sc_fe(2));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 3.9, dr: 0.1, ..Default::default()
        }).unwrap();
        let sq_res = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 1.0, q_max: 10.0, dq: 0.5, ..Default::default()
        });
        let meta = sq_res.meta_lines(&gr_res).join("\n");
        assert!(meta.contains("[g(r) parameters]") && meta.contains("[S(q) parameters]"));
        assert!(meta.contains("q_max") && meta.contains("dr"));
    }
}
