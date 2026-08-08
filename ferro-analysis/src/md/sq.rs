//! Structure factor S(q) calculation and output.
//!
//! Workflow: compute g(r) first (`calc_gr`), then call `calc_sq_from_gr` to obtain S(q).
//! Output: two files — `.gr` (written by the gr module) and `.sq` (written here).
//!
//! Formula (Faber-Ziman, from code2/dump2sq.c CalcSq):
//!   S_ij(q) = 1 + (4πρ/q) Σ_r r[g_ij(r)−1] sin(qr) Δr

use super::gr::{sorted_keys, GrResult, VERSION};
use super::scattering_data::{form_factor_xrd, neutron_bcoh};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::io::{BufWriter, Write};

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
    /// q step size \[Å⁻¹\] (default: 0.05)
    pub dq: f64,
    /// Scattering factor weighting for weighted total S(q) (default: None)
    pub weighting: SqWeighting,
}

impl Default for SqParams {
    fn default() -> Self {
        SqParams { q_min: 0.1, q_max: 25.0, dq: 0.05, weighting: SqWeighting::None }
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
    let mut sq_map: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for label in &canonical {
        let gr_vals = &gr.gr[label];
        let sq_vals: Vec<f64> = q_vals.par_iter().map(|&qi| {
            if qi.abs() < 1e-10 { return 1.0; }
            let prefactor = pi4 * rho / qi;
            let integral: f64 = gr.r.iter().zip(gr_vals.iter())
                .map(|(&ri, &gri)| ri * (gri - 1.0) * (qi * ri).sin() * dr)
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

/// Write S(q) data to a tab-separated text file (`.sq`).
///
/// The header records both g(r) parameters (used as input) and S(q) parameters.
///
/// Columns: `q[Ang^-1]`, then `total_xrd` / `total_neutron` where computed, then per
/// pair a `_sq` / `_xrd` / `_neutron` triple ordered by (Z, label).
///
/// - `pair = None` → every canonical pair.
/// - `pair = Some((a, b))` → only that pair. Used with label-resolved partials, where
///   the pair count grows quadratically in the number of site labels and a full table
///   would run to hundreds of columns.
pub fn write_sq(
    gr: &GrResult,
    sq: &SqResult,
    path: &str,
    pair: Option<(&str, &str)>,
) -> std::io::Result<()> {
    let keys: Vec<String> = match pair {
        Some((a, b)) => {
            let ab = format!("{a}-{b}");
            let ba = format!("{b}-{a}");
            let key = if sq.sq.contains_key(&ab) {
                ab
            } else if sq.sq.contains_key(&ba) {
                ba
            } else {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    format!("pair '{ab}' not present in the trajectory"),
                ));
            };
            vec![key]
        }
        None => sorted_keys(&sq.sq),
    };

    let mut w = BufWriter::new(std::fs::File::create(path)?);

    // 文件头
    writeln!(w, "# ferro v{}", VERSION)?;
    writeln!(w, "# Structure Factor S(q) [computed from g(r) via Fourier sine transform]")?;
    writeln!(w, "# {}", "-".repeat(60))?;
    // g(r) 计算参数
    writeln!(w, "# [g(r) parameters]")?;
    writeln!(w, "# r_min   = {} Ang", gr.params.r_min)?;
    writeln!(w, "# r_max   = {} Ang", gr.params.r_max)?;
    writeln!(w, "# dr      = {} Ang", gr.params.dr)?;
    writeln!(w, "# frames  = {}", gr.n_frames)?;
    writeln!(w, "# volume  = {:.3} Ang^3", gr.avg_volume)?;
    writeln!(w, "# density = {:.6e} Ang^-3", gr.rho)?;
    writeln!(w, "# atoms:")?;
    for elem in &gr.elements {
        let count = gr.element_counts.get(elem).copied().unwrap_or(0);
        writeln!(w, "#   {:<10}: {}", elem, count)?;
    }
    // S(q) 参数
    writeln!(w, "# [S(q) parameters]")?;
    writeln!(w, "# q_min   = {} Ang^-1", sq.params.q_min)?;
    writeln!(w, "# q_max   = {} Ang^-1", sq.params.q_max)?;
    writeln!(w, "# dq      = {} Ang^-1", sq.params.dq)?;
    writeln!(w, "# partials: <A>-<B>_sq unweighted; _xrd / _neutron are w_ij(q)*S_ij(q),")?;
    writeln!(w, "#           which sum over pairs to total_xrd / total_neutron")?;
    writeln!(w, "# {}", "-".repeat(60))?;

    // 列标题
    write!(w, "# q[Ang^-1]")?;
    if sq.total_xrd.is_some() { write!(w, "\ttotal_xrd")?; }
    if sq.total_neutron.is_some() { write!(w, "\ttotal_neutron")?; }
    for k in &keys {
        write!(w, "\t{k}_sq")?;
        if sq.total_xrd.is_some() { write!(w, "\t{k}_xrd")?; }
        if sq.total_neutron.is_some() { write!(w, "\t{k}_neutron")?; }
    }
    writeln!(w)?;

    // 数据
    for i in 0..sq.q.len() {
        write!(w, "{:.6e}", sq.q[i])?;
        if let Some(v) = &sq.total_xrd { write!(w, "\t{:.6e}", v[i])?; }
        if let Some(v) = &sq.total_neutron { write!(w, "\t{:.6e}", v[i])?; }
        for k in &keys {
            write!(w, "\t{:.6e}", sq.sq.get(k).map(|v| v[i]).unwrap_or(0.0))?;
            if sq.total_xrd.is_some() {
                write!(w, "\t{:.6e}", sq.sq_xrd.get(k).map(|v| v[i]).unwrap_or(0.0))?;
            }
            if sq.total_neutron.is_some() {
                write!(w, "\t{:.6e}", sq.sq_neutron.get(k).map(|v| v[i]).unwrap_or(0.0))?;
            }
        }
        writeln!(w)?;
    }
    Ok(())
}

// ─── 测试 ────────────────────────────────────────────────────────────────────

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
    fn test_write_sq_pair_selection_columns() {
        use std::io::Read;
        let traj = Trajectory::from_frame(make_multi_crystal(4, &["O", "Si"]));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 5.9, dr: 0.05, ..Default::default()
        }).unwrap();
        let sq = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 1.0, q_max: 10.0, dq: 0.5, weighting: SqWeighting::Both,
        });

        // 指定单对：q + 2 条 total + 该对 3 列 = 6 列，参数顺序无关
        let p1 = "/tmp/test_ferro_pair.sq";
        write_sq(&gr_res, &sq, p1, Some(("Si", "O"))).unwrap();
        let mut s = String::new();
        std::fs::File::open(p1).unwrap().read_to_string(&mut s).unwrap();
        assert!(s.contains("# q[Ang^-1]\ttotal_xrd\ttotal_neutron\tO-Si_sq\tO-Si_xrd\tO-Si_neutron\n"));
        let row = s.lines().find(|l| !l.starts_with('#')).unwrap();
        assert_eq!(row.split('\t').count(), 6);

        // 全部对：q + 2 条 total + 3 对 × 3 列 = 12 列
        let p2 = "/tmp/test_ferro_all.sq";
        write_sq(&gr_res, &sq, p2, None).unwrap();
        let mut s2 = String::new();
        std::fs::File::open(p2).unwrap().read_to_string(&mut s2).unwrap();
        let row2 = s2.lines().find(|l| !l.starts_with('#')).unwrap();
        assert_eq!(row2.split('\t').count(), 12);
        assert!(!s2.contains("\ttotal\t") && !s2.trim_end().ends_with("\ttotal"));

        assert!(write_sq(&gr_res, &sq, "/tmp/test_ferro_bad.sq", Some(("Zn", "O"))).is_err());
    }

    #[test]
    fn test_write_sq() {
        use std::io::Read;
        let traj = Trajectory::from_frame(make_sc_fe(2));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 3.9, dr: 0.1, ..Default::default()
        }).unwrap();
        let sq_res = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 1.0, q_max: 10.0, dq: 0.5, ..Default::default()
        });
        let sq_path = "/tmp/test_ferro.sq";
        write_sq(&gr_res, &sq_res, sq_path, None).expect("write_sq failed");

        let mut content = String::new();
        std::fs::File::open(sq_path).unwrap().read_to_string(&mut content).unwrap();
        assert!(content.starts_with("# ferro v"));
        assert!(content.contains("# q[Ang^-1]"));
        assert!(content.contains("Fourier sine transform"));
    }
}
