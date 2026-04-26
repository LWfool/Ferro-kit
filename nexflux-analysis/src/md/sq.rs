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
/// Key semantics match `GrResult.gr`:
/// - `"El1-El2"` — partial S_ij(q) with Z(El1) ≤ Z(El2)
/// - `"total"`   — S(q) computed from the total g(r)
#[derive(Debug, Clone)]
pub struct SqResult {
    /// q values \[Å⁻¹\]
    pub q: Vec<f64>,
    /// Partial and total S(q): same pair keys as `GrResult.gr`
    pub sq: BTreeMap<String, Vec<f64>>,
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
    use nexflux_core::data::elements::by_symbol;
    if let Some(e) = by_symbol(symbol) { return e.atomic_number as usize; }
    let b = symbol.as_bytes();
    if b.len() >= 2 {
        let two = [b[0].to_ascii_uppercase(), b[1].to_ascii_lowercase()];
        if let Ok(s) = std::str::from_utf8(&two) {
            if let Some(e) = by_symbol(s) { return e.atomic_number as usize; }
        }
    }
    if let Some(&first) = b.first() {
        let one = [first.to_ascii_uppercase()];
        if let Ok(s) = std::str::from_utf8(&one) {
            if let Some(e) = by_symbol(s) { return e.atomic_number as usize; }
        }
    }
    0
}

/// Compute S(q) from a pre-computed g(r) result via Fourier sine transform.
///
/// Formula (code2 `CalcSq`):
///   S_ij(q) = 1 + (4πρ/q) Σ_r  r [g_ij(r) − 1] sin(qr) Δr
///
/// Applies to all pairs present in `gr.gr`, including `"total"`.
///
/// When `params.weighting` is `Xrd`, `Neutron`, or `Both`, the result also
/// contains weighted total S(q) keys `"total_xrd"` and/or `"total_neutron"`:
///
/// S_weighted(q) = Σ_{i≤j} w_ij(q) · S_ij(q)
///
/// XRD weights:  w_ij = (2−δᵢⱼ)·cᵢcⱼfᵢ(q)fⱼ(q) / [Σₖ cₖfₖ(q)]²
/// Neutron weights: same formula with fᵢ → bcoh_i (q-independent)
pub fn calc_sq_from_gr(gr: &GrResult, params: &SqParams) -> SqResult {
    let n_q = ((params.q_max - params.q_min) / params.dq).floor() as usize + 1;
    let q_vals: Vec<f64> = (0..n_q)
        .map(|i| params.q_min + i as f64 * params.dq)
        .collect();

    let rho = gr.rho;
    let dr = gr.params.dr;
    let pi4 = 4.0 * std::f64::consts::PI;

    let mut sq_map: BTreeMap<String, Vec<f64>> = BTreeMap::new();

    // ── 计算各 partial S_ij(q) ────────────────────────────────────────────────
    for (label, gr_vals) in &gr.gr {
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

    // ── 散射因子加权总 S(q) ──────────────────────────────────────────────────
    let want_xrd = matches!(params.weighting, SqWeighting::Xrd | SqWeighting::Both);
    let want_neu = matches!(params.weighting, SqWeighting::Neutron | SqWeighting::Both);

    if want_xrd || want_neu {
        // 元素浓度 cᵢ = Nᵢ / N_total
        let n_total: usize = gr.element_counts.values().sum();
        if n_total > 0 {
            let elems: Vec<(&str, f64, usize)> = gr.elements.iter()
                .map(|e| {
                    let c = gr.element_counts.get(e).copied().unwrap_or(0) as f64
                        / n_total as f64;
                    let z = elem_z_local(e);
                    (e.as_str(), c, z)
                })
                .collect();

            // Neutron 权重与 q 无关，预先计算
            let neu_weights: Option<BTreeMap<String, f64>> = if want_neu {
                Some(build_neutron_weights(&elems, &gr.gr))
            } else {
                None
            };

            let (xrd_total, neu_total): (Option<Vec<f64>>, Option<Vec<f64>>) = if want_xrd && want_neu {
                let pairs: Vec<(f64, f64)> = (0..n_q).into_par_iter().map(|qi| {
                    let xw = build_xrd_weights_at_q(&elems, q_vals[qi], &gr.gr);
                    let xv = weighted_total_indexed(&xw, &sq_map, qi);
                    let nv = weighted_total_indexed(neu_weights.as_ref().unwrap(), &sq_map, qi);
                    (xv, nv)
                }).collect();
                let (x, n) = pairs.into_iter().unzip();
                (Some(x), Some(n))
            } else if want_xrd {
                let v: Vec<f64> = (0..n_q).into_par_iter().map(|qi| {
                    let xw = build_xrd_weights_at_q(&elems, q_vals[qi], &gr.gr);
                    weighted_total_indexed(&xw, &sq_map, qi)
                }).collect();
                (Some(v), None)
            } else {
                let v: Vec<f64> = (0..n_q).into_par_iter().map(|qi| {
                    weighted_total_indexed(neu_weights.as_ref().unwrap(), &sq_map, qi)
                }).collect();
                (None, Some(v))
            };

            if let Some(v) = xrd_total { sq_map.insert("total_xrd".to_string(), v); }
            if let Some(v) = neu_total { sq_map.insert("total_neutron".to_string(), v); }
        }
    }

    SqResult { q: q_vals, sq: sq_map, params: params.clone(), rho }
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
    for (ia, (ea, ca, _)) in elems.iter().enumerate() {
        for (ib, (eb, cb, _)) in elems.iter().enumerate() {
            // 仅处理对称对 key（Z(a) ≤ Z(b)，或 a==b）
            let key = canonical_pair_key(ea, eb, elems[ia].2, elems[ib].2);
            if !gr_gr.contains_key(&key) { continue; }
            let factor = if ia == ib { 1.0 } else { 2.0 };
            *w.entry(key).or_insert(0.0) +=
                factor * ca * cb * factors[ia] * factors[ib] / denom;
        }
    }
    w
}

/// Compute the canonical symmetric pair key "El1-El2" with Z(El1) ≤ Z(El2).
fn canonical_pair_key(ea: &str, eb: &str, za: usize, zb: usize) -> String {
    if za < zb || (za == zb && ea <= eb) {
        format!("{}-{}", ea, eb)
    } else {
        format!("{}-{}", eb, ea)
    }
}

/// Evaluate Σ_pairs w_pair · S_pair[qi] for a single q index.
fn weighted_total_indexed(
    weights: &BTreeMap<String, f64>,
    sq_map: &BTreeMap<String, Vec<f64>>,
    qi: usize,
) -> f64 {
    weights.iter()
        .filter_map(|(key, &w)| {
            sq_map.get(key).map(|v| w * v[qi])
        })
        .sum()
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

/// Write S(q) data to a tab-separated text file (`.sq`).
///
/// The header records both g(r) parameters (used as input) and S(q) parameters.
/// Columns follow the same periodic-table ordering as the `.gr` file.
pub fn write_sq(gr: &GrResult, sq: &SqResult, path: &str) -> std::io::Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(path)?);

    // 文件头
    writeln!(w, "# nexflux v{}", VERSION)?;
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
        writeln!(w, "#   {:<4}: {}", elem, count)?;
    }
    // S(q) 参数
    writeln!(w, "# [S(q) parameters]")?;
    writeln!(w, "# q_min   = {} Ang^-1", sq.params.q_min)?;
    writeln!(w, "# q_max   = {} Ang^-1", sq.params.q_max)?;
    writeln!(w, "# dq      = {} Ang^-1", sq.params.dq)?;
    writeln!(w, "# {}", "-".repeat(60))?;

    // 列标题
    let keys = sorted_keys(&sq.sq);
    write!(w, "# q[Ang^-1]")?;
    for k in &keys { write!(w, "\t{}", k)?; }
    writeln!(w)?;

    // 数据
    let n = sq.q.len();
    for i in 0..n {
        write!(w, "{:.6e}", sq.q[i])?;
        for k in &keys {
            let v = sq.sq.get(k).map(|v| v[i]).unwrap_or(0.0);
            write!(w, "\t{:.6e}", v)?;
        }
        writeln!(w)?;
    }
    Ok(())
}

// ─── 测试 ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::gr::{GrParams, calc_gr};
    use nexflux_core::{Atom, Cell, Frame, Trajectory};
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
            r_min: 0.1, r_max: 3.9, dr: 0.01, r_cut: 3.0,
        }).unwrap();
        let sq_res = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 0.5, q_max: 10.0, dq: 0.1, ..Default::default()
        });
        let expected = ((10.0_f64 - 0.5) / 0.1).floor() as usize + 1;
        assert_eq!(sq_res.q.len(), expected);
    }

    #[test]
    fn test_sq_pair_keys_match_gr() {
        // S(q) 的 pair key 应与 g(r) 一致
        let traj = Trajectory::from_frame(make_sc_fe(3));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 3.9, dr: 0.01, r_cut: 3.0,
        }).unwrap();
        let sq_res = calc_sq_from_gr(&gr_res, &SqParams::default());

        for key in gr_res.gr.keys() {
            assert!(sq_res.sq.contains_key(key), "missing key {} in sq", key);
        }
    }

    #[test]
    fn test_sq_large_q_tail_near_one() {
        let traj = Trajectory::from_frame(make_sc_fe(4));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 5.5, dr: 0.01, r_cut: 3.0,
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
    fn test_xrd_weighted_adds_total_xrd_key() {
        let traj = Trajectory::from_frame(make_sc_fe(3));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 3.9, dr: 0.01, r_cut: 3.0,
        }).unwrap();
        let sq_res = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 1.0, q_max: 10.0, dq: 0.5,
            weighting: SqWeighting::Xrd,
        });
        assert!(sq_res.sq.contains_key("total_xrd"), "should contain total_xrd key");
        // For single-element system total_xrd ≈ total (weights are all 1)
        let xrd = &sq_res.sq["total_xrd"];
        let tot = &sq_res.sq["total"];
        let max_diff = xrd.iter().zip(tot.iter()).map(|(a,b)| (a-b).abs()).fold(0.0_f64, f64::max);
        assert!(max_diff < 0.01, "single-element XRD total should match unweighted total, diff={:.4}", max_diff);
    }

    #[test]
    fn test_neutron_weighted_adds_total_neutron_key() {
        let traj = Trajectory::from_frame(make_sc_fe(3));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 3.9, dr: 0.01, r_cut: 3.0,
        }).unwrap();
        let sq_res = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 1.0, q_max: 10.0, dq: 0.5,
            weighting: SqWeighting::Neutron,
        });
        assert!(sq_res.sq.contains_key("total_neutron"), "should contain total_neutron key");
        let neu = &sq_res.sq["total_neutron"];
        let tot = &sq_res.sq["total"];
        let max_diff = neu.iter().zip(tot.iter()).map(|(a,b)| (a-b).abs()).fold(0.0_f64, f64::max);
        assert!(max_diff < 0.01, "single-element neutron total should match unweighted total, diff={:.4}", max_diff);
    }

    #[test]
    fn test_write_sq() {
        use std::io::Read;
        let traj = Trajectory::from_frame(make_sc_fe(2));
        let gr_res = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 3.9, dr: 0.1, r_cut: 3.0,
        }).unwrap();
        let sq_res = calc_sq_from_gr(&gr_res, &SqParams {
            q_min: 1.0, q_max: 10.0, dq: 0.5, ..Default::default()
        });
        let sq_path = "/tmp/test_nexflux.sq";
        write_sq(&gr_res, &sq_res, sq_path).expect("write_sq failed");

        let mut content = String::new();
        std::fs::File::open(sq_path).unwrap().read_to_string(&mut content).unwrap();
        assert!(content.starts_with("# nexflux v"));
        assert!(content.contains("# q[Ang^-1]"));
        assert!(content.contains("Fourier sine transform"));
    }
}
