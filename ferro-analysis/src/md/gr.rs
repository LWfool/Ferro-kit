//! Radial distribution function g(r) calculation and output.
//!
//! Design reference:
//!   - Parameter style follows code1/gr.c (GrParams fields map to DR/RMAX/RMIN/RCUT macros)
//!   - Multi-component partial g(r) and directed CN follow code2/dump2sq.c CalcCn_pp + CalcGr
//!
//! Column ordering: sorted by atomic number (periodic table order).
//!
//! Two key semantics:
//!   - `gr`  — symmetric pair `"El1-El2"` with Z(El1) ≤ Z(El2)
//!   - `cn`  — **directed pair** `"center-neighbor"`; for A≠B, `"A-B"` and `"B-A"` are independent with different CN values

use ferro_core::Trajectory;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::io::{BufWriter, Write};

/// ferro package version (from Cargo.toml)
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// ─── 内部辅助 ────────────────────────────────────────────────────────────────

/// Look up the atomic number for an element symbol, including pseudo-element labels (e.g. "Onb", "P1", "Fe2").
///
/// Matching strategy (tried in order):
///   1. Exact match ("Fe" → 26, "O" → 8).
///   2. First two bytes (uppercase + lowercase) as a chemical symbol ("Fe1" → try "Fe" → 26).
///   3. First byte (uppercase letter) as a single-character element ("Onb" → "O" → 8, "P1" → "P" → 15).
///   4. Unrecognised → 255 (sorted to end, with string secondary ordering).
pub(super) fn elem_z(symbol: &str) -> u8 {
    use ferro_core::data::elements::by_symbol;
    if let Some(e) = by_symbol(symbol) { return e.atomic_number; }
    let b = symbol.as_bytes();
    if b.len() >= 2 && b[0].is_ascii_uppercase() && b[1].is_ascii_lowercase() {
        if let Some(e) = by_symbol(&symbol[..2]) { return e.atomic_number; }
    }
    if !b.is_empty() && b[0].is_ascii_uppercase() {
        if let Some(e) = by_symbol(&symbol[..1]) { return e.atomic_number; }
    }
    255
}

/// Split a `"El1-El2"` key into `(El1, El2)`.
fn split_pair(key: &str) -> (&str, &str) {
    key.split_once('-').unwrap_or((key, ""))
}

/// Sort map keys by (Z(left), left_label, Z(right), right_label), placing "total" last.
///
/// The secondary string sort ensures deterministic ordering when multiple pseudo-element
/// labels share the same base atomic number (e.g., "Ob" and "Onb" both derived from O).
pub fn sorted_keys(map: &BTreeMap<String, Vec<f64>>) -> Vec<String> {
    let mut keys: Vec<String> = map.keys()
        .filter(|k| k.as_str() != "total")
        .cloned()
        .collect();
    keys.sort_by(|ka, kb| {
        let (a1, b1) = split_pair(ka);
        let (a2, b2) = split_pair(kb);
        (elem_z(a1), a1, elem_z(b1), b1)
            .cmp(&(elem_z(a2), a2, elem_z(b2), b2))
    });
    keys.push("total".to_string());
    keys
}

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Parameters for g(r) calculation.
#[derive(Debug, Clone)]
pub struct GrParams {
    /// Minimum distance \[Å\] (default: 0.005)
    pub r_min: f64,
    /// Maximum distance \[Å\]; must satisfy r_max < L_min/2 (default: 10.005)
    pub r_max: f64,
    /// Distance bin width \[Å\] (default: 0.01)
    pub dr: f64,
    /// Cut-off radius \[Å\] for short-range bond distance statistics (default: 2.3)
    pub r_cut: f64,
}

impl Default for GrParams {
    fn default() -> Self {
        GrParams { r_min: 0.005, r_max: 10.005, dr: 0.01, r_cut: 2.3 }
    }
}

impl GrParams {
    /// Determine `r_max` automatically as half the shortest cell vector in the first frame.
    pub fn with_auto_rmax(traj: &Trajectory) -> Self {
        let r_max = traj
            .first()
            .and_then(|f| f.cell.as_ref())
            .map(|cell| {
                let l = cell.lengths();
                l[0].min(l[1]).min(l[2]) * 0.5
            })
            .unwrap_or(10.005);
        GrParams { r_max, ..Default::default() }
    }
}

// ─── 结果结构体 ──────────────────────────────────────────────────────────────

/// Short-range bond distance statistics within `r_cut`.
#[derive(Debug, Clone)]
pub struct PairStats {
    /// Mean bond distance \[Å\]
    pub mean_dist: f64,
    /// Standard deviation \[Å\]
    pub std_dist: f64,
    /// Total number of counted bonds across all frames
    pub count: usize,
}

/// Result of a g(r) / CN(r) calculation.
///
/// # Key conventions
/// - `gr` — symmetric pairs `"El1-El2"` with Z(El1) ≤ Z(El2), plus `"total"`.
/// - `cn` — **directed** cumulative coordination numbers:
///   key `"A-B"` = average number of B atoms within r from each A atom.
///   For A ≠ B, both `"A-B"` and `"B-A"` are present (and generally unequal).
#[derive(Debug, Clone)]
pub struct GrResult {
    /// Bin-centre r values \[Å\]
    pub r: Vec<f64>,
    /// Symmetric partial g(r): key `"El1-El2"` (Z(El1)≤Z(El2)) and `"total"`
    pub gr: BTreeMap<String, Vec<f64>>,
    /// Directed cumulative CN(r): key `"center-neighbor"` sorted by (Z_center, Z_neighbor)
    pub cn: BTreeMap<String, Vec<f64>>,
    /// Bond distance statistics within r_cut (symmetric key `"El1-El2"`)
    pub pair_stats: BTreeMap<String, PairStats>,
    /// Element types present, sorted by atomic number
    pub elements: Vec<String>,
    /// Atom count per element type (from the first frame)
    pub element_counts: BTreeMap<String, usize>,
    pub n_frames: usize,
    /// Average box volume \[Å³\]
    pub avg_volume: f64,
    /// Total number density N/V \[Å⁻³\]
    pub rho: f64,
    pub params: GrParams,
}

// ─── 计算函数 ────────────────────────────────────────────────────────────────

/// Compute partial and total g(r) and directed CN(r) for all element pairs.
///
/// This is the `--all` mode: all element pairs are calculated simultaneously.
/// Requires periodic cells in every frame; uses the minimum-image convention.
/// Returns `None` when the trajectory is empty or any frame has no cell.
pub fn calc_gr(traj: &Trajectory, params: &GrParams) -> Option<GrResult> {
    if traj.frames.is_empty() { return None; }

    let n_bins = ((params.r_max - params.r_min) / params.dr).floor() as usize;
    if n_bins == 0 { return None; }

    // ── 元素列表，按原子序数升序 ─────────────────────────────────────────
    let first_frame = traj.frames.first()?;
    let mut elem_set = std::collections::HashSet::new();
    for a in &first_frame.atoms { elem_set.insert(a.element.clone()); }
    let mut elements: Vec<String> = elem_set.into_iter().collect();
    elements.sort_by_key(|e| elem_z(e));
    let n_types = elements.len();

    let elem_idx: std::collections::HashMap<&str, usize> = elements
        .iter().enumerate().map(|(i, e)| (e.as_str(), i)).collect();

    // ── 正则对列表 (ti ≤ tj) ─────────────────────────────────────────────
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..n_types {
        for j in i..n_types { pairs.push((i, j)); }
    }
    let n_pairs = pairs.len();

    // O(1) 对索引查表
    let mut pair_lookup = vec![vec![0usize; n_types]; n_types];
    for (pidx, &(i, j)) in pairs.iter().enumerate() {
        pair_lookup[i][j] = pidx;
        pair_lookup[j][i] = pidx;
    }

    // 正则对标签（g(r) 文件列名）
    let sym_labels: Vec<String> = pairs.iter()
        .map(|&(i, j)| format!("{}-{}", elements[i], elements[j]))
        .collect();

    // 先校验所有帧都有 cell（并行路径内无法使用 ?）
    if traj.frames.iter().any(|f| f.cell.is_none()) { return None; }

    // 并行逐帧计数：每个线程独立维护 (hist, cut_dists, volume)，最后 reduce 合并
    let init = || (
        vec![vec![0.0f64; n_bins]; n_pairs + 1],
        vec![Vec::<f64>::new(); n_pairs],
        0.0f64,
    );
    let (hist, cut_dists, total_volume) = traj.frames.par_iter()
        .fold(init, |(mut h, mut c, mut v), frame| {
            let cell = frame.cell.as_ref().unwrap();
            v += cell.volume();
            let n_atoms = frame.atoms.len();
            for i in 0..n_atoms {
                for j in (i + 1)..n_atoms {
                    let ai = &frame.atoms[i];
                    let aj = &frame.atoms[j];
                    let diff = cell.minimum_image(aj.position - ai.position);
                    let r_val = diff.norm();
                    if r_val < params.r_min || r_val >= params.r_max { continue; }
                    let bin = ((r_val - params.r_min) / params.dr) as usize;
                    if bin >= n_bins { continue; }
                    let Some(&ti) = elem_idx.get(ai.element.as_str()) else { continue; };
                    let Some(&tj) = elem_idx.get(aj.element.as_str()) else { continue; };
                    let pidx = pair_lookup[ti.min(tj)][ti.max(tj)];
                    h[pidx][bin] += 1.0;
                    h[n_pairs][bin] += 1.0;
                    if r_val < params.r_cut { c[pidx].push(r_val); }
                }
            }
            (h, c, v)
        })
        .reduce(init, |(mut ha, mut ca, va), (hb, cb, vb)| {
            for p in 0..ha.len() {
                for b in 0..ha[p].len() { ha[p][b] += hb[p][b]; }
            }
            for (dst, src) in ca.iter_mut().zip(cb) { dst.extend(src); }
            (ha, ca, va + vb)
        });

    // ── normalisation ───────────────────────────────────────────────────────
    let n_frames = traj.frames.len();
    let avg_volume = total_volume / n_frames as f64;
    if avg_volume <= 0.0 { return None; }

    // atom counts per element from first frame
    let mut type_counts = vec![0.0f64; n_types];
    let mut element_counts: BTreeMap<String, usize> = BTreeMap::new();
    for atom in &first_frame.atoms {
        if let Some(&ti) = elem_idx.get(atom.element.as_str()) {
            type_counts[ti] += 1.0;
            *element_counts.entry(atom.element.clone()).or_insert(0) += 1;
        }
    }
    let n_total = first_frame.atoms.len() as f64;
    let rho = n_total / avg_volume;

    let r_centers: Vec<f64> = (0..n_bins)
        .map(|i| params.r_min + (i as f64 + 0.5) * params.dr)
        .collect();

    let pi4 = 4.0 * std::f64::consts::PI;
    let mut gr_map: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut cn_map: BTreeMap<String, Vec<f64>> = BTreeMap::new();

    // ── partial g(r) + directed CN ──────────────────────────────────────
    // Normalisation formula (code2 CalcGr):
    //   Same type (A=A): g = 2·cr / (4πr²Δr·(N_A-1)/V·N_A·steps)
    //   Different type (A-B): g = cr  / (4πr²Δr·N_B/V·N_A·steps)
    //
    // Directed CN:
    //   CN("A-B") = Σ_{bins} [count_AB / N_A / steps] * multiplicity
    //   CN("B-A") = Σ_{bins} [count_AB / N_B / steps]  (heterogeneous pairs only)
    for (pidx, &(ti, tj)) in pairs.iter().enumerate() {
        let n_a = type_counts[ti];
        let n_b = type_counts[tj];
        if n_a < 1.0 || n_b < 1.0 { continue; }

        let same_type = ti == tj;
        let (ni, n_center, n_neighbor) = if same_type {
            (2.0f64, n_a, n_a - 1.0)
        } else {
            (1.0f64, n_a, n_b)
        };

        // 对称 g(r)
        let mut gr_vec = vec![0.0f64; n_bins];
        for (i, &r_c) in r_centers.iter().enumerate() {
            let bunbo = pi4 * r_c * r_c * params.dr * n_neighbor / avg_volume * n_center;
            if bunbo > 0.0 {
                gr_vec[i] = ni * hist[pidx][i] / (bunbo * n_frames as f64);
            }
        }
        gr_map.insert(sym_labels[pidx].clone(), gr_vec);

        // 有向 CN("A-B")：center=A（低 Z），neighbor=B（高 Z 或同种）
        {
            let mut cn_vec = vec![0.0f64; n_bins];
            let mut running = 0.0f64;
            for i in 0..n_bins {
                running += ni * hist[pidx][i] / (n_center * n_frames as f64);
                cn_vec[i] = running;
            }
            cn_map.insert(sym_labels[pidx].clone(), cn_vec);
        }

        // 有向 CN("B-A")：center=B（高 Z），neighbor=A — 仅对异种对
        if !same_type {
            let label_ba = format!("{}-{}", elements[tj], elements[ti]);
            let mut cn_vec = vec![0.0f64; n_bins];
            let mut running = 0.0f64;
            for i in 0..n_bins {
                running += hist[pidx][i] / (n_b * n_frames as f64);
                cn_vec[i] = running;
            }
            cn_map.insert(label_ba, cn_vec);
        }
    }

    // total（全原子视为同种）
    {
        let mut gr_t = vec![0.0f64; n_bins];
        let mut cn_t = vec![0.0f64; n_bins];
        let mut running = 0.0f64;
        for (i, &r_c) in r_centers.iter().enumerate() {
            let bunbo = pi4 * r_c * r_c * params.dr * (n_total - 1.0) / avg_volume * n_total;
            if bunbo > 0.0 {
                gr_t[i] = 2.0 * hist[n_pairs][i] / (bunbo * n_frames as f64);
            }
            running += 2.0 * hist[n_pairs][i] / (n_total * n_frames as f64);
            cn_t[i] = running;
        }
        gr_map.insert("total".to_string(), gr_t);
        cn_map.insert("total".to_string(), cn_t);
    }

    // 键长统计（正则对标签）
    let mut pair_stats: BTreeMap<String, PairStats> = BTreeMap::new();
    for (pidx, dists) in cut_dists.iter().enumerate() {
        if dists.is_empty() { continue; }
        let n = dists.len() as f64;
        let mean = dists.iter().sum::<f64>() / n;
        let var = dists.iter().map(|&d| (d - mean).powi(2)).sum::<f64>() / n;
        pair_stats.insert(
            sym_labels[pidx].clone(),
            PairStats { mean_dist: mean, std_dist: var.sqrt(), count: dists.len() },
        );
    }

    Some(GrResult {
        r: r_centers,
        gr: gr_map,
        cn: cn_map,
        pair_stats,
        elements,
        element_counts,
        n_frames,
        avg_volume,
        rho,
        params: params.clone(),
    })
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

/// Write header lines common to both .gr and .cn files.
fn write_header(
    w: &mut impl Write,
    result: &GrResult,
    analysis_name: &str,
) -> std::io::Result<()> {
    writeln!(w, "# ferro v{}", VERSION)?;
    writeln!(w, "# {}", analysis_name)?;
    writeln!(w, "# {}", "-".repeat(60))?;
    writeln!(w, "# r_min   = {} Ang", result.params.r_min)?;
    writeln!(w, "# r_max   = {} Ang", result.params.r_max)?;
    writeln!(w, "# dr      = {} Ang", result.params.dr)?;
    writeln!(w, "# r_cut   = {} Ang", result.params.r_cut)?;
    writeln!(w, "# frames  = {}", result.n_frames)?;
    writeln!(w, "# volume  = {:.3} Ang^3", result.avg_volume)?;
    writeln!(w, "# density = {:.6e} Ang^-3", result.rho)?;
    writeln!(w, "# atoms:")?;
    for elem in &result.elements {
        let count = result.element_counts.get(elem).copied().unwrap_or(0);
        writeln!(w, "#   {:<4}: {}", elem, count)?;
    }
    writeln!(w, "# {}", "-".repeat(60))?;
    Ok(())
}

/// Write g(r) data to a tab-separated text file (`.gr`).
///
/// Columns: `r[Ang]`, then symmetric pair columns in periodic-table order
/// (`El1-El2` with Z(El1) ≤ Z(El2)), and finally `total`.
pub fn write_gr(result: &GrResult, path: &str) -> std::io::Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    write_header(&mut w, result, "Radial Distribution Function g(r)")?;

    let keys = sorted_keys(&result.gr);
    write!(w, "# r[Ang]")?;
    for k in &keys { write!(w, "\t{}", k)?; }
    writeln!(w)?;

    let n = result.r.len();
    for i in 0..n {
        write!(w, "{:.6e}", result.r[i])?;
        for k in &keys {
            let v = result.gr.get(k).map(|v| v[i]).unwrap_or(0.0);
            write!(w, "\t{:.6e}", v)?;
        }
        writeln!(w)?;
    }
    Ok(())
}

/// Write CN(r) data to a tab-separated text file (`.cn`).
///
/// Columns: `r[Ang]`, then **directed** pair columns in (Z_center, Z_neighbor) order.
/// Key `"A-B"` = cumulative average number of B atoms within r from each A atom.
/// For A ≠ B, both `"A-B"` and `"B-A"` columns are present.
pub fn write_cn(result: &GrResult, path: &str) -> std::io::Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    write_header(&mut w, result, "Coordination Number CN(r) [directed: center-neighbor]")?;

    let keys = sorted_keys(&result.cn);
    write!(w, "# r[Ang]")?;
    for k in &keys { write!(w, "\t{}", k)?; }
    writeln!(w)?;

    let n = result.r.len();
    for i in 0..n {
        write!(w, "{:.6e}", result.r[i])?;
        for k in &keys {
            let v = result.cn.get(k).map(|v| v[i]).unwrap_or(0.0);
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
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    /// 构造简单立方 Fe 超胞（n×n×n，a=2.87 Å）
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

    /// 构造两组分系统（O 和 Si 交替，a=3.0 Å）
    fn make_o_si_crystal(n: usize) -> Frame {
        let a = 3.0_f64;
        let side = n as f64 * a;
        let cell = Cell::from_lengths_angles(side, side, side, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let pos = Vector3::new(i as f64 * a, j as f64 * a, k as f64 * a);
                    let elem = if (i + j + k) % 2 == 0 { "O" } else { "Si" };
                    frame.add_atom(Atom::new(elem, pos));
                }
            }
        }
        frame
    }

    #[test]
    fn test_elements_sorted_by_z() {
        // Si(Z=14) 出现在第一帧第一个，O(Z=8) 在第二个
        // 排序后应为 [O, Si]
        let cell = Cell::from_lengths_angles(9.0, 9.0, 9.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Si", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("O",  Vector3::new(1.5, 0.0, 0.0)));
        frame.add_atom(Atom::new("O",  Vector3::new(3.0, 0.0, 0.0)));
        let traj = Trajectory::from_frame(frame);
        let params = GrParams { r_min: 0.1, r_max: 3.5, dr: 0.01, r_cut: 2.0 };
        let res = calc_gr(&traj, &params).unwrap();
        assert_eq!(res.elements, vec!["O".to_string(), "Si".to_string()]);
    }

    #[test]
    fn test_gr_symmetric_pair_keys_canonical_order() {
        // g(r) key 应为 "O-Si"（O 原子序数 8 < Si 14），不是 "Si-O"
        let traj = Trajectory::from_frame(make_o_si_crystal(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.01, r_cut: 2.0 };
        let res = calc_gr(&traj, &params).unwrap();
        assert!(res.gr.contains_key("O-Si"), "expected O-Si key in gr");
        assert!(!res.gr.contains_key("Si-O"), "Si-O should NOT be in gr (asymmetric)");
    }

    #[test]
    fn test_cn_directed_both_directions() {
        // CN 文件中 O-Si 和 Si-O 均应存在，且值不同（N_O ≠ N_Si）
        let traj = Trajectory::from_frame(make_o_si_crystal(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.01, r_cut: 2.0 };
        let res = calc_gr(&traj, &params).unwrap();
        assert!(res.cn.contains_key("O-Si"), "O-Si missing from cn");
        assert!(res.cn.contains_key("Si-O"), "Si-O missing from cn");

        // 若 N_O = N_Si，CN 相等；若不等则不同——此晶体 O 和 Si 各占一半，N_O = N_Si
        // 所以相等，但键名应该都存在
    }

    #[test]
    fn test_cn_directed_unequal_when_n_differs() {
        // P(Z=15) 和 O(Z=8) 数量不同时，CN("O-P") ≠ CN("P-O")
        let cell = Cell::from_lengths_angles(15.0, 15.0, 15.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        // 1 个 P + 4 个 O
        frame.add_atom(Atom::new("P", Vector3::new(0.0, 0.0, 0.0)));
        for i in 1..=4 {
            frame.add_atom(Atom::new("O", Vector3::new(i as f64 * 3.0, 0.0, 0.0)));
        }
        let traj = Trajectory::from_frame(frame);
        let params = GrParams { r_min: 0.1, r_max: 7.0, dr: 0.05, r_cut: 13.0 };
        let res = calc_gr(&traj, &params).unwrap();

        let cn_op = res.cn["O-P"].last().copied().unwrap_or(0.0); // CN of P around each O
        let cn_po = res.cn["P-O"].last().copied().unwrap_or(0.0); // CN of O around each P
        // N_P=1, N_O=4 → CN_PO/CN_OP ≈ N_O/N_P = 4
        assert!(
            cn_po > cn_op * 2.0,
            "CN(P→O)={:.2} should be >> CN(O→P)={:.2} since N_O/N_P=4",
            cn_po, cn_op
        );
    }

    #[test]
    fn test_gr_first_peak_sc_fe() {
        let traj = Trajectory::from_frame(make_sc_fe(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.01, r_cut: 3.0 };
        let res = calc_gr(&traj, &params).unwrap();

        let fe_gr = &res.gr["Fe-Fe"];
        let (peak_bin, _) = fe_gr.iter().enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();
        let peak_r = res.r[peak_bin];
        assert!(
            (peak_r - 2.87).abs() < 0.02,
            "first peak at {:.3} Å, expected ~2.87 Å",
            peak_r
        );
    }

    #[test]
    fn test_cn_sc_fe_first_shell() {
        let traj = Trajectory::from_frame(make_sc_fe(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.01, r_cut: 3.0 };
        let res = calc_gr(&traj, &params).unwrap();
        let cn_last = *res.cn["Fe-Fe"].last().unwrap();
        assert!((cn_last - 6.0).abs() < 0.3, "CN(3.9) = {:.2}, expected ~6", cn_last);
    }

    #[test]
    fn test_sorted_keys_order() {
        // sorted_keys 应按 (Z_left, Z_right) 排，"total" 最后
        let traj = Trajectory::from_frame(make_o_si_crystal(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.01, r_cut: 2.0 };
        let res = calc_gr(&traj, &params).unwrap();

        let gr_keys = sorted_keys(&res.gr);
        // 期望顺序：O-O, O-Si, Si-Si, total
        let non_total: Vec<_> = gr_keys.iter().filter(|k| k.as_str() != "total").collect();
        assert_eq!(non_total[0].as_str(), "O-O");
        assert_eq!(non_total[1].as_str(), "O-Si");
        assert_eq!(non_total[2].as_str(), "Si-Si");
        assert_eq!(gr_keys.last().unwrap().as_str(), "total");
    }

    #[test]
    fn test_write_gr_and_cn() {
        use std::io::Read;
        let traj = Trajectory::from_frame(make_sc_fe(2));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.1, r_cut: 3.0 };
        let res = calc_gr(&traj, &params).unwrap();

        let gr_path = "/tmp/test_ferro.gr";
        let cn_path = "/tmp/test_ferro.cn";
        write_gr(&res, gr_path).expect("write_gr failed");
        write_cn(&res, cn_path).expect("write_cn failed");

        // 验证文件存在且有 # 开头的行
        let mut gr_content = String::new();
        std::fs::File::open(gr_path).unwrap().read_to_string(&mut gr_content).unwrap();
        assert!(gr_content.starts_with("# ferro v"));
        assert!(gr_content.contains("# r[Ang]"));

        let mut cn_content = String::new();
        std::fs::File::open(cn_path).unwrap().read_to_string(&mut cn_content).unwrap();
        assert!(cn_content.contains("directed"));
    }

    #[test]
    fn test_elem_z_pseudo_elements() {
        // elem_z 应能正确识别伪元素标签并返回基础元素的原子序数
        assert_eq!(elem_z("O"),   8,  "exact match");
        assert_eq!(elem_z("Ob"),  8,  "bridging oxygen");
        assert_eq!(elem_z("Onb"), 8,  "non-bridging oxygen");
        assert_eq!(elem_z("P"),   15, "exact match");
        assert_eq!(elem_z("P0"),  15, "P with 0 BO");
        assert_eq!(elem_z("P1"),  15, "P with 1 BO");
        assert_eq!(elem_z("P3"),  15, "P with 3 BO");
        assert_eq!(elem_z("Fe"),  26, "exact match");
        assert_eq!(elem_z("Fe1"), 26, "Fe site 1");
        assert_eq!(elem_z("Fe2"), 26, "Fe site 2");
        assert_eq!(elem_z("Na1"), 11, "Na site 1");
        assert_eq!(elem_z("??"),  255, "totally unknown");
    }

    #[test]
    fn test_pseudo_element_gr_key_ordering() {
        // 伪元素轨迹：Ob(Z=8) 和 P1(Z=15)，gr key 应为 "Ob-P1"（低 Z 在前）
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("P1", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Ob", Vector3::new(1.6, 0.0, 0.0)));
        frame.add_atom(Atom::new("Onb", Vector3::new(0.0, 1.6, 0.0)));
        let traj = Trajectory::from_frame(frame);
        let params = GrParams { r_min: 0.1, r_max: 4.5, dr: 0.05, r_cut: 2.0 };
        let res = calc_gr(&traj, &params).unwrap();

        // 应有 "Ob-P1" 和 "Onb-P1" 键（Z(O)=8 < Z(P)=15，低 Z 在前）
        assert!(res.gr.contains_key("Ob-P1"),  "Ob-P1 missing");
        assert!(res.gr.contains_key("Onb-P1"), "Onb-P1 missing");
        assert!(!res.gr.contains_key("P1-Ob"),  "P1-Ob should not exist (non-canonical)");

        // sorted_keys 应使 Ob 系列在 Onb 系列之前（字典序 "Ob" < "Onb"）
        let keys = sorted_keys(&res.gr);
        let ob_pos  = keys.iter().position(|k| k == "Ob-P1").unwrap();
        let onb_pos = keys.iter().position(|k| k == "Onb-P1").unwrap();
        assert!(ob_pos < onb_pos, "Ob-P1 should sort before Onb-P1");
    }
}
