//! Radial distribution function g(r) calculation and output.
//!
//! Design reference:
//!   - Parameter style follows code1/gr.c (GrParams fields map to DR/RMAX/RMIN macros)
//!   - Multi-component partial g(r) and directed CN follow code2/dump2sq.c CalcCn_pp + CalcGr
//!
//! Column ordering: sorted by (atomic number, label string), i.e. periodic-table order
//! with a deterministic tie-break for site labels sharing an element.
//!
//! Two key semantics, both keyed by the **ordered** pair `"A-B"`:
//!   - `gr` — symmetric: `"A-B"` and `"B-A"` hold pointwise identical values
//!   - `cn` — directed: `"A-B"` is the average number of B around each A, so
//!     `"A-B"` and `"B-A"` generally differ

use ferro_core::Trajectory;
use ferro_core::error::ChemError;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::io::{BufWriter, Write};
use super::angle::CellList;

/// ferro package version (from Cargo.toml)
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

// ─── 内部辅助 ────────────────────────────────────────────────────────────────

/// Look up the atomic number for an element symbol or site label (e.g. "O_b_P_P", "P_0").
///
/// Matching strategy (tried in order):
///   1. Exact match ("Fe" → 26, "O" → 8).
///   2. First two bytes (uppercase + lowercase) as a chemical symbol ("Zn_f" → "Zn" → 30).
///   3. First byte (uppercase letter) as a single-character element ("O_b_P_P" → "O" → 8).
///   4. Unrecognised → 255 (sorted to end, with string secondary ordering).
pub(super) fn elem_z(symbol: &str) -> u8 {
    ferro_core::data::elements::symbol_to_z(symbol)
}

/// How atoms are grouped into the "types" that partial g(r) / CN are resolved over.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum GroupBy {
    /// Group by `Atom::element` — the chemical element (default).
    #[default]
    Element,
    /// Group by `Atom::label` (site type), falling back to `element` where no label is set.
    Label,
}

/// The grouping key of a single atom under the selected mode.
pub(super) fn group_key(atom: &ferro_core::Atom, by: GroupBy) -> &str {
    match by {
        GroupBy::Element => atom.element.as_str(),
        GroupBy::Label => atom.label.as_deref().unwrap_or(atom.element.as_str()),
    }
}

/// Collect the distinct grouping keys of a frame, ordered by (Z, string).
///
/// The string tie-break is what makes column order reproducible: site labels sharing an
/// element (`O_f` / `O_b_P_P`, `P_0` … `P_3`) all map to the same Z, so sorting by Z alone
/// would leave their relative order at the mercy of hash iteration order.
pub(super) fn sorted_types(frame: &ferro_core::Frame, by: GroupBy) -> Vec<String> {
    let mut set = std::collections::HashSet::new();
    for a in &frame.atoms {
        set.insert(group_key(a, by).to_string());
    }
    let mut types: Vec<String> = set.into_iter().collect();
    types.sort_by(|a, b| (elem_z(a), a.as_str()).cmp(&(elem_z(b), b.as_str())));
    types
}

/// Split a `"A-B"` key into `(A, B)`.
fn split_pair(key: &str) -> (&str, &str) {
    key.split_once('-').unwrap_or((key, ""))
}

/// Sort pair keys by (Z(left), left, Z(right), right).
///
/// Ordering by the left member first groups all columns sharing a centre type together.
pub fn sorted_keys(map: &BTreeMap<String, Vec<f64>>) -> Vec<String> {
    let mut keys: Vec<String> = map.keys().cloned().collect();
    keys.sort_by(|ka, kb| {
        let (a1, b1) = split_pair(ka);
        let (a2, b2) = split_pair(kb);
        (elem_z(a1), a1, elem_z(b1), b1)
            .cmp(&(elem_z(a2), a2, elem_z(b2), b2))
    });
    keys
}

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Parameters for g(r) calculation.
#[derive(Debug, Clone)]
pub struct GrParams {
    /// Minimum distance \[Å\] (default: 0.005)
    pub r_min: f64,
    /// Maximum distance \[Å\]. Clamped internally to half the smallest interplanar
    /// spacing across all frames — the minimum-image upper bound (default: 10.005)
    pub r_max: f64,
    /// Distance bin width \[Å\] (default: 0.01)
    pub dr: f64,
    /// Whether partials are resolved over elements or site labels (default: `Element`)
    pub group_by: GroupBy,
}

impl Default for GrParams {
    fn default() -> Self {
        GrParams { r_min: 0.005, r_max: 10.005, dr: 0.01, group_by: GroupBy::default() }
    }
}

impl GrParams {
    /// Determine `r_max` automatically as the minimum-image cutoff of the first frame.
    pub fn with_auto_rmax(traj: &Trajectory) -> Self {
        let r_max = traj
            .first()
            .and_then(|f| f.cell.as_ref())
            .and_then(|cell| cell.minimum_image_cutoff().ok())
            .unwrap_or(10.005);
        GrParams { r_max, ..Default::default() }
    }
}

// ─── 结果结构体 ──────────────────────────────────────────────────────────────

/// Result of a g(r) / CN(r) calculation.
///
/// Both maps are keyed by the **ordered** pair `"A-B"` and contain all n² pairs.
/// - `gr` — symmetric, so `gr["A-B"]` and `gr["B-A"]` are pointwise identical.
/// - `cn` — directed: `cn["A-B"]` is the cumulative average number of B atoms within r
///   of each A atom, so `cn["A-B"]` and `cn["B-A"]` differ whenever N_A ≠ N_B.
#[derive(Debug, Clone)]
pub struct GrResult {
    /// Bin-centre r values \[Å\]
    pub r: Vec<f64>,
    /// Partial g(r), key `"A-B"` (all n² ordered pairs)
    pub gr: BTreeMap<String, Vec<f64>>,
    /// Directed cumulative CN(r), key `"centre-neighbour"` (all n² ordered pairs)
    pub cn: BTreeMap<String, Vec<f64>>,
    /// Grouping keys present (elements or site labels), sorted by (Z, string)
    pub elements: Vec<String>,
    /// Atom count per grouping key (from the first frame)
    pub element_counts: BTreeMap<String, usize>,
    /// How `elements` was derived
    pub group_by: GroupBy,
    pub n_frames: usize,
    /// Average box volume \[Å³\]
    pub avg_volume: f64,
    /// Total number density N/V \[Å⁻³\]
    pub rho: f64,
    /// Parameters actually used (`r_max` reflects the applied clamp)
    pub params: GrParams,
}

// ─── 计算函数 ────────────────────────────────────────────────────────────────

/// Compute partial g(r) and directed CN(r) for every ordered pair of types.
///
/// Requires periodic cells in every frame; uses the minimum-image convention.
/// Returns `Err` when the trajectory is empty, bin count is zero, or any frame has no cell.
pub fn calc_gr(traj: &Trajectory, params: &GrParams) -> ferro_core::Result<GrResult> {
    if traj.frames.is_empty() {
        return Err(ChemError::ValidationError("trajectory is empty".into()));
    }

    // 最小镜像有效上界 = 所有帧最小面间距的一半。r_max 超过它时 CellList 每轴退化为
    // 单个 cell，最小镜像在该距离以上失效，会让 g(r) 尾部被错误压向 0。用面间距而非
    // 边长：非正交晶胞 d < L，用边长会高估可用截断。
    let mut mic_max = f64::INFINITY;
    for frame in &traj.frames {
        let cell = frame.cell.as_ref().ok_or_else(|| {
            ChemError::ValidationError("all frames must have a periodic cell".into())
        })?;
        mic_max = mic_max.min(cell.minimum_image_cutoff()?);
    }
    let r_max = if mic_max.is_finite() { params.r_max.min(mic_max) } else { params.r_max };

    let n_bins = ((r_max - params.r_min) / params.dr).floor() as usize;
    if n_bins == 0 {
        return Err(ChemError::ValidationError("r range produces zero bins".into()));
    }

    // ── 分组类型列表（元素或位点标签），按 (Z, 字符串) 升序 ──────────────
    let first_frame = traj.frames.first()
        .ok_or_else(|| ChemError::ValidationError("empty trajectory".into()))?;
    let by = params.group_by;
    let types = sorted_types(first_frame, by);
    let n_types = types.len();

    let type_idx: std::collections::HashMap<&str, usize> = types
        .iter().enumerate().map(|(i, e)| (e.as_str(), i)).collect();

    // ── 无序对列表 (ti ≤ tj)，直方图按无序对累加，输出时展开为 n² 个有序对 ──
    let mut upairs: Vec<(usize, usize)> = Vec::new();
    for i in 0..n_types {
        for j in i..n_types { upairs.push((i, j)); }
    }
    let n_upairs = upairs.len();

    // O(1) 无序对索引查表
    let mut pair_lookup = vec![vec![0usize; n_types]; n_types];
    for (pidx, &(i, j)) in upairs.iter().enumerate() {
        pair_lookup[i][j] = pidx;
        pair_lookup[j][i] = pidx;
    }

    // 并行逐帧计数：每个线程独立维护 (hist, volume)，最后 reduce 合并
    // 使用 linked-cell list 将帧内原子对遍历从 O(N²) 降至 O(N)
    let init = || (vec![vec![0.0f64; n_bins]; n_upairs], 0.0f64);
    let (hist, total_volume) = traj.frames.par_iter()
        .fold(init, |(mut h, mut v), frame| {
            let cell = frame.cell.as_ref().unwrap();
            v += cell.volume();
            let n_atoms = frame.atoms.len();
            let cl = CellList::build(frame, cell, r_max);
            for i in 0..n_atoms {
                let Some(&ti) = type_idx.get(group_key(&frame.atoms[i], by)) else { continue; };
                for (j, r_val, _) in cl.neighbors_of(i, r_max) {
                    if j <= i { continue; }
                    if r_val < params.r_min { continue; }
                    let bin = ((r_val - params.r_min) / params.dr) as usize;
                    if bin >= n_bins { continue; }
                    let Some(&tj) = type_idx.get(group_key(&frame.atoms[j], by)) else { continue; };
                    h[pair_lookup[ti][tj]][bin] += 1.0;
                }
            }
            (h, v)
        })
        .reduce(init, |(mut ha, va), (hb, vb)| {
            for p in 0..ha.len() {
                for b in 0..ha[p].len() { ha[p][b] += hb[p][b]; }
            }
            (ha, va + vb)
        });

    // ── normalisation ───────────────────────────────────────────────────────
    let n_frames = traj.frames.len();
    let avg_volume = total_volume / n_frames as f64;
    if avg_volume <= 0.0 {
        return Err(ChemError::ValidationError("average cell volume is zero or negative".into()));
    }

    // atom counts per type from first frame
    let mut type_counts = vec![0.0f64; n_types];
    let mut element_counts: BTreeMap<String, usize> = BTreeMap::new();
    for atom in &first_frame.atoms {
        if let Some(&ti) = type_idx.get(group_key(atom, by)) {
            type_counts[ti] += 1.0;
            *element_counts.entry(types[ti].clone()).or_insert(0) += 1;
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

    // ── 展开为 n² 个有序对 ───────────────────────────────────────────────
    // g(r) 归一化（code2 CalcGr）：
    //   同种 (A=A): g = 2·cr / (4πr²Δr·(N_A−1)/V·N_A·steps)
    //   异种 (A-B): g = cr  / (4πr²Δr·N_B/V·N_A·steps)      —— 对调 A,B 数值不变
    // 有向 CN：
    //   CN(A→B) = Σ_bins ni·count_AB / (N_A·steps)，同种 ni=2，异种 ni=1
    for ti in 0..n_types {
        for tj in 0..n_types {
            let n_a = type_counts[ti];
            let n_b = type_counts[tj];
            if n_a < 1.0 || n_b < 1.0 { continue; }

            let pidx = pair_lookup[ti][tj];
            let same_type = ti == tj;
            let (ni, n_neighbor) = if same_type { (2.0f64, n_a - 1.0) } else { (1.0f64, n_b) };
            let key = format!("{}-{}", types[ti], types[tj]);

            let mut gr_vec = vec![0.0f64; n_bins];
            let mut cn_vec = vec![0.0f64; n_bins];
            let mut running = 0.0f64;
            for (i, &r_c) in r_centers.iter().enumerate() {
                let bunbo = pi4 * r_c * r_c * params.dr * n_neighbor / avg_volume * n_a;
                if bunbo > 0.0 {
                    gr_vec[i] = ni * hist[pidx][i] / (bunbo * n_frames as f64);
                }
                running += ni * hist[pidx][i] / (n_a * n_frames as f64);
                cn_vec[i] = running;
            }
            gr_map.insert(key.clone(), gr_vec);
            cn_map.insert(key, cn_vec);
        }
    }

    Ok(GrResult {
        r: r_centers,
        gr: gr_map,
        cn: cn_map,
        elements: types,
        element_counts,
        group_by: by,
        n_frames,
        avg_volume,
        rho,
        params: GrParams { r_max, ..params.clone() },
    })
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

/// Write header lines shared by the `.dat` output.
fn write_header(w: &mut impl Write, result: &GrResult) -> std::io::Result<()> {
    let what = match result.group_by {
        GroupBy::Element => "element",
        GroupBy::Label => "label",
    };
    writeln!(w, "# ferro v{}", VERSION)?;
    writeln!(w, "# Radial Distribution Function g(r) and Coordination Number CN(r)")?;
    writeln!(w, "# {}", "-".repeat(60))?;
    writeln!(w, "# r_min   = {} Ang", result.params.r_min)?;
    writeln!(w, "# r_max   = {} Ang", result.params.r_max)?;
    writeln!(w, "# dr      = {} Ang", result.params.dr)?;
    writeln!(w, "# frames  = {}", result.n_frames)?;
    writeln!(w, "# volume  = {:.3} Ang^3", result.avg_volume)?;
    writeln!(w, "# density = {:.6e} Ang^-3", result.rho)?;
    writeln!(w, "# grouped by {what}:")?;
    for elem in &result.elements {
        let count = result.element_counts.get(elem).copied().unwrap_or(0);
        writeln!(w, "#   {:<10}: {}", elem, count)?;
    }
    writeln!(w, "# pairs: <A>-<B>_gr is symmetric, <A>-<B>_cn is directed (centre=A, neighbour=B)")?;
    writeln!(w, "# {}", "-".repeat(60))?;
    Ok(())
}

/// Write g(r) and CN(r) to a single tab-separated text file.
///
/// - `pair = Some((a, b))` → three columns: `r[Ang]`, `a-b_gr`, `a-b_cn`.
/// - `pair = None` → wide table: `r[Ang]` then every ordered pair as an adjacent
///   `_gr` / `_cn` column duplet, ordered by (Z_centre, centre, Z_neighbour, neighbour)
///   so all columns sharing a centre stay together.
pub fn write_gr(
    result: &GrResult,
    path: &str,
    pair: Option<(&str, &str)>,
) -> std::io::Result<()> {
    let keys: Vec<String> = match pair {
        Some((a, b)) => {
            let key = format!("{a}-{b}");
            if !result.gr.contains_key(&key) {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    format!("pair '{key}' not present in the trajectory"),
                ));
            }
            vec![key]
        }
        None => sorted_keys(&result.gr),
    };

    let mut w = BufWriter::new(std::fs::File::create(path)?);
    write_header(&mut w, result)?;

    write!(w, "# r[Ang]")?;
    for k in &keys { write!(w, "\t{k}_gr\t{k}_cn")?; }
    writeln!(w)?;

    for i in 0..result.r.len() {
        write!(w, "{:.6e}", result.r[i])?;
        for k in &keys {
            let g = result.gr.get(k).map(|v| v[i]).unwrap_or(0.0);
            let c = result.cn.get(k).map(|v| v[i]).unwrap_or(0.0);
            write!(w, "\t{g:.6e}\t{c:.6e}")?;
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
        // Si(Z=14) 出现在第一帧第一个，O(Z=8) 在第二个；排序后应为 [O, Si]
        let cell = Cell::from_lengths_angles(9.0, 9.0, 9.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Si", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("O",  Vector3::new(1.5, 0.0, 0.0)));
        frame.add_atom(Atom::new("O",  Vector3::new(3.0, 0.0, 0.0)));
        let traj = Trajectory::from_frame(frame);
        let params = GrParams { r_min: 0.1, r_max: 3.5, dr: 0.01, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();
        assert_eq!(res.elements, vec!["O".to_string(), "Si".to_string()]);
    }

    #[test]
    fn test_all_ordered_pairs_present() {
        // n 个类型 → n² 个有序对，gr 与 cn 键集合完全一致
        let traj = Trajectory::from_frame(make_o_si_crystal(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.01, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();

        let n = res.elements.len();
        assert_eq!(res.gr.len(), n * n, "gr should hold n² ordered pairs");
        assert_eq!(res.cn.len(), n * n, "cn should hold n² ordered pairs");
        for a in &res.elements {
            for b in &res.elements {
                let k = format!("{a}-{b}");
                assert!(res.gr.contains_key(&k), "gr missing {k}");
                assert!(res.cn.contains_key(&k), "cn missing {k}");
            }
        }
        assert!(!res.gr.contains_key("total"), "total must not be emitted");
    }

    #[test]
    fn test_gr_is_symmetric_across_mirror_keys() {
        // g(r) 无方向：A-B 与 B-A 应逐点相同
        let traj = Trajectory::from_frame(make_o_si_crystal(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.01, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();
        for (x, y) in res.gr["O-Si"].iter().zip(res.gr["Si-O"].iter()) {
            assert!((x - y).abs() < 1e-12, "g(O-Si) != g(Si-O): {x} vs {y}");
        }
    }

    #[test]
    fn test_cn_directed_ratio_matches_count_ratio() {
        // P(Z=15) 1 个、O(Z=8) 4 个 → CN(P→O)/CN(O→P) == N_O/N_P == 4
        let cell = Cell::from_lengths_angles(15.0, 15.0, 15.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("P", Vector3::new(0.0, 0.0, 0.0)));
        for i in 1..=4 {
            frame.add_atom(Atom::new("O", Vector3::new(i as f64 * 3.0, 0.0, 0.0)));
        }
        let traj = Trajectory::from_frame(frame);
        let params = GrParams { r_min: 0.1, r_max: 7.0, dr: 0.05, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();

        let cn_op = *res.cn["O-P"].last().unwrap(); // 每个 O 周围的 P 数
        let cn_po = *res.cn["P-O"].last().unwrap(); // 每个 P 周围的 O 数
        assert!(cn_op > 0.0, "CN(O→P) should be non-zero");
        assert!(
            (cn_po / cn_op - 4.0).abs() < 1e-9,
            "CN(P→O)/CN(O→P) = {:.4}, expected N_O/N_P = 4",
            cn_po / cn_op
        );
    }

    #[test]
    fn test_gr_first_peak_sc_fe() {
        let traj = Trajectory::from_frame(make_sc_fe(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.01, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();

        let fe_gr = &res.gr["Fe-Fe"];
        let (peak_bin, _) = fe_gr.iter().enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();
        let peak_r = res.r[peak_bin];
        assert!(
            (peak_r - 2.87).abs() < 0.02,
            "first peak at {peak_r:.3} Å, expected ~2.87 Å"
        );
    }

    #[test]
    fn test_cn_sc_fe_first_shell() {
        let traj = Trajectory::from_frame(make_sc_fe(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.01, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();
        let cn_last = *res.cn["Fe-Fe"].last().unwrap();
        assert!((cn_last - 6.0).abs() < 0.3, "CN(3.9) = {cn_last:.2}, expected ~6");
    }

    #[test]
    fn test_sorted_keys_group_by_centre() {
        // 列顺序按 (Z_centre, centre, Z_neighbour, neighbour)，同一中心的列连在一起
        let traj = Trajectory::from_frame(make_o_si_crystal(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.01, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();
        assert_eq!(
            sorted_keys(&res.gr),
            vec!["O-O", "O-Si", "Si-O", "Si-Si"]
        );
    }

    #[test]
    fn test_write_gr_pair_three_columns() {
        use std::io::Read;
        let traj = Trajectory::from_frame(make_o_si_crystal(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.1, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();

        let path = "/tmp/test_ferro_pair.dat";
        write_gr(&res, path, Some(("Si", "O"))).expect("write_gr failed");

        let mut s = String::new();
        std::fs::File::open(path).unwrap().read_to_string(&mut s).unwrap();
        assert!(s.starts_with("# ferro v"));
        assert!(s.contains("# r[Ang]\tSi-O_gr\tSi-O_cn\n"), "header columns wrong");
        assert!(!s.contains("r_cut"), "r_cut should be gone from the header");
        assert!(!s.contains("pair stats"), "pair stats block should be gone");
        assert!(!s.contains("total"), "no total column or row");
        // 数据行恰好 3 列
        let first_data = s.lines().find(|l| !l.starts_with('#')).unwrap();
        assert_eq!(first_data.split('\t').count(), 3);
    }

    #[test]
    fn test_write_gr_wide_table_column_count() {
        use std::io::Read;
        let traj = Trajectory::from_frame(make_o_si_crystal(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.1, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();

        let path = "/tmp/test_ferro_wide.dat";
        write_gr(&res, path, None).expect("write_gr failed");

        let mut s = String::new();
        std::fs::File::open(path).unwrap().read_to_string(&mut s).unwrap();
        // 2 个类型 → 4 个有序对 → r + 4×2 = 9 列
        let first_data = s.lines().find(|l| !l.starts_with('#')).unwrap();
        assert_eq!(first_data.split('\t').count(), 9);
        assert!(s.contains("O-O_gr\tO-O_cn\tO-Si_gr\tO-Si_cn\tSi-O_gr"));
    }

    #[test]
    fn test_write_gr_unknown_pair_errors() {
        let traj = Trajectory::from_frame(make_sc_fe(2));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.1, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();
        let err = write_gr(&res, "/tmp/test_ferro_missing.dat", Some(("Zn", "O")));
        assert!(err.is_err(), "unknown pair should be rejected, not silently zero-filled");
    }

    #[test]
    fn test_elem_z_site_labels() {
        assert_eq!(elem_z("O"), 8, "exact match");
        assert_eq!(elem_z("O_f"), 8, "free oxygen label");
        assert_eq!(elem_z("O_b_P_P"), 8, "bridging oxygen label");
        assert_eq!(elem_z("P"), 15, "exact match");
        assert_eq!(elem_z("P_0"), 15, "P with 0 BO");
        assert_eq!(elem_z("Zn_f"), 30, "two-letter symbol prefix");
        assert_eq!(elem_z("??"), 255, "totally unknown");
    }

    #[test]
    fn test_gr_rmax_clamped_to_minimum_image_cutoff() {
        // 盒子 6 Å → MIC 上界 = 3.0。请求 r_max=10 时若不 clamp，CellList 每轴退化为
        // 1 个 cell，最小镜像在 r > 3 失效，g(r) 尾部会被错误压向 0。
        let cell = Cell::from_lengths_angles(6.0, 6.0, 6.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(3.0, 0.0, 0.0)));
        let traj = Trajectory::from_frame(frame);
        let params = GrParams { r_min: 0.1, r_max: 10.0, dr: 0.01, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();

        assert!(res.r.last().copied().unwrap() <= 3.0 + 1e-9);
        assert!(res.params.r_max <= 3.0 + 1e-9);
    }

    #[test]
    fn test_gr_rmax_clamp_uses_spacing_not_length_for_skewed_cell() {
        // 三斜晶胞：面间距 < 边长，clamp 必须用面间距。
        // 用边长会得到 min(L)/2 = 5.0，用面间距则明显更小。
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 60.0, 70.0, 80.0).unwrap();
        let mic = cell.minimum_image_cutoff().unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(2.5, 0.0, 0.0)));
        let traj = Trajectory::from_frame(frame);
        let params = GrParams { r_min: 0.1, r_max: 10.0, dr: 0.01, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();

        assert!(mic < 5.0, "sanity: skewed cell spacing/2 = {mic:.4} should be < min(L)/2 = 5");
        assert!(
            (res.params.r_max - mic).abs() < 1e-9,
            "r_max clamped to {:.6}, expected interplanar bound {mic:.6}",
            res.params.r_max
        );
    }

    #[test]
    fn test_gr_rmax_not_clamped_when_within_bound() {
        let traj = Trajectory::from_frame(make_sc_fe(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.01, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();
        assert!((res.params.r_max - 3.9).abs() < 1e-12);
    }

    #[test]
    fn test_group_by_label_resolves_sites() {
        // 同一元素的两个位点标签应分成独立类型；未设 label 的原子回退到 element
        let cell = Cell::from_lengths_angles(12.0, 12.0, 12.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        let mut push = |elem: &str, label: Option<&str>, x: f64| {
            let mut a = Atom::new(elem, Vector3::new(x, 0.0, 0.0));
            a.label = label.map(|s| s.to_string());
            frame.add_atom(a);
        };
        push("P", Some("P_0"), 0.0);
        push("O", Some("O_f"), 1.6);
        push("O", Some("O_b_P_P"), 3.2);
        push("Zn", None, 4.8);
        let traj = Trajectory::from_frame(frame);

        let by_elem = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 5.9, dr: 0.05, group_by: GroupBy::Element,
        }).unwrap();
        assert_eq!(by_elem.elements, vec!["O", "P", "Zn"]);

        let by_label = calc_gr(&traj, &GrParams {
            r_min: 0.1, r_max: 5.9, dr: 0.05, group_by: GroupBy::Label,
        }).unwrap();
        // (Z, 字符串) 排序：O(8) 的两个标签在前（"O_b_P_P" < "O_f"），P(15)，Zn(30)
        assert_eq!(by_label.elements, vec!["O_b_P_P", "O_f", "P_0", "Zn"]);
        assert!(by_label.gr.contains_key("P_0-O_b_P_P"));
        assert_eq!(by_label.gr.len(), 16);
    }

    #[test]
    fn test_same_z_label_order_is_deterministic() {
        // 同 Z 的多个标签必须每次给出相同顺序（否则下游按列位置取数会错位）
        let cell = Cell::from_lengths_angles(20.0, 20.0, 20.0, 90.0, 90.0, 90.0).unwrap();
        let labels = ["O_f", "O_n_P", "O_b_P_P", "O_x"];
        let mut frame = Frame::with_cell(cell, [true; 3]);
        for (i, l) in labels.iter().enumerate() {
            let mut a = Atom::new("O", Vector3::new(i as f64 * 2.0, 0.0, 0.0));
            a.label = Some(l.to_string());
            frame.add_atom(a);
        }
        let traj = Trajectory::from_frame(frame);
        let params = GrParams {
            r_min: 0.1, r_max: 9.5, dr: 0.05, group_by: GroupBy::Label,
        };
        let expected = vec!["O_b_P_P", "O_f", "O_n_P", "O_x"];
        for _ in 0..8 {
            let res = calc_gr(&traj, &params).unwrap();
            assert_eq!(res.elements, expected, "label ordering must be reproducible");
        }
    }
}
