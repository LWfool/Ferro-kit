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

use ferro_core::{Table, Trajectory};
use ferro_core::error::ChemError;
use rayon::prelude::*;
use std::collections::BTreeMap;
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
    /// Minimum distance \[Å\] (default: 0.001)
    pub r_min: f64,
    /// Maximum distance \[Å\]. Clamped internally to half the smallest interplanar
    /// spacing across all frames — the minimum-image upper bound (default: 10.005)
    pub r_max: f64,
    /// Distance bin width \[Å\] (default: 0.002)
    pub dr: f64,
    /// Whether partials are resolved over elements or site labels (default: `Element`)
    pub group_by: GroupBy,
}

impl Default for GrParams {
    fn default() -> Self {
        GrParams { r_min: 0.001, r_max: 10.005, dr: 0.002, group_by: GroupBy::default() }
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
    /// Per-frame average of `ρ_f · g_f(r)`, key `"A-B"` (all n² ordered pairs).
    ///
    /// This is the quantity S(q) actually needs, and it is *not* `rho * gr` under NPT:
    /// `g_f ∝ V_f` while `ρ_f ∝ 1/V_f`, so the box volume cancels out of the product and
    /// `⟨ρ_f·g_f⟩` reduces to a plain (volume-free) pair count. Keeping it separate is what
    /// lets `calc_sq_from_gr` reproduce a strictly per-frame transform from time-averaged
    /// data alone. Under NVT it degenerates to exactly `rho * gr`.
    pub rho_g: BTreeMap<String, Vec<f64>>,
    /// Grouping keys present (elements or site labels), sorted by (Z, string)
    pub elements: Vec<String>,
    /// Atom count per grouping key (from the first frame)
    pub element_counts: BTreeMap<String, usize>,
    /// How `elements` was derived
    pub group_by: GroupBy,
    pub n_frames: usize,
    /// Mean box volume ⟨V⟩ \[Å³\] — reported only; normalisation is per-frame
    pub avg_volume: f64,
    /// Population standard deviation of the box volume \[Å³\] (0 for NVT)
    pub volume_std: f64,
    /// Mean number density ⟨ρ_f⟩ = N·⟨1/V⟩ \[Å⁻³\].
    ///
    /// Note this is the average of the per-frame densities, **not** `N/⟨V⟩` — the two
    /// differ by an O(σ_V²/⟨V⟩²) term whenever the box fluctuates.
    pub rho: f64,
    /// Parameters actually used (`r_max` reflects the applied clamp)
    pub params: GrParams,
}

// ─── 计算函数 ────────────────────────────────────────────────────────────────

/// Per-thread fold accumulator for [`calc_gr`].
struct Acc {
    /// `Σ_f count` — volume-free counts, feeding CN(r) and `⟨ρ_f·g_f⟩`
    plain: Vec<Vec<f64>>,
    /// `Σ_f count · V_f` — volume-weighted counts, feeding g(r)
    vol_w: Vec<Vec<f64>>,
    /// Scratch buffer holding one frame's counts; reused across frames
    cr: Vec<Vec<f64>>,
}

impl Acc {
    fn new(n_upairs: usize, n_bins: usize) -> Self {
        Acc {
            plain: vec![vec![0.0f64; n_bins]; n_upairs],
            vol_w: vec![vec![0.0f64; n_bins]; n_upairs],
            cr:    vec![vec![0.0f64; n_bins]; n_upairs],
        }
    }

    fn merge(mut a: Self, b: Self) -> Self {
        for p in 0..a.plain.len() {
            for i in 0..a.plain[p].len() {
                a.plain[p][i] += b.plain[p][i];
                a.vol_w[p][i] += b.vol_w[p][i];
            }
        }
        a
    }
}

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
    let mut volumes: Vec<f64> = Vec::with_capacity(traj.frames.len());
    for frame in &traj.frames {
        let cell = frame.cell.as_ref().ok_or_else(|| {
            ChemError::ValidationError("all frames must have a periodic cell".into())
        })?;
        mic_max = mic_max.min(cell.minimum_image_cutoff()?);
        volumes.push(cell.volume());
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

    // ── 粒子数守恒校验 ───────────────────────────────────────────────────────
    // 归一化把 N_A / N_B / N 当作常数提出求和号，这只有在粒子数逐帧不变时成立。
    // 不校验的话，把两个不同体系拼成的轨迹喂进来会静默算出垃圾：下面主循环的
    // `else { continue }` 会把首帧没见过的类型悄悄丢掉，结果看着正常却完全错误。
    let count_types = |frame: &ferro_core::Frame| -> Vec<usize> {
        let mut c = vec![0usize; n_types];
        for atom in &frame.atoms {
            if let Some(&ti) = type_idx.get(group_key(atom, by)) { c[ti] += 1; }
        }
        c
    };
    let first_counts = count_types(first_frame);
    for (fi, frame) in traj.frames.iter().enumerate().skip(1) {
        if frame.atoms.len() != first_frame.atoms.len() {
            return Err(ChemError::ValidationError(format!(
                "atom count changes across the trajectory (frame 0 has {}, frame {} has {}); \
                 g(r) normalisation assumes a fixed particle number",
                first_frame.atoms.len(), fi, frame.atoms.len()
            )));
        }
        if count_types(frame) != first_counts {
            return Err(ChemError::ValidationError(format!(
                "per-type atom counts change at frame {fi}; \
                 g(r) normalisation assumes a fixed particle number per type"
            )));
        }
    }

    let type_counts: Vec<f64> = first_counts.iter().map(|&c| c as f64).collect();
    let element_counts: BTreeMap<String, usize> = types
        .iter().cloned().zip(first_counts.iter().copied())
        .filter(|&(_, c)| c > 0)
        .collect();
    let n_total = first_frame.atoms.len() as f64;

    // 并行逐帧计数：每个线程独立维护两份直方图与三个体积标量，最后 reduce 合并
    // 使用 linked-cell list 将帧内原子对遍历从 O(N²) 降至 O(N)
    //
    // 为什么是两份直方图：g(r) 要的是 ⟨hist·V_f⟩（逐帧归一化后再时间平均，见 code1/gr.c
    // 的「1ステップ毎にgrの計算」），而 CN(r) 与 ⟨ρ_f·g_f⟩ 要的是不带体积权重的纯计数。
    // NVT 下两者互为常数倍，NPT 下则相差一个 Cov(hist, V) 项。
    let init = || Acc::new(n_upairs, n_bins);
    let acc = traj.frames.par_iter()
        .fold(init, |mut acc, frame| {
            let cell = frame.cell.as_ref().unwrap();
            let vol = cell.volume();
            let n_atoms = frame.atoms.len();
            let cl = CellList::build(frame, cell, r_max);
            // 帧内纯计数缓冲（对应 code1 的 cr / code2 的 cn_），帧末一次性加权归并 ——
            // 比在每次命中处累加两份直方图便宜得多（命中数 ≫ bin 数）
            for row in acc.cr.iter_mut() { row.fill(0.0); }
            for i in 0..n_atoms {
                let Some(&ti) = type_idx.get(group_key(&frame.atoms[i], by)) else { continue; };
                for (j, r_val, _) in cl.neighbors_of(i, r_max) {
                    if j <= i { continue; }
                    if r_val < params.r_min { continue; }
                    let bin = ((r_val - params.r_min) / params.dr) as usize;
                    if bin >= n_bins { continue; }
                    let Some(&tj) = type_idx.get(group_key(&frame.atoms[j], by)) else { continue; };
                    acc.cr[pair_lookup[ti][tj]][bin] += 1.0;
                }
            }
            for p in 0..n_upairs {
                for b in 0..n_bins {
                    let c = acc.cr[p][b];
                    if c != 0.0 {
                        acc.plain[p][b] += c;
                        acc.vol_w[p][b] += c * vol;
                    }
                }
            }
            acc
        })
        .reduce(init, Acc::merge);

    // ── normalisation ───────────────────────────────────────────────────────
    let n_frames = traj.frames.len();
    let nf = n_frames as f64;
    let avg_volume = volumes.iter().sum::<f64>() / nf;
    if avg_volume <= 0.0 {
        return Err(ChemError::ValidationError("average cell volume is zero or negative".into()));
    }
    // 两遍算法：⟨V²⟩−⟨V⟩² 在 V≈6e4 时抵消误差会放大到 1e-3 量级，定容轨迹也报不出 0
    let volume_std = (volumes.iter().map(|v| (v - avg_volume).powi(2)).sum::<f64>() / nf).sqrt();
    // ⟨ρ_f⟩ = N·⟨1/V⟩ —— 逐帧密度的平均，不是 N/⟨V⟩
    let rho = n_total * volumes.iter().map(|v| 1.0 / v).sum::<f64>() / nf;

    let r_centers: Vec<f64> = (0..n_bins)
        .map(|i| params.r_min + (i as f64 + 0.5) * params.dr)
        .collect();

    let pi4 = 4.0 * std::f64::consts::PI;
    let mut gr_map: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut cn_map: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut rho_g_map: BTreeMap<String, Vec<f64>> = BTreeMap::new();

    // ── 展开为 n² 个有序对 ───────────────────────────────────────────────
    // 逐帧归一化后时间平均（code1/gr.c:139-146、code2/dump2sq.c CalcGr）：
    //   同种 (A=A): g_f = 2·cr_f·V_f / (4πr²Δr·N_A·(N_A−1))
    //   异种 (A-B): g_f =   cr_f·V_f / (4πr²Δr·N_A·N_B)     —— 对调 A,B 数值不变
    //   g = ⟨g_f⟩ —— NPT 下 ≠ ⟨cr⟩·⟨V⟩/…，两者相差一个 Cov(cr, V) 项
    // 有向 CN（与体积无关）：
    //   CN(A→B) = Σ_bins ni·⟨count_AB⟩ / N_A，同种 ni=2，异种 ni=1
    // S(q) 用的 ⟨ρ_f·g_f⟩：ρ_f = N/V_f 与 g_f ∝ V_f 相乘后体积消去，只剩纯计数
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
            let mut rho_g_vec = vec![0.0f64; n_bins];
            let mut running = 0.0f64;
            for (i, &r_c) in r_centers.iter().enumerate() {
                // 几何 + 计数因子，不含体积
                let bunbo = pi4 * r_c * r_c * params.dr * n_neighbor * n_a;
                if bunbo > 0.0 {
                    gr_vec[i] = ni * acc.vol_w[pidx][i] / (bunbo * nf);
                    rho_g_vec[i] = ni * n_total * acc.plain[pidx][i] / (bunbo * nf);
                }
                running += ni * acc.plain[pidx][i] / (n_a * nf);
                cn_vec[i] = running;
            }
            gr_map.insert(key.clone(), gr_vec);
            cn_map.insert(key.clone(), cn_vec);
            rho_g_map.insert(key, rho_g_vec);
        }
    }

    Ok(GrResult {
        r: r_centers,
        gr: gr_map,
        cn: cn_map,
        rho_g: rho_g_map,
        elements: types,
        element_counts,
        group_by: by,
        n_frames,
        avg_volume,
        volume_std,
        rho,
        params: GrParams { r_max, ..params.clone() },
    })
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

/// Write header lines shared by the `.dat` output.
impl GrResult {
    /// Projects the result into the long-format table the writers consume.
    ///
    /// One row per `(r, pair)`: columns `r, center, neighbor, gr, cn`.
    ///
    /// Long rather than wide because **different trajectories can have different
    /// element sets** — a Zn-P-O run and an Al-P-O run stack without any column
    /// alignment, whereas a pair-per-column layout would need a column union with
    /// holes. Splitting the pair into `center`/`neighbor` also writes CN's
    /// directedness into the table structure instead of leaving it to a doc comment.
    ///
    /// `pair = Some((a, b))` restricts the output to that ordered pair.
    /// The caller adds the `file` column when stacking several inputs
    /// (see `ferro_core::Table::concat_union`).
    pub fn to_tables(&self, pair: Option<(&str, &str)>) -> Result<Vec<(String, Table)>, ChemError> {
        let keys: Vec<String> = match pair {
            Some((a, b)) => {
                let key = format!("{a}-{b}");
                if !self.gr.contains_key(&key) {
                    return Err(ChemError::ValidationError(format!(
                        "pair '{key}' not present in the trajectory"
                    )));
                }
                vec![key]
            }
            None => sorted_keys(&self.gr),
        };

        let n = self.r.len();
        let rows = n * keys.len();
        let mut r_col        = Vec::with_capacity(rows);
        let mut centre_col   = Vec::with_capacity(rows);
        let mut neighbour_col = Vec::with_capacity(rows);
        let mut gr_col       = Vec::with_capacity(rows);
        let mut cn_col       = Vec::with_capacity(rows);

        for key in &keys {
            let (centre, neighbour) = split_pair(key);
            let g = self.gr.get(key);
            let c = self.cn.get(key);
            for i in 0..n {
                r_col.push(self.r[i]);
                centre_col.push(centre.to_string());
                neighbour_col.push(neighbour.to_string());
                gr_col.push(g.map(|v| v[i]).unwrap_or(f64::NAN));
                cn_col.push(c.map(|v| v[i]).unwrap_or(f64::NAN));
            }
        }

        let mut t = Table::new();
        t.push_num("r", r_col)
            .push_text("center", centre_col)
            .push_text("neighbor", neighbour_col)
            .push_num("gr", gr_col)
            .push_num("cn", cn_col);
        Ok(vec![("gr".to_string(), t)])
    }

    /// The analysis parameters, as `key = value` strings for the comment block.
    ///
    /// **Shared across a batch only.** Anything that varies per input — the composition,
    /// the clamped `r_max` — belongs in the per-input summary instead, where it cannot
    /// be mistaken for a global fact. Deliberately not machine-readable; anything a
    /// script must parse belongs in a column.
    pub fn meta_lines(&self) -> Vec<String> {
        let what = match self.group_by {
            GroupBy::Element => "element",
            GroupBy::Label => "label",
        };
        vec![
            format!("r_min   = {} Ang", self.params.r_min),
            format!("r_max   = {} Ang (requested; clamped per input, see [inputs])", self.params.r_max),
            format!("dr      = {} Ang", self.params.dr),
            format!("grouped by {what}"),
            "gr is symmetric; cn is directed (center -> neighbor)".to_string(),
        ]
    }

    /// Per-input composition, e.g. `O:1314 P:372 Zn:186` — a summary column, not a
    /// header line, because it differs from input to input.
    pub fn composition(&self) -> String {
        self.elements
            .iter()
            .map(|e| format!("{e}:{}", self.element_counts.get(e).copied().unwrap_or(0)))
            .collect::<Vec<_>>()
            .join(" ")
    }
}

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

    /// 构造合成 NPT 轨迹：同一 O/Si 晶体逐帧按 `scales` 缩放，晶胞与坐标同比变化，
    /// 故每帧都是合法晶体而体积逐帧不同。用于检验体积归一化口径。
    fn make_npt_traj(n: usize, scales: &[f64]) -> Trajectory {
        let base = make_o_si_crystal(n);
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

    /// 逐帧归一化后时间平均 —— 即 code1/code2 字面意义上的做法，作为对拍基准。
    ///
    /// 单帧调用 `calc_gr` 恰好就是「用该帧体积归一化」，所以对每帧各算一次再取
    /// 算术平均，就是参考实现的逐帧循环。
    fn literal_per_frame(traj: &Trajectory, gp: &GrParams) -> BTreeMap<String, Vec<f64>> {
        let mut acc: BTreeMap<String, Vec<f64>> = BTreeMap::new();
        for frame in &traj.frames {
            let single = Trajectory::from_frame(frame.clone());
            let g = calc_gr(&single, gp).unwrap();
            for (k, v) in &g.gr {
                let e = acc.entry(k.clone()).or_insert_with(|| vec![0.0; v.len()]);
                for (a, b) in e.iter_mut().zip(v.iter()) { *a += b; }
            }
        }
        let nf = traj.frames.len() as f64;
        for v in acc.values_mut() { for x in v.iter_mut() { *x /= nf; } }
        acc
    }

    #[test]
    fn test_gr_matches_literal_per_frame_normalisation() {
        // g(r) 必须等于「逐帧归一化再平均」，而不是「先平均计数再乘 ⟨V⟩」。
        // r_max 取得足够小，使各帧都不触发最小镜像 clamp，否则单帧与全轨迹的
        // bin 数会不一致而无法逐点比较。
        let traj = make_npt_traj(4, &[1.0, 1.08, 0.94, 1.12, 0.90]);
        let gp = GrParams { r_min: 0.1, r_max: 5.0, dr: 0.02, ..Default::default() };

        let folded = calc_gr(&traj, &gp).unwrap();
        let literal = literal_per_frame(&traj, &gp);

        assert_eq!(folded.gr.len(), literal.len());
        for (k, v) in &folded.gr {
            let l = &literal[k];
            let d = v.iter().zip(l.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            assert!(d < 1e-12, "pair {k}: folded vs literal per-frame differs by {d:.3e}");
        }
    }

    #[test]
    fn test_gr_differs_from_average_volume_normalisation_under_npt() {
        // 反向断言：若仍用 ⟨hist⟩·⟨V⟩ 归一化，结果与逐帧法**不同**。
        // 没有这条，上面的对拍测试在两种实现恰好等价时也会通过，失去意义。
        let traj = make_npt_traj(4, &[1.0, 1.08, 0.94, 1.12, 0.90]);
        let gp = GrParams { r_min: 0.1, r_max: 5.0, dr: 0.02, ..Default::default() };
        let folded = calc_gr(&traj, &gp).unwrap();

        // 用平均体积重算（旧口径）：g_old = g_new · ⟨V⟩·Σhist / Σ(hist·V)
        // 直接从 rho_g 反推 Σhist 的比例即可 —— rho_g ∝ Σhist、gr ∝ Σ(hist·V)
        let key = folded.gr.keys().next().unwrap().clone();
        let gr_v = &folded.gr[key.as_str()];
        let rg_v = &folded.rho_g[key.as_str()];
        // 若体积恒定则 gr/rho_g 处处为常数 1/rho；NPT 下该比值随 bin 漂移
        let ratios: Vec<f64> = gr_v.iter().zip(rg_v.iter())
            .filter(|(g, r)| **g > 1e-9 && **r > 1e-9)
            .map(|(g, r)| g / r)
            .collect();
        assert!(ratios.len() > 2, "need several populated bins for this check");
        let spread = ratios.iter().cloned().fold(f64::MIN, f64::max)
            / ratios.iter().cloned().fold(f64::MAX, f64::min) - 1.0;
        assert!(spread > 1e-6,
            "under NPT g(r) and <rho*g> must not stay proportional (spread {spread:.3e}); \
             if they do, the volume weighting was lost");
    }

    #[test]
    fn test_rho_g_degenerates_to_rho_times_gr_under_nvt() {
        // NVT（体积恒定）下 ⟨ρ_f·g_f⟩ 必须精确等于 ρ·g —— 这正是本次改动
        // 在定容轨迹上零影响的原因。
        let traj = make_npt_traj(4, &[1.0, 1.0, 1.0]);
        let gp = GrParams { r_min: 0.1, r_max: 5.0, dr: 0.02, ..Default::default() };
        let g = calc_gr(&traj, &gp).unwrap();
        assert!(g.volume_std / g.avg_volume < 1e-14,
            "constant-volume trajectory must report a negligible spread");
        for (k, v) in &g.gr {
            for (i, &gv) in v.iter().enumerate() {
                let expect = g.rho * gv;
                assert!((g.rho_g[k][i] - expect).abs() <= 1e-9 * (1.0 + expect.abs()),
                    "pair {k} bin {i}: rho_g {} != rho*gr {}", g.rho_g[k][i], expect);
            }
        }
    }

    #[test]
    fn test_volume_statistics_reported() {
        let traj = make_npt_traj(3, &[1.0, 2.0]);
        let gp = GrParams { r_min: 0.1, r_max: 4.0, dr: 0.05, ..Default::default() };
        let g = calc_gr(&traj, &gp).unwrap();
        let v1 = 9.0_f64.powi(3);
        let v2 = 18.0_f64.powi(3);
        assert!((g.avg_volume - (v1 + v2) / 2.0).abs() < 1e-6);
        assert!((g.volume_std - (v2 - v1) / 2.0).abs() < 1e-6);
        // rho = N·⟨1/V⟩，不是 N/⟨V⟩
        let n = 27.0_f64;
        assert!((g.rho - n * (1.0 / v1 + 1.0 / v2) / 2.0).abs() < 1e-12);
        assert!((g.rho - n / g.avg_volume).abs() > 1e-9, "rho must not be N/<V>");
    }

    #[test]
    fn test_changing_atom_count_is_rejected() {
        let f1 = make_o_si_crystal(3);
        let mut f2 = f1.clone();
        f2.add_atom(Atom::new("O", Vector3::new(1.5, 1.5, 1.5)));
        let mut traj = Trajectory::from_frame(f1);
        traj.add_frame(f2);
        let err = calc_gr(&traj, &GrParams::default()).unwrap_err();
        assert!(format!("{err}").contains("atom count changes"), "got: {err}");
    }

    #[test]
    fn test_changing_per_type_counts_is_rejected() {
        // 总原子数不变但组成变了 —— 这类轨迹拼接最容易被静默接受
        let f1 = make_o_si_crystal(3);
        let mut f2 = f1.clone();
        f2.atoms[0].element = "Si".to_string();
        let mut traj = Trajectory::from_frame(f1);
        traj.add_frame(f2);
        let err = calc_gr(&traj, &GrParams::default()).unwrap_err();
        assert!(format!("{err}").contains("per-type atom counts change"), "got: {err}");
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
    fn test_to_tables_pair_long_format() {
        let traj = Trajectory::from_frame(make_o_si_crystal(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.1, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();

        let tables = res.to_tables(Some(("Si", "O"))).unwrap();
        assert_eq!(tables.len(), 1);
        let (name, t) = &tables[0];
        assert_eq!(name, "gr");
        assert_eq!(t.names(), vec!["r", "center", "neighbor", "gr", "cn"]);
        assert_eq!(t.n_rows(), res.r.len(), "点名一对 → 行数 = r 网格长度");
        assert!(t.validate().is_ok());

        // center/neighbor 按调用方给的顺序写死,不按 Z 重排
        match t.column("center").unwrap() {
            ferro_core::Column::Text(v) => assert!(v.iter().all(|x| x == "Si")),
            _ => panic!("center must be text"),
        }
    }

    #[test]
    fn test_to_tables_all_pairs_stacks_rows_not_columns() {
        let traj = Trajectory::from_frame(make_o_si_crystal(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.1, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();

        let tables = res.to_tables(None).unwrap();
        let (_, t) = &tables[0];
        // 2 个类型 → 4 个有序对;长表下列数恒为 5,行数才随配对数增长
        assert_eq!(t.n_cols(), 5, "长表列数与体系组成无关");
        assert_eq!(t.n_rows(), 4 * res.r.len());
        assert!(t.validate().is_ok());
    }

    #[test]
    fn test_to_tables_unknown_pair_errors() {
        let traj = Trajectory::from_frame(make_sc_fe(2));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.1, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();
        assert!(
            res.to_tables(Some(("Zn", "O"))).is_err(),
            "unknown pair should be rejected, not silently zero-filled"
        );
    }

    #[test]
    fn test_meta_lines_report_clamped_rmax_and_composition() {
        let traj = Trajectory::from_frame(make_o_si_crystal(3));
        let params = GrParams { r_min: 0.1, r_max: 3.9, dr: 0.1, ..Default::default() };
        let res = calc_gr(&traj, &params).unwrap();
        let meta = res.meta_lines().join("\n");
        assert!(meta.contains("r_max"));
        assert!(meta.contains("clamped per input"), "r_max 是逐文件 clamp 的,头里要说清");
        assert!(!meta.contains("Si:"), "逐文件的组成不能混进共享参数块");
        assert!(res.composition().contains("Si:"), "组成走 composition()");
        assert!(!meta.contains("total"), "no total anywhere");
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
