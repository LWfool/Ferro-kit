//! Glass network trajectory statistics: bridging-ligand speciation (Qn for P),
//! ligand type distribution, coordination numbers.
//!
//! Atom-level classification is delegated to `ferro_core::classify_frame`.
//! This module accumulates per-frame [`AtomType`]s into time-averaged distributions.
//!
//! # What "fraction ± sd" means
//!
//! Every distribution row reports the **mean over frames of that frame's fraction**,
//! and `sd` is the sample standard deviation (ddof = 1) of the same per-frame
//! fraction.  A bin absent from a frame contributes 0 to both, not a gap.
//!
//! `sd` is a spread, **not a standard error**: consecutive MD frames are strongly
//! correlated, so it neither shrinks like `1/√N` nor estimates a physical
//! fluctuation.  Read it as "how much this number moved between snapshots".
//!
//! # Usage
//! ```ignore
//! use ferro_core::TypeParams;
//! use std::collections::BTreeMap;
//!
//! let mut cutoffs = BTreeMap::new();
//! cutoffs.insert(("P".into(), "O".into()), 2.3);
//! let params = TypeParams { cutoffs, modifier_cutoffs: BTreeMap::new() };
//! let result = calc_network(&traj, &params).unwrap();
//! ```

pub use ferro_core::{AtomType, CutoffTable, TypeParams};

use ferro_core::{classify_frame, Frame, Table, Trajectory};
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};

/// Ligand distribution key: `(sort rank, label, partner elements)`.
///
/// The rank comes from [`AtomType::class_rank`], so a `BTreeMap` keyed on this is
/// sorted by construction — no post-hoc `sort_by` that re-derives the order from the
/// label text.  Partners are kept only for the 0/1/2-former cases; a tricluster
/// collapses to one row per label, its pairs being reported by the linkage table.
type OxyKey = (u8, String, Vec<String>);

fn oxy_key(t: &AtomType) -> OxyKey {
    let partners = match t {
        AtomType::Ligand { partners, .. } if partners.len() <= 2 => partners.clone(),
        _ => Vec::new(),
    };
    (t.class_rank(), t.label(), partners)
}

// ─── 逐帧比例的一阶与二阶矩 ────────────────────────────────────────────────────

/// Pooled count plus a running mean/variance of the per-frame fraction.
///
/// Memory is O(number of bins), not O(bins × frames): no frame's histogram is kept.
/// Welford rather than `Σf` / `Σf²` — the textbook two-sum form computes the variance
/// of a constant series as the difference of two nearly equal numbers, which reports
/// `2e-10` where the answer is exactly `0`, and a spurious wobble in an `sd` column
/// is exactly the kind of thing a reader would take for physics.
///
/// Frames where the bin does not occur contribute `f = 0`.  They are **not** pushed
/// (the bin may not exist yet); [`Moments::finish`] folds them in once the total frame
/// count is known.
#[derive(Debug, Clone, Copy, Default)]
struct Moments {
    /// 逐帧计数之和（整数，精确）
    count: usize,
    /// 该 bin 出现过的帧数
    seen: usize,
    /// 逐帧比例的均值
    mean: f64,
    /// Σ(f − mean)²
    m2: f64,
}

impl Moments {
    fn push(&mut self, count: usize, frac: f64) {
        self.count += count;
        self.seen += 1;
        let delta = frac - self.mean;
        self.mean += delta / self.seen as f64;
        self.m2 += delta * (frac - self.mean);
    }

    /// Chan/Golub/LeVeque pairwise combination — needed because rayon reduces
    /// partial accumulators rather than folding one sequence.
    fn merge(&mut self, o: &Moments) {
        if o.seen == 0 { return; }
        if self.seen == 0 { *self = *o; return; }
        let (na, nb) = (self.seen as f64, o.seen as f64);
        let n = na + nb;
        let delta = o.mean - self.mean;
        self.count += o.count;
        self.m2 += o.m2 + delta * delta * na * nb / n;
        self.mean += delta * nb / n;
        self.seen += o.seen;
    }

    /// Folds in the frames where this bin did not occur (each contributing `f = 0`),
    /// then reports the mean and sample standard deviation over **all** frames.
    fn finish(&self, n_frames: usize) -> Bin {
        if n_frames == 0 {
            return Bin { count: self.count, fraction: f64::NAN, sd: f64::NAN };
        }
        let (mut mean, mut m2) = (self.mean, self.m2);
        let missing = n_frames.saturating_sub(self.seen);
        if missing > 0 {
            // 缺席组：均值 0、m2 0，套用同一条合并公式
            let (na, nb) = (self.seen as f64, missing as f64);
            let n = na + nb;
            let delta = -mean;
            m2 += delta * delta * na * nb / n;
            mean += delta * nb / n;
        }
        let sd = if n_frames < 2 {
            f64::NAN
        } else {
            (m2 / (n_frames as f64 - 1.0)).max(0.0).sqrt()
        };
        Bin { count: self.count, fraction: mean, sd }
    }
}

/// One row of a distribution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bin {
    /// Occurrences summed over every frame.
    pub count: usize,
    /// Mean over frames of the per-frame fraction.
    pub fraction: f64,
    /// Sample standard deviation (ddof = 1) of the per-frame fraction; `NaN` for a
    /// single frame.  A spread between snapshots, not a standard error.
    pub sd: f64,
}

// ─── 结果结构体 ───────────────────────────────────────────────────────────────

/// 网络分析时间平均结果。
#[derive(Debug, Clone)]
pub struct NetworkResult {
    /// former_elem → 桥接配体数分布（P 的即 Qn），按数值升序
    pub bridge_dist: HashMap<String, Vec<(u32, Bin)>>,
    /// former_elem → 平均桥接配体数
    pub mean_bridge: HashMap<String, f64>,
    /// 元素 → 总配位数分布（所有配体类型之和）。**含修饰子** ——
    /// 修饰子没有桥接数，配位数就是描述它的全部
    pub cn_dist: HashMap<String, Vec<(u32, Bin)>>,
    /// 元素 → 平均总配位数（含修饰子）
    pub mean_cn: HashMap<String, f64>,
    /// 配体类型分布：`(label, partner_elements, bin)`，按 `AtomType::class_rank` 排序
    pub oxy_dist: Vec<(String, Vec<String>, Bin)>,
    /// 参与统计的帧数
    pub n_frames: usize,
    /// 每帧原子数
    pub n_atoms: usize,
    /// 本次分析的参数（截断表 / 修饰子），供 `meta_lines` 使用
    pub params: TypeParams,
}

impl NetworkResult {
    /// Batch-wide parameters for the `#` comment block.
    ///
    /// Only what is shared by every input goes here — the cutoff table and the
    /// modifier list.  Per-input quantities (frames, atoms, composition) belong in
    /// the CLI's `[inputs]` summary, or they masquerade as global facts.
    pub fn meta_lines(&self) -> Vec<String> {
        let fmt = |t: &CutoffTable| -> String {
            t.iter()
                .map(|((a, b), c)| format!("{a}-{b}={c}"))
                .collect::<Vec<_>>()
                .join("  ")
        };
        let mut v = vec![format!("cutoffs   : {}", fmt(&self.params.cutoffs))];
        if !self.params.modifier_cutoffs.is_empty() {
            v.push(format!("modifiers : {}", fmt(&self.params.modifier_cutoffs)));
        }
        v.push("fraction  : mean over frames; sd = sample std (ddof=1) of the".to_string());
        v.push("            per-frame fraction — a spread, not a standard error".to_string());
        v
    }

    /// Long-format projections, one per product granularity.
    ///
    /// | table | one row per | columns |
    /// |---|---|---|
    /// | `bridge` | former × bridging count | `former, n_bridge, count, fraction, sd` |
    /// | `oxy` | ligand type × partner pair | `type, former_a, former_b, count, fraction, sd` |
    /// | `cn` | element × coordination number | `element, cn, count, fraction, sd` |
    /// | `mean` | element | `element, mean_n_bridge, mean_cn` |
    ///
    /// The three distributions have different value-column semantics, so they stay
    /// three tables — merging them would make `groupby` on the value column
    /// meaningless.  The means are one row per element rather than per value, which
    /// is a different granularity again, hence a fourth table instead of a column
    /// repeated down every row.
    ///
    /// `n_bridge` is deliberately not called `qn`: the count is defined for every
    /// former, but the Qn notation is a tetrahedral-former convention.  For P it is
    /// Qn; for Al, read the `cn` table instead.
    pub fn to_tables(&self) -> Vec<(String, Table)> {
        let mut out = Vec::with_capacity(4);

        let mut formers: Vec<&String> = self.bridge_dist.keys().collect();
        formers.sort();
        let mut elems: Vec<&String> = self.cn_dist.keys().collect();
        elems.sort();

        // bridge：只有形成子有桥接数
        let (mut f_col, mut v_col) = (Vec::new(), Vec::new());
        let mut bins = Vec::new();
        for former in &formers {
            for (v, b) in &self.bridge_dist[*former] {
                f_col.push((*former).clone());
                v_col.push(*v as f64);
                bins.push(*b);
            }
        }
        let mut t = Table::new();
        t.push_text("former", f_col).push_num("n_bridge", v_col);
        push_bins(&mut t, &bins);
        out.push(("bridge".to_string(), t));

        // oxy：伙伴元素以数据列给出,不编码进标签
        let (mut ty, mut fa, mut fb) = (Vec::new(), Vec::new(), Vec::new());
        let mut bins = Vec::new();
        for (label, partners, b) in &self.oxy_dist {
            ty.push(label.clone());
            fa.push(partners.first().cloned().unwrap_or_default());
            fb.push(partners.get(1).cloned().unwrap_or_default());
            bins.push(*b);
        }
        let mut t = Table::new();
        t.push_text("type", ty).push_text("former_a", fa).push_text("former_b", fb);
        push_bins(&mut t, &bins);
        out.push(("oxy".to_string(), t));

        // cn：形成子与修饰子共用一张表
        let (mut e_col, mut v_col) = (Vec::new(), Vec::new());
        let mut bins = Vec::new();
        for elem in &elems {
            for (v, b) in &self.cn_dist[*elem] {
                e_col.push((*elem).clone());
                v_col.push(*v as f64);
                bins.push(*b);
            }
        }
        let mut t = Table::new();
        t.push_text("element", e_col).push_num("cn", v_col);
        push_bins(&mut t, &bins);
        out.push(("cn".to_string(), t));

        // mean：一元素一行。修饰子没有桥接数,mean_n_bridge 留空(NaN)而不是补 0；
        // 轨迹里不存在的元素在 finalize 里已被丢掉,这里不会出现
        let (mut e_col, mut mb_col, mut mc_col) = (Vec::new(), Vec::new(), Vec::new());
        for elem in &elems {
            e_col.push((*elem).clone());
            mb_col.push(self.mean_bridge.get(*elem).copied().unwrap_or(f64::NAN));
            mc_col.push(self.mean_cn.get(*elem).copied().unwrap_or(f64::NAN));
        }
        let mut t = Table::new();
        t.push_text("element", e_col)
            .push_num("mean_n_bridge", mb_col)
            .push_num("mean_cn", mc_col);
        out.push(("mean".to_string(), t));

        out
    }
}

fn push_bins(t: &mut Table, bins: &[Bin]) {
    t.push_num("count", bins.iter().map(|b| b.count as f64).collect())
        .push_num("fraction", bins.iter().map(|b| b.fraction).collect())
        .push_num("sd", bins.iter().map(|b| b.sd).collect());
}

// ─── 逐帧中间数据 ─────────────────────────────────────────────────────────────

/// One frame's histograms with the denominators that turn them into fractions.
struct FrameData {
    /// 元素 → (桥接数直方图, 该元素原子数)
    bridge: HashMap<String, (HashMap<u32, usize>, usize)>,
    /// 元素 → (配位数直方图, 该元素原子数)，形成子与修饰子共用
    cn: HashMap<String, (HashMap<u32, usize>, usize)>,
    /// 配体类型 → 计数
    oxy: HashMap<OxyKey, usize>,
    /// 配体原子总数（oxy 的分母）
    n_ligand: usize,
}

// ─── 顶层入口 ─────────────────────────────────────────────────────────────────

/// 对整条轨迹执行网络统计分析。要求每帧有 Cell（PBC）。
///
/// 返回 `None` 如果轨迹为空、没有形成子截断、或所有帧都缺少 Cell。
pub fn calc_network(traj: &Trajectory, params: &TypeParams) -> Option<NetworkResult> {
    if traj.frames.is_empty() || params.cutoffs.is_empty() { return None; }

    let usable: Vec<(&Frame, &ferro_core::Cell)> = traj.frames.iter()
        .filter_map(|frame| frame.cell.as_ref().map(|cell| (frame, cell)))
        .collect();
    if usable.is_empty() { return None; }

    let n_frames = usable.len();
    let n_atoms = usable[0].0.atoms.len();

    let acc = usable
        .par_iter()
        .filter_map(|(frame, cell)| compute_frame(frame, cell, params))
        .fold(Accumulator::default, |mut acc, fd| { acc.push(&fd); acc })
        .reduce(Accumulator::default, |mut a, b| { a.merge(b); a });

    Some(acc.finalize(n_frames, n_atoms, params.clone()))
}

// ─── 单帧计算 ─────────────────────────────────────────────────────────────────

fn compute_frame(
    frame: &Frame,
    cell: &ferro_core::Cell,
    params: &TypeParams,
) -> Option<FrameData> {
    let types = classify_frame(frame, cell, params);

    let mut bridge: HashMap<String, (HashMap<u32, usize>, usize)> = HashMap::new();
    let mut cn: HashMap<String, (HashMap<u32, usize>, usize)> = HashMap::new();
    let mut oxy: HashMap<OxyKey, usize> = HashMap::new();
    let mut n_ligand = 0usize;

    let ligands = params.ligands();

    for (idx, t) in types.iter().enumerate() {
        match t {
            AtomType::Former { elem, bridging, cn: c } => {
                let b = bridge.entry(elem.clone()).or_default();
                *b.0.entry(*bridging).or_insert(0) += 1;
                b.1 += 1;
                let e = cn.entry(elem.clone()).or_default();
                *e.0.entry(*c).or_insert(0) += 1;
                e.1 += 1;
            }
            AtomType::Modifier { elem, cn: c } => {
                let e = cn.entry(elem.clone()).or_default();
                *e.0.entry(*c).or_insert(0) += 1;
                e.1 += 1;
            }
            AtomType::Ligand { .. } => {
                *oxy.entry(oxy_key(t)).or_insert(0) += 1;
                n_ligand += 1;
            }
            // 配体元素的原子若被别的角色覆盖（元素同时是配体与形成子）已在上面计入
            AtomType::Other { elem } => {
                debug_assert!(!ligands.contains(elem), "ligand atom {idx} left unclassified");
            }
        }
    }

    Some(FrameData { bridge, cn, oxy, n_ligand })
}

// ─── 跨帧累加器 ───────────────────────────────────────────────────────────────

#[derive(Default)]
struct Accumulator {
    /// 元素 → { 桥接数 → 矩 }
    bridge: HashMap<String, HashMap<u32, Moments>>,
    /// 元素 → { 配位数 → 矩 }
    cn: HashMap<String, HashMap<u32, Moments>>,
    /// 配体类型 → 矩（BTreeMap：键含 rank，天然有序）
    oxy: BTreeMap<OxyKey, Moments>,
    /// 元素 → (Σ值·计数, Σ计数)，用于池化均值
    bridge_sum: HashMap<String, (f64, usize)>,
    cn_sum: HashMap<String, (f64, usize)>,
}

/// Folds one frame's histogram into `dst`, converting counts to fractions with
/// `total` as the denominator, and updating the pooled sum for the mean.
fn push_hist(
    dst: &mut HashMap<String, HashMap<u32, Moments>>,
    sums: &mut HashMap<String, (f64, usize)>,
    elem: &str,
    hist: &HashMap<u32, usize>,
    total: usize,
) {
    if total == 0 { return; }
    let slot = dst.entry(elem.to_string()).or_default();
    let sum = sums.entry(elem.to_string()).or_insert((0.0, 0));
    for (&v, &c) in hist {
        slot.entry(v).or_default().push(c, c as f64 / total as f64);
        sum.0 += v as f64 * c as f64;
        sum.1 += c;
    }
}

impl Accumulator {
    fn push(&mut self, fd: &FrameData) {
        for (elem, (hist, total)) in &fd.bridge {
            push_hist(&mut self.bridge, &mut self.bridge_sum, elem, hist, *total);
        }
        for (elem, (hist, total)) in &fd.cn {
            push_hist(&mut self.cn, &mut self.cn_sum, elem, hist, *total);
        }
        if fd.n_ligand > 0 {
            for (key, &c) in &fd.oxy {
                self.oxy.entry(key.clone()).or_default()
                    .push(c, c as f64 / fd.n_ligand as f64);
            }
        }
    }

    fn merge(&mut self, other: Self) {
        let merge_nested = |dst: &mut HashMap<String, HashMap<u32, Moments>>,
                            src: HashMap<String, HashMap<u32, Moments>>| {
            for (k, inner) in src {
                let m = dst.entry(k).or_default();
                for (v, mo) in inner { m.entry(v).or_default().merge(&mo); }
            }
        };
        let merge_sums = |dst: &mut HashMap<String, (f64, usize)>,
                          src: HashMap<String, (f64, usize)>| {
            for (k, (s, n)) in src {
                let e = dst.entry(k).or_insert((0.0, 0));
                e.0 += s;
                e.1 += n;
            }
        };
        merge_nested(&mut self.bridge, other.bridge);
        merge_nested(&mut self.cn, other.cn);
        merge_sums(&mut self.bridge_sum, other.bridge_sum);
        merge_sums(&mut self.cn_sum, other.cn_sum);
        for (key, mo) in other.oxy {
            self.oxy.entry(key).or_default().merge(&mo);
        }
    }

    fn finalize(self, n_frames: usize, n_atoms: usize, params: TypeParams) -> NetworkResult {
        let to_dist = |m: &HashMap<u32, Moments>| -> Vec<(u32, Bin)> {
            let mut rows: Vec<(u32, Bin)> =
                m.iter().map(|(&v, mo)| (v, mo.finish(n_frames))).collect();
            rows.sort_by_key(|r| r.0);
            rows
        };
        let mean_of = |sums: &HashMap<String, (f64, usize)>, k: &String| -> f64 {
            match sums.get(k) {
                Some(&(s, n)) if n > 0 => s / n as f64,
                _ => f64::NAN,
            }
        };

        // 只保留轨迹里真正出现过的元素。参数里点名但不存在的元素若留下,
        // 均值会算成 0.0,把「不存在」伪装成「平均值为零」;长表下缺席就是不出行
        let bridge_dist: HashMap<_, _> = self.bridge.iter()
            .map(|(e, m)| (e.clone(), to_dist(m))).collect();
        let mean_bridge: HashMap<_, _> = self.bridge.keys()
            .map(|e| (e.clone(), mean_of(&self.bridge_sum, e))).collect();
        let cn_dist: HashMap<_, _> = self.cn.iter()
            .map(|(e, m)| (e.clone(), to_dist(m))).collect();
        let mean_cn: HashMap<_, _> = self.cn.keys()
            .map(|e| (e.clone(), mean_of(&self.cn_sum, e))).collect();

        let oxy_dist: Vec<(String, Vec<String>, Bin)> = self.oxy.iter()
            .map(|((_, lbl, partners), mo)| {
                (lbl.clone(), partners.clone(), mo.finish(n_frames))
            })
            .collect();

        NetworkResult {
            bridge_dist, mean_bridge, cn_dist, mean_cn, oxy_dist,
            n_frames, n_atoms, params,
        }
    }
}

// ─── 测试 ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame, Trajectory, TypeParams};
    use nalgebra::{Matrix3, Vector3};
    use std::collections::BTreeMap;

    fn atom(elem: &str, x: f64, y: f64, z: f64) -> Atom {
        Atom { element: elem.to_string(), position: Vector3::new(x, y, z),
               label: None, mass: None, magmom: None, charge: None }
    }

    fn cube(atoms: Vec<Atom>) -> Frame {
        let cell = Cell::from_matrix(Matrix3::from_diagonal(&Vector3::new(20.0, 20.0, 20.0)));
        Frame { atoms, cell: Some(cell), ..Frame::default() }
    }

    fn traj_of(frames: Vec<Frame>) -> Trajectory {
        Trajectory { frames, metadata: Default::default() }
    }

    fn make_params(p_o: f64) -> TypeParams {
        let mut c = BTreeMap::new();
        c.insert(("P".into(), "O".into()), p_o);
        TypeParams { cutoffs: c, modifier_cutoffs: BTreeMap::new() }
    }

    /// 取配体类型的 label → count 映射
    fn oxy_counts(res: &NetworkResult) -> HashMap<&str, usize> {
        let mut m: HashMap<&str, usize> = HashMap::new();
        for (l, _, b) in &res.oxy_dist { *m.entry(l.as_str()).or_insert(0) += b.count; }
        m
    }

    /// P–O–P 体系：两个 Q1 的 P，一个桥氧
    #[test]
    fn test_q1_system() {
        let frame = cube(vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0), // 桥氧
            atom("P",  3.2, 0.0, 0.0),
            atom("O", -1.6, 0.0, 0.0), // NBO(P1)
            atom("O",  4.8, 0.0, 0.0), // NBO(P2)
        ]);
        let res = calc_network(&traj_of(vec![frame]), &make_params(2.3)).unwrap();

        let bp = &res.bridge_dist["P"];
        assert_eq!(bp.len(), 1);
        assert_eq!(bp[0].0, 1);          // 桥接数 = 1
        assert_eq!(bp[0].1.count, 2);    // 两个 P

        let oxy = oxy_counts(&res);
        assert_eq!(oxy["O_b"], 1);
        assert_eq!(oxy["O_n"], 2);
    }

    /// 孤立 PO4（Q0）
    #[test]
    fn test_q0_system() {
        let frame = cube(vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0),
            atom("O", -1.6, 0.0, 0.0),
            atom("O",  0.0, 1.6, 0.0),
            atom("O",  0.0,-1.6, 0.0),
        ]);
        let res = calc_network(&traj_of(vec![frame]), &make_params(2.3)).unwrap();
        assert_eq!(res.bridge_dist["P"][0].0, 0);
        assert!(res.oxy_dist.iter().all(|(l, _, _)| l == "O_n"));
    }

    /// 异质桥（P–O–Al）与同质桥（P–O–P）共用标签 `O_b`，但伙伴元素以
    /// **数据列**区分 —— 标签合并，信息不丢。
    #[test]
    fn test_bridge_partners_are_data_not_label() {
        let frame = cube(vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0), // P–O–Al
            atom("Al", 3.2, 0.0, 0.0),
            atom("P",  0.0, 5.0, 0.0),
            atom("O",  1.6, 5.0, 0.0), // P–O–P
            atom("P",  3.2, 5.0, 0.0),
        ]);
        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("P".to_string(), "O".to_string()), 2.3);
        cutoffs.insert(("Al".to_string(), "O".to_string()), 2.3);
        let params = TypeParams { cutoffs, modifier_cutoffs: BTreeMap::new() };
        let res = calc_network(&traj_of(vec![frame]), &params).unwrap();

        // 标签都是 O_b
        assert!(res.oxy_dist.iter().all(|(l, _, _)| l == "O_b"));
        // 但按伙伴分成两行
        let mut rows: Vec<(&Vec<String>, usize)> =
            res.oxy_dist.iter().map(|(_, p, b)| (p, b.count)).collect();
        rows.sort();
        assert_eq!(rows.len(), 2, "P-O-Al 与 P-O-P 必须是两行: {rows:?}");
        assert_eq!(rows[0].0, &vec!["Al".to_string(), "P".to_string()]);
        assert_eq!(rows[1].0, &vec!["P".to_string(), "P".to_string()]);
    }

    /// 修饰子只进 CN 表，且**不影响配体分类** —— 挨着修饰子的非桥氧仍是 `O_n`。
    /// 这正是 `--modifier` 存在的理由。
    #[test]
    fn test_modifier_counts_cn_only() {
        let frame = cube(vec![
            atom("P",   0.0, 0.0, 0.0),
            atom("O",   1.6, 0.0, 0.0),
            atom("O",  -1.6, 0.0, 0.0),
            atom("Zn",  1.6, 2.0, 0.0), // 距一个 O 2.0 Å，距另一个 3.8 Å
        ]);
        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("P".to_string(), "O".to_string()), 2.3);
        let mut modifier_cutoffs = BTreeMap::new();
        modifier_cutoffs.insert(("Zn".to_string(), "O".to_string()), 2.6);
        let params = TypeParams { cutoffs, modifier_cutoffs };
        let res = calc_network(&traj_of(vec![frame]), &params).unwrap();

        assert_eq!(res.cn_dist["Zn"][0].0, 1, "Zn 配位数 1");
        assert_eq!(res.cn_dist["Zn"][0].1.count, 1);
        assert!(!res.bridge_dist.contains_key("Zn"), "修饰子没有桥接数");
        assert!(res.oxy_dist.iter().all(|(l, _, _)| l == "O_n"),
                "Zn 不能把非桥氧变成桥氧");
    }

    /// 参数里点名但轨迹里没有的元素不能出现在结果里 —— 否则均值会算成 0.0，
    /// 把「这个体系没有 Al」写成「Al 的平均配位数是 0」。
    #[test]
    fn test_absent_element_produces_no_rows() {
        let frame = cube(vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0),
            atom("O", -1.6, 0.0, 0.0),
        ]);
        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("P".to_string(), "O".to_string()), 2.3);
        cutoffs.insert(("Al".to_string(), "O".to_string()), 2.3);
        let mut modifier_cutoffs = BTreeMap::new();
        modifier_cutoffs.insert(("Zn".to_string(), "O".to_string()), 2.6);
        let params = TypeParams { cutoffs, modifier_cutoffs };
        let res = calc_network(&traj_of(vec![frame]), &params).unwrap();

        for absent in ["Al", "Zn"] {
            assert!(!res.bridge_dist.contains_key(absent));
            assert!(!res.cn_dist.contains_key(absent));
            assert!(!res.mean_bridge.contains_key(absent));
            assert!(!res.mean_cn.contains_key(absent));
        }
        let tables = res.to_tables();
        let mean = &tables.iter().find(|(n, _)| n == "mean").unwrap().1;
        assert_eq!(mean.n_rows(), 1, "只有 P 一行");
    }

    /// 恒定体系的 sd 必须恰好为 0，单帧必须是 NaN。
    #[test]
    fn test_sd_is_zero_for_a_constant_system() {
        let make = || cube(vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0),
            atom("O", -1.6, 0.0, 0.0),
        ]);
        let params = make_params(2.3);

        let one = calc_network(&traj_of(vec![make()]), &params).unwrap();
        assert!(one.oxy_dist[0].2.sd.is_nan(), "单帧无法给出样本标准差");

        let three = calc_network(&traj_of(vec![make(), make(), make()]), &params).unwrap();
        let bin = three.oxy_dist[0].2;
        assert_eq!(bin.count, 6, "3 帧 × 2 个 O");
        assert!((bin.fraction - 1.0).abs() < 1e-12);
        assert!(bin.sd < 1e-12, "三帧完全相同,sd 必须为 0,实得 {}", bin.sd);
    }

    /// 逐帧比例在变时，fraction 是逐帧比例的平均，sd 是它们的样本标准差。
    /// 手算对照：两帧的 O_b 比例分别是 1/3 与 0 → mean 1/6，sd = 1/(3√2)。
    #[test]
    fn test_fraction_and_sd_are_per_frame_statistics() {
        // 帧 1：P–O–P 桥 + 两个 NBO  → 3 个 O，其中 1 个桥氧
        let f1 = cube(vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0),
            atom("P",  3.2, 0.0, 0.0),
            atom("O", -1.6, 0.0, 0.0),
            atom("O",  4.8, 0.0, 0.0),
        ]);
        // 帧 2：桥氧被拉开 → 3 个 O 全是 NBO
        let f2 = cube(vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0),
            atom("P",  6.0, 0.0, 0.0),
            atom("O", -1.6, 0.0, 0.0),
            atom("O",  7.6, 0.0, 0.0),
        ]);
        let res = calc_network(&traj_of(vec![f1, f2]), &make_params(2.3)).unwrap();

        let ob = res.oxy_dist.iter().find(|(l, _, _)| l == "O_b").unwrap().2;
        assert_eq!(ob.count, 1, "只有帧 1 有一个桥氧");
        assert!((ob.fraction - 1.0 / 6.0).abs() < 1e-12,
                "(1/3 + 0)/2 = 1/6，实得 {}", ob.fraction);
        let expect_sd = (1.0 / 3.0f64) / 2.0f64.sqrt();
        assert!((ob.sd - expect_sd).abs() < 1e-12,
                "ddof=1 的样本标准差应为 1/(3√2)={expect_sd}，实得 {}", ob.sd);
    }

    /// Welford 对恒定序列必须给出**精确** 0。两遍求和的 `Σf² − N·mean²`
    /// 在这里会返回 ~2e-10 的抵消噪声，而 sd 列里的假抖动会被当成物理。
    #[test]
    fn test_sd_of_a_constant_bin_is_exactly_zero() {
        let mut m = Moments::default();
        // 比例恒为 25/1860 —— 取自参考轨迹里 P 的 Q0 档
        let frac = 25.0 / 1860.0;
        for _ in 0..5 { m.push(25, frac); }
        let bin = m.finish(5);
        assert_eq!(bin.sd, 0.0, "恒定序列的样本标准差必须精确为 0,实得 {}", bin.sd);
        assert_eq!(bin.count, 125);

        // 并行归并路径同样精确
        let (mut a, mut b) = (Moments::default(), Moments::default());
        for _ in 0..2 { a.push(25, frac); }
        for _ in 0..3 { b.push(25, frac); }
        a.merge(&b);
        assert_eq!(a.finish(5).sd, 0.0, "归并后仍须精确为 0");
    }

    /// 缺席帧按 f = 0 计入，而不是被当作缺失值跳过。
    #[test]
    fn test_absent_frames_count_as_zero() {
        let mut m = Moments::default();
        m.push(1, 1.0);           // 4 帧里只有 1 帧出现,比例 1.0
        let bin = m.finish(4);
        assert!((bin.fraction - 0.25).abs() < 1e-15, "(1+0+0+0)/4 = 0.25");
        // 样本标准差: 值为 [1,0,0,0], mean=0.25, Σ(x-μ)²=0.75, /3 → 0.25 → sd=0.5
        assert!((bin.sd - 0.5).abs() < 1e-15, "实得 {}", bin.sd);
    }
}
