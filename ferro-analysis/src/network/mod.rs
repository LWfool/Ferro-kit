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
//! let params = TypeParams::new(cutoffs, BTreeMap::new());
//! let result = calc_network(&traj, &params).unwrap();
//! ```

pub use ferro_core::{AtomType, CutoffTable, TypeParams};

use ferro_core::{classify_frame_detailed, Frame, Table, Trajectory};
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};

/// Qn partner key: the Qn value, plus how many of those bridges lead to each partner
/// element.  `(2, {Al: 1, P: 1})` is the Q²(1Al) of the NMR literature.
///
/// Keeping the decomposition per *atom* rather than per bridge is what makes it
/// irrecoverable from the linkage table: two P–O–Al bridges could come from one P
/// with `m_Al = 2` or from two P atoms with `m_Al = 1` each.
pub type BridgeKey = (u32, BTreeMap<String, u32>);

/// One end of a bridge: `(element, homopolar connection count, coordination number)`.
///
/// Both numbers are kept for both ends.  The reference Python stores one number per
/// element — Qn for P, CN for Al — which answers "what Qn of P does a 4-coordinate Al
/// bond to" but not "what is the bridging count of that 4-coordinate Al".
pub type SiteState = (String, u32, u32);

/// Linkage key: the two ends in canonical order, the **bridging ligand's element**,
/// and how many formers that ligand joins (2 for a true bridge, ≥3 for a tricluster).
///
/// The ligand element is part of the key because a system with two ligand species
/// (`--Al-O=2.4 --Al-F=2.1`) has genuinely different bridges: merging Al–O–P with
/// Al–F–P would report a count no experiment could reproduce.
pub type LinkKey = (SiteState, SiteState, String, usize);

/// Ligand distribution key: `(element, sort rank, label, partner elements)`.
///
/// Element first so a two-ligand system groups by species; the rank then comes from
/// [`AtomType::class_rank`], so a `BTreeMap` keyed on this is sorted by construction —
/// no post-hoc `sort_by` that re-derives the order from the label text.  Partners are
/// kept only for the 0/1/2-former cases; a tricluster collapses to one row per label,
/// its pairs being reported by the linkage table.
type OxyKey = (String, u8, String, Vec<String>);

/// Ligand type key: `(class rank, label)`.  The rank leads so a `BTreeMap` orders
/// `_f` < `_n` < `_b` < `_t` by construction rather than alphabetically by label.
type LigandKey = (u8, String);

/// Per-frame ligand histogram: ligand element → (type histogram, that element's atom
/// count).  The count is the denominator, and it is **per element** — see
/// [`NetworkResult::ligand_dist`].
type LigandHist = HashMap<String, (HashMap<LigandKey, usize>, usize)>;

fn oxy_key(t: &AtomType) -> OxyKey {
    let partners = match t {
        AtomType::Ligand { partners, .. } if partners.len() <= 2 => partners.clone(),
        _ => Vec::new(),
    };
    (t.element().to_string(), t.class_rank(), t.label(), partners)
}

// ─── 两套标签词汇 ─────────────────────────────────────────────────────────────
//
// **原子词汇** `P_2` / `Al_4`：指一个原子。用于导出轨迹（必须能被 dump reader 按首个
// 下划线拆回 element）与 `linkage` 表（它描述的是原子之间的连接）。
//
// **单元词汇** `P-Q2` / `Al_4`：指一类结构单元。用于 `qn` / `qn_partner` /
// `composition` 三张分布表——它们数的是「有多少个 Q2 单元」，不是「哪个原子连着哪个」。
//
// 两者对非 Qn 形成子恰好相同（都是 `Al_4`），因为 Al 没有单元记号可用。

/// Atom vocabulary: routed through [`AtomType::label`] so the tables and the exported
/// trajectory cannot disagree about what `Al_4` means.
fn former_label(elem: &str, qn: Option<u32>, cn: u32) -> String {
    AtomType::Former {
        elem: elem.to_string(),
        qn,
        n_bo: qn.unwrap_or(0),
        cn,
        bridges_to: BTreeMap::new(),
    }
    .label()
}

/// Unit vocabulary: `P-Q2` for a Qn former, the atom label otherwise.
///
/// The element is kept in front because `Q2` alone is ambiguous the moment a run has
/// two Qn formers (B + P, Si + Al): the rows would collide and no adjacent column
/// could separate two *different* elements sharing a Qn.
fn species_label(elem: &str, qn: Option<u32>, cn: u32) -> String {
    match qn {
        Some(n) => format!("{elem}-Q{n}"),
        None => former_label(elem, None, cn),
    }
}

/// Label of one bridge end, given the run's Qn convention.  Atom vocabulary — a bridge
/// joins two *atoms*, whereas Qn names a unit that contains several.
fn site_label(params: &TypeParams, (elem, n, cn): &SiteState) -> String {
    former_label(elem, params.is_qn_former(elem).then_some(*n), *cn)
}

/// Ligand row label: the formers it joins, with the ligand's own type in the middle.
///
/// `Al-O_b-P`, `P-O_n`, `O_f`, `O_t`.  Missing ends are left off rather than filled
/// with a placeholder — a free ligand has no formers, and a tricluster has no single
/// pair (its pairs are in the linkage table).
fn ligand_row_label(ty: &str, partners: &[String]) -> String {
    match partners {
        [a, b] => format!("{a}-{ty}-{b}"),
        [a] => format!("{a}-{ty}"),
        _ => ty.to_string(),
    }
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
    /// former_elem → Qn 分布，按数值升序。**边际分布** —— 主产物，
    /// 一行一个 (形成子, Qn)。
    ///
    /// 只含 [`TypeParams::qn_formers`] 里的元素。Al 之类由配位数刻画的形成子
    /// 不在此表：Qn 是四面体形成子的记号，给 Al 报一个 Qn 等于发明一个
    /// 文献里不存在的量
    pub qn_dist: HashMap<String, Vec<(u32, Bin)>>,
    /// former_elem → Qn × 伙伴元素的**联合分布**（Q^n(mAl) 记号）。
    /// 与 `qn_dist` 是两个粒度：后者是前者对伙伴维度的边际
    pub qn_partner_dist: HashMap<String, Vec<(BridgeKey, Bin)>>,
    /// former_elem → 平均 Qn（同样只含 Qn 形成子）
    pub mean_qn: HashMap<String, f64>,
    /// 元素 → 总配位数分布（所有配体类型之和）。**含修饰子** ——
    /// 修饰子没有桥接数，配位数就是描述它的全部
    pub cn_dist: HashMap<String, Vec<(u32, Bin)>>,
    /// 元素 → 平均总配位数（含修饰子）
    pub mean_cn: HashMap<String, f64>,
    /// Former element → mean number of **bridging oxygens** per atom.
    ///
    /// A third quantity next to `mean_qn` (homopolar connections) and `mean_cn`
    /// (all ligands): it counts bridging ligands, so a tricluster adds one here but
    /// two to `qn + Σm`.  Reported alongside `mean_qn` because the two coincide in
    /// a corner-sharing, tricluster-free network and diverge exactly where the
    /// literature's conventional Qn analysis stops applying.
    pub mean_n_bo: HashMap<String, f64>,
    /// 配体类型分布：`(label, partner_elements, bin)`，按 (元素, `class_rank`) 排序。
    /// `fraction` 的分母是**该配体元素**的原子数
    pub oxy_dist: Vec<(String, Vec<String>, Bin)>,
    /// 配体元素 → 类型分布（`O_f`/`O_n`/`O_b`/`O_t`），已对伙伴维度聚合。
    /// 与 `oxy_dist` 是两个粒度，`sd` 各自累加
    pub ligand_dist: HashMap<String, Vec<(String, Bin)>>,
    /// 桥联统计：每个桥联配体贡献一条「两端位点状态」的观测，
    /// 只存规范半边（两端按 `SiteState` 排序，小的在前）
    pub linkage: Vec<(LinkKey, Bin)>,
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
    /// | table | one row per | key columns |
    /// |---|---|---|
    /// | `qn` | Qn former × Qn | `label, former, qn` |
    /// | `qn_partner` | Qn former × Qn × partner split | `label, former, qn, m_<X>…` |
    /// | `ligand_type` | ligand type × partner pair | `label, former_a, former_b` |
    /// | `coordination` | element × coordination number | `element, cn` |
    /// | `average` | element | `element, mean_qn, mean_cn` |
    /// | `linkage` | bridge × both site states | `linkage, ligand, elem_a, n_bridge_a, cn_a, …` |
    ///
    /// Every distribution table then carries `count, fraction, sd`.
    ///
    /// **Why six tables.**  `qn` is the marginal of `qn_partner` over the partner
    /// dimension, and is a separate table rather than a `groupby` because the plain
    /// Qn distribution is a primary product: it has to be readable by opening the
    /// file.  `sd` genuinely has to be re-accumulated rather than summed across
    /// partner rows — the variance of a sum is not the sum of variances when the
    /// terms are correlated.  The three distributions have different value-column
    /// semantics, so merging them would make a `groupby` on the value column
    /// meaningless, and `average` is one row per element rather than per value,
    /// a different granularity again.
    ///
    /// **`qn` and `coordination` are different quantities.**  `qn` counts only the
    /// bridging ligands and exists only for [`TypeParams::qn_formers`];
    /// `coordination` counts every ligand inside the cutoff and exists for every
    /// former and modifier.  A former with no non-bridging ligand makes them equal,
    /// which is a fact about that system, not a definition.
    ///
    /// The `qn` and `qn_partner` tables are **omitted entirely** when no former is a
    /// Qn element: a header-only file reads as "measured, and it was zero".
    pub fn to_tables(&self) -> Vec<(String, Table)> {
        let mut out = Vec::with_capacity(6);

        let mut formers: Vec<&String> = self.qn_dist.keys().collect();
        formers.sort();
        let mut elems: Vec<&String> = self.cn_dist.keys().collect();
        elems.sort();

        // qn：Qn 的**边际**分布 —— 打开文件就能读的那张
        if !formers.is_empty() {
            let (mut l_col, mut f_col, mut v_col) = (Vec::new(), Vec::new(), Vec::new());
            let mut bins = Vec::new();
            for former in &formers {
                for (n, b) in &self.qn_dist[*former] {
                    l_col.push(species_label(former, Some(*n), 0));
                    f_col.push((*former).clone());
                    v_col.push(*n as f64);
                    bins.push(*b);
                }
            }
            let mut t = Table::new();
            t.push_text("label", l_col).push_text("former", f_col).push_num("qn", v_col);
            push_bins(&mut t, &bins);
            for line in [
                "table   : Qn speciation — one row per (former, Qn).",
                "          n = HOMOPOLAR connections only (P-O-P for a P), following",
                "          the literature's Q^n_m / P^n_(mAl,xB): the heteropolar ones",
                "          are the m values in network_qn_partner.csv, and the TOTAL",
                "          bridge count is n + Sum(m), not n.  A P-Q1 row therefore",
                "          still contains phosphorus with any number of P-O-Al bridges",
                "label   : structural unit, <element>-Q<n>.  NOTE this is the UNIT",
                "          vocabulary; the exported trajectory labels the same site as",
                "          P_2 (atom vocabulary), which is what -x/-y select on",
                "former  : network former element",
                "qn      : homopolar connections on the site (the literature's n)",
                "count   : occurrences summed over all frames",
                "fraction: mean over frames of that frame's fraction of this former",
                "sd      : sample std (ddof=1) of the per-frame fraction; a spread",
                "          between snapshots, NOT a standard error — MD frames correlate",
                "note    : formers described by coordination number instead (Al, …) are",
                "          absent here by design; read network_coordination.csv for them",
            ] { t.meta_line(line); }
            out.push(("qn".to_string(), t));
        }

        // qn_partner：Qn × 伙伴元素的**联合**分布，文献的 Q^n(mAl) 记号
        if !self.qn_partner_dist.is_empty() {
            let partner_elems: Vec<String> = {
                let mut s: Vec<String> = self.qn_partner_dist.values()
                    .flat_map(|rows| rows.iter().flat_map(|((_, to), _)| to.keys().cloned()))
                    .collect();
                s.sort();
                s.dedup();
                s
            };
            let mut pf: Vec<&String> = self.qn_partner_dist.keys().collect();
            pf.sort();
            let (mut l_col, mut f_col, mut v_col) = (Vec::new(), Vec::new(), Vec::new());
            let mut m_cols: Vec<Vec<f64>> = vec![Vec::new(); partner_elems.len()];
            let mut bins = Vec::new();
            for former in &pf {
                for ((n, to), b) in &self.qn_partner_dist[*former] {
                    l_col.push(species_label(former, Some(*n), 0));
                    f_col.push((*former).clone());
                    v_col.push(*n as f64);
                    for (i, p) in partner_elems.iter().enumerate() {
                        m_cols[i].push(to.get(p).copied().unwrap_or(0) as f64);
                    }
                    bins.push(*b);
                }
            }
            let mut t = Table::new();
            t.push_text("label", l_col).push_text("former", f_col).push_num("qn", v_col);
            for (p, col) in partner_elems.iter().zip(m_cols) {
                t.push_num(format!("m_{p}"), col);
            }
            push_bins(&mut t, &bins);
            for line in [
                "table   : Q^n_m — the Qn distribution split by partner element.",
                "          network_qn.csv is this table's marginal over the m_ columns.",
                "label   : structural unit, <element>-Q<n>; repeats across the rows",
                "          that share a Qn.  The partner split is kept in columns",
                "          rather than folded into the label: with two heteropolar",
                "          partners the label would grow to P-Q1(2Al,1B), and these",
                "          columns are what you filter on anyway",
                "qn      : homopolar connections on the site (the literature's n)",
                "m_<X>   : HETEROPOLAR connections to element X (the literature's m).",
                "          The former's own element has no column: it would repeat qn",
                "          exactly.  Total bridging connections = qn + Sum(m_).",
                "          A ligand shared by 3 formers connects this site to TWO",
                "          partners and contributes 2, so qn + Sum(m_) can exceed the",
                "          number of bridging oxygens — see n_bo in the composition",
                "          table for that count",
                "count/fraction/sd: as in network_qn.csv, but re-accumulated at this",
                "          granularity — sd cannot be summed across partner rows,",
                "          the terms are correlated",
            ] { t.meta_line(line); }
            out.push(("qn_partner".to_string(), t));
        }

        // ligand_type：伙伴元素以数据列给出,不编码进标签
        let (mut lbl, mut ty) = (Vec::new(), Vec::new());
        let (mut fa, mut fb) = (Vec::new(), Vec::new());
        let mut bins = Vec::new();
        for (label, partners, b) in &self.oxy_dist {
            lbl.push(ligand_row_label(label, partners));
            ty.push(label.clone());
            fa.push(partners.first().cloned().unwrap_or_default());
            fb.push(partners.get(1).cloned().unwrap_or_default());
            bins.push(*b);
        }
        let mut t = Table::new();
        t.push_text("label", lbl).push_text("type", ty)
            .push_text("former_a", fa).push_text("former_b", fb);
        push_bins(&mut t, &bins);
        for line in [
            "table   : ligand speciation — how many formers each ligand atom touches.",
            "label   : the whole row read as one string, <former_a>-<type>-<former_b>.",
            "          Missing ends are left off, not padded: O_f, P-O_n, Al-O_b-P, O_t",
            "type    : <element>_f free (0 formers)   <element>_n non-bridging (1)",
            "          <element>_b bridging (2)       <element>_t tricluster (3+)",
            "former_a: the formers it bridges, kept as their own columns so",
            "former_b: former_a=='Al' and former_b=='Al' stays a plain query.",
            "          Both empty for _f, only former_a set for _n.  A tricluster",
            "          leaves both empty: it has no single pair, its pairs are",
            "          enumerated in network_linkage.csv",
            "fraction: denominator is that LIGAND ELEMENT's atom count — the",
            "          literature's BO-fraction convention.  In a two-ligand system",
            "          (--Al-O and --Al-F) O_b is a fraction of O, not of O+F",
        ] { t.meta_line(line); }
        out.push(("ligand_type".to_string(), t));

        // coordination：形成子与修饰子共用一张表
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
        for line in [
            "table   : coordination number distribution, one row per (element, cn).",
            "element : every former AND every modifier — a modifier has no Qn, so",
            "          this table is the whole of what describes it",
            "cn      : ligands within the cutoff, counting bridging and non-bridging",
            "          alike.  NOT the same quantity as network_qn.csv's qn, which",
            "          counts bridging ligands only; the two coincide only for a",
            "          former that happens to carry no non-bridging ligand",
            "fraction: denominator is that element's atom count, so it sums to 1",
            "          within each element",
        ] { t.meta_line(line); }
        out.push(("coordination".to_string(), t));

        // composition：一物种一行,分母恒为**该元素**的原子数。
        // 每个元素只出现一种刻画:Qn 形成子出 Qn,其余形成子与修饰子出配位数,
        // 配体出类型。故「每个元素的 fraction 求和为 1」是一条可核对的恒等式
        let (mut l_col, mut e_col) = (Vec::new(), Vec::new());
        let mut bins = Vec::new();
        for former in &formers {
            for (n, b) in &self.qn_dist[*former] {
                l_col.push(species_label(former, Some(*n), 0));
                e_col.push((*former).clone());
                bins.push(*b);
            }
        }
        for elem in &elems {
            // Qn 形成子已由上一段给出,不再按配位数重复一遍
            if self.qn_dist.contains_key(*elem) { continue; }
            for (v, b) in &self.cn_dist[*elem] {
                l_col.push(species_label(elem, None, *v));
                e_col.push((*elem).clone());
                bins.push(*b);
            }
        }
        let mut ligs: Vec<&String> = self.ligand_dist.keys().collect();
        ligs.sort();
        for elem in ligs {
            for (label, b) in &self.ligand_dist[elem] {
                l_col.push(label.clone());
                e_col.push(elem.clone());
                bins.push(*b);
            }
        }
        let mut t = Table::new();
        t.push_text("label", l_col).push_text("element", e_col);
        push_bins(&mut t, &bins);
        for line in [
            "table   : structural composition — one row per species, at a glance.",
            "          A digest of the other tables, not a new measurement.",
            "label   : P-Q2 (Qn unit)  Al_4 (4-coordinate former)  O_b (bridging",
            "          ligand)  Zn_4 (4-coordinate modifier)",
            "element : the element the species belongs to; group on it",
            "fraction: denominator is that ELEMENT's atom count.  Q2 as a fraction of",
            "          all P, Al_4 as a fraction of all Al, O_b as a fraction of all O",
            "          — so each element's rows sum to 1, which is worth checking",
            "note    : every element is characterised ONE way — Qn formers by Qn,",
            "          other formers and modifiers by coordination number, ligands by",
            "          type.  P therefore has no cn rows here; see",
            "          network_coordination.csv for those",
            "sd      : re-accumulated per frame at THIS granularity, never summed from",
            "          the source tables — O_b aggregates three partner rows, and the",
            "          variance of a sum is not the sum of variances",
            "means   : mean Qn and mean CN are in the [inputs] block above, one line",
            "          per input; they are per-input scalars, not rows",
        ] { t.meta_line(line); }
        out.push(("composition".to_string(), t));

        // linkage：一行一种「配体 × 两端位点状态」的组合
        let mut txt: [Vec<String>; 4] = Default::default();
        let mut nums: [Vec<f64>; 5] = Default::default();
        let mut bins = Vec::new();
        for ((a, b_site, ligand, n_formers), bin) in &self.linkage {
            // 展示列走 AtomType::label,与导出轨迹用同一套词汇
            txt[0].push(format!("{}-{ligand}-{}",
                                site_label(&self.params, a), site_label(&self.params, b_site)));
            txt[1].push(ligand.clone());
            txt[2].push(a.0.clone());
            txt[3].push(b_site.0.clone());
            nums[0].push(a.1 as f64);
            nums[1].push(a.2 as f64);
            nums[2].push(b_site.1 as f64);
            nums[3].push(b_site.2 as f64);
            nums[4].push(*n_formers as f64);
            bins.push(*bin);
        }
        let [c_link, c_lig, c_a, c_b] = txt;
        let [n_a, cn_a, n_b, cn_b, nf] = nums;
        let mut t = Table::new();
        t.push_text("linkage", c_link)
            .push_text("ligand", c_lig)
            .push_text("elem_a", c_a)
            .push_num("qn_a", n_a)
            .push_num("cn_a", cn_a)
            .push_text("elem_b", c_b)
            .push_num("qn_b", n_b)
            .push_num("cn_b", cn_b)
            .push_num("n_formers", nf);
        push_bins(&mut t, &bins);
        for line in [
            "table   : bridge connectivity — one row per (ligand, state of each end).",
            "linkage : human-readable form, <label>-<ligand>-<label>.  The digit is",
            "          the Qn for a Qn former and the COORDINATION NUMBER for any",
            "          other former, which is the literature's own convention",
            "          (Q2 vs Al[4]).  Filter on the numeric columns, not this text",
            "ligand  : element of the bridging atom.  Al-O-P and Al-F-P are different",
            "          bridges and never share a row",
            "qn_*    : HOMOPOLAR connections on that end (P-O-P for a P), the n of",
            "          the literature's Q^n_m.  Defined for every former, including",
            "          those not reported with a Qn speciation (Al-O-Al for an Al)",
            "cn_*    : coordination number of that end",
            "n_formers: 2 for a true bridge.  A ligand shared by n formers is a",
            "          tricluster and contributes C(n,2) rows, all tagged with its n",
            "note    : each pair is stored ONCE, canonically ordered, because a bridge",
            "          has no direction — so a row sum is not a site's total",
            "          involvement.  Both ends carry both numbers, so \"what is the",
            "          homopolar count of that 4-coordinate Al\" is answerable",
        ] { t.meta_line(line); }
        out.push(("linkage".to_string(), t));

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
    /// Qn 元素 → (Qn 直方图, 该元素原子数) —— 边际。非 Qn 形成子不在内
    qn: HashMap<String, (HashMap<u32, usize>, usize)>,
    /// Qn 元素 → (Qn×伙伴分解直方图, 该元素原子数) —— 联合
    qn_partner: HashMap<String, (HashMap<BridgeKey, usize>, usize)>,
    /// 元素 → (配位数直方图, 该元素原子数)，形成子与修饰子共用
    cn: HashMap<String, (HashMap<u32, usize>, usize)>,
    /// 配体类型 → 计数
    oxy: HashMap<OxyKey, usize>,
    /// 配体元素 → (类型直方图, 该元素原子数) —— composition 的配体部分,
    /// 且是 oxy 的分母来源。逐元素而非全体配体:「桥氧占全部氧的比例」是文献的
    /// BO fraction 口径,「桥氧占 O+F 的比例」不对应任何常用量
    ligand: LigandHist,
    /// 元素 → (桥氧个数直方图, 该元素原子数)，只为池化 mean_n_bo
    n_bo: HashMap<String, (HashMap<u32, usize>, usize)>,
    /// 桥联组合 → 计数
    linkage: HashMap<LinkKey, usize>,
    /// 本帧的桥联观测总数（linkage 的分母）
    n_links: usize,
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
    let ft = classify_frame_detailed(frame, cell, params);
    let types = &ft.types;

    let mut bridge: HashMap<String, (HashMap<u32, usize>, usize)> = HashMap::new();
    let mut partner: HashMap<String, (HashMap<BridgeKey, usize>, usize)> = HashMap::new();
    let mut cn: HashMap<String, (HashMap<u32, usize>, usize)> = HashMap::new();
    let mut oxy: HashMap<OxyKey, usize> = HashMap::new();
    let mut ligand: LigandHist = HashMap::new();
    let mut n_bo_hist: HashMap<String, (HashMap<u32, usize>, usize)> = HashMap::new();

    let ligands = params.ligands();

    for (idx, t) in types.iter().enumerate() {
        match t {
            // `qn` 是 Some 才有 Qn 分布可言；非 Qn 形成子（Al）只进配位数表。
            // 判据取自类型自身而不是重查 params，因为类型携带的正是它被分类时
            // 所用的那套约定
            AtomType::Former { elem, qn: q, cn: c, bridges_to, n_bo } => {
                let nb = n_bo_hist.entry(elem.clone()).or_default();
                *nb.0.entry(*n_bo).or_insert(0) += 1;
                nb.1 += 1;
                if let Some(n) = q {
                    let b = bridge.entry(elem.clone()).or_default();
                    *b.0.entry(*n).or_insert(0) += 1;
                    b.1 += 1;
                    let p = partner.entry(elem.clone()).or_default();
                    // 只留异核连接：同元素那一列恒等于 n 本身,是重复信息。
                    // 文献 Q^n_mAl,xB 也没给自己留下标位
                    let mut hetero = bridges_to.clone();
                    hetero.remove(elem);
                    *p.0.entry((*n, hetero)).or_insert(0) += 1;
                    p.1 += 1;
                }
                let e = cn.entry(elem.clone()).or_default();
                *e.0.entry(*c).or_insert(0) += 1;
                e.1 += 1;
            }
            AtomType::Modifier { elem, cn: c } => {
                let e = cn.entry(elem.clone()).or_default();
                *e.0.entry(*c).or_insert(0) += 1;
                e.1 += 1;
            }
            AtomType::Ligand { elem, .. } => {
                *oxy.entry(oxy_key(t)).or_insert(0) += 1;
                let l = ligand.entry(elem.clone()).or_default();
                *l.0.entry((t.class_rank(), t.label())).or_insert(0) += 1;
                l.1 += 1;
            }
            // 配体元素的原子若被别的角色覆盖（元素同时是配体与形成子）已在上面计入
            AtomType::Other { elem } => {
                debug_assert!(!ligands.contains(elem), "ligand atom {idx} left unclassified");
            }
        }
    }

    // 桥联：每个连 ≥2 个形成子的配体，把两端形成子的位点状态两两记下。
    // 三配位配体贡献 C(3,2)=3 条记录，靠 n_formers 列区分,不静默丢掉
    let site = |i: usize| -> Option<SiteState> {
        match &types[i] {
            AtomType::Former { elem, cn, bridges_to, .. } => {
                // 同元素连接数 = 文献的 n。对非 Qn 形成子(Al)同样有定义
                // (Al-O-Al 数),故 linkage 两端都能给出这个量
                let n = bridges_to.get(elem).copied().unwrap_or(0);
                Some((elem.clone(), n, *cn))
            }
            _ => None,
        }
    };
    let mut linkage: HashMap<LinkKey, usize> = HashMap::new();
    let mut n_links = 0usize;
    for (&lig_idx, formers) in &ft.ligand_formers {
        if formers.len() < 2 { continue; }
        let n_formers = formers.len();
        // 桥中间那个原子的元素。多配体体系里 Al-O-P 与 Al-F-P 是两种桥,
        // 合并计数会报出一个实验无法复现的数
        let ligand_elem = &frame.atoms[lig_idx].element;
        for i in 0..n_formers {
            for j in (i + 1)..n_formers {
                let (Some(a), Some(b)) = (site(formers[i]), site(formers[j])) else { continue };
                // 规范半边：桥联无方向，A-O-B 与 B-O-A 是同一条观测
                let (a, b) = if a <= b { (a, b) } else { (b, a) };
                *linkage.entry((a, b, ligand_elem.clone(), n_formers)).or_insert(0) += 1;
                n_links += 1;
            }
        }
    }

    Some(FrameData { qn: bridge, qn_partner: partner, cn, oxy, ligand,
                     n_bo: n_bo_hist, linkage, n_links })
}

// ─── 跨帧累加器 ───────────────────────────────────────────────────────────────

#[derive(Default)]
struct Accumulator {
    /// Qn 元素 → { Qn → 矩 }（边际）
    qn: HashMap<String, BTreeMap<u32, Moments>>,
    /// Qn 元素 → { (Qn, 伙伴分解) → 矩 }（联合）。
    /// 必须与边际分开累加：和的方差不等于方差的和
    qn_partner: HashMap<String, BTreeMap<BridgeKey, Moments>>,
    /// 元素 → { 配位数 → 矩 }
    cn: HashMap<String, BTreeMap<u32, Moments>>,
    /// 配体类型 × 伙伴对 → 矩（键含元素与 rank，天然有序）
    oxy: BTreeMap<OxyKey, Moments>,
    /// 配体元素 → { 类型标签 → 矩 }。**独立累加**,不是把 oxy 的行相加 ——
    /// 相关项之和的方差不等于方差之和,`O_b` 的 sd 只能逐帧重算
    ligand: HashMap<String, BTreeMap<LigandKey, Moments>>,
    /// 桥联组合 → 矩
    linkage: BTreeMap<LinkKey, Moments>,
    /// 元素 → (Σ值·计数, Σ计数)，用于池化均值
    qn_sum: HashMap<String, (f64, usize)>,
    cn_sum: HashMap<String, (f64, usize)>,
    /// 桥氧个数的池化均值。与 qn_sum 是两个量：qn 只数同元素连接，
    /// n_bo 数桥氧个数，三簇氧下 n_bo < qn + Σm
    n_bo_sum: HashMap<String, (f64, usize)>,
}

/// Folds one frame's histogram into `dst`, converting counts to fractions with
/// `total` as the denominator, and updating the pooled sum via `value_of`.
fn push_hist<K: Ord + Clone>(
    dst: &mut HashMap<String, BTreeMap<K, Moments>>,
    sums: &mut HashMap<String, (f64, usize)>,
    elem: &str,
    hist: &HashMap<K, usize>,
    total: usize,
    value_of: impl Fn(&K) -> f64,
) {
    if total == 0 { return; }
    let slot = dst.entry(elem.to_string()).or_default();
    let sum = sums.entry(elem.to_string()).or_insert((0.0, 0));
    for (k, &c) in hist {
        slot.entry(k.clone()).or_default().push(c, c as f64 / total as f64);
        sum.0 += value_of(k) * c as f64;
        sum.1 += c;
    }
}

fn merge_keyed<K: Ord>(dst: &mut BTreeMap<K, Moments>, src: BTreeMap<K, Moments>) {
    for (k, mo) in src { dst.entry(k).or_default().merge(&mo); }
}

impl Accumulator {
    fn push(&mut self, fd: &FrameData) {
        for (elem, (hist, total)) in &fd.qn {
            push_hist(&mut self.qn, &mut self.qn_sum, elem, hist, *total,
                      |n| *n as f64);
        }
        // 池化均值只从边际取一次,故这里传一个丢弃用的累加器
        let mut ignored = HashMap::new();
        for (elem, (hist, total)) in &fd.qn_partner {
            push_hist(&mut self.qn_partner, &mut ignored, elem, hist, *total, |_| 0.0);
        }
        for (elem, (hist, total)) in &fd.cn {
            push_hist(&mut self.cn, &mut self.cn_sum, elem, hist, *total, |v| *v as f64);
        }
        for (elem, (hist, total)) in &fd.ligand {
            push_hist(&mut self.ligand, &mut ignored, elem, hist, *total, |_| 0.0);
        }
        // oxy 的分母是该配体元素的原子数,与 ligand 一致
        for (key, &c) in &fd.oxy {
            let Some((_, total)) = fd.ligand.get(&key.0) else { continue };
            if *total == 0 { continue; }
            self.oxy.entry(key.clone()).or_default().push(c, c as f64 / *total as f64);
        }
        for (elem, (hist, _)) in &fd.n_bo {
            let e = self.n_bo_sum.entry(elem.clone()).or_insert((0.0, 0));
            for (&v, &c) in hist {
                e.0 += v as f64 * c as f64;
                e.1 += c;
            }
        }
        if fd.n_links > 0 {
            for (key, &c) in &fd.linkage {
                self.linkage.entry(key.clone()).or_default()
                    .push(c, c as f64 / fd.n_links as f64);
            }
        }
    }

    fn merge(&mut self, other: Self) {
        for (k, inner) in other.qn {
            merge_keyed(self.qn.entry(k).or_default(), inner);
        }
        for (k, inner) in other.qn_partner {
            merge_keyed(self.qn_partner.entry(k).or_default(), inner);
        }
        for (k, inner) in other.cn {
            merge_keyed(self.cn.entry(k).or_default(), inner);
        }
        merge_keyed(&mut self.oxy, other.oxy);
        for (k, inner) in other.ligand {
            merge_keyed(self.ligand.entry(k).or_default(), inner);
        }
        merge_keyed(&mut self.linkage, other.linkage);
        for (dst, src) in [(&mut self.qn_sum, other.qn_sum),
                           (&mut self.cn_sum, other.cn_sum),
                           (&mut self.n_bo_sum, other.n_bo_sum)] {
            for (k, (s, n)) in src {
                let e = dst.entry(k).or_insert((0.0, 0));
                e.0 += s;
                e.1 += n;
            }
        }
    }

    fn finalize(self, n_frames: usize, n_atoms: usize, params: TypeParams) -> NetworkResult {
        let to_dist = |m: &BTreeMap<u32, Moments>| -> Vec<(u32, Bin)> {
            m.iter().map(|(&v, mo)| (v, mo.finish(n_frames))).collect()
        };
        let mean_of = |sums: &HashMap<String, (f64, usize)>, k: &String| -> f64 {
            match sums.get(k) {
                Some(&(s, n)) if n > 0 => s / n as f64,
                _ => f64::NAN,
            }
        };

        // 只保留轨迹里真正出现过的元素。参数里点名但不存在的元素若留下,
        // 均值会算成 0.0,把「不存在」伪装成「平均值为零」;长表下缺席就是不出行
        let qn_dist: HashMap<_, _> = self.qn.iter()
            .map(|(e, m)| (e.clone(), to_dist(m))).collect();
        let qn_partner_dist: HashMap<_, _> = self.qn_partner.iter()
            .map(|(e, m)| {
                let rows = m.iter().map(|(k, mo)| (k.clone(), mo.finish(n_frames))).collect();
                (e.clone(), rows)
            })
            .collect();
        let mean_qn: HashMap<_, _> = self.qn.keys()
            .map(|e| (e.clone(), mean_of(&self.qn_sum, e))).collect();
        let cn_dist: HashMap<_, _> = self.cn.iter()
            .map(|(e, m)| (e.clone(), to_dist(m))).collect();
        let mean_cn: HashMap<_, _> = self.cn.keys()
            .map(|e| (e.clone(), mean_of(&self.cn_sum, e))).collect();
        // 键取自 n_bo_sum 而非 qn:非 Qn 形成子(Al)也有桥氧数
        let mean_n_bo: HashMap<_, _> = self.n_bo_sum.keys()
            .map(|e| (e.clone(), mean_of(&self.n_bo_sum, e))).collect();

        let ligand_dist: HashMap<_, _> = self.ligand.iter()
            .map(|(e, m)| {
                let rows: Vec<(String, Bin)> = m.iter()
                    .map(|((_, lbl), mo)| (lbl.clone(), mo.finish(n_frames)))
                    .collect();
                (e.clone(), rows)
            })
            .collect();

        let oxy_dist: Vec<(String, Vec<String>, Bin)> = self.oxy.iter()
            .map(|((_, _, lbl, partners), mo)| {
                (lbl.clone(), partners.clone(), mo.finish(n_frames))
            })
            .collect();

        let linkage: Vec<(LinkKey, Bin)> = self.linkage.iter()
            .map(|(k, mo)| (k.clone(), mo.finish(n_frames)))
            .collect();

        NetworkResult {
            qn_dist, qn_partner_dist, mean_qn, cn_dist, mean_cn, mean_n_bo, oxy_dist, ligand_dist,
            linkage, n_frames, n_atoms, params,
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
        TypeParams::new(c, BTreeMap::new())
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

        let bp = &res.qn_dist["P"];
        assert_eq!(bp.len(), 1);
        assert_eq!(bp[0].0, 1);                       // 桥接数 = 1
        assert_eq!(bp[0].1.count, 2);                 // 两个 P
        // 联合分布另成一表。该桥通向 P，是**同元素**连接，已计入 n；
        // m 只装异核伙伴，故这里为空（P 自己没有 m_P 列）
        let pp = &res.qn_partner_dist["P"];
        assert_eq!(pp[0].0.0, 1);
        assert!(pp[0].0.1.is_empty(), "同元素连接进 n，不再重复成 m_P 列: {:?}", pp[0].0.1);

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
        assert_eq!(res.qn_dist["P"][0].0, 0);
        assert!(res.qn_partner_dist["P"][0].0.1.is_empty(), "Q0 没有任何桥");
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
        let params = TypeParams::new(cutoffs, BTreeMap::new());
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
        let params = TypeParams::new(cutoffs, modifier_cutoffs);
        let res = calc_network(&traj_of(vec![frame]), &params).unwrap();

        assert_eq!(res.cn_dist["Zn"][0].0, 1, "Zn 配位数 1");
        assert_eq!(res.cn_dist["Zn"][0].1.count, 1);
        assert!(!res.qn_dist.contains_key("Zn"), "修饰子没有桥接数");
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
        let params = TypeParams::new(cutoffs, modifier_cutoffs);
        let res = calc_network(&traj_of(vec![frame]), &params).unwrap();

        for absent in ["Al", "Zn"] {
            assert!(!res.qn_dist.contains_key(absent));
            assert!(!res.cn_dist.contains_key(absent));
            assert!(!res.mean_qn.contains_key(absent));
            assert!(!res.mean_cn.contains_key(absent));
        }
        let tables = res.to_tables();
        let comp = &tables.iter().find(|(n, _)| n == "composition").unwrap().1;
        let c = comp.column("element").unwrap();
        let elems: Vec<String> = (0..c.len()).map(|i| c.cell(i)).collect();
        assert!(!elems.contains(&"Al".to_string()) && !elems.contains(&"Zn".to_string()),
                "参数里点名但轨迹里没有的元素不该出行: {elems:?}");
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

    /// 桥联表：两端各带元素 / 桥接数 / 配位数，且只存规范半边。
    #[test]
    fn test_linkage_records_both_ends_once() {
        // Al–O–P 一座桥，两端各带一个非桥氧
        let frame = cube(vec![
            atom("Al", 0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0), // 桥
            atom("P",  3.2, 0.0, 0.0),
            atom("O",  0.0, 1.6, 0.0), // Al 的非桥氧
            atom("O",  3.2, 1.6, 0.0), // P 的非桥氧
        ]);
        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("P".to_string(), "O".to_string()), 2.3);
        cutoffs.insert(("Al".to_string(), "O".to_string()), 2.3);
        let params = TypeParams::new(cutoffs, BTreeMap::new());
        let res = calc_network(&traj_of(vec![frame]), &params).unwrap();

        assert_eq!(res.linkage.len(), 1, "一座桥只存一行: {:?}", res.linkage);
        let ((a, b, ligand, n_formers), bin) = &res.linkage[0];
        assert_eq!(*n_formers, 2, "普通桥氧连两个形成子");
        assert_eq!(ligand, "O", "桥中间那个原子的元素进键");
        assert_eq!(bin.count, 1);
        // 规范半边：Al < P，故 Al 在前。两端都带 (元素, 同元素连接数, 配位数)。
        // 这座桥是 Al-O-P：两端各自的同元素连接数都是 0（没有 Al-O-Al，
        // 也没有 P-O-P），配位数都是 2
        assert_eq!(*a, ("Al".to_string(), 0, 2));
        assert_eq!(*b, ("P".to_string(), 0, 2));
    }

    /// 三配位氧不是一对，展开成 C(3,2)=3 条记录，靠 `n_formers` 列区分。
    /// 参考 Python 把它们整体丢进一个计数器，Al 连的是谁就此丢失。
    #[test]
    fn test_tricluster_expands_into_three_pairs() {
        // 一个 O 同时被 3 个 P 包围（正三角形中心）
        let r = 1.8;
        let frame = cube(vec![
            atom("O", 0.0, 0.0, 0.0),
            atom("P", r, 0.0, 0.0),
            atom("P", -r / 2.0,  r * 0.8660254, 0.0),
            atom("P", -r / 2.0, -r * 0.8660254, 0.0),
        ]);
        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("P".to_string(), "O".to_string()), 2.3);
        let params = TypeParams::new(cutoffs, BTreeMap::new());
        let res = calc_network(&traj_of(vec![frame]), &params).unwrap();

        assert_eq!(res.oxy_dist.len(), 1);
        assert_eq!(res.oxy_dist[0].0, "O_t");
        // 三个 Al 位点状态完全相同 → 三条记录合成一行,计数 3
        let total: usize = res.linkage.iter().map(|(_, b)| b.count).sum();
        assert_eq!(total, 3, "C(3,2) = 3 对: {:?}", res.linkage);
        assert!(res.linkage.iter().all(|((_, _, _, n), _)| *n == 3),
                "n_formers 必须标出这是三配位配体");

        // 三簇配体贡献的是**连接**：该 P 通过这一个氧连上了 2 个 P，故 n = 2,
        // 而桥氧只有 1 个。这正是 n 与 n_bo 分叉的地方 —— 传统 Qn 分析假设
        // 「无三键氧」,分叉一出现该假设就不成立了
        let bp = &res.qn_dist["P"];
        assert_eq!(bp.len(), 1);
        assert_eq!(bp[0].0, 2, "三簇氧把该 P 连上了 2 个 P，n = 2");
        assert_eq!(res.mean_n_bo["P"], 1.0, "但桥氧只有 1 个：一个氧、两个连接");
        let (key, _) = &res.qn_partner_dist["P"][0];
        assert!(key.1.is_empty(), "全是同元素连接，无异核伙伴: {key:?}");
    }

    /// 简单 Qn 分布必须**直接可读**，而不是要靠 groupby 从伙伴分解里聚出来。
    /// 0.2.1 曾把两者合成一张表（理由是「groupby 能退回去」），结果是打开
    /// 文件根本看不到 P 的 Qn 占比 —— 它散在 14 行里。
    #[test]
    fn test_qn_table_is_the_plain_distribution() {
        // 两个 P 各连 1 个桥氧，一个通向 P、一个通向 Al
        let frame = cube(vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0),
            atom("Al", 3.2, 0.0, 0.0),
            atom("P",  0.0, 5.0, 0.0),
            atom("O",  1.6, 5.0, 0.0),
            atom("P",  3.2, 5.0, 0.0),
        ]);
        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("P".to_string(), "O".to_string()), 2.3);
        cutoffs.insert(("Al".to_string(), "O".to_string()), 2.3);
        let params = TypeParams::new(cutoffs, BTreeMap::new());
        let res = calc_network(&traj_of(vec![frame]), &params).unwrap();

        // 边际：三个 P 各有 1 个桥氧，但 n 只数同元素连接 —— 连向 Al 的那个
        // P 是 Q0，互相连的两个 P 是 Q1。这正是新旧口径的分界：旧口径三个
        // 都是「1 个桥」合成一行，文献口径分成两行
        let bp = &res.qn_dist["P"];
        assert_eq!(bp.len(), 2, "n=0 与 n=1 各一行: {bp:?}");
        assert_eq!((bp[0].0, bp[0].1.count), (0, 1), "连向 Al 的那个 P 是 Q0");
        assert_eq!((bp[1].0, bp[1].1.count), (1, 2), "互连的两个 P 是 Q1");

        // 联合：Q0 那个带 m_Al=1，两个 Q1 无异核伙伴
        let pp = &res.qn_partner_dist["P"];
        assert_eq!(pp.len(), 2, "联合分布按伙伴拆开: {pp:?}");
        assert_eq!(pp.iter().map(|(_, b)| b.count).sum::<usize>(), 3);
        assert_eq!(pp[0].0.1.get("Al"), Some(&1), "Q0 的那个 P 连着一个 Al: {pp:?}");

        // 两张表在 count 上必须闭合
        let tables = res.to_tables();
        let names: Vec<&str> = tables.iter().map(|(n, _)| n.as_str()).collect();
        assert!(names.contains(&"qn") && names.contains(&"qn_partner"),
                "两张表都要出: {names:?}");
    }

    // ─── Qn 元素的收窄 ────────────────────────────────────────────────────────

    /// Al 由配位数刻画，不由 Qn。它必须从 qn / qn_partner 两张表里**整体消失**，
    /// 同时保留在 coordination / ligand_type / linkage —— 否则 Al-O-Al 无从统计。
    #[test]
    fn test_non_qn_former_leaves_the_qn_tables_but_stays_elsewhere() {
        let frame = cube(vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0),
            atom("Al", 3.2, 0.0, 0.0),
        ]);
        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("P".to_string(), "O".to_string()), 2.3);
        cutoffs.insert(("Al".to_string(), "O".to_string()), 2.3);
        let params = TypeParams::new(cutoffs, BTreeMap::new());
        let res = calc_network(&traj_of(vec![frame]), &params).unwrap();

        assert!(!res.qn_dist.contains_key("Al"), "Al 不该有 Qn 分布");
        assert!(!res.qn_partner_dist.contains_key("Al"));
        assert!(!res.mean_qn.contains_key("Al"));
        assert!(res.qn_dist.contains_key("P"));
        // 但配位数、桥联、氧分类照旧含 Al
        assert!(res.cn_dist.contains_key("Al"));
        assert!(res.mean_cn.contains_key("Al"));
        assert_eq!(res.linkage.len(), 1, "P-O-Al 这座桥必须还在");
        assert_eq!(res.oxy_dist[0].1, vec!["Al".to_string(), "P".to_string()]);

        let tables = res.to_tables();
        let col = |name: &str, c: &str| -> Vec<String> {
            let t = &tables.iter().find(|(n, _)| n == name).unwrap().1;
            let c = t.column(c).unwrap();
            (0..c.len()).map(|i| c.cell(i)).collect()
        };
        assert!(!col("qn", "former").contains(&"Al".to_string()));
        assert!(col("coordination", "element").contains(&"Al".to_string()));

        // composition:P 走单元词汇,Al 走原子词汇 Al_1,两者并存一表。
        // 这个 P 的唯一桥氧通向 Al,是异核连接,故 n=0 → P-Q0
        let labels = col("composition", "label");
        assert!(labels.contains(&"P-Q0".to_string()), "{labels:?}");
        assert!(labels.contains(&"Al_1".to_string()), "{labels:?}");
        // 每个元素只出现一种刻画:P 不该再有配位数行
        assert_eq!(labels.iter().filter(|l| l.starts_with("P")).count(), 1,
                   "P 只应有 Qn 行,不该按配位数再来一遍: {labels:?}");

        // ligand_type 的 label 合并成一行可读的连接式
        assert_eq!(col("ligand_type", "label"), vec!["Al-O_b-P".to_string()]);
        assert_eq!(col("ligand_type", "type"), vec!["O_b".to_string()]);
    }

    /// 没有任何 Qn 形成子时，两张 Qn 表**整个不出**。只有表头的文件读起来
    /// 像「测了，Qn 是零」，而实情是压根没测。
    #[test]
    fn test_qn_tables_are_omitted_when_no_former_has_a_qn() {
        let frame = cube(vec![
            atom("Al", 0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0),
            atom("Al", 3.2, 0.0, 0.0),
        ]);
        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("Al".to_string(), "O".to_string()), 2.3);
        let res = calc_network(&traj_of(vec![frame]),
                               &TypeParams::new(cutoffs, BTreeMap::new())).unwrap();

        let tables = res.to_tables();
        let names: Vec<&str> = tables.iter().map(|(n, _)| n.as_str()).collect();
        assert!(!names.contains(&"qn"), "不该出空的 qn 表: {names:?}");
        assert!(!names.contains(&"qn_partner"));
        assert!(names.contains(&"coordination") && names.contains(&"linkage"));
    }

    /// `--qn` 覆盖默认列表：同一构型,Al 进 Qn 表而 P 退出。
    #[test]
    fn test_qn_override_moves_which_element_gets_a_qn_table() {
        let frame = cube(vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0),
            atom("Al", 3.2, 0.0, 0.0),
        ]);
        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("P".to_string(), "O".to_string()), 2.3);
        cutoffs.insert(("Al".to_string(), "O".to_string()), 2.3);
        let params = TypeParams::new(cutoffs, BTreeMap::new())
            .with_qn_elements(["Al".to_string()]);
        let res = calc_network(&traj_of(vec![frame]), &params).unwrap();

        assert!(res.qn_dist.contains_key("Al"));
        assert!(!res.qn_dist.contains_key("P"));
    }

    /// 配体元素进 linkage 的键：Al-O-P 与 Al-F-P 是两种桥，合并会报出一个
    /// 实验无法复现的数。这份体系只有 O，所以这条只能靠构造双配体来验。
    #[test]
    fn test_two_ligand_species_do_not_share_a_linkage_row() {
        let frame = cube(vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0),
            atom("Al", 3.2, 0.0, 0.0),
            atom("P",  0.0, 6.0, 0.0),
            atom("F",  1.6, 6.0, 0.0),
            atom("Al", 3.2, 6.0, 0.0),
        ]);
        let mut cutoffs = BTreeMap::new();
        for (f, l) in [("P", "O"), ("Al", "O"), ("P", "F"), ("Al", "F")] {
            cutoffs.insert((f.to_string(), l.to_string()), 2.3);
        }
        let res = calc_network(&traj_of(vec![frame]),
                               &TypeParams::new(cutoffs, BTreeMap::new())).unwrap();

        assert_eq!(res.linkage.len(), 2, "两种配体两行: {:?}", res.linkage);
        let mut ligands: Vec<&str> = res.linkage.iter()
            .map(|((_, _, l, _), _)| l.as_str()).collect();
        ligands.sort_unstable();
        assert_eq!(ligands, vec!["F", "O"]);

        // 展示列把配体摆在中间,且两端数字按各自约定:P 报 Qn,Al 报配位数
        let tables = res.to_tables();
        let t = &tables.iter().find(|(n, _)| n == "linkage").unwrap().1;
        let c = t.column("linkage").unwrap();
        let mut rows: Vec<String> = (0..c.len()).map(|i| c.cell(i)).collect();
        rows.sort();
        // P 的数字是同元素连接数：这两座桥都通向 Al，故 P 是 Q0
        assert_eq!(rows, vec!["Al_1-F-P_0".to_string(), "Al_1-O-P_0".to_string()]);
    }

    // ─── composition 表与配体分母 ──────────────────────────────────────────────

    /// composition 表的核心恒等式：**每个元素的 fraction 求和为 1**。
    /// 分母恒为该元素的原子数——Q2 占全部 P，Al_4 占全部 Al，O_b 占全部 O。
    #[test]
    fn test_composition_fractions_sum_to_one_within_each_element() {
        let frame = cube(vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0),
            atom("Al", 3.2, 0.0, 0.0),
            atom("O",  0.0, 1.6, 0.0),   // 只连 P 的非桥氧
            atom("O",  6.0, 6.0, 6.0),   // 游离氧
            atom("Zn", 9.0, 9.0, 9.0),
        ]);
        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("P".to_string(), "O".to_string()), 2.3);
        cutoffs.insert(("Al".to_string(), "O".to_string()), 2.3);
        let mut modifier_cutoffs = BTreeMap::new();
        modifier_cutoffs.insert(("Zn".to_string(), "O".to_string()), 2.6);
        let res = calc_network(&traj_of(vec![frame]),
                               &TypeParams::new(cutoffs, modifier_cutoffs)).unwrap();

        let tables = res.to_tables();
        let t = &tables.iter().find(|(n, _)| n == "composition").unwrap().1;
        let (e, f) = (t.column("element").unwrap(), t.column("fraction").unwrap());
        let mut sums: std::collections::HashMap<String, f64> = HashMap::new();
        for i in 0..e.len() {
            *sums.entry(e.cell(i)).or_insert(0.0) += f.cell(i).parse::<f64>().unwrap();
        }
        assert!(sums.contains_key("P") && sums.contains_key("Al")
                    && sums.contains_key("O") && sums.contains_key("Zn"),
                "形成子、配体、修饰子都要在表里: {sums:?}");
        // 容差取 CSV 自己的精度:`Column::cell` 按 {:.6e} 渲染,1/3 写成 3.333333e-1,
        // 三个相加就是 0.9999999。恒等式在 f64 上是精确的,下面顺带验一遍
        for (elem, sum) in &sums {
            assert!((sum - 1.0).abs() < 1e-6, "{elem} 的 fraction 求和为 {sum}, 应为 1");
        }
        for (elem, rows) in &res.ligand_dist {
            let s: f64 = rows.iter().map(|(_, b)| b.fraction).sum();
            assert!((s - 1.0).abs() < 1e-12, "{elem} 在 f64 上必须精确求和为 1, 实得 {s}");
        }
    }

    /// 配体 `fraction` 的分母是**该配体元素**的原子数，不是全体配体原子。
    /// 「桥氧占全部氧的比例」是文献的 BO fraction；「占 O+F 的比例」不对应任何常用量。
    #[test]
    fn test_ligand_fraction_denominator_is_per_element() {
        // 2 个 O（1 桥 1 非桥）+ 2 个 F（都非桥）。若分母取全体配体，
        // O_b 会是 1/4 = 0.25；正确答案是 1/2 = 0.5
        let frame = cube(vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0),
            atom("Al", 3.2, 0.0, 0.0),
            atom("O",  0.0, 1.6, 0.0),
            atom("F", -1.6, 0.0, 0.0),
            atom("F",  0.0, -1.6, 0.0),
        ]);
        let mut cutoffs = BTreeMap::new();
        for (fm, l) in [("P", "O"), ("Al", "O"), ("P", "F"), ("Al", "F")] {
            cutoffs.insert((fm.to_string(), l.to_string()), 2.3);
        }
        let res = calc_network(&traj_of(vec![frame]),
                               &TypeParams::new(cutoffs, BTreeMap::new())).unwrap();

        let ob = res.oxy_dist.iter().find(|(l, _, _)| l == "O_b").unwrap();
        assert!((ob.2.fraction - 0.5).abs() < 1e-12,
                "O_b 应占全部 O 的 1/2, 实得 {} (取全体配体分母会给 0.25)", ob.2.fraction);
        let fn_ = res.oxy_dist.iter().find(|(l, _, _)| l == "F_n").unwrap();
        assert!((fn_.2.fraction - 1.0).abs() < 1e-12,
                "两个 F 都是非桥, 应为 1.0, 实得 {}", fn_.2.fraction);
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
