//! Network atom type classification: per-frame labeling for glass network analysis.
//!
//! Shared by `ferro-structure` (single-frame typing / cluster finding) and
//! `ferro-analysis` (trajectory statistics) without creating cross-crate dependencies.
//!
//! # What classification produces
//!
//! [`classify_frame`] returns one [`AtomType`] per atom — a *structured* value that
//! keeps the numbers it computed (bridging count, coordination number, the partner
//! elements of a ligand).  Rendering a type as text happens in exactly one place,
//! [`AtomType::label`]; nothing downstream parses a label back into a number.
//!
//! # Label scheme
//!
//! Every label is `<element>_<suffix>`, split at the **first** underscore — the same
//! convention the LAMMPS dump reader uses, so a labelled trajectory round-trips and
//! its labels can be selected with `-x/-y/-z`.
//!
//! | Atom role | Label |
//! |-----------|-------|
//! | Qn former | `P_0`, `P_3`, `Si_4`, … (digit = Qn, i.e. bridging ligands) |
//! | Other former | `Al_4`, `Al_5`, … (digit = **coordination number**) |
//! | Free ligand (0 NF) | `O_f` |
//! | Non-bridging ligand (1 NF) | `O_n` |
//! | Bridging ligand (2 NF) | `O_b` |
//! | Tricluster ligand (≥3 NF) | `O_t` |
//! | Modifier | element symbol unchanged (`Zn`) |
//! | Other atoms | element symbol unchanged |
//!
//! **A former's digit is not always the same quantity.**  Qn is a tetrahedral-former
//! convention, so it is shown only for the elements in [`crate::data::qn_elements`];
//! every other former shows its coordination number, which is what the literature
//! quotes for it (`Al[4]`/`Al[5]`/`Al[6]`).  The two coincide for a former with no
//! non-bridging ligand — true of Al in aluminophosphates — so a system where they
//! agree proves nothing about which one is being displayed.  See
//! [`AtomType::site_digit`].
//!
//! Ligand labels carry no partner suffix: which formers a bridge joins is a
//! *statistic*, reported through [`AtomType::Ligand::partners`] as data columns,
//! not encoded into the label text.
//!
//! Modifiers carry no role suffix either.  The former `Zn_f`/`Zn_t`/`Zn_b`/`X`
//! scheme binned by non-bridging-ligand count 0/1/2/≥3, which for a modifier of
//! ordinary coordination puts almost everything in the ≥3 catch-all (97 % of Zn in
//! the reference trajectory), so it carried no resolving power.  Modifiers are now
//! described by their coordination number instead.

use crate::data::qn_elements;
use crate::{Cell, Frame};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

/// Cutoff table: `(element_A, element_B)` → distance [Å]
pub type CutoffTable = BTreeMap<(String, String), f64>;

/// Parameters for network atom type classification.
#[derive(Debug, Clone)]
pub struct TypeParams {
    /// `(former_elem, ligand_elem)` → cutoff [Å]
    pub cutoffs: CutoffTable,
    /// `(modifier_elem, ligand_elem)` → cutoff [Å].  Empty → no modifier classification.
    pub modifier_cutoffs: CutoffTable,
    /// Formers described by a Qn speciation; everything else is described by its
    /// coordination number.  See [`crate::data::qn_elements`].
    ///
    /// Held here rather than looked up from the static table at render time because
    /// `--qn` overrides it per run: a classification must carry the convention it was
    /// made under, or a label means different things in two files of the same batch.
    pub qn_elements: BTreeSet<String>,
}

impl TypeParams {
    /// Cutoff tables with the default Qn element set.
    pub fn new(cutoffs: CutoffTable, modifier_cutoffs: CutoffTable) -> Self {
        Self { cutoffs, modifier_cutoffs, qn_elements: qn_elements::default_qn_set() }
    }

    /// Replace the Qn element set (`ferro net --qn`).
    pub fn with_qn_elements(mut self, elems: impl IntoIterator<Item = String>) -> Self {
        self.qn_elements = elems.into_iter().collect();
        self
    }

    /// Whether this former is reported with a Qn speciation.
    pub fn is_qn_former(&self, elem: &str) -> bool {
        self.qn_elements.contains(elem)
    }

    /// Former elements reported with a Qn speciation, sorted — the row set of the
    /// `qn` and `qn_partner` tables.  May be empty.
    pub fn qn_formers(&self) -> Vec<String> {
        self.formers().into_iter().filter(|e| self.is_qn_former(e)).collect()
    }

    /// Unique former elements, sorted.
    pub fn formers(&self) -> Vec<String> {
        unique_keys_left(&self.cutoffs)
    }

    /// Unique ligand elements (from former cutoffs), sorted.
    pub fn ligands(&self) -> Vec<String> {
        unique_keys_right(&self.cutoffs)
    }

    /// Unique modifier elements, sorted.
    pub fn modifiers(&self) -> Vec<String> {
        unique_keys_left(&self.modifier_cutoffs)
    }

    /// Former–ligand cutoff, if defined.
    pub fn cutoff(&self, former: &str, ligand: &str) -> Option<f64> {
        self.cutoffs.get(&(former.to_string(), ligand.to_string())).copied()
    }

    /// All `(former_elem, cutoff)` pairs for a given ligand element.
    pub fn formers_for_ligand(&self, ligand: &str) -> Vec<(String, f64)> {
        self.cutoffs.iter()
            .filter(|((_, l), _)| l == ligand)
            .map(|((f, _), &c)| (f.clone(), c))
            .collect()
    }
}

// ─── Structured classification result ─────────────────────────────────────────

/// Structural role of one atom, with the quantities the classification computed.
///
/// Consumers read the fields directly.  [`label`](Self::label) is the single point
/// where a type becomes text — changing the label scheme touches that method only.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum AtomType {
    /// Network former.
    ///
    /// | field | meaning |
    /// |---|---|
    /// | `qn` | `Some(n)` for a Qn former, `None` otherwise — the flag *and* the value |
    /// | `n_bo` | bridging ligands: ligands shared with at least one other former |
    /// | `cn` | **all** ligands within cutoff, bridging or not |
    /// | `bridges_to` | connections to each partner element (a tricluster counts twice) |
    ///
    /// `qn` is the **homopolar** connection count — P–O–P for a P, i.e.
    /// `bridges_to[elem]`.  This is the literature's `n` in `Q^n_m` / `P^n_(mAl,xB)`;
    /// the heteropolar connections are the `m` values and live in `bridges_to` under
    /// their own elements.  The total bridging-oxygen count is `n_bo`, which is a
    /// *third* number: it equals `n + Σm` only when no tricluster is present.
    ///
    /// `n_bo` and `cn` are different quantities and must not be conflated: `cn`
    /// counts every ligand inside the cutoff, `n_bo` only those a second former
    /// also touches.  They happen to coincide for a former carrying no non-bridging
    /// ligand, which is true of Al in aluminophosphate glasses and false in general.
    ///
    /// `bridges_to` counts **connections**, not bridging ligands: a ligand shared
    /// with three formers connects this one to two partners and contributes two.
    /// So `Σ bridges_to ≥ bridging`, the excess being triclusters.  Counting
    /// bridging ligands instead would leave `Σ bridges_to` short of the literature's
    /// `n + Σm` identity, which is what `Q^n_m` / `P^n_(mAl,xB)` add up.
    Former {
        elem: String,
        qn: Option<u32>,
        n_bo: u32,
        cn: u32,
        bridges_to: BTreeMap<String, u32>,
    },
    /// Ligand (oxygen).  `partners` = elements of the formers bonded to it, sorted;
    /// its length is the classification (0 free, 1 non-bridging, 2 bridging,
    /// ≥3 tricluster).
    Ligand { elem: String, partners: Vec<String> },
    /// Network modifier: counted for coordination, but excluded from the bridging
    /// count and from ligand classification.  `cn` = ligands within cutoff.
    Modifier { elem: String, cn: u32 },
    /// Any atom none of the cutoff tables mentions.
    Other { elem: String },
}

impl AtomType {
    /// The chemical element, whatever the role.
    pub fn element(&self) -> &str {
        match self {
            AtomType::Former { elem, .. }
            | AtomType::Ligand { elem, .. }
            | AtomType::Modifier { elem, .. }
            | AtomType::Other { elem } => elem,
        }
    }

    /// The digit a former's label carries: its Qn if it has one, its coordination
    /// number otherwise.
    ///
    /// Two different quantities behind one slot, which is the literature's own
    /// convention — `Q²` and `Al[4]` are read by different rules and nobody confuses
    /// them.  The alternative, showing the bridging count for every former, puts a
    /// number in the Al slot that no aluminophosphate paper quotes and that happens to
    /// equal the coordination number whenever Al carries no non-bridging ligand — a
    /// coincidence of the system, not a definition, so the error would be invisible.
    pub fn site_digit(&self) -> Option<u32> {
        match self {
            AtomType::Former { qn: Some(n), .. } => Some(*n),
            AtomType::Former { qn: None, cn, .. } => Some(*cn),
            _ => None,
        }
    }

    /// Render as a site label.  **The only place a type becomes text.**
    pub fn label(&self) -> String {
        match self {
            AtomType::Former { elem, qn, cn, .. } => {
                format!("{elem}_{}", qn.unwrap_or(*cn))
            }
            AtomType::Ligand { elem, partners } => match partners.len() {
                0 => format!("{elem}_f"),
                1 => format!("{elem}_n"),
                2 => format!("{elem}_b"),
                _ => format!("{elem}_t"),
            },
            AtomType::Modifier { elem, .. } | AtomType::Other { elem } => elem.clone(),
        }
    }

    /// Sort rank *within one distribution table* (all rows share a role).
    ///
    /// Ligand: `_f` < `_n` < `_b` < `_t`.  Former: by the digit its label shows, so
    /// the sort order and the printed order cannot disagree.
    pub fn class_rank(&self) -> u8 {
        match self {
            AtomType::Former { qn, cn, .. } => qn.unwrap_or(*cn).min(u8::MAX as u32) as u8,
            AtomType::Ligand { partners, .. } => partners.len().min(3) as u8,
            AtomType::Modifier { .. } | AtomType::Other { .. } => 0,
        }
    }

    /// Sort rank when every role is printed in *one* table:
    /// formers, then ligands (free → non-bridging → bridging → tricluster),
    /// then modifiers, then everything else.
    pub fn display_rank(&self) -> u8 {
        match self {
            AtomType::Former { .. } => 0,
            AtomType::Ligand { partners, .. } => 1 + partners.len().min(3) as u8,
            AtomType::Modifier { .. } => 5,
            AtomType::Other { .. } => 6,
        }
    }

    /// True for a ligand bridging **exactly two** formers.
    ///
    /// Not a connectivity test.  A tricluster joins three formers and is the most
    /// strongly connected node in the network, yet it is `false` here — use
    /// `partners.len() >= 2` (equivalently [`crate::LigandKind::Bridging`]) when the
    /// question is whether the ligand links formers at all.  Reach for this method
    /// only when a true two-ended bridge is what is meant, such as separating the
    /// pairwise `linkage` rows from the tricluster ones.
    pub fn is_bridging(&self) -> bool {
        matches!(self, AtomType::Ligand { partners, .. } if partners.len() == 2)
    }
}

// ─── Public entry point ───────────────────────────────────────────────────────

/// Classification of one frame, plus the adjacency it was derived from.
///
/// [`AtomType`] describes each atom on its own; statistics about *pairs* of formers
/// (which site bridges to which) need the graph as well, and recomputing the
/// neighbour search to recover it would repeat the most expensive part of the work.
#[derive(Debug, Clone)]
pub struct FrameTypes {
    /// One type per atom, same order as `frame.atoms`.
    pub types: Vec<AtomType>,
    /// Ligand atom index → indices of the former atoms bonded to it, ascending.
    /// Every ligand atom has an entry, including those bonded to no former.
    pub ligand_formers: HashMap<usize, Vec<usize>>,
}

/// Classify every atom in one frame.  Returns one type per atom (same order as `frame.atoms`).
pub fn classify_frame(frame: &Frame, cell: &Cell, params: &TypeParams) -> Vec<AtomType> {
    classify_frame_detailed(frame, cell, params).types
}

/// [`classify_frame`] plus the ligand→former adjacency, for pairwise statistics.
pub fn classify_frame_detailed(
    frame: &Frame,
    cell: &Cell,
    params: &TypeParams,
) -> FrameTypes {
    // 按元素建立索引
    let elem_map = build_elem_map(frame);

    // 1. 配体原子的 NF 邻居表（ligand_idx → Vec<(former_elem, former_idx)>）
    let nf_map = build_nf_map(frame, cell, params, &elem_map);

    // 2. 组装：默认 Other，再依次被配体 / 形成子 / 修饰子覆盖（顺序同旧实现）
    let mut types: Vec<AtomType> = frame.atoms.iter()
        .map(|a| AtomType::Other { elem: a.element.clone() })
        .collect();

    let mut ligand_formers: HashMap<usize, Vec<usize>> = HashMap::with_capacity(nf_map.len());
    for (&idx, nf) in &nf_map {
        let mut partners: Vec<String> = nf.iter().map(|(e, _)| e.clone()).collect();
        partners.sort_unstable();
        types[idx] = AtomType::Ligand { elem: frame.atoms[idx].element.clone(), partners };

        let mut idxs: Vec<usize> = nf.iter().map(|(_, i)| *i).collect();
        idxs.sort_unstable();
        ligand_formers.insert(idx, idxs);
    }
    for (idx, t) in classify_formers(frame, cell, params, &elem_map, &nf_map) {
        types[idx] = t;
    }
    for (idx, t) in classify_modifiers_inner(frame, cell, params, &elem_map) {
        types[idx] = t;
    }
    FrameTypes { types, ligand_formers }
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

fn build_elem_map(frame: &Frame) -> HashMap<String, Vec<usize>> {
    let mut map: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, atom) in frame.atoms.iter().enumerate() {
        map.entry(atom.element.clone()).or_default().push(idx);
    }
    map
}

/// `ligand_idx` → `Vec<(former_elem, former_idx)>` — includes ligand atoms with 0 NF neighbors.
fn build_nf_map(
    frame: &Frame,
    cell: &Cell,
    params: &TypeParams,
    elem_map: &HashMap<String, Vec<usize>>,
) -> HashMap<usize, Vec<(String, usize)>> {
    let mut nf_map: HashMap<usize, Vec<(String, usize)>> = HashMap::new();

    // 预先为所有配体原子建立空条目
    for ligand in params.ligands() {
        if let Some(idxs) = elem_map.get(&ligand) {
            for &idx in idxs { nf_map.entry(idx).or_default(); }
        }
    }

    for ligand in params.ligands() {
        let ligand_idxs = match elem_map.get(&ligand) {
            Some(v) => v.clone(),
            None => continue,
        };
        let former_pairs = params.formers_for_ligand(&ligand);

        for &la_idx in &ligand_idxs {
            let la_pos = frame.atoms[la_idx].position;
            for (former, cutoff) in &former_pairs {
                let c2 = cutoff * cutoff;
                let Some(former_idxs) = elem_map.get(former) else { continue };
                for &fa_idx in former_idxs {
                    if fa_idx == la_idx { continue; }
                    let diff = cell.minimum_image(frame.atoms[fa_idx].position - la_pos)
                        .expect("cell must be non-singular");
                    if diff.norm_squared() < c2 {
                        nf_map.entry(la_idx).or_default().push((former.clone(), fa_idx));
                    }
                }
            }
        }
    }
    nf_map
}

/// Former classification.  One pass yields both the bridging count (Qn) and the
/// total coordination number — the analysis layer used to re-scan the same pairs
/// for `cn` alone.
fn classify_formers(
    frame: &Frame,
    cell: &Cell,
    params: &TypeParams,
    elem_map: &HashMap<String, Vec<usize>>,
    nf_map: &HashMap<usize, Vec<(String, usize)>>,
) -> Vec<(usize, AtomType)> {
    let mut result: Vec<(usize, AtomType)> = Vec::new();

    for former_elem in params.formers() {
        let Some(former_idxs) = elem_map.get(&former_elem) else { continue };

        for &fa_idx in former_idxs {
            let fa_pos = frame.atoms[fa_idx].position;
            let mut n_bo = 0u32;
            let mut cn = 0u32;
            let mut bridges_to: BTreeMap<String, u32> = BTreeMap::new();

            for ligand_elem in params.ligands() {
                let Some(&cutoff) = params.cutoffs.get(&(former_elem.clone(), ligand_elem.clone()))
                else { continue };
                let c2 = cutoff * cutoff;
                let Some(ligand_idxs) = elem_map.get(&ligand_elem) else { continue };

                for &la_idx in ligand_idxs {
                    if la_idx == fa_idx { continue; }
                    let diff = cell.minimum_image(frame.atoms[la_idx].position - fa_pos)
                        .expect("cell must be non-singular");
                    if diff.norm_squared() >= c2 { continue; }
                    cn += 1;

                    // 桥接判断：该配体有 ≥2 个 NF 邻居
                    let nf = nf_map.get(&la_idx).map(|v| v.as_slice()).unwrap_or(&[]);
                    if nf.len() < 2 { continue; }
                    n_bo += 1;

                    // 伙伴分解数的是**连接**，不是桥氧：该配体上除自己之外的
                    // 每个形成子各算一个连接。三簇配体 P-O(-Al)(-Al) 因此给出
                    // m_Al += 2 —— 文献 (Q^n_m / P^n_mAl,xB) 数的正是
                    // 「connections with phosphate」，而一个三簇氧确实把该 P
                    // 连上了两个 Al。按桥氧数记则 Σm 会亏空，无法与 n 相加
                    for (partner, i) in nf {
                        if *i == fa_idx { continue; }
                        *bridges_to.entry(partner.clone()).or_insert(0) += 1;
                    }
                }
            }
            // n 只数**同元素**连接（P-O-P）,异核连接留在 bridges_to 里当 m。
            // 这是文献 Q^n_m 的口径:总桥数由 n 与各 m 相加得到,而不是 n 本身。
            // 约定在分类时定死而非渲染时查表:`--qn` 能覆盖它,一个类型必须携带
            // 自己是在哪套约定下产生的
            let qn = params.is_qn_former(&former_elem)
                .then(|| bridges_to.get(&former_elem).copied().unwrap_or(0));
            result.push((fa_idx, AtomType::Former {
                elem: former_elem.clone(), qn, n_bo, cn, bridges_to,
            }));
        }
    }
    result
}

/// Modifier classification: coordination number only.  A modifier never enters the
/// ligands' NF list, so it neither carries a bridging count nor affects how a ligand
/// is classified — that separation is the whole point of `--modifier`.
fn classify_modifiers_inner(
    frame: &Frame,
    cell: &Cell,
    params: &TypeParams,
    elem_map: &HashMap<String, Vec<usize>>,
) -> Vec<(usize, AtomType)> {
    let mut result: Vec<(usize, AtomType)> = Vec::new();

    for mod_elem in params.modifiers() {
        let Some(mod_idxs) = elem_map.get(&mod_elem) else { continue };

        for &ma_idx in mod_idxs {
            let ma_pos = frame.atoms[ma_idx].position;
            let mut cn = 0u32;

            for ((m, ligand), &cutoff) in &params.modifier_cutoffs {
                if m != &mod_elem { continue; }
                let c2 = cutoff * cutoff;
                let Some(ligand_idxs) = elem_map.get(ligand) else { continue };
                for &la_idx in ligand_idxs {
                    if la_idx == ma_idx { continue; }
                    let diff = cell.minimum_image(frame.atoms[la_idx].position - ma_pos)
                        .expect("cell must be non-singular");
                    if diff.norm_squared() < c2 { cn += 1; }
                }
            }

            result.push((ma_idx, AtomType::Modifier { elem: mod_elem.clone(), cn }));
        }
    }
    result
}

// ─── Utility ─────────────────────────────────────────────────────────────────

fn unique_keys_left(table: &CutoffTable) -> Vec<String> {
    let mut v: Vec<String> = table.keys().map(|(a, _)| a.clone())
        .collect::<HashSet<_>>().into_iter().collect();
    v.sort();
    v
}

fn unique_keys_right(table: &CutoffTable) -> Vec<String> {
    let mut v: Vec<String> = table.keys().map(|(_, b)| b.clone())
        .collect::<HashSet<_>>().into_iter().collect();
    v.sort();
    v
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Qn former: the label digit is the bridging count.
    fn former(elem: &str, n_bo: u32, cn: u32) -> AtomType {
        AtomType::Former {
            elem: elem.to_string(), qn: Some(n_bo), n_bo, cn,
            bridges_to: BTreeMap::new(),
        }
    }

    /// Non-Qn former: the label digit is the coordination number.
    fn cn_former(elem: &str, n_bo: u32, cn: u32) -> AtomType {
        AtomType::Former {
            elem: elem.to_string(), qn: None, n_bo, cn, bridges_to: BTreeMap::new(),
        }
    }

    #[test]
    fn test_every_label_splits_at_first_underscore() {
        // 全部标签必须形如 <元素>_<后缀>（或裸元素），否则 LAMMPS dump 读回来
        // 拆不出正确的 element，`-a P` 就选不中导出的结构
        let cases: Vec<(AtomType, &str, &str)> = vec![
            (former("P", 0, 4), "P_0", "P"),
            (former("P", 3, 4), "P_3", "P"),
            (cn_former("Al", 2, 5), "Al_5", "Al"),
            (AtomType::Ligand { elem: "O".into(), partners: vec![] }, "O_f", "O"),
            (AtomType::Ligand { elem: "O".into(), partners: vec!["P".into()] }, "O_n", "O"),
            (AtomType::Ligand {
                elem: "O".into(), partners: vec!["Al".into(), "P".into()],
            }, "O_b", "O"),
            (AtomType::Ligand {
                elem: "O".into(), partners: vec!["Al".into(), "Al".into(), "P".into()],
            }, "O_t", "O"),
            (AtomType::Modifier { elem: "Zn".into(), cn: 5 }, "Zn", "Zn"),
            (AtomType::Other { elem: "Ar".into() }, "Ar", "Ar"),
        ];
        for (t, expect, elem) in cases {
            let label = t.label();
            assert_eq!(label, expect);
            assert_eq!(label.split('_').next().unwrap(), elem,
                       "{label} 的首段必须是元素符号");
            assert_eq!(t.element(), elem);
        }
    }

    #[test]
    fn test_ligand_label_carries_no_partner_suffix() {
        // 伙伴元素是统计量,不进标签 —— 不同伙伴的桥氧共用一个标签
        let lig = |p: &[&str]| AtomType::Ligand {
            elem: "O".into(),
            partners: p.iter().map(|s| s.to_string()).collect(),
        };
        assert_eq!(lig(&["P", "P"]).label(), lig(&["Al", "P"]).label());
        assert_eq!(lig(&["Al", "Al"]).label(), "O_b");
        // 但 partners 本身仍在类型里,统计层能拿到
        let AtomType::Ligand { partners, .. } = lig(&["Al", "P"]) else { unreachable!() };
        assert_eq!(partners, vec!["Al".to_string(), "P".to_string()]);
    }

    #[test]
    fn test_class_rank_orders_ligands() {
        let lig = |n: usize| AtomType::Ligand {
            elem: "O".into(),
            partners: vec!["P".to_string(); n],
        };
        // _f < _n < _b < _t，且 ≥3 全部落在 _t 一档
        assert!(lig(0).class_rank() < lig(1).class_rank());
        assert!(lig(1).class_rank() < lig(2).class_rank());
        assert!(lig(2).class_rank() < lig(3).class_rank());
        assert_eq!(lig(3).class_rank(), lig(5).class_rank());
        assert_eq!(lig(4).label(), "O_t");
    }

    #[test]
    fn test_display_rank_orders_all_roles() {
        let former = former("P", 2, 4);
        let lig = |n: usize| AtomType::Ligand {
            elem: "O".into(),
            partners: vec!["P".to_string(); n],
        };
        let modif = AtomType::Modifier { elem: "Zn".into(), cn: 4 };
        let other = AtomType::Other { elem: "Ar".into() };

        let ranks = [
            former.display_rank(),
            lig(0).display_rank(), lig(1).display_rank(),
            lig(2).display_rank(), lig(3).display_rank(),
            modif.display_rank(),
            other.display_rank(),
        ];
        assert!(ranks.windows(2).all(|w| w[0] < w[1]), "混排顺序必须严格递增: {ranks:?}");

        // 旧方案里过配位氧与过配位修饰子都叫 "X" 而合并成一行；
        // 现在 O_t 与 Zn 是两个标签,不再塌缩
        assert_ne!(lig(3).label(), modif.label());
    }

    #[test]
    fn test_is_bridging_is_exactly_two_partners() {
        let lig = |n: usize| AtomType::Ligand {
            elem: "O".into(),
            partners: vec!["P".to_string(); n],
        };
        assert!(!lig(0).is_bridging());
        assert!(!lig(1).is_bridging());
        assert!(lig(2).is_bridging());
        // ≥3 渲染为 X，旧实现的 `starts_with("Ob_")` 同样不匹配
        assert!(!lig(3).is_bridging());
    }

    // ─── Qn 与配位数的分叉 ────────────────────────────────────────────────────

    /// A frame in a 20 Å cubic box; positions in Å.
    fn frame_of(atoms: &[(&str, f64, f64, f64)]) -> (Frame, Cell) {
        use crate::Atom;
        use nalgebra::Vector3;
        let cell = Cell::from_lengths_angles(20.0, 20.0, 20.0, 90.0, 90.0, 90.0).unwrap();
        let mut f = Frame::new();
        f.atoms = atoms.iter()
            .map(|(e, x, y, z)| Atom::new(*e, Vector3::new(*x, *y, *z)))
            .collect();
        f.cell = Some(cell.clone());
        f.pbc = [true; 3];
        (f, cell)
    }

    fn params_po_alo() -> TypeParams {
        let mut c = CutoffTable::new();
        c.insert(("P".into(), "O".into()), 2.0);
        c.insert(("Al".into(), "O".into()), 2.0);
        TypeParams::new(c, CutoffTable::new())
    }

    #[test]
    fn test_non_qn_former_labels_by_coordination_not_bridging() {
        // 盲区回归:参考轨迹里的 Al 恰好不带非桥氧,于是 bridging == cn,
        // 两种口径的分叉从未被触发过。这里造一个 Al 带非桥氧的构型把它钉死。
        //
        //   O0 —— 同时挨着 Al 和 P  → 桥氧,计入 Al 的 bridging 与 cn
        //   O1 —— 只挨着 Al        → 非桥氧,只计入 Al 的 cn
        let (frame, cell) = frame_of(&[
            ("Al", 0.0, 0.0, 0.0),
            ("O",  1.5, 0.0, 0.0),   // Al-O-P 桥
            ("P",  3.0, 0.0, 0.0),
            ("O",  0.0, 1.5, 0.0),   // 只连 Al 的非桥氧
        ]);
        let types = classify_frame(&frame, &cell, &params_po_alo());

        let AtomType::Former { qn, n_bo, cn, .. } = &types[0] else {
            panic!("Al 应被分类为形成子, got {:?}", types[0])
        };
        assert_eq!(*qn, None, "Al 不是 Qn 元素");
        assert_eq!(*n_bo, 1, "只有一个氧同时连着 P");
        assert_eq!(*cn, 2, "cn 数的是截断内全部氧,含非桥氧");
        assert_ne!(*n_bo, *cn, "构型没造对:这个测试要的正是两者分叉");
        assert_eq!(types[0].label(), "Al_2", "非 Qn 形成子按配位数标注");
        assert_eq!(types[0].site_digit(), Some(2));

        // 同一帧里的 P 走另一套:数字是 Qn。这个 P 唯一的桥氧通向 Al,
        // 是**异核**连接 —— n 只数 P-O-P,故 n=0 而桥氧数 n_bo=1
        let AtomType::Former { qn, n_bo, cn, bridges_to, .. } = &types[2] else {
            panic!("P 应被分类为形成子")
        };
        assert_eq!(*qn, Some(0), "唯一的桥通向 Al,不是 P-O-P");
        assert_eq!((*n_bo, *cn), (1, 1));
        assert_eq!(bridges_to.get("Al"), Some(&1), "那一个连接记在 m_Al 上");
        assert_eq!(types[2].label(), "P_0");

        assert_eq!(types[1].label(), "O_b");
        assert_eq!(types[3].label(), "O_n");
    }

    #[test]
    fn test_qn_override_switches_which_digit_a_label_shows() {
        // `--qn` 改的是约定,而约定必须随分类走 —— 同一构型换个约定,
        // Al 的标签要从配位数变成桥接数
        let (frame, cell) = frame_of(&[
            ("Al", 0.0, 0.0, 0.0),
            ("O",  1.5, 0.0, 0.0),
            ("P",  3.0, 0.0, 0.0),
            ("O",  0.0, 1.5, 0.0),
        ]);
        let params = params_po_alo().with_qn_elements(["Al".to_string()]);
        let types = classify_frame(&frame, &cell, &params);

        // Al 现在报 Qn = 同元素连接数。这个构型里没有 Al-O-Al,故为 0 ——
        // 与它的配位数 2 不同,正好证明约定确实切换了
        assert_eq!(types[0].label(), "Al_0", "Al 报 Qn(同元素连接)=0,而非配位数 2");
        assert_eq!(types[2].label(), "P_1", "P 被移出 Qn 列表,改报配位数 1");
        assert_eq!(params.qn_formers(), vec!["Al".to_string()]);
    }

    #[test]
    fn test_default_qn_set_excludes_al() {
        let p = params_po_alo();
        assert!(p.is_qn_former("P"));
        assert!(!p.is_qn_former("Al"));
        assert_eq!(p.formers(), vec!["Al".to_string(), "P".to_string()]);
        assert_eq!(p.qn_formers(), vec!["P".to_string()]);
    }
}
