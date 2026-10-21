//! Network atom type classification: per-frame labeling for glass network analysis.
//!
//! Used by `ferro net` (statistics and the labelled trajectory) and by the dataset
//! filter's coordination criteria.
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
//! convention, so it is shown only for the elements in [`super::qn_elements`];
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

use super::qn_elements;
use crate::neighbor::NeighborList;
use ferro_core::error::Result;
use ferro_core::{Frame, PeriodicCell};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

/// Cutoff table: `(element_A, element_B)` → distance \[Å\]
pub type CutoffTable = BTreeMap<(String, String), f64>;

/// Parameters for network atom type classification.
#[derive(Debug, Clone)]
pub struct TypeParams {
    /// `(former_elem, ligand_elem)` → cutoff \[Å\]
    pub cutoffs: CutoffTable,
    /// `(modifier_elem, ligand_elem)` → cutoff \[Å\].  Empty → no modifier classification.
    pub modifier_cutoffs: CutoffTable,
    /// Formers described by a Qn speciation; everything else is described by its
    /// coordination number.  See `super::qn_elements`.
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
            AtomType::Former { elem, .. } => {
                format!("{elem}_{}", self.site_digit().expect("every former has a site digit"))
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
            AtomType::Former { .. } => {
                self.site_digit().expect("every former has a site digit").min(u8::MAX as u32) as u8
            }
            AtomType::Ligand { partners, .. } => partners.len().min(3) as u8,
            AtomType::Modifier { .. } | AtomType::Other { .. } => 0,
        }
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
///
/// `Err` when the frame has no 3-D periodic cell, the cell is singular (e.g. a slab
/// written with `c = 0`), or a cutoff exceeds its minimum-image bound.
pub fn classify_frame(frame: &Frame, params: &TypeParams) -> Result<Vec<AtomType>> {
    Ok(classify_frame_detailed(frame, params)?.types)
}

/// [`classify_frame`] plus the ligand→former adjacency, for pairwise statistics.
pub fn classify_frame_detailed(frame: &Frame, params: &TypeParams) -> Result<FrameTypes> {
    let cell = PeriodicCell::of(frame)?;
    // One neighbour list at the longest cutoff serves every pair; each pair is then
    // measured again against its own cutoff
    let rcut = params.cutoffs.values().chain(params.modifier_cutoffs.values())
        .copied().fold(0.0, f64::max);
    let nl = if rcut > 0.0 {
        Some(NeighborList::build(frame.positions(), &cell, rcut)?)
    } else {
        None
    };
    let near = |i: usize| -> Vec<usize> {
        nl.as_ref().map_or_else(Vec::new, |nl| nl.neighbors_of(i).into_iter().map(|(j, _, _)| j).collect())
    };
    let d2 = |from: usize, to: usize| {
        cell.minimum_image(frame.positions()[to] - frame.positions()[from]).norm_squared()
    };

    // Atom indices by element
    let elem_map = build_elem_map(frame);

    // 1. Each ligand atom's former neighbours (ligand_idx → Vec<(former_elem, former_idx)>)
    let nf_map = build_nf_map(frame, params, &elem_map, &near, &d2);

    // 2. Assemble: Other by default, then overwritten by ligand / former / modifier (the old order)
    let mut types: Vec<AtomType> = frame.atoms()
        .map(|a| AtomType::Other { elem: a.element.to_string() })
        .collect();

    let mut ligand_formers: HashMap<usize, Vec<usize>> = HashMap::with_capacity(nf_map.len());
    for (&idx, nf) in &nf_map {
        let mut partners: Vec<String> = nf.iter().map(|(e, _)| e.clone()).collect();
        partners.sort_unstable();
        types[idx] = AtomType::Ligand { elem: frame.element(idx).to_string(), partners };

        let mut idxs: Vec<usize> = nf.iter().map(|(_, i)| *i).collect();
        idxs.sort_unstable();
        ligand_formers.insert(idx, idxs);
    }
    for (idx, t) in classify_formers(params, &elem_map, &nf_map) {
        types[idx] = t;
    }
    for (idx, t) in classify_modifiers_inner(frame, params, &elem_map, &near, &d2) {
        types[idx] = t;
    }
    Ok(FrameTypes { types, ligand_formers })
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

fn build_elem_map(frame: &Frame) -> HashMap<String, Vec<usize>> {
    let mut map: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, atom) in frame.atoms().enumerate() {
        map.entry(atom.element.to_string()).or_default().push(idx);
    }
    map
}

/// `ligand_idx` → `Vec<(former_elem, former_idx)>` — includes ligand atoms with 0 NF neighbors.
///
/// `near(i)` lists the atoms near `i` in ascending index order, so each former's
/// atoms are visited in the same order as a scan over all of them.
fn build_nf_map(
    frame: &Frame,
    params: &TypeParams,
    elem_map: &HashMap<String, Vec<usize>>,
    near: &dyn Fn(usize) -> Vec<usize>,
    d2: &dyn Fn(usize, usize) -> f64,
) -> HashMap<usize, Vec<(String, usize)>> {
    let mut nf_map: HashMap<usize, Vec<(String, usize)>> = HashMap::new();

    // Every ligand atom gets an entry, bonded or not
    for ligand in params.ligands() {
        if let Some(idxs) = elem_map.get(&ligand) {
            for &idx in idxs { nf_map.entry(idx).or_default(); }
        }
    }

    for ligand in params.ligands() {
        let Some(ligand_idxs) = elem_map.get(&ligand) else { continue };
        let former_pairs = params.formers_for_ligand(&ligand);

        for &la_idx in ligand_idxs {
            let around = near(la_idx);
            for (former, cutoff) in &former_pairs {
                let c2 = cutoff * cutoff;
                for &fa_idx in &around {
                    if frame.element(fa_idx) != *former { continue; }
                    if d2(la_idx, fa_idx) < c2 {
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
///
/// The former–ligand bonds are read off `nf_map`, which holds every one of them:
/// it is built over the same cutoff table, and the minimum-image distance is
/// exactly symmetric, so measuring each pair again from the former's side (as this
/// used to) finds the same bonds.
fn classify_formers(
    params: &TypeParams,
    elem_map: &HashMap<String, Vec<usize>>,
    nf_map: &HashMap<usize, Vec<(String, usize)>>,
) -> Vec<(usize, AtomType)> {
    // former atom → the ligand atoms bonded to it
    let mut bonded: HashMap<usize, Vec<usize>> = HashMap::new();
    for (&la_idx, nf) in nf_map {
        for &(_, fa_idx) in nf {
            bonded.entry(fa_idx).or_default().push(la_idx);
        }
    }

    let mut result: Vec<(usize, AtomType)> = Vec::new();
    for former_elem in params.formers() {
        let Some(former_idxs) = elem_map.get(&former_elem) else { continue };

        for &fa_idx in former_idxs {
            let mut n_bo = 0u32;
            let mut cn = 0u32;
            let mut bridges_to: BTreeMap<String, u32> = BTreeMap::new();

            for &la_idx in bonded.get(&fa_idx).map(|v| v.as_slice()).unwrap_or(&[]) {
                cn += 1;

                // Bridging: the ligand has >= 2 former neighbours
                let nf = &nf_map[&la_idx];
                if nf.len() < 2 { continue; }
                n_bo += 1;

                // The partner breakdown counts *connections*, not bridging ligands: every
                // other former on the ligand is one connection. A tricluster P-O(-Al)(-Al)
                // therefore gives m_Al += 2 — the literature's Q^n_m / P^n_mAl,xB counts
                // exactly the "connections with phosphate", and one tricluster oxygen does
                // link this P to two Al. Counting bridging ligands instead would leave Σm
                // short, and it could no longer be added to n
                for (partner, i) in nf {
                    if *i == fa_idx { continue; }
                    *bridges_to.entry(partner.clone()).or_insert(0) += 1;
                }
            }
            // n counts only *same-element* connections (P-O-P); hetero connections stay in
            // bridges_to as m. That is the literature's Q^n_m: the total number of bridges
            // is n plus every m, not n alone. The convention is fixed at classification
            // time rather than looked up at render time: `--qn` can override it, so a type
            // must carry the convention it was produced under

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
    params: &TypeParams,
    elem_map: &HashMap<String, Vec<usize>>,
    near: &dyn Fn(usize) -> Vec<usize>,
    d2: &dyn Fn(usize, usize) -> f64,
) -> Vec<(usize, AtomType)> {
    let mut result: Vec<(usize, AtomType)> = Vec::new();

    for mod_elem in params.modifiers() {
        let Some(mod_idxs) = elem_map.get(&mod_elem) else { continue };

        for &ma_idx in mod_idxs {
            let around = near(ma_idx);
            let mut cn = 0u32;

            for ((m, ligand), &cutoff) in &params.modifier_cutoffs {
                if m != &mod_elem { continue; }
                let c2 = cutoff * cutoff;
                for &la_idx in &around {
                    if frame.element(la_idx) != *ligand { continue; }
                    if d2(ma_idx, la_idx) < c2 { cn += 1; }
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
    use ferro_core::Cell;

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
        // The old scheme called an over-coordinated oxygen and an over-coordinated
        // modifier both "X" and merged them into one row; O_t and Zn are now distinct
        let modif = AtomType::Modifier { elem: "Zn".into(), cn: 4 };
        assert_ne!(lig(3).label(), modif.label());
    }

    // ─── Qn 与配位数的分叉 ────────────────────────────────────────────────────

    /// A frame in a 20 Å cubic box; positions in Å.
    fn frame_of(atoms: &[(&str, f64, f64, f64)]) -> Frame {
        use ferro_core::Atom;
        use nalgebra::Vector3;
        let cell = Cell::from_lengths_angles(20.0, 20.0, 20.0, 90.0, 90.0, 90.0).unwrap();
        let mut f = Frame::new();
        f.set_atoms(atoms.iter()
            .map(|(e, x, y, z)| Atom::new(*e, Vector3::new(*x, *y, *z))));
        f.cell = Some(cell);
        f.pbc = [true; 3];
        f
    }

    fn params_po_alo() -> TypeParams {
        let mut c = CutoffTable::new();
        c.insert(("P".into(), "O".into()), 2.0);
        c.insert(("Al".into(), "O".into()), 2.0);
        TypeParams::new(c, CutoffTable::new())
    }

    /// 审查 D-M17：奇异胞（c=0 的 slab）曾在最小镜像处 panic，公开函数应报错
    #[test]
    fn test_singular_cell_is_an_error_not_a_panic() {
        let mut frame = frame_of(&[("P", 0.0, 0.0, 0.0), ("O", 1.5, 0.0, 0.0)]);
        frame.cell = Some(Cell::from_matrix(nalgebra::Matrix3::from_diagonal(&nalgebra::Vector3::new(10.0, 10.0, 0.0))));
        let err = classify_frame(&frame, &params_po_alo()).unwrap_err().to_string();
        assert!(err.contains("singular"), "实际 {err}");
    }

    #[test]
    fn test_non_qn_former_labels_by_coordination_not_bridging() {
        // 盲区回归:参考轨迹里的 Al 恰好不带非桥氧,于是 bridging == cn,
        // 两种口径的分叉从未被触发过。这里造一个 Al 带非桥氧的构型把它钉死。
        //
        //   O0 —— 同时挨着 Al 和 P  → 桥氧,计入 Al 的 bridging 与 cn
        //   O1 —— 只挨着 Al        → 非桥氧,只计入 Al 的 cn
        let frame = frame_of(&[
            ("Al", 0.0, 0.0, 0.0),
            ("O",  1.5, 0.0, 0.0),   // Al-O-P 桥
            ("P",  3.0, 0.0, 0.0),
            ("O",  0.0, 1.5, 0.0),   // 只连 Al 的非桥氧
        ]);
        let types = classify_frame(&frame, &params_po_alo()).unwrap();

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

    /// The label, the sort rank and `site_digit` read one number, Qn former or not
    #[test]
    fn test_label_and_rank_show_the_site_digit() {
        let former = |qn, cn| AtomType::Former {
            elem: "Al".into(), qn, n_bo: 1, cn, bridges_to: BTreeMap::new(),
        };
        for t in [former(Some(3), 5), former(None, 5)] {
            let d = t.site_digit().unwrap();
            assert_eq!(t.label(), format!("Al_{d}"));
            assert_eq!(t.class_rank() as u32, d);
        }
        assert_eq!(former(Some(3), 5).site_digit(), Some(3));
        assert_eq!(former(None, 5).site_digit(), Some(5));
    }

    #[test]
    fn test_qn_override_switches_which_digit_a_label_shows() {
        // `--qn` 改的是约定,而约定必须随分类走 —— 同一构型换个约定,
        // Al 的标签要从配位数变成桥接数
        let frame = frame_of(&[
            ("Al", 0.0, 0.0, 0.0),
            ("O",  1.5, 0.0, 0.0),
            ("P",  3.0, 0.0, 0.0),
            ("O",  0.0, 1.5, 0.0),
        ]);
        let params = params_po_alo().with_qn_elements(["Al".to_string()]);
        let types = classify_frame(&frame, &params).unwrap();

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
