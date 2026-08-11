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
//! | Network former | `P_0`, `P_3`, `Al_2`, … (digit = bridging ligands) |
//! | Free ligand (0 NF) | `O_f` |
//! | Non-bridging ligand (1 NF) | `O_n` |
//! | Bridging ligand (2 NF) | `O_b` |
//! | Tricluster ligand (≥3 NF) | `O_t` |
//! | Modifier | element symbol unchanged (`Zn`) |
//! | Other atoms | element symbol unchanged |
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

use crate::{Cell, Frame};
use std::collections::{BTreeMap, HashMap, HashSet};

/// Cutoff table: `(element_A, element_B)` → distance [Å]
pub type CutoffTable = BTreeMap<(String, String), f64>;

/// Parameters for network atom type classification.
#[derive(Debug, Clone)]
pub struct TypeParams {
    /// `(former_elem, ligand_elem)` → cutoff [Å]
    pub cutoffs: CutoffTable,
    /// `(modifier_elem, ligand_elem)` → cutoff [Å].  Empty → no modifier classification.
    pub modifier_cutoffs: CutoffTable,
}

impl TypeParams {
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
    /// Network former.  `bridging` = number of bridging ligands (the Qn value for
    /// elements where Qn applies); `cn` = total ligands within cutoff.
    Former { elem: String, bridging: u32, cn: u32 },
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

    /// Render as a site label.  **The only place a type becomes text.**
    pub fn label(&self) -> String {
        match self {
            AtomType::Former { elem, bridging, .. } => format!("{elem}_{bridging}"),
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
    /// Ligand: `_f` < `_n` < `_b` < `_t`.
    pub fn class_rank(&self) -> u8 {
        match self {
            AtomType::Former { bridging, .. } => (*bridging).min(u8::MAX as u32) as u8,
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

    /// True for a ligand bridging exactly two formers.
    pub fn is_bridging(&self) -> bool {
        matches!(self, AtomType::Ligand { partners, .. } if partners.len() == 2)
    }
}

// ─── Public entry point ───────────────────────────────────────────────────────

/// Classify every atom in one frame.  Returns one type per atom (same order as `frame.atoms`).
pub fn classify_frame(frame: &Frame, cell: &Cell, params: &TypeParams) -> Vec<AtomType> {
    // 按元素建立索引
    let elem_map = build_elem_map(frame);

    // 1. 配体原子的 NF 邻居表（ligand_idx → Vec<(former_elem, former_idx)>）
    let nf_map = build_nf_map(frame, cell, params, &elem_map);

    // 2. 组装：默认 Other，再依次被配体 / 形成子 / 修饰子覆盖（顺序同旧实现）
    let mut types: Vec<AtomType> = frame.atoms.iter()
        .map(|a| AtomType::Other { elem: a.element.clone() })
        .collect();

    for (&idx, nf) in &nf_map {
        let mut partners: Vec<String> = nf.iter().map(|(e, _)| e.clone()).collect();
        partners.sort_unstable();
        types[idx] = AtomType::Ligand { elem: frame.atoms[idx].element.clone(), partners };
    }
    for (idx, t) in classify_formers(frame, cell, params, &elem_map, &nf_map) {
        types[idx] = t;
    }
    for (idx, t) in classify_modifiers_inner(frame, cell, params, &elem_map) {
        types[idx] = t;
    }
    types
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
            let mut bridging = 0u32;
            let mut cn = 0u32;

            for ligand_elem in params.ligands() {
                let Some(&cutoff) = params.cutoffs.get(&(former_elem.clone(), ligand_elem.clone()))
                else { continue };
                let c2 = cutoff * cutoff;
                let Some(ligand_idxs) = elem_map.get(&ligand_elem) else { continue };

                for &la_idx in ligand_idxs {
                    if la_idx == fa_idx { continue; }
                    let diff = cell.minimum_image(frame.atoms[la_idx].position - fa_pos)
                        .expect("cell must be non-singular");
                    if diff.norm_squared() < c2 {
                        cn += 1;
                        // 桥接判断：该配体有 ≥2 个 NF 邻居
                        let nf_cnt = nf_map.get(&la_idx).map(|v| v.len()).unwrap_or(0);
                        if nf_cnt >= 2 { bridging += 1; }
                    }
                }
            }
            result.push((fa_idx, AtomType::Former {
                elem: former_elem.clone(), bridging, cn,
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

    #[test]
    fn test_every_label_splits_at_first_underscore() {
        // 全部标签必须形如 <元素>_<后缀>（或裸元素），否则 LAMMPS dump 读回来
        // 拆不出正确的 element，`-a P` 就选不中导出的结构
        let cases: Vec<(AtomType, &str, &str)> = vec![
            (AtomType::Former { elem: "P".into(), bridging: 0, cn: 4 }, "P_0", "P"),
            (AtomType::Former { elem: "P".into(), bridging: 3, cn: 4 }, "P_3", "P"),
            (AtomType::Former { elem: "Al".into(), bridging: 2, cn: 5 }, "Al_2", "Al"),
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
        let former = AtomType::Former { elem: "P".into(), bridging: 2, cn: 4 };
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
}
