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
//! | Atom role | Label |
//! |-----------|-------|
//! | Network former (Qn) | `P0`, `P1`, `Al2`, … |
//! | Free oxygen (0 NF) | `Of` |
//! | Non-bridging oxygen (1 NF) | `On_P`, `On_Al`, … |
//! | Bridging oxygen (2 NF) | `Ob_Al_P`, `Ob_P_P`, … (alphabetical) |
//! | Over-bridging oxygen (≥3 NF) | `X` |
//! | Modifier – free | `Zn_f`, `Na_f`, … |
//! | Modifier – terminal (1 NBO) | `Zn_t`, … |
//! | Modifier – bridging (2 NBO) | `Zn_b`, … |
//! | Modifier – over (≥3 NBO) | `X` |
//! | Other atoms | element symbol unchanged |

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
    /// its length is the classification (0 free, 1 non-bridging, 2 bridging, ≥3 over).
    Ligand { elem: String, partners: Vec<String> },
    /// Network modifier.  `nbo` = non-bridging ligands within reach; `cn` = total
    /// ligands within cutoff.
    Modifier { elem: String, nbo: u32, cn: u32 },
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
            AtomType::Former { elem, bridging, .. } => format!("{elem}{bridging}"),
            AtomType::Ligand { partners, .. } => match partners.len() {
                0 => "Of".to_string(),
                1 => format!("On_{}", partners[0]),
                2 => format!("Ob_{}_{}", partners[0], partners[1]),
                _ => "X".to_string(),
            },
            AtomType::Modifier { elem, nbo, .. } => match nbo {
                0 => format!("{elem}_f"),
                1 => format!("{elem}_t"),
                2 => format!("{elem}_b"),
                _ => "X".to_string(),
            },
            AtomType::Other { elem } => elem.clone(),
        }
    }

    /// Sort rank *within one distribution table* (all rows share a role).
    ///
    /// Oxygen: `Of` < `On_*` < `Ob_*` < `X`.  Modifier: `_f` < `_t` < `_b` < `X`.
    pub fn class_rank(&self) -> u8 {
        match self {
            AtomType::Former { bridging, .. } => (*bridging).min(u8::MAX as u32) as u8,
            AtomType::Ligand { partners, .. } => partners.len().min(3) as u8,
            AtomType::Modifier { nbo, .. } => (*nbo).min(3) as u8,
            AtomType::Other { .. } => 0,
        }
    }

    /// Sort rank when every role is printed in *one* table:
    /// formers, then oxygen (free → non-bridging → bridging), then modifiers,
    /// then the over-coordinated `X` bucket, then everything else.
    pub fn display_rank(&self) -> u8 {
        match self {
            AtomType::Former { .. } => 0,
            AtomType::Ligand { partners, .. } => match partners.len() {
                0 => 1,
                1 => 2,
                2 => 3,
                _ => 5,
            },
            AtomType::Modifier { nbo, .. } => if *nbo <= 2 { 4 } else { 5 },
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

    // 2. NBO 集合（恰好 1 个 NF 邻居的配体原子）
    let nbo_set: HashSet<usize> = nf_map.iter()
        .filter(|(_, nf)| nf.len() == 1)
        .map(|(&idx, _)| idx)
        .collect();

    // 3. 组装：默认 Other，再依次被配体 / 形成子 / 修饰子覆盖（顺序同旧实现）
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
    for (idx, t) in classify_modifiers_inner(frame, cell, params, &elem_map, &nbo_set) {
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

fn classify_modifiers_inner(
    frame: &Frame,
    cell: &Cell,
    params: &TypeParams,
    elem_map: &HashMap<String, Vec<usize>>,
    nbo_set: &HashSet<usize>,
) -> Vec<(usize, AtomType)> {
    let mut result: Vec<(usize, AtomType)> = Vec::new();

    for mod_elem in params.modifiers() {
        let Some(mod_idxs) = elem_map.get(&mod_elem) else { continue };
        // NBO 计数用该修饰子的最大截断半径（与逐配体截断的 cn 口径不同,
        // 但 nbo 只服务于 `Zn_f`/`Zn_t`/`Zn_b` 角色标签,该标签即将退役）
        let max_cut = params.modifier_cutoffs.iter()
            .filter(|((m, _), _)| *m == mod_elem)
            .map(|(_, &c)| c)
            .fold(0.0_f64, f64::max);
        let max_cut2 = max_cut * max_cut;

        for &ma_idx in mod_idxs {
            let ma_pos = frame.atoms[ma_idx].position;
            let nbo_count = nbo_set.iter()
                .filter(|&&nbo_idx| {
                    if nbo_idx == ma_idx { return false; }
                    let diff = cell.minimum_image(frame.atoms[nbo_idx].position - ma_pos)
                        .expect("cell must be non-singular");
                    diff.norm_squared() < max_cut2
                })
                .count() as u32;

            // 总配位数：逐配体截断
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

            result.push((ma_idx, AtomType::Modifier {
                elem: mod_elem.clone(), nbo: nbo_count, cn,
            }));
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
    fn test_labels_match_legacy_scheme() {
        let f = |b| AtomType::Former { elem: "P".into(), bridging: b, cn: 4 };
        assert_eq!(f(0).label(), "P0");
        assert_eq!(f(3).label(), "P3");

        let lig = |p: &[&str]| AtomType::Ligand {
            elem: "O".into(),
            partners: p.iter().map(|s| s.to_string()).collect(),
        };
        assert_eq!(lig(&[]).label(), "Of");
        assert_eq!(lig(&["P"]).label(), "On_P");
        assert_eq!(lig(&["Al", "P"]).label(), "Ob_Al_P");
        assert_eq!(lig(&["Al", "P", "P"]).label(), "X");

        let m = |n| AtomType::Modifier { elem: "Zn".into(), nbo: n, cn: 5 };
        assert_eq!(m(0).label(), "Zn_f");
        assert_eq!(m(1).label(), "Zn_t");
        assert_eq!(m(2).label(), "Zn_b");
        assert_eq!(m(3).label(), "X");

        assert_eq!(AtomType::Other { elem: "Ar".into() }.label(), "Ar");
    }

    #[test]
    fn test_class_rank_orders_oxygen_and_modifier() {
        let lig = |n: usize| AtomType::Ligand {
            elem: "O".into(),
            partners: vec!["P".to_string(); n],
        };
        // Of < On_* < Ob_* < X，且 ≥3 全部落在同一档
        assert!(lig(0).class_rank() < lig(1).class_rank());
        assert!(lig(1).class_rank() < lig(2).class_rank());
        assert!(lig(2).class_rank() < lig(3).class_rank());
        assert_eq!(lig(3).class_rank(), lig(5).class_rank());

        let m = |n| AtomType::Modifier { elem: "Zn".into(), nbo: n, cn: 4 };
        assert!(m(0).class_rank() < m(1).class_rank());
        assert!(m(2).class_rank() < m(3).class_rank());
        assert_eq!(m(3).class_rank(), m(9).class_rank());
    }

    #[test]
    fn test_display_rank_puts_both_x_buckets_together() {
        // 过配位氧与过配位修饰子都渲染成 "X"，混排时必须落在同一档，
        // 否则同名两行会在类型统计表里分开出现
        let over_o = AtomType::Ligand {
            elem: "O".into(),
            partners: vec!["P".into(), "P".into(), "P".into()],
        };
        let over_m = AtomType::Modifier { elem: "Zn".into(), nbo: 4, cn: 6 };
        assert_eq!(over_o.label(), over_m.label());
        assert_eq!(over_o.display_rank(), over_m.display_rank());

        let former = AtomType::Former { elem: "P".into(), bridging: 2, cn: 4 };
        let free_o = AtomType::Ligand { elem: "O".into(), partners: vec![] };
        let modif = AtomType::Modifier { elem: "Zn".into(), nbo: 1, cn: 4 };
        let other = AtomType::Other { elem: "Ar".into() };
        assert!(former.display_rank() < free_o.display_rank());
        assert!(free_o.display_rank() < modif.display_rank());
        assert!(modif.display_rank() < over_o.display_rank());
        assert!(over_o.display_rank() < other.display_rank());
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
