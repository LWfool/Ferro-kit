//! Network atom type classification: per-frame labeling for glass network analysis.
//!
//! Shared by `ferro-structure` (single-frame typing / cluster finding) and
//! `ferro-analysis` (trajectory statistics) without creating cross-crate dependencies.
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

// ─── Public entry point ───────────────────────────────────────────────────────

/// Classify every atom in one frame.  Returns one label per atom (same order as `frame.atoms`).
///
/// Requires `frame.cell` to be `Some` (PBC).  Returns `None` if cell is missing.
pub fn classify_frame(frame: &Frame, cell: &Cell, params: &TypeParams) -> Vec<String> {
    // 按元素建立索引
    let elem_map = build_elem_map(frame);

    // 1. 配体原子的 NF 邻居表（ligand_idx → Vec<(former_elem, former_idx)>）
    let nf_map = build_nf_map(frame, cell, params, &elem_map);

    // 2. 配体标签
    let oxy_labels: HashMap<usize, String> = nf_map.iter()
        .map(|(&idx, nf)| (idx, oxygen_label(nf)))
        .collect();

    // 3. NBO 集合（恰好 1 个 NF 邻居的配体原子）
    let nbo_set: HashSet<usize> = nf_map.iter()
        .filter(|(_, nf)| nf.len() == 1)
        .map(|(&idx, _)| idx)
        .collect();

    // 4. 形成子标签（Qn = 桥氧数）
    let former_labels = classify_formers(frame, cell, params, &elem_map, &nf_map);

    // 5. 修饰子标签
    let modifier_labels = classify_modifiers_inner(frame, cell, params, &elem_map, &nbo_set);

    // 6. 组装
    let mut labels: Vec<String> = frame.atoms.iter().map(|a| a.element.clone()).collect();
    for (idx, lbl) in oxy_labels      { labels[idx] = lbl; }
    for (idx, lbl) in former_labels   { labels[idx] = lbl; }
    for (idx, lbl) in modifier_labels { labels[idx] = lbl; }
    labels
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

fn classify_formers(
    frame: &Frame,
    cell: &Cell,
    params: &TypeParams,
    elem_map: &HashMap<String, Vec<usize>>,
    nf_map: &HashMap<usize, Vec<(String, usize)>>,
) -> HashMap<usize, String> {
    let mut result: HashMap<usize, String> = HashMap::new();

    for former_elem in params.formers() {
        let Some(former_idxs) = elem_map.get(&former_elem) else { continue };

        for &fa_idx in former_idxs {
            let fa_pos = frame.atoms[fa_idx].position;
            let mut bridging = 0u32;

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
                        // 桥接判断：该配体有 ≥2 个 NF 邻居
                        let nf_cnt = nf_map.get(&la_idx).map(|v| v.len()).unwrap_or(0);
                        if nf_cnt >= 2 { bridging += 1; }
                    }
                }
            }
            result.insert(fa_idx, format!("{former_elem}{bridging}"));
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
) -> HashMap<usize, String> {
    let mut result: HashMap<usize, String> = HashMap::new();

    for mod_elem in params.modifiers() {
        let Some(mod_idxs) = elem_map.get(&mod_elem) else { continue };
        // 取该修饰子的最大截断半径
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
            result.insert(ma_idx, modifier_label(&mod_elem, nbo_count));
        }
    }
    result
}

/// 配体标签：`Of` / `On_X` / `Ob_X_Y` / `X`
fn oxygen_label(nf: &[(String, usize)]) -> String {
    let mut elems: Vec<&str> = nf.iter().map(|(e, _)| e.as_str()).collect();
    elems.sort_unstable();
    match elems.len() {
        0 => "Of".to_string(),
        1 => format!("On_{}", elems[0]),
        2 => format!("Ob_{}_{}", elems[0], elems[1]),
        _ => "X".to_string(),
    }
}

/// 修饰子标签：`{Elem}_f` / `{Elem}_t` / `{Elem}_b` / `X`
fn modifier_label(elem: &str, nbo_count: u32) -> String {
    match nbo_count {
        0 => format!("{elem}_f"),
        1 => format!("{elem}_t"),
        2 => format!("{elem}_b"),
        _ => "X".to_string(),
    }
}

// ─── Label ordering helpers (re-exported for sorting in CLI/analysis) ─────────

/// 配体标签排序键：`Of` < `On_*` < `Ob_*` < `X` < other
pub fn oxygen_label_order(label: &str) -> u8 {
    if label == "Of"           { 0 }
    else if label.starts_with("On_") { 1 }
    else if label.starts_with("Ob_") { 2 }
    else if label == "X"       { 3 }
    else                       { 4 }
}

/// 修饰子标签排序键：`_f` < `_t` < `_b` < `X`
pub fn modifier_label_order(label: &str) -> u8 {
    if label.ends_with("_f")   { 0 }
    else if label.ends_with("_t") { 1 }
    else if label.ends_with("_b") { 2 }
    else                       { 3 } // X
}

/// Qn 形成子标签排序键：提取尾部数字，无数字排最后。
pub fn former_label_order(label: &str) -> u32 {
    // "P0" → 0, "P3" → 3, "Al12" → 12
    label.chars()
        .skip_while(|c| c.is_alphabetic())
        .collect::<String>()
        .parse()
        .unwrap_or(u32::MAX)
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
