//! Ligand atom classification: FO / NBO(X) / BO(X-Y) / OBO(X-Y-Z).

use ferro_core::{Cell, Frame};
use std::collections::HashMap;
use super::NetworkParams;

/// 构建配体原子的 NF 邻居表。
///
/// 返回：`ligand_atom_index` → `Vec<(former_elem, former_atom_index)>`
/// 仅包含截断半径内的网络形成体（NF）邻居。
pub fn build_ligand_nf_map(
    frame: &Frame,
    cell: &Cell,
    params: &NetworkParams,
) -> HashMap<usize, Vec<(String, usize)>> {
    // 按元素建立原子索引列表
    let mut elem_atoms: HashMap<&str, Vec<usize>> = HashMap::new();
    for (idx, atom) in frame.atoms.iter().enumerate() {
        elem_atoms.entry(atom.element.as_str()).or_default().push(idx);
    }

    let mut nf_map: HashMap<usize, Vec<(String, usize)>> = HashMap::new();

    for ligand in params.ligands() {
        let ligand_indices = match elem_atoms.get(ligand.as_str()) {
            Some(v) => v,
            None => continue,
        };

        // 该配体涉及的所有 (former, cutoff)
        let former_cutoffs = params.formers_for_ligand(&ligand);

        for &la_idx in ligand_indices {
            let la_pos = frame.atoms[la_idx].position;
            let entry = nf_map.entry(la_idx).or_default();

            for (former, cutoff) in &former_cutoffs {
                let cutoff2 = cutoff * cutoff;
                let former_indices = match elem_atoms.get(former.as_str()) {
                    Some(v) => v,
                    None => continue,
                };
                for &fa_idx in former_indices {
                    if fa_idx == la_idx { continue; }
                    let diff = cell.minimum_image(frame.atoms[fa_idx].position - la_pos);
                    if diff.norm_squared() < cutoff2 {
                        entry.push((former.clone(), fa_idx));
                    }
                }
            }
        }
    }
    nf_map
}

/// 根据 NF 邻居表对每个配体原子打上分类标签。
///
/// 分类规则：
/// - 0 个 NF 邻居 → `FO`
/// - 1 个 NF 邻居 → `NBO(X)`（X = former 元素符号）
/// - 2 个 NF 邻居 → `BO(X-Y)`（X, Y 字母序排列）
/// - ≥3 个 NF 邻居 → `OBO(X-Y-Z)`（字母序排列）
///
/// 返回：`ligand_elem → Vec<label>`（长度 = 该配体类型原子数）
pub fn classify_ligands(
    frame: &Frame,
    ligand_nf_map: &HashMap<usize, Vec<(String, usize)>>,
    params: &NetworkParams,
) -> HashMap<String, Vec<String>> {
    let mut result: HashMap<String, Vec<String>> = HashMap::new();

    for ligand in params.ligands() {
        let labels = result.entry(ligand.clone()).or_default();
        for (idx, atom) in frame.atoms.iter().enumerate() {
            if atom.element != ligand { continue; }
            let nf_neighbors = ligand_nf_map.get(&idx).map(|v| v.as_slice()).unwrap_or(&[]);
            labels.push(make_label(nf_neighbors));
        }
    }
    result
}

/// 根据 NF 邻居列表生成分类标签字符串。
pub fn make_label(nf_neighbors: &[(String, usize)]) -> String {
    let mut former_elems: Vec<&str> = nf_neighbors.iter().map(|(e, _)| e.as_str()).collect();
    former_elems.sort_unstable();

    match former_elems.len() {
        0 => "FO".to_string(),
        1 => format!("NBO({})", former_elems[0]),
        2 => format!("BO({}-{})", former_elems[0], former_elems[1]),
        _ => format!("OBO({})", former_elems.join("-")),
    }
}
