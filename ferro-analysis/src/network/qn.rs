//! Qn speciation calculation.
//!
//! Qn label definition:
//! - n = number of ligands bonded to the center (former) atom that are of "bridging" type (BO or OBO)
//! - mX = number of those n bridging ligands that bridge through a heterogeneous former X (X ≠ center)
//!
//! Label format:
//! - `Q0`
//! - `Q3(0Al)` — 3 bridging ligands, 0 bridging through Al (all are P-O-P)
//! - `Q3(2Al)` — 3 bridging ligands, 2 of which are P-O-Al bridges
//! - Multiple heterogeneous formers: `Q4(1Al,2Si)` (alphabetical order)

use ferro_core::Frame;
use std::collections::HashMap;
use super::{NetworkParams, ligand_class::make_label};

/// 为每个 former 原子计算 Qn 标签。
///
/// 返回：`former_elem → Vec<qn_label>`（长度 = 该 former 元素原子数）
pub fn calc_qn(
    frame: &Frame,
    neighbor_map: &HashMap<(String, String), Vec<Vec<usize>>>,
    ligand_nf_map: &HashMap<usize, Vec<(String, usize)>>,
    params: &NetworkParams,
) -> HashMap<String, Vec<String>> {
    let formers = params.formers();
    let mut result: HashMap<String, Vec<String>> = HashMap::new();

    // 按 former 元素建立全局原子索引列表（与 neighbor_map 中相同的顺序）
    let mut former_local_to_global: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, atom) in frame.atoms.iter().enumerate() {
        if formers.contains(&atom.element) {
            former_local_to_global.entry(atom.element.clone()).or_default().push(idx);
        }
    }

    for former in &formers {
        let global_indices = match former_local_to_global.get(former) {
            Some(v) => v,
            None => {
                result.insert(former.clone(), Vec::new());
                continue;
            }
        };

        let labels: Vec<String> = global_indices.iter().enumerate().map(|(local_i, &_global_i)| {
            // 收集该 former 原子在各配体类型下的所有邻居配体 atom indices
            let mut all_ligand_neighbors: Vec<usize> = Vec::new();
            for ligand in params.ligands() {
                let key = (former.clone(), ligand);
                if let Some(pair_neighbors) = neighbor_map.get(&key) {
                    if let Some(atom_neighbors) = pair_neighbors.get(local_i) {
                        all_ligand_neighbors.extend_from_slice(atom_neighbors);
                    }
                }
            }

            // 对每个配体邻居，查询其分类标签
            let mut n_bridging: u32 = 0;
            // 统计通过各异元素 X (X ≠ former) 桥接的数量
            let mut hetero_counts: HashMap<String, u32> = HashMap::new();

            for la_idx in all_ligand_neighbors {
                let nf_neighbors = ligand_nf_map.get(&la_idx).map(|v| v.as_slice()).unwrap_or(&[]);
                let label = make_label(nf_neighbors);

                // BO 或 OBO → 桥接配体
                if label.starts_with("BO") || label.starts_with("OBO") {
                    n_bridging += 1;
                    // 每个桥接配体原子，对每种异元素 former 只计一次
                    // （避免 OBO 氧有多个 Zn 邻居时重复计数）
                    let mut seen: std::collections::HashSet<&str> = std::collections::HashSet::new();
                    for (nf_elem, _) in nf_neighbors {
                        if nf_elem != former && seen.insert(nf_elem.as_str()) {
                            *hetero_counts.entry(nf_elem.clone()).or_insert(0) += 1;
                        }
                    }
                }
            }

            make_qn_label(n_bridging, &hetero_counts)
        }).collect();

        result.insert(former.clone(), labels);
    }
    result
}

/// 生成 Qn 标签字符串。
///
/// - 无异元素桥接：`Q3`
/// - 有异元素：`Q3(2Al)` 或 `Q4(1Al,1Si)` （元素字母序）
fn make_qn_label(n: u32, hetero_counts: &HashMap<String, u32>) -> String {
    if hetero_counts.is_empty() {
        return format!("Q{n}");
    }
    // 按元素字母序排列
    let mut parts: Vec<String> = hetero_counts.iter()
        .map(|(elem, &cnt)| format!("{cnt}{elem}"))
        .collect();
    parts.sort();
    format!("Q{n}({})", parts.join(","))
}

/// Qn 标签排序键：Q0 → 0, Q1 → 1, Q2 → 2, …
pub fn qn_label_order(label: &str) -> u32 {
    label.chars().nth(1).and_then(|c| c.to_digit(10)).unwrap_or(99)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_make_qn_label_no_hetero() {
        let h: HashMap<String, u32> = HashMap::new();
        assert_eq!(make_qn_label(0, &h), "Q0");
        assert_eq!(make_qn_label(3, &h), "Q3");
    }

    #[test]
    fn test_make_qn_label_with_hetero() {
        let mut h = HashMap::new();
        h.insert("Al".to_string(), 2u32);
        assert_eq!(make_qn_label(3, &h), "Q3(2Al)");
    }

    #[test]
    fn test_make_qn_label_two_hetero() {
        let mut h = HashMap::new();
        h.insert("Si".to_string(), 1u32);
        h.insert("Al".to_string(), 1u32);
        let label = make_qn_label(4, &h);
        // 字母序：Al before Si
        assert_eq!(label, "Q4(1Al,1Si)");
    }
}
