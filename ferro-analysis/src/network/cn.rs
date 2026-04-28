//! Per-atom coordination number (CN) calculation.

use ferro_core::{Cell, Frame};
use std::collections::HashMap;
use super::NetworkParams;

/// 构建邻居索引表。
///
/// 返回：`(former_elem, ligand_elem)` → 每个 former 原子的邻居 ligand 原子索引列表。
/// 外层 Vec 长度 = 该 former 元素的原子数；内层 Vec 包含邻居的全局原子索引。
pub fn build_neighbor_map(
    frame: &Frame,
    cell: &Cell,
    params: &NetworkParams,
) -> HashMap<(String, String), Vec<Vec<usize>>> {
    // 按元素建立原子索引列表
    let mut elem_atoms: HashMap<&str, Vec<usize>> = HashMap::new();
    for (idx, atom) in frame.atoms.iter().enumerate() {
        elem_atoms.entry(atom.element.as_str()).or_default().push(idx);
    }

    let mut result: HashMap<(String, String), Vec<Vec<usize>>> = HashMap::new();

    for ((former, ligand), &cutoff) in &params.cutoffs {
        let former_indices = match elem_atoms.get(former.as_str()) {
            Some(v) => v,
            None => {
                // 帧中没有该元素：插入空列表
                result.insert((former.clone(), ligand.clone()), Vec::new());
                continue;
            }
        };
        let ligand_indices = match elem_atoms.get(ligand.as_str()) {
            Some(v) => v,
            None => {
                result.insert(
                    (former.clone(), ligand.clone()),
                    vec![Vec::new(); former_indices.len()],
                );
                continue;
            }
        };

        let cutoff2 = cutoff * cutoff;
        let mut pair_neighbors: Vec<Vec<usize>> = vec![Vec::new(); former_indices.len()];

        for (fi, &fa_idx) in former_indices.iter().enumerate() {
            let fa_pos = frame.atoms[fa_idx].position;
            for &la_idx in ligand_indices {
                // 排除 former == ligand 时的自身
                if la_idx == fa_idx { continue; }
                let diff = cell.minimum_image(frame.atoms[la_idx].position - fa_pos);
                if diff.norm_squared() < cutoff2 {
                    pair_neighbors[fi].push(la_idx);
                }
            }
        }
        result.insert((former.clone(), ligand.clone()), pair_neighbors);
    }
    result
}

/// 从邻居表提取每个 former 原子对各配体的 CN 值。
pub fn collect_cn(
    neighbor_map: &HashMap<(String, String), Vec<Vec<usize>>>,
) -> HashMap<(String, String), Vec<u32>> {
    neighbor_map.iter()
        .map(|(pair, neighbors)| {
            let cns = neighbors.iter().map(|v| v.len() as u32).collect();
            (pair.clone(), cns)
        })
        .collect()
}

/// 计算每个 former 原子的总 CN（跨所有配体类型求和）。
///
/// 对同一 former 的多个配体 pair 的 CN 按原子位置相加。
pub fn collect_cn_total(
    neighbor_map: &HashMap<(String, String), Vec<Vec<usize>>>,
    params: &NetworkParams,
) -> HashMap<String, Vec<u32>> {
    let mut total: HashMap<String, Vec<u32>> = HashMap::new();

    for former in params.formers() {
        // 收集该 former 的所有配体 pair 的 CN 列表
        let mut base: Option<Vec<u32>> = None;
        for ligand in params.ligands() {
            let key = (former.clone(), ligand);
            if let Some(neighbors) = neighbor_map.get(&key) {
                let cns: Vec<u32> = neighbors.iter().map(|v| v.len() as u32).collect();
                match base {
                    None => base = Some(cns),
                    Some(ref mut b) => {
                        for (bi, ci) in b.iter_mut().zip(cns) { *bi += ci; }
                    }
                }
            }
        }
        if let Some(v) = base {
            total.insert(former, v);
        }
    }
    total
}
