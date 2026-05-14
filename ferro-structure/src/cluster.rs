//! Network cluster identification.
//!
//! A "cluster" is a connected component of network-former atoms linked through
//! bridging oxygen atoms (oxygen with ≥2 NF neighbors).  Used by `cube_sdf`
//! and downstream analysis that needs per-cluster statistics.

use ferro_core::{classify_frame, Frame, TypeParams};
use std::collections::HashMap;

/// Result of cluster identification for one frame.
#[derive(Debug, Clone)]
pub struct ClusterResult {
    /// Per-atom cluster ID.  `None` for atoms that are not network formers.
    pub cluster_id: Vec<Option<usize>>,
    /// Total number of clusters found.
    pub n_clusters: usize,
}

impl ClusterResult {
    /// Atom indices belonging to cluster `id`.
    pub fn members(&self, id: usize) -> Vec<usize> {
        self.cluster_id.iter().enumerate()
            .filter_map(|(i, c)| if *c == Some(id) { Some(i) } else { None })
            .collect()
    }

    /// All cluster IDs present.
    pub fn ids(&self) -> Vec<usize> {
        (0..self.n_clusters).collect()
    }
}

/// Identify connected clusters of network-former atoms in one frame.
///
/// Two former atoms are in the same cluster when they share at least one
/// bridging oxygen (a ligand atom bonded to ≥2 formers).
///
/// Returns `None` if `frame.cell` is missing.
pub fn find_clusters(frame: &Frame, params: &TypeParams) -> Option<ClusterResult> {
    let cell = frame.cell.as_ref()?;
    let labels = classify_frame(frame, cell, params);

    // 1. 找到所有形成子原子的索引
    let formers = params.formers();
    let former_set: std::collections::HashSet<&str> =
        formers.iter().map(|s| s.as_str()).collect();
    let former_indices: Vec<usize> = frame.atoms.iter().enumerate()
        .filter(|(_, a)| former_set.contains(a.element.as_str()))
        .map(|(i, _)| i)
        .collect();

    // 2. Union-Find
    let n = former_indices.len();
    let mut parent: Vec<usize> = (0..n).collect();

    // 建立 atom_idx → local_idx 映射
    let local: HashMap<usize, usize> = former_indices.iter().enumerate()
        .map(|(li, &ai)| (ai, li))
        .collect();

    // 3. 找桥氧：标签以 "Ob_" 开头
    //    对每个桥氧，找其 NF 邻居，合并这些 NF 所属的连通分量
    let mut elem_map: HashMap<&str, Vec<usize>> = HashMap::new();
    for (idx, atom) in frame.atoms.iter().enumerate() {
        elem_map.entry(atom.element.as_str()).or_default().push(idx);
    }

    for (la_idx, label) in labels.iter().enumerate() {
        if !label.starts_with("Ob_") { continue; }
        // 收集该桥氧在截断内的所有形成子邻居
        let la_pos = frame.atoms[la_idx].position;
        let mut nf_locals: Vec<usize> = Vec::new();

        for former_elem in &formers {
            let Some(&cutoff) = params.cutoffs.keys()
                .find(|(f, _)| f == former_elem)
                .and_then(|k| params.cutoffs.get(k))
            else { continue };
            let c2 = cutoff * cutoff;
            let Some(fa_idxs) = elem_map.get(former_elem.as_str()) else { continue };
            for &fa_idx in fa_idxs {
                if fa_idx == la_idx { continue; }
                let diff = cell.minimum_image(frame.atoms[fa_idx].position - la_pos)
                    .expect("cell must be non-singular");
                if diff.norm_squared() < c2 {
                    if let Some(&li) = local.get(&fa_idx) {
                        nf_locals.push(li);
                    }
                }
            }
        }
        // 合并所有相邻形成子
        for i in 1..nf_locals.len() {
            union(&mut parent, nf_locals[0], nf_locals[i]);
        }
    }

    // 4. 压缩根节点，生成连续 cluster ID
    let mut root_to_id: HashMap<usize, usize> = HashMap::new();
    let mut next_id = 0usize;
    let local_cluster: Vec<usize> = (0..n).map(|li| {
        let r = find(&mut parent, li);
        *root_to_id.entry(r).or_insert_with(|| { let i = next_id; next_id += 1; i })
    }).collect();

    // 5. 组装结果（per-atom）
    let mut cluster_id: Vec<Option<usize>> = vec![None; frame.atoms.len()];
    for (li, &ai) in former_indices.iter().enumerate() {
        cluster_id[ai] = Some(local_cluster[li]);
    }

    Some(ClusterResult { cluster_id, n_clusters: next_id })
}

// ─── Union-Find ───────────────────────────────────────────────────────────────

fn find(parent: &mut Vec<usize>, x: usize) -> usize {
    if parent[x] != x { parent[x] = find(parent, parent[x]); }
    parent[x]
}

fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
    let ra = find(parent, a);
    let rb = find(parent, b);
    if ra != rb { parent[ra] = rb; }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame, TypeParams};
    use nalgebra::{Matrix3, Vector3};
    use std::collections::BTreeMap;

    fn make_params(p_o: f64) -> TypeParams {
        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("P".into(), "O".into()), p_o);
        TypeParams { cutoffs, modifier_cutoffs: BTreeMap::new() }
    }

    fn atom(elem: &str, x: f64, y: f64, z: f64) -> Atom {
        Atom { element: elem.to_string(), position: Vector3::new(x, y, z),
               label: None, mass: None, magmom: None, charge: None }
    }

    /// 两个 P 原子共享一个桥氧 → 同一团簇
    #[test]
    fn test_two_p_bridged() {
        // P1 – O – P2，盒子 20 Å
        let atoms = vec![
            atom("P", 0.0, 0.0, 0.0),
            atom("O", 1.6, 0.0, 0.0),
            atom("P", 3.2, 0.0, 0.0),
        ];
        let cell = Cell::from_matrix(Matrix3::from_diagonal(&Vector3::new(20.0, 20.0, 20.0)));
        let frame = Frame { atoms, cell: Some(cell), ..Frame::default() };
        let params = make_params(2.3);
        let res = find_clusters(&frame, &params).unwrap();
        // P1 和 P2 在同一团簇
        assert_eq!(res.cluster_id[0], res.cluster_id[2]);
        assert_eq!(res.n_clusters, 1);
    }

    /// 两个孤立的 PO4（无共享桥氧）→ 两个团簇
    #[test]
    fn test_two_isolated_po4() {
        // P1 完全被 NBO 包围，P2 也是，不共享桥氧
        let atoms = vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0), // NBO(P1)
            atom("P", 10.0, 0.0, 0.0),
            atom("O", 11.6, 0.0, 0.0), // NBO(P2)
        ];
        let cell = Cell::from_matrix(Matrix3::from_diagonal(&Vector3::new(30.0, 30.0, 30.0)));
        let frame = Frame { atoms, cell: Some(cell), ..Frame::default() };
        let params = make_params(2.3);
        let res = find_clusters(&frame, &params).unwrap();
        // P1 和 P2 在不同团簇
        assert_ne!(res.cluster_id[0], res.cluster_id[2]);
        assert_eq!(res.n_clusters, 2);
    }
}
