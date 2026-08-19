//! Network cluster identification.
//!
//! A "cluster" is a connected component of network-former atoms linked through
//! bridging oxygen atoms (oxygen with ≥2 NF neighbors).  Used by `cube_sdf`
//! and downstream analysis that needs per-cluster statistics.

use ferro_core::{classify_frame, connected_components, AtomType, Frame, TypeParams};
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
    let types = classify_frame(frame, cell, params);

    // 1. 找到所有形成子原子的索引
    let formers = params.formers();
    let former_set: std::collections::HashSet<&str> =
        formers.iter().map(|s| s.as_str()).collect();
    let former_indices: Vec<usize> = frame.atoms.iter().enumerate()
        .filter(|(_, a)| former_set.contains(a.element.as_str()))
        .map(|(i, _)| i)
        .collect();

    // 2. atom_idx → local_idx 映射
    let n = former_indices.len();
    let local: HashMap<usize, usize> = former_indices.iter().enumerate()
        .map(|(li, &ai)| (ai, li))
        .collect();

    // 3. 找桥氧：恰好连接两个形成子的配体
    //    对每个桥氧，找其 NF 邻居，合并这些 NF 所属的连通分量
    let mut elem_map: HashMap<&str, Vec<usize>> = HashMap::new();
    for (idx, atom) in frame.atoms.iter().enumerate() {
        elem_map.entry(atom.element.as_str()).or_default().push(idx);
    }

    let mut edges: Vec<(usize, usize)> = Vec::new();
    for (la_idx, atom_type) in types.iter().enumerate() {
        // 只有配体角色参与连边。桥联判据是「≥2 个形成子」，在下面按实际邻居数
        // 判定 —— 不能用 `AtomType::is_bridging()`，它是「**恰好**两个」，会把
        // 三簇配体（连 3 个形成子，网络里连得最紧的节点）当成完全不连通
        let AtomType::Ligand { elem: lig_elem, .. } = atom_type else { continue };
        // 收集该配体在截断内的所有形成子邻居
        let la_pos = frame.atoms[la_idx].position;
        let mut nf_locals: Vec<usize> = Vec::new();

        for former_elem in &formers {
            // 截断按 (形成子, **该配体元素**) 取。双配体体系（`--Al-O` + `--Al-F`）
            // 下 Al-O 与 Al-F 是两个不同的键长，取错会把一种键判成全断或全连
            let Some(cutoff) = params.cutoff(former_elem, lig_elem) else { continue };
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
        if nf_locals.len() < 2 { continue; }
        // 该桥联配体的相邻形成子两两连通（链式即可；三簇配体给 2 条边连通 3 个）
        for i in 1..nf_locals.len() {
            edges.push((nf_locals[0], nf_locals[i]));
        }
    }

    // 4. 连通分量（复用 ferro-core 并查集，分量 ID 按首见根确定）
    let (local_cluster, next_id) = connected_components(n, &edges);

    // 5. 组装结果（per-atom）
    let mut cluster_id: Vec<Option<usize>> = vec![None; frame.atoms.len()];
    for (li, &ai) in former_indices.iter().enumerate() {
        cluster_id[ai] = Some(local_cluster[li]);
    }

    Some(ClusterResult { cluster_id, n_clusters: next_id })
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
        TypeParams::new(cutoffs, BTreeMap::new())
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

    /// 三簇氧（连 3 个形成子）必须连通它的三个形成子。
    /// 模块 doc 与 `ferro_core::LigandKind` 都定义桥联为「≥2 个形成子」。
    #[test]
    fn test_tricluster_connects_three_formers() {
        // 一个 O 被 3 个 P 包围（正三角形中心）
        let r = 1.6;
        let atoms = vec![
            atom("O", 0.0, 0.0, 0.0),
            atom("P", r, 0.0, 0.0),
            atom("P", -r / 2.0,  r * 0.8660254, 0.0),
            atom("P", -r / 2.0, -r * 0.8660254, 0.0),
        ];
        let cell = Cell::from_matrix(Matrix3::from_diagonal(&Vector3::new(30.0, 30.0, 30.0)));
        let frame = Frame { atoms, cell: Some(cell), ..Frame::default() };
        let params = make_params(2.3);
        let res = find_clusters(&frame, &params).unwrap();
        assert_eq!(res.n_clusters, 1, "三簇氧把三个 P 连成一个团簇");
        assert_eq!(res.cluster_id[1], res.cluster_id[2]);
        assert_eq!(res.cluster_id[1], res.cluster_id[3]);
    }

    /// 双配体体系：截断必须按 (形成子, **该配体元素**) 取。
    /// `BTreeMap` 里 ("P","F") 排在 ("P","O") 之前，按形成子元素取「第一个」
    /// 截断会拿 P-F 的 1.4 Å 去判 P-O 的 1.6 Å 键，桥氧整体判不出来。
    #[test]
    fn test_cutoff_is_per_ligand_element() {
        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("P".to_string(), "O".to_string()), 2.3);
        cutoffs.insert(("P".to_string(), "F".to_string()), 1.4);
        let params = TypeParams::new(cutoffs, BTreeMap::new());

        let atoms = vec![
            atom("P", 0.0, 0.0, 0.0),
            atom("O", 1.6, 0.0, 0.0), // 桥氧：P-O = 1.6 Å，在 2.3 内、在 1.4 外
            atom("P", 3.2, 0.0, 0.0),
            atom("F", 0.0, 9.0, 0.0), // 远处的 F，只为让 F 进入配体元素集
        ];
        let cell = Cell::from_matrix(Matrix3::from_diagonal(&Vector3::new(30.0, 30.0, 30.0)));
        let frame = Frame { atoms, cell: Some(cell), ..Frame::default() };
        let res = find_clusters(&frame, &params).unwrap();
        assert_eq!(res.n_clusters, 1, "P-O 桥必须按 P-O 的截断判定");
        assert_eq!(res.cluster_id[0], res.cluster_id[2]);
    }
}
