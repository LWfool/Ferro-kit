//! 网络连通性原语：形成子–配体邻接图 + 连通分量。
//!
//! 玻璃 / 网络体系中，「形成子」（network former，如 P/Si）通过「桥联配体」
//! （bridging ligand，如连接 ≥2 个形成子的 O）相互连接成团簇。本模块提供
//! 计算该图的共享原语，供 `ferro-structure::find_clusters` 与
//! `ferro-analysis::cube_sdf` 等复用，避免重复实现连通分量逻辑。

use crate::{Cell, Frame};

/// 配体按所连形成子数量分类。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LigandKind {
    /// 不连任何形成子（free / 自由氧）
    Free,
    /// 仅连 1 个形成子（non-bridging / 非桥氧）
    NonBridging,
    /// 连 ≥2 个形成子（bridging / 桥氧）
    Bridging,
}

impl LigandKind {
    /// 按形成子邻居数分类。
    pub fn from_former_count(n: usize) -> Self {
        match n {
            0 => LigandKind::Free,
            1 => LigandKind::NonBridging,
            _ => LigandKind::Bridging,
        }
    }

    /// 输出排序权重：桥氧(0) < 非桥氧(1) < 自由氧(2)。
    pub fn sort_rank(&self) -> u8 {
        match self {
            LigandKind::Bridging => 0,
            LigandKind::NonBridging => 1,
            LigandKind::Free => 2,
        }
    }
}

/// 单 (former, ligand) 对的网络邻接图。所有 `*_local` 索引为该图内的
/// 局部下标；`former_global` / `ligand_global` 将其映射回 `frame.atoms` 全局下标。
#[derive(Debug, Clone)]
pub struct NetworkGraph {
    /// former 局部下标 → 全局原子下标
    pub former_global: Vec<usize>,
    /// ligand 局部下标 → 全局原子下标
    pub ligand_global: Vec<usize>,
    /// former_local → 截断内的 [ligand_local]
    pub f_l_adj: Vec<Vec<usize>>,
    /// ligand_local → 截断内的 [former_local]
    pub l_f_adj: Vec<Vec<usize>>,
    /// 每个 ligand 的类型（Free/NonBridging/Bridging）
    pub ligand_kind: Vec<LigandKind>,
    /// 每个 former 的个人 Qn（所连桥联配体数）
    pub former_qn: Vec<u8>,
    /// 连通分量，每个元素为一组 former_local 下标（共享桥联配体即连通）
    pub components: Vec<Vec<usize>>,
}

/// 构建单 (former, ligand) 对的网络图。
///
/// 两个形成子若共享至少一个桥联配体（连 ≥2 形成子的配体）则属同一连通分量。
/// `cutoff` 为形成子–配体键长截断 \[Å\]。无形成子或无配体时返回空图
/// （`components` 为空）。
pub fn build_network_graph(
    frame: &Frame,
    cell: &Cell,
    former: &str,
    ligand: &str,
    cutoff: f64,
) -> NetworkGraph {
    let mut former_global: Vec<usize> = Vec::new();
    let mut ligand_global: Vec<usize> = Vec::new();
    for (i, atom) in frame.atoms.iter().enumerate() {
        if atom.element == former {
            former_global.push(i);
        } else if atom.element == ligand {
            ligand_global.push(i);
        }
    }

    let nf = former_global.len();
    let nl = ligand_global.len();
    let mut f_l_adj: Vec<Vec<usize>> = vec![Vec::new(); nf];
    let mut l_f_adj: Vec<Vec<usize>> = vec![Vec::new(); nl];

    if nf == 0 || nl == 0 {
        return NetworkGraph {
            former_global,
            ligand_global,
            f_l_adj,
            l_f_adj,
            ligand_kind: vec![LigandKind::Free; nl],
            former_qn: vec![0; nf],
            components: Vec::new(),
        };
    }

    let cut2 = cutoff * cutoff;
    for (fi, &fa) in former_global.iter().enumerate() {
        let fpos = frame.atoms[fa].position;
        for (li, &la) in ligand_global.iter().enumerate() {
            let diff = cell
                .minimum_image(frame.atoms[la].position - fpos)
                .expect("cell must be non-singular");
            if diff.norm_squared() < cut2 {
                f_l_adj[fi].push(li);
                l_f_adj[li].push(fi);
            }
        }
    }

    let ligand_kind: Vec<LigandKind> = l_f_adj
        .iter()
        .map(|fs| LigandKind::from_former_count(fs.len()))
        .collect();

    let former_qn: Vec<u8> = f_l_adj
        .iter()
        .map(|ls| {
            ls.iter()
                .filter(|&&li| ligand_kind[li] == LigandKind::Bridging)
                .count() as u8
        })
        .collect();

    // 边：每个桥联配体把其形成子邻居两两相连（链式即可触发并查集合并）
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for (li, kind) in ligand_kind.iter().enumerate() {
        if *kind == LigandKind::Bridging {
            let fs = &l_f_adj[li];
            for w in fs.windows(2) {
                edges.push((w[0], w[1]));
            }
        }
    }

    let (comp_id, n_comp) = connected_components(nf, &edges);
    let mut components: Vec<Vec<usize>> = vec![Vec::new(); n_comp];
    for (fi, &cid) in comp_id.iter().enumerate() {
        components[cid].push(fi);
    }

    NetworkGraph {
        former_global,
        ligand_global,
        f_l_adj,
        l_f_adj,
        ligand_kind,
        former_qn,
        components,
    }
}

/// 通用并查集连通分量。
///
/// 返回 `(component_id_per_node, n_components)`。分量 ID 按节点 `0..n`
/// 遍历时根节点首次出现的顺序分配，结果确定。
pub fn connected_components(n: usize, edges: &[(usize, usize)]) -> (Vec<usize>, usize) {
    let mut parent: Vec<usize> = (0..n).collect();
    for &(a, b) in edges {
        union(&mut parent, a, b);
    }
    let mut root_to_id: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut next_id = 0usize;
    let comp_id: Vec<usize> = (0..n)
        .map(|x| {
            let r = find(&mut parent, x);
            *root_to_id.entry(r).or_insert_with(|| {
                let i = next_id;
                next_id += 1;
                i
            })
        })
        .collect();
    (comp_id, next_id)
}

fn find(parent: &mut Vec<usize>, x: usize) -> usize {
    if parent[x] != x {
        parent[x] = find(parent, parent[x]);
    }
    parent[x]
}

fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
    let ra = find(parent, a);
    let rb = find(parent, b);
    if ra != rb {
        parent[ra] = rb;
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Atom, Cell, Frame};
    use nalgebra::{Matrix3, Vector3};

    fn atom(el: &str, x: f64, y: f64, z: f64) -> Atom {
        Atom::new(el, Vector3::new(x, y, z))
    }

    fn cubic(frame_atoms: Vec<Atom>, l: f64) -> (Frame, Cell) {
        let cell = Cell::from_matrix(Matrix3::from_diagonal(&Vector3::new(l, l, l)));
        let frame = Frame { atoms: frame_atoms, cell: Some(cell.clone()), ..Frame::default() };
        (frame, cell)
    }

    #[test]
    fn connected_components_basic() {
        // 0-1, 2 孤立, 3-4-5 链
        let (id, n) = connected_components(6, &[(0, 1), (3, 4), (4, 5)]);
        assert_eq!(n, 3);
        assert_eq!(id[0], id[1]);
        assert_ne!(id[0], id[2]);
        assert_eq!(id[3], id[5]);
        assert_ne!(id[2], id[3]);
    }

    #[test]
    fn two_formers_bridged_one_cluster() {
        // P – O – P，共享桥氧 → 1 个团簇
        let (frame, cell) = cubic(vec![
            atom("P", 0.0, 0.0, 0.0),
            atom("O", 1.6, 0.0, 0.0),
            atom("P", 3.2, 0.0, 0.0),
        ], 20.0);
        let g = build_network_graph(&frame, &cell, "P", "O", 2.3);
        assert_eq!(g.components.len(), 1);
        assert_eq!(g.components[0].len(), 2);
        assert_eq!(g.ligand_kind[0], LigandKind::Bridging);
        assert_eq!(g.former_qn, vec![1, 1]);
    }

    #[test]
    fn isolated_units_two_clusters() {
        // 两个孤立 P，各自一个非桥氧 → 2 个团簇
        let (frame, cell) = cubic(vec![
            atom("P", 0.0, 0.0, 0.0),
            atom("O", 1.6, 0.0, 0.0),
            atom("P", 10.0, 0.0, 0.0),
            atom("O", 11.6, 0.0, 0.0),
        ], 30.0);
        let g = build_network_graph(&frame, &cell, "P", "O", 2.3);
        assert_eq!(g.components.len(), 2);
        assert_eq!(g.ligand_kind[0], LigandKind::NonBridging);
        assert_eq!(g.former_qn, vec![0, 0]);
    }

    #[test]
    fn empty_when_no_ligand() {
        let (frame, cell) = cubic(vec![atom("P", 0.0, 0.0, 0.0)], 20.0);
        let g = build_network_graph(&frame, &cell, "P", "O", 2.3);
        assert!(g.components.is_empty());
        assert_eq!(g.former_global.len(), 1);
    }
}
