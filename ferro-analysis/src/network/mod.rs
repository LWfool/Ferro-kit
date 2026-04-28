//! Glass network analysis: per-atom CN, ligand classification (FO/NBO/BO/OBO), Qn speciation.
//!
//! # Usage
//! ```ignore
//! let mut cutoffs = HashMap::new();
//! cutoffs.insert(("P".into(), "O".into()), 2.3);
//! cutoffs.insert(("P".into(), "F".into()), 2.1);
//! let params = NetworkParams { cutoffs };
//! let result = calc_network(&traj, &params).unwrap();
//! ```

mod cn;
mod ligand_class;
mod qn;

pub use qn::qn_label_order;

use ferro_core::{Cell, Frame, Trajectory};
use std::collections::HashMap;

/// (network_former_element, ligand_element) → cutoff radius [Å]
pub type CutoffTable = HashMap<(String, String), f64>;

/// Input parameters for network analysis.
#[derive(Debug, Clone)]
pub struct NetworkParams {
    pub cutoffs: CutoffTable,
}

impl NetworkParams {
    pub fn formers(&self) -> Vec<String> {
        let mut v: Vec<String> = self.cutoffs.keys()
            .map(|(f, _)| f.clone())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        v.sort();
        v
    }

    pub fn ligands(&self) -> Vec<String> {
        let mut v: Vec<String> = self.cutoffs.keys()
            .map(|(_, l)| l.clone())
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();
        v.sort();
        v
    }

    pub fn cutoff(&self, former: &str, ligand: &str) -> Option<f64> {
        self.cutoffs.get(&(former.to_string(), ligand.to_string())).copied()
    }

    /// 查询指定配体 l 对应的所有 (former, cutoff) 对
    pub fn formers_for_ligand(&self, ligand: &str) -> Vec<(String, f64)> {
        self.cutoffs.iter()
            .filter(|((_, l), _)| l == ligand)
            .map(|((f, _), &c)| (f.clone(), c))
            .collect()
    }
}

// ─── 结果结构体 ───────────────────────────────────────────────────────────────

/// 最终输出结果（所有帧累加后归一化）
#[derive(Debug, Clone)]
pub struct NetworkResult {
    /// (former, ligand) → CN 分布行：(cn_value, count, fraction)
    pub cn_dist: HashMap<(String, String), Vec<(u32, usize, f64)>>,
    /// former → total CN 分布（所有配体之和）
    pub cn_total: HashMap<String, Vec<(u32, usize, f64)>>,
    /// (former, ligand) → 平均 CN
    pub mean_cn: HashMap<(String, String), f64>,
    /// former → 平均总 CN
    pub mean_cn_total: HashMap<String, f64>,
    /// ligand → 分类分布：(class_label, count, fraction)，按标签字母序
    pub ligand_classes: HashMap<String, Vec<(String, usize, f64)>>,
    /// former → Qn 物种分布：(species_label, count, fraction)，按 Qn 序
    pub qn_species: HashMap<String, Vec<(String, usize, f64)>>,
}

// ─── 逐帧中间数据 ─────────────────────────────────────────────────────────────

/// 单帧分析的原始原子级数据（无需对外暴露）
struct FrameData {
    /// (former, ligand) → 每个 former 原子的 CN 值（长度 = 该 former 类型原子数）
    cn_by_pair: HashMap<(String, String), Vec<u32>>,
    /// former → 每个 former 原子的总 CN（所有配体合计）
    cn_total_by_former: HashMap<String, Vec<u32>>,
    /// ligand → 每个配体原子的分类标签（长度 = 该配体类型原子数）
    ligand_labels: HashMap<String, Vec<String>>,
    /// former → 每个 former 原子的 Qn 标签
    qn_labels: HashMap<String, Vec<String>>,
}

// ─── 顶层入口 ─────────────────────────────────────────────────────────────────

/// 对整条轨迹执行网络分析，要求每帧有 Cell（PBC）。
pub fn calc_network(traj: &Trajectory, params: &NetworkParams) -> Option<NetworkResult> {
    if traj.frames.is_empty() { return None; }
    if traj.frames.iter().any(|f| f.cell.is_none()) { return None; }
    if params.cutoffs.is_empty() { return None; }

    let mut acc = Accumulator::new(params);
    for frame in &traj.frames {
        let cell = frame.cell.as_ref().unwrap();
        if let Some(fd) = compute_frame(frame, cell, params) {
            acc.push(&fd);
        }
    }
    Some(acc.finalize())
}

// ─── 单帧计算 ─────────────────────────────────────────────────────────────────

fn compute_frame(frame: &Frame, cell: &Cell, params: &NetworkParams) -> Option<FrameData> {
    // 1. 对每个 former 原子，按配体类型建立邻居索引表
    //    neighbor_map: (former, ligand) → Vec<Vec<usize>>
    //    外层 Vec 长度 = 该 former 元素的原子数，内层 = 该原子的配体邻居 atom index
    let neighbor_map = cn::build_neighbor_map(frame, cell, params);

    // 2. 对每个配体原子，建立其 NF 邻居列表（用于配体分类和 Qn）
    //    ligand_nf_map: atom_index → Vec<(former_elem, former_atom_index)>
    let ligand_nf_map = ligand_class::build_ligand_nf_map(frame, cell, params);

    // 3. CN 统计
    let cn_by_pair = cn::collect_cn(&neighbor_map);

    // 4. 总 CN（每个 former 原子跨所有配体类型的 CN 之和）
    let cn_total_by_former = cn::collect_cn_total(&neighbor_map, params);

    // 5. 配体分类
    let ligand_labels = ligand_class::classify_ligands(frame, &ligand_nf_map, params);

    // 6. Qn 计算
    let qn_labels = qn::calc_qn(frame, &neighbor_map, &ligand_nf_map, params);

    Some(FrameData { cn_by_pair, cn_total_by_former, ligand_labels, qn_labels })
}

// ─── 跨帧累加器 ───────────────────────────────────────────────────────────────

struct Accumulator {
    /// (former, ligand) → { cn_value → count }
    cn_pair: HashMap<(String, String), HashMap<u32, usize>>,
    /// former → { cn_value → count }
    cn_total: HashMap<String, HashMap<u32, usize>>,
    /// ligand → { class_label → count }
    lig_class: HashMap<String, HashMap<String, usize>>,
    /// former → { qn_label → count }
    qn: HashMap<String, HashMap<String, usize>>,
}

impl Accumulator {
    fn new(params: &NetworkParams) -> Self {
        let cn_pair = params.cutoffs.keys().map(|k| (k.clone(), HashMap::new())).collect();
        let cn_total = params.formers().into_iter().map(|f| (f, HashMap::new())).collect();
        let lig_class = params.ligands().into_iter().map(|l| (l, HashMap::new())).collect();
        let qn = params.formers().into_iter().map(|f| (f, HashMap::new())).collect();
        Accumulator { cn_pair, cn_total, lig_class, qn }
    }

    fn push(&mut self, fd: &FrameData) {
        for (pair, cns) in &fd.cn_by_pair {
            let m = self.cn_pair.entry(pair.clone()).or_default();
            for &v in cns { *m.entry(v).or_insert(0) += 1; }
        }
        for (former, cns) in &fd.cn_total_by_former {
            let m = self.cn_total.entry(former.clone()).or_default();
            for &v in cns { *m.entry(v).or_insert(0) += 1; }
        }
        for (lig, labels) in &fd.ligand_labels {
            let m = self.lig_class.entry(lig.clone()).or_default();
            for l in labels { *m.entry(l.clone()).or_insert(0) += 1; }
        }
        for (former, labels) in &fd.qn_labels {
            let m = self.qn.entry(former.clone()).or_default();
            for l in labels { *m.entry(l.clone()).or_insert(0) += 1; }
        }
    }

    fn finalize(self) -> NetworkResult {
        let to_dist = |counts: &HashMap<u32, usize>| -> Vec<(u32, usize, f64)> {
            let total: usize = counts.values().sum();
            let mut rows: Vec<_> = counts.iter()
                .map(|(&cn, &c)| (cn, c, if total > 0 { c as f64 / total as f64 } else { 0.0 }))
                .collect();
            rows.sort_by_key(|r| r.0);
            rows
        };
        let mean_of = |counts: &HashMap<u32, usize>| -> f64 {
            let total: usize = counts.values().sum();
            if total == 0 { return 0.0; }
            counts.iter().map(|(&cn, &c)| cn as f64 * c as f64).sum::<f64>() / total as f64
        };

        let cn_dist: HashMap<_, _> = self.cn_pair.iter()
            .map(|(k, v)| (k.clone(), to_dist(v)))
            .collect();
        let mean_cn: HashMap<_, _> = self.cn_pair.iter()
            .map(|(k, v)| (k.clone(), mean_of(v)))
            .collect();

        let cn_total: HashMap<_, _> = self.cn_total.iter()
            .map(|(k, v)| (k.clone(), to_dist(v)))
            .collect();
        let mean_cn_total: HashMap<_, _> = self.cn_total.iter()
            .map(|(k, v)| (k.clone(), mean_of(v)))
            .collect();

        let ligand_classes: HashMap<_, _> = self.lig_class.iter()
            .map(|(lig, counts)| {
                let total: usize = counts.values().sum();
                let mut rows: Vec<(String, usize, f64)> = counts.iter()
                    .map(|(l, &c)| (l.clone(), c, if total > 0 { c as f64 / total as f64 } else { 0.0 }))
                    .collect();
                rows.sort_by(|a, b| ligand_class_order(&a.0).cmp(&ligand_class_order(&b.0)).then(a.0.cmp(&b.0)));
                (lig.clone(), rows)
            })
            .collect();

        let qn_species: HashMap<_, _> = self.qn.iter()
            .map(|(former, counts)| {
                let total: usize = counts.values().sum();
                let mut rows: Vec<(String, usize, f64)> = counts.iter()
                    .map(|(l, &c)| (l.clone(), c, if total > 0 { c as f64 / total as f64 } else { 0.0 }))
                    .collect();
                rows.sort_by(|a, b| {
                    qn_label_order(&a.0).cmp(&qn_label_order(&b.0)).then(a.0.cmp(&b.0))
                });
                (former.clone(), rows)
            })
            .collect();

        NetworkResult { cn_dist, cn_total, mean_cn, mean_cn_total, ligand_classes, qn_species }
    }
}

/// FO < NBO < BO < OBO 排序键
fn ligand_class_order(label: &str) -> u8 {
    if label.starts_with("FO")  { 0 }
    else if label.starts_with("NBO") { 1 }
    else if label.starts_with("BO")  { 2 }
    else                             { 3 } // OBO
}
