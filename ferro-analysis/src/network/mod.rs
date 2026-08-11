//! Glass network trajectory statistics: Qn speciation, oxygen type distribution, CN distribution.
//!
//! Atom-level classification is delegated to `ferro_core::classify_frame`.
//! This module accumulates per-frame labels into time-averaged distributions.
//!
//! # Usage
//! ```ignore
//! use ferro_core::TypeParams;
//! use std::collections::BTreeMap;
//!
//! let mut cutoffs = BTreeMap::new();
//! cutoffs.insert(("P".into(), "O".into()), 2.3);
//! let params = TypeParams { cutoffs, modifier_cutoffs: BTreeMap::new() };
//! let result = calc_network(&traj, &params).unwrap();
//! ```

pub use ferro_core::{AtomType, CutoffTable, TypeParams};

use ferro_core::{classify_frame, Frame, Trajectory};
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};

/// Distribution key: `(sort rank, rendered label)`.
///
/// The rank comes from [`AtomType::class_rank`], so a `BTreeMap` keyed on this is
/// sorted by construction — no post-hoc `sort_by` that re-derives the order from
/// the label text.  Labels that collide (`X` covers both an over-coordinated
/// oxygen and an over-coordinated modifier) share one entry, as they always have.
type LabelKey = (u8, String);

fn label_key(t: &AtomType) -> LabelKey {
    (t.class_rank(), t.label())
}

// ─── 结果结构体 ───────────────────────────────────────────────────────────────

/// 网络分析时间平均结果。
#[derive(Debug, Clone)]
pub struct NetworkResult {
    /// former_elem → Qn 分布：`(qn_value, count, fraction)`，按 qn 升序
    pub qn_dist: HashMap<String, Vec<(u32, usize, f64)>>,
    /// former_elem → 平均 Qn
    pub mean_qn: HashMap<String, f64>,
    /// 元素 → 总 CN 分布（所有配体类型之和）：`(cn_value, count, fraction)`。
    /// **含修饰子** —— 修饰子没有桥氧数，配位数就是描述它的全部。
    pub cn_dist: HashMap<String, Vec<(u32, usize, f64)>>,
    /// 元素 → 平均总 CN（含修饰子）
    pub mean_cn: HashMap<String, f64>,
    /// 配体类型分布：`(type_label, count, fraction)`，按 `AtomType::class_rank` 排序
    pub oxy_dist: Vec<(String, usize, f64)>,
}

// ─── 逐帧中间数据 ─────────────────────────────────────────────────────────────

struct FrameData {
    /// former_elem → Vec<(bridging, cn)>，长度 = 该元素的原子数
    former_stats: HashMap<String, Vec<(u32, u32)>>,
    /// modifier_elem → Vec<cn>，长度 = 该元素的原子数
    modifier_cn: HashMap<String, Vec<u32>>,
    /// 各配体标签的计数
    oxy_counts: HashMap<LabelKey, usize>,
}

// ─── 顶层入口 ─────────────────────────────────────────────────────────────────

/// 对整条轨迹执行网络统计分析。要求每帧有 Cell（PBC）。
///
/// 返回 `None` 如果轨迹为空或所有帧都缺少 Cell。
pub fn calc_network(traj: &Trajectory, params: &TypeParams) -> Option<NetworkResult> {
    if traj.frames.is_empty() || params.cutoffs.is_empty() { return None; }

    let acc = traj.frames
        .par_iter()
        .filter_map(|frame| frame.cell.as_ref().map(|cell| (frame, cell)))
        .filter_map(|(frame, cell)| compute_frame(frame, cell, params))
        .fold(
            || Accumulator::new(params),
            |mut acc, fd| { acc.push(&fd); acc },
        )
        .reduce(
            || Accumulator::new(params),
            |mut a, b| { a.merge(b); a },
        );

    Some(acc.finalize())
}

// ─── 单帧计算 ─────────────────────────────────────────────────────────────────

fn compute_frame(
    frame: &Frame,
    cell: &ferro_core::Cell,
    params: &TypeParams,
) -> Option<FrameData> {
    // 1. 每个原子的结构化类型
    let types = classify_frame(frame, cell, params);

    // 2. 建立元素索引
    let mut elem_atoms: HashMap<&str, Vec<usize>> = HashMap::new();
    for (idx, atom) in frame.atoms.iter().enumerate() {
        elem_atoms.entry(atom.element.as_str()).or_default().push(idx);
    }

    let formers = params.formers();
    let ligands = params.ligands();
    let modifiers = params.modifiers();

    // 3. 形成子统计：桥氧数与总 CN 都由 classify_frame 一并算出，
    //    此处不再重扫一遍邻居
    let mut former_stats: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    for former_elem in &formers {
        let Some(fa_idxs) = elem_atoms.get(former_elem.as_str()) else { continue };
        let stats: Vec<(u32, u32)> = fa_idxs.iter()
            .filter_map(|&fa_idx| match &types[fa_idx] {
                AtomType::Former { bridging, cn, .. } => Some((*bridging, *cn)),
                _ => None,
            })
            .collect();
        former_stats.insert(former_elem.clone(), stats);
    }

    // 4. 修饰子：只有配位数
    let mut modifier_cn: HashMap<String, Vec<u32>> = HashMap::new();
    for mod_elem in &modifiers {
        let Some(ma_idxs) = elem_atoms.get(mod_elem.as_str()) else { continue };
        let cns: Vec<u32> = ma_idxs.iter()
            .filter_map(|&ma_idx| match &types[ma_idx] {
                AtomType::Modifier { cn, .. } => Some(*cn),
                _ => None,
            })
            .collect();
        modifier_cn.insert(mod_elem.clone(), cns);
    }

    // 5. 统计配体标签
    let mut oxy_counts: HashMap<LabelKey, usize> = HashMap::new();
    for ligand_elem in &ligands {
        let Some(la_idxs) = elem_atoms.get(ligand_elem.as_str()) else { continue };
        for &la_idx in la_idxs {
            *oxy_counts.entry(label_key(&types[la_idx])).or_insert(0) += 1;
        }
    }

    Some(FrameData { former_stats, modifier_cn, oxy_counts })
}

// ─── 跨帧累加器 ───────────────────────────────────────────────────────────────

struct Accumulator {
    /// former_elem → { bridging → count }
    qn: HashMap<String, HashMap<u32, usize>>,
    /// 元素 → { cn → count }（形成子 + 修饰子）
    cn: HashMap<String, HashMap<u32, usize>>,
    /// 配体标签 → count（BTreeMap：键含 rank，天然有序）
    oxy: BTreeMap<LabelKey, usize>,
}

impl Accumulator {
    fn new(params: &TypeParams) -> Self {
        let qn = params.formers().into_iter().map(|f| (f, HashMap::new())).collect();
        let cn = params.formers().into_iter().chain(params.modifiers())
            .map(|e| (e, HashMap::new())).collect();
        Accumulator { qn, cn, oxy: BTreeMap::new() }
    }

    fn push(&mut self, fd: &FrameData) {
        for (former, stats) in &fd.former_stats {
            let qm = self.qn.entry(former.clone()).or_default();
            let cm = self.cn.entry(former.clone()).or_default();
            for &(bridging, cn) in stats {
                *qm.entry(bridging).or_insert(0) += 1;
                *cm.entry(cn).or_insert(0) += 1;
            }
        }
        for (mod_elem, cns) in &fd.modifier_cn {
            let cm = self.cn.entry(mod_elem.clone()).or_default();
            for &cn in cns { *cm.entry(cn).or_insert(0) += 1; }
        }
        for (key, &count) in &fd.oxy_counts {
            *self.oxy.entry(key.clone()).or_insert(0) += count;
        }
    }

    fn merge(&mut self, other: Self) {
        for (k, inner) in other.qn {
            let m = self.qn.entry(k).or_default();
            for (q, c) in inner { *m.entry(q).or_insert(0) += c; }
        }
        for (k, inner) in other.cn {
            let m = self.cn.entry(k).or_default();
            for (c, cnt) in inner { *m.entry(c).or_insert(0) += cnt; }
        }
        for (label, c) in other.oxy {
            *self.oxy.entry(label).or_insert(0) += c;
        }
    }

    fn finalize(self) -> NetworkResult {
        let mean_of = |counts: &HashMap<u32, usize>| -> f64 {
            let total: usize = counts.values().sum();
            if total == 0 { return 0.0; }
            counts.iter().map(|(&v, &c)| v as f64 * c as f64).sum::<f64>() / total as f64
        };
        let to_dist = |counts: &HashMap<u32, usize>| -> Vec<(u32, usize, f64)> {
            let total: usize = counts.values().sum();
            let mut rows: Vec<_> = counts.iter()
                .map(|(&v, &c)| (v, c, if total > 0 { c as f64 / total as f64 } else { 0.0 }))
                .collect();
            rows.sort_by_key(|r| r.0);
            rows
        };

        let qn_dist: HashMap<_, _> = self.qn.iter()
            .map(|(f, m)| (f.clone(), to_dist(m))).collect();
        let mean_qn: HashMap<_, _> = self.qn.iter()
            .map(|(f, m)| (f.clone(), mean_of(m))).collect();
        let cn_dist: HashMap<_, _> = self.cn.iter()
            .map(|(f, m)| (f.clone(), to_dist(m))).collect();
        let mean_cn: HashMap<_, _> = self.cn.iter()
            .map(|(f, m)| (f.clone(), mean_of(m))).collect();

        // 配体类型分布（BTreeMap 的键是 (class_rank, label)，已按 _f < _n < _b < _t 排好）
        let oxy_total: usize = self.oxy.values().sum();
        let oxy_dist: Vec<(String, usize, f64)> = self.oxy.iter()
            .map(|((_, lbl), &c)| (
                lbl.clone(), c,
                if oxy_total > 0 { c as f64 / oxy_total as f64 } else { 0.0 },
            ))
            .collect();

        NetworkResult { qn_dist, mean_qn, cn_dist, mean_cn, oxy_dist }
    }
}

// ─── 测试 ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame, Trajectory, TypeParams};
    use nalgebra::{Matrix3, Vector3};
    use std::collections::BTreeMap;

    fn atom(elem: &str, x: f64, y: f64, z: f64) -> Atom {
        Atom { element: elem.to_string(), position: Vector3::new(x, y, z),
               label: None, mass: None, magmom: None, charge: None }
    }

    fn make_params(p_o: f64) -> TypeParams {
        let mut c = BTreeMap::new();
        c.insert(("P".into(), "O".into()), p_o);
        TypeParams { cutoffs: c, modifier_cutoffs: BTreeMap::new() }
    }

    /// P–O–P 体系：一个 Q1 P，一个 Q1 P，桥氧一个
    #[test]
    fn test_q1_system() {
        // P1 – Ob – P2，两个 NBO（各自一个），一个桥氧
        //  positions: P1@0, Ob@1.6, P2@3.2, On1@-1.6, On2@4.8, On3@0(y=1.6), On4@3.2(y=1.6)
        let atoms = vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0), // 桥氧
            atom("P",  3.2, 0.0, 0.0),
            atom("O", -1.6, 0.0, 0.0), // NBO(P1)
            atom("O",  4.8, 0.0, 0.0), // NBO(P2)
        ];
        let cell = Cell::from_matrix(Matrix3::from_diagonal(&Vector3::new(20.0, 20.0, 20.0)));
        let frame = Frame { atoms, cell: Some(cell), ..Frame::default() };
        let traj = Trajectory { frames: vec![frame], metadata: Default::default() };
        let params = make_params(2.3);
        let res = calc_network(&traj, &params).unwrap();

        // 两个 P 都是 Q1
        let qn_p = &res.qn_dist["P"];
        assert_eq!(qn_p.len(), 1);
        assert_eq!(qn_p[0].0, 1); // Qn=1
        assert_eq!(qn_p[0].1, 2); // 2 个原子

        // 氧分布：1 桥氧(O_b) + 2 非桥氧(O_n)
        let oxy: HashMap<&str, usize> = res.oxy_dist.iter()
            .map(|(l, c, _)| (l.as_str(), *c)).collect();
        assert_eq!(oxy["O_b"], 1);
        assert_eq!(oxy["O_n"], 2);
    }

    /// 孤立 PO4（Q0）
    #[test]
    fn test_q0_system() {
        let atoms = vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0),
            atom("O", -1.6, 0.0, 0.0),
            atom("O",  0.0, 1.6, 0.0),
            atom("O",  0.0,-1.6, 0.0),
        ];
        let cell = Cell::from_matrix(Matrix3::from_diagonal(&Vector3::new(20.0, 20.0, 20.0)));
        let frame = Frame { atoms, cell: Some(cell), ..Frame::default() };
        let traj = Trajectory { frames: vec![frame], metadata: Default::default() };
        let params = make_params(2.3);
        let res = calc_network(&traj, &params).unwrap();

        let qn_p = &res.qn_dist["P"];
        assert_eq!(qn_p[0].0, 0); // Qn=0

        // 所有氧都是非桥氧
        assert!(res.oxy_dist.iter().all(|(l, _, _)| l == "O_n"));
    }

    /// 异质桥（P–O–Al）与同质桥（P–O–P）共用标签 `O_b`，合并成一行。
    /// 伙伴元素是统计量，不再编码进标签 —— 这是标签方案里唯一有信息损失的地方，
    /// 由 linkage 表补回。
    #[test]
    fn test_hetero_and_homo_bridges_share_one_label() {
        let atoms = vec![
            atom("P",  0.0, 0.0, 0.0),
            atom("O",  1.6, 0.0, 0.0), // P–O–Al 桥
            atom("Al", 3.2, 0.0, 0.0),
            atom("P",  0.0, 5.0, 0.0),
            atom("O",  1.6, 5.0, 0.0), // P–O–P 桥
            atom("P",  3.2, 5.0, 0.0),
        ];
        let cell = Cell::from_matrix(Matrix3::from_diagonal(&Vector3::new(20.0, 20.0, 20.0)));
        let frame = Frame { atoms, cell: Some(cell), ..Frame::default() };
        let traj = Trajectory { frames: vec![frame], metadata: Default::default() };

        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("P".to_string(), "O".to_string()), 2.3);
        cutoffs.insert(("Al".to_string(), "O".to_string()), 2.3);
        let params = TypeParams { cutoffs, modifier_cutoffs: BTreeMap::new() };
        let res = calc_network(&traj, &params).unwrap();

        assert_eq!(res.oxy_dist.len(), 1, "两种桥氧必须并成一行: {:?}", res.oxy_dist);
        assert_eq!(res.oxy_dist[0].0, "O_b");
        assert_eq!(res.oxy_dist[0].1, 2);
    }

    /// 修饰子只进 CN 表，且**不影响配体分类** —— 挨着修饰子的非桥氧仍是 `O_n`，
    /// 不会因为多了个 Zn 邻居而变成桥氧。这正是 `--modifier` 存在的理由。
    #[test]
    fn test_modifier_counts_cn_only() {
        let atoms = vec![
            atom("P",   0.0, 0.0, 0.0),
            atom("O",   1.6, 0.0, 0.0),
            atom("O",  -1.6, 0.0, 0.0),
            atom("Zn",  1.6, 2.0, 0.0), // 距 O(1.6,0,0) 2.0 Å，距另一个 O 3.8 Å
        ];
        let cell = Cell::from_matrix(Matrix3::from_diagonal(&Vector3::new(20.0, 20.0, 20.0)));
        let frame = Frame { atoms, cell: Some(cell), ..Frame::default() };
        let traj = Trajectory { frames: vec![frame], metadata: Default::default() };

        let mut cutoffs = BTreeMap::new();
        cutoffs.insert(("P".to_string(), "O".to_string()), 2.3);
        let mut modifier_cutoffs = BTreeMap::new();
        modifier_cutoffs.insert(("Zn".to_string(), "O".to_string()), 2.6);
        let params = TypeParams { cutoffs, modifier_cutoffs };
        let res = calc_network(&traj, &params).unwrap();

        assert_eq!(res.cn_dist["Zn"], vec![(1, 1, 1.0)], "Zn 进 CN 表");
        assert!(!res.qn_dist.contains_key("Zn"), "Zn 不进 Qn 表 —— 修饰子没有桥氧数");
        assert!(res.oxy_dist.iter().all(|(l, _, _)| l == "O_n"),
                "Zn 不能把非桥氧变成桥氧: {:?}", res.oxy_dist);
        assert_eq!(res.qn_dist["P"][0].0, 0, "P 仍是 Q0");
    }
}
