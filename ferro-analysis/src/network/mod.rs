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

pub use ferro_core::{
    TypeParams, CutoffTable,
    oxygen_label_order, modifier_label_order, former_label_order,
};

use ferro_core::{classify_frame, Frame, Trajectory};
use rayon::prelude::*;
use std::collections::HashMap;

// ─── 结果结构体 ───────────────────────────────────────────────────────────────

/// 网络分析时间平均结果。
#[derive(Debug, Clone)]
pub struct NetworkResult {
    /// former_elem → Qn 分布：`(qn_value, count, fraction)`，按 qn 升序
    pub qn_dist: HashMap<String, Vec<(u32, usize, f64)>>,
    /// former_elem → 平均 Qn
    pub mean_qn: HashMap<String, f64>,
    /// former_elem → 总 CN 分布（所有配体类型之和）：`(cn_value, count, fraction)`
    pub cn_dist: HashMap<String, Vec<(u32, usize, f64)>>,
    /// former_elem → 平均总 CN
    pub mean_cn: HashMap<String, f64>,
    /// 氧类型分布：`(type_label, count, fraction)`，按 oxygen_label_order 排序
    pub oxy_dist: Vec<(String, usize, f64)>,
    /// modifier_elem → 角色分布：`(role_label, count, fraction)`，按 modifier_label_order 排序
    pub modifier_dist: HashMap<String, Vec<(String, usize, f64)>>,
}

// ─── 逐帧中间数据 ─────────────────────────────────────────────────────────────

struct FrameData {
    /// former_elem → Vec<(qn_value, cn_value)>，长度 = 该元素的原子数
    former_stats: HashMap<String, Vec<(u32, u32)>>,
    /// 各氧标签的计数：label → count
    oxy_counts: HashMap<String, usize>,
    /// modifier_elem → Vec<role_label>
    modifier_labels: HashMap<String, Vec<String>>,
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
    // 1. 获取每个原子的类型标签
    let labels = classify_frame(frame, cell, params);

    // 2. 建立元素索引（用于 CN 计算）
    let mut elem_atoms: HashMap<&str, Vec<usize>> = HashMap::new();
    for (idx, atom) in frame.atoms.iter().enumerate() {
        elem_atoms.entry(atom.element.as_str()).or_default().push(idx);
    }

    let formers = params.formers();
    let ligands = params.ligands();
    let modifiers = params.modifiers();

    // 3. 提取形成子统计（Qn + 总 CN）
    let mut former_stats: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    for former_elem in &formers {
        let Some(fa_idxs) = elem_atoms.get(former_elem.as_str()) else { continue };
        let stats: Vec<(u32, u32)> = fa_idxs.iter().map(|&fa_idx| {
            // 从标签提取 Qn（标签格式 "P0"、"Al3" 等）
            let qn = extract_qn(&labels[fa_idx]);

            // 计算总 CN（所有配体类型中截断内的邻居数之和）
            let total_cn: u32 = ligands.iter().map(|ligand_elem| {
                let Some(&cutoff) = params.cutoffs.get(
                    &(former_elem.clone(), ligand_elem.clone())
                ) else { return 0u32 };
                let c2 = cutoff * cutoff;
                let fa_pos = frame.atoms[fa_idx].position;
                let Some(la_idxs) = elem_atoms.get(ligand_elem.as_str()) else { return 0u32 };
                la_idxs.iter().filter(|&&la_idx| {
                    if la_idx == fa_idx { return false; }
                    let diff = cell.minimum_image(frame.atoms[la_idx].position - fa_pos)
                        .expect("cell is non-singular");
                    diff.norm_squared() < c2
                }).count() as u32
            }).sum();

            (qn, total_cn)
        }).collect();
        former_stats.insert(former_elem.clone(), stats);
    }

    // 4. 统计氧标签
    let mut oxy_counts: HashMap<String, usize> = HashMap::new();
    for ligand_elem in &ligands {
        let Some(la_idxs) = elem_atoms.get(ligand_elem.as_str()) else { continue };
        for &la_idx in la_idxs {
            *oxy_counts.entry(labels[la_idx].clone()).or_insert(0) += 1;
        }
    }

    // 5. 提取修饰子标签
    let mut modifier_labels: HashMap<String, Vec<String>> = HashMap::new();
    for mod_elem in &modifiers {
        let Some(ma_idxs) = elem_atoms.get(mod_elem.as_str()) else { continue };
        let role_labels: Vec<String> = ma_idxs.iter()
            .map(|&ma_idx| labels[ma_idx].clone())
            .collect();
        modifier_labels.insert(mod_elem.clone(), role_labels);
    }

    Some(FrameData { former_stats, oxy_counts, modifier_labels })
}

/// 从形成子标签提取 Qn 数字（"P3" → 3，"Al12" → 12，解析失败 → 0）
fn extract_qn(label: &str) -> u32 {
    label.chars()
        .skip_while(|c| c.is_alphabetic())
        .collect::<String>()
        .parse()
        .unwrap_or(0)
}

// ─── 跨帧累加器 ───────────────────────────────────────────────────────────────

struct Accumulator {
    /// former_elem → { qn → count }
    qn: HashMap<String, HashMap<u32, usize>>,
    /// former_elem → { cn → count }
    cn: HashMap<String, HashMap<u32, usize>>,
    /// oxygen_label → count
    oxy: HashMap<String, usize>,
    /// modifier_elem → { role_label → count }
    modifier: HashMap<String, HashMap<String, usize>>,
}

impl Accumulator {
    fn new(params: &TypeParams) -> Self {
        let qn = params.formers().into_iter().map(|f| (f, HashMap::new())).collect();
        let cn = params.formers().into_iter().map(|f| (f, HashMap::new())).collect();
        let oxy = HashMap::new();
        let modifier = params.modifiers().into_iter().map(|m| (m, HashMap::new())).collect();
        Accumulator { qn, cn, oxy, modifier }
    }

    fn push(&mut self, fd: &FrameData) {
        for (former, stats) in &fd.former_stats {
            let qm = self.qn.entry(former.clone()).or_default();
            let cm = self.cn.entry(former.clone()).or_default();
            for &(qn, cn) in stats {
                *qm.entry(qn).or_insert(0) += 1;
                *cm.entry(cn).or_insert(0) += 1;
            }
        }
        for (label, &count) in &fd.oxy_counts {
            *self.oxy.entry(label.clone()).or_insert(0) += count;
        }
        for (mod_elem, labels) in &fd.modifier_labels {
            let mm = self.modifier.entry(mod_elem.clone()).or_default();
            for lbl in labels { *mm.entry(lbl.clone()).or_insert(0) += 1; }
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
        for (k, inner) in other.modifier {
            let m = self.modifier.entry(k).or_default();
            for (lbl, c) in inner { *m.entry(lbl).or_insert(0) += c; }
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

        // 氧类型分布（排序：Of < On_* < Ob_* < X）
        let oxy_total: usize = self.oxy.values().sum();
        let mut oxy_dist: Vec<(String, usize, f64)> = self.oxy.iter()
            .map(|(lbl, &c)| (
                lbl.clone(), c,
                if oxy_total > 0 { c as f64 / oxy_total as f64 } else { 0.0 },
            ))
            .collect();
        oxy_dist.sort_by(|a, b| {
            oxygen_label_order(&a.0).cmp(&oxygen_label_order(&b.0)).then(a.0.cmp(&b.0))
        });

        // 修饰子分布
        let modifier_dist: HashMap<_, _> = self.modifier.iter()
            .map(|(mod_elem, counts)| {
                let total: usize = counts.values().sum();
                let mut rows: Vec<(String, usize, f64)> = counts.iter()
                    .map(|(lbl, &c)| (
                        lbl.clone(), c,
                        if total > 0 { c as f64 / total as f64 } else { 0.0 },
                    ))
                    .collect();
                rows.sort_by_key(|r| modifier_label_order(&r.0));
                (mod_elem.clone(), rows)
            })
            .collect();

        NetworkResult { qn_dist, mean_qn, cn_dist, mean_cn, oxy_dist, modifier_dist }
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

        // 氧分布：1 桥氧(Ob_P_P) + 2 NBO(On_P)
        let oxy: HashMap<&str, usize> = res.oxy_dist.iter()
            .map(|(l, c, _)| (l.as_str(), *c)).collect();
        assert_eq!(oxy["Ob_P_P"], 1);
        assert_eq!(oxy["On_P"], 2);
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

        // 所有氧都是 NBO
        assert!(res.oxy_dist.iter().all(|(l, _, _)| l == "On_P"));
    }
}
