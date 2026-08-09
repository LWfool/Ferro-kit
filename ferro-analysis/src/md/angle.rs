//! Bond angle distribution calculation and output.
//!
//! Three-atom angle A-B-C with B as the center atom.
//! Automatically enumerates all chemical element triplets.
//!
//! Canonical key `"ElemA-ElemCenter-ElemC"`:
//!   - B is always the center; endpoints sorted by atomic number, Z(A) ≤ Z(C).
//!   - Symmetric pairs (e.g. Si-O-Si) are not double-counted: enumeration requires A_idx < C_idx.
//!
//! `rcut_ab`: cutoff from end atom A to center; `rcut_bc`: from end atom C to center.
//! Which physical end is "A" comes from `AngleParams::ends` — the order the caller wrote
//! `-a`/`-c` — falling back to the canonical (Z, symbol) order when the triplet was not
//! named. When both end atoms are the same type, `min(rcut_ab, rcut_bc)` is used for both.
//!
//! Counting: each geometric angle once (a PO₄ tetrahedron gives 6 O-P-O angles, not 12).
//! `code1/dump2analysis` enumerates ordered end pairs instead, so its histogram is exactly
//! twice this one whenever both ends are the same type; mean/std/peak positions are equal.
//!
//! Parallelism: per-frame `par_iter().fold().reduce()`, same pattern as gr.rs.
//! Algorithm reference: code1/angle.c (`EstimateAngle`).

use nalgebra::{Matrix3, Vector3};
use rayon::prelude::*;
use std::collections::BTreeMap;
use ferro_core::{Table, Trajectory};
use super::gr::{GroupBy, elem_z, group_key, sorted_types};

// ─── 内部辅助 ─────────────────────────────────────────────────────────────────

/// Construct a canonical triplet key `"ElemA-ElemCenter-ElemC"` with Z(A) ≤ Z(C).
/// When Z values are equal, lexicographic order is used to guarantee deterministic ordering for pseudo-element labels.
pub(crate) fn canonical_triplet(a: &str, b: &str, c: &str) -> String {
    // 先比较 Z，Z 相同时按字符串字典序决定顺序
    let swap = match elem_z(a).cmp(&elem_z(c)) {
        std::cmp::Ordering::Greater => true,
        std::cmp::Ordering::Equal   => a > c,
        std::cmp::Ordering::Less    => false,
    };
    if swap { format!("{}-{}-{}", c, b, a) } else { format!("{}-{}-{}", a, b, c) }
}

/// Sort triplet keys by (Z_center, center_label, Z_left, left_label, Z_right, right_label).
/// Secondary string ordering ensures deterministic ordering for pseudo-element labels sharing the same Z.
fn sort_triplet_keys(map: &BTreeMap<String, Vec<u64>>) -> Vec<String> {
    let mut keys: Vec<String> = map.keys().cloned().collect();
    keys.sort_by(|ka, kb| {
        // 格式 "A-B-C"：split_once('-') 得 (A, "B-C")，再 split_once('-') 得 (B, C)
        // 使用 String 避免生命周期问题
        let parse = |k: &str| -> (u8, String, u8, String, u8, String) {
            if let Some((a, rest)) = k.split_once('-') {
                if let Some((b, c)) = rest.split_once('-') {
                    return (elem_z(b), b.to_string(), elem_z(a), a.to_string(),
                            elem_z(c), c.to_string());
                }
            }
            (255, k.to_string(), 255, String::new(), 255, String::new())
        };
        parse(ka).cmp(&parse(kb))
    });
    keys
}

// ─── linked-cell list ────────────────────────────────────────────────────────

pub(super) struct CellList {
    nx: usize,
    ny: usize,
    nz: usize,
    /// cells[ix + iy*nx + iz*nx*ny] = atom indices in that cell
    cells: Vec<Vec<usize>>,
    /// 各原子分数坐标（各分量 ∈ [0, 1)）
    frac: Vec<[f64; 3]>,
    /// cell.matrix.transpose()，用于分数差 → 笛卡尔向量
    mat_t: Matrix3<f64>,
}

impl CellList {
    pub(super) fn build(frame: &ferro_core::Frame, cell: &ferro_core::Cell, max_rcut: f64) -> Self {
        let mat_t = cell.matrix.transpose();
        // (mat_t)^{-1} 的第 i 行是倒格矢 i，其模 = 1/d_i（面间距）
        // 每轴 cell 数 = floor(d_i / max_rcut).max(1)，对三斜晶胞完全正确
        let mat_t_inv = mat_t.try_inverse().unwrap_or(mat_t);
        let nx = ((1.0 / (max_rcut * mat_t_inv.row(0).norm())).floor() as usize).max(1);
        let ny = ((1.0 / (max_rcut * mat_t_inv.row(1).norm())).floor() as usize).max(1);
        let nz = ((1.0 / (max_rcut * mat_t_inv.row(2).norm())).floor() as usize).max(1);

        let mut cells = vec![Vec::new(); nx * ny * nz];

        let frac: Vec<[f64; 3]> = frame.atoms.iter().enumerate().map(|(i, a)| {
            let f = mat_t_inv * a.position;
            let fx = f.x.rem_euclid(1.0);
            let fy = f.y.rem_euclid(1.0);
            let fz = f.z.rem_euclid(1.0);
            let ix = ((fx * nx as f64) as usize).min(nx - 1);
            let iy = ((fy * ny as f64) as usize).min(ny - 1);
            let iz = ((fz * nz as f64) as usize).min(nz - 1);
            cells[ix + iy * nx + iz * nx * ny].push(i);
            [fx, fy, fz]
        }).collect();

        CellList { nx, ny, nz, cells, frac, mat_t }
    }

    /// 返回 b_idx 在 max_rcut 内的所有邻居：(atom_idx, dist, [dx, dy, dz] Å)
    pub(super) fn neighbors_of(&self, b_idx: usize, max_rcut: f64) -> Vec<(usize, f64, [f64; 3])> {
        let fb = self.frac[b_idx];
        let radius2 = max_rcut * max_rcut;

        let cx = (fb[0] * self.nx as f64) as i64;
        let cy = (fb[1] * self.ny as f64) as i64;
        let cz = (fb[2] * self.nz as f64) as i64;

        let mut result = Vec::new();
        // 最多 27 个相邻 cell；盒子过小时同一 cell 会被多次映射，用线性扫描去重
        let mut seen: Vec<usize> = Vec::with_capacity(27);

        for dz in -1i64..=1 {
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let ix = ((cx + dx).rem_euclid(self.nx as i64)) as usize;
                    let iy = ((cy + dy).rem_euclid(self.ny as i64)) as usize;
                    let iz = ((cz + dz).rem_euclid(self.nz as i64)) as usize;
                    let flat = ix + iy * self.nx + iz * self.nx * self.ny;
                    if seen.contains(&flat) { continue; }
                    seen.push(flat);

                    for &other in &self.cells[flat] {
                        if other == b_idx { continue; }
                        let fo = self.frac[other];
                        let mut df = Vector3::new(
                            fo[0] - fb[0],
                            fo[1] - fb[1],
                            fo[2] - fb[2],
                        );
                        // 分数空间最小镜像，等价于 cell.minimum_image 但省去矩阵求逆
                        df.x -= df.x.round();
                        df.y -= df.y.round();
                        df.z -= df.z.round();
                        let cart = self.mat_t * df;
                        let dist2 = cart.norm_squared();
                        if dist2 > 0.0 && dist2 < radius2 {
                            result.push((other, dist2.sqrt(), [cart.x, cart.y, cart.z]));
                        }
                    }
                }
            }
        }
        result
    }
}

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Parameters for bond angle distribution calculation.
#[derive(Debug, Clone)]
pub struct AngleParams {
    /// Distance cutoff from end atom A to the center B \[Å\] (default: 2.3).
    ///
    /// Which physical end counts as "A" is decided by `ends`: when the caller names the
    /// triplet, A is the end they wrote first; otherwise A is the lower-Z end.
    pub r_cut_ab: f64,
    /// Distance cutoff from end atom C to the center B \[Å\] (default: 2.3). See `r_cut_ab`.
    pub r_cut_bc: f64,
    /// Lower edge of the histogram \[degrees\] (default: 0.0)
    pub angle_min: f64,
    /// Upper edge of the histogram \[degrees\] (default: 180.0)
    pub angle_max: f64,
    /// Histogram bin width \[degrees\] (default: 0.1, same as code1 180/1800 split)
    pub d_angle: f64,
    /// Whether triplets are resolved over elements or site labels (default: `Element`)
    pub group_by: GroupBy,
    /// End-atom order exactly as the caller wrote it — `(A, C)` from `-a`/`-c` (or `-x`/`-z`).
    ///
    /// `r_cut_ab` is applied to A and `r_cut_bc` to C. Without this the two cutoffs would
    /// be handed out by the canonical (Z, symbol) order, which silently ignores the order
    /// the caller typed: `-a Zn -b O -c P --r-cut-ab 2.5 --r-cut-bc 2.0` would give 2.5 to
    /// P, because Z(P)=15 < Z(Zn)=30.
    ///
    /// `None` (the default, and the only option when scanning every triplet at once) falls
    /// back to the canonical order. Ignored when both ends are the same type — see
    /// `calc_angle`.
    pub ends: Option<(String, String)>,
}

impl Default for AngleParams {
    fn default() -> Self {
        AngleParams {
            r_cut_ab: 2.3,
            r_cut_bc: 2.3,
            angle_min: 0.0,
            angle_max: 180.0,
            d_angle: 0.1,
            group_by: GroupBy::default(),
            ends: None,
        }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Per-triplet angle statistics.
#[derive(Debug, Clone)]
pub struct AngleStats {
    /// Mean angle \[degrees\]
    pub mean: f64,
    /// Standard deviation \[degrees\]
    pub std: f64,
    /// Total number of angle measurements
    pub count: u64,
}

/// Result of a bond angle distribution calculation.
///
/// Key format: `"ElemA-ElemCenter-ElemC"` with Z(ElemA) ≤ Z(ElemC).
/// The histogram stores raw integer counts accumulated over all frames.
#[derive(Debug, Clone)]
pub struct AngleResult {
    /// Bin-centre angles \[degrees\]
    pub angle: Vec<f64>,
    /// Raw count histogram per triplet key
    pub hist: BTreeMap<String, Vec<u64>>,
    /// Per-triplet statistics (mean, std, count)
    pub stats: BTreeMap<String, AngleStats>,
    pub n_frames: usize,
    pub params: AngleParams,
    /// Types present (elements or site labels), sorted by (Z, string)
    pub elements: Vec<String>,
}

// ─── 计算 ────────────────────────────────────────────────────────────────────

/// Compute bond angle distribution for all element triplets (A-B-C, B = center).
///
/// All possible center elements B and endpoint pairs (A, C) are considered
/// automatically; no manual selection is needed.
///
/// Returns `None` if:
/// - The trajectory is empty
/// - Any frame is missing a cell
/// - No valid angle triplets are found
///
/// # Canonical key and counting convention
/// For a given B center, its neighbors are enumerated as unordered pairs
/// `(α, β)` with `α_global_index < β_global_index` to avoid double-counting.
/// The canonical key uses Z(ElemA) ≤ Z(ElemC).
///
/// Each geometric angle is therefore counted **once**: a PO₄ tetrahedron contributes 6
/// O-P-O angles, not 12. `code1/dump2analysis -m angle` instead enumerates ordered end
/// pairs and counts (O₁,P,O₂) and (O₂,P,O₁) separately, so its histogram is exactly
/// twice this one whenever the two ends are the same type (for A ≠ C the two agree).
/// The factor cancels out of `mean`, `std` and the peak positions.
///
/// # Cutoff assignment
/// `r_cut_ab` goes to end A and `r_cut_bc` to end C, where A/C are taken from
/// `params.ends` when the caller named the triplet, and from the canonical (Z, symbol)
/// order otherwise. When both ends are the same type, which one is "A" would depend on
/// the neighbour enumeration order, so `min(r_cut_ab, r_cut_bc)` is applied to both —
/// the only assignment that gives a reproducible result.
pub fn calc_angle(traj: &Trajectory, params: &AngleParams) -> Option<AngleResult> {
    if traj.frames.is_empty() { return None; }

    // 校验所有帧都有 cell
    if traj.frames.iter().any(|f| f.cell.is_none()) { return None; }

    if params.angle_max <= params.angle_min || params.d_angle <= 0.0 { return None; }
    let n_bins = ((params.angle_max - params.angle_min) / params.d_angle).ceil() as usize;
    if n_bins == 0 { return None; }
    let n_frames = traj.frames.len();

    // 类型列表（元素或位点标签），按 (Z, 字符串) 排序：字符串二级比较使同 Z 的
    // 多个标签（O_f / O_b_P_P）顺序可复现，而非随 HashSet 迭代序漂移。
    let first_frame = traj.frames.first()?;
    let by = params.group_by;
    let elements = sorted_types(first_frame, by);

    let max_rcut = params.r_cut_ab.max(params.r_cut_bc);

    // 并行逐帧计数：每个线程独立维护局部直方图，最后 reduce 合并
    let init = || BTreeMap::<String, Vec<u64>>::new();

    let total_hist: BTreeMap<String, Vec<u64>> = traj.frames.par_iter()
        .fold(init, |mut h, frame| {
            let cell = frame.cell.as_ref().unwrap();
            let n = frame.atoms.len();

            // 每帧构建 linked-cell list，将邻居搜索从 O(N²) 降至 O(N)
            let cl = CellList::build(frame, cell, max_rcut);

            for b_idx in 0..n {
                let b_elem = group_key(&frame.atoms[b_idx], by);

                // 只搜索 27 个相邻 cell，平均 ~30 个原子而非全部 N 个
                let neighbors = cl.neighbors_of(b_idx, max_rcut);

                // 枚举所有无序邻居对 (ai, ci)，ai < ci，避免重复计数
                let nn = neighbors.len();
                for ai in 0..nn {
                    for ci in (ai + 1)..nn {
                        let (a_idx, a_dist, a_vec) = neighbors[ai];
                        let (c_idx, c_dist, c_vec) = neighbors[ci];

                        let a_elem = group_key(&frame.atoms[a_idx], by);
                        let c_elem = group_key(&frame.atoms[c_idx], by);

                        // 规范排序：低 Z 端 → rcut_ab，高 Z 端 → rcut_bc
                        // 相同元素两端使用 min(rcut_ab, rcut_bc)
                        // 与 canonical_triplet 用同一套 (Z, 字符串) 比较，
                        // 否则同 Z 不同标签的两端谁拿 r_cut_ab 取决于邻居枚举顺序
                        let a_first = (elem_z(a_elem), a_elem) <= (elem_z(c_elem), c_elem);
                        let (lo_elem, lo_dist, lo_vec, hi_elem, hi_dist, hi_vec) =
                            if a_first {
                                (a_elem, a_dist, a_vec, c_elem, c_dist, c_vec)
                            } else {
                                (c_elem, c_dist, c_vec, a_elem, a_dist, a_vec)
                            };

                        let (rcut_lo, rcut_hi) = if lo_elem == hi_elem {
                            // 两端同类型时谁被排成 lo 取决于邻居枚举顺序，
                            // 只有取两者较小值才能给出可复现的结果
                            let m = params.r_cut_ab.min(params.r_cut_bc);
                            (m, m)
                        } else {
                            match &params.ends {
                                // 调用方点名了三元组：--r-cut-ab 归它写在前面的那一端。
                                // 仅当用户的 A 落在 hi 端时才需要互换。
                                Some((ea, ec))
                                    if hi_elem == ea.as_str() && lo_elem == ec.as_str() =>
                                {
                                    (params.r_cut_bc, params.r_cut_ab)
                                }
                                _ => (params.r_cut_ab, params.r_cut_bc),
                            }
                        };

                        if lo_dist >= rcut_lo || hi_dist >= rcut_hi { continue; }

                        // 计算夹角 A-B-C（vec_BA · vec_BC）
                        let dot = lo_vec[0]*hi_vec[0] + lo_vec[1]*hi_vec[1] + lo_vec[2]*hi_vec[2];
                        let cos_a = (dot / (lo_dist * hi_dist)).clamp(-1.0, 1.0);
                        let angle_deg = cos_a.acos().to_degrees();
                        // 直方图区间可由调用方收窄，区间外的角丢弃。上界取闭区间：
                        // 完全共线的三元组 acos(-1) 恰为 180.0，半开区间会把它整个丢掉
                        if angle_deg < params.angle_min || angle_deg > params.angle_max {
                            continue;
                        }
                        let bin = (((angle_deg - params.angle_min) / params.d_angle) as usize)
                            .min(n_bins - 1);

                        // lo_elem 已是低 Z 端，直接调用 canonical_triplet 确保一致性
                        let key = canonical_triplet(lo_elem, b_elem, hi_elem);
                        h.entry(key).or_insert_with(|| vec![0u64; n_bins])[bin] += 1;
                    }
                }
            }
            h
        })
        .reduce(init, |mut a, b| {
            for (key, b_hist) in b {
                let entry = a.entry(key).or_insert_with(|| vec![0u64; n_bins]);
                for (x, y) in entry.iter_mut().zip(b_hist.iter()) { *x += y; }
            }
            a
        });

    if total_hist.is_empty() { return None; }

    // 角度轴（bin 中心）
    let angle: Vec<f64> = (0..n_bins)
        .map(|i| params.angle_min + (i as f64 + 0.5) * params.d_angle)
        .collect();

    // 从直方图计算 mean 和 std（加权平均）
    let mut stats: BTreeMap<String, AngleStats> = BTreeMap::new();
    for (key, bins) in &total_hist {
        let total_count: u64 = bins.iter().sum();
        if total_count == 0 { continue; }
        let tc = total_count as f64;
        let mean: f64 = angle.iter().zip(bins.iter())
            .map(|(&a, &c)| a * c as f64).sum::<f64>() / tc;
        let var: f64 = angle.iter().zip(bins.iter())
            .map(|(&a, &c)| (a - mean).powi(2) * c as f64).sum::<f64>() / tc;
        stats.insert(key.clone(), AngleStats { mean, std: var.sqrt(), count: total_count });
    }

    Some(AngleResult { angle, hist: total_hist, stats, n_frames, params: params.clone(), elements })
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

/// Write bond angle distribution to a tab-separated text file (`.angle`).
///
/// Each column (after the angle axis) contains the raw count histogram for one
impl AngleResult {
    /// Projects the histogram into the long-format table the writers consume.
    ///
    /// One row per `(angle, triplet)`: columns `angle, end_a, center, end_c, count, p`.
    /// Same reasoning as `GrResult::to_tables` — the triplet becomes data, not column
    /// names, so trajectories with different compositions stack without alignment.
    ///
    /// Both a raw `count` and a normalised `p` are emitted. `count` is kept because the
    /// integer histogram is what `scripts/compare_angle.py` checks against
    /// dump2analysis bin for bin; `p` is what any plot actually wants.
    pub fn to_tables(&self) -> Vec<(String, Table)> {
        let keys = sort_triplet_keys(&self.hist);
        let n = self.angle.len();
        let rows = n * keys.len();

        let mut angle_col = Vec::with_capacity(rows);
        let mut a_col     = Vec::with_capacity(rows);
        let mut b_col     = Vec::with_capacity(rows);
        let mut c_col     = Vec::with_capacity(rows);
        let mut cnt_col   = Vec::with_capacity(rows);
        let mut p_col     = Vec::with_capacity(rows);

        for key in &keys {
            let (ea, centre, ec) = split_triplet(key);
            let h = self.hist.get(key);
            let total: f64 = h.map(|v| v.iter().sum::<u64>() as f64).unwrap_or(0.0);
            for i in 0..n {
                let c = h.map(|v| v[i]).unwrap_or(0);
                angle_col.push(self.angle[i]);
                a_col.push(ea.to_string());
                b_col.push(centre.to_string());
                c_col.push(ec.to_string());
                cnt_col.push(c as f64);
                p_col.push(if total > 0.0 { c as f64 / total } else { f64::NAN });
            }
        }

        let mut t = Table::new();
        t.push_num("angle", angle_col)
            .push_text("end_a", a_col)
            .push_text("center", b_col)
            .push_text("end_c", c_col)
            .push_num("count", cnt_col)
            .push_num("p", p_col);
        vec![("angle".to_string(), t)]
    }

    /// Parameter block plus the per-triplet mean/std/count summary.
    pub fn meta_lines(&self) -> Vec<String> {
        let (end_a, end_c) = match &self.params.ends {
            Some((a, c)) => (format!("end A = {a}"), format!("end C = {c}")),
            None => ("low-Z end".to_string(), "high-Z end".to_string()),
        };
        let mut v = vec![
            format!("r_cut_ab = {} Ang  ({end_a} to center)", self.params.r_cut_ab),
            format!("r_cut_bc = {} Ang  ({end_c} to center)", self.params.r_cut_bc),
            format!("angle    = {} .. {} deg", self.params.angle_min, self.params.angle_max),
            format!("d_angle  = {} deg", self.params.d_angle),
            "counting : each geometric angle once  \
             (code1/dump2analysis counts A-B-C twice when both ends are the same type)"
                .to_string(),
            format!("grouped by {}:", match self.params.group_by {
                GroupBy::Element => "element",
                GroupBy::Label => "label",
            }),
        ];
        for elem in &self.elements { v.push(format!("  {elem}")); }
        v.push("[statistics]".to_string());
        for key in sort_triplet_keys(&self.hist) {
            if let Some(s) = self.stats.get(&key) {
                v.push(format!(
                    "{key:<14}: mean={:7.3}  std={:6.3}  count={}",
                    s.mean, s.std, s.count
                ));
            }
        }
        v
    }
}

/// Split an `"A-B-C"` triplet key into its three type names.
fn split_triplet(key: &str) -> (&str, &str, &str) {
    match key.split_once('-').and_then(|(a, rest)| rest.split_once('-').map(|(b, c)| (a, b, c))) {
        Some(t) => t,
        None => (key, "", ""),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    /// 简单 O-Si-O 三原子分子（非周期），验证角度计算
    /// Si 在原点，两个 O 分别在 (2, 0, 0) 和 (0, 2, 0) → 角度 90°
    fn make_90deg_angle() -> Trajectory {
        // 用大盒子模拟非周期性（cell 存在但 pbc 等同于真空）
        let cell = Cell::from_lengths_angles(20.0, 20.0, 20.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Si", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("O",  Vector3::new(1.6, 0.0, 0.0)));  // O-Si 键长 1.6 Å
        frame.add_atom(Atom::new("O",  Vector3::new(0.0, 1.6, 0.0)));
        Trajectory::from_frame(frame)
    }

    /// 线形分子：中心原子和两端原子在一条线上 → 角度 180°
    fn make_180deg_angle() -> Trajectory {
        let cell = Cell::from_lengths_angles(20.0, 20.0, 20.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Si", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("O",  Vector3::new( 1.6, 0.0, 0.0)));
        frame.add_atom(Atom::new("O",  Vector3::new(-1.6, 0.0, 0.0)));
        Trajectory::from_frame(frame)
    }

    #[test]
    fn test_angle_90deg() {
        let traj = make_90deg_angle();
        let params = AngleParams { r_cut_ab: 2.0, r_cut_bc: 2.0, d_angle: 1.0, ..Default::default() };
        let res = calc_angle(&traj, &params).unwrap();

        // 键 "O-Si-O"：Z(O)=8 < Z(Si)=14，规范 key 为 "O-Si-O"
        let key = "O-Si-O";
        assert!(res.hist.contains_key(key), "missing key {}", key);

        let stats = &res.stats[key];
        assert!((stats.mean - 90.0).abs() < 1.0,
            "expected ~90°, got {:.2}°", stats.mean);
        assert_eq!(stats.count, 1, "should have exactly 1 angle pair");
    }

    #[test]
    fn test_angle_180deg() {
        let traj = make_180deg_angle();
        let params = AngleParams { r_cut_ab: 2.0, r_cut_bc: 2.0, d_angle: 1.0, ..Default::default() };
        let res = calc_angle(&traj, &params).unwrap();

        let stats = &res.stats["O-Si-O"];
        assert!((stats.mean - 180.0).abs() < 1.0,
            "expected ~180°, got {:.2}°", stats.mean);
        assert_eq!(stats.count, 1);
    }

    #[test]
    fn test_canonical_key_order() {
        // 无论 a, c 哪个先放，规范 key 应以低 Z 在前
        // O(Z=8), Si(Z=14)：canonical = "O-Si-O" 或 "Si-O-Si"，不是 "Si-O-O" 或 "O-Si-Si"
        assert_eq!(canonical_triplet("Si", "O", "Si"), "Si-O-Si");
        assert_eq!(canonical_triplet("Si", "O", "O"), "O-O-Si");  // Z(O)=8 < Z(Si)=14
        assert_eq!(canonical_triplet("O",  "Si", "O"), "O-Si-O");
        assert_eq!(canonical_triplet("P",  "O", "Si"), "Si-O-P"); // Z(Si)=14 < Z(P)=15
    }

    #[test]
    fn test_angle_no_double_counting() {
        // 90° 系统只有 1 个有效 (A, B, C) 三元组（A_idx < C_idx 约束）
        let traj = make_90deg_angle();
        let params = AngleParams::default();
        let res = calc_angle(&traj, &params).unwrap();
        let total: u64 = res.hist["O-Si-O"].iter().sum();
        assert_eq!(total, 1, "should count exactly 1 unordered O-O pair around Si");
    }

    #[test]
    fn test_rcut_filters_distant_atoms() {
        // Si 在原点，一个 O 在 1.6 Å（近），另一个 O 在 4.0 Å（远）
        // 用 rcut=2.3 Å 时，远端 O 不在截断内，不应形成角度
        let cell = Cell::from_lengths_angles(20.0, 20.0, 20.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Si", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("O",  Vector3::new(1.6, 0.0, 0.0)));
        frame.add_atom(Atom::new("O",  Vector3::new(0.0, 4.0, 0.0))); // 超出截断
        let traj = Trajectory::from_frame(frame);

        let params = AngleParams { r_cut_ab: 2.3, r_cut_bc: 2.3, d_angle: 1.0, ..Default::default() };
        let res = calc_angle(&traj, &params);

        // 只有 1 个 O 在截断内，无法形成三元组 → 结果为 None 或 O-Si-O 计数为 0
        let total: u64 = res.as_ref()
            .and_then(|r| r.hist.get("O-Si-O"))
            .map(|h| h.iter().sum())
            .unwrap_or(0);
        assert_eq!(total, 0, "far O should be excluded by rcut");
    }

    #[test]
    fn test_sort_triplet_keys_order() {
        // 排序应以 (Z_center, Z_left, Z_right) 为键
        // O-O-O (Z_center=8), O-O-Si (Z_center=8), O-Si-O (Z_center=14)
        let mut map: BTreeMap<String, Vec<u64>> = BTreeMap::new();
        map.insert("O-Si-O".to_string(), vec![0]);
        map.insert("O-O-O".to_string(),  vec![0]);
        map.insert("O-O-Si".to_string(), vec![0]);
        let keys = sort_triplet_keys(&map);
        // 期望：O-O-O, O-O-Si (Z_center=8), O-Si-O (Z_center=14)
        assert_eq!(keys[0].as_str(), "O-O-O");
        assert_eq!(keys[1].as_str(), "O-O-Si");
        assert_eq!(keys[2].as_str(), "O-Si-O");
    }

    #[test]
    fn test_to_tables_long_format_splits_triplet() {
        let traj = make_90deg_angle();
        let res = calc_angle(&traj, &AngleParams::default()).unwrap();
        let (name, t) = res.to_tables().remove(0);
        assert_eq!(name, "angle");
        assert_eq!(t.names(), vec!["angle", "end_a", "center", "end_c", "count", "p"]);
        assert_eq!(t.n_rows(), res.angle.len() * res.hist.len());
        assert!(t.validate().is_ok());
    }

    #[test]
    fn test_to_tables_probability_normalises_per_triplet() {
        let traj = make_90deg_angle();
        let res = calc_angle(&traj, &AngleParams::default()).unwrap();
        let (_, t) = res.to_tables().remove(0);
        let (counts, ps) = match (t.column("count").unwrap(), t.column("p").unwrap()) {
            (ferro_core::Column::Num(c), ferro_core::Column::Num(p)) => (c.clone(), p.clone()),
            _ => panic!("count/p must be numeric"),
        };
        // 每个三元组段内 p 之和为 1（有计数时）
        let n = res.angle.len();
        for chunk in 0..res.hist.len() {
            let lo = chunk * n;
            let total: f64 = counts[lo..lo + n].iter().sum();
            if total > 0.0 {
                let psum: f64 = ps[lo..lo + n].iter().sum();
                assert!((psum - 1.0).abs() < 1e-12, "p 未按三元组归一: {psum}");
            }
        }
    }

    #[test]
    fn test_meta_lines_keep_stats_and_counting_convention() {
        let traj = make_90deg_angle();
        let res = calc_angle(&traj, &AngleParams::default()).unwrap();
        let meta = res.meta_lines().join("\n");
        assert!(meta.contains("[statistics]"));
        assert!(meta.contains("mean="));
        assert!(meta.contains("each geometric angle once"), "计数约定差异必须留在头里");
    }

    #[test]
    fn test_angle_group_by_label_resolves_sites() {
        // 同一元素的两个位点标签应产生独立的三元组 key
        let cell = Cell::from_lengths_angles(20.0, 20.0, 20.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        let mut push = |elem: &str, label: Option<&str>, x: f64, y: f64| {
            let mut a = Atom::new(elem, Vector3::new(x, y, 0.0));
            a.label = label.map(|s| s.to_string());
            frame.add_atom(a);
        };
        // O_f — P — O_b_P_P，两条键都在 2.3 Å 内
        push("P", Some("P_0"), 0.0, 0.0);
        push("O", Some("O_f"), 1.6, 0.0);
        push("O", Some("O_b_P_P"), 0.0, 1.6);
        let traj = Trajectory::from_frame(frame);

        let by_elem = calc_angle(&traj, &AngleParams {
            group_by: GroupBy::Element, ..Default::default()
        }).unwrap();
        assert!(by_elem.hist.contains_key("O-P-O"), "element mode key");

        let by_label = calc_angle(&traj, &AngleParams {
            group_by: GroupBy::Label, ..Default::default()
        }).unwrap();
        assert!(
            by_label.hist.contains_key("O_b_P_P-P_0-O_f"),
            "label mode keys: {:?}", by_label.hist.keys().collect::<Vec<_>>()
        );
        assert_eq!(by_label.elements, vec!["O_b_P_P", "O_f", "P_0"]);
    }

    #[test]
    fn test_angle_same_z_label_order_is_deterministic() {
        // 同 Z 的多个标签必须每次给出相同的类型顺序与三元组 key
        let cell = Cell::from_lengths_angles(20.0, 20.0, 20.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        let mut push = |elem: &str, label: &str, x: f64, y: f64| {
            let mut a = Atom::new(elem, Vector3::new(x, y, 0.0));
            a.label = Some(label.to_string());
            frame.add_atom(a);
        };
        push("P", "P_0", 0.0, 0.0);
        push("O", "O_x", 1.6, 0.0);
        push("O", "O_f", 0.0, 1.6);
        push("O", "O_b_P_P", -1.6, 0.0);
        let traj = Trajectory::from_frame(frame);
        let params = AngleParams { group_by: GroupBy::Label, ..Default::default() };

        let expected_types = vec!["O_b_P_P", "O_f", "O_x", "P_0"];
        let first = calc_angle(&traj, &params).unwrap();
        let expected_keys: Vec<String> = first.hist.keys().cloned().collect();
        for _ in 0..8 {
            let res = calc_angle(&traj, &params).unwrap();
            assert_eq!(res.elements, expected_types);
            assert_eq!(res.hist.keys().cloned().collect::<Vec<_>>(), expected_keys);
        }
    }

    // ── Q6: cutoff 按调用方给的端原子顺序分配 ────────────────────────────────

    /// Zn-O-P：中心 O，Zn 在 2.4 Å、P 在 2.0 Å。
    ///
    /// Z(P)=15 < Z(Zn)=30，所以规范顺序里 P 是 lo 端、Zn 是 hi 端 —— 与用户写的
    /// `-a Zn -c P` 正好相反。取 r_cut_ab=2.5 / r_cut_bc=2.2 时两种分配给出
    /// 完全不同的结果，足以把这个行为钉死。
    fn make_zn_o_p() -> Trajectory {
        let cell = Cell::from_lengths_angles(20.0, 20.0, 20.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("O",  Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Zn", Vector3::new(2.4, 0.0, 0.0)));
        frame.add_atom(Atom::new("P",  Vector3::new(0.0, 2.0, 0.0)));
        Trajectory::from_frame(frame)
    }

    #[test]
    fn test_cutoff_follows_caller_end_order() {
        let traj = make_zn_o_p();
        let base = AngleParams { r_cut_ab: 2.5, r_cut_bc: 2.2, d_angle: 1.0, ..Default::default() };

        // 用户写 -a Zn -c P：Zn 拿 2.5（2.4 < 2.5 ✓），P 拿 2.2（2.0 < 2.2 ✓）→ 计入
        let named = AngleParams {
            ends: Some(("Zn".to_string(), "P".to_string())), ..base.clone()
        };
        let res = calc_angle(&traj, &named).unwrap();
        assert_eq!(res.stats["P-O-Zn"].count, 1, "端原子顺序应由 ends 决定");

        // 反过来写 -a P -c Zn：P 拿 2.5（✓），Zn 拿 2.2（2.4 > 2.2 ✗）→ 落选
        let flipped = AngleParams {
            ends: Some(("P".to_string(), "Zn".to_string())), ..base.clone()
        };
        assert!(calc_angle(&traj, &flipped).is_none(),
            "调换 -a/-c 应当改变 cutoff 归属");

        // 未点名三元组时退回规范 (Z, 符号) 顺序：lo=P 拿 2.5，hi=Zn 拿 2.2 → 落选
        assert!(calc_angle(&traj, &base).is_none(),
            "ends=None 时应沿用规范顺序");
    }

    #[test]
    fn test_symmetric_ends_ignore_order() {
        // 两端同类型时 ends 不生效，始终取 min(r_cut_ab, r_cut_bc)：
        // 谁被排成 lo 端取决于邻居枚举顺序，否则结果不可复现
        let traj = make_90deg_angle();
        let p = AngleParams {
            r_cut_ab: 2.0, r_cut_bc: 1.5, d_angle: 1.0,
            ends: Some(("O".to_string(), "O".to_string())), ..Default::default()
        };
        // 两个 O 都在 1.6 Å，min(2.0, 1.5) = 1.5 → 全部落选
        assert!(calc_angle(&traj, &p).is_none());
    }

    // ── Q3: 直方图上下界可调 ─────────────────────────────────────────────────

    #[test]
    fn test_angle_window_is_configurable() {
        let traj = make_90deg_angle();   // 唯一一个角是 90°

        // 收窄到 [80, 100)：角度轴从 80 起、bin 中心带 angle_min 偏移
        let inside = AngleParams {
            angle_min: 80.0, angle_max: 100.0, d_angle: 1.0,
            r_cut_ab: 2.0, r_cut_bc: 2.0, ..Default::default()
        };
        let res = calc_angle(&traj, &inside).unwrap();
        assert_eq!(res.angle.len(), 20);
        assert!((res.angle[0] - 80.5).abs() < 1e-12, "bin 中心应为 angle_min + (i+0.5)*d");
        assert_eq!(res.stats["O-Si-O"].count, 1);
        assert!((res.stats["O-Si-O"].mean - 90.0).abs() < 1.0);

        // 窗口挪开后该角被丢弃，没有任何三元组留下
        let outside = AngleParams {
            angle_min: 100.0, angle_max: 180.0, d_angle: 1.0,
            r_cut_ab: 2.0, r_cut_bc: 2.0, ..Default::default()
        };
        assert!(calc_angle(&traj, &outside).is_none());
    }

    #[test]
    fn test_180deg_stays_in_last_bin() {
        // acos(-1) 恰好是 180.0；上界若按半开区间处理会把共线三元组整个丢掉
        let traj = make_180deg_angle();
        let p = AngleParams { r_cut_ab: 2.0, r_cut_bc: 2.0, d_angle: 1.0, ..Default::default() };
        let res = calc_angle(&traj, &p).unwrap();
        let bins = &res.hist["O-Si-O"];
        assert_eq!(bins[bins.len() - 1], 1, "180° 应落在最后一个 bin");
    }

    #[test]
    fn test_invalid_window_rejected() {
        let traj = make_90deg_angle();
        for p in [
            AngleParams { angle_min: 100.0, angle_max: 100.0, ..Default::default() },
            AngleParams { angle_min: 100.0, angle_max: 50.0, ..Default::default() },
            AngleParams { d_angle: 0.0, ..Default::default() },
        ] {
            assert!(calc_angle(&traj, &p).is_none());
        }
    }
}
