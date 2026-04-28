//! Bond angle distribution calculation and output.
//!
//! Three-atom angle A-B-C with B as the center atom.
//! Automatically enumerates all chemical element triplets.
//!
//! Canonical key `"ElemA-ElemCenter-ElemC"`:
//!   - B is always the center; endpoints sorted by atomic number, Z(A) ≤ Z(C).
//!   - Symmetric pairs (e.g. Si-O-Si) are not double-counted: enumeration requires A_idx < C_idx.
//!
//! `rcut_ab`: cutoff from the lower-Z end atom (A) to center.
//! `rcut_bc`: cutoff from the higher-Z end atom (C) to center.
//! When both end atoms are the same element, `min(rcut_ab, rcut_bc)` is used.
//!
//! Parallelism: per-frame `par_iter().fold().reduce()`, same pattern as gr.rs.
//! Algorithm reference: code1/angle.c (`EstimateAngle`).

use nalgebra::{Matrix3, Vector3};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::io::{BufWriter, Write};
use ferro_core::Trajectory;
use super::gr::{VERSION, elem_z};

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

struct CellList {
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
    fn build(frame: &ferro_core::Frame, cell: &ferro_core::Cell, max_rcut: f64) -> Self {
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
    fn neighbors_of(&self, b_idx: usize, max_rcut: f64) -> Vec<(usize, f64, [f64; 3])> {
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
    /// Distance cutoff for the lower-Z end atom (A) to center (B) \[Å\] (default: 2.3)
    pub r_cut_ab: f64,
    /// Distance cutoff for the higher-Z end atom (C) to center (B) \[Å\] (default: 2.3)
    pub r_cut_bc: f64,
    /// Histogram bin width \[degrees\] (default: 0.1, same as code1 180/1800 split)
    pub d_angle: f64,
}

impl Default for AngleParams {
    fn default() -> Self {
        AngleParams { r_cut_ab: 2.3, r_cut_bc: 2.3, d_angle: 0.1 }
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
    /// Element types present, sorted by atomic number
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
/// The canonical key uses Z(ElemA) ≤ Z(ElemC); when both ends are the same
/// element the cutoff applied is `min(r_cut_ab, r_cut_bc)`.
pub fn calc_angle(traj: &Trajectory, params: &AngleParams) -> Option<AngleResult> {
    if traj.frames.is_empty() { return None; }

    // 校验所有帧都有 cell
    if traj.frames.iter().any(|f| f.cell.is_none()) { return None; }

    let n_bins = (180.0 / params.d_angle).ceil() as usize;
    let n_frames = traj.frames.len();

    // 元素列表（按 Z 排序，用于文件头）
    let first_frame = traj.frames.first()?;
    let mut elem_set = std::collections::HashSet::new();
    for a in &first_frame.atoms { elem_set.insert(a.element.clone()); }
    let mut elements: Vec<String> = elem_set.into_iter().collect();
    elements.sort_by_key(|e| elem_z(e));

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
                let b_elem = frame.atoms[b_idx].element.as_str();

                // 只搜索 27 个相邻 cell，平均 ~30 个原子而非全部 N 个
                let neighbors = cl.neighbors_of(b_idx, max_rcut);

                // 枚举所有无序邻居对 (ai, ci)，ai < ci，避免重复计数
                let nn = neighbors.len();
                for ai in 0..nn {
                    for ci in (ai + 1)..nn {
                        let (a_idx, a_dist, a_vec) = neighbors[ai];
                        let (c_idx, c_dist, c_vec) = neighbors[ci];

                        let a_elem = frame.atoms[a_idx].element.as_str();
                        let c_elem = frame.atoms[c_idx].element.as_str();

                        // 规范排序：低 Z 端 → rcut_ab，高 Z 端 → rcut_bc
                        // 相同元素两端使用 min(rcut_ab, rcut_bc)
                        let (lo_elem, lo_dist, lo_vec, hi_elem, hi_dist, hi_vec) =
                            if elem_z(a_elem) <= elem_z(c_elem) {
                                (a_elem, a_dist, a_vec, c_elem, c_dist, c_vec)
                            } else {
                                (c_elem, c_dist, c_vec, a_elem, a_dist, a_vec)
                            };

                        let (rcut_lo, rcut_hi) = if lo_elem == hi_elem {
                            let m = params.r_cut_ab.min(params.r_cut_bc);
                            (m, m)
                        } else {
                            (params.r_cut_ab, params.r_cut_bc)
                        };

                        if lo_dist >= rcut_lo || hi_dist >= rcut_hi { continue; }

                        // 计算夹角 A-B-C（vec_BA · vec_BC）
                        let dot = lo_vec[0]*hi_vec[0] + lo_vec[1]*hi_vec[1] + lo_vec[2]*hi_vec[2];
                        let cos_a = (dot / (lo_dist * hi_dist)).clamp(-1.0, 1.0);
                        let angle_deg = cos_a.acos().to_degrees();
                        let bin = ((angle_deg / params.d_angle) as usize).min(n_bins - 1);

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
        .map(|i| (i as f64 + 0.5) * params.d_angle)
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
/// triplet.  The header reports mean and standard deviation per triplet.
pub fn write_angle(result: &AngleResult, path: &str) -> std::io::Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(path)?);

    writeln!(w, "# ferro v{}", VERSION)?;
    writeln!(w, "# Bond Angle Distribution A-B-C  (B = center)")?;
    writeln!(w, "# {}", "-".repeat(60))?;
    writeln!(w, "# r_cut_ab = {} Ang  (low-Z end to center)", result.params.r_cut_ab)?;
    writeln!(w, "# r_cut_bc = {} Ang  (high-Z end to center)", result.params.r_cut_bc)?;
    writeln!(w, "# d_angle  = {} deg", result.params.d_angle)?;
    writeln!(w, "# frames   = {}", result.n_frames)?;
    writeln!(w, "# atoms:")?;
    for elem in &result.elements { writeln!(w, "#   {}", elem)?; }
    writeln!(w, "# [statistics]")?;
    let keys = sort_triplet_keys(&result.hist);
    for key in &keys {
        if let Some(s) = result.stats.get(key) {
            writeln!(w, "# {:<14}: mean={:7.3}  std={:6.3}  count={}",
                key, s.mean, s.std, s.count)?;
        }
    }
    writeln!(w, "# {}", "-".repeat(60))?;
    write!(w, "# angle[deg]")?;
    for k in &keys { write!(w, "\t{}", k)?; }
    writeln!(w)?;

    let n = result.angle.len();
    for i in 0..n {
        write!(w, "{:.4}", result.angle[i])?;
        for k in &keys {
            let v = result.hist.get(k).map(|h| h[i]).unwrap_or(0);
            write!(w, "\t{}", v)?;
        }
        writeln!(w)?;
    }
    Ok(())
}

// ─── 测试 ────────────────────────────────────────────────────────────────────

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
        let params = AngleParams { r_cut_ab: 2.0, r_cut_bc: 2.0, d_angle: 1.0 };
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
        let params = AngleParams { r_cut_ab: 2.0, r_cut_bc: 2.0, d_angle: 1.0 };
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

        let params = AngleParams { r_cut_ab: 2.3, r_cut_bc: 2.3, d_angle: 1.0 };
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
    fn test_write_angle() {
        use std::io::Read;
        let traj = make_90deg_angle();
        let res = calc_angle(&traj, &AngleParams::default()).unwrap();
        let path = "/tmp/test_ferro.angle";
        write_angle(&res, path).expect("write_angle failed");

        let mut content = String::new();
        std::fs::File::open(path).unwrap().read_to_string(&mut content).unwrap();
        assert!(content.starts_with("# ferro v"));
        assert!(content.contains("# angle[deg]"));
        assert!(content.contains("center"));
    }
}
