//! Charge-density cluster SDF (ChgSDF) — averages electron density over aligned Qn clusters.
//!
//! For each input `(Frame, ChargeGrid)` pair (one per QE pp.x cube file):
//!   1. Identify Qn clusters with the same logic as [`calc_cluster_sdf`].
//!   2. Extract a cubic sub-grid of the charge density centered on the cluster anchor.
//!   3. Apply the Kabsch rotation (same as atom SDF) to the sub-grid via pull trilinear interpolation.
//!   4. Accumulate; divide by cluster count to obtain the average density.
//!
//! Output: one [`CubeData`] per signature family, written as a Gaussian cube file.
//! Values are in the `ChargeGrid` convention (`ρ_phys × V_cell`); normalise by `V_cell` for `e/Å³`.
//!
//! See `dev/bader.md` and `cube_sdf.rs` for background.

use std::collections::HashMap;

use ferro_core::{Atom, ChargeGrid, CubeData, Frame};
use nalgebra::{Matrix3, Vector3};

use crate::md::cube_sdf::{
    ClusterSdfParams, ClusterSnapshot, align_to_reference, process_frame,
};

// ─── 公开参数 ─────────────────────────────────────────────────────────────────

/// Parameters for charge-density cluster SDF analysis.
#[derive(Debug, Clone)]
pub struct ChgSdfParams {
    /// Network-former element symbol, e.g. `"P"`
    pub former: String,
    /// Bridging-ligand element symbol, e.g. `"O"`
    pub ligand: String,
    /// Target cluster Qn level (0/1/2/3)
    pub target_qn: u8,
    /// Former–ligand cutoff radius \[Å\]
    pub former_ligand_cutoff: f64,
    /// Modifier cation element (e.g. `Some("Zn")`); `None` excludes modifiers
    pub modifier: Option<String>,
    /// Modifier–ligand cutoff radius \[Å\]
    pub modifier_cutoff: f64,
    /// Sub-grid boundary margin \[Å\] around the cluster anchor
    pub padding: f64,
    /// RMSD threshold above which a warning is printed \[Å\]
    pub rmsd_warn_threshold: f64,
}

impl Default for ChgSdfParams {
    fn default() -> Self {
        Self {
            former: "P".into(),
            ligand: "O".into(),
            target_qn: 2,
            former_ligand_cutoff: 2.4,
            modifier: None,
            modifier_cutoff: 2.8,
            padding: 6.0,
            rmsd_warn_threshold: 0.5,
        }
    }
}

// ─── 公开结果 ─────────────────────────────────────────────────────────────────

/// Alignment quality statistics.
#[derive(Debug, Clone, Default)]
pub struct ChgRmsdStats {
    pub mean: f64,
    pub max: f64,
    pub n_warned: usize,
}

/// Averaged charge density for one cluster-signature family.
#[derive(Debug)]
pub struct ChgSdfFamily {
    /// Signature string matching [`ClusterFamily`](crate::md::cube_sdf::ClusterFamily)
    pub signature: String,
    /// Averaged charge density in a Cartesian sub-grid centered on the cluster anchor
    pub cube: CubeData,
    /// Number of clusters accumulated (including reference)
    pub n_clusters: usize,
    /// Alignment RMSD statistics
    pub rmsd_stats: ChgRmsdStats,
}

/// Return type of [`calc_chg_sdf`].
pub struct ChgSdfResult {
    /// Signature → ChgSdfFamily
    pub families: HashMap<String, ChgSdfFamily>,
    /// Total frames processed (frames with a periodic cell)
    pub n_frames: usize,
    /// Total clusters across all families
    pub n_clusters_total: usize,
}

// ─── 内部累加器 ───────────────────────────────────────────────────────────────

struct ChgAcc {
    reference: ClusterSnapshot,
    /// Accumulated charge density (sum over all clusters)
    grid_sum: Vec<f64>,
    n: usize,
    n_clusters: usize,
    rmsd_sum: f64,
    rmsd_max: f64,
    n_warned: usize,
}

impl ChgAcc {
    fn new(reference: ClusterSnapshot, initial_sub: Vec<f64>, n: usize) -> Self {
        Self {
            reference,
            grid_sum: initial_sub,
            n,
            n_clusters: 1,
            rmsd_sum: 0.0,
            rmsd_max: 0.0,
            n_warned: 0,
        }
    }

    fn push(&mut self, rotated: Vec<f64>, rmsd: f64, threshold: f64) {
        for (a, s) in self.grid_sum.iter_mut().zip(&rotated) {
            *a += s;
        }
        self.n_clusters += 1;
        self.rmsd_sum += rmsd;
        if rmsd > self.rmsd_max { self.rmsd_max = rmsd; }
        if rmsd > threshold {
            self.n_warned += 1;
            eprintln!(
                "[ChgSDF] RMSD = {rmsd:.3} Å > 阈值 {threshold:.3} Å，\
                 团簇电子密度偏差可能较大，请人工甄别。"
            );
        }
    }
}

// ─── 主入口 ───────────────────────────────────────────────────────────────────

/// Compute averaged charge-density cluster SDF from a list of (Frame, ChargeGrid) pairs.
///
/// Each pair is one QE pp.x cube file (read via [`ferro_io::read_cube_as_chg`]).
/// Returns `None` if no target-Qn cluster is found.
pub fn calc_chg_sdf(
    pairs: &[(Frame, ChargeGrid)],
    params: &ChgSdfParams,
) -> Option<ChgSdfResult> {
    if pairs.is_empty() { return None; }

    // 从第一个 ChargeGrid 确定输出格分辨率
    let first_chg = &pairs[0].1;
    let voxel_size = (0..3)
        .map(|i| first_chg.lat2car.column(i).norm())
        .fold(f64::MAX, f64::min);
    let half_n = (params.padding / voxel_size).ceil() as usize;
    let n = 2 * half_n + 1;

    // ClusterSdfParams wrapper for process_frame
    let sdf_params = ClusterSdfParams {
        former: params.former.clone(),
        ligand: params.ligand.clone(),
        target_qn: params.target_qn,
        former_ligand_cutoff: params.former_ligand_cutoff,
        modifier: params.modifier.clone(),
        modifier_cutoff: params.modifier_cutoff,
        grid_res: 0.1,
        sigma: 1.5,
        padding: params.padding,
        rmsd_warn_threshold: params.rmsd_warn_threshold,
    };

    let mut accumulators: HashMap<String, ChgAcc> = HashMap::new();
    let mut n_frames = 0usize;

    for (frame, chg) in pairs {
        let cell = match frame.cell.as_ref() {
            Some(c) => c,
            None => continue,
        };
        n_frames += 1;

        for mut snapshot in process_frame(frame, cell, &sdf_params) {
            let sig = snapshot.signature();
            let center = snapshot.anchor_global_pos;

            if let Some(acc) = accumulators.get_mut(&sig) {
                let (rmsd, rot) = align_to_reference(&mut snapshot, &acc.reference, params.target_qn);
                let sub = extract_subgrid(chg, center, half_n, voxel_size);
                let rotated = rotate_grid(&sub, n, rot);
                acc.push(rotated, rmsd, params.rmsd_warn_threshold);
            } else {
                let sub = extract_subgrid(chg, center, half_n, voxel_size);
                accumulators.insert(sig, ChgAcc::new(snapshot, sub, n));
            }
        }
    }

    if accumulators.is_empty() { return None; }

    let n_clusters_total = accumulators.values().map(|a| a.n_clusters).sum();
    let families = accumulators
        .into_iter()
        .map(|(sig, acc)| {
            let family = build_family(acc, voxel_size, half_n);
            (sig, family)
        })
        .collect();

    Some(ChgSdfResult { families, n_frames, n_clusters_total })
}

// ─── 子格提取 ─────────────────────────────────────────────────────────────────

/// 从 ChargeGrid 提取以 `center`（Å）为中心的立方子格，保持原始 cube 分辨率。
///
/// 输出索引约定与 CubeData 相同：`idx = ix*(n*n) + iy*n + iz`（X 最慢，Z 最快）。
fn extract_subgrid(
    chg: &ChargeGrid,
    center: Vector3<f64>,
    half_n: usize,
    voxel_size: f64,
) -> Vec<f64> {
    let n = 2 * half_n + 1;
    let mut sub = vec![0.0_f64; n * n * n];
    let half = half_n as f64;

    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let r_cart = center + Vector3::new(
                    (ix as f64 - half) * voxel_size,
                    (iy as f64 - half) * voxel_size,
                    (iz as f64 - half) * voxel_size,
                );
                let r_lat = chg.car2lat * r_cart;
                sub[ix * n * n + iy * n + iz] = trilinear_interp_pbc(chg, r_lat);
            }
        }
    }
    sub
}

/// 对 ChargeGrid 做 PBC 三线性插值，输入为格点坐标（可为分数值）。
fn trilinear_interp_pbc(chg: &ChargeGrid, r_lat: Vector3<f64>) -> f64 {
    let [n1, n2, n3] = chg.shape;

    // PBC 折叠至 [0, Ni)
    let lx = r_lat.x.rem_euclid(n1 as f64);
    let ly = r_lat.y.rem_euclid(n2 as f64);
    let lz = r_lat.z.rem_euclid(n3 as f64);

    let x0 = lx.floor() as usize;
    let y0 = ly.floor() as usize;
    let z0 = lz.floor() as usize;
    let x1 = (x0 + 1) % n1;
    let y1 = (y0 + 1) % n2;
    let z1 = (z0 + 1) % n3;
    let fx = lx - lx.floor();
    let fy = ly - ly.floor();
    let fz = lz - lz.floor();

    let idx = |i: usize, j: usize, k: usize| chg.rho[i + n1 * (j + n2 * k)];

    let c00 = idx(x0, y0, z0) * (1.0 - fx) + idx(x1, y0, z0) * fx;
    let c10 = idx(x0, y1, z0) * (1.0 - fx) + idx(x1, y1, z0) * fx;
    let c01 = idx(x0, y0, z1) * (1.0 - fx) + idx(x1, y0, z1) * fx;
    let c11 = idx(x0, y1, z1) * (1.0 - fx) + idx(x1, y1, z1) * fx;

    let c0 = c00 * (1.0 - fy) + c10 * fy;
    let c1 = c01 * (1.0 - fy) + c11 * fy;

    c0 * (1.0 - fz) + c1 * fz
}

// ─── 子格旋转 ─────────────────────────────────────────────────────────────────

/// 用 pull 插值旋转立方子格：对输出格每个体素，逆变换回源格坐标取值。
///
/// `rot` 满足 `rot * mobile_atom ≈ reference_atom`，即把移动帧对齐到参考帧的旋转。
/// Pull: `output[r] = input[rot^T * r]`（R^T = R^{-1} 对正交旋转矩阵）。
fn rotate_grid(sub: &[f64], n: usize, rot: Matrix3<f64>) -> Vec<f64> {
    let rot_inv = rot.transpose();
    let half = (n as f64 - 1.0) * 0.5;
    let mut out = vec![0.0_f64; n * n * n];

    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let r_out = Vector3::new(
                    ix as f64 - half,
                    iy as f64 - half,
                    iz as f64 - half,
                );
                let r_src = rot_inv * r_out + Vector3::repeat(half);
                out[ix * n * n + iy * n + iz] = trilinear_interp_subgrid(sub, n, r_src);
            }
        }
    }
    out
}

/// 在子格（无 PBC）中三线性插值，越界返回 0.0。
fn trilinear_interp_subgrid(sub: &[f64], n: usize, pos: Vector3<f64>) -> f64 {
    let ni = n as isize;
    let x0 = pos.x.floor() as isize;
    let y0 = pos.y.floor() as isize;
    let z0 = pos.z.floor() as isize;
    let fx = pos.x - pos.x.floor();
    let fy = pos.y - pos.y.floor();
    let fz = pos.z - pos.z.floor();

    let get = |ix: isize, iy: isize, iz: isize| -> f64 {
        if ix < 0 || iy < 0 || iz < 0 || ix >= ni || iy >= ni || iz >= ni {
            return 0.0;
        }
        sub[ix as usize * n * n + iy as usize * n + iz as usize]
    };

    let c00 = get(x0, y0, z0) * (1.0 - fx) + get(x0 + 1, y0, z0) * fx;
    let c10 = get(x0, y0 + 1, z0) * (1.0 - fx) + get(x0 + 1, y0 + 1, z0) * fx;
    let c01 = get(x0, y0, z0 + 1) * (1.0 - fx) + get(x0 + 1, y0, z0 + 1) * fx;
    let c11 = get(x0, y0 + 1, z0 + 1) * (1.0 - fx) + get(x0 + 1, y0 + 1, z0 + 1) * fx;

    let c0 = c00 * (1.0 - fy) + c10 * fy;
    let c1 = c01 * (1.0 - fy) + c11 * fy;

    c0 * (1.0 - fz) + c1 * fz
}

// ─── 结果构建 ─────────────────────────────────────────────────────────────────

fn build_family(acc: ChgAcc, voxel_size: f64, half_n: usize) -> ChgSdfFamily {
    let sig = acc.reference.signature();
    let n = acc.n;
    let n_clusters = acc.n_clusters;

    // 平均电荷密度
    let avg: Vec<f64> = acc.grid_sum.iter().map(|v| v / n_clusters as f64).collect();

    // 输出格原点在局部坐标系（锚原子在 [half_n, half_n, half_n] 即原点）
    let origin = Vector3::repeat(-(half_n as f64) * voxel_size);
    let spacing = Matrix3::from_diagonal(&Vector3::repeat(voxel_size));

    // cube 文件头：参考团簇原子（局部坐标，锚在原点）
    let ref_frame = build_ref_frame(&acc.reference, origin);

    let rmsd_stats = ChgRmsdStats {
        mean: if n_clusters > 1 { acc.rmsd_sum / (n_clusters - 1) as f64 } else { 0.0 },
        max: acc.rmsd_max,
        n_warned: acc.n_warned,
    };

    let cube = CubeData {
        frame: ref_frame,
        data: avg,
        shape: [n, n, n],
        origin,
        spacing,
    };

    ChgSdfFamily { signature: sig, cube, n_clusters, rmsd_stats }
}

/// 构建参考团簇的 Frame（用于 cube 文件原子头），坐标为局部笛卡尔（锚在原点）。
fn build_ref_frame(snapshot: &ClusterSnapshot, _origin: Vector3<f64>) -> Frame {
    use crate::md::cube_sdf::type_to_element;
    let mut frame = Frame::new();
    for (t, &pos) in snapshot.types.iter().zip(&snapshot.positions) {
        frame.add_atom(Atom::new(type_to_element(t), pos));
    }
    frame
}

// ─── 测试 ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::Cell;
    use nalgebra::Matrix3;

    fn make_cubic_chg(n: usize, cell_a: f64) -> (Frame, ChargeGrid) {
        let cell = Cell::from_lengths_angles(cell_a, cell_a, cell_a, 90.0, 90.0, 90.0).unwrap();
        let vol = cell.volume();
        // 均匀密度 rho_stored = 1.0 * vol
        let rho = vec![vol; n * n * n];
        let chg = ChargeGrid::new(rho, [n, n, n], &cell);
        let frame = Frame::with_cell(cell, [true; 3]);
        (frame, chg)
    }

    #[test]
    fn test_trilinear_interp_pbc_uniform() {
        // 均匀密度：任意位置插值均等于常数
        let (_, chg) = make_cubic_chg(10, 10.0);
        let val = chg.rho[0];
        for _ in 0..5 {
            let r_lat = Vector3::new(3.7, 2.1, 8.9);
            let interp = trilinear_interp_pbc(&chg, r_lat);
            assert!((interp - val).abs() < 1e-10, "uniform density interpolation failed");
        }
    }

    #[test]
    fn test_extract_subgrid_uniform() {
        // 均匀密度子格：所有体素应等于常数
        let (_, chg) = make_cubic_chg(20, 10.0);
        let expected = chg.rho[0];
        let voxel_size = chg.lat2car.column(0).norm();
        let center = Vector3::new(5.0, 5.0, 5.0);
        let sub = extract_subgrid(&chg, center, 3, voxel_size);
        for (i, &v) in sub.iter().enumerate() {
            assert!((v - expected).abs() < 1e-8, "sub[{i}] = {v} != {expected}");
        }
    }

    #[test]
    fn test_rotate_grid_identity() {
        // 恒等旋转：旋转后子格与原子格完全相同
        let n = 5;
        let sub: Vec<f64> = (0..n * n * n).map(|i| i as f64).collect();
        let rotated = rotate_grid(&sub, n, Matrix3::identity());
        for (a, b) in sub.iter().zip(&rotated) {
            assert!((a - b).abs() < 1e-8, "identity rotation changed grid values");
        }
    }

    #[test]
    fn test_rotate_grid_180z() {
        // 绕 Z 轴旋转 180°：中心体素不变
        let n = 5;
        let mut sub = vec![0.0_f64; n * n * n];
        // 仅中心体素非零
        let mid = n / 2;
        sub[mid * n * n + mid * n + mid] = 1.0;

        let rot = Matrix3::new(
            -1.0, 0.0, 0.0,
             0.0,-1.0, 0.0,
             0.0, 0.0, 1.0,
        );
        let rotated = rotate_grid(&sub, n, rot);
        let center_val = rotated[mid * n * n + mid * n + mid];
        assert!(center_val > 0.9, "center voxel should survive 180° rotation, got {center_val}");
    }
}
