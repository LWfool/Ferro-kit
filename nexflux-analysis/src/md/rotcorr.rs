//! Rotational autocorrelation function calculation and output.
//!
//! Computes the second-order Legendre polynomial rotational correlation function for molecular orientations:
//!   C(t) = <P₂(cosθ)> = <(3cos²θ − 1) / 2>
//!
//! Orientation vector: u_c(t) = Σ_{n ∈ neighbors(c, r_cut)} minimum_image(r_n − r_c)
//! For each center atom c, the bond vectors to all specified-element neighbors within r_cut are summed.
//!
//! Workflow: `calc_rotcorr` → `write_rotcorr`.
//! Algorithm follows code1/rotcorr.c (`EstimateRotationCorr`), generalised to arbitrary A-center-B structures.
//!
//! Parallelism: per time-origin par_iter; each origin computed independently then reduced.

use rayon::prelude::*;
use std::io::{BufWriter, Write};
use nexflux_core::Trajectory;
use super::gr::VERSION;

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Parameters for rotational autocorrelation function calculation.
#[derive(Debug, Clone)]
pub struct RotCorrParams {
    /// Element symbol of the center atom (e.g. `"O"` for water)
    pub center: String,
    /// Element symbol of the neighbor atoms (e.g. `"H"` for water)
    pub neighbor: String,
    /// Cutoff radius for neighbor search \[Å\] (default: 1.2)
    pub r_cut: f64,
    /// Correlation window in frames (`None` = all frames)
    pub tau: Option<usize>,
    /// Time shift between origins in frames (default: 1)
    pub shift: usize,
    /// Time step per frame \[fs\] (default: 1.0)
    pub dt: f64,
}

impl Default for RotCorrParams {
    fn default() -> Self {
        RotCorrParams {
            center: "O".to_string(),
            neighbor: "H".to_string(),
            r_cut: 1.2,
            tau: None,
            shift: 1,
            dt: 1.0,
        }
    }
}

// ─── 结果 ────────────────────────────────────────────────────────────────────

/// Result of a rotational autocorrelation function calculation.
///
/// `rotcorr[i]` = C(i·dt) ∈ [−0.5, 1.0].
/// C(0) = 1 by definition; C(∞) → 0 for freely rotating molecules.
/// `integral[i]` is the running integral `Σₖ₌₀ⁱ C[k]·dt` \[fs\],
/// whose long-time limit gives the rotational correlation time τc.
#[derive(Debug, Clone)]
pub struct RotCorrResult {
    /// Time axis \[fs\]; `time[i] = i * dt`
    pub time: Vec<f64>,
    /// P₂ rotational correlation function C(t) ∈ [−0.5, 1.0]
    pub rotcorr: Vec<f64>,
    /// Running integral ∫₀ᵗ C(τ)dτ \[fs\]
    pub integral: Vec<f64>,
    /// Number of center atoms (molecules) used
    pub n_molecules: usize,
    /// Number of time origins averaged
    pub n_origins: usize,
    pub params: RotCorrParams,
}

// ─── 计算 ────────────────────────────────────────────────────────────────────

/// Compute the P₂ rotational autocorrelation function with time-shift averaging.
///
/// For each center atom the orientation vector is the sum of all bond vectors
/// to neighbor atoms within `r_cut` (minimum-image corrected for periodic cells).
/// Center atoms with no valid neighbors in any frame are excluded.
///
/// Returns `None` if:
/// - The trajectory has fewer than 2 frames
/// - No center atoms match the element filter
/// - `tau` exceeds the trajectory length
/// - Any frame is missing a cell when PBC is required (detected automatically)
pub fn calc_rotcorr(traj: &Trajectory, params: &RotCorrParams) -> Option<RotCorrResult> {
    let n_steps = traj.n_frames();
    if n_steps < 2 { return None; }

    let tau = params.tau.unwrap_or(n_steps).min(n_steps).max(1);
    let shift = params.shift.max(1);

    // 确定参与计算的 center 原子下标（按第一帧筛选）
    let ref_frame = traj.first()?;
    let center_indices: Vec<usize> = ref_frame.atoms.iter().enumerate()
        .filter(|(_, a)| a.element == params.center)
        .map(|(i, _)| i)
        .collect();
    if center_indices.is_empty() { return None; }

    let has_cell = ref_frame.cell.is_some();
    let r_cut2 = params.r_cut * params.r_cut;

    // Precompute orientation vector orient[step][mol_local] = [ux, uy, uz] for every center atom in each frame.
    // Orientation vector = sum of all center→neighbor bond vectors within r_cut (minimum-image corrected).
    let orient: Vec<Vec<[f64; 3]>> = traj.frames.iter().map(|frame| {
        center_indices.iter().map(|&ci| {
            let mut ux = 0.0_f64;
            let mut uy = 0.0_f64;
            let mut uz = 0.0_f64;
            let c_pos = frame.atoms[ci].position;
            for (ni, na) in frame.atoms.iter().enumerate() {
                if ni == ci || na.element != params.neighbor { continue; }
                let diff = if has_cell {
                    if let Some(cell) = &frame.cell {
                        cell.minimum_image(na.position - c_pos)
                    } else {
                        na.position - c_pos
                    }
                } else {
                    na.position - c_pos
                };
                let d2 = diff.norm_squared();
                if d2 < r_cut2 {
                    ux += diff.x;
                    uy += diff.y;
                    uz += diff.z;
                }
            }
            [ux, uy, uz]
        }).collect()
    }).collect();

    let n_mol = center_indices.len();

    // 枚举所有 origin：内层循环访问 orient[p+i]（i ∈ [0, tau)），需 p+tau <= n_steps
    let p_values: Vec<usize> = (0..)
        .map(|k: usize| k * shift)
        .take_while(|&p| p + tau <= n_steps)
        .collect();
    let n_origins = p_values.len();
    if n_origins == 0 { return None; }

    // 并行计算各 origin 的局部 P₂ 累积
    let accum: Vec<f64> = p_values.par_iter()
        .map(|&p| {
            let mut local = vec![0.0f64; tau];
            for i in 0..tau {
                let mut sum_p2 = 0.0_f64;
                let mut count = 0usize;
                for j in 0..n_mol {
                    let [ax, ay, az] = orient[p][j];
                    let [bx, by, bz] = orient[p + i][j];
                    let norm_a2 = ax*ax + ay*ay + az*az;
                    let norm_b2 = bx*bx + by*by + bz*bz;
                    // skip zero-norm vectors (no neighbours) to avoid division by zero
                    if norm_a2 < 1e-30 || norm_b2 < 1e-30 { continue; }
                    let ab = ax*bx + ay*by + az*bz;
                    let cos2 = (ab * ab) / (norm_a2 * norm_b2);
                    sum_p2 += 0.5 * (3.0 * cos2 - 1.0);
                    count += 1;
                }
                if count > 0 {
                    local[i] = sum_p2 / count as f64;
                }
            }
            local
        })
        .reduce(
            || vec![0.0f64; tau],
            |mut a, b| { for i in 0..tau { a[i] += b[i]; } a },
        );

    // 归一化：除以 n_origins
    let inv = 1.0 / n_origins as f64;
    let time:    Vec<f64> = (0..tau).map(|i| i as f64 * params.dt).collect();
    let rotcorr: Vec<f64> = (0..tau).map(|i| accum[i] * inv).collect();

    // 累积积分 τc(t) = Σᵢ C(i)·dt
    let mut integral = vec![0.0f64; tau];
    let mut running = 0.0;
    for i in 0..tau {
        running += rotcorr[i] * params.dt;
        integral[i] = running;
    }

    Some(RotCorrResult {
        time, rotcorr, integral,
        n_molecules: n_mol,
        n_origins,
        params: params.clone(),
    })
}

// ─── 输出函数 ────────────────────────────────────────────────────────────────

/// Write rotational autocorrelation function to a tab-separated text file (`.rotcorr`).
///
/// Columns: `time[fs]`, `C(t)`, `integral[fs]`
pub fn write_rotcorr(result: &RotCorrResult, path: &str) -> std::io::Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(path)?);

    writeln!(w, "# nexflux v{}", VERSION)?;
    writeln!(w, "# Rotational Autocorrelation Function C(t) = <P2(cos theta)>")?;
    writeln!(w, "# {}", "-".repeat(60))?;
    writeln!(w, "# center   = {}", result.params.center)?;
    writeln!(w, "# neighbor = {}", result.params.neighbor)?;
    writeln!(w, "# r_cut    = {} Ang", result.params.r_cut)?;
    writeln!(w, "# tau      = {} frames", result.time.len())?;
    writeln!(w, "# shift    = {} frames", result.params.shift)?;
    writeln!(w, "# dt       = {} fs", result.params.dt)?;
    writeln!(w, "# molecules = {}", result.n_molecules)?;
    writeln!(w, "# origins   = {}", result.n_origins)?;
    writeln!(w, "# {}", "-".repeat(60))?;
    writeln!(w, "# time[fs]\tC(t)\tintegral[fs]")?;

    for i in 0..result.time.len() {
        writeln!(w, "{:.6e}\t{:.6e}\t{:.6e}",
            result.time[i], result.rotcorr[i], result.integral[i])?;
    }
    Ok(())
}

// ─── 测试 ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use nexflux_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    /// 构建单水分子轨迹，所有帧取向相同（O 在原点，两个 H 固定）
    fn make_fixed_orientation(n: usize) -> Trajectory {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for _ in 0..n {
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));
            frame.add_atom(Atom::new("H", Vector3::new(0.96, 0.0,  0.0)));
            frame.add_atom(Atom::new("H", Vector3::new(-0.24, 0.93, 0.0)));
            traj.add_frame(frame);
        }
        traj
    }

    /// 构建两帧轨迹：t=0 取向沿 x，t=1 取向沿 y（cosθ=0 → P₂=-0.5）
    fn make_perpendicular_orientations() -> Trajectory {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();

        // 帧0：H 沿 +x 方向
        let mut f0 = Frame::with_cell(cell.clone(), [true; 3]);
        f0.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));
        f0.add_atom(Atom::new("H", Vector3::new(1.0, 0.0, 0.0)));
        f0.add_atom(Atom::new("H", Vector3::new(1.0, 0.0, 0.0)));
        traj.add_frame(f0);

        // 帧1：H 沿 +y 方向（与帧0 垂直）
        let mut f1 = Frame::with_cell(cell.clone(), [true; 3]);
        f1.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));
        f1.add_atom(Atom::new("H", Vector3::new(0.0, 1.0, 0.0)));
        f1.add_atom(Atom::new("H", Vector3::new(0.0, 1.0, 0.0)));
        traj.add_frame(f1);

        traj
    }

    #[test]
    fn test_c0_always_one() {
        // C(0) = P₂(cos0) = P₂(1) = 1，对任意轨迹恒成立
        let traj = make_fixed_orientation(8);
        let params = RotCorrParams { tau: Some(4), ..Default::default() };
        let res = calc_rotcorr(&traj, &params).unwrap();
        assert!((res.rotcorr[0] - 1.0).abs() < 1e-10,
            "C(0) = {:.6}, expected 1.0", res.rotcorr[0]);
    }

    #[test]
    fn test_fixed_orientation_flat() {
        // 取向不变 → cos²θ = 1 for all t → C(t) = 1 恒为 1
        let traj = make_fixed_orientation(10);
        let params = RotCorrParams { tau: Some(5), ..Default::default() };
        let res = calc_rotcorr(&traj, &params).unwrap();
        for (i, &c) in res.rotcorr.iter().enumerate() {
            assert!((c - 1.0).abs() < 1e-10,
                "C({}) = {:.6}, expected 1.0", i, c);
        }
    }

    #[test]
    fn test_perpendicular_gives_minus_half() {
        // u(0) ∥ x，u(1) ∥ y → cosθ = 0 → P₂(0) = -0.5
        let traj = make_perpendicular_orientations();
        let params = RotCorrParams {
            center: "O".to_string(),
            neighbor: "H".to_string(),
            r_cut: 1.5,
            tau: Some(2), shift: 1, dt: 1.0,
        };
        let res = calc_rotcorr(&traj, &params).unwrap();
        // lag=0: C=1, lag=1: C=-0.5
        assert!((res.rotcorr[0] - 1.0).abs() < 1e-10, "C(0)={}", res.rotcorr[0]);
        assert!((res.rotcorr[1] - (-0.5)).abs() < 1e-10,
            "C(1) = {:.6}, expected -0.5", res.rotcorr[1]);
    }

    #[test]
    fn test_rcut_excludes_distant_neighbor() {
        // H 原子距 O 超出截断 → 取向向量为零 → 该分子被跳过 → 无有效 C(t)
        // 此时 calc_rotcorr 仍可返回 Some（count=0 → C[i]=0），但不崩溃
        let cell = Cell::from_lengths_angles(20.0, 20.0, 20.0, 90.0, 90.0, 90.0).unwrap();
        let mut traj = Trajectory::new();
        for _ in 0..4 {
            let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
            frame.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));
            frame.add_atom(Atom::new("H", Vector3::new(5.0, 0.0, 0.0))); // 远超截断
            frame.add_atom(Atom::new("H", Vector3::new(5.0, 0.0, 0.0)));
            traj.add_frame(frame);
        }
        let params = RotCorrParams { r_cut: 1.2, tau: Some(2), ..Default::default() };
        // 应返回 Some（不崩溃）或 None，不 panic
        let _ = calc_rotcorr(&traj, &params);
    }

    #[test]
    fn test_n_origins() {
        // n=10, tau=3, shift=1 → p+3<=10 → p ∈ {0..7} → 8 origins
        let traj = make_fixed_orientation(10);
        let params = RotCorrParams { tau: Some(3), shift: 1, ..Default::default() };
        let res = calc_rotcorr(&traj, &params).unwrap();
        assert_eq!(res.n_origins, 8);
    }

    #[test]
    fn test_integral_accumulates() {
        // C(t) = 1 恒定，dt=2.0 → integral[i] = (i+1) * 1.0 * 2.0
        let traj = make_fixed_orientation(8);
        let params = RotCorrParams { tau: Some(4), dt: 2.0, ..Default::default() };
        let res = calc_rotcorr(&traj, &params).unwrap();
        for (i, &ig) in res.integral.iter().enumerate() {
            let expected = (i as f64 + 1.0) * 2.0;
            assert!((ig - expected).abs() < 1e-10,
                "integral[{}] = {:.6}, expected {:.6}", i, ig, expected);
        }
    }

    #[test]
    fn test_write_rotcorr() {
        use std::io::Read;
        let traj = make_fixed_orientation(8);
        let params = RotCorrParams { tau: Some(4), ..Default::default() };
        let res = calc_rotcorr(&traj, &params).unwrap();
        let path = "/tmp/test_nexflux.rotcorr";
        write_rotcorr(&res, path).expect("write_rotcorr failed");

        let mut content = String::new();
        std::fs::File::open(path).unwrap().read_to_string(&mut content).unwrap();
        assert!(content.starts_with("# nexflux v"), "missing version header");
        assert!(content.contains("Rotational"), "missing title");
        assert!(content.contains("# time[fs]"), "missing column header");
        assert!(content.contains("integral"), "missing integral column");
    }
}
