//! 参数取值范围的检查，供各 `XxxParams::validate` 共用（审查 M3）。
//!
//! 一律写成「合法则放行，否则报错」：NaN 与任何数比较都为假，落到报错分支，
//! `--dr nan` 也会被拦下；反过来写 `if v <= 0.0 { 报错 }` 会把 NaN 放过去。名字用 CLI 的参数拼法（`r-max`），
//! 与 `correlate::resolve_max_lag` 的 `max-lag` 一致。

use ferro_core::error::ChemError;
use ferro_core::{PeriodicCell, Trajectory};

type Check = ferro_core::Result<()>;

fn fail(msg: String) -> Check {
    Err(ChemError::ValidationError(msg))
}

/// [`Trajectory::check_same_atoms`], saying why a time-correlation analysis needs it.
pub(crate) fn same_atoms(traj: &Trajectory) -> Check {
    traj.check_same_atoms().map_err(|e| match e {
        ChemError::ValidationError(m) => ChemError::ValidationError(format!(
            "{m}; this analysis follows atoms by index and needs the same atoms in every frame"
        )),
        e => e,
    })
}

/// 步长、截断、时间步：有限且 > 0（`inf` 的步长会得到 0 个 bin）
pub(crate) fn positive(name: &str, v: f64) -> Check {
    if v > 0.0 && v.is_finite() { return Ok(()); }
    fail(format!("{name} must be a finite number > 0, got {v}"))
}

/// 有限且 >= 0（`--q-max inf` 之类在下游得到 0 个或无穷多个 bin）
pub(crate) fn non_negative(name: &str, v: f64) -> Check {
    if v >= 0.0 && v.is_finite() { return Ok(()); }
    fail(format!("{name} must be a finite number >= 0, got {v}"))
}

/// `lo < hi`，两端都有限
pub(crate) fn ordered(lo_name: &str, lo: f64, hi_name: &str, hi: f64) -> Check {
    if lo < hi && lo.is_finite() && hi.is_finite() { return Ok(()); }
    if !hi.is_finite() { return fail(format!("{hi_name} must be a finite number, got {hi}")); }
    if !lo.is_finite() { return fail(format!("{lo_name} must be a finite number, got {lo}")); }
    fail(format!("{lo_name} ({lo}) must be < {hi_name} ({hi})"))
}

/// 网格点数、间隔帧数这类计数
pub(crate) fn at_least_one(name: &str, n: usize) -> Check {
    if n == 0 { return fail(format!("{name} must be >= 1, got 0")); }
    Ok(())
}

/// 可省略的计数：给了就必须 >= 1（上限取决于帧数，由 `resolve_max_lag` 逐文件查）
pub(crate) fn at_least_one_if_given(name: &str, n: Option<usize>) -> Check {
    n.map_or(Ok(()), |n| at_least_one(name, n))
}

/// `[lo, hi)` 至少容得下一个宽 `step` 的 bin
pub(crate) fn holds_a_bin(range: &str, lo: f64, hi: f64, step: &str, w: f64) -> Check {
    if (hi - lo) / w >= 1.0 { return Ok(()); }
    fail(format!("{range} [{lo}, {hi}] is narrower than one {step} = {w} bin"))
}

/// Every frame's cell as a [`PeriodicCell`], `None` for a frame that has no cell.
///
/// A frame that carries a cell must be 3-D periodic with a finite, invertible
/// matrix; anything else is an error naming the frame. The analyses only know bulk
/// periodicity, and a slab or a molecule in a box used to be measured as one
/// (review A-7). Built once at the entry, the cells make the minimum image in the
/// parallel loops infallible.
///
/// With `cutoff = Some((name, r))`, `r` must not exceed the minimum-image bound
/// (half the smallest interplanar spacing): past it only the nearest image is seen
/// and neighbours in farther images are silently lost. Every frame is checked, not
/// just frame 0 — an NPT box shrinks, and a late frame past the bound is just as
/// wrong. The tightest frame is reported so the user learns the usable limit.
///
/// gr does not pass a cutoff: it clamps `r_max` to the bound, as the manual says.
pub(crate) fn periodic_cells(
    traj: &Trajectory,
    cutoff: Option<(&str, f64)>,
) -> ferro_core::Result<Vec<Option<PeriodicCell>>> {
    let cells = traj.frames.iter().enumerate().map(|(i, f)| {
        f.cell.as_ref().map(|c| PeriodicCell::new(c, f.pbc).map_err(|e| match e {
            ChemError::ValidationError(m) => ChemError::ValidationError(format!("frame {i}: {m}")),
            e => e,
        })).transpose()
    }).collect::<ferro_core::Result<Vec<_>>>()?;
    if let Some((name, r)) = cutoff {
        let tightest = cells.iter().enumerate()
            .filter_map(|(i, c)| c.as_ref().map(|c| (i, c.mic_bound())))
            .reduce(|a, b| if b.1 < a.1 { b } else { a });
        if let Some((i, bound)) = tightest {
            if r > bound {
                fail(format!(
                    "{name} {r:.3} A exceeds the minimum-image bound {bound:.3} A of the cell (frame {i})"
                ))?;
            }
        }
    }
    Ok(cells)
}

/// [`periodic_cells`] for an analysis that needs a cell in every frame.
pub(crate) fn all_periodic_cells(
    traj: &Trajectory,
    cutoff: Option<(&str, f64)>,
) -> ferro_core::Result<Vec<PeriodicCell>> {
    periodic_cells(traj, cutoff)?.into_iter().enumerate().map(|(i, c)| c.ok_or_else(|| {
        ChemError::ValidationError(format!("all frames must have a periodic cell (frame {i} has none)"))
    })).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nan_is_rejected_everywhere() {
        assert!(positive("dr", f64::NAN).is_err());
        assert!(non_negative("r-min", f64::NAN).is_err());
        assert!(ordered("r-min", 0.0, "r-max", f64::NAN).is_err());
        assert!(holds_a_bin("r", 0.0, 5.0, "dr", f64::NAN).is_err());
    }

    #[test]
    fn test_boundaries() {
        assert!(positive("dr", 0.0).is_err());
        assert!(positive("dr", 1e-300).is_ok());
        assert!(positive("dr", f64::INFINITY).is_err());
        assert!(non_negative("r-min", 0.0).is_ok());
        assert!(non_negative("r-min", -1e-12).is_err());
        // 审查 B-2：以前只有 positive 查 is_finite，这两个放行 +inf
        assert!(non_negative("r-min", f64::INFINITY).is_err());
        assert!(ordered("a", 0.0, "b", f64::INFINITY).is_err());
        assert!(ordered("a", f64::NEG_INFINITY, "b", 0.0).is_err());
        assert!(ordered("a", 1.0, "b", 1.0).is_err());
        assert!(ordered("a", 1.0, "b", 1.5).is_ok());
        assert!(at_least_one("nx", 0).is_err());
        assert!(at_least_one_if_given("max-lag", None).is_ok());
        assert!(at_least_one_if_given("max-lag", Some(0)).is_err());
        assert!(holds_a_bin("r", 0.0, 0.5, "dr", 1.0).is_err());
        assert!(holds_a_bin("r", 0.0, 1.0, "dr", 1.0).is_ok());
    }

    /// 每个参数结构体：默认值合法，每个受检字段各改一个坏值都被拦下。
    /// 盯的是「字段有没有接进 validate」，判据本身由上面两个测试管
    #[test]
    fn test_every_params_struct_checks_each_field() {
        use crate::dft::ChgSdfParams;
        use crate::md::*;
        let nan = f64::NAN;
        macro_rules! case {
            ($ty:ident, $( $field:ident : $bad:expr ),+ $(,)?) => {{
                assert!($ty::default().validate().is_ok(), "{} 的默认值应合法", stringify!($ty));
                $(
                    let p = $ty { $field: $bad, ..$ty::default() };
                    assert!(p.validate().is_err(), "{}.{} = {:?} 应被拦下",
                            stringify!($ty), stringify!($field), $bad);
                )+
            }};
        }
        case!(GrParams, r_min: -0.1, r_max: 0.0, dr: 0.0, dr: nan, dr: 20.0);
        let inf = f64::INFINITY;
        case!(SqParams, q_min: -1.0, q_min: 0.0, q_max: 0.0, q_max: inf, dq: 0.0, dq: nan);
        case!(MsdParams, dt: 0.0, dt: nan, max_lag: Some(0),
              fit_range: Some((0.8, 0.2)), fit_range: Some((0.0, 1.5)), fit_range: Some((nan, 0.5)));
        case!(AngleParams, r_cut_ab: 0.0, r_cut_bc: nan, angle_min: -1.0, angle_max: 0.0,
              angle_max: 181.0, d_angle: 0.0, d_angle: f64::INFINITY);
        case!(VacfParams, dt: -1.0, max_lag: Some(0));
        case!(RotCorrParams, r_cut: 0.0, dt: nan, max_lag: Some(0));
        case!(BondLifeParams, r_bond: 0.0, r_bond: nan, r_break: Some(1.0), dt: 0.0, max_lag: Some(0));
        case!(VanHoveParams, dt: 0.0, tau: Some(0), shift: 0, r_min: -1.0, dr: 0.0, r_max: 0.0,
              r_max: inf);
        case!(CubeDensityParams, nx: 0, ny: 0, nz: 0);
        case!(CubeJumpParams, nx: 0, ny: 0, nz: 0, tau: 0, threshold: 0.0, threshold: nan,
              threshold: inf);
        case!(CubeRadiusParams, nx: 0, ny: 0, nz: 0, radius: 0.0, radius: nan);
        case!(ClusterSdfParams, former_ligand_cutoff: 0.0, modifier_cutoff: nan, grid_res: 0.0,
              sigma: -1.0, sigma: inf, padding: -1.0, padding: inf, rmsd_warn_threshold: nan);
        case!(ChgSdfParams, former_ligand_cutoff: 0.0, modifier_cutoff: -1.0, padding: nan,
              rmsd_warn_threshold: -0.5);
    }

    /// 审查 B-1：ASE 的二维材料约定 c=0（pbc TTF）给出奇异胞。这几处在闭包里
    /// `expect` 可逆，以前直接 panic、打断整个批次（gr / msd 等同一文件是报错跳过）
    fn singular_traj() -> Trajectory {
        use ferro_core::{Atom, Cell, Frame};
        use nalgebra::{Matrix3, Vector3};
        let cell = Cell::from_matrix(Matrix3::new(10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0));
        let mut traj = Trajectory::new();
        for k in 0..3 {
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            f.add_atom(Atom::new("P", Vector3::new(5.0 + 0.01 * k as f64, 5.0, 5.0)));
            for d in [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]] {
                f.add_atom(Atom::new("O", Vector3::new(5.0 + 0.9 * d[0], 5.0 + 0.9 * d[1], 5.0 + 0.9 * d[2])));
            }
            traj.add_frame(f);
        }
        traj
    }

    #[test]
    fn test_singular_cells_are_an_error_not_a_panic() {
        use crate::md::*;
        let traj = singular_traj();
        let msg = calc_vanhove(&traj, &VanHoveParams::default()).err().map(|e| e.to_string());
        assert!(msg.as_deref().is_some_and(|m| m.contains("singular") && m.contains("frame 0")),
            "vanhove 应报奇异胞并点名帧号：{msg:?}");
        let msg = calc_cube_jump(&traj, &CubeJumpParams { tau: 1, ..Default::default() }).err().map(|e| e.to_string());
        assert!(msg.as_deref().is_some_and(|m| m.contains("singular")), "cube jump: {msg:?}");
        let params = ClusterSdfParams {
            former: "P".into(), ligand: "O".into(), target_qn: 0, ..Default::default()
        };
        let msg = calc_cluster_sdf(&traj, &params).err().map(|e| e.to_string());
        assert!(msg.as_deref().is_some_and(|m| m.contains("singular")), "map sdf 应报奇异胞：{msg:?}");
    }

    /// 四个按截断找邻居的分析都要拦下超过最小镜像上界的截断，并点名最紧的帧。
    /// 盒子 10 Å → 8 Å（上界 5 → 4 Å），模拟 NPT 后段收缩：只查第 0 帧会放过 4.5 Å
    #[test]
    fn test_cutoffs_past_the_minimum_image_bound_are_rejected() {
        use crate::md::*;
        use crate::network::{CutoffTable, TypeParams};
        use ferro_core::{Atom, Cell, Frame};
        use nalgebra::Vector3;
        let frame = |l: f64| {
            let cell = Cell::from_lengths_angles(l, l, l, 90.0, 90.0, 90.0).unwrap();
            let mut f = Frame::with_cell(cell, [true; 3]);
            for (el, x) in [("Si", 0.0), ("O", 1.6), ("H", 2.5), ("O", 3.2)] {
                f.add_atom(Atom::new(el, Vector3::new(x, 0.0, 0.0)));
            }
            f
        };
        let mut traj = Trajectory::new();
        traj.add_frame(frame(10.0));
        traj.add_frame(frame(8.0));

        let net = |r: f64| {
            let mut c = CutoffTable::new();
            c.insert(("Si".into(), "O".into()), r);
            crate::calc_network(&traj, &TypeParams::new(c, CutoffTable::new())).map(|_| ())
        };
        let run = |r: f64| -> Vec<(&str, ferro_core::Result<()>)> {
            vec![
                ("angle", calc_angle(&traj, &AngleParams { r_cut_bc: r, ..Default::default() }).map(|_| ())),
                ("bondlife", calc_bondlife(&traj, &BondLifeParams {
                    r_bond: 1.0, r_break: Some(r), ..Default::default() }).map(|_| ())),
                ("rotcorr", calc_rotcorr(&traj, &RotCorrParams { r_cut: r, ..Default::default() }).map(|_| ())),
                ("network", net(r)),
                ("map sdf", calc_cluster_sdf(&traj, &ClusterSdfParams {
                    former: "Si".into(), ligand: "O".into(), former_ligand_cutoff: r, ..Default::default()
                }).map(|_| ())),
                // Review D-A5: the scan reaches past the selection's own cutoff
                ("cutoff scan", crate::ml::diagnostics::cutoff_scan(&traj, "Si", "O", 4, &[r - 0.6, r], 10)
                    .map(|_| ())),
            ]
        };
        for (who, res) in run(4.5) {
            let msg = res.err().map(|e| e.to_string()).unwrap_or_default();
            assert!(msg.contains("minimum-image") && msg.contains("frame 1"), "{who} 应拦下 4.5 Å：{msg:?}");
        }
        for (who, res) in run(3.9) {
            let msg = res.err().map(|e| e.to_string()).unwrap_or_default();
            assert!(!msg.contains("minimum-image"), "{who} 不应拦下 3.9 Å：{msg}");
        }
    }

    /// Review A-7: a frame that carries a box but is not 3-D periodic (a molecule in
    /// a box, a slab) used to be measured as a bulk. Two O atoms 2.0 A apart in a
    /// 3 A box: gr saw them at 1.0 A through an image that does not exist.
    #[test]
    fn frames_that_are_not_3d_periodic_are_refused() {
        use crate::md::*;
        use crate::ml::filter::{filter_frames, FilterParams};
        use crate::network::{CutoffTable, TypeParams};
        use ferro_core::{Atom, Cell, Frame};
        use nalgebra::Vector3;
        for pbc in [[false; 3], [true, true, false]] {
            let mut traj = Trajectory::new();
            for k in 0..3 {
                let cell = Cell::from_lengths_angles(3.0, 3.0, 3.0, 90.0, 90.0, 90.0).unwrap();
                let mut f = Frame::with_cell(cell, pbc);
                let x = 0.01 * k as f64;
                f.add_atom(Atom::new("O", Vector3::new(0.5 + x, 1.5, 1.5)));
                f.add_atom(Atom::new("O", Vector3::new(2.5, 1.5, 1.5)));
                f.add_atom(Atom::new("P", Vector3::new(1.5, 0.5, 1.5)));
                traj.add_frame(f);
            }
            let mut c = CutoffTable::new();
            c.insert(("P".into(), "O".into()), 1.2);
            let res: Vec<(&str, ferro_core::Result<()>)> = vec![
                ("gr", calc_gr(&traj, &GrParams::default()).map(|_| ())),
                ("angle", calc_angle(&traj, &AngleParams { r_cut_ab: 1.2, r_cut_bc: 1.2, ..Default::default() }).map(|_| ())),
                ("msd", calc_msd(&traj, &MsdParams::default()).map(|_| ())),
                ("rotcorr", calc_rotcorr(&traj, &RotCorrParams {
                    center: "P".into(), neighbor: "O".into(), r_cut: 1.2, ..Default::default() }).map(|_| ())),
                ("bondlife", calc_bondlife(&traj, &BondLifeParams {
                    center: "O".into(), neighbor: "O".into(), r_bond: 1.2, ..Default::default() }).map(|_| ())),
                ("vanhove", calc_vanhove(&traj, &VanHoveParams::default()).map(|_| ())),
                ("map radius", calc_cube_radius(&traj, &CubeRadiusParams::default()).map(|_| ())),
                ("map sdf", calc_cluster_sdf(&traj, &ClusterSdfParams {
                    former: "P".into(), ligand: "O".into(), target_qn: 0, ..Default::default() }).map(|_| ())),
                ("network", crate::calc_network(&traj, &TypeParams::new(c, CutoffTable::new())).map(|_| ())),
                ("filter", filter_frames(&traj, &FilterParams { oo_min: 1.2, ..Default::default() }).map(|_| ())),
                ("map density", calc_cube_density(&traj, &CubeDensityParams::default()).map(|_| ())),
                ("cube jump", calc_cube_jump(&traj, &CubeJumpParams { tau: 1, ..Default::default() }).map(|_| ())),
            ];
            for (who, r) in res {
                let msg = r.err().map(|e| e.to_string()).unwrap_or_default();
                assert!(msg.contains("frame 0") && msg.contains("only 3-D periodic"), "{who} {pbc:?}: {msg:?}");
            }
        }
    }

    /// The five time-correlation analyses check the frame spacing themselves
    /// (review A-5): a restart written twice is refused whoever calls them.
    #[test]
    fn time_correlations_refuse_a_duplicated_step() {
        use crate::md::*;
        use ferro_core::{Atom, Cell, Frame};
        use nalgebra::Vector3;
        let mut traj = Trajectory::new();
        for (k, step) in [0, 10, 10, 20].into_iter().enumerate() {
            let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
            let mut f = Frame::with_cell(cell, [true; 3]);
            f.add_atom(Atom::new("P", Vector3::new(5.0 + 0.01 * k as f64, 5.0, 5.0)));
            f.add_atom(Atom::new("O", Vector3::new(6.5, 5.0, 5.0)));
            f.set_velocities(Some(vec![Vector3::zeros(); 2])).unwrap();
            f.step = Some(step);
            traj.add_frame(f);
        }
        let res: Vec<(&str, ferro_core::Result<()>)> = vec![
            ("msd", calc_msd(&traj, &MsdParams::default()).map(|_| ())),
            ("vacf", calc_vacf(&traj, &VacfParams::default()).map(|_| ())),
            ("rotcorr", calc_rotcorr(&traj, &RotCorrParams {
                center: "P".into(), neighbor: "O".into(), r_cut: 2.0, ..Default::default() }).map(|_| ())),
            ("bondlife", calc_bondlife(&traj, &BondLifeParams {
                center: "P".into(), neighbor: "O".into(), r_bond: 2.0, ..Default::default() }).map(|_| ())),
            ("vanhove", calc_vanhove(&traj, &VanHoveParams::default()).map(|_| ())),
        ];
        for (who, r) in res {
            let msg = r.err().map(|e| e.to_string()).unwrap_or_default();
            assert!(msg.contains("steps 10 and 10"), "{who}: {msg:?}");
        }
    }
}
