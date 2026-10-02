//! 参数取值范围的检查，供各 `XxxParams::validate` 共用（审查 M3）。
//!
//! 一律写成「合法则放行，否则报错」：NaN 与任何数比较都为假，落到报错分支，
//! `--dr nan` 也会被拦下；反过来写 `if v <= 0.0 { 报错 }` 会把 NaN 放过去。名字用 CLI 的参数拼法（`r-max`），
//! 与 `correlate::resolve_max_lag` 的 `max-lag` 一致。

use ferro_core::error::ChemError;

type Check = ferro_core::Result<()>;

fn fail(msg: String) -> Check {
    Err(ChemError::ValidationError(msg))
}

/// 步长、截断、时间步：有限且 > 0（`inf` 的步长会得到 0 个 bin）
pub(crate) fn positive(name: &str, v: f64) -> Check {
    if v > 0.0 && v.is_finite() { return Ok(()); }
    fail(format!("{name} must be a finite number > 0, got {v}"))
}

pub(crate) fn non_negative(name: &str, v: f64) -> Check {
    if v >= 0.0 { return Ok(()); }
    fail(format!("{name} must be >= 0, got {v}"))
}

/// `lo < hi`
pub(crate) fn ordered(lo_name: &str, lo: f64, hi_name: &str, hi: f64) -> Check {
    if lo < hi { return Ok(()); }
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
        case!(SqParams, q_min: -1.0, q_max: 0.0, dq: 0.0, dq: nan);
        case!(MsdParams, dt: 0.0, dt: nan, max_lag: Some(0),
              fit_range: Some((0.8, 0.2)), fit_range: Some((0.0, 1.5)), fit_range: Some((nan, 0.5)));
        case!(AngleParams, r_cut_ab: 0.0, r_cut_bc: nan, angle_min: -1.0, angle_max: 0.0,
              angle_max: 181.0, d_angle: 0.0, d_angle: f64::INFINITY);
        case!(VacfParams, dt: -1.0, max_lag: Some(0));
        case!(RotCorrParams, r_cut: 0.0, dt: nan, max_lag: Some(0));
        case!(BondLifeParams, r_bond: 0.0, r_bond: nan, r_break: Some(1.0), dt: 0.0, max_lag: Some(0));
        case!(VanHoveParams, dt: 0.0, tau: Some(0), shift: 0, r_min: -1.0, dr: 0.0, r_max: 0.0);
        case!(CubeDensityParams, nx: 0, ny: 0, nz: 0);
        case!(CubeRadiusParams, nx: 0, ny: 0, nz: 0, radius: 0.0, radius: nan);
        case!(ClusterSdfParams, former_ligand_cutoff: 0.0, modifier_cutoff: nan, grid_res: 0.0,
              sigma: -1.0, padding: -1.0, rmsd_warn_threshold: nan);
        case!(ChgSdfParams, former_ligand_cutoff: 0.0, modifier_cutoff: -1.0, padding: nan,
              rmsd_warn_threshold: -0.5);
    }
}
