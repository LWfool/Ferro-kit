//! 分析绑定（g(r) / MSD）。返回 `dict[str, list[float]]`。

use std::collections::HashMap;

use ferro_analysis::{calc_gr, calc_msd, GrParams, MsdParams};
use pyo3::prelude::*;

use crate::pyerr;
use crate::types::PyTrajectory;

/// 径向分布函数 g(r) 与配位数 CN(r)。
///
/// 返回字典：`"r"` → bin 中心；`"gr:<El1-El2>"` / `"gr:total"` → 偏 g(r)；
/// `"cn:<center-neighbor>"` → 有向累积配位数。`r_max=None` 时按首帧最短
/// 晶胞向量的一半自动确定。
#[pyfunction]
#[pyo3(signature = (traj, r_max=None, dr=0.01, r_cut=2.3, r_min=0.005))]
fn gr(
    traj: &PyTrajectory,
    r_max: Option<f64>,
    dr: f64,
    r_cut: f64,
    r_min: f64,
) -> PyResult<HashMap<String, Vec<f64>>> {
    let r_max = r_max.unwrap_or_else(|| GrParams::with_auto_rmax(&traj.inner).r_max);
    let params = GrParams { r_min, r_max, dr, r_cut };
    let res = calc_gr(&traj.inner, &params).map_err(pyerr)?;

    let mut out: HashMap<String, Vec<f64>> = HashMap::new();
    out.insert("r".to_string(), res.r);
    for (k, v) in res.gr {
        out.insert(format!("gr:{k}"), v);
    }
    for (k, v) in res.cn {
        out.insert(format!("cn:{k}"), v);
    }
    Ok(out)
}

/// 均方位移 MSD(t)（时间原点平均，NPT 安全）。
///
/// 返回字典：`"time"` \[fs\]、`"msd"`（总）、`"msd_a"/"msd_b"/"msd_c"`
/// （周期体系沿晶轴，非周期体系沿 x/y/z）。
#[pyfunction]
#[pyo3(signature = (traj, dt=1.0, shift=1, tau=None, elements=None))]
fn msd(
    traj: &PyTrajectory,
    dt: f64,
    shift: usize,
    tau: Option<usize>,
    elements: Option<Vec<String>>,
) -> PyResult<HashMap<String, Vec<f64>>> {
    let params = MsdParams { tau, shift, dt, elements };
    let res = calc_msd(&traj.inner, &params).map_err(pyerr)?;

    let mut out: HashMap<String, Vec<f64>> = HashMap::new();
    out.insert("time".to_string(), res.time);
    out.insert("msd".to_string(), res.msd);
    out.insert("msd_a".to_string(), res.msd_a);
    out.insert("msd_b".to_string(), res.msd_b);
    out.insert("msd_c".to_string(), res.msd_c);
    Ok(out)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(gr, m)?)?;
    m.add_function(wrap_pyfunction!(msd, m)?)?;
    Ok(())
}
