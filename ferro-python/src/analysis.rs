//! 分析绑定（g(r) / MSD）。返回 `dict[str, list[float]]`。

use std::collections::HashMap;

use ferro_analysis::{calc_gr, calc_msd, GrParams, GroupBy, MsdParams};
use pyo3::prelude::*;

use crate::pyerr;
use crate::types::PyTrajectory;

/// 单个配对的径向分布函数与有向配位数。
///
/// `a` 为中心原子，`b` 为近邻原子，顺序有意义：`cn` 是每个 `a` 周围的 `b` 数。
/// `by="element"`（默认）按元素分组，`by="label"` 按位点标签分组。
///
/// 返回字典：`"r"` → bin 中心；`"<a>-<b>_gr"`；`"<a>-<b>_cn"`。
/// `r_max=None` 时按首帧的最小镜像上界（最小面间距的一半）自动确定。
#[pyfunction]
#[pyo3(signature = (traj, a, b, by="element", r_max=None, dr=0.01, r_min=0.005))]
fn gr_pair(
    traj: &PyTrajectory,
    a: &str,
    b: &str,
    by: &str,
    r_max: Option<f64>,
    dr: f64,
    r_min: f64,
) -> PyResult<HashMap<String, Vec<f64>>> {
    let res = run_gr(traj, by, r_max, dr, r_min)?;
    let key = format!("{a}-{b}");
    let g = res.gr.get(&key).ok_or_else(|| {
        pyo3::exceptions::PyKeyError::new_err(format!(
            "pair '{key}' not present; available types: {:?}", res.elements
        ))
    })?;
    let cn = &res.cn[&key];

    let mut out: HashMap<String, Vec<f64>> = HashMap::new();
    out.insert("r".to_string(), res.r.clone());
    out.insert(format!("{key}_gr"), g.clone());
    out.insert(format!("{key}_cn"), cn.clone());
    Ok(out)
}

/// 全部有序配对的 g(r) 与 CN(r)（一次遍历算出所有对）。
///
/// 返回字典：`"r"`，以及每个有序对的 `"<A>-<B>_gr"` / `"<A>-<B>_cn"`。
/// n 个类型给出 n² 个配对；`gr` 对称（`A-B` 与 `B-A` 相同），`cn` 有向。
#[pyfunction]
#[pyo3(signature = (traj, by="element", r_max=None, dr=0.01, r_min=0.005))]
fn gr_all(
    traj: &PyTrajectory,
    by: &str,
    r_max: Option<f64>,
    dr: f64,
    r_min: f64,
) -> PyResult<HashMap<String, Vec<f64>>> {
    let res = run_gr(traj, by, r_max, dr, r_min)?;

    let mut out: HashMap<String, Vec<f64>> = HashMap::new();
    out.insert("r".to_string(), res.r);
    for (k, v) in res.gr {
        out.insert(format!("{k}_gr"), v);
    }
    for (k, v) in res.cn {
        out.insert(format!("{k}_cn"), v);
    }
    Ok(out)
}

/// `gr_pair` / `gr_all` 共用的计算入口。
fn run_gr(
    traj: &PyTrajectory,
    by: &str,
    r_max: Option<f64>,
    dr: f64,
    r_min: f64,
) -> PyResult<ferro_analysis::GrResult> {
    let group_by = match by {
        "element" => GroupBy::Element,
        "label" => GroupBy::Label,
        other => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "by must be 'element' or 'label', got '{other}'"
            )))
        }
    };
    let r_max = r_max.unwrap_or_else(|| GrParams::with_auto_rmax(&traj.inner).r_max);
    let params = GrParams { r_min, r_max, dr, group_by };
    calc_gr(&traj.inner, &params).map_err(pyerr)
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
    let params = MsdParams { tau, shift, dt, elements, fit_range: None };
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
    m.add_function(wrap_pyfunction!(gr_pair, m)?)?;
    m.add_function(wrap_pyfunction!(gr_all, m)?)?;
    m.add_function(wrap_pyfunction!(msd, m)?)?;
    Ok(())
}
