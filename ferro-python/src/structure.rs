//! 结构操作绑定（超胞 / 真空层 / 合并）。

use ferro_core::Trajectory;
use ferro_structure::{add_vacuum, make_supercell, merge_frames};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::pyerr;
use crate::types::PyTrajectory;

/// 对轨迹每一帧构建 nx×ny×nz 超胞，返回新轨迹。
#[pyfunction]
fn supercell(traj: &PyTrajectory, nx: usize, ny: usize, nz: usize) -> PyResult<PyTrajectory> {
    let mut out = traj.inner.clone();
    for f in out.frames.iter_mut() {
        *f = make_supercell(f, nx, ny, nz).map_err(pyerr)?;
    }
    Ok(PyTrajectory::new(out))
}

/// 沿 `axis`（"x"/"y"/"z"）为每一帧添加 `thickness` Å 真空层。
#[pyfunction]
fn add_vacuum_layer(
    traj: &PyTrajectory,
    axis: &str,
    thickness: f64,
) -> PyResult<PyTrajectory> {
    let mut out = traj.inner.clone();
    for f in out.frames.iter_mut() {
        *f = add_vacuum(f, axis, thickness).map_err(pyerr)?;
    }
    Ok(PyTrajectory::new(out))
}

/// 沿 `axis` 以 `gap` Å 间隙合并两个体系的首帧，返回单帧轨迹。
#[pyfunction]
fn merge(
    a: &PyTrajectory,
    b: &PyTrajectory,
    axis: &str,
    gap: f64,
) -> PyResult<PyTrajectory> {
    let fa = a
        .inner
        .first()
        .ok_or_else(|| PyValueError::new_err("first trajectory is empty"))?;
    let fb = b
        .inner
        .first()
        .ok_or_else(|| PyValueError::new_err("second trajectory is empty"))?;
    let merged = merge_frames(fa, fb, axis, gap).map_err(pyerr)?;
    Ok(PyTrajectory::new(Trajectory::from_frame(merged)))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(supercell, m)?)?;
    m.add_function(wrap_pyfunction!(add_vacuum_layer, m)?)?;
    m.add_function(wrap_pyfunction!(merge, m)?)?;
    Ok(())
}
