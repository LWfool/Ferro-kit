//! Structure operation bindings (supercell / vacuum layer / merge / mixed box).

use ferro_core::{matrix3_from_row_major, Cell, Trajectory};
use ferro_structure::{
    add_vacuum, build_box as build_mixed_box, check_supercell_total, make_supercell, merge_frames, Component,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::pyerr;
use crate::types::PyTrajectory;

/// 对轨迹每一帧构建 nx×ny×nz 超胞，返回新轨迹。
#[pyfunction]
fn supercell(traj: &PyTrajectory, nx: usize, ny: usize, nz: usize) -> PyResult<PyTrajectory> {
    // Every frame is held at once: cap the whole trajectory before the copy (review fb7d8c3-1)
    check_supercell_total(traj.inner.frames.iter().map(|f| f.n_atoms()), [nx, ny, nz]).map_err(pyerr)?;
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
///
/// 两块都不拉伸、只整体平移；`gap` 沿界面法向量（同 `add_vacuum`），细节见 `merge_frames`。
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

/// Pack `components`, `[(name or formula, count), ...]`, at random into a cubic box of
/// `density` g/cm³; a one-frame trajectory. Defaults are those of `ferro build box`.
#[pyfunction]
#[pyo3(signature = (components, density, min_dist = 1.5, relax_steps = 500, seed = 666))]
fn build_box(
    components: Vec<(String, usize)>,
    density: f64,
    min_dist: f64,
    relax_steps: usize,
    seed: u64,
) -> PyResult<PyTrajectory> {
    let comps: Vec<Component> =
        components.into_iter().map(|(compound, n_molecules)| Component { compound, n_molecules }).collect();
    let frame = build_mixed_box(&comps, density, min_dist, relax_steps, seed).map_err(pyerr)?;
    Ok(PyTrajectory::new(Trajectory::from_frame(frame)))
}

/// Interplanar spacings (d_a, d_b, d_c) [Å] of a cell given as rows a, b, c [Å]: the distance
/// between the lattice planes normal to each reciprocal vector, `1 / |b_i|` (2π not included).
/// The k-point mesh of `ferro.jobs` is built on it; a singular cell is an error.
#[pyfunction]
fn interplanar_spacings(cell: [[f64; 3]; 3]) -> PyResult<(f64, f64, f64)> {
    let [a, b, c] = cell;
    let m = matrix3_from_row_major(&[a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]]);
    let [da, db, dc] = Cell::from_matrix(m).interplanar_spacings().map_err(pyerr)?;
    Ok((da, db, dc))
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(interplanar_spacings, m)?)?;
    m.add_function(wrap_pyfunction!(supercell, m)?)?;
    m.add_function(wrap_pyfunction!(add_vacuum_layer, m)?)?;
    m.add_function(wrap_pyfunction!(merge, m)?)?;
    m.add_function(wrap_pyfunction!(build_box, m)?)?;
    Ok(())
}
