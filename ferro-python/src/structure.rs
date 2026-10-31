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
    // Each frame is built anew, so the input frames are not deep-copied first (review R8-10)
    let frames = traj.inner.frames.iter()
        .map(|f| make_supercell(f, nx, ny, nz))
        .collect::<ferro_core::Result<Vec<_>>>().map_err(pyerr)?;
    Ok(PyTrajectory::new(Trajectory { frames, metadata: traj.inner.metadata.clone() }))
}

/// 沿 `axis`（"x"/"y"/"z"）为每一帧添加 `thickness` Å 真空层。
#[pyfunction]
fn add_vacuum_layer(
    traj: &PyTrajectory,
    axis: &str,
    thickness: f64,
) -> PyResult<PyTrajectory> {
    // Each frame is built anew, so the input frames are not deep-copied first (review R8-10)
    let frames = traj.inner.frames.iter()
        .map(|f| add_vacuum(f, axis, thickness))
        .collect::<ferro_core::Result<Vec<_>>>().map_err(pyerr)?;
    Ok(PyTrajectory::new(Trajectory { frames, metadata: traj.inner.metadata.clone() }))
}

/// Stacks two one-frame trajectories along `axis` with a `gap` Å, giving a one-frame trajectory.
///
/// Neither block is strained, only translated; `gap` is along the interface normal (as in
/// `add_vacuum`), see `merge_frames`. A trajectory of several frames is refused, as `ferro build
/// merge` refuses it: which frame to take is the caller's to say (review 2026-10-30 A6-M-5).
#[pyfunction]
fn merge(
    a: &PyTrajectory,
    b: &PyTrajectory,
    axis: &str,
    gap: f64,
) -> PyResult<PyTrajectory> {
    fn one<'a>(t: &'a PyTrajectory, which: &str) -> PyResult<&'a ferro_core::Frame> {
        match t.inner.frames.as_slice() {
            [f] => Ok(f),
            frames => Err(PyValueError::new_err(format!(
                "{which} trajectory holds {} frames; merge stacks two structures — take one first \
                 with .subset([i])",
                frames.len()
            ))),
        }
    }
    let (fa, fb) = (one(a, "first")?, one(b, "second")?);
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
