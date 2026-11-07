//! Structure operation bindings (supercell / vacuum layer / merge / mixed box).

use ferro_core::{matrix3_from_row_major, Cell, Trajectory};
use ferro_structure::{
    add_vacuum, build_box as build_mixed_box, make_supercell, merge_frames, Component,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::pyerr;
use crate::types::PyTrajectory;

/// The one frame of `t`. A structure operation builds a starting structure, so a trajectory is
/// refused rather than worked frame by frame, and which frame to take is the caller's to say,
/// as `ferro build` refuses it (user, 2026-11-07; merge refused it before, review 2026-10-30 A6-M-5).
fn one<'a>(t: &'a PyTrajectory, which: &str, op: &str) -> PyResult<&'a ferro_core::Frame> {
    match t.inner.frames.as_slice() {
        [f] => Ok(f),
        frames => Err(PyValueError::new_err(format!(
            "{which} holds {} frames; {op} works on one structure — take one first with .subset([i])",
            frames.len()
        ))),
    }
}

/// The `nx` × `ny` × `nz` supercell of a one-frame trajectory, as a one-frame trajectory.
#[pyfunction]
fn supercell(traj: &PyTrajectory, nx: usize, ny: usize, nz: usize) -> PyResult<PyTrajectory> {
    let made = make_supercell(one(traj, "the trajectory", "supercell")?, nx, ny, nz).map_err(pyerr)?;
    Ok(PyTrajectory::new(Trajectory { frames: vec![made], metadata: traj.inner.metadata.clone() }))
}

/// A `thickness` Å vacuum layer along `axis` ("x"/"y"/"z") of a one-frame trajectory, as a
/// one-frame trajectory.
#[pyfunction]
fn add_vacuum_layer(
    traj: &PyTrajectory,
    axis: &str,
    thickness: f64,
) -> PyResult<PyTrajectory> {
    let made = add_vacuum(one(traj, "the trajectory", "add_vacuum_layer")?, axis, thickness).map_err(pyerr)?;
    Ok(PyTrajectory::new(Trajectory { frames: vec![made], metadata: traj.inner.metadata.clone() }))
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
    let (fa, fb) = (one(a, "the first trajectory", "merge")?, one(b, "the second trajectory", "merge")?);
    let merged = merge_frames(fa, fb, axis, gap).map_err(pyerr)?;
    Ok(PyTrajectory::new(Trajectory::from_frame(merged)))
}

// The seed is written out so `help(ferro.build_box)` shows it (a path shows as `seed=...`);
// this fails to compile if it drifts from the CLI's (review 2026-11-06 D3)
const _: () = assert!(ferro_analysis::ml::dataset::DEFAULT_SEED == 666);

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
