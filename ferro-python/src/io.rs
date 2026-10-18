//! File read / write bindings, through the format table of `ferro_io::format`
//! (the same one `ferro convert` uses).

use std::path::Path;

use ferro_io::{AtomStyle, Format, LammpsUnits, ReadOptions, read_aimd_with_stats, TableFormat};
use ferro_core::Table;
use pyo3::prelude::*;

use crate::pyerr;
use crate::types::PyTrajectory;

/// Raises each reader warning as a Python `UserWarning`, the channel a Python caller
/// filters and captures; the Rust readers never print.
fn emit_warnings(py: Python<'_>, warnings: &[String]) -> PyResult<()> {
    let category = py.get_type::<pyo3::exceptions::PyUserWarning>();
    for w in warnings {
        let msg = std::ffi::CString::new(w.as_str()).map_err(pyerr)?;
        PyErr::warn(py, &category, &msg, 1)?;
    }
    Ok(())
}

fn parse_units(units: Option<&str>) -> PyResult<Option<LammpsUnits>> {
    units.map(|u| u.parse::<LammpsUnits>().map_err(pyerr)).transpose()
}

/// Reads a structure / trajectory file, recognising its format as `ferro convert` does.
///
/// Recognised by name (any case): xyz, extxyz, pdb, cif, POSCAR/CONTCAR/vasp/pos,
/// in/qe (Quantum ESPRESSO input), inp/restart (CP2K), lammpstrj/dump (LAMMPS dump),
/// data/lmp/lammps (LAMMPS data). A file whose name says nothing (OUTCAR, vasprun.xml,
/// a CP2K / pw.x `.out`) is read by its program banner, dropping unconverged frames;
/// a directory holding `type.raw` is a DeePMD system. What the reader noticed (a
/// dropped frame, a site-label split) arrives as `UserWarning`s.
///
/// A LAMMPS dump with velocity / force columns needs `units` ("real" / "metal"): the
/// dump does not record them, and they are never guessed. LAMMPS data needs
/// `atom_style` ("atomic" / "charge" / "full"), never guessed from the comment or the
/// column count. A dump without an element column takes `type_map`
/// (`{1: "Zr", 2: "O"}`), the same function as the CLI's `--type-map`.
#[pyfunction]
#[pyo3(signature = (path, units = None, atom_style = None, type_map = None))]
fn read(
    py: Python<'_>,
    path: &str,
    units: Option<&str>,
    atom_style: Option<&str>,
    type_map: Option<std::collections::BTreeMap<u32, String>>,
) -> PyResult<PyTrajectory> {
    let p = Path::new(path);
    let by_name = Format::from_name(p);
    // Checked before reading, in Python's own words: the parameters are what is wrong
    if type_map.is_some() && by_name != Some(Format::LammpsDump) {
        return Err(pyerr(format!("type_map applies to LAMMPS dump input only, not {path}")));
    }
    if atom_style.is_none() && by_name == Some(Format::LammpsData) {
        return Err(pyerr(format!(
            "{path} is a LAMMPS data file: pass atom_style='atomic', 'charge' or 'full'"
        )));
    }
    let options = ReadOptions {
        units: parse_units(units)?,
        atom_style: atom_style.map(|s| s.parse::<AtomStyle>().map_err(pyerr)).transpose()?,
        type_map: type_map.map(|m| m.into_iter().collect()),
    };
    // Reading is the long part; other Python threads run meanwhile
    let (traj, report) = py.detach(|| ferro_io::read(p, &options)).map_err(pyerr)?;
    emit_warnings(py, &report.warnings)?;
    Ok(PyTrajectory::new(traj))
}

/// Writes a trajectory in the format its name declares.
///
/// xyz, extxyz, pdb, lammpstrj/dump hold a whole trajectory; cif, POSCAR/CONTCAR/vasp/pos
/// and data/lmp/lammps hold one structure, so several frames go to one file each, named
/// as `ferro convert` names them (`conf.cif` -> `conf_0000.cif`, ...). The inputs and
/// program outputs are read-only. A trajectory with velocities / forces written to a
/// dump needs `units` ("real" / "metal").
#[pyfunction]
#[pyo3(signature = (traj, path, units = None))]
fn write(py: Python<'_>, traj: &PyTrajectory, path: &str, units: Option<&str>) -> PyResult<()> {
    let units = parse_units(units)?;
    py.detach(|| ferro_io::write(&traj.inner, Path::new(path), units)).map_err(pyerr)?;
    Ok(())
}

/// Reads a calculation output (VASP OUTCAR / vasprun.xml, CP2K single-point and MD
/// output, pw.x output), identifying the program by content, as `ferro dataset collect`
/// does. Returns `(Trajectory, stats)`: frames carry energy, forces and stress; frames
/// dropped for an unconverged SCF and the like are not in the trajectory, their counts
/// are in the stats (`scf_failed`, `incomplete`, ...). `version_note` says why the
/// program version deserves a second look (None when it does not), and `step_range` is the
/// `(first, last)` step number the file claims.
#[pyfunction]
fn read_output(
    py: Python<'_>,
    path: &str,
) -> PyResult<(PyTrajectory, std::collections::BTreeMap<String, Py<PyAny>>)> {
    let (traj, s) = py.detach(|| read_aimd_with_stats(Path::new(path))).map_err(pyerr)?;
    {
        let mut stats = std::collections::BTreeMap::new();
        let mut put = |k: &str, v: Py<PyAny>| { stats.insert(k.to_string(), v); };
        put("format", s.format.name().into_pyobject(py)?.into_any().unbind());
        for (k, v) in [("steps", s.n_steps), ("kept", s.n_kept), ("scf_failed", s.n_scf_failed),
                       ("incomplete", s.n_incomplete), ("bad_composition", s.n_bad_composition),
                       ("no_cell", s.n_no_cell), ("layout_drift", s.n_layout_drift),
                       ("restarts", s.n_restarts)] {
            put(k, v.into_pyobject(py)?.into_any().unbind());
        }
        put("version", s.version.into_pyobject(py)?.into_any().unbind());
        put("version_note", s.version_note.into_pyobject(py)?.into_any().unbind());
        put("step_range", s.steps.into_pyobject(py)?.into_any().unbind());

        Ok((PyTrajectory::new(traj), stats))
    }
}

/// 写一张表，格式与 ferro 其他分析产物相同（`# ` 开头的元数据行、表头、数据；NaN 写成空格）。
/// `columns` 为 `[(列名, 数值列表或字符串列表)]`，各列等长。
#[pyfunction]
#[pyo3(signature = (path, columns, meta = Vec::new()))]
fn write_table(path: &str, columns: Vec<(String, Bound<'_, PyAny>)>, meta: Vec<String>) -> PyResult<()> {
    let mut table = Table::new();
    for line in meta {
        table.meta_line(line);
    }
    for (name, values) in columns {
        if let Ok(v) = values.extract::<Vec<f64>>() {
            table.push_num(name, v);
        } else if let Ok(v) = values.extract::<Vec<String>>() {
            table.push_text(name, v);
        } else {
            return Err(pyerr(format!("column {name}: values must be all numbers or all strings")));
        }
    }
    ferro_io::write_table(&table, Path::new(path), TableFormat::Csv).map_err(pyerr)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read, m)?)?;
    m.add_function(wrap_pyfunction!(write, m)?)?;
    m.add_function(wrap_pyfunction!(read_output, m)?)?;
    m.add_function(wrap_pyfunction!(write_table, m)?)?;
    Ok(())
}
