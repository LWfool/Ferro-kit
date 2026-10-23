//! File read / write bindings, through the format table of `ferro_io::format`
//! (the same one `ferro convert` uses).

use std::path::Path;

use ferro_io::{AtomStyle, Format, LammpsUnits, ReadOptions, read_aimd, TableFormat};
use ferro_core::Table;
use pyo3::prelude::*;

use crate::pyerr;
use crate::types::PyTrajectory;

/// Raises each reader notice as a `ferro.FerroWarning` / `ferro.FerroNote`, the channel a
/// Python caller filters and captures; the Rust readers never print.
///
/// The text starts with `source[:line]: `, so the same notice from two files is two warnings:
/// Python shows a warning once per (text, category, caller line), and without the path a
/// loop over several files showed only the first (scan E-12). `source`, `line` and `detail`
/// ride on the instance for ferro-py's `_show_warning`, which renders the CLI's layout.
fn emit_warnings(py: Python<'_>, source: &str, notices: &[ferro_io::Notice]) -> PyResult<()> {
    let warn = py.import("warnings")?.getattr("warn")?;
    for n in notices {
        let category = match n.level {
            ferro_io::Level::Note => py.get_type::<crate::FerroNote>(),
            ferro_io::Level::Warning => py.get_type::<crate::FerroWarning>(),
        };
        let place = match n.line {
            Some(l) => format!("{source}:{l}"),
            None => source.to_string(),
        };
        let w = category.call1((format!("{place}: {}", n.message),))?;
        w.setattr("source", source)?;
        w.setattr("line", n.line)?;
        w.setattr("detail", &n.message)?;
        // stacklevel 1 is the Python line that called into ferro: there is no Python frame here
        warn.call1((w, py.None(), 1))?;
    }
    Ok(())
}

fn parse_units(units: Option<&str>) -> PyResult<Option<LammpsUnits>> {
    units.map(|u| u.parse::<LammpsUnits>().map_err(pyerr)).transpose()
}

fn parse_atom_style(atom_style: Option<&str>) -> PyResult<Option<AtomStyle>> {
    atom_style.map(|s| s.parse::<AtomStyle>().map_err(pyerr)).transpose()
}

/// Checks `units` / `atom_style` as `read` does, for `ferro-py` to call before reading the
/// first file: a typo is reported at once, not after a long trajectory has been read.
#[pyfunction]
#[pyo3(name = "_check_read_options", signature = (units = None, atom_style = None))]
fn check_read_options(units: Option<&str>, atom_style: Option<&str>) -> PyResult<()> {
    parse_units(units)?;
    parse_atom_style(atom_style)?;
    Ok(())
}

/// The frame-selection checks of `Trajectory.frame_indices` that need no trajectory.
#[pyfunction]
#[pyo3(name = "_check_selection", signature = (start = None, end = None, stride = None, number = None))]
fn check_selection(
    start: Option<usize>,
    end: Option<usize>,
    stride: Option<usize>,
    number: Option<usize>,
) -> PyResult<()> {
    crate::types::check_selection(start, end, stride, number)
}

/// Reads a structure / trajectory file, recognising its format as `ferro convert` does.
///
/// Recognised by name (any case): xyz, extxyz, pdb, cif, POSCAR/CONTCAR/vasp/pos,
/// in/qe (Quantum ESPRESSO input), inp/restart (CP2K), lammpstrj/dump (LAMMPS dump),
/// data/lmp/lammps (LAMMPS data). A file whose name says nothing (OUTCAR, vasprun.xml,
/// a CP2K / pw.x `.out`) is read by its program banner, dropping unconverged frames;
/// a directory holding `type.raw` is a DeePMD system. What the reader noticed (a
/// dropped frame, a site-label split) arrives as `ferro.FerroWarning` / `ferro.FerroNote`
/// (both `UserWarning`s), the path at the front of the text.
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
        atom_style: parse_atom_style(atom_style)?,
        type_map: type_map.map(|m| m.into_iter().collect()),
    };
    // Reading is the long part; other Python threads run meanwhile
    let (traj, report) = py.detach(|| ferro_io::read(p, &options)).map_err(pyerr)?;
    emit_warnings(py, path, &report.warnings)?;
    Ok(PyTrajectory::new(traj))
}

/// Writes a trajectory in the format its name declares.
///
/// xyz, extxyz, pdb, lammpstrj/dump hold a whole trajectory; cif, POSCAR/CONTCAR/vasp/pos
/// and data/lmp/lammps hold one structure, so several frames go to one file each,
/// numbered by frame (`conf.cif` -> `conf_0000.cif`, ...). The inputs and
/// program outputs are read-only. A trajectory with velocities / forces written to a
/// dump needs `units` ("real" / "metal"). A charge / magnetic-moment column that only
/// some atoms of a frame carry is an error, never filled with 0; `drop_partial=True`
/// leaves it out (LAMMPS data, whose charge column is mandatory, then writes 0 for all).
#[pyfunction]
#[pyo3(signature = (traj, path, units = None, drop_partial = false))]
fn write(py: Python<'_>, traj: &PyTrajectory, path: &str, units: Option<&str>, drop_partial: bool) -> PyResult<()> {
    let opts = ferro_io::WriteOptions { units: parse_units(units)?, drop_partial };
    py.detach(|| ferro_io::write(&traj.inner, Path::new(path), &opts)).map_err(pyerr)?;
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
    let (traj, s) = py.detach(|| read_aimd(Path::new(path))).map_err(pyerr)?;
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

/// Write a table in the format of ferro's other analysis products (`# ` metadata lines, a
/// header, the data; NaN is written as an empty field). `columns` is
/// `[(name, list of numbers or of strings)]`, all of the same length.
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
    m.add_function(wrap_pyfunction!(check_read_options, m)?)?;
    m.add_function(wrap_pyfunction!(check_selection, m)?)?;
    Ok(())
}
