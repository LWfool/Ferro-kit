//! 文件读写绑定（按扩展名分派）。

use std::path::Path;

use ferro_io::{
    read_cif, read_contcar, read_cp2k_inp, read_cp2k_restart, read_extxyz, read_lammps_data,
    read_lammps_dump, read_pdb, read_poscar, read_qe_input, read_xyz, write_cif, write_extxyz,
    write_lammps_data, write_lammps_dump, write_pdb, write_poscar, write_qe_input, write_xyz,
    LammpsUnits,
};
use pyo3::prelude::*;

use crate::pyerr;
use crate::types::PyTrajectory;

/// 按文件名/扩展名归一化的格式标识。
fn detect(path: &str) -> String {
    let p = Path::new(path);
    let stem = p
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_uppercase();
    if stem.starts_with("POSCAR") {
        return "poscar".into();
    }
    if stem.starts_with("CONTCAR") {
        return "contcar".into();
    }
    p.extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_lowercase()
}

/// 读取结构 / 轨迹文件，按扩展名自动识别格式。
///
/// 支持：xyz, extxyz, pdb, cif, POSCAR/CONTCAR, in/qe (Quantum ESPRESSO),
/// inp/restart (CP2K), lammpstrj/dump/lammps (LAMMPS dump),
/// data/lmp (LAMMPS data)。LAMMPS dump 默认 real 单位，`metal_units=True`
/// 切换为 metal 单位。
#[pyfunction]
#[pyo3(signature = (path, metal_units = false))]
fn read(path: &str, metal_units: bool) -> PyResult<PyTrajectory> {
    let units = if metal_units {
        LammpsUnits::Metal
    } else {
        LammpsUnits::Real
    };
    let p = std::path::Path::new(path);
    let traj = match detect(path).as_str() {
        "xyz" => read_xyz(p).map_err(pyerr)?,
        "extxyz" => read_extxyz(p).map_err(pyerr)?,
        "pdb" => read_pdb(p).map_err(pyerr)?,
        "cif" => read_cif(p).map_err(pyerr)?,
        "poscar" => read_poscar(p).map_err(pyerr)?,
        "contcar" => read_contcar(p).map_err(pyerr)?,
        "in" | "qe" => read_qe_input(p).map_err(pyerr)?,
        "inp" => read_cp2k_inp(p).map_err(pyerr)?,
        "restart" => read_cp2k_restart(p).map_err(pyerr)?,
        "lammpstrj" | "dump" | "lammps" => read_lammps_dump(p, units).map_err(pyerr)?,
        "data" | "lmp" => read_lammps_data(p).map_err(pyerr)?,
        other => {
            return Err(pyerr(format!(
                "unsupported input format '{other}' for path '{path}'"
            )))
        }
    };
    Ok(PyTrajectory::new(traj))
}

/// 写出轨迹，按扩展名选择格式。
///
/// 支持：xyz, extxyz, pdb, cif, POSCAR, in/qe, data/lmp,
/// lammpstrj/dump。
#[pyfunction]
#[pyo3(signature = (traj, path, metal_units = false))]
fn write(traj: &PyTrajectory, path: &str, metal_units: bool) -> PyResult<()> {
    let units = if metal_units {
        LammpsUnits::Metal
    } else {
        LammpsUnits::Real
    };
    let t = &traj.inner;
    let p = std::path::Path::new(path);
    match detect(path).as_str() {
        "xyz" => write_xyz(t, p).map_err(pyerr),
        "extxyz" => write_extxyz(t, p).map_err(pyerr),
        "pdb" => write_pdb(t, p).map_err(pyerr),
        "cif" => write_cif(t, p).map_err(pyerr),
        "poscar" => write_poscar(t, p).map_err(pyerr),
        "in" | "qe" => write_qe_input(t, p).map_err(pyerr),
        "data" | "lmp" => write_lammps_data(t, p).map_err(pyerr),
        "lammpstrj" | "dump" => write_lammps_dump(t, p, units).map_err(pyerr),
        other => Err(pyerr(format!(
            "unsupported output format '{other}' for path '{path}'"
        ))),
    }
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read, m)?)?;
    m.add_function(wrap_pyfunction!(write, m)?)?;
    Ok(())
}
