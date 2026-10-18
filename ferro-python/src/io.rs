//! 文件读写绑定（按扩展名分派）。

use std::path::Path;

use ferro_io::{
    apply_type_map, is_extxyz, read_cif, read_contcar, read_cp2k_inp, read_cp2k_restart, read_extxyz, read_lammps_data, AtomStyle,
    read_lammps_dump, read_pdb, read_poscar, read_qe_input, read_xyz, write_cif, write_extxyz,
    write_lammps_data, write_lammps_dump, write_pdb, write_poscar, write_xyz,
    index_width, indexed_path,
    LammpsUnits, read_aimd_with_stats, TableFormat,
};
use ferro_core::Table;
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

fn parse_units(units: Option<&str>) -> PyResult<Option<LammpsUnits>> {
    units.map(|u| u.parse::<LammpsUnits>().map_err(pyerr)).transpose()
}

/// 读取结构 / 轨迹文件，按扩展名自动识别格式。
///
/// 支持：xyz, extxyz, pdb, cif, POSCAR/CONTCAR/vasp/pos, in/qe (Quantum ESPRESSO),
/// inp/restart (CP2K), lammpstrj/dump (LAMMPS dump),
/// data/lmp/lammps (LAMMPS data)。LAMMPS dump 含速度/力列时必须给 `units`（"real" /
/// "metal"），dump 不记录单位，不猜。LAMMPS data 必须给 `atom_style`（"atomic" / "charge" /
/// "full"），不按注释或列数猜。无 element 列的 dump 用 `type_map`（`{1: "Zr", 2: "O"}`）
/// 把类型号换成元素，与 CLI 的 `--type-map` 同一个函数。
#[pyfunction]
#[pyo3(signature = (path, units = None, atom_style = None, type_map = None))]
fn read(
    path: &str,
    units: Option<&str>,
    atom_style: Option<&str>,
    type_map: Option<std::collections::BTreeMap<u32, String>>,
) -> PyResult<PyTrajectory> {
    let units = parse_units(units)?;
    let p = std::path::Path::new(path);
    let format = detect(path);
    let is_dump = matches!(format.as_str(), "lammpstrj" | "dump");
    if type_map.is_some() && !is_dump {
        return Err(pyerr(format!("type_map applies to LAMMPS dump input only, not {path}")));
    }
    let mut traj = match format.as_str() {
        // .xyz 也常装 extxyz（ASE、GPUMD），与 CLI 的 io_dispatch 同一个判据
        "xyz" if is_extxyz(p).map_err(pyerr)? => read_extxyz(p).map_err(pyerr)?,
        "xyz" => read_xyz(p).map_err(pyerr)?,
        "extxyz" => read_extxyz(p).map_err(pyerr)?,
        "pdb" => read_pdb(p).map_err(pyerr)?,
        "cif" => read_cif(p).map_err(pyerr)?,
        "poscar" | "vasp" | "pos" => read_poscar(p).map_err(pyerr)?,
        "contcar" => read_contcar(p).map_err(pyerr)?,
        "in" | "qe" => read_qe_input(p).map_err(pyerr)?,
        "inp" => read_cp2k_inp(p).map_err(pyerr)?,
        "restart" => read_cp2k_restart(p).map_err(pyerr)?,
        "lammpstrj" | "dump" => read_lammps_dump(p, units).map_err(pyerr)?,
        "data" | "lmp" | "lammps" => {
            let style = atom_style.ok_or_else(|| pyerr(format!(
                "{path} is a LAMMPS data file: pass atom_style='atomic', 'charge' or 'full'"
            )))?;
            read_lammps_data(p, style.parse::<AtomStyle>().map_err(pyerr)?).map_err(pyerr)?
        }
        other => {
            return Err(pyerr(format!(
                "unsupported input format '{other}' for path '{path}'"
            )))
        }
    };
    if let Some(map) = type_map {
        let map: Vec<(u32, String)> = map.into_iter().collect();
        apply_type_map(&mut traj, &map).map_err(pyerr)?;
    }
    Ok(PyTrajectory::new(traj))
}

/// 写出轨迹，按扩展名选择格式。
///
/// 支持：xyz, extxyz, pdb, cif, POSCAR/CONTCAR/vasp/pos, data/lmp/lammps, lammpstrj/dump（QE 输入只读）。
/// 带速度/力的轨迹写 dump 时必须给 `units`（"real" / "metal"）。
#[pyfunction]
#[pyo3(signature = (traj, path, units = None))]
fn write(traj: &PyTrajectory, path: &str, units: Option<&str>) -> PyResult<()> {
    let units = parse_units(units)?;
    let t = &traj.inner;
    let p = std::path::Path::new(path);
    // POSCAR, CIF and LAMMPS data hold one structure: several frames go to one file
    // each, named as `ferro convert` names them (`conf.cif` -> `conf_0000.cif`, ...)
    let single: fn(&ferro_core::Trajectory, &Path) -> PyResult<()> = match detect(path).as_str() {
        "xyz" => return write_xyz(t, p).map_err(pyerr),
        "extxyz" => return write_extxyz(t, p).map_err(pyerr),
        "pdb" => return write_pdb(t, p).map_err(pyerr),
        "lammpstrj" | "dump" => return write_lammps_dump(t, p, units).map_err(pyerr),
        "cif" => |t, p| write_cif(t, p).map_err(pyerr),
        "poscar" | "contcar" | "vasp" | "pos" => |t, p| write_poscar(t, p).map_err(pyerr),
        "data" | "lmp" | "lammps" => |t, p| write_lammps_data(t, p).map_err(pyerr),
        other => return Err(pyerr(format!(
            "unsupported output format '{other}' for path '{path}'"
        ))),
    };
    if t.n_frames() <= 1 {
        return single(t, p);
    }
    let width = index_width(t.n_frames() - 1);
    for (i, frame) in t.frames.iter().enumerate() {
        let one = ferro_core::Trajectory { frames: vec![frame.clone()], metadata: t.metadata.clone() };
        single(&one, &indexed_path(p, i, width))?;
    }
    Ok(())
}

/// 读计算输出（VASP OUTCAR / vasprun.xml、CP2K 单点与 MD 输出），按文件内容识别程序，
/// 同 `ferro dataset collect`。返回 `(Trajectory, 统计)`：帧带能量、力、应力；SCF 未收敛等
/// 被丢弃的帧不在轨迹里，个数在统计里（`scf_failed`、`incomplete`…），调用方据此报原因。
#[pyfunction]
fn read_output(path: &str) -> PyResult<(PyTrajectory, std::collections::BTreeMap<String, Py<PyAny>>)> {
    let (traj, s) = read_aimd_with_stats(Path::new(path)).map_err(pyerr)?;
    Python::attach(|py| {
        let mut stats = std::collections::BTreeMap::new();
        let mut put = |k: &str, v: Py<PyAny>| { stats.insert(k.to_string(), v); };
        put("format", s.format.name().into_pyobject(py)?.into_any().unbind());
        for (k, v) in [("steps", s.n_steps), ("kept", s.n_kept), ("scf_failed", s.n_scf_failed),
                       ("incomplete", s.n_incomplete), ("bad_composition", s.n_bad_composition),
                       ("no_cell", s.n_no_cell)] {
            put(k, v.into_pyobject(py)?.into_any().unbind());
        }
        put("version", s.version.into_pyobject(py)?.into_any().unbind());
        Ok((PyTrajectory::new(traj), stats))
    })
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
