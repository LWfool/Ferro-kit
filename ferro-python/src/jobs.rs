//! `ferro.jobs`（纯 Python 的输入生成）要用的 Rust 侧接口：结构块字符串、
//! QE namelist 的拆分与写回、原子排序。格式规范只在 ferro-io 一份，这里只做转换。

use ferro_analysis::ml::merge;
use ferro_io::{Namelist, QeInputParts};
use pyo3::exceptions::PyIndexError;
use pyo3::prelude::*;

use crate::pyerr;
use crate::types::PyTrajectory;

fn frame(traj: &PyTrajectory, index: usize) -> PyResult<&ferro_core::Frame> {
    traj.inner
        .frame(index)
        .ok_or_else(|| PyIndexError::new_err(format!("frame index {index} out of range")))
}

/// 逐帧按 (原子序数, 符号) 重排原子（同元素成块、块按周期表序，同元素内保持原相对顺序），
/// 与 `ferro convert` 对所有产物的排法相同。逐帧各取各的置换。
#[pyfunction]
fn sort_atoms(traj: &PyTrajectory) -> PyTrajectory {
    let mut out = traj.inner.clone();
    out.frames = traj
        .inner
        .frames
        .iter()
        .map(|f| merge::sort_atoms(&ferro_core::Trajectory::from_frame(f.clone())).frames.remove(0))
        .collect();
    PyTrajectory::new(out)
}

/// 指定帧的 POSCAR 全文（`comment` 为首行）。
#[pyfunction]
#[pyo3(signature = (traj, frame_index = 0, comment = "ferro"))]
fn poscar_string(traj: &PyTrajectory, frame_index: usize, comment: &str) -> PyResult<String> {
    ferro_io::poscar_string(frame(traj, frame_index)?, comment).map_err(pyerr)
}

/// 指定帧的 QE 三张结构卡片；`pseudos` 为 元素 → 赝势文件名，须覆盖全部元素。
#[pyfunction]
#[pyo3(signature = (traj, pseudos, frame_index = 0))]
fn qe_structure_cards(
    traj: &PyTrajectory,
    pseudos: std::collections::BTreeMap<String, String>,
    frame_index: usize,
) -> PyResult<String> {
    let pseudos: Vec<(String, String)> = pseudos.into_iter().collect();
    ferro_io::qe_structure_cards(frame(traj, frame_index)?, &pseudos).map_err(pyerr)
}

/// 指定帧的 CP2K `&CELL` 与 `&COORD` 两段。
#[pyfunction]
#[pyo3(signature = (traj, frame_index = 0))]
fn cp2k_structure_sections(traj: &PyTrajectory, frame_index: usize) -> PyResult<String> {
    ferro_io::cp2k_structure_sections(frame(traj, frame_index)?).map_err(pyerr)
}

/// Python 侧的一个 namelist：`(名, [(键, 值原文), …])`
type PyNamelist = (String, Vec<(String, String)>);

/// pw.x 输入 → `([(namelist 名, [(键, 值原文), …]), …], 其余卡片原文)`。
#[pyfunction]
fn parse_qe_input_parts(text: &str) -> PyResult<(Vec<PyNamelist>, String)> {
    let QeInputParts { namelists, cards } = ferro_io::parse_qe_input_parts(text).map_err(pyerr)?;
    Ok((namelists.into_iter().map(|n| (n.name, n.entries)).collect(), cards))
}

/// `parse_qe_input_parts` 第一项的逆：写回 namelist 文本（值原样写出）。
#[pyfunction]
fn write_namelists(namelists: Vec<PyNamelist>) -> String {
    let nl: Vec<Namelist> = namelists
        .into_iter()
        .map(|(name, entries)| Namelist { name, entries })
        .collect();
    ferro_io::write_namelists(&nl)
}

pub(crate) fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sort_atoms, m)?)?;
    m.add_function(wrap_pyfunction!(poscar_string, m)?)?;
    m.add_function(wrap_pyfunction!(qe_structure_cards, m)?)?;
    m.add_function(wrap_pyfunction!(cp2k_structure_sections, m)?)?;
    m.add_function(wrap_pyfunction!(parse_qe_input_parts, m)?)?;
    m.add_function(wrap_pyfunction!(write_namelists, m)?)?;
    Ok(())
}
