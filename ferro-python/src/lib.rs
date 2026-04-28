//! ChemTools Python 绑定
//! 
//! 将 Rust 库导出为 Python 模块

use pyo3::prelude::*;

/// 读取 XYZ 文件
#[pyfunction]
fn read_xyz(path: String) -> PyResult<String> {
    use ferro_io::read_xyz;
    
    let mol = read_xyz(&path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    
    Ok(format!(
        "分子包含 {} 个原子，总质量 {:.2} amu",
        mol.atom_count(),
        mol.total_mass()
    ))
}

/// 写入 XYZ 文件
#[pyfunction]
fn write_xyz(input_path: String, output_path: String) -> PyResult<()> {
    use ferro_io::{read_xyz, write_xyz};
    
    let mol = read_xyz(&input_path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    
    write_xyz(&mol, &output_path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    
    Ok(())
}

/// 计算分子质心
#[pyfunction]
fn center_of_mass(path: String) -> PyResult<(f64, f64, f64)> {
    use ferro_io::read_xyz;
    
    let mol = read_xyz(&path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    
    let com = mol.center_of_mass();
    Ok((com.x, com.y, com.z))
}

/// Python 模块定义
#[pymodule]
fn ferro(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_xyz, m)?)?;
    m.add_function(wrap_pyfunction!(write_xyz, m)?)?;
    m.add_function(wrap_pyfunction!(center_of_mass, m)?)?;
    
    Ok(())
}
