//! ferro — Python 绑定（PyO3）。
//!
//! 纯 Rust 库（ferro-core / io / structure / analysis）对 Python 零感知；
//! 本 crate 是唯一的 PyO3 胶水层。用 maturin 构建：`maturin develop`。
//!
//! ```python
//! import ferro
//! t = ferro.read("traj.lammpstrj", metal_units=True)
//! sc = ferro.supercell(t, 2, 2, 1)
//! g  = ferro.gr(t, r_max=10.0, dr=0.02)   # dict[str, list[float]]
//! d  = ferro.msd(t, dt=2.0, elements=["Li"])
//! ```

use pyo3::prelude::*;

mod analysis;
mod io;
mod structure;
mod types;

/// 把任意可显示错误（anyhow::Error / ferro_core::ChemError 等）转为 Python 异常。
pub(crate) fn pyerr<E: std::fmt::Display>(e: E) -> PyErr {
    pyo3::exceptions::PyRuntimeError::new_err(e.to_string())
}

#[pymodule]
fn ferro(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<types::PyTrajectory>()?;
    io::register(m)?;
    structure::register(m)?;
    analysis::register(m)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
