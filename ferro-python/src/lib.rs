//! ferro — Python 绑定（PyO3）。
//!
//! 纯 Rust 库（ferro-core / io / structure / analysis）对 Python 零感知；
//! 本 crate 是唯一的 PyO3 胶水层。用 maturin 构建：`maturin develop`。
//!
//! ```python
//! import ferro
//! t = ferro.read("traj.lammpstrj", units="metal")
//! sc = ferro.supercell(t, 2, 2, 1)
//! g  = ferro.gr(t, r_max=10.0, dr=0.02)   # dict[str, list[float]]
//! d  = ferro.msd(t, dt=2.0, elements=["Li"])
//! ```

use pyo3::prelude::*;

mod analysis;
mod doc;
mod io;
mod jobs;
mod structure;
mod types;

pyo3::create_exception!(
    ferro,
    FerroError,
    pyo3::exceptions::PyRuntimeError,
    "An error reported by ferro's Rust libraries (a file that cannot be read, a missing cell, ...). \
     A subclass of RuntimeError, so `except RuntimeError` written before it existed still catches it."
);

/// Turns any displayable error (anyhow::Error, ferro_core::ChemError, ...) into a `ferro.FerroError`:
/// its own class, so a caller catches what ferro reports without swallowing Python's own RuntimeErrors.
pub(crate) fn pyerr<E: std::fmt::Display>(e: E) -> PyErr {
    // `{:#}` keeps anyhow's whole context chain; to_string() keeps only the outer "parsing FILE"
    FerroError::new_err(format!("{e:#}"))
}

#[pymodule]
fn _ferro(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<types::PyTrajectory>()?;
    m.add("FerroError", m.py().get_type::<FerroError>())?;
    io::register(m)?;
    jobs::register(m)?;
    structure::register(m)?;
    analysis::register(m)?;
    doc::register(m)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
