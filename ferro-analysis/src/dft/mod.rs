//! DFT post-processing analysis modules.
//!
//! Current implementations:
//!   - [`bader`] — Bader charge analysis (Henkelman method)

pub mod bader;
pub mod bader_grid;
pub mod bader_weight;

pub use bader::{BaderAnalyzer, BaderMethod, BaderParams, BaderResult};
