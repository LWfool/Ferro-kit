//! DFT post-processing analysis modules.
//!
//! Current implementations:
//!   - [`bader`] — Bader charge analysis (Henkelman method)
//!   - [`chg_sdf`] — Averaged charge-density SDF over aligned Qn clusters

pub mod bader;
pub mod bader_grid;
pub mod bader_weight;
pub mod chg_sdf;

pub use bader::{BaderAnalyzer, BaderMethod, BaderParams, BaderResult};
pub use chg_sdf::{ChgSdfParams, ChgSdfResult, ChgSdfFamily, ChgRmsdStats, calc_chg_sdf};
