//! Sub-modules for MD trajectory post-processing analysis.
//!
//! Current implementations:
//!   - [`gr`]       — radial distribution function g(r) and coordination number CN(r)
//!   - [`sq`]       — structure factor S(q) via Fourier transform of g(r)
//!   - [`msd`]      — mean squared displacement (time-shift averaging, NPT-safe)
//!   - [`angle`]    — bond angle distribution
//!   - [`vanhove`]  — van Hove self-correlation function Gs(r, τ)
//!
pub mod gr;
pub mod sq;
pub mod msd;
pub mod angle;
pub mod vanhove;
pub mod vacf;
pub mod rotcorr;
pub mod cube_density;
pub mod cube_jump;
pub mod cube_radius;
pub mod cube_sdf;
pub mod scattering_data;

pub use gr::{GrParams, GrResult, PairStats, calc_gr, write_gr, write_cn};
pub use sq::{SqParams, SqResult, SqWeighting, calc_sq_from_gr, write_sq};
pub use msd::{MsdParams, MsdResult, calc_msd, write_msd};
pub use angle::{AngleParams, AngleResult, AngleStats, calc_angle, write_angle};
pub use vanhove::{VanHoveParams, VanHoveResult, calc_vanhove, write_vanhove};
pub use vacf::{VacfParams, VacfResult, calc_vacf, write_vacf};
pub use rotcorr::{RotCorrParams, RotCorrResult, calc_rotcorr, write_rotcorr};
pub use cube_density::{CubeMode, CubeDensityParams, CubeDensityResult, calc_cube_density};
pub use cube_jump::{JumpPosition, CubeJumpParams, CubeJumpResult, calc_cube_jump};
pub use cube_radius::{CubeRadiusParams, CubeRadiusResult, calc_cube_radius};
pub use cube_sdf::{ClusterSdfParams, ClusterFamily, RmsdStats, ClusterSdfResult, calc_cluster_sdf};
