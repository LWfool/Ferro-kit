//! ferro-analysis — 后处理分析模块

mod check;
mod neighbor;
// One public path per item (review A-12): the analyses are reached from the crate
// root; `ml` stays a module because its parts are an API of their own
mod md;
mod network;
mod dft;
pub mod ml;

pub use network::{NetworkResult, calc_network, calc_network_typed, AtomType, CutoffTable, TypeParams};
pub use dft::{BaderAnalyzer, BaderMethod, BaderParams, BaderResult,
              ChgSdfParams, ChgSdfResult, ChgSdfFamily, ChgRmsdStats, calc_chg_sdf};
pub use md::{
    GrParams, GrResult, GroupBy, calc_gr,
    SqParams, SqResult, SqWeighting, calc_sq_from_gr,
    MsdParams, MsdResult, calc_msd,
    AngleParams, AngleResult, AngleStats, calc_angle,
    VanHoveParams, VanHoveResult, calc_vanhove,
    VacfParams, VacfResult, calc_vacf,
    Legendre, RotCorrParams, RotCorrResult, RotVector, calc_rotcorr,
    BondLifeParams, BondLifeResult, calc_bondlife,
    CubeMode, CubeDensityParams, CubeDensityResult, calc_cube_density,
    CubeRadiusParams, CubeRadiusResult, calc_cube_radius,
    ClusterSdfParams, ClusterFamily, RmsdStats, ClusterSdfResult, calc_cluster_sdf,
    ShellCutoff, first_shell_cutoffs,
    JumpPosition, CubeJumpParams, CubeJumpResult, calc_cube_jump,
};

/// ferro-analysis version (from Cargo.toml)
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
