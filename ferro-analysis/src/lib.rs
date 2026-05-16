//! ferro-analysis — 后处理分析模块

pub mod geometry;
pub mod trajectory_analysis;
pub mod md;
pub mod network;
pub mod dft;

pub use geometry::*;
pub use trajectory_analysis::*;
pub use network::{
    NetworkResult, calc_network,
    oxygen_label_order, modifier_label_order, former_label_order,
};
pub use ferro_core::{TypeParams, CutoffTable};
pub use dft::{BaderAnalyzer, BaderMethod, BaderParams, BaderResult,
              ChgSdfParams, ChgSdfResult, ChgSdfFamily, ChgRmsdStats, calc_chg_sdf};
pub use md::{
    GrParams, GrResult, PairStats, calc_gr, write_gr, write_cn,
    SqParams, SqResult, SqWeighting, calc_sq_from_gr, write_sq,
    MsdParams, MsdResult, calc_msd, write_msd,
    AngleParams, AngleResult, AngleStats, calc_angle, write_angle,
    VanHoveParams, VanHoveResult, calc_vanhove, write_vanhove,
    VacfParams, VacfResult, calc_vacf, write_vacf,
    RotCorrParams, RotCorrResult, calc_rotcorr, write_rotcorr,
    CubeMode, CubeDensityParams, CubeDensityResult, calc_cube_density,
    CubeRadiusParams, CubeRadiusResult, calc_cube_radius,
    ClusterSdfParams, ClusterFamily, RmsdStats, ClusterSdfResult, calc_cluster_sdf,
};
