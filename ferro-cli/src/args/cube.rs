use clap::ValueEnum;
use ferro_analysis::CubeMode;

#[derive(ValueEnum, Clone, Debug)]
pub enum CubeCliMode {
    Density,
    Velocity,
    Force,
    Radius,
    Sdf,
}

impl From<CubeCliMode> for CubeMode {
    fn from(m: CubeCliMode) -> Self {
        match m {
            CubeCliMode::Density  => CubeMode::Density,
            CubeCliMode::Velocity => CubeMode::Velocity,
            CubeCliMode::Force    => CubeMode::Force,
            // 以下两个走独立执行路径，不经此转换
            CubeCliMode::Radius   => CubeMode::Density,
            CubeCliMode::Sdf      => CubeMode::Density,
        }
    }
}
