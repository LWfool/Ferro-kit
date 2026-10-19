use clap::ValueEnum;
use ferro_analysis::CubeMode;

#[derive(ValueEnum, Clone, Debug)]
pub enum CubeCliMode {
    Density,
    Velocity,
    Force,
}

impl From<CubeCliMode> for CubeMode {
    fn from(m: CubeCliMode) -> Self {
        match m {
            CubeCliMode::Density  => CubeMode::Density,
            CubeCliMode::Velocity => CubeMode::Velocity,
            CubeCliMode::Force    => CubeMode::Force,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every CLI mode maps to the library mode of the same name; no placeholder arm
    #[test]
    fn each_mode_maps_to_its_namesake() {
        for m in CubeCliMode::value_variants() {
            let name = m.to_possible_value().unwrap().get_name().to_string();
            assert_eq!(format!("{:?}", CubeMode::from(m.clone())).to_lowercase(), name);
        }
        assert_eq!(CubeCliMode::value_variants().len(), 3);
    }
}
