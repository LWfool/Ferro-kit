use clap::ValueEnum;
use ferro_analysis::SqWeighting;

#[derive(ValueEnum, Clone, Debug)]
pub enum TrajMode {
    Gr,
    Sq,
    Msd,
    Angle,
}

/// CLI-side weighting enum; converts into `SqWeighting`.
#[derive(ValueEnum, Clone, Debug, Default)]
pub enum SqWeightingCli {
    None,
    Xrd,
    Neutron,
    #[default]
    Both,
}

impl From<SqWeightingCli> for SqWeighting {
    fn from(w: SqWeightingCli) -> Self {
        match w {
            SqWeightingCli::None    => SqWeighting::None,
            SqWeightingCli::Xrd     => SqWeighting::Xrd,
            SqWeightingCli::Neutron => SqWeighting::Neutron,
            SqWeightingCli::Both    => SqWeighting::Both,
        }
    }
}
