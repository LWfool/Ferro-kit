use clap::ValueEnum;

#[derive(ValueEnum, Clone, Debug)]
pub enum CorrMode {
    Vacf,
    Rotcorr,
    Vanhove,
}
