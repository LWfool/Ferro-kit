use clap::Args;
use std::path::PathBuf;

#[derive(Args, Clone, Debug)]
pub struct CommonArgs {
    /// Input file (format auto-detected from extension / filename)
    #[arg(short, long)]
    pub input: Option<PathBuf>,

    /// Output file
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Use only the last N frames (skip equilibration)
    #[arg(long)]
    pub last_n: Option<usize>,

    /// Parallel threads (default: all CPU cores)
    #[arg(long)]
    pub ncore: Option<usize>,
}
