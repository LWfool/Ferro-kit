use anyhow::{anyhow, Result};
use clap::Parser;
use ferro::io_dispatch::read_trajectory;
use ferro_io::LammpsUnits;
use ferro_workflow::GaussianJobBuilder;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "fe-job", about = "Generate QC software input files")]
struct Cli {
    /// Input structure file (format auto-detected)
    #[arg(short, long)]
    input: PathBuf,

    /// Target software: gaussian | gromacs
    #[arg(short, long)]
    software: String,

    /// DFT method (e.g. B3LYP, PBE)
    #[arg(short, long)]
    method: Option<String>,

    /// Basis set (e.g. 6-31G*)
    #[arg(short, long)]
    basis: Option<String>,

    /// Output file
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Use LAMMPS metal units for dump files (velocities Å/ps, forces eV/Å)
    #[arg(long)]
    metal_units: bool,
}

fn main() -> Result<()> {
    let args = Cli::parse();
    let units = if args.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let traj = read_trajectory(&args.input, units)?;
    let frame = traj.frames.into_iter().next()
        .ok_or_else(|| anyhow!("No frames in input file"))?;

    match args.software.to_lowercase().as_str() {
        "gaussian" => {
            let mut builder = GaussianJobBuilder::new(frame);
            if let Some(m) = args.method {
                builder.method = m;
            }
            if let Some(b) = args.basis {
                builder.basis_set = b;
            }
            let content = builder.build()?;
            let out = args.output.as_deref()
                .map(|p| p.to_string_lossy().into_owned())
                .unwrap_or_else(|| "job.gjf".to_string());
            std::fs::write(&out, content)?;
            println!("Gaussian input written to: {out}");
        }
        other => anyhow::bail!("Unsupported software: {other}  (supported: gaussian)"),
    }

    Ok(())
}
