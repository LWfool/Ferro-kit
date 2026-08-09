use anyhow::{bail, Result};
use clap::Args;
use crate::io_dispatch::{read_trajectory, supported_formats, write_trajectory};
use ferro_io::LammpsUnits;
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct ConvertCmd {
    /// Input file (format auto-detected)
    #[arg(short, long)]
    pub input: PathBuf,

    /// Output file (format auto-detected)
    #[arg(short, long)]
    pub output: PathBuf,

    /// Use LAMMPS metal units for dump files (velocities Å/ps, forces eV/Å)
    #[arg(long)]
    pub metal_units: bool,
}

pub fn run(args: &ConvertCmd) -> Result<()> {

    if args.input == args.output {
        bail!("Input and output paths are identical");
    }

    let units = if args.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let traj = read_trajectory(&args.input, units)?;
    let n = traj.frames.len();
    write_trajectory(&traj, &args.output, units)?;

    println!(
        "Converted {} ({} frame{}) -> {}",
        args.input.display(),
        n,
        if n == 1 { "" } else { "s" },
        args.output.display()
    );
    Ok(())
}

#[allow(dead_code)]
fn _formats() {
    println!("{}", supported_formats());
}
