use anyhow::{anyhow, Result};
use clap::Parser;
use ferro::{
    args::corr::CorrMode,
    help::print_corr_help,
    io_dispatch::read_trajectory,
};
use ferro_io::LammpsUnits;
use ferro_analysis::{
    calc_rotcorr, calc_vacf, calc_vanhove,
    write_rotcorr, write_vacf, write_vanhove,
    RotCorrParams, VacfParams, VanHoveParams,
};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(
    name = "mol-corr",
    about = "Correlation functions: vacf | rotcorr | vanhove  (run without -i for mode help)"
)]
struct Cli {
    /// Analysis mode
    #[arg(short = 'm', long, value_enum)]
    mode: CorrMode,

    /// Input trajectory file (omit to show mode-specific help)
    #[arg(short, long)]
    input: Option<PathBuf>,

    /// Output file
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Use only the last N frames
    #[arg(long)]
    last_n: Option<usize>,

    // ── shared time params ──────────────────────────────────────────────────

    /// Timestep between frames [fs]
    #[arg(long, default_value = "1.0")]
    dt: f64,

    /// Time-origin stride
    #[arg(long, default_value = "1")]
    shift: usize,

    /// Lag time in frames (default: half trajectory)
    #[arg(long)]
    tau: Option<usize>,

    // ── vacf / vanhove ──────────────────────────────────────────────────────

    /// [vacf/vanhove] Element filter, e.g. Fe,O
    #[arg(long, value_delimiter = ',')]
    elements: Option<Vec<String>>,

    // ── rotcorr ─────────────────────────────────────────────────────────────

    /// [rotcorr] Central atom element
    #[arg(long)]
    center: Option<String>,

    /// [rotcorr] Neighbor atom element
    #[arg(long)]
    neighbor: Option<String>,

    /// [rotcorr] Bond search cutoff [Å]
    #[arg(long, default_value = "1.2")]
    r_cut: f64,

    // ── vanhove ─────────────────────────────────────────────────────────────

    /// [vanhove] Max displacement [Å]
    #[arg(long, default_value = "10.0")]
    r_max: f64,

    /// [vanhove] Bin width [Å]
    #[arg(long, default_value = "0.01")]
    dr: f64,

    /// Use LAMMPS metal units for dump files (velocities Å/ps, forces eV/Å)
    #[arg(long)]
    metal_units: bool,
}

fn main() -> Result<()> {
    let args = Cli::parse();

    let input = match &args.input {
        Some(p) => p.clone(),
        None => {
            print_corr_help(&args.mode);
            return Ok(());
        }
    };

    let units = if args.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let mut traj = read_trajectory(&input, units)?;
    if let Some(n) = args.last_n {
        traj = traj.tail(n);
    }

    match args.mode {
        CorrMode::Vacf    => run_vacf(&args, &traj)?,
        CorrMode::Rotcorr => run_rotcorr(&args, &traj)?,
        CorrMode::Vanhove => run_vanhove(&args, &traj)?,
    }

    Ok(())
}

fn run_vacf(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let params = VacfParams {
        dt: args.dt,
        shift: args.shift,
        tau: args.tau,
        elements: args.elements.clone(),
        ..VacfParams::default()
    };
    let result = calc_vacf(traj, &params)
        .ok_or_else(|| anyhow!("VACF calc failed (missing velocities or empty trajectory?)"))?;

    let out = args.output.as_deref().unwrap_or(Path::new("vacf.dat"));
    write_vacf(&result, out.to_str().unwrap_or("vacf.dat"))?;
    println!("VACF -> {}", out.display());
    Ok(())
}

fn run_rotcorr(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let center = args.center.clone()
        .ok_or_else(|| anyhow!("--center is required for rotcorr (run without -i to see help)"))?;
    let neighbor = args.neighbor.clone()
        .ok_or_else(|| anyhow!("--neighbor is required for rotcorr (run without -i to see help)"))?;

    let params = RotCorrParams {
        center,
        neighbor,
        r_cut: args.r_cut,
        dt: args.dt,
        shift: args.shift,
        tau: args.tau,
    };
    let result = calc_rotcorr(traj, &params)
        .ok_or_else(|| anyhow!("RotCorr calc failed (no matching atom pairs found?)"))?;

    let out = args.output.as_deref().unwrap_or(Path::new("rotcorr.dat"));
    write_rotcorr(&result, out.to_str().unwrap_or("rotcorr.dat"))?;
    println!("RotCorr -> {}", out.display());
    Ok(())
}

fn run_vanhove(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let params = VanHoveParams {
        tau: args.tau,
        dt: args.dt,
        shift: args.shift,
        r_max: args.r_max,
        dr: args.dr,
        elements: args.elements.clone(),
        ..VanHoveParams::default()
    };
    let result = calc_vanhove(traj, &params)
        .ok_or_else(|| anyhow!("VanHove calc failed (trajectory too short?)"))?;

    let out = args.output.as_deref().unwrap_or(Path::new("vanhove.dat"));
    write_vanhove(&result, out.to_str().unwrap_or("vanhove.dat"))?;
    println!("VanHove -> {}", out.display());
    Ok(())
}
