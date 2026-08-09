use anyhow::{anyhow, Result};
use clap::Parser;
use ferro::{
    args::corr::CorrMode,
    batch::{self, Summary},
    help::{print_corr_help, print_fe_corr_overview},
    io_dispatch::read_trajectory_tail,
};
use ferro_io::LammpsUnits;
use ferro_analysis::{
    calc_rotcorr, calc_vacf, calc_vanhove,
    RotCorrParams, RotCorrResult, VacfParams, VacfResult, VanHoveParams, VanHoveResult,
};
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "fe-corr",
    about = "Correlation functions  (vacf | rotcorr | vanhove)",
    disable_help_flag = true,
)]
struct Cli {
    /// Analysis mode; omit to see overview
    #[arg(short = 'm', long, value_enum)]
    mode: Option<CorrMode>,

    /// Show help: overview when -m is absent, mode-specific when -m is given
    #[arg(short = 'h', long = "help", action = clap::ArgAction::SetTrue)]
    help: bool,

    /// Input trajectory file(s); glob patterns allowed (quote them)
    #[arg(short, long, num_args = 1.., value_name = "FILE")]
    input: Vec<PathBuf>,

    /// Output name suffix: results go to <mode>_<suffix>.csv
    #[arg(short, long, value_name = "SUFFIX")]
    output: Option<String>,

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

    // 无 -m → 概览
    let mode = match args.mode.clone() {
        Some(m) => m,
        None => {
            print_fe_corr_overview();
            return Ok(());
        }
    };

    // -h 或无 -i → 模式专属帮助
    if args.help || args.input.is_empty() {
        print_corr_help(&mode);
        return Ok(());
    }

    let inputs = batch::expand_inputs(&args.input)?;
    println!("Inputs: {} file(s)", inputs.len());

    let failures = match mode {
        CorrMode::Vacf    => run_vacf(&args, &inputs)?,
        CorrMode::Rotcorr => run_rotcorr(&args, &inputs)?,
        CorrMode::Vanhove => run_vanhove(&args, &inputs)?,
    };

    if failures > 0 {
        eprintln!("{failures} input(s) failed; see the [inputs] block in the output file");
        std::process::exit(1);
    }
    Ok(())
}

/// 读一条轨迹并按 --last-n 裁剪。
fn load(args: &Cli, path: &std::path::Path) -> Result<ferro_core::Trajectory> {
    let units = if args.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    read_trajectory_tail(path, units, args.last_n)
}

fn run_vacf(args: &Cli, inputs: &[PathBuf]) -> Result<usize> {
    let params = VacfParams {
        dt: args.dt,
        shift: args.shift,
        tau: args.tau,
        elements: args.elements.clone(),
    };

    let (results, failures) = batch::map_inputs(inputs, |p| {
        let traj = load(args, p)?;
        calc_vacf(&traj, &params)
            .ok_or_else(|| anyhow!("VACF calc failed (missing velocities or empty trajectory?)"))
    });
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }

    let tables = batch::stack(&results, |r: &VacfResult| Ok(r.to_tables()))?;
    let mut summary = Summary::new(&["origins"]);
    for (path, r) in &results {
        summary.ok(batch::label_of(path), r.time.len(), r.n_atoms, &[r.n_origins as f64]);
    }
    summary.failed(&failures);

    batch::write_all(
        "vacf",
        "Velocity Autocorrelation Function (VACF)",
        &results[0].1.meta_lines(),
        &summary.into_table(),
        tables,
        args.output.as_deref(),
    )?;
    Ok(failures.len())
}

fn run_rotcorr(args: &Cli, inputs: &[PathBuf]) -> Result<usize> {
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

    let (results, failures) = batch::map_inputs(inputs, |p| {
        let traj = load(args, p)?;
        calc_rotcorr(&traj, &params)
            .ok_or_else(|| anyhow!("RotCorr calc failed (no matching atom pairs found?)"))
    });
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }

    let tables = batch::stack(&results, |r: &RotCorrResult| Ok(r.to_tables()))?;
    let mut summary = Summary::new(&["origins"]);
    for (path, r) in &results {
        summary.ok(batch::label_of(path), r.time.len(), r.n_molecules, &[r.n_origins as f64]);
    }
    summary.failed(&failures);

    batch::write_all(
        "rotcorr",
        "Rotational Autocorrelation Function C(t) = <P2(cos theta)>",
        &results[0].1.meta_lines(),
        &summary.into_table(),
        tables,
        args.output.as_deref(),
    )?;
    Ok(failures.len())
}

fn run_vanhove(args: &Cli, inputs: &[PathBuf]) -> Result<usize> {
    let params = VanHoveParams {
        tau: args.tau,
        dt: args.dt,
        shift: args.shift,
        r_max: args.r_max,
        dr: args.dr,
        elements: args.elements.clone(),
        ..VanHoveParams::default()
    };

    let (results, failures) = batch::map_inputs(inputs, |p| {
        let traj = load(args, p)?;
        calc_vanhove(&traj, &params)
            .ok_or_else(|| anyhow!("VanHove calc failed (trajectory too short?)"))
    });
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }

    let tables = batch::stack(&results, |r: &VanHoveResult| Ok(r.to_tables()))?;
    // tau 逐文件相同,但 time = tau*dt 与 origins 值得横向看一眼
    let mut summary = Summary::new(&["tau_frames", "time_fs", "origins"]);
    for (path, r) in &results {
        summary.ok(
            batch::label_of(path),
            r.r.len(),
            r.n_atoms,
            &[r.tau_frames as f64, r.time, r.n_origins as f64],
        );
    }
    summary.failed(&failures);

    batch::write_all(
        "vanhove",
        "van Hove Self-Correlation Function Gs(r, tau)",
        &results[0].1.meta_lines(),
        &summary.into_table(),
        tables,
        args.output.as_deref(),
    )?;
    Ok(failures.len())
}
