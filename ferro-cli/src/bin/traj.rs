use anyhow::{anyhow, Result};
use clap::Parser;
use ferro::{
    args::traj::{SqWeightingCli, TrajMode},
    help::print_traj_help,
    io_dispatch::read_trajectory,
};
use ferro_io::LammpsUnits;
use ferro_analysis::{
    calc_angle, calc_gr, calc_msd, calc_sq_from_gr,
    write_angle, write_cn, write_gr, write_msd, write_sq,
    AngleParams, GrParams, MsdParams, SqParams,
};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(
    name = "mol-traj",
    about = "Structural analysis: gr | sq | msd | angle  (run without -i for mode help)"
)]
struct Cli {
    /// Analysis mode
    #[arg(short = 'm', long, value_enum)]
    mode: TrajMode,

    /// Input trajectory file (omit to show mode-specific help)
    #[arg(short, long)]
    input: Option<PathBuf>,

    /// Output file
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Use only the last N frames
    #[arg(long)]
    last_n: Option<usize>,

    /// Parallel threads (default: all cores)
    #[arg(long)]
    ncore: Option<usize>,

    // ── gr / sq shared ──────────────────────────────────────────────────────

    /// [gr/sq] Max cutoff radius [Å]
    #[arg(long, default_value = "10.005")]
    r_max: f64,

    /// [gr/sq] Histogram bin width [Å]
    #[arg(long, default_value = "0.01")]
    dr: f64,

    /// [gr] CN integration cutoff [Å]
    #[arg(long, default_value = "2.3")]
    r_cut: f64,

    // ── sq ──────────────────────────────────────────────────────────────────

    /// [sq] Max q [Å⁻¹]
    #[arg(long, default_value = "25.0")]
    q_max: f64,

    /// [sq] q bin width [Å⁻¹]
    #[arg(long, default_value = "0.05")]
    dq: f64,

    /// [sq] Scattering weighting scheme
    #[arg(long, value_enum, default_value = "both")]
    weighting: SqWeightingCli,

    // ── msd ─────────────────────────────────────────────────────────────────

    /// [msd] Timestep between frames [fs]
    #[arg(long, default_value = "1.0")]
    dt: f64,

    /// [msd] Time-origin stride
    #[arg(long, default_value = "1")]
    shift: usize,

    /// [msd] Track only these elements, e.g. Fe,O
    #[arg(long, value_delimiter = ',')]
    elements: Option<Vec<String>>,

    // ── angle ───────────────────────────────────────────────────────────────

    /// [angle] A-to-B distance cutoff [Å]
    #[arg(long, default_value = "2.3")]
    r_cut_ab: f64,

    /// [angle] C-to-B distance cutoff [Å]
    #[arg(long, default_value = "2.3")]
    r_cut_bc: f64,

    /// [angle] Histogram bin width [°]
    #[arg(long, default_value = "0.1")]
    d_angle: f64,

    /// Use LAMMPS metal units for dump files (velocities Å/ps, forces eV/Å)
    #[arg(long)]
    metal_units: bool,
}

fn main() -> Result<()> {
    let args = Cli::parse();

    // 没有输入文件时显示模式专属帮助
    let input = match &args.input {
        Some(p) => p.clone(),
        None => {
            print_traj_help(&args.mode);
            return Ok(());
        }
    };

    if let Some(n) = args.ncore {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
    }

    let units = if args.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let mut traj = read_trajectory(&input, units)?;
    if let Some(n) = args.last_n {
        traj = traj.tail(n);
    }

    match args.mode {
        TrajMode::Gr => run_gr(&args, &traj)?,
        TrajMode::Sq => run_sq(&args, &traj)?,
        TrajMode::Msd => run_msd(&args, &traj)?,
        TrajMode::Angle => run_angle(&args, &traj)?,
    }

    Ok(())
}

fn run_gr(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let params = GrParams {
        r_max: args.r_max,
        dr: args.dr,
        r_cut: args.r_cut,
        ..GrParams::default()
    };
    let result = calc_gr(traj, &params).ok_or_else(|| anyhow!("GR calc failed (empty trajectory or missing cell)"))?;

    let out = args.output.as_deref().unwrap_or(Path::new("gr.dat"));
    let out_str = out.to_str().unwrap_or("gr.dat");
    write_gr(&result, out_str)?;

    // 同时输出配位数文件
    let cn_path = sibling_path(out, "_cn");
    write_cn(&result, &cn_path)?;

    println!("GR  -> {out_str}");
    println!("CN  -> {cn_path}");
    Ok(())
}

fn run_sq(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let gr_params = GrParams {
        r_max: args.r_max,
        dr: args.dr,
        r_cut: args.r_cut,
        ..GrParams::default()
    };
    let gr = calc_gr(traj, &gr_params).ok_or_else(|| anyhow!("GR calc failed (needed for SQ)"))?;

    let sq_params = SqParams {
        q_max: args.q_max,
        dq: args.dq,
        weighting: args.weighting.clone().into(),
        ..SqParams::default()
    };
    let sq = calc_sq_from_gr(&gr, &sq_params);

    let out = args.output.as_deref().unwrap_or(Path::new("sq.dat"));
    write_sq(&gr, &sq, out.to_str().unwrap_or("sq.dat"))?;
    println!("SQ  -> {}", out.display());
    Ok(())
}

fn run_msd(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let params = MsdParams {
        dt: args.dt,
        shift: args.shift,
        elements: args.elements.clone(),
        ..MsdParams::default()
    };
    let result = calc_msd(traj, &params).ok_or_else(|| anyhow!("MSD calc failed (trajectory too short?)"))?;

    let out = args.output.as_deref().unwrap_or(Path::new("msd.dat"));
    write_msd(&result, out.to_str().unwrap_or("msd.dat"))?;
    println!("MSD -> {}", out.display());
    Ok(())
}

fn run_angle(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let params = AngleParams {
        r_cut_ab: args.r_cut_ab,
        r_cut_bc: args.r_cut_bc,
        d_angle: args.d_angle,
    };
    let result = calc_angle(traj, &params).ok_or_else(|| anyhow!("Angle calc failed (empty trajectory?)"))?;

    let out = args.output.as_deref().unwrap_or(Path::new("angle.dat"));
    write_angle(&result, out.to_str().unwrap_or("angle.dat"))?;
    println!("Angle -> {}", out.display());
    Ok(())
}

/// Build a sibling output path: `gr.dat` + `_cn` → `gr_cn.dat`
fn sibling_path(base: &Path, suffix: &str) -> String {
    let stem = base.file_stem().and_then(|s| s.to_str()).unwrap_or("gr");
    let ext  = base.extension().and_then(|e| e.to_str()).unwrap_or("dat");
    let dir  = base.parent().map(|p| p.to_str().unwrap_or("")).unwrap_or("");
    if dir.is_empty() {
        format!("{stem}{suffix}.{ext}")
    } else {
        format!("{dir}/{stem}{suffix}.{ext}")
    }
}
