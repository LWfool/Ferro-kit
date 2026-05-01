use anyhow::{anyhow, Result};
use clap::Parser;
use ferro::{
    args::traj::{SqWeightingCli, TrajMode},
    help::{print_fe_traj_overview, print_traj_help},
    io_dispatch::read_trajectory,
    plot::{open_plot, plot_angle, plot_gr, plot_sq},
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
    name = "fe-traj",
    about = "Trajectory structural analysis  (gr | sq | msd | angle)",
    disable_help_flag = true,
)]
struct Cli {
    /// Analysis mode  (gr | sq | msd | angle); omit to see overview
    #[arg(short = 'm', long, value_enum)]
    mode: Option<TrajMode>,

    /// Show help: overview when -m is absent, mode-specific when -m is given
    #[arg(short = 'h', long = "help", action = clap::ArgAction::SetTrue)]
    help: bool,

    /// Input trajectory file
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

    // ── pair / triplet filter ────────────────────────────────────────────────

    /// [gr] Element A of the pair (e.g. O); [angle] end atom A — requires -b
    #[arg(short = 'a', long)]
    atom_a: Option<String>,

    /// [gr] Element B of the pair (e.g. P); [angle] center atom B — requires -a
    #[arg(short = 'b', long)]
    atom_b: Option<String>,

    /// [angle] End atom C (e.g. O); corresponds to --r-cut-bc — requires -a -b
    #[arg(short = 'c', long)]
    atom_c: Option<String>,

    // ── angle ───────────────────────────────────────────────────────────────

    /// [angle] Cutoff for A-to-center-B bond [Å]; A is the atom given by -a
    #[arg(long, default_value = "2.3")]
    r_cut_ab: f64,

    /// [angle] Cutoff for center-B-to-C bond [Å]; C is the atom given by -c
    #[arg(long, default_value = "2.3")]
    r_cut_bc: f64,

    /// [angle] Histogram bin width [°]
    #[arg(long, default_value = "0.1")]
    d_angle: f64,

    /// Use LAMMPS metal units for dump files (velocities Å/ps, forces eV/Å)
    #[arg(long)]
    metal_units: bool,

    /// Generate a PNG plot and open it after calculation
    #[arg(long)]
    plot: bool,
}

fn main() -> Result<()> {
    let args = Cli::parse();

    // 无 -m → 概览
    let mode = match args.mode.clone() {
        Some(m) => m,
        None => {
            print_fe_traj_overview();
            return Ok(());
        }
    };

    // -h 或无 -i → 模式专属帮助
    if args.help || args.input.is_none() {
        print_traj_help(&mode);
        return Ok(());
    }

    let input = args.input.as_ref().unwrap().clone();

    if let Some(n) = args.ncore {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
    }

    let units = if args.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let mut traj = read_trajectory(&input, units)?;
    if let Some(n) = args.last_n {
        traj = traj.tail(n);
    }

    match mode {
        TrajMode::Gr    => run_gr(&args, &traj)?,
        TrajMode::Sq    => run_sq(&args, &traj)?,
        TrajMode::Msd   => run_msd(&args, &traj)?,
        TrajMode::Angle => run_angle(&args, &traj)?,
    }

    Ok(())
}

fn run_gr(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    match (&args.atom_a, &args.atom_b) {
        (Some(_), None) | (None, Some(_)) =>
            return Err(anyhow!("GR pair filter: -a and -b must be specified together")),
        _ => {}
    }

    let params = GrParams {
        r_max: args.r_max,
        dr: args.dr,
        r_cut: args.r_cut,
        ..GrParams::default()
    };
    let mut result = calc_gr(traj, &params)
        .ok_or_else(|| anyhow!("GR calc failed (empty trajectory or missing cell)"))?;

    // 指定原子对时过滤输出列，保留 total
    if let (Some(a), Some(b)) = (&args.atom_a, &args.atom_b) {
        let ab = format!("{a}-{b}");
        let ba = format!("{b}-{a}");
        let keep_sym = |k: &str| k == "total" || k == ab || k == ba;
        result.gr.retain(|k, _| keep_sym(k));
        result.cn.retain(|k, _| keep_sym(k));
        result.pair_stats.retain(|k, _| keep_sym(k));
        if result.gr.len() <= 1 {
            return Err(anyhow!("pair '{ab}' not found in trajectory (check element symbols)"));
        }
    }

    let out = args.output.as_deref().unwrap_or(Path::new("gr.dat"));
    let out_str = out.to_str().unwrap_or("gr.dat");
    write_gr(&result, out_str)?;

    // 同时输出配位数文件
    let cn_path = sibling_path(out, "_cn");
    write_cn(&result, &cn_path)?;

    println!("GR  -> {out_str}");
    println!("CN  -> {cn_path}");

    if args.plot {
        let png = plot_gr(&result, out_str)?;
        println!("Plot-> {png}");
        open_plot(&png);
    }
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
    let out_str = out.to_str().unwrap_or("sq.dat");
    write_sq(&gr, &sq, out_str)?;
    println!("SQ  -> {out_str}");

    if args.plot {
        let png = plot_sq(&sq, out_str)?;
        println!("Plot-> {png}");
        open_plot(&png);
    }
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
    let n_filter = [&args.atom_a, &args.atom_b, &args.atom_c].iter().filter(|x| x.is_some()).count();
    if n_filter > 0 && n_filter < 3 {
        return Err(anyhow!("Angle filter: -a (end A), -b (center B), -c (end C) must all be specified together"));
    }

    let params = AngleParams {
        r_cut_ab: args.r_cut_ab,
        r_cut_bc: args.r_cut_bc,
        d_angle: args.d_angle,
    };
    let mut result = calc_angle(traj, &params)
        .ok_or_else(|| anyhow!("Angle calc failed (empty trajectory?)"))?;

    // 指定三元组时过滤输出列（key 格式 "A-Center-C"，端原子按 Z 排序，两种顺序均检查）
    if let (Some(a), Some(b), Some(c)) = (&args.atom_a, &args.atom_b, &args.atom_c) {
        let key1 = format!("{a}-{b}-{c}");
        let key2 = format!("{c}-{b}-{a}");
        result.hist.retain(|k, _| k == &key1 || k == &key2);
        result.stats.retain(|k, _| k == &key1 || k == &key2);
        if result.hist.is_empty() {
            return Err(anyhow!(
                "triplet '{key1}' not found (check element symbols and that -b is the center atom)"
            ));
        }
    }

    let out = args.output.as_deref().unwrap_or(Path::new("angle.dat"));
    let out_str = out.to_str().unwrap_or("angle.dat");
    write_angle(&result, out_str)?;
    println!("Angle -> {out_str}");

    if args.plot {
        let png = plot_angle(&result, out_str)?;
        println!("Plot -> {png}");
        open_plot(&png);
    }
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
