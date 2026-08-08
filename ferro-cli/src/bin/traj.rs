use anyhow::{anyhow, Result};
use clap::Parser;
use ferro::{
    args::traj::{SqWeightingCli, TrajMode},
    help::{print_fe_traj_overview, print_traj_help},
    io_dispatch::read_trajectory,
    plot::{open_plot, plot_angle, plot_gr, plot_msd, plot_sq},
};
use ferro_io::LammpsUnits;
use ferro_analysis::{
    calc_angle, calc_gr, calc_msd, calc_sq_from_gr,
    write_angle, write_gr, write_msd, write_sq,
    AngleParams, GrParams, GroupBy, MsdParams, SqParams,
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

    /// [msd] Linear-fit window as trajectory fractions FMIN,FMAX (e.g. 0.3,0.8) -> self-diffusion D
    #[arg(long, value_delimiter = ',')]
    fit_range: Option<Vec<f64>>,

    // ── pair / triplet filter ────────────────────────────────────────────────

    /// [gr/sq] Centre element (e.g. P); [angle] end atom A — by element, requires -b
    #[arg(short = 'a', long)]
    atom_a: Option<String>,

    /// [gr/sq] Neighbour element (e.g. O); [angle] centre atom B — by element, requires -a
    #[arg(short = 'b', long)]
    atom_b: Option<String>,

    /// [angle] End atom C by element; corresponds to --r-cut-bc — requires -a -b
    #[arg(short = 'c', long)]
    atom_c: Option<String>,

    /// [gr/sq] Centre site label (e.g. P_0); [angle] end atom A — by label, requires -y
    #[arg(short = 'x', long)]
    label_x: Option<String>,

    /// [gr/sq] Neighbour site label (e.g. O_b_P_P); [angle] centre atom B — by label, requires -x
    #[arg(short = 'y', long)]
    label_y: Option<String>,

    /// [angle] End atom C by site label; corresponds to --r-cut-bc — requires -x -y
    #[arg(short = 'z', long)]
    label_z: Option<String>,

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

/// 解析选择参数：`-a/-b/-c` 按 element，`-x/-y/-z` 按 label，两组互斥。
///
/// 返回分组模式与三个槽位（angle 用满三个，gr/sq 只用前两个）。
fn resolve_selection(args: &Cli) -> Result<(GroupBy, [Option<String>; 3])> {
    let by_elem = [&args.atom_a, &args.atom_b, &args.atom_c];
    let by_label = [&args.label_x, &args.label_y, &args.label_z];
    let n_elem = by_elem.iter().filter(|x| x.is_some()).count();
    let n_label = by_label.iter().filter(|x| x.is_some()).count();

    if n_elem > 0 && n_label > 0 {
        return Err(anyhow!(
            "-a/-b/-c (by element) and -x/-y/-z (by label) are mutually exclusive; use one group"
        ));
    }
    let slots = if n_label > 0 { by_label } else { by_elem };
    let mode = if n_label > 0 { GroupBy::Label } else { GroupBy::Element };
    Ok((mode, [slots[0].clone(), slots[1].clone(), slots[2].clone()]))
}

/// gr / sq 的配对选择：要么两个都不给，要么中心与近邻都给。
fn resolve_pair(args: &Cli) -> Result<(GroupBy, Option<(String, String)>)> {
    let (mode, slots) = resolve_selection(args)?;
    if slots[2].is_some() {
        return Err(anyhow!("-c / -z applies to -m angle only"));
    }
    match (&slots[0], &slots[1]) {
        (Some(a), Some(b)) => Ok((mode, Some((a.clone(), b.clone())))),
        (None, None) => Ok((mode, None)),
        _ => Err(anyhow!(
            "pair filter needs both members: -a CENTRE -b NEIGHBOUR (by element) \
             or -x CENTRE -y NEIGHBOUR (by label)"
        )),
    }
}

fn run_gr(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let (group_by, pair) = resolve_pair(args)?;

    let params = GrParams {
        r_max: args.r_max,
        dr: args.dr,
        group_by,
        ..GrParams::default()
    };
    let result = calc_gr(traj, &params)?;

    let out = args.output.as_deref().unwrap_or(Path::new("gr.dat"));
    let out_str = out.to_str().unwrap_or("gr.dat");
    let pair_ref = pair.as_ref().map(|(a, b)| (a.as_str(), b.as_str()));
    write_gr(&result, out_str, pair_ref)?;
    println!("GR/CN -> {out_str}");

    if args.plot {
        match pair_ref {
            Some((a, b)) => {
                let png = plot_gr(&result, out_str, a, b)?;
                println!("Plot  -> {png}");
                open_plot(&png);
            }
            None => println!(
                "Note  : --plot needs a pair; add -a CENTRE -b NEIGHBOUR (or -x/-y) to plot one"
            ),
        }
    }
    Ok(())
}

fn run_sq(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let (group_by, pair) = resolve_pair(args)?;

    let gr_params = GrParams {
        r_max: args.r_max,
        dr: args.dr,
        group_by,
        ..GrParams::default()
    };
    let gr = calc_gr(traj, &gr_params)?;

    let sq_params = SqParams {
        q_max: args.q_max,
        dq: args.dq,
        weighting: args.weighting.clone().into(),
        ..SqParams::default()
    };
    let sq = calc_sq_from_gr(&gr, &sq_params);

    let out = args.output.as_deref().unwrap_or(Path::new("sq.dat"));
    let out_str = out.to_str().unwrap_or("sq.dat");
    let pair_ref = pair.as_ref().map(|(a, b)| (a.as_str(), b.as_str()));
    write_sq(&gr, &sq, out_str, pair_ref)?;
    println!("SQ  -> {out_str}");

    if args.plot {
        let png = plot_sq(&sq, out_str)?;
        println!("Plot-> {png}");
        open_plot(&png);
    }
    Ok(())
}

fn run_msd(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let fit_range = match &args.fit_range {
        None => None,
        Some(v) if v.len() == 2 => Some((v[0], v[1])),
        Some(v) => return Err(anyhow!(
            "--fit-range expects exactly two comma-separated fractions, e.g. 0.3,0.8 (got {} value(s))",
            v.len()
        )),
    };

    let params = MsdParams {
        dt: args.dt,
        shift: args.shift,
        elements: args.elements.clone(),
        fit_range,
        ..MsdParams::default()
    };
    let result = calc_msd(traj, &params)?;

    let out = args.output.as_deref().unwrap_or(Path::new("msd.dat"));
    let out_str = out.to_str().unwrap_or("msd.dat");
    write_msd(&result, out_str)?;
    println!("MSD -> {out_str}");

    if let Some(f) = &result.fit {
        println!(
            "Fit  range [{:.2}, {:.2}]  ->  t in [{:.1}, {:.1}] fs  ({} pts)",
            f.frac_lo, f.frac_hi, f.t_lo, f.t_hi, f.n_points
        );
        println!(
            "D (total) = {:.6e} Ang^2/fs = {:.6e} cm^2/s = {:.6e} m^2/s  (R^2={:.4})",
            f.d_ang2_per_fs,
            f.d_ang2_per_fs * 0.1,
            f.d_ang2_per_fs * 1e-5,
            f.r2
        );
    }

    if args.plot {
        let png = plot_msd(&result, out_str)?;
        println!("Plot -> {png}");
        open_plot(&png);
    }
    Ok(())
}

fn run_angle(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let (group_by, slots) = resolve_selection(args)?;
    let n_given = slots.iter().filter(|x| x.is_some()).count();
    if n_given > 0 && n_given < 3 {
        return Err(anyhow!(
            "angle filter needs all three: -a END_A -b CENTRE -c END_C (by element) \
             or -x END_A -y CENTRE -z END_C (by label)"
        ));
    }

    let params = AngleParams {
        r_cut_ab: args.r_cut_ab,
        r_cut_bc: args.r_cut_bc,
        d_angle: args.d_angle,
        group_by,
    };
    let mut result = calc_angle(traj, &params)
        .ok_or_else(|| anyhow!("Angle calc failed (empty trajectory?)"))?;

    // 指定三元组时过滤输出列（key 格式 "A-Centre-C"，端原子已规范排序，两种顺序均检查）
    if let [Some(a), Some(b), Some(c)] = &slots {
        let key1 = format!("{a}-{b}-{c}");
        let key2 = format!("{c}-{b}-{a}");
        result.hist.retain(|k, _| k == &key1 || k == &key2);
        result.stats.retain(|k, _| k == &key1 || k == &key2);
        if result.hist.is_empty() {
            return Err(anyhow!(
                "triplet '{key1}' not found (check the symbols and that the centre is the middle one)"
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
