use anyhow::{anyhow, Result};
use clap::Parser;
use ferro::{
    args::traj::{SqWeightingCli, TrajMode},
    help::{print_fe_traj_overview, print_traj_help},
    io_dispatch::read_trajectory_tail,
    plot::{open_plot, plot_angle, plot_gr, plot_msd, plot_sq},
};
use ferro_io::LammpsUnits;
use ferro_analysis::{
    calc_angle, calc_gr, calc_msd, calc_sq_from_gr,
    AngleParams, AngleResult, GrParams, GrResult, GroupBy, MsdParams, MsdResult,
    SqParams, SqResult,
};
use ferro::batch::{self, Summary};
use std::path::PathBuf;

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

    /// Input trajectory file(s); glob patterns allowed (quote them)
    #[arg(short, long, num_args = 1.., value_name = "FILE")]
    input: Vec<PathBuf>,

    /// Output name suffix: results go to <mode>[_<table>]_<suffix>.csv
    #[arg(short, long, value_name = "SUFFIX")]
    output: Option<String>,

    /// Use only the last N frames
    #[arg(long)]
    last_n: Option<usize>,

    /// Parallel threads (default: all cores)
    #[arg(long)]
    ncore: Option<usize>,

    // ── gr / sq shared ──────────────────────────────────────────────────────

    /// [gr/sq] Min cutoff radius [Å]
    #[arg(long, default_value = "0.001")]
    r_min: f64,

    /// [gr/sq] Max cutoff radius [Å]
    #[arg(long, default_value = "10.005")]
    r_max: f64,

    /// [gr/sq] Histogram bin width [Å]
    #[arg(long, default_value = "0.002")]
    dr: f64,

    // ── sq ──────────────────────────────────────────────────────────────────

    /// [sq] Min q [Å⁻¹]
    #[arg(long, default_value = "0.1")]
    q_min: f64,

    /// [sq] Max q [Å⁻¹]
    #[arg(long, default_value = "25.0")]
    q_max: f64,

    /// [sq] q bin width [Å⁻¹]
    #[arg(long, default_value = "0.02")]
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

    /// [angle] Histogram lower edge [°]
    #[arg(long, default_value = "0.0")]
    angle_min: f64,

    /// [angle] Histogram upper edge [°]
    #[arg(long, default_value = "180.0")]
    angle_max: f64,

    /// [angle] Histogram bin width [°]
    #[arg(long, default_value = "0.1")]
    d_angle: f64,

    /// Use LAMMPS metal units for dump files (velocities Å/ps, forces eV/Å)
    #[arg(long)]
    metal_units: bool,

    /// Generate a PNG plot next to the data file
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
    if args.help || args.input.is_empty() {
        print_traj_help(&mode);
        return Ok(());
    }

    if let Some(n) = args.ncore {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
    }

    // 索引在前:模式串拼错、文件名写错要在读第一个文件之前就报出来
    let inputs = batch::expand_inputs(&args.input)?;
    println!("Inputs: {} file(s)", inputs.len());

    let failures = match mode {
        TrajMode::Gr    => run_gr(&args, &inputs)?,
        TrajMode::Sq    => run_sq(&args, &inputs)?,
        TrajMode::Msd   => run_msd(&args, &inputs)?,
        TrajMode::Angle => run_angle(&args, &inputs)?,
    };

    if failures > 0 {
        // 退出码非零:否则 shell 里 `fe-traj ... && next-step` 会把批内失败当成功
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

fn run_gr(args: &Cli, inputs: &[PathBuf]) -> Result<usize> {
    let (group_by, pair) = resolve_pair(args)?;
    let params = GrParams {
        r_min: args.r_min,
        r_max: args.r_max,
        dr: args.dr,
        group_by,
    };

    let (results, failures) = batch::map_inputs(inputs, |p| {
        let traj = load(args, p)?;
        Ok(calc_gr(&traj, &params)?)
    });
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }

    let pair_ref = pair.as_ref().map(|(a, b)| (a.as_str(), b.as_str()));
    let tables = batch::stack(&results, |r: &GrResult| Ok(r.to_tables(pair_ref)?))?;

    // r_max 逐文件 clamp 到各自盒子的最小面间距,值可能不同,所以进清单
    let mut summary = Summary::new(&["volume", "volume_std", "r_max"]);
    for (path, r) in &results {
        let atoms: usize = r.element_counts.values().sum();
        summary.ok(
            batch::label_of(path),
            r.n_frames,
            atoms,
            &[r.avg_volume, r.volume_std, r.params.r_max],
        );
        summary.note("composition", r.composition());
    }
    summary.failed(&failures);

    let meta_params = results[0].1.meta_lines();
    let data_path = batch::write_all(
        "gr",
        "Radial Distribution Function g(r) and Coordination Number CN(r)",
        &meta_params,
        &summary.into_table(),
        tables,
        args.output.as_deref(),
    )?;

    if args.plot {
        match pair_ref {
            Some((a, b)) => {
                let refs: Vec<(String, &GrResult)> =
                    results.iter().map(|(p, r)| (batch::label_of(p), r)).collect();
                let png = plot_gr(&refs, &data_path, a, b)?;
                println!("Plot    -> {png}");
                if results.len() == 1 {
                    open_plot(&png);
                }
            }
            None => println!(
                "Note  : --plot needs a pair; add -a CENTRE -b NEIGHBOUR (or -x/-y) to plot one"
            ),
        }
    }
    Ok(failures.len())
}

fn run_sq(args: &Cli, inputs: &[PathBuf]) -> Result<usize> {
    let (group_by, pair) = resolve_pair(args)?;
    let gr_params = GrParams {
        r_min: args.r_min,
        r_max: args.r_max,
        dr: args.dr,
        group_by,
    };
    let sq_params = SqParams {
        q_min: args.q_min,
        q_max: args.q_max,
        dq: args.dq,
        weighting: args.weighting.clone().into(),
    };

    let (results, failures) = batch::map_inputs(inputs, |p| {
        let traj = load(args, p)?;
        let gr = calc_gr(&traj, &gr_params)?;
        let sq = calc_sq_from_gr(&gr, &sq_params);
        Ok((gr, sq))
    });
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }

    let pair_ref = pair.as_ref().map(|(a, b)| (a.as_str(), b.as_str()));
    let tables = batch::stack(&results, |(gr, sq): &(GrResult, SqResult)| {
        Ok(sq.to_tables(gr, pair_ref)?)
    })?;

    let mut summary = Summary::new(&["volume", "volume_std", "r_max"]);
    for (path, (gr, _)) in &results {
        let atoms: usize = gr.element_counts.values().sum();
        summary.ok(
            batch::label_of(path),
            gr.n_frames,
            atoms,
            &[gr.avg_volume, gr.volume_std, gr.params.r_max],
        );
        summary.note("composition", gr.composition());
    }
    summary.failed(&failures);

    let meta_params = results[0].1 .1.meta_lines(&results[0].1 .0);
    let data_path = batch::write_all(
        "sq",
        "Structure Factor S(q) [computed from g(r) via Fourier sine transform]",
        &meta_params,
        &summary.into_table(),
        tables,
        args.output.as_deref(),
    )?;

    if args.plot {
        let refs: Vec<(String, &SqResult)> =
            results.iter().map(|(p, (_, sq))| (batch::label_of(p), sq)).collect();
        let png = plot_sq(&refs, &data_path)?;
        println!("Plot    -> {png}");
        if results.len() == 1 {
            open_plot(&png);
        }
    }
    Ok(failures.len())
}

fn run_msd(args: &Cli, inputs: &[PathBuf]) -> Result<usize> {
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

    let (results, failures) = batch::map_inputs(inputs, |p| {
        let traj = load(args, p)?;
        Ok(calc_msd(&traj, &params)?)
    });
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }

    let tables = batch::stack(&results, |r: &MsdResult| Ok(r.to_tables()))?;

    // D 与 R² 逐文件不同,是这个模式里最该横向比的两个数
    let mut summary = Summary::new(&["origins", "d_ang2_per_fs", "r2"]);
    for (path, r) in &results {
        let (d, r2) = match &r.fit {
            Some(f) => (f.d_ang2_per_fs, f.r2),
            None => (f64::NAN, f64::NAN),
        };
        summary.ok(batch::label_of(path), r.time.len(), r.n_atoms, &[r.n_origins as f64, d, r2]);
    }
    summary.failed(&failures);

    for (path, r) in &results {
        if let Some(f) = &r.fit {
            println!(
                "{}: D = {:.6e} Ang^2/fs = {:.6e} cm^2/s = {:.6e} m^2/s  (R^2={:.4})",
                batch::label_of(path),
                f.d_ang2_per_fs,
                f.d_ang2_per_fs * 0.1,
                f.d_ang2_per_fs * 1e-5,
                f.r2
            );
        }
    }

    let meta_params = results[0].1.meta_lines();
    let data_path = batch::write_all(
        "msd",
        "Mean Squared Displacement (MSD)",
        &meta_params,
        &summary.into_table(),
        tables,
        args.output.as_deref(),
    )?;

    if args.plot {
        let refs: Vec<(String, &MsdResult)> =
            results.iter().map(|(p, r)| (batch::label_of(p), r)).collect();
        let png = plot_msd(&refs, &data_path)?;
        println!("Plot    -> {png}");
        if results.len() == 1 {
            open_plot(&png);
        }
    }
    Ok(failures.len())
}

fn run_angle(args: &Cli, inputs: &[PathBuf]) -> Result<usize> {
    let (group_by, slots) = resolve_selection(args)?;
    let n_given = slots.iter().filter(|x| x.is_some()).count();
    if n_given > 0 && n_given < 3 {
        return Err(anyhow!(
            "angle filter needs all three: -a END_A -b CENTRE -c END_C (by element) \
             or -x END_A -y CENTRE -z END_C (by label)"
        ));
    }

    // 点名三元组时把端原子顺序原样传下去：--r-cut-ab 归 -a/-x 写的那一端，
    // --r-cut-bc 归 -c/-z 的那一端，不再按原子序数派发（见 AngleParams::ends）
    let ends = match (&slots[0], &slots[2]) {
        (Some(a), Some(c)) => Some((a.clone(), c.clone())),
        _ => None,
    };
    let params = AngleParams {
        r_cut_ab: args.r_cut_ab,
        r_cut_bc: args.r_cut_bc,
        angle_min: args.angle_min,
        angle_max: args.angle_max,
        d_angle: args.d_angle,
        group_by,
        ends,
    };

    let triplet_keys = slots.iter().all(|s| s.is_some()).then(|| {
        let (a, b, c) = (slots[0].as_ref().unwrap(), slots[1].as_ref().unwrap(), slots[2].as_ref().unwrap());
        (format!("{a}-{b}-{c}"), format!("{c}-{b}-{a}"))
    });

    let (results, failures) = batch::map_inputs(inputs, |p| {
        let traj = load(args, p)?;
        let mut result = calc_angle(&traj, &params)
            .ok_or_else(|| anyhow!("Angle calc failed (empty trajectory?)"))?;
        // 指定三元组时过滤输出（端原子已规范排序，两种顺序均检查）
        if let Some((key1, key2)) = &triplet_keys {
            result.hist.retain(|k, _| k == key1 || k == key2);
            result.stats.retain(|k, _| k == key1 || k == key2);
            if result.hist.is_empty() {
                return Err(anyhow!(
                    "triplet '{key1}' not found (check the symbols and that the centre is the middle one)"
                ));
            }
        }
        Ok(result)
    });
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }

    let tables = batch::stack(&results, |r: &AngleResult| Ok(r.to_tables()))?;

    let mut summary = Summary::new(&["triplets"]);
    for (path, r) in &results {
        let atoms = r.elements.len();
        summary.ok(batch::label_of(path), r.n_frames, atoms, &[r.hist.len() as f64]);
    }
    summary.failed(&failures);

    let meta_params = results[0].1.meta_lines();
    let data_path = batch::write_all(
        "angle",
        "Bond Angle Distribution A-B-C  (B = center)",
        &meta_params,
        &summary.into_table(),
        tables,
        args.output.as_deref(),
    )?;

    if args.plot {
        match &triplet_keys {
            Some((key1, key2)) => {
                let refs: Vec<(String, &AngleResult)> =
                    results.iter().map(|(p, r)| (batch::label_of(p), r)).collect();
                // 端原子可能被规范排序反转,两个 key 都试
                let key = if refs.iter().any(|(_, r)| r.hist.contains_key(key1)) { key1 } else { key2 };
                let png = plot_angle(&refs, &data_path, key)?;
                println!("Plot    -> {png}");
                if results.len() == 1 {
                    open_plot(&png);
                }
            }
            None => println!(
                "Note  : --plot needs a triplet; add -a END_A -b CENTRE -c END_C (or -x/-y/-z)"
            ),
        }
    }
    Ok(failures.len())
}
