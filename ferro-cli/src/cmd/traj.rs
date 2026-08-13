//! `ferro traj` — the seven trajectory analyses that share one output pipeline.
//!
//! Grouped together because they produce the **same product**: a long- or wide-format
//! CSV with a `file` column plus an optional quick-look PNG. The old `fe-traj` /
//! `fe-corr` split tried to separate "structural" from "dynamic", but `msd` is a time
//! correlation and sat on the structural side, so the line never held.
//!
//! Each subcommand owns its own argument struct. That is the point of the subcommand
//! layout: `--dt` means one thing for `msd` and another for `vacf`, and they no longer
//! have to share a field.

use anyhow::{anyhow, Result};
use clap::{Args, Subcommand};
use ferro_analysis::{
    calc_angle, calc_gr, calc_msd, calc_rotcorr, calc_sq_from_gr, calc_vacf, calc_vanhove,
    AngleParams, AngleResult, GroupBy, GrParams, GrResult, MsdParams, MsdResult, RotCorrParams,
    RotCorrResult, SqParams, SqResult, VacfParams, VacfResult, VanHoveParams, VanHoveResult,
};

use crate::args::common::{CommonArgs, SelectArgs};
use crate::args::traj::SqWeightingCli;
use crate::batch::{self, Summary};
use crate::help;
use crate::plot::{open_plot, plot_angle, plot_gr, plot_msd, plot_sq};

#[derive(Subcommand, Debug)]
pub enum TrajCmd {
    /// Radial distribution function g(r) and coordination number CN(r)
    Gr(GrCmd),
    /// Structure factor S(q) via Fourier transform of g(r)
    Sq(SqCmd),
    /// Mean square displacement MSD(t) and self-diffusion coefficient
    Msd(MsdCmd),
    /// Bond angle distribution P(θ) for A-B-C triplets
    Angle(AngleCmd),
    /// Velocity autocorrelation function and Green-Kubo diffusion
    Vacf(VacfCmd),
    /// Rotational correlation C₂(t) for molecular bond vectors
    Rotcorr(RotcorrCmd),
    /// Van Hove self-correlation Gs(r, τ)
    Vanhove(VanhoveCmd),
}

// ─── 参数 ────────────────────────────────────────────────────────────────────

/// Options shared by g(r) and the S(q) that is derived from it.
#[derive(Args, Debug)]
pub struct GrKnobs {
    /// Min cutoff radius [Å]
    #[arg(long, default_value = "0.001")]
    pub r_min: f64,

    /// Max cutoff radius [Å]; clamped to half the smallest interplanar spacing
    #[arg(long, default_value = "10.005")]
    pub r_max: f64,

    /// Histogram bin width [Å]
    #[arg(long, default_value = "0.002")]
    pub dr: f64,
}

impl GrKnobs {
    fn params(&self, group_by: ferro_analysis::GroupBy) -> GrParams {
        GrParams { r_min: self.r_min, r_max: self.r_max, dr: self.dr, group_by }
    }
}

#[derive(Args, Debug)]
pub struct GrCmd {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub select: SelectArgs,
    #[command(flatten)]
    pub knobs: GrKnobs,
    /// Write a PNG next to the data file (needs a pair)
    #[arg(long)]
    pub plot: bool,
}

#[derive(Args, Debug)]
/// No `SelectArgs`: S(q) has no type selection.
///
/// The primary product is the pair of weighted totals; the partials are a diagnostic
/// decomposition that sums back to them (`Σ w_ij·S_ij == total`), so filtering to one
/// pair hides the very closure they exist to show. `-x/-y` used to double as the only
/// way to reach label-resolved partials — dropped with it, since a site label rarely
/// carries enough atoms for its partial to show a signal.
pub struct SqCmd {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub knobs: GrKnobs,

    /// Min q [Å⁻¹]
    #[arg(long, default_value = "0.1")]
    pub q_min: f64,
    /// Max q [Å⁻¹]
    #[arg(long, default_value = "25.0")]
    pub q_max: f64,
    /// q bin width [Å⁻¹]
    #[arg(long, default_value = "0.02")]
    pub dq: f64,
    /// Scattering weighting scheme
    #[arg(long, value_enum, default_value = "both")]
    pub weighting: SqWeightingCli,
    /// Write a PNG next to the data file
    #[arg(long)]
    pub plot: bool,
}

#[derive(Args, Debug)]
pub struct MsdCmd {
    #[command(flatten)]
    pub common: CommonArgs,

    /// Timestep between frames [fs]
    #[arg(long, default_value = "1.0")]
    pub dt: f64,
    /// Time-origin stride
    #[arg(long, default_value = "1")]
    pub shift: usize,
    /// Track only these elements, e.g. Fe,O
    #[arg(long, value_delimiter = ',')]
    pub elements: Option<Vec<String>>,
    /// Linear-fit window as trajectory fractions FMIN,FMAX (e.g. 0.3,0.8) -> D
    #[arg(long, value_delimiter = ',')]
    pub fit_range: Option<Vec<f64>>,
    /// Write a PNG next to the data file
    #[arg(long)]
    pub plot: bool,
}

#[derive(Args, Debug)]
pub struct AngleCmd {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub select: SelectArgs,

    /// Cutoff for the A-to-centre bond [Å]; A is the atom given by -a / -x
    #[arg(long, default_value = "2.3")]
    pub r_cut_ab: f64,
    /// Cutoff for the centre-to-C bond [Å]; C is the atom given by -c / -z
    #[arg(long, default_value = "2.3")]
    pub r_cut_bc: f64,
    /// Histogram lower edge [°]
    #[arg(long, default_value = "0.0")]
    pub angle_min: f64,
    /// Histogram upper edge [°]
    #[arg(long, default_value = "180.0")]
    pub angle_max: f64,
    /// Histogram bin width [°]
    #[arg(long, default_value = "0.1")]
    pub d_angle: f64,
    /// Write a PNG next to the data file (needs a triplet)
    #[arg(long)]
    pub plot: bool,
}

/// Time-axis options shared by the three correlation functions.
#[derive(Args, Debug)]
pub struct TimeKnobs {
    /// Timestep between frames [fs]
    #[arg(long, default_value = "1.0")]
    pub dt: f64,
    /// Time-origin stride
    #[arg(long, default_value = "1")]
    pub shift: usize,
    /// Lag time in frames (default: half the trajectory)
    #[arg(long)]
    pub tau: Option<usize>,
}

#[derive(Args, Debug)]
pub struct VacfCmd {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub time: TimeKnobs,
    /// Element filter, e.g. Fe,O
    #[arg(long, value_delimiter = ',')]
    pub elements: Option<Vec<String>>,
}

#[derive(Args, Debug)]
pub struct RotcorrCmd {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub time: TimeKnobs,
    /// Central atom element
    #[arg(long)]
    pub center: Option<String>,
    /// Neighbour atom element
    #[arg(long)]
    pub neighbor: Option<String>,
    /// Bond search cutoff [Å]
    #[arg(long, default_value = "1.2")]
    pub r_cut: f64,
}

#[derive(Args, Debug)]
pub struct VanhoveCmd {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub time: TimeKnobs,
    /// Element filter, e.g. Fe,O
    #[arg(long, value_delimiter = ',')]
    pub elements: Option<Vec<String>>,
    /// Max displacement [Å]
    #[arg(long, default_value = "10.0")]
    pub r_max: f64,
    /// Bin width [Å]
    #[arg(long, default_value = "0.01")]
    pub dr: f64,
}

// ─── 分派 ────────────────────────────────────────────────────────────────────

/// Runs one trajectory analysis. Returns the number of inputs that failed.
///
/// Every arm follows the same four steps — build params (validating before any file is
/// read), map over the inputs, stack the per-input tables, write and plot.
pub fn run(cmd: &TrajCmd) -> Result<usize> {
    match cmd {
        TrajCmd::Gr(c)      => run_gr(c),
        TrajCmd::Sq(c)      => run_sq(c),
        TrajCmd::Msd(c)     => run_msd(c),
        TrajCmd::Angle(c)   => run_angle(c),
        TrajCmd::Vacf(c)    => run_vacf(c),
        TrajCmd::Rotcorr(c) => run_rotcorr(c),
        TrajCmd::Vanhove(c) => run_vanhove(c),
    }
}

/// Whether a subcommand was invoked bare (no `-i`), meaning "show me the help".
pub fn wants_help(cmd: &TrajCmd) -> bool {
    let common = match cmd {
        TrajCmd::Gr(c)      => &c.common,
        TrajCmd::Sq(c)      => &c.common,
        TrajCmd::Msd(c)     => &c.common,
        TrajCmd::Angle(c)   => &c.common,
        TrajCmd::Vacf(c)    => &c.common,
        TrajCmd::Rotcorr(c) => &c.common,
        TrajCmd::Vanhove(c) => &c.common,
    };
    common.input.is_empty()
}

pub fn print_help(cmd: &TrajCmd) {
    match cmd {
        TrajCmd::Gr(_)      => help::print_gr(),
        TrajCmd::Sq(_)      => help::print_sq(),
        TrajCmd::Msd(_)     => help::print_msd(),
        TrajCmd::Angle(_)   => help::print_angle(),
        TrajCmd::Vacf(_)    => help::print_vacf(),
        TrajCmd::Rotcorr(_) => help::print_rotcorr(),
        TrajCmd::Vanhove(_) => help::print_vanhove(),
    }
}

// ─── 各分析 ──────────────────────────────────────────────────────────────────

fn run_gr(c: &GrCmd) -> Result<usize> {
    let (group_by, pair) = c.select.resolve_pair()?;
    let params = c.knobs.params(group_by);
    // label 进文件名,故在读第一个文件之前校验并建目录
    let out = c.common.out(Some(match &pair {
        Some((a, b)) => batch::file_label(&[a, b])?,
        None => batch::file_label::<&str>(&[])?,
    }));
    out.prepare()?;
    let inputs = batch::expand_inputs(&c.common.input)?;
    c.common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());

    let (results, failures) = batch::map_inputs(&inputs, |p| {
        let traj = c.common.load(p)?;
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

    let data_path = batch::write_all(
        "gr",
        "Radial Distribution Function g(r) and Coordination Number CN(r)",
        &results[0].1.meta_lines(),
        &summary.into_table(),
        tables,
        &out,
    )?;

    if c.plot {
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

fn run_sq(c: &SqCmd) -> Result<usize> {
    // S(q) 恒输出全部 partial 与两条 total:主产物是那两条 total,partial 是能加回
    // total 的诊断分解,只留一对反而看不出闭合。故这里没有类型选择,也没有 label 段
    let gr_params = c.knobs.params(GroupBy::Element);
    let out = c.common.out(None);
    out.prepare()?;
    let sq_params = SqParams {
        q_min: c.q_min,
        q_max: c.q_max,
        dq: c.dq,
        weighting: c.weighting.clone().into(),
    };
    let inputs = batch::expand_inputs(&c.common.input)?;
    c.common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());

    let (results, failures) = batch::map_inputs(&inputs, |p| {
        let traj = c.common.load(p)?;
        let gr = calc_gr(&traj, &gr_params)?;
        let sq = calc_sq_from_gr(&gr, &sq_params);
        Ok((gr, sq))
    });
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }

    let tables =
        batch::stack(&results, |(gr, sq): &(GrResult, SqResult)| Ok(sq.to_tables(gr)?))?;

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

    let data_path = batch::write_all(
        "sq",
        "Structure Factor S(q) [computed from g(r) via Fourier sine transform]",
        &results[0].1 .1.meta_lines(&results[0].1 .0),
        &summary.into_table(),
        tables,
        &out,
    )?;

    if c.plot {
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

fn run_msd(c: &MsdCmd) -> Result<usize> {
    let fit_range = match &c.fit_range {
        None => None,
        Some(v) if v.len() == 2 => Some((v[0], v[1])),
        Some(v) => return Err(anyhow!(
            "--fit-range expects exactly two comma-separated fractions, e.g. 0.3,0.8 (got {} value(s))",
            v.len()
        )),
    };
    let params = MsdParams {
        dt: c.dt,
        shift: c.shift,
        elements: c.elements.clone(),
        fit_range,
        ..MsdParams::default()
    };
    let out = c.common.out(Some(batch::set_label(c.elements.as_ref())?));
    out.prepare()?;
    let inputs = batch::expand_inputs(&c.common.input)?;
    c.common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());

    let (results, failures) = batch::map_inputs(&inputs, |p| {
        let traj = c.common.load(p)?;
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

    let data_path = batch::write_all(
        "msd",
        "Mean Squared Displacement (MSD)",
        &results[0].1.meta_lines(),
        &summary.into_table(),
        tables,
        &out,
    )?;

    if c.plot {
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

fn run_angle(c: &AngleCmd) -> Result<usize> {
    let (group_by, slots) = c.select.resolve_triplet()?;

    // 点名三元组时把端原子顺序原样传下去：--r-cut-ab 归 -a/-x 写的那一端，
    // --r-cut-bc 归 -c/-z 的那一端，不再按原子序数派发（见 AngleParams::ends）
    let ends = match (&slots[0], &slots[2]) {
        (Some(a), Some(c2)) => Some((a.clone(), c2.clone())),
        _ => None,
    };
    let params = AngleParams {
        r_cut_ab: c.r_cut_ab,
        r_cut_bc: c.r_cut_bc,
        angle_min: c.angle_min,
        angle_max: c.angle_max,
        d_angle: c.d_angle,
        group_by,
        ends,
    };
    let out = c.common.out(Some(match slots.iter().all(|s| s.is_some()) {
        true => batch::file_label(&[
            slots[0].as_ref().unwrap(),
            slots[1].as_ref().unwrap(),
            slots[2].as_ref().unwrap(),
        ])?,
        false => batch::file_label::<&str>(&[])?,
    }));
    out.prepare()?;

    let triplet_keys = slots.iter().all(|s| s.is_some()).then(|| {
        let (a, b, cc) = (
            slots[0].as_ref().unwrap(),
            slots[1].as_ref().unwrap(),
            slots[2].as_ref().unwrap(),
        );
        (format!("{a}-{b}-{cc}"), format!("{cc}-{b}-{a}"))
    });

    let inputs = batch::expand_inputs(&c.common.input)?;
    c.common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());

    let (results, failures) = batch::map_inputs(&inputs, |p| {
        let traj = c.common.load(p)?;
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
        summary.ok(batch::label_of(path), r.n_frames, r.elements.len(), &[r.hist.len() as f64]);
    }
    summary.failed(&failures);

    let data_path = batch::write_all(
        "angle",
        "Bond Angle Distribution A-B-C  (B = center)",
        &results[0].1.meta_lines(),
        &summary.into_table(),
        tables,
        &out,
    )?;

    if c.plot {
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

fn run_vacf(c: &VacfCmd) -> Result<usize> {
    let params = VacfParams {
        dt: c.time.dt,
        shift: c.time.shift,
        tau: c.time.tau,
        elements: c.elements.clone(),
    };
    let out = c.common.out(Some(batch::set_label(c.elements.as_ref())?));
    out.prepare()?;
    let inputs = batch::expand_inputs(&c.common.input)?;
    c.common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());

    let (results, failures) = batch::map_inputs(&inputs, |p| {
        let traj = c.common.load(p)?;
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
        &out,
    )?;
    Ok(failures.len())
}

fn run_rotcorr(c: &RotcorrCmd) -> Result<usize> {
    let center = c.center.clone()
        .ok_or_else(|| anyhow!("--center is required for rotcorr (run without -i to see help)"))?;
    let neighbor = c.neighbor.clone()
        .ok_or_else(|| anyhow!("--neighbor is required for rotcorr (run without -i to see help)"))?;

    // 两个参数都是必填的(上面已 bail),所以 rotcorr 恒有 label,走不到 "all"
    let out = c.common.out(Some(batch::file_label(&[&center, &neighbor])?));
    out.prepare()?;

    let params = RotCorrParams {
        center,
        neighbor,
        r_cut: c.r_cut,
        dt: c.time.dt,
        shift: c.time.shift,
        tau: c.time.tau,
    };
    let inputs = batch::expand_inputs(&c.common.input)?;
    c.common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());

    let (results, failures) = batch::map_inputs(&inputs, |p| {
        let traj = c.common.load(p)?;
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
        &out,
    )?;
    Ok(failures.len())
}

fn run_vanhove(c: &VanhoveCmd) -> Result<usize> {
    let params = VanHoveParams {
        tau: c.time.tau,
        dt: c.time.dt,
        shift: c.time.shift,
        r_max: c.r_max,
        dr: c.dr,
        elements: c.elements.clone(),
        ..VanHoveParams::default()
    };
    let out = c.common.out(Some(batch::set_label(c.elements.as_ref())?));
    out.prepare()?;
    let inputs = batch::expand_inputs(&c.common.input)?;
    c.common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());

    let (results, failures) = batch::map_inputs(&inputs, |p| {
        let traj = c.common.load(p)?;
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
        &out,
    )?;
    Ok(failures.len())
}
