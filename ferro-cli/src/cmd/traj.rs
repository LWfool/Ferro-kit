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
    AngleParams, AngleResult, GroupBy, GrParams, GrResult, Legendre, MsdParams, MsdResult, RotCorrParams,
    RotCorrResult, SqParams, SqResult, VacfParams, VacfResult, VanHoveParams, VanHoveResult,
};

use crate::args::common::{CommonArgs, SelectArgs};
use crate::args::traj::{RotVectorCli, SqWeightingCli};
use crate::batch::{self, Summary};
use crate::help;
use ferro_core::Trajectory;
use std::path::PathBuf;

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
}

#[derive(Args, Debug)]
pub struct MsdCmd {
    #[command(flatten)]
    pub common: CommonArgs,

    /// Timestep between frames [fs]
    #[arg(long, default_value = "1.0")]
    pub dt: f64,
    /// Longest lag in frames (default: half the trajectory); every lag uses all time origins
    #[arg(long)]
    pub max_lag: Option<usize>,
    /// Track only these elements, e.g. Fe,O
    #[arg(long, value_delimiter = ',')]
    pub elements: Option<Vec<String>>,
    /// Linear-fit window as trajectory fractions FMIN,FMAX (e.g. 0.3,0.8) -> D
    #[arg(long, value_delimiter = ',')]
    pub fit_range: Option<Vec<f64>>,
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
}

/// Time-axis options of the correlation functions: every lag uses all time origins.
#[derive(Args, Debug)]
pub struct LagKnobs {
    /// Timestep between frames [fs]
    #[arg(long, default_value = "1.0")]
    pub dt: f64,
    /// Longest lag in frames (default: half the trajectory); every lag uses all time origins
    #[arg(long)]
    pub max_lag: Option<usize>,
}

/// Time-axis options of `vanhove`: one fixed lag, origins every `shift` frames.
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
    pub time: LagKnobs,
    /// Element filter, e.g. Fe,O
    #[arg(long, value_delimiter = ',')]
    pub elements: Option<Vec<String>>,
}

#[derive(Args, Debug)]
pub struct RotcorrCmd {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub time: LagKnobs,
    /// Central atom element
    #[arg(long)]
    pub center: Option<String>,
    /// Neighbour atom element
    #[arg(long)]
    pub neighbor: Option<String>,
    /// Bond search cutoff [Å]
    #[arg(long, default_value = "1.2")]
    pub r_cut: f64,
    /// Orientation vector: sum of the centre's bonds, or each bond fixed in frame 0
    #[arg(long, value_enum, default_value_t = RotVectorCli::Sum)]
    pub vector: RotVectorCli,
    /// Legendre order of the correlation: 1 or 2
    #[arg(long, default_value = "2", value_parser = clap::value_parser!(u8).range(1..=2))]
    pub legendre: u8,
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
/// read), map over the inputs, stack the per-input tables, write.
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

/// What [`drive`] hands back: one result per input that parsed, the inputs that did not,
/// and the prepared output location.
type Driven<T> = (Vec<(PathBuf, T)>, Vec<batch::Failure>, batch::Output);

/// The half of the pipeline all seven analyses share: prepare the output directory,
/// expand the inputs, start the thread pool, run one analysis per file, and refuse an
/// empty result set.
///
/// It stops there on purpose.  The second half — the `Summary` columns, the product name
/// and title — differs in every one of the seven, and threading those through as more
/// closures would cost more to read than the eight lines it saves.
/// `cmd/map.rs::drive` draws the line at the same place.
///
/// The output directory is built **before the first file is read**, because the label
/// goes into the file name and a bad selection should fail immediately rather than after
/// a long batch.
fn drive<T>(
    common: &CommonArgs,
    label: Option<String>,
    calc: impl Fn(&Trajectory) -> Result<T>,
) -> Result<Driven<T>> {
    let out = common.out(label);
    out.prepare()?;
    let inputs = batch::expand_inputs(&common.input)?;
    common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());

    let (results, failures) = batch::map_inputs(&inputs, |p| calc(&common.load(p)?));
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }
    Ok((results, failures, out))
}

fn run_gr(c: &GrCmd) -> Result<usize> {
    let (group_by, pair) = c.select.resolve_pair()?;
    let params = c.knobs.params(group_by);
    let label = match &pair {
        Some((a, b)) => batch::file_label(&[a, b])?,
        None => batch::file_label::<&str>(&[])?,
    };
    let (results, failures, out) =
        drive(&c.common, Some(label), |traj| Ok(calc_gr(traj, &params)?))?;

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

    batch::write_all(
        "gr",
        "Radial Distribution Function g(r) and Coordination Number CN(r)",
        &results[0].1.meta_lines(),
        &summary.into_table(),
        tables,
        &out,
    )?;
    Ok(failures.len())
}

fn run_sq(c: &SqCmd) -> Result<usize> {
    // S(q) 恒输出全部 partial 与两条 total:主产物是那两条 total,partial 是能加回
    // total 的诊断分解,只留一对反而看不出闭合。故这里没有类型选择,也没有 label 段
    let gr_params = c.knobs.params(GroupBy::Element);
    let sq_params = SqParams {
        q_min: c.q_min,
        q_max: c.q_max,
        dq: c.dq,
        weighting: c.weighting.clone().into(),
    };
    let (results, failures, out) = drive(&c.common, None, |traj| {
        let gr = calc_gr(traj, &gr_params)?;
        let sq = calc_sq_from_gr(&gr, &sq_params);
        Ok((gr, sq))
    })?;

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

    batch::write_all(
        "sq",
        "Structure Factor S(q) [computed from g(r) via Fourier sine transform]",
        &results[0].1 .1.meta_lines(&results[0].1 .0),
        &summary.into_table(),
        tables,
        &out,
    )?;
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
        max_lag: c.max_lag,
        elements: c.elements.clone(),
        fit_range,
    };
    let label = batch::set_label(c.elements.as_ref())?;
    let (results, failures, out) =
        drive(&c.common, Some(label), |traj| Ok(calc_msd(traj, &params)?))?;

    let tables = batch::stack(&results, |r: &MsdResult| Ok(r.to_tables()))?;

    // 拟合的一切（窗口换算成的 t、斜率、截距、D、误差、R²）都逐文件不同，
    // 只能进 [inputs]；头部共享区只放比例这类全批一致的参数
    let fit_cols = ["t_lo", "t_hi", "points", "slope", "intercept", "d_ang2_per_fs", "d_err", "r2"];
    // max_lag 默认随帧数取半，逐文件不同；min_origins 是最长 lag 的原点数
    let mut cols = vec!["max_lag", "min_origins"];
    if params.fit_range.is_some() {
        cols.extend(fit_cols);
    }
    let mut summary = Summary::new(&cols);
    for (path, r) in &results {
        let mut vals = vec![(r.time.len() - 1) as f64, r.min_origins as f64];
        if let Some(f) = &r.fit {
            vals.extend([
                f.t_lo, f.t_hi, f.n_points as f64, f.slope, f.intercept,
                f.d_ang2_per_fs, f.d_err, f.r2,
            ]);
        }
        summary.ok(batch::label_of(path), r.n_frames, r.n_atoms, &vals);
        summary.note("species", r.elements.join(" "));
    }
    summary.failed(&failures);

    for (path, r) in &results {
        if let Some(f) = &r.fit {
            println!(
                "{}: D = {:.6e} ± {:.1e} Ang^2/fs = {:.6e} cm^2/s = {:.6e} m^2/s  (R^2={:.4})",
                batch::label_of(path),
                f.d_ang2_per_fs,
                f.d_err,
                f.d_ang2_per_fs * 0.1,
                f.d_ang2_per_fs * 1e-5,
                f.r2
            );
        }
    }

    batch::write_all(
        "msd",
        "Mean Squared Displacement (MSD)",
        &results[0].1.meta_lines(),
        &summary.into_table(),
        tables,
        &out,
    )?;
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
    let label = match slots.iter().all(|s| s.is_some()) {
        true => batch::file_label(&[
            slots[0].as_ref().unwrap(),
            slots[1].as_ref().unwrap(),
            slots[2].as_ref().unwrap(),
        ])?,
        false => batch::file_label::<&str>(&[])?,
    };

    let triplet_keys = slots.iter().all(|s| s.is_some()).then(|| {
        let (a, b, cc) = (
            slots[0].as_ref().unwrap(),
            slots[1].as_ref().unwrap(),
            slots[2].as_ref().unwrap(),
        );
        (format!("{a}-{b}-{cc}"), format!("{cc}-{b}-{a}"))
    });

    let (results, failures, out) = drive(&c.common, Some(label), |traj| {
        let mut result = calc_angle(traj, &params)
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
    })?;

    let tables = batch::stack(&results, |r: &AngleResult| Ok(r.to_tables()))?;

    let mut summary = Summary::new(&["triplets"]);
    for (path, r) in &results {
        summary.ok(batch::label_of(path), r.n_frames, r.elements.len(), &[r.hist.len() as f64]);
    }
    summary.failed(&failures);

    batch::write_all(
        "angle",
        "Bond Angle Distribution A-B-C  (B = center)",
        &results[0].1.meta_lines(),
        &summary.into_table(),
        tables,
        &out,
    )?;
    Ok(failures.len())
}

fn run_vacf(c: &VacfCmd) -> Result<usize> {
    let params = VacfParams {
        dt: c.time.dt,
        max_lag: c.time.max_lag,
        elements: c.elements.clone(),
    };
    let label = batch::set_label(c.elements.as_ref())?;
    let (results, failures, out) =
        drive(&c.common, Some(label), |traj| Ok(calc_vacf(traj, &params)?))?;

    let tables = batch::stack(&results, |r: &VacfResult| Ok(r.to_tables()))?;
    // Green-Kubo 积分的末值逐文件不同，放进清单方便横向比（D 要看 diffusion 列走平处）
    let mut summary = Summary::new(&["max_lag", "min_origins", "diffusion_end"]);
    for (path, r) in &results {
        summary.ok(batch::label_of(path), r.n_frames, r.n_atoms, &[
            (r.time.len() - 1) as f64, r.min_origins as f64, *r.diffusion.last().unwrap(),
        ]);
        summary.note("species", r.elements.join(" "));
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
    let label = batch::file_label(&[&center, &neighbor])?;

    let params = RotCorrParams {
        center,
        neighbor,
        r_cut: c.r_cut,
        dt: c.time.dt,
        max_lag: c.time.max_lag,
        vector: c.vector.into(),
        legendre: if c.legendre == 1 { Legendre::P1 } else { Legendre::P2 },
    };
    let (results, failures, out) =
        drive(&c.common, Some(label), |traj| Ok(calc_rotcorr(traj, &params)?))?;

    let tables = batch::stack(&results, |r: &RotCorrResult| Ok(r.to_tables()))?;
    // valid_fraction：有取向向量的 (分子, 帧) 占比 —— 偏低说明 r_cut 抓不稳邻居
    // atoms = 中心原子数；units = 参与相关的单元（sum 模式同中心数，bond 模式为键数）
    let mut summary = Summary::new(&["units", "max_lag", "min_origins", "valid_fraction"]);
    for (path, r) in &results {
        summary.ok(batch::label_of(path), r.n_frames, r.n_centers, &[
            r.n_units as f64, (r.time.len() - 1) as f64, r.min_origins as f64, r.valid_fraction,
        ]);
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
    let label = batch::set_label(c.elements.as_ref())?;
    let (results, failures, out) = drive(&c.common, Some(label), |traj| {
        calc_vanhove(traj, &params)
            .ok_or_else(|| anyhow!("VanHove calc failed (trajectory too short?)"))
    })?;

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
