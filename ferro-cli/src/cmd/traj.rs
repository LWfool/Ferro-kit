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
    calc_angle, calc_bondlife, calc_gr, calc_msd, calc_rotcorr, calc_sq_from_gr, calc_vacf, calc_vanhove,
    AngleParams, AngleResult, BondLifeParams, BondLifeResult, GroupBy, GrParams, GrResult, Legendre, MsdParams, MsdResult, RotCorrParams,
    RotCorrResult, SqParams, SqResult, VacfParams, VacfResult, VanHoveParams, VanHoveResult,
};

use crate::args::common::{CommonArgs, SelectArgs};
use crate::args::traj::{RotVectorCli, SqWeightingCli};
use crate::batch::{self, Summary};
use crate::help;
use ferro_core::Trajectory;

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
    /// Bond lifetimes (intermittent / continuous) and bond formation / breaking events
    Bondlife(BondlifeCmd),
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

    /// Time between stored frames [fs] = MD timestep × dump interval; required
    #[arg(long)]
    pub dt: Option<f64>,
    /// Longest lag in frames (default: half the trajectory); every lag uses all time origins
    #[arg(long)]
    pub max_lag: Option<usize>,
    /// Track only these elements, e.g. Fe,O
    #[arg(long, value_delimiter = ',')]
    pub elements: Option<Vec<String>>,
    /// Linear-fit window as fractions of the lag axis FMIN,FMAX (e.g. 0.3,0.8) -> D
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
    /// Time between stored frames [fs] = MD timestep × dump interval; required
    #[arg(long)]
    pub dt: Option<f64>,
    /// Longest lag in frames (default: half the trajectory); every lag uses all time origins
    #[arg(long)]
    pub max_lag: Option<usize>,
}

/// Time-axis options of `vanhove`: one fixed lag, origins every `shift` frames.
#[derive(Args, Debug)]
pub struct TimeKnobs {
    /// Time between stored frames [fs] = MD timestep × dump interval; required
    #[arg(long)]
    pub dt: Option<f64>,
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
pub struct BondlifeCmd {
    #[command(flatten)]
    pub common: CommonArgs,
    #[command(flatten)]
    pub time: LagKnobs,
    /// Central atom element
    #[arg(long)]
    pub center: Option<String>,
    /// Neighbour atom element (may equal --center)
    #[arg(long)]
    pub neighbor: Option<String>,
    /// A free pair bonds at r <= r-bond [Å] (first g(r) minimum)
    #[arg(long)]
    pub r_bond: Option<f64>,
    /// A bonded pair stays bonded while r <= r-break [Å] (default: --r-bond)
    #[arg(long)]
    pub r_break: Option<f64>,
    /// Breaks of at most this many frames are filled (continuous function and events)
    #[arg(long, default_value = "0")]
    pub intermittency: usize,
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
        TrajCmd::Bondlife(c) => run_bondlife(c),
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
        TrajCmd::Bondlife(c) => &c.common,
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
        TrajCmd::Bondlife(_) => help::print_bondlife(),
    }
}

// ─── 各分析 ──────────────────────────────────────────────────────────────────

/// What [`drive`] hands back: one result per input that parsed, the inputs that did not,
/// and the prepared output location.
type Driven<T> = (Vec<(batch::Input, T)>, Vec<batch::Failure>, batch::Output);

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
    // 参数（含 --atom-style）在建目录、读第一个文件之前查完
    let inputs = common.inputs()?;
    let out = common.out(label);
    out.prepare()?;
    common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());

    let (results, failures) = batch::map_inputs(&inputs, |inp| calc(&common.load(&inp.path)?));
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }
    Ok((results, failures, out))
}

/// 时间相关分析的前提：帧等间隔、不重复。轨迹带步号（LAMMPS dump、CP2K、OUTCAR）时
/// 逐帧核对；没有步号的格式无从查起，只能信 `--dt`。
///
/// 重启拼接的 dump 常把重启点写两次（同一个 TIMESTEP），或中途改过 dump 间隔；
/// 全原点平均把每一对帧都当成 `lag × dt` 相隔，这两种情形都会静默算错。
fn check_frame_spacing(traj: &Trajectory) -> Result<()> {
    let steps: Option<Vec<i64>> = traj.frames.iter().map(|f| f.step).collect();
    let Some(steps) = steps else { return Ok(()) };
    let Some(first_gap) = steps.get(1).map(|s1| s1 - steps[0]) else { return Ok(()) };
    for (k, w) in steps.windows(2).enumerate() {
        let gap = w[1] - w[0];
        if gap <= 0 {
            return Err(anyhow!(
                "frames {k} and {} have steps {} and {}: duplicated or out-of-order frames \
                 (a restart written twice?); remove them before a time-correlation analysis",
                k + 1, w[0], w[1]
            ));
        }
        if gap != first_gap {
            return Err(anyhow!(
                "frames are not evenly spaced: step gap {first_gap} at the start, {gap} between \
                 frames {k} and {} (steps {} → {}); --dt assumes one constant spacing",
                k + 1, w[0], w[1]
            ));
        }
    }
    Ok(())
}

fn run_gr(c: &GrCmd) -> Result<usize> {
    let (group_by, pair) = c.select.resolve_pair()?;
    let params = c.knobs.params(group_by);
    // 取值范围在建目录、读第一个文件之前查完（审查 M3）
    params.validate()?;
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
    for (input, r) in &results {
        let atoms: usize = r.element_counts.values().sum();
        summary.ok(
            input.label.clone(),
            r.n_frames,
            atoms,
            &[r.avg_volume, r.volume_std, r.r_max_used],
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
    // 取值范围在建目录、读第一个文件之前查完（审查 M3）
    gr_params.validate()?;
    sq_params.validate()?;
    let (results, failures, out) = drive(&c.common, None, |traj| {
        let gr = calc_gr(traj, &gr_params)?;
        let sq = calc_sq_from_gr(&gr, &sq_params)?;
        Ok((gr, sq))
    })?;

    let tables =
        batch::stack(&results, |(gr, sq): &(GrResult, SqResult)| Ok(sq.to_tables(gr)?))?;

    let mut summary = Summary::new(&["volume", "volume_std", "r_max", "q_trunc"]);
    for (input, (gr, _)) in &results {
        let atoms: usize = gr.element_counts.values().sum();
        // g(r) 截断在 r_max，q ≲ 2π/r_max 的 S(q) 由截断振荡主导（审查 B-3）。
        // r_max 逐输入截到最小镜像上界，故此值逐输入、进 [inputs] 而不进共享头部
        let q_trunc = 2.0 * std::f64::consts::PI / gr.r_max_used;
        summary.ok(
            input.label.clone(),
            gr.n_frames,
            atoms,
            &[gr.avg_volume, gr.volume_std, gr.r_max_used, q_trunc],
        );
        summary.note("composition", gr.composition());
        if c.q_min < q_trunc {
            eprintln!(
                "        warning: {}: S(q) below q_trunc = 2π/r_max = {:.3} Å⁻¹ (r_max = {:.3} Å) \
                 is dominated by the truncation of g(r); --q-min is {}",
                input.label, q_trunc, gr.r_max_used, c.q_min
            );
        }
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
        dt: c.dt.ok_or_else(|| anyhow!("--dt is required for msd (run without -i to see help)"))?,
        max_lag: c.max_lag,
        elements: c.elements.clone(),
        fit_range,
    };
    // 取值范围在建目录、读第一个文件之前查完（审查 M3）
    params.validate()?;
    let label = batch::set_label(c.elements.as_ref())?;
    let (results, failures, out) =
        drive(&c.common, Some(label), |traj| { check_frame_spacing(traj)?; Ok(calc_msd(traj, &params)?) })?;

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
    for (input, r) in &results {
        let mut vals = vec![(r.time.len() - 1) as f64, r.min_origins as f64];
        if let Some(f) = &r.fit {
            vals.extend([
                f.t_lo, f.t_hi, f.n_points as f64, f.slope, f.intercept,
                f.d_ang2_per_fs, f.d_err, f.r2,
            ]);
        }
        summary.ok(input.label.clone(), r.n_frames, r.n_atoms, &vals);
        summary.note("species", r.elements.join(" "));
    }
    summary.failed(&failures);

    for (input, r) in &results {
        if let Some(f) = &r.fit {
            println!(
                "{}: D = {:.6e} ± {:.1e} Ang^2/fs = {:.6e} cm^2/s = {:.6e} m^2/s  (R^2={:.4})",
                input.label.clone(),
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
    // 取值范围在建目录、读第一个文件之前查完（审查 M3）
    params.validate()?;
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
        let mut result = calc_angle(traj, &params)?;
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
    for (input, r) in &results {
        summary.ok(input.label.clone(), r.n_frames, r.n_atoms, &[r.hist.len() as f64]);
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
        dt: c.time.dt.ok_or_else(|| anyhow!("--dt is required for vacf (run without -i to see help)"))?,
        max_lag: c.time.max_lag,
        elements: c.elements.clone(),
    };
    // 取值范围在建目录、读第一个文件之前查完（审查 M3）
    params.validate()?;
    let label = batch::set_label(c.elements.as_ref())?;
    let (results, failures, out) =
        drive(&c.common, Some(label), |traj| { check_frame_spacing(traj)?; Ok(calc_vacf(traj, &params)?) })?;

    let tables = batch::stack(&results, |r: &VacfResult| Ok(r.to_tables()))?;
    // Green-Kubo 积分的末值逐文件不同，放进清单方便横向比（D 要看 diffusion 列走平处）
    let mut summary = Summary::new(&["max_lag", "min_origins", "diffusion_end"]);
    for (input, r) in &results {
        summary.ok(input.label.clone(), r.n_frames, r.n_atoms, &[
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
        dt: c.time.dt.ok_or_else(|| anyhow!("--dt is required for rotcorr (run without -i to see help)"))?,
        max_lag: c.time.max_lag,
        vector: c.vector.into(),
        legendre: if c.legendre == 1 { Legendre::P1 } else { Legendre::P2 },
    };
    // 取值范围在建目录、读第一个文件之前查完（审查 M3）
    params.validate()?;
    let (results, failures, out) =
        drive(&c.common, Some(label), |traj| { check_frame_spacing(traj)?; Ok(calc_rotcorr(traj, &params)?) })?;

    let tables = batch::stack(&results, |r: &RotCorrResult| Ok(r.to_tables()))?;
    // valid_fraction：有取向向量的 (分子, 帧) 占比 —— 偏低说明 r_cut 抓不稳邻居
    // atoms = 中心原子数；units = 参与相关的单元（sum 模式同中心数，bond 模式为键数）
    let mut summary = Summary::new(&["units", "max_lag", "min_origins", "valid_fraction"]);
    for (input, r) in &results {
        summary.ok(input.label.clone(), r.n_frames, r.n_centers, &[
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

fn run_bondlife(c: &BondlifeCmd) -> Result<usize> {
    let center = c.center.clone()
        .ok_or_else(|| anyhow!("--center is required for bondlife (run without -i to see help)"))?;
    let neighbor = c.neighbor.clone()
        .ok_or_else(|| anyhow!("--neighbor is required for bondlife (run without -i to see help)"))?;
    let r_bond = c.r_bond
        .ok_or_else(|| anyhow!("--r-bond is required: take the first minimum of the {center}-{neighbor} g(r)"))?;
    let label = batch::file_label(&[&center, &neighbor])?;
    let params = BondLifeParams {
        center, neighbor, r_bond,
        r_break: c.r_break,
        intermittency: c.intermittency,
        max_lag: c.time.max_lag,
        dt: c.time.dt.ok_or_else(|| anyhow!("--dt is required for bondlife (run without -i to see help)"))?,
    };
    // 取值范围（含 r-break >= r-bond）在建目录、读第一个文件之前查完
    params.validate()?;
    let (results, failures, out) =
        drive(&c.common, Some(label), |traj| { check_frame_spacing(traj)?; Ok(calc_bondlife(traj, &params)?) })?;

    let tables = batch::stack(&results, |r: &BondLifeResult| Ok(r.to_tables()))?;
    let mut summary = Summary::new(&[
        "candidates", "mean_bonds", "max_lag", "min_origins",
        "tau_int_integral", "tau_int_1e", "tau_cont_integral", "tau_cont_1e",
    ]);
    for (input, r) in &results {
        let (ii, i1, ci, c1) = r.taus();
        summary.ok(input.label.clone(), r.n_frames, r.n_centers, &[
            r.n_candidates as f64, r.mean_bonds(), (r.time.len() - 1) as f64, r.min_origins as f64,
            ii, i1, ci, c1,
        ]);
    }
    summary.failed(&failures);

    batch::write_all(
        "bondlife",
        "Bond lifetimes: intermittent C_I(t), continuous S_C(t); per-frame bond events",
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
        dt: c.time.dt.ok_or_else(|| anyhow!("--dt is required for vanhove (run without -i to see help)"))?,
        shift: c.time.shift,
        r_max: c.r_max,
        dr: c.dr,
        elements: c.elements.clone(),
        ..VanHoveParams::default()
    };
    // 取值范围在建目录、读第一个文件之前查完（审查 M3）
    params.validate()?;
    let label = batch::set_label(c.elements.as_ref())?;
    let (results, failures, out) = drive(&c.common, Some(label), |traj| {
        check_frame_spacing(traj)?;
        Ok(calc_vanhove(traj, &params)?)
    })?;

    let tables = batch::stack(&results, |r: &VanHoveResult| Ok(r.to_tables()))?;
    // tau 逐文件相同,但 time = tau*dt 与 origins 值得横向看一眼；
    // outside_fraction 是落在 [r_min, r_max) 外、没画进 p_r 的位移比例，逐文件不同
    let mut summary = Summary::new(&["tau_frames", "time_fs", "origins", "outside_fraction"]);
    for (input, r) in &results {
        summary.ok(
            input.label.clone(),
            r.n_frames,
            r.n_atoms,
            &[r.tau_frames as f64, r.time, r.n_origins as f64, r.outside_fraction],
        );
        if r.outside_fraction > 0.01 {
            eprintln!(
                "        warning: {}: {:.1}% of displacements exceed r_max = {} Å and are not in p_r; \
                 raise --r-max",
                input.label, r.outside_fraction * 100.0, c.r_max
            );
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[derive(Parser)]
    struct Wrap {
        #[command(subcommand)]
        cmd: TrajCmd,
    }

    /// 坏参数必须在读文件、建目录之前报错（审查 M3）。输入文件不存在：若先读文件，
    /// 报的会是「找不到输入」而不是参数错误；`-o` 目录不存在且给了 --mkdir：若先建
    fn frames_with_steps(steps: &[i64]) -> Trajectory {
        let mut t = Trajectory::new();
        for &s in steps {
            let mut f = ferro_core::Frame::new();
            f.step = Some(s);
            t.add_frame(f);
        }
        t
    }

    #[test]
    fn test_frame_spacing() {
        assert!(check_frame_spacing(&frames_with_steps(&[0, 100, 200, 300])).is_ok());
        assert!(check_frame_spacing(&frames_with_steps(&[0])).is_ok(), "单帧无间隔可查");
        assert!(check_frame_spacing(&Trajectory::from_frame(ferro_core::Frame::new())).is_ok(), "无步号的格式只能信 --dt");
        let dup = format!("{:#}", check_frame_spacing(&frames_with_steps(&[0, 100, 100, 200])).unwrap_err());
        assert!(dup.contains("frames 1 and 2") && dup.contains("duplicated"), "重启点写两次应点名帧，实际 {dup}");
        let uneven = format!("{:#}", check_frame_spacing(&frames_with_steps(&[0, 100, 200, 400])).unwrap_err());
        assert!(uneven.contains("not evenly spaced") && uneven.contains("frames 2 and 3"), "改过 dump 间隔应点名帧，实际 {uneven}");
    }

    /// 目录，测试结束时它就在
    #[test]
    fn test_bad_values_fail_before_any_file_or_directory() {
        let out = std::env::temp_dir().join("ferro_m3_traj_must_not_exist");
        let _ = std::fs::remove_dir_all(&out);
        for (args, want) in [
            ("gr --dr 0", "dr must be"),
            ("gr --dr nan", "dr must be"),
            ("gr --r-min 5 --r-max 3", "r-min (5) must be < r-max (3)"),
            ("gr --r-min=-1", "r-min must be a finite number >= 0"),
            ("sq --q-max inf", "q-max must be a finite number"),
            ("sq --q-min 0", "q-min must be a finite number > 0"),
            ("sq --dq 0", "dq must be"),
            ("sq --q-min 5 --q-max 1", "q-min (5) must be < q-max (1)"),
            ("msd --dt 0", "dt must be"),
            ("msd --dt 1 --max-lag 0", "max-lag must be >= 1"),
            ("angle -a O -b P -c O --d-angle 0", "d-angle must be"),
            ("angle -a O -b P -c O --angle-max 200", "angle-max must be <= 180"),
            ("vacf --dt=-1", "dt must be"),
            ("rotcorr --dt 1 --center P --neighbor O --r-cut 0", "r-cut must be"),
            ("bondlife --dt 1 --center P --neighbor O --r-bond 2 --r-break 1", "r-break (1) must be >= r-bond (2)"),
            ("vanhove --dt 1 --dr 0", "dr must be"),
            ("vanhove --dt 1 --shift 0", "shift must be >= 1"),
            ("msd", "--dt is required"),
            ("vacf", "--dt is required"),
            ("rotcorr --center P --neighbor O", "--dt is required"),
            ("bondlife --center P --neighbor O --r-bond 2", "--dt is required"),
            ("vanhove", "--dt is required"),
        ] {
            let mut argv = vec!["traj"];
            argv.extend(args.split_whitespace());
            argv.extend(["-i", "no_such_input.lammpstrj", "-o", out.to_str().unwrap(), "--mkdir"]);
            let w = Wrap::try_parse_from(&argv).unwrap_or_else(|e| panic!("{args}: {e}"));
            let err = format!("{:#}", run(&w.cmd).expect_err(args));
            assert!(err.contains(want), "{args}：应报「{want}」，实际 {err}");
            assert!(!out.exists(), "{args}：参数错误时不该建 -o 目录");
        }
    }

    /// 不带 -i 显示帮助页的前提是 clap 先放行：任何必填参数都会让 clap 在
    /// `wants_help` 之前报错退出，帮助页就永远到不了。必填项一律在 run 里查
    #[test]
    fn test_bare_subcommand_reaches_help_page() {
        for sub in ["gr", "sq", "msd", "angle", "vacf", "rotcorr", "vanhove", "bondlife"] {
            let w = Wrap::try_parse_from(["traj", sub])
                .unwrap_or_else(|e| panic!("traj {sub}：不带参数应能解析以显示帮助页，实际 {e}"));
            assert!(wants_help(&w.cmd), "traj {sub}：不带 -i 应显示帮助页");
        }
    }
}
