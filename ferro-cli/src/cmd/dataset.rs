//! `ferro dataset` — build and maintain machine-learning training sets.
//!
//! Three steps, deliberately separate commands rather than one pipeline flag,
//! because the first one is expensive and its output is what gets backed up:
//!
//! - `collect` — AIMD output → DeePMD system directories (this file)
//! - `filter`  — quality selection on an existing dataset (not implemented yet)
//! - `merge`   — combine same-composition datasets, resize sets (not yet)
//!
//! `collect` writes one system directory per input file. Merging inputs into a
//! single system is NOT done here: a DeePMD system holds one composition, and
//! deciding which inputs belong together is `merge`'s job.

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::{Args, Subcommand};

use crate::batch::expand_inputs;
use ferro_analysis::ml::diagnostics::{
    coordination_table, count_histogram, cutoff_scan, distribution_table, pooled_coordination,
    scan_table,
};
use ferro_analysis::ml::merge::{
    composition_key, group_name, shuffle_order, sort_atoms, DEFAULT_SEED,
};
use ferro_analysis::ml::{filter_frames, first_shell_cutoff, FilterParams, FilterResult};
use ferro_core::units::{convert_pressure, PressureUnit};
use ferro_io::{
    read_cp2k_out_with_stats, read_deepmd_npy_with_warnings, write_deepmd_npy,
    write_deepmd_npy_bounds, write_deepmd_npy_sets, Cp2kOutStats,
};

#[derive(Subcommand, Debug)]
pub enum DatasetCmd {
    /// Extract AIMD output into DeePMD system directories
    Collect(CollectCmd),
    /// Drop low-quality frames from existing datasets
    Filter(FilterCmd),
    /// Combine datasets of the same composition
    Merge(MergeCmd),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum MergeMode {
    /// Concatenate everything of one composition, shuffle, then cut sets
    Shuffle,
    /// One set per input system; no mixing, no shuffling, no resizing
    BySource,
}

#[derive(Args, Debug)]
pub struct MergeCmd {
    /// System directories to combine (glob patterns allowed)
    #[arg(short, long, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Output root; one directory per composition is created under it
    #[arg(short, long, value_name = "DIR")]
    pub outdir: Option<PathBuf>,

    /// How frames from different sources are laid out         [default: shuffle]
    #[arg(long, value_enum, default_value_t = MergeMode::Shuffle)]
    pub mode: MergeMode,

    /// Shuffle seed; ignored by --mode by-source                  [default: 666]
    #[arg(long, value_name = "N")]
    pub seed: Option<u64>,

    /// Frames per output set; 0 keeps everything in one set. Only --mode
    /// shuffle uses it — by-source gives each system its own set    [400]
    #[arg(long, value_name = "N")]
    pub set_size: Option<usize>,

    /// Force this suffix on output directories; default inherits a shared one
    #[arg(long, value_name = "EXT")]
    pub suffix: Option<String>,

    /// Allow writing into an existing non-empty output directory
    #[arg(long)]
    pub overwrite: bool,
}

#[derive(Args, Debug)]
pub struct FilterCmd {
    /// System directories, or a directory holding them (searched recursively)
    #[arg(short, long, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Output root; each system is rebuilt under its path relative to -i.
    /// Omit for a read-only run that only reports.
    #[arg(short, long, value_name = "DIR")]
    pub outdir: Option<PathBuf>,

    /// Drop frames whose largest force magnitude exceeds this, eV/A; 0 = off
    #[arg(short = 'f', long, value_name = "EV_PER_A", default_value_t = 20.0)]
    pub f_max: f64,

    /// Drop frames whose largest |stress component| exceeds this, GPa; 0 = off
    #[arg(short = 's', long, value_name = "GPA", default_value_t = 10.0)]
    pub s_max: f64,

    /// First surviving frame to take (0-based, inclusive)          [default: 0]
    #[arg(long, value_name = "N")]
    pub start: Option<usize>,

    /// Last surviving frame to take (0-based, INCLUSIVE)        [default: last]
    #[arg(long, value_name = "N")]
    pub end: Option<usize>,

    /// Take every Nth surviving frame
    #[arg(long, value_name = "N")]
    pub stride: Option<usize>,

    /// Take this many surviving frames, spread evenly
    #[arg(short = 'N', long, value_name = "N", conflicts_with = "stride")]
    pub number: Option<usize>,

    /// Drop frames whose smallest O-O distance is below this, A.
    /// Bare --oo-min uses 2.0; omit the flag to switch the criterion off
    #[arg(long, num_args = 0..=1, default_missing_value = "2.0", value_name = "DMIN")]
    pub oo_min: Option<f64>,

    /// Keep only frames holding a 6-coordinated Al. Bare --al6 takes the cutoff
    /// from the Al-O RDF; give a number to set it by hand
    #[arg(long, num_args = 0..=1, default_missing_value = "auto", value_name = "RCUT")]
    pub al6: Option<String>,

    /// Frames per output set; 0 keeps everything in one set    [default: 400]
    #[arg(long, value_name = "N", default_value_t = 400)]
    pub set_size: usize,

    /// Allow writing into an existing non-empty output directory
    #[arg(long)]
    pub overwrite: bool,
}

#[derive(Args, Debug)]
pub struct CollectCmd {
    /// AIMD output files (glob patterns allowed; omit to print the full help)
    #[arg(short, long, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Directory to hold the system directories                 [default: .]
    #[arg(short, long, value_name = "DIR")]
    pub outdir: Option<PathBuf>,
}

/// True when `ferro dataset collect` was typed with no input.
pub fn wants_help(cmd: &DatasetCmd) -> bool {
    match cmd {
        DatasetCmd::Collect(c) => c.input.is_empty(),
        DatasetCmd::Filter(c) => c.input.is_empty(),
        DatasetCmd::Merge(c) => c.input.is_empty(),
    }
}

pub fn print_help(cmd: &DatasetCmd) {
    match cmd {
        DatasetCmd::Collect(_) => crate::help::print_dataset_collect(),
        DatasetCmd::Filter(_) => crate::help::print_dataset_filter(),
        DatasetCmd::Merge(_) => crate::help::print_dataset_merge(),
    }
}

/// Returns the number of inputs that failed, for the process exit code.
pub fn run(cmd: &DatasetCmd) -> Result<usize> {
    match cmd {
        DatasetCmd::Collect(c) => run_collect(c),
        DatasetCmd::Filter(c) => run_filter(c),
        DatasetCmd::Merge(c) => run_merge(c),
    }
}

fn run_collect(args: &CollectCmd) -> Result<usize> {
    let inputs = expand_inputs(&args.input)?;
    let root = args.outdir.clone().unwrap_or_else(|| PathBuf::from("."));
    // 与其余命令一致：路径问题在读第一个文件之前就暴露，而不是跑完才发现写不出去
    std::fs::create_dir_all(&root)
        .with_context(|| format!("cannot create {}", root.display()))?;

    let names = system_names(&inputs)?;
    let mut failures = 0usize;

    for (path, name) in inputs.iter().zip(&names) {
        let dir = root.join(name);
        match collect_one(path, &dir) {
            Ok(stats) => report(path, &dir, &stats),
            Err(e) => {
                eprintln!("SKIP {}: {e:#}", path.display());
                failures += 1;
            }
        }
    }

    if failures > 0 {
        eprintln!("\n{failures} of {} input(s) failed", inputs.len());
    }
    Ok(failures)
}

fn collect_one(path: &Path, dir: &Path) -> Result<Cp2kOutStats> {
    let s = path.to_string_lossy().to_string();
    let (traj, stats) = read_cp2k_out_with_stats(&s)?;
    write_deepmd_npy(&traj, dir)?;
    Ok(stats)
}

fn report(path: &Path, dir: &Path, st: &Cp2kOutStats) {
    println!("{} -> {}", path.display(), dir.display());
    println!(
        "  {} step(s), {} kept, {} dropped",
        st.n_steps, st.n_kept, st.n_dropped()
    );
    if st.n_dropped() > 0 {
        // 丢帧从不静默：5000 帧里丢掉 3000 说明 SCF 设置有问题，用户得当场知道
        println!(
            "    SCF not converged {} | incomplete block {} | composition {}",
            st.n_scf_failed, st.n_incomplete, st.n_bad_composition
        );
    }
    if st.n_restarts > 0 {
        println!("  {} restart(s) concatenated", st.n_restarts);
    }
    if st.n_layout_drift > 0 {
        println!(
            "  WARNING: {} frame(s) print their blocks at a different offset than the first;\n           extra output may be interleaved — check a few frames by hand",
            st.n_layout_drift
        );
    }
}

/// One directory name per input, disambiguated when the file stems collide.
///
/// CP2K logs are routinely all named the same thing (`total.out`) and told
/// apart by their directory, so a bare stem would make `run1/total.out` and
/// `run2/total.out` overwrite each other. Colliding stems fall back to
/// `<parent>_<stem>`; if that still collides, it is an error rather than a
/// silent overwrite.
fn system_names(inputs: &[PathBuf]) -> Result<Vec<String>> {
    let stem = |p: &Path| {
        p.file_stem()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "system".to_string())
    };
    let mut count: HashMap<String, usize> = HashMap::new();
    for p in inputs {
        *count.entry(stem(p)).or_default() += 1;
    }

    let mut names = Vec::with_capacity(inputs.len());
    let mut seen: HashMap<String, PathBuf> = HashMap::new();
    for p in inputs {
        let s = stem(p);
        let name = if count[&s] == 1 {
            s
        } else {
            let parent = p
                .parent()
                .and_then(|d| d.file_name())
                .map(|d| d.to_string_lossy().to_string())
                .unwrap_or_default();
            if parent.is_empty() { s } else { format!("{parent}_{s}") }
        };
        if let Some(first) = seen.insert(name.clone(), p.clone()) {
            bail!(
                "{} and {} would both write the system directory `{name}`; rename one or pass them separately",
                first.display(), p.display()
            );
        }
        names.push(name);
    }
    Ok(names)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unique_stems_are_used_as_is() {
        let inputs = vec![PathBuf::from("a/x.out"), PathBuf::from("b/y.out")];
        assert_eq!(system_names(&inputs).unwrap(), vec!["x", "y"]);
    }

    #[test]
    fn colliding_stems_fall_back_to_the_parent_directory() {
        let inputs = vec![PathBuf::from("run1/total.out"), PathBuf::from("run2/total.out")];
        assert_eq!(system_names(&inputs).unwrap(), vec!["run1_total", "run2_total"]);
    }

    #[test]
    fn a_shared_split_suffix_is_inherited_and_a_mixed_one_is_not() {
        let same = vec![PathBuf::from("a/x.train"), PathBuf::from("b/y.train")];
        assert_eq!(shared_suffix(&same).as_deref(), Some(".train"));
        let mixed = vec![PathBuf::from("a/x.train"), PathBuf::from("b/y.test")];
        assert_eq!(shared_suffix(&mixed), None);
        let bare = vec![PathBuf::from("a/sys.001"), PathBuf::from("b/sys.002")];
        assert_eq!(shared_suffix(&bare), None);
    }

    #[test]
    fn a_remaining_collision_is_an_error() {
        let inputs = vec![PathBuf::from("r/total.out"), PathBuf::from("x/r/total.out")];
        assert!(system_names(&inputs).is_err());
    }
}

// ── filter ───────────────────────────────────────────────────────────────────

fn run_filter(args: &FilterCmd) -> Result<usize> {
    // 参数级错误在读第一个数据集之前暴露
    if args.f_max < 0.0 || args.s_max < 0.0 {
        bail!("thresholds cannot be negative (0 switches the criterion off)");
    }
    let manual_rcut = match args.al6.as_deref() {
        None | Some("auto") => None,
        Some(v) => Some(
            v.parse::<f64>()
                .with_context(|| format!("--al6 expects a cutoff in Angstrom or nothing, got `{v}`"))?,
        ),
    };
    if let Some(r) = manual_rcut {
        if r <= 0.0 {
            bail!("--al6 cutoff must be positive");
        }
    }
    if args.oo_min.is_some_and(|v| v <= 0.0) {
        bail!("--oo-min must be positive (omit the flag to switch the criterion off)");
    }

    let params = FilterParams {
        f_max: args.f_max,
        // CLI 收 GPa（用起来顺手），内部一律 eV/Å³
        s_max: convert_pressure(args.s_max, PressureUnit::GPa, PressureUnit::EVPerAng3),
        start: args.start.unwrap_or(0),
        end: args.end,
        stride: args.stride.unwrap_or(1),
        number: args.number,
        oo_min: args.oo_min.unwrap_or(0.0),
        // 每个 system 各算各的，此处只放手动值
        al6_rcut: manual_rcut,
    };

    let roots = crate::batch::expand_dirs(&args.input)?;
    let mut jobs: Vec<(PathBuf, PathBuf)> = Vec::new(); // (system, 相对路径)
    for root in &roots {
        for sys in find_systems(root)? {
            let rel = sys.strip_prefix(root).unwrap_or(Path::new(""));
            let rel = if rel.as_os_str().is_empty() {
                PathBuf::from(root.file_name().unwrap_or(root.as_os_str()))
            } else {
                rel.to_path_buf()
            };
            jobs.push((sys, rel));
        }
    }
    if jobs.is_empty() {
        bail!("no DeePMD system (a directory holding type.raw) found under the given paths");
    }

    if let Some(out) = &args.outdir {
        std::fs::create_dir_all(out)
            .with_context(|| format!("cannot create {}", out.display()))?;
    } else {
        println!("(read-only: no -o given, nothing will be written)\n");
    }

    let mut failures = 0usize;
    let mut auto_rcuts: Vec<f64> = Vec::new();
    for (sys, rel) in &jobs {
        match filter_one(sys, rel, args, &params) {
            Ok(rcut) => auto_rcuts.extend(rcut),
            Err(e) => {
                eprintln!("SKIP {}: {e:#}", sys.display());
                failures += 1;
            }
        }
    }
    // 自动截断参与了删帧决定，不能是个看不见的数；多 system 时报均值与范围
    if auto_rcuts.len() > 1 {
        let mean = auto_rcuts.iter().sum::<f64>() / auto_rcuts.len() as f64;
        let lo = auto_rcuts.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = auto_rcuts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        println!(
            "Al-O cutoff over {} system(s): mean {mean:.3} A  (range {lo:.3} - {hi:.3})",
            auto_rcuts.len()
        );
    }
    if failures > 0 {
        eprintln!("\n{failures} of {} system(s) failed", jobs.len());
    }
    Ok(failures)
}

/// Returns the automatically derived Al-O cutoff, when one was derived.
fn filter_one(
    sys: &Path,
    rel: &Path,
    args: &FilterCmd,
    params: &FilterParams,
) -> Result<Option<f64>> {
    let (traj, warnings) = read_deepmd_npy_with_warnings(sys)?;
    for w in &warnings {
        eprintln!("WARNING: {w}");
    }

    // --al6 不带值：从这个 system 自己的 Al-O RDF 取第一壳层的外沿。
    // 逐 system 各算各的 —— 成分不同，壳层位置本来就不同
    let mut params = params.clone();
    let mut derived = None;
    if args.al6.as_deref() == Some("auto") {
        let shell = first_shell_cutoff(&traj, "Al", "O")
            .map_err(|e| anyhow::anyhow!("{e}"))?
            .context("no Al-O pair in this system, so --al6 has no cutoff to derive")?;
        println!(
            "  Al-O first shell: peak {:.2} A (g={:.1}), cutoff {:.2} A (g={:.3})",
            shell.peak_r, shell.peak_g, shell.min_r, shell.depth
        );
        if shell.depth > 0.5 {
            println!(
                "  WARNING: that minimum is shallow (g={:.2}); Al-O may have no clear shell here",
                shell.depth
            );
        }
        params.al6_rcut = Some(shell.min_r);
        derived = Some(shell.min_r);
    }

    let result = filter_frames(&traj, &params)?;
    println!("{}", sys.display());
    print_report(&result);

    let Some(out_root) = &args.outdir else {
        // 诊断只在只读模式算：它比筛选本身贵，而写出时人已经定好参数了
        print_diagnostics(&traj, &result, &params);
        println!();
        return Ok(derived);
    };
    if result.keep.is_empty() {
        bail!("every frame was dropped; nothing to write");
    }

    let dest = out_root.join(rel);
    if !args.overwrite && dest.exists() && std::fs::read_dir(&dest)?.next().is_some() {
        bail!("{} exists and is not empty (pass --overwrite)", dest.display());
    }
    // 筛过的轨迹替换原轨迹，原数据集不动
    let kept = traj.subset(&result.keep);
    write_deepmd_npy_sets(&kept, &dest, args.set_size)?;
    println!("  -> {}\n", dest.display());
    Ok(derived)
}

fn print_report(r: &FilterResult) {
    for line in r.meta_lines() {
        println!("  {line}");
    }
    for (name, table) in r.to_tables() {
        println!("  [{name}]");
        for line in table.to_comment_lines() {
            println!("    {line}");
        }
    }
}

/// Directories holding a `type.raw`, searched depth-first and not descended into.
fn find_systems(root: &Path) -> Result<Vec<PathBuf>> {
    if !root.is_dir() {
        bail!("{} is not a directory", root.display());
    }
    if root.join("type.raw").exists() {
        return Ok(vec![root.to_path_buf()]);
    }
    let mut out = Vec::new();
    let mut stack = vec![root.to_path_buf()];
    while let Some(dir) = stack.pop() {
        for entry in std::fs::read_dir(&dir)
            .with_context(|| format!("cannot list {}", dir.display()))?
            .filter_map(|e| e.ok())
        {
            let p = entry.path();
            if !p.is_dir() {
                continue;
            }
            if p.join("type.raw").exists() {
                out.push(p); // 认作 system 就不再往下走
            } else {
                stack.push(p);
            }
        }
    }
    out.sort();
    Ok(out)
}

/// The tables that answer "should I be filtering, and at what value".
///
/// Only printed in read-only mode. A selection whose outcome swings with its
/// cutoff is chosen by the cutoff rather than by the structure, and a minimum
/// distance drawn from a smooth distribution has no outliers to remove — neither
/// is visible from the funnel alone.
fn print_diagnostics(
    traj: &ferro_core::Trajectory,
    r: &FilterResult,
    params: &FilterParams,
) {
    if params.oo_min > 0.0 {
        let v: Vec<f64> = r.verdicts.iter().filter_map(|x| x.min_oo).collect();
        println!("  [min_oo distribution]");
        print_table(&distribution_table("min d(O-O) [A]", &v, 16));
    }

    let Some(rcut) = params.al6_rcut else { return };

    let n6: Vec<usize> = r.verdicts.iter().filter_map(|x| x.n_al6).collect();
    if !n6.is_empty() {
        println!("  [Al6 per frame]");
        print_table(&count_histogram("n_al6", &n6));
    }

    let mut cut = std::collections::BTreeMap::new();
    cut.insert(("Al".to_string(), "O".to_string()), rcut);
    let tp = ferro_core::TypeParams::new(cut, Default::default());
    let hist = pooled_coordination(traj, &tp, "Al");
    if !hist.is_empty() {
        println!("  [Al coordination at rcut = {rcut:.2} A]");
        print_table(&coordination_table(&hist));
    }

    // 以当前截断为中心扫一圈：陡不陡才是这张表要说的事
    let rcuts: Vec<f64> = (-3..=3).map(|k| rcut + k as f64 * 0.1).filter(|v| *v > 0.0).collect();
    let scan = cutoff_scan(traj, "Al", "O", 6, &rcuts, 200);
    println!("  [rcut sensitivity]");
    print_table(&scan_table(&scan));
}

fn print_table(t: &ferro_core::Table) {
    for line in t.to_comment_lines() {
        println!("    {line}");
    }
}

// ── merge ────────────────────────────────────────────────────────────────────

/// dpgen / dpdata split suffixes an output directory may inherit.
const SPLIT_SUFFIXES: [&str; 3] = [".train", ".test", ".valid"];

/// Frames per set when `--mode shuffle` is not told otherwise.
const DEFAULT_SET_SIZE: usize = 400;

fn run_merge(args: &MergeCmd) -> Result<usize> {
    let Some(out_root) = &args.outdir else {
        bail!("merge needs an output directory (-o DIR)");
    };
    let roots = crate::batch::expand_dirs(&args.input)?;
    let mut systems: Vec<PathBuf> = Vec::new();
    for root in &roots {
        systems.extend(find_systems(root)?);
    }
    if systems.is_empty() {
        bail!("no DeePMD system (a directory holding type.raw) found under the given paths");
    }
    std::fs::create_dir_all(out_root)
        .with_context(|| format!("cannot create {}", out_root.display()))?;

    // 分组不看目录名 —— init.011 这类名字说明不了里面装的是什么。
    // 按逐原子的元素序列（规范序）分组，成分相同才合并
    let mut groups: BTreeMap<Vec<String>, Vec<(PathBuf, ferro_core::Trajectory)>> =
        BTreeMap::new();
    let mut failures = 0usize;
    for sys in &systems {
        match read_deepmd_npy_with_warnings(sys) {
            Ok((traj, warns)) => {
                for w in warns {
                    eprintln!("WARNING: {w}");
                }
                groups.entry(composition_key(&traj)).or_default().push((sys.clone(), traj));
            }
            Err(e) => {
                eprintln!("SKIP {}: {e:#}", sys.display());
                failures += 1;
            }
        }
    }

    for (_, members) in groups {
        if let Err(e) = merge_group(&members, out_root, args) {
            eprintln!("SKIP group: {e:#}");
            failures += 1;
        }
    }
    if failures > 0 {
        eprintln!("\n{failures} failure(s)");
    }
    Ok(failures)
}

fn merge_group(
    members: &[(PathBuf, ferro_core::Trajectory)],
    out_root: &Path,
    args: &MergeCmd,
) -> Result<()> {
    // 各 system 的原子排列与 type_map 顺序都可能不同；统一到规范序，
    // 逐原子数据跟着同一个置换走。DP 对原子编号置换不变，改的是记法不是物理
    let sorted: Vec<(PathBuf, ferro_core::Trajectory)> = members
        .iter()
        .map(|(p, t)| (p.clone(), sort_atoms(t)))
        .collect();

    let name = group_name(&sorted[0].1);
    let suffix = args
        .suffix
        .clone()
        .or_else(|| shared_suffix(&sorted.iter().map(|(p, _)| p.clone()).collect::<Vec<_>>()))
        .unwrap_or_default();
    let dest = out_root.join(format!("{name}{suffix}"));

    if !args.overwrite && dest.exists() && std::fs::read_dir(&dest)?.next().is_some() {
        bail!("{} exists and is not empty (pass --overwrite)", dest.display());
    }

    let mut all = ferro_core::Trajectory::new();
    let mut source_spans: Vec<(PathBuf, usize, usize)> = Vec::new();
    for (p, t) in &sorted {
        let lo = all.frames.len();
        all.frames.extend(t.frames.iter().cloned());
        source_spans.push((p.clone(), lo, all.frames.len()));
    }
    all.metadata = sorted[0].1.metadata.clone();

    println!(
        "{name}{suffix}: {} system(s), {} frames",
        sorted.len(),
        all.n_frames()
    );
    for (p, lo, hi) in &source_spans {
        println!("  {:5} frames  {}", hi - lo, p.display());
    }

    match args.mode {
        MergeMode::Shuffle => {
            let seed = args.seed.unwrap_or(DEFAULT_SEED);
            let order = shuffle_order(all.n_frames(), seed);
            let mixed = all.subset(&order);
            write_deepmd_npy_sets(&mixed, &dest, args.set_size.unwrap_or(DEFAULT_SET_SIZE))?;
            println!("  shuffled with seed {seed} -> {}", dest.display());
        }
        MergeMode::BySource => {
            // 一个 system 一个 set：不混合、不打乱、也不重切。set 与来源
            // 一一对应，于是每个 set 就是一个条件的完整留出集
            if args.set_size.is_some() {
                println!("  note: --set-size does not apply to --mode by-source (one set per system)");
            }
            let bounds: Vec<(usize, usize)> =
                source_spans.iter().map(|(_, lo, hi)| (*lo, *hi)).collect();
            write_deepmd_npy_bounds(&all, &dest, &bounds)?;
            let mut txt = String::from("# set  frames  source\n");
            for (i, (p, lo, hi)) in source_spans.iter().enumerate() {
                txt.push_str(&format!("set.{i:03}  {}  {}\n", hi - lo, p.display()));
            }
            std::fs::write(dest.join("sets_source.txt"), txt)?;
            println!(
                "  {} set(s), one per system -> {}",
                bounds.len(),
                dest.display()
            );
        }
    }
    Ok(())
}

/// The split suffix every input shares, if they all share one.
fn shared_suffix(paths: &[PathBuf]) -> Option<String> {
    let suffix_of = |p: &PathBuf| -> Option<String> {
        let name = p.file_name()?.to_str()?;
        SPLIT_SUFFIXES
            .iter()
            .find(|s| name.ends_with(**s))
            .map(|s| s.to_string())
    };
    let first = suffix_of(&paths[0])?;
    paths
        .iter()
        .all(|p| suffix_of(p).as_deref() == Some(first.as_str()))
        .then_some(first)
}
