//! `ferro dataset` — build and maintain machine-learning training sets.
//!
//! Three steps, deliberately separate commands rather than one pipeline flag,
//! because the first one is expensive and its output is what gets backed up:
//!
//! - `collect` — AIMD output → DeePMD system directories (this file)
//! - `filter`  — quality selection on an existing dataset (not implemented yet)
//! - `merge`   — combine same-composition datasets, resize sets (not yet)
//!
//! `collect` writes one system directory per input DIRECTORY: the `.out` files
//! sitting together are the restart segments of one run, so putting them back
//! together is restoring a trajectory, not merging datasets. That is the line
//! between the two commands — `collect` reassembles the pieces of ONE run,
//! `merge` combines DIFFERENT runs of the same composition.

use std::collections::BTreeMap;
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
use ferro_core::Trajectory;
use ferro_io::{
    read_cp2k_out_with_stats, read_deepmd_npy_with_warnings, write_deepmd_npy,
    write_deepmd_npy_bounds, write_deepmd_npy_sets, write_extxyz_with, Cp2kOutStats,
    StressKey,
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

/// What `filter` / `merge` write out.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, clap::ValueEnum)]
pub enum OutType {
    /// DeePMD system directory (`type.raw` + `set.NNN/*.npy`)
    #[default]
    Deepmd,
    /// One extxyz file per system, stress carried as `virial=` (GPUMD / NEP)
    Nep,
    /// One extxyz file per system, stress carried as `stress=` (ASE convention)
    Extxyz,
}

impl OutType {
    /// The stress key, or `None` for the DeePMD path.
    fn stress_key(self) -> Option<StressKey> {
        match self {
            Self::Deepmd => None,
            Self::Nep => Some(StressKey::Virial),
            Self::Extxyz => Some(StressKey::Stress),
        }
    }
}

/// The three split parts, in output order.
const SPLIT_PARTS: [(&str, &str); 3] =
    [("train", ".train"), ("valid", ".valid"), ("test", ".test")];

/// Frame-level train/valid/test split.
///
/// Membership is drawn from a shuffled order — taking the tail as a test set
/// would hand it the end of the trajectory, which is one contiguous stretch of
/// a single state. The indices of each part are then sorted back into
/// trajectory order, so only membership is random and the output stays
/// byte-identical across runs with the same seed.
///
/// Note what this cannot fix: frames of one MD run are correlated, so a
/// frame-level test set still shares its neighbourhood with the training set
/// and reads optimistic. A split across whole systems is the honest estimate;
/// this one is the convenient one.
#[derive(Clone, Copy, Debug)]
struct Split {
    valid: f64,
    test: f64,
    seed: u64,
}

impl Split {
    /// Parses `--ratio`: `train:valid:test`, or `train:test` with two fields.
    ///
    /// The numbers are weights, not fractions — `8:1:1` and `80:10:10` are the
    /// same split — so no one has to make them sum to one. Two fields mean
    /// train:test because that is the pair NEP asks for (`train.xyz` +
    /// `test.xyz`) and the order the phrase "train/test split" already implies;
    /// every run prints the three parts by name, so a misread shows up on the
    /// first line of output rather than in a silently mislabelled dataset.
    fn parse(spec: &str, seed: u64) -> Result<Self> {
        let fields: Vec<&str> = spec.split(':').map(|f| f.trim()).collect();
        let bad = || {
            anyhow::anyhow!(
                "--ratio must look like 8:1:1 (train:valid:test) or 9:1 \
                 (train:test), got {spec:?}"
            )
        };
        if !(2..=3).contains(&fields.len()) {
            return Err(bad());
        }
        let mut w = Vec::with_capacity(3);
        for f in &fields {
            let v: f64 = f.parse().map_err(|_| bad())?;
            if !v.is_finite() || v < 0.0 {
                bail!("--ratio takes non-negative numbers, got {f:?}");
            }
            w.push(v);
        }
        let (train, valid, test) = match w.len() {
            2 => (w[0], 0.0, w[1]),
            _ => (w[0], w[1], w[2]),
        };
        let total = train + valid + test;
        if total <= 0.0 {
            bail!("--ratio is all zeros, so nothing would be written");
        }
        if valid <= 0.0 && test <= 0.0 {
            bail!("--ratio {spec:?} holds out nothing; omit --ratio instead");
        }
        Ok(Split { valid: valid / total, test: test / total, seed })
    }

    /// The split off state, used when `--ratio` is absent.
    fn off(seed: u64) -> Self {
        Split { valid: 0.0, test: 0.0, seed }
    }

    fn is_off(&self) -> bool {
        self.valid <= 0.0 && self.test <= 0.0
    }

    /// `[train, valid, test]` frame indices, each in ascending order.
    fn parts(&self, n: usize, who: &Path) -> Result<[Vec<usize>; 3]> {
        if self.is_off() {
            return Ok([(0..n).collect(), Vec::new(), Vec::new()]);
        }
        // 比例向上取到至少 1 帧:给了比例却拿到 0 帧,是静默地没有验证集
        let take = |r: f64| -> usize {
            if r <= 0.0 { 0 } else { ((n as f64 * r).round() as usize).max(1) }
        };
        let (n_valid, n_test) = (take(self.valid), take(self.test));
        if n_valid + n_test >= n {
            bail!(
                "{}: {n} frame(s) cannot give {n_valid} validation + {n_test} test \\
                 and still leave a training set",
                who.display()
            );
        }
        let order = shuffle_order(n, self.seed);
        let mut valid: Vec<usize> = order[..n_valid].to_vec();
        let mut test: Vec<usize> = order[n_valid..n_valid + n_test].to_vec();
        let mut train: Vec<usize> = order[n_valid + n_test..].to_vec();
        for v in [&mut train, &mut valid, &mut test] {
            v.sort_unstable();
        }
        Ok([train, valid, test])
    }
}

/// Appends a split suffix to the last component of a path.
fn with_suffix(p: &Path, suffix: &str) -> PathBuf {
    if suffix.is_empty() {
        return p.to_path_buf();
    }
    let name = p.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
    p.with_file_name(format!("{name}{suffix}"))
}

/// Refuses to write over an existing dataset unless told to.
fn ensure_writable(dest: &Path, overwrite: bool) -> Result<()> {
    if overwrite {
        return Ok(());
    }
    let occupied = if dest.is_dir() {
        std::fs::read_dir(dest)?.next().is_some()
    } else {
        dest.exists()
    };
    if occupied {
        bail!("{} exists and is not empty (pass --overwrite)", dest.display());
    }
    Ok(())
}

/// Writes one trajectory in the requested form; returns the path written.
///
/// `set_size` only reaches the DeePMD path — an extxyz file has no sets.
fn write_as(
    traj: &Trajectory,
    base: &Path,
    ty: OutType,
    set_size: usize,
    overwrite: bool,
) -> Result<PathBuf> {
    match ty.stress_key() {
        None => {
            ensure_writable(base, overwrite)?;
            write_deepmd_npy_sets(traj, base, set_size)?;
            Ok(base.to_path_buf())
        }
        Some(key) => {
            // 不用 with_extension:base 常带 .train 这类后缀,那会被它当扩展名换掉
            let path = PathBuf::from(format!("{}.xyz", base.display()));
            ensure_writable(&path, overwrite)?;
            if let Some(parent) = path.parent() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("cannot create {}", parent.display()))?;
            }
            write_extxyz_with(traj, &path.to_string_lossy(), key)?;
            Ok(path)
        }
    }
}

/// Writes every non-empty part of a split and reports each line.
fn write_split(
    traj: &Trajectory,
    base: &Path,
    split: &Split,
    ty: OutType,
    set_size: usize,
    overwrite: bool,
) -> Result<()> {
    let parts = split.parts(traj.n_frames(), base)?;
    for (idx, (label, suffix)) in SPLIT_PARTS.iter().enumerate() {
        if parts[idx].is_empty() {
            continue;
        }
        let suffix = if split.is_off() { "" } else { *suffix };
        let sub = traj.subset(&parts[idx]);
        let dest = write_as(&sub, &with_suffix(base, suffix), ty, set_size, overwrite)?;
        if split.is_off() {
            println!("  -> {}", dest.display());
        } else {
            println!("  {:<5} {:5} frames -> {}", label, parts[idx].len(), dest.display());
        }
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum MergeMode {
    /// Concatenate everything of one composition, shuffle, then cut sets
    Shuffle,
    /// Never mix systems; set boundaries fall on system edges
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

    /// Frames per output set; 0 keeps everything in one set    [default: 400]
    #[arg(long, value_name = "N", default_value_t = DEFAULT_SET_SIZE)]
    pub set_size: usize,

    /// Force this suffix on output directories; default inherits a shared one
    #[arg(long, value_name = "EXT")]
    pub suffix: Option<String>,

    /// Allow writing into an existing non-empty output directory
    #[arg(long)]
    pub overwrite: bool,

    /// What to write                                       [default: deepmd]
    #[arg(long = "type", value_enum, default_value_t = OutType::Deepmd)]
    pub out_type: OutType,

    /// Split each group train:valid:test, e.g. 8:1:1 (or 9:1 for train:test)
    #[arg(long, value_name = "A:B:C")]
    pub ratio: Option<String>,
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

    /// Shuffle the kept frames before writing, after every criterion has run
    #[arg(long)]
    pub shuffle: bool,

    /// Seed for --shuffle                                        [default: 666]
    #[arg(long, value_name = "N")]
    pub seed: Option<u64>,

    /// Frames per output set; 0 keeps everything in one set    [default: 400]
    #[arg(long, value_name = "N", default_value_t = 400)]
    pub set_size: usize,

    /// Allow writing into an existing non-empty output directory
    #[arg(long)]
    pub overwrite: bool,

    /// What to write                                       [default: deepmd]
    #[arg(long = "type", value_enum, default_value_t = OutType::Deepmd)]
    pub out_type: OutType,

    /// Split each system train:valid:test, e.g. 8:1:1 (or 9:1 for train:test)
    #[arg(long, value_name = "A:B:C")]
    pub ratio: Option<String>,
}

#[derive(Args, Debug)]
pub struct CollectCmd {
    /// AIMD output files (glob patterns allowed; omit to print the full help)
    #[arg(short, long, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Output root; the directory tree under -i is rebuilt inside it
    #[arg(short, long, value_name = "DIR")]
    pub outdir: Option<PathBuf>,

    /// Allow writing into an existing non-empty output directory
    #[arg(long)]
    pub overwrite: bool,
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
    // -o 是必填而不是默认 `.`：产物是一棵目录树，默认落在 cwd 会把 npy
    // 撒进正在工作的目录。merge 也是必填，filter 的缺省有「只读」这个明确语义
    let Some(root) = &args.outdir else {
        bail!("collect needs an output directory (-o DIR)");
    };
    let inputs = expand_inputs(&args.input)?;
    // 与其余命令一致：路径问题在读第一个文件之前就暴露，而不是跑完才发现写不出去
    std::fs::create_dir_all(root)
        .with_context(|| format!("cannot create {}", root.display()))?;

    let groups = group_by_directory(&inputs);
    let mut failures = 0usize;
    let mut skipped: Vec<PathBuf> = Vec::new();

    for group in &groups {
        let dest = root.join(&group.rel);
        match collect_group(group, &dest, args.overwrite, &mut skipped) {
            Ok(()) => {}
            Err(e) => {
                eprintln!("SKIP {}: {e:#}", group.dir.display());
                failures += 1;
            }
        }
    }

    if !skipped.is_empty() {
        // 合并语义下坏文件不毒化整个 system，但 system 目录看着是正常的，
        // 帧数少了却无从察觉 —— 所以这份清单要在最后再说一遍
        eprintln!("\n{} file(s) skipped and NOT in any system:", skipped.len());
        for p in &skipped {
            eprintln!("  {}", p.display());
        }
        failures += skipped.len();
    }
    if failures > 0 {
        eprintln!("\n{failures} failure(s)");
    }
    Ok(failures)
}

/// The `.out` files of one directory, which become one system.
struct Group {
    /// The directory itself, as the user wrote it — for messages.
    dir: PathBuf,
    /// Where the system goes under `-o`; empty when there is only one group.
    rel: PathBuf,
    files: Vec<PathBuf>,
}

/// One group per directory, named by the path below the shared ancestor.
///
/// The files of a directory are the restart segments of one run, so they become
/// one system rather than one each. Naming keeps the directory structure instead
/// of flattening it with separators: `-i /s/a/md/x.out /s/b/md/x.out` gives
/// `a/md` and `b/md`, and the file stem never enters the name at all.
///
/// With a single group the shared ancestor is the whole path, so `rel` is empty
/// and the system is written into `-o` itself — there is nothing to tell apart.
fn group_by_directory(inputs: &[PathBuf]) -> Vec<Group> {
    // 分组键取规范化路径，`a/x.out` 与 `./a/y.out` 才落进同一组；
    // 显示与命名仍用规范化后的路径，两者一致
    let key_of = |p: &Path| -> PathBuf {
        let dir = p.parent().unwrap_or(Path::new("."));
        std::fs::canonicalize(dir).unwrap_or_else(|_| dir.to_path_buf())
    };

    let mut order: Vec<PathBuf> = Vec::new();
    let mut by_dir: BTreeMap<PathBuf, (PathBuf, Vec<PathBuf>)> = BTreeMap::new();
    for p in inputs {
        let k = key_of(p);
        if !by_dir.contains_key(&k) {
            order.push(k.clone());
        }
        // 显示路径取用户写下的那一个（规范化后是绝对路径，刷屏且认不出）
        let as_written = p.parent().unwrap_or(Path::new(".")).to_path_buf();
        by_dir
            .entry(k)
            .or_insert_with(|| (as_written, Vec::new()))
            .1
            .push(p.clone());
    }
    order.sort();

    let ancestor = common_ancestor(&order);
    order
        .into_iter()
        .map(|key| {
            let rel = key.strip_prefix(&ancestor).unwrap_or(Path::new("")).to_path_buf();
            let (dir, files) = by_dir.remove(&key).unwrap_or_default();
            Group { dir, rel, files }
        })
        .collect()
}

/// The longest path prefix every input shares, component by component.
///
/// A shared prefix carries no distinguishing information by definition, so what
/// is left after stripping it is exactly what tells the systems apart.
fn common_ancestor(dirs: &[PathBuf]) -> PathBuf {
    let Some(first) = dirs.first() else {
        return PathBuf::new();
    };
    let mut prefix: Vec<_> = first.components().collect();
    for d in &dirs[1..] {
        let comps: Vec<_> = d.components().collect();
        let keep = prefix
            .iter()
            .zip(&comps)
            .take_while(|(a, b)| a == b)
            .count();
        prefix.truncate(keep);
    }
    prefix.iter().collect()
}

/// Reads every file of a group, concatenates them, and writes one system.
fn collect_group(
    group: &Group,
    dest: &Path,
    overwrite: bool,
    skipped: &mut Vec<PathBuf>,
) -> Result<()> {
    if !overwrite && dest.exists() && std::fs::read_dir(dest)?.next().is_some() {
        bail!("{} exists and is not empty (pass --overwrite)", dest.display());
    }

    // 先全部读进来，坏文件跳过而不毒化整个 system —— 与 reader 对坏帧的态度一致
    let mut parts: Vec<(PathBuf, Trajectory, Cp2kOutStats)> = Vec::new();
    for path in &group.files {
        match read_cp2k_out_with_stats(&path.to_string_lossy()) {
            Ok((traj, stats)) => parts.push((path.clone(), traj, stats)),
            Err(e) => {
                eprintln!("SKIP {}: {e:#}", path.display());
                skipped.push(path.clone());
            }
        }
    }
    if parts.is_empty() {
        bail!("no usable file in this directory");
    }

    // 按首个 step 号排序，文件内保持原序。全局逐帧排序看着更彻底，但重启
    // 若从 0 重新计数就会把两段真实轨迹交错洗牌，比不排序更糟；这里最坏
    // 情况退化成「按文件名拼」，不比原来差
    parts.sort_by(|a, b| {
        let ka = (a.2.steps.map(|(s, _)| s), a.0.clone());
        let kb = (b.2.steps.map(|(s, _)| s), b.0.clone());
        ka.cmp(&kb)
    });

    // 一个 system 的 type.raw 只写一次，故各文件的原子序列必须逐项相同。
    // 不一致是「把两个体系放进了一个目录」这个人的错误，不是数据的问题 ——
    // 当作坏帧丢掉会把它渲染成完全不同的一件事
    let reference = symbols_of(&parts[0].1);
    for (path, traj, _) in &parts[1..] {
        let here = symbols_of(traj);
        if here != reference {
            bail!(
                "{} and {} hold different compositions ({} vs {}); \
                 a system holds one composition, so put them in separate directories",
                parts[0].0.display(),
                path.display(),
                formula_of(&reference),
                formula_of(&here),
            );
        }
    }

    let mut all = Trajectory::new();
    all.metadata = parts[0].1.metadata.clone();
    for (_, traj, _) in &parts {
        all.frames.extend(traj.frames.iter().cloned());
    }

    report_group(dest, &parts, all.n_frames());
    write_deepmd_npy(&all, dest)?;
    Ok(())
}

fn symbols_of(traj: &Trajectory) -> Vec<String> {
    traj.frames
        .first()
        .map(|f| f.symbols().into_iter().map(|s| s.to_string()).collect())
        .unwrap_or_default()
}

/// `Al32O64Zn16` from a per-atom element sequence, for the mismatch message.
fn formula_of(symbols: &[String]) -> String {
    let mut count: BTreeMap<&str, usize> = BTreeMap::new();
    for s in symbols {
        *count.entry(s.as_str()).or_default() += 1;
    }
    count.iter().map(|(el, n)| format!("{el}{n}")).collect()
}

/// One section per system, one line per source file.
///
/// The step span is printed because the overlap of a restart is otherwise
/// invisible: frames are concatenated without de-duplication, on the grounds
/// that a restart re-runs at most a few steps and identical positions give
/// identical energies. That premise is checkable only if the spans are shown.
fn report_group(dest: &Path, parts: &[(PathBuf, Trajectory, Cp2kOutStats)], n_frames: usize) {
    println!("{}  ({} file(s), {n_frames} frames)", dest.display(), parts.len());
    // 列宽按本组实际路径算：写死的宽度装不下真实的 CP2K 目录名，span 列会错开
    let w = parts.iter().map(|(p, _, _)| p.display().to_string().len()).max().unwrap_or(0);
    for (path, _, st) in parts {
        let span = match st.steps {
            Some((a, b)) => format!("steps {a}-{b}"),
            None => "steps ?".to_string(),
        };
        println!(
            "  {:<w$}  {span:<16} {} kept, {} dropped",
            path.display().to_string(),
            st.n_kept,
            st.n_dropped()
        );
        if st.n_dropped() > 0 {
            // 丢帧从不静默：5000 帧里丢掉 3000 说明 SCF 设置有问题，用户得当场知道
            println!(
                "    SCF not converged {} | incomplete block {} | composition {}",
                st.n_scf_failed, st.n_incomplete, st.n_bad_composition
            );
        }
        if st.n_restarts > 0 {
            println!("    {} restart(s) concatenated within this file", st.n_restarts);
        }
        if st.n_layout_drift > 0 {
            println!(
                "    WARNING: {} frame(s) print their blocks at a different offset than the first;\n             extra output may be interleaved — check a few frames by hand",
                st.n_layout_drift
            );
        }
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
    let split = filter_split(args)?;
    if args.seed.is_some() && !args.shuffle && split.is_off() {
        bail!("--seed only means something with --shuffle or --ratio");
    }
    if !split.is_off() && args.outdir.is_none() {
        bail!("--ratio needs an output directory (-o DIR); a read-only run writes nothing");
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
        shuffle: args.shuffle.then(|| args.seed.unwrap_or(DEFAULT_SEED)),
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
    // 报告按表名分组堆叠，行标签用**相对路径**而不是目录名 ——
    // 嵌套结构下 a/md 与 b/md 的叶子名相同，堆起来就分不出是谁
    let mut order: Vec<String> = Vec::new();
    let mut groups: Vec<Vec<(String, ferro_core::Table)>> = Vec::new();
    let mut summary = crate::batch::Summary::new(&["frames_out"]);

    for (sys, rel) in &jobs {
        let label = rel.display().to_string();
        match filter_one(sys, rel, args, &params) {
            Ok(one) => {
                auto_rcuts.extend(one.rcut);
                summary.ok(label.clone(), one.n_input, one.n_atoms, &[one.n_kept as f64]);
                for (name, table) in one.tables {
                    match order.iter().position(|n| *n == name) {
                        Some(i) => groups[i].push((label.clone(), table)),
                        None => {
                            order.push(name);
                            groups.push(vec![(label.clone(), table)]);
                        }
                    }
                }
            }
            Err(e) => {
                eprintln!("SKIP {}: {e:#}", sys.display());
                summary.failed_one(label.clone(), format!("{e:#}"));
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

    // 报告经 write_table 落盘（分析产物的唯一出口），与 traj / net 同一条路。
    // 平铺在 -o 根下而不是塞进子目录：expand_dirs 只收 is_dir(),所以 csv 会
    // 被后续 `merge -i clean/*` 自动滤掉,而一个 report/ 子目录反倒会被收进去
    if let Some(out) = &args.outdir {
        let mut tables = Vec::with_capacity(order.len());
        for (name, parts) in order.into_iter().zip(groups) {
            let merged = ferro_core::Table::concat_union("system", parts)
                .map_err(|e| anyhow::anyhow!(e))?;
            tables.push((name, merged));
        }
        let outp = crate::batch::Output {
            dir: Some(out.clone()),
            label: None,
            suffix: None,
        };
        println!();
        crate::batch::write_all(
            "filter",
            "ferro dataset filter — frame selection report",
            &report_params(args, &params),
            &summary.into_table_named("system"),
            tables,
            &outp,
        )?;
    }

    if failures > 0 {
        eprintln!("\n{failures} of {} system(s) failed", jobs.len());
    }
    Ok(failures)
}

/// The shared parameter block at the top of every report file.
fn report_params(args: &FilterCmd, params: &FilterParams) -> Vec<String> {
    let off = |v: f64| if v > 0.0 { format!("{v}") } else { "off".to_string() };
    let mut v = vec![
        format!("f_max     = {} eV/Ang", off(args.f_max)),
        format!("s_max     = {} GPa", off(args.s_max)),
        format!("oo_min    = {} Ang", off(params.oo_min)),
        match args.al6.as_deref() {
            None => "al6       = off".to_string(),
            Some("auto") => "al6       = on (cutoff from the Al-O RDF, per system)".to_string(),
            Some(v) => format!("al6       = on (cutoff {v} Ang)"),
        },
        format!("set_size  = {}", args.set_size),
    ];
    if let Some(seed) = params.shuffle {
        v.push(format!("shuffle   = yes (seed {seed})"));
    }
    v
}

/// What one system contributed: its report tables and its frame counts.
struct FilterOne {
    /// The automatically derived Al-O cutoff, when one was derived.
    rcut: Option<f64>,
    tables: Vec<(String, ferro_core::Table)>,
    n_input: usize,
    n_kept: usize,
    n_atoms: usize,
}

fn filter_one(
    sys: &Path,
    rel: &Path,
    args: &FilterCmd,
    params: &FilterParams,
) -> Result<FilterOne> {
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

    // 诊断表恒算：实测 1110 帧 / 302 原子挂钟时间与不算时相同（rayon 跑满），
    // 而它是选阈值的依据，只在只读模式算就等于永远落不了盘
    let diagnostics = diagnostic_tables(&traj, &result, &params);
    let mut tables = result.to_tables();
    tables.extend(diagnostics.iter().cloned());

    if args.outdir.is_none() {
        // 只读模式：四张诊断表也打出来，但一个字不落盘
        for (name, table) in &diagnostics {
            println!("  [{name}]");
            print_table(table);
        }
        println!();
    }

    let one = FilterOne {
        rcut: derived,
        tables,
        n_input: result.n_input,
        n_kept: result.keep.len(),
        n_atoms: traj.frames.first().map(|f| f.n_atoms()).unwrap_or(0),
    };

    let Some(out_root) = &args.outdir else {
        return Ok(one);
    };
    if result.keep.is_empty() {
        bail!("every frame was dropped; nothing to write");
    }

    // 筛过的轨迹替换原轨迹，原数据集不动
    let kept = traj.subset(&result.keep);
    let base = out_root.join(rel);
    write_split(
        &kept,
        &base,
        &filter_split(args)?,
        args.out_type,
        args.set_size,
        args.overwrite,
    )?;
    println!();
    Ok(one)
}

/// The split `filter` was asked for; the seed is shared with `--shuffle`.
fn filter_split(args: &FilterCmd) -> Result<Split> {
    let seed = args.seed.unwrap_or(DEFAULT_SEED);
    match &args.ratio {
        None => Ok(Split::off(seed)),
        Some(spec) => Split::parse(spec, seed),
    }
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
/// A selection whose outcome swings with its cutoff is chosen by the cutoff
/// rather than by the structure, and a minimum distance drawn from a smooth
/// distribution has no outliers to remove — neither is visible from the funnel
/// alone. Only the criteria that are switched on contribute a table.
fn diagnostic_tables(
    traj: &ferro_core::Trajectory,
    r: &FilterResult,
    params: &FilterParams,
) -> Vec<(String, ferro_core::Table)> {
    let mut out = Vec::new();
    if params.oo_min > 0.0 {
        let v: Vec<f64> = r.verdicts.iter().filter_map(|x| x.min_oo).collect();
        out.push(("min_oo".to_string(), distribution_table("min d(O-O) [A]", &v, 16)));
    }

    let Some(rcut) = params.al6_rcut else { return out };

    let n6: Vec<usize> = r.verdicts.iter().filter_map(|x| x.n_al6).collect();
    if !n6.is_empty() {
        out.push(("al6".to_string(), count_histogram("n_al6", &n6)));
    }

    let mut cut = std::collections::BTreeMap::new();
    cut.insert(("Al".to_string(), "O".to_string()), rcut);
    let tp = ferro_core::TypeParams::new(cut, Default::default());
    let hist = pooled_coordination(traj, &tp, "Al");
    if !hist.is_empty() {
        let mut t = coordination_table(&hist);
        t.meta_line(format!("Al coordination at rcut = {rcut:.2} A"));
        out.push(("al_cn".to_string(), t));
    }

    // 以当前截断为中心扫一圈：陡不陡才是这张表要说的事
    let rcuts: Vec<f64> = (-3..=3).map(|k| rcut + k as f64 * 0.1).filter(|v| *v > 0.0).collect();
    let scan = cutoff_scan(traj, "Al", "O", 6, &rcuts, 200);
    out.push(("rcut_scan".to_string(), scan_table(&scan)));
    out
}

fn print_table(t: &ferro_core::Table) {
    for line in t.to_comment_lines() {
        println!("    {line}");
    }
}

// ── merge ────────────────────────────────────────────────────────────────────

/// dpgen / dpdata split suffixes an output directory may inherit.
const SPLIT_SUFFIXES: [&str; 3] = [".train", ".test", ".valid"];

/// Frames per set unless told otherwise.
const DEFAULT_SET_SIZE: usize = 400;

/// `[lo, hi)` set spans within ONE system, remainder spread evenly.
///
/// Splitting happens inside a system so no set ever straddles two of them. The
/// remainder is spread rather than left at the end: 500 frames at 400 gives
/// 250 + 250, not 400 + 100 — the lopsided pair is worse for both training
/// balance and for using a set as a validation split.
fn set_spans(n: usize, set_size: usize) -> Vec<(usize, usize)> {
    if set_size == 0 || n <= set_size {
        return vec![(0, n)];
    }
    let n_sets = n.div_ceil(set_size);
    let base = n / n_sets;
    let extra = n % n_sets;
    let mut out = Vec::with_capacity(n_sets);
    let mut lo = 0;
    for i in 0..n_sets {
        let take = base + usize::from(i < extra);
        out.push((lo, lo + take));
        lo += take;
    }
    out
}

fn run_merge(args: &MergeCmd) -> Result<usize> {
    let Some(out_root) = &args.outdir else {
        bail!("merge needs an output directory (-o DIR)");
    };
    let split = merge_split(args)?;
    if !split.is_off() {
        // by-source 的全部意义是 set 边界落在 system 边界上,每个 set 出自单一
        // 条件;帧级随机划分正好把这条打碎
        if args.mode == MergeMode::BySource {
            bail!(
                "--mode by-source keeps set boundaries on system edges, which a \
                 frame-level split would break. Use --mode shuffle, or split the \
                 systems with `ferro dataset filter`"
            );
        }
        // --suffix 与划分后缀是同一个位置的两个主张
        if args.suffix.is_some() {
            bail!("--suffix and --ratio both name the output suffix; pick one");
        }
    }
    // extxyz 没有 set 的概念,by-source 的边界无处安放 —— 与其写出一个丢了边界
    // 信息的文件,不如在读第一个 system 之前就说不支持
    if args.mode == MergeMode::BySource && args.out_type != OutType::Deepmd {
        bail!(
            "--mode by-source carries set boundaries, which an extxyz file has \
             nowhere to put (it has no sets). Use --mode shuffle, or --type deepmd"
        );
    }
    let roots = crate::batch::expand_dirs(&args.input)?;
    let mut systems: Vec<PathBuf> = Vec::new();
    for root in &roots {
        systems.extend(find_systems(root)?);
    }
    if systems.is_empty() {
        bail!("no DeePMD system (a directory holding type.raw) found under the given paths");
    }
    // 已经是 xxx.train 的输入再划一次,会得到 xxx.train.test 这种自相矛盾的名字。
    // 只看路径,故在读第一条轨迹之前就能判
    if !split.is_off() {
        if let Some(had) = systems.iter().find_map(|p| {
            let name = p.file_name()?.to_str()?;
            SPLIT_SUFFIXES.iter().find(|s| name.ends_with(**s)).map(|s| (p.clone(), *s))
        }) {
            bail!(
                "{} already carries the split suffix `{}`; splitting an already-split \
                 dataset would produce names like `X.train.test`. Merge these without \
                 --ratio, or split the unsplit sources",
                had.0.display(), had.1
            );
        }
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
    let inherited =
        shared_suffix(&sorted.iter().map(|(p, _)| p.clone()).collect::<Vec<_>>());
    let split = merge_split(args)?;
    let suffix = args.suffix.clone().or(inherited).unwrap_or_default();
    let dest = out_root.join(format!("{name}{suffix}"));

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
            println!("  shuffled with seed {seed}");
            write_split(
                &mixed,
                &dest,
                &split,
                args.out_type,
                args.set_size,
                args.overwrite,
            )?;
        }
        MergeMode::BySource => {
            // 不混合、不打乱：每个 system 独立切 set，边界落在 system 边界上，
            // 于是每个 set 仍出自单一条件；system 内部的余数均分，避免
            // 400 + 100 这种一大一小
            let mut bounds: Vec<(usize, usize)> = Vec::new();
            let mut record: Vec<(usize, usize, PathBuf)> = Vec::new();
            for (p, lo, hi) in &source_spans {
                for (a, b) in set_spans(hi - lo, args.set_size) {
                    record.push((bounds.len(), b - a, p.clone()));
                    bounds.push((lo + a, lo + b));
                }
            }
            ensure_writable(&dest, args.overwrite)?;
            write_deepmd_npy_bounds(&all, &dest, &bounds)?;
            let mut txt = String::from("# set  frames  source\n");
            for (i, n, p) in &record {
                txt.push_str(&format!("set.{i:03}  {n}  {}\n", p.display()));
            }
            std::fs::write(dest.join("sets_source.txt"), txt)?;
            println!(
                "  {} set(s), boundaries kept on system edges -> {}",
                bounds.len(),
                dest.display()
            );
        }
    }
    Ok(())
}

/// The split `merge` was asked for; the seed is shared with `--mode shuffle`.
fn merge_split(args: &MergeCmd) -> Result<Split> {
    let seed = args.seed.unwrap_or(DEFAULT_SEED);
    match &args.ratio {
        None => Ok(Split::off(seed)),
        Some(spec) => Split::parse(spec, seed),
    }
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

// ── tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn rels(inputs: &[&str]) -> Vec<String> {
        let paths: Vec<PathBuf> = inputs.iter().map(PathBuf::from).collect();
        group_by_directory(&paths)
            .iter()
            .map(|g| g.rel.display().to_string())
            .collect()
    }

    #[test]
    fn the_shared_ancestor_is_stripped_and_the_rest_kept_nested() {
        // 公共祖先 /s 剥掉，其余层级原样保留 —— 不用分隔符压平
        assert_eq!(rels(&["/s/a/md/x.out", "/s/b/md/x.out"]), vec!["a/md", "b/md"]);
        assert_eq!(rels(&["run1/total.out", "run2/total.out"]), vec!["run1", "run2"]);
    }

    #[test]
    fn a_single_directory_leaves_an_empty_name() {
        // 只有一组时公共祖先就是整条路径，产物直接写进 -o 本身
        assert_eq!(rels(&["/s/run1/a.out", "/s/run1/b.out"]), vec![""]);
    }

    #[test]
    fn files_of_one_directory_become_one_group() {
        let paths: Vec<PathBuf> = ["/s/run1/a.out", "/s/run1/b.out", "/s/run2/c.out"]
            .iter()
            .map(PathBuf::from)
            .collect();
        let groups = group_by_directory(&paths);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].files.len(), 2);
        assert_eq!(groups[1].files.len(), 1);
    }

    #[test]
    fn uneven_depth_is_kept_as_given() {
        // 输入本来就不齐，产物忠实反映；find_systems 与 merge 都是递归的
        assert_eq!(rels(&["/s/a/total.out", "/s/b/md/total.out"]), vec!["a", "b/md"]);
    }

    #[test]
    fn the_common_ancestor_stops_at_the_first_difference() {
        let dirs = vec![PathBuf::from("/s/a/md"), PathBuf::from("/s/a/opt")];
        assert_eq!(common_ancestor(&dirs), PathBuf::from("/s/a"));
        let disjoint = vec![PathBuf::from("/x/a"), PathBuf::from("/y/b")];
        assert_eq!(common_ancestor(&disjoint), PathBuf::from("/"));
    }

    #[test]
    fn set_spans_split_within_a_system_and_spread_the_remainder() {
        // 500 帧 / 400：均分成 250+250，而不是 400+100
        assert_eq!(set_spans(500, 400), vec![(0, 250), (250, 500)]);
        assert_eq!(
            set_spans(2000, 400),
            vec![(0, 400), (400, 800), (800, 1200), (1200, 1600), (1600, 2000)]
        );
        assert_eq!(set_spans(120, 400), vec![(0, 120)]);
        assert_eq!(set_spans(120, 0), vec![(0, 120)]);
        // 余数摊开：1000 / 400 -> 3 个 set，334+333+333
        assert_eq!(set_spans(1000, 400), vec![(0, 334), (334, 667), (667, 1000)]);
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
    fn a_composition_is_rendered_for_the_mismatch_message() {
        let syms: Vec<String> = ["O", "Al", "O", "Zn", "O"].iter().map(|s| s.to_string()).collect();
        assert_eq!(formula_of(&syms), "Al1O3Zn1");
    }
    fn split_of(spec: &str) -> Split {
        Split::parse(spec, DEFAULT_SEED).unwrap()
    }

    #[test]
    fn a_split_partitions_every_frame_exactly_once() {
        let [tr, va, te] = split_of("7:2:1").parts(20, Path::new("x")).unwrap();
        assert_eq!((tr.len(), va.len(), te.len()), (14, 4, 2));
        let mut all: Vec<usize> = tr.iter().chain(&va).chain(&te).copied().collect();
        all.sort_unstable();
        assert_eq!(all, (0..20).collect::<Vec<_>>());
    }

    #[test]
    fn each_part_keeps_trajectory_order() {
        // 成员是随机抽的,但每部分内部按帧序排列 —— 同 seed 下产物逐字节可复现
        let [tr, va, te] = split_of("2:1:1").parts(40, Path::new("x")).unwrap();
        for part in [&tr, &va, &te] {
            assert!(part.windows(2).all(|w| w[0] < w[1]), "{part:?}");
        }
    }

    #[test]
    fn a_split_is_not_the_tail_of_the_trajectory() {
        // 直接切尾巴的话 test 会全是最后几帧;抽样必须先打乱
        let [_, _, te] = split_of("3:0:1").parts(40, Path::new("x")).unwrap();
        assert!(te.iter().any(|&i| i < 30), "test set looks like a tail: {te:?}");
    }

    #[test]
    fn a_ratio_too_small_to_reach_one_frame_still_gets_one() {
        // 给了比例却拿到 0 帧,等于静默地没有验证集
        let [tr, va, _] = split_of("99:1:0").parts(20, Path::new("x")).unwrap();
        assert_eq!(va.len(), 1);
        assert_eq!(tr.len(), 19);
    }

    #[test]
    fn a_split_that_leaves_no_training_frames_is_an_error() {
        let e = split_of("0.001:1:1").parts(4, Path::new("sysA")).unwrap_err();
        assert!(format!("{e:#}").contains("sysA"), "{e:#}");
    }

    #[test]
    fn a_ratio_is_weights_not_fractions() {
        // 8:1:1 与 80:10:10 是同一个划分 —— 用户不必凑成和为 1
        let a = split_of("8:1:1").parts(100, Path::new("x")).unwrap();
        let b = split_of("80:10:10").parts(100, Path::new("x")).unwrap();
        assert_eq!(a[0].len(), b[0].len());
        assert_eq!((a[0].len(), a[1].len(), a[2].len()), (80, 10, 10));
    }

    #[test]
    fn two_fields_mean_train_and_test() {
        // 9:1 的 1 是 test,不是 valid;每次运行都按名字打出各部分帧数,
        // 解释错了第一行输出就看得见
        let [tr, va, te] = split_of("9:1").parts(100, Path::new("x")).unwrap();
        assert_eq!((tr.len(), va.len(), te.len()), (90, 0, 10));
    }

    #[test]
    fn a_malformed_ratio_is_rejected_with_the_expected_shape() {
        for bad in ["", "8", "8:1:1:1", "8:x:1", "8:-1:1", "0:0:0", "10:0:0"] {
            let e = Split::parse(bad, DEFAULT_SEED).unwrap_err();
            let msg = format!("{e:#}");
            assert!(
                msg.contains("--ratio"),
                "{bad:?} produced an unhelpful message: {msg}"
            );
        }
        assert!(Split::parse("8:1:1", DEFAULT_SEED).is_ok());
        assert!(Split::parse(" 8 : 1 : 1 ", DEFAULT_SEED).is_ok());
        assert!(Split::parse("0.8:0.1:0.1", DEFAULT_SEED).is_ok());
    }

    #[test]
    fn a_split_suffix_lands_on_the_last_component_only() {
        // sys.train 这种名字里的点不是扩展名,不能被 with_extension 吃掉
        assert_eq!(with_suffix(Path::new("a/b/sysA"), ".train"),
                   PathBuf::from("a/b/sysA.train"));
        assert_eq!(with_suffix(Path::new("a/b/sysA"), ""), PathBuf::from("a/b/sysA"));
    }

    #[test]
    fn the_stress_key_follows_the_output_type() {
        assert_eq!(OutType::Deepmd.stress_key(), None);
        assert_eq!(OutType::Nep.stress_key(), Some(StressKey::Virial));
        assert_eq!(OutType::Extxyz.stress_key(), Some(StressKey::Stress));
    }

}
