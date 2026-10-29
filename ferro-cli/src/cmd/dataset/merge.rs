//! `ferro dataset merge`: combine same-composition datasets, resize sets.

use super::*;

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
    pub output: Option<PathBuf>,

    /// Create the output root without asking (required when there is no terminal)
    #[arg(long)]
    pub mkdir: bool,

    /// How frames from different sources are laid out
    #[arg(long, value_enum, default_value_t = MergeMode::Shuffle)]
    pub mode: MergeMode,

    /// Shuffle seed; by-source shuffles within each system         [default: 666]
    #[arg(long, value_name = "N")]
    pub seed: Option<u64>,

    /// Frames per output set; 0 keeps everything in one set
    #[arg(long, value_name = "N", default_value_t = DEFAULT_SET_SIZE)]
    pub set_size: usize,

    /// Force this suffix on output directories; default inherits a shared one
    // 不叫 --suffix：全仓 -s / --suffix 是产物名末段的批次标记，这里是目录名后缀（审查 E2-2）
    #[arg(long = "dir-suffix", value_name = "EXT")]
    pub dir_suffix: Option<String>,

    /// Allow writing into an existing non-empty output directory
    #[arg(long)]
    pub overwrite: bool,

    /// What to write
    #[arg(long = "type", value_enum, default_value_t = OutType::Deepmd)]
    pub out_type: OutType,

    /// Split each group train:valid:test, e.g. 8:1:1 (or 9:1 for train:test)
    #[arg(long, value_name = "A:B:C")]
    pub ratio: Option<String>,
}

pub(super) fn run_merge(args: &MergeCmd) -> Result<usize> {
    let Some(out_root) = &args.output else {
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
        // --dir-suffix 与划分后缀是同一个位置的两个主张
        if args.dir_suffix.is_some() {
            bail!("--dir-suffix and --ratio both name the output suffix; pick one");
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
    // The suffix is glued onto a directory name: `/` or `..` in it would leave `-o`
    if let Some(s) = &args.dir_suffix {
        let ok = s.strip_prefix('.').is_some_and(|rest| !rest.is_empty() && rest.chars().all(crate::batch::label_char_ok));
        if !ok {
            bail!(
                "invalid --dir-suffix {s:?}: it is appended to a directory name, so it must start \
                 with '.' followed by letters, digits, '_', '+' or '-' (e.g. .train)"
            );
        }
    }
    let roots = crate::batch::expand_dirs(&args.input)?;
    let systems = find_systems(&roots)?;
    if systems.is_empty() {
        bail!("no DeePMD system (a directory holding type.raw) found under the given paths");
    }
    if !split.is_off() {
        refuse_held_out_inputs(&systems)?;
    }
    crate::outpath::ensure_dir(out_root, args.mkdir)?;

    // 分组不看目录名 —— init.011 这类名字说明不了里面装的是什么。
    // 按逐原子的元素序列（规范序）分组，成分相同才合并
    let mut groups: BTreeMap<Vec<String>, Vec<(PathBuf, ferro_core::Trajectory)>> =
        BTreeMap::new();
    let mut failures = 0usize;
    for sys in &systems {
        match read_deepmd_npy(sys) {
            Ok((traj, warns)) => {
                crate::report::notices(sys, &warns);
                groups.entry(composition_key(&traj)).or_default().push((sys.clone(), traj));
            }
            Err(e) => {
                crate::report::error(sys.display(), format!("{e:#}; skipped"));
                failures += 1;
            }
        }
    }

    for (_, members) in groups {
        if let Err(e) = merge_group(&members, out_root, args) {
            crate::report::error(crate::report::PROGRAM, format!("{e:#}; group skipped"));
            failures += 1;
        }
    }
    if failures > 0 {
        crate::report::error(crate::report::PROGRAM, format!("{failures} failure(s)"));
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
        .map(|(p, t)| (p.clone(), t.canonically_sorted()))
        .collect();

    let name = group_name(&sorted[0].1);
    let split = merge_split(args)?;
    // 划分时后缀由 write_split 给三部分各自命名,这里不能先贴一个
    let suffix = if split.is_off() {
        match args.dir_suffix.clone() {
            Some(forced) => forced,
            None => match suffix_survey(&sorted.iter().map(|(p, _)| p.clone()).collect::<Vec<_>>()) {
                SuffixSurvey::Bare => TRAIN_SUFFIX.to_string(),
                SuffixSurvey::Shared(s) => s,
                SuffixSurvey::Mixed(seen) => {
                    let lines: Vec<String> = seen
                        .iter()
                        .map(|(p, s)| format!("  {} -> {s}", p.display()))
                        .collect();
                    bail!(
                        "this composition mixes split parts, so the merged set belongs to \
                         no one of them:\n{}\nMerge one part at a time, or name the result \
                         with --dir-suffix",
                        lines.join("\n")
                    );
                }
            },
        }
    } else {
        String::new()
    };
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
            // Never mixed across systems: each is shuffled on its own and cut into
            // its own sets, so set boundaries fall on system edges and every set
            // still comes from one condition. Without the shuffle a set would be a
            // run of consecutive MD frames. The remainder is spread within each
            // system, avoiding a lopsided 400 + 100.
            let seed = args.seed.unwrap_or(DEFAULT_SEED);
            let spans: Vec<(usize, usize)> = source_spans.iter().map(|(_, lo, hi)| (*lo, *hi)).collect();
            let all = all.subset(&shuffle_within(&spans, seed));
            println!("  shuffled within each system with seed {seed}");
            let mut bounds: Vec<(usize, usize)> = Vec::new();
            let mut record: Vec<(usize, usize, PathBuf)> = Vec::new();
            for (p, lo, hi) in &source_spans {
                for (a, b) in set_bounds(hi - lo, args.set_size) {
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
        Some(spec) => Ok(Split::parse(spec, seed)?),
    }
}

/// What the inputs say about which part of a split they are.
enum SuffixSurvey {
    /// No input carries a split suffix — the output is plain training data.
    Bare,
    /// Every input carries the same one; the output inherits it.
    Shared(String),
    /// Different parts in one group. Merging a training set into a test set is a
    /// glob written too wide far more often than it is a plan.
    Mixed(Vec<(PathBuf, String)>),
}

/// Reads the split suffixes off the input paths.
fn suffix_survey(paths: &[PathBuf]) -> SuffixSurvey {
    let seen: Vec<(PathBuf, Option<&str>)> = paths
        .iter()
        .map(|p| {
            let s = p.file_name().and_then(|n| n.to_str()).and_then(split_suffix_of);
            (p.clone(), s)
        })
        .collect();
    if seen.iter().all(|(_, s)| s.is_none()) {
        return SuffixSurvey::Bare;
    }
    let first = seen[0].1;
    if seen.iter().all(|(_, s)| *s == first) {
        // 全体一致且非 None
        return SuffixSurvey::Shared(first.unwrap_or(TRAIN_SUFFIX).to_string());
    }
    SuffixSurvey::Mixed(
        seen.into_iter()
            .map(|(p, s)| (p, s.unwrap_or("(none)").to_string()))
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::test_ids::*;

    /// merge 打乱、重排原子、切 set 之后，每帧的能量、盒子、应力、坐标与力仍属同一帧同一原子
    /// `--dir-suffix` is glued onto a directory name: one with `/` or `..` used to write
    /// outside `-o`. Refused before any input is read or any directory created.
    #[test]
    fn a_merge_suffix_that_could_leave_the_output_is_refused_up_front() {
        #[derive(clap::Parser)]
        struct Wrap {
            #[command(flatten)]
            cmd: MergeCmd,
        }
        use clap::Parser;
        let out = std::env::temp_dir().join(format!("ferro_merge_suffix_{}", std::process::id()));
        for bad in ["/../../escaped", ".a/b", "train", "."] {
            let argv = ["merge", "-i", "no/such/dir", "-o", out.to_str().unwrap(), "--mkdir", "--dir-suffix", bad];
            let err = run_merge(&Wrap::try_parse_from(argv).unwrap().cmd).unwrap_err().to_string();
            assert!(err.contains("--dir-suffix"), "{bad}: {err}");
        }
        assert!(!out.exists(), "nothing may be created for a refused suffix");
        for good in [".train", ".valid", ".run_2"] {
            let argv = ["merge", "-i", "no/such/dir", "-o", out.to_str().unwrap(), "--mkdir", "--dir-suffix", good];
            let err = run_merge(&Wrap::try_parse_from(argv).unwrap().cmd).unwrap_err().to_string();
            assert!(!err.contains("--dir-suffix"), "{good} is a valid suffix: {err}");
        }
    }

    /// Review E2-2: `-s` is the batch tag everywhere else, so `filter -s 5` no longer
    /// means 5 GPa, and merge's directory suffix is no longer called `--suffix`.
    #[test]
    fn the_old_dataset_spellings_are_refused() {
        use clap::Parser;
        #[derive(clap::Parser)]
        struct F {
            #[command(flatten)]
            cmd: FilterCmd,
        }
        #[derive(clap::Parser)]
        struct M {
            #[command(flatten)]
            cmd: MergeCmd,
        }
        assert!(F::try_parse_from(["filter", "-i", "x", "-s", "5"]).is_err());
        assert!(F::try_parse_from(["filter", "-i", "x", "--s-max", "5"]).is_ok());
        assert!(M::try_parse_from(["merge", "-i", "x", "-o", "y", "--suffix", ".train"]).is_err());
        assert!(M::try_parse_from(["merge", "-i", "x", "-o", "y", "--dir-suffix", ".train"]).is_ok());
    }

    /// `--mode by-source` shuffles each system with `--seed` and never moves a frame
    /// across systems; it used to ignore the seed and keep MD order.
    #[test]
    fn merge_by_source_shuffles_within_each_system_only() {
        #[derive(clap::Parser)]
        struct Wrap {
            #[command(flatten)]
            cmd: MergeCmd,
        }
        use clap::Parser;
        let root = std::env::temp_dir().join(format!("ferro_merge_by_source_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        ferro_io::write_deepmd_npy(&id_system(&[0, 1, 2, 3, 4, 5], &[0, 1, 2]), &root.join("in/a")).unwrap();
        ferro_io::write_deepmd_npy(&id_system(&[10, 11, 12, 13, 14], &[1, 2, 0]), &root.join("in/b")).unwrap();
        let (a, b) = (root.join("in/a"), root.join("in/b"));
        let run = |out: &str, seed: &str| {
            let out = root.join(out);
            let argv = [
                "merge", "-i", a.to_str().unwrap(), b.to_str().unwrap(), "-o", out.to_str().unwrap(),
                "--mkdir", "--mode", "by-source", "--set-size", "2", "--seed", seed,
            ];
            assert_eq!(run_merge(&Wrap::try_parse_from(argv).unwrap().cmd).unwrap(), 0);
            check_ids(&out)
        };
        let got = run("out1", "42");
        assert_eq!(got, run("out2", "42"), "same seed, different order");
        assert_ne!(got, run("out3", "43"), "the seed is ignored");
        let (from_a, from_b) = got.split_at(6);
        assert_ne!(from_a, [0, 1, 2, 3, 4, 5], "system a kept MD order");
        let mut sorted = from_a.to_vec();
        sorted.sort();
        assert_eq!(sorted, [0, 1, 2, 3, 4, 5], "a frame crossed a system edge");
        let mut sorted = from_b.to_vec();
        sorted.sort();
        assert_eq!(sorted, [10, 11, 12, 13, 14]);
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn merge_shuffle_keeps_every_label_with_its_frame_and_atom() {
        #[derive(clap::Parser)]
        struct Wrap {
            #[command(flatten)]
            cmd: MergeCmd,
        }
        use clap::Parser;
        let root = std::env::temp_dir().join(format!("ferro_merge_ids_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        // 两个来源原子顺序不同：b 的 Si 在最前
        ferro_io::write_deepmd_npy(&id_system(&[0, 1, 2, 3, 4], &[0, 1, 2]), &root.join("in/a")).unwrap();
        ferro_io::write_deepmd_npy(&id_system(&[10, 11, 12, 13], &[1, 2, 0]), &root.join("in/b")).unwrap();
        let (a, b, out) = (root.join("in/a"), root.join("in/b"), root.join("out"));
        let argv = [
            "merge", "-i", a.to_str().unwrap(), b.to_str().unwrap(), "-o", out.to_str().unwrap(),
            "--mkdir", "--set-size", "3",
        ];
        let w = Wrap::try_parse_from(argv).unwrap();
        assert_eq!(run_merge(&w.cmd).unwrap(), 0);

        let got = check_ids(&out);
        let input_order = vec![0, 1, 2, 3, 4, 10, 11, 12, 13];
        assert_ne!(got, input_order, "默认 shuffle 却没有打乱");
        let mut sorted = got.clone();
        sorted.sort();
        assert_eq!(sorted, input_order, "每帧应恰好出现一次");
        let _ = std::fs::remove_dir_all(&root);
    }

    /// Review R7-1: `-i in -i in/a` found system `a` under both roots and merged it twice.
    #[test]
    fn overlapping_roots_merge_each_system_once() {
        #[derive(clap::Parser)]
        struct Wrap {
            #[command(flatten)]
            cmd: MergeCmd,
        }
        use clap::Parser;
        let root = std::env::temp_dir().join(format!("ferro_merge_overlap_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        ferro_io::write_deepmd_npy(&id_system(&[0, 1, 2, 3, 4], &[0, 1, 2]), &root.join("in/a")).unwrap();
        ferro_io::write_deepmd_npy(&id_system(&[10, 11, 12, 13], &[1, 2, 0]), &root.join("in/b")).unwrap();
        let (inp, a, out) = (root.join("in"), root.join("in/a"), root.join("out"));
        let argv = [
            "merge", "-i", inp.to_str().unwrap(), a.to_str().unwrap(), "-o", out.to_str().unwrap(),
            "--mkdir", "--set-size", "3",
        ];
        assert_eq!(run_merge(&Wrap::try_parse_from(argv).unwrap().cmd).unwrap(), 0);
        let mut got = check_ids(&out);
        got.sort();
        assert_eq!(got, [0, 1, 2, 3, 4, 10, 11, 12, 13], "a frame was merged twice");
        let _ = std::fs::remove_dir_all(&root);
    }

    /// 后缀继承:一致就继承,没有就是训练数据,混着就说不清这份产物是哪一部分。
    #[test]
    fn the_survey_separates_shared_bare_and_mixed_inputs() {
        let same = vec![PathBuf::from("a/x.test"), PathBuf::from("b/y.test")];
        assert!(matches!(suffix_survey(&same), SuffixSurvey::Shared(s) if s == ".test"));
        let bare = vec![PathBuf::from("a/sys.001"), PathBuf::from("b/sys.002")];
        assert!(matches!(suffix_survey(&bare), SuffixSurvey::Bare));
        let mixed = vec![PathBuf::from("a/x.train"), PathBuf::from("b/y.test")];
        assert!(matches!(suffix_survey(&mixed), SuffixSurvey::Mixed(_)));
        // 一半带后缀一半不带也是混合 —— 默认的 .train 会把无名的那半也断言成训练集
        let half = vec![PathBuf::from("a/x.train"), PathBuf::from("b/sys.002")];
        assert!(matches!(suffix_survey(&half), SuffixSurvey::Mixed(_)));
    }

    /// `.db` is stripped, and it is not a split suffix.
    #[test]
    fn the_collect_suffix_is_stripped_and_is_not_a_split_suffix() {
        assert_eq!(without_collect_suffix(Path::new("o/run1.db")), PathBuf::from("o/run1"));
        assert_eq!(without_collect_suffix(Path::new("o/run1")), PathBuf::from("o/run1"));
        // 不进 SPLIT_SUFFIXES：否则「.valid 不许再划分」那条守卫会连 .db 一起拒掉
        assert_eq!(split_suffix_of("run1.db"), None);
        // 一批全 .db 判为「都没有划分后缀」→ .train，而不是混合后缀报错
        assert!(matches!(
            suffix_survey(&[PathBuf::from("a.db"), PathBuf::from("b.db")]),
            SuffixSurvey::Bare
        ));
    }
}
