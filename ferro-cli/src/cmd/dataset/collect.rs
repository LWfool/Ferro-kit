//! `ferro dataset collect`: AIMD output -> DeePMD system directories.

use super::*;

#[derive(Args, Debug)]
pub struct CollectCmd {
    /// AIMD output files (glob patterns allowed; PJM job logs x0-<job>.<id>.out are
    /// skipped; omit to print the full help)
    #[arg(short, long, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Output root; default is beside each input, as <AIMD dir>.db
    #[arg(short, long, value_name = "DIR")]
    pub output: Option<PathBuf>,

    /// Create the output root without asking (required when there is no terminal)
    #[arg(long)]
    pub mkdir: bool,

    /// Allow writing into an existing non-empty output directory
    #[arg(long)]
    pub overwrite: bool,

    /// What to write
    #[arg(long = "type", value_enum, default_value_t = CollectType::Deepmd)]
    pub out_type: CollectType,

    /// Read the inputs as this format instead of trusting their banner
    #[arg(long, value_enum, value_name = "FMT")]
    pub format: Option<InputFormat>,
}

/// The AIMD layout a file is to be read as, when the user names it.
///
/// The five values mirror [`AimdFormat`] one for one. VASP is two values rather
/// than one `vasp/md`: OUTCAR and vasprun.xml record the same frames but decide
/// convergence differently, so they drop a different number of them.
#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
pub enum InputFormat {
    #[value(name = "cp2k/md")]
    Cp2kMd,
    #[value(name = "cp2k/sp")]
    Cp2kSp,
    #[value(name = "vasp/outcar")]
    VaspOutcar,
    #[value(name = "vasp/xml")]
    VaspXml,
    #[value(name = "qe/pw")]
    QePw,
}

impl From<InputFormat> for AimdFormat {
    fn from(f: InputFormat) -> Self {
        match f {
            InputFormat::Cp2kMd => AimdFormat::Cp2kMd,
            InputFormat::Cp2kSp => AimdFormat::Cp2kSp,
            InputFormat::VaspOutcar => AimdFormat::VaspOutcar,
            InputFormat::VaspXml => AimdFormat::VaspXml,
            InputFormat::QePw => AimdFormat::QeOut,
        }
    }
}

pub(super) fn run_collect(args: &CollectCmd) -> Result<usize> {
    let inspect = args.out_type == CollectType::Inspect;
    // -o 的含义是「数据集写到哪」，而 --type inspect 一个 npy 都不写。静默忽略
    // 会让用户回头去那个目录找数据集 —— 参数级错误在读第一个文件之前失败
    if inspect && args.output.is_some() {
        bail!(
            "--type inspect writes no dataset; its three files always go to <AIMD dir>/{}. Drop -o, or drop --type inspect",
            super::inspect::DIR_NAME
        );
    }
    let inputs = drop_pjm_logs(expand_inputs(&args.input)?)?;
    // 与其余命令一致：路径问题在读第一个文件之前就暴露，而不是跑完才发现写不出去
    if let Some(root) = &args.output {
        crate::outpath::ensure_dir(root, args.mkdir)?;
    }

    let fmt = args.format.map(AimdFormat::from);

    let groups = group_by_directory(&inputs);
    let mut failures = 0usize;
    let mut skipped: Vec<PathBuf> = Vec::new();

    for group in &groups {
        let dest = if inspect {
            group.dir.join(super::inspect::DIR_NAME)
        } else {
            collect_dest(args.output.as_deref(), group)
        };
        match collect_group(group, &dest, args.overwrite, inspect, fmt, &mut skipped) {
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

/// Removes PJM job logs (`x0-<job>.<id>.out`) before grouping.
///
/// `-i '*.out'` matches them beside the real output, and they are empty or hold
/// scheduler messages — each one used to count as a failure and turn the exit
/// code to 1, which hid the files that were really broken.
fn drop_pjm_logs(inputs: Vec<PathBuf>) -> Result<Vec<PathBuf>> {
    let total = inputs.len();
    let (logs, kept): (Vec<PathBuf>, Vec<PathBuf>) = inputs.into_iter().partition(|p| {
        p.file_name()
            .and_then(|n| n.to_str())
            .is_some_and(is_pjm_log)
    });
    if logs.is_empty() {
        return Ok(kept);
    }
    // 在分组之前剔除：只剩 PJM 日志的目录不成组，也就不会报 no usable file
    if kept.is_empty() {
        bail!("all {total} input(s) are PJM job logs (x<N>-<job>.<id>.out); nothing to collect");
    }
    eprintln!(
        "NOTE: excluded {} PJM job log(s) matching x<N>-<job>.<id>.out",
        logs.len()
    );
    Ok(kept)
}

/// `x<digits>-<jobname>.<jobid digits>.out`, the stdout/stderr name PJM gives.
///
/// The jobid part is what keeps a user's own `x1-run.out` out of the match.
fn is_pjm_log(name: &str) -> bool {
    let all_digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let Some(rest) = name.strip_prefix('x') else {
        return false;
    };
    let Some((n, rest)) = rest.split_once('-') else {
        return false;
    };
    let Some(rest) = rest.strip_suffix(".out") else {
        return false;
    };
    // jobname 自身可含 `.` 与 `-`，故 jobid 取最后一个 `.` 之后
    let Some((job, id)) = rest.rsplit_once('.') else {
        return false;
    };
    all_digits(n) && !job.is_empty() && all_digits(id)
}

/// Where one group's system directory goes.
///
/// Without `-o` the system lands **beside** the AIMD directory it came from —
/// `/data/md/*.out` gives `/data/md.db`. That is the shape dpgen and the
/// reference script use, and it keeps a collected dataset next to the run that
/// produced it. A default of `.` was rejected earlier for scattering npy files
/// into the working directory; a default that follows the input cannot do that,
/// which is why `-o` is no longer required.
///
/// With `-o` the directory tree below the shared ancestor is rebuilt inside it;
/// a single group, whose path below the ancestor is empty, goes in under its
/// directory's own name, so `-o` is the folder the systems go into whatever
/// their number. Either way the name carries `.db`, so a downstream `filter -i`
/// sees one shape.
fn collect_dest(root: Option<&Path>, group: &Group) -> PathBuf {
    let with_suffix = |p: &Path| -> PathBuf {
        let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("system");
        // 已经带着 .db 的不叠加（重跑同一条命令时 -i 可能指到上一轮的产物旁）
        if name.ends_with(COLLECT_SUFFIX) {
            p.to_path_buf()
        } else {
            p.with_file_name(format!("{name}{COLLECT_SUFFIX}"))
        }
    };
    // The path as the user wrote it, not canonicalized: the absolute path fills
    // the screen and is hard to recognise. `with_file_name` swaps the last
    // component, so beside it is exactly where the system lands. Only `.` / `..`
    // (`-i *.out` in the current directory) have no last name and fall back.
    let dir = if group.dir.file_name().is_some() {
        group.dir.clone()
    } else {
        group.dir.canonicalize().unwrap_or_else(|_| group.dir.clone())
    };
    match root {
        Some(r) if group.rel.as_os_str().is_empty() => {
            with_suffix(&r.join(dir.file_name().unwrap_or("system".as_ref())))
        }
        Some(r) => with_suffix(&r.join(&group.rel)),
        None => with_suffix(&dir),
    }
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
/// and `collect_dest` names the system after its directory.
fn group_by_directory(inputs: &[PathBuf]) -> Vec<Group> {
    // 分组键取规范化路径，`a/x.out` 与 `./a/y.out` 才落进同一组；
    // 显示与命名仍用规范化后的路径，两者一致
    // A bare `a.out` has parent `Some("")`, which canonicalize rejects: read it as `.`.
    let dir_of = |p: &Path| -> PathBuf {
        match p.parent() {
            Some(d) if !d.as_os_str().is_empty() => d.to_path_buf(),
            _ => PathBuf::from("."),
        }
    };
    let key_of = |p: &Path| -> PathBuf {
        let dir = dir_of(p);
        std::fs::canonicalize(&dir).unwrap_or(dir)
    };

    let mut order: Vec<PathBuf> = Vec::new();
    let mut by_dir: BTreeMap<PathBuf, (PathBuf, Vec<PathBuf>)> = BTreeMap::new();
    for p in inputs {
        let k = key_of(p);
        if !by_dir.contains_key(&k) {
            order.push(k.clone());
        }
        // 显示路径取用户写下的那一个（规范化后是绝对路径，刷屏且认不出）
        let as_written = dir_of(p);
        by_dir
            .entry(k)
            .or_insert_with(|| (as_written, Vec::new()))
            .1
            .push(p.clone());
    }
    order.sort();

    let mut ancestor = common_ancestor(&order);
    // A group sitting at the ancestor itself would get an empty name and land in
    // `<-o>.db`, beside the `<-o>/` holding its siblings. One level up, every
    // group keeps a name and all of them nest under `-o`.
    if order.len() > 1 && order.contains(&ancestor) {
        if let Some(up) = ancestor.parent() {
            ancestor = up.to_path_buf();
        }
    }
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
    inspect: bool,
    format: Option<AimdFormat>,
    skipped: &mut Vec<PathBuf>,
) -> Result<()> {
    // inspect 那条路自己查，因为它的目录是 ferro_inspect/ 而不是 system 目录。
    // 复用 ensure_writable 而不是手抄（审查 E2-11）：手抄的那份漏了「dest 是普通文件」，
    // 报出来的是 read_dir 的 io 错误
    if !inspect {
        ensure_writable(dest, overwrite)?;
    }

    // 先全部读进来，坏文件跳过而不毒化整个 system —— 与 reader 对坏帧的态度一致。
    //
    // Each file goes into canonical order as soon as it is read (a permutation per
    // frame): two files written independently have their own atom orders, and one
    // file does not promise a constant order from frame to frame either
    let mut parts: Vec<(PathBuf, Trajectory, AimdStats)> = Vec::new();
    // 嗅探结果单独留一份：--format 一给，stats.format 就恒等于它，下面那条
    // 「同目录两种格式」的守卫会失去比较的对象。嗅探本来就要跑（它出告警），
    // 守卫白拿
    let mut sniffed: Vec<(PathBuf, AimdFormat)> = Vec::new();
    for path in &group.files {
        let banner = ferro_io::sniff(path);
        if let Ok(f) = &banner {
            sniffed.push((path.clone(), *f));
        }
        let fmt = match (format, &banner) {
            // 用户说了算，但手滑写错与有意覆盖长得一样，所以说一声
            (Some(asked), Ok(seen)) if *seen != asked => {
                eprintln!(
                    "NOTE: {} looks like {} but --format says {}; reading as {}",
                    path.display(),
                    seen.name(),
                    asked.name(),
                    asked.name()
                );
                asked
            }
            (Some(asked), _) => asked,
            (None, _) => match banner {
                Ok(f) => f,
                Err(e) => {
                    // The banner check is ferro-io's, and knows nothing of this flag
                    eprintln!("SKIP {}: {e:#} (--format reads it as a given program)", path.display());
                    skipped.push(path.clone());
                    continue;
                }
            },
        };
        match read_aimd_as(path, fmt) {
            Ok((traj, stats)) => parts.push((path.clone(), traj.canonically_sorted(), stats)),
            Err(e) => {
                eprintln!("SKIP {}: {e:#}", path.display());
                skipped.push(path.clone());
            }
        }
    }
    if parts.is_empty() {
        bail!("no usable file in this directory");
    }

    // A real VASP run directory holds OUTCAR and vasprun.xml side by side, recording the
    // same frames. collect's rule is "files of one directory = segments of one run",
    // and following it would double the frame count silently: same composition, both
    // files readable, no other symptom
    if let Some((first, fmt0)) = sniffed.first() {
        if let Some((path, other)) = sniffed.iter().find(|(_, f)| f != fmt0) {
            bail!(
                "{} is {} but {} is {}. A run directory holds both, and they record \
                 the same frames — concatenating them would double the dataset. \
                 Narrow -i to one of the two",
                first.display(), fmt0.name(), path.display(), other.name(),
            );
        }
    }

    // Sort by first step, check one composition, concatenate (ferro-analysis)
    let (all, sources) = stitch_segments(&mut parts, |(path, traj, stats)| {
        (path.as_path(), traj, stats.steps.map(|(s, _)| s))
    })?;

    report_group(dest, &parts, all.n_frames());
    report_kind_labels(&all);
    if inspect {
        let name = group
            .dir
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| "run".to_string());
        for p in super::inspect::write_all(&all, parts[0].2.format, dest, &name, &sources, overwrite)? {
            println!("  -> {}", p.display());
        }
    } else {
        write_deepmd_npy(&all, dest)?;
    }
    Ok(())
}

/// One section per system, one line per source file.
///
/// The step span is printed because the overlap of a restart is otherwise
/// invisible: frames are concatenated without de-duplication, on the grounds
/// that a restart re-runs at most a few steps and identical positions give
/// identical energies. That premise is checkable only if the spans are shown.
fn report_group(dest: &Path, parts: &[(PathBuf, Trajectory, AimdStats)], n_frames: usize) {
    let fmt = parts[0].2.format;
    println!(
        "{}  ({} file(s), {n_frames} frames, {})",
        dest.display(), parts.len(), fmt.name()
    );
    // 列宽按本组实际路径算：写死的宽度装不下真实的 CP2K 目录名，span 列会错开
    let w = parts.iter().map(|(p, _, _)| p.display().to_string().len()).max().unwrap_or(0);
    for (path, _, st) in parts {
        let span = match (st.steps, st.format) {
            (Some((a, b)), _) => format!("steps {a}-{b}"),
            // 单点没有步号,"steps ?" 会被读成"步号丢了",而其实是没有步号这回事
            (None, AimdFormat::Cp2kSp | AimdFormat::QeOut) => format!("{} single point(s)", st.n_steps),
            (None, _) => "steps ?".to_string(),
        };
        println!(
            "  {:<w$}  {span:<16} {} kept, {} dropped",
            path.display().to_string(),
            st.n_kept,
            st.n_dropped()
        );
        if st.n_dropped() > 0 {
            // 丢帧从不静默：5000 帧里丢掉 3000 说明 SCF 设置有问题，用户得当场知道。
            // 判据随格式变（OUTCAR 读 VASP 自己的结论，vasprun 只能数 SCF 步），
            // 所以把它一并打出来 —— 否则同一次运行换个来源、丢帧数不同会没人说得清
            println!("    convergence rule: {}", st.format.convergence_rule());
            println!(
                "    SCF not converged {} | incomplete block {} | no cell line {} | composition {}",
                st.n_scf_failed, st.n_incomplete, st.n_no_cell, st.n_bad_composition
            );
            if st.n_no_cell > 0 {
                println!("    (no cell line: print &MOTION/&PRINT/&CELL to __STD_OUT__, see `ferro doc dataset collect`)");
            }
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
        // 版本提示由 reader 判定、这里打印:库自己 print 就没人关得掉,
        // 而 collect 一次要读几十个文件
        if let Some(note) = &st.version_note {
            println!("    NOTE: {note}");
        }
    }
}

/// Prints the kind names CP2K carried in, once per system.
///
/// `type_map.raw` is built from the ELEMENT, so two kinds sharing one element
/// (`Fe1`/`Fe2`, a spin guess) become one training type. That is the right
/// default for a DeePMD model, but it throws away a distinction the person
/// deliberately made in the CP2K input, so it is said out loud rather than
/// left to be discovered in the type map.
fn report_kind_labels(traj: &Trajectory) {
    let Some(first) = traj.first() else { return };
    let mut seen: BTreeMap<(&str, &str), usize> = BTreeMap::new();
    for a in first.atoms() {
        if let Some(l) = a.label {
            *seen.entry((l, a.element)).or_default() += 1;
        }
    }
    if seen.is_empty() {
        return;
    }
    let pairs: Vec<String> = seen
        .iter()
        .map(|((label, element), n)| format!("{label} -> {element} ({n})"))
        .collect();
    println!("    kind names kept as labels, type_map uses the element: {}", pairs.join("  "));
}

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

    /// One CP2K single point of the same two atoms, written in a chosen order.
    ///
    /// `si_first` flips only the order the atoms appear in — the physics is
    /// identical, which is what makes the assertion below sharp.
    fn sp_text(si_first: bool) -> String {
        let head = " CP2K| version string:                                       CP2K version 2025.2\n\
                    \x20GLOBAL| Run type                                                   ENERGY_FORCE\n\
                    \x20CELL| Volume [angstrom^3]:                                           125.000000\n\
                    \x20CELL| Vector a [angstrom]:       5.000     0.000     0.000   |a| =     5.000000\n\
                    \x20CELL| Vector b [angstrom]:       0.000     5.000     0.000   |b| =     5.000000\n\
                    \x20CELL| Vector c [angstrom]:       0.000     0.000     5.000   |c| =     5.000000\n\
                    \x20CELL| Periodicity                                                           XYZ\n\n\
                    \x20MODULE QUICKSTEP: ATOMIC COORDINATES IN ANGSTROM\n\n\
                    \x20\x20 Atom Kind Element         X             Y             Z       Z(eff)     Mass\n";
        let si = "Si   14      0.000000      0.000000      0.000000   4.0000  28.0855";
        let ox = "O     8      1.000000      0.000000      0.000000   6.0000  15.9994";
        let (first, second) = if si_first { (si, ox) } else { (ox, si) };
        // 力跟着各自的原子走：Si 受 +x，O 受 -x
        let (f1, f2) = if si_first { (1.0, -1.0) } else { (-1.0, 1.0) };
        format!(
            "{head}      1    1 {first}\n      2    2 {second}\n\n\
             \x20*** SCF run converged in     4 steps ***\n\n\
             \x20ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]            -10.000000000000000\n\n\
             \x20FORCES| Atomic forces [hartree/bohr]\n\
             \x20FORCES|   Atom     x               y               z               |f|\n\
             \x20FORCES|      1 {f1:.8E}  0.00000000E+00  0.00000000E+00   1.00000000E-02\n\
             \x20FORCES|      2 {f2:.8E}  0.00000000E+00  0.00000000E+00   1.00000000E-02\n\
             \x20FORCES| Sum     0.00000000E+00  0.00000000E+00  0.00000000E+00\n",
            f1 = f1 * 1e-2,
            f2 = f2 * 1e-2,
        )
    }

    /// A plain file where the system directory should go is refused by name, not with
    /// the raw io error of listing it (review E2-11).
    #[test]
    fn a_file_in_the_way_of_the_system_directory_is_refused() {
        let dir = std::env::temp_dir().join("ferro_collect_file_in_the_way");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("a.out"), sp_text(true)).unwrap();
        let dest = dir.join("out.db");
        std::fs::write(&dest, "occupied").unwrap();
        let group = Group { dir: dir.clone(), rel: PathBuf::new(), files: vec![dir.join("a.out")] };
        let err = collect_group(&group, &dest, false, false, None, &mut Vec::new()).unwrap_err();
        assert!(format!("{err:#}").contains("--overwrite"), "{err:#}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Two files of one composition written in different atom orders are one system.
    ///
    /// CP2K writes the atoms in whatever order the input listed them, so two single
    /// points of the same material routinely disagree. Before the canonical sort this
    /// was an error message holding two identical formulas.
    #[test]
    fn files_whose_atom_order_differs_collect_into_one_system() {
        let dir = std::env::temp_dir().join("ferro_collect_atom_order");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("a.out"), sp_text(true)).unwrap();
        std::fs::write(dir.join("b.out"), sp_text(false)).unwrap();

        let group = Group {
            dir: dir.clone(),
            rel: PathBuf::new(),
            files: vec![dir.join("a.out"), dir.join("b.out")],
        };
        let dest = dir.join("out.db");
        let mut skipped = Vec::new();
        collect_group(&group, &dest, false, false, None, &mut skipped).unwrap();
        assert!(skipped.is_empty());

        let (traj, _) = read_deepmd_npy_with_warnings(&dest).unwrap();
        assert_eq!(traj.n_frames(), 2);
        // 规范序是 (Z, 符号)：O(8) 在前，Si(14) 在后 —— 与写出的顺序无关
        let elems: Vec<&str> = traj.frames[0].atoms().map(|a| a.element).collect();
        assert_eq!(elems, vec!["O", "Si"]);

        // 两份文件描述的是同一个构型，重排之后两帧必须逐位相同
        let (f0, f1) = (&traj.frames[0], &traj.frames[1]);
        for i in 0..2 {
            assert!(
                (f0.positions()[i] - f1.positions()[i]).norm() < 1e-12,
                "atom {i}: {:?} vs {:?}",
                f0.positions()[i],
                f1.positions()[i]
            );
            let (a, b) = (f0.forces().unwrap()[i], f1.forces().unwrap()[i]);
            assert!((a - b).norm() < 1e-12, "force {i}: {a:?} vs {b:?}");
        }
        // 力确实跟着原子走了，而不是两帧都成了零
        assert!(f0.forces().unwrap()[0].x < 0.0, "O 受 -x");
        assert!(f0.forces().unwrap()[1].x > 0.0, "Si 受 +x");

        let _ = std::fs::remove_dir_all(&dir);
    }

    fn collect_dir(name: &str, files: &[(&str, String)]) -> (PathBuf, Result<Trajectory>) {
        let dir = std::env::temp_dir().join(name);
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        for (file, text) in files {
            std::fs::write(dir.join(file), text).unwrap();
        }
        let group = Group {
            dir: dir.clone(),
            rel: PathBuf::new(),
            files: files.iter().map(|(f, _)| dir.join(f)).collect(),
        };
        let dest = dir.join("out.db");
        let mut skipped = Vec::new();
        let read = collect_group(&group, &dest, false, false, None, &mut skipped)
            .and_then(|_| Ok(read_deepmd_npy_with_warnings(&dest)?.0));
        (dir, read)
    }

    /// 一批单点里有一个没收敛：reader 对「一帧不剩」报错，collect 把这份文件 SKIP 掉，
    /// 不能让它以空组成去和别的文件比、使整个目录失败。排在第一个的就是它时也一样
    #[test]
    fn a_file_whose_frames_were_all_dropped_does_not_fail_the_directory() {
        let failed = sp_text(true).replace("SCF run converged", "SCF run NOT converged");
        let (dir, read) = collect_dir("ferro_collect_all_dropped",
            &[("a.out", failed), ("b.out", sp_text(true)), ("c.out", sp_text(false))]);
        assert_eq!(read.unwrap().n_frames(), 2);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// QE pw.x 实算输出：单点按文件成帧（没收敛的那份被 SKIP），md 逐步成帧，
    /// virial = 应力 × 体积（同号，正 = 受压）
    #[test]
    fn quantum_espresso_outputs_collect() {
        let fixture = |n: &str| std::fs::read_to_string(format!("../tests/{n}")).unwrap();
        let (dir, read) = collect_dir("ferro_collect_qe_sp",
            &[("a.out", fixture("qe75_sp_noconv.out")), ("b.out", fixture("qe75_sp.out"))]);
        let sp = read.unwrap();
        assert_eq!(sp.n_frames(), 1);
        assert!((sp.frames[0].energy.unwrap() - -3567.569297968025).abs() < 1e-8);
        let _ = std::fs::remove_dir_all(&dir);

        let (dir, read) = collect_dir("ferro_collect_qe_md", &[("pw.out", fixture("qe75_md.out"))]);
        let md = read.unwrap();
        assert_eq!(md.n_frames(), 5);
        let direct = ferro_io::read_qe_out(Path::new("../tests/qe75_md.out")).unwrap();
        for (a, b) in md.frames.iter().zip(&direct.frames) {
            assert!((a.energy.unwrap() - b.energy.unwrap()).abs() < 1e-9);
            let ds = (a.stress.unwrap() - b.stress.unwrap()).abs().max();
            assert!(ds < 1e-12, "virial.npy 读回的应力与 reader 一致（差 {ds}）");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn pjm_log_names_are_recognised_and_nothing_else() {
        assert!(is_pjm_log("x0-job.123.out"));
        assert!(is_pjm_log("x1-a.b-c.99.out"), "jobname 可含 . 与 -");
        assert!(!is_pjm_log("x1-run.out"), "无 jobid 的是用户自己的文件");
        assert!(!is_pjm_log("x0-.123.out"), "jobname 不能为空");
        assert!(!is_pjm_log("x0-job.12a.out"), "jobid 必须全是数字");
        assert!(!is_pjm_log("total.out"));
    }

    #[test]
    fn a_directory_of_only_pjm_logs_forms_no_group() {
        let paths: Vec<PathBuf> = ["/s/run1/total.out", "/s/run1/x0-j.1.out", "/s/run2/x1-j.2.out"]
            .iter()
            .map(PathBuf::from)
            .collect();
        let kept = drop_pjm_logs(paths).unwrap();
        assert_eq!(kept, vec![PathBuf::from("/s/run1/total.out")]);
        assert_eq!(group_by_directory(&kept).len(), 1, "run2 只有 PJM 日志，不应成组");
    }

    #[test]
    fn all_inputs_being_pjm_logs_is_an_error() {
        let paths = vec![PathBuf::from("x0-j.1.out"), PathBuf::from("x1-j.1.out")];
        let err = drop_pjm_logs(paths).unwrap_err().to_string();
        assert!(err.contains("all 2 input(s)"), "报错应点明全部被剔除：{err}");
    }

    #[test]
    fn a_single_directory_leaves_an_empty_name() {
        // One group: the shared ancestor is the whole path; collect_dest names it
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
    fn a_group_sitting_at_the_ancestor_is_named_by_its_directory() {
        // Otherwise its rel is empty and it lands in `<-o>.db`, outside the
        // `<-o>/` that holds its sibling groups.
        assert_eq!(rels(&["/s/md/x.out", "/s/md/sub/y.out"]), vec!["md", "md/sub"]);
    }

    #[test]
    fn a_bare_file_name_groups_under_the_current_directory() {
        // `Path::new("a.out").parent()` is `Some("")`, not `None`. cargo runs
        // tests from the crate directory, so `.` canonicalises to `ferro-cli`.
        assert_eq!(rels(&["a.out", "src/b.out"]), vec!["ferro-cli", "ferro-cli/src"]);
        let groups = group_by_directory(&[PathBuf::from("a.out")]);
        let dest = collect_dest(None, &groups[0]);
        assert_eq!(dest.file_name().unwrap(), "ferro-cli.db", "got {}", dest.display());
    }

    #[test]
    fn the_common_ancestor_stops_at_the_first_difference() {
        let dirs = vec![PathBuf::from("/s/a/md"), PathBuf::from("/s/a/opt")];
        assert_eq!(common_ancestor(&dirs), PathBuf::from("/s/a"));
        let disjoint = vec![PathBuf::from("/x/a"), PathBuf::from("/y/b")];
        assert_eq!(common_ancestor(&disjoint), PathBuf::from("/"));
    }

    /// The default destination sits beside the AIMD directory and carries `.db`.
    #[test]
    fn collect_defaults_beside_the_input_and_never_doubles_the_suffix() {
        let g = |dir: &str, rel: &str| Group {
            dir: PathBuf::from(dir),
            rel: PathBuf::from(rel),
            files: vec![],
        };
        // 不给 -o：同级，加 .db
        assert_eq!(collect_dest(None, &g("data/md", "")), PathBuf::from("data/md.db"));
        // 给了 -o：重建目录结构，末级同样加 .db
        assert_eq!(
            collect_dest(Some(Path::new("out")), &g("data/md", "a/md")),
            PathBuf::from("out/a/md.db")
        );
        // An empty rel is a single group: it still goes inside -o, under its own name.
        // It used to become `out.db` beside an empty `out/`.
        assert_eq!(
            collect_dest(Some(Path::new("out")), &g("data/md", "")),
            PathBuf::from("out/md.db")
        );
        // 已带 .db 的不叠加
        assert_eq!(collect_dest(None, &g("data/md.db", "")), PathBuf::from("data/md.db"));
    }
}
