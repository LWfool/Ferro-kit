//! `ferro dataset filter`: quality selection on an existing dataset.

use super::*;

#[derive(Args, Debug)]
pub struct FilterCmd {
    /// System directories, or a directory holding them (searched recursively)
    #[arg(short, long, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Output root; each system is rebuilt under its path relative to -i.
    /// Omit for a read-only run that only reports.
    #[arg(short, long, value_name = "DIR")]
    pub output: Option<PathBuf>,

    /// Create the output root without asking (required when there is no terminal)
    #[arg(long)]
    pub mkdir: bool,

    // Option 而非 default_value_t（审查 D-4）：要分清「用户给了」与「默认值」——
    // 默认值遇到缺该标签的 system 只关掉这一条，显式给了才报错
    /// Drop frames whose largest force magnitude exceeds this, eV/A; 0 = off   [default: 20]
    #[arg(short = 'f', long, value_name = "EV_PER_A")]
    pub f_max: Option<f64>,

    /// Drop frames whose largest |stress component| exceeds this, GPa; 0 = off  [default: 10]
    // 没有短名：全仓 -s 是批次标记，`filter -s 5` 曾被读成 5 GPa（审查 E2-2）
    #[arg(long, value_name = "GPA")]
    pub s_max: Option<f64>,

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

    /// Seed for --shuffle and for --ratio                        [default: 666]
    #[arg(long, value_name = "N")]
    pub seed: Option<u64>,

    /// Frames per output set; 0 keeps everything in one set
    #[arg(long, value_name = "N", default_value_t = DEFAULT_SET_SIZE)]
    pub set_size: usize,

    /// Allow writing into an existing non-empty output directory
    #[arg(long)]
    pub overwrite: bool,

    /// What to write
    #[arg(long = "type", value_enum, default_value_t = OutType::Deepmd)]
    pub out_type: OutType,

    /// Split each system train:valid:test, e.g. 8:1:1 (or 9:1 for train:test)
    #[arg(long, value_name = "A:B:C")]
    pub ratio: Option<String>,
}

impl FilterCmd {
    /// `-f` 的实际值，没给取默认
    fn f_max(&self) -> f64 { self.f_max.unwrap_or(20.0) }
    /// `-s` 的实际值（GPa），没给取默认
    fn s_max(&self) -> f64 { self.s_max.unwrap_or(10.0) }
}

pub(super) fn run_filter(args: &FilterCmd) -> Result<usize> {
    // Parameter errors surface before the first system is read; the value ranges are
    // FilterParams::validate's alone (review C-2)
    let manual_rcut = match args.al6.as_deref() {
        None | Some("auto") => None,
        Some(v) => Some(
            v.parse::<f64>()
                .with_context(|| format!("--al6 expects a cutoff in Angstrom or nothing, got `{v}`"))?,
        ),
    };
    let split = filter_split(args)?;
    if args.seed.is_some() && !args.shuffle && split.is_off() {
        bail!("--seed only means something with --shuffle or --ratio");
    }
    if !split.is_off() && args.output.is_none() {
        bail!("--ratio needs an output directory (-o DIR); a read-only run writes nothing");
    }

    let params = FilterParams {
        f_max: args.f_max(),
        // CLI 收 GPa（用起来顺手），内部一律 eV/Å³
        s_max: convert_pressure(args.s_max(), PressureUnit::GPa, PressureUnit::EVPerAng3),
        start: args.start.unwrap_or(0),
        end: args.end,
        stride: args.stride.unwrap_or(1),
        number: args.number,
        oo_min: args.oo_min.unwrap_or(0.0),
        shuffle: args.shuffle.then(|| args.seed.unwrap_or(DEFAULT_SEED)),
        // 每个 system 各算各的，此处只放手动值
        al6_rcut: manual_rcut,
    };
    params.validate()?;

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
    // merge 侧同一条守卫:只看路径,故在读第一条轨迹之前就判得出
    if !split.is_off() {
        refuse_held_out_inputs(&jobs.iter().map(|(sys, _)| sys.clone()).collect::<Vec<_>>())?;
    }

    if let Some(out) = &args.output {
        crate::outpath::ensure_dir(out, args.mkdir)?;
    } else {
        println!("(read-only: no -o given, nothing will be written)\n");
    }

    let mut failures = 0usize;
    let mut auto_rcuts: Vec<f64> = Vec::new();
    // 报告按表名分组堆叠，行标签用**相对路径**而不是目录名 ——
    // 嵌套结构下 a/md 与 b/md 的叶子名相同，堆起来就分不出是谁
    let mut order: Vec<String> = Vec::new();
    let mut groups: Vec<Vec<(String, ferro_core::Table)>> = Vec::new();
    // 逐 system 的统计（O–O 分位、实际所用 Al–O 截断）进 [inputs] 的列，
    // 不进表头：表头堆叠后只剩第一个 system 的，会冒充整批
    let mut columns = vec!["frames_out"];
    if params.oo_min > 0.0 {
        columns.extend(OO_COLUMNS);
    }
    if args.al6.is_some() {
        columns.push("al6_rcut");
    }
    let mut summary = crate::batch::Summary::new(&columns);

    for (sys, rel) in &jobs {
        let label = rel.display().to_string();
        match filter_one(sys, rel, args, &params) {
            Ok(one) => {
                auto_rcuts.extend(one.rcut);
                let mut values = vec![one.n_kept as f64];
                if params.oo_min > 0.0 {
                    values.extend(one.min_oo.unwrap_or([f64::NAN; 6]));
                }
                if args.al6.is_some() {
                    values.push(one.rcut.or(params.al6_rcut).unwrap_or(f64::NAN));
                }
                summary.ok(label.clone(), one.n_input, one.n_atoms, &values);
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
                crate::report::error(sys.display(), format!("{e:#}; skipped"));
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
    if let Some(out) = &args.output {
        let mut tables = Vec::with_capacity(order.len());
        for (name, parts) in order.into_iter().zip(groups) {
            let merged = ferro_core::Table::concat_union("system", parts)?;
            tables.push((name, merged));
        }
        let outp = crate::batch::Output {
            mkdir: true,
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
        crate::report::error(crate::report::PROGRAM, format!("{failures} of {} system(s) failed", jobs.len()));
    }
    Ok(failures)
}

/// The shared parameter block at the top of every report file.
fn report_params(args: &FilterCmd, params: &FilterParams) -> Vec<String> {
    let off = |v: f64| if v > 0.0 { format!("{v}") } else { "off".to_string() };
    let mut v = vec![
        format!("f_max     = {} eV/Ang", off(args.f_max())),
        format!("s_max     = {} GPa", off(args.s_max())),
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

/// `[inputs]` columns for the quantiles of the per-frame minimum O–O distance,
/// in the order [`quantiles`] returns them.
const OO_COLUMNS: [&str; 6] =
    ["min_oo_min", "min_oo_p1", "min_oo_p50", "min_oo_p99", "min_oo_max", "min_oo_mean"];

/// What one system contributed: its report tables and its frame counts.
struct FilterOne {
    /// The automatically derived Al-O cutoff, when one was derived.
    rcut: Option<f64>,
    /// Quantiles of the per-frame minimum O–O distance, when `--oo-min` is on.
    min_oo: Option<[f64; 6]>,
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
    let (traj, warnings) = read_deepmd_npy(sys)?;
    crate::report::notices(sys, &warnings);

    // --al6 不带值：从这个 system 自己的 Al-O RDF 取第一壳层的外沿。
    // 逐 system 各算各的 —— 成分不同，壳层位置本来就不同
    let mut params = params.clone();
    // 审查 D-4：-f / -s 默认开着，无 virial 的数据集很常见。用户没给时，缺该标签的 system
    // 只关掉这一条并说明；显式给了才交给 filter_frames 按「给了阈值却缺标签」报错
    if args.s_max.is_none() && params.s_max > 0.0 && traj.frames.iter().any(|f| f.stress.is_none()) {
        params.s_max = 0.0;
        crate::report::note(sys.display(), "no virial; the default stress criterion is off for it");
    }
    if args.f_max.is_none() && params.f_max > 0.0 && traj.frames.iter().any(|f| f.forces().is_none()) {
        params.f_max = 0.0;
        crate::report::note(sys.display(), "no forces; the default force criterion is off for it");
    }
    let mut derived = None;
    if args.al6.as_deref() == Some("auto") {
        let shell = first_shell_cutoffs(&traj, &[("Al", "O")])
            .map_err(|e| anyhow::anyhow!("{e}"))?[0]
            .context("no Al-O pair in this system, so --al6 has no cutoff to derive")?;
        println!(
            "  Al-O first shell: peak {:.2} A (g={:.1}), cutoff {:.2} A (g={:.3})",
            shell.peak_r, shell.peak_g, shell.min_r, shell.depth
        );
        if shell.depth > 0.5 {
            crate::report::warn(sys.display(), format!(
                "the Al-O g(r) minimum at {:.2} A is shallow (g={:.2}); Al-O may have no clear shell here",
                shell.min_r, shell.depth
            ));
        }
        params.al6_rcut = Some(shell.min_r);
        derived = Some(shell.min_r);
    }

    let result = filter_frames(&traj, &params)?;
    println!("{}", sys.display());
    print_report(&result);

    // 诊断表恒算：实测 1110 帧 / 302 原子挂钟时间与不算时相同（rayon 跑满），
    // 而它是选阈值的依据，只在只读模式算就等于永远落不了盘
    let diagnostics = filter_diagnostics(&traj, &result, &params)?;
    let mut tables = result.to_tables();
    tables.extend(diagnostics.iter().cloned());
    let min_oo = if params.oo_min > 0.0 {
        quantiles(&result.verdicts.iter().filter_map(|x| x.min_oo).collect::<Vec<_>>())
    } else {
        None
    };

    if args.output.is_none() {
        // 只读模式：四张诊断表也打出来，但一个字不落盘
        for (name, table) in &diagnostics {
            println!("  [{name}]");
            if let (Some(q), "min_oo") = (min_oo, name.as_str()) {
                println!(
                    "    min d(O-O) [A]: min {:.3}  p1 {:.3}  p50 {:.3}  p99 {:.3}  max {:.3}  mean {:.3}",
                    q[0], q[1], q[2], q[3], q[4], q[5]
                );
            }
            print_table(table);
        }
        println!();
    }

    let one = FilterOne {
        rcut: derived,
        min_oo,
        tables,
        n_input: result.n_input,
        n_kept: result.keep.len(),
        n_atoms: traj.frames.first().map(|f| f.n_atoms()).unwrap_or(0),
    };

    let Some(out_root) = &args.output else {
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
        Some(spec) => Ok(Split::parse(spec, seed)?),
    }
}

fn print_report(r: &FilterResult) {
    for line in r.meta_lines() {
        println!("  {line}");
    }
    for (name, table) in r.to_tables() {
        println!("  [{name}]");
        print_table(&table);
    }
}

fn print_table(t: &ferro_core::Table) {
    for line in t.to_comment_lines() {
        println!("    {line}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::test_ids::*;

    /// 审查 C-D4 连带：`-f nan` / `--oo-min nan` 曾静默关掉判据，`--al6 nan` 要删光帧才失败
    #[test]
    fn non_finite_thresholds_fail_before_any_input_is_read() {
        #[derive(clap::Parser)]
        struct Wrap {
            #[command(flatten)]
            cmd: FilterCmd,
        }
        use clap::Parser;
        for (args, want) in [
            ("-f nan", "f-max must be a finite number"),
            ("-f inf", "f-max must be a finite number"),
            ("--s-max nan", "s-max must be a finite number"),
            ("--oo-min nan", "oo-min must be a finite number"),
            ("--oo-min inf", "oo-min must be a finite number"),
            ("--al6 nan", "al6 must be a finite number"),
            ("--al6 inf", "al6 must be a finite number"),
        ] {
            let mut argv = vec!["filter", "-i", "no_such_system_dir"];
            argv.extend(args.split_whitespace());
            let w = Wrap::try_parse_from(&argv).unwrap_or_else(|e| panic!("{args}: {e}"));
            let err = format!("{:#}", run_filter(&w.cmd).expect_err(args));
            assert!(err.contains(want), "{args}：应报「{want}」，实际 {err}");
        }
    }

    /// `--stride 0`, `-N 0` and `--start` past `--end` fail before any input is read, as
    /// they do in convert (review D-1)
    #[test]
    fn bad_frame_selection_fails_before_any_input_is_read() {
        #[derive(clap::Parser)]
        struct Wrap {
            #[command(flatten)]
            cmd: FilterCmd,
        }
        use clap::Parser;
        for (args, want) in [
            ("--stride 0", "stride must be >= 1"),
            ("-N 0", "number must be >= 1"),
            ("--start 3 --end 1", "start 3 is past end 1"),
        ] {
            let mut argv = vec!["filter", "-i", "no_such_system_dir"];
            argv.extend(args.split_whitespace());
            let w = Wrap::try_parse_from(&argv).unwrap_or_else(|e| panic!("{args}: {e}"));
            let err = format!("{:#}", run_filter(&w.cmd).expect_err(args));
            assert!(err.contains(want), "{args}: expected `{want}`, got {err}");
        }
    }

    /// 审查 C-D6：多 system 的报告曾只留第一个 system 的表头，
    /// 它的 O–O 分位数与 Al 截断摆在头部冒充整批
    #[test]
    fn per_system_statistics_go_to_inputs_not_the_table_header() {
        use ferro_core::{Atom, Cell, Frame};
        use nalgebra::{Matrix3, Vector3};
        #[derive(clap::Parser)]
        struct Wrap {
            #[command(flatten)]
            cmd: FilterCmd,
        }
        use clap::Parser;

        let root = std::env::temp_dir().join(format!("ferro_cd6_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        // 两个 system 各一个 AlO6 八面体，Al–O 分别 1.8 / 2.0 Å，故最近 O–O 为 d·√2；
        // 要有 Al6：--al6 删的是不含六配位 Al 的帧
        for (name, d) in [("a", 1.8), ("b", 2.0)] {
            let cell = Cell::from_matrix(Matrix3::from_diagonal_element(10.0));
            let frames = (0..3)
                .map(|_| {
                    let mut f = Frame::with_cell(cell.clone(), [true; 3]);
                    let c = Vector3::new(5.0, 5.0, 5.0);
                    f.set_atoms(vec![Atom::new("Al", c)]);
                    for k in 0..3 {
                        for sgn in [1.0, -1.0] {
                            let mut v = Vector3::zeros();
                            v[k] = sgn * d;
                            f.add_atom(Atom::new("O", c + v));
                        }
                    }
                    f.energy = Some(-1.0);
                    f.set_forces(Some(vec![Vector3::zeros(); 7])).unwrap();
                    f.stress = Some(Matrix3::zeros());
                    f
                })
                .collect();
            let t = Trajectory { frames, metadata: Default::default() };
            ferro_io::write_deepmd_npy(&t, &root.join("in").join(name)).unwrap();
        }
        let (inp, out) = (root.join("in"), root.join("out"));
        let argv = [
            "filter", "-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap(),
            "--mkdir", "--oo-min", "1.0", "--al6", "2.3",
        ];
        let w = Wrap::try_parse_from(argv).unwrap();
        assert_eq!(run_filter(&w.cmd).unwrap(), 0);

        let text = std::fs::read_to_string(out.join("filter_min_oo.csv")).unwrap();
        assert!(!text.contains("min d(O-O)"), "分位数不应留在表头：\n{text}");
        let text = std::fs::read_to_string(out.join("filter_al_cn.csv")).unwrap();
        assert!(!text.contains("rcut ="), "截断不应留在表头：\n{text}");

        // [inputs] 两行各带自己的分位数与截断
        let row = |sys: &str| -> String {
            text.lines()
                .find(|l| l.split_whitespace().nth(1) == Some(sys))
                .unwrap_or_else(|| panic!("[inputs] 缺 {sys} 行：\n{text}"))
                .to_string()
        };
        let header = text.lines().find(|l| l.contains("min_oo_p50")).expect("[inputs] 缺 min_oo_p50 列");
        let col = |name: &str| header.split_whitespace().position(|c| c == name).unwrap();
        for (sys, d) in [("a", 1.8f64), ("b", 2.0)] {
            let r = row(sys);
            let cells: Vec<&str> = r.split_whitespace().collect();
            let p50: f64 = cells[col("min_oo_p50")].parse().unwrap();
            assert!((p50 - d * 2f64.sqrt()).abs() < 1e-3, "{sys} 的 p50：{r}");
            assert_eq!(cells[col("al6_rcut")].parse::<f64>().unwrap(), 2.3, "{sys} 的截断：{r}");
        }
        let _ = std::fs::remove_dir_all(&root);
    }

    /// 审查 D-4：`-s` 默认开着 10 GPa，无 virial 的数据集按默认参数整批失败，报错却说「你给了阈值」
    #[test]
    fn a_default_threshold_is_switched_off_for_a_system_without_that_label() {
        use ferro_core::{Atom, Cell, Frame};
        use nalgebra::{Matrix3, Vector3};
        #[derive(clap::Parser)]
        struct Wrap {
            #[command(flatten)]
            cmd: FilterCmd,
        }
        use clap::Parser;

        let root = std::env::temp_dir().join(format!("ferro_d4_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let cell = Cell::from_matrix(Matrix3::from_diagonal_element(10.0));
        let frames = (0..3)
            .map(|_| {
                let mut f = Frame::with_cell(cell.clone(), [true; 3]);
                f.set_atoms(vec![Atom::new("O", Vector3::new(1.0, 1.0, 1.0)), Atom::new("H", Vector3::new(2.0, 1.0, 1.0))]);
                f.energy = Some(-1.0);
                f.set_forces(Some(vec![Vector3::zeros(); 2])).unwrap();
                f
            })
            .collect();
        let sys = root.join("novirial");
        ferro_io::write_deepmd_npy(&Trajectory { frames, metadata: Default::default() }, &sys).unwrap();
        let run = |extra: &[&str]| {
            let mut argv = vec!["filter", "-i", sys.to_str().unwrap()];
            argv.extend_from_slice(extra);
            run_filter(&Wrap::try_parse_from(argv).unwrap().cmd)
        };
        // 默认参数：应力判据对这个 system 关掉，照常完成
        assert_eq!(run(&[]).unwrap(), 0, "默认参数下无 virial 的 system 不应失败");
        // 显式给了 -s：仍按「给了阈值却缺标签」处理
        assert!(!matches!(run(&["--s-max", "10"]), Ok(0)), "显式 --s-max 时应失败");
        let _ = std::fs::remove_dir_all(&root);
    }

    /// filter 删掉的正是坏帧；留下的帧打乱、切 set 之后各量仍属同一帧同一原子
    #[test]
    fn filter_drops_exactly_the_bad_frames_and_keeps_the_rest_intact() {
        #[derive(clap::Parser)]
        struct Wrap {
            #[command(flatten)]
            cmd: FilterCmd,
        }
        use clap::Parser;
        let root = std::env::temp_dir().join(format!("ferro_filter_ids_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let mut t = id_system(&[0, 1, 2, 3, 4, 5], &[1, 0, 2]);
        // 帧 1、4 的力超过默认的 20 eV/Å
        for g in [1, 4] {
            for fv in t.frames[g].forces_mut().unwrap() {
                fv.z = 50.0;
            }
        }
        ferro_io::write_deepmd_npy(&t, &root.join("in/sys")).unwrap();
        let (inp, out) = (root.join("in"), root.join("out"));
        let argv = [
            "filter", "-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap(),
            "--mkdir", "--shuffle", "--set-size", "2",
        ];
        let w = Wrap::try_parse_from(argv).unwrap();
        assert_eq!(run_filter(&w.cmd).unwrap(), 0);

        let mut got = check_ids(&out);
        got.sort();
        assert_eq!(got, vec![0, 2, 3, 5], "应恰好删掉帧 1 与 4");
        let _ = std::fs::remove_dir_all(&root);
    }
}
