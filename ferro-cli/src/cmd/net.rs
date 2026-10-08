//! `ferro net` — glass network topology.
//!
//! A leaf command, not a subcommand group. The old `net qn` / `net type` split was
//! never two analyses: both classify every atom of every frame the same way, and
//! `type` merely wrote the classification out instead of summarising it. So the
//! export is a **flag**, `--export-traj`, and the statistics always run.
//!
//! Cutoffs are given as `--<Former>-<Ligand>=<Å>` (e.g. `--P-O=2.3`), which clap cannot
//! model as a fixed flag set: the element pair is part of the flag name. `main` strips
//! those out of argv before clap parses and hands them here.
//!
//! A cutoff may also be `auto` (`--P-O=auto`): the first minimum of that pair's g(r)
//! behind its first peak, taken from **each input separately**. The shell edge moves
//! with composition, and one shared hard cutoff either cuts into the shell of one
//! composition or reaches past it in another, which shifts the coordination fractions.

use anyhow::{anyhow, bail, Result};
use clap::{Args, ValueEnum};
use ferro_analysis::md::{first_shell_cutoffs, ShellCutoff};
use ferro_analysis::{calc_network, NetworkResult};
use ferro_core::{Trajectory, TypeParams};
use ferro_io::{write_extxyz, write_lammps_dump, LammpsUnits};
use ferro_structure::{apply_type_labels, classify_trajectory, fold_labels};
use std::collections::BTreeMap;

use crate::args::common::CommonArgs;
use crate::batch::{self, Summary};

#[derive(Args, Debug)]
pub struct NetCmd {
    #[command(flatten)]
    pub common: CommonArgs,

    /// Modifier elements, comma separated (e.g. Zn,Na). They count towards
    /// coordination numbers but take no part in bridging counts or ligand types
    #[arg(long)]
    pub modifier: Option<String>,

    /// Formers reported as a Qn speciation, comma separated. REPLACES the default
    /// list (B,P,Si); every other former is described by its coordination number
    #[arg(long)]
    pub qn: Option<String>,

    /// Also write the classified trajectory, one file per input:
    /// <input stem>_types[_<suffix>].<ext>. Defaults to lammpstrj
    #[arg(long, value_enum, num_args = 0..=1, default_missing_value = "lammpstrj")]
    pub export_traj: Option<ExportFormat>,
}

/// Where a classified trajectory can go, and how the label survives the trip.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum ExportFormat {
    /// One name column only, so the label replaces the element; the dump reader
    /// splits `<element>_<suffix>` back apart, so nothing downstream sees a new column
    Lammpstrj,
    /// Self-describing columns, so the label gets its own `label:S:1` and
    /// `species` stays a pure element symbol — lossless in both directions
    Extxyz,
}

impl ExportFormat {
    fn ext(self) -> &'static str {
        match self {
            ExportFormat::Lammpstrj => "lammpstrj",
            ExportFormat::Extxyz => "extxyz",
        }
    }
}

/// One input's result and the `auto` shells it was classified with (`None`: pair absent).
type Resolved = (NetworkResult, Vec<Option<ShellCutoff>>);

pub fn wants_help(cmd: &NetCmd) -> bool {
    cmd.common.input.is_empty()
}

pub fn print_help() {
    println!("{}", HELP_EXTRA);
}

/// Runs the network analysis. `pair_args` are the `--P-O=2.3`-style cutoffs `main`
/// pulled out of argv.
pub fn run(cmd: &NetCmd, pair_args: &[String]) -> Result<usize> {
    // 参数错误在读第一个文件之前就失败
    if pair_args.is_empty() {
        bail!("No pair cutoffs specified. Use --Former-Ligand=cutoff, e.g. --P-O=2.3");
    }
    let (params, auto) = build_params(pair_args, cmd.modifier.as_deref(), cmd.qn.as_deref())?;
    if params.cutoffs.is_empty() {
        bail!("Every cutoff names a modifier element; at least one former is required");
    }

    // 参数（含 --atom-style）在建目录、读第一个文件之前查完
    let inputs = cmd.common.inputs()?;
    // net 没有类型选择,故无 label 段;-o 的目录在读第一个文件之前建好
    let out = cmd.common.out(None);
    out.prepare()?;

    cmd.common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());
    print_label_scheme(&params);

    let (results, failures) = batch::map_inputs(&inputs, |inp| {
        let traj = cmd.common.load(&inp.path)?;
        let (params, shells) = resolve_auto(&params, &auto, &traj)?;
        let result = calc_network(&traj, &params)?;
        if let Some(fmt) = cmd.export_traj {
            export_labelled(&traj, &params, &inp.label, &out, fmt, cmd.common.read.units)?;
        }
        Ok((result, shells))
    });
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }

    let tables = batch::stack(&results, |(r, _): &Resolved| {
        Ok(r.to_tables())
    })?;
    note_missing_qn_tables(&params);
    warn_edge_sharing(&results);

    let mut summary = Summary::new(&[]);
    for (input, (r, shells)) in &results {
        summary.ok(input.label.clone(), r.n_frames, r.n_atoms, &[]);
        // auto 截断逐输入各取各的，故值进 [inputs] 而不是共享头部；g_min 是
        // 极小处的 g(r)，读者据此判断那个极小是不是一道干净的空隙
        if !auto.is_empty() {
            let pairs = || auto.iter().zip(shells);
            let show = |f: fn(&ShellCutoff) -> String| pairs()
                .map(|((a, b), s)| format!("{a}-{b}={}", s.as_ref().map_or("-".into(), f)))
                .collect::<Vec<_>>().join(" ");
            summary.note("cutoff", show(|s| format!("{:.3}", s.min_r)));
            summary.note("g_min", show(|s| format!("{:.3}", s.depth)));
        }
        summary.note("mean_qn", fmt_means(&r.mean_qn, &r.qn_dist));
        // mean_n_bo 与 mean_qn 并列：前者是桥氧个数(旧口径的那个数)，后者只数
        // 同元素连接。两者在无三簇氧的纯磷酸盐里相等，混合体系里分叉，摆在
        // 一行让读者一眼看出这条轨迹的异核桥有多少
        summary.note("mean_n_bo", fmt_means(&r.mean_n_bo, &r.cn_dist));
        summary.note("mean_cn", fmt_means(&r.mean_cn, &r.cn_dist));
    }
    summary.failed(&failures);

    batch::write_all(
        "network",
        "Glass Network Topology — Qn speciation, ligand types, coordination numbers",
        &results[0].1.0.meta_lines(&auto),
        &summary.into_table(),
        tables,
        &out,
    )?;

    Ok(failures.len())
}

/// Prints the label scheme once per run, with this run's elements filled in.
///
/// It lives here rather than in the help text because it is **derived from the
/// arguments**: which element gets a Qn and which gets a coordination number depends
/// on `--qn` and on the cutoffs given, so a static block in `--help` would have to
/// describe every case in the abstract. Printed once, above the results, it says what
/// the labels in *this* run's files mean.
fn print_label_scheme(params: &TypeParams) {
    let join = |v: Vec<String>| -> String {
        if v.is_empty() { "-".to_string() } else { v.join(",") }
    };
    let qn: Vec<String> = params.qn_formers();
    let cn_formers: Vec<String> = params.formers().into_iter()
        .filter(|e| !params.is_qn_former(e)).collect();

    println!("Labels:");
    println!("  {:<10} <elem>_<Qn>   digit = HOMOPOLAR connections (P-O-P), the n",
             join(qn.clone()));
    if !qn.is_empty() {
        println!("  {:<10}               of Q^n_m; P-O-Al etc. are the m_ columns", "");
    }
    if !cn_formers.is_empty() {
        println!("  {:<10} <elem>_<CN>   digit = COORDINATION number, not Qn",
                 join(cn_formers));
    }
    println!("  {:<10} _f free  _n non-bridging  _b bridging  _t tricluster",
             join(params.ligands()));
    if !params.modifiers().is_empty() {
        println!("  {:<10} bare element symbol, no role suffix",
                 join(params.modifiers()));
    }
}

/// Warns when the Qn tables will not be written, and why.
///
/// A run whose formers are all coordination-described produces four tables instead of
/// six. Silence would read as a bug; an empty `network_qn.csv` would read as "measured,
/// and the answer was zero".
fn note_missing_qn_tables(params: &TypeParams) {
    if !params.qn_formers().is_empty() { return; }
    println!(
        "        note: no former is a Qn element ({} given, default list is B,P,Si), \n\
         \x20             so network_qn.csv and network_qn_partner.csv are not written.\n\
         \x20             Use --qn <ELEM> to report one of them as a Qn speciation.",
        params.formers().join(",")
    );
}

/// `Zn=3.98 P=4.00` — per-input, so it belongs in the `[inputs]` block rather than
/// the shared parameter block.
/// Warns when formers share two or more ligands (edge-sharing polyhedra).
///
/// The conventional Qn analysis assumes an all-corner-sharing network.  Edge
/// sharing is genuinely rare in phosphates, so the likelier cause is a cutoff
/// that reaches into the second shell — say so rather than leaving the user to
/// wonder why one neighbour produced two bridges.
fn warn_edge_sharing(results: &[(batch::Input, Resolved)]) {
    for (input, (r, _)) in results {
        if r.n_edge_sharing == 0 { continue; }
        eprintln!(
            "warning: {}: {} former pair(s) share 2+ ligands (edge-sharing).\n\
             \x20        Qn assumes corner sharing, so one neighbour here yields two bridges.\n\
             \x20        Edge sharing is very rare in phosphates — check the cutoffs first.",
            input.label.clone(), r.n_edge_sharing);
    }
}

fn fmt_means<T>(
    means: &std::collections::HashMap<String, f64>,
    present: &std::collections::HashMap<String, T>,
) -> String {
    let mut elems: Vec<&String> = present.keys().collect();
    elems.sort();
    elems.iter()
        .filter_map(|e| means.get(*e).map(|v| format!("{e}={v:.2}")))
        .collect::<Vec<_>>()
        .join(" ")
}

// ─── 标注轨迹导出 ─────────────────────────────────────────────────────────────

/// Writes the classified trajectory as `<input stem>_types[_<suffix>].<ext>`.
///
/// The name carries the input stem because this is a **one product per input**
/// product, like `ferro map`'s cubes — a fixed output path would make the second
/// input overwrite the first.
///
/// The labels live in `Atom::label` throughout; only the LAMMPS-dump branch folds
/// them into the element column, and only here — `ferro convert` keeps writing
/// clean element symbols.
fn export_labelled(
    traj: &Trajectory,
    params: &TypeParams,
    label: &str,
    out: &batch::Output,
    fmt: ExportFormat,
    units: Option<LammpsUnits>,
) -> Result<()> {
    let per_frame = classify_trajectory(traj, params);
    if per_frame.len() != traj.frames.len() {
        bail!(
            "cannot export: {} of {} frames have no cell",
            traj.frames.len() - per_frame.len(),
            traj.frames.len()
        );
    }

    let mut skipped = 0usize;
    let frames = traj.frames.iter().zip(&per_frame)
        .map(|(frame, types)| {
            let labels: Vec<String> = types.iter().map(|t| t.label()).collect();
            let labelled = apply_type_labels(frame, &labels);
            match fmt {
                ExportFormat::Lammpstrj => {
                    let (folded, n) = fold_labels(&labelled);
                    skipped += n;
                    folded
                }
                ExportFormat::Extxyz => labelled,
            }
        })
        .collect();
    let out_traj = Trajectory { frames, metadata: traj.metadata.clone() };

    let stem = label;
    let ext = fmt.ext();
    let name = match out.suffix.as_deref().filter(|s| !s.is_empty()) {
        Some(s) => format!("{stem}_types_{s}.{ext}"),
        None => format!("{stem}_types.{ext}"),
    };
    let path = out.join(&name);
    match fmt {
        // 同 convert：沿用读入时的 --units；没给而帧里有速度 / 力时由 writer 报错要求给出
        ExportFormat::Lammpstrj => write_lammps_dump(&out_traj, &path, units)?,
        ExportFormat::Extxyz => write_extxyz(&out_traj, &path)?,
    }
    if skipped > 0 {
        println!(
            "        note: {skipped} label(s) not of the form <element>_<suffix>; \
             wrote the element instead"
        );
    }
    println!("        traj -> {}", path.display());
    Ok(())
}

// ─── Pair 参数解析 ────────────────────────────────────────────────────────────

/// Splits `--Former-Ligand=cutoff` arguments out of argv; the rest goes to clap.
pub fn split_pair_args(all: &[String]) -> (Vec<String>, Vec<String>) {
    let mut pairs = Vec::new();
    let mut clap  = Vec::new();
    for arg in all {
        if is_pair_arg(arg) { pairs.push(arg.clone()); } else { clap.push(arg.clone()); }
    }
    (pairs, clap)
}

fn is_pair_arg(s: &str) -> bool {
    if !s.starts_with("--") { return false; }
    let inner = &s[2..];
    inner.starts_with(|c: char| c.is_ascii_uppercase()) && inner.contains('=')
}

fn split_elems(s: Option<&str>) -> Vec<String> {
    s.map(|s| s.split(',').map(|e| e.trim().to_string()).filter(|e| !e.is_empty()).collect())
        .unwrap_or_default()
}

/// Builds the run's parameters and lists the `auto` pairs.
///
/// An `auto` pair holds `NaN` in the returned table until [`resolve_auto`] fills it
/// from each input; the table is needed before any input is read, to validate the
/// element roles and print the label scheme, and those only look at the keys.
fn build_params(
    pair_args: &[String],
    modifier: Option<&str>,
    qn: Option<&str>,
) -> Result<(TypeParams, Vec<(String, String)>)> {
    let modifier_elems: std::collections::HashSet<String> =
        split_elems(modifier).into_iter().collect();

    let mut cutoffs = BTreeMap::new();
    let mut modifier_cutoffs = BTreeMap::new();
    let mut auto = Vec::new();
    for ((elem, ligand), cutoff) in parse_pairs(pair_args)? {
        let cutoff = cutoff.unwrap_or_else(|| {
            auto.push((elem.clone(), ligand.clone()));
            f64::NAN
        });
        if modifier_elems.contains(&elem) {
            modifier_cutoffs.insert((elem, ligand), cutoff);
        } else {
            cutoffs.insert((elem, ligand), cutoff);
        }
    }

    // 点名了修饰子却没给它截断:静默当成形成子会让氧的分类整体错位
    for m in &modifier_elems {
        if !modifier_cutoffs.keys().any(|(e, _)| e == m) {
            bail!("--modifier names {m} but no --{m}-<Ligand>=<cutoff> was given");
        }
    }

    let mut params = TypeParams::new(cutoffs, modifier_cutoffs);
    // --qn 替换默认列表而不是叠加:「本体系里 B 不当形成子」是真实需求,
    // 追加式标志没法把默认项摘出去
    if let Some(list) = qn {
        let elems = split_elems(Some(list));
        if elems.is_empty() { bail!("--qn was given an empty element list"); }
        // 两条静默失败的路都堵死,理由同 --modifier:错的分类不会报错,只会给出
        // 一份看着正常的错数据
        for e in &elems {
            if modifier_elems.contains(e) {
                bail!("--qn names {e}, but --modifier already claims it; \
                       a modifier has no bridging count, so it can have no Qn");
            }
            if !params.cutoffs.keys().any(|(f, _)| f == e) {
                bail!("--qn names {e} but no --{e}-<Ligand>=<cutoff> was given, \
                       so {e} is not a former in this run");
            }
        }
        params = params.with_qn_elements(elems);
    }
    Ok((params, auto))
}

/// Fills the `auto` cutoffs from this input's own g(r) and prints what it chose.
///
/// Each input gets its own value: the first-shell edge moves with composition, so one
/// shared number would cut into the shell of some inputs and reach past it in others.
fn resolve_auto(
    template: &TypeParams,
    auto: &[(String, String)],
    traj: &Trajectory,
) -> Result<(TypeParams, Vec<Option<ShellCutoff>>)> {
    if auto.is_empty() {
        return Ok((template.clone(), Vec::new()));
    }
    let pairs: Vec<(&str, &str)> = auto.iter().map(|(a, b)| (a.as_str(), b.as_str())).collect();
    let found = first_shell_cutoffs(traj, &pairs)?;

    let present = traj.frames.first().map(|f| f.unique_elements()).unwrap_or_default();
    let mut params = template.clone();
    let mut shells = Vec::with_capacity(auto.len());
    for ((a, b), shell) in auto.iter().zip(found) {
        let key = (a.clone(), b.clone());
        let slot = params.cutoffs.get_mut(&key)
            .or_else(|| params.modifier_cutoffs.get_mut(&key))
            .expect("auto 对来自 build_params，必在两张表之一");
        // 这一对元素在本输入里不存在（按成分扫描时常见，如不含 Al 的样品）：与给了
        // 定值截断时一样照常分析。没有这对原子，截断取什么都不影响结果；取 0 而不
        // 留 NaN，是为了不让它进 calc_network 的最小镜像上界检查
        if !present.contains(a) || !present.contains(b) {
            println!("        {a}-{b} cutoff -     (no {a}-{b} pair in this input)");
            *slot = 0.0;
            shells.push(None);
            continue;
        }
        let s = shell.ok_or_else(|| anyhow!(
            "--{a}-{b}=auto: g(r) of {a}-{b} has no first peak above 1 in this input \
             ({a} and {b} are not bonded); give the cutoff as a number"))?;
        println!(
            "        {a}-{b} cutoff {:.3} A  (first shell peak {:.2} A, g(min) = {:.3})",
            s.min_r, s.peak_r, s.depth
        );
        // 阈值与 dataset filter --al6 的同一条告警相同（cmd/dataset.rs filter_one）
        if s.depth > 0.5 {
            eprintln!(
                "warning: {a}-{b}: the g(r) minimum at {:.2} A is shallow (g = {:.2}); \
                 there is no clean gap behind the first shell, so the coordination \
                 numbers depend on where inside the trough the cutoff falls",
                s.min_r, s.depth
            );
        }
        *slot = s.min_r;
        shells.push(Some(s));
    }
    Ok((params, shells))
}

/// `None` is `auto`: the cutoff is derived from each input's g(r).
fn parse_pairs(pair_args: &[String]) -> Result<BTreeMap<(String, String), Option<f64>>> {
    let mut map = BTreeMap::new();
    for arg in pair_args {
        let inner = arg.trim_start_matches('-');
        let (pair, cutoff_str) = inner.split_once('=')
            .ok_or_else(|| anyhow!("Invalid pair argument (missing '='): {arg}"))?;
        let (former, ligand) = pair.split_once('-')
            .ok_or_else(|| anyhow!("Invalid pair argument (missing '-'): {arg}"))?;
        let cutoff = if cutoff_str == "auto" {
            None
        } else {
            let c: f64 = cutoff_str.parse()
                .map_err(|_| anyhow!("Invalid cutoff value in '{arg}' (a number or auto)"))?;
            if c <= 0.0 { bail!("Cutoff must be positive, got {c} in '{arg}'"); }
            Some(c)
        };
        map.insert((former.to_string(), ligand.to_string()), cutoff);
    }
    Ok(map)
}

const HELP_EXTRA: &str = "\
ferro net — Glass network topology

  Classifies every atom as former, ligand or modifier from the pair cutoffs you
  give, and reports six tables. At least one pair is required.

Parameters:
  --<Former>-<Ligand>=<cutoff>
                        Pair cutoff [Å], e.g. --P-O=2.4 --Al-F=2.1, or auto:
                        the first minimum of that pair's g(r), taken from each
                        input separately. At least one is required
  -i, --input  FILE...  Input trajectory files; glob patterns allowed (quote them)
  -o, --output DIR      Write every product here, tables and --export-traj alike;
                        --mkdir creates it without asking
  -s, --suffix SUFFIX   Output name suffix: network_<table>_<suffix>.csv
      --last-n N        Use only the last N frames (skip equilibration)
      --ncore N         Parallel threads                            [all cores]
      --atom-style S    LAMMPS data input: atomic | charge | full (required)
      --units U         Dump vx/fx units: real | metal (required if present)
      --modifier E,E    Elements counted for coordination only: no bridging
                        count, no part in ligand classification. Give each a
                        cutoff too
      --qn E,E          Formers reported as a Qn speciation. REPLACES the
                        default list (B,P,Si); every other former is described
                        by its coordination number
      --export-traj [FMT]
                        Also write the classified trajectory, one file per
                        input: <input stem>_types[_<suffix>].<ext>
                        FMT is lammpstrj (default) or extxyz

Output:
  network_<table>[_<suffix>].csv, six stacked tables with a `file` column:
  composition qn qn_partner ligand_type coordination linkage
  (the two qn tables only when some former is a Qn element)
  With auto, each input's cutoffs and g(min) are in the [inputs] block

Examples:
  ferro net -i traj.lammpstrj --P-O=2.4
  ferro net -i traj.lammpstrj --P-O=2.4 --Al-O=2.4 --Zn-O=2.6 --modifier Zn
  ferro net -i 'runs/*/prod.lammpstrj' --P-O=2.4 -o scan --export-traj
  ferro net -i 'runs/*/prod.lammpstrj' --P-O=auto --Al-O=auto -o scan

Full documentation:  ferro doc net";


#[cfg(test)]
mod tests {
    use super::*;

    fn args(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn test_split_pair_args_keeps_the_rest_for_clap() {
        let (pairs, rest) = split_pair_args(&args(&[
            "ferro", "net", "-i", "a.dump", "--P-O=2.4", "--modifier", "Zn", "--Zn-O=2.6",
        ]));
        assert_eq!(pairs, args(&["--P-O=2.4", "--Zn-O=2.6"]));
        assert_eq!(rest, args(&["ferro", "net", "-i", "a.dump", "--modifier", "Zn"]));
    }

    #[test]
    fn test_modifier_cutoffs_are_routed_out_of_the_former_table() {
        let p = build_params(&args(&["--P-O=2.4", "--Zn-O=2.6"]), Some("Zn"), None).unwrap().0;
        assert_eq!(p.formers(), vec!["P".to_string()]);
        assert_eq!(p.modifiers(), vec!["Zn".to_string()]);
    }

    #[test]
    fn test_modifier_without_its_cutoff_is_rejected() {
        // 静默通过的话 Zn 会被当成形成子,氧的分类整体错位
        let err = build_params(&args(&["--P-O=2.4"]), Some("Zn"), None).unwrap_err();
        assert!(err.to_string().contains("--modifier names Zn"), "{err}");
    }

    #[test]
    fn test_qn_replaces_the_default_list_rather_than_adding_to_it() {
        let p = build_params(&args(&["--P-O=2.4", "--Al-O=2.4"]), None, Some("Al")).unwrap().0;
        assert_eq!(p.qn_formers(), vec!["Al".to_string()], "P 必须被替换掉,不是叠加");
        // 不给 --qn 时走默认表:P 有 Qn,Al 没有
        let d = build_params(&args(&["--P-O=2.4", "--Al-O=2.4"]), None, None).unwrap().0;
        assert_eq!(d.qn_formers(), vec!["P".to_string()]);
    }

    #[test]
    fn test_qn_rejects_a_non_former_and_a_modifier() {
        // 点名了不是形成子的元素:静默忽略会让人以为报了 Qn 而其实没有
        let err = build_params(&args(&["--P-O=2.4"]), None, Some("Ti")).unwrap_err();
        assert!(err.to_string().contains("not a former"), "{err}");
        // 点名了修饰子:修饰子没有桥接数,给它 Qn 是自相矛盾
        let err = build_params(&args(&["--P-O=2.4", "--Zn-O=2.6"]), Some("Zn"), Some("Zn"))
            .unwrap_err();
        assert!(err.to_string().contains("--modifier already claims it"), "{err}");
    }

    #[test]
    fn test_bad_cutoffs_are_rejected() {
        assert!(parse_pairs(&args(&["--P-O=abc"])).is_err());
        assert!(parse_pairs(&args(&["--P-O=-1.0"])).is_err());
        assert!(parse_pairs(&args(&["--PO=2.4"])).is_err());
        assert!(parse_pairs(&args(&["--P-O=AUTO"])).is_err(), "只认小写 auto");
        assert_eq!(parse_pairs(&args(&["--P-O=auto"])).unwrap()
                   [&("P".to_string(), "O".to_string())], None);
    }

    #[test]
    fn test_auto_pairs_are_listed_and_hold_nan_until_resolved() {
        let (p, auto) = build_params(&args(&["--P-O=auto", "--Al-O=2.4", "--Zn-O=auto"]),
                                     Some("Zn"), None).unwrap();
        assert_eq!(auto, vec![("P".to_string(), "O".to_string()),
                              ("Zn".to_string(), "O".to_string())]);
        assert!(p.cutoffs[&("P".to_string(), "O".to_string())].is_nan());
        assert!(p.modifier_cutoffs[&("Zn".to_string(), "O".to_string())].is_nan());
        assert_eq!(p.cutoffs[&("Al".to_string(), "O".to_string())], 2.4);
    }

    /// 周期盒里一圈 P 各带 4 个 O（P-O 1.5 Å），其余空着：P-O 第一壳层外是干净空隙
    fn phosphate_box() -> Trajectory {
        use ferro_core::{Atom, Cell, Frame};
        use nalgebra::{Matrix3, Vector3};
        let cell = Cell::from_matrix(Matrix3::identity() * 20.0);
        let mut f = Frame::with_cell(cell, [true; 3]);
        let d = 1.5 / 3f64.sqrt();
        for c in [[3.0, 3.0, 3.0], [3.0, 13.0, 8.0], [13.0, 6.0, 15.0], [10.0, 15.0, 4.0]] {
            let c = Vector3::new(c[0], c[1], c[2]);
            f.atoms.push(Atom::new("P", c));
            for s in [[1.0, 1.0, 1.0], [1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0]] {
                f.atoms.push(Atom::new("O", c + Vector3::new(s[0], s[1], s[2]) * d));
            }
        }
        Trajectory { frames: vec![f], metadata: Default::default() }
    }

    // 审查 D-S4：--export-traj lammpstrj 原先恒以 real 写速度，无视 --units
    #[test]
    fn test_export_lammpstrj_follows_units() {
        let mut traj = phosphate_box();
        let n = traj.frames[0].atoms.len();
        let v = nalgebra::Vector3::new(0.001, 0.0, 0.0); // 内部 Å/fs = metal 的 1 Å/ps
        traj.frames[0].velocities = Some(vec![v; n]);
        let (p, _) = build_params(&args(&["--P-O=2.0"]), None, None).unwrap();
        let dir = std::env::temp_dir().join("ferro_net_export_units");
        let out = batch::Output { dir: Some(dir.clone()), label: None, suffix: None, mkdir: true };
        out.prepare().unwrap();

        export_labelled(&traj, &p, "box", &out, ExportFormat::Lammpstrj, Some(LammpsUnits::Metal)).unwrap();
        let path = dir.join("box_types.lammpstrj");
        let back = ferro_io::read_lammps_dump(&path, Some(LammpsUnits::Metal)).unwrap();
        let got = back.frames[0].velocities.as_ref().unwrap()[0];
        assert!((got - v).norm() < 1e-12, "按 metal 写出再按 metal 读回应不变，得到 {got:?}");

        let err = export_labelled(&traj, &p, "box", &out, ExportFormat::Lammpstrj, None).unwrap_err();
        assert!(format!("{err:#}").contains("--units"), "没给 --units 且有速度应报错：{err:#}");
    }

    #[test]
    fn test_resolve_auto_fills_the_cutoff_from_the_gap_behind_the_first_shell() {
        let (p, auto) = build_params(&args(&["--P-O=auto"]), None, None).unwrap();
        let (p, shells) = resolve_auto(&p, &auto, &phosphate_box()).unwrap();
        let rc = p.cutoffs[&("P".to_string(), "O".to_string())];
        assert!(rc > 1.6 && rc < 4.6, "截断应落在 1.5 Å 壳层之后的空隙里，得到 {rc}");
        let s = shells[0].expect("P-O 有第一壳层");
        assert_eq!(s.min_r, rc);
        assert_eq!(s.depth, 0.0);
    }

    #[test]
    fn test_resolve_auto_skips_a_pair_absent_from_the_input() {
        // 不含 Al 的样品照常分析：与给了定值截断时的行为一致，而不是整个文件失败
        let (p, auto) = build_params(&args(&["--P-O=auto", "--Al-O=auto"]), None, None).unwrap();
        let (p, shells) = resolve_auto(&p, &auto, &phosphate_box()).unwrap();
        assert!(shells[0].is_none(), "Al-O 不存在，应记为 None");
        assert_eq!(p.cutoffs[&("Al".to_string(), "O".to_string())], 0.0);
        assert!(p.cutoffs.values().all(|c| c.is_finite()), "NaN 不能留到 calc_network");
        calc_network(&phosphate_box(), &p).expect("缺 Al 的输入应能照常分析");
    }
}
