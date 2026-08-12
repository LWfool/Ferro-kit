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

use anyhow::{anyhow, bail, Result};
use clap::{Args, ValueEnum};
use ferro_analysis::{calc_network, NetworkResult};
use ferro_core::{Trajectory, TypeParams};
use ferro_io::{write_extxyz, write_lammps_dump};
use ferro_structure::{apply_type_labels, classify_trajectory, fold_labels};
use std::collections::BTreeMap;
use std::path::Path;

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
    let params = build_params(pair_args, cmd.modifier.as_deref(), cmd.qn.as_deref())?;
    if params.cutoffs.is_empty() {
        bail!("Every cutoff names a modifier element; at least one former is required");
    }

    let inputs = batch::expand_inputs(&cmd.common.input)?;
    cmd.common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());
    print_label_scheme(&params);

    let (results, failures) = batch::map_inputs(&inputs, |p| {
        let traj = cmd.common.load(p)?;
        let result = calc_network(&traj, &params).ok_or_else(|| {
            anyhow!("no usable frame (every frame is missing a cell; PBC required)")
        })?;
        if let Some(fmt) = cmd.export_traj {
            export_labelled(&traj, &params, p, cmd.common.suffix(), fmt)?;
        }
        Ok(result)
    });
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }

    let tables = batch::stack(&results, |r: &NetworkResult| Ok(r.to_tables()))?;
    note_missing_qn_tables(&params);

    let mut summary = Summary::new(&[]);
    for (path, r) in &results {
        summary.ok(batch::label_of(path), r.n_frames, r.n_atoms, &[]);
        summary.note("mean_qn", fmt_means(&r.mean_qn, &r.qn_dist));
        summary.note("mean_cn", fmt_means(&r.mean_cn, &r.cn_dist));
    }
    summary.failed(&failures);

    batch::write_all(
        "network",
        "Glass Network Topology — Qn speciation, ligand types, coordination numbers",
        &results[0].1.meta_lines(),
        &summary.into_table(),
        tables,
        cmd.common.suffix(),
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
    println!("  {:<10} <elem>_<Qn>   digit = bridging ligands (Qn)", join(qn));
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
    input: &Path,
    suffix: Option<&str>,
    fmt: ExportFormat,
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

    let stem = batch::label_of(input);
    let ext = fmt.ext();
    let path = match suffix {
        Some(s) if !s.is_empty() => format!("{stem}_types_{s}.{ext}"),
        _ => format!("{stem}_types.{ext}"),
    };
    match fmt {
        ExportFormat::Lammpstrj => write_lammps_dump(&out_traj, &path, ferro_io::LammpsUnits::Real)?,
        ExportFormat::Extxyz => write_extxyz(&out_traj, &path)?,
    }
    if skipped > 0 {
        println!(
            "        note: {skipped} label(s) not of the form <element>_<suffix>; \
             wrote the element instead"
        );
    }
    println!("        traj -> {path}");
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

fn build_params(
    pair_args: &[String],
    modifier: Option<&str>,
    qn: Option<&str>,
) -> Result<TypeParams> {
    let modifier_elems: std::collections::HashSet<String> =
        split_elems(modifier).into_iter().collect();

    let mut cutoffs = BTreeMap::new();
    let mut modifier_cutoffs = BTreeMap::new();
    for ((elem, ligand), cutoff) in parse_pairs(pair_args)? {
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
    Ok(params)
}

fn parse_pairs(pair_args: &[String]) -> Result<BTreeMap<(String, String), f64>> {
    let mut map = BTreeMap::new();
    for arg in pair_args {
        let inner = arg.trim_start_matches('-');
        let (pair, cutoff_str) = inner.split_once('=')
            .ok_or_else(|| anyhow!("Invalid pair argument (missing '='): {arg}"))?;
        let (former, ligand) = pair.split_once('-')
            .ok_or_else(|| anyhow!("Invalid pair argument (missing '-'): {arg}"))?;
        let cutoff: f64 = cutoff_str.parse()
            .map_err(|_| anyhow!("Invalid cutoff value in '{arg}'"))?;
        if cutoff <= 0.0 { bail!("Cutoff must be positive, got {cutoff} in '{arg}'"); }
        map.insert((former.to_string(), ligand.to_string()), cutoff);
    }
    Ok(map)
}

const HELP_EXTRA: &str = "\
ferro net — Glass network topology

USAGE:
  ferro net -i <FILE>... --<Former>-<Ligand>=<cutoff> [OPTIONS]

PAIR ARGUMENTS (required, at least one):
  --P-O=2.4            P former, O ligand, cutoff 2.4 Å
  --Al-O=2.4 --Al-F=2.1
  The element pair lives in the flag name, so these are stripped from argv before
  clap parses; everything else follows the usual flag rules.

OPTIONS:
  -i, --input  FILE...  Input trajectory files; glob patterns allowed (quote them)
  -o, --output SUFFIX   Output name suffix: network_<table>_<suffix>.csv
      --last-n N        Use only the last N frames (skip equilibration)
      --ncore N         Parallel threads (default: all cores)
      --metal-units     LAMMPS metal units; only affects --export-traj extxyz
      --modifier E,E    Elements counted for coordination only: no bridging count,
                        no part in ligand classification. Give each a cutoff too
      --qn E,E          Formers reported as a Qn speciation. REPLACES the default
                        list (B,P,Si); every other former is described by its
                        coordination number
      --export-traj [FMT]
                        Also write the classified trajectory, one file per input:
                        <input stem>_types[_<suffix>].<ext>
                        FMT is lammpstrj (default) or extxyz

LABELS:
  Printed with this run's elements filled in when the analysis starts, and
  repeated in the header of every file that has a label column.

  NOTE: an exported trajectory can be selected by element (-a/-b/-c) over any
  number of frames, but by label (-x/-y/-z) only from a SINGLE frame — g(r)
  requires a fixed particle count per type and labels change as the run evolves.
  Use --last-n 1, or select by element.

OUTPUT — stacked CSVs, each with a `file` column. Every file carries a `#` header
describing its own columns; `pandas.read_csv(comment='#')` drops it.

  network_qn.csv            Qn speciation — the plain distribution, readable as-is
  network_qn_partner.csv    the same, split by partner element: Q^n(mAl)
  network_ligand_type.csv   ligand speciation: free / non-bridging / bridging /
                            tricluster, with the formers each one joins
  network_coordination.csv  coordination number distribution, formers + modifiers
  network_average.csv       per-element mean Qn and mean coordination number
  network_linkage.csv       bridge connectivity, both ends and the ligand between

  The first two are omitted when no former is a Qn element; a reason is printed.

EXAMPLES:
  ferro net -i traj.lammpstrj --P-O=2.4
  ferro net -i traj.lammpstrj --P-O=2.4 --Al-O=2.4 --Zn-O=2.6 --modifier Zn
  ferro net -i 'runs/*/prod.lammpstrj' --P-O=2.4 -o scan
  ferro net -i traj.lammpstrj --P-O=2.4 --last-n 50 --export-traj
  ferro net -i traj.lammpstrj --Al-O=2.4 --Si-O=2.0 --qn Si,Al";


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
        let p = build_params(&args(&["--P-O=2.4", "--Zn-O=2.6"]), Some("Zn"), None).unwrap();
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
        let p = build_params(&args(&["--P-O=2.4", "--Al-O=2.4"]), None, Some("Al")).unwrap();
        assert_eq!(p.qn_formers(), vec!["Al".to_string()], "P 必须被替换掉,不是叠加");
        // 不给 --qn 时走默认表:P 有 Qn,Al 没有
        let d = build_params(&args(&["--P-O=2.4", "--Al-O=2.4"]), None, None).unwrap();
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
    }
}
