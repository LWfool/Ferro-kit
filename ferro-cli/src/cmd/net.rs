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
use clap::Args;
use ferro_analysis::{calc_network, NetworkResult};
use ferro_core::{Trajectory, TypeParams};
use ferro_io::write_lammps_dump;
use ferro_structure::{apply_type_labels, classify_trajectory};
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

    /// Also write the classified trajectory to <input>_types[_<suffix>].lammpstrj
    #[arg(long)]
    pub export_traj: bool,
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
    let params = build_params(pair_args, cmd.modifier.as_deref())?;
    if params.cutoffs.is_empty() {
        bail!("Every cutoff names a modifier element; at least one former is required");
    }

    let inputs = batch::expand_inputs(&cmd.common.input)?;
    cmd.common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());

    let (results, failures) = batch::map_inputs(&inputs, |p| {
        let traj = cmd.common.load(p)?;
        let result = calc_network(&traj, &params).ok_or_else(|| {
            anyhow!("no usable frame (every frame is missing a cell; PBC required)")
        })?;
        if cmd.export_traj {
            export_labelled(&traj, &params, p, cmd.common.suffix())?;
        }
        Ok(result)
    });
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }

    let tables = batch::stack(&results, |r: &NetworkResult| Ok(r.to_tables()))?;

    let mut summary = Summary::new(&[]);
    for (path, r) in &results {
        summary.ok(batch::label_of(path), r.n_frames, r.n_atoms, &[]);
        summary.note("mean_n_bridge", fmt_means(&r.mean_bridge, &r.bridge_dist));
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

/// Writes the classified trajectory as `<input stem>_types[_<suffix>].lammpstrj`.
///
/// The name carries the input stem because this is a **one product per input**
/// product, like `ferro map`'s cubes — a fixed output path would make the second
/// input overwrite the first.
fn export_labelled(
    traj: &Trajectory,
    params: &TypeParams,
    input: &Path,
    suffix: Option<&str>,
) -> Result<()> {
    let per_frame = classify_trajectory(traj, params);
    if per_frame.len() != traj.frames.len() {
        bail!(
            "cannot export: {} of {} frames have no cell",
            traj.frames.len() - per_frame.len(),
            traj.frames.len()
        );
    }

    let frames = traj.frames.iter().zip(&per_frame)
        .map(|(frame, types)| {
            let labels: Vec<String> = types.iter().map(|t| t.label()).collect();
            apply_type_labels(frame, &labels)
        })
        .collect();
    let labelled = Trajectory { frames, metadata: traj.metadata.clone() };

    let stem = batch::label_of(input);
    let path = match suffix {
        Some(s) if !s.is_empty() => format!("{stem}_types_{s}.lammpstrj"),
        _ => format!("{stem}_types.lammpstrj"),
    };
    write_lammps_dump(&labelled, &path, ferro_io::LammpsUnits::Real)?;
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

fn build_params(pair_args: &[String], modifier: Option<&str>) -> Result<TypeParams> {
    let modifier_elems: std::collections::HashSet<String> = modifier
        .map(|s| s.split(',').map(|e| e.trim().to_string()).collect())
        .unwrap_or_default();

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
    Ok(TypeParams { cutoffs, modifier_cutoffs })
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
      --metal-units     LAMMPS metal units; only affects velocities/forces, which
                        this analysis does not read
      --modifier E,E    Elements counted for coordination only: they take no part
                        in bridging counts or ligand classification.
                        Supply each one's cutoff too, e.g. --Zn-O=2.6
      --export-traj     Also write the classified trajectory, one file per input:
                        <input stem>_types[_<suffix>].lammpstrj

LABELS:
  Former      P_0 P_1 P_2 …    digit = number of bridging ligands (Qn for P)
  Free        O_f              bonded to no former
  Non-bridge  O_n              bonded to one former
  Bridge      O_b              bonded to two formers
  Tricluster  O_t              bonded to three or more
  Modifier    Zn               element symbol, no role suffix
  Every label splits at the first underscore into element + label, so an exported
  trajectory can be selected either way:
    ferro traj gr -i run_types.lammpstrj -a P -b O     (by element)
    ferro traj gr -i run_types.lammpstrj -x P_3 -y O_b (by label)

OUTPUT — stacked CSVs, each with a `file` column:
  network_bridge.csv  file, former, n_bridge, count, fraction, sd
  network_oxy.csv     file, type, former_a, former_b, count, fraction, sd
  network_cn.csv      file, element, cn, count, fraction, sd
  network_mean.csv    file, element, mean_n_bridge, mean_cn

  n_bridge is the number of bridging ligands. For P that is Qn; Al has no Qn in
  the literature, so the column is named for what it counts and Al is described
  by the cn table instead.

  The oxy table keeps the partner elements as data columns, not in the label:
  P-O-P and P-O-Al are both labelled O_b but are separate rows (former_a/_b).
  A tricluster (O_t) leaves them empty — its pairs are in the linkage table.

  fraction is the mean over frames of that frame's fraction; sd is the sample
  standard deviation (ddof=1) of the same quantity. Consecutive MD frames are
  correlated, so sd is a spread between snapshots, NOT a standard error.

EXAMPLES:
  ferro net -i traj.lammpstrj --P-O=2.4
  ferro net -i traj.lammpstrj --P-O=2.4 --Al-O=2.4 --Zn-O=2.6 --modifier Zn
  ferro net -i 'runs/*/prod.lammpstrj' --P-O=2.4 -o scan
  ferro net -i traj.lammpstrj --P-O=2.4 --last-n 50 --export-traj";

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
        let p = build_params(&args(&["--P-O=2.4", "--Zn-O=2.6"]), Some("Zn")).unwrap();
        assert_eq!(p.formers(), vec!["P".to_string()]);
        assert_eq!(p.modifiers(), vec!["Zn".to_string()]);
    }

    #[test]
    fn test_modifier_without_its_cutoff_is_rejected() {
        // 静默通过的话 Zn 会被当成形成子,氧的分类整体错位
        let err = build_params(&args(&["--P-O=2.4"]), Some("Zn")).unwrap_err();
        assert!(err.to_string().contains("--modifier names Zn"), "{err}");
    }

    #[test]
    fn test_bad_cutoffs_are_rejected() {
        assert!(parse_pairs(&args(&["--P-O=abc"])).is_err());
        assert!(parse_pairs(&args(&["--P-O=-1.0"])).is_err());
        assert!(parse_pairs(&args(&["--PO=2.4"])).is_err());
    }
}
