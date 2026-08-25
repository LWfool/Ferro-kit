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

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::{Args, Subcommand};

use crate::batch::expand_inputs;
use ferro_io::{read_cp2k_out_with_stats, write_deepmd_npy, Cp2kOutStats};

#[derive(Subcommand, Debug)]
pub enum DatasetCmd {
    /// Extract AIMD output into DeePMD system directories
    Collect(CollectCmd),
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
    }
}

pub fn print_help(cmd: &DatasetCmd) {
    match cmd {
        DatasetCmd::Collect(_) => crate::help::print_dataset_collect(),
    }
}

/// Returns the number of inputs that failed, for the process exit code.
pub fn run(cmd: &DatasetCmd) -> Result<usize> {
    match cmd {
        DatasetCmd::Collect(c) => run_collect(c),
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
    fn a_remaining_collision_is_an_error() {
        let inputs = vec![PathBuf::from("r/total.out"), PathBuf::from("x/r/total.out")];
        assert!(system_names(&inputs).is_err());
    }
}
