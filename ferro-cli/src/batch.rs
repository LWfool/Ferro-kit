//! Multi-input driving: index files, run one analysis per file, stack the results.
//!
//! Everything here is **generic over the result type**. `batch` knows about paths,
//! loops, failures and file names; it never names `GrResult` or any other analysis
//! type, and the analysis crates never learn that batching exists.
//!
//! There is exactly one code path: a single input is the `N = 1` case of `N`, not a
//! special mode. Dispatching on the number of inputs would make the shape of the
//! output depend on how many files a glob happened to match that day.

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use anyhow::{bail, Result};
use ferro_core::Table;
use ferro_io::write_table;
use ferro_io::writers::TableFormat;

/// One input that could not be analysed. Collected rather than fatal.
pub struct Failure {
    pub path: PathBuf,
    pub reason: String,
}

/// The label used for an input in the `file` column and in messages.
///
/// The file stem, not the full path: `runs/700K/prod.lammpstrj` → `prod`. Legends and
/// `groupby("file")` both want something short.
pub fn label_of(path: &Path) -> String {
    path.file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string_lossy().into_owned())
}

/// Expands glob patterns and literal paths into an ordered, deduplicated file list.
///
/// Patterns keep the order they were given on the command line; matches within one
/// pattern are sorted lexicographically. `{a,b}` braces are **not** supported — quote
/// them for the shell instead, since shell-expanded paths pass through here untouched.
///
/// A pattern matching nothing is an error. Silently analysing zero files is the
/// hardest failure mode to notice, especially inside a script.
pub fn expand_inputs(patterns: &[PathBuf]) -> Result<Vec<PathBuf>> {
    if patterns.is_empty() {
        bail!("no input given (-i FILE [FILE ...], glob patterns allowed)");
    }

    let mut out: Vec<PathBuf> = Vec::new();
    let mut seen: HashSet<PathBuf> = HashSet::new();
    let mut duplicates = 0usize;

    for pat in patterns {
        let pat_str = pat.to_string_lossy();
        let mut matched: Vec<PathBuf> = Vec::new();

        match glob::glob(&pat_str) {
            Ok(paths) => {
                for entry in paths.flatten() {
                    if entry.is_file() {
                        matched.push(entry);
                    }
                }
            }
            // 不是合法模式(例如路径里有未闭合的 `[`)时按字面路径处理
            Err(_) if pat.is_file() => matched.push(pat.clone()),
            Err(e) => bail!("bad input pattern '{pat_str}': {e}"),
        }

        if matched.is_empty() {
            bail!("no file matches '{pat_str}'");
        }
        matched.sort();

        for m in matched {
            let key = std::fs::canonicalize(&m).unwrap_or_else(|_| m.clone());
            if seen.insert(key) {
                out.push(m);
            } else {
                duplicates += 1;
            }
        }
    }

    if duplicates > 0 {
        println!("Note  : {duplicates} duplicate input(s) skipped");
    }
    Ok(out)
}

/// Runs `f` over every input, keeping successes and collecting failures.
///
/// Serial by design: one trajectory at a time is the memory peak, and each analysis
/// already parallelises internally over frames (`--ncore` keeps its meaning). Results
/// are small enough to all stay resident; trajectories are dropped as we go.
pub fn map_inputs<T>(
    inputs: &[PathBuf],
    f: impl Fn(&Path) -> Result<T>,
) -> (Vec<(PathBuf, T)>, Vec<Failure>) {
    let mut ok = Vec::new();
    let mut failed = Vec::new();

    for (i, path) in inputs.iter().enumerate() {
        println!("[{}/{}] {}", i + 1, inputs.len(), path.display());
        match f(path) {
            Ok(v) => ok.push((path.clone(), v)),
            Err(e) => {
                // 跳过而非中止:跑一晚上的批量不该因为第 7 个文件丢掉前 6 个的结果
                println!("        skipped: {e:#}");
                failed.push(Failure { path: path.clone(), reason: format!("{e:#}") });
            }
        }
    }
    (ok, failed)
}

/// Stacks each input's tables into one table per table name, adding the `file` column.
///
/// Inputs contributing different columns (different element sets) are unioned, with
/// the gaps left empty — see [`ferro_core::Table::concat_union`].
pub fn stack<T>(
    results: &[(PathBuf, T)],
    to_tables: impl Fn(&T) -> Result<Vec<(String, Table)>>,
) -> Result<Vec<(String, Table)>> {
    let mut order: Vec<String> = Vec::new();
    let mut groups: Vec<Vec<(String, Table)>> = Vec::new();

    for (path, res) in results {
        let label = label_of(path);
        for (name, table) in to_tables(res)? {
            match order.iter().position(|n| *n == name) {
                Some(i) => groups[i].push((label.clone(), table)),
                None => {
                    order.push(name);
                    groups.push(vec![(label.clone(), table)]);
                }
            }
        }
    }

    let mut out = Vec::with_capacity(order.len());
    for (name, parts) in order.into_iter().zip(groups) {
        let merged = Table::concat_union("file", parts).map_err(|e| anyhow::anyhow!(e))?;
        out.push((name, merged));
    }
    Ok(out)
}

/// Output file name: `<mode>[_<table>]_<suffix>.csv`.
///
/// `-o` supplies the suffix, so one batch's products sort together under `ls <mode>_*`.
/// The table name is dropped when it merely repeats the mode (the single-table case).
pub fn out_path(mode: &str, table: &str, suffix: Option<&str>) -> String {
    let stem = if table == mode { mode.to_string() } else { format!("{mode}_{table}") };
    match suffix {
        Some(s) if !s.is_empty() => format!("{stem}_{s}.csv"),
        _ => format!("{stem}.csv"),
    }
}

/// Builds the `#` comment block written above the data.
///
/// Version, title, the analysis parameters, then a per-input summary table listing
/// every input — successes and failures alike, so "which files went in" has a single
/// source of truth. Not machine-readable by design: `read_csv(comment="#")` drops it.
pub fn meta_block(title: &str, params: &[String], summary: &Table) -> Vec<String> {
    let mut v = vec![
        format!("ferro v{}", env!("CARGO_PKG_VERSION")),
        title.to_string(),
        "-".repeat(60),
    ];
    v.extend(params.iter().cloned());
    if summary.n_rows() > 0 {
        v.push("[inputs]".to_string());
        v.extend(summary.to_comment_lines());
    }
    v.push("-".repeat(60));
    v
}

/// Writes each stacked table, all sharing one comment block. Returns the first path
/// (the one a plot is named after).
pub fn write_all(
    mode: &str,
    title: &str,
    params: &[String],
    summary: &Table,
    tables: Vec<(String, Table)>,
    suffix: Option<&str>,
) -> Result<String> {
    let meta = meta_block(title, params, summary);
    let mut first = String::new();
    for (name, mut table) in tables {
        table.meta = meta.clone();
        let path = out_path(mode, &name, suffix);
        write_table(&table, &path, TableFormat::Csv)?;
        println!("{:<8}-> {path}", name.to_uppercase());
        if first.is_empty() {
            first = path;
        }
    }
    Ok(first)
}

/// One row per input: statistics for the ones that worked, the reason for the ones
/// that did not. Lands in the `[inputs]` block so "which files went in" has a single
/// source of truth.
///
/// Values are stored pre-formatted as text. This block is read by humans, not by
/// `read_csv`, and `{:.6e}` turns "5 frames" into `5.000000e0`.
pub struct Summary {
    file: Vec<String>,
    status: Vec<String>,
    extra: Vec<(String, Vec<String>)>,
}

/// Compact rendering for the summary block: integers stay integers, everything else
/// gets four significant digits. Non-finite values render as `-`.
pub fn fmt_compact(v: f64) -> String {
    if !v.is_finite() {
        "-".to_string()
    } else if v.fract() == 0.0 && v.abs() < 1e9 {
        format!("{v:.0}")
    } else if v.abs() >= 1e5 || (v != 0.0 && v.abs() < 1e-3) {
        format!("{v:.4e}")
    } else {
        format!("{v:.4}")
    }
}

impl Summary {
    /// `extra_names` are the mode-specific columns, in order, appended after `atoms`.
    pub fn new(extra_names: &[&str]) -> Self {
        let mut extra: Vec<(String, Vec<String>)> =
            vec![("frames".into(), Vec::new()), ("atoms".into(), Vec::new())];
        extra.extend(extra_names.iter().map(|n| (n.to_string(), Vec::new())));
        Summary { file: Vec::new(), status: Vec::new(), extra }
    }

    pub fn ok(&mut self, file: String, frames: usize, atoms: usize, extra: &[f64]) {
        self.file.push(file);
        self.status.push("ok".into());
        self.extra[0].1.push(frames.to_string());
        self.extra[1].1.push(atoms.to_string());
        for (slot, v) in self.extra[2..].iter_mut().zip(extra) {
            slot.1.push(fmt_compact(*v));
        }
    }

    /// Adds a free-text column value to the row just pushed by [`Summary::ok`].
    ///
    /// For things that are per-input but not numeric — a composition, a clamped
    /// setting — which must not sit in the shared parameter block pretending to be
    /// global.
    pub fn note(&mut self, name: &str, value: String) {
        let rows = self.file.len();
        let slot = match self.extra.iter_mut().find(|(n, _)| n == name) {
            Some(s) => s,
            None => {
                self.extra.push((name.to_string(), vec!["-".to_string(); rows - 1]));
                self.extra.last_mut().unwrap()
            }
        };
        while slot.1.len() < rows - 1 {
            slot.1.push("-".to_string());
        }
        slot.1.push(value);
    }

    pub fn failed(&mut self, failures: &[Failure]) {
        for f in failures {
            self.file.push(label_of(&f.path));
            self.status.push(f.reason.clone());
            for slot in self.extra.iter_mut() {
                slot.1.push("-".to_string());
            }
        }
    }

    pub fn into_table(mut self) -> Table {
        let rows = self.file.len();
        let mut t = Table::new();
        t.push_text("file", self.file);
        for (name, mut values) in std::mem::take(&mut self.extra) {
            // note() 只填了部分行时补齐,否则 validate 会拒绝这张表
            while values.len() < rows {
                values.push("-".to_string());
            }
            t.push_text(name, values);
        }
        t.push_text("status", self.status);
        t
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../tests")
    }

    #[test]
    fn test_expand_rejects_zero_match() {
        let err = expand_inputs(&[fixture_dir().join("no_such_*.lammpstrj")]).unwrap_err();
        assert!(err.to_string().contains("no file matches"), "{err}");
    }

    #[test]
    fn test_expand_rejects_empty_input_list() {
        assert!(expand_inputs(&[]).is_err());
    }

    #[test]
    fn test_expand_glob_sorted_and_deduplicated() {
        let dir = fixture_dir();
        let pattern = dir.join("*.lammpstrj");
        let by_glob = expand_inputs(std::slice::from_ref(&pattern)).unwrap();
        assert!(by_glob.len() >= 2, "fixtures should provide several trajectories");
        let mut sorted = by_glob.clone();
        sorted.sort();
        assert_eq!(by_glob, sorted, "matches within one pattern are lexicographic");

        // 同一个文件经模式与字面路径各来一次 → 只留一份
        let dup = expand_inputs(&[pattern, by_glob[0].clone()]).unwrap();
        assert_eq!(dup.len(), by_glob.len(), "duplicates must be dropped");
    }

    #[test]
    fn test_label_is_the_stem() {
        assert_eq!(label_of(Path::new("runs/700K/prod.lammpstrj")), "prod");
    }

    #[test]
    fn test_map_inputs_skips_failures_and_keeps_the_rest() {
        let inputs = vec![PathBuf::from("a"), PathBuf::from("b"), PathBuf::from("c")];
        let (ok, failed) = map_inputs(&inputs, |p| {
            if p == Path::new("b") { bail!("boom") } else { Ok(p.to_string_lossy().into_owned()) }
        });
        assert_eq!(ok.len(), 2);
        assert_eq!(failed.len(), 1);
        assert_eq!(label_of(&failed[0].path), "b");
        assert!(failed[0].reason.contains("boom"));
    }

    #[test]
    fn test_stack_groups_by_table_name() {
        let results = vec![
            (PathBuf::from("a.dump"), 1.0f64),
            (PathBuf::from("b.dump"), 2.0f64),
        ];
        let stacked = stack(&results, |v| {
            let mut t = Table::new();
            t.push_num("x", vec![*v]);
            let mut u = Table::new();
            u.push_num("y", vec![*v * 10.0]);
            Ok(vec![("main".into(), t), ("aux".into(), u)])
        })
        .unwrap();

        assert_eq!(stacked.len(), 2);
        assert_eq!(stacked[0].0, "main", "表名保持首见顺序");
        assert_eq!(stacked[0].1.names(), vec!["file", "x"]);
        assert_eq!(stacked[0].1.n_rows(), 2, "两个输入各贡献一行");
    }

    #[test]
    fn test_out_path_naming() {
        assert_eq!(out_path("gr", "gr", Some("run1")), "gr_run1.csv");
        assert_eq!(out_path("gr", "gr", None), "gr.csv");
        assert_eq!(out_path("network", "qn", Some("run1")), "network_qn_run1.csv");
        assert_eq!(out_path("network", "qn", None), "network_qn.csv");
    }

    #[test]
    fn test_meta_block_lists_inputs() {
        let mut summary = Table::new();
        summary.push_text("file", vec!["a".into(), "b".into()])
            .push_text("status", vec!["ok".into(), "read error".into()]);
        let block = meta_block("g(r)", &["dr = 0.002".into()], &summary).join("\n");
        assert!(block.contains("ferro v"));
        assert!(block.contains("dr = 0.002"));
        assert!(block.contains("[inputs]"));
        assert!(block.contains("read error"), "失败的输入也要留在清单里");
    }
}
