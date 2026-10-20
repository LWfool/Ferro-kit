//! Multi-input driving: index files, run one analysis per file, stack the results.
//!
//! Everything here is **generic over the result type**. `batch` knows about paths,
//! loops, failures and file names; it never names `GrResult` or any other analysis
//! type, and the analysis crates never learn that batching exists.
//!
//! There is exactly one code path: a single input is the `N = 1` case of `N`, not a
//! special mode. Dispatching on the number of inputs would make the shape of the
//! output depend on how many files a glob happened to match that day.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Result};
use ferro_core::{Table, Trajectory};
use ferro_io::write_table;
use ferro_io::TableFormat;

/// One input and the label it carries into the `file` column, `[inputs]` and the
/// file names of per-input products. Labels are unique within a batch, see
/// [`unique_labels`].
pub struct Input {
    pub path: PathBuf,
    pub label: String,
}

/// One input that could not be analysed. Collected rather than fatal.
pub struct Failure {
    pub path: PathBuf,
    pub label: String,
    pub reason: String,
    /// Position in the input list, so [`Summary::failed`] can put the row back where
    /// the screen's `[k/N]` showed it.
    pub index: usize,
}

/// The label one input would get on its own: the file stem, not the full path
/// (`runs/700K/prod.lammpstrj` → `prod`). Legends and `groupby("file")` both want
/// something short. Within a batch use [`unique_labels`], which renames repeats.
fn label_of(path: &Path) -> String {
    path.file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string_lossy().into_owned())
}

/// Labels for a whole batch: the file stem, made unique.
///
/// Inputs in different directories often share a name (`runs/700K/prod.lammpstrj`,
/// `runs/900K/prod.lammpstrj`); with the bare stem the second `map` cube overwrites
/// the first and the stacked csv's `file` column cannot tell them apart. The second
/// and later holders of a stem become `<stem>_2`, `<stem>_3`, … in input order, and
/// each rename is printed. Not `<parent>_<stem>`: that makes product names long
/// (user decision 2026-10-02).
pub fn unique_labels(inputs: &[PathBuf]) -> Vec<String> {
    let stems: Vec<String> = inputs.iter().map(|p| label_of(p)).collect();
    // 撞名按小写比：macOS / Windows 的文件系统不分大小写，Prod 与 prod 的产物会互相
    // 覆盖（审查 E-12）。显示仍用原样的 stem
    // 改名不能撞上别的输入原本就叫的名字（prod、prod、prod_2 三个输入）
    let mut used: HashSet<String> = stems.iter().map(|s| s.to_lowercase()).collect();
    let mut count: HashMap<String, usize> = HashMap::new();
    let mut out = Vec::with_capacity(stems.len());
    for (path, stem) in inputs.iter().zip(&stems) {
        let n = count.entry(stem.to_lowercase()).or_insert(0);
        *n += 1;
        if *n == 1 {
            out.push(stem.clone());
            continue;
        }
        let mut k = *n;
        let label = loop {
            let cand = format!("{stem}_{k}");
            if used.insert(cand.to_lowercase()) { break cand; }
            k += 1;
        };
        crate::report::note(path.display(),
                            format!("shares the name '{stem}' with an earlier input, labelled '{label}'"));
        out.push(label);
    }
    out
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
    expand_with(patterns, |p| p.is_file(), "FILE", "file")
}

/// [`expand_inputs`] for the commands that read trajectories (`traj`, `map`, `net`):
/// a DeePMD system directory is one trajectory there, read like any file
/// (`ferro_io::Format::Deepmd`). `dataset collect` keeps to files, where a
/// directory is not an input.
pub fn expand_trajectories(patterns: &[PathBuf]) -> Result<Vec<PathBuf>> {
    expand_with(patterns, |p| p.is_file() || ferro_io::is_deepmd_system(p), "FILE", "file")
}

/// The directory form of [`expand_inputs`], for commands whose input is a
/// directory rather than a file — a DeePMD system is a directory, not a file.
pub fn expand_dirs(patterns: &[PathBuf]) -> Result<Vec<PathBuf>> {
    expand_with(patterns, |p| p.is_dir(), "DIR", "directory")
}

fn expand_with(
    patterns: &[PathBuf],
    keep: fn(&std::path::Path) -> bool,
    arg: &str,
    noun: &str,
) -> Result<Vec<PathBuf>> {
    if patterns.is_empty() {
        bail!("no input given (-i {arg} [{arg} ...], glob patterns allowed)");
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
                    if keep(&entry) {
                        matched.push(entry);
                    }
                }
            }
            // 不是合法模式(例如路径里有未闭合的 `[`)时按字面路径处理
            Err(_) if keep(pat) => matched.push(pat.clone()),
            Err(e) => bail!("bad input pattern '{pat_str}': {e}"),
        }

        if matched.is_empty() {
            bail!("no {noun} matches '{pat_str}'");
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
        crate::report::note(crate::report::PROGRAM, format!("{duplicates} duplicate input(s) skipped"));
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
    f: impl Fn(&Input) -> Result<T>,
) -> (Vec<(Input, T)>, Vec<Failure>) {
    let mut ok = Vec::new();
    let mut failed = Vec::new();

    // 标签在读第一个文件之前一次定好，撞名的改名提示也在这里打
    let labels = unique_labels(inputs);
    for (i, (path, label)) in inputs.iter().zip(labels).enumerate() {
        println!("[{}/{}] {}", i + 1, inputs.len(), path.display());
        let input = Input { path: path.clone(), label };
        match f(&input) {
            Ok(v) => ok.push((input, v)),
            Err(e) => {
                // 跳过而非中止:跑一晚上的批量不该因为第 7 个文件丢掉前 6 个的结果。
                // On stderr, as dataset's are (review E2-8)
                crate::report::error(input.path.display(), format!("{e:#}; skipped"));
                failed.push(Failure {
                    path: input.path,
                    label: input.label,
                    reason: format!("{e:#}"),
                    index: i,
                });
            }
        }
    }
    (ok, failed)
}

/// What [`drive`] hands back: one result per input that worked, the inputs that did
/// not, and the prepared output location.
pub type Driven<T> = (Vec<(Input, T)>, Vec<Failure>, Output);

/// The half of the pipeline every trajectory command shares (`traj`, `map`, `net`):
/// check the inputs, build `-o`, start the thread pool, run `calc` on each file, and
/// refuse an empty result set.
///
/// It stops there on purpose. The second half — the `Summary` columns, the product
/// name and title — differs in every caller, and threading those through as more
/// closures would cost more to read than the lines it saves.
///
/// The output directory is built **before the first file is read**, because the label
/// goes into the file name and a bad selection should fail immediately rather than
/// after a long batch.
pub fn drive<T>(
    common: &crate::args::common::CommonArgs,
    label: Option<String>,
    calc: impl Fn(&Input, &Trajectory) -> Result<T>,
) -> Result<Driven<T>> {
    // 参数（含 --atom-style）在建目录、读第一个文件之前查完
    let inputs = common.inputs()?;
    let out = common.out(label);
    out.prepare()?;
    common.init_threads();
    println!("Inputs: {} file(s)", inputs.len());

    let (results, failures) = map_inputs(&inputs, |inp| calc(inp, &common.load(&inp.path)?));
    if results.is_empty() {
        return Err(anyhow!("every input failed; nothing to write"));
    }
    Ok((results, failures, out))
}

/// Stacks each input's tables into one table per table name, adding the `file` column.
///
/// Inputs contributing different columns (different element sets) are unioned, with
/// the gaps left empty — see [`ferro_core::Table::concat_union`].
pub fn stack<T>(
    results: &[(Input, T)],
    to_tables: impl Fn(&T) -> Result<Vec<(String, Table)>>,
) -> Result<Vec<(String, Table)>> {
    let mut order: Vec<String> = Vec::new();
    let mut groups: Vec<Vec<(String, Table)>> = Vec::new();

    for (input, res) in results {
        let label = &input.label;
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
        let merged = Table::concat_union("file", parts)?;
        out.push((name, merged));
    }
    Ok(out)
}

/// Characters allowed in a file-name label. Everything a chemical symbol or a site
/// label can legitimately contain, and nothing that can redirect a path.
pub(crate) fn label_char_ok(c: char) -> bool {
    c.is_ascii_alphanumeric() || matches!(c, '_' | '+' | '-')
}

/// Checks the `-s` batch tag, which becomes part of a file name, before any file is
/// read: `-s 'a/../../x'` used to write outside `-o`.
pub fn check_suffix(suffix: Option<&str>) -> Result<()> {
    if let Some(bad) = suffix.and_then(|s| s.chars().find(|c| !label_char_ok(*c))) {
        bail!(
            "invalid character {bad:?} in -s: the batch tag goes into the output file \
             name, so only letters, digits, '_', '+' and '-' are accepted"
        );
    }
    Ok(())
}

/// Joins the selected types into the label that goes in the file name.
///
/// `["P", "O"] -> "P-O"`, `[] -> "all"`. The parts are what the caller wrote, so
/// `-a P -b O` and `-a O -b P` land in different files — that is correct, `g(r)` is
/// symmetric but `CN(r)` is directed. Callers whose selection is an unordered set
/// (`--elements`) sort and dedup before calling, or the same data would be written
/// under two names.
///
/// Validated **before** any file is read: these strings become a path component, and
/// a `/` in `-a` would otherwise silently write somewhere else.
pub fn file_label<S: AsRef<str>>(parts: &[S]) -> Result<String> {
    if parts.is_empty() {
        return Ok("all".to_string());
    }
    for p in parts {
        let p = p.as_ref();
        if p.is_empty() {
            bail!("empty type selection: a selected element or label cannot be blank");
        }
        if let Some(bad) = p.chars().find(|c| !label_char_ok(*c)) {
            bail!(
                "invalid character {bad:?} in '{p}': a selected element or site label \
                 goes into the output file name, so only letters, digits, '_', '+' and \
                 '-' are accepted"
            );
        }
    }
    Ok(parts.iter().map(|p| p.as_ref()).collect::<Vec<_>>().join("-"))
}

/// Sorts and dedups an unordered selection before it becomes a label.
///
/// `--elements` is a set: `O,P` and `P,O` select the same atoms and must not produce
/// two files holding identical data.
pub fn set_label(parts: Option<&Vec<String>>) -> Result<String> {
    match parts {
        Some(v) => {
            let mut v: Vec<&str> = v.iter().map(|s| s.as_str()).collect();
            v.sort_unstable();
            v.dedup();
            file_label(&v)
        }
        None => Ok("all".to_string()),
    }
}

/// Where a batch's products go and what they are called.
///
/// Collected into one struct because the three fields always travel together and
/// `write_all` would otherwise take eight positional arguments.
#[derive(Clone, Debug, Default)]
pub struct Output {
    /// `-o`, the output directory; `None` = the current directory.
    pub dir: Option<PathBuf>,
    /// What was analysed (`P-O`, `O-P-O`, `all`). `None` for modes with no selection.
    pub label: Option<String>,
    /// `-s`, the batch tag.
    pub suffix: Option<String>,
    /// `--mkdir`: create `dir` without asking.
    pub mkdir: bool,
}

impl Output {
    /// Creates `-o` if it does not exist. Call once, **before the first input is
    /// read**, so a bad path fails alongside the other parameter errors rather than
    /// after a long analysis.
    pub fn prepare(&self) -> Result<()> {
        match &self.dir {
            Some(d) => crate::outpath::ensure_dir(d, self.mkdir),
            None => Ok(()),
        }
    }

    /// Places a finished file name inside `-o`.
    pub fn join(&self, name: &str) -> PathBuf {
        match &self.dir {
            Some(d) => d.join(name),
            None => PathBuf::from(name),
        }
    }
}

/// The one naming rule for every product:
/// `<mode>[_<table>][_<selection>][_<input>][_<suffix>].<ext>`.
///
/// `segments` are mode, table, selection and input in that order; empty ones are
/// skipped. The `-s` tag always comes last, so `ls *_<tag>.*` lists one batch, and
/// the selection precedes the input, so `ls gr_P-O_*` lists one pair across batches.
/// A per-input product (`map`, `net --export-traj`, `bader`) always carries its input
/// label, even for a single input: whether a name has that segment must not depend on
/// how many files a glob matched (user decision 2026-10-18).
pub fn product_name(segments: &[&str], suffix: Option<&str>, ext: &str) -> String {
    let mut parts: Vec<&str> = segments.iter().copied().filter(|s| !s.is_empty()).collect();
    if let Some(s) = suffix.filter(|s| !s.is_empty()) {
        parts.push(s);
    }
    format!("{}.{ext}", parts.join("_"))
}

/// A stacked table's path: `<-o dir>/<mode>[_<table>][_<selection>][_<suffix>].csv`.
///
/// The table name is dropped when it merely repeats the mode (the single-table case).
pub fn out_path(mode: &str, table: &str, out: &Output) -> PathBuf {
    let table = if table == mode { "" } else { table };
    let label = out.label.as_deref().unwrap_or("");
    out.join(&product_name(&[mode, table, label], out.suffix.as_deref(), "csv"))
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

/// Writes each stacked table, all sharing one comment block.
///
/// The shared block is **prepended to**, not substituted for, whatever the analysis
/// already put in `Table::meta`. A per-table description is the only place a column's
/// meaning can be stated at the point of use — help text scrolls away, a manual is
/// somewhere else, but the `#` header travels with the file. Overwriting here would
/// silently discard it.
pub fn write_all(
    mode: &str,
    title: &str,
    params: &[String],
    summary: &Table,
    tables: Vec<(String, Table)>,
    out: &Output,
) -> Result<()> {
    let meta = meta_block(title, params, summary);
    for (name, mut table) in tables {
        if table.meta.is_empty() {
            table.meta = meta.clone();
        } else {
            let own = std::mem::take(&mut table.meta);
            table.meta = meta.clone();
            table.meta.extend(own);
            table.meta.push("-".repeat(60));
        }
        let path = out_path(mode, &name, out);
        write_table(&table, &path, TableFormat::Csv)?;
        println!("{:<12} -> {}", name.to_uppercase(), path.display());
    }
    Ok(())
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

    /// One failed row whose label is given directly rather than derived from a path.
    ///
    /// `failed` labels rows with the file stem, which is right when the inputs are
    /// files. A DeePMD system is a directory, and its label is the path relative to
    /// `-i` — nested systems `a/md` and `b/md` share a stem and would collide.
    pub fn failed_one(&mut self, file: String, reason: String) {
        self.file.push(file);
        self.status.push(reason);
        for slot in self.extra.iter_mut() {
            slot.1.push("-".to_string());
        }
    }

    /// Adds the failed inputs, each at its place in the input order.
    ///
    /// Call after every [`Summary::ok`] row, which go in input order: inserting the
    /// failures by ascending index then reproduces the order of the screen's `[k/N]`
    /// (审查 E2-7: they used to all sit at the end).
    pub fn failed(&mut self, failures: &[Failure]) {
        let rows = self.file.len();
        for slot in self.extra.iter_mut() {
            // note() 只填了部分行时先补齐，插入位置才对得上
            slot.1.resize(rows, "-".to_string());
        }
        let mut sorted: Vec<&Failure> = failures.iter().collect();
        sorted.sort_by_key(|f| f.index);
        for f in sorted {
            let at = f.index.min(self.file.len());
            self.file.insert(at, f.label.clone());
            self.status.insert(at, f.reason.clone());
            for slot in self.extra.iter_mut() {
                slot.1.insert(at, "-".to_string());
            }
        }
    }

    pub fn into_table(self) -> Table {
        self.into_table_named("file")
    }

    /// [`Summary::into_table`] with the row-label column named.
    ///
    /// `file` is right for the commands whose inputs are files. `ferro dataset`
    /// takes system directories, and calling their column `file` would contradict
    /// the `system` column of the tables the same report carries.
    pub fn into_table_named(mut self, label: &str) -> Table {
        let rows = self.file.len();
        let mut t = Table::new();
        t.push_text(label, self.file);
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
    fn test_unique_labels_rename_repeats_with_a_counter() {
        let p = |v: &[&str]| v.iter().map(PathBuf::from).collect::<Vec<_>>();
        // 不撞名时就是 stem，与此前逐字节相同
        assert_eq!(unique_labels(&p(&["a/x.dump", "b/y.dump"])), ["x", "y"]);
        // 不同目录下的同名文件：第二个起按出现顺序加 _2、_3
        assert_eq!(unique_labels(&p(&["700K/prod.dump", "900K/prod.dump", "1100K/prod.dump"])),
                   ["prod", "prod_2", "prod_3"]);
        // 改名不能撞上别的输入原本的名字
        assert_eq!(unique_labels(&p(&["a/prod.dump", "prod_2.dump", "b/prod.dump"])),
                   ["prod", "prod_2", "prod_3"]);
        // 审查 E-12：macOS / Windows 的文件系统不分大小写，Prod 与 prod 的产物会互相覆盖
        assert_eq!(unique_labels(&p(&["a/Prod.dump", "b/prod.dump"])), ["Prod", "prod_2"]);
        assert_eq!(unique_labels(&p(&["a/prod.dump", "PROD_2.dump", "b/prod.dump"])),
                   ["prod", "PROD_2", "prod_3"]);
    }

    #[test]
    fn test_map_inputs_skips_failures_and_keeps_the_rest() {
        let inputs = vec![PathBuf::from("a"), PathBuf::from("b"), PathBuf::from("c")];
        let (ok, failed) = map_inputs(&inputs, |inp| {
            if inp.path == Path::new("b") { bail!("boom") } else { Ok(inp.label.clone()) }
        });
        assert_eq!(ok.len(), 2);
        assert_eq!(failed.len(), 1);
        assert_eq!(failed[0].label, "b");
        assert!(failed[0].reason.contains("boom"));
    }

    /// 审查 E2-7：失败行要落在它的输入位置上，与屏幕的 `[k/N]` 对得上，
    /// 而不是全堆在 `[inputs]` 末尾；note() 只填了成功行的列也要随之对齐
    #[test]
    fn test_summary_keeps_failures_in_input_order() {
        let inputs: Vec<PathBuf> = ["a", "b", "c", "d"].iter().map(PathBuf::from).collect();
        let (ok, failed) = map_inputs(&inputs, |inp| {
            if inp.label == "a" || inp.label == "c" { bail!("boom") } else { Ok(()) }
        });
        let mut s = Summary::new(&[]);
        for (inp, _) in &ok {
            s.ok(inp.label.clone(), 1, 1, &[]);
            s.note("species", format!("sp-{}", inp.label));
        }
        s.failed(&failed);
        let t = s.into_table();
        let col = |name: &str| match t.column(name) {
            Some(ferro_core::Column::Text(v)) => v.clone(),
            other => panic!("{name}: {other:?}"),
        };
        assert_eq!(col("file"), ["a", "b", "c", "d"]);
        assert_eq!(col("species"), ["-", "sp-b", "-", "sp-d"]);
        assert_eq!(col("status")[1], "ok");
        assert!(col("status")[2].contains("boom"));
    }

    #[test]
    fn test_stack_groups_by_table_name() {
        let input = |p: &str, l: &str| Input { path: PathBuf::from(p), label: l.into() };
        let results = vec![(input("a.dump", "a"), 1.0f64), (input("b.dump", "b"), 2.0f64)];
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

    fn out(label: Option<&str>, suffix: Option<&str>) -> Output {
        Output {
            dir: None,
            label: label.map(str::to_string),
            suffix: suffix.map(str::to_string),
            mkdir: false,
        }
    }

    /// One rule for every product; the `-s` tag is always last and empty segments
    /// leave no stray `_`.
    #[test]
    fn test_product_name_order() {
        assert_eq!(product_name(&["density", "prod"], None, "cube"), "density_prod.cube");
        assert_eq!(product_name(&["density", "prod"], Some("700K"), "cube"), "density_prod_700K.cube");
        assert_eq!(product_name(&["sdf", "", "Ob", "prod"], Some(""), "cube"), "sdf_Ob_prod.cube");
        assert_eq!(product_name(&["sdf", "fam1", "Ob", "prod"], Some("x"), "cube"), "sdf_fam1_Ob_prod_x.cube");
        assert_eq!(product_name(&["bader", "ACF", "CHGCAR"], None, "dat"), "bader_ACF_CHGCAR.dat");
    }

    #[test]
    fn test_out_path_naming() {
        let p = |o: Output| out_path("gr", "gr", &o).to_string_lossy().into_owned();
        assert_eq!(p(out(None, None)), "gr.csv");
        assert_eq!(p(out(None, Some("run1"))), "gr_run1.csv");
        assert_eq!(p(out(Some("P-O"), None)), "gr_P-O.csv");
        assert_eq!(p(out(Some("P-O"), Some("run1"))), "gr_P-O_run1.csv");
        // label 在 suffix 之前:`ls gr_P-O_*` 能列出同一对在各批次的结果
        assert_eq!(p(out(Some("all"), Some("run1"))), "gr_all_run1.csv");

        let q = |o: Output| out_path("network", "qn", &o).to_string_lossy().into_owned();
        assert_eq!(q(out(None, None)), "network_qn.csv", "net 无 label 段");
        assert_eq!(q(out(None, Some("run1"))), "network_qn_run1.csv");
    }

    #[test]
    fn test_out_path_honours_the_output_dir() {
        let o = Output {
            dir: Some(PathBuf::from("results/700K")),
            label: Some("P-O".into()),
            suffix: None,
            mkdir: false,
        };
        assert_eq!(out_path("gr", "gr", &o), PathBuf::from("results/700K/gr_P-O.csv"));
    }

    #[test]
    fn test_file_label_joins_and_defaults_to_all() {
        assert_eq!(file_label::<&str>(&[]).unwrap(), "all");
        assert_eq!(file_label(&["P", "O"]).unwrap(), "P-O");
        assert_eq!(file_label(&["O", "P", "O"]).unwrap(), "O-P-O");
        // 位点标签自带下划线,原样拼接不转义
        assert_eq!(file_label(&["P_3", "O_b"]).unwrap(), "P_3-O_b");
        // 顺序照写:gr 的 CN 有向,-a P -b O 与 -a O -b P 是两份不同的数据
        assert_eq!(file_label(&["O", "P"]).unwrap(), "O-P");
    }

    #[test]
    fn test_file_label_rejects_path_characters() {
        // 这个串会成为路径的一段,放行等于让 -a 决定写到哪里
        for bad in ["P/2", "../etc", "P O", "P.O", ""] {
            assert!(file_label(&[bad]).is_err(), "{bad:?} 应被拒绝");
        }
    }

    #[test]
    fn test_set_label_sorts_and_dedups() {
        let v = |s: &[&str]| Some(s.iter().map(|x| x.to_string()).collect::<Vec<_>>());
        // --elements 是集合,两种写法必须落到同一个文件名
        assert_eq!(set_label(v(&["O", "P"]).as_ref()).unwrap(), "O-P");
        assert_eq!(set_label(v(&["P", "O"]).as_ref()).unwrap(), "O-P");
        assert_eq!(set_label(v(&["P", "O", "P"]).as_ref()).unwrap(), "O-P");
        assert_eq!(set_label(None).unwrap(), "all");
    }

    #[test]
    fn test_output_prepare_creates_missing_dir() {
        let base = std::env::temp_dir().join(format!("ferro_outdir_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&base);
        let o = Output { dir: Some(base.join("deep/nested")), mkdir: true, ..Default::default() };
        o.prepare().unwrap();
        assert!(o.dir.as_ref().unwrap().is_dir(), "缺失目录应被创建");
        o.prepare().unwrap(); // 已存在时是 no-op
        let _ = std::fs::remove_dir_all(&base);
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
