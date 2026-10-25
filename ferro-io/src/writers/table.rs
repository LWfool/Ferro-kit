//! Writer for [`ferro_core::Table`] — the single exit point for every analysis result.
//!
//! Replaces the nine per-analysis `write_*` functions that used to live in
//! `ferro-analysis`. Analysis code now only builds a `Table`; deciding what the bytes
//! look like happens here.

use std::path::Path;
use ferro_core::Table;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{bail, Context, Result};

/// Output container for a [`Table`].
///
/// `Csv` is the only format the analysis binaries emit. The enum exists so adding a
/// second container is a new match arm rather than a new call path through the CLI.
/// A one-variant enum was raised in the 2026-10-11 review; kept (user's call, 2026-10-19).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum TableFormat {
    #[default]
    Csv,
}

/// Writes `table` to `path`.
///
/// Layout: the metadata block first, one `# `-prefixed line each, then the header row
/// of column names, then the data. Missing numeric values render as empty fields
/// (see [`ferro_core::Column::cell`]), so `pandas.read_csv` reads them back as `NaN`.
///
/// The table is validated first — a ragged table fails before the file is created.
pub fn write_table(table: &Table, path: &Path, format: TableFormat) -> Result<()> {
    let path_ = path.display();
    if let Err(e) = table.validate() {
        bail!("refusing to write {path_}: {e}");
    }
    match format {
        TableFormat::Csv => write_csv(table, path),
    }
}

fn write_csv(table: &Table, path: &Path) -> Result<()> {
    let path_ = path.display();
    let file = File::create(path).context(format!("cannot create {path_}"))?;
    let mut w = BufWriter::new(file);

    // A meta line holding a newline is written as several, each behind its own `#`
    // (review R2-6)
    for line in table.meta.iter().flat_map(|l| l.split('\n')) {
        if line.is_empty() {
            writeln!(w, "#")?;
        } else {
            writeln!(w, "# {line}")?;
        }
    }

    let names: Vec<String> = table.names().into_iter().map(escape).collect();
    writeln!(w, "{}", names.join(","))?;

    for row in 0..table.n_rows() {
        let cells: Vec<String> = table.cols.iter().map(|(_, c)| escape(&c.cell(row))).collect();
        writeln!(w, "{}", cells.join(","))?;
    }

    w.flush()?;
    Ok(())
}

/// Quotes a field only when it would otherwise break the CSV structure.
///
/// Element symbols and site labels never need it; the escape exists so a stray comma
/// in a file name or a failure message cannot shift every following column. `#` is
/// quoted because the manual reads these files with `read_csv(comment="#")`, which
/// would cut the row there, and `\r` because it ends a line as `\n` does (review R2-5).
fn escape(s: &str) -> String {
    if s.contains([',', '"', '\n', '\r', '#']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::Table;

    fn tmp(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(name)
    }

    fn sample() -> Table {
        let mut t = Table::new();
        t.meta_line("ferro test")
            .meta_line("")
            .push_text("file", vec!["a".into(), "b".into()])
            .push_num("r", vec![1.0, f64::NAN]);
        t
    }

    #[test]
    fn test_csv_layout() {
        let p = tmp("ferro_table_layout.csv");
        write_table(&sample(), &p, TableFormat::Csv).unwrap();
        let s = std::fs::read_to_string(&p).unwrap();
        let lines: Vec<&str> = s.lines().collect();

        assert_eq!(lines[0], "# ferro test");
        assert_eq!(lines[1], "#", "空 meta 行只写 # 不留尾随空格");
        assert_eq!(lines[2], "file,r", "表头紧跟 meta 块");
        assert_eq!(lines[3], format!("a,{:.6e}", 1.0));
        assert_eq!(lines[4], "b,", "NaN 写成空字段");
        let _ = std::fs::remove_file(&p);
    }

    #[test]
    fn test_ragged_table_is_refused_before_creating_file() {
        let p = tmp("ferro_table_ragged.csv");
        let _ = std::fs::remove_file(&p);
        let mut t = sample();
        t.push_num("bad", vec![1.0]);
        assert!(write_table(&t, &p, TableFormat::Csv).is_err());
        assert!(!std::path::Path::new(&p).exists(), "失败时不应留下半个文件");
    }

    #[test]
    fn test_comma_in_field_is_quoted() {
        let p = tmp("ferro_table_quote.csv");
        let mut t = Table::new();
        t.push_text("status", vec!["failed: pair 'P-O' missing, skipped".into()]);
        write_table(&t, &p, TableFormat::Csv).unwrap();
        let s = std::fs::read_to_string(&p).unwrap();
        assert!(s.lines().nth(1).unwrap().starts_with('"'));
        let _ = std::fs::remove_file(&p);
    }

    /// The manual reads these files with `read_csv(comment="#")`: an unquoted `#` cut the
    /// rest of the row off, and a bare `\r` ended it. A meta line holding a newline left its
    /// second half outside the comment block; a column name was never escaped at all.
    #[test]
    fn test_hash_carriage_return_and_newlines_cannot_break_the_layout() {
        let p = tmp("ferro_table_hash.csv");
        let mut t = Table::new();
        t.meta_line("first\nsecond")
            .push_text("a,b", vec!["run#2".into()])
            .push_text("c", vec!["x\ry".into()]);
        write_table(&t, &p, TableFormat::Csv).unwrap();
        let s = std::fs::read_to_string(&p).unwrap();
        assert_eq!(s, "# first\n# second\n\"a,b\",c\n\"run#2\",\"x\ry\"\n");
        let _ = std::fs::remove_file(&p);
    }
}
