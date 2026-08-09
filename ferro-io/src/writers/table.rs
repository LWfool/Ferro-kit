//! Writer for [`ferro_core::Table`] — the single exit point for every analysis result.
//!
//! Replaces the nine per-analysis `write_*` functions that used to live in
//! `ferro-analysis`. Analysis code now only builds a `Table`; deciding what the bytes
//! look like happens here.

use ferro_core::Table;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{bail, Context, Result};

/// Output container for a [`Table`].
///
/// `Csv` is the only format the analysis binaries emit. The enum exists so adding a
/// second container is a new match arm rather than a new call path through the CLI.
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
pub fn write_table(table: &Table, path: &str, format: TableFormat) -> Result<()> {
    if let Err(e) = table.validate() {
        bail!("refusing to write {path}: {e}");
    }
    match format {
        TableFormat::Csv => write_csv(table, path),
    }
}

fn write_csv(table: &Table, path: &str) -> Result<()> {
    let file = File::create(path).context(format!("cannot create {path}"))?;
    let mut w = BufWriter::new(file);

    for line in &table.meta {
        if line.is_empty() {
            writeln!(w, "#")?;
        } else {
            writeln!(w, "# {line}")?;
        }
    }

    writeln!(w, "{}", table.names().join(","))?;

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
/// in a file name or a failure message cannot shift every following column.
fn escape(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::Table;

    fn tmp(name: &str) -> String {
        std::env::temp_dir().join(name).to_string_lossy().into_owned()
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
}
