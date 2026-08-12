//! `Table` — the neutral carrier for analysis results on their way to disk.
//!
//! `ferro-analysis` produces it (`XxxResult::to_tables()`), `ferro-io` consumes it
//! (`write_table`). Neither crate needs to name a type owned by the other, so the
//! layering rule in `CLAUDE.md` ("middle layers must not depend on each other") holds
//! while results still reach the filesystem. See the shared-type criterion there:
//! a type belongs here iff two or more middle layers need to name it.
//!
//! Shape conventions live with the producers, not here. A `Table` is just an ordered
//! list of equal-length named columns plus a free-form metadata block.

use std::fmt::Write as _;

/// One column of a [`Table`].
///
/// New variants (integers, booleans, dates) can be added here; both the comment-block
/// renderer below and the writers in `ferro-io` match exhaustively, so the compiler
/// will point at every place that needs updating.
#[derive(Clone, Debug, PartialEq)]
pub enum Column {
    Num(Vec<f64>),
    Text(Vec<String>),
}

impl Column {
    pub fn len(&self) -> usize {
        match self {
            Column::Num(v)  => v.len(),
            Column::Text(v) => v.len(),
        }
    }

    pub fn is_empty(&self) -> bool { self.len() == 0 }

    /// Renders cell `i` using the project-wide output convention.
    ///
    /// Numbers use `{:.6e}`; **NaN renders as an empty field** so that a missing value
    /// stays missing (pandas reads it back as `NaN`) instead of being silently invented.
    /// That distinction is the whole reason column-union tables are allowed at all.
    pub fn cell(&self, i: usize) -> String {
        match self {
            Column::Num(v) => match v.get(i) {
                Some(x) if x.is_nan() => String::new(),
                Some(x) => format!("{x:.6e}"),
                None    => String::new(),
            },
            Column::Text(v) => v.get(i).cloned().unwrap_or_default(),
        }
    }
}

/// An ordered set of equal-length named columns, plus free-form metadata lines.
///
/// The metadata block is written above the data as `#`-prefixed comments. It is
/// deliberately *not* machine-readable: `pandas.read_csv(comment="#")` drops it.
/// Anything a downstream script must parse belongs in a column.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Table {
    pub meta: Vec<String>,
    pub cols: Vec<(String, Column)>,
}

impl Table {
    pub fn new() -> Self { Self::default() }

    /// Appends a metadata line (written later as `# <line>`).
    pub fn meta_line(&mut self, line: impl Into<String>) -> &mut Self {
        self.meta.push(line.into());
        self
    }

    pub fn push_num(&mut self, name: impl Into<String>, values: Vec<f64>) -> &mut Self {
        self.cols.push((name.into(), Column::Num(values)));
        self
    }

    pub fn push_text(&mut self, name: impl Into<String>, values: Vec<String>) -> &mut Self {
        self.cols.push((name.into(), Column::Text(values)));
        self
    }

    pub fn n_cols(&self) -> usize { self.cols.len() }

    /// Row count, taken from the first column (0 when the table has no columns).
    pub fn n_rows(&self) -> usize {
        self.cols.first().map(|(_, c)| c.len()).unwrap_or(0)
    }

    pub fn names(&self) -> Vec<&str> {
        self.cols.iter().map(|(n, _)| n.as_str()).collect()
    }

    pub fn column(&self, name: &str) -> Option<&Column> {
        self.cols.iter().find(|(n, _)| n == name).map(|(_, c)| c)
    }

    /// Checks that every column has the same length.
    ///
    /// Writers call this before touching the filesystem, so a ragged table fails loudly
    /// instead of producing a half-written file.
    pub fn validate(&self) -> Result<(), String> {
        let n = self.n_rows();
        for (name, col) in &self.cols {
            if col.len() != n {
                return Err(format!(
                    "column '{name}' has {} rows, expected {n} (from column '{}')",
                    col.len(),
                    self.cols[0].0
                ));
            }
        }
        Ok(())
    }

    /// Stacks per-input tables into one, prepending a `label` column.
    ///
    /// Columns are unioned in first-seen order. A part missing a column contributes
    /// empty cells there (`NaN` for numbers) — never a fabricated value. This is what
    /// lets trajectories with **different element sets** share one wide table: a
    /// Zn-P-O run simply has no `Al-O_sq` values, and that reads back as missing
    /// rather than as zero.
    ///
    /// Errors if two parts disagree on a column's kind, which would otherwise write
    /// text into a numeric column.
    pub fn concat_union(label: &str, parts: Vec<(String, Table)>) -> Result<Table, String> {
        let mut names: Vec<String> = Vec::new();
        let mut is_text: Vec<bool> = Vec::new();
        for (_, t) in &parts {
            for (name, col) in &t.cols {
                let text = matches!(col, Column::Text(_));
                match names.iter().position(|n| n == name) {
                    Some(i) if is_text[i] != text => {
                        return Err(format!("column '{name}' is numeric in one input and text in another"));
                    }
                    Some(_) => {}
                    None => {
                        names.push(name.clone());
                        is_text.push(text);
                    }
                }
            }
        }

        let mut labels: Vec<String> = Vec::new();
        let mut nums: Vec<Vec<f64>>    = vec![Vec::new(); names.len()];
        let mut texts: Vec<Vec<String>> = vec![Vec::new(); names.len()];

        for (tag, t) in &parts {
            t.validate()?;
            let rows = t.n_rows();
            labels.extend(std::iter::repeat_n(tag.clone(), rows));
            for (ci, name) in names.iter().enumerate() {
                match t.column(name) {
                    Some(Column::Num(v))  => nums[ci].extend_from_slice(v),
                    Some(Column::Text(v)) => texts[ci].extend_from_slice(v),
                    // 该输入没有这一列:补空,不补零
                    None if is_text[ci] => texts[ci].extend(std::iter::repeat_n(String::new(), rows)),
                    None                => nums[ci].extend(std::iter::repeat_n(f64::NAN, rows)),
                }
            }
        }

        let mut out = Table::new();
        // 表的自述属于「这类表」,不属于某个输入 —— 各部分的说明是同一段静态文本,
        // 取第一份即可;拼接会把同一段话重复 N 遍
        out.meta = parts.first().map(|(_, t)| t.meta.clone()).unwrap_or_default();
        out.push_text(label, labels);
        for (ci, name) in names.into_iter().enumerate() {
            if is_text[ci] {
                out.push_text(name, std::mem::take(&mut texts[ci]));
            } else {
                out.push_num(name, std::mem::take(&mut nums[ci]));
            }
        }
        Ok(out)
    }

    /// Renders the whole table as space-aligned plain-text lines, header included.
    ///
    /// Used to embed a small table (the per-input metadata summary) inside another
    /// table's comment block. Returns lines **without** the `#` prefix; the writer adds it.
    pub fn to_comment_lines(&self) -> Vec<String> {
        if self.cols.is_empty() {
            return Vec::new();
        }
        let n = self.n_rows();
        let cells: Vec<Vec<String>> = self
            .cols
            .iter()
            .map(|(name, col)| {
                let mut v = Vec::with_capacity(n + 1);
                v.push(name.clone());
                for i in 0..n { v.push(col.cell(i)); }
                v
            })
            .collect();

        let widths: Vec<usize> = cells
            .iter()
            .map(|c| c.iter().map(|s| s.chars().count()).max().unwrap_or(0))
            .collect();

        (0..=n)
            .map(|row| {
                let mut line = String::new();
                for (ci, col) in cells.iter().enumerate() {
                    if ci > 0 { line.push_str("  "); }
                    let _ = write!(line, "{:<width$}", col[row], width = widths[ci]);
                }
                line.trim_end().to_string()
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Table {
        let mut t = Table::new();
        t.meta_line("ferro test")
            .push_text("file", vec!["a".into(), "b".into()])
            .push_num("r", vec![1.0, 2.0]);
        t
    }

    #[test]
    fn test_shape_accessors() {
        let t = sample();
        assert_eq!(t.n_cols(), 2);
        assert_eq!(t.n_rows(), 2);
        assert_eq!(t.names(), vec!["file", "r"]);
        assert!(t.column("r").is_some());
        assert!(t.column("missing").is_none());
        assert!(t.validate().is_ok());
    }

    #[test]
    fn test_validate_catches_ragged_columns() {
        let mut t = sample();
        t.push_num("extra", vec![1.0]);
        let err = t.validate().unwrap_err();
        assert!(err.contains("'extra'"), "error should name the bad column: {err}");
    }

    #[test]
    fn test_nan_renders_as_empty_field() {
        // 缺失值必须留空,不能写成 NaN 字面量,更不能被填成 0
        let c = Column::Num(vec![f64::NAN, 1.5]);
        assert_eq!(c.cell(0), "");
        assert_eq!(c.cell(1), format!("{:.6e}", 1.5));
    }

    #[test]
    fn test_num_uses_six_digit_scientific() {
        let c = Column::Num(vec![1234.5678]);
        assert_eq!(c.cell(0), "1.234568e3");
    }

    #[test]
    fn test_out_of_range_index_is_empty_not_panic() {
        let c = Column::Text(vec!["x".into()]);
        assert_eq!(c.cell(5), "");
    }

    #[test]
    fn test_comment_lines_align_columns() {
        let mut t = Table::new();
        t.push_text("file", vec!["short".into(), "a_much_longer_name".into()])
            .push_text("status", vec!["ok".into(), "failed".into()]);
        let lines = t.to_comment_lines();
        assert_eq!(lines.len(), 3, "header + 2 rows");
        assert!(lines[0].starts_with("file"));
        // 头行与数据行的第二列起始位置一致
        let pos = |s: &str, pat: &str| s.find(pat).unwrap();
        assert_eq!(pos(&lines[0], "status"), pos(&lines[1], "ok"));
        assert_eq!(pos(&lines[1], "ok"), pos(&lines[2], "failed"));
    }

    fn part(tag: &str, cols: &[(&str, Vec<f64>)]) -> (String, Table) {
        let mut t = Table::new();
        for (n, v) in cols { t.push_num(*n, v.clone()); }
        (tag.to_string(), t)
    }

    #[test]
    fn test_concat_prepends_label_and_stacks_rows() {
        let out = Table::concat_union(
            "file",
            vec![part("a", &[("r", vec![1.0, 2.0])]), part("b", &[("r", vec![3.0])])],
        )
        .unwrap();
        assert_eq!(out.names(), vec!["file", "r"]);
        assert_eq!(out.n_rows(), 3);
        match out.column("file").unwrap() {
            Column::Text(v) => assert_eq!(v, &vec!["a".to_string(), "a".into(), "b".into()]),
            _ => panic!("label column must be text"),
        }
    }

    #[test]
    fn test_concat_unions_columns_and_leaves_missing_empty() {
        // Zn-P-O 与 Al-P-O 同批:各自缺对方的配对列
        let out = Table::concat_union(
            "file",
            vec![
                part("znpo", &[("q", vec![1.0]), ("Zn-O_sq", vec![0.5])]),
                part("alpo", &[("q", vec![1.0]), ("Al-O_sq", vec![0.7])]),
            ],
        )
        .unwrap();
        assert_eq!(out.names(), vec!["file", "q", "Zn-O_sq", "Al-O_sq"],
                   "列按首见顺序取并集");
        let zn = out.column("Zn-O_sq").unwrap();
        assert_eq!(zn.cell(0), format!("{:.6e}", 0.5));
        assert_eq!(zn.cell(1), "", "缺失必须留空,不能是 0");
    }

    #[test]
    fn test_concat_rejects_kind_conflict() {
        let mut a = Table::new();
        a.push_num("x", vec![1.0]);
        let mut b = Table::new();
        b.push_text("x", vec!["one".into()]);
        let err = Table::concat_union("file", vec![("a".into(), a), ("b".into(), b)]).unwrap_err();
        assert!(err.contains("'x'"), "{err}");
    }

    #[test]
    fn test_concat_of_single_part_is_the_n_equals_one_case() {
        // N=1 不是特例:结构与多文件完全一致,只是 file 列取值唯一
        let out = Table::concat_union("file", vec![part("only", &[("r", vec![1.0, 2.0])])]).unwrap();
        assert_eq!(out.names(), vec!["file", "r"]);
        assert_eq!(out.n_rows(), 2);
    }

    #[test]
    fn test_empty_table_has_no_comment_lines() {
        assert!(Table::new().to_comment_lines().is_empty());
        assert_eq!(Table::new().n_rows(), 0);
    }
}
