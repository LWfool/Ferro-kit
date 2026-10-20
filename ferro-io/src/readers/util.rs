//! Helpers shared by more than two readers.

use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, ensure, Context, Result};

/// The leading floats of `line`, stopping at the first token that is not a number.
///
/// **Stopping** is the point: header lines carry a trailing comment or label that must not
/// abort the parse.  `readers/vasp_outcar.rs` needs the opposite (skip non-numbers and keep
/// going, for the `in kB` label at the *start* of a line) and therefore keeps its own —
/// merging the two would silently change how both parse.
pub(crate) fn floats(line: &str, min: usize) -> Result<Vec<f64>> {
    let v: Vec<f64> = line
        .split_whitespace()
        .map_while(|s| s.parse::<f64>().ok())
        .collect();
    ensure!(v.len() >= min, "expected ≥{min} floats on line {line:?}, got {}", v.len());
    Ok(v)
}

/// A file ending inside its last frame is an MD run cut off: that frame is dropped with a
/// warning and the frames before it stand. When it was the only frame nothing is left, and
/// reading on would hand back an empty trajectory with no error, so that is an error.
/// Shared by the xyz, extxyz, PDB and LAMMPS dump readers.
pub(crate) fn drop_truncated(warnings: &mut Vec<crate::Notice>, k: usize, line: usize, why: &str) -> Result<()> {
    ensure!(k > 0, "the only frame (line {line}) is incomplete: {why}");
    warnings.push(crate::Notice::warning(Some(line), format!("last frame (frame {k}) is incomplete and was dropped: {why}")));
    Ok(())
}

/// A species column entry: an element symbol, kept as written, or an atomic number, which
/// becomes its symbol (as ASE reads both). `8` used to be the element name. Shared by the xyz
/// and extxyz readers; extxyz's `Z` column, a number by definition, does not come here.
pub(crate) fn species_element(token: &str) -> Result<String> {
    match token.parse::<u8>() {
        Ok(z) => Ok(ferro_core::data::elements::by_number(z)
            .ok_or_else(|| anyhow!("no element has atomic number {z}"))?
            .symbol.to_string()),
        Err(_) => Ok(token.to_string()),
    }
}

/// The element a site label names: its leading letters, capitalised as a symbol
/// (`fe2` -> `Fe`, `O1` -> `O`, `CA` -> `Ca`); `X` when it starts with no letter.
///
/// Not checked against the element table — the callers keep a label such as `Ow`
/// as written and leave validation to the analyses. Shared by the CP2K input, QE
/// input / output and CIF readers.
pub(crate) fn element_from_label(label: &str) -> String {
    let letters: String = label.chars().take_while(|c| c.is_alphabetic()).collect();
    if letters.is_empty() { "X".to_string() } else { capitalize_symbol(&letters) }
}

/// An element symbol in the case the element table uses: `FE` / `fe` → `Fe`. Shared by
/// [`element_from_label`] and the PDB and LAMMPS data readers (scan B, one case rule).
pub(crate) fn capitalize_symbol(s: &str) -> String {
    let mut ch = s.chars();
    match ch.next() {
        None => String::new(),
        Some(first) => first.to_uppercase().chain(ch.flat_map(char::to_lowercase)).collect(),
    }
}

/// Reads a whole text file.
///
/// Bytes that are not UTF-8 are an error naming the file, line and column of the
/// first one (review C-8): a parse past them would be guessing, and "cannot open"
/// for a file that opened fine sends the user looking for the wrong problem.
pub(crate) fn read_text(path: &Path) -> Result<String> {
    let bytes = std::fs::read(path).with_context(|| format!("cannot open {}", path.display()))?;
    String::from_utf8(bytes).map_err(|e| {
        let at = e.utf8_error().valid_up_to();
        let b = e.as_bytes();
        let line_start = b[..at].iter().rposition(|&c| c == b'\n').map_or(0, |p| p + 1);
        let line = 1 + b[..at].iter().filter(|&&c| c == b'\n').count();
        not_utf8(path, line, at - line_start + 1, b[at])
    })
}

/// The lines of a text file, read one at a time, ending as `str::lines` ends them
/// (`\n` or `\r\n` removed). A line that is not UTF-8 is the same error as in
/// [`read_text`], so the streaming readers report the same position.
pub(crate) fn text_lines(path: &Path) -> Result<TextLines> {
    let file = std::fs::File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    Ok(TextLines { reader: BufReader::with_capacity(1 << 16, file), path: path.to_path_buf(), n: 0, width: 0 })
}

pub(crate) struct TextLines {
    reader: BufReader<std::fs::File>,
    path: PathBuf,
    n: usize,
    /// Length of the previous line: sizing the next buffer from it saves the
    /// reallocations that dominate reading many short lines
    width: usize,
}

impl Iterator for TextLines {
    type Item = Result<String>;

    fn next(&mut self) -> Option<Result<String>> {
        let mut buf = Vec::with_capacity(self.width + 8);
        match self.reader.read_until(b'\n', &mut buf) {
            Ok(0) => return None,
            Ok(_) => {}
            Err(e) => return Some(Err(anyhow::Error::new(e).context(format!("reading {}", self.path.display())))),
        }
        self.n += 1;
        self.width = buf.len();
        if buf.last() == Some(&b'\n') {
            buf.pop();
            if buf.last() == Some(&b'\r') {
                buf.pop();
            }
        }
        Some(String::from_utf8(buf).map_err(|e| {
            let at = e.utf8_error().valid_up_to();
            not_utf8(&self.path, self.n, at + 1, e.as_bytes()[at])
        }))
    }
}

fn not_utf8(path: &Path, line: usize, column: usize, byte: u8) -> anyhow::Error {
    anyhow!(
        "{}: line {line}, column {column}: byte 0x{byte:02x} is not UTF-8 text \
         (a binary file, or another program's output mixed into this one?)",
        path.display()
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One bad byte is an error at its own line and column, from both the whole-file
    /// and the line-by-line path; good text reads unchanged, `\r\n` included.
    #[test]
    fn non_utf8_names_the_line_and_column() {
        let path = std::env::temp_dir().join("ferro_util_not_utf8.txt");
        std::fs::write(&path, b"first\r\nsec\xffond\nthird\n").unwrap();
        for err in [
            read_text(&path).unwrap_err(),
            text_lines(&path).unwrap().collect::<Result<Vec<_>>>().unwrap_err(),
        ] {
            let msg = err.to_string();
            assert!(msg.contains("line 2, column 4: byte 0xff"), "{msg}");
        }

        std::fs::write(&path, "a\r\nb\nc").unwrap();
        assert_eq!(read_text(&path).unwrap().lines().collect::<Vec<_>>(), ["a", "b", "c"]);
        assert_eq!(text_lines(&path).unwrap().collect::<Result<Vec<_>>>().unwrap(), ["a", "b", "c"]);
        let _ = std::fs::remove_file(&path);

        let missing = read_text(Path::new("/no/such/ferro/file")).unwrap_err().to_string();
        assert!(missing.starts_with("cannot open"), "{missing}");
    }

    #[test]
    fn element_from_label_takes_the_leading_letters() {
        for (label, want) in [("Fe1a", "Fe"), ("H12", "H"), ("CA", "Ca"), ("fe2", "Fe"),
                              ("O", "O"), ("Ow", "Ow"), ("1O", "X"), ("", "X")] {
            assert_eq!(element_from_label(label), want, "{label}");
        }
    }
}
