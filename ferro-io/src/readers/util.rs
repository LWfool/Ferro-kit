//! Helpers shared by more than two readers.

use anyhow::{ensure, Result};

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
