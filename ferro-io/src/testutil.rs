//! Shared helpers for the reader/writer round-trip tests.
//!
//! Every reader test needs the same thing: put a literal into a file, hand the path to
//! the parser.  That was written out thirteen times, in three shapes, before this module
//! existed — the two shapes that appear more than once live here.

use std::path::PathBuf;

/// Writes `content` under the OS temp dir and returns the path as a `String`.
///
/// The `String` shape is what most reader signatures still take (see the `&str` → `&Path`
/// item in `dev/plan.md`); [`write_tmp`] returns the `PathBuf` the rest expect.
pub(crate) fn write_tmp_str(name: &str, content: &str) -> String {
    write_tmp(name, content).to_str().unwrap().to_string()
}

/// Writes `content` under the OS temp dir and returns the path.
pub(crate) fn write_tmp(name: &str, content: &str) -> PathBuf {
    let path = std::env::temp_dir().join(name);
    std::fs::write(&path, content).unwrap();
    path
}
