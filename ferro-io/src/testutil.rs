//! Shared helpers for the reader/writer round-trip tests.
//!
//! Every reader test needs the same thing: put a literal into a file, hand the path to
//! the parser.  That was written out thirteen times, in three shapes, before this module
//! existed.

use std::path::PathBuf;

/// Writes `content` under the OS temp dir and returns the path.
pub(crate) fn write_tmp(name: &str, content: &str) -> PathBuf {
    let path = std::env::temp_dir().join(name);
    std::fs::write(&path, content).unwrap();
    path
}
