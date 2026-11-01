//! Shared helpers for the reader/writer round-trip tests.
//!
//! Every reader test needs the same thing: put a literal into a file, hand the path to
//! the parser.  That was written out thirteen times, in three shapes, before this module
//! existed.

use std::path::PathBuf;

/// `name` in this test process's own directory under the OS temp dir. Tests used to share
/// fixed names in the temp dir itself, so two `cargo test` runs at once could overwrite each
/// other's files (review 2026-10-30 A1-L14); a per-process directory keeps them apart.
pub(crate) fn tmp_path(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("ferro-io-tests-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir.join(name)
}

/// Writes `content` to [`tmp_path`]`(name)` and returns the path.
pub(crate) fn write_tmp(name: &str, content: &str) -> PathBuf {
    let path = tmp_path(name);
    std::fs::write(&path, content).unwrap();
    path
}

/// Two-atom BCC Fe in its cubic cell (a = 2.87 A): the frame the writer tests round-trip.
pub(crate) fn bcc_fe() -> ferro_core::Frame {
    use ferro_core::{Atom, Cell, Frame};
    use nalgebra::Vector3;
    let cell = Cell::from_lengths_angles(2.87, 2.87, 2.87, 90.0, 90.0, 90.0).unwrap();
    let mut frame = Frame::with_cell(cell, [true; 3]);
    frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
    frame.add_atom(Atom::new("Fe", Vector3::new(1.435, 1.435, 1.435)));
    frame
}
