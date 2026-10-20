//! Every message ferro prints outside its products: one line on stderr,
//! `<source>[:<line>]: <category>: <message>`, Python's warning layout.
//!
//! The source is the input path as the user typed it, or [`PROGRAM`] for a message
//! about no file in particular. The categories are `FerroNote` (worth knowing),
//! `FerroWarning` (may be wrong) and `FerroError` (failed). stdout carries results
//! and progress only, so a pipe or a redirect never mixes the two.
//!
//! ferro-py renders the same layout in Python (`ferro/cli.py`, `_show_warning`): two
//! implementations of one format, too small to share across the language boundary.

use std::fmt::Display;
use std::path::Path;

use ferro_io::{Level, Notice};

/// The source of a message about no file in particular.
pub const PROGRAM: &str = "ferro";

fn emit(source: &dyn Display, line: Option<usize>, category: &str, message: &dyn Display) {
    match line {
        Some(n) => eprintln!("{source}:{n}: {category}: {message}"),
        None => eprintln!("{source}: {category}: {message}"),
    }
}

pub fn note(source: impl Display, message: impl Display) {
    emit(&source, None, &Level::Note.to_string(), &message);
}

pub fn warn(source: impl Display, message: impl Display) {
    emit(&source, None, &Level::Warning.to_string(), &message);
}

pub fn error(source: impl Display, message: impl Display) {
    emit(&source, None, "FerroError", &message);
}

/// What a reader noticed about `path`.
pub fn notices(path: &Path, notices: &[Notice]) {
    for n in notices {
        emit(&path.display(), n.line, &n.level.to_string(), &n.message);
    }
}
