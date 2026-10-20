//! What a reader noticed beyond the frames themselves.
//!
//! The readers never print: they hand these back, and the caller renders them as
//! `<source>[:<line>]: <level>: <message>` (the CLI's `report.rs`, ferro-py's
//! `showwarning`). The source is the caller's to supply — it knows the path as the
//! user typed it — so a notice carries only the line, when there is one.

use std::fmt;

/// How much a notice matters. Errors are not notices: a reader that cannot go on returns `Err`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Level {
    /// Information worth knowing; nothing looks wrong
    Note,
    /// Something may be wrong with the input or with how it was read
    Warning,
}

impl fmt::Display for Level {
    /// The category name every channel shows: `FerroNote`, `FerroWarning`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Level::Note => "FerroNote",
            Level::Warning => "FerroWarning",
        })
    }
}

/// One message from a reader.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Notice {
    pub level: Level,
    /// 1-based line of the file the message is about, when there is one
    pub line: Option<usize>,
    pub message: String,
}

impl Notice {
    pub fn note(line: Option<usize>, message: impl Into<String>) -> Self {
        Notice { level: Level::Note, line, message: message.into() }
    }

    pub fn warning(line: Option<usize>, message: impl Into<String>) -> Self {
        Notice { level: Level::Warning, line, message: message.into() }
    }
}
