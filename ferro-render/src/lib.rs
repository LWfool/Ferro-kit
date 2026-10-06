//! Terminal rendering of the manual's markdown, for `ferro doc`.
//!
//! Used only when stdout is a terminal: rendered output is escape codes and
//! box-drawing characters, which in a redirected file are noise, so
//! `ferro doc net > net.md` still gets the source as written.
//!
//! A crate of its own so that the rest of the CLI recompiles without it: the
//! rendering style rarely changes. [`markdown`] lays out the page, [`latex`]
//! turns formulas into Unicode; nothing here knows about ferro's data.
//!
//! What a renderer cannot convert — a table without an alignment row, a math
//! environment, a TeX command not known here — comes out as written, so new
//! manual content can fail to be rendered but never be rendered wrong.

mod latex;
mod markdown;

pub use markdown::{render, terminal_width, Style};
