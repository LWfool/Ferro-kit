// The doc comments on the clap structs double as `--help` text, which is what users
// read: `[Å]` and `<suffix>` are right there, and escaping them for rustdoc would put
// the backslashes into `--help`.  Real links elsewhere in this crate lose the check too;
// that is the price of one text serving both readers.
#![allow(rustdoc::broken_intra_doc_links, rustdoc::invalid_html_tags)]

pub mod args;
pub mod batch;
pub mod cmd;
pub mod doc;
pub mod help;
pub mod io_dispatch;
pub mod outpath;
