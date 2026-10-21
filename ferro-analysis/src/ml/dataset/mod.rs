//! ML dataset preparation: the computation behind `ferro dataset filter / merge`
//! and the inspect tables.
//!
//! Reading the dataset and writing the result belong to `ferro-io`; this module
//! only decides. Every public item has one path, `ml::dataset::<item>`; the
//! submodules are visible only inside the crate.

pub(crate) mod diagnostics;
pub(crate) mod filter;
mod geometry;
pub(crate) mod merge;
mod split;

pub use diagnostics::{filter_diagnostics, frame_table, quantiles};
pub use filter::{filter_frames, Criterion, FilterParams, FilterResult, FrameVerdict};
pub use merge::{composition_key, group_name, shuffle_order, shuffle_within, stitch_segments, DEFAULT_SEED};
pub use split::Split;
