//! Machine-learning dataset preparation.
//!
//! Sits beside `md`, `network` and `dft`: pure computation over a
//! [`ferro_core::Trajectory`], no filesystem. Reading the dataset and writing
//! the result belong to `ferro-io`; this module only decides.

pub mod filter;
pub mod diagnostics;
pub mod geometry;
pub mod merge;
pub mod split;

pub use filter::{filter_frames, Criterion, FilterParams, FilterResult, FrameVerdict};
pub use merge::{composition_key, group_name, shuffle_order, stitch_segments};
pub use split::Split;
pub use geometry::{coordination_histogram, count_with_coordination, min_pair_distance};
