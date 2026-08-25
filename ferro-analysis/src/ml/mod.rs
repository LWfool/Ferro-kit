//! Machine-learning dataset preparation.
//!
//! Sits beside `md`, `network` and `dft`: pure computation over a
//! [`ferro_core::Trajectory`], no filesystem. Reading the dataset and writing
//! the result belong to `ferro-io`; this module only decides.

pub mod filter;

pub use filter::{filter_frames, Criterion, FilterParams, FilterResult, FrameVerdict};
