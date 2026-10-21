//! Machine learning, one namespace per task.
//!
//! Sits beside `md`, `network` and `dft`: pure computation over a
//! [`ferro_core::Trajectory`], no filesystem. Each task is a submodule with its own
//! API; `dataset` (what `ferro dataset` runs) is the only one so far.

pub mod dataset;
