//! Subcommand implementations.
//!
//! Grouping follows **what a command produces**, not which crate implements it:
//! `traj` emits stacked tables, `map` emits 3-D grid files, `net` emits topology
//! statistics, and `convert` / `info` / `job` / `bader` are one-offs at the top level.

pub mod bader;
pub mod convert;
pub mod info;
pub mod job;
pub mod map;
pub mod net;
pub mod traj;
