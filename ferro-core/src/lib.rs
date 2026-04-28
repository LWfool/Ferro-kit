//! ferro-core — core data structures and utilities.
//!
//! Type hierarchy:
//! ```text
//! Trajectory (Vec<Frame>)
//!   └── Frame  (atoms + cell + pbc + computed results)
//!         └── Atom (element + position + optional properties)
//!               └── data::elements (static element table, used for mass lookup)
//! ```

pub mod atom;
pub mod cell;
pub mod frame;
pub mod trajectory;
pub mod cube_data;
pub mod data;
pub mod units;
pub mod error;

// top-level re-exports for downstream crates
pub use atom::Atom;
pub use cell::Cell;
pub use frame::Frame;
pub use trajectory::{Trajectory, TrajectoryMetadata};
pub use cube_data::CubeData;
pub use error::{ChemError, Result};
