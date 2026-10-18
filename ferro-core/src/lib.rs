//! ferro-core — core data structures and utilities.
//!
//! Type hierarchy:
//! ```text
//! Trajectory (Vec<Frame>)
//!   └── Frame  (atoms + cell + pbc + computed results)
//!         └── Atom (element + position + optional properties)
//!               └── data::elements (static element table, used for mass lookup)
//! ```

pub mod array_order;
pub mod atom;
pub mod cell;
pub mod periodic;
pub mod charge_grid;
pub mod frame;
pub mod trajectory;
pub mod cube_data;
pub mod table;
pub mod data;
pub mod units;
pub mod error;

// top-level re-exports for downstream crates
pub use array_order::{matrix3_from_row_major, matrix3_row_major};
pub use atom::{Atom, AtomRef, Site};
pub use cell::Cell;
pub use periodic::PeriodicCell;
pub use charge_grid::ChargeGrid;
pub use frame::Frame;
pub use trajectory::{select_range, spread_range, Trajectory, TrajectoryMetadata};
pub use cube_data::CubeData;
pub use table::{Column, Table};
pub use error::{ChemError, Result};
