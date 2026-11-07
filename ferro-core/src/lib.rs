//! ferro-core — core data structures and utilities.
//!
//! Type hierarchy:
//! ```text
//! Trajectory (Vec<Frame>)
//!   └── Frame  (cell + pbc + per-frame results: energy, stress, temperature, step)
//!         ├── Arc<Vec<Site>>  what each atom is (element, label, mass, magmom),
//!         │                   one table shared by consecutive frames
//!         └── per-atom arrays positions, charges, forces, velocities, the frame's own
//! Atom     an owned atom, what a frame is built from (Frame::add_atom / set_atoms)
//! AtomRef  the view of one atom a frame hands out (Frame::atom / atoms)
//! data::elements  static element table (symbol, Z, mass)
//! ```

mod array_order;
mod atom;
mod cell;
mod periodic;
mod charge_grid;
mod frame;
mod trajectory;
mod cube_data;
mod table;
pub mod data;
pub mod units;
mod error;

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
