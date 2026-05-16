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
pub mod charge_grid;
pub mod frame;
pub mod trajectory;
pub mod cube_data;
pub mod data;
pub mod units;
pub mod error;
pub mod network_type;
pub mod cluster;
pub mod spin;

// top-level re-exports for downstream crates
pub use atom::Atom;
pub use cell::Cell;
pub use charge_grid::ChargeGrid;
pub use frame::Frame;
pub use trajectory::{Trajectory, TrajectoryMetadata};
pub use cube_data::CubeData;
pub use error::{ChemError, Result};
pub use spin::{
    assign_oxidation_states, guess_spin, parity_min_multiplicity,
    total_electron_count, SpinGuess, SpinMethod,
};
pub use network_type::{
    TypeParams, CutoffTable, classify_frame,
    oxygen_label_order, modifier_label_order, former_label_order,
};
pub use cluster::{
    build_network_graph, connected_components, LigandKind, NetworkGraph,
};
