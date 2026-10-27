pub mod box_builder;
mod compounds;
pub mod merge;
pub mod supercell;
pub mod typing;
pub mod vacuum;

pub use box_builder::{build_box, estimate_box_length, Component};
pub use merge::merge_frames;
pub use supercell::{check_supercell_total, find_supercell_dims, make_supercell, MAX_SUPERCELL_ATOMS};
pub use typing::apply_type_labels;
pub use vacuum::add_vacuum;
