//! Single-frame atom type classification for structural export and downstream use.
//!
//! Wraps `ferro_core::classify_frame` and adds helpers for applying labels to
//! a `Frame` (so it can be written out with type labels in place of element symbols).

use ferro_core::{classify_frame, AtomType, Frame, Trajectory, TypeParams};

/// Classify every atom in every frame of a trajectory.
///
/// Returns `Vec<Vec<AtomType>>` — outer index = frame, inner index = atom.
/// Frames without a cell are skipped, so the outer index is not a frame index
/// when the trajectory mixes periodic and non-periodic frames.
pub fn classify_trajectory(traj: &Trajectory, params: &TypeParams) -> Vec<Vec<AtomType>> {
    traj.frames.iter().filter_map(|f| {
        f.cell.as_ref().map(|cell| classify_frame(f, cell, params))
    }).collect()
}

/// Return a copy of `frame` with atom `element` fields replaced by type labels.
///
/// Atoms not covered by `params` (e.g. non-network elements) keep their original element.
pub fn apply_type_labels(frame: &Frame, labels: &[String]) -> Frame {
    let mut out = frame.clone();
    for (atom, label) in out.atoms.iter_mut().zip(labels.iter()) {
        atom.element = label.clone();
    }
    out
}
