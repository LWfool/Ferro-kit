//! Single-frame atom type classification for structural export and downstream use.
//!
//! Wraps `ferro_core::classify_frame` and adds helpers for applying labels to
//! a `Frame` (so it can be written out with type labels in place of element symbols).

use ferro_core::error::Result;
use ferro_core::{classify_frame, AtomType, Frame, Trajectory, TypeParams};

/// Classify every atom in every frame of a trajectory.
///
/// Returns `Vec<Vec<AtomType>>` — outer index = frame, inner index = atom.
/// Frames without a cell are skipped, so the outer index is not a frame index
/// when the trajectory mixes periodic and non-periodic frames. A singular cell
/// is an error.
pub fn classify_trajectory(traj: &Trajectory, params: &TypeParams) -> Result<Vec<Vec<AtomType>>> {
    traj.frames.iter().filter_map(|f| {
        f.cell.as_ref().map(|_| classify_frame(f, params))
    }).collect()
}

/// Return a copy of `frame` with the type labels written to `Atom::label`.
///
/// `Atom::element` is left alone: a `Frame` whose element field is not an element
/// cannot be used for anything else (its mass, its scattering factor, `-a`/`-b`
/// selection all break), and the labelled frame stays usable in the same process —
/// which the script/REPL mode will depend on.
///
/// Formats that have nowhere to put a second per-atom string fold the label into
/// their element column at **write** time; see `ferro_io::dump_name`.
pub fn apply_type_labels(frame: &Frame, labels: &[String]) -> Frame {
    let mut out = frame.clone();
    for (atom, label) in out.atoms.iter_mut().zip(labels.iter()) {
        atom.label = Some(label.clone());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::Atom;
    use nalgebra::Vector3;

    fn frame_of(specs: &[(&str, Option<&str>)]) -> Frame {
        let mut f = Frame::new();
        for (elem, label) in specs {
            let mut a = Atom::new(*elem, Vector3::zeros());
            a.label = label.map(str::to_string);
            f.add_atom(a);
        }
        f
    }

    #[test]
    fn test_apply_writes_label_and_leaves_element_alone() {
        let f = frame_of(&[("P", None), ("O", None)]);
        let out = apply_type_labels(&f, &["P_3".into(), "O_b".into()]);
        assert_eq!(out.atoms[0].element, "P", "element 必须仍是元素");
        assert_eq!(out.atoms[0].label.as_deref(), Some("P_3"));
        assert_eq!(out.atoms[1].element, "O");
        assert_eq!(out.atoms[1].label.as_deref(), Some("O_b"));
    }
}
