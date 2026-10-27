//! Applying per-atom type labels to a `Frame`, so it can be written out with them.
//!
//! The labels come from the network classification in `ferro-analysis`; this crate
//! only puts them on the atoms.

use ferro_core::error::{ChemError, Result};
use ferro_core::Frame;

/// Return a copy of `frame` with the type labels written to `Site::label`.
///
/// `Site::element` is left alone: a `Frame` whose element field is not an element
/// cannot be used for anything else (its mass, its scattering factor, `-a`/`-b`
/// selection all break), and the labelled frame stays usable in the same process.
///
/// Formats that have nowhere to put a second per-atom string fold the label into
/// their element column at **write** time; see `ferro_io::dump_name`.
///
/// `labels` is either empty (the frame is returned unchanged) or one per atom;
/// any other length is a `ValidationError` rather than a silent partial labelling.
pub fn apply_type_labels(frame: &Frame, labels: &[String]) -> Result<Frame> {
    let mut out = frame.clone();
    if labels.is_empty() {
        return Ok(out);
    }
    if labels.len() != frame.n_atoms() {
        return Err(ChemError::ValidationError(format!(
            "{} type labels for {} atoms", labels.len(), frame.n_atoms()
        )));
    }
    for (atom, label) in out.sites_mut().iter_mut().zip(labels.iter()) {
        atom.label = Some(label.clone());
    }
    Ok(out)
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
        let out = apply_type_labels(&f, &["P_3".into(), "O_b".into()]).unwrap();
        assert_eq!(out.element(0), "P", "element 必须仍是元素");
        assert_eq!(out.atom(0).label, Some("P_3"));
        assert_eq!(out.element(1), "O");
        assert_eq!(out.atom(1).label, Some("O_b"));
    }

    /// Review R6-5: `zip` used to stop at the shorter side, labelling only some atoms.
    #[test]
    fn labels_are_none_or_one_per_atom() {
        let f = frame_of(&[("P", Some("P_x")), ("O", None)]);
        let out = apply_type_labels(&f, &[]).unwrap();
        assert_eq!(out.atom(0).label, Some("P_x"));
        assert_eq!(out.atom(1).label, None);
        for labels in [vec!["P_3".to_string()], vec!["P_3".into(), "O_b".into(), "O_t".into()]] {
            let e = apply_type_labels(&f, &labels).unwrap_err().to_string();
            assert!(e.contains(&format!("{} type labels for 2 atoms", labels.len())), "{e}");
        }
    }
}
