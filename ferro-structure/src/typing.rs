//! Applying per-atom type labels to a `Frame`, so it can be written out with them.
//!
//! The labels come from the network classification in `ferro-analysis`; this crate
//! only puts them on the atoms.

use ferro_core::Frame;

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
    for (atom, label) in out.sites_mut().iter_mut().zip(labels.iter()) {
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
        assert_eq!(out.element(0), "P", "element 必须仍是元素");
        assert_eq!(out.atom(0).label, Some("P_3"));
        assert_eq!(out.element(1), "O");
        assert_eq!(out.atom(1).label, Some("O_b"));
    }
}
