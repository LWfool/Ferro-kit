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

/// Return a copy of `frame` with the type labels written to `Atom::label`.
///
/// `Atom::element` is left alone: a `Frame` whose element field is not an element
/// cannot be used for anything else (its mass, its scattering factor, `-a`/`-b`
/// selection all break), and the labelled frame stays usable in the same process —
/// which the script/REPL mode will depend on.
///
/// Formats that have nowhere to put a second per-atom string fold the label into
/// their element column at **write** time; see `ferro_structure::fold_labels`.
pub fn apply_type_labels(frame: &Frame, labels: &[String]) -> Frame {
    let mut out = frame.clone();
    for (atom, label) in out.atoms.iter_mut().zip(labels.iter()) {
        atom.label = Some(label.clone());
    }
    out
}

/// Return a copy of `frame` with `Atom::label` folded into `Atom::element`.
///
/// For formats with a single per-atom name column (LAMMPS dump), where an extra
/// column would break every downstream tool that already parses the file. The dump
/// reader splits `<element>_<suffix>` back apart, so the round trip is exact.
///
/// A label that does not start with `<element>_` is **not** folded: CIF/CP2K/QE site
/// names (`O1`, `Fe1`) would otherwise land in the element column and read back as
/// the bogus element `O1`. Returns the number of atoms left unfolded for that reason.
pub fn fold_labels(frame: &Frame) -> (Frame, usize) {
    let mut out = frame.clone();
    let mut skipped = 0usize;
    for atom in out.atoms.iter_mut() {
        let Some(label) = atom.label.clone() else { continue };
        if label == atom.element {
            continue;
        }
        if label.starts_with(&format!("{}_", atom.element)) {
            atom.element = label;
        } else {
            skipped += 1;
        }
    }
    (out, skipped)
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

    #[test]
    fn test_fold_only_touches_conforming_labels() {
        // CIF / CP2K / QE 的位点名(O1、Fe1)不合 <元素>_<后缀> 约定,
        // 折进 element 列会读回成 "O1" 这个不存在的元素
        let f = frame_of(&[
            ("P",  Some("P_3")),   // 合规 → 折叠
            ("O",  Some("O1")),    // CIF 位点名 → 不折
            ("Fe", Some("Fe1")),   // CP2K kind 名 → 不折
            ("Zn", Some("Zn")),    // 与元素相同 → 无需折叠,也不算跳过
            ("Ar", None),          // 无标签
        ]);
        let (out, skipped) = fold_labels(&f);
        assert_eq!(out.atoms[0].element, "P_3");
        assert_eq!(out.atoms[1].element, "O");
        assert_eq!(out.atoms[2].element, "Fe");
        assert_eq!(out.atoms[3].element, "Zn");
        assert_eq!(out.atoms[4].element, "Ar");
        assert_eq!(skipped, 2, "两个不合约定的标签被跳过并计数");
    }
}
