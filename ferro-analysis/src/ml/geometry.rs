//! Per-frame geometric quantities used by the dataset filter.
//!
//! Both quantities are periodic and both are computed under the minimum-image
//! convention, which is exact only while the distance of interest stays below
//! half the smallest interplanar spacing. The callers check that bound and say
//! so rather than returning a quietly wrong number.

use ferro_core::{AtomType, Frame, PeriodicCell, TypeParams};

/// Smallest distance between two atoms of the given elements, minimum image.
///
/// `None` when the frame holds fewer than the two atoms the pair needs. When
/// `a == b` a distance is measured between distinct atoms only — the zero of an
/// atom with itself is not a contact, while its periodic images are, and those
/// are far away by construction (they sit at the box widths).
pub fn min_pair_distance(frame: &Frame, cell: &PeriodicCell, a: &str, b: &str) -> Option<f64> {
    let ia: Vec<usize> = idx_of(frame, a);
    let ib: Vec<usize> = if a == b { ia.clone() } else { idx_of(frame, b) };
    if ia.is_empty() || ib.is_empty() || (a == b && ia.len() < 2) {
        return None;
    }
    let mut best = f64::INFINITY;
    for (n, &i) in ia.iter().enumerate() {
        let pi = frame.positions()[i];
        // 同元素时只走上三角，距离对称
        let rest = if a == b { &ib[n + 1..] } else { &ib[..] };
        for &j in rest {
            let d = cell.minimum_image(frame.positions()[j] - pi);
            let r2 = d.norm_squared();
            if r2 < best {
                best = r2;
            }
        }
    }
    (best.is_finite()).then(|| best.sqrt())
}

fn idx_of(frame: &Frame, elem: &str) -> Vec<usize> {
    frame
        .atoms()
        .enumerate()
        .filter_map(|(i, a)| (a.element == elem).then_some(i))
        .collect()
}

/// How many atoms of `elem` have exactly `cn` neighbours, via the network classifier.
///
/// Goes through [`ferro_core::classify_frame`] rather than counting distances
/// here, so the answer is the same number `ferro net` reports — one definition
/// of coordination for the whole project. The coordination is read from the
/// `cn` FIELD of [`AtomType::Former`]; never parse it back out of `label()`.
pub fn count_with_coordination(
    frame: &Frame,
    params: &TypeParams,
    elem: &str,
    cn: u32,
) -> ferro_core::Result<usize> {
    Ok(ferro_core::classify_frame(frame, params)?
        .iter()
        .filter(|t| matches!(t, AtomType::Former { elem: e, cn: c, .. } if e == elem && *c == cn))
        .count())
}

/// Distribution of coordination numbers for one element, `cn -> count`.
///
/// Feeds the read-only diagnostics: a selection on "coordination == 6" is only
/// as meaningful as the spread of that distribution, and the spread is what
/// tells you whether the cutoff sits on a shoulder or in a gap.
pub fn coordination_histogram(
    frame: &Frame,
    params: &TypeParams,
    elem: &str,
) -> ferro_core::Result<Vec<(u32, usize)>> {
    let mut hist: std::collections::BTreeMap<u32, usize> = Default::default();
    for t in ferro_core::classify_frame(frame, params)? {
        if let AtomType::Former { elem: e, cn, .. } = t {
            if e == elem {
                *hist.entry(cn).or_default() += 1;
            }
        }
    }
    Ok(hist.into_iter().collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell};
    use nalgebra::{Matrix3, Vector3};

    fn cubic(a: f64, atoms: Vec<(&str, [f64; 3])>) -> (Frame, PeriodicCell) {
        let cell = Cell::from_matrix(Matrix3::identity() * a);
        let mut f = Frame::with_cell(cell, [true; 3]);
        f.set_atoms(atoms
            .into_iter()
            .map(|(e, p)| Atom::new(e, Vector3::new(p[0], p[1], p[2])))
            );
        let pc = PeriodicCell::of(&f).unwrap();
        (f, pc)
    }

    #[test]
    fn same_element_pair_skips_the_atom_itself() {
        let (f, c) = cubic(10.0, vec![("O", [0.0; 3]), ("O", [2.5, 0.0, 0.0])]);
        assert!((min_pair_distance(&f, &c, "O", "O").unwrap() - 2.5).abs() < 1e-12);
    }

    #[test]
    fn the_periodic_image_wins_when_it_is_closer() {
        // 盒子 10 Å，两个 O 相隔 9 Å —— 镜像只有 1 Å
        let (f, c) = cubic(10.0, vec![("O", [0.5, 0.0, 0.0]), ("O", [9.5, 0.0, 0.0])]);
        assert!((min_pair_distance(&f, &c, "O", "O").unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn a_missing_element_gives_none() {
        let (f, c) = cubic(10.0, vec![("O", [0.0; 3]), ("O", [2.5, 0.0, 0.0])]);
        assert!(min_pair_distance(&f, &c, "Al", "O").is_none());
        // 只有一个 O 时同元素对不存在
        let (f1, c1) = cubic(10.0, vec![("O", [0.0; 3])]);
        assert!(min_pair_distance(&f1, &c1, "O", "O").is_none());
    }

    #[test]
    fn coordination_counts_match_the_network_classifier() {
        // 一个 Al 被 6 个 O 围住，恰好 6 配位
        let mut atoms = vec![("Al", [5.0, 5.0, 5.0])];
        for d in [
            [1.9, 0.0, 0.0], [-1.9, 0.0, 0.0], [0.0, 1.9, 0.0],
            [0.0, -1.9, 0.0], [0.0, 0.0, 1.9], [0.0, 0.0, -1.9],
        ] {
            atoms.push(("O", [5.0 + d[0], 5.0 + d[1], 5.0 + d[2]]));
        }
        let (f, _) = cubic(20.0, atoms);
        let mut cut = std::collections::BTreeMap::new();
        cut.insert(("Al".to_string(), "O".to_string()), 2.4);
        let p = TypeParams::new(cut, Default::default());
        assert_eq!(count_with_coordination(&f, &p, "Al", 6).unwrap(), 1);
        assert_eq!(count_with_coordination(&f, &p, "Al", 4).unwrap(), 0);
        assert_eq!(coordination_histogram(&f, &p, "Al").unwrap(), vec![(6, 1)]);
    }
}
