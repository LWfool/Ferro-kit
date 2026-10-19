//! Combining datasets of the same composition.
//!
//! Grouping does NOT go by directory name — `init.011` and
//! `Al2ZnO4_1.000_0300K.train` say nothing reliable about what is inside. It
//! goes by the per-atom element sequence, so two systems merge exactly when
//! they hold the same atoms.
//!
//! # Atom order
//!
//! Systems of one composition may still list their atoms in different orders
//! and carry different `type_map` orders. Merging sorts every system's atoms
//! into one canonical order ([`Trajectory::canonically_sorted`], in `ferro-core`
//! because `convert` and the Python jobs sort too) and permutes the per-atom
//! arrays along with them.
//! A DP model is invariant under atom renumbering, so this changes notation
//! rather than physics — the same thing dpdata does with `sort_atom_names` plus
//! `sort_atom_types`.
//!
//! The canonical order is **(Z, symbol)**, matching the `type_map` that
//! `collect` writes and the grouping order `ferro traj gr` uses. dpdata sorts
//! alphabetically instead; both are self-describing through `type_map.raw`, and
//! being consistent inside ferro matters more than matching dpdata's choice.

use std::collections::BTreeMap;
use std::path::Path;

use ferro_core::data::elements::canonical_order;
use ferro_core::error::ChemError;
use ferro_core::Trajectory;
use rand::seq::SliceRandom;
use rand::SeedableRng;

/// Default shuffle seed, so a run without `--seed` is still reproducible.
pub const DEFAULT_SEED: u64 = 666;

/// The key two systems must share to be merged: their sorted element sequence.
pub fn composition_key(traj: &Trajectory) -> Vec<String> {
    let mut v: Vec<String> = traj
        .frames
        .first()
        .map(|f| f.atoms().map(|a| a.element.to_string()).collect())
        .unwrap_or_default();
    let order = canonical_order(&v);
    v = order.into_iter().map(|i| v[i].clone()).collect();
    v
}

/// `<natoms>_<formula>`, e.g. `7_Al2O4Zn` — the output directory name.
///
/// The formula is NOT reduced: the subscripts are the actual atom counts, so
/// `Al96O192Zn48` is 336 atoms. The count prefix is there to read the scale at a
/// glance and to make `ls` group systems of equal size together. Elements run
/// alphabetically, which is how a formula is normally read, independent of the
/// (Z, symbol) order used for `type_map`.
///
/// A count of 1 drops its subscript (`Al2O4Zn`), unlike `formula_of`
/// (`Al2O4Zn1`).  The two are deliberately different, not an oversight: this one is read
/// as a name in an `ls` listing, that one is read stacked against a second formula in a
/// mismatch message, where the padded columns are what make the difference visible.
pub fn group_name(traj: &Trajectory) -> String {
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    if let Some(f) = traj.frames.first() {
        for a in f.atoms() {
            *counts.entry(a.element.to_string()).or_default() += 1;
        }
    }
    let n: usize = counts.values().sum();
    let formula: String = counts
        .iter()
        .map(|(e, c)| if *c == 1 { e.clone() } else { format!("{e}{c}") })
        .collect();
    format!("{n}_{formula}")
}

/// Puts the restart segments of one run back together, as `collect` does.
///
/// `parts` are the files of one run directory; `view` gives each one's path, frames
/// and first MD step. The parts are sorted in place by first step, then path, with
/// each file keeping its own frame order: a global sort by step looks more
/// thorough, but a restart that counts from 0 again would interleave two real
/// trajectories, which is worse than not sorting; the worst case here is
/// concatenating by file name. Every part must hold the composition of the first —
/// different atom sequences are two systems put in one directory, a person's
/// mistake rather than bad data, and dropping them as bad frames would report a
/// different thing altogether. The sequences compared are canonical (each file is
/// sorted as it is read), so any difference left is one of counts.
///
/// Returns the concatenated trajectory, with the first part's metadata, and the
/// file name each frame came from: the seams between segments are what deserve a
/// look. Frames are not de-duplicated; the step numbers make an overlap visible.
pub fn stitch_segments<T>(
    parts: &mut [T],
    view: impl Fn(&T) -> (&Path, &Trajectory, Option<i64>),
) -> ferro_core::Result<(Trajectory, Vec<String>)> {
    parts.sort_by(|a, b| {
        let ((pa, _, sa), (pb, _, sb)) = (view(a), view(b));
        (sa, pa).cmp(&(sb, pb))
    });
    let Some(first) = parts.first() else {
        return Ok((Trajectory::new(), Vec::new()));
    };
    let (first_path, first_traj, _) = view(first);
    let reference = symbols_of(first_traj);
    for part in &parts[1..] {
        let (path, traj, _) = view(part);
        let here = symbols_of(traj);
        if here != reference {
            return Err(ChemError::ValidationError(format!(
                "{} and {} hold different compositions ({} vs {}; {}); \
                 a system holds one composition, so put them in separate directories",
                first_path.display(),
                path.display(),
                formula_of(&reference),
                formula_of(&here),
                composition_diff(&reference, &here),
            )));
        }
    }

    let mut all = Trajectory::new();
    all.metadata = first_traj.metadata.clone();
    let mut sources: Vec<String> = Vec::new();
    for part in parts.iter() {
        let (path, traj, _) = view(part);
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.display().to_string());
        sources.extend(std::iter::repeat_n(name, traj.n_frames()));
        all.frames.extend(traj.frames.iter().cloned());
    }
    Ok((all, sources))
}

fn symbols_of(traj: &Trajectory) -> Vec<String> {
    traj.frames
        .first()
        .map(|f| f.symbols().into_iter().map(|s| s.to_string()).collect())
        .unwrap_or_default()
}

/// `Al32O64Zn16` from a per-atom element sequence, for the mismatch message.
///
/// Counts of 1 keep their subscript here, unlike [`group_name`] — the two
/// formulas in a mismatch message are read one above the other, and dropping the `1`s
/// misaligns the element that differs.  Directory names have the opposite need.
fn formula_of(symbols: &[String]) -> String {
    let mut count: BTreeMap<&str, usize> = BTreeMap::new();
    for s in symbols {
        *count.entry(s.as_str()).or_default() += 1;
    }
    count.iter().map(|(el, n)| format!("{el}{n}")).collect()
}

/// The elements whose counts differ, as `Al 2 vs 3`.
///
/// Two formulas printed side by side answer "are they different" but not "where":
/// `Al2O303P100Zn50 vs Al2O303P100Zn50` used to be a possible message, back when
/// the comparison ran on unsorted sequences and the formulas could not show it.
/// The sequences are canonical now, so every remaining difference is a count.
fn composition_diff(a: &[String], b: &[String]) -> String {
    let count = |v: &[String]| {
        let mut m: BTreeMap<String, usize> = BTreeMap::new();
        for s in v {
            *m.entry(s.clone()).or_default() += 1;
        }
        m
    };
    let (ca, cb) = (count(a), count(b));
    let diffs: Vec<String> = ca
        .keys()
        .chain(cb.keys())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .filter_map(|el| {
            let (na, nb) = (ca.get(el).copied().unwrap_or(0), cb.get(el).copied().unwrap_or(0));
            (na != nb).then(|| format!("{el} {na} vs {nb}"))
        })
        .collect();
    if diffs.is_empty() {
        "same counts, different atom order".to_string()
    } else {
        format!("differing: {}", diffs.join(", "))
    }
}

/// A deterministic shuffle of `0..n`.
pub fn shuffle_order(n: usize, seed: u64) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..n).collect();
    // 不用 StdRng：rand 声明它的算法可能随版本更换，同 seed 的划分就复现不了
    let mut rng = chacha20::ChaCha12Rng::seed_from_u64(seed);
    idx.shuffle(&mut rng);
    idx
}

/// A deterministic shuffle of each `lo..hi` span on its own: a frame moves only
/// within its span. One generator serves the spans in turn, so two spans of
/// equal length do not get the same permutation.
pub fn shuffle_within(spans: &[(usize, usize)], seed: u64) -> Vec<usize> {
    let mut rng = chacha20::ChaCha12Rng::seed_from_u64(seed);
    let mut order = Vec::new();
    for &(lo, hi) in spans {
        let mut idx: Vec<usize> = (lo..hi).collect();
        idx.shuffle(&mut rng);
        order.extend(idx);
    }
    order
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame};
    use nalgebra::{Matrix3, Vector3};

    fn traj_of(elements: &[&str]) -> Trajectory {
        let cell = Cell::from_matrix(Matrix3::identity() * 10.0);
        let mut f = Frame::with_cell(cell, [true; 3]);
        f.set_atoms(elements
            .iter()
            .enumerate()
            .map(|(i, e)| Atom::new(*e, Vector3::new(i as f64, 0.0, 0.0)))
            );
        f.set_forces(Some(
            (0..elements.len())
                .map(|i| Vector3::new(i as f64 * 0.1, 0.0, 0.0))
                .collect(),
        )).unwrap();
        Trajectory { frames: vec![f], metadata: Default::default() }
    }

    #[test]
    fn canonical_order_is_by_z_then_symbol() {
        // O(8) < Al(13) < Zn(30)
        let t = traj_of(&["Zn", "Al", "O", "Al"]);
        assert_eq!(composition_key(&t), vec!["O", "Al", "Al", "Zn"]);
    }

    /// Different atom orders of one composition must land in the same group.
    #[test]
    fn different_orders_of_one_composition_share_a_key() {
        assert_eq!(
            composition_key(&traj_of(&["Zn", "Al", "O", "O"])),
            composition_key(&traj_of(&["O", "Zn", "O", "Al"]))
        );
        assert_ne!(
            composition_key(&traj_of(&["Zn", "Al", "O", "O"])),
            composition_key(&traj_of(&["Zn", "Al", "O"]))
        );
    }

    #[test]
    fn a_composition_is_rendered_for_the_mismatch_message() {
        let syms: Vec<String> = ["O", "Al", "O", "Zn", "O"].iter().map(|s| s.to_string()).collect();
        assert_eq!(formula_of(&syms), "Al1O3Zn1");
    }

    /// Segments go in first-step order, each keeping its own frames, with the file
    /// each frame came from; a second composition is refused, naming both files.
    #[test]
    fn segments_are_stitched_by_first_step_and_must_share_a_composition() {
        let seg = |x: f64| {
            let mut t = traj_of(&["O", "Al"]);
            t.frames[0].positions_mut()[0].x = x;
            t
        };
        let mut parts = vec![
            (Path::new("run/b.out").to_path_buf(), seg(2.0), Some(200)),
            (Path::new("run/a.out").to_path_buf(), seg(1.0), Some(100)),
            (Path::new("run/c.out").to_path_buf(), seg(0.0), None),
        ];
        fn view(p: &(std::path::PathBuf, Trajectory, Option<i64>)) -> (&Path, &Trajectory, Option<i64>) {
            (p.0.as_path(), &p.1, p.2)
        }
        let (all, sources) = stitch_segments(&mut parts, view).unwrap();
        let xs: Vec<f64> = all.frames.iter().map(|f| f.positions()[0].x).collect();
        assert_eq!(xs, [0.0, 1.0, 2.0], "no step sorts first, then by first step");
        assert_eq!(sources, ["c.out", "a.out", "b.out"]);

        parts.push((Path::new("run/d.out").to_path_buf(), traj_of(&["O", "Al", "Al"]), Some(300)));
        let msg = stitch_segments(&mut parts, view).unwrap_err().to_string();
        assert!(msg.contains("run/d.out") && msg.contains("Al 1 vs 2"), "{msg}");
    }

    #[test]
    fn group_name_counts_atoms_and_spells_the_formula_alphabetically() {
        assert_eq!(group_name(&traj_of(&["Zn", "Al", "O", "O", "Al", "O", "O"])), "7_Al2O4Zn");
        assert_eq!(group_name(&traj_of(&["O", "H", "H"])), "3_H2O");
    }

    #[test]
    fn shuffle_is_a_permutation_and_reproducible() {
        let a = shuffle_order(50, 42);
        let b = shuffle_order(50, 42);
        assert_eq!(a, b);
        assert_ne!(a, shuffle_order(50, 43));
        let mut s = a.clone();
        s.sort();
        assert_eq!(s, (0..50).collect::<Vec<_>>());
    }

    #[test]
    fn shuffle_within_keeps_every_index_in_its_span() {
        let spans = [(0, 6), (6, 12), (12, 15)];
        let a = shuffle_within(&spans, 42);
        assert_eq!(a, shuffle_within(&spans, 42));
        assert_ne!(a, shuffle_within(&spans, 43));
        for &(lo, hi) in &spans {
            let mut part = a[lo..hi].to_vec();
            part.sort();
            assert_eq!(part, (lo..hi).collect::<Vec<_>>());
        }
        let shape = |s: &[usize], lo: usize| s.iter().map(|i| i - lo).collect::<Vec<_>>();
        assert_ne!(shape(&a[0..6], 0), shape(&a[6..12], 6), "equal spans, same permutation");
    }

    #[test]
    fn shuffle_order_is_pinned_across_dependency_upgrades() {
        // 期望值取自换 ChaCha12Rng 之前的 StdRng（rand 0.10.2），换后应逐位相同。
        // 依赖升级改了随机源或洗牌算法时这里会红：同 seed 的历史划分将无法复现，
        // 那是需要用户决定的破坏性改动，不能让它静默发生
        assert_eq!(shuffle_order(20, 666),
            [6, 16, 15, 18, 13, 0, 14, 19, 9, 17, 7, 8, 5, 4, 2, 12, 3, 10, 11, 1]);
        assert_eq!(shuffle_order(20, 42),
            [17, 16, 13, 2, 1, 5, 15, 14, 19, 0, 18, 6, 8, 3, 7, 10, 11, 12, 9, 4]);
    }
}
