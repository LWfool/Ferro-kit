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

use ferro_core::data::elements::canonical_order;
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
/// A count of 1 drops its subscript (`Al2O4Zn`), unlike `cmd/dataset.rs`'s `formula_of`
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

/// A deterministic shuffle of `0..n`.
pub fn shuffle_order(n: usize, seed: u64) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..n).collect();
    // 不用 StdRng：rand 声明它的算法可能随版本更换，同 seed 的划分就复现不了
    let mut rng = chacha20::ChaCha12Rng::seed_from_u64(seed);
    idx.shuffle(&mut rng);
    idx
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
