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
//! into one canonical order and permutes the per-atom arrays along with them.
//! A DP model is invariant under atom renumbering, so this changes notation
//! rather than physics — the same thing dpdata does with `sort_atom_names` plus
//! `sort_atom_types`.
//!
//! The canonical order is **(Z, symbol)**, matching the `type_map` that
//! `collect` writes and the grouping order `ferro traj gr` uses. dpdata sorts
//! alphabetically instead; both are self-describing through `type_map.raw`, and
//! being consistent inside ferro matters more than matching dpdata's choice.

use std::collections::BTreeMap;

use ferro_core::data::elements::symbol_to_z;
use ferro_core::Trajectory;
use rand::seq::SliceRandom;
use rand::SeedableRng;

/// Default shuffle seed, so a run without `--seed` is still reproducible.
pub const DEFAULT_SEED: u64 = 666;

/// The permutation putting a frame's atoms into canonical (Z, symbol) order.
///
/// Stable: atoms of one element keep their relative order, so a system already
/// in canonical order is left untouched rather than reshuffled.
pub fn canonical_order(elements: &[String]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..elements.len()).collect();
    idx.sort_by(|&a, &b| {
        let (za, zb) = (symbol_to_z(&elements[a]), symbol_to_z(&elements[b]));
        za.cmp(&zb)
            .then_with(|| elements[a].cmp(&elements[b]))
            .then_with(|| a.cmp(&b))
    });
    idx
}

/// Reorders every frame's atoms — and the per-atom arrays — into canonical order.
///
/// The permutation comes from frame 0 and is applied to every frame, so this is
/// only correct for a trajectory whose atom order is constant. Files that were
/// written independently must each be sorted on their own, *before* being
/// concatenated.
pub fn sort_atoms(traj: &Trajectory) -> Trajectory {
    let Some(first) = traj.frames.first() else {
        return traj.clone();
    };
    let elements: Vec<String> = first.atoms.iter().map(|a| a.element.clone()).collect();
    let order = canonical_order(&elements);
    if order.iter().enumerate().all(|(i, &j)| i == j) {
        return traj.clone();
    }
    // bonds 存的是原子编号，跟着置换走要的是反向映射：old -> new
    let mut new_of = vec![0usize; order.len()];
    for (new, &old) in order.iter().enumerate() {
        new_of[old] = new;
    }

    let mut out = traj.clone();
    for f in &mut out.frames {
        let atoms: Vec<_> = order.iter().map(|&i| f.atoms[i].clone()).collect();
        f.atoms = atoms;
        // 逐原子的旁路数组必须跟着同一个置换走；能量/盒子/应力与原子编号无关。
        // Frame 上共三个：forces、velocities 与 bonds —— 漏掉任何一个都是静默的
        if let Some(forces) = &f.forces {
            f.forces = Some(order.iter().map(|&i| forces[i]).collect());
        }
        if let Some(velocities) = &f.velocities {
            f.velocities = Some(order.iter().map(|&i| velocities[i]).collect());
        }
        if let Some(bonds) = &f.bonds {
            f.bonds = Some(bonds.iter().map(|&(a, b)| (new_of[a], new_of[b])).collect());
        }
    }
    out
}

/// The key two systems must share to be merged: their sorted element sequence.
pub fn composition_key(traj: &Trajectory) -> Vec<String> {
    let mut v: Vec<String> = traj
        .frames
        .first()
        .map(|f| f.atoms.iter().map(|a| a.element.clone()).collect())
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
        for a in &f.atoms {
            *counts.entry(a.element.clone()).or_default() += 1;
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
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
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
        f.atoms = elements
            .iter()
            .enumerate()
            .map(|(i, e)| Atom::new(*e, Vector3::new(i as f64, 0.0, 0.0)))
            .collect();
        f.forces = Some(
            (0..elements.len())
                .map(|i| Vector3::new(i as f64 * 0.1, 0.0, 0.0))
                .collect(),
        );
        Trajectory { frames: vec![f], metadata: Default::default() }
    }

    #[test]
    fn canonical_order_is_by_z_then_symbol() {
        // O(8) < Al(13) < Zn(30)
        let t = traj_of(&["Zn", "Al", "O", "Al"]);
        assert_eq!(composition_key(&t), vec!["O", "Al", "Al", "Zn"]);
    }

    #[test]
    fn sorting_carries_the_per_atom_arrays_along() {
        let t = traj_of(&["Zn", "O"]);
        let s = sort_atoms(&t);
        let f = &s.frames[0];
        assert_eq!(f.atoms[0].element, "O");
        // O 原来在索引 1，位置与力都必须跟过来
        assert_eq!(f.atoms[0].position.x, 1.0);
        assert!((f.forces.as_ref().unwrap()[0].x - 0.1).abs() < 1e-12);
        assert_eq!(f.atoms[1].element, "Zn");
        assert_eq!(f.atoms[1].position.x, 0.0);
    }

    /// Velocities and bonds are per-atom too, and were once left behind.
    ///
    /// A DeePMD system carries neither, so `merge` never showed the bug; the
    /// diagnostics `collect` exports (a dump with velocities) do.
    #[test]
    fn sorting_carries_velocities_and_bonds_too() {
        let mut t = traj_of(&["Zn", "O", "Al"]);
        let f = &mut t.frames[0];
        f.velocities = Some(
            (0..3).map(|i| Vector3::new(i as f64 * 10.0, 0.0, 0.0)).collect(),
        );
        // Zn(0)-O(1) 与 O(1)-Al(2)；排序后 O->0, Al->1, Zn->2
        f.bonds = Some(vec![(0, 1), (1, 2)]);

        let s = sort_atoms(&t);
        let f = &s.frames[0];
        let v = f.velocities.as_ref().unwrap();
        assert_eq!(f.atoms[0].element, "O");
        assert!((v[0].x - 10.0).abs() < 1e-12, "O 的速度必须跟过来");
        assert!((v[2].x - 0.0).abs() < 1e-12, "Zn 的速度也是");
        // 键连的是原子，不是位置：两端都要重编号。端点的先后不动 —— 规范成
        // a < b 是另一件事，这里只做重编号
        assert_eq!(f.bonds.clone().unwrap(), vec![(2, 0), (0, 1)]);
    }

    #[test]
    fn an_already_sorted_system_is_untouched() {
        let t = traj_of(&["O", "O", "Al"]);
        let s = sort_atoms(&t);
        let before: Vec<f64> = t.frames[0].atoms.iter().map(|a| a.position.x).collect();
        let after: Vec<f64> = s.frames[0].atoms.iter().map(|a| a.position.x).collect();
        assert_eq!(before, after);
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
}
