//! Which network formers are conventionally described by a **Qn speciation**.
//!
//! Qn is a tetrahedral-former convention: Q<sup>n</sup> counts the bridging ligands on
//! a site whose coordination is otherwise fixed, so the single digit says everything
//! about how that site is connected.  Al breaks the premise — its coordination varies
//! (Al[4]/Al[5]/Al[6]) and *that* is what the literature reports for it, not a Qn.
//! Listing Al here would produce a `Q^n` column that no aluminophosphate paper quotes.
//!
//! This is a **convention**, not a property of the atom, so it lives in its own table
//! rather than as a field on [`ferro_core::data::elements`] — nothing about boron's mass or
//! electron configuration implies it gets a Qn.
//!
//! The list is a *default*.  `ferro net --qn <ELEM,…>` replaces it wholesale, because
//! whether an element acts as a Qn former is a property of the system being modelled:
//! in an aluminosilicate one may well want Al reported as Q<sup>n</sup>(mSi), and in a
//! borophosphate one may want B left out.

use std::collections::BTreeSet;

/// Formers reported with a Qn speciation by default, sorted.
///
/// Deliberately short.  An element earns a place here only when the Qn notation for it
/// is standard in the glass literature *and* its coordination is effectively fixed —
/// the two conditions that make a single digit sufficient.
pub const QN_ELEMENTS: &[&str] = &["B", "P", "Si"];

/// [`QN_ELEMENTS`] as an owned set, for use as a [`super::TypeParams`] default.
pub fn default_qn_set() -> BTreeSet<String> {
    QN_ELEMENTS.iter().map(|s| s.to_string()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_al_is_not_a_qn_element() {
        // Al is described by its coordination, not by Qn. Should this break,
        // network_qn.csv gains an Al row no paper reports, and the Al label turns
        // from a coordination number into a bridging count
        assert!(!QN_ELEMENTS.contains(&"Al"));
        for e in ["P", "Si", "B"] {
            assert!(QN_ELEMENTS.contains(&e), "{e}");
        }
    }

    #[test]
    fn test_list_is_sorted_and_unique() {
        // BTreeSet 的构造依赖不了它,但排序过的常量表读起来才能一眼看全
        let mut sorted = QN_ELEMENTS.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted, QN_ELEMENTS);
        assert_eq!(default_qn_set().len(), QN_ELEMENTS.len());
    }
}
