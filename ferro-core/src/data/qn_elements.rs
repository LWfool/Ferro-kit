//! Which network formers are conventionally described by a **Qn speciation**.
//!
//! Qn is a tetrahedral-former convention: Q<sup>n</sup> counts the bridging ligands on
//! a site whose coordination is otherwise fixed, so the single digit says everything
//! about how that site is connected.  Al breaks the premise — its coordination varies
//! (Al[4]/Al[5]/Al[6]) and *that* is what the literature reports for it, not a Qn.
//! Listing Al here would produce a `Q^n` column that no aluminophosphate paper quotes.
//!
//! This is a **convention**, not a property of the atom, so it lives in its own table
//! rather than as a field on [`crate::data::elements`] — nothing about boron's mass or
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

/// Whether `elem` is reported with a Qn speciation by default.
pub fn has_qn(elem: &str) -> bool {
    QN_ELEMENTS.contains(&elem)
}

/// [`QN_ELEMENTS`] as an owned set, for use as a [`crate::TypeParams`] default.
pub fn default_qn_set() -> BTreeSet<String> {
    QN_ELEMENTS.iter().map(|s| s.to_string()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_al_is_not_a_qn_element() {
        // Al 由配位数刻画,不由 Qn。这条一旦破了,network_qn.csv 会多出
        // 一个文献里不存在的 Al 行,而 Al 的 label 会从配位数变成桥接数
        assert!(!has_qn("Al"));
        assert!(has_qn("P"));
        assert!(has_qn("Si"));
        assert!(has_qn("B"));
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
