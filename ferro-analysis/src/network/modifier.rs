//! Modifier cation role classification: Free / T / B / M.
//!
//! Each modifier atom is classified by how many NBO (non-bridging ligand)
//! atoms it coordinates within the given cutoff radius:
//!
//! - `Free` — 0 NBO neighbors
//! - `T` (terminal)  — 1 NBO neighbor
//! - `B` (bridging)  — 2 NBO neighbors
//! - `M` (multi)     — ≥3 NBO neighbors

use ferro_core::{Cell, Frame};
use std::collections::{HashMap, HashSet};

/// Build the set of NBO atom indices from the ligand NF map.
///
/// NBO: any ligand atom bonded to exactly 1 network former.
pub fn build_nbo_set(ligand_nf_map: &HashMap<usize, Vec<(String, usize)>>) -> HashSet<usize> {
    ligand_nf_map
        .iter()
        .filter(|(_, nf)| nf.len() == 1)
        .map(|(&idx, _)| idx)
        .collect()
}

/// Classify modifier atoms by their NBO coordination count.
///
/// Returns a `Vec<label>` with one entry per modifier atom in frame order.
pub fn classify_modifiers(
    frame: &Frame,
    cell: &Cell,
    modifier_elem: &str,
    modifier_cutoff: f64,
    nbo_set: &HashSet<usize>,
) -> Vec<String> {
    let cutoff2 = modifier_cutoff * modifier_cutoff;

    frame
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.element == modifier_elem)
        .map(|(mi, ma)| {
            let ma_pos = ma.position;
            let nbo_count = nbo_set
                .iter()
                .filter(|&&nbo_idx| {
                    if nbo_idx == mi {
                        return false;
                    }
                    let diff = cell
                        .minimum_image(frame.atoms[nbo_idx].position - ma_pos)
                        .expect("cell is non-singular");
                    diff.norm_squared() < cutoff2
                })
                .count() as u32;
            role_label(nbo_count)
        })
        .collect()
}

fn role_label(nbo_count: u32) -> String {
    match nbo_count {
        0 => "Free".to_string(),
        1 => "T".to_string(),
        2 => "B".to_string(),
        _ => "M".to_string(),
    }
}

/// Ordering key for modifier role labels: Free < T < B < M.
pub fn modifier_role_order(label: &str) -> u8 {
    match label {
        "Free" => 0,
        "T" => 1,
        "B" => 2,
        "M" => 3,
        _ => 4,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_role_label() {
        assert_eq!(role_label(0), "Free");
        assert_eq!(role_label(1), "T");
        assert_eq!(role_label(2), "B");
        assert_eq!(role_label(3), "M");
        assert_eq!(role_label(5), "M");
    }

    #[test]
    fn test_modifier_role_order() {
        assert!(modifier_role_order("Free") < modifier_role_order("T"));
        assert!(modifier_role_order("T") < modifier_role_order("B"));
        assert!(modifier_role_order("B") < modifier_role_order("M"));
    }

    #[test]
    fn test_build_nbo_set() {
        let mut map: HashMap<usize, Vec<(String, usize)>> = HashMap::new();
        // atom 0: FO (no formers)
        map.insert(0, vec![]);
        // atom 1: NBO (one former)
        map.insert(1, vec![("P".to_string(), 10)]);
        // atom 2: BO (two formers)
        map.insert(2, vec![("P".to_string(), 10), ("P".to_string(), 20)]);

        let nbo = build_nbo_set(&map);
        assert!(nbo.contains(&1));
        assert!(!nbo.contains(&0));
        assert!(!nbo.contains(&2));
    }
}
