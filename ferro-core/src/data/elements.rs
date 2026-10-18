//! Static element table: atomic number, symbol and atomic mass (amu) of
//! elements 1-86 (H to Rn).

pub struct ElementData {
    pub atomic_number: u8,
    pub symbol: &'static str,
    pub atomic_mass: f64,
}

/// 元素表，下标 0 = H（atomic_number - 1 = 数组下标）
pub static ELEMENTS: &[ElementData] = &[
    // 第一周期
    ElementData { atomic_number: 1,  symbol: "H",     atomic_mass: 1.008 },
    ElementData { atomic_number: 2,  symbol: "He",       atomic_mass: 4.0026 },
    // 第二周期
    ElementData { atomic_number: 3,  symbol: "Li",      atomic_mass: 6.941 },
    ElementData { atomic_number: 4,  symbol: "Be",    atomic_mass: 9.0122 },
    ElementData { atomic_number: 5,  symbol: "B",        atomic_mass: 10.811 },
    ElementData { atomic_number: 6,  symbol: "C",       atomic_mass: 12.011 },
    ElementData { atomic_number: 7,  symbol: "N",     atomic_mass: 14.007 },
    ElementData { atomic_number: 8,  symbol: "O",       atomic_mass: 15.999 },
    ElementData { atomic_number: 9,  symbol: "F",     atomic_mass: 18.998 },
    ElementData { atomic_number: 10, symbol: "Ne",         atomic_mass: 20.180 },
    // 第三周期
    ElementData { atomic_number: 11, symbol: "Na",       atomic_mass: 22.990 },
    ElementData { atomic_number: 12, symbol: "Mg",    atomic_mass: 24.305 },
    ElementData { atomic_number: 13, symbol: "Al",     atomic_mass: 26.982 },
    ElementData { atomic_number: 14, symbol: "Si",      atomic_mass: 28.086 },
    ElementData { atomic_number: 15, symbol: "P",   atomic_mass: 30.974 },
    ElementData { atomic_number: 16, symbol: "S",       atomic_mass: 32.060 },
    ElementData { atomic_number: 17, symbol: "Cl",     atomic_mass: 35.450 },
    ElementData { atomic_number: 18, symbol: "Ar",        atomic_mass: 39.948 },
    // 第四周期
    ElementData { atomic_number: 19, symbol: "K",    atomic_mass: 39.098 },
    ElementData { atomic_number: 20, symbol: "Ca",      atomic_mass: 40.078 },
    ElementData { atomic_number: 21, symbol: "Sc",     atomic_mass: 44.956 },
    ElementData { atomic_number: 22, symbol: "Ti",     atomic_mass: 47.867 },
    ElementData { atomic_number: 23, symbol: "V",     atomic_mass: 50.942 },
    ElementData { atomic_number: 24, symbol: "Cr",     atomic_mass: 51.996 },
    ElementData { atomic_number: 25, symbol: "Mn",    atomic_mass: 54.938 },
    ElementData { atomic_number: 26, symbol: "Fe",         atomic_mass: 55.845 },
    ElementData { atomic_number: 27, symbol: "Co",       atomic_mass: 58.933 },
    ElementData { atomic_number: 28, symbol: "Ni",       atomic_mass: 58.693 },
    ElementData { atomic_number: 29, symbol: "Cu",       atomic_mass: 63.546 },
    ElementData { atomic_number: 30, symbol: "Zn",         atomic_mass: 65.380 },
    ElementData { atomic_number: 31, symbol: "Ga",      atomic_mass: 69.723 },
    ElementData { atomic_number: 32, symbol: "Ge",    atomic_mass: 72.630 },
    ElementData { atomic_number: 33, symbol: "As",      atomic_mass: 74.922 },
    ElementData { atomic_number: 34, symbol: "Se",     atomic_mass: 78.971 },
    ElementData { atomic_number: 35, symbol: "Br",      atomic_mass: 79.904 },
    ElementData { atomic_number: 36, symbol: "Kr",      atomic_mass: 83.798 },
    // 第五周期
    ElementData { atomic_number: 37, symbol: "Rb",     atomic_mass: 85.468 },
    ElementData { atomic_number: 38, symbol: "Sr",    atomic_mass: 87.620 },
    ElementData { atomic_number: 39, symbol: "Y",      atomic_mass: 88.906 },
    ElementData { atomic_number: 40, symbol: "Zr",    atomic_mass: 91.224 },
    ElementData { atomic_number: 41, symbol: "Nb",      atomic_mass: 92.906 },
    ElementData { atomic_number: 42, symbol: "Mo",   atomic_mass: 95.950 },
    ElementData { atomic_number: 43, symbol: "Tc",   atomic_mass: 98.000 },
    ElementData { atomic_number: 44, symbol: "Ru",    atomic_mass: 101.07 },
    ElementData { atomic_number: 45, symbol: "Rh",      atomic_mass: 102.91 },
    ElementData { atomic_number: 46, symbol: "Pd",    atomic_mass: 106.42 },
    ElementData { atomic_number: 47, symbol: "Ag",       atomic_mass: 107.87 },
    ElementData { atomic_number: 48, symbol: "Cd",      atomic_mass: 112.41 },
    ElementData { atomic_number: 49, symbol: "In",       atomic_mass: 114.82 },
    ElementData { atomic_number: 50, symbol: "Sn",          atomic_mass: 118.71 },
    ElementData { atomic_number: 51, symbol: "Sb",     atomic_mass: 121.76 },
    ElementData { atomic_number: 52, symbol: "Te",    atomic_mass: 127.60 },
    ElementData { atomic_number: 53, symbol: "I",       atomic_mass: 126.90 },
    ElementData { atomic_number: 54, symbol: "Xe",        atomic_mass: 131.29 },
    // 第六周期
    ElementData { atomic_number: 55, symbol: "Cs",       atomic_mass: 132.91 },
    ElementData { atomic_number: 56, symbol: "Ba",       atomic_mass: 137.33 },
    ElementData { atomic_number: 57, symbol: "La",    atomic_mass: 138.91 },
    ElementData { atomic_number: 58, symbol: "Ce",       atomic_mass: 140.12 },
    ElementData { atomic_number: 59, symbol: "Pr", atomic_mass: 140.91 },
    ElementData { atomic_number: 60, symbol: "Nd",    atomic_mass: 144.24 },
    ElementData { atomic_number: 61, symbol: "Pm",   atomic_mass: 145.00 },
    ElementData { atomic_number: 62, symbol: "Sm",     atomic_mass: 150.36 },
    ElementData { atomic_number: 63, symbol: "Eu",     atomic_mass: 151.96 },
    ElementData { atomic_number: 64, symbol: "Gd",   atomic_mass: 157.25 },
    ElementData { atomic_number: 65, symbol: "Tb",      atomic_mass: 158.93 },
    ElementData { atomic_number: 66, symbol: "Dy",   atomic_mass: 162.50 },
    ElementData { atomic_number: 67, symbol: "Ho",      atomic_mass: 164.93 },
    ElementData { atomic_number: 68, symbol: "Er",       atomic_mass: 167.26 },
    ElementData { atomic_number: 69, symbol: "Tm",      atomic_mass: 168.93 },
    ElementData { atomic_number: 70, symbol: "Yb",    atomic_mass: 173.04 },
    ElementData { atomic_number: 71, symbol: "Lu",     atomic_mass: 174.97 },
    ElementData { atomic_number: 72, symbol: "Hf",      atomic_mass: 178.49 },
    ElementData { atomic_number: 73, symbol: "Ta",     atomic_mass: 180.95 },
    ElementData { atomic_number: 74, symbol: "W",     atomic_mass: 183.84 },
    ElementData { atomic_number: 75, symbol: "Re",      atomic_mass: 186.21 },
    ElementData { atomic_number: 76, symbol: "Os",       atomic_mass: 190.23 },
    ElementData { atomic_number: 77, symbol: "Ir",      atomic_mass: 192.22 },
    ElementData { atomic_number: 78, symbol: "Pt",     atomic_mass: 195.08 },
    ElementData { atomic_number: 79, symbol: "Au",         atomic_mass: 196.97 },
    ElementData { atomic_number: 80, symbol: "Hg",      atomic_mass: 200.59 },
    ElementData { atomic_number: 81, symbol: "Tl",     atomic_mass: 204.38 },
    ElementData { atomic_number: 82, symbol: "Pb",         atomic_mass: 207.20 },
    ElementData { atomic_number: 83, symbol: "Bi",      atomic_mass: 208.98 },
    ElementData { atomic_number: 84, symbol: "Po",     atomic_mass: 209.00 },
    ElementData { atomic_number: 85, symbol: "At",     atomic_mass: 210.00 },
    ElementData { atomic_number: 86, symbol: "Rn",        atomic_mass: 222.00 },
];

/// 按元素符号查找，O(n)
pub fn by_symbol(symbol: &str) -> Option<&'static ElementData> {
    ELEMENTS.iter().find(|e| e.symbol == symbol)
}

/// 按原子序数查找，O(1)（atomic_number 从 1 开始）
pub fn by_number(n: u8) -> Option<&'static ElementData> {
    if n == 0 || n as usize > ELEMENTS.len() {
        return None;
    }
    Some(&ELEMENTS[(n - 1) as usize])
}

/// The canonical element order: by atomic number, then by the string.
///
/// The one rule behind the DeePMD `type_map` (`ferro-io`), the partial g(r)
/// column order and the atom order `merge` sorts into (review B12): a dataset
/// written by `collect` and one sorted by `merge` must agree on it, so it lives
/// in one place. The string tiebreak orders site labels sharing an element
/// (`P_3` < `P_4`) and keeps unknown symbols (Z = 255) deterministic.
pub fn canonical_cmp(a: &str, b: &str) -> std::cmp::Ordering {
    (symbol_to_z(a), a).cmp(&(symbol_to_z(b), b))
}

/// The permutation putting `elements` into canonical order ([`canonical_cmp`]).
///
/// Stable: atoms of one element keep their relative order, so a system already
/// in canonical order is left untouched rather than reshuffled.
pub fn canonical_order(elements: &[String]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..elements.len()).collect();
    idx.sort_by(|&a, &b| canonical_cmp(&elements[a], &elements[b]).then_with(|| a.cmp(&b)));
    idx
}

/// Look up the atomic number for a symbol, with pseudo-element label fallback.
///
/// Matching strategy (tried in order):
///   1. Exact match ("Fe" → 26, "O" → 8).
///   2. First two bytes (uppercase + lowercase) as a chemical symbol ("Fe1" → "Fe" → 26).
///   3. First byte (uppercase letter) as a single-character element ("Onb" → "O" → 8).
///   4. Unrecognised → 255.
pub fn symbol_to_z(symbol: &str) -> u8 {
    if let Some(e) = by_symbol(symbol) { return e.atomic_number; }
    let b = symbol.as_bytes();
    if b.len() >= 2 && b[0].is_ascii_uppercase() && b[1].is_ascii_lowercase() {
        if let Some(e) = by_symbol(&symbol[..2]) { return e.atomic_number; }
    }
    if !b.is_empty() && b[0].is_ascii_uppercase() {
        if let Some(e) = by_symbol(&symbol[..1]) { return e.atomic_number; }
    }
    255
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_by_symbol() {
        let fe = by_symbol("Fe").unwrap();
        assert_eq!(fe.atomic_number, 26);
        assert!((fe.atomic_mass - 55.845).abs() < 1e-3);
    }

    #[test]
    fn test_by_number() {
        let o = by_number(8).unwrap();
        assert_eq!(o.symbol, "O");
    }

    #[test]
    fn test_index_consistency() {
        for (i, e) in ELEMENTS.iter().enumerate() {
            assert_eq!(e.atomic_number as usize, i + 1, "Element {} index mismatch", e.symbol);
        }
    }
}
