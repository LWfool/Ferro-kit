//! 元素静态数据表
//!
//! 包含原子序数 1-86 (H 到 Rn) 的基本属性。
//! 内部单位：质量 amu，电负性 Pauling 标度。

pub struct ElementData {
    pub atomic_number: u8,
    pub symbol: &'static str,
    pub name: &'static str,
    pub atomic_mass: f64,
    /// 常见氧化态（从小到大排列）
    pub common_oxidation_states: &'static [i8],
    /// Pauling 电负性，稀有气体为 None
    pub electronegativity: Option<f64>,
}

/// 元素表，下标 0 = H（atomic_number - 1 = 数组下标）
pub static ELEMENTS: &[ElementData] = &[
    // 第一周期
    ElementData { atomic_number: 1,  symbol: "H",  name: "Hydrogen",     atomic_mass: 1.008,    common_oxidation_states: &[-1, 1],                 electronegativity: Some(2.20) },
    ElementData { atomic_number: 2,  symbol: "He", name: "Helium",       atomic_mass: 4.0026,   common_oxidation_states: &[0],                     electronegativity: None },
    // 第二周期
    ElementData { atomic_number: 3,  symbol: "Li", name: "Lithium",      atomic_mass: 6.941,    common_oxidation_states: &[1],                     electronegativity: Some(0.98) },
    ElementData { atomic_number: 4,  symbol: "Be", name: "Beryllium",    atomic_mass: 9.0122,   common_oxidation_states: &[2],                     electronegativity: Some(1.57) },
    ElementData { atomic_number: 5,  symbol: "B",  name: "Boron",        atomic_mass: 10.811,   common_oxidation_states: &[3],                     electronegativity: Some(2.04) },
    ElementData { atomic_number: 6,  symbol: "C",  name: "Carbon",       atomic_mass: 12.011,   common_oxidation_states: &[-4, -3, -2, -1, 0, 1, 2, 3, 4], electronegativity: Some(2.55) },
    ElementData { atomic_number: 7,  symbol: "N",  name: "Nitrogen",     atomic_mass: 14.007,   common_oxidation_states: &[-3, -2, -1, 1, 2, 3, 4, 5], electronegativity: Some(3.04) },
    ElementData { atomic_number: 8,  symbol: "O",  name: "Oxygen",       atomic_mass: 15.999,   common_oxidation_states: &[-2, -1, 0, 1, 2],       electronegativity: Some(3.44) },
    ElementData { atomic_number: 9,  symbol: "F",  name: "Fluorine",     atomic_mass: 18.998,   common_oxidation_states: &[-1],                    electronegativity: Some(3.98) },
    ElementData { atomic_number: 10, symbol: "Ne", name: "Neon",         atomic_mass: 20.180,   common_oxidation_states: &[0],                     electronegativity: None },
    // 第三周期
    ElementData { atomic_number: 11, symbol: "Na", name: "Sodium",       atomic_mass: 22.990,   common_oxidation_states: &[1],                     electronegativity: Some(0.93) },
    ElementData { atomic_number: 12, symbol: "Mg", name: "Magnesium",    atomic_mass: 24.305,   common_oxidation_states: &[2],                     electronegativity: Some(1.31) },
    ElementData { atomic_number: 13, symbol: "Al", name: "Aluminum",     atomic_mass: 26.982,   common_oxidation_states: &[3],                     electronegativity: Some(1.61) },
    ElementData { atomic_number: 14, symbol: "Si", name: "Silicon",      atomic_mass: 28.086,   common_oxidation_states: &[-4, 4],                 electronegativity: Some(1.90) },
    ElementData { atomic_number: 15, symbol: "P",  name: "Phosphorus",   atomic_mass: 30.974,   common_oxidation_states: &[-3, 3, 5],              electronegativity: Some(2.19) },
    ElementData { atomic_number: 16, symbol: "S",  name: "Sulfur",       atomic_mass: 32.060,   common_oxidation_states: &[-2, 2, 4, 6],           electronegativity: Some(2.58) },
    ElementData { atomic_number: 17, symbol: "Cl", name: "Chlorine",     atomic_mass: 35.450,   common_oxidation_states: &[-1, 1, 3, 5, 7],        electronegativity: Some(3.16) },
    ElementData { atomic_number: 18, symbol: "Ar", name: "Argon",        atomic_mass: 39.948,   common_oxidation_states: &[0],                     electronegativity: None },
    // 第四周期
    ElementData { atomic_number: 19, symbol: "K",  name: "Potassium",    atomic_mass: 39.098,   common_oxidation_states: &[1],                     electronegativity: Some(0.82) },
    ElementData { atomic_number: 20, symbol: "Ca", name: "Calcium",      atomic_mass: 40.078,   common_oxidation_states: &[2],                     electronegativity: Some(1.00) },
    ElementData { atomic_number: 21, symbol: "Sc", name: "Scandium",     atomic_mass: 44.956,   common_oxidation_states: &[3],                     electronegativity: Some(1.36) },
    ElementData { atomic_number: 22, symbol: "Ti", name: "Titanium",     atomic_mass: 47.867,   common_oxidation_states: &[2, 3, 4],               electronegativity: Some(1.54) },
    ElementData { atomic_number: 23, symbol: "V",  name: "Vanadium",     atomic_mass: 50.942,   common_oxidation_states: &[2, 3, 4, 5],            electronegativity: Some(1.63) },
    ElementData { atomic_number: 24, symbol: "Cr", name: "Chromium",     atomic_mass: 51.996,   common_oxidation_states: &[2, 3, 6],               electronegativity: Some(1.66) },
    ElementData { atomic_number: 25, symbol: "Mn", name: "Manganese",    atomic_mass: 54.938,   common_oxidation_states: &[2, 3, 4, 6, 7],         electronegativity: Some(1.55) },
    ElementData { atomic_number: 26, symbol: "Fe", name: "Iron",         atomic_mass: 55.845,   common_oxidation_states: &[2, 3],                  electronegativity: Some(1.83) },
    ElementData { atomic_number: 27, symbol: "Co", name: "Cobalt",       atomic_mass: 58.933,   common_oxidation_states: &[2, 3],                  electronegativity: Some(1.88) },
    ElementData { atomic_number: 28, symbol: "Ni", name: "Nickel",       atomic_mass: 58.693,   common_oxidation_states: &[2, 3],                  electronegativity: Some(1.91) },
    ElementData { atomic_number: 29, symbol: "Cu", name: "Copper",       atomic_mass: 63.546,   common_oxidation_states: &[1, 2],                  electronegativity: Some(1.90) },
    ElementData { atomic_number: 30, symbol: "Zn", name: "Zinc",         atomic_mass: 65.380,   common_oxidation_states: &[2],                     electronegativity: Some(1.65) },
    ElementData { atomic_number: 31, symbol: "Ga", name: "Gallium",      atomic_mass: 69.723,   common_oxidation_states: &[3],                     electronegativity: Some(1.81) },
    ElementData { atomic_number: 32, symbol: "Ge", name: "Germanium",    atomic_mass: 72.630,   common_oxidation_states: &[-4, 2, 4],              electronegativity: Some(2.01) },
    ElementData { atomic_number: 33, symbol: "As", name: "Arsenic",      atomic_mass: 74.922,   common_oxidation_states: &[-3, 3, 5],              electronegativity: Some(2.18) },
    ElementData { atomic_number: 34, symbol: "Se", name: "Selenium",     atomic_mass: 78.971,   common_oxidation_states: &[-2, 2, 4, 6],           electronegativity: Some(2.55) },
    ElementData { atomic_number: 35, symbol: "Br", name: "Bromine",      atomic_mass: 79.904,   common_oxidation_states: &[-1, 1, 3, 5],           electronegativity: Some(2.96) },
    ElementData { atomic_number: 36, symbol: "Kr", name: "Krypton",      atomic_mass: 83.798,   common_oxidation_states: &[0],                     electronegativity: None },
    // 第五周期
    ElementData { atomic_number: 37, symbol: "Rb", name: "Rubidium",     atomic_mass: 85.468,   common_oxidation_states: &[1],                     electronegativity: Some(0.82) },
    ElementData { atomic_number: 38, symbol: "Sr", name: "Strontium",    atomic_mass: 87.620,   common_oxidation_states: &[2],                     electronegativity: Some(0.95) },
    ElementData { atomic_number: 39, symbol: "Y",  name: "Yttrium",      atomic_mass: 88.906,   common_oxidation_states: &[3],                     electronegativity: Some(1.22) },
    ElementData { atomic_number: 40, symbol: "Zr", name: "Zirconium",    atomic_mass: 91.224,   common_oxidation_states: &[4],                     electronegativity: Some(1.33) },
    ElementData { atomic_number: 41, symbol: "Nb", name: "Niobium",      atomic_mass: 92.906,   common_oxidation_states: &[3, 5],                  electronegativity: Some(1.60) },
    ElementData { atomic_number: 42, symbol: "Mo", name: "Molybdenum",   atomic_mass: 95.950,   common_oxidation_states: &[2, 3, 4, 6],            electronegativity: Some(2.16) },
    ElementData { atomic_number: 43, symbol: "Tc", name: "Technetium",   atomic_mass: 98.000,   common_oxidation_states: &[4, 7],                  electronegativity: Some(1.90) },
    ElementData { atomic_number: 44, symbol: "Ru", name: "Ruthenium",    atomic_mass: 101.07,   common_oxidation_states: &[2, 3, 4, 8],            electronegativity: Some(2.20) },
    ElementData { atomic_number: 45, symbol: "Rh", name: "Rhodium",      atomic_mass: 102.91,   common_oxidation_states: &[3],                     electronegativity: Some(2.28) },
    ElementData { atomic_number: 46, symbol: "Pd", name: "Palladium",    atomic_mass: 106.42,   common_oxidation_states: &[2, 4],                  electronegativity: Some(2.20) },
    ElementData { atomic_number: 47, symbol: "Ag", name: "Silver",       atomic_mass: 107.87,   common_oxidation_states: &[1],                     electronegativity: Some(1.93) },
    ElementData { atomic_number: 48, symbol: "Cd", name: "Cadmium",      atomic_mass: 112.41,   common_oxidation_states: &[2],                     electronegativity: Some(1.69) },
    ElementData { atomic_number: 49, symbol: "In", name: "Indium",       atomic_mass: 114.82,   common_oxidation_states: &[3],                     electronegativity: Some(1.78) },
    ElementData { atomic_number: 50, symbol: "Sn", name: "Tin",          atomic_mass: 118.71,   common_oxidation_states: &[2, 4],                  electronegativity: Some(1.96) },
    ElementData { atomic_number: 51, symbol: "Sb", name: "Antimony",     atomic_mass: 121.76,   common_oxidation_states: &[-3, 3, 5],              electronegativity: Some(2.05) },
    ElementData { atomic_number: 52, symbol: "Te", name: "Tellurium",    atomic_mass: 127.60,   common_oxidation_states: &[-2, 2, 4, 6],           electronegativity: Some(2.10) },
    ElementData { atomic_number: 53, symbol: "I",  name: "Iodine",       atomic_mass: 126.90,   common_oxidation_states: &[-1, 1, 3, 5, 7],        electronegativity: Some(2.66) },
    ElementData { atomic_number: 54, symbol: "Xe", name: "Xenon",        atomic_mass: 131.29,   common_oxidation_states: &[0],                     electronegativity: None },
    // 第六周期
    ElementData { atomic_number: 55, symbol: "Cs", name: "Cesium",       atomic_mass: 132.91,   common_oxidation_states: &[1],                     electronegativity: Some(0.79) },
    ElementData { atomic_number: 56, symbol: "Ba", name: "Barium",       atomic_mass: 137.33,   common_oxidation_states: &[2],                     electronegativity: Some(0.89) },
    ElementData { atomic_number: 57, symbol: "La", name: "Lanthanum",    atomic_mass: 138.91,   common_oxidation_states: &[3],                     electronegativity: Some(1.10) },
    ElementData { atomic_number: 58, symbol: "Ce", name: "Cerium",       atomic_mass: 140.12,   common_oxidation_states: &[3, 4],                  electronegativity: Some(1.12) },
    ElementData { atomic_number: 59, symbol: "Pr", name: "Praseodymium", atomic_mass: 140.91,   common_oxidation_states: &[3, 4],                  electronegativity: Some(1.13) },
    ElementData { atomic_number: 60, symbol: "Nd", name: "Neodymium",    atomic_mass: 144.24,   common_oxidation_states: &[3],                     electronegativity: Some(1.14) },
    ElementData { atomic_number: 61, symbol: "Pm", name: "Promethium",   atomic_mass: 145.00,   common_oxidation_states: &[3],                     electronegativity: Some(1.13) },
    ElementData { atomic_number: 62, symbol: "Sm", name: "Samarium",     atomic_mass: 150.36,   common_oxidation_states: &[2, 3],                  electronegativity: Some(1.17) },
    ElementData { atomic_number: 63, symbol: "Eu", name: "Europium",     atomic_mass: 151.96,   common_oxidation_states: &[2, 3],                  electronegativity: Some(1.20) },
    ElementData { atomic_number: 64, symbol: "Gd", name: "Gadolinium",   atomic_mass: 157.25,   common_oxidation_states: &[3],                     electronegativity: Some(1.20) },
    ElementData { atomic_number: 65, symbol: "Tb", name: "Terbium",      atomic_mass: 158.93,   common_oxidation_states: &[3, 4],                  electronegativity: Some(1.10) },
    ElementData { atomic_number: 66, symbol: "Dy", name: "Dysprosium",   atomic_mass: 162.50,   common_oxidation_states: &[3],                     electronegativity: Some(1.22) },
    ElementData { atomic_number: 67, symbol: "Ho", name: "Holmium",      atomic_mass: 164.93,   common_oxidation_states: &[3],                     electronegativity: Some(1.23) },
    ElementData { atomic_number: 68, symbol: "Er", name: "Erbium",       atomic_mass: 167.26,   common_oxidation_states: &[3],                     electronegativity: Some(1.24) },
    ElementData { atomic_number: 69, symbol: "Tm", name: "Thulium",      atomic_mass: 168.93,   common_oxidation_states: &[3],                     electronegativity: Some(1.25) },
    ElementData { atomic_number: 70, symbol: "Yb", name: "Ytterbium",    atomic_mass: 173.04,   common_oxidation_states: &[2, 3],                  electronegativity: Some(1.10) },
    ElementData { atomic_number: 71, symbol: "Lu", name: "Lutetium",     atomic_mass: 174.97,   common_oxidation_states: &[3],                     electronegativity: Some(1.27) },
    ElementData { atomic_number: 72, symbol: "Hf", name: "Hafnium",      atomic_mass: 178.49,   common_oxidation_states: &[4],                     electronegativity: Some(1.30) },
    ElementData { atomic_number: 73, symbol: "Ta", name: "Tantalum",     atomic_mass: 180.95,   common_oxidation_states: &[5],                     electronegativity: Some(1.50) },
    ElementData { atomic_number: 74, symbol: "W",  name: "Tungsten",     atomic_mass: 183.84,   common_oxidation_states: &[4, 6],                  electronegativity: Some(2.36) },
    ElementData { atomic_number: 75, symbol: "Re", name: "Rhenium",      atomic_mass: 186.21,   common_oxidation_states: &[4, 7],                  electronegativity: Some(1.90) },
    ElementData { atomic_number: 76, symbol: "Os", name: "Osmium",       atomic_mass: 190.23,   common_oxidation_states: &[4, 8],                  electronegativity: Some(2.20) },
    ElementData { atomic_number: 77, symbol: "Ir", name: "Iridium",      atomic_mass: 192.22,   common_oxidation_states: &[3, 4],                  electronegativity: Some(2.20) },
    ElementData { atomic_number: 78, symbol: "Pt", name: "Platinum",     atomic_mass: 195.08,   common_oxidation_states: &[2, 4],                  electronegativity: Some(2.28) },
    ElementData { atomic_number: 79, symbol: "Au", name: "Gold",         atomic_mass: 196.97,   common_oxidation_states: &[1, 3],                  electronegativity: Some(2.54) },
    ElementData { atomic_number: 80, symbol: "Hg", name: "Mercury",      atomic_mass: 200.59,   common_oxidation_states: &[1, 2],                  electronegativity: Some(2.00) },
    ElementData { atomic_number: 81, symbol: "Tl", name: "Thallium",     atomic_mass: 204.38,   common_oxidation_states: &[1, 3],                  electronegativity: Some(1.62) },
    ElementData { atomic_number: 82, symbol: "Pb", name: "Lead",         atomic_mass: 207.20,   common_oxidation_states: &[2, 4],                  electronegativity: Some(2.33) },
    ElementData { atomic_number: 83, symbol: "Bi", name: "Bismuth",      atomic_mass: 208.98,   common_oxidation_states: &[3, 5],                  electronegativity: Some(2.02) },
    ElementData { atomic_number: 84, symbol: "Po", name: "Polonium",     atomic_mass: 209.00,   common_oxidation_states: &[2, 4],                  electronegativity: Some(2.00) },
    ElementData { atomic_number: 85, symbol: "At", name: "Astatine",     atomic_mass: 210.00,   common_oxidation_states: &[-1, 1],                 electronegativity: Some(2.20) },
    ElementData { atomic_number: 86, symbol: "Rn", name: "Radon",        atomic_mass: 222.00,   common_oxidation_states: &[0],                     electronegativity: None },
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
