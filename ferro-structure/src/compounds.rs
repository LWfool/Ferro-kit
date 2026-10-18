//! Reference data for common compounds.
//!
//! Used by `box_builder` to turn molecule counts into a mass, from which the box size
//! follows at the density the caller gives.
//! Molecular mass unit: g/mol.

pub struct CompoundData {
    /// 常用名称（小写，用于查找）
    pub name: &'static str,
    /// 化学式
    pub formula: &'static str,
    /// 相对分子质量 (g/mol)
    pub molecular_mass: f64,
}

pub static COMPOUNDS: &[CompoundData] = &[
    // 水和醇类
    CompoundData { name: "water",        formula: "H2O",       molecular_mass: 18.015 },
    CompoundData { name: "methanol",     formula: "CH3OH",     molecular_mass: 32.042 },
    CompoundData { name: "ethanol",      formula: "C2H5OH",    molecular_mass: 46.069 },
    CompoundData { name: "isopropanol",  formula: "C3H7OH",    molecular_mass: 60.096 },
    CompoundData { name: "glycerol",     formula: "C3H8O3",    molecular_mass: 92.094 },
    // 酮和醛类
    CompoundData { name: "acetone",      formula: "C3H6O",     molecular_mass: 58.080 },
    CompoundData { name: "acetaldehyde", formula: "C2H4O",     molecular_mass: 44.053 },
    // 腈类
    CompoundData { name: "acetonitrile", formula: "CH3CN",     molecular_mass: 41.053 },
    // 亚砜类
    CompoundData { name: "dmso",         formula: "C2H6OS",    molecular_mass: 78.133 },
    CompoundData { name: "dmf",          formula: "C3H7NO",    molecular_mass: 73.094 },
    // 醚类
    CompoundData { name: "thf",          formula: "C4H8O",     molecular_mass: 72.107 },
    CompoundData { name: "diethylether", formula: "C4H10O",    molecular_mass: 74.123 },
    CompoundData { name: "dme",          formula: "C4H10O2",   molecular_mass: 90.122 },
    // 卤代烃
    CompoundData { name: "chloroform",   formula: "CHCl3",     molecular_mass: 119.38 },
    CompoundData { name: "dcm",          formula: "CH2Cl2",    molecular_mass: 84.933 },
    CompoundData { name: "ccl4",         formula: "CCl4",      molecular_mass: 153.82 },
    // 芳烃
    CompoundData { name: "benzene",      formula: "C6H6",      molecular_mass: 78.114 },
    CompoundData { name: "toluene",      formula: "C7H8",      molecular_mass: 92.141 },
    CompoundData { name: "xylene",       formula: "C8H10",     molecular_mass: 106.17 },
    // 烷烃
    CompoundData { name: "hexane",       formula: "C6H14",     molecular_mass: 86.178 },
    CompoundData { name: "cyclohexane",  formula: "C6H12",     molecular_mass: 84.162 },
    CompoundData { name: "octane",       formula: "C8H18",     molecular_mass: 114.23 },
    // 酸类
    CompoundData { name: "aceticacid",   formula: "C2H4O2",    molecular_mass: 60.052 },
    CompoundData { name: "formicacid",   formula: "CH2O2",     molecular_mass: 46.026 },
    // 常见离子液体前体 / 其他
    CompoundData { name: "ethyleneglycol", formula: "C2H6O2",  molecular_mass: 62.068 },
];

/// 按名称或化学式查找（大小写不敏感）
pub fn find(query: &str) -> Option<&'static CompoundData> {
    let q = query.to_lowercase();
    COMPOUNDS.iter().find(|c| {
        c.name == q.as_str() || c.formula.to_lowercase() == q.as_str()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_water() {
        let w = find("water").unwrap();
        assert_eq!(w.formula, "H2O");
        assert!((w.molecular_mass - 18.015).abs() < 1e-3);
    }

    #[test]
    fn test_find_by_formula() {
        let c = find("CHCl3").unwrap();
        assert_eq!(c.name, "chloroform");
    }
}
