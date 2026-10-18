//! Reference data for common compounds.
//!
//! Used by `ferro-structure` to estimate the initial simulation box size from molecule counts and density.
//! Density unit: g/cm³ (standard conditions, 25 °C); `None` for gases.
//! Molecular mass unit: g/mol.

pub struct CompoundData {
    /// 常用名称（小写，用于查找）
    pub name: &'static str,
    /// 化学式
    pub formula: &'static str,
    /// 相对分子质量 (g/mol)
    pub molecular_mass: f64,
    /// 密度 (g/cm³)，气体为 None
    pub density: Option<f64>,
    /// CAS 号
    pub cas: Option<&'static str>,
}

pub static COMPOUNDS: &[CompoundData] = &[
    // 水和醇类
    CompoundData { name: "water",        formula: "H2O",       molecular_mass: 18.015,  density: Some(0.997),  cas: Some("7732-18-5")  },
    CompoundData { name: "methanol",     formula: "CH3OH",     molecular_mass: 32.042,  density: Some(0.791),  cas: Some("67-56-1")    },
    CompoundData { name: "ethanol",      formula: "C2H5OH",    molecular_mass: 46.069,  density: Some(0.789),  cas: Some("64-17-5")    },
    CompoundData { name: "isopropanol",  formula: "C3H7OH",    molecular_mass: 60.096,  density: Some(0.786),  cas: Some("67-63-0")    },
    CompoundData { name: "glycerol",     formula: "C3H8O3",    molecular_mass: 92.094,  density: Some(1.261),  cas: Some("56-81-5")    },
    // 酮和醛类
    CompoundData { name: "acetone",      formula: "C3H6O",     molecular_mass: 58.080,  density: Some(0.791),  cas: Some("67-64-1")    },
    CompoundData { name: "acetaldehyde", formula: "C2H4O",     molecular_mass: 44.053,  density: Some(0.788),  cas: Some("75-07-0")    },
    // 腈类
    CompoundData { name: "acetonitrile", formula: "CH3CN",     molecular_mass: 41.053,  density: Some(0.786),  cas: Some("75-05-8")    },
    // 亚砜类
    CompoundData { name: "dmso",         formula: "C2H6OS",    molecular_mass: 78.133,  density: Some(1.100),  cas: Some("67-68-5")    },
    CompoundData { name: "dmf",          formula: "C3H7NO",    molecular_mass: 73.094,  density: Some(0.948),  cas: Some("68-12-2")    },
    // 醚类
    CompoundData { name: "thf",          formula: "C4H8O",     molecular_mass: 72.107,  density: Some(0.889),  cas: Some("109-99-9")   },
    CompoundData { name: "diethylether", formula: "C4H10O",    molecular_mass: 74.123,  density: Some(0.713),  cas: Some("60-29-7")    },
    CompoundData { name: "dme",          formula: "C4H10O2",   molecular_mass: 90.122,  density: Some(0.867),  cas: Some("110-71-4")   },
    // 卤代烃
    CompoundData { name: "chloroform",   formula: "CHCl3",     molecular_mass: 119.38,  density: Some(1.489),  cas: Some("67-66-3")    },
    CompoundData { name: "dcm",          formula: "CH2Cl2",    molecular_mass: 84.933,  density: Some(1.325),  cas: Some("75-09-2")    },
    CompoundData { name: "ccl4",         formula: "CCl4",      molecular_mass: 153.82,  density: Some(1.594),  cas: Some("56-23-5")    },
    // 芳烃
    CompoundData { name: "benzene",      formula: "C6H6",      molecular_mass: 78.114,  density: Some(0.879),  cas: Some("71-43-2")    },
    CompoundData { name: "toluene",      formula: "C7H8",      molecular_mass: 92.141,  density: Some(0.867),  cas: Some("108-88-3")   },
    CompoundData { name: "xylene",       formula: "C8H10",     molecular_mass: 106.17,  density: Some(0.864),  cas: Some("1330-20-7")  },
    // 烷烃
    CompoundData { name: "hexane",       formula: "C6H14",     molecular_mass: 86.178,  density: Some(0.659),  cas: Some("110-54-3")   },
    CompoundData { name: "cyclohexane",  formula: "C6H12",     molecular_mass: 84.162,  density: Some(0.779),  cas: Some("110-82-7")   },
    CompoundData { name: "octane",       formula: "C8H18",     molecular_mass: 114.23,  density: Some(0.703),  cas: Some("111-65-9")   },
    // 酸类
    CompoundData { name: "aceticacid",   formula: "C2H4O2",    molecular_mass: 60.052,  density: Some(1.049),  cas: Some("64-19-7")    },
    CompoundData { name: "formicacid",   formula: "CH2O2",     molecular_mass: 46.026,  density: Some(1.220),  cas: Some("64-18-6")    },
    // 常见离子液体前体 / 其他
    CompoundData { name: "ethyleneglycol", formula: "C2H6O2",  molecular_mass: 62.068,  density: Some(1.113),  cas: Some("107-21-1")   },
];

/// 按名称或化学式查找（大小写不敏感）
pub fn find(query: &str) -> Option<&'static CompoundData> {
    let q = query.to_lowercase();
    COMPOUNDS.iter().find(|c| {
        c.name == q.as_str() || c.formula.to_lowercase() == q.as_str()
    })
}

/// 返回所有有密度数据的化合物
pub fn with_density() -> impl Iterator<Item = &'static CompoundData> {
    COMPOUNDS.iter().filter(|c| c.density.is_some())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_water() {
        let w = find("water").unwrap();
        assert_eq!(w.formula, "H2O");
        assert!((w.molecular_mass - 18.015).abs() < 1e-3);
        assert!((w.density.unwrap() - 0.997).abs() < 1e-3);
    }

    #[test]
    fn test_find_by_formula() {
        let c = find("CHCl3").unwrap();
        assert_eq!(c.name, "chloroform");
    }
}
