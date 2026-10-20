//! Unit conversion utilities.
//!
//! Internal standard units (DeePMD-kit / VASP convention):
//! - Length: Å
//! - Energy: eV
//! - Force: eV/Å
//! - Stress: eV/Å³
//! - Time: fs
//! - Mass: amu
//! - Charge: e
//! - Temperature: K

// ── 物理常数 ──────────────────────────────────────────────────────────────────

/// Bohr 半径，单位 Å
pub const BOHR_TO_ANG: f64 = 0.529_177_210_9;
pub const ANG_TO_BOHR: f64 = 1.0 / BOHR_TO_ANG;

/// Hartree → eV
pub const HARTREE_TO_EV: f64 = 27.211_386_245_988;
pub const EV_TO_HARTREE: f64 = 1.0 / HARTREE_TO_EV;

/// eV → kcal/mol, CODATA 2018: e·N_A / 4184 J (thermochemical calorie), both exact.
/// Was 23.0605419, 2.6e-7 off; its source is unknown, no CODATA set gives it (2014 gives
/// 23.0605480, 2010 23.0605489).
pub const EV_TO_KCAL_MOL: f64 = 23.060_547_830_619_03;
pub const KCAL_MOL_TO_EV: f64 = 1.0 / EV_TO_KCAL_MOL;

/// eV → kJ/mol
pub const EV_TO_KJ_MOL: f64 = 96.485_332_1;
pub const KJ_MOL_TO_EV: f64 = 1.0 / EV_TO_KJ_MOL;

/// eV → cm⁻¹（波数）
pub const EV_TO_WAVENUMBER: f64 = 8_065.544_0;
pub const WAVENUMBER_TO_EV: f64 = 1.0 / EV_TO_WAVENUMBER;

/// eV/Å³ → GPa
pub const EV_ANG3_TO_GPA: f64 = 160.217_663_4;
pub const GPA_TO_EV_ANG3: f64 = 1.0 / EV_ANG3_TO_GPA;

/// GPa → bar（CP2K 的 `STRESS|` 默认单位）
pub const GPA_TO_BAR: f64 = 1e4;
pub const BAR_TO_GPA: f64 = 1e-4;

/// GPa → kBar
pub const GPA_TO_KBAR: f64 = 10.0;
pub const KBAR_TO_GPA: f64 = 0.1;

/// Avogadro 常数
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// 质量密度：amu/Å³ → g/cm³。
///
/// 1 amu = 1/N_A g 且 1 Å³ = 1e-24 cm³，两者之比即 `1e24 / AVOGADRO`。
pub const AMU_ANG3_TO_G_CM3: f64 = 1e24 / AVOGADRO;

/// Diffusion coefficient: Å²/fs → cm²/s (1 Å² = 1e-16 cm², 1 fs = 1e-15 s).
pub const ANG2_FS_TO_CM2_S: f64 = 0.1;
/// Diffusion coefficient: Å²/fs → m²/s (1 Å² = 1e-20 m², 1 fs = 1e-15 s).
pub const ANG2_FS_TO_M2_S: f64 = 1e-5;

/// Boltzmann 常数（eV/K）。
///
/// CODATA 2018 下 k_B = 1.380649e-23 J/K 与 e = 1.602176634e-19 C 都是精确值，
/// 故这个比值也是精确的。用于从离子动能反算瞬时温度：`T = 2·E_kin/(3N·k_B)`。
pub const BOLTZMANN_EV_K: f64 = 8.617_333_262e-5;

/// 1 fs 折合多少个 ASE 时间单位（Å·√(amu/eV)），即 `ase.units.fs`。
///
/// ASE 的速度是 Å/该单位、`momenta` 是 amu·Å/该单位，换成 Å/fs 要乘这个数。取 ASE 默认的
/// CODATA 2014 下的值（ASE 3.29 实测），不用本文件的 2018 常数自己算 —— 要与 ASE 写出的
/// `momenta` 逐位对上（两版差 1e-9 量级）
pub const ASE_FS: f64 = 0.098_226_947_884_640_63;

// ── 长度 ──────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LengthUnit {
    Angstrom,
    Bohr,
    Nanometer,
    Picometer,
}

impl LengthUnit {
    /// The unit a program's output or input names, case-insensitively: `angstrom`,
    /// `bohr`. What a name means is fixed; where it sits in a line is each reader's.
    pub fn from_name(name: &str) -> Option<Self> {
        match name.to_ascii_lowercase().as_str() {
            "angstrom" => Some(Self::Angstrom),
            "bohr" => Some(Self::Bohr),
            _ => None,
        }
    }

    fn to_angstrom(self, v: f64) -> f64 {
        match self {
            Self::Angstrom  => v,
            Self::Bohr      => v * BOHR_TO_ANG,
            Self::Nanometer => v * 10.0,
            Self::Picometer => v * 0.01,
        }
    }
    fn ang_into(self, v: f64) -> f64 {
        match self {
            Self::Angstrom  => v,
            Self::Bohr      => v * ANG_TO_BOHR,
            Self::Nanometer => v * 0.1,
            Self::Picometer => v * 100.0,
        }
    }
}

pub fn convert_length(value: f64, from: LengthUnit, to: LengthUnit) -> f64 {
    to.ang_into(from.to_angstrom(value))
}

// ── 能量 ──────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyUnit {
    EV,
    Hartree,
    KcalPerMol,
    KJPerMol,
    Wavenumber,
}

impl EnergyUnit {
    /// The unit a program's output names, case-insensitively: `hartree` (also
    /// `a.u.` / `au`, the atomic unit of energy, as CP2K wrote it before 2025), `ev`.
    pub fn from_name(name: &str) -> Option<Self> {
        match name.to_ascii_lowercase().as_str() {
            "hartree" | "a.u." | "au" => Some(Self::Hartree),
            "ev" => Some(Self::EV),
            _ => None,
        }
    }

    fn to_ev(self, v: f64) -> f64 {
        match self {
            Self::EV         => v,
            Self::Hartree    => v * HARTREE_TO_EV,
            Self::KcalPerMol => v * KCAL_MOL_TO_EV,
            Self::KJPerMol   => v * KJ_MOL_TO_EV,
            Self::Wavenumber => v * WAVENUMBER_TO_EV,
        }
    }
    fn ev_into(self, v: f64) -> f64 {
        match self {
            Self::EV         => v,
            Self::Hartree    => v * EV_TO_HARTREE,
            Self::KcalPerMol => v * EV_TO_KCAL_MOL,
            Self::KJPerMol   => v * EV_TO_KJ_MOL,
            Self::Wavenumber => v * EV_TO_WAVENUMBER,
        }
    }
}

pub fn convert_energy(value: f64, from: EnergyUnit, to: EnergyUnit) -> f64 {
    to.ev_into(from.to_ev(value))
}

// ── 压强/应力 ─────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PressureUnit {
    /// eV/Å³（内部单位）
    EVPerAng3,
    GPa,
    Kbar,
    /// CP2K 的 `STRESS|` 默认单位
    Bar,
}

impl PressureUnit {
    /// The unit a program's output names, case-insensitively: `bar`, `gpa`, `kbar`.
    pub fn from_name(name: &str) -> Option<Self> {
        match name.to_ascii_lowercase().as_str() {
            "bar" => Some(Self::Bar),
            "gpa" => Some(Self::GPa),
            "kbar" => Some(Self::Kbar),
            _ => None,
        }
    }

    fn to_ev_ang3(self, v: f64) -> f64 {
        match self {
            Self::EVPerAng3 => v,
            Self::GPa       => v * GPA_TO_EV_ANG3,
            Self::Kbar      => v * KBAR_TO_GPA * GPA_TO_EV_ANG3,
            Self::Bar       => v * BAR_TO_GPA * GPA_TO_EV_ANG3,
        }
    }
    fn ev_ang3_into(self, v: f64) -> f64 {
        match self {
            Self::EVPerAng3 => v,
            Self::GPa       => v * EV_ANG3_TO_GPA,
            Self::Kbar      => v * EV_ANG3_TO_GPA * GPA_TO_KBAR,
            Self::Bar       => v * EV_ANG3_TO_GPA * GPA_TO_BAR,
        }
    }
}

pub fn convert_pressure(value: f64, from: PressureUnit, to: PressureUnit) -> f64 {
    to.ev_ang3_into(from.to_ev_ang3(value))
}

// ── 时间 ──────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeUnit {
    Femtosecond,
    Picosecond,
}

impl TimeUnit {
    fn to_fs(self, v: f64) -> f64 {
        match self {
            Self::Femtosecond => v,
            Self::Picosecond  => v * 1000.0,
        }
    }
    fn fs_into(self, v: f64) -> f64 {
        match self {
            Self::Femtosecond => v,
            Self::Picosecond  => v * 0.001,
        }
    }
}

pub fn convert_time(value: f64, from: TimeUnit, to: TimeUnit) -> f64 {
    to.fs_into(from.to_fs(value))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_density_constant_against_known_solids() {
        // 用两个手册值反查：NaCl 立方胞 a=5.6402 Å，Z=4（4 个 NaCl 单元），
        // 实验密度 2.165 g/cm3
        let m_nacl = 4.0 * (22.990 + 35.453);
        let v_nacl = 5.6402_f64.powi(3);
        let rho = m_nacl / v_nacl * AMU_ANG3_TO_G_CM3;
        assert!((rho - 2.165).abs() < 5e-3, "NaCl rho = {rho}");

        // 金刚石 a=3.5670 Å，Z=8 个 C，实验密度 3.515 g/cm3
        let rho_c = 8.0 * 12.011 / 3.5670_f64.powi(3) * AMU_ANG3_TO_G_CM3;
        assert!((rho_c - 3.515).abs() < 5e-3, "diamond rho = {rho_c}");
    }

    #[test]
    fn test_length() {
        let bohr = 1.0;
        let ang = convert_length(bohr, LengthUnit::Bohr, LengthUnit::Angstrom);
        assert!((ang - BOHR_TO_ANG).abs() < 1e-10);

        let back = convert_length(ang, LengthUnit::Angstrom, LengthUnit::Bohr);
        assert!((back - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_energy_hartree_ev() {
        let ev = convert_energy(1.0, EnergyUnit::Hartree, EnergyUnit::EV);
        assert!((ev - HARTREE_TO_EV).abs() < 1e-6);
    }

    #[test]
    fn test_pressure_roundtrip() {
        let gpa = 10.0;
        let internal = convert_pressure(gpa, PressureUnit::GPa, PressureUnit::EVPerAng3);
        let back = convert_pressure(internal, PressureUnit::EVPerAng3, PressureUnit::GPa);
        assert!((back - gpa).abs() < 1e-9);
    }

    #[test]
    fn test_time() {
        let ps = 1.0;
        let fs = convert_time(ps, TimeUnit::Picosecond, TimeUnit::Femtosecond);
        assert!((fs - 1000.0).abs() < 1e-10);
    }

    /// The unit names program outputs write, case-insensitively.
    #[test]
    fn unit_names() {
        assert_eq!(LengthUnit::from_name("Bohr"), Some(LengthUnit::Bohr));
        assert_eq!(LengthUnit::from_name("angstrom"), Some(LengthUnit::Angstrom));
        assert_eq!(LengthUnit::from_name("alat"), None);
        for name in ["hartree", "a.u.", "AU"] {
            assert_eq!(EnergyUnit::from_name(name), Some(EnergyUnit::Hartree), "{name}");
        }
        assert_eq!(EnergyUnit::from_name("eV"), Some(EnergyUnit::EV));
        assert_eq!(EnergyUnit::from_name("Ry"), None);
        assert_eq!(PressureUnit::from_name("GPa"), Some(PressureUnit::GPa));
        assert_eq!(PressureUnit::from_name("bar"), Some(PressureUnit::Bar));
        assert_eq!(PressureUnit::from_name("kbar"), Some(PressureUnit::Kbar));
        assert_eq!(PressureUnit::from_name("atm"), None);
    }

    /// CODATA 2018: e and N_A are exact, the thermochemical calorie is 4.184 J.
    #[test]
    fn kcal_per_mol_is_codata_2018() {
        let want = 1.602_176_634e-19 * AVOGADRO / 4184.0;
        assert!((EV_TO_KCAL_MOL - want).abs() < 1e-12, "{EV_TO_KCAL_MOL} vs {want}");
    }
}
