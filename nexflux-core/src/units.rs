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
pub const HARTREE_TO_EV: f64 = 27.211_396_132;
pub const EV_TO_HARTREE: f64 = 1.0 / HARTREE_TO_EV;

/// eV → kcal/mol
pub const EV_TO_KCAL_MOL: f64 = 23.060_541_9;
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

/// GPa → kBar
pub const GPA_TO_KBAR: f64 = 10.0;
pub const KBAR_TO_GPA: f64 = 0.1;

/// Avogadro 常数
pub const AVOGADRO: f64 = 6.022_140_76e23;

// ── 长度 ──────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LengthUnit {
    Angstrom,
    Bohr,
    Nanometer,
    Picometer,
}

impl LengthUnit {
    fn to_angstrom(self, v: f64) -> f64 {
        match self {
            Self::Angstrom  => v,
            Self::Bohr      => v * BOHR_TO_ANG,
            Self::Nanometer => v * 10.0,
            Self::Picometer => v * 0.01,
        }
    }
    fn from_angstrom(self, v: f64) -> f64 {
        match self {
            Self::Angstrom  => v,
            Self::Bohr      => v * ANG_TO_BOHR,
            Self::Nanometer => v * 0.1,
            Self::Picometer => v * 100.0,
        }
    }
}

pub fn convert_length(value: f64, from: LengthUnit, to: LengthUnit) -> f64 {
    to.from_angstrom(from.to_angstrom(value))
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
    fn to_ev(self, v: f64) -> f64 {
        match self {
            Self::EV         => v,
            Self::Hartree    => v * HARTREE_TO_EV,
            Self::KcalPerMol => v * KCAL_MOL_TO_EV,
            Self::KJPerMol   => v * KJ_MOL_TO_EV,
            Self::Wavenumber => v * WAVENUMBER_TO_EV,
        }
    }
    fn from_ev(self, v: f64) -> f64 {
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
    to.from_ev(from.to_ev(value))
}

// ── 压强/应力 ─────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PressureUnit {
    /// eV/Å³（内部单位）
    EVPerAng3,
    GPa,
    Kbar,
}

impl PressureUnit {
    fn to_ev_ang3(self, v: f64) -> f64 {
        match self {
            Self::EVPerAng3 => v,
            Self::GPa       => v * GPA_TO_EV_ANG3,
            Self::Kbar      => v * KBAR_TO_GPA * GPA_TO_EV_ANG3,
        }
    }
    fn from_ev_ang3(self, v: f64) -> f64 {
        match self {
            Self::EVPerAng3 => v,
            Self::GPa       => v * EV_ANG3_TO_GPA,
            Self::Kbar      => v * EV_ANG3_TO_GPA * GPA_TO_KBAR,
        }
    }
}

pub fn convert_pressure(value: f64, from: PressureUnit, to: PressureUnit) -> f64 {
    to.from_ev_ang3(from.to_ev_ang3(value))
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
    fn from_fs(self, v: f64) -> f64 {
        match self {
            Self::Femtosecond => v,
            Self::Picosecond  => v * 0.001,
        }
    }
}

pub fn convert_time(value: f64, from: TimeUnit, to: TimeUnit) -> f64 {
    to.from_fs(from.to_fs(value))
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
