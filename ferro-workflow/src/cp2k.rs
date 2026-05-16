use ferro_core::{guess_spin, Frame};
use anyhow::Result;
use std::fmt::Write;

use crate::cp2k_basis_db::{self, DbFunc};

// ── Task ────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Default)]
pub enum Cp2kTask {
    #[default]
    Energy,
    ForceEnergy,
    GeoOpt,
    CellOpt,
    MD,
    Frequency,
}

impl Cp2kTask {
    fn run_type(&self) -> &'static str {
        match self {
            Cp2kTask::Energy      => "ENERGY",
            Cp2kTask::ForceEnergy => "ENERGY_FORCE",
            Cp2kTask::GeoOpt      => "GEO_OPT",
            Cp2kTask::CellOpt     => "CELL_OPT",
            Cp2kTask::MD          => "MD",
            Cp2kTask::Frequency   => "VIBRATIONAL_ANALYSIS",
        }
    }
}

// ── Functional ──────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Default)]
pub enum Cp2kFunctional {
    #[default]
    PBE,
    BLYP,
    PBE0,
    B3LYP,
    RevPBE,
    PBEsol,
    SCAN,
    R2SCAN,
    HSE06,
    Custom(String),
}

// ── Basis set ───────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Default)]
pub enum Cp2kBasis {
    #[default]
    DzvpMoloptSr,   // DZVP-MOLOPT-SR-GTH (default, good all-round)
    TzvpMolopt,     // TZVP-MOLOPT-GTH
    Tzv2pMolopt,    // TZV2P-MOLOPT-GTH
    DzvpGth,        // DZVP-GTH (older GTH style)
    TzvpGth,        // TZVP-GTH
    PobDzvp,        // pob-DZVP (all-electron, periodic)
    PobTzvp,        // pob-TZVP (all-electron, periodic)
    Custom(String),
}

impl Cp2kBasis {
    /// Family prefix used for precise database matching (`None` = not DB-backed).
    fn db_prefix(&self) -> Option<&'static str> {
        match self {
            Cp2kBasis::DzvpMoloptSr => Some("DZVP-MOLOPT-SR"),
            Cp2kBasis::TzvpMolopt   => Some("TZVP-MOLOPT"),
            Cp2kBasis::Tzv2pMolopt  => Some("TZV2P-MOLOPT"),
            Cp2kBasis::PobDzvp      => Some("pob-DZVP"),
            Cp2kBasis::PobTzvp      => Some("pob-TZVP"),
            Cp2kBasis::DzvpGth | Cp2kBasis::TzvpGth | Cp2kBasis::Custom(_) => None,
        }
    }

    fn is_all_electron(&self) -> bool {
        matches!(self, Cp2kBasis::PobDzvp | Cp2kBasis::PobTzvp)
    }
}

// ── Dispersion ──────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Default)]
pub enum Cp2kDispersion {
    #[default]
    None,
    D3,
    D3BJ,
}

// ── SCF solver ──────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Default)]
pub enum Cp2kScf {
    #[default]
    Diag,   // Diagonalisation + Broyden (robust, good for metals)
    OT,     // Orbital Transform (better for insulators/band-gap systems)
}

// ── PBC ─────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Default)]
pub enum Cp2kPbc {
    #[default]
    XYZ,
    Z,
    None,
}

impl Cp2kPbc {
    fn as_str(&self) -> &'static str {
        match self {
            Cp2kPbc::XYZ  => "XYZ",
            Cp2kPbc::Z    => "Z",
            Cp2kPbc::None => "NONE",
        }
    }
}

// ── Output options ──────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Default)]
pub enum Cp2kCubePrint {
    #[default]
    None,
    Density,
    ELF,
    HartreePot,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub enum Cp2kAtomCharge {
    #[default]
    None,
    Mulliken,
    Hirshfeld,
    HirshfeldI,
}

// ── MD settings ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Default)]
pub enum Cp2kThermostat {
    #[default]
    CSVR,       // Canonical sampling — robust default
    NoseHoover,
    Langevin,
    None,       // NVE
}

#[derive(Debug, Clone)]
pub struct Cp2kMdParams {
    pub steps: u32,
    pub timestep: f64,      // fs
    pub temperature: f64,   // K
    pub thermostat: Cp2kThermostat,
    pub traj_freq: u32,     // write every N steps
    pub barostat: bool,     // NPT
}

impl Default for Cp2kMdParams {
    fn default() -> Self {
        Self {
            steps: 10000,
            timestep: 1.0,
            temperature: 298.15,
            thermostat: Cp2kThermostat::CSVR,
            traj_freq: 100,
            barostat: false,
        }
    }
}

// ── Builder ─────────────────────────────────────────────────────────────────

pub struct Cp2kJobBuilder {
    pub frame: Frame,
    pub project: String,
    pub task: Cp2kTask,
    pub functional: Cp2kFunctional,
    pub basis: Cp2kBasis,
    pub dispersion: Cp2kDispersion,
    pub scf: Cp2kScf,
    pub pbc: Cp2kPbc,
    pub kpoints: Option<[u32; 3]>,
    pub cutoff: u32,
    pub rel_cutoff: u32,
    pub smear: bool,
    pub cube_print: Cp2kCubePrint,
    pub atom_charge: Cp2kAtomCharge,
    pub export_molden: bool,
    pub md: Cp2kMdParams,
    /// 自动从结构推断自旋多重度（magmom→氧化态→奇偶）。
    /// 关闭后使用 `frame.multiplicity`（供 CLI 显式 `--multiplicity` 覆盖）。
    pub auto_spin: bool,
}

impl Cp2kJobBuilder {
    pub fn new(frame: Frame) -> Self {
        let pbc = if frame.cell.is_some() {
            match frame.pbc {
                [true, true, true]   => Cp2kPbc::XYZ,
                [false, false, true] => Cp2kPbc::Z,
                _                    => Cp2kPbc::XYZ,
            }
        } else {
            Cp2kPbc::None
        };
        Self {
            project: "ferro".to_string(),
            task: Cp2kTask::Energy,
            functional: Cp2kFunctional::PBE,
            basis: Cp2kBasis::DzvpMoloptSr,
            dispersion: Cp2kDispersion::None,
            scf: Cp2kScf::Diag,
            pbc,
            kpoints: None,
            cutoff: 400,
            rel_cutoff: 50,
            smear: false,
            cube_print: Cp2kCubePrint::None,
            atom_charge: Cp2kAtomCharge::None,
            export_molden: false,
            md: Cp2kMdParams::default(),
            auto_spin: true,
            frame,
        }
    }

    /// 解析有效电荷与自旋多重度。
    ///
    /// 返回 `(charge, multiplicity, needs_uks)`。`auto_spin` 开启时用
    /// [`ferro_core::guess_spin`] 推断；关闭时取 `frame.multiplicity`。
    /// 多重度 > 1 时启用非限制性 Kohn-Sham (UKS)。
    fn resolved_spin(&self) -> (i32, u32, bool) {
        let mult = if self.auto_spin {
            guess_spin(&self.frame).multiplicity
        } else {
            self.frame.multiplicity
        };
        (self.frame.charge, mult, mult > 1)
    }

    pub fn build(&self) -> Result<String> {
        let mut out = String::new();
        writeln!(out, "# Generated by Ferro")?;
        writeln!(out)?;
        self.write_global(&mut out)?;
        self.write_force_eval(&mut out)?;
        self.write_motion(&mut out)?;
        Ok(out)
    }

    // ── sections ──────────────────────────────────────────────────────────

    fn write_global(&self, out: &mut String) -> Result<()> {
        writeln!(out, "&GLOBAL")?;
        writeln!(out, "  PROJECT {}", self.project)?;
        writeln!(out, "  PRINT_LEVEL LOW")?;
        writeln!(out, "  RUN_TYPE {}", self.task.run_type())?;
        writeln!(out, "&END GLOBAL")?;
        writeln!(out)?;
        Ok(())
    }

    fn write_force_eval(&self, out: &mut String) -> Result<()> {
        writeln!(out, "&FORCE_EVAL")?;
        writeln!(out, "  METHOD Quickstep")?;
        writeln!(out)?;
        self.write_subsys(out)?;
        self.write_dft(out)?;
        if self.task == Cp2kTask::ForceEnergy {
            writeln!(out, "  &PRINT")?;
            writeln!(out, "    &FORCES ON")?;
            writeln!(out, "    &END FORCES")?;
            writeln!(out, "  &END PRINT")?;
        }
        if matches!(self.task, Cp2kTask::CellOpt) {
            writeln!(out, "  STRESS_TENSOR ANALYTICAL")?;
        }
        writeln!(out, "&END FORCE_EVAL")?;
        writeln!(out)?;
        Ok(())
    }

    fn write_subsys(&self, out: &mut String) -> Result<()> {
        writeln!(out, "  &SUBSYS")?;

        // &CELL
        writeln!(out, "    &CELL")?;
        if let Some(cell) = &self.frame.cell {
            let m = &cell.matrix;
            writeln!(out, "      A  {:14.8}  {:14.8}  {:14.8}", m[(0,0)], m[(0,1)], m[(0,2)])?;
            writeln!(out, "      B  {:14.8}  {:14.8}  {:14.8}", m[(1,0)], m[(1,1)], m[(1,2)])?;
            writeln!(out, "      C  {:14.8}  {:14.8}  {:14.8}", m[(2,0)], m[(2,1)], m[(2,2)])?;
        }
        writeln!(out, "      PERIODIC {}", self.pbc.as_str())?;
        writeln!(out, "    &END CELL")?;
        writeln!(out)?;

        // &COORD
        writeln!(out, "    &COORD")?;
        for atom in &self.frame.atoms {
            writeln!(out, "      {:4}  {:14.8}  {:14.8}  {:14.8}",
                atom.element, atom.position.x, atom.position.y, atom.position.z)?;
        }
        writeln!(out, "    &END COORD")?;
        writeln!(out)?;

        // &KIND per unique element
        let elements = self.frame.unique_elements();
        for elem in &elements {
            let (bname, q) = self.basis_name_for(elem);
            let pname = self.potential_name_for(elem, q);
            writeln!(out, "    &KIND {elem}")?;
            writeln!(out, "      ELEMENT {elem}")?;
            writeln!(out, "      BASIS_SET {bname}")?;
            writeln!(out, "      POTENTIAL {pname}")?;
            writeln!(out, "    &END KIND")?;
        }

        writeln!(out, "  &END SUBSYS")?;
        writeln!(out)?;
        Ok(())
    }

    fn write_dft(&self, out: &mut String) -> Result<()> {
        writeln!(out, "  &DFT")?;
        for bf in self.basis_files() {
            writeln!(out, "    BASIS_SET_FILE_NAME  {bf}")?;
        }
        for pf in self.potential_files() {
            writeln!(out, "    POTENTIAL_FILE_NAME  {pf}")?;
        }
        writeln!(out)?;
        let (charge, mult, needs_uks) = self.resolved_spin();
        writeln!(out, "    CHARGE {charge}")?;
        writeln!(out, "    MULTIPLICITY {mult}")?;
        if needs_uks {
            writeln!(out, "    UKS")?;
        }
        writeln!(out)?;

        // &QS
        writeln!(out, "    &QS")?;
        writeln!(out, "      EPS_DEFAULT 1E-10")?;
        if self.task == Cp2kTask::MD {
            writeln!(out, "      EXTRAPOLATION ASPC")?;
            writeln!(out, "      EXTRAPOLATION_ORDER 3")?;
        }
        writeln!(out, "    &END QS")?;
        writeln!(out)?;

        // &MGRID
        writeln!(out, "    &MGRID")?;
        writeln!(out, "      CUTOFF {}", self.cutoff)?;
        writeln!(out, "      REL_CUTOFF {}", self.rel_cutoff)?;
        if matches!(self.basis, Cp2kBasis::TzvpMolopt | Cp2kBasis::Tzv2pMolopt) {
            writeln!(out, "      NGRIDS 5")?;
        }
        writeln!(out, "    &END MGRID")?;
        writeln!(out)?;

        // &SCF
        let eps_scf = match self.task {
            Cp2kTask::MD        => "1E-6",
            Cp2kTask::Frequency => "1E-7",
            _                   => "5E-7",
        };
        writeln!(out, "    &SCF")?;
        match self.scf {
            Cp2kScf::Diag => {
                writeln!(out, "      MAX_SCF 128")?;
                writeln!(out, "      EPS_SCF {eps_scf}")?;
                writeln!(out, "#     SCF_GUESS RESTART")?;
                writeln!(out, "      &DIAGONALIZATION")?;
                writeln!(out, "        ALGORITHM STANDARD")?;
                writeln!(out, "      &END DIAGONALIZATION")?;
                writeln!(out, "      &MIXING")?;
                writeln!(out, "        METHOD BROYDEN_MIXING")?;
                writeln!(out, "        ALPHA 0.4")?;
                writeln!(out, "        NBROYDEN 8")?;
                writeln!(out, "      &END MIXING")?;
                if self.smear {
                    writeln!(out, "      &SMEAR")?;
                    writeln!(out, "        METHOD FERMI_DIRAC")?;
                    writeln!(out, "        ELECTRONIC_TEMPERATURE 300")?;
                    writeln!(out, "      &END SMEAR")?;
                }
            }
            Cp2kScf::OT => {
                writeln!(out, "      MAX_SCF 25")?;
                writeln!(out, "      EPS_SCF {eps_scf}")?;
                writeln!(out, "#     SCF_GUESS RESTART")?;
                let prec = if self.frame.atoms.len() < 300 {
                    "FULL_ALL"
                } else {
                    "FULL_KINETIC"
                };
                writeln!(out, "      &OT")?;
                writeln!(out, "        PRECONDITIONER {prec}")?;
                writeln!(out, "        MINIMIZER DIIS")?;
                writeln!(out, "      &END OT")?;
                writeln!(out, "      &OUTER_SCF")?;
                writeln!(out, "        MAX_SCF 20")?;
                writeln!(out, "        EPS_SCF {eps_scf}")?;
                writeln!(out, "      &END OUTER_SCF")?;
            }
        }
        writeln!(out, "      &PRINT")?;
        if matches!(self.task, Cp2kTask::MD | Cp2kTask::Frequency) {
            writeln!(out, "        &RESTART OFF")?;
            writeln!(out, "        &END RESTART")?;
        } else {
            writeln!(out, "        &RESTART")?;
            writeln!(out, "          BACKUP_COPIES 0")?;
            writeln!(out, "        &END RESTART")?;
        }
        writeln!(out, "      &END PRINT")?;
        writeln!(out, "    &END SCF")?;
        writeln!(out)?;

        // &XC
        writeln!(out, "    &XC")?;
        self.write_xc_functional(out)?;
        self.write_dispersion(out)?;
        writeln!(out, "    &END XC")?;
        writeln!(out)?;

        // &POISSON
        writeln!(out, "    &POISSON")?;
        writeln!(out, "      PERIODIC {}", self.pbc.as_str())?;
        writeln!(out, "    &END POISSON")?;

        // &KPOINTS
        if let Some([k1, k2, k3]) = self.kpoints {
            writeln!(out)?;
            writeln!(out, "    &KPOINTS")?;
            writeln!(out, "      SCHEME MONKHORST-PACK {k1} {k2} {k3}")?;
            writeln!(out, "    &END KPOINTS")?;
        }

        // DFT/PRINT (cube, charges, molden)
        let has_print = !matches!(self.cube_print, Cp2kCubePrint::None)
            || !matches!(self.atom_charge, Cp2kAtomCharge::None)
            || self.export_molden;
        if has_print {
            writeln!(out)?;
            writeln!(out, "    &PRINT")?;
            match &self.cube_print {
                Cp2kCubePrint::Density    => {
                    writeln!(out, "      &E_DENSITY_CUBE")?;
                    writeln!(out, "        STRIDE 1")?;
                    writeln!(out, "      &END E_DENSITY_CUBE")?;
                }
                Cp2kCubePrint::ELF        => {
                    writeln!(out, "      &ELF_CUBE")?;
                    writeln!(out, "        STRIDE 1")?;
                    writeln!(out, "      &END ELF_CUBE")?;
                }
                Cp2kCubePrint::HartreePot => {
                    writeln!(out, "      &V_HARTREE_CUBE")?;
                    writeln!(out, "        STRIDE 1")?;
                    writeln!(out, "      &END V_HARTREE_CUBE")?;
                }
                Cp2kCubePrint::None => {}
            }
            match &self.atom_charge {
                Cp2kAtomCharge::Mulliken => {
                    writeln!(out, "      &MULLIKEN")?;
                    writeln!(out, "        PRINT_ALL F")?;
                    writeln!(out, "      &END MULLIKEN")?;
                }
                Cp2kAtomCharge::Hirshfeld => {
                    writeln!(out, "      &HIRSHFELD")?;
                    writeln!(out, "        SHAPE_FUNCTION DENSITY")?;
                    writeln!(out, "      &END HIRSHFELD")?;
                }
                Cp2kAtomCharge::HirshfeldI => {
                    writeln!(out, "      &HIRSHFELD")?;
                    writeln!(out, "        SHAPE_FUNCTION DENSITY")?;
                    writeln!(out, "        SELF_CONSISTENT T")?;
                    writeln!(out, "      &END HIRSHFELD")?;
                }
                Cp2kAtomCharge::None => {}
            }
            if self.export_molden {
                writeln!(out, "      &MO_MOLDEN")?;
                writeln!(out, "        NDIGITS 9")?;
                writeln!(out, "      &END MO_MOLDEN")?;
            }
            writeln!(out, "    &END PRINT")?;
        }

        writeln!(out, "  &END DFT")?;
        writeln!(out)?;
        Ok(())
    }

    fn write_xc_functional(&self, out: &mut String) -> Result<()> {
        match &self.functional {
            Cp2kFunctional::PBE => {
                writeln!(out, "      &XC_FUNCTIONAL PBE")?;
                writeln!(out, "      &END XC_FUNCTIONAL")?;
            }
            Cp2kFunctional::BLYP => {
                writeln!(out, "      &XC_FUNCTIONAL BLYP")?;
                writeln!(out, "      &END XC_FUNCTIONAL")?;
            }
            Cp2kFunctional::RevPBE => {
                writeln!(out, "      &XC_FUNCTIONAL")?;
                writeln!(out, "        &PBE")?;
                writeln!(out, "          PARAMETRIZATION REVPBE")?;
                writeln!(out, "        &END PBE")?;
                writeln!(out, "      &END XC_FUNCTIONAL")?;
            }
            Cp2kFunctional::PBEsol => {
                writeln!(out, "      &XC_FUNCTIONAL")?;
                writeln!(out, "        &PBE")?;
                writeln!(out, "          PARAMETRIZATION PBESOL")?;
                writeln!(out, "        &END PBE")?;
                writeln!(out, "      &END XC_FUNCTIONAL")?;
            }
            Cp2kFunctional::PBE0 => {
                writeln!(out, "      &XC_FUNCTIONAL")?;
                writeln!(out, "        &PBE")?;
                writeln!(out, "          SCALE_X 0.75")?;
                writeln!(out, "          SCALE_C 1.0")?;
                writeln!(out, "        &END PBE")?;
                writeln!(out, "      &END XC_FUNCTIONAL")?;
                writeln!(out, "      &HF")?;
                writeln!(out, "        FRACTION 0.25")?;
                writeln!(out, "        &SCREENING")?;
                writeln!(out, "          EPS_SCHWARZ 1E-7")?;
                writeln!(out, "          SCREEN_ON_INITIAL_P F")?;
                writeln!(out, "        &END SCREENING")?;
                writeln!(out, "      &END HF")?;
            }
            Cp2kFunctional::B3LYP => {
                writeln!(out, "      &XC_FUNCTIONAL")?;
                writeln!(out, "        &LYP")?;
                writeln!(out, "          SCALE_C 0.81")?;
                writeln!(out, "        &END")?;
                writeln!(out, "        &BECKE88")?;
                writeln!(out, "          SCALE_X 0.72")?;
                writeln!(out, "        &END")?;
                writeln!(out, "        &VWN")?;
                writeln!(out, "          FUNCTIONAL_TYPE VWN3")?;
                writeln!(out, "          SCALE_C 0.19")?;
                writeln!(out, "        &END")?;
                writeln!(out, "        &XALPHA")?;
                writeln!(out, "          SCALE_X 0.08")?;
                writeln!(out, "        &END")?;
                writeln!(out, "      &END XC_FUNCTIONAL")?;
                writeln!(out, "      &HF")?;
                writeln!(out, "        FRACTION 0.20")?;
                writeln!(out, "        &SCREENING")?;
                writeln!(out, "          EPS_SCHWARZ 1E-7")?;
                writeln!(out, "          SCREEN_ON_INITIAL_P F")?;
                writeln!(out, "        &END SCREENING")?;
                writeln!(out, "      &END HF")?;
            }
            Cp2kFunctional::HSE06 => {
                writeln!(out, "      &XC_FUNCTIONAL")?;
                writeln!(out, "        &XWPBE")?;
                writeln!(out, "          SCALE_X -0.25")?;
                writeln!(out, "          SCALE_X0 1.0")?;
                writeln!(out, "          OMEGA 0.11")?;
                writeln!(out, "        &END XWPBE")?;
                writeln!(out, "        &PBE")?;
                writeln!(out, "          SCALE_X 0.0")?;
                writeln!(out, "          SCALE_C 1.0")?;
                writeln!(out, "        &END PBE")?;
                writeln!(out, "      &END XC_FUNCTIONAL")?;
                writeln!(out, "      &HF")?;
                writeln!(out, "        FRACTION 0.25")?;
                writeln!(out, "        &SCREENING")?;
                writeln!(out, "          EPS_SCHWARZ 1E-7")?;
                writeln!(out, "          SCREEN_ON_INITIAL_P F")?;
                writeln!(out, "        &END SCREENING")?;
                writeln!(out, "      &END HF")?;
            }
            Cp2kFunctional::SCAN => {
                writeln!(out, "      &XC_FUNCTIONAL")?;
                writeln!(out, "        &MGGA_X_SCAN")?;
                writeln!(out, "        &END MGGA_X_SCAN")?;
                writeln!(out, "        &MGGA_C_SCAN")?;
                writeln!(out, "        &END MGGA_C_SCAN")?;
                writeln!(out, "      &END XC_FUNCTIONAL")?;
            }
            Cp2kFunctional::R2SCAN => {
                writeln!(out, "      &XC_FUNCTIONAL")?;
                writeln!(out, "        &MGGA_X_R2SCAN")?;
                writeln!(out, "        &END MGGA_X_R2SCAN")?;
                writeln!(out, "        &MGGA_C_R2SCAN")?;
                writeln!(out, "        &END MGGA_C_R2SCAN")?;
                writeln!(out, "      &END XC_FUNCTIONAL")?;
            }
            Cp2kFunctional::Custom(name) => {
                writeln!(out, "      &XC_FUNCTIONAL {name}")?;
                writeln!(out, "      &END XC_FUNCTIONAL")?;
            }
        }
        Ok(())
    }

    fn write_dispersion(&self, out: &mut String) -> Result<()> {
        if self.dispersion == Cp2kDispersion::None {
            return Ok(());
        }
        let func_name = match &self.functional {
            Cp2kFunctional::PBE        => "PBE",
            Cp2kFunctional::BLYP       => "BLYP",
            Cp2kFunctional::PBE0       => "PBE0",
            Cp2kFunctional::B3LYP      => "B3LYP",
            Cp2kFunctional::RevPBE     => "REVPBE",
            Cp2kFunctional::PBEsol     => "PBESOL",
            Cp2kFunctional::SCAN       => "SCAN",
            Cp2kFunctional::R2SCAN     => "r2SCAN",
            Cp2kFunctional::HSE06      => "HSE06",
            Cp2kFunctional::Custom(s)  => s.as_str(),
        };
        writeln!(out, "      &VDW_POTENTIAL")?;
        writeln!(out, "        POTENTIAL_TYPE PAIR_POTENTIAL")?;
        writeln!(out, "        &PAIR_POTENTIAL")?;
        writeln!(out, "          PARAMETER_FILE_NAME dftd3.dat")?;
        match self.dispersion {
            Cp2kDispersion::D3   => writeln!(out, "          TYPE DFTD3")?,
            Cp2kDispersion::D3BJ => writeln!(out, "          TYPE DFTD3(BJ)")?,
            Cp2kDispersion::None => {}
        }
        writeln!(out, "          REFERENCE_FUNCTIONAL {func_name}")?;
        writeln!(out, "        &END PAIR_POTENTIAL")?;
        writeln!(out, "      &END VDW_POTENTIAL")?;
        Ok(())
    }

    fn write_motion(&self, out: &mut String) -> Result<()> {
        match &self.task {
            Cp2kTask::GeoOpt => {
                writeln!(out, "&MOTION")?;
                writeln!(out, "  &GEO_OPT")?;
                writeln!(out, "    TYPE MINIMIZATION")?;
                writeln!(out, "    OPTIMIZER BFGS")?;
                writeln!(out, "    MAX_ITER 500")?;
                writeln!(out, "    MAX_DR    3E-3")?;
                writeln!(out, "    RMS_DR    1.5E-3")?;
                writeln!(out, "    MAX_FORCE 4.5E-4")?;
                writeln!(out, "    RMS_FORCE 3E-4")?;
                writeln!(out, "    &BFGS")?;
                writeln!(out, "      TRUST_RADIUS 0.2")?;
                writeln!(out, "    &END BFGS")?;
                writeln!(out, "  &END GEO_OPT")?;
                writeln!(out, "  &PRINT")?;
                writeln!(out, "    &TRAJECTORY")?;
                writeln!(out, "      FORMAT XYZ")?;
                writeln!(out, "      &EACH")?;
                writeln!(out, "        GEO_OPT 1")?;
                writeln!(out, "      &END EACH")?;
                writeln!(out, "    &END TRAJECTORY")?;
                writeln!(out, "  &END PRINT")?;
                writeln!(out, "&END MOTION")?;
            }
            Cp2kTask::CellOpt => {
                writeln!(out, "&MOTION")?;
                writeln!(out, "  &CELL_OPT")?;
                writeln!(out, "    TYPE DIRECT_CELL_OPT")?;
                writeln!(out, "    OPTIMIZER BFGS")?;
                writeln!(out, "    MAX_ITER 400")?;
                writeln!(out, "    MAX_DR    3E-3")?;
                writeln!(out, "    RMS_DR    1.5E-3")?;
                writeln!(out, "    MAX_FORCE 4.5E-4")?;
                writeln!(out, "    RMS_FORCE 3E-4")?;
                writeln!(out, "    PRESSURE_TOLERANCE 100")?;
                writeln!(out, "    KEEP_ANGLES F")?;
                writeln!(out, "    KEEP_SYMMETRY F")?;
                writeln!(out, "    &BFGS")?;
                writeln!(out, "      TRUST_RADIUS 0.2")?;
                writeln!(out, "    &END BFGS")?;
                writeln!(out, "  &END CELL_OPT")?;
                writeln!(out, "  &PRINT")?;
                writeln!(out, "    &TRAJECTORY")?;
                writeln!(out, "      FORMAT XYZ")?;
                writeln!(out, "      &EACH")?;
                writeln!(out, "        CELL_OPT 1")?;
                writeln!(out, "      &END EACH")?;
                writeln!(out, "    &END TRAJECTORY")?;
                writeln!(out, "  &END PRINT")?;
                writeln!(out, "&END MOTION")?;
            }
            Cp2kTask::MD => {
                let md = &self.md;
                let ensemble = if md.barostat { "NPT_F" } else { "NVT" };
                writeln!(out, "&MOTION")?;
                writeln!(out, "  &MD")?;
                writeln!(out, "    ENSEMBLE {ensemble}")?;
                writeln!(out, "    STEPS {}", md.steps)?;
                writeln!(out, "    TIMESTEP {:.4} # fs", md.timestep)?;
                writeln!(out, "    TEMPERATURE {:.2} # K", md.temperature)?;
                match md.thermostat {
                    Cp2kThermostat::CSVR => {
                        writeln!(out, "    &THERMOSTAT")?;
                        writeln!(out, "      TYPE CSVR")?;
                        writeln!(out, "      &CSVR")?;
                        writeln!(out, "        TIMECON_CSVR 100 # fs")?;
                        writeln!(out, "      &END CSVR")?;
                        writeln!(out, "    &END THERMOSTAT")?;
                    }
                    Cp2kThermostat::NoseHoover => {
                        writeln!(out, "    &THERMOSTAT")?;
                        writeln!(out, "      TYPE NOSE")?;
                        writeln!(out, "      &NOSE")?;
                        writeln!(out, "        TIMECON 1000 # fs")?;
                        writeln!(out, "      &END NOSE")?;
                        writeln!(out, "    &END THERMOSTAT")?;
                    }
                    Cp2kThermostat::Langevin => {
                        writeln!(out, "    &THERMOSTAT")?;
                        writeln!(out, "      TYPE LANGEVIN")?;
                        writeln!(out, "      &LANGEVIN")?;
                        writeln!(out, "        GAMMA 0.01 # 1/fs")?;
                        writeln!(out, "      &END LANGEVIN")?;
                        writeln!(out, "    &END THERMOSTAT")?;
                    }
                    Cp2kThermostat::None => {}
                }
                if md.barostat {
                    writeln!(out, "    &BAROSTAT")?;
                    writeln!(out, "      PRESSURE 1.01325E+05 # bar")?;
                    writeln!(out, "      TIMECON 1000 # fs")?;
                    writeln!(out, "    &END BAROSTAT")?;
                }
                writeln!(out, "  &END MD")?;
                writeln!(out, "  &PRINT")?;
                writeln!(out, "    &TRAJECTORY")?;
                writeln!(out, "      FORMAT XYZ")?;
                writeln!(out, "      &EACH")?;
                writeln!(out, "        MD {}", md.traj_freq)?;
                writeln!(out, "      &END EACH")?;
                writeln!(out, "    &END TRAJECTORY")?;
                writeln!(out, "    &VELOCITIES")?;
                writeln!(out, "      &EACH")?;
                writeln!(out, "        MD {}", md.traj_freq)?;
                writeln!(out, "      &END EACH")?;
                writeln!(out, "    &END VELOCITIES")?;
                writeln!(out, "    &FORCES")?;
                writeln!(out, "      &EACH")?;
                writeln!(out, "        MD {}", md.traj_freq)?;
                writeln!(out, "      &END EACH")?;
                writeln!(out, "    &END FORCES")?;
                writeln!(out, "    &RESTART_HISTORY OFF")?;
                writeln!(out, "    &END RESTART_HISTORY")?;
                writeln!(out, "    &RESTART")?;
                writeln!(out, "      BACKUP_COPIES 1")?;
                writeln!(out, "    &END RESTART")?;
                writeln!(out, "  &END PRINT")?;
                writeln!(out, "&END MOTION")?;
            }
            Cp2kTask::Frequency => {
                writeln!(out, "&VIBRATIONAL_ANALYSIS")?;
                writeln!(out, "  NPROC_REP 1")?;
                writeln!(out, "  FULLY_PERIODIC F")?;
                writeln!(out, "&END VIBRATIONAL_ANALYSIS")?;
            }
            _ => {} // Energy, ForceEnergy: no MOTION block
        }
        Ok(())
    }

    // ── helpers ───────────────────────────────────────────────────────────

    /// Database functional class for the current basis / functional.
    fn db_func(&self) -> DbFunc {
        if self.basis.is_all_electron() {
            return DbFunc::AllElectron;
        }
        match self.functional {
            Cp2kFunctional::SCAN | Cp2kFunctional::R2SCAN => DbFunc::Scan,
            _ => DbFunc::Pbe, // PBE / revPBE / PBEsol / BLYP / hybrids → GGA-class
        }
    }

    /// Precise per-element basis-set name. Falls back to the heuristic
    /// `gth_nval` q-guess only when the database has no entry.
    fn basis_name_for(&self, elem: &str) -> (String, Option<u8>) {
        if let Cp2kBasis::Custom(s) = &self.basis {
            return (s.clone(), None);
        }
        if let Some(prefix) = self.basis.db_prefix() {
            if let Some((name, q)) = cp2k_basis_db::basis(elem, prefix, self.db_func()) {
                return (name.to_string(), q);
            }
        }
        // fallback: legacy heuristic
        match &self.basis {
            Cp2kBasis::DzvpMoloptSr => (format!("DZVP-MOLOPT-SR-GTH-q{}", gth_nval(elem)), Some(gth_nval(elem))),
            Cp2kBasis::TzvpMolopt   => (format!("TZVP-MOLOPT-GTH-q{}",  gth_nval(elem)), Some(gth_nval(elem))),
            Cp2kBasis::Tzv2pMolopt  => (format!("TZV2P-MOLOPT-GTH-q{}", gth_nval(elem)), Some(gth_nval(elem))),
            Cp2kBasis::PobDzvp      => ("pob-DZVP".to_string(), None),
            Cp2kBasis::PobTzvp      => ("pob-TZVP".to_string(), None),
            Cp2kBasis::DzvpGth      => ("DZVP-GTH".to_string(), None),
            Cp2kBasis::TzvpGth      => ("TZVP-GTH".to_string(), None),
            Cp2kBasis::Custom(s)    => (s.clone(), None),
        }
    }

    /// Precise per-element pseudopotential name matching the basis valence `q`.
    fn potential_name_for(&self, elem: &str, q: Option<u8>) -> String {
        let func = self.db_func();
        if let Some(name) = cp2k_basis_db::potential(elem, func, q) {
            return name.to_string();
        }
        // fallback: legacy generic names
        match func {
            DbFunc::AllElectron => "ALLELECTRON".to_string(),
            _ => match &self.functional {
                Cp2kFunctional::BLYP | Cp2kFunctional::B3LYP => "GTH-BLYP".to_string(),
                _ => format!("GTH-PBE-q{}", gth_nval(elem)),
            },
        }
    }

    /// BASIS_SET_FILE_NAME entries to emit (CP2K searches all listed files).
    fn basis_files(&self) -> &'static [&'static str] {
        match &self.basis {
            Cp2kBasis::PobDzvp | Cp2kBasis::PobTzvp =>
                &["BASIS_pob", "BASIS_MOLOPT_UZH"],
            Cp2kBasis::DzvpGth | Cp2kBasis::TzvpGth =>
                &["GTH_BASIS_SETS", "BASIS_MOLOPT"],
            _ =>
                &["BASIS_MOLOPT", "BASIS_MOLOPT_UCL", "BASIS_MOLOPT_UZH"],
        }
    }

    /// POTENTIAL_FILE_NAME entries to emit.
    fn potential_files(&self) -> &'static [&'static str] {
        &["POTENTIAL", "POTENTIAL_UZH"]
    }
}

// ── GTH valence electron count (fallback) ────────────────────────────────────

/// Fallback valence-electron count `q`, derived from the basis database rather
/// than a parallel hand-maintained table.
///
/// Returns the `q` of the element's default `DZVP-MOLOPT-SR` basis (the same
/// resolver the builder uses for that family — this encodes the conventional
/// valence choice, e.g. Zn → 12 rather than the deep-semicore 20).  Falls back
/// to any Pbe basis entry's `q`, or a conservative 4 for elements absent from
/// the database.  Only reached when the precise lookup misses.
fn gth_nval(element: &str) -> u8 {
    use cp2k_basis_db::{DbFunc, DbKind, DB};

    if let Some((_, Some(q))) =
        cp2k_basis_db::basis(element, "DZVP-MOLOPT-SR", DbFunc::Pbe)
    {
        return q;
    }
    DB.iter()
        .filter(|e| {
            e.element == element && e.kind == DbKind::Basis && e.func == DbFunc::Pbe
        })
        .filter_map(|e| e.q)
        .min()
        .unwrap_or(4)
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Frame};
    use nalgebra::Vector3;

    fn water_frame() -> Frame {
        Frame {
            atoms: vec![
                Atom::new("O", Vector3::new(0.0,  0.0,  0.0)),
                Atom::new("H", Vector3::new(0.96, 0.0,  0.0)),
                Atom::new("H", Vector3::new(-0.24, 0.93, 0.0)),
            ],
            ..Default::default()
        }
    }

    #[test]
    fn test_energy_molecular() {
        let mut b = Cp2kJobBuilder::new(water_frame());
        b.pbc = Cp2kPbc::None;
        let inp = b.build().unwrap();
        assert!(inp.contains("RUN_TYPE ENERGY"));
        assert!(inp.contains("PERIODIC NONE"));
        assert!(inp.contains("&KIND O"));
        assert!(inp.contains("&KIND H"));
        assert!(inp.contains("DZVP-MOLOPT-SR-GTH-q6")); // O has 6 val elec
        assert!(inp.contains("DZVP-MOLOPT-SR-GTH-q1")); // H has 1 val elec
        assert!(inp.contains("&XC_FUNCTIONAL PBE"));
    }

    #[test]
    fn test_geo_opt_with_dispersion() {
        let mut b = Cp2kJobBuilder::new(water_frame());
        b.pbc = Cp2kPbc::None;
        b.task = Cp2kTask::GeoOpt;
        b.dispersion = Cp2kDispersion::D3BJ;
        let inp = b.build().unwrap();
        assert!(inp.contains("RUN_TYPE GEO_OPT"));
        assert!(inp.contains("&GEO_OPT"));
        assert!(inp.contains("TYPE DFTD3(BJ)"));
        assert!(inp.contains("REFERENCE_FUNCTIONAL PBE"));
    }

    #[test]
    fn test_md_nvt() {
        let mut b = Cp2kJobBuilder::new(water_frame());
        b.pbc = Cp2kPbc::None;
        b.task = Cp2kTask::MD;
        b.md.steps = 5000;
        b.md.temperature = 1000.0;
        let inp = b.build().unwrap();
        assert!(inp.contains("RUN_TYPE MD"));
        assert!(inp.contains("ENSEMBLE NVT"));
        assert!(inp.contains("STEPS 5000"));
        assert!(inp.contains("TEMPERATURE 1000.00"));
        assert!(inp.contains("TYPE CSVR"));
    }

    #[test]
    fn test_pbe0_with_ot() {
        let mut b = Cp2kJobBuilder::new(water_frame());
        b.pbc = Cp2kPbc::None;
        b.functional = Cp2kFunctional::PBE0;
        b.scf = Cp2kScf::OT;
        let inp = b.build().unwrap();
        assert!(inp.contains("FRACTION 0.25"));
        assert!(inp.contains("&OT"));
        assert!(inp.contains("FULL_ALL")); // small system
    }

    #[test]
    fn test_gth_nval() {
        assert_eq!(gth_nval("O"),  6);
        assert_eq!(gth_nval("Si"), 4);
        assert_eq!(gth_nval("Fe"), 16);
        assert_eq!(gth_nval("Zn"), 12);
        assert_eq!(gth_nval("P"),  5);
    }

    fn cell_frame(elems: &[(&str, usize)]) -> Frame {
        use ferro_core::Cell;
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut f = Frame::with_cell(cell, [true; 3]);
        let mut x = 0.0;
        for &(s, n) in elems {
            for _ in 0..n {
                f.add_atom(Atom::new(s, Vector3::new(x, 0.0, 0.0)));
                x += 1.6;
            }
        }
        f
    }

    #[test]
    fn test_auto_spin_fe2o3_high_spin() {
        // Fe³⁺ d⁵ 高自旋 → 5×2 = 10 未成对 → 多重度 11，启用 UKS
        let b = Cp2kJobBuilder::new(cell_frame(&[("Fe", 2), ("O", 3)]));
        assert!(b.auto_spin, "默认应开启 auto_spin");
        let inp = b.build().unwrap();
        assert!(inp.contains("MULTIPLICITY 11"), "Fe2O3 应为多重度 11");
        assert!(inp.contains("\n    UKS\n"), "开壳层应启用 UKS");
    }

    #[test]
    fn test_auto_spin_znp2o6_singlet() {
        // 全闭壳 → 多重度 1，无 UKS
        let b = Cp2kJobBuilder::new(cell_frame(&[("Zn", 1), ("P", 2), ("O", 6)]));
        let inp = b.build().unwrap();
        assert!(inp.contains("MULTIPLICITY 1"));
        assert!(!inp.contains("\n    UKS\n"), "闭壳层不应有 UKS");
    }

    #[test]
    fn test_auto_spin_off_respects_frame() {
        // 关闭 auto_spin → 使用 frame.multiplicity（CLI --multiplicity 覆盖路径）
        let mut f = cell_frame(&[("Fe", 2), ("O", 3)]);
        f.multiplicity = 1;
        let mut b = Cp2kJobBuilder::new(f);
        b.auto_spin = false;
        let inp = b.build().unwrap();
        assert!(inp.contains("MULTIPLICITY 1"), "应尊重手动 multiplicity");
    }

    #[test]
    fn test_db_precise_basis_and_potential() {
        // Fe 经数据库应得精确 q16 基组与匹配赝势，而非启发式猜测名
        let mut b = Cp2kJobBuilder::new(cell_frame(&[("Fe", 1), ("O", 1)]));
        b.auto_spin = false;
        let inp = b.build().unwrap();
        assert!(inp.contains("BASIS_SET DZVP-MOLOPT-SR-GTH-q16"), "Fe 精确基组");
        assert!(inp.contains("BASIS_SET DZVP-MOLOPT-SR-GTH-q6"),  "O 精确基组");
        assert!(inp.contains("POTENTIAL GTH-PBE-q16"), "Fe 赝势 q 应与基组一致");
        assert!(inp.contains("BASIS_SET_FILE_NAME  BASIS_MOLOPT"));
        assert!(inp.contains("POTENTIAL_FILE_NAME  POTENTIAL_UZH"));
    }

    #[test]
    fn test_scan_functional_uses_scan_basis() {
        let mut b = Cp2kJobBuilder::new(cell_frame(&[("O", 1)]));
        b.auto_spin = false;
        b.functional = Cp2kFunctional::SCAN;
        b.basis = Cp2kBasis::TzvpMolopt;
        let inp = b.build().unwrap();
        assert!(inp.contains("BASIS_SET TZVP-MOLOPT-SCAN-GTH-q6"), "SCAN 基组");
    }

    #[test]
    fn test_pob_all_electron() {
        let mut b = Cp2kJobBuilder::new(cell_frame(&[("Si", 1), ("O", 2)]));
        b.auto_spin = false;
        b.basis = Cp2kBasis::PobTzvp;
        let inp = b.build().unwrap();
        assert!(inp.contains("BASIS_SET pob-TZVP"), "全电子 pob 基组");
        assert!(inp.contains("POTENTIAL ALL"), "全电子应使用 ALL/ALLELECTRON");
        assert!(inp.contains("BASIS_SET_FILE_NAME  BASIS_pob"));
    }
}
