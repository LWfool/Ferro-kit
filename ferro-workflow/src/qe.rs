//! Quantum ESPRESSO `pw.x` 输入文件构建器。
//!
//! 与 [`crate::cp2k::Cp2kJobBuilder`] 平行的工作流层 builder：支持 scf/nscf/
//! bands/relax/vc-relax/md/vc-md 任务，泛函、平面波截断、k 点、展宽、
//! 自旋（复用 [`ferro_core::guess_spin`]，与 CP2K 路径一致）等。
//!
//! 这是 **计算任务** 用途（`fe-job -s qe`）。纯结构格式转换（`fe-convert`，
//! 可与读取器往返）请用 `ferro_io::write_qe_input`。两者分属不能互相依赖的
//! 中间层，QE 卡片格式各写一份，修改其一须同步另一处。

use ferro_core::{guess_spin, Frame};
use anyhow::Result;
use std::fmt::Write;

// ── Task ────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Default)]
pub enum QeTask {
    #[default]
    Scf,
    Nscf,
    Bands,
    Relax,    // 原子弛豫
    VcRelax,  // 变胞弛豫
    Md,       // Born-Oppenheimer MD
    VcMd,     // 变胞 MD
}

impl QeTask {
    fn calculation(&self) -> &'static str {
        match self {
            QeTask::Scf     => "scf",
            QeTask::Nscf    => "nscf",
            QeTask::Bands   => "bands",
            QeTask::Relax   => "relax",
            QeTask::VcRelax => "vc-relax",
            QeTask::Md      => "md",
            QeTask::VcMd    => "vc-md",
        }
    }
    fn needs_ions(&self) -> bool {
        matches!(self, QeTask::Relax | QeTask::VcRelax | QeTask::Md | QeTask::VcMd)
    }
    fn needs_cell(&self) -> bool {
        matches!(self, QeTask::VcRelax | QeTask::VcMd)
    }
    fn is_md(&self) -> bool {
        matches!(self, QeTask::Md | QeTask::VcMd)
    }
}

// ── Functional ──────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Default)]
pub enum QeFunctional {
    #[default]
    PBE,
    PBEsol,
    RevPBE,
    BLYP,
    SCAN,
    R2SCAN,
    PBE0,
    HSE06,
    Custom(String),
}

impl QeFunctional {
    /// `input_dft` 取值；`None` 表示由赝势自带泛函决定（如普通 PBE）。
    fn input_dft(&self) -> Option<&str> {
        match self {
            QeFunctional::PBE     => None,
            QeFunctional::PBEsol  => Some("PBESOL"),
            QeFunctional::RevPBE  => Some("REVPBE"),
            QeFunctional::BLYP    => Some("BLYP"),
            QeFunctional::SCAN    => Some("SCAN"),
            QeFunctional::R2SCAN  => Some("R2SCAN"),
            QeFunctional::PBE0    => Some("PBE0"),
            QeFunctional::HSE06   => Some("HSE"),
            QeFunctional::Custom(s) => Some(s.as_str()),
        }
    }
}

// ── Smearing ────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Default)]
pub enum QeSmearing {
    #[default]
    None,                 // 绝缘体：固定占据
    Gaussian,
    MethfesselPaxton,     // 金属推荐
    MarzariVanderbilt,    // cold smearing
    FermiDirac,
}

impl QeSmearing {
    fn as_str(&self) -> &'static str {
        match self {
            QeSmearing::None              => "fixed",
            QeSmearing::Gaussian          => "gaussian",
            QeSmearing::MethfesselPaxton  => "mp",
            QeSmearing::MarzariVanderbilt => "mv",
            QeSmearing::FermiDirac        => "fd",
        }
    }
}

// ── MD 参数 ─────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct QeMdParams {
    pub steps: u32,
    pub dt_au: f64,        // QE 时间步，单位 Rydberg a.u.（≈ 0.048 fs）
    pub temperature: f64,  // K
}

impl Default for QeMdParams {
    fn default() -> Self {
        Self { steps: 1000, dt_au: 20.0, temperature: 300.0 }
    }
}

// ── Builder ─────────────────────────────────────────────────────────────────

pub struct QeJobBuilder {
    pub frame: Frame,
    pub prefix: String,
    pub task: QeTask,
    pub functional: QeFunctional,
    pub ecutwfc: f64,           // Ry
    pub ecutrho: Option<f64>,   // Ry；None → QE 默认 4×ecutwfc
    pub kpoints: Option<[u32; 3]>,   // None → Gamma 点
    pub kpoints_shift: [u32; 3],
    pub smearing: QeSmearing,
    pub degauss: f64,           // Ry
    pub conv_thr: f64,
    pub mixing_beta: f64,
    pub pseudo_dir: String,
    pub pseudo_suffix: String,  // 赝势文件后缀，默认 ".UPF"
    /// 自动由结构推断自旋（nspin/tot_magnetization）；关闭则用 frame.multiplicity。
    pub auto_spin: bool,
    pub md: QeMdParams,
}

impl QeJobBuilder {
    pub fn new(frame: Frame) -> Self {
        Self {
            prefix: "ferro".to_string(),
            task: QeTask::Scf,
            functional: QeFunctional::PBE,
            ecutwfc: 50.0,
            ecutrho: None,
            kpoints: None,
            kpoints_shift: [0, 0, 0],
            smearing: QeSmearing::None,
            degauss: 0.01,
            conv_thr: 1e-8,
            mixing_beta: 0.4,
            pseudo_dir: "./pseudo".to_string(),
            pseudo_suffix: ".UPF".to_string(),
            auto_spin: true,
            md: QeMdParams::default(),
            frame,
        }
    }

    /// 解析有效自旋：返回 `(nspin, tot_magnetization)`。
    fn resolved_spin(&self) -> (u32, u32) {
        let mult = if self.auto_spin {
            guess_spin(&self.frame).multiplicity
        } else {
            self.frame.multiplicity
        };
        if mult > 1 { (2, mult - 1) } else { (1, 0) }
    }

    pub fn build(&self) -> Result<String> {
        let mut out = String::new();
        let elems = self.frame.unique_elements();
        let nat = self.frame.atoms.len();
        let ntyp = elems.len();
        let (nspin, tot_mag) = self.resolved_spin();

        // &CONTROL
        writeln!(out, "&CONTROL")?;
        writeln!(out, "  calculation = '{}',", self.task.calculation())?;
        writeln!(out, "  prefix = '{}',", self.prefix)?;
        writeln!(out, "  pseudo_dir = '{}',", self.pseudo_dir)?;
        writeln!(out, "  outdir = './out',")?;
        if self.task.needs_ions() {
            writeln!(out, "  tprnfor = .true.,")?;
            writeln!(out, "  tstress = .true.,")?;
        }
        if self.task.is_md() {
            writeln!(out, "  nstep = {},", self.md.steps)?;
            writeln!(out, "  dt = {:.4},", self.md.dt_au)?;
        } else if self.task.needs_ions() {
            writeln!(out, "  nstep = 200,")?;
        }
        writeln!(out, "/")?;
        writeln!(out)?;

        // &SYSTEM
        writeln!(out, "&SYSTEM")?;
        writeln!(out, "  ibrav = 0,")?;
        writeln!(out, "  nat = {nat},")?;
        writeln!(out, "  ntyp = {ntyp},")?;
        writeln!(out, "  ecutwfc = {:.1},", self.ecutwfc)?;
        let ecutrho = self.ecutrho.unwrap_or(self.ecutwfc * 4.0);
        writeln!(out, "  ecutrho = {ecutrho:.1},")?;
        if self.frame.charge != 0 {
            writeln!(out, "  tot_charge = {},", self.frame.charge)?;
        }
        if nspin == 2 {
            writeln!(out, "  nspin = 2,")?;
            writeln!(out, "  tot_magnetization = {tot_mag},")?;
        }
        if self.smearing != QeSmearing::None {
            writeln!(out, "  occupations = 'smearing',")?;
            writeln!(out, "  smearing = '{}',", self.smearing.as_str())?;
            writeln!(out, "  degauss = {:.4},", self.degauss)?;
        }
        if let Some(dft) = self.functional.input_dft() {
            writeln!(out, "  input_dft = '{dft}',")?;
        }
        writeln!(out, "/")?;
        writeln!(out)?;

        // &ELECTRONS
        writeln!(out, "&ELECTRONS")?;
        writeln!(out, "  conv_thr = {:.1e},", self.conv_thr)?;
        writeln!(out, "  mixing_beta = {:.2},", self.mixing_beta)?;
        writeln!(out, "  diagonalization = 'david',")?;
        writeln!(out, "/")?;
        writeln!(out)?;

        // &IONS
        if self.task.needs_ions() {
            writeln!(out, "&IONS")?;
            if self.task.is_md() {
                writeln!(out, "  ion_dynamics = 'verlet',")?;
                writeln!(out, "  ion_temperature = 'rescaling',")?;
                writeln!(out, "  tempw = {:.1},", self.md.temperature)?;
            } else {
                writeln!(out, "  ion_dynamics = 'bfgs',")?;
            }
            writeln!(out, "/")?;
            writeln!(out)?;
        }

        // &CELL
        if self.task.needs_cell() {
            writeln!(out, "&CELL")?;
            writeln!(out, "  cell_dynamics = '{}',",
                if self.task.is_md() { "pr" } else { "bfgs" })?;
            writeln!(out, "  press = 0.0,")?;
            writeln!(out, "  press_conv_thr = 0.5,")?;
            writeln!(out, "/")?;
            writeln!(out)?;
        }

        // ATOMIC_SPECIES
        writeln!(out, "ATOMIC_SPECIES")?;
        for el in &elems {
            let mass = ferro_core::data::elements::by_symbol(el)
                .map(|e| e.atomic_mass)
                .unwrap_or(1.0);
            writeln!(out, "  {el}  {mass:.4}  {el}{}", self.pseudo_suffix)?;
        }
        writeln!(out)?;

        // CELL_PARAMETERS
        if let Some(cell) = &self.frame.cell {
            writeln!(out, "CELL_PARAMETERS angstrom")?;
            for i in 0..3 {
                let r = cell.matrix.row(i);
                writeln!(out, "  {:14.10}  {:14.10}  {:14.10}", r[0], r[1], r[2])?;
            }
            writeln!(out)?;
        }

        // ATOMIC_POSITIONS
        writeln!(out, "ATOMIC_POSITIONS angstrom")?;
        for atom in &self.frame.atoms {
            writeln!(out, "  {:4}  {:14.10}  {:14.10}  {:14.10}",
                atom.element, atom.position.x, atom.position.y, atom.position.z)?;
        }
        writeln!(out)?;

        // K_POINTS
        match self.kpoints {
            None => {
                writeln!(out, "K_POINTS gamma")?;
            }
            Some([k1, k2, k3]) => {
                writeln!(out, "K_POINTS automatic")?;
                let [s1, s2, s3] = self.kpoints_shift;
                writeln!(out, "  {k1} {k2} {k3} {s1} {s2} {s3}")?;
            }
        }

        Ok(out)
    }
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame};
    use nalgebra::Vector3;

    fn fe_bcc() -> Frame {
        let cell = Cell::from_lengths_angles(2.87, 2.87, 2.87, 90.0, 90.0, 90.0).unwrap();
        let mut f = Frame::with_cell(cell, [true; 3]);
        f.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        f.add_atom(Atom::new("Fe", Vector3::new(1.435, 1.435, 1.435)));
        f
    }

    fn water() -> Frame {
        let mut f = Frame::new();
        f.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));
        f.add_atom(Atom::new("H", Vector3::new(0.96, 0.0, 0.0)));
        f.add_atom(Atom::new("H", Vector3::new(-0.24, 0.93, 0.0)));
        f
    }

    #[test]
    fn test_scf_basic() {
        let b = QeJobBuilder::new(water());
        let inp = b.build().unwrap();
        assert!(inp.contains("calculation = 'scf'"));
        assert!(inp.contains("ibrav = 0"));
        assert!(inp.contains("nat = 3"));
        assert!(inp.contains("ntyp = 2"));
        assert!(inp.contains("ATOMIC_SPECIES"));
        assert!(inp.contains("O.UPF"));
        assert!(inp.contains("K_POINTS gamma"));
    }

    #[test]
    fn test_vc_relax_has_ions_and_cell() {
        let mut b = QeJobBuilder::new(fe_bcc());
        b.task = QeTask::VcRelax;
        b.auto_spin = false;
        let inp = b.build().unwrap();
        assert!(inp.contains("calculation = 'vc-relax'"));
        assert!(inp.contains("&IONS"));
        assert!(inp.contains("&CELL"));
        assert!(inp.contains("tprnfor = .true."));
        assert!(inp.contains("CELL_PARAMETERS angstrom"));
    }

    #[test]
    fn test_kpoints_automatic() {
        let mut b = QeJobBuilder::new(fe_bcc());
        b.kpoints = Some([4, 4, 4]);
        b.kpoints_shift = [1, 1, 1];
        b.auto_spin = false;
        let inp = b.build().unwrap();
        assert!(inp.contains("K_POINTS automatic"));
        assert!(inp.contains("4 4 4 1 1 1"));
    }

    #[test]
    fn test_auto_spin_fe_metal() {
        // bcc Fe（单质）→ 氧化态法不适用，回退奇偶；2 个 Fe → 52 电子偶 → 单重
        // 显式测自旋通道：手动多重度
        let mut f = fe_bcc();
        f.multiplicity = 5; // 强制开壳层
        let mut b = QeJobBuilder::new(f);
        b.auto_spin = false;
        let inp = b.build().unwrap();
        assert!(inp.contains("nspin = 2"));
        assert!(inp.contains("tot_magnetization = 4"));
    }

    #[test]
    fn test_scan_and_smearing() {
        let mut b = QeJobBuilder::new(fe_bcc());
        b.functional = QeFunctional::SCAN;
        b.smearing = QeSmearing::MethfesselPaxton;
        b.auto_spin = false;
        let inp = b.build().unwrap();
        assert!(inp.contains("input_dft = 'SCAN'"));
        assert!(inp.contains("occupations = 'smearing'"));
        assert!(inp.contains("smearing = 'mp'"));
    }

    #[test]
    fn test_md_block() {
        let mut b = QeJobBuilder::new(fe_bcc());
        b.task = QeTask::Md;
        b.md.steps = 500;
        b.md.temperature = 800.0;
        b.auto_spin = false;
        let inp = b.build().unwrap();
        assert!(inp.contains("calculation = 'md'"));
        assert!(inp.contains("nstep = 500"));
        assert!(inp.contains("ion_dynamics = 'verlet'"));
        assert!(inp.contains("tempw = 800.0"));
    }
}
