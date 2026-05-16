use anyhow::{anyhow, bail, Result};
use clap::Parser;
use ferro::help::{print_fe_job_overview, print_job_help};
use ferro::io_dispatch::read_trajectory;
use ferro_core::{guess_spin, SpinMethod};
use ferro_io::LammpsUnits;
use ferro_workflow::{
    GaussianJobBuilder,
    Cp2kJobBuilder, Cp2kTask, Cp2kFunctional, Cp2kBasis, Cp2kDispersion,
    Cp2kScf, Cp2kPbc, Cp2kCubePrint, Cp2kAtomCharge, Cp2kThermostat,
    QeJobBuilder, QeTask, QeFunctional, QeSmearing,
};
use std::path::PathBuf;

#[derive(Parser)]
#[command(
    name = "fe-job",
    about = "Generate QC software input files  (gaussian | cp2k)",
    disable_help_flag = true,
)]
struct Cli {
    /// Show help: overview when -s is absent, software-specific when -s is given
    #[arg(short = 'h', long = "help", action = clap::ArgAction::SetTrue)]
    help: bool,

    /// Input structure file (omit to show software-specific help)
    #[arg(short, long)]
    input: Option<PathBuf>,

    /// Target software: gaussian | cp2k  (omit to see overview)
    #[arg(short, long)]
    software: Option<String>,

    /// Output file (default: job.gjf / job.inp)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Use LAMMPS metal units for dump files
    #[arg(long)]
    metal_units: bool,

    // ── Charge / spin (shared) ───────────────────────────────────────────────
    /// Override total system charge
    #[arg(long)]
    charge: Option<i32>,

    /// Override spin multiplicity 2S+1 (1=singlet, 2=doublet, …)
    #[arg(long)]
    multiplicity: Option<u32>,

    /// Auto-guess multiplicity from structure (magmom → oxidation state → parity)
    #[arg(long)]
    auto_spin: bool,

    // ── Gaussian ─────────────────────────────────────────────────────────────
    /// DFT method  [gaussian]
    #[arg(short = 'm', long)]
    method: Option<String>,

    /// Basis set  [gaussian]
    #[arg(short = 'b', long)]
    basis: Option<String>,

    // ── CP2K shared ──────────────────────────────────────────────────────────
    /// Task type  [cp2k]  energy|force|geo-opt|cell-opt|md|freq
    #[arg(long, default_value = "energy")]
    task: String,

    /// DFT functional  [cp2k]  pbe|blyp|pbe0|b3lyp|revpbe|pbesol|scan|r2scan|hse06
    #[arg(long, default_value = "pbe")]
    functional: String,

    /// Basis set  [cp2k]  dzvp-molopt-sr|tzvp-molopt|tzv2p-molopt|dzvp-gth|tzvp-gth
    #[arg(long, default_value = "dzvp-molopt-sr")]
    cp2k_basis: String,

    /// Dispersion correction  [cp2k]  none|d3|d3bj
    #[arg(long, default_value = "none")]
    dispersion: String,

    /// SCF solver  [cp2k]  diag|ot
    #[arg(long, default_value = "diag")]
    scf: String,

    /// PBC direction  [cp2k]  xyz|z|none  (auto-detected from cell if omitted)
    #[arg(long)]
    pbc: Option<String>,

    /// k-point mesh  [cp2k]  e.g. --kpoints 2 2 2
    #[arg(long, num_args = 3, value_names = ["K1","K2","K3"])]
    kpoints: Option<Vec<u32>>,

    /// Plane-wave cutoff [Ry]  [cp2k]
    #[arg(long, default_value = "400")]
    cutoff: u32,

    /// Relative cutoff [Ry]  [cp2k]
    #[arg(long, default_value = "50")]
    rel_cutoff: u32,

    /// Enable Fermi-Dirac smearing  [cp2k]
    #[arg(long)]
    smear: bool,

    // ── CP2K output ───────────────────────────────────────────────────────────
    /// Atomic charge scheme  [cp2k]  none|mulliken|hirshfeld|hirshfeld-i
    #[arg(long, default_value = "none")]
    atom_charge: String,

    /// Export cube file  [cp2k]  none|density|elf|hartree
    #[arg(long, default_value = "none")]
    cube: String,

    /// Export Molden file  [cp2k]
    #[arg(long)]
    molden: bool,

    /// Project name  [cp2k]
    #[arg(long, default_value = "ferro")]
    project: String,

    // ── CP2K MD ───────────────────────────────────────────────────────────────
    /// MD steps  [cp2k]
    #[arg(long, default_value = "10000")]
    md_steps: u32,

    /// MD timestep [fs]  [cp2k]
    #[arg(long, default_value = "1.0")]
    md_timestep: f64,

    /// MD temperature [K]  [cp2k]
    #[arg(long, default_value = "298.15")]
    temperature: f64,

    /// Thermostat  [cp2k]  csvr|nose|langevin|none
    #[arg(long, default_value = "csvr")]
    thermostat: String,

    /// MD trajectory write frequency  [cp2k]
    #[arg(long, default_value = "100")]
    traj_freq: u32,

    /// Enable NPT barostat  [cp2k]
    #[arg(long)]
    barostat: bool,

    // ── Quantum ESPRESSO ─────────────────────────────────────────────────────
    /// Calculation type  [qe]  scf|nscf|bands|relax|vc-relax|md|vc-md
    #[arg(long, default_value = "scf")]
    qe_task: String,

    /// DFT functional  [qe]  pbe|pbesol|revpbe|blyp|scan|r2scan|pbe0|hse06
    #[arg(long, default_value = "pbe")]
    qe_functional: String,

    /// Plane-wave cutoff ecutwfc [Ry]  [qe]
    #[arg(long, default_value = "50")]
    ecutwfc: f64,

    /// Smearing  [qe]  none|gaussian|mp|mv|fd
    #[arg(long, default_value = "none")]
    smearing: String,

    /// Pseudopotential directory  [qe]
    #[arg(long, default_value = "./pseudo")]
    pseudo_dir: String,
}

fn main() -> Result<()> {
    let args = Cli::parse();

    // 无 -s → 概览
    let software = match args.software.as_deref() {
        Some(s) => s.to_lowercase(),
        None => {
            print_fe_job_overview();
            return Ok(());
        }
    };

    // -h 或无 -i → 软件专属帮助
    if args.help || args.input.is_none() {
        print_job_help(&software);
        return Ok(());
    }

    let input = args.input.as_ref().unwrap();
    let units = if args.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let traj = read_trajectory(input, units)?;
    let mut frame = traj.frames.into_iter().next()
        .ok_or_else(|| anyhow!("No frames in input file"))?;

    // 电荷覆盖（先于自旋推断，因 guess_spin 依赖电荷）
    if let Some(c) = args.charge {
        frame.charge = c;
    }

    // 自旋多重度：显式 --multiplicity > --auto-spin > 文件值
    if let Some(m) = args.multiplicity {
        frame.multiplicity = m;
    } else if args.auto_spin {
        let g = guess_spin(&frame);
        let method = match g.method {
            SpinMethod::Magmom => "magmom 求和",
            SpinMethod::OxidationState => "氧化态 + Hund 规则",
            SpinMethod::Parity => "电子数奇偶下限",
        };
        println!(
            "Auto-spin: multiplicity = {} (未成对电子 {}, 方法: {})",
            g.multiplicity, g.n_unpaired, method
        );
        if let Some(ox) = &g.oxidation_states {
            let s: Vec<String> = ox.iter().map(|(e, o)| format!("{e}{o:+}")).collect();
            println!("           氧化态: {}", s.join(" "));
        }
        for w in &g.warnings {
            eprintln!("[warn] {w}");
        }
        frame.multiplicity = g.multiplicity;
    }

    match software.as_str() {
        "gaussian" => {
            let mut builder = GaussianJobBuilder::new(frame);
            if let Some(m) = args.method { builder.method = m; }
            if let Some(b) = args.basis  { builder.basis_set = b; }
            let content = builder.build()?;
            let out = args.output.as_deref()
                .map(|p| p.to_string_lossy().into_owned())
                .unwrap_or_else(|| "job.gjf".to_string());
            std::fs::write(&out, content)?;
            println!("Gaussian input written to: {out}");
        }
        "cp2k" => {
            let mut builder = Cp2kJobBuilder::new(frame);
            builder.project = args.project.clone();
            // 显式 --multiplicity 时关闭 builder 自动推断，尊重手动值；
            // 否则（--auto-spin 或缺省）由 builder 经 guess_spin 推断。
            builder.auto_spin = args.multiplicity.is_none();

            builder.task = match args.task.as_str() {
                "energy"   => Cp2kTask::Energy,
                "force"    => Cp2kTask::ForceEnergy,
                "geo-opt"  => Cp2kTask::GeoOpt,
                "cell-opt" => Cp2kTask::CellOpt,
                "md"       => Cp2kTask::MD,
                "freq"     => Cp2kTask::Frequency,
                other      => bail!("Unknown task: {other}  (energy|force|geo-opt|cell-opt|md|freq)"),
            };

            builder.functional = match args.functional.as_str() {
                "pbe"     => Cp2kFunctional::PBE,
                "blyp"    => Cp2kFunctional::BLYP,
                "pbe0"    => Cp2kFunctional::PBE0,
                "b3lyp"   => Cp2kFunctional::B3LYP,
                "revpbe"  => Cp2kFunctional::RevPBE,
                "pbesol"  => Cp2kFunctional::PBEsol,
                "scan"    => Cp2kFunctional::SCAN,
                "r2scan"  => Cp2kFunctional::R2SCAN,
                "hse06"   => Cp2kFunctional::HSE06,
                other     => Cp2kFunctional::Custom(other.to_string()),
            };

            builder.basis = match args.cp2k_basis.as_str() {
                "dzvp-molopt-sr" => Cp2kBasis::DzvpMoloptSr,
                "tzvp-molopt"    => Cp2kBasis::TzvpMolopt,
                "tzv2p-molopt"   => Cp2kBasis::Tzv2pMolopt,
                "dzvp-gth"       => Cp2kBasis::DzvpGth,
                "tzvp-gth"       => Cp2kBasis::TzvpGth,
                "pob-dzvp"       => Cp2kBasis::PobDzvp,
                "pob-tzvp"       => Cp2kBasis::PobTzvp,
                other            => Cp2kBasis::Custom(other.to_string()),
            };

            builder.dispersion = match args.dispersion.as_str() {
                "none" => Cp2kDispersion::None,
                "d3"   => Cp2kDispersion::D3,
                "d3bj" => Cp2kDispersion::D3BJ,
                other  => bail!("Unknown dispersion: {other}  (none|d3|d3bj)"),
            };

            builder.scf = match args.scf.as_str() {
                "diag" => Cp2kScf::Diag,
                "ot"   => Cp2kScf::OT,
                other  => bail!("Unknown scf: {other}  (diag|ot)"),
            };

            if let Some(pbc) = &args.pbc {
                builder.pbc = match pbc.to_lowercase().as_str() {
                    "xyz"  => Cp2kPbc::XYZ,
                    "z"    => Cp2kPbc::Z,
                    "none" => Cp2kPbc::None,
                    other  => bail!("Unknown pbc: {other}  (xyz|z|none)"),
                };
            }

            if let Some(kp) = args.kpoints {
                builder.kpoints = Some([kp[0], kp[1], kp[2]]);
            }

            builder.cutoff     = args.cutoff;
            builder.rel_cutoff = args.rel_cutoff;
            builder.smear      = args.smear;

            builder.atom_charge = match args.atom_charge.as_str() {
                "none"        => Cp2kAtomCharge::None,
                "mulliken"    => Cp2kAtomCharge::Mulliken,
                "hirshfeld"   => Cp2kAtomCharge::Hirshfeld,
                "hirshfeld-i" => Cp2kAtomCharge::HirshfeldI,
                other         => bail!("Unknown atom-charge: {other}"),
            };

            builder.cube_print = match args.cube.as_str() {
                "none"    => Cp2kCubePrint::None,
                "density" => Cp2kCubePrint::Density,
                "elf"     => Cp2kCubePrint::ELF,
                "hartree" => Cp2kCubePrint::HartreePot,
                other     => bail!("Unknown cube: {other}  (none|density|elf|hartree)"),
            };

            builder.export_molden   = args.molden;
            builder.md.steps        = args.md_steps;
            builder.md.timestep     = args.md_timestep;
            builder.md.temperature  = args.temperature;
            builder.md.traj_freq    = args.traj_freq;
            builder.md.barostat     = args.barostat;
            builder.md.thermostat   = match args.thermostat.as_str() {
                "csvr"     => Cp2kThermostat::CSVR,
                "nose"     => Cp2kThermostat::NoseHoover,
                "langevin" => Cp2kThermostat::Langevin,
                "none"     => Cp2kThermostat::None,
                other      => bail!("Unknown thermostat: {other}  (csvr|nose|langevin|none)"),
            };

            let content = builder.build()?;
            let out = args.output.as_deref()
                .map(|p| p.to_string_lossy().into_owned())
                .unwrap_or_else(|| "job.inp".to_string());
            std::fs::write(&out, content)?;
            println!("CP2K input written to: {out}");
        }
        "qe" | "espresso" | "pwscf" => {
            let mut builder = QeJobBuilder::new(frame);
            builder.auto_spin = args.multiplicity.is_none();
            builder.ecutwfc = args.ecutwfc;
            builder.pseudo_dir = args.pseudo_dir.clone();
            builder.task = match args.qe_task.as_str() {
                "scf"      => QeTask::Scf,
                "nscf"     => QeTask::Nscf,
                "bands"    => QeTask::Bands,
                "relax"    => QeTask::Relax,
                "vc-relax" => QeTask::VcRelax,
                "md"       => QeTask::Md,
                "vc-md"    => QeTask::VcMd,
                other      => bail!("Unknown qe-task: {other}  (scf|nscf|bands|relax|vc-relax|md|vc-md)"),
            };
            builder.functional = match args.qe_functional.as_str() {
                "pbe"    => QeFunctional::PBE,
                "pbesol" => QeFunctional::PBEsol,
                "revpbe" => QeFunctional::RevPBE,
                "blyp"   => QeFunctional::BLYP,
                "scan"   => QeFunctional::SCAN,
                "r2scan" => QeFunctional::R2SCAN,
                "pbe0"   => QeFunctional::PBE0,
                "hse06"  => QeFunctional::HSE06,
                other    => QeFunctional::Custom(other.to_uppercase()),
            };
            builder.smearing = match args.smearing.as_str() {
                "none"     => QeSmearing::None,
                "gaussian" => QeSmearing::Gaussian,
                "mp"       => QeSmearing::MethfesselPaxton,
                "mv"       => QeSmearing::MarzariVanderbilt,
                "fd"       => QeSmearing::FermiDirac,
                other      => bail!("Unknown smearing: {other}  (none|gaussian|mp|mv|fd)"),
            };
            if let Some(kp) = args.kpoints {
                builder.kpoints = Some([kp[0], kp[1], kp[2]]);
            }
            builder.md.steps = args.md_steps;
            builder.md.temperature = args.temperature;

            let content = builder.build()?;
            let out = args.output.as_deref()
                .map(|p| p.to_string_lossy().into_owned())
                .unwrap_or_else(|| "pw.in".to_string());
            std::fs::write(&out, content)?;
            println!("Quantum ESPRESSO input written to: {out}");
        }
        other => bail!("Unsupported software: {other}  (gaussian | cp2k | qe)"),
    }

    Ok(())
}
