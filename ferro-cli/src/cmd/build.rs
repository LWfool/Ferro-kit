//! `ferro build` — make a new structure from existing ones (or from nothing, `box`).
//!
//! One product per input, `<-o dir>/<op>_<input>[_<suffix>].<ext>`, in the input's
//! format unless `--format` says otherwise. The structure work is `ferro-structure`'s;
//! this module only reads, names and writes.

use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::{Args, Subcommand};

use crate::args::common::ReadArgs;
use crate::batch;
use crate::io_dispatch::read_trajectory;
use ferro_core::Trajectory;
use ferro_io::{Format, WriteOptions};
use ferro_structure::{add_vacuum, build_box, find_supercell_dims, make_supercell, merge_frames, Component};

#[derive(Subcommand, Debug)]
pub enum BuildCmd {
    /// Replicate a cell nx × ny × nz times, given or sized from --min-atoms / --min-length
    Supercell(SupercellCmd),
    /// Add a vacuum layer along one lattice vector
    Vacuum(VacuumCmd),
    /// Stack two structures along one lattice vector, with a gap between them
    Merge(MergeCmd),
    /// Pack molecules at random into a cubic box of a given density
    Box(BoxCmd),
}

/// The lattice vector a layer or a stack goes along.
const AXES: [&str; 3] = ["x", "y", "z"];

/// Where and how the products are written; shared by every `build` command.
#[derive(Args, Debug, Default)]
pub struct OutArgs {
    /// Directory to write the products into (default: current dir)
    #[arg(short, long, value_name = "DIR")]
    pub output: Option<PathBuf>,

    /// Batch tag, the last segment of every product name: <op>_<input>_<suffix>.<ext>
    #[arg(short, long, value_name = "SUFFIX")]
    pub suffix: Option<String>,

    /// Create the output directory without asking (required when there is no terminal)
    #[arg(long)]
    pub mkdir: bool,

    /// Output format as an extension: vasp, cif, extxyz, data… (default: the input's)
    #[arg(long, value_name = "EXT")]
    pub format: Option<String>,
}

impl OutArgs {
    fn out(&self) -> batch::Output {
        batch::Output { dir: self.output.clone(), label: None, suffix: self.suffix.clone(), mkdir: self.mkdir }
    }

    /// The extension of the product made from `input`: `--format`, else the input's own
    /// (`POSCAR` / `CONTCAR` → `vasp`). Checked to be writable, so a read-only input
    /// format (OUTCAR, a QE input) fails before anything is read.
    fn ext(&self, input: &Path) -> Result<String> {
        if let Some(f) = &self.format {
            let ext = f.trim_start_matches('.').to_lowercase();
            Format::for_extension(&ext).with_context(|| format!("--format {f}"))?;
            return Ok(ext);
        }
        match Format::from_name(input) {
            Some(Format::Poscar) => Ok("vasp".to_string()),
            Some(f) if f.can_write() => {
                Ok(input.extension().map(|e| e.to_string_lossy().to_lowercase()).unwrap_or_default())
            }
            _ => bail!(
                "{}: ferro does not write this input's format, so the product needs --format EXT \
                 (vasp, cif, extxyz, ...)",
                input.display()
            ),
        }
    }
}

#[derive(Args, Debug)]
pub struct SupercellCmd {
    /// Input structure file(s); glob patterns allowed — quote them
    #[arg(short, long, num_args = 1.., value_name = "FILE")]
    pub input: Vec<PathBuf>,

    /// Replicas along a, b, c
    #[arg(long, num_args = 3, value_names = ["NX", "NY", "NZ"],
          conflicts_with_all = ["min_atoms", "min_length"])]
    pub dims: Option<Vec<usize>>,

    /// Smallest supercell with at least this many atoms
    #[arg(long, value_name = "N")]
    pub min_atoms: Option<usize>,

    /// Smallest supercell whose interplanar spacings all reach this length [Å]
    #[arg(long, value_name = "LEN")]
    pub min_length: Option<f64>,

    /// Print the replicas and atom counts, write nothing
    #[arg(long)]
    pub dry_run: bool,

    #[command(flatten)]
    pub out: OutArgs,

    #[command(flatten)]
    pub read: ReadArgs,
}

#[derive(Args, Debug)]
pub struct VacuumCmd {
    /// Input structure file(s); glob patterns allowed — quote them
    #[arg(short, long, num_args = 1.., value_name = "FILE")]
    pub input: Vec<PathBuf>,

    /// Lattice vector to lengthen: x = a, y = b, z = c
    #[arg(long, default_value = "z", value_parser = clap::builder::PossibleValuesParser::new(AXES))]
    pub axis: String,

    /// Vacuum thickness [Å], measured perpendicular to the other two vectors
    #[arg(long, value_name = "LEN")]
    pub thickness: Option<f64>,

    #[command(flatten)]
    pub out: OutArgs,

    #[command(flatten)]
    pub read: ReadArgs,
}

#[derive(Args, Debug)]
pub struct MergeCmd {
    /// The two structures, bottom first
    #[arg(short, long, num_args = 2, value_names = ["A", "B"])]
    pub input: Vec<PathBuf>,

    /// Lattice vector to stack along: x = a, y = b, z = c
    #[arg(long, default_value = "z", value_parser = clap::builder::PossibleValuesParser::new(AXES))]
    pub axis: String,

    /// Gap between the two blocks [Å], measured along the interface normal
    #[arg(long, value_name = "LEN")]
    pub gap: Option<f64>,

    #[command(flatten)]
    pub out: OutArgs,

    #[command(flatten)]
    pub read: ReadArgs,
}

#[derive(Args, Debug)]
pub struct BoxCmd {
    /// One component, NAME:COUNT — a formula or a known name (H2O:100, 'Ca3(PO4)2:10'); repeat for more
    #[arg(short, long, value_name = "SPEC", value_parser = parse_component)]
    pub component: Vec<Component>,

    /// Mass density [g/cm³]; sets the cubic box length
    #[arg(long, value_name = "RHO")]
    pub density: Option<f64>,

    /// Closest allowed distance between two atoms [Å]
    #[arg(long, value_name = "LEN", default_value_t = 1.5)]
    pub min_dist: f64,

    /// Soft-core relaxation steps that push overlapping atoms apart
    #[arg(long, value_name = "N", default_value_t = 500)]
    pub relax_steps: usize,

    /// Random seed; the same seed gives the same box, atom for atom
    // 与 `dataset --seed` 同一默认值
    #[arg(long, default_value_t = 666)]
    pub seed: u64,

    #[command(flatten)]
    pub out: OutArgs,
}

fn parse_component(s: &str) -> std::result::Result<Component, String> {
    let (name, count) = s.rsplit_once(':').ok_or_else(|| format!("'{s}' is not NAME:COUNT (e.g. H2O:100)"))?;
    let n_molecules = count.parse::<usize>().map_err(|_| format!("'{count}' in '{s}' is not a count"))?;
    if name.is_empty() || n_molecules == 0 {
        return Err(format!("'{s}' needs a name and a count of at least 1"));
    }
    Ok(Component { compound: name.to_string(), n_molecules })
}

pub fn wants_help(cmd: &BuildCmd) -> bool {
    match cmd {
        BuildCmd::Supercell(c) => c.input.is_empty(),
        BuildCmd::Vacuum(c) => c.input.is_empty(),
        BuildCmd::Merge(c) => c.input.is_empty(),
        BuildCmd::Box(c) => c.component.is_empty(),
    }
}

pub fn print_help(cmd: &BuildCmd) {
    match cmd {
        BuildCmd::Supercell(_) => crate::help::print_build_supercell(),
        BuildCmd::Vacuum(_) => crate::help::print_build_vacuum(),
        BuildCmd::Merge(_) => crate::help::print_build_merge(),
        BuildCmd::Box(_) => crate::help::print_build_box(),
    }
}

/// Returns the number of inputs that failed, for the process exit code.
pub fn run(cmd: &BuildCmd) -> Result<usize> {
    match cmd {
        BuildCmd::Supercell(c) => run_supercell(c),
        BuildCmd::Vacuum(c) => run_vacuum(c),
        BuildCmd::Merge(c) => run_merge(c).map(|_| 0),
        BuildCmd::Box(c) => run_box(c).map(|_| 0),
    }
}

/// Writes one product, `<-o>/<segments…>[_<suffix>].<ext>`, and returns the line to print.
fn write_product(out_args: &OutArgs, segments: &[&str], ext: &str, made: &Trajectory, units: Option<ferro_io::LammpsUnits>) -> Result<String> {
    let out = out_args.out();
    let path = out.join(&batch::product_name(segments, out.suffix.as_deref(), ext));
    let opts = WriteOptions { units, drop_partial: false };
    let written = ferro_io::write(&made.canonically_sorted(), &path, &opts)?;
    Ok(match written.as_slice() {
        [one] => one.display().to_string(),
        many => format!("{} files  {} ... {}", many.len(), many[0].display(), many[many.len() - 1].display()),
    })
}

/// Reads each input, hands it to `make`, and writes what comes back as
/// `<op>_<input>[_<suffix>].<ext>`, atoms in (Z, symbol) order as `convert` writes
/// them. `make` returns `None` to write nothing (`--dry-run`) and a note for the
/// input's line.
fn per_input(
    op: &str,
    inputs: &[PathBuf],
    read: &ReadArgs,
    out_args: &OutArgs,
    make: impl Fn(&Trajectory) -> Result<(Option<Trajectory>, String)>,
) -> Result<usize> {
    // 参数级错误（-s、目标格式）在建目录、读第一个文件之前查完
    batch::check_suffix(out_args.suffix.as_deref())?;
    let inputs = batch::expand_trajectories(inputs)?;
    for p in &inputs {
        out_args.ext(p)?;
    }
    out_args.out().prepare()?;

    let (ok, failed) = batch::map_inputs(&inputs, |inp| {
        let traj = read_trajectory(&inp.path, read)?;
        let (made, note) = make(&traj)?;
        let Some(made) = made else {
            println!("  {note}");
            return Ok(());
        };
        let ext = out_args.ext(&inp.path)?;
        println!("  {note} -> {}", write_product(out_args, &[op, &inp.label], &ext, &made, read.units)?);
        Ok(())
    });
    if ok.is_empty() {
        bail!("every input failed; nothing to write");
    }
    Ok(failed.len())
}

fn run_supercell(c: &SupercellCmd) -> Result<usize> {
    if c.dims.is_none() && c.min_atoms.is_none() && c.min_length.is_none() {
        bail!("give the replicas (--dims NX NY NZ) or a size to reach (--min-atoms N, --min-length LEN)");
    }
    if let Some(d) = &c.dims {
        if d.contains(&0) {
            bail!("--dims takes three numbers of at least 1, got {d:?}");
        }
    }
    if let Some(l) = c.min_length {
        if !(l >= 0.0 && l.is_finite()) {
            bail!("--min-length must be a finite length >= 0, got {l}");
        }
    }

    per_input("supercell", &c.input, &c.read, &c.out, |traj| {
        let Some(first) = traj.first() else { bail!("no frames") };
        // 自动模式按第 0 帧定倍数，整条轨迹用同一组：每帧原子数不同的轨迹写不出去
        let (dims, how) = match &c.dims {
            Some(d) => ([d[0], d[1], d[2]], "given"),
            None => {
                let Some(cell) = &first.cell else { bail!("the structure has no cell") };
                let dims = find_supercell_dims(
                    cell, first.n_atoms(), c.min_length.unwrap_or(0.0), c.min_atoms.unwrap_or(0))?;
                (dims, "auto")
            }
        };
        let [nx, ny, nz] = dims;
        let n = first.n_atoms();
        let frames = traj.frames.iter().map(|f| make_supercell(f, nx, ny, nz)).collect::<Result<Vec<_>, _>>()?;
        // 面间距才是最小镜像截断的上界（r_cut ≤ d/2），斜胞上它比边长短
        let d = frames[0].cell.as_ref().map(|cl| cl.interplanar_spacings()).transpose()?.unwrap_or([0.0; 3]);
        let note = format!(
            "{nx}x{ny}x{nz} ({how}): {n} -> {} atoms, spacings {:.3} {:.3} {:.3} Å{}",
            n * nx * ny * nz, d[0], d[1], d[2],
            if dims == [1, 1, 1] { ", already large enough" } else { "" }
        );
        let made = Trajectory { frames, metadata: traj.metadata.clone() };
        Ok((if c.dry_run { None } else { Some(made) }, note))
    })
}

fn run_vacuum(c: &VacuumCmd) -> Result<usize> {
    let Some(t) = c.thickness else { bail!("give the vacuum thickness: --thickness LEN") };
    // add_vacuum 也查，但它在读完文件之后才查到
    if !(t > 0.0 && t.is_finite()) {
        bail!("--thickness must be a finite length > 0, got {t}");
    }
    per_input("vacuum", &c.input, &c.read, &c.out, |traj| {
        let frames = traj.frames.iter().map(|f| add_vacuum(f, &c.axis, t)).collect::<Result<Vec<_>, _>>()?;
        let note = format!("{t} Å of vacuum along {}, {} frame(s)", c.axis, frames.len());
        Ok((Some(Trajectory { frames, metadata: traj.metadata.clone() }), note))
    })
}

/// Not a batch: two inputs make one product, so a failure fails the run.
fn run_merge(c: &MergeCmd) -> Result<()> {
    let [a, b] = c.input.as_slice() else {
        bail!("merge takes exactly two structures: -i A B, got {}", c.input.len());
    };
    let Some(gap) = c.gap else { bail!("give the gap between the two blocks: --gap LEN (0 allowed)") };
    if !(gap >= 0.0 && gap.is_finite()) {
        bail!("--gap must be a finite length >= 0, got {gap}");
    }
    batch::check_suffix(c.out.suffix.as_deref())?;
    let ext = c.out.ext(a)?;
    let labels = batch::unique_labels(&c.input);
    c.out.out().prepare()?;

    // 只拼两个结构：多帧时取哪一帧该由用户说（convert --start），不默认取第 0 帧
    let one = |p: &Path| -> Result<ferro_core::Frame> {
        let traj = read_trajectory(p, &c.read)?;
        match traj.frames.len() {
            1 => Ok(traj.frames.into_iter().next().expect("one frame")),
            n => bail!(
                "{} holds {n} frames; merge stacks two structures — take one first with \
                 `ferro convert -i {} -o DIR --format extxyz --start K --end K`",
                p.display(), p.display()
            ),
        }
    };
    let (fa, fb) = (one(a)?, one(b)?);
    let merged = merge_frames(&fa, &fb, &c.axis, gap)?;
    let note = format!("{} + {} atoms along {}, gap {gap} Å", fa.n_atoms(), fb.n_atoms(), c.axis);
    let line = write_product(&c.out, &["merge", &labels[0], &labels[1]], &ext, &Trajectory::from_frame(merged), c.read.units)?;
    println!("Merged {note} -> {line}");
    Ok(())
}

fn run_box(c: &BoxCmd) -> Result<()> {
    let Some(rho) = c.density else { bail!("give the mass density: --density RHO (g/cm³)") };
    batch::check_suffix(c.out.suffix.as_deref())?;
    // 没有输入可沿用格式，默认写 POSCAR
    let ext = match &c.out.format {
        Some(_) => c.out.ext(Path::new(""))?,
        None => "vasp".to_string(),
    };
    // 先建盒再建目录：密度、最小间距、化学式的错误都由 build_box 报，不该先留下一个空目录
    let frame = build_box(&c.component, rho, c.min_dist, c.relax_steps, c.seed)?;
    c.out.out().prepare()?;
    let side = frame.cell.as_ref().map_or(0.0, |cl| cl.matrix[(0, 0)]);
    let n = frame.n_atoms();
    let line = write_product(&c.out, &["box"], &ext, &Trajectory::from_frame(frame), None)?;
    println!("Box {side:.3} Å, {n} atoms, {rho} g/cm³ -> {line}");
    Ok(())
}
