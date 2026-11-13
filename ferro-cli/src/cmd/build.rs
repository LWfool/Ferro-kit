//! `ferro build` — make a new structure from existing ones (or from nothing, `box`).
//!
//! One input, one frame, one product: `<-o dir>/<op>_<input>[_<suffix>].<ext>` (merge takes
//! two inputs and names both), in the input's format unless `--format-out` says otherwise. The structure work is `ferro-structure`'s;
//! this module only reads, names and writes.

use std::path::{Path, PathBuf};

use anyhow::{bail, Result};
use clap::{Args, Subcommand};

use crate::args::common::ReadArgs;
use crate::batch;
use crate::io_dispatch::{check_read_args, format_out, read_trajectory};
use ferro_core::{Frame, Trajectory};
use ferro_io::{Format, WriteOptions};
use ferro_structure::{
    add_vacuum, build_box, find_supercell_dims, make_supercell, merge_frames, Component, MAX_SUPERCELL_ATOMS,
};

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
    pub format_out: Option<String>,
}

impl OutArgs {
    fn out(&self) -> batch::Output {
        batch::Output { dir: self.output.clone(), label: None, suffix: self.suffix.clone(), mkdir: self.mkdir }
    }

    /// The extension of the product made from `input`: `--format-out`, else that of the format
    /// the input is read as (`--format-in`, else what `Format::detect` finds, so an extended XYZ
    /// `.xyz` gives `.extxyz`). The input's own extension stays when its name declares that
    /// same format (`.data` stays `.data`); `POSCAR` / `CONTCAR` give `vasp`. Checked to be
    /// writable, so a read-only input format (OUTCAR, a QE input) fails before anything is read.
    fn ext(&self, input: &Path, read: &ReadArgs) -> Result<String> {
        if let Some(f) = &self.format_out {
            return Ok(format_out(f)?.0);
        }
        // detect looks inside: on an input that cannot be opened fall back on its name, so the
        // failure is reported by the read with its own message, not here
        let format = read.format_in
            .or_else(|| Format::detect(input).ok())
            .or_else(|| Format::from_name(input));
        match format.and_then(|f| Some((f, f.extension()?))) {
            Some((f, _)) if f != Format::Poscar && Format::from_name(input) == Some(f) => {
                Ok(input.extension().map(|e| e.to_string_lossy().to_lowercase()).unwrap_or_default())
            }
            Some((_, ext)) => Ok(ext.to_string()),
            None => bail!(
                "{}: ferro does not write this input's format, so the product needs --format-out EXT \
                 (vasp, cif, extxyz, ...)",
                input.display()
            ),
        }
    }
}

#[derive(Args, Debug)]
pub struct SupercellCmd {
    /// Input structure file, one frame
    // Any count reaches exact_inputs, which names the files when a shell glob gave several
    #[arg(short, long, num_args = 1.., value_name = "FILE")]
    pub input: Vec<PathBuf>,

    /// Replicas along a, b, c
    #[arg(long, num_args = 3, value_names = ["NX", "NY", "NZ"],
          value_parser = at_least_one,
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
    /// Input structure file, one frame
    // Any count reaches exact_inputs, which names the files when a shell glob gave several
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
    #[arg(short, long, num_args = 1.., value_names = ["A", "B"])]
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
    // The same default as `dataset --seed`, by name (review 2026-10-30 A5-L5)
    #[arg(long, default_value_t = ferro_analysis::ml::dataset::DEFAULT_SEED)]
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

/// A `--dims` replica count: clap refuses 0 as a usage error (review 2026-11-06 D15); its own
/// range parser would print the bound as `1..18446744073709551615`
fn at_least_one(v: &str) -> std::result::Result<usize, String> {
    match v.parse::<usize>() {
        Ok(n) if n >= 1 => Ok(n),
        _ => Err(format!("expects a whole number >= 1, got `{v}`")),
    }
}

pub fn run(cmd: &BuildCmd) -> Result<()> {
    match cmd {
        BuildCmd::Supercell(c) => run_supercell(c),
        BuildCmd::Vacuum(c) => run_vacuum(c),
        BuildCmd::Merge(c) => run_merge(c),
        BuildCmd::Box(c) => run_box(c),
    }
}

/// Writes one product, `<-o>/<segments…>[_<suffix>].<ext>`, and returns the line to print.
fn write_product(out_args: &OutArgs, segments: &[&str], ext: &str, made: &Trajectory, units: Option<ferro_io::LammpsUnits>) -> Result<String> {
    let out = out_args.out();
    let path = out.join(&batch::product_name(segments, out.suffix.as_deref(), ext));
    let opts = WriteOptions { units, drop_partial: false };
    // One frame always lands in `path` itself: `write` numbers files only for several frames
    let report = ferro_io::write(&made.canonically_sorted(), &path, &opts)?;
    crate::report::notices(&path, &report.warnings);
    Ok(path.display().to_string())
}

/// Reads `path` as one structure. A structure operation builds a starting structure, so a
/// trajectory is refused rather than worked frame by frame, and which frame to take is the
/// user's to say (user, 2026-11-07; merge refused it before). The one-frame trajectory is
/// returned whole: the writers take the file's title from its metadata.
fn read_structure(path: &Path, read: &ReadArgs, op: &str) -> Result<Trajectory> {
    let traj = read_trajectory(path, read)?;
    if traj.n_frames() != 1 {
        bail!(
            "{} holds {} frames; {op} works on one structure — take one first with \
             `ferro convert -i {} -o DIR --format-out extxyz --start K --end K`",
            path.display(), traj.n_frames(), path.display()
        );
    }
    Ok(traj)
}

/// Reads the one input of `op` as one structure, hands it to `make`, and writes what comes
/// back as `<op>_<input>[_<suffix>].<ext>`, atoms in (Z, symbol) order as `convert` writes
/// them. `make` returns the structure and a note to print; `dry_run` prints the note only,
/// so it needs neither an output format nor the `-o` directory (review 2026-11-10 M8).
fn one_input(
    op: &str,
    inputs: &[PathBuf],
    read: &ReadArgs,
    out_args: &OutArgs,
    dry_run: bool,
    make: impl Fn(&Frame) -> Result<(Frame, String)>,
) -> Result<()> {
    let input = &batch::exact_inputs(&format!("build {op}"), inputs, 1)?[0];
    // Parameter errors (-s, read arguments, output format) are checked before any directory is made or file read
    batch::check_suffix(out_args.suffix.as_deref())?;
    check_read_args(input, read)?;
    let ext = (!dry_run).then(|| out_args.ext(input, read)).transpose()?;
    let label = &batch::unique_labels(inputs)[0];
    if ext.is_some() {
        out_args.out().prepare()?;
    }

    let traj = read_structure(input, read, op)?;
    let (made, note) = make(&traj.frames[0])?;
    let Some(ext) = ext else {
        println!("{note}");
        return Ok(());
    };
    let made = Trajectory { frames: vec![made], metadata: traj.metadata };
    println!("{note} -> {}", write_product(out_args, &[op, label], &ext, &made, read.units)?);
    Ok(())
}

fn run_supercell(c: &SupercellCmd) -> Result<()> {
    if c.dims.is_none() && c.min_atoms.is_none() && c.min_length.is_none() {
        bail!("give the replicas (--dims NX NY NZ) or a size to reach (--min-atoms N, --min-length LEN)");
    }
    if let Some(l) = c.min_length {
        if !(l >= 0.0 && l.is_finite()) {
            bail!("--min-length must be a finite number >= 0, got {l}");
        }
    }

    // find_supercell_dims refuses it too, but only after the input is read
    if let Some(n) = c.min_atoms {
        if n > MAX_SUPERCELL_ATOMS {
            bail!("--min-atoms must be at most {MAX_SUPERCELL_ATOMS}, got {n}");
        }
    }

    one_input("supercell", &c.input, &c.read, &c.out, c.dry_run, |first| {
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
        let first_sc = make_supercell(first, nx, ny, nz)?;
        // The interplanar spacing bounds a minimum-image cutoff (rcut <= d/2); on an oblique cell it is shorter than the edge
        let d = first_sc.cell.as_ref().expect("make_supercell keeps the cell").interplanar_spacings()?;
        let note = format!(
            "{nx}x{ny}x{nz} ({how}): {n} -> {} atoms, spacings {:.3} {:.3} {:.3} Å{}",
            n * nx * ny * nz, d[0], d[1], d[2],
            if dims == [1, 1, 1] { ", already large enough" } else { "" }
        );
        Ok((first_sc, note))
    })
}

fn run_vacuum(c: &VacuumCmd) -> Result<()> {
    let Some(t) = c.thickness else { bail!("give the vacuum thickness: --thickness LEN") };
    // add_vacuum 也查，但它在读完文件之后才查到
    if !(t > 0.0 && t.is_finite()) {
        bail!("--thickness must be a finite number > 0, got {t}");
    }
    one_input("vacuum", &c.input, &c.read, &c.out, false, |frame| {
        Ok((add_vacuum(frame, &c.axis, t)?, format!("{t} Å of vacuum along {}", c.axis)))
    })
}

/// Not a batch: two inputs make one product, so a failure fails the run.
fn run_merge(c: &MergeCmd) -> Result<()> {
    // None goes to the help page (review 2026-10-30 A5-L7); any count but two is refused here
    let inputs = batch::exact_inputs("build merge", &c.input, 2)?;
    let (a, b) = (&inputs[0], &inputs[1]);
    let Some(gap) = c.gap else { bail!("give the gap between the two blocks: --gap LEN (0 allowed)") };
    if !(gap >= 0.0 && gap.is_finite()) {
        bail!("--gap must be a finite number >= 0, got {gap}");
    }
    batch::check_suffix(c.out.suffix.as_deref())?;
    // One --type-map reaches both inputs: it goes to the dump among them, so a dump without
    // element names merges with another format (review 2026-11-10 B6-6). Neither a dump,
    // check_read_args refuses it as for any other command
    let is_dump = |p: &Path| c.read.format_in.or_else(|| Format::from_name(p)) == Some(Format::LammpsDump);
    let read_for = |p: &Path| {
        let mut r = c.read.clone();
        if is_dump(a) || is_dump(b) {
            r.type_map = r.type_map.filter(|_| is_dump(p));
        }
        r
    };
    let (read_a, read_b) = (read_for(a), read_for(b));
    check_read_args(a, &read_a)?;
    check_read_args(b, &read_b)?;
    let ext = c.out.ext(a, &c.read)?;
    let labels = batch::unique_labels(&c.input);
    c.out.out().prepare()?;

    let one = |p: &Path, r: &ReadArgs| -> Result<Frame> { Ok(read_structure(p, r, "merge")?.frames.remove(0)) };
    let (fa, fb) = (one(a, &read_a)?, one(b, &read_b)?);
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
    let ext = match &c.out.format_out {
        Some(f) => format_out(f)?.0,
        None => "vasp".to_string(),
    };
    // 先建盒再建目录：密度、最小间距、化学式的错误都由 build_box 报，不该先留下一个空目录
    let frame = build_box(&c.component, rho, c.min_dist, c.relax_steps, c.seed)?;
    c.out.out().prepare()?;
    let side = frame.cell.as_ref().expect("build_box makes a cubic cell").matrix[(0, 0)];
    let n = frame.n_atoms();
    let line = write_product(&c.out, &["box"], &ext, &Trajectory::from_frame(frame), None)?;
    println!("Box {side:.3} Å, {n} atoms, {rho} g/cm³ -> {line}");
    Ok(())
}
