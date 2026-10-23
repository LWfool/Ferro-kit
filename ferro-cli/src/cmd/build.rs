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
use ferro_structure::{find_supercell_dims, make_supercell};

#[derive(Subcommand, Debug)]
pub enum BuildCmd {
    /// Replicate a cell nx × ny × nz times, given or sized from --min-atoms / --min-length
    Supercell(SupercellCmd),
}

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
            Format::for_output(Path::new(&format!("x.{ext}"))).with_context(|| format!("--format {f}"))?;
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

pub fn wants_help(cmd: &BuildCmd) -> bool {
    match cmd {
        BuildCmd::Supercell(c) => c.input.is_empty(),
    }
}

pub fn print_help(cmd: &BuildCmd) {
    match cmd {
        BuildCmd::Supercell(_) => crate::help::print_build_supercell(),
    }
}

/// Returns the number of inputs that failed, for the process exit code.
pub fn run(cmd: &BuildCmd) -> Result<usize> {
    match cmd {
        BuildCmd::Supercell(c) => run_supercell(c),
    }
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
    let out = out_args.out();
    out.prepare()?;

    let opts = WriteOptions { units: read.units, drop_partial: false };
    let (ok, failed) = batch::map_inputs(&inputs, |inp| {
        let traj = read_trajectory(&inp.path, read)?;
        let (made, note) = make(&traj)?;
        let Some(made) = made else {
            println!("  {note}");
            return Ok(());
        };
        let ext = out_args.ext(&inp.path)?;
        let path = out.join(&batch::product_name(&[op, &inp.label], out.suffix.as_deref(), &ext));
        let written = ferro_io::write(&made.canonically_sorted(), &path, &opts)?;
        match written.as_slice() {
            [one] => println!("  {note} -> {}", one.display()),
            many => println!("  {note} -> {} files  {} ... {}",
                             many.len(), many[0].display(), many[many.len() - 1].display()),
        }
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
