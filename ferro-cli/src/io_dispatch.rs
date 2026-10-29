//! The CLI side of reading: the format itself is recognised by `ferro_io::Format`;
//! this file adds what is particular to the command line — errors that name the
//! flags, the reader warnings printed to stderr, and the help-page format matrix.

use std::path::Path;
use anyhow::{bail, Result};
use ferro_core::Trajectory;
use ferro_io::{Format, ReadOptions};

use crate::args::common::ReadArgs;

/// Fails when the reading options do not fit `path`: a LAMMPS data file with no
/// `--atom-style`, or `--type-map` on anything but a LAMMPS dump.
///
/// The `Atoms` layout is not recoverable from the file: the `# style` comment is
/// optional and charge/molecular share a column count, so it is never guessed.
/// Batch commands call this on every input before reading the first one.
pub fn check_read_args(path: &Path, opts: &ReadArgs) -> Result<()> {
    let format = Format::from_name(path);
    if opts.type_map.is_some() && format != Some(Format::LammpsDump) {
        bail!("--type-map applies to LAMMPS dump input only, not {}", path.display());
    }
    if opts.atom_style.is_none() && format == Some(Format::LammpsData) {
        bail!(
            "{} is a LAMMPS data file: give its atom style with --atom-style atomic|charge|full \
             (Ferro never guesses it from the Atoms comment or the column count)",
            path.display()
        );
    }
    Ok(())
}

/// Reads any supported file or DeePMD directory, printing what the reader noticed.
pub fn read_trajectory(path: &Path, opts: &ReadArgs) -> Result<Trajectory> {
    check_read_args(path, opts)?;
    let options = ReadOptions {
        units: opts.units,
        atom_style: opts.atom_style,
        type_map: opts.type_map.as_ref().map(|m| m.0.clone()),
    };
    let (traj, report) = ferro_io::read(path, &options)?;
    crate::report::notices(path, &report.warnings);
    Ok(traj)
}

/// Returns a human-readable read/write matrix of the supported formats.
///
/// The rows come from `Format::ALL`, so the matrix cannot drift from what is read
/// and written. Read and write are NOT symmetric — the inputs and the program
/// outputs are read-only — so the two directions get their own column.
pub fn supported_formats() -> String {
    let mut rows = String::new();
    for f in Format::ALL {
        let frames = match f {
            _ if !f.can_write() => "-",
            Format::Pdb => "all (MODEL records)",
            _ if f.multi_frame() => "all",
            _ => "one per file",
        };
        let write = if f.can_write() { "y" } else { "-" };
        rows.push_str(&format!("  {:<14}{:<31}{:<7}{:<7}{frames}\n", f.name(), f.detected_by(), "y", write));
    }
    format!(
        "  Format        Detected by                    Read   Write  Frames on write
  ------------------------------------------------------------------
{rows}
  Format is taken from the file NAME first: an extension (any case), or a
  POSCAR/CONTCAR prefix for VASP files that have none.
  * A .xyz whose comment line (line 2) declares Lattice= or Properties= is
  read as extended XYZ, as ASE and GPUMD write it. Writing .xyz still gives
  plain XYZ; use .extxyz to keep the cell.
  A name that says nothing (OUTCAR, vasprun.xml, a .out log) is read by its
  program banner. Frames whose SCF did not converge are dropped and counted,
  as `ferro dataset collect` does. A directory holding type.raw is a DeePMD
  system. --format takes an extension from a row with Write = y.

  `-` under Write means read-only: an input or a program output can be
  converted FROM, not TO. To generate a CP2K or QE input use
  `ferro-py job -s cp2k|qe`, which writes a full run setup rather than bare
  coordinates.

  `one per file` means the format holds one structure: 500 frames of
  traj.dump with --format cif give traj_0000.cif ... traj_0499.cif. Use
  --stride / --number to thin them out."
    )
}


#[cfg(test)]
mod tests {
    use super::*;

    /// The reading rules the format table does not own: flag-named errors that fail
    /// before any file is read. Format recognition itself is tested in `ferro_io::format`.
    #[test]
    fn flag_errors_name_the_flag() {
        for name in ["t.data", "t.lmp", "T.LAMMPS"] {
            let err = check_read_args(Path::new(name), &ReadArgs::default()).unwrap_err().to_string();
            assert!(err.contains("--atom-style"), "{name}: {err}");
        }
        assert!(check_read_args(Path::new("t.lammpstrj"), &ReadArgs::default()).is_ok());
        assert!(check_read_args(Path::new("t.xyz"), &ReadArgs::default()).is_ok());

        let opts = ReadArgs { type_map: Some(crate::args::common::TypeMap(vec![(1, "P".into())])), ..Default::default() };
        assert!(check_read_args(Path::new("t.lammpstrj"), &opts).is_ok());
        let err = read_trajectory(Path::new("t.extxyz"), &opts).unwrap_err().to_string();
        assert!(err.contains("--type-map applies to LAMMPS dump input only"), "{err}");
    }

    #[test]
    fn the_matrix_lists_every_format_once() {
        let table = supported_formats();
        for f in Format::ALL {
            let rows = table.lines().filter(|l| l.trim_start().starts_with(f.name())
                && l.contains(f.detected_by())).count();
            assert_eq!(rows, 1, "{f:?}\n{table}");
        }
    }
}
