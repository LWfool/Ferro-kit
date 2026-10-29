//! The CLI side of reading: the format itself is recognised by `ferro_io::Format`;
//! this file adds what is particular to the command line — errors that name the
//! flags, the reader warnings printed to stderr, and the help-page format matrix.

use std::path::Path;
use anyhow::{bail, Result};
use ferro_core::Trajectory;
use ferro_io::{Format, ReadOptions};

use crate::args::common::ReadArgs;

/// Fails when the reading options do not fit `path`: a LAMMPS data file with no
/// `--atom-style`, or `--type-map` on anything but a LAMMPS dump. The format is
/// `--format-in`'s when given, else the one the name declares.
///
/// The `Atoms` layout is not recoverable from the file: the `# style` comment is
/// optional and charge/molecular share a column count, so it is never guessed.
/// Batch commands call this on every input before reading the first one.
pub fn check_read_args(path: &Path, opts: &ReadArgs) -> Result<()> {
    let format = opts.format_in.or_else(|| Format::from_name(path));
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
    let (traj, report) = match opts.format_in {
        Some(f) => ferro_io::read_as(path, f, &options)?,
        None => ferro_io::read(path, &options)?,
    };
    crate::report::notices(path, &report.warnings);
    Ok(traj)
}

/// Returns the format matrix of `ferro convert`'s page: how each format is
/// recognised, its `--format-in` name(s), and what writing it gives.
///
/// The rows come from `Format::ALL`, so the matrix cannot drift from what is read
/// and written. Every format is read; the inputs and the program outputs are not
/// written, so writing gets its own column.
pub fn supported_formats() -> String {
    let mut rows = String::new();
    for f in Format::ALL {
        // CP2K's two logs share a row and a name, and have a label each
        let labels: Vec<&str> = Format::LABELS.iter()
            .filter(|(_, g)| g.name() == f.name()).map(|(l, _)| *l).collect();
        let write = match f {
            _ if !f.can_write() => "-",
            _ if f.multi_frame() => "all",
            _ => "one per file",
        };
        rows.push_str(&format!("  {:<14}{:<28}{:<18}{write}\n", f.name(), f.detected_by(), labels.join("  ")));
    }
    format!(
        "  Format        Detected by                 --format-in       Write
  ------------------------------------------------------------------------
{rows}
  Without --format-in the format comes from the file NAME: an extension in
  any case, or a POSCAR/CONTCAR prefix. Then from the content:
  * a .xyz whose line 2 declares Lattice= or Properties= is extended XYZ;
    writing .xyz still gives plain XYZ, use extxyz to keep the cell
  - a name that says nothing (OUTCAR, vasprun.xml, .out) is read by its
    program banner; frames whose SCF did not converge are dropped
  - a directory holding type.raw is a DeePMD system
  Write: `all` = one file holds the trajectory; `one per file` = one file
  per frame (thin with --stride / --number); `-` = read-only, generate a
  CP2K / QE input with `ferro-py job -s cp2k|qe`."
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

        // --format-in overrides the name in both directions
        let data = ReadArgs { format_in: Some(Format::LammpsData), ..Default::default() };
        assert!(check_read_args(Path::new("t.txt"), &data).unwrap_err().to_string().contains("--atom-style"));
        let xyz = ReadArgs { format_in: Some(Format::Xyz), ..Default::default() };
        assert!(check_read_args(Path::new("t.data"), &xyz).is_ok());
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
