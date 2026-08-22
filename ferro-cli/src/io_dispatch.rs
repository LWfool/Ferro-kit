use std::path::Path;
use anyhow::{bail, Result};
use ferro_core::Trajectory;
use ferro_io::{self, LammpsUnits, *};

pub fn read_trajectory(path: &Path, lammps_units: LammpsUnits) -> Result<Trajectory> {
    let s = path.to_str().unwrap_or_default();
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or_default();
    let upper = name.to_uppercase();

    if upper.starts_with("POSCAR") {
        return read_poscar(s);
    }
    if upper.starts_with("CONTCAR") {
        return read_contcar(s);
    }

    match path.extension().and_then(|e| e.to_str()) {
        Some("xyz")                      => Ok(read_xyz(s)?),
        Some("pdb")                      => Ok(read_pdb(s)?),
        Some("cif")                      => Ok(read_cif(s)?),
        Some("extxyz")                   => Ok(read_extxyz(s)?),
        Some("lammps") | Some("data") | Some("lmp") => Ok(read_lammps_data(s)?),
        Some("dump") | Some("lammpstrj")             => Ok(read_lammps_dump(s, lammps_units)?),
        Some("inp")                      => Ok(read_cp2k_inp(s)?),
        Some("restart")                  => Ok(read_cp2k_restart(s)?),
        Some("in") | Some("qe")          => Ok(read_qe_input(s)?),
        Some(ext) => bail!("Unsupported input format: .{ext}"),
        None      => bail!("Cannot determine format (no extension): {s}"),
    }
}

/// Reads a trajectory and keeps only its last `n` frames when `last_n` is given.
///
/// The `--last-n` skip-equilibration step every analysis binary performs, in one place
/// so the batch loop stays a one-liner in each of them.
pub fn read_trajectory_tail(
    path: &Path,
    lammps_units: LammpsUnits,
    last_n: Option<usize>,
) -> Result<Trajectory> {
    let mut traj = read_trajectory(path, lammps_units)?;
    if let Some(n) = last_n {
        traj = traj.tail(n);
    }
    Ok(traj)
}

pub fn write_trajectory(traj: &Trajectory, path: &Path, lammps_units: LammpsUnits) -> Result<()> {
    let s = path.to_str().unwrap_or_default();
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or_default();
    let upper = name.to_uppercase();

    if upper.starts_with("POSCAR") || upper.starts_with("CONTCAR") {
        return write_poscar(traj, s);
    }

    match path.extension().and_then(|e| e.to_str()) {
        Some("xyz")                      => Ok(write_xyz(traj, s)?),
        Some("pdb")                      => Ok(write_pdb(traj, s)?),
        Some("cif")                      => Ok(write_cif(traj, s)?),
        Some("extxyz")                   => Ok(write_extxyz(traj, s)?),
        Some("lammps") | Some("data") | Some("lmp") => Ok(write_lammps_data(traj, s)?),
        Some("dump") | Some("lammpstrj") => Ok(write_lammps_dump(traj, s, lammps_units)?),
        Some("in") | Some("qe")          => Ok(write_qe_input(traj, s)?),
        Some(ext) => bail!("Unsupported output format: .{ext}"),
        None      => bail!("Cannot determine format (no extension): {s}"),
    }
}

/// Whether the format named by `path` can hold more than one frame on write.
///
/// The `Frames on write` column of [`supported_formats`] in code form: XYZ, extxyz,
/// PDB, CIF and LAMMPS dump carry a whole trajectory, while POSCAR, LAMMPS data and
/// QE input hold a single structure and would silently keep only frame 0. Callers
/// writing several frames use this to decide between one file and one file per frame.
///
/// Unknown extensions answer `true` so the write itself produces the "Unsupported
/// output format" error, rather than this function turning it into a pile of
/// per-frame failures.
pub fn holds_multiple_frames(path: &Path) -> bool {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or_default();
    let upper = name.to_uppercase();
    if upper.starts_with("POSCAR") || upper.starts_with("CONTCAR") {
        return false;
    }
    !matches!(
        path.extension().and_then(|e| e.to_str()),
        Some("lammps") | Some("data") | Some("lmp") | Some("in") | Some("qe")
    )
}

/// Returns a human-readable read/write matrix of the supported formats.
///
/// Read and write are NOT symmetric — the CP2K inputs are read-only — so the two
/// directions get their own column rather than one flat list.
pub fn supported_formats() -> &'static str {
    "  Format        Detected by                  Read   Write  Frames on write
  ----------------------------------------------------------------
  XYZ           .xyz                         y      y      all
  extended XYZ  .extxyz                      y      y      all
  PDB           .pdb                         y      y      all (MODEL records)
  CIF           .cif                         y      y      all (data blocks)
  LAMMPS dump   .dump  .lammpstrj            y      y      all
  VASP          POSCAR* / CONTCAR* (prefix)  y      y      FIRST only
  LAMMPS data   .lammps  .data  .lmp         y      y      FIRST only
  QE (pw.x)     .in  .qe                     y      y      FIRST only
  CP2K input    .inp                         y      -      -
  CP2K restart  .restart                     y      -      -

  Format is taken from the file NAME, never from a flag: extension for most,
  a POSCAR/CONTCAR prefix for VASP (case-insensitive). Writing to a CONTCAR
  name emits POSCAR-format content.

  `-` under Write means read-only: a CP2K .inp can be converted FROM, not TO.
  To generate CP2K input use `ferro job -s cp2k`, which writes a full run
  setup rather than bare coordinates.

  `FIRST only` means the format holds one structure: a 500-frame trajectory
  written to POSCAR gives you frame 0 and no warning."
}


#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame};
    use nalgebra::Vector3;

    fn one_atom_traj() -> Trajectory {
        let mut frame = Frame::new();
        frame.cell = Some(Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap());
        frame.pbc = [true; 3];
        frame.atoms.push(Atom::new("Si", Vector3::new(0.0, 0.0, 0.0)));
        Trajectory { frames: vec![frame], metadata: Default::default() }
    }

    /// `supported_formats()` 的 Write 列写着 `-` 的两个格式，dispatch 必须真的拒绝。
    /// 表与 match 分支是两处手写的事实，漂了就是文档在说谎。
    #[test]
    fn test_read_only_formats_are_rejected_on_write() {
        let traj = one_atom_traj();
        let dir = std::env::temp_dir();
        for ext in ["inp", "restart"] {
            let path = dir.join(format!("ferro_dispatch_test.{ext}"));
            let err = write_trajectory(&traj, &path, LammpsUnits::Real).unwrap_err();
            assert!(
                err.to_string().contains("Unsupported output format"),
                ".{ext} should be read-only, got: {err}"
            );
            assert!(!path.exists(), "拒绝写的格式不该留下空文件");
        }
    }

    #[test]
    fn test_format_table_marks_cp2k_read_only() {
        let table = supported_formats();
        for line in table.lines() {
            let Some(rest) = line.trim().strip_prefix("CP2K") else { continue };
            // 该行形如 `CP2K input    .inp    y    -    -`
            let cols: Vec<&str> = rest.split_whitespace().collect();
            assert_eq!(cols[2], "y", "CP2K 应可读: {line}");
            assert_eq!(cols[3], "-", "CP2K 应不可写: {line}");
        }
    }

    /// `holds_multiple_frames` 与格式表的 `Frames on write` 列是两处手写的同一事实
    #[test]
    fn test_multi_frame_capability_agrees_with_the_table() {
        for (name, expected) in [
            ("t.xyz", true), ("t.extxyz", true), ("t.pdb", true), ("t.cif", true),
            ("t.dump", true), ("t.lammpstrj", true),
            ("POSCAR", false), ("CONTCAR", false), ("poscar", false),
            ("run1_POSCAR", true),   // 前缀匹配，不是包含匹配
            ("t.lmp", false), ("t.data", false), ("t.lammps", false),
            ("t.in", false), ("t.qe", false),
        ] {
            assert_eq!(
                holds_multiple_frames(Path::new(name)),
                expected,
                "{name} 与 supported_formats() 的 Frames 列不一致"
            );
        }
    }

    /// 未知扩展名答 true，好让 write_trajectory 报一次「不支持的格式」，
    /// 而不是被拆成 N 次逐帧失败
    #[test]
    fn test_unknown_extension_defers_to_the_writer_error() {
        assert!(holds_multiple_frames(Path::new("t.nosuchfmt")));
    }

    #[test]
    fn test_unknown_extension_is_rejected_both_ways() {
        let traj = one_atom_traj();
        let path = std::env::temp_dir().join("ferro_dispatch_test.nosuchfmt");
        assert!(write_trajectory(&traj, &path, LammpsUnits::Real).is_err());
        assert!(read_trajectory(&path, LammpsUnits::Real).is_err());
    }
}
