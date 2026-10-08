use std::path::Path;
use anyhow::{bail, Result};
use ferro_core::Trajectory;
use ferro_io::{self, LammpsUnits, *};

use crate::args::common::ReadArgs;

fn is_lammps_data(path: &Path) -> bool {
    matches!(path.extension().and_then(|e| e.to_str()), Some("lammps") | Some("data") | Some("lmp"))
}

/// Fails when `path` is a LAMMPS data file and no `--atom-style` was given.
///
/// The `Atoms` layout is not recoverable from the file: the `# style` comment is
/// optional and charge/molecular share a column count, so it is never guessed.
/// Batch commands call this on every input before reading the first one.
pub fn check_atom_style(path: &Path, style: Option<AtomStyle>) -> Result<Option<AtomStyle>> {
    if style.is_none() && is_lammps_data(path) {
        bail!(
            "{} is a LAMMPS data file: give its atom style with --atom-style atomic|charge|full \
             (Ferro never guesses it from the Atoms comment or the column count)",
            path.display()
        );
    }
    Ok(style)
}

pub fn read_trajectory(path: &Path, opts: &ReadArgs) -> Result<Trajectory> {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or_default();
    let upper = name.to_uppercase();

    if upper.starts_with("POSCAR") {
        return read_poscar(path);
    }
    if upper.starts_with("CONTCAR") {
        return read_contcar(path);
    }

    match path.extension().and_then(|e| e.to_str()) {
        // .xyz 也常装 extxyz（ASE、GPUMD、本仓 dataset --type nep），按第 2 行判
        Some("xyz") if is_extxyz(path)?  => Ok(read_extxyz(path)?),
        Some("xyz")                      => Ok(read_xyz(path)?),
        Some("pdb")                      => Ok(read_pdb(path)?),
        Some("cif")                      => Ok(read_cif(path)?),
        Some("extxyz")                   => Ok(read_extxyz(path)?),
        Some("vasp") | Some("pos")       => read_poscar(path),
        Some("lammps") | Some("data") | Some("lmp") => {
            let style = check_atom_style(path, opts.atom_style)?.expect("check_atom_style 已拒绝 None");
            Ok(read_lammps_data(path, style)?)
        }
        Some("dump") | Some("lammpstrj")             => Ok(read_lammps_dump(path, opts.units)?),
        Some("inp")                      => Ok(read_cp2k_inp(path)?),
        Some("restart")                  => Ok(read_cp2k_restart(path)?),
        Some("in") | Some("qe")          => Ok(read_qe_input(path)?),
        Some(ext) => bail!("Unsupported input format: .{ext}"),
        None      => bail!("Cannot determine format (no extension): {}", path.display()),
    }
}

/// Reads a trajectory and keeps only its last `n` frames when `last_n` is given.
///
/// The `--last-n` skip-equilibration step every analysis binary performs, in one place
/// so the batch loop stays a one-liner in each of them.
pub fn read_trajectory_tail(
    path: &Path,
    opts: &ReadArgs,
    last_n: Option<usize>,
) -> Result<Trajectory> {
    let mut traj = read_trajectory(path, opts)?;
    if let Some(n) = last_n {
        traj = traj.tail(n);
    }
    Ok(traj)
}

/// Write-side format, resolved from the file name alone.
pub enum OutFormat { Poscar, Xyz, Pdb, Cif, Extxyz, LammpsData, LammpsDump, Qe }

/// Resolves the write format without touching the disk, so a caller can reject a bad
/// `-o` before reading a multi-GB input (`convert`); `write_trajectory` goes through it.
pub fn out_format(path: &Path) -> Result<OutFormat> {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or_default();
    let upper = name.to_uppercase();

    if upper.starts_with("POSCAR") || upper.starts_with("CONTCAR") {
        return Ok(OutFormat::Poscar);
    }

    Ok(match path.extension().and_then(|e| e.to_str()) {
        Some("xyz")                      => OutFormat::Xyz,
        Some("pdb")                      => OutFormat::Pdb,
        Some("cif")                      => OutFormat::Cif,
        Some("extxyz")                   => OutFormat::Extxyz,
        Some("vasp") | Some("pos")       => OutFormat::Poscar,
        Some("lammps") | Some("data") | Some("lmp") => OutFormat::LammpsData,
        Some("dump") | Some("lammpstrj") => OutFormat::LammpsDump,
        Some("in") | Some("qe")          => OutFormat::Qe,
        Some(ext) => bail!("Unsupported output format: .{ext}"),
        None      => bail!("Cannot determine format (no extension): {}", path.display()),
    })
}

pub fn write_trajectory(traj: &Trajectory, path: &Path, lammps_units: Option<LammpsUnits>) -> Result<()> {
    match out_format(path)? {
        OutFormat::Poscar     => write_poscar(traj, path),
        OutFormat::Xyz        => Ok(write_xyz(traj, path)?),
        OutFormat::Pdb        => Ok(write_pdb(traj, path)?),
        OutFormat::Cif        => Ok(write_cif(traj, path)?),
        OutFormat::Extxyz     => Ok(write_extxyz(traj, path)?),
        OutFormat::LammpsData => Ok(write_lammps_data(traj, path)?),
        OutFormat::LammpsDump => Ok(write_lammps_dump(traj, path, lammps_units)?),
        OutFormat::Qe         => Ok(write_qe_input(traj, path)?),
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
            | Some("vasp") | Some("pos")
    )
}

/// Returns a human-readable read/write matrix of the supported formats.
///
/// Read and write are NOT symmetric — the CP2K inputs are read-only — so the two
/// directions get their own column rather than one flat list.
pub fn supported_formats() -> &'static str {
    "  Format        Detected by                    Read   Write  Frames on write
  ------------------------------------------------------------------
  XYZ           .xyz                           y      y      all
  extended XYZ  .extxyz  .xyz*                 y      y      all
  PDB           .pdb                           y      y      all (MODEL records)
  CIF           .cif                           y      y      all (data blocks)
  LAMMPS dump   .dump  .lammpstrj              y      y      all
  VASP          .vasp  .pos  POSCAR*/CONTCAR*  y      y      FIRST only
  LAMMPS data   .lammps  .data  .lmp           y      y      FIRST only
  QE (pw.x)     .in  .qe                       y      y      FIRST only
  CP2K input    .inp                           y      -      -
  CP2K restart  .restart                       y      -      -

  Format is taken from the file NAME, never from a flag: an extension, or a
  POSCAR/CONTCAR prefix (case-insensitive) for VASP files that have none.
  * The one look inside: a .xyz whose comment line (line 2) declares Lattice=
  or Properties= is read as extended XYZ, as ASE and GPUMD write it. Writing
  .xyz still gives plain XYZ; use .extxyz to keep the cell.
  Writing to a CONTCAR name emits POSCAR-format content.

  `-` under Write means read-only: a CP2K .inp can be converted FROM, not TO.
  To generate CP2K input use `ferro job -s cp2k`, which writes a full run
  setup rather than bare coordinates.

  `FIRST only` means the format holds one structure: a 500-frame trajectory
  written to POSCAR gives you frame 0 and no warning. Use --stride / --number
  to write one file per frame instead."
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

    /// `.xyz` 里装的 extxyz 必须走 extxyz reader，否则胞、能量、力、速度全部静默丢掉。
    #[test]
    fn test_xyz_extension_holding_extxyz() {
        let dir = std::env::temp_dir();
        let read = |name: &str, text: &str| {
            let p = dir.join(name);
            std::fs::write(&p, text).unwrap();
            read_trajectory(&p, &ReadArgs::default()).unwrap()
        };

        // ASE 写法：Lattice + Properties，带力、速度、能量
        let t = read("dispatch_ase.xyz", "1\nLattice=\"5 0 0 0 5 0 0 0 5\" \
            Properties=species:S:1:pos:R:3:forces:R:3:velocities:R:3 energy=-1.5 pbc=\"T T T\"\n\
            Si 0 0 0 0.1 0.2 0.3 0.01 0.02 0.03\n");
        let f = &t.frames[0];
        assert!(f.cell.is_some(), "Lattice 应被读成晶胞");
        assert!(f.forces.is_some() && f.velocities.is_some(), "力与速度不应丢");
        assert_eq!(f.energy, Some(-1.5));

        // GPUMD 的 train.xyz 常写小写键
        let t = read("dispatch_gpumd.xyz", "1\nlattice=\"5 0 0 0 5 0 0 0 5\" \
            properties=species:S:1:pos:R:3:force:R:3 energy=-2.0\nSi 0 0 0 0.1 0.2 0.3\n");
        assert!(t.frames[0].cell.is_some() && t.frames[0].forces.is_some(), "小写键同样是 extxyz");

        // ferro 自己 dataset --type nep 写出的 .xyz：自产自读必须无损
        let p = dir.join("dispatch_self.xyz");
        let mut src = one_atom_traj();
        src.frames[0].forces = Some(vec![Vector3::new(0.1, -0.2, 0.3)]);
        ferro_io::write_extxyz(&src, &p).unwrap();
        let back = read_trajectory(&p, &ReadArgs::default()).unwrap();
        assert!(back.frames[0].cell.is_some() && back.frames[0].forces.is_some(), "自产 extxyz 读回应无损");

        // 纯 XYZ（含 CP2K 轨迹那种带 = 的注释）仍走 xyz reader
        let t = read("dispatch_plain.xyz", "2\ni = 0, time = 0.000, E = -1.0\nO 0 0 0\nH 0 0 1\n");
        assert_eq!(t.frames[0].n_atoms(), 2);
        assert!(t.frames[0].cell.is_none());
    }

    /// `supported_formats()` 的 Write 列写着 `-` 的两个格式，dispatch 必须真的拒绝。
    /// 表与 match 分支是两处手写的事实，漂了就是文档在说谎。
    #[test]
    fn test_read_only_formats_are_rejected_on_write() {
        let traj = one_atom_traj();
        let dir = std::env::temp_dir();
        for ext in ["inp", "restart"] {
            let path = dir.join(format!("ferro_dispatch_test.{ext}"));
            let err = write_trajectory(&traj, &path, Some(LammpsUnits::Real)).unwrap_err();
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
            ("t.vasp", false), ("t.pos", false), ("conf_0000.vasp", false),
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
    /// `.vasp` / `.pos` 与 POSCAR 前缀走同一对 reader/writer
    #[test]
    fn test_vasp_extensions_round_trip_like_poscar() {
        let traj = one_atom_traj();
        let dir = std::env::temp_dir();
        for name in ["ferro_vasp_test.vasp", "ferro_vasp_test.pos"] {
            let path = dir.join(name);
            write_trajectory(&traj, &path, Some(LammpsUnits::Real)).unwrap();
            let back = read_trajectory(&path, &ReadArgs::default()).unwrap();
            assert_eq!(back.n_frames(), 1, "{name}");
            assert_eq!(back.frames[0].atoms[0].element, "Si", "{name}");
            std::fs::remove_file(&path).ok();
        }
    }

    #[test]
    fn test_unknown_extension_defers_to_the_writer_error() {
        assert!(holds_multiple_frames(Path::new("t.nosuchfmt")));
    }

    #[test]
    fn test_unknown_extension_is_rejected_both_ways() {
        let traj = one_atom_traj();
        let path = std::env::temp_dir().join("ferro_dispatch_test.nosuchfmt");
        assert!(write_trajectory(&traj, &path, Some(LammpsUnits::Real)).is_err());
        assert!(read_trajectory(&path, &ReadArgs::default()).is_err());
    }

    /// LAMMPS data 没给 --atom-style 就报错（三种扩展名都算），给了才读；
    /// 别的格式不受影响
    #[test]
    fn test_lammps_data_needs_an_explicit_atom_style() {
        let traj = one_atom_traj();
        for name in ["ferro_style_test.data", "ferro_style_test.lmp", "ferro_style_test.lammps"] {
            let path = std::env::temp_dir().join(name);
            write_trajectory(&traj, &path, Some(LammpsUnits::Real)).unwrap();
            let err = read_trajectory(&path, &ReadArgs::default()).unwrap_err().to_string();
            assert!(err.contains("--atom-style"), "{name}：{err}");
            let full = ReadArgs { atom_style: Some(AtomStyle::Full), ..ReadArgs::default() };
            let back = read_trajectory(&path, &full).unwrap();
            assert_eq!(back.frames[0].atoms[0].element, "Si", "{name}");
            std::fs::remove_file(&path).ok();
        }
        assert!(check_atom_style(Path::new("t.lammpstrj"), None).is_ok());
        assert!(check_atom_style(Path::new("t.xyz"), None).is_ok());
    }
}
