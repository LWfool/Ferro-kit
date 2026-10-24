//! The one table of file formats: how each is recognised, whether it can be
//! written, and whether it holds more than one frame.
//!
//! The CLI and the Python bindings both read and write through [`read`] and
//! [`write`], so a format added here is reachable from both, and the help-page
//! format matrix is generated from [`Format::ALL`] rather than written by hand.

use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use ferro_core::{ChargeGrid, Frame, Trajectory};
use crate::Notice;

use crate::readers::{
    apply_type_map, is_extxyz, read_aimd_as, read_chgcar, read_cif, read_cp2k_inp, read_cp2k_restart, read_cube_as_chg,
    read_deepmd_npy, read_extxyz, read_lammps_data, read_lammps_dump,
    read_pdb, read_poscar, read_qe_input, read_xyz, sniff, AimdFormat, AimdStats,
    AtomStyle, LammpsUnits,
};
use crate::writers::{
    index_width, indexed_path, write_cif, write_extxyz_with, write_lammps_data,
    write_lammps_dump, StressKey,
    write_pdb, write_poscar, write_xyz,
};

/// A file format ferro reads, and for some of them writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Xyz,
    Extxyz,
    Pdb,
    LammpsDump,
    Cif,
    Poscar,
    LammpsData,
    QeInput,
    Cp2kInput,
    Cp2kRestart,
    /// A program's output log, recognised by its banner rather than its name
    Aimd(AimdFormat),
    /// A DeePMD system directory (`type.raw` + `set.*/`)
    Deepmd,
}

impl Format {
    /// Every format, in the order the help page lists them.
    ///
    /// CP2K's single-point and MD logs share one row: they share a name and a
    /// banner, and only the run type inside tells them apart.
    pub const ALL: [Format; 15] = [
        Format::Xyz, Format::Extxyz, Format::Pdb, Format::LammpsDump,
        Format::Cif, Format::Poscar, Format::LammpsData,
        Format::QeInput, Format::Cp2kInput, Format::Cp2kRestart,
        Format::Aimd(AimdFormat::VaspOutcar), Format::Aimd(AimdFormat::VaspXml),
        Format::Aimd(AimdFormat::Cp2kMd), Format::Aimd(AimdFormat::QeOut), Format::Deepmd,
    ];

    /// Short name, for the format matrix and messages.
    pub fn name(self) -> &'static str {
        match self {
            Format::Xyz => "XYZ",
            Format::Extxyz => "extended XYZ",
            Format::Pdb => "PDB",
            Format::LammpsDump => "LAMMPS dump",
            Format::Cif => "CIF",
            Format::Poscar => "VASP",
            Format::LammpsData => "LAMMPS data",
            Format::QeInput => "QE (pw.x)",
            Format::Cp2kInput => "CP2K input",
            Format::Cp2kRestart => "CP2K restart",
            Format::Aimd(AimdFormat::VaspOutcar) => "VASP OUTCAR",
            Format::Aimd(AimdFormat::VaspXml) => "vasprun.xml",
            Format::Aimd(AimdFormat::Cp2kMd | AimdFormat::Cp2kSp) => "CP2K output",
            Format::Aimd(AimdFormat::QeOut) => "QE output",
            Format::Deepmd => "DeePMD",
        }
    }

    /// What recognises the format, for the format matrix.
    pub fn detected_by(self) -> &'static str {
        match self {
            Format::Xyz => ".xyz",
            Format::Extxyz => ".extxyz  .xyz*",
            Format::Pdb => ".pdb",
            Format::LammpsDump => ".dump  .lammpstrj",
            Format::Cif => ".cif",
            Format::Poscar => ".vasp .poscar .pos POSCAR*",
            Format::LammpsData => ".lammps  .data  .lmp",
            Format::QeInput => ".in  .qe",
            Format::Cp2kInput => ".inp",
            Format::Cp2kRestart => ".restart",
            Format::Aimd(AimdFormat::VaspOutcar) => "banner (OUTCAR)",
            Format::Aimd(AimdFormat::VaspXml) => "banner (.xml)",
            Format::Aimd(AimdFormat::Cp2kMd | AimdFormat::Cp2kSp) => "banner (.out)",
            Format::Aimd(AimdFormat::QeOut) => "banner (.out)",
            Format::Deepmd => "directory with type.raw",
        }
    }

    /// Whether ferro writes this format. The inputs and the program outputs are
    /// read-only: a full CP2K / QE input is a job generator's work, not a converter's.
    pub fn can_write(self) -> bool {
        matches!(
            self,
            Format::Xyz | Format::Extxyz | Format::Pdb | Format::LammpsDump
                | Format::Cif | Format::Poscar | Format::LammpsData
        )
    }

    /// Whether one file holds a whole trajectory. CIF can technically carry several
    /// `data_` blocks, but VESTA / Materials Studio show only the first, so it holds
    /// one structure like POSCAR (user decision, 2026-10-18).
    pub fn multi_frame(self) -> bool {
        matches!(self, Format::Xyz | Format::Extxyz | Format::Pdb | Format::LammpsDump)
    }

    /// The format a file NAME declares, without touching the disk; `None` when the
    /// name says nothing (an OUTCAR, a `.out`, an unknown extension).
    ///
    /// A POSCAR / CONTCAR prefix comes first, for VASP files with no extension;
    /// extensions are compared case-insensitively.
    pub fn from_name(path: &Path) -> Option<Format> {
        let name = path.file_name().and_then(|n| n.to_str()).unwrap_or_default();
        let upper = name.to_uppercase();
        if upper.starts_with("POSCAR") || upper.starts_with("CONTCAR") {
            return Some(Format::Poscar);
        }
        let ext = path.extension().and_then(|e| e.to_str())?.to_lowercase();
        Some(match ext.as_str() {
            "xyz" => Format::Xyz,
            "extxyz" => Format::Extxyz,
            "pdb" => Format::Pdb,
            "dump" | "lammpstrj" => Format::LammpsDump,
            "cif" => Format::Cif,
            "vasp" | "poscar" | "pos" => Format::Poscar,
            "lammps" | "data" | "lmp" => Format::LammpsData,
            "in" | "qe" => Format::QeInput,
            "inp" => Format::Cp2kInput,
            "restart" => Format::Cp2kRestart,
            _ => return None,
        })
    }

    /// The format to read `path` as: the name first, the content as fallback.
    ///
    /// Two looks inside: a `.xyz` whose comment line declares `Lattice=` or
    /// `Properties=` is extended XYZ (ASE and GPUMD write it that way), and a file
    /// whose name says nothing is identified by its program banner, because VASP
    /// writes `OUTCAR` with no extension and `.out` belongs to no one program.
    pub fn detect(path: &Path) -> Result<Format> {
        if path.is_dir() {
            if is_deepmd_system(path) {
                return Ok(Format::Deepmd);
            }
            bail!("{} is a directory but not a DeePMD system (no type.raw)", path.display());
        }
        match Format::from_name(path) {
            Some(Format::Xyz) if is_extxyz(path)? => Ok(Format::Extxyz),
            Some(f) => Ok(f),
            None => {
                let sniffed = sniff(path).map(Format::Aimd);
                // Program outputs carry no extension or a generic one: a failure there
                // is the banner's to explain. Any other extension is simply unknown
                let ext = path.extension().map(|e| e.to_string_lossy().to_lowercase());
                match ext.as_deref() {
                    None | Some("out" | "xml" | "log") => sniffed,
                    Some(ext) => sniffed.with_context(|| format!("Unsupported input format: .{ext}")),
                }
            }
        }
    }

    /// The format to write `path` as, from its name alone, so a caller can refuse a
    /// bad output path before reading a multi-GB input.
    pub fn for_output(path: &Path) -> Result<Format> {
        let ext = || path.extension().map(|e| e.to_string_lossy().to_string()).unwrap_or_default();
        match Format::from_name(path) {
            Some(f) if f.can_write() => Ok(f),
            Some(Format::QeInput) => bail!(
                "Unsupported output format: .{} (QE input is read-only here; \
                 generate a pw.x input with `ferro-py job -s qe`)", ext()
            ),
            Some(_) => bail!("Unsupported output format: .{} (read-only)", ext()),
            None if path.extension().is_some() => bail!("Unsupported output format: .{}", ext()),
            None => bail!("Cannot determine format (no extension): {}", path.display()),
        }
    }

    /// The format a user names to write, as an extension (`cif`, `.vasp`, `EXTXYZ`), with
    /// the checks of [`Format::for_output`]. For the CLI's `--format` and Python's
    /// `write(format=)`: the product path is not asked, because [`Format::from_name`] takes
    /// a `POSCAR_x.cif` for a POSCAR.
    pub fn for_extension(ext: &str) -> Result<Format> {
        Format::for_output(Path::new(&format!("x.{}", ext.trim_start_matches('.').to_lowercase())))
    }
}

/// Reads a charge-density grid: a `.cube` file (Gaussian / QE pp.x), or a VASP CHGCAR
/// for any other name (VASP writes `CHGCAR`, `AECCAR0`, `CHG` with no extension).
///
/// The grid formats are not rows of [`Format`]: they carry one structure plus a
/// volumetric field and feed `bader` only, never [`read`] / `convert`. The extension
/// test lives here all the same, so format recognition stays in this one file.
pub fn read_charge_grid(path: &Path) -> Result<(Frame, ChargeGrid)> {
    let is_cube = path.extension().and_then(|e| e.to_str()).is_some_and(|e| e.eq_ignore_ascii_case("cube"));
    if is_cube { read_cube_as_chg(path) } else { read_chgcar(path) }
}

/// Whether `dir` looks like a DeePMD system: a directory holding `type.raw`.
pub fn is_deepmd_system(dir: &Path) -> bool {
    dir.join("type.raw").is_file()
}

/// Reading choices a file cannot carry.
#[derive(Debug, Clone, Default)]
pub struct ReadOptions {
    /// LAMMPS units of a dump's velocity / force columns; never guessed
    pub units: Option<LammpsUnits>,
    /// LAMMPS data atom style; never guessed
    pub atom_style: Option<AtomStyle>,
    /// LAMMPS dump without an element column: type number -> element symbol
    pub type_map: Option<Vec<(u32, String)>>,
}

/// What a read noticed beyond the frames themselves.
#[derive(Debug, Clone, Default)]
pub struct ReadReport {
    /// Messages for the caller to show; the readers never print
    pub warnings: Vec<Notice>,
    /// Drop counts and version of a program output, when one was read
    pub aimd: Option<AimdStats>,
}

/// Reads any supported file or DeePMD directory, recognising its format.
///
/// Program outputs drop frames the way `dataset collect` does (unconverged SCF,
/// an incomplete last step); the counts are in the report, and a file with no
/// usable frame is an error.
pub fn read(path: &Path, opts: &ReadOptions) -> Result<(Trajectory, ReadReport)> {
    read_as(path, Format::detect(path)?, opts)
}

/// Reads `path` as the given format.
pub fn read_as(path: &Path, format: Format, opts: &ReadOptions) -> Result<(Trajectory, ReadReport)> {
    if opts.type_map.is_some() && format != Format::LammpsDump {
        bail!("a type map applies to LAMMPS dump input only, not {}", path.display());
    }
    let mut report = ReadReport::default();
    let mut warned = |(traj, w): (Trajectory, Vec<Notice>)| {
        report.warnings.extend(w);
        traj
    };
    let traj = match format {
        Format::Xyz => warned(read_xyz(path)?),
        Format::Extxyz => warned(read_extxyz(path)?),
        Format::Pdb => warned(read_pdb(path)?),
        Format::Cif => read_cif(path)?,
        Format::Poscar => warned(read_poscar(path)?),
        Format::LammpsData => {
            let Some(style) = opts.atom_style else {
                bail!(
                    "{} is a LAMMPS data file: its atom style (atomic, charge or full) must be \
                     given; ferro never guesses it from the Atoms comment or the column count",
                    path.display()
                );
            };
            warned(read_lammps_data(path, style)?)
        }
        Format::LammpsDump => {
            let mut traj = warned(read_lammps_dump(path, opts.units)?);
            if let Some(map) = &opts.type_map {
                apply_type_map(&mut traj, map)?;
            }
            traj
        }
        Format::QeInput => read_qe_input(path)?,
        Format::Cp2kInput => read_cp2k_inp(path)?,
        Format::Cp2kRestart => read_cp2k_restart(path)?,
        Format::Deepmd => warned(read_deepmd_npy(path)?),
        Format::Aimd(fmt) => {
            let (traj, stats) = read_aimd_as(path, fmt)?;
            report.warnings.extend(aimd_warnings(&stats));
            report.aimd = Some(stats);
            traj
        }
    };
    let mut traj = traj;
    // Every reader then stores one site table for a run of identical frames (review B15)
    traj.share_sites();
    Ok((traj, report))
}

/// The drop counts and version note of a program output, as messages.
fn aimd_warnings(stats: &AimdStats) -> Vec<Notice> {
    let mut out = Vec::new();
    if stats.n_dropped() > 0 {
        out.push(Notice::warning(None, format!(
            "{} of {} step(s) dropped ({} SCF not converged ({}), {} incomplete, \
             {} composition, {} no cell line)",
            stats.n_dropped(), stats.n_steps, stats.n_scf_failed, stats.format.convergence_rule(),
            stats.n_incomplete, stats.n_bad_composition, stats.n_no_cell
        )));
    }
    if let Some(note) = &stats.version_note {
        out.push(Notice::note(None, note.clone()));
    }
    out
}

/// Writing choices a trajectory cannot carry.
#[derive(Debug, Clone, Copy, Default)]
pub struct WriteOptions {
    /// LAMMPS units of a dump's velocity / force columns; never guessed
    pub units: Option<LammpsUnits>,
    /// Leave out a charge / magnetic-moment column that only some atoms of a frame
    /// carry (LAMMPS data: write every charge as 0), instead of refusing the frame
    pub drop_partial: bool,
}

/// Writes `traj` to `path` in the format its name declares, returning the files
/// written. A format that holds one structure gets one file per frame,
/// `conf.cif` -> `conf_0000.cif`, `conf_0001.cif`, ...
pub fn write(traj: &Trajectory, path: &Path, opts: &WriteOptions) -> Result<Vec<PathBuf>> {
    let numbers: Vec<usize> = (0..traj.n_frames()).collect();
    write_numbered(traj, path, Format::for_output(path)?, opts, &numbers)
}

/// [`write()`] in a given `format`, naming the per-frame files by `numbers` instead of
/// by position — `ferro convert` numbers them by the frame they came from in the input.
///
/// The format is not read off `path`: `convert` names its product after the input, and
/// [`Format::from_name`] takes a `POSCAR_x.extxyz` for a POSCAR whatever `--format` said.
pub fn write_numbered(
    traj: &Trajectory,
    path: &Path,
    format: Format,
    opts: &WriteOptions,
    numbers: &[usize],
) -> Result<Vec<PathBuf>> {
    if format.multi_frame() || traj.n_frames() <= 1 {
        write_one(traj, path, format, opts)?;
        return Ok(vec![path.to_path_buf()]);
    }
    assert_eq!(numbers.len(), traj.n_frames(), "one number per frame");
    let width = index_width(numbers.iter().copied().max().unwrap_or(0));
    let mut written = Vec::with_capacity(numbers.len());
    for (&n, frame) in numbers.iter().zip(&traj.frames) {
        let one = Trajectory { frames: vec![frame.clone()], metadata: traj.metadata.clone() };
        let file = indexed_path(path, n, width);
        write_one(&one, &file, format, opts)
            .with_context(|| format!("writing frame {n} to {}", file.display()))?;
        written.push(file);
    }
    Ok(written)
}

fn write_one(traj: &Trajectory, path: &Path, format: Format, opts: &WriteOptions) -> Result<()> {
    match format {
        Format::Xyz => write_xyz(traj, path),
        Format::Extxyz => write_extxyz_with(traj, path, StressKey::Stress, opts.drop_partial),
        Format::Pdb => write_pdb(traj, path),
        Format::LammpsDump => write_lammps_dump(traj, path, opts.units, opts.drop_partial),
        Format::Cif => write_cif(traj, path),
        Format::Poscar => write_poscar(traj, path),
        Format::LammpsData => write_lammps_data(traj, path, opts.drop_partial),
        other => bail!("Unsupported output format: {}", other.name()),
    }
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
        frame.add_atom(Atom::new("Si", Vector3::new(0.0, 0.0, 0.0)));
        Trajectory { frames: vec![frame], metadata: Default::default() }
    }

    fn tmp_dir(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("ferro_format_{name}"));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn read_plain(path: &Path) -> Result<Trajectory> {
        Ok(read(path, &ReadOptions::default())?.0)
    }

    /// An extxyz inside a `.xyz` must go to the extxyz reader, or the cell, energy,
    /// forces and velocities are all dropped without a word.
    #[test]
    fn xyz_extension_holding_extxyz() {
        let dir = tmp_dir("xyz");
        let read = |name: &str, text: &str| {
            let p = dir.join(name);
            std::fs::write(&p, text).unwrap();
            read_plain(&p).unwrap()
        };

        // ASE: Lattice + Properties, with forces, velocities and energy
        let t = read("ase.xyz", "1\nLattice=\"5 0 0 0 5 0 0 0 5\" \
            Properties=species:S:1:pos:R:3:forces:R:3:velocities:R:3 energy=-1.5 pbc=\"T T T\"\n\
            Si 0 0 0 0.1 0.2 0.3 0.01 0.02 0.03\n");
        let f = &t.frames[0];
        assert!(f.cell.is_some(), "Lattice should become the cell");
        assert!(f.forces().is_some() && f.velocities().is_some(), "forces and velocities must survive");
        assert_eq!(f.energy, Some(-1.5));

        // GPUMD's train.xyz often uses lower-case keys
        let t = read("gpumd.xyz", "1\nlattice=\"5 0 0 0 5 0 0 0 5\" \
            properties=species:S:1:pos:R:3:force:R:3 energy=-2.0\nSi 0 0 0 0.1 0.2 0.3\n");
        assert!(t.frames[0].cell.is_some() && t.frames[0].forces().is_some(), "lower-case keys are extxyz too");

        // The .xyz ferro's own `dataset --type nep` writes reads back losslessly
        let p = dir.join("self.xyz");
        let mut src = one_atom_traj();
        src.frames[0].set_forces(Some(vec![Vector3::new(0.1, -0.2, 0.3)])).unwrap();
        crate::writers::write_extxyz(&src, &p).unwrap();
        let back = read_plain(&p).unwrap();
        assert!(back.frames[0].cell.is_some() && back.frames[0].forces().is_some(), "own extxyz must read back");

        // Plain XYZ, including a CP2K trajectory comment with `=` in it, stays XYZ
        let t = read("plain.xyz", "2\ni = 0, time = 0.000, E = -1.0\nO 0 0 0\nH 0 0 1\n");
        assert_eq!(t.frames[0].n_atoms(), 2);
        assert!(t.frames[0].cell.is_none());
    }

    /// The read-only formats are refused on write, and leave no empty file behind.
    #[test]
    fn read_only_formats_are_rejected_on_write() {
        let traj = one_atom_traj();
        let dir = tmp_dir("ro");
        for ext in ["inp", "restart", "in", "qe", "out", "xml"] {
            let path = dir.join(format!("t.{ext}"));
            let err = write(&traj, &path, &WriteOptions { units: Some(LammpsUnits::Real), ..Default::default() }).unwrap_err();
            assert!(err.to_string().contains("Unsupported output format"), ".{ext}: {err}");
            assert!(!path.exists(), ".{ext} must not leave a file");
        }
    }

    #[test]
    fn multi_frame_follows_the_name() {
        for (name, multi) in [
            ("t.xyz", true), ("t.extxyz", true), ("t.pdb", true), ("t.cif", false),
            ("t.dump", true), ("t.lammpstrj", true),
            ("POSCAR", false), ("CONTCAR", false), ("poscar", false),
            ("t.vasp", false), ("t.poscar", false), ("t.pos", false), ("conf_0000.vasp", false),
            ("t.lmp", false), ("t.data", false), ("t.lammps", false),
            ("T.XYZ", true), ("T.CIF", false),
        ] {
            assert_eq!(Format::for_output(Path::new(name)).unwrap().multi_frame(), multi, "{name}");
        }
        // A prefix match, not a substring match
        assert!(Format::for_output(Path::new("run1_POSCAR")).is_err());
    }

    /// Several frames to a one-structure format: one file per frame, named by position
    /// (`write`) or by the numbers given (`write_numbered`).
    #[test]
    fn single_structure_formats_split_into_numbered_files() {
        let mut traj = one_atom_traj();
        traj.frames.push(traj.frames[0].clone());
        traj.frames.push(traj.frames[0].clone());
        let dir = tmp_dir("split");
        let made = write(&traj, &dir.join("conf.cif"), &WriteOptions::default()).unwrap();
        assert_eq!(made, ["conf_0000.cif", "conf_0001.cif", "conf_0002.cif"].map(|n| dir.join(n)));
        let made = write_numbered(&traj, &dir.join("POSCAR"), Format::Poscar, &WriteOptions::default(), &[0, 50, 12345]).unwrap();
        assert_eq!(made, ["POSCAR_00000", "POSCAR_00050", "POSCAR_12345"].map(|n| dir.join(n)));
        assert!(made.iter().all(|p| p.is_file()));
        let one = write(&traj, &dir.join("all.extxyz"), &WriteOptions::default()).unwrap();
        assert_eq!(read_plain(&one[0]).unwrap().n_frames(), 3);
    }

    /// `.vasp` / `.pos` go through the same reader and writer as a POSCAR name,
    /// and an upper-case extension is the same format.
    #[test]
    fn vasp_extensions_and_case_round_trip() {
        let traj = one_atom_traj();
        let dir = tmp_dir("vasp");
        for name in ["t.vasp", "t.poscar", "t.pos", "T.VASP", "t.XYZ"] {
            let path = dir.join(name);
            write(&traj, &path, &WriteOptions { units: Some(LammpsUnits::Real), ..Default::default() }).unwrap();
            let back = read_plain(&path).unwrap();
            assert_eq!(back.n_frames(), 1, "{name}");
            assert_eq!(back.frames[0].element(0), "Si", "{name}");
        }
    }

    #[test]
    fn unknown_extension_is_rejected_both_ways() {
        let traj = one_atom_traj();
        let path = tmp_dir("unknown").join("t.nosuchfmt");
        assert!(write(&traj, &path, &WriteOptions::default()).unwrap_err().to_string().contains("Unsupported output format"));
        std::fs::write(&path, "hello\n").unwrap();
        let err = format!("{:#}", read_plain(&path).unwrap_err());
        assert!(err.contains("Unsupported input format: .nosuchfmt"), "{err}");
    }

    /// A LAMMPS data file is not read without its atom style, whatever the extension;
    /// a type map on anything but a dump is refused rather than ignored.
    #[test]
    fn lammps_data_needs_an_atom_style() {
        let traj = one_atom_traj();
        let dir = tmp_dir("data");
        for name in ["t.data", "t.lmp", "t.lammps"] {
            let path = dir.join(name);
            write(&traj, &path, &WriteOptions::default()).unwrap();
            let err = read_plain(&path).unwrap_err().to_string();
            assert!(err.contains("atom style"), "{name}: {err}");
            let full = ReadOptions { atom_style: Some(AtomStyle::Full), ..Default::default() };
            assert_eq!(read(&path, &full).unwrap().0.frames[0].element(0), "Si", "{name}");
            let mapped = ReadOptions { type_map: Some(vec![(1, "P".into())]), ..full };
            assert!(read(&path, &mapped).unwrap_err().to_string().contains("LAMMPS dump input only"));
        }
    }

    /// A stray byte in a program log (MPI noise sent through `2>&1`) is an error at
    /// its line, from the banner check and from the reader alike (review C-8; user
    /// decision: strict, not replaced).
    #[test]
    fn a_non_utf8_byte_in_a_log_names_its_line() {
        let mut bytes = std::fs::read("../tests/cp2k_sp_v61.out").unwrap();
        let at = bytes.iter().enumerate().filter(|(_, &c)| c == b'\n').nth(199).unwrap().0 + 1;
        bytes.insert(at, 0xff);
        let path = tmp_dir("utf8").join("bad.out");
        std::fs::write(&path, &bytes).unwrap();
        let err = format!("{:#}", read_plain(&path).unwrap_err());
        assert!(err.contains("line 201, column 1: byte 0xff"), "{err}");
        let err = format!("{:#}", read_as(&path, Format::Aimd(AimdFormat::Cp2kSp), &ReadOptions::default()).unwrap_err());
        assert!(err.contains("line 201, column 1: byte 0xff"), "{err}");
    }

    /// Program outputs and DeePMD systems are found by content and by directory, so
    /// `convert` / `info` / `traj` read them like any structure file (review C-2).
    #[test]
    fn program_outputs_and_deepmd_systems_are_detected() {
        let fx = Path::new("../tests");
        for (name, want) in [
            ("vasp_OUTCAR_2frames", Format::Aimd(AimdFormat::VaspOutcar)),
            ("vasp651_sp/frame_0011/OUTCAR", Format::Aimd(AimdFormat::VaspOutcar)),
            ("vasp651_sp/frame_0011/vasprun.xml", Format::Aimd(AimdFormat::VaspXml)),
            ("cp2k_sp_v61.out", Format::Aimd(AimdFormat::Cp2kSp)),
        ] {
            let path = fx.join(name);
            assert_eq!(Format::detect(&path).unwrap(), want, "{name}");
            let (traj, report) = read(&path, &ReadOptions::default()).unwrap();
            assert!(traj.n_frames() >= 1 && report.aimd.is_some(), "{name}");
        }

        let dir = tmp_dir("deepmd");
        let mut traj = one_atom_traj();
        traj.frames[0].energy = Some(-1.0);
        traj.frames[0].set_forces(Some(vec![Vector3::zeros()])).unwrap();
        crate::write_deepmd_npy(&traj, &dir).unwrap();
        assert_eq!(Format::detect(&dir).unwrap(), Format::Deepmd);
        assert_eq!(read_plain(&dir).unwrap().frames[0].element(0), "Si");
        let err = Format::detect(dir.join("set.000").as_path()).unwrap_err().to_string();
        assert!(err.contains("not a DeePMD system"), "{err}");
    }

    /// The name picks the reader: `.cube` in any case goes to the cube reader, anything
    /// else to CHGCAR. The same CHGCAR text under a `.cube` name fails as a cube (review D-10)
    #[test]
    fn read_charge_grid_dispatches_on_the_extension() {
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("../tests/CHGCAR_2atoms");
        let text = std::fs::read_to_string(&fixture).unwrap();
        let (frame, chg) = read_charge_grid(&fixture).unwrap();
        assert_eq!((frame.n_atoms(), chg.shape), (2, [8, 8, 8]));
        for name in ["ferro_grid_dispatch.cube", "ferro_grid_dispatch.CUBE"] {
            let as_cube = crate::testutil::write_tmp(name, &text);
            assert!(read_charge_grid(&as_cube).is_err(), "{name} must go to the cube reader");
        }
        let as_chgcar = crate::testutil::write_tmp("ferro_grid_dispatch_AECCAR0", &text);
        assert_eq!(read_charge_grid(&as_chgcar).unwrap().1.shape, [8, 8, 8]);
    }
}
