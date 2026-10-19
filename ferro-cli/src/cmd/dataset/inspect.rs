//! Diagnostic products for one collected AIMD run.
//!
//! `ferro dataset collect --type inspect` writes three files into a
//! `ferro_inspect/` folder **inside** the AIMD directory, so each run keeps its
//! own diagnostics next to the output that produced them — CP2K logs are very
//! often all called `total.out` and are told apart only by their directory.
//!
//! The three are a trajectory to look at, the last structure to carry on from,
//! and one table of the per-frame scalars that say whether the run was healthy.
//! Nothing here is part of a dataset: no npy is written on this path.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use ferro_core::table::Table;
use ferro_analysis::ml::diagnostics::frame_table;
use ferro_core::Trajectory;
use ferro_io::{write_lammps_data, write_lammps_dump, write_table, AimdFormat, LammpsUnits, TableFormat};

/// The folder the three products go into, inside the AIMD directory.
///
/// Prefixed because it is created in the user's own data directory rather than
/// under `-o`: a bare `inspect/` gives nobody a way to tell, three months later,
/// what put it there or whether it is safe to delete.
pub const DIR_NAME: &str = "ferro_inspect";

/// Writes the trajectory, the last structure and the per-frame table.
///
/// `sources` names the file each frame came from, one entry per frame; `collect`
/// stitches restart segments, and that column plus the step numbers are what make
/// a seam — or a doubled run — visible. `format` is the reader the frames came
/// through.
pub fn write_all(
    traj: &Trajectory,
    format: AimdFormat,
    dir: &Path,
    name: &str,
    sources: &[String],
    overwrite: bool,
) -> Result<Vec<PathBuf>> {
    if !overwrite && dir.exists() && std::fs::read_dir(dir)?.next().is_some() {
        anyhow::bail!("{} exists and is not empty (pass --overwrite)", dir.display());
    }
    std::fs::create_dir_all(dir)
        .with_context(|| format!("cannot create {}", dir.display()))?;

    let trj = dir.join(format!("{name}.lammpstrj"));
    write_lammps_dump(traj, &trj, Some(LammpsUnits::Metal), false)?;

    // 末帧而不是第 0 帧：第 0 帧是喂进 CP2K 的初始构型，用户手上本来就有；
    // 末帧才是这次 AIMD 产出的新东西，接着跑或做下一步 DFT 用的是它
    let data = dir.join(format!("{name}.data"));
    let last = Trajectory::from_frame(
        traj.last().context("the trajectory has no frame")?.clone(),
    );
    write_lammps_data(&last, &data, false)?;

    let info = dir.join(format!("{name}_info.csv"));
    write_table(&info_table(traj, format, sources), &info, TableFormat::Csv)?;

    Ok(vec![trj, data, info])
}

/// The per-frame table, plus the note a vasprun.xml run needs.
///
/// vasprun.xml has no temperature field, so the column is derived; without saying
/// so, the same run read from another file gives numbers nobody can explain.
/// Decided by the reader the frames came through, not by a metadata string, which
/// no reader ever set to "vasprun.xml" (review B4).
fn info_table(traj: &Trajectory, format: AimdFormat, sources: &[String]) -> Table {
    let (mut t, warnings) = frame_table(traj, sources);
    for w in &warnings {
        eprintln!("WARNING: {w}");
    }
    if format == AimdFormat::VaspXml && traj.frames.iter().any(|f| f.temperature.is_some()) {
        t.meta_line(
            "note: temperature is derived from the ionic kinetic energy \
             (T = 2*E_kin/(3N*k_B)); vasprun.xml does not print it"
                .to_string(),
        );
    }
    t
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame};
    use nalgebra::Vector3;

    fn traj_of(temps: &[Option<f64>], steps: &[Option<i64>]) -> Trajectory {
        let mut t = Trajectory::new();
        for (i, (temp, step)) in temps.iter().zip(steps).enumerate() {
            let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
            let mut f = Frame::with_cell(cell, [true; 3]);
            f.add_atom(Atom::new("Si", Vector3::new(0.0, 0.0, 0.0)));
            f.add_atom(Atom::new("O", Vector3::new(1.0, 1.0, 1.0)));
            f.energy = Some(-10.0 - i as f64);
            f.temperature = *temp;
            f.step = *step;
            t.add_frame(f);
        }
        t
    }

    /// The derived-temperature note is shown for vasprun.xml and only for it.
    #[test]
    fn the_vasprun_temperature_note_follows_the_reader() {
        let traj = traj_of(&[Some(300.0), Some(310.0)], &[None, None]);
        let note = |f| info_table(&traj, f, &["a".into(), "a".into()]).meta.join("\n");
        assert!(note(AimdFormat::VaspXml).contains("derived from the ionic kinetic energy"));
        assert!(!note(AimdFormat::VaspOutcar).contains("derived"));
    }
}
