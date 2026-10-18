//! File writer modules.

pub mod xyz;
pub mod pdb;
pub mod cif;
pub mod vasp;
pub mod extxyz;
pub mod lammps_data;
pub mod lammps_dump;
pub mod qe;
pub mod cp2k;
pub mod cube;
pub mod table;
pub mod deepmd;

pub use xyz::write_xyz;
pub use pdb::write_pdb;
pub use cif::write_cif;
pub use vasp::{write_poscar, poscar_string};
pub use extxyz::{write_extxyz, write_extxyz_with, StressKey};
pub use lammps_data::write_lammps_data;
pub use lammps_dump::{dump_name, write_lammps_dump};
pub use qe::qe_structure_cards;
pub use cp2k::cp2k_structure_sections;
pub use cube::write_cube;
pub use table::{write_table, TableFormat};
pub use deepmd::{write_deepmd_npy, write_deepmd_npy_bounds, write_deepmd_npy_sets};

use ferro_core::{Frame, Trajectory};
use std::path::{Path, PathBuf};

/// The one frame a single-structure format (POSCAR, CIF, LAMMPS data) holds.
///
/// Several frames are an error rather than "write the first": the rest would vanish
/// with a success status. Callers wanting one file per frame split first, naming the
/// files with [`indexed_path`].
pub(crate) fn single_frame<'a>(trajectory: &'a Trajectory, format: &str) -> anyhow::Result<&'a Frame> {
    match trajectory.n_frames() {
        0 => anyhow::bail!("trajectory is empty"),
        1 => Ok(&trajectory.frames[0]),
        n => anyhow::bail!(
            "{format} holds one structure but the trajectory has {n} frames; \
             write one file per frame or select a frame first"
        ),
    }
}

/// Zero-padded width for frame indices, at least 4 so `ls` sorts the products in
/// frame order.
pub fn index_width(max_index: usize) -> usize {
    max_index.to_string().len().max(4)
}

/// Inserts `_<frame index>` before the extension: `out.vasp` -> `out_0010.vasp`,
/// and `POSCAR` -> `POSCAR_0010`.
///
/// The number goes before the extension so the product keeps the extension its
/// format is detected by; extension-less VASP names keep their POSCAR prefix, which
/// is what detects them.
pub fn indexed_path(base: &Path, frame_idx: usize, width: usize) -> PathBuf {
    let stem = base.file_stem().and_then(|s| s.to_str()).unwrap_or("frame");
    let mut name = format!("{stem}_{frame_idx:0width$}");
    if let Some(ext) = base.extension().and_then(|e| e.to_str()) {
        name.push('.');
        name.push_str(ext);
    }
    base.with_file_name(name)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    /// POSCAR, CIF and LAMMPS data hold one structure: several frames are refused
    /// (before any file is created) instead of silently writing the first one.
    #[test]
    fn single_structure_writers_refuse_several_frames() {
        let cell = Cell::from_lengths_angles(3.0, 3.0, 3.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::zeros()));
        let mut traj = Trajectory::from_frame(frame.clone());
        traj.add_frame(frame);

        let dir = std::env::temp_dir().join(format!("ferro_single_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        type Writer = fn(&Trajectory, &std::path::Path) -> anyhow::Result<()>;
        let writers: [(&str, Writer); 3] = [
            ("x.vasp", write_poscar),
            ("x.cif", write_cif),
            ("x.data", write_lammps_data),
        ];
        for (name, write) in writers {
            let path = dir.join(name);
            let err = write(&traj, &path).unwrap_err().to_string();
            assert!(err.contains("2 frames"), "{name}: {err}");
            assert!(!path.exists(), "{name}: a refused write must not leave a file");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn indexed_path_keeps_the_name_its_format_is_detected_by() {
        // Extension formats: the index goes before the extension, which must survive
        assert_eq!(indexed_path(Path::new("out.lmp"), 10, 4), PathBuf::from("out_0010.lmp"));
        // POSCAR has no extension and is detected by its prefix: the index is appended
        assert_eq!(indexed_path(Path::new("POSCAR"), 2, 4), PathBuf::from("POSCAR_0002"));
        // The directory part is kept as written
        assert_eq!(
            indexed_path(Path::new("runs/700K/conf.lmp"), 7, 4),
            PathBuf::from("runs/700K/conf_0007.lmp")
        );
    }

    /// Indices are padded to at least 4 digits, or `ls` puts frame 10 before frame 2
    #[test]
    fn index_width_pads_for_ls_ordering() {
        assert_eq!(index_width(9), 4);
        assert_eq!(index_width(9999), 4);
        assert_eq!(index_width(10000), 5);
        let w = index_width(12345);
        assert_eq!(indexed_path(Path::new("a.xyz"), 7, w), PathBuf::from("a_00007.xyz"));
    }
}
