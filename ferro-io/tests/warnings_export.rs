//! Reader warnings reach a library caller through `ferro_io::read` (`ReadReport::warnings`).
//! Scan B-13 exported the four `_with_warnings` readers for this; P3 took the per-format
//! readers off the crate root again, so the same checks now go through `read`.

use ferro_io::{read, write_pdb, AtomStyle, Level, ReadOptions};
use ferro_core::{Atom, Frame, Trajectory};
use nalgebra::Vector3;
use std::path::PathBuf;

fn tmp(name: &str, text: &str) -> PathBuf {
    let p = std::env::temp_dir().join(name);
    std::fs::write(&p, text).unwrap();
    p
}

#[test]
fn xyz_truncated_last_frame_is_reported() {
    let p = tmp("ferro_export_trunc.xyz", "2\nc\nO 0 0 0\nH 1 0 0\n2\nc\nO 0 0 0\n");
    let (traj, r) = read(&p, &ReadOptions::default()).unwrap();
    let w = r.warnings;
    assert_eq!(traj.n_frames(), 1);
    assert!(w.iter().any(|m| m.level == Level::Warning && m.line == Some(5)
        && m.message.contains("1 of 2 atoms")), "{w:?}");
}

#[test]
fn extxyz_truncated_last_frame_is_reported() {
    let block = "2\nProperties=species:S:1:pos:R:3\nO 0 0 0\nH 1 0 0\n";
    let p = tmp("ferro_export_trunc.extxyz",
                &format!("{block}2\nProperties=species:S:1:pos:R:3\nO 0 0 0\n"));
    let (traj, r) = read(&p, &ReadOptions::default()).unwrap();
    let w = r.warnings;
    assert_eq!(traj.n_frames(), 1);
    assert!(w.iter().any(|m| m.level == Level::Warning && m.line == Some(5)
        && m.message.contains("1 of 2 atoms")), "{w:?}");
}

#[test]
fn pdb_last_block_without_end_is_reported() {
    let mut f = Frame::new();
    f.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.0)));
    let p = std::env::temp_dir().join("ferro_export_noend.pdb");
    write_pdb(&Trajectory::from_frame(f), &p).unwrap();
    let one = std::fs::read_to_string(&p).unwrap();
    let p = tmp("ferro_export_noend.pdb", &format!("{one}{}", one.trim_end_matches("END\n")));
    let (traj, r) = read(&p, &ReadOptions::default()).unwrap();
    let w = r.warnings;
    assert_eq!(traj.n_frames(), 2);
    assert!(w.iter().any(|m| m.level == Level::Warning && m.message.contains("no END")), "{w:?}");
}

#[test]
fn lammps_data_masses_comment_is_reported() {
    let text = "test\n\n1 atoms\n1 atom types\n\n0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\n\
                Masses\n\n1 15.999 # Ow\n\nAtoms # atomic\n\n1 1 0 0 0\n";
    let p = tmp("ferro_export_masses.data", text);
    let opts = ReadOptions { atom_style: Some(AtomStyle::Atomic), ..Default::default() };
    let (traj, r) = read(&p, &opts).unwrap();
    let w = r.warnings;
    assert_eq!(traj.frames[0].atom(0).element, "O");
    assert!(w.iter().any(|m| m.line == Some(12) && m.message.contains("\"Ow\"")), "{w:?}");
}
