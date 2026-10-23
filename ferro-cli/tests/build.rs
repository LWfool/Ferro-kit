//! End-to-end cover for `ferro build`: the products are read back with `ferro_io::read`.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use ferro_core::Trajectory;

fn ferro(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_ferro")).args(args).output().unwrap()
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("ferro_build_{name}"));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn read(path: &Path) -> Trajectory {
    ferro_io::read(path, &ferro_io::ReadOptions::default()).unwrap().0
}

fn ok(out: &Output) {
    assert_eq!(out.status.code(), Some(0), "{}", String::from_utf8_lossy(&out.stderr));
}

/// Three atoms in a triclinic box, two frames.
const TRICLINIC: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../tests/triclinic_2frames.lammpstrj");
const OUTCAR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../tests/vasp651_sp/frame_0011/OUTCAR");

#[test]
fn supercell_with_given_dims_replicates_every_frame() {
    let d = scratch("sc_dims");
    let o = d.to_str().unwrap();
    ok(&ferro(&["build", "supercell", "-i", TRICLINIC, "--dims", "2", "1", "3",
                "--format", "extxyz", "-o", o, "-s", "t"]));
    let src = read(Path::new(TRICLINIC));
    let made = read(&d.join("supercell_triclinic_2frames_t.extxyz"));
    assert_eq!(made.n_frames(), 2);
    for (a, b) in src.frames.iter().zip(&made.frames) {
        assert_eq!(b.n_atoms(), a.n_atoms() * 6);
        let (ca, cb) = (a.cell.as_ref().unwrap().matrix, b.cell.as_ref().unwrap().matrix);
        for (row, n) in [(0, 2.0), (1, 1.0), (2, 3.0)] {
            assert!((cb.row(row) - ca.row(row) * n).norm() < 1e-9, "row {row}");
        }
    }
}

/// `--min-length` is reached by every interplanar spacing; the triclinic box is
/// where that differs from the vector length.
#[test]
fn supercell_auto_reaches_min_length_and_min_atoms() {
    let d = scratch("sc_auto");
    ok(&ferro(&["build", "supercell", "-i", TRICLINIC, "--min-length", "25", "--min-atoms", "100",
                "--format", "extxyz", "-o", d.to_str().unwrap()]));
    let made = read(&d.join("supercell_triclinic_2frames.extxyz"));
    let f = &made.frames[0];
    assert!(f.n_atoms() >= 100);
    for s in f.cell.as_ref().unwrap().interplanar_spacings().unwrap() {
        assert!(s >= 25.0, "spacing {s}");
    }
}

#[test]
fn supercell_dry_run_writes_nothing() {
    let d = scratch("sc_dry");
    let out = ferro(&["build", "supercell", "-i", TRICLINIC, "--dims", "2", "2", "2",
                      "--dry-run", "-o", d.to_str().unwrap()]);
    ok(&out);
    assert!(String::from_utf8_lossy(&out.stdout).contains("2x2x2 (given): 3 -> 24 atoms"));
    assert_eq!(std::fs::read_dir(&d).unwrap().count(), 0);
}

/// An OUTCAR cannot be written back, and that is known before it is read.
#[test]
fn a_read_only_input_format_needs_format() {
    let d = scratch("sc_readonly").join("new");
    let out = ferro(&["build", "supercell", "-i", OUTCAR, "--dims", "2", "1", "1",
                      "-o", d.to_str().unwrap(), "--mkdir"]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("needs --format EXT"));
    assert!(!d.exists(), "-o was created for a run that could not write");
}
