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
                "--format-out", "extxyz", "-o", o, "-s", "t"]));
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

/// The product takes the format the input was read as (`--format-in`, or an extended XYZ
/// inside a `.xyz`), not the one its name says (review 2026-10-30 H6). An extended XYZ
/// `.xyz` used to come out as plain XYZ, its cell lost; `--format-in` was ignored.
#[test]
fn supercell_writes_the_format_the_input_was_read_as() {
    let d = scratch("sc_format_in");
    let ext = "2\nLattice=\"4 0 0 0 4 0 0 0 4\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n\
               Na 0 0 0\nCl 2 2 2\n";
    for name in ["a.xyz", "b.txt", "c.lmp"] {
        std::fs::write(d.join(name), ext).unwrap();
    }
    let o = d.join("out");
    let o = o.to_str().unwrap();
    let run = |name: &str, extra: &[&str]| {
        let input = d.join(name);
        let mut args = vec!["build", "supercell", "-i", input.to_str().unwrap(), "--dims", "2", "1", "1",
                            "-o", o, "--mkdir"];
        args.extend_from_slice(extra);
        ok(&ferro(&args));
    };
    run("a.xyz", &[]);
    run("b.txt", &["--format-in", "extxyz"]);
    run("c.lmp", &["--format-in", "extxyz"]);
    for stem in ["a", "b", "c"] {
        let made = read(&Path::new(o).join(format!("supercell_{stem}.extxyz")));
        let f = &made.frames[0];
        assert_eq!(f.n_atoms(), 4);
        assert!((f.cell.as_ref().expect("cell lost").matrix[(0, 0)] - 8.0).abs() < 1e-9);
    }
}

/// `--min-length` is reached by every interplanar spacing; the triclinic box is
/// where that differs from the vector length.
#[test]
fn supercell_auto_reaches_min_length_and_min_atoms() {
    let d = scratch("sc_auto");
    ok(&ferro(&["build", "supercell", "-i", TRICLINIC, "--min-length", "25", "--min-atoms", "100",
                "--format-out", "extxyz", "-o", d.to_str().unwrap()]));
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

/// Each frame stays under the 1e8 cap (3 atoms x 2e7), the two-frame trajectory does
/// not; it is refused before any supercell is built (review fb7d8c3-1).
#[test]
fn supercell_over_the_trajectory_total_is_refused() {
    let d = scratch("sc_total");
    let out = ferro(&["build", "supercell", "-i", TRICLINIC, "--dims", "200", "200", "500",
                      "--format-out", "extxyz", "-o", d.to_str().unwrap()]);
    assert_ne!(out.status.code(), Some(0));
    assert!(String::from_utf8_lossy(&out.stderr).contains("120000000 atoms in total"),
            "{}", String::from_utf8_lossy(&out.stderr));
    assert_eq!(std::fs::read_dir(&d).unwrap().count(), 0);
}

/// A parameter error, so it fails before the (missing) input is read.
#[test]
fn supercell_min_atoms_over_the_cap_fails_before_reading() {
    let out = ferro(&["build", "supercell", "-i", "no_such_file.cif", "--min-atoms", "200000000"]);
    let err = String::from_utf8_lossy(&out.stderr);
    assert_ne!(out.status.code(), Some(0));
    assert!(err.contains("--min-atoms must be at most") && !err.contains("no_such_file"), "{err}");
}

/// The spacing along the vacuum axis grows by exactly the thickness; the others stay.
#[test]
fn vacuum_grows_the_spacing_along_its_axis() {
    let d = scratch("vac");
    ok(&ferro(&["build", "vacuum", "-i", TRICLINIC, "--axis", "y", "--thickness", "7.5",
                "--format-out", "extxyz", "-o", d.to_str().unwrap()]));
    let src = read(Path::new(TRICLINIC));
    let made = read(&d.join("vacuum_triclinic_2frames.extxyz"));
    assert_eq!(made.n_frames(), 2);
    for (a, b) in src.frames.iter().zip(&made.frames) {
        let sa = a.cell.as_ref().unwrap().interplanar_spacings().unwrap();
        let sb = b.cell.as_ref().unwrap().interplanar_spacings().unwrap();
        assert!((sb[1] - sa[1] - 7.5).abs() < 1e-9 && (sb[0] - sa[0]).abs() < 1e-9 && (sb[2] - sa[2]).abs() < 1e-9,
                "{sa:?} -> {sb:?}");
    }
}

/// A box built twice with one seed is the same file; merging it with itself doubles it.
#[test]
fn box_is_reproducible_and_merge_stacks_two() {
    let d = scratch("box");
    let o = d.to_str().unwrap();
    for s in ["a", "b"] {
        ok(&ferro(&["build", "box", "-c", "H2O:30", "-c", "ZnO:5", "--density", "1.2", "--seed", "7",
                    "-o", o, "-s", s]));
    }
    let (a, b) = (d.join("box_a.vasp"), d.join("box_b.vasp"));
    assert_eq!(std::fs::read(&a).unwrap(), std::fs::read(&b).unwrap());
    assert_eq!(read(&a).frames[0].n_atoms(), 30 * 3 + 5 * 2);

    ok(&ferro(&["build", "merge", "-i", a.to_str().unwrap(), b.to_str().unwrap(), "--gap", "2", "-o", o]));
    let merged = read(&d.join("merge_box_a_box_b.vasp"));
    assert_eq!(merged.frames[0].n_atoms(), 2 * 100);
}

/// merge does not pick a frame of a trajectory on the user's behalf.
#[test]
fn merge_refuses_a_trajectory() {
    let d = scratch("merge_traj");
    let out = ferro(&["build", "merge", "-i", TRICLINIC, TRICLINIC, "--gap", "1", "-o", d.to_str().unwrap()]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("holds 2 frames"));
}

/// An OUTCAR cannot be written back, and that is known before it is read.
#[test]
fn a_read_only_input_format_needs_format_out() {
    let d = scratch("sc_readonly").join("new");
    let out = ferro(&["build", "supercell", "-i", OUTCAR, "--dims", "2", "1", "1",
                      "-o", d.to_str().unwrap(), "--mkdir"]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("needs --format-out EXT"));
    assert!(!d.exists(), "-o was created for a run that could not write");
}

/// The reader options are checked for every input before the first is read, as traj / map /
/// net do: a LAMMPS data file without --atom-style used to be found only when its turn came,
/// after the dump before it had been written (review 2026-10-30 A5-M2).
#[test]
fn reader_options_are_checked_before_any_input_is_read() {
    let d = scratch("read_args_first");
    let data = d.join("b.data");
    std::fs::write(&data, "x\n\n1 atoms\n1 atom types\n0 4 xlo xhi\n0 4 ylo yhi\n0 4 zlo zhi\n\n\
                           Masses\n\n1 22.99\n\nAtoms # atomic\n\n1 1 0 0 0\n").unwrap();
    let o = d.join("out");
    let out = ferro(&["build", "supercell", "-i", TRICLINIC, data.to_str().unwrap(),
                      "--dims", "2", "1", "1", "-o", o.to_str().unwrap(), "--mkdir"]);
    assert_ne!(out.status.code(), Some(0));
    assert!(String::from_utf8_lossy(&out.stderr).contains("--atom-style"));
    assert!(!o.exists(), "nothing may be written");
    let out = ferro(&["build", "merge", "-i", TRICLINIC, data.to_str().unwrap(), "--gap", "1",
                      "-o", o.to_str().unwrap(), "--mkdir"]);
    assert!(String::from_utf8_lossy(&out.stderr).contains("--atom-style"));
    assert!(!o.exists(), "nothing may be written");
}
