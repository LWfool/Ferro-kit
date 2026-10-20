//! Every message outside the products is one stderr line, `<source>[:<line>]: <category>: <message>`
//! (see `src/report.rs`); stdout keeps results and progress only.

use std::path::PathBuf;
use std::process::{Command, Output};

fn ferro(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_ferro")).args(args).output().unwrap()
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("ferro_messages_{name}"));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

const NPT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../tests/43Z43P15A_NPT_5.lammpstrj");

/// The batch driver's notes and per-input failures, and the run's failure count: stderr only,
/// so stdout holds nothing but progress.
#[test]
fn batch_notes_and_failures() {
    let d = scratch("batch");
    std::fs::create_dir_all(d.join("sub")).unwrap();
    let twin = d.join("sub/43Z43P15A_NPT_5.lammpstrj");
    std::fs::copy(NPT, &twin).unwrap();
    let broken = d.join("broken.lammpstrj");
    std::fs::write(&broken, "garbage\n").unwrap();
    let out = ferro(&["traj", "gr", "-i", NPT, NPT, twin.to_str().unwrap(), broken.to_str().unwrap(),
                      "-o", d.join("o").to_str().unwrap(), "--mkdir"]);
    assert_eq!(out.status.code(), Some(1));
    let err = text(&out.stderr);
    let lines: Vec<&str> = err.lines().collect();
    assert_eq!(lines[0], "ferro: FerroNote: 1 duplicate input(s) skipped", "{err}");
    assert_eq!(lines[1], format!("{}: FerroNote: shares the name '43Z43P15A_NPT_5' with an earlier input, \
                                  labelled '43Z43P15A_NPT_5_2'", twin.display()), "{err}");
    assert!(lines[2].starts_with(&format!("{}: FerroError: ", broken.display())) && lines[2].ends_with("; skipped"),
            "{err}");
    assert_eq!(*lines.last().unwrap(), "ferro: FerroError: 1 input(s) failed; see the messages above");
    assert!(!text(&out.stdout).contains("Note"), "notes left stdout");
}

/// A fatal error is one line, not anyhow's `Error:` + `Caused by:` block.
#[test]
fn fatal_error_is_one_line() {
    let out = ferro(&["traj", "gr", "-i", "no_such_file.lammpstrj"]);
    assert_eq!(out.status.code(), Some(1));
    assert_eq!(text(&out.stderr), "ferro: FerroError: no file matches 'no_such_file.lammpstrj'\n");
}

/// A reader's notice names the file as typed and the line it is about.
#[test]
fn reader_notice_carries_source_and_line() {
    let d = scratch("reader");
    let input = d.join("cut.xyz");
    std::fs::write(&input, "2\nc\nO 0 0 0\nH 1 0 0\n2\nc\nO 0 0 0\n").unwrap();
    let out = ferro(&["convert", "-i", input.to_str().unwrap(), "-o", d.join("o.xyz").to_str().unwrap()]);
    assert!(out.status.success());
    assert_eq!(
        text(&out.stderr),
        format!("{}:5: FerroWarning: last frame (frame 1) is incomplete and was dropped: 1 of 2 atoms\n",
                input.display())
    );
}
