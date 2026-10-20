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
