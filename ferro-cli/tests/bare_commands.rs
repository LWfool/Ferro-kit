//! A command typed alone prints its page; options without an input are an error.
//!
//! Both used to exit 0, so `ferro traj gr --dr 0.1 && next-step` ran `next-step`
//! after printing a help page nobody asked for.

use std::process::Command;

fn ferro(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_ferro")).args(args).output().unwrap()
}

#[test]
fn the_command_name_alone_prints_its_page() {
    for args in [&["traj", "gr"][..], &["map", "chg-sdf"], &["dataset", "merge"], &["net"], &["info"]] {
        let out = ferro(args);
        assert!(out.status.success(), "{args:?}");
        assert!(!out.stdout.is_empty(), "{args:?} printed no page");
    }
}

#[test]
fn options_without_an_input_fail() {
    for args in [
        &["traj", "gr", "--dr", "0.1"][..],
        &["map", "chg-sdf", "--qn", "3"],
        &["dataset", "merge", "--seed", "1"],
        &["net", "--P-O=2.3"],
        &["convert", "-o", "out"],
    ] {
        let out = ferro(args);
        assert!(!out.status.success(), "{args:?} exited 0");
        let err = String::from_utf8_lossy(&out.stderr);
        assert!(err.contains("no input"), "{args:?}: {err}");
    }
}
