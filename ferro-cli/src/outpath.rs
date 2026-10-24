//! Creating output directories, with the user's consent.
//!
//! Every command takes its destination as `-o`, a path whose parent may not exist yet.
//! Creating it silently turns a typo into a stray directory tree that looks like a
//! successful run, so the terminal asks first. Scripts have no one to ask, hence
//! `--mkdir`: non-interactive runs must say up front that creating is fine.
//!
//! It is called **before the first input is read**, so a bad path fails next to the
//! other parameter errors rather than after an hour of analysis. `convert` is the
//! exception: it creates its one directory without asking, so it waits until there is
//! something to write and a failed read leaves no empty directory.

use std::io::{BufRead, IsTerminal, Write};
use std::path::Path;

use anyhow::{bail, Result};

/// Makes sure `dir` is a usable directory, asking before creating it.
pub fn ensure_dir(dir: &Path, mkdir: bool) -> Result<()> {
    if dir.as_os_str().is_empty() {
        return Ok(());
    }
    if dir.exists() {
        if !dir.is_dir() {
            bail!("-o '{}' exists and is not a directory", dir.display());
        }
        return Ok(());
    }
    if !mkdir && !approved(dir)? {
        bail!("no output directory: '{}' was not created", dir.display());
    }
    std::fs::create_dir_all(dir)
        .map_err(|e| anyhow::anyhow!("cannot create '{}': {e}", dir.display()))?;
    println!("Created: {}", dir.display());
    Ok(())
}

/// Asks on the terminal; refuses to guess when there is no terminal to ask.
///
/// The prompt goes to stderr so that `ferro ... > out.txt` still shows it, and the
/// answer defaults to no — the question exists to catch a mistyped path, and a
/// reflexive Enter should not create one.
fn approved(dir: &Path) -> Result<bool> {
    if !std::io::stdin().is_terminal() {
        bail!(
            "output directory '{}' does not exist, and there is no terminal to ask. \
             Pass --mkdir to create it.",
            dir.display()
        );
    }
    eprint!("Directory '{}' does not exist. Create it? [y/N] ", dir.display());
    std::io::stderr().flush().ok();
    let mut answer = String::new();
    std::io::stdin().lock().read_line(&mut answer)?;
    Ok(matches!(answer.trim(), "y" | "Y" | "yes" | "Yes"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("ferro_outpath_{}_{name}", std::process::id()))
    }

    #[test]
    fn an_existing_directory_needs_no_permission() {
        let d = tmp("existing");
        std::fs::create_dir_all(&d).unwrap();
        ensure_dir(&d, false).unwrap();
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn mkdir_creates_the_whole_chain() {
        let d = tmp("a").join("b").join("c");
        let _ = std::fs::remove_dir_all(tmp("a"));
        ensure_dir(&d, true).unwrap();
        assert!(d.is_dir());
        let _ = std::fs::remove_dir_all(tmp("a"));
    }

    /// A file where a directory belongs is a mistake worth naming, not something to
    /// discover halfway through writing the products.
    #[test]
    fn a_file_in_the_directorys_place_is_refused() {
        let f = tmp("not_a_dir");
        std::fs::write(&f, "x").unwrap();
        let err = ensure_dir(&f, true).unwrap_err().to_string();
        assert!(err.contains("is not a directory"), "{err}");
        let _ = std::fs::remove_file(&f);
    }

    /// The test runner has no terminal, so this exercises the non-interactive path:
    /// without `--mkdir` the run stops instead of guessing.
    #[test]
    fn without_a_terminal_and_without_mkdir_it_refuses() {
        let d = tmp("never_created");
        let _ = std::fs::remove_dir_all(&d);
        let err = ensure_dir(&d, false).unwrap_err().to_string();
        assert!(err.contains("--mkdir"), "{err}");
        assert!(!d.exists());
    }
}
