//! End-to-end cover for `ferro bader`: CHGCAR in, three reports out.
//!
//! The unit tests in `ferro-analysis` build their grid in memory, so nothing used to
//! exercise the CHGCAR reader, the analyzer and the report rendering together. The
//! fixture is a synthetic 8×8×8 grid with two Gaussian peaks — small enough to read
//! (9 KB) and to run under every method.

use std::path::Path;

use ferro_analysis::{BaderAnalyzer, BaderMethod};
use ferro_io::read_chgcar;

const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../tests/CHGCAR_2atoms");

/// Charge is conserved: what the grid holds ends up on the atoms or in the vacuum.
///
/// The grid stores `ρ × V_cell` (VASP convention), so the electron count is the mean
/// of the stored values — that mean is the number every method has to reproduce.
#[test]
fn every_method_conserves_the_grid_charge() {
    let (frame, chg) = read_chgcar(Path::new(FIXTURE)).unwrap();
    let expected: f64 = chg.rho.iter().sum::<f64>() / chg.nrho as f64;

    for method in [BaderMethod::OnGrid, BaderMethod::NearGrid, BaderMethod::Weight] {
        let (_, chg) = read_chgcar(Path::new(FIXTURE)).unwrap();
        let result = BaderAnalyzer::new(chg, frame.clone()).method(method).run();
        let got: f64 = result.ionchg.iter().sum();
        assert!(
            (got - expected).abs() < 1e-6,
            "{method:?}: {got} e on the atoms, grid holds {expected} e"
        );
        assert!(result.nvols >= frame.n_atoms(), "{method:?}: {} volumes", result.nvols);
    }
}

/// The vacuum bucket is empty on a fixture whose density never drops below the
/// threshold — every voxel belongs to a basin.
///
/// `BaderMethod::Weight` is left out: it reports the **last Bader volume's** charge as
/// the vacuum charge (`bader_weight.rs:265` indexes `volchg[nvols]` while that array is
/// 1-indexed there), so its ACF footer prints a `Total` 1.5× the real electron count on
/// this fixture. Pre-existing, see `dev/issues.md`.
#[test]
fn the_grid_methods_leave_the_vacuum_empty() {
    let (frame, _) = read_chgcar(Path::new(FIXTURE)).unwrap();
    for method in [BaderMethod::OnGrid, BaderMethod::NearGrid] {
        let (_, chg) = read_chgcar(Path::new(FIXTURE)).unwrap();
        let result = BaderAnalyzer::new(chg, frame.clone()).method(method).run();
        assert_eq!(result.vacchg, 0.0, "{method:?}");
    }
}

/// The three reports are what external tools parse, so their shape is part of the
/// contract: one ACF row per atom, one BCF row per Bader volume, one AVF row per atom.
#[test]
fn the_three_reports_keep_their_shape() {
    let (frame, chg) = read_chgcar(Path::new(FIXTURE)).unwrap();
    let result = BaderAnalyzer::new(chg, frame.clone()).method(BaderMethod::NearGrid).run();

    let rows = |text: &str| text.lines().filter(|l| l.starts_with(' ')).count();

    let acf = result.acf_text(&frame);
    assert!(acf.starts_with("# ACF.dat"), "{acf}");
    assert_eq!(rows(&acf), frame.n_atoms());
    assert!(acf.contains("# Vacuum charge:") && acf.contains("# Total:"));

    let bcf = result.bcf_text();
    assert!(bcf.starts_with("# BCF.dat"), "{bcf}");
    assert_eq!(rows(&bcf), result.nvols);

    let avf = result.avf_text();
    assert!(avf.starts_with("# AVF.dat"), "{avf}");
    assert_eq!(rows(&avf), frame.n_atoms());
}
