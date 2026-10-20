use anyhow::{bail, Result};
use clap::Args;
use crate::args::common::ReadArgs;
use crate::io_dispatch::read_trajectory;
use ferro_core::units::AMU_ANG3_TO_G_CM3;
use ferro_core::Frame;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

#[derive(Args, Debug)]
pub struct InfoCmd {
    /// Input file (format auto-detected; omit to show what this reports)
    #[arg(short, long)]
    pub input: Option<PathBuf>,

    #[command(flatten)]
    pub read: ReadArgs,
}

/// True when `ferro info` was typed with no input: print the help page rather
/// than clap's bare "required argument" error.
pub fn wants_help(args: &InfoCmd) -> bool {
    args.input.is_none()
}

/// Returns the number of failures: 1 when a density could not be computed, so the
/// exit code says so while the rest is still printed (review A-3).
pub fn run(args: &InfoCmd) -> Result<usize> {
    // input 为空由 main 分派到帮助页，这里的 bail 只是防御
    let Some(input) = &args.input else {
        bail!("info needs an input file: -i <FILE>");
    };

    let traj = read_trajectory(input, &args.read)?;
    let n = traj.frames.len();

    println!("File:   {}", input.display());
    println!("Frames: {n}");

    let mut ok = true;
    if let Some(frame) = traj.frames.first() {
        ok &= print_frame_info(frame, 0, input);
    }
    if n > 1 {
        if let Some(frame) = traj.frames.last() {
            ok &= print_frame_info(frame, n - 1, input);
        }
    }

    Ok(usize::from(!ok))
}

/// Total mass of a frame in amu, or the atoms with no mass counted per symbol.
///
/// An atom has no mass when its element is missing from the element table and no
/// mass was given (review A-3: it used to count as 1 amu, dragging the density down
/// with no other sign). `collect --type inspect` lists the same atoms in
/// `ferro_analysis::ml::diagnostics::frame_table`; two places, not merged (R6,
/// review E2-12).
fn frame_mass_amu(frame: &Frame) -> std::result::Result<f64, BTreeMap<&str, usize>> {
    let mut missing: BTreeMap<&str, usize> = BTreeMap::new();
    let mut total = 0.0;
    for a in frame.atoms() {
        match a.effective_mass() {
            Ok(m) => total += m,
            Err(_) => *missing.entry(a.element).or_insert(0) += 1,
        }
    }
    if missing.is_empty() { Ok(total) } else { Err(missing) }
}

/// Prints one frame; `false` when its density could not be computed.
fn print_frame_info(frame: &Frame, idx: usize, source: &Path) -> bool {
    let mut ok = true;
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for a in frame.atoms() {
        *counts.entry(a.element).or_insert(0) += 1;
    }
    let composition: Vec<String> = counts.iter().map(|(e, n)| format!("{e}: {n}")).collect();

    println!("\nFrame {idx}:");
    println!("  Atoms:  {}  ({})", frame.n_atoms(), composition.join(", "));

    if let Some(cell) = &frame.cell {
        let [a, b, c] = cell.lengths();
        let [al, be, ga] = cell.angles();
        println!("  Cell:   a={a:.4} b={b:.4} c={c:.4} Å   α={al:.2} β={be:.2} γ={ga:.2}°");

        let volume = cell.volume();
        println!("  Volume: {volume:.4} Å³");

        // No volume or no mass, no density: the line is left out rather than filled
        // with a placeholder posing as measured
        if volume > 0.0 {
            match frame_mass_amu(frame) {
                Ok(mass) => println!("  Density: {:.4} g/cm³", mass / volume * AMU_ANG3_TO_G_CM3),
                Err(missing) => {
                    let n_atoms: usize = missing.values().sum();
                    let listed: Vec<String> =
                        missing.iter().map(|(e, n)| format!("{e}×{n}")).collect();
                    crate::report::error(source.display(), format!(
                        "frame {idx}: no mass for {n_atoms} atom(s) ({}): not in the element \
                         table and no mass given; density not computed",
                        listed.join(", ")
                    ));
                    ok = false;
                }
            }
        }
    } else {
        println!("  Cell:   none (non-periodic)");
    }

    let pbc = frame.pbc;
    println!("  PBC:    [{}, {}, {}]", pbc[0], pbc[1], pbc[2]);

    if let Some(e) = frame.energy {
        println!("  Energy: {e:.6} eV");
    }
    if frame.forces().is_some() {
        println!("  Forces: yes");
    }
    if frame.velocities().is_some() {
        println!("  Velocities: yes");
    }
    ok
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell};
    use nalgebra::Vector3;

    /// 构造一个立方胞，按给定元素各放一个原子（位置无关紧要，密度只用总质量与体积）
    fn cubic_frame(a: f64, elements: &[&str]) -> Frame {
        let mut frame = Frame::new();
        frame.cell = Some(Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap());
        frame.pbc = [true; 3];
        for e in elements {
            frame.add_atom(Atom::new(*e, Vector3::new(0.0, 0.0, 0.0)));
        }
        frame
    }

    fn density_g_cm3(frame: &Frame) -> f64 {
        let mass = frame_mass_amu(frame).unwrap();
        mass / frame.cell.as_ref().unwrap().volume() * AMU_ANG3_TO_G_CM3
    }

    #[test]
    fn test_density_matches_handbook_values() {
        // NaCl: 立方 a=5.6402 Å，Z=4（4 个 NaCl 单元），实验密度 2.165 g/cm3
        let nacl = cubic_frame(5.6402, &["Na", "Cl", "Na", "Cl", "Na", "Cl", "Na", "Cl"]);
        assert!((density_g_cm3(&nacl) - 2.165).abs() < 5e-3);

        // 金刚石: a=3.5670 Å，Z=8 个 C，实验密度 3.515 g/cm3
        let diamond = cubic_frame(3.5670, &["C"; 8]);
        assert!((density_g_cm3(&diamond) - 3.515).abs() < 5e-3);
    }

    #[test]
    fn test_explicit_mass_overrides_table_and_is_not_a_fallback() {
        let mut frame = cubic_frame(10.0, &["Si"]);
        frame.sites_mut()[0].mass = Some(100.0);
        let mass = frame_mass_amu(&frame).expect("an explicit mass is the user's, not a missing one");
        assert!((mass - 100.0).abs() < 1e-12);
    }

    // Review A-3: an unknown symbol used to count as 1 amu, a density too low and nothing else
    #[test]
    fn test_unknown_element_has_no_mass_and_is_counted() {
        let frame = cubic_frame(10.0, &["Si", "Xx", "Xx"]);
        let missing = frame_mass_amu(&frame).unwrap_err();
        assert_eq!(missing.len(), 1);
        assert_eq!(missing["Xx"], 2);
    }

    #[test]
    fn test_wants_help_only_without_input() {
        let bare = InfoCmd { input: None, read: ReadArgs::default() };
        assert!(wants_help(&bare));
        let with_input = InfoCmd { input: Some(PathBuf::from("a.xyz")), read: ReadArgs::default() };
        assert!(!wants_help(&with_input));
    }
}
