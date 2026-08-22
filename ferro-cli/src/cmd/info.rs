use anyhow::{bail, Result};
use clap::Args;
use crate::io_dispatch::read_trajectory;
use ferro_core::data::elements::by_symbol;
use ferro_core::units::AMU_ANG3_TO_G_CM3;
use ferro_core::Frame;
use ferro_io::LammpsUnits;
use std::collections::BTreeMap;
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct InfoCmd {
    /// Input file (format auto-detected; omit to show what this reports)
    #[arg(short, long)]
    pub input: Option<PathBuf>,

    /// Use LAMMPS metal units for dump files (velocities Å/ps, forces eV/Å)
    #[arg(long)]
    pub metal_units: bool,
}

/// True when `ferro info` was typed with no input: print the help page rather
/// than clap's bare "required argument" error.
pub fn wants_help(args: &InfoCmd) -> bool {
    args.input.is_none()
}

pub fn run(args: &InfoCmd) -> Result<()> {
    // input 为空由 main 分派到帮助页，这里的 bail 只是防御
    let Some(input) = &args.input else {
        bail!("info needs an input file: -i <FILE>");
    };

    let units = if args.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let traj = read_trajectory(input, units)?;
    let n = traj.frames.len();

    println!("File:   {}", input.display());
    println!("Frames: {n}");

    if let Some(frame) = traj.frames.first() {
        print_frame_info(frame, 0);
    }
    if n > 1 {
        if let Some(frame) = traj.frames.last() {
            print_frame_info(frame, n - 1);
        }
    }

    Ok(())
}

/// Total mass of a frame in amu, plus the atoms whose element is missing from the
/// element table.
///
/// `Atom::effective_mass` falls back to 1 amu for an unknown symbol, which drags
/// the density down with no other visible sign — so the fallbacks are counted per
/// symbol here and reported next to the density.
fn frame_mass_amu(frame: &Frame) -> (f64, BTreeMap<&str, usize>) {
    let mut fallbacks: BTreeMap<&str, usize> = BTreeMap::new();
    let mut total = 0.0;
    for a in &frame.atoms {
        total += a.effective_mass();
        // 只有既无显式质量、元素表又查不到的才是回退：显式质量是用户给的，不算未知
        if a.mass.is_none() && by_symbol(&a.element).is_none() {
            *fallbacks.entry(a.element.as_str()).or_insert(0) += 1;
        }
    }
    (total, fallbacks)
}

fn print_frame_info(frame: &Frame, idx: usize) {
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for a in &frame.atoms {
        *counts.entry(a.element.as_str()).or_insert(0) += 1;
    }
    let composition: Vec<String> = counts.iter().map(|(e, n)| format!("{e}: {n}")).collect();

    println!("\nFrame {idx}:");
    println!("  Atoms:  {}  ({})", frame.atoms.len(), composition.join(", "));

    if let Some(cell) = &frame.cell {
        let [a, b, c] = cell.lengths();
        let [al, be, ga] = cell.angles();
        println!("  Cell:   a={a:.4} b={b:.4} c={c:.4} Å   α={al:.2} β={be:.2} γ={ga:.2}°");

        let volume = cell.volume();
        println!("  Volume: {volume:.4} Å³");

        // 无体积就没有密度：宁可不打这一行，也不写占位符冒充「测到了」
        if volume > 0.0 {
            let (mass, fallbacks) = frame_mass_amu(frame);
            println!("  Density: {:.4} g/cm³", mass / volume * AMU_ANG3_TO_G_CM3);
            if !fallbacks.is_empty() {
                let n_atoms: usize = fallbacks.values().sum();
                let listed: Vec<String> =
                    fallbacks.iter().map(|(e, n)| format!("{e}×{n}")).collect();
                println!(
                    "           WARNING: {n_atoms} atom(s) not in the element table \
                     ({}) counted as 1 amu — the density is too low",
                    listed.join(", ")
                );
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
    if frame.forces.is_some() {
        println!("  Forces: yes");
    }
    if frame.velocities.is_some() {
        println!("  Velocities: yes");
    }
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
            frame.atoms.push(Atom::new(*e, Vector3::new(0.0, 0.0, 0.0)));
        }
        frame
    }

    fn density_g_cm3(frame: &Frame) -> f64 {
        let (mass, _) = frame_mass_amu(frame);
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
        frame.atoms[0].mass = Some(100.0);
        let (mass, fallbacks) = frame_mass_amu(&frame);
        assert!((mass - 100.0).abs() < 1e-12);
        assert!(fallbacks.is_empty(), "显式质量是用户给的，不该算未知元素");
    }

    #[test]
    fn test_unknown_element_falls_back_to_one_amu_and_is_counted() {
        // 未知符号在 effective_mass 里静默变成 1 amu，密度会偏低而无其他征兆，
        // 所以必须数出来告警 —— 这个测试钉住的正是那条告警的触发条件
        let frame = cubic_frame(10.0, &["Si", "Xx", "Xx"]);
        let (mass, fallbacks) = frame_mass_amu(&frame);
        assert!((mass - (28.085 + 2.0)).abs() < 1e-3);
        assert_eq!(fallbacks.len(), 1);
        assert_eq!(fallbacks["Xx"], 2);
    }

    #[test]
    fn test_wants_help_only_without_input() {
        let bare = InfoCmd { input: None, metal_units: false };
        assert!(wants_help(&bare));
        let with_input = InfoCmd { input: Some(PathBuf::from("a.xyz")), metal_units: false };
        assert!(!wants_help(&with_input));
    }
}
