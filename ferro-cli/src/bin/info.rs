use anyhow::Result;
use clap::Parser;
use ferro::io_dispatch::read_trajectory;
use ferro_core::Frame;
use ferro_io::LammpsUnits;
use std::collections::BTreeMap;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "mol-info", about = "Display structure/trajectory information")]
struct Cli {
    /// Input file (format auto-detected)
    #[arg(short, long)]
    input: PathBuf,

    /// Use LAMMPS metal units for dump files (velocities Å/ps, forces eV/Å)
    #[arg(long)]
    metal_units: bool,
}

fn main() -> Result<()> {
    let args = Cli::parse();
    let units = if args.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let traj = read_trajectory(&args.input, units)?;
    let n = traj.frames.len();

    println!("File:   {}", args.input.display());
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
        println!("  Volume: {:.4} Å³", cell.volume());
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
