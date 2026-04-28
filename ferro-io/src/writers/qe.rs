use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};

/// 写 Quantum ESPRESSO pw.x 输入文件（ibrav=0，ATOMIC_POSITIONS angstrom）。
pub fn write_qe_input(trajectory: &Trajectory, path: &str) -> Result<()> {
    let frame = trajectory.first().context("trajectory is empty")?;

    let file = File::create(path).with_context(|| format!("cannot create {path}"))?;
    let mut w = BufWriter::new(file);

    let prefix = trajectory.metadata.source.as_deref().unwrap_or("ferro");

    // Collect unique element labels in first-appearance order
    let mut elem_order: Vec<&str> = Vec::new();
    for atom in &frame.atoms {
        if !elem_order.contains(&atom.element.as_str()) {
            elem_order.push(atom.element.as_str());
        }
    }
    let ntyp = elem_order.len();
    let nat = frame.n_atoms();

    // &CONTROL
    writeln!(w, "&CONTROL")?;
    writeln!(w, "  calculation = 'scf',")?;
    writeln!(w, "  prefix = '{prefix}',")?;
    writeln!(w, "/")?;
    writeln!(w)?;

    // &SYSTEM
    writeln!(w, "&SYSTEM")?;
    writeln!(w, "  ibrav = 0,")?;
    writeln!(w, "  nat = {nat},")?;
    writeln!(w, "  ntyp = {ntyp},")?;
    let charge = frame.charge;
    if charge != 0 {
        writeln!(w, "  tot_charge = {charge},")?;
    }
    let mult = frame.multiplicity;
    if mult != 1 {
        writeln!(w, "  tot_magnetization = {},", mult - 1)?;
    }
    writeln!(w, "/")?;
    writeln!(w)?;

    // &ELECTRONS (empty placeholder)
    writeln!(w, "&ELECTRONS")?;
    writeln!(w, "/")?;
    writeln!(w)?;

    // ATOMIC_SPECIES
    writeln!(w, "ATOMIC_SPECIES")?;
    for elem in &elem_order {
        let mass = ferro_core::data::elements::by_symbol(elem)
            .map(|e| e.atomic_mass)
            .unwrap_or(1.0);
        writeln!(w, "  {elem}  {mass:.4}  {elem}.UPF")?;
    }
    writeln!(w)?;

    // CELL_PARAMETERS angstrom
    if let Some(cell) = &frame.cell {
        writeln!(w, "CELL_PARAMETERS {{angstrom}}")?;
        for i in 0..3 {
            let r = cell.matrix.row(i);
            writeln!(w, "  {:.10}  {:.10}  {:.10}", r[0], r[1], r[2])?;
        }
        writeln!(w)?;
    }

    // ATOMIC_POSITIONS angstrom
    writeln!(w, "ATOMIC_POSITIONS {{angstrom}}")?;
    for atom in &frame.atoms {
        writeln!(w, "  {}  {:.10}  {:.10}  {:.10}",
            atom.element,
            atom.position.x, atom.position.y, atom.position.z)?;
    }
    writeln!(w)?;

    // K_POINTS gamma (placeholder)
    writeln!(w, "K_POINTS {{gamma}}")?;

    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::qe::read_qe_input;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    fn bcc_traj() -> Trajectory {
        let cell = Cell::from_lengths_angles(2.87, 2.87, 2.87, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(1.435, 1.435, 1.435)));
        Trajectory::from_frame(frame)
    }

    #[test]
    fn test_roundtrip() {
        let path = std::env::temp_dir().join("bcc_rt.qe");
        let p = path.to_str().unwrap();
        let orig = bcc_traj();
        write_qe_input(&orig, p).unwrap();

        let loaded = read_qe_input(p).unwrap();
        let f = loaded.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
        assert!((f.atom(1).position.x - 1.435).abs() < 1e-6);
        let [a, ..] = f.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-4);
    }
}
