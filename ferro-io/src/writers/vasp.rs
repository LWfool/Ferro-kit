use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};

/// 写 VASP5 POSCAR 格式，坐标使用 Direct（分数坐标），原子按元素分组。
pub fn write_poscar(trajectory: &Trajectory, path: &str) -> Result<()> {
    let frame = trajectory.first().context("trajectory is empty")?;
    let cell = frame.cell.as_ref().context("frame has no cell (POSCAR requires periodic frame)")?;

    let file = File::create(path).with_context(|| format!("cannot create {path}"))?;
    let mut w = BufWriter::new(file);

    // Comment
    writeln!(w, "{}", trajectory.metadata.source.as_deref().unwrap_or("ferro"))?;
    // Scale
    writeln!(w, "  1.00000000000000")?;
    // Lattice vectors (row = lattice vector)
    for i in 0..3 {
        let r = cell.matrix.row(i);
        writeln!(w, "   {:>20.16}  {:>20.16}  {:>20.16}", r[0], r[1], r[2])?;
    }

    // Collect elements in first-appearance order
    let mut elem_order: Vec<&str> = Vec::new();
    for atom in &frame.atoms {
        if !elem_order.contains(&atom.element.as_str()) {
            elem_order.push(atom.element.as_str());
        }
    }
    let counts: Vec<usize> = elem_order.iter()
        .map(|e| frame.atoms.iter().filter(|a| a.element == *e).count())
        .collect();

    // Element and count lines (VASP5)
    writeln!(w, "  {}", elem_order.join("   "))?;
    writeln!(w, "  {}", counts.iter().map(|n| n.to_string()).collect::<Vec<_>>().join("   "))?;

    // Coordinates — atoms grouped by element (required by POSCAR format)
    writeln!(w, "Direct")?;
    for elem in &elem_order {
        for atom in frame.atoms.iter().filter(|a| a.element == *elem) {
            let f = cell.cartesian_to_fractional(atom.position);
            writeln!(w, "  {:>18.16}  {:>18.16}  {:>18.16}", f.x, f.y, f.z)?;
        }
    }

    // Optional velocities
    if let Some(vels) = &frame.velocities {
        writeln!(w)?;
        for vel in vels {
            writeln!(w, "  {:>18.16}  {:>18.16}  {:>18.16}", vel.x, vel.y, vel.z)?;
        }
    }

    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::vasp::read_poscar;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    fn bcc_traj() -> Trajectory {
        let cell = Cell::from_lengths_angles(2.87, 2.87, 2.87, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(1.435, 1.435, 1.435)));
        let mut traj = Trajectory::from_frame(frame);
        traj.metadata.source = Some("BCC Fe".to_string());
        traj
    }

    #[test]
    fn test_roundtrip() {
        let path = std::env::temp_dir().join("bcc_rt.poscar");
        let p = path.to_str().unwrap();
        let orig = bcc_traj();
        write_poscar(&orig, p).unwrap();

        let loaded = read_poscar(p).unwrap();
        let f = loaded.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
        let [a, ..] = f.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-10);
    }
}
