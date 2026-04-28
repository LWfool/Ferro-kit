use std::collections::HashMap;
use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};

/// 写 LAMMPS data 文件，atom_style full，real 单位（长度 Å）。
/// 周期性帧写晶格信息；非周期性帧写最小包围盒。
pub fn write_lammps_data(trajectory: &Trajectory, path: &str) -> Result<()> {
    let frame = trajectory.first().context("trajectory is empty")?;

    let file = File::create(path).with_context(|| format!("cannot create {path}"))?;
    let mut w = BufWriter::new(file);

    let comment = trajectory.metadata.source.as_deref().unwrap_or("LAMMPS data file written by ferro");
    writeln!(w, "{comment}")?;
    writeln!(w)?;

    let n = frame.n_atoms();

    // Collect unique elements in first-appearance order → type IDs
    let mut elem_order: Vec<&str> = Vec::new();
    for atom in &frame.atoms {
        if !elem_order.contains(&atom.element.as_str()) {
            elem_order.push(atom.element.as_str());
        }
    }
    let n_types = elem_order.len();
    let elem_to_type: HashMap<&str, usize> = elem_order.iter()
        .enumerate()
        .map(|(i, e)| (*e, i + 1))
        .collect();


    writeln!(w, "{n} atoms")?;
    writeln!(w, "{n_types} atom types")?;
    writeln!(w)?;

    // Box bounds
    let (lx, ly, lz, xy, xz, yz) = match &frame.cell {
        Some(cell) => cell_to_lammps(cell),
        None => {
            // Non-periodic: compute bounding box
            let (minx, maxx, miny, maxy, minz, maxz) = bounding_box(frame);
            (maxx - minx, maxy - miny, maxz - minz, 0.0, 0.0, 0.0)
        }
    };

    let is_triclinic = xy != 0.0 || xz != 0.0 || yz != 0.0;

    writeln!(w, "{:.10} {:.10} xlo xhi", 0.0, lx)?;
    writeln!(w, "{:.10} {:.10} ylo yhi", 0.0, ly)?;
    writeln!(w, "{:.10} {:.10} zlo zhi", 0.0, lz)?;
    if is_triclinic {
        writeln!(w, "{:.10} {:.10} {:.10} xy xz yz", xy, xz, yz)?;
    }
    writeln!(w)?;

    // Masses section
    writeln!(w, "Masses")?;
    writeln!(w)?;
    for (i, elem) in elem_order.iter().enumerate() {
        let mass = ferro_core::data::elements::by_symbol(elem)
            .map(|e| e.atomic_mass)
            .unwrap_or(1.0);
        writeln!(w, "{} {:.4}  # {elem}", i + 1, mass)?;
    }
    writeln!(w)?;

    // Atoms section (full style)
    writeln!(w, "Atoms # full")?;
    writeln!(w)?;

    // Build LAMMPS cell for coordinate transformation
    let lammps_cell = lammps_cell_matrix(lx, ly, lz, xy, xz, yz);

    for (i, atom) in frame.atoms.iter().enumerate() {
        let tp = elem_to_type[atom.element.as_str()];
        let q = atom.charge.unwrap_or(0.0);

        // Transform position to LAMMPS coordinate frame
        let pos = match &frame.cell {
            Some(orig_cell) => {
                let frac = orig_cell.cartesian_to_fractional(atom.position);
                lammps_cell.fractional_to_cartesian(frac)
            }
            None => atom.position,
        };

        // id mol-id type charge x y z
        writeln!(w, "{} 1 {tp} {q:.6} {:.10} {:.10} {:.10}",
            i + 1, pos.x, pos.y, pos.z)?;
    }

    w.flush()?;
    Ok(())
}

/// Convert Cell to LAMMPS parameters: returns (lx, ly, lz, xy, xz, yz)
fn cell_to_lammps(cell: &ferro_core::Cell) -> (f64, f64, f64, f64, f64, f64) {
    let [a, b, c] = cell.lengths();
    let [alpha, beta, gamma] = cell.angles();
    let (al, be, ga) = (alpha.to_radians(), beta.to_radians(), gamma.to_radians());

    let lx = a;
    let xy = b * ga.cos();
    let xz = c * be.cos();
    let ly = (b * b - xy * xy).max(0.0).sqrt();
    let yz = if ly > 1e-10 { (b * c * al.cos() - xy * xz) / ly } else { 0.0 };
    let lz = (c * c - xz * xz - yz * yz).max(0.0).sqrt();

    (lx, ly, lz, xy, xz, yz)
}

fn lammps_cell_matrix(lx: f64, ly: f64, lz: f64, xy: f64, xz: f64, yz: f64) -> ferro_core::Cell {
    use nalgebra::Matrix3;
    ferro_core::Cell::from_matrix(Matrix3::new(
        lx,  0.0, 0.0,
        xy,  ly,  0.0,
        xz,  yz,  lz,
    ))
}

fn bounding_box(frame: &ferro_core::Frame) -> (f64, f64, f64, f64, f64, f64) {
    let mut mn = [f64::MAX; 3];
    let mut mx = [f64::MIN; 3];
    for a in &frame.atoms {
        mn[0] = mn[0].min(a.position.x);
        mn[1] = mn[1].min(a.position.y);
        mn[2] = mn[2].min(a.position.z);
        mx[0] = mx[0].max(a.position.x);
        mx[1] = mx[1].max(a.position.y);
        mx[2] = mx[2].max(a.position.z);
    }
    (mn[0], mx[0], mn[1], mx[1], mn[2], mx[2])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::lammps_data::read_lammps_data;
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
        let path = std::env::temp_dir().join("bcc_rt.lammps");
        let p = path.to_str().unwrap();
        write_lammps_data(&bcc_traj(), p).unwrap();

        let loaded = read_lammps_data(p).unwrap();
        let f = loaded.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
        let [a, ..] = f.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-4);
    }
}
