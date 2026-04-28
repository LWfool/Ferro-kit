use crate::readers::lammps_dump::LammpsUnits;
use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};

const EV_TO_KCAL: f64 = 1.0 / 0.04336410; // eV/Å → kcal/(mol·Å)

/// 写 LAMMPS dump 文件。
/// 包含列：id type element x y z [vx vy vz] [fx fy fz] [q]
pub fn write_lammps_dump(trajectory: &Trajectory, path: &str, units: LammpsUnits) -> Result<()> {
    let file = File::create(path).with_context(|| format!("cannot create {path}"))?;
    let mut w = BufWriter::new(file);

    for (ts, frame) in trajectory.frames.iter().enumerate() {
        let n = frame.n_atoms();

        writeln!(w, "ITEM: TIMESTEP")?;
        writeln!(w, "{ts}")?;
        writeln!(w, "ITEM: NUMBER OF ATOMS")?;
        writeln!(w, "{n}")?;

        // BOX BOUNDS
        let (lx, ly, lz, xy, xz, yz) = match &frame.cell {
            Some(cell) => cell_to_lammps(cell),
            None => {
                let (mnx, mxx, mny, mxy, mnz, mxz) = bounding_box(frame);
                (mxx - mnx, mxy - mny, mxz - mnz, 0.0, 0.0, 0.0)
            }
        };

        let is_triclinic = xy != 0.0 || xz != 0.0 || yz != 0.0;
        if is_triclinic {
            writeln!(w, "ITEM: BOX BOUNDS xy xz yz pp pp pp")?;
            writeln!(w, "{:.10} {:.10} {:.10}", 0.0, lx, xy)?;
            writeln!(w, "{:.10} {:.10} {:.10}", 0.0, ly, xz)?;
            writeln!(w, "{:.10} {:.10} {:.10}", 0.0, lz, yz)?;
        } else {
            writeln!(w, "ITEM: BOX BOUNDS pp pp pp")?;
            writeln!(w, "{:.10} {:.10}", 0.0, lx)?;
            writeln!(w, "{:.10} {:.10}", 0.0, ly)?;
            writeln!(w, "{:.10} {:.10}", 0.0, lz)?;
        }

        // ATOMS header
        let has_vel = frame.velocities.is_some();
        let has_force = frame.forces.is_some();
        let has_charge = frame.atoms.iter().any(|a| a.charge.is_some());

        let mut header = "ITEM: ATOMS id type element x y z".to_string();
        if has_vel { header.push_str(" vx vy vz"); }
        if has_force { header.push_str(" fx fy fz"); }
        if has_charge { header.push_str(" q"); }
        writeln!(w, "{header}")?;

        // Atom type assignment (element → integer type)
        let mut elem_types: Vec<&str> = Vec::new();
        for atom in &frame.atoms {
            if !elem_types.contains(&atom.element.as_str()) {
                elem_types.push(atom.element.as_str());
            }
        }

        let lammps_cell = lammps_cell_matrix(lx, ly, lz, xy, xz, yz);

        for (i, atom) in frame.atoms.iter().enumerate() {
            let tp = elem_types.iter().position(|e| *e == atom.element).unwrap_or(0) + 1;

            // Transform to LAMMPS frame
            let pos = match &frame.cell {
                Some(orig) => {
                    let frac = orig.cartesian_to_fractional(atom.position);
                    lammps_cell.fractional_to_cartesian(frac)
                }
                None => atom.position,
            };

            let mut line = format!("{} {tp} {} {:.10} {:.10} {:.10}",
                i + 1, atom.element, pos.x, pos.y, pos.z);

            if has_vel {
                let v = frame.velocities.as_ref()
                    .and_then(|vv| vv.get(i))
                    .copied()
                    .unwrap_or_default();
                // real: Å/fs (no-op), metal: Å/fs → Å/ps (×1000)
                let vscale = match units {
                    LammpsUnits::Real  => 1.0,
                    LammpsUnits::Metal => 1000.0,
                };
                line.push_str(&format!(" {:.10} {:.10} {:.10}", v.x * vscale, v.y * vscale, v.z * vscale));
            }
            if has_force {
                let f = frame.forces.as_ref()
                    .and_then(|ff| ff.get(i))
                    .copied()
                    .unwrap_or_default();
                // real: eV/Å → kcal/(mol·Å), metal: eV/Å (no-op)
                let fscale = match units {
                    LammpsUnits::Real  => EV_TO_KCAL,
                    LammpsUnits::Metal => 1.0,
                };
                line.push_str(&format!(" {:.10} {:.10} {:.10}", f.x * fscale, f.y * fscale, f.z * fscale));
            }
            if has_charge {
                line.push_str(&format!(" {:.6}", atom.charge.unwrap_or(0.0)));
            }

            writeln!(w, "{line}")?;
        }
    }

    w.flush()?;
    Ok(())
}

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
    ferro_core::Cell::from_matrix(Matrix3::new(lx, 0.0, 0.0, xy, ly, 0.0, xz, yz, lz))
}

fn bounding_box(frame: &ferro_core::Frame) -> (f64, f64, f64, f64, f64, f64) {
    let mut mn = [f64::MAX; 3]; let mut mx = [f64::MIN; 3];
    for a in &frame.atoms {
        mn[0]=mn[0].min(a.position.x); mn[1]=mn[1].min(a.position.y); mn[2]=mn[2].min(a.position.z);
        mx[0]=mx[0].max(a.position.x); mx[1]=mx[1].max(a.position.y); mx[2]=mx[2].max(a.position.z);
    }
    (mn[0],mx[0],mn[1],mx[1],mn[2],mx[2])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::lammps_dump::read_lammps_dump;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    fn bcc_traj() -> Trajectory {
        let cell = Cell::from_lengths_angles(2.87, 2.87, 2.87, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(1.435, 1.435, 1.435)));
        let mut traj = Trajectory::new();
        traj.add_frame(frame.clone());
        traj.add_frame(frame);
        traj
    }

    #[test]
    fn test_roundtrip() {
        use crate::readers::lammps_dump::LammpsUnits;
        let path = std::env::temp_dir().join("bcc_rt.dump");
        let p = path.to_str().unwrap();
        let orig = bcc_traj();
        write_lammps_dump(&orig, p, LammpsUnits::Real).unwrap();

        let loaded = read_lammps_dump(p, LammpsUnits::Real).unwrap();
        assert_eq!(loaded.n_frames(), 2);
        let f = loaded.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
    }
}
