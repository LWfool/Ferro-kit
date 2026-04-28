use std::collections::HashMap;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{Context, Result};

const KCAL_TO_EV: f64 = 0.04336410;

/// LAMMPS unit system for dump files.
///
/// Positions are always Å in both systems. Differences:
/// - `Real`:  velocities Å/fs, forces kcal/(mol·Å)
/// - `Metal`: velocities Å/ps, forces eV/Å
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum LammpsUnits {
    #[default]
    Real,
    Metal,
}

pub fn read_lammps_dump(path: &str, units: LammpsUnits) -> Result<Trajectory> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path}"))?;
    parse_lammps_dump(&content, units).with_context(|| format!("parsing {path}"))
}

fn parse_lammps_dump(content: &str, units: LammpsUnits) -> Result<Trajectory> {
    let lines: Vec<&str> = content.lines().collect();
    let mut i = 0;
    let mut traj = Trajectory::new();

    while i < lines.len() {
        // Find ITEM: TIMESTEP
        if !lines[i].trim().starts_with("ITEM: TIMESTEP") { i += 1; continue; }
        i += 1;
        let _timestep: i64 = lines.get(i)
            .and_then(|l| l.trim().parse().ok())
            .unwrap_or(0);
        i += 1;

        // NUMBER OF ATOMS
        while i < lines.len() && !lines[i].contains("NUMBER OF ATOMS") { i += 1; }
        i += 1;
        let n: usize = lines.get(i)
            .and_then(|l| l.trim().parse().ok())
            .unwrap_or(0);
        i += 1;

        // BOX BOUNDS
        while i < lines.len() && !lines[i].contains("BOX BOUNDS") { i += 1; }
        let bounds_header = lines[i];
        let is_triclinic = bounds_header.contains("xy");
        i += 1;

        let mut lo = [0.0_f64; 3];
        let mut hi = [0.0_f64; 3];
        let mut tilt = [0.0_f64; 3]; // xy, xz, yz

        for dim in 0..3 {
            let line = lines.get(i).unwrap_or(&"");
            i += 1;
            let vals: Vec<f64> = line.split_whitespace()
                .map_while(|s| s.parse().ok())
                .collect();
            lo[dim] = *vals.first().unwrap_or(&0.0);
            hi[dim] = *vals.get(1).unwrap_or(&0.0);
            if vals.len() >= 3 {
                tilt[dim] = vals[2]; // xy on dim 0, xz on dim 1, yz on dim 2
            }
        }

        // User request C: also detect triclinic by value count (≥3 values per line)
        // already handled above: vals.len() >= 3 sets tilt

        let lx = hi[0] - lo[0];
        let ly = hi[1] - lo[1];
        let lz = hi[2] - lo[2];
        let (xy, xz, yz) = (tilt[0], tilt[1], tilt[2]);

        let cell = if is_triclinic || xy != 0.0 || xz != 0.0 || yz != 0.0 {
            Cell::from_matrix(Matrix3::new(
                lx,  0.0, 0.0,
                xy,  ly,  0.0,
                xz,  yz,  lz,
            ))
        } else {
            Cell::from_matrix(Matrix3::new(
                lx,  0.0, 0.0,
                0.0, ly,  0.0,
                0.0, 0.0, lz,
            ))
        };

        // ITEM: ATOMS col1 col2 ...
        while i < lines.len() && !lines[i].contains("ITEM: ATOMS") { i += 1; }
        let atoms_header = lines[i];
        i += 1;

        let col_names: Vec<&str> = atoms_header
            .trim_start_matches("ITEM: ATOMS")
            .split_whitespace()
            .collect();
        let col: HashMap<&str, usize> = col_names.iter()
            .enumerate()
            .map(|(idx, &name)| (name, idx))
            .collect();

        let get_col = |name: &str| col.get(name).copied();

        let mut atoms_raw: Vec<(usize, Atom, Option<Vector3<f64>>, Option<Vector3<f64>>)> = Vec::new();

        for _ in 0..n {
            let line = lines.get(i).unwrap_or(&"");
            i += 1;
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.is_empty() { continue; }

            let atom_id: usize = get_col("id")
                .and_then(|c| parts.get(c))
                .and_then(|s| s.parse().ok())
                .unwrap_or(0);

            let tp: usize = get_col("type")
                .and_then(|c| parts.get(c))
                .and_then(|s| s.parse().ok())
                .unwrap_or(1);

            let element = get_col("element")
                .and_then(|c| parts.get(c))
                .map(|s| s.to_string())
                .unwrap_or_else(|| format!("X{tp}"));

            // Position — try x/y/z first, then xs/ys/zs (scaled), then xu/yu/zu (unwrapped)
            let pos = if let (Some(cx), Some(cy), Some(cz)) =
                (get_col("x"), get_col("y"), get_col("z"))
            {
                let x: f64 = parts.get(cx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let y: f64 = parts.get(cy).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let z: f64 = parts.get(cz).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                Vector3::new(x, y, z)
            } else if let (Some(cx), Some(cy), Some(cz)) =
                (get_col("xs"), get_col("ys"), get_col("zs"))
            {
                // Scaled [0,1) → Cartesian via cell
                let sx: f64 = parts.get(cx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let sy: f64 = parts.get(cy).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let sz: f64 = parts.get(cz).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                cell.fractional_to_cartesian(Vector3::new(sx, sy, sz))
            } else if let (Some(cx), Some(cy), Some(cz)) =
                (get_col("xu"), get_col("yu"), get_col("zu"))
            {
                let x: f64 = parts.get(cx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let y: f64 = parts.get(cy).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let z: f64 = parts.get(cz).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                Vector3::new(x, y, z)
            } else {
                Vector3::zeros()
            };

            let mut atom = Atom::new(element, pos);
            // Charge
            if let Some(c) = get_col("q") {
                atom.charge = parts.get(c).and_then(|s| s.parse().ok());
            }

            // Velocity: real Å/fs (internal), metal Å/ps → ×1e-3
            let vel = if let (Some(vx), Some(vy), Some(vz)) =
                (get_col("vx"), get_col("vy"), get_col("vz"))
            {
                let x: f64 = parts.get(vx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let y: f64 = parts.get(vy).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let z: f64 = parts.get(vz).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let vscale = match units {
                    LammpsUnits::Real  => 1.0,
                    LammpsUnits::Metal => 1e-3,
                };
                Some(Vector3::new(x, y, z) * vscale)
            } else { None };

            // Force: real kcal/(mol·Å) → eV/Å, metal eV/Å already
            let fscale = match units {
                LammpsUnits::Real  => KCAL_TO_EV,
                LammpsUnits::Metal => 1.0,
            };
            let force = if let (Some(fx), Some(fy), Some(fz)) =
                (get_col("fx"), get_col("fy"), get_col("fz"))
            {
                let x: f64 = parts.get(fx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let y: f64 = parts.get(fy).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                let z: f64 = parts.get(fz).and_then(|s| s.parse().ok()).unwrap_or(0.0);
                Some(Vector3::new(x * fscale, y * fscale, z * fscale))
            } else { None };

            atoms_raw.push((atom_id, atom, vel, force));
        }

        // Sort by atom id
        atoms_raw.sort_by_key(|(id, _, _, _)| *id);

        let mut frame = Frame::with_cell(cell, [true; 3]);
        let mut all_vels = Vec::new();
        let mut all_forces = Vec::new();
        let mut has_vel = false;
        let mut has_force = false;

        for (_, atom, vel, force) in atoms_raw {
            frame.add_atom(atom);
            if let Some(v) = vel { all_vels.push(v); has_vel = true; }
            else { all_vels.push(Vector3::zeros()); }
            if let Some(f) = force { all_forces.push(f); has_force = true; }
            else { all_forces.push(Vector3::zeros()); }
        }

        if has_vel { frame.velocities = Some(all_vels); }
        if has_force { frame.forces = Some(all_forces); }

        traj.add_frame(frame);
    }

    Ok(traj)
}

#[cfg(test)]
mod tests {
    use super::*;

    const DUMP_ORTHO: &str = "ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 2.87
0 2.87
0 2.87
ITEM: ATOMS id type element x y z
1 1 Fe 0.0   0.0   0.0
2 1 Fe 1.435 1.435 1.435
ITEM: TIMESTEP
10
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 2.87
0 2.87
0 2.87
ITEM: ATOMS id type element x y z
1 1 Fe 0.01  0.0   0.0
2 1 Fe 1.445 1.435 1.435
";

    fn tmp(n: &str, c: &str) -> String {
        let p = std::env::temp_dir().join(n);
        std::fs::write(&p, c).unwrap();
        p.to_str().unwrap().to_string()
    }

    const DUMP_VEL: &str = "ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0 5.0
0 5.0
0 5.0
ITEM: ATOMS id type element x y z vx vy vz
1 1 Fe 0.0 0.0 0.0 2.0 0.0 0.0
";

    #[test]
    fn test_multiframe() {
        let traj = read_lammps_dump(&tmp("bcc.dump", DUMP_ORTHO), LammpsUnits::Real).unwrap();
        assert_eq!(traj.n_frames(), 2);
        assert_eq!(traj.first().unwrap().n_atoms(), 2);
        assert_eq!(traj.first().unwrap().atom(0).element, "Fe");
    }

    #[test]
    fn test_box_bounds() {
        let traj = read_lammps_dump(&tmp("box.dump", DUMP_ORTHO), LammpsUnits::Real).unwrap();
        let [a, ..] = traj.first().unwrap().cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-6);
    }

    #[test]
    fn test_velocity_real_units() {
        let traj = read_lammps_dump(&tmp("vel_real.dump", DUMP_VEL), LammpsUnits::Real).unwrap();
        let vx = traj.first().unwrap().velocities.as_ref().unwrap()[0].x;
        assert!((vx - 2.0).abs() < 1e-10, "real units: vx should be 2.0 Å/fs, got {vx}");
    }

    #[test]
    fn test_velocity_metal_units() {
        let traj = read_lammps_dump(&tmp("vel_metal.dump", DUMP_VEL), LammpsUnits::Metal).unwrap();
        let vx = traj.first().unwrap().velocities.as_ref().unwrap()[0].x;
        // 2.0 Å/ps × 1e-3 = 0.002 Å/fs
        assert!((vx - 0.002).abs() < 1e-12, "metal units: vx should be 0.002 Å/fs, got {vx}");
    }
}
