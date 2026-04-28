use std::collections::HashMap;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{Context, Result};

pub fn read_lammps_data(path: &str) -> Result<Trajectory> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path}"))?;
    parse_lammps_data(&content).with_context(|| format!("parsing {path}"))
}

fn parse_lammps_data(content: &str) -> Result<Trajectory> {
    let lines: Vec<&str> = content.lines().collect();
    let mut i = 0;

    let skip = |l: &str| -> bool { l.trim().is_empty() || l.trim().starts_with('#') };

    // ── Header ────────────────────────────────────────────────────────────────
    // First non-blank line is the comment
    while i < lines.len() && skip(lines[i]) { i += 1; }
    let comment = lines[i].trim().to_string();
    i += 1;

    // Box bounds
    let mut xlo = 0.0_f64; let mut xhi = 0.0_f64;
    let mut ylo = 0.0_f64; let mut yhi = 0.0_f64;
    let mut zlo = 0.0_f64; let mut zhi = 0.0_f64;
    let mut xy = 0.0_f64; let mut xz = 0.0_f64; let mut yz = 0.0_f64;
    let mut is_triclinic = false;

    while i < lines.len() {
        let line = strip_comment(lines[i]);
        if line.is_empty() { i += 1; continue; }

        if line.ends_with("atoms") && !line.contains("atom types") {
            // count stored but not needed — we read atoms greedily
        } else if line.ends_with("atom types") {
            // type count stored but not needed
        } else if line.contains("xlo xhi") {
            let v = floats(line, 2)?;
            xlo = v[0]; xhi = v[1];
        } else if line.contains("ylo yhi") {
            let v = floats(line, 2)?;
            ylo = v[0]; yhi = v[1];
        } else if line.contains("zlo zhi") {
            let v = floats(line, 2)?;
            zlo = v[0]; zhi = v[1];
        } else if line.contains("xy xz yz") {
            let v = floats(line, 3)?;
            xy = v[0]; xz = v[1]; yz = v[2];
            is_triclinic = true;
        } else if line == "Masses" || line == "Atoms" || line == "Atoms # full"
                || line == "Atoms # atomic" || line == "Atoms # charge"
                || line.starts_with("Masses") || line.starts_with("Atoms") {
            break;
        }
        i += 1;
    }

    // Build cell
    let lx = xhi - xlo;
    let ly = yhi - ylo;
    let lz = zhi - zlo;
    let cell = if is_triclinic {
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

    // ── Masses section ────────────────────────────────────────────────────────
    let mut type_element: HashMap<usize, String> = HashMap::new();
    while i < lines.len() {
        let raw = lines[i];
        let line = strip_comment(raw).trim();
        if line == "Masses" { i += 1; continue; }
        if line.is_empty() { i += 1; continue; }

        // Check if this line is a section header (Atoms, Velocities, ...)
        if is_section_header(line) { break; }

        // Masses line: type_id mass  # element_comment
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let Ok(tid) = parts[0].parse::<usize>() {
                let mass: f64 = parts[1].parse().unwrap_or(1.0);
                // Try comment for element symbol: "# Fe" or "# Fe ..."
                let elem = raw.find('#')
                    .and_then(|p| {
                        raw[p+1..].trim().split_whitespace().next().map(|s| s.to_string())
                    })
                    .unwrap_or_else(|| element_from_mass(mass).to_string());
                type_element.insert(tid, elem);
            }
        }
        i += 1;
    }

    // ── Atoms section ─────────────────────────────────────────────────────────
    // Detect atom_style from section header comment
    // "Atoms # full" → full, "Atoms # atomic" → atomic, "Atoms # charge" → charge
    let mut style = AtomStyle::Full; // default per user preference

    while i < lines.len() {
        let line = strip_comment(lines[i]).trim();
        if line.is_empty() { i += 1; continue; }
        if line.starts_with("Atoms") {
            if lines[i].contains("atomic") { style = AtomStyle::Atomic; }
            else if lines[i].contains("charge") { style = AtomStyle::Charge; }
            i += 1; break;
        }
        i += 1;
    }

    let mut atoms_raw: Vec<(usize, usize, f64, f64, f64, f64)> = Vec::new(); // (id, type, charge, x, y, z)

    while i < lines.len() {
        let line = strip_comment(lines[i]).trim();
        i += 1;
        if line.is_empty() { continue; }
        if is_section_header(line) { break; }

        let parts: Vec<&str> = line.split_whitespace().collect();
        match style {
            AtomStyle::Atomic => {
                // id type x y z [ix iy iz]
                if parts.len() < 5 { continue; }
                let id: usize = parts[0].parse().unwrap_or(0);
                let tp: usize = parts[1].parse().unwrap_or(1);
                let x: f64 = parts[2].parse().unwrap_or(0.0);
                let y: f64 = parts[3].parse().unwrap_or(0.0);
                let z: f64 = parts[4].parse().unwrap_or(0.0);
                atoms_raw.push((id, tp, 0.0, x, y, z));
            }
            AtomStyle::Charge => {
                // id type charge x y z [ix iy iz]
                if parts.len() < 6 { continue; }
                let id: usize = parts[0].parse().unwrap_or(0);
                let tp: usize = parts[1].parse().unwrap_or(1);
                let q: f64 = parts[2].parse().unwrap_or(0.0);
                let x: f64 = parts[3].parse().unwrap_or(0.0);
                let y: f64 = parts[4].parse().unwrap_or(0.0);
                let z: f64 = parts[5].parse().unwrap_or(0.0);
                atoms_raw.push((id, tp, q, x, y, z));
            }
            AtomStyle::Full => {
                // id mol-id type charge x y z [ix iy iz]
                if parts.len() < 7 { continue; }
                let id: usize = parts[0].parse().unwrap_or(0);
                let tp: usize = parts[2].parse().unwrap_or(1);
                let q: f64 = parts[3].parse().unwrap_or(0.0);
                let x: f64 = parts[4].parse().unwrap_or(0.0);
                let y: f64 = parts[5].parse().unwrap_or(0.0);
                let z: f64 = parts[6].parse().unwrap_or(0.0);
                atoms_raw.push((id, tp, q, x, y, z));
            }
        }
    }

    // Sort by atom id
    atoms_raw.sort_by_key(|(id, _, _, _, _, _)| *id);

    let mut frame = Frame::with_cell(cell, [true; 3]);
    for (_, tp, q, x, y, z) in atoms_raw {
        let elem = type_element.get(&tp)
            .cloned()
            .unwrap_or_else(|| format!("X{tp}"));
        let mut atom = Atom::new(elem, Vector3::new(x, y, z));
        if q != 0.0 { atom.charge = Some(q); }
        frame.add_atom(atom);
    }

    let mut traj = Trajectory::from_frame(frame);
    if !comment.is_empty() { traj.metadata.source = Some(comment); }
    Ok(traj)
}

#[derive(Clone, Copy)]
enum AtomStyle { Atomic, Charge, Full }

fn strip_comment(l: &str) -> &str {
    l.split('#').next().unwrap_or("").trim()
}

fn is_section_header(line: &str) -> bool {
    matches!(line, "Masses" | "Atoms" | "Velocities" | "Bonds" | "Angles"
        | "Dihedrals" | "Impropers" | "Pair Coeffs" | "Bond Coeffs"
        | "Angle Coeffs")
    || line.starts_with("Atoms #")
}

fn floats(line: &str, min: usize) -> Result<Vec<f64>> {
    let v: Vec<f64> = line.split_whitespace()
        .map_while(|s| s.parse::<f64>().ok())
        .collect();
    anyhow::ensure!(v.len() >= min, "expected ≥{min} floats in {line:?}");
    Ok(v)
}

fn element_from_mass(mass: f64) -> &'static str {
    ferro_core::data::elements::ELEMENTS
        .iter()
        .filter(|e| (e.atomic_mass - mass).abs() < 0.5)
        .min_by(|a, b| {
            (a.atomic_mass - mass).abs()
                .partial_cmp(&(b.atomic_mass - mass).abs())
                .unwrap()
        })
        .map(|e| e.symbol)
        .unwrap_or("X")
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER_FULL: &str = "LAMMPS data file

6 atoms
2 atom types

0.0 20.0 xlo xhi
0.0 20.0 ylo yhi
0.0 20.0 zlo zhi

Masses

1 15.999  # O
2 1.008   # H

Atoms # full

1 1 1 -0.834 0.000  0.000  0.000
2 1 2  0.417 0.758  0.587  0.000
3 1 2  0.417 -0.758 0.587  0.000
4 2 1 -0.834 10.000 0.000  0.000
5 2 2  0.417 10.758 0.587  0.000
6 2 2  0.417 9.242  0.587  0.000
";

    fn tmp(n: &str, c: &str) -> String {
        let p = std::env::temp_dir().join(n);
        std::fs::write(&p, c).unwrap();
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn test_water_full() {
        let traj = read_lammps_data(&tmp("water.lammps", WATER_FULL)).unwrap();
        let frame = traj.first().unwrap();
        assert_eq!(frame.n_atoms(), 6);
        assert_eq!(frame.atom(0).element, "O");
        assert_eq!(frame.atom(1).element, "H");
        assert!((frame.atom(0).charge.unwrap() - (-0.834)).abs() < 1e-6);
        let [a, ..] = frame.cell.as_ref().unwrap().lengths();
        assert!((a - 20.0).abs() < 1e-6);
    }
}
