use std::collections::HashMap;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{bail, Context, Result};

pub fn read_extxyz(path: &str) -> Result<Trajectory> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path}"))?;
    parse_extxyz(&content).with_context(|| format!("parsing {path}"))
}

fn parse_extxyz(content: &str) -> Result<Trajectory> {
    let mut lines = content.lines().peekable();
    let mut traj = Trajectory::new();

    loop {
        // Skip blank lines between frames
        while lines.peek().map(|l| l.trim().is_empty()).unwrap_or(false) {
            lines.next();
        }
        let count_line = match lines.next() {
            None => break,
            Some(l) => l.trim(),
        };
        if count_line.is_empty() { break; }

        let n: usize = count_line.parse()
            .with_context(|| format!("expected atom count, got {count_line:?}"))?;

        let comment = lines.next().context("missing comment line")?;
        let kv = parse_comment(comment);

        // Cell and pbc
        let (cell, pbc) = match kv.get("lattice") {
            Some(lat) => {
                let c = parse_lattice(lat).context("invalid Lattice")?;
                let pbc = kv.get("pbc").map(|s| parse_pbc(s)).unwrap_or([true; 3]);
                (Some(c), pbc)
            }
            None => (None, [false; 3]),
        };

        // Scalar properties in comment
        let energy: Option<f64> = kv.get("energy").and_then(|s| s.parse().ok());
        let stress: Option<Matrix3<f64>> = kv.get("stress")
            .or_else(|| kv.get("virial"))
            .and_then(|s| parse_matrix9(s));

        // Properties column spec
        let props = kv.get("properties").map(|s| parse_properties(s))
            .unwrap_or_else(|| vec![
                ("species".to_string(), 'S', 1),
                ("pos".to_string(), 'R', 3),
            ]);
        let ncols: usize = props.iter().map(|(_, _, c)| c).sum();

        // Build column → field map
        let mut col_idx: HashMap<&str, usize> = HashMap::new();
        let mut offset = 0;
        for (name, _, count) in &props {
            col_idx.insert(name.as_str(), offset);
            offset += count;
        }

        let find = |name: &str| col_idx.get(name).copied();

        let mut frame = match cell {
            Some(c) => Frame::with_cell(c, pbc),
            None => Frame::new(),
        };
        frame.energy = energy;
        frame.stress = stress;

        let mut all_forces: Vec<Vector3<f64>> = Vec::new();
        let mut all_vels: Vec<Vector3<f64>> = Vec::new();

        for i in 0..n {
            let line = lines.next()
                .with_context(|| format!("missing atom line {i}"))?;
            let cols: Vec<&str> = line.split_whitespace().collect();
            if cols.len() < ncols {
                bail!("atom line {i}: expected {ncols} columns, got {}", cols.len());
            }

            let element = find("species")
                .map(|c| cols[c].to_string())
                .unwrap_or_else(|| "X".to_string());

            let pos = if let Some(c) = find("pos") {
                let x: f64 = cols[c].parse().context("invalid pos.x")?;
                let y: f64 = cols[c+1].parse().context("invalid pos.y")?;
                let z: f64 = cols[c+2].parse().context("invalid pos.z")?;
                Vector3::new(x, y, z)
            } else {
                Vector3::zeros()
            };

            let mut atom = Atom::new(element, pos);
            if let Some(c) = find("charges") {
                atom.charge = cols[c].parse().ok();
            }
            if let Some(c) = find("masses") {
                atom.mass = cols[c].parse().ok();
            }
            if let Some(c) = find("magmoms") {
                atom.magmom = cols[c].parse().ok();
            }
            frame.add_atom(atom);

            if let Some(c) = find("forces") {
                let fx: f64 = cols[c].parse().unwrap_or(0.0);
                let fy: f64 = cols[c+1].parse().unwrap_or(0.0);
                let fz: f64 = cols[c+2].parse().unwrap_or(0.0);
                all_forces.push(Vector3::new(fx, fy, fz));
            }
            if let Some(c) = find("velocities").or_else(|| find("momenta")) {
                let vx: f64 = cols[c].parse().unwrap_or(0.0);
                let vy: f64 = cols[c+1].parse().unwrap_or(0.0);
                let vz: f64 = cols[c+2].parse().unwrap_or(0.0);
                all_vels.push(Vector3::new(vx, vy, vz));
            }
        }

        if all_forces.len() == n { frame.forces = Some(all_forces); }
        if all_vels.len() == n { frame.velocities = Some(all_vels); }

        if traj.n_frames() == 0 {
            if let Some(src) = kv.get("config_type").or_else(|| kv.get("comment")) {
                traj.metadata.source = Some(src.clone());
            }
        }
        traj.add_frame(frame);
    }

    Ok(traj)
}

// ─── Comment line parser ──────────────────────────────────────────────────────

fn parse_comment(line: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    let mut s = line.trim();
    while !s.is_empty() {
        // Find next '='
        let eq = match s.find('=') {
            Some(p) => p,
            None => break,
        };
        let key = s[..eq].trim().to_lowercase();
        s = &s[eq + 1..];
        let (val, rest) = read_value(s);
        if !key.is_empty() { map.insert(key, val); }
        s = rest.trim_start();
    }
    map
}

fn read_value(s: &str) -> (String, &str) {
    let s = s.trim_start();
    if let Some(q) = s.chars().next().filter(|&c| c == '"' || c == '\'') {
        let inner = &s[1..];
        // Find closing quote (not escaped)
        let end = inner.find(q).unwrap_or(inner.len());
        (inner[..end].to_string(), &inner[end + 1..])
    } else {
        let end = s.find(char::is_whitespace).unwrap_or(s.len());
        (s[..end].to_string(), &s[end..])
    }
}

fn parse_lattice(s: &str) -> Result<Cell> {
    let v: Vec<f64> = s.split_whitespace()
        .map(|x| x.parse::<f64>().context("float"))
        .collect::<Result<_>>()?;
    anyhow::ensure!(v.len() == 9, "Lattice must have 9 values");
    Ok(Cell::from_matrix(Matrix3::new(
        v[0], v[1], v[2],
        v[3], v[4], v[5],
        v[6], v[7], v[8],
    )))
}

fn parse_pbc(s: &str) -> [bool; 3] {
    let mut pbc = [false; 3];
    for (i, tok) in s.split_whitespace().take(3).enumerate() {
        pbc[i] = matches!(tok.to_uppercase().as_str(), "T" | "TRUE" | "1");
    }
    pbc
}

fn parse_matrix9(s: &str) -> Option<Matrix3<f64>> {
    let v: Vec<f64> = s.split_whitespace()
        .filter_map(|x| x.parse().ok())
        .collect();
    if v.len() == 9 {
        Some(Matrix3::new(v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8]))
    } else {
        None
    }
}

fn parse_properties(spec: &str) -> Vec<(String, char, usize)> {
    let parts: Vec<&str> = spec.split(':').collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i + 2 < parts.len() {
        let name = parts[i].to_lowercase();
        let tc = parts[i+1].chars().next().unwrap_or('S').to_ascii_uppercase();
        let count: usize = parts[i+2].parse().unwrap_or(1);
        out.push((name, tc, count));
        i += 3;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const SINGLE: &str = r#"2
Lattice="2.87 0.0 0.0 0.0 2.87 0.0 0.0 0.0 2.87" Properties=species:S:1:pos:R:3 pbc="T T T" energy=-17.43
Fe 0.0    0.0    0.0
Fe 1.435  1.435  1.435
"#;

    const WITH_FORCES: &str = r#"2
Lattice="5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0" Properties=species:S:1:pos:R:3:forces:R:3 pbc="T T T"
Fe 0.0 0.0 0.0 0.1 0.0 0.0
Fe 2.5 2.5 2.5 -0.1 0.0 0.0
"#;

    fn tmp(name: &str, c: &str) -> String {
        let p = std::env::temp_dir().join(name);
        std::fs::write(&p, c).unwrap();
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn test_basic() {
        let traj = read_extxyz(&tmp("bcc.extxyz", SINGLE)).unwrap();
        assert_eq!(traj.n_frames(), 1);
        let f = traj.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
        assert!((f.energy.unwrap() - (-17.43)).abs() < 1e-10);
        let [a, ..] = f.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-6);
    }

    #[test]
    fn test_forces() {
        let traj = read_extxyz(&tmp("forces.extxyz", WITH_FORCES)).unwrap();
        let f = traj.first().unwrap();
        let forces = f.forces.as_ref().unwrap();
        assert!((forces[0].x - 0.1).abs() < 1e-10);
        assert!((forces[1].x - (-0.1)).abs() < 1e-10);
    }
}
