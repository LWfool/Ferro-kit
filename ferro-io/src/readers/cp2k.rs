use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{bail, Context, Result};

pub fn read_cp2k_inp(path: &str) -> Result<Trajectory> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path}"))?;
    parse_cp2k(&content).with_context(|| format!("parsing {path}"))
}

/// CP2K restart 与 inp 格式相同，共用同一解析器。
pub fn read_cp2k_restart(path: &str) -> Result<Trajectory> {
    read_cp2k_inp(path)
}

const BOHR: f64 = 0.52917721; // Bohr → Å

fn parse_cp2k(content: &str) -> Result<Trajectory> {
    let lines: Vec<&str> = content.lines().collect();

    // Build a stripped-comment line iterator
    let stripped: Vec<&str> = lines.iter()
        .map(|l| {
            let l = l.trim();
            // CP2K comments start with ! or #
            if let Some(p) = l.find('!').or_else(|| l.find('#')) { &l[..p] } else { l }
        })
        .collect();

    // ── Find &SUBSYS ──────────────────────────────────────────────────────────
    let subsys_start = stripped.iter().position(|l| {
        let low = l.to_lowercase();
        low.starts_with("&subsys")
    }).context("&SUBSYS section not found")?;

    let subsys_end = find_end(&stripped, subsys_start)
        .context("&END SUBSYS not found")?;

    let subsys: &[&str] = &stripped[subsys_start..=subsys_end];

    // ── &KIND label → element map ─────────────────────────────────────────────
    let mut kind_map: std::collections::HashMap<String, String> = Default::default();
    {
        let mut j = 0;
        while j < subsys.len() {
            let low = subsys[j].to_lowercase();
            if low.starts_with("&kind") {
                let label = low["&kind".len()..].trim().to_string();
                let end = find_end(subsys, j).unwrap_or(subsys.len() - 1);
                for l in &subsys[j..=end] {
                    let kv: Vec<&str> = l.split_whitespace().collect();
                    if kv.first().map(|s| s.to_lowercase()) == Some("element".to_string()) {
                        if let Some(elem) = kv.get(1) {
                            kind_map.insert(label.clone(), elem.to_string());
                        }
                    }
                }
            }
            j += 1;
        }
    }

    // ── &CELL ─────────────────────────────────────────────────────────────────
    let cell_start = subsys.iter().position(|l| l.to_lowercase().starts_with("&cell"))
        .context("&CELL not found")?;
    let cell_end = find_end(subsys, cell_start).context("&END CELL not found")?;
    let cell_section: &[&str] = &subsys[cell_start..=cell_end];

    let cell = parse_cell_section(cell_section)?;

    // ── &COORD ────────────────────────────────────────────────────────────────
    let coord_start = subsys.iter().position(|l| l.to_lowercase().starts_with("&coord"))
        .context("&COORD not found")?;
    let coord_end = find_end(subsys, coord_start).context("&END COORD not found")?;
    let coord_section: &[&str] = &subsys[coord_start..=coord_end];

    let mut unit_bohr = false;
    let mut scaled = false;
    for l in coord_section {
        let low = l.to_lowercase();
        if low.starts_with("unit") && low.contains("bohr") { unit_bohr = true; }
        if low.starts_with("scaled") && (low.contains("true") || low.contains(".true.")) {
            scaled = true;
        }
    }

    let mut frame = Frame::with_cell(cell.clone(), [true; 3]);

    for l in &coord_section[1..coord_section.len().saturating_sub(1)] {
        let l = l.trim();
        if l.is_empty() || l.to_lowercase().starts_with("unit")
            || l.to_lowercase().starts_with("scaled") { continue; }

        let parts: Vec<&str> = l.split_whitespace().collect();
        if parts.len() < 4 { continue; }

        let label = parts[0].to_string();
        let elem = kind_map.get(&label.to_lowercase())
            .cloned()
            .unwrap_or_else(|| extract_element(&label));

        let x: f64 = parts[1].parse().unwrap_or(0.0);
        let y: f64 = parts[2].parse().unwrap_or(0.0);
        let z: f64 = parts[3].parse().unwrap_or(0.0);

        let frac = Vector3::new(x, y, z);
        let cart_raw = Vector3::new(
            if unit_bohr { x * BOHR } else { x },
            if unit_bohr { y * BOHR } else { y },
            if unit_bohr { z * BOHR } else { z },
        );

        let position = if scaled {
            cell.fractional_to_cartesian(frac)
        } else {
            cart_raw
        };

        let mut atom = Atom::new(elem, position);
        atom.label = Some(label);

        // Optional velocity (columns 5-7 in restart)
        if parts.len() >= 7 {
            if let (Ok(vx), Ok(vy), Ok(vz)) = (
                parts[4].parse::<f64>(),
                parts[5].parse::<f64>(),
                parts[6].parse::<f64>(),
            ) {
                // CP2K velocity units: Bohr/fs → Å/fs
                if frame.velocities.is_none() {
                    frame.velocities = Some(Vec::new());
                }
                frame.velocities.as_mut().unwrap()
                    .push(Vector3::new(vx * BOHR, vy * BOHR, vz * BOHR));
            }
        }

        frame.add_atom(atom);
    }

    // Sanity check velocities length
    if let Some(vels) = &frame.velocities {
        if vels.len() != frame.n_atoms() {
            frame.velocities = None;
        }
    }

    Ok(Trajectory::from_frame(frame))
}

fn parse_cell_section(section: &[&str]) -> Result<Cell> {
    // Supports:
    //   A x y z / B x y z / C x y z  (explicit vectors)
    //   ABC a b c  (diagonal only, Angstrom)
    //   ALPHA_BETA_GAMMA α β γ (angles, paired with ABC)
    let mut vecs: [Option<[f64; 3]>; 3] = [None; 3];
    let mut abc: Option<[f64; 3]> = None;
    let mut angles: [f64; 3] = [90.0; 3];
    let mut unit_bohr = false;

    for l in section {
        let low = l.to_lowercase();
        let parts: Vec<&str> = l.split_whitespace().collect();
        if low.starts_with("unit") && low.contains("bohr") { unit_bohr = true; continue; }

        if low.starts_with("a ") || low.starts_with("a\t") {
            if let Some(v) = parse_vec3(&parts[1..]) { vecs[0] = Some(v); }
        } else if low.starts_with("b ") || low.starts_with("b\t") {
            if let Some(v) = parse_vec3(&parts[1..]) { vecs[1] = Some(v); }
        } else if low.starts_with("c ") || low.starts_with("c\t") {
            if let Some(v) = parse_vec3(&parts[1..]) { vecs[2] = Some(v); }
        } else if low.starts_with("abc ") {
            abc = parse_vec3(&parts[1..]).map(|v| v);
        } else if low.starts_with("alpha_beta_gamma") {
            if let Some(v) = parse_vec3(&parts[1..]) { angles = v; }
        }
    }

    let scale = if unit_bohr { BOHR } else { 1.0 };

    if let (Some(a), Some(b), Some(c)) = (vecs[0], vecs[1], vecs[2]) {
        Ok(Cell::from_matrix(Matrix3::new(
            a[0]*scale, a[1]*scale, a[2]*scale,
            b[0]*scale, b[1]*scale, b[2]*scale,
            c[0]*scale, c[1]*scale, c[2]*scale,
        )))
    } else if let Some([a, b, c]) = abc {
        Cell::from_lengths_angles(a*scale, b*scale, c*scale, angles[0], angles[1], angles[2])
            .map_err(|e| anyhow::anyhow!("{e}"))
    } else {
        bail!("cannot parse &CELL: need A/B/C vectors or ABC lengths")
    }
}

fn find_end(lines: &[&str], start: usize) -> Option<usize> {
    let _section_name: String = lines[start].to_lowercase()
        .trim_start_matches('&')
        .split_whitespace().next().unwrap_or("").to_string();
    let mut depth = 0usize;
    for (j, l) in lines[start..].iter().enumerate() {
        let low = l.to_lowercase().trim().to_string();
        if low.starts_with('&') && !low.starts_with("&end") {
            depth += 1;
        } else if low.starts_with("&end") {
            depth = depth.saturating_sub(1);
            if depth == 0 { return Some(start + j); }
        }
    }
    None
}

fn parse_vec3(parts: &[&str]) -> Option<[f64; 3]> {
    if parts.len() < 3 { return None; }
    Some([
        parts[0].parse().ok()?,
        parts[1].parse().ok()?,
        parts[2].parse().ok()?,
    ])
}

fn extract_element(label: &str) -> String {
    let alpha: String = label.chars().take_while(|c| c.is_alphabetic()).collect();
    if alpha.is_empty() { return "X".to_string(); }
    let mut it = alpha.chars();
    let first = it.next().unwrap().to_uppercase().to_string();
    let rest: String = it.collect::<String>().to_lowercase();
    format!("{first}{rest}")
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER_INP: &str = "
&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    &CELL
      ABC 10.0 10.0 10.0
    &END CELL
    &COORD
      O  0.000  0.000  0.119
      H  0.000  0.763 -0.477
      H  0.000 -0.763 -0.477
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
";

    const FE_INP: &str = "
&FORCE_EVAL
  &SUBSYS
    &CELL
      A 2.87 0.0  0.0
      B 0.0  2.87 0.0
      C 0.0  0.0  2.87
    &END CELL
    &COORD
      Fe1  0.0   0.0   0.0
      Fe2  1.435 1.435 1.435
    &END COORD
    &KIND Fe1
      ELEMENT Fe
    &END KIND
    &KIND Fe2
      ELEMENT Fe
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
";

    fn tmp(n: &str, c: &str) -> String {
        let p = std::env::temp_dir().join(n);
        std::fs::write(&p, c).unwrap();
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn test_water_abc() {
        let traj = read_cp2k_inp(&tmp("water.inp", WATER_INP)).unwrap();
        let f = traj.first().unwrap();
        assert_eq!(f.n_atoms(), 3);
        assert_eq!(f.atom(0).element, "O");
    }

    #[test]
    fn test_fe_bcc_explicit_vectors() {
        let traj = read_cp2k_inp(&tmp("fe.inp", FE_INP)).unwrap();
        let f = traj.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
        let [a, ..] = f.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-4);
    }
}
