use std::collections::HashMap;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{bail, Context, Result};

const BOHR: f64 = 0.52917721;

pub fn read_qe_input(path: &str) -> Result<Trajectory> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path}"))?;
    parse_qe(&content).with_context(|| format!("parsing {path}"))
}

fn parse_qe(content: &str) -> Result<Trajectory> {
    let lines: Vec<&str> = content.lines().collect();

    // Strip Fortran comments (lines starting with !)
    let stripped: Vec<&str> = lines.iter()
        .map(|l| {
            let t = l.trim();
            if let Some(p) = t.find('!') { &t[..p] } else { t }
        })
        .collect();

    // ── &SYSTEM namelist ──────────────────────────────────────────────────────
    let sys = collect_namelist(&stripped, "system");
    let ibrav: i32 = sys.get("ibrav").and_then(|s| s.parse().ok()).unwrap_or(0);
    if ibrav != 0 {
        bail!("ibrav={ibrav} is not supported; use ibrav=0 with explicit CELL_PARAMETERS");
    }
    let _ntyp: usize = sys.get("ntyp").and_then(|s| s.parse().ok()).unwrap_or(0);

    // starting_magnetization(n) → type_index → magmom
    let mut type_magmom: HashMap<usize, f64> = HashMap::new();
    for (k, v) in &sys {
        if k.starts_with("starting_magnetization(") {
            let idx: usize = k.trim_start_matches("starting_magnetization(")
                .trim_end_matches(')')
                .parse().unwrap_or(0);
            if let Ok(m) = v.parse::<f64>() { type_magmom.insert(idx, m); }
        }
    }

    // ── ATOMIC_SPECIES card ───────────────────────────────────────────────────
    // label → element (first line word treated as label, second as mass (ignored), third as pseudo)
    let mut species: HashMap<String, (String, usize)> = HashMap::new(); // label → (element, type_index)
    if let Some(start) = find_card(&stripped, "ATOMIC_SPECIES") {
        let mut tidx = 1usize;
        for l in stripped[start+1..].iter() {
            let l = l.trim();
            if l.is_empty() || is_card_or_namelist(l) { break; }
            let parts: Vec<&str> = l.split_whitespace().collect();
            if parts.len() >= 2 {
                let label = parts[0].to_string();
                // Try to extract element: first alphabetic chars of label
                let elem = extract_element(&label);
                species.insert(label, (elem, tidx));
                tidx += 1;
            }
        }
    }

    // ── CELL_PARAMETERS card ──────────────────────────────────────────────────
    let (cell_start, cell_unit) = match find_card_with_unit(&stripped, "CELL_PARAMETERS") {
        Some(r) => r,
        None => bail!("CELL_PARAMETERS card not found"),
    };
    let cell_scale = match cell_unit.as_str() {
        "bohr" => BOHR,
        _ => 1.0, // angstrom (default)
    };

    let mut cell_vecs: [[f64; 3]; 3] = [[0.0; 3]; 3];
    let mut count = 0;
    for l in stripped[cell_start+1..].iter() {
        if count >= 3 { break; }
        let l = l.trim();
        if l.is_empty() || is_card_or_namelist(l) { break; }
        let parts: Vec<f64> = l.split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();
        if parts.len() >= 3 {
            cell_vecs[count] = [parts[0]*cell_scale, parts[1]*cell_scale, parts[2]*cell_scale];
            count += 1;
        }
    }
    anyhow::ensure!(count == 3, "CELL_PARAMETERS must have 3 vectors");
    let cell = Cell::from_matrix(Matrix3::new(
        cell_vecs[0][0], cell_vecs[0][1], cell_vecs[0][2],
        cell_vecs[1][0], cell_vecs[1][1], cell_vecs[1][2],
        cell_vecs[2][0], cell_vecs[2][1], cell_vecs[2][2],
    ));

    // ── ATOMIC_POSITIONS card ─────────────────────────────────────────────────
    let (pos_start, pos_unit) = match find_card_with_unit(&stripped, "ATOMIC_POSITIONS") {
        Some(r) => r,
        None => bail!("ATOMIC_POSITIONS card not found"),
    };

    anyhow::ensure!(
        !matches!(pos_unit.as_str(), "alat"),
        "ATOMIC_POSITIONS alat units are not supported; use angstrom, bohr, or crystal"
    );

    let mut frame = Frame::with_cell(cell.clone(), [true; 3]);

    for l in stripped[pos_start+1..].iter() {
        let l = l.trim();
        if l.is_empty() || is_card_or_namelist(l) { break; }
        let parts: Vec<&str> = l.split_whitespace().collect();
        if parts.len() < 4 { continue; }

        let label = parts[0].to_string();
        let (elem, tidx) = species.get(&label)
            .cloned()
            .unwrap_or_else(|| (extract_element(&label), 0));

        let x: f64 = parts[1].parse().unwrap_or(0.0);
        let y: f64 = parts[2].parse().unwrap_or(0.0);
        let z: f64 = parts[3].parse().unwrap_or(0.0);

        let position = match pos_unit.as_str() {
            "bohr" => Vector3::new(x * BOHR, y * BOHR, z * BOHR),
            "crystal" => cell.fractional_to_cartesian(Vector3::new(x, y, z)),
            _ => Vector3::new(x, y, z), // angstrom
        };

        let mut atom = Atom::new(elem, position);
        atom.label = Some(label);
        if tidx > 0 { atom.magmom = type_magmom.get(&tidx).copied(); }
        frame.add_atom(atom);
    }

    Ok(Trajectory::from_frame(frame))
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

fn collect_namelist(lines: &[&str], name: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    let target = format!("&{name}");
    let start = match lines.iter().position(|l| l.to_lowercase().trim() == target) {
        Some(p) => p + 1,
        None => return map,
    };
    for l in &lines[start..] {
        let l = l.trim();
        if l == "/" || l.starts_with('&') { break; }
        // Parse "key = value," pairs (may be multiple per line)
        for pair in l.split(',') {
            let pair = pair.trim();
            if let Some(eq) = pair.find('=') {
                let key = pair[..eq].trim().to_lowercase()
                    .replace(' ', "");
                let val = pair[eq+1..].trim()
                    .trim_matches('\'')
                    .trim_matches('"')
                    .to_string();
                if !key.is_empty() { map.insert(key, val); }
            }
        }
    }
    map
}

fn find_card(lines: &[&str], card: &str) -> Option<usize> {
    lines.iter().position(|l| {
        l.to_uppercase().split_whitespace().next() == Some(card)
    })
}

fn find_card_with_unit(lines: &[&str], card: &str) -> Option<(usize, String)> {
    lines.iter().position(|l| {
        l.to_uppercase().split_whitespace().next() == Some(card)
    }).map(|pos| {
        let l = lines[pos].to_lowercase();
        let unit = l.find('{')
            .and_then(|s| l.find('}').map(|e| l[s+1..e].trim().to_string()))
            .unwrap_or_else(|| "angstrom".to_string());
        (pos, unit)
    })
}

fn is_card_or_namelist(l: &str) -> bool {
    let up = l.to_uppercase();
    matches!(
        up.split_whitespace().next().unwrap_or(""),
        "ATOMIC_POSITIONS" | "ATOMIC_SPECIES" | "CELL_PARAMETERS" |
        "K_POINTS" | "CONSTRAINTS" | "ATOMIC_FORCES"
    ) || up.trim().starts_with('&')
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

    const WATER_QE: &str = "
&CONTROL
  calculation = 'scf',
  prefix = 'water',
/
&SYSTEM
  ibrav = 0,
  nat = 3,
  ntyp = 2,
  ecutwfc = 30.0,
/
&ELECTRONS
/
ATOMIC_SPECIES
O  15.999  O.pbesol.UPF
H  1.008   H.pbesol.UPF

ATOMIC_POSITIONS {angstrom}
O  0.000  0.000  0.119
H  0.000  0.763 -0.477
H  0.000 -0.763 -0.477

CELL_PARAMETERS {angstrom}
10.0  0.0  0.0
0.0  10.0  0.0
0.0   0.0 10.0

K_POINTS {gamma}
";

    const CRYSTAL_POS: &str = "
&SYSTEM
  ibrav = 0,
  nat = 2,
  ntyp = 1,
/
ATOMIC_SPECIES
Fe  55.845  Fe.UPF

ATOMIC_POSITIONS {crystal}
Fe  0.0  0.0  0.0
Fe  0.5  0.5  0.5

CELL_PARAMETERS {angstrom}
2.87  0.0   0.0
0.0   2.87  0.0
0.0   0.0   2.87
";

    fn tmp(n: &str, c: &str) -> String {
        let p = std::env::temp_dir().join(n);
        std::fs::write(&p, c).unwrap();
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn test_water() {
        let traj = read_qe_input(&tmp("water.qe", WATER_QE)).unwrap();
        let f = traj.first().unwrap();
        assert_eq!(f.n_atoms(), 3);
        assert_eq!(f.atom(0).element, "O");
        assert_eq!(f.atom(1).element, "H");
    }

    #[test]
    fn test_crystal_coords() {
        let traj = read_qe_input(&tmp("bcc.qe", CRYSTAL_POS)).unwrap();
        let f = traj.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        // Fe at (0,0,0) should be at Cartesian (0,0,0)
        assert!(f.atom(0).position.norm() < 1e-6);
        // Fe at (0.5,0.5,0.5) fractional = (1.435,1.435,1.435) Cartesian
        assert!((f.atom(1).position.x - 1.435).abs() < 1e-3);
    }
}
