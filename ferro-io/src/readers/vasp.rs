use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{ensure, Context, Result};

pub fn read_poscar(path: &str) -> Result<Trajectory> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path}"))?;
    parse_poscar(&content).with_context(|| format!("parsing {path}"))
}

pub fn read_contcar(path: &str) -> Result<Trajectory> {
    read_poscar(path)
}

fn parse_poscar(content: &str) -> Result<Trajectory> {
    let mut lines = content.lines();
    let mut next = |what: &str| -> Result<&str> {
        lines.next().with_context(|| format!("unexpected EOF before {what}"))
    };

    let comment = next("comment")?.trim().to_string();

    let scale: f64 = next("scaling factor")?
        .trim().parse().context("invalid scaling factor")?;
    ensure!(scale > 0.0, "negative scaling factor (volume-based) is not supported");

    // Lattice vectors — rows of the Cell matrix
    let mut m = [[0.0_f64; 3]; 3];
    for (i, row) in m.iter_mut().enumerate() {
        let v = floats(next(&format!("lattice vector {i}"))?, 3)
            .with_context(|| format!("invalid lattice vector {i}"))?;
        *row = [v[0] * scale, v[1] * scale, v[2] * scale];
    }
    let cell = Cell::from_matrix(Matrix3::new(
        m[0][0], m[0][1], m[0][2],
        m[1][0], m[1][1], m[1][2],
        m[2][0], m[2][1], m[2][2],
    ));

    // VASP4 vs VASP5 detection
    let line5 = next("element/count line")?.trim();
    let (elements, counts): (Vec<String>, Vec<usize>) =
        if line5.split_whitespace().all(|t| t.parse::<u64>().is_ok()) {
            // VASP4: counts only
            let c: Vec<usize> = line5.split_whitespace().map(|s| s.parse().unwrap()).collect();
            let e = (1..=c.len()).map(|i| format!("X{i}")).collect();
            (e, c)
        } else {
            // VASP5: element symbols + counts on next line
            let e: Vec<String> = line5.split_whitespace().map(|s| s.to_string()).collect();
            let c: Vec<usize> = next("atom counts")?
                .split_whitespace()
                .map(|s| s.parse().context("invalid count"))
                .collect::<Result<_>>()?;
            ensure!(e.len() == c.len(), "element/count length mismatch");
            (e, c)
        };

    // Optional "Selective dynamics"
    let mut coord_type = next("coordinate type")?;
    if coord_type.trim().to_lowercase().starts_with('s') {
        coord_type = next("coordinate type after Selective dynamics")?;
    }
    let is_direct = coord_type.trim().to_lowercase().starts_with('d');

    // Atom positions
    let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
    let total: usize = counts.iter().sum();
    for (elem, &count) in elements.iter().zip(counts.iter()) {
        for _ in 0..count {
            let v = floats(next("coordinate")?, 3).context("invalid coordinate")?;
            let pos = if is_direct {
                cell.fractional_to_cartesian(Vector3::new(v[0], v[1], v[2]))
            } else {
                Vector3::new(v[0] * scale, v[1] * scale, v[2] * scale)
            };
            frame.add_atom(Atom::new(elem.as_str(), pos));
        }
    }

    // Optional velocity block (same count, no T/F flags)
    let rem: Vec<&str> = lines
        .map(|l| l.trim())
        .filter(|l| !l.is_empty())
        .collect();
    if rem.len() >= total {
        let vels: Vec<Vector3<f64>> = rem.iter().take(total)
            .filter_map(|l| floats(l, 3).ok())
            .map(|v| Vector3::new(v[0], v[1], v[2]))
            .collect();
        if vels.len() == total {
            frame.velocities = Some(vels);
        }
    }

    let mut traj = Trajectory::from_frame(frame);
    if !comment.is_empty() { traj.metadata.source = Some(comment); }
    Ok(traj)
}

fn floats(line: &str, min: usize) -> Result<Vec<f64>> {
    let v: Vec<f64> = line.split_whitespace()
        .map_while(|s| s.parse::<f64>().ok())
        .collect();
    ensure!(v.len() >= min, "expected ≥{min} floats, got {}", v.len());
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    const BCC_FE: &str = "BCC Fe
  1.00000000
     2.87000000   0.00000000   0.00000000
     0.00000000   2.87000000   0.00000000
     0.00000000   0.00000000   2.87000000
   Fe
   2
Direct
  0.0000000  0.0000000  0.0000000
  0.5000000  0.5000000  0.5000000
";

    const SCALE_POSCAR: &str = "Scaled
  2.00000000
     1.435 0.0 0.0
     0.0   1.435 0.0
     0.0   0.0   1.435
Fe
1
Direct
  0.0 0.0 0.0
";

    fn tmp(name: &str, content: &str) -> String {
        let p = std::env::temp_dir().join(name);
        std::fs::write(&p, content).unwrap();
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn test_bcc_fe() {
        let traj = read_poscar(&tmp("bcc.poscar", BCC_FE)).unwrap();
        let frame = traj.first().unwrap();
        assert_eq!(frame.n_atoms(), 2);
        assert_eq!(frame.atom(0).element, "Fe");
        let [a, ..] = frame.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-6);
    }

    #[test]
    fn test_scaling_factor() {
        let traj = read_poscar(&tmp("scale.poscar", SCALE_POSCAR)).unwrap();
        let [a, ..] = traj.first().unwrap().cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-6);
    }
}
