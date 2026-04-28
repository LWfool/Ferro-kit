use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::Vector3;
use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::{Context, Result};

/// 读取 PDB 文件，支持多模型（MODEL/ENDMDL 记录）。
pub fn read_pdb(path: &str) -> Result<Trajectory> {
    let file = File::open(path).context(format!("cannot open {path}"))?;
    let reader = BufReader::new(file);

    let mut traj = Trajectory::new();
    let mut current = Frame::new();
    let mut has_model = false;
    let mut cell: Option<Cell> = None;

    for line in reader.lines() {
        let line = line.context("read error")?;
        let tag = &line[..line.len().min(6)];

        match tag {
            "HEADER" => {
                let source = line[line.len().min(10)..].trim();
                if !source.is_empty() {
                    traj.metadata.source = Some(source.to_string());
                }
            }
            "CRYST1" => {
                cell = parse_cryst1(&line);
            }
            "MODEL " => {
                has_model = true;
            }
            "ENDMDL" => {
                let mut frame = std::mem::replace(&mut current, Frame::new());
                if let Some(c) = &cell {
                    frame.cell = Some(c.clone());
                    frame.pbc = [true; 3];
                }
                traj.add_frame(frame);
            }
            "ATOM  " | "HETATM" => {
                if let Some(atom) = parse_atom_record(&line) {
                    current.add_atom(atom);
                }
            }
            _ => {}
        }
    }

    // 没有 MODEL 记录时，视为单帧
    if !has_model && current.n_atoms() > 0 {
        if let Some(c) = cell {
            current.cell = Some(c);
            current.pbc = [true; 3];
        }
        traj.add_frame(current);
    }

    Ok(traj)
}

/// CRYST1 format: CRYST1  a  b  c  α  β  γ  sgroup  z
fn parse_cryst1(line: &str) -> Option<Cell> {
    if line.len() < 54 { return None; }
    let a:  f64 = line[6..15].trim().parse().ok()?;
    let b:  f64 = line[15..24].trim().parse().ok()?;
    let c:  f64 = line[24..33].trim().parse().ok()?;
    let al: f64 = line[33..40].trim().parse().ok()?;
    let be: f64 = line[40..47].trim().parse().ok()?;
    let ga: f64 = line[47..54].trim().parse().ok()?;
    Cell::from_lengths_angles(a, b, c, al, be, ga).ok()
}

fn parse_atom_record(line: &str) -> Option<Atom> {
    if line.len() < 54 {
        return None;
    }
    let element = if line.len() >= 78 {
        line[76..78].trim().to_string()
    } else {
        // 从原子名称字段推断（取第一个非空字符）
        line[12..16]
            .trim()
            .chars()
            .next()
            .map(|c| c.to_string())
            .unwrap_or_else(|| "X".to_string())
    };
    let x: f64 = line[30..38].trim().parse().ok()?;
    let y: f64 = line[38..46].trim().parse().ok()?;
    let z: f64 = line[46..54].trim().parse().ok()?;
    Some(Atom::new(element, Vector3::new(x, y, z)))
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER_PDB: &str = "\
HEADER    Water molecule
ATOM      1  O   UNK A   1       0.000   0.000   0.119  1.00  0.00           O
ATOM      2  H   UNK A   1       0.000   0.763  -0.477  1.00  0.00           H
ATOM      3  H   UNK A   1       0.000  -0.763  -0.477  1.00  0.00           H
END
";

    const MULTI_PDB: &str = "\
MODEL        1
ATOM      1  C   UNK A   1       0.000   0.000   0.000  1.00  0.00           C
ENDMDL
MODEL        2
ATOM      1  C   UNK A   1       1.400   0.000   0.000  1.00  0.00           C
ENDMDL
";

    fn write_tmp(name: &str, content: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(name);
        std::fs::write(&path, content).unwrap();
        path
    }

    #[test]
    fn test_single_frame() {
        let path = write_tmp("test_water.pdb", WATER_PDB);
        let traj = read_pdb(path.to_str().unwrap()).unwrap();
        assert_eq!(traj.n_frames(), 1);
        let frame = traj.first().unwrap();
        assert_eq!(frame.n_atoms(), 3);
        assert_eq!(frame.atom(0).element, "O");
        assert_eq!(traj.metadata.source.as_deref(), Some("Water molecule"));
    }

    #[test]
    fn test_multi_model() {
        let path = write_tmp("test_multi.pdb", MULTI_PDB);
        let traj = read_pdb(path.to_str().unwrap()).unwrap();
        assert_eq!(traj.n_frames(), 2);
    }
}
