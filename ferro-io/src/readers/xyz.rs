use ferro_core::{Atom, Frame, Trajectory};
use nalgebra::Vector3;
use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::{Context, Result};

/// 读取 XYZ 文件，支持多帧（多个连续 block）。
pub fn read_xyz(path: &str) -> Result<Trajectory> {
    let file = File::open(path).context(format!("cannot open {path}"))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let mut traj = Trajectory::new();

    loop {
        // 寻找下一个帧的原子数行，跳过空行
        let count_line = loop {
            match lines.next() {
                None => return Ok(traj),
                Some(l) => {
                    let l = l.context("read error")?;
                    let trimmed = l.trim().to_string();
                    if !trimmed.is_empty() {
                        break trimmed;
                    }
                }
            }
        };

        let n_atoms: usize = count_line
            .parse()
            .context(format!("expected atom count, got: {count_line:?}"))?;

        let comment = lines
            .next()
            .context("missing comment line")?
            .context("read error")?;

        let mut frame = Frame::new();
        for i in 0..n_atoms {
            let line = lines
                .next()
                .context(format!("missing atom line {}", i + 1))?
                .context("read error")?;
            let parts: Vec<&str> = line.split_whitespace().collect();
            anyhow::ensure!(
                parts.len() >= 4,
                "atom line {}: expected element + 3 coords, got {:?}",
                i + 1,
                line
            );
            let x: f64 = parts[1].parse().context("invalid x")?;
            let y: f64 = parts[2].parse().context("invalid y")?;
            let z: f64 = parts[3].parse().context("invalid z")?;
            frame.add_atom(Atom::new(parts[0], Vector3::new(x, y, z)));
        }

        if traj.n_frames() == 0 {
            let s = comment.trim();
            if !s.is_empty() {
                traj.metadata.source = Some(s.to_string());
            }
        }

        traj.add_frame(frame);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER_XYZ: &str = "\
3
Water molecule
O     0.000000    0.000000    0.119262
H     0.000000    0.763239   -0.477047
H     0.000000   -0.763239   -0.477047
";

    const MULTI_XYZ: &str = "\
2
frame 0
C  0.0  0.0  0.0
C  1.4  0.0  0.0

2
frame 1
C  0.0  0.0  0.1
C  1.4  0.0  0.1
";

    fn write_tmp(name: &str, content: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(name);
        std::fs::write(&path, content).unwrap();
        path
    }

    #[test]
    fn test_single_frame() {
        let path = write_tmp("test_water.xyz", WATER_XYZ);
        let traj = read_xyz(path.to_str().unwrap()).unwrap();
        assert_eq!(traj.n_frames(), 1);
        let frame = traj.first().unwrap();
        assert_eq!(frame.n_atoms(), 3);
        assert_eq!(frame.atom(0).element, "O");
        assert_eq!(traj.metadata.source.as_deref(), Some("Water molecule"));
    }

    #[test]
    fn test_multi_frame() {
        let path = write_tmp("test_multi.xyz", MULTI_XYZ);
        let traj = read_xyz(path.to_str().unwrap()).unwrap();
        assert_eq!(traj.n_frames(), 2);
        assert_eq!(traj.frame(0).unwrap().n_atoms(), 2);
        assert_eq!(traj.frame(1).unwrap().n_atoms(), 2);
    }
}
