use std::path::Path;
use ferro_core::{Atom, Frame, Trajectory};
use nalgebra::Vector3;
use anyhow::{bail, Context, Result};
use ferro_core::data::elements::by_number;

/// 读取 XYZ 文件，支持多帧（多个连续 block）。
pub fn read_xyz(path: &Path) -> Result<Trajectory> {
    Ok(read_xyz_with_warnings(path)?.0)
}

/// [`read_xyz`], returning what the reader noticed: an incomplete last frame is
/// dropped and named, as the LAMMPS dump reader does.
pub fn read_xyz_with_warnings(path: &Path) -> Result<(Trajectory, Vec<String>)> {
    let mut lines = super::util::text_lines(path)?.enumerate().map(|(i, l)| (i + 1, l)).peekable();
    let mut traj = Trajectory::new();
    let mut warnings = Vec::new();
    // A file ending inside its last frame is an MD run cut off: drop that frame, keep the others
    // (as the dump reader does)
    let truncated = |k: usize, start: usize, why: String| {
        format!("warning: last frame (frame {k}, line {start}) is incomplete and was dropped: {why}")
    };

    'frames: loop {
        // 寻找下一个帧的原子数行，跳过空行
        let (count_no, count_line) = loop {
            match lines.next() {
                None => break 'frames,
                Some((no, l)) => {
                    let trimmed = l.context("read error")?.trim().to_string();
                    if !trimmed.is_empty() {
                        break (no, trimmed);
                    }
                }
            }
        };

        let n_atoms: usize = count_line
            .parse()
            .context(format!("expected atom count, got: {count_line:?}"))?;

        let k = traj.n_frames();
        let Some((_, comment)) = lines.next() else {
            warnings.push(truncated(k, count_no, "no comment line".into()));
            break;
        };
        let comment = comment.context("read error")?;

        let mut frame = Frame::new();
        for i in 0..n_atoms {
            let Some((line_no, line)) = lines.next() else {
                warnings.push(truncated(k, count_no, format!("{i} of {n_atoms} atoms")));
                break 'frames;
            };
            let line = line.context("read error")?;
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 4 {
                // The file's last line written only in part is a cut-off too
                if lines.peek().is_none() {
                    warnings.push(truncated(k, count_no, format!("{i} of {n_atoms} atoms")));
                    break 'frames;
                }
                bail!("line {line_no}: expected element + 3 coords, got {line:?}");
            }
            let x: f64 = parts[1].parse().context("invalid x")?;
            let y: f64 = parts[2].parse().context("invalid y")?;
            let z: f64 = parts[3].parse().context("invalid z")?;
            // An atomic number in the first column becomes its symbol (as in ASE); `8` used to be the
            // element name
            let element = match parts[0].parse::<u8>() {
                Ok(z) => by_number(z)
                    .with_context(|| format!("line {line_no}: no element has atomic number {z}"))?
                    .symbol.to_string(),
                Err(_) => parts[0].to_string(),
            };
            frame.add_atom(Atom::new(element, Vector3::new(x, y, z)));
        }

        if traj.n_frames() == 0 {
            let s = comment.trim();
            if !s.is_empty() {
                traj.metadata.title = Some(s.to_string());
            }
        }

        traj.add_frame(frame);
    }
    Ok((traj, warnings))
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

    use crate::testutil::write_tmp;

    #[test]
    fn test_single_frame() {
        let path = write_tmp("test_water.xyz", WATER_XYZ);
        let traj = read_xyz(&path).unwrap();
        assert_eq!(traj.n_frames(), 1);
        let frame = traj.first().unwrap();
        assert_eq!(frame.n_atoms(), 3);
        assert_eq!(frame.atom(0).element, "O");
        assert_eq!(traj.metadata.title.as_deref(), Some("Water molecule"));
    }

    #[test]
    fn test_multi_frame() {
        let path = write_tmp("test_multi.xyz", MULTI_XYZ);
        let traj = read_xyz(&path).unwrap();
        assert_eq!(traj.n_frames(), 2);
        assert_eq!(traj.frame(0).unwrap().n_atoms(), 2);
        assert_eq!(traj.frame(1).unwrap().n_atoms(), 2);
    }

    /// MD cut off mid-frame: the earlier frames stay, the cut one is dropped and named.
    #[test]
    fn test_truncated_last_frame_is_dropped() {
        for (name, tail) in [
            ("trunc_atoms.xyz", "2\nframe 1\nC 0 0 0.1\n"),       // 1 of 2 atom lines
            ("trunc_half.xyz", "2\nframe 1\nC 0 0 0.1\nC 1.4"),   // last line cut
            ("trunc_comment.xyz", "2\n"),                         // only the count line
        ] {
            let path = write_tmp(name, &format!("2\nframe 0\nC 0 0 0\nC 1.4 0 0\n{tail}"));
            let (traj, warnings) = read_xyz_with_warnings(&path).unwrap();
            assert_eq!(traj.n_frames(), 1, "{name}");
            assert_eq!(warnings.len(), 1, "{name}");
            assert!(warnings[0].contains("frame 1, line 5"), "{name}: {}", warnings[0]);
        }
        // An incomplete middle frame is still an error
        let path = write_tmp("trunc_middle.xyz", "2\nf0\nC 0 0 0\n2\nf1\nC 0 0 0\nC 1 0 0\n");
        assert!(read_xyz(&path).is_err());
    }

    /// An atomic number in the first column is an element, as ASE reads it.
    #[test]
    fn test_atomic_number_column() {
        let path = write_tmp("znum.xyz", "2\n\n8 0 0 0\n1 0.96 0 0\n");
        let frame = read_xyz(&path).unwrap().frames.remove(0);
        assert_eq!(frame.atom(0).element, "O");
        assert_eq!(frame.atom(1).element, "H");
        let path = write_tmp("znum_bad.xyz", "1\n\n200 0 0 0\n");
        assert!(read_xyz(&path).unwrap_err().to_string().contains("atomic number 200"));
    }
}
