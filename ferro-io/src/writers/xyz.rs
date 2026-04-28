use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};

/// 将轨迹写入 XYZ 文件。多帧轨迹写为连续的多个 XYZ block。
pub fn write_xyz(trajectory: &Trajectory, path: &str) -> Result<()> {
    let file = File::create(path).context(format!("cannot create {path}"))?;
    let mut writer = BufWriter::new(file);

    for (i, frame) in trajectory.frames.iter().enumerate() {
        writeln!(writer, "{}", frame.n_atoms())?;
        let comment = if i == 0 {
            trajectory.metadata.source.as_deref().unwrap_or("")
        } else {
            ""
        };
        writeln!(writer, "{comment}")?;
        for atom in &frame.atoms {
            writeln!(
                writer,
                "{:<2}  {:>12.6}  {:>12.6}  {:>12.6}",
                atom.element, atom.position.x, atom.position.y, atom.position.z
            )?;
        }
    }

    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::xyz::read_xyz;
    use ferro_core::{Atom, Frame};
    use nalgebra::Vector3;

    fn make_traj() -> Trajectory {
        let mut frame = Frame::new();
        frame.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.119)));
        frame.add_atom(Atom::new("H", Vector3::new(0.0, 0.763, -0.477)));
        frame.add_atom(Atom::new("H", Vector3::new(0.0, -0.763, -0.477)));
        let mut traj = Trajectory::from_frame(frame);
        traj.metadata.source = Some("Water molecule".to_string());
        traj
    }

    #[test]
    fn test_roundtrip() {
        let path = std::env::temp_dir().join("roundtrip_water.xyz");
        let path_str = path.to_str().unwrap();
        let original = make_traj();
        write_xyz(&original, path_str).unwrap();

        let loaded = read_xyz(path_str).unwrap();
        assert_eq!(loaded.n_frames(), 1);
        assert_eq!(loaded.first().unwrap().n_atoms(), 3);
        assert_eq!(loaded.metadata.source.as_deref(), Some("Water molecule"));
    }
}
