use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};

/// 将轨迹写入 PDB 文件。多帧使用 MODEL/ENDMDL 记录。
pub fn write_pdb(trajectory: &Trajectory, path: &str) -> Result<()> {
    let file = File::create(path).context(format!("cannot create {path}"))?;
    let mut writer = BufWriter::new(file);

    if let Some(source) = &trajectory.metadata.source {
        writeln!(writer, "HEADER    {source}")?;
    }

    let multi = trajectory.n_frames() > 1;

    // CRYST1 — write cell from first frame that has one
    if let Some(cell) = trajectory.frames.iter().find_map(|f| f.cell.as_ref()) {
        let [a, b, c] = cell.lengths();
        let [al, be, ga] = cell.angles();
        writeln!(writer, "CRYST1{:>9.3}{:>9.3}{:>9.3}{:>7.2}{:>7.2}{:>7.2} P 1           1", a, b, c, al, be, ga)?;
    }

    for (model_idx, frame) in trajectory.frames.iter().enumerate() {
        if multi {
            writeln!(writer, "MODEL     {:>4}", model_idx + 1)?;
        }
        for (i, atom) in frame.atoms.iter().enumerate() {
            writeln!(
                writer,
                "{:<6}{:>5} {:^4} {:3} {:1}{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
                "ATOM",
                i + 1,
                atom.element,
                "UNK",
                "A",
                1,
                atom.position.x,
                atom.position.y,
                atom.position.z,
                1.00,
                0.00,
                atom.element,
            )?;
        }
        if multi {
            writeln!(writer, "ENDMDL")?;
        }
    }

    writeln!(writer, "END")?;
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::pdb::read_pdb;
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
        let path = std::env::temp_dir().join("roundtrip_water.pdb");
        let path_str = path.to_str().unwrap();
        let original = make_traj();
        write_pdb(&original, path_str).unwrap();

        let loaded = read_pdb(path_str).unwrap();
        assert_eq!(loaded.n_frames(), 1);
        assert_eq!(loaded.first().unwrap().n_atoms(), 3);
    }
}
