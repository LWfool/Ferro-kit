use std::path::Path;
use ferro_core::{Cell, Trajectory};
use nalgebra::Matrix3;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};

/// 将轨迹写入 PDB 文件。多帧使用 MODEL/ENDMDL 记录。
pub fn write_pdb(trajectory: &Trajectory, path: &Path) -> Result<()> {
    let path_ = path.display();
    let file = File::create(path).context(format!("cannot create {path_}"))?;
    let mut writer = BufWriter::new(file);

    if let Some(source) = &trajectory.metadata.source {
        writeln!(writer, "HEADER    {source}")?;
    }

    let multi = trajectory.n_frames() > 1;

    for (model_idx, frame) in trajectory.frames.iter().enumerate() {
        // 逐帧写 CRYST1（NPT 每帧胞不同），放在 MODEL 之前，与 ASE 同序。
        // CRYST1 只存 a b c α β γ，reader 按 a 沿 x、b 在 xy 平面重建胞（同 ASE
        // standard_form），故坐标须经分数坐标映到该取向，否则非标准取向的胞读回后错位
        let mut to_std = Matrix3::identity();
        if let Some(cell) = &frame.cell {
            let [a, b, c] = cell.lengths();
            let [al, be, ga] = cell.angles();
            let m = &cell.matrix;
            // 已是标准取向（LAMMPS 胞恒是）就不变换：M⁻¹·S 的 1e-16 噪声会把恰在
            // x.xxx5 的坐标推过 3 位小数的舍入点
            let standard = m[(0, 1)] == 0.0 && m[(0, 2)] == 0.0 && m[(1, 2)] == 0.0
                && m[(0, 0)] > 0.0 && m[(1, 1)] > 0.0 && m[(2, 2)] > 0.0;
            if !standard {
                let frame_no = model_idx + 1;
                let std_cell = Cell::from_lengths_angles(a, b, c, al, be, ga)
                    .context(format!("frame {frame_no}: cell cannot be written as CRYST1"))?;
                let inv = m.try_inverse()
                    .with_context(|| format!("frame {frame_no}: cell is singular"))?;
                // 行向量约定：p_std = p · M⁻¹ · S，列向量下即 (M⁻¹ S)ᵀ p
                to_std = (inv * std_cell.matrix).transpose();
            }
            writeln!(writer, "CRYST1{:>9.3}{:>9.3}{:>9.3}{:>7.2}{:>7.2}{:>7.2} P 1           1", a, b, c, al, be, ga)?;
        }
        if multi {
            writeln!(writer, "MODEL     {:>4}", model_idx + 1)?;
        }
        for (i, atom) in frame.atoms.iter().enumerate() {
            let p = to_std * atom.position;
            writeln!(
                writer,
                "{:<6}{:>5} {:^4} {:3} {:1}{:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
                "ATOM",
                i + 1,
                atom.element,
                "UNK",
                "A",
                1,
                p.x,
                p.y,
                p.z,
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
        let path_str = &path;
        let original = make_traj();
        write_pdb(&original, path_str).unwrap();

        let loaded = read_pdb(path_str).unwrap();
        assert_eq!(loaded.n_frames(), 1);
        assert_eq!(loaded.first().unwrap().n_atoms(), 3);
    }

    fn cubic_frame(l: f64) -> Frame {
        let mut frame = Frame::new();
        frame.add_atom(Atom::new("O", Vector3::new(1.0, 2.0, 3.0)));
        frame.cell = Some(Cell::from_matrix(Matrix3::from_diagonal_element(l)));
        frame.pbc = [true; 3];
        frame
    }

    // 审查 A2：NPT 轨迹每帧胞不同，只写第 0 帧的 CRYST1 会让读回的所有帧共用一个胞
    #[test]
    fn test_cryst1_per_model() {
        let traj = Trajectory { frames: vec![cubic_frame(10.0), cubic_frame(11.0)], ..Default::default() };
        let path = std::env::temp_dir().join("ferro_pdb_cryst1_per_model.pdb");
        write_pdb(&traj, &path).unwrap();

        let text = std::fs::read_to_string(&path).unwrap();
        let tags: Vec<&str> = text.lines()
            .map(|l| &l[..l.len().min(6)])
            .filter(|t| *t == "CRYST1" || *t == "MODEL ")
            .collect();
        assert_eq!(tags, ["CRYST1", "MODEL ", "CRYST1", "MODEL "], "每个 MODEL 前应有本帧的 CRYST1");

        let loaded = read_pdb(&path).unwrap();
        for (f, l) in loaded.frames.iter().zip([10.0_f64, 11.0]) {
            let v = f.cell.as_ref().unwrap().volume();
            assert!((v - l.powi(3)).abs() / l.powi(3) < 1e-4, "读回体积 {v}，应为 {}", l.powi(3));
        }
    }

    // CRYST1 只存 a b c α β γ，reader 按 a 沿 x、b 在 xy 平面重建胞；
    // 非标准取向的胞若不先旋转坐标，读回后坐标与胞对不上
    #[test]
    fn test_rotated_cell_keeps_fractional() {
        // a 沿 y、b 沿 -x：右手系，但不是标准取向
        let m = Matrix3::new(
            0.0, 5.0, 0.0,
            -6.0, 0.0, 0.0,
            0.0, 0.0, 7.0,
        );
        let cell = Cell::from_matrix(m);
        let fracs = [Vector3::new(0.1, 0.2, 0.3), Vector3::new(0.7, 0.4, 0.9)];
        let mut frame = Frame::new();
        for s in &fracs {
            frame.add_atom(Atom::new("O", cell.fractional_to_cartesian(*s)));
        }
        frame.cell = Some(cell);
        frame.pbc = [true; 3];

        let path = std::env::temp_dir().join("ferro_pdb_rotated_cell.pdb");
        write_pdb(&Trajectory::from_frame(frame), &path).unwrap();

        let loaded = read_pdb(&path).unwrap();
        let f = loaded.first().unwrap();
        let c = f.cell.as_ref().unwrap();
        for (atom, s) in f.atoms.iter().zip(&fracs) {
            let got = c.cartesian_to_fractional(atom.position).unwrap();
            assert!((got - s).norm() < 1e-3, "读回分数坐标 {got:?}，应为 {s:?}");
        }
    }
}
