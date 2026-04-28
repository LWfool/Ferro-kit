use std::collections::HashMap;
use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};

/// 将轨迹写入 CIF 文件。每帧作为独立的 `data_` block 写出。
///
/// 周期性帧：用分数坐标（`_atom_site_fract_*`），写 P1 对称性。
/// 非周期性帧：用 Cartesian 坐标（`_atom_site_Cartn_*`），省略晶格参数。
pub fn write_cif(trajectory: &Trajectory, path: &str) -> Result<()> {
    let file = File::create(path).with_context(|| format!("cannot create {path}"))?;
    let mut w = BufWriter::new(file);

    for (idx, frame) in trajectory.frames.iter().enumerate() {
        // 生成 data block 名称
        let block_name = trajectory.metadata.source
            .as_deref()
            .map(|s| sanitize_block_name(s))
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| format!("structure_{}", idx + 1));

        // 多帧时给 block 名称加序号避免重复
        let block_name = if trajectory.n_frames() > 1 {
            format!("{block_name}_{}", idx + 1)
        } else {
            block_name
        };

        writeln!(w, "data_{block_name}")?;
        writeln!(w)?;

        if let Some(cell) = &frame.cell {
            let [a, b, c] = cell.lengths();
            let [alpha, beta, gamma] = cell.angles();

            writeln!(w, "_cell_length_a         {a:.6}")?;
            writeln!(w, "_cell_length_b         {b:.6}")?;
            writeln!(w, "_cell_length_c         {c:.6}")?;
            writeln!(w, "_cell_angle_alpha       {alpha:.4}")?;
            writeln!(w, "_cell_angle_beta        {beta:.4}")?;
            writeln!(w, "_cell_angle_gamma       {gamma:.4}")?;
            writeln!(w)?;
            writeln!(w, "_space_group_name_H-M   'P 1'")?;
            writeln!(w, "_space_group_IT_number  1")?;
            writeln!(w)?;

            writeln!(w, "loop_")?;
            writeln!(w, "_atom_site_label")?;
            writeln!(w, "_atom_site_type_symbol")?;
            writeln!(w, "_atom_site_fract_x")?;
            writeln!(w, "_atom_site_fract_y")?;
            writeln!(w, "_atom_site_fract_z")?;

            let mut counters: HashMap<String, usize> = HashMap::new();
            for atom in &frame.atoms {
                let label = if let Some(l) = &atom.label {
                    l.clone()
                } else {
                    let n = counters.entry(atom.element.clone()).or_insert(0);
                    *n += 1;
                    format!("{}{n}", atom.element)
                };
                let frac = cell.cartesian_to_fractional(atom.position);
                writeln!(
                    w,
                    "{label:<8} {:<4} {:>10.6} {:>10.6} {:>10.6}",
                    atom.element, frac.x, frac.y, frac.z
                )?;
            }
        } else {
            // 非周期性：Cartesian 坐标
            writeln!(w, "loop_")?;
            writeln!(w, "_atom_site_label")?;
            writeln!(w, "_atom_site_type_symbol")?;
            writeln!(w, "_atom_site_Cartn_x")?;
            writeln!(w, "_atom_site_Cartn_y")?;
            writeln!(w, "_atom_site_Cartn_z")?;

            let mut counters: HashMap<String, usize> = HashMap::new();
            for atom in &frame.atoms {
                let label = if let Some(l) = &atom.label {
                    l.clone()
                } else {
                    let n = counters.entry(atom.element.clone()).or_insert(0);
                    *n += 1;
                    format!("{}{n}", atom.element)
                };
                writeln!(
                    w,
                    "{label:<8} {:<4} {:>10.6} {:>10.6} {:>10.6}",
                    atom.element,
                    atom.position.x,
                    atom.position.y,
                    atom.position.z
                )?;
            }
        }

        writeln!(w)?;
    }

    w.flush()?;
    Ok(())
}

fn sanitize_block_name(s: &str) -> String {
    s.chars()
        .map(|c| if c.is_alphanumeric() || c == '_' || c == '-' { c } else { '_' })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::cif::read_cif;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    fn make_bcc_frame() -> Frame {
        let cell = Cell::from_lengths_angles(2.87, 2.87, 2.87, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(1.435, 1.435, 1.435)));
        frame
    }

    #[test]
    fn test_roundtrip_periodic() {
        let path = std::env::temp_dir().join("bcc_fe_rt.cif");
        let path_str = path.to_str().unwrap();

        let mut traj = Trajectory::from_frame(make_bcc_frame());
        traj.metadata.source = Some("BCC_Fe".to_string());
        write_cif(&traj, path_str).unwrap();

        let loaded = read_cif(path_str).unwrap();
        assert_eq!(loaded.n_frames(), 1);
        let frame = loaded.first().unwrap();
        assert_eq!(frame.n_atoms(), 2);
        assert_eq!(frame.atom(0).element, "Fe");

        let [a, ..] = frame.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-4);
    }

    #[test]
    fn test_roundtrip_nonperiodic() {
        let path = std::env::temp_dir().join("water_rt.cif");
        let path_str = path.to_str().unwrap();

        let mut frame = Frame::new();
        frame.add_atom(Atom::new("O", Vector3::new(0.0, 0.0, 0.119)));
        frame.add_atom(Atom::new("H", Vector3::new(0.0, 0.763, -0.477)));
        frame.add_atom(Atom::new("H", Vector3::new(0.0, -0.763, -0.477)));
        let traj = Trajectory::from_frame(frame);
        write_cif(&traj, path_str).unwrap();

        // 非周期性 CIF 不能 roundtrip via read_cif (无 cell)，只验证写入不报错
        let content = std::fs::read_to_string(path_str).unwrap();
        assert!(content.contains("_atom_site_Cartn_x"));
        assert!(content.contains("O1"));
    }

    #[test]
    fn test_multi_frame_block_names() {
        let path = std::env::temp_dir().join("multi_frame.cif");
        let path_str = path.to_str().unwrap();

        let mut traj = Trajectory::new();
        traj.add_frame(make_bcc_frame());
        traj.add_frame(make_bcc_frame());
        traj.metadata.source = Some("Fe".to_string());
        write_cif(&traj, path_str).unwrap();

        let content = std::fs::read_to_string(path_str).unwrap();
        assert!(content.contains("data_Fe_1"));
        assert!(content.contains("data_Fe_2"));
    }
}
