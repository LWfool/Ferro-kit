use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};

/// 写 extxyz 格式，多帧轨迹写为连续 block。
pub fn write_extxyz(trajectory: &Trajectory, path: &str) -> Result<()> {
    let file = File::create(path).with_context(|| format!("cannot create {path}"))?;
    let mut w = BufWriter::new(file);

    for frame in &trajectory.frames {
        // Line 1: atom count
        writeln!(w, "{}", frame.n_atoms())?;

        // Line 2: key=value comment
        let mut parts: Vec<String> = Vec::new();

        // Lattice (row-major: 9 values = rows a, b, c)
        if let Some(cell) = &frame.cell {
            let m = &cell.matrix;
            parts.push(format!(
                "Lattice=\"{} {} {} {} {} {} {} {} {}\"",
                fmt(m[(0,0)]), fmt(m[(0,1)]), fmt(m[(0,2)]),
                fmt(m[(1,0)]), fmt(m[(1,1)]), fmt(m[(1,2)]),
                fmt(m[(2,0)]), fmt(m[(2,1)]), fmt(m[(2,2)]),
            ));
            let pbc = frame.pbc;
            parts.push(format!("pbc=\"{} {} {}\"",
                bool_str(pbc[0]), bool_str(pbc[1]), bool_str(pbc[2])));
        }

        // Properties spec
        let mut prop_spec = "species:S:1:pos:R:3".to_string();
        if frame.forces.is_some() { prop_spec.push_str(":forces:R:3"); }
        if frame.velocities.is_some() { prop_spec.push_str(":velocities:R:3"); }
        let has_charge = frame.atoms.iter().any(|a| a.charge.is_some());
        if has_charge { prop_spec.push_str(":charges:R:1"); }
        let has_magmom = frame.atoms.iter().any(|a| a.magmom.is_some());
        if has_magmom { prop_spec.push_str(":magmoms:R:1"); }
        // extxyz 的列是自描述的,故位点标签走自己的一列,species 保持纯元素。
        // 这是与 LAMMPS dump 的区别:那边没地方放第二个名字,只能折进 element 列
        let has_label = frame.atoms.iter().any(|a| a.label.is_some());
        if has_label { prop_spec.push_str(":label:S:1"); }
        parts.push(format!("Properties={prop_spec}"));

        if let Some(e) = frame.energy { parts.push(format!("energy={}", fmt(e))); }
        if let Some(s) = &frame.stress {
            parts.push(format!(
                "stress=\"{} {} {} {} {} {} {} {} {}\"",
                fmt(s[(0,0)]), fmt(s[(0,1)]), fmt(s[(0,2)]),
                fmt(s[(1,0)]), fmt(s[(1,1)]), fmt(s[(1,2)]),
                fmt(s[(2,0)]), fmt(s[(2,1)]), fmt(s[(2,2)]),
            ));
        }

        writeln!(w, "{}", parts.join(" "))?;

        // Atom lines
        let dummy_forces = vec![];
        let dummy_vels = vec![];
        let forces = frame.forces.as_deref().unwrap_or(&dummy_forces);
        let vels = frame.velocities.as_deref().unwrap_or(&dummy_vels);

        for (i, atom) in frame.atoms.iter().enumerate() {
            let mut line = format!("{} {} {} {}",
                atom.element,
                fmt(atom.position.x), fmt(atom.position.y), fmt(atom.position.z));

            if let Some(f) = forces.get(i) {
                line.push_str(&format!(" {} {} {}", fmt(f.x), fmt(f.y), fmt(f.z)));
            }
            if let Some(v) = vels.get(i) {
                line.push_str(&format!(" {} {} {}", fmt(v.x), fmt(v.y), fmt(v.z)));
            }
            if has_charge {
                line.push_str(&format!(" {}", fmt(atom.charge.unwrap_or(0.0))));
            }
            if has_magmom {
                line.push_str(&format!(" {}", fmt(atom.magmom.unwrap_or(0.0))));
            }
            if has_label {
                // 无标签的原子回退为元素符号:S 列不能留空,否则列数对不齐
                line.push_str(&format!(" {}", atom.label.as_deref().unwrap_or(&atom.element)));
            }
            writeln!(w, "{line}")?;
        }
    }

    w.flush()?;
    Ok(())
}

fn fmt(v: f64) -> String { format!("{v:.10}") }
fn bool_str(b: bool) -> &'static str { if b { "T" } else { "F" } }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::extxyz::read_extxyz;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    fn bcc_traj() -> Trajectory {
        let cell = Cell::from_lengths_angles(2.87, 2.87, 2.87, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(1.435, 1.435, 1.435)));
        frame.energy = Some(-17.43);
        frame.forces = Some(vec![
            Vector3::new(0.1, 0.0, 0.0),
            Vector3::new(-0.1, 0.0, 0.0),
        ]);
        Trajectory::from_frame(frame)
    }

    /// 位点标签走自己的 `label:S:1` 列,species 保持纯元素 —— 两个方向都无损。
    #[test]
    fn test_label_column_roundtrips_without_touching_species() {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        for (elem, label) in [("P", "P_3"), ("O", "O_b"), ("Zn", "Zn")] {
            let mut a = Atom::new(elem, Vector3::new(0.0, 0.0, 0.0));
            a.label = Some(label.to_string());
            frame.add_atom(a);
        }
        let traj = Trajectory::from_frame(frame);

        let path = std::env::temp_dir().join("label_col.extxyz");
        let p = path.to_str().unwrap();
        write_extxyz(&traj, p).unwrap();

        let text = std::fs::read_to_string(p).unwrap();
        assert!(text.contains("label:S:1"), "Properties 必须声明 label 列:\n{text}");
        let atom_line = text.lines().nth(2).unwrap();
        assert!(atom_line.starts_with("P "), "species 列仍是纯元素: {atom_line}");
        assert!(atom_line.ends_with("P_3"), "标签在自己的列里: {atom_line}");

        let back = read_extxyz(p).unwrap();
        let atoms = &back.frames[0].atoms;
        assert_eq!(atoms[0].element, "P");
        assert_eq!(atoms[0].label.as_deref(), Some("P_3"));
        assert_eq!(atoms[1].element, "O");
        assert_eq!(atoms[1].label.as_deref(), Some("O_b"));
        assert_eq!(atoms[2].label.as_deref(), Some("Zn"));
    }

    #[test]
    fn test_roundtrip() {
        let path = std::env::temp_dir().join("bcc_rt.extxyz");
        let p = path.to_str().unwrap();
        let orig = bcc_traj();
        write_extxyz(&orig, p).unwrap();

        let loaded = read_extxyz(p).unwrap();
        let f = loaded.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
        assert!((f.energy.unwrap() - (-17.43)).abs() < 1e-6);
        let forces = f.forces.as_ref().unwrap();
        assert!((forces[0].x - 0.1).abs() < 1e-6);
    }
}
