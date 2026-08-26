use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{bail, Context, Result};

/// Which key carries the stress tensor in the comment line.
///
/// Both are legal extxyz and GPUMD reads either, but they are different
/// quantities: `stress` is eV/Å³ with ASE's sign (positive = tension), while
/// `virial` is eV with the sign [`ferro_core::Frame::stress`] already uses
/// (positive = compression). GPUMD prefers `virial` when a file carries both,
/// which is why the NEP path writes that one.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum StressKey {
    /// `stress="..."`, eV/Å³, negated on the way out (ASE convention).
    #[default]
    Stress,
    /// `virial="..."`, eV, `stress * V` with no sign change.
    Virial,
}

/// 写 extxyz 格式，多帧轨迹写为连续 block。
pub fn write_extxyz(trajectory: &Trajectory, path: &str) -> Result<()> {
    write_extxyz_with(trajectory, path, StressKey::Stress)
}

/// [`write_extxyz`] with a choice of stress key; see [`StressKey`].
pub fn write_extxyz_with(
    trajectory: &Trajectory,
    path: &str,
    stress_key: StressKey,
) -> Result<()> {
    let file = File::create(path).with_context(|| format!("cannot create {path}"))?;
    let mut w = BufWriter::new(file);

    for (fi, frame) in trajectory.frames.iter().enumerate() {
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
        // 恒只写一个键 —— 两个键就是两处可能互相矛盾的事实,读侧为此专门做了
        // 交叉校验,没有理由自己生产这种文件。
        //   stress= 是 ASE 约定(正 = 拉伸),故变号
        //   virial= 是 eV 且正 = 压缩,故乘体积、不变号
        if let Some(s) = &frame.stress {
            let (key, t) = match stress_key {
                StressKey::Stress => ("stress", -s),
                StressKey::Virial => {
                    let cell = frame.cell.as_ref().with_context(|| format!(
                        "frame {fi} has a stress but no cell; virial= is eV and needs \
                         the volume. Write stress= instead"
                    ))?;
                    let vol = cell.volume();
                    if vol.abs() < 1e-12 {
                        bail!("cannot write virial= for a zero-volume cell");
                    }
                    ("virial", s * vol)
                }
            };
            parts.push(format!(
                "{key}=\"{} {} {} {} {} {} {} {} {}\"",
                fmt(t[(0,0)]), fmt(t[(0,1)]), fmt(t[(0,2)]),
                fmt(t[(1,0)]), fmt(t[(1,1)]), fmt(t[(1,2)]),
                fmt(t[(2,0)]), fmt(t[(2,1)]), fmt(t[(2,2)]),
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

// -0.0 与 0.0 在这里没有区别,但 "-0.0000000000" 读起来像个 bug
fn fmt(v: f64) -> String {
    let v = if v == 0.0 { 0.0 } else { v };
    format!("{v:.10}")
}
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
    #[test]
    fn test_stress_written_in_ase_sign() {
        // 锚点是 ASE 3.29.0 对同一张量写出的文本,不是本 writer 自己的往返 ——
        // 写侧与读侧同时漏掉变号时,往返测试照样通过。
        //   σ_ferro(正 = 压缩) = -[[0.01,0.002,0.003],[...]]
        //   ASE 对应写出的 stress= 必须是 +0.01 0.002 0.003 ...
        use nalgebra::Matrix3;
        let mut traj = bcc_traj();
        traj.frames[0].stress = Some(-Matrix3::new(
            0.01, 0.002, 0.003,
            0.002, 0.02, 0.004,
            0.003, 0.004, 0.03));
        let path = std::env::temp_dir().join("stress_sign.extxyz");
        write_extxyz(&traj, path.to_str().unwrap()).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let line = text.lines().nth(1).unwrap();
        assert!(
            line.contains("stress=\"0.0100000000 0.0020000000 0.0030000000 \
0.0020000000 0.0200000000 0.0040000000 0.0030000000 0.0040000000 0.0300000000\""),
            "{line}");
        // virial= 不写:两个键 = 两处可能互相矛盾的事实
        assert!(!line.contains("virial="), "{line}");
    }

    #[test]
    fn test_virial_key_is_stress_times_volume() {
        // NEP 侧要 virial=:eV,正 = 压缩,故乘体积、不变号。锚点仍是外部换算式
        // (dpdata 1.0.2 的 virials = -V * stress_ase),不是本 writer 的往返。
        use nalgebra::Matrix3;
        let mut traj = bcc_traj();
        let vol = traj.frames[0].cell.as_ref().unwrap().volume();
        let sigma = Matrix3::new(0.01, 0.002, 0.003,
                                 0.002, 0.02, 0.004,
                                 0.003, 0.004, 0.03);
        traj.frames[0].stress = Some(-sigma);       // σ_ferro = -σ_ase
        let path = std::env::temp_dir().join("virial_key.extxyz");
        write_extxyz_with(&traj, path.to_str().unwrap(), StressKey::Virial).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        let line = text.lines().nth(1).unwrap();
        assert!(line.contains("virial="), "{line}");
        assert!(!line.contains("stress="), "{line}");

        // 读回:reader 对 virial= 除体积、不变号,应回到 σ_ferro
        let back = read_extxyz(path.to_str().unwrap()).unwrap();
        let s = back.first().unwrap().stress.unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert!((s[(i, j)] + sigma[(i, j)]).abs() < 1e-9);
            }
        }
        // 且文本里的数就是 -V*σ_ase
        let first: f64 = line.split("virial=\"").nth(1).unwrap()
            .split_whitespace().next().unwrap().parse().unwrap();
        assert!((first - (-vol * 0.01)).abs() < 1e-6, "{first}");
    }

    #[test]
    fn test_virial_without_cell_is_an_error() {
        use nalgebra::Matrix3;
        let mut traj = bcc_traj();
        traj.frames[0].cell = None;
        traj.frames[0].stress = Some(Matrix3::identity());
        let path = std::env::temp_dir().join("virial_nocell.extxyz");
        let e = write_extxyz_with(&traj, path.to_str().unwrap(), StressKey::Virial)
            .unwrap_err();
        assert!(format!("{e:#}").contains("no cell"), "{e:#}");
    }

}
