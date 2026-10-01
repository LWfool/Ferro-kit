use std::path::Path;
use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{bail, Context, Result};

/// 写 VASP5 POSCAR 格式，坐标使用 Direct（分数坐标），原子按元素分组。
pub fn write_poscar(trajectory: &Trajectory, path: &Path) -> Result<()> {
    let path_ = path.display();
    let frame = trajectory.first().context("trajectory is empty")?;
    let cell = frame.cell.as_ref().context("frame has no cell (POSCAR requires periodic frame)")?;

    if let Some(vels) = &frame.velocities {
        if vels.len() != frame.n_atoms() {
            bail!("{} velocities for {} atoms", vels.len(), frame.n_atoms());
        }
    }

    let file = File::create(path).with_context(|| format!("cannot create {path_}"))?;
    let mut w = BufWriter::new(file);

    // Comment
    writeln!(w, "{}", trajectory.metadata.source.as_deref().unwrap_or("ferro"))?;
    // Scale
    writeln!(w, "  1.00000000000000")?;
    // Lattice vectors (row = lattice vector)
    for i in 0..3 {
        let r = cell.matrix.row(i);
        writeln!(w, "   {:>20.16}  {:>20.16}  {:>20.16}", r[0], r[1], r[2])?;
    }

    // Unique elements in first-appearance order
    let elem_order = frame.unique_elements();
    let counts: Vec<usize> = elem_order.iter()
        .map(|e| frame.count_element(e))
        .collect();

    // Element and count lines (VASP5)
    writeln!(w, "  {}", elem_order.join("   "))?;
    writeln!(w, "  {}", counts.iter().map(|n| n.to_string()).collect::<Vec<_>>().join("   "))?;

    // 原子按元素分组写出（与计数行对应）。分组后的下标序只算一次，坐标与速度都按它
    // 遍历 —— 速度块必须与坐标块逐行对应（VASP 文档：velocities of each ion defined
    // in the "Ion positions" section），此前速度按原序写，元素交错时对错了原子
    let order: Vec<usize> = elem_order.iter()
        .flat_map(|elem| frame.atoms.iter().enumerate()
            .filter(move |(_, a)| a.element == *elem)
            .map(|(i, _)| i))
        .collect();

    writeln!(w, "Direct")?;
    for &i in &order {
        let f = cell.cartesian_to_fractional(frame.atoms[i].position)?;
        writeln!(w, "  {:>18.16}  {:>18.16}  {:>18.16}", f.x, f.y, f.z)?;
    }

    // Optional velocities（空行 = Cartesian，Å/fs）
    if let Some(vels) = &frame.velocities {
        writeln!(w)?;
        for &i in &order {
            let vel = vels[i];
            writeln!(w, "  {:>18.16}  {:>18.16}  {:>18.16}", vel.x, vel.y, vel.z)?;
        }
    }

    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::vasp::read_poscar;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    fn bcc_traj() -> Trajectory {
        let cell = Cell::from_lengths_angles(2.87, 2.87, 2.87, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(1.435, 1.435, 1.435)));
        let mut traj = Trajectory::from_frame(frame);
        traj.metadata.source = Some("BCC Fe".to_string());
        traj
    }

    #[test]
    fn test_roundtrip() {
        let path = std::env::temp_dir().join("bcc_rt.poscar");
        let p = &path;
        let orig = bcc_traj();
        write_poscar(&orig, p).unwrap();

        let loaded = read_poscar(p).unwrap();
        let f = loaded.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
        let [a, ..] = f.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-10);
    }

    /// 元素交错时坐标按元素分组写出，速度必须跟着同一个顺序走。
    ///
    /// 断言写出的文本而非读回：用速度的 x 分量给每个原子编号（等于它的 x 坐标 / 10），
    /// 坐标块第 k 行与速度块第 k 行必须是同一个原子
    #[test]
    fn test_interleaved_elements_keep_velocity_with_its_atom() {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        for (el, x) in [("O", 1.0), ("Si", 2.0), ("O", 3.0), ("Si", 4.0)] {
            frame.add_atom(Atom::new(el, Vector3::new(x, 0.0, 0.0)));
        }
        frame.velocities = Some((1..=4).map(|k| Vector3::new(k as f64 / 10.0, 0.0, 0.0)).collect());
        let path = std::env::temp_dir().join("interleaved_vel.poscar");
        write_poscar(&Trajectory::from_frame(frame), &path).unwrap();

        let text = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let first_x = |l: &str| l.split_whitespace().next().unwrap().parse::<f64>().unwrap();
        let start = lines.iter().position(|l| *l == "Direct").unwrap() + 1;
        let pos_x: Vec<f64> = lines[start..start + 4].iter().map(|l| first_x(l)).collect();
        let vel_x: Vec<f64> = lines[start + 5..start + 9].iter().map(|l| first_x(l)).collect();
        // 分组后是 O(x=1) O(x=3) Si(x=2) Si(x=4)，分数坐标 = x/10，速度 = x/10
        assert_eq!(pos_x, vec![0.1, 0.3, 0.2, 0.4], "坐标应按元素分组");
        assert_eq!(vel_x, pos_x, "速度块第 k 行应与坐标块第 k 行是同一个原子");
    }
}
