use std::path::Path;
use crate::readers::lammps_dump::LammpsUnits;
use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};

const EV_TO_KCAL: f64 = 1.0 / 0.04336410; // eV/Å → kcal/(mol·Å)

/// 写 LAMMPS dump 文件。
/// 包含列：id type element x y z [vx vy vz] [fx fy fz] [q]
pub fn write_lammps_dump(trajectory: &Trajectory, path: &Path, units: LammpsUnits) -> Result<()> {
    let path_ = path.display();
    let file = File::create(path).with_context(|| format!("cannot create {path_}"))?;
    let mut w = BufWriter::new(file);

    // type 编号在**整条轨迹**上确定一次。逐帧重建会让编号跟着「该帧碰巧先出现
    // 哪个名字」走 —— 元素只有几种时看不出来,但写位点标签时(某帧恰好没有 P_4)
    // 同一个 type 会在不同帧指向不同的东西
    let mut elem_types: Vec<&str> = Vec::new();
    for frame in &trajectory.frames {
        for atom in &frame.atoms {
            if !elem_types.contains(&atom.element.as_str()) {
                elem_types.push(atom.element.as_str());
            }
        }
    }

    for (ts, frame) in trajectory.frames.iter().enumerate() {
        let n = frame.n_atoms();

        writeln!(w, "ITEM: TIMESTEP")?;
        writeln!(w, "{ts}")?;
        writeln!(w, "ITEM: NUMBER OF ATOMS")?;
        writeln!(w, "{n}")?;

        // BOX BOUNDS
        let (lx, ly, lz, xy, xz, yz) = match &frame.cell {
            Some(cell) => cell_to_lammps(cell),
            None => {
                let (mnx, mxx, mny, mxy, mnz, mxz) = bounding_box(frame);
                (mxx - mnx, mxy - mny, mxz - mnz, 0.0, 0.0, 0.0)
            }
        };

        let is_triclinic = xy != 0.0 || xz != 0.0 || yz != 0.0;
        if is_triclinic {
            // 三斜时 dump 的六个数是 *_bound（倾斜后的外接盒），不是 xlo/xhi ——
            // 倾斜向量会把盒子探出 [0, l) 之外，bound 把那部分包进来。写成 xlo/xhi
            // 的话 reader 端 hi-lo 得到的边长偏小，而读写两侧一致地错时往返测试
            // 全绿（与 extxyz 应力符号同一个陷阱）
            let (xlo_b, xhi_b) = (min4(0.0, xy, xz, xy + xz), lx + max4(0.0, xy, xz, xy + xz));
            let (ylo_b, yhi_b) = (yz.min(0.0), ly + yz.max(0.0));
            writeln!(w, "ITEM: BOX BOUNDS xy xz yz pp pp pp")?;
            writeln!(w, "{xlo_b:.10} {xhi_b:.10} {xy:.10}")?;
            writeln!(w, "{ylo_b:.10} {yhi_b:.10} {xz:.10}")?;
            writeln!(w, "{:.10} {:.10} {yz:.10}", 0.0, lz)?;
        } else {
            writeln!(w, "ITEM: BOX BOUNDS pp pp pp")?;
            writeln!(w, "{:.10} {:.10}", 0.0, lx)?;
            writeln!(w, "{:.10} {:.10}", 0.0, ly)?;
            writeln!(w, "{:.10} {:.10}", 0.0, lz)?;
        }

        // ATOMS header
        let has_vel = frame.velocities.is_some();
        let has_force = frame.forces.is_some();
        let has_charge = frame.atoms.iter().any(|a| a.charge.is_some());

        let mut header = "ITEM: ATOMS id type element x y z".to_string();
        if has_vel { header.push_str(" vx vy vz"); }
        if has_force { header.push_str(" fx fy fz"); }
        if has_charge { header.push_str(" q"); }
        writeln!(w, "{header}")?;

        let lammps_cell = lammps_cell_matrix(lx, ly, lz, xy, xz, yz);
        // 原胞 → LAMMPS 规范胞（下三角）的线性变换，即 ASE `Prism.rot_mat` 那个旋转。
        // 坐标、速度、力**必须用同一个变换**：只转坐标的话，非下三角胞下 v、f 还在
        // 原坐标系里，与坐标对不上（ASE 的 lammpsdata / lammpsrun 三者同转）
        let to_lammps = |v: nalgebra::Vector3<f64>| -> Result<nalgebra::Vector3<f64>> {
            Ok(match &frame.cell {
                Some(orig) => lammps_cell.fractional_to_cartesian(orig.cartesian_to_fractional(v)?),
                None => v,
            })
        };

        for (i, atom) in frame.atoms.iter().enumerate() {
            let tp = elem_types.iter().position(|e| *e == atom.element).unwrap_or(0) + 1;

            let pos = to_lammps(atom.position)?;

            let mut line = format!("{} {tp} {} {:.10} {:.10} {:.10}",
                i + 1, atom.element, pos.x, pos.y, pos.z);

            if has_vel {
                let v = frame.velocities.as_ref()
                    .and_then(|vv| vv.get(i))
                    .copied()
                    .unwrap_or_default();
                let v = to_lammps(v)?;
                // real: Å/fs (no-op), metal: Å/fs → Å/ps (×1000)
                let vscale = match units {
                    LammpsUnits::Real  => 1.0,
                    LammpsUnits::Metal => 1000.0,
                };
                line.push_str(&format!(" {:.10} {:.10} {:.10}", v.x * vscale, v.y * vscale, v.z * vscale));
            }
            if has_force {
                let f = frame.forces.as_ref()
                    .and_then(|ff| ff.get(i))
                    .copied()
                    .unwrap_or_default();
                let f = to_lammps(f)?;
                // real: eV/Å → kcal/(mol·Å), metal: eV/Å (no-op)
                let fscale = match units {
                    LammpsUnits::Real  => EV_TO_KCAL,
                    LammpsUnits::Metal => 1.0,
                };
                line.push_str(&format!(" {:.10} {:.10} {:.10}", f.x * fscale, f.y * fscale, f.z * fscale));
            }
            if has_charge {
                line.push_str(&format!(" {:.6}", atom.charge.unwrap_or(0.0)));
            }

            writeln!(w, "{line}")?;
        }
    }

    w.flush()?;
    Ok(())
}

fn min4(a: f64, b: f64, c: f64, d: f64) -> f64 { a.min(b).min(c).min(d) }
fn max4(a: f64, b: f64, c: f64, d: f64) -> f64 { a.max(b).max(c).max(d) }

fn cell_to_lammps(cell: &ferro_core::Cell) -> (f64, f64, f64, f64, f64, f64) {
    let [a, b, c] = cell.lengths();
    let [alpha, beta, gamma] = cell.angles();
    let (al, be, ga) = (alpha.to_radians(), beta.to_radians(), gamma.to_radians());
    let lx = a;
    let xy = b * ga.cos();
    let xz = c * be.cos();
    let ly = (b * b - xy * xy).max(0.0).sqrt();
    let yz = if ly > 1e-10 { (b * c * al.cos() - xy * xz) / ly } else { 0.0 };
    let lz = (c * c - xz * xz - yz * yz).max(0.0).sqrt();
    (lx, ly, lz, xy, xz, yz)
}

fn lammps_cell_matrix(lx: f64, ly: f64, lz: f64, xy: f64, xz: f64, yz: f64) -> ferro_core::Cell {
    use nalgebra::Matrix3;
    ferro_core::Cell::from_matrix(Matrix3::new(lx, 0.0, 0.0, xy, ly, 0.0, xz, yz, lz))
}

fn bounding_box(frame: &ferro_core::Frame) -> (f64, f64, f64, f64, f64, f64) {
    let mut mn = [f64::MAX; 3]; let mut mx = [f64::MIN; 3];
    for a in &frame.atoms {
        mn[0]=mn[0].min(a.position.x); mn[1]=mn[1].min(a.position.y); mn[2]=mn[2].min(a.position.z);
        mx[0]=mx[0].max(a.position.x); mx[1]=mx[1].max(a.position.y); mx[2]=mx[2].max(a.position.z);
    }
    (mn[0],mx[0],mn[1],mx[1],mn[2],mx[2])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::lammps_dump::read_lammps_dump;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    fn bcc_traj() -> Trajectory {
        let cell = Cell::from_lengths_angles(2.87, 2.87, 2.87, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(1.435, 1.435, 1.435)));
        let mut traj = Trajectory::new();
        traj.add_frame(frame.clone());
        traj.add_frame(frame);
        traj
    }

    /// type 编号在整条轨迹上确定一次。逐帧重建时,某帧碰巧缺一个名字
    /// (标注轨迹里很常见:某帧没有 P_4)会让后续所有编号错位。
    #[test]
    fn test_type_numbering_is_stable_across_frames() {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mk = |elems: &[&str]| {
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            for (i, e) in elems.iter().enumerate() {
                f.add_atom(Atom::new(*e, Vector3::new(i as f64, 0.0, 0.0)));
            }
            f
        };
        // 帧 0 有 P_4,帧 1 没有 —— 逐帧编号会让帧 1 的 O_b 抢到 type 2
        let traj = Trajectory { frames: vec![mk(&["P_2", "P_4", "O_b"]), mk(&["P_2", "O_b"])],
                                metadata: Default::default() };

        let path = std::env::temp_dir().join("type_stable.lammpstrj");
        let p = &path;
        write_lammps_dump(&traj, p, LammpsUnits::Real).unwrap();

        let text = std::fs::read_to_string(p).unwrap();
        let mut seen: std::collections::HashMap<&str, &str> = std::collections::HashMap::new();
        for line in text.lines() {
            let c: Vec<&str> = line.split_whitespace().collect();
            if c.len() == 6 && c[0].parse::<usize>().is_ok() {
                if let Some(prev) = seen.insert(c[2], c[1]) {
                    assert_eq!(prev, c[1], "{} 的 type 编号在帧间变了", c[2]);
                }
            }
        }
        assert_eq!(seen.len(), 3, "三种名字各自一个稳定编号: {seen:?}");
    }

    #[test]
    fn test_roundtrip() {
        use crate::readers::lammps_dump::LammpsUnits;
        let path = std::env::temp_dir().join("bcc_rt.dump");
        let p = &path;
        let orig = bcc_traj();
        write_lammps_dump(&orig, p, LammpsUnits::Real).unwrap();

        let loaded = read_lammps_dump(p, LammpsUnits::Real).unwrap();
        assert_eq!(loaded.n_frames(), 2);
        let f = loaded.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
    }

    /// 坐标转进 LAMMPS 规范胞时，速度与力必须跟着同一个变换走。
    ///
    /// 期望值取自**独立实现**：ASE 3.29.0 的 `Prism(cell).vector_to_lammps`
    /// （`ase/calculators/lammps/coordinatetransform.py`），断言写出的文本而非读回，
    /// 理由同下一个测试。`rotated` 是 a 沿 y 的纯旋转胞；`general` 三条边都不沿轴
    #[test]
    fn test_non_lower_triangular_cell_rotates_velocity_and_force() {
        type V = [f64; 3];
        // (名字, 晶胞行, 坐标, 速度 Å/fs, 力 eV/Å, ASE 给出的 r, v, f)
        type Case = (&'static str, [V; 3], V, V, V, [V; 3]);
        let cases: [Case; 2] = [
            ("rotated",
             [[0.0, 4.0, 0.0], [-4.0, 0.0, 0.0], [0.0, 0.0, 4.0]],
             [1.0, 2.0, 3.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0],
             [[2.0, -1.0, 3.0], [0.0, -0.5, 0.0], [0.0, -1.0, 0.0]]),
            ("general",
             [[3.1, 0.4, -0.2], [0.7, 4.2, 0.5], [-0.3, 0.9, 5.3]],
             [1.2, -0.4, 2.5], [0.3, -0.7, 0.2], [-1.1, 0.6, 0.9],
             [[0.976982810926, -0.184550872348, 2.619436115402],
              [0.194758011328, -0.695518253611, 0.313566063083],
              [-1.069572685164, 0.845496016435, 0.721907720796]]),
        ];
        for (name, rows, r, v, f, want) in cases {
            let m = nalgebra::Matrix3::from_fn(|i, j| rows[i][j]);
            let mut frame = Frame::with_cell(Cell::from_matrix(m), [true; 3]);
            frame.add_atom(Atom::new("Si", Vector3::from(r)));
            frame.velocities = Some(vec![Vector3::from(v)]);
            frame.forces = Some(vec![Vector3::from(f)]);
            let path = std::env::temp_dir().join(format!("rot_{name}.dump"));
            // metal：力不换算，速度 Å/fs → Å/ps 乘 1000
            write_lammps_dump(&Trajectory::from_frame(frame), &path, LammpsUnits::Metal).unwrap();

            let text = std::fs::read_to_string(&path).unwrap();
            let cols: Vec<f64> = text.lines().last().unwrap()
                .split_whitespace().skip(3).map(|s| s.parse().unwrap()).collect();
            let got = [
                [cols[0], cols[1], cols[2]],
                [cols[3] / 1000.0, cols[4] / 1000.0, cols[5] / 1000.0],
                [cols[6], cols[7], cols[8]],
            ];
            for (q, (g, w)) in ["坐标", "速度", "力"].iter().zip(got.iter().zip(want.iter())) {
                for k in 0..3 {
                    assert!((g[k] - w[k]).abs() < 1e-8, "{name} 的{q}第 {k} 分量：ferro {}，ASE {}", g[k], w[k]);
                }
            }
        }
    }

    /// 三斜时写出的六个数是 `*_bound`，不是 xlo/xhi。
    ///
    /// Asserted on the literal text rather than through a read-back: the reader
    /// applies the inverse of whatever the writer does, so a round-trip stays green
    /// even when both sides share the same wrong convention — the same trap as the
    /// extxyz stress sign. Expected values follow the LAMMPS dump spec,
    /// xlo_bound = xlo + MIN(0,xy,xz,xy+xz) and xhi_bound = xhi + MAX(...).
    #[test]
    fn test_triclinic_writes_bound_not_xlo_xhi() {
        // lx/ly/lz = 10/12/14, xy/xz/yz = 2/-3/1
        let cell = Cell::from_matrix(nalgebra::Matrix3::new(
            10.0, 0.0, 0.0,
            2.0, 12.0, 0.0,
            -3.0, 1.0, 14.0,
        ));
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Si", Vector3::new(0.5, 0.5, 0.5)));
        let mut traj = Trajectory::new();
        traj.add_frame(frame);
        let path = std::env::temp_dir().join("tri_w.dump");
        write_lammps_dump(&traj, &path, crate::readers::lammps_dump::LammpsUnits::Metal).unwrap();

        let text = std::fs::read_to_string(&path).unwrap();
        let box_lines: Vec<&str> = text
            .lines()
            .skip_while(|l| !l.starts_with("ITEM: BOX BOUNDS"))
            .skip(1)
            .take(3)
            .collect();
        let nums: Vec<Vec<f64>> = box_lines
            .iter()
            .map(|l| l.split_whitespace().map(|s| s.parse().unwrap()).collect())
            .collect();
        // xlo_b = min(0,2,-3,-1) = -3 ; xhi_b = 10 + max(0,2,-3,-1) = 12
        let want = [[-3.0, 12.0, 2.0], [0.0, 13.0, -3.0], [0.0, 14.0, 1.0]];
        for (i, row) in want.iter().enumerate() {
            for (j, w) in row.iter().enumerate() {
                assert!(
                    (nums[i][j] - w).abs() < 1e-9,
                    "box line {i} field {j}: got {}, want {w} (bound, not xlo/xhi)",
                    nums[i][j]
                );
            }
        }
    }
}
