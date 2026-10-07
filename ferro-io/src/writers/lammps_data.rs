use std::path::Path;
use std::collections::HashMap;
use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};

/// 无胞帧的盒子取包围盒每侧外扩这么多（Å），坐标原样写（dump 与 data 都允许 lo ≠ 0）。
/// 余量无物理意义，只保证 lo < hi（平面分子某维厚度为 0）、原子不落在 hi 边界上被
/// `f` 边界判丢。与 `writers/lammps_dump.rs` 同值
const CELLLESS_PAD: f64 = 1.0;

/// 写 LAMMPS data 文件，atom_style full，real 单位（长度 Å）。
/// 周期性帧写晶格信息；非周期性帧写最小包围盒。
pub fn write_lammps_data(trajectory: &Trajectory, path: &Path) -> Result<()> {
    let path_ = path.display();
    let frame = trajectory.first().context("trajectory is empty")?;

    let file = File::create(path).with_context(|| format!("cannot create {path_}"))?;
    let mut w = BufWriter::new(file);

    let comment = trajectory.metadata.source.as_deref().unwrap_or("LAMMPS data file written by ferro");
    writeln!(w, "{comment}")?;
    writeln!(w)?;

    let n = frame.n_atoms();

    // Unique elements in first-appearance order → type IDs
    let elem_order = frame.unique_elements();
    let n_types = elem_order.len();
    let elem_to_type: HashMap<&str, usize> = elem_order.iter()
        .enumerate()
        .map(|(i, e)| (e.as_str(), i + 1))
        .collect();


    writeln!(w, "{n} atoms")?;
    writeln!(w, "{n_types} atom types")?;
    writeln!(w)?;

    // Box bounds
    let (lo, (lx, ly, lz, xy, xz, yz)) = match &frame.cell {
        Some(cell) => ([0.0; 3], cell_to_lammps(cell)),
        None => {
            // Non-periodic: compute bounding box
            let (minx, maxx, miny, maxy, minz, maxz) = bounding_box(frame);
            let w = 2.0 * CELLLESS_PAD;
            ([minx - CELLLESS_PAD, miny - CELLLESS_PAD, minz - CELLLESS_PAD],
             (maxx - minx + w, maxy - miny + w, maxz - minz + w, 0.0, 0.0, 0.0))
        }
    };

    let is_triclinic = xy != 0.0 || xz != 0.0 || yz != 0.0;

    writeln!(w, "{:.10} {:.10} xlo xhi", lo[0], lo[0] + lx)?;
    writeln!(w, "{:.10} {:.10} ylo yhi", lo[1], lo[1] + ly)?;
    writeln!(w, "{:.10} {:.10} zlo zhi", lo[2], lo[2] + lz)?;
    if is_triclinic {
        writeln!(w, "{:.10} {:.10} {:.10} xy xz yz", xy, xz, yz)?;
    }
    writeln!(w)?;

    // Masses section
    writeln!(w, "Masses")?;
    writeln!(w)?;
    for (i, elem) in elem_order.iter().enumerate() {
        let mass = ferro_core::data::elements::by_symbol(elem)
            .map(|e| e.atomic_mass)
            .unwrap_or(1.0);
        writeln!(w, "{} {:.4}  # {elem}", i + 1, mass)?;
    }
    writeln!(w)?;

    // Atoms section (full style)
    writeln!(w, "Atoms # full")?;
    writeln!(w)?;

    // Build LAMMPS cell for coordinate transformation
    let lammps_cell = lammps_cell_matrix(lx, ly, lz, xy, xz, yz);

    for (i, atom) in frame.atoms.iter().enumerate() {
        let tp = elem_to_type[atom.element.as_str()];
        let q = atom.charge.unwrap_or(0.0);

        // Transform position to LAMMPS coordinate frame
        let pos = match &frame.cell {
            Some(orig_cell) => {
                let frac = orig_cell.cartesian_to_fractional(atom.position)?;
                lammps_cell.fractional_to_cartesian(frac)
            }
            None => atom.position,
        };

        // id mol-id type charge x y z
        writeln!(w, "{} 1 {tp} {q:.6} {:.10} {:.10} {:.10}",
            i + 1, pos.x, pos.y, pos.z)?;
    }

    w.flush()?;
    Ok(())
}

// cell_to_lammps / lammps_cell_matrix / bounding_box 与 `writers/lammps_dump.rs` 逐字两份，有意不合
// （R6：只有两个使用点，合并只多一个跨文件依赖）。改一处必须同改另一处
/// Convert Cell to LAMMPS parameters: returns (lx, ly, lz, xy, xz, yz)
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
    // 90° 的 cos 残留（~1e-16·L）会让正交盒被判成三斜；相对最大边长 1e-10 以下的倾斜置 0
    let tol = 1e-10 * lx.max(ly).max(lz);
    let snap = |t: f64| if t.abs() < tol { 0.0 } else { t };
    let (xy, xz, yz) = (snap(xy), snap(xz), snap(yz));

    (lx, ly, lz, xy, xz, yz)
}

fn lammps_cell_matrix(lx: f64, ly: f64, lz: f64, xy: f64, xz: f64, yz: f64) -> ferro_core::Cell {
    use nalgebra::Matrix3;
    ferro_core::Cell::from_matrix(Matrix3::new(
        lx,  0.0, 0.0,
        xy,  ly,  0.0,
        xz,  yz,  lz,
    ))
}

fn bounding_box(frame: &ferro_core::Frame) -> (f64, f64, f64, f64, f64, f64) {
    let mut mn = [f64::MAX; 3];
    let mut mx = [f64::MIN; 3];
    for a in &frame.atoms {
        mn[0] = mn[0].min(a.position.x);
        mn[1] = mn[1].min(a.position.y);
        mn[2] = mn[2].min(a.position.z);
        mx[0] = mx[0].max(a.position.x);
        mx[1] = mx[1].max(a.position.y);
        mx[2] = mx[2].max(a.position.z);
    }
    (mn[0], mx[0], mn[1], mx[1], mn[2], mx[2])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::lammps_data::{read_lammps_data, AtomStyle};
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    fn bcc_traj() -> Trajectory {
        let cell = Cell::from_lengths_angles(2.87, 2.87, 2.87, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(1.435, 1.435, 1.435)));
        Trajectory::from_frame(frame)
    }

    #[test]
    fn test_roundtrip() {
        let path = std::env::temp_dir().join("bcc_rt.lammps");
        let p = &path;
        write_lammps_data(&bcc_traj(), p).unwrap();

        let loaded = read_lammps_data(p, AtomStyle::Full).unwrap();
        let f = loaded.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
        let [a, ..] = f.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-4);
    }

    // 审查 A3：无胞帧以前按包围盒尺寸写成 0..L、坐标不平移，原子落在盒外
    #[test]
    fn test_cellless_frame_box_contains_atoms() {
        let mut f = Frame::new();
        f.add_atom(Atom::new("O", Vector3::new(1.0, 2.0, 3.0)));
        f.add_atom(Atom::new("H", Vector3::new(1.76, 2.59, 3.0)));
        f.add_atom(Atom::new("H", Vector3::new(0.24, 2.59, 3.0)));
        let path = std::env::temp_dir().join("ferro_data_cellless.lammps");
        write_lammps_data(&Trajectory::from_frame(f.clone()), &path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        for (lo, hi, tag) in [(-0.76, 2.76, "xlo xhi"), (1.0, 3.59, "ylo yhi"), (2.0, 4.0, "zlo zhi")] {
            let line = text.lines().find(|l| l.ends_with(tag)).unwrap();
            let v: Vec<f64> = line.split_whitespace().take(2).map(|x| x.parse().unwrap()).collect();
            assert!((v[0] - lo).abs() < 1e-9 && (v[1] - hi).abs() < 1e-9, "{tag} 应为 {lo} {hi}，实得 {line}");
        }
        let back = read_lammps_data(&path, AtomStyle::Full).unwrap().frames.remove(0);
        for (x, y) in back.atoms.iter().zip(&f.atoms) {
            assert!((x.position - y.position).norm() < 1e-9, "坐标应原样写出");
        }
    }

    // 审查 A14：90° 的 cos 残留让正交盒写出 xy xz yz 行
    #[test]
    fn test_orthogonal_cell_from_angles_is_not_triclinic() {
        let path = std::env::temp_dir().join("ferro_data_ortho.lammps");
        write_lammps_data(&bcc_traj(), &path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(!text.contains("xy xz yz"), "正交盒不应写倾斜量：\n{text}");
    }
}
