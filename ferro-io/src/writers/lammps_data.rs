use std::path::Path;
use std::collections::HashMap;
use ferro_core::Trajectory;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{bail, Context, Result};

/// 无胞帧的盒子取包围盒每侧外扩这么多（Å），坐标原样写（dump 与 data 都允许 lo ≠ 0）。
/// 余量无物理意义，只保证 lo < hi（平面分子某维厚度为 0）、原子不落在 hi 边界上被
/// `f` 边界判丢。与 `writers/lammps_dump.rs` 同值
const CELLLESS_PAD: f64 = 1.0;

/// 写 LAMMPS data 文件，atom_style full，real 单位（长度 Å）。
/// 周期性帧写晶格信息；非周期性帧写最小包围盒。
///
/// Every atom needs a value in the full style's charge column: with no charges at all it is 0
/// (the structure carries none), with charges on some atoms only it is an error, and with
/// `drop_partial` every charge is written as 0 (see `partial_column`).
pub fn write_lammps_data(trajectory: &Trajectory, path: &Path, drop_partial: bool) -> Result<()> {
    let path_ = path.display();
    let frame = super::single_frame(trajectory, "LAMMPS data")?;

    let n = frame.n_atoms();
    let with_charge = frame.atoms().filter(|a| a.charge.is_some()).count();
    let use_charge = super::partial_column("charge", with_charge, n, 0, drop_partial)?;

    // Unique elements in first-appearance order → type IDs
    let elem_order = frame.unique_elements();
    let n_types = elem_order.len();
    let elem_to_type: HashMap<&str, usize> = elem_order.iter()
        .enumerate()
        .map(|(i, e)| (e.as_str(), i + 1))
        .collect();

    // The atom's own mass (a data file's Masses value, an isotope), not the table's.
    // One mass per type: atoms of one element that disagree cannot be written.
    // Checked before the file is created, so a refusal leaves nothing half-written.
    let type_masses = elem_order.iter().map(|elem| {
        // An element with no mass is an error (review A-3; it used to be written as 1 amu)
        let masses = frame.atoms().filter(|a| a.element == *elem)
            .map(|a| a.effective_mass()).collect::<ferro_core::Result<Vec<f64>>>()?;
        let mass = masses[0]; // elem_order lists the frame's own elements
        match masses.iter().copied().find(|m| (m - mass).abs() > 1e-6) {
            Some(other) => bail!("atoms of element {elem} carry different masses ({mass} and {other}); \
                                  LAMMPS data gives one mass per atom type"),
            None => Ok(mass),
        }
    }).collect::<Result<Vec<f64>>>()?;

    let file = File::create(path).with_context(|| format!("cannot create {path_}"))?;
    let mut w = BufWriter::new(file);

    let comment = trajectory.metadata.title.as_deref().unwrap_or("LAMMPS data file written by ferro");
    writeln!(w, "{comment}")?;
    writeln!(w)?;

    writeln!(w, "{n} atoms")?;
    writeln!(w, "{n_types} atom types")?;
    writeln!(w)?;

    // Box bounds
    let (lo, (lx, ly, lz, xy, xz, yz)) = match &frame.cell {
        Some(cell) => (cell.origin.into(), cell_to_lammps(cell)),
        None => {
            // Non-periodic: compute bounding box
            let (minx, maxx, miny, maxy, minz, maxz) = bounding_box(frame);
            let w = 2.0 * CELLLESS_PAD;
            ([minx - CELLLESS_PAD, miny - CELLLESS_PAD, minz - CELLLESS_PAD],
             (maxx - minx + w, maxy - miny + w, maxz - minz + w, 0.0, 0.0, 0.0))
        }
    };

    // The box on a reduced basis; the positions below go through the unreduced cell
    let (bxy, bxz, byz) = super::reduce_tilts(lx, ly, xy, xz, yz);
    let is_triclinic = bxy != 0.0 || bxz != 0.0 || byz != 0.0;

    writeln!(w, "{:.10} {:.10} xlo xhi", lo[0], lo[0] + lx)?;
    writeln!(w, "{:.10} {:.10} ylo yhi", lo[1], lo[1] + ly)?;
    writeln!(w, "{:.10} {:.10} zlo zhi", lo[2], lo[2] + lz)?;
    if is_triclinic {
        writeln!(w, "{bxy:.10} {bxz:.10} {byz:.10} xy xz yz")?;
    }
    writeln!(w)?;

    // Masses section
    writeln!(w, "Masses")?;
    writeln!(w)?;
    for (i, (elem, mass)) in elem_order.iter().zip(&type_masses).enumerate() {
        // Shortest exact form: `{:.4}` cut an explicit isotope mass (review R2-2)
        writeln!(w, "{} {mass}  # {elem}", i + 1)?;
    }
    writeln!(w)?;

    // Atoms section (full style)
    writeln!(w, "Atoms # full")?;
    writeln!(w)?;

    // Build LAMMPS cell for coordinate transformation
    let lammps_cell = lammps_cell_matrix(lx, ly, lz, xy, xz, yz);

    for (i, atom) in frame.atoms().enumerate() {
        let tp = elem_to_type[atom.element];
        let q = if use_charge { atom.charge.unwrap_or(0.0) } else { 0.0 };

        // Transform position to LAMMPS coordinate frame
        let pos = match &frame.cell {
            // 旋转绕盒子原点做：非规范取向的胞带非零原点时，原点那段平移不能跟着转
            Some(orig_cell) => {
                let frac = orig_cell.cartesian_to_fractional(atom.position - orig_cell.origin)?;
                lammps_cell.fractional_to_cartesian(frac) + orig_cell.origin
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
    for a in frame.atoms() {
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
    use ferro_core::{Atom, Frame, Trajectory};
    use nalgebra::Vector3;

    fn bcc_traj() -> Trajectory {
        Trajectory::from_frame(crate::testutil::bcc_fe())
    }

    #[test]
    fn test_roundtrip() {
        let path = std::env::temp_dir().join("bcc_rt.lammps");
        let p = &path;
        write_lammps_data(&bcc_traj(), p, false).unwrap();

        let loaded = read_lammps_data(p, AtomStyle::Full).unwrap().0;
        let f = loaded.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
        let [a, ..] = f.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-4);
    }

    // Review A-3: the Masses section used to say 1.0 for an element with no mass
    #[test]
    fn test_element_with_no_mass_is_an_error() {
        let mut traj = bcc_traj();
        traj.frames[0].add_atom(Atom::new("Xx", Vector3::new(0.5, 0.5, 0.5)));
        let path = std::env::temp_dir().join("ferro_data_nomass.lammps");
        let _ = std::fs::remove_file(&path);
        let e = write_lammps_data(&traj, &path, false).unwrap_err();
        assert!(format!("{e:#}").contains("no mass for element 'Xx'"), "{e:#}");
        assert!(!path.exists(), "a refused write must not leave a file");
    }

    // 审查 A3：无胞帧以前按包围盒尺寸写成 0..L、坐标不平移，原子落在盒外
    #[test]
    fn test_cellless_frame_box_contains_atoms() {
        let mut f = Frame::new();
        f.add_atom(Atom::new("O", Vector3::new(1.0, 2.0, 3.0)));
        f.add_atom(Atom::new("H", Vector3::new(1.76, 2.59, 3.0)));
        f.add_atom(Atom::new("H", Vector3::new(0.24, 2.59, 3.0)));
        let path = std::env::temp_dir().join("ferro_data_cellless.lammps");
        write_lammps_data(&Trajectory::from_frame(f.clone()), &path, false).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        for (lo, hi, tag) in [(-0.76, 2.76, "xlo xhi"), (1.0, 3.59, "ylo yhi"), (2.0, 4.0, "zlo zhi")] {
            let line = text.lines().find(|l| l.ends_with(tag)).unwrap();
            let v: Vec<f64> = line.split_whitespace().take(2).map(|x| x.parse().unwrap()).collect();
            assert!((v[0] - lo).abs() < 1e-9 && (v[1] - hi).abs() < 1e-9, "{tag} 应为 {lo} {hi}，实得 {line}");
        }
        let back = read_lammps_data(&path, AtomStyle::Full).unwrap().0.frames.remove(0);
        for (x, y) in back.atoms().zip(f.atoms()) {
            assert!((x.position - y.position).norm() < 1e-9, "坐标应原样写出");
        }
    }

    // 审查 A14：90° 的 cos 残留让正交盒写出 xy xz yz 行
    #[test]
    fn test_orthogonal_cell_from_angles_is_not_triclinic() {
        let path = std::env::temp_dir().join("ferro_data_ortho.lammps");
        write_lammps_data(&bcc_traj(), &path, false).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(!text.contains("xy xz yz"), "正交盒不应写倾斜量：\n{text}");
    }

    // 审查 A-12：data 往返同样丢 lo
    /// An explicit isotope mass survives the round trip; `{:.4}` used to cut deuterium's
    /// 2.014101778 to 2.0141.
    #[test]
    fn test_isotope_mass_is_written_in_full() {
        let text = "t\n\n1 atoms\n1 atom types\n\n0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\n\
                    Masses\n\n1 2.014101778  # H\n\nAtoms # full\n\n1 1 1 0.0 1.0 1.0 1.0\n";
        let src = std::env::temp_dir().join("ferro_data_isotope_in.lammps");
        let dst = std::env::temp_dir().join("ferro_data_isotope_out.lammps");
        std::fs::write(&src, text).unwrap();
        write_lammps_data(&read_lammps_data(&src, AtomStyle::Full).unwrap().0, &dst, false).unwrap();
        let back = read_lammps_data(&dst, AtomStyle::Full).unwrap().0;
        assert_eq!(back.frames[0].atom(0).effective_mass().unwrap(), 2.014101778);
    }

    #[test]
    fn test_box_origin_survives_roundtrip() {
        let text = "test\n\n1 atoms\n1 atom types\n\n1.5 11.5 xlo xhi\n-2.0 8.0 ylo yhi\n3.0 13.0 zlo zhi\n\
                    2.0 -1.0 0.5 xy xz yz\n\nMasses\n\n1 15.999  # O\n\nAtoms # full\n\n1 1 1 0.0 2.0 -1.0 12.0\n";
        let src = std::env::temp_dir().join("ferro_data_origin_in.lammps");
        let dst = std::env::temp_dir().join("ferro_data_origin_out.lammps");
        std::fs::write(&src, text).unwrap();
        let a = read_lammps_data(&src, AtomStyle::Full).unwrap().0;
        write_lammps_data(&a, &dst, false).unwrap();
        let b = read_lammps_data(&dst, AtomStyle::Full).unwrap().0;
        let (fa, fb) = (&a.frames[0], &b.frames[0]);
        let (ca, cb) = (fa.cell.as_ref().unwrap(), fb.cell.as_ref().unwrap());
        assert!((ca.origin - Vector3::new(1.5, -2.0, 3.0)).norm() < 1e-12, "原点应为 lo，实得 {:?}", ca.origin);
        assert!((ca.origin - cb.origin).norm() < 1e-9, "原点应往返不变：{:?} 变成 {:?}", ca.origin, cb.origin);
        assert!((ca.matrix - cb.matrix).norm() < 1e-9, "晶格应往返不变");
        assert!((fa.positions()[0] - fb.positions()[0]).norm() < 1e-9, "坐标应往返不变");
    }
}
