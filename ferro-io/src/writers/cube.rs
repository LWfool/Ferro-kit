//! Gaussian cube 格式写入器

use std::io::{BufWriter, Write};
use anyhow::{Context, Result};
use ferro_core::CubeData;
use ferro_core::units::ANG_TO_BOHR;
use ferro_core::data::elements::by_symbol;

pub fn write_cube(path: &str, cube: &CubeData) -> Result<()> {
    let file = std::fs::File::create(path)
        .with_context(|| format!("cannot create {path}"))?;
    let mut w = BufWriter::new(file);
    write_cube_to(&mut w, cube).context("writing cube data")
}

fn write_cube_to<W: Write>(w: &mut W, cube: &CubeData) -> Result<()> {
    let frame = &cube.frame;
    let n_atoms = frame.n_atoms();
    let (nx, ny, nz) = cube.shape();

    // 2行注释
    writeln!(w, "Cube file written by ferro")?;
    writeln!(w, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z")?;

    // 原点（Å → Bohr）
    let ox = cube.origin.x * ANG_TO_BOHR;
    let oy = cube.origin.y * ANG_TO_BOHR;
    let oz = cube.origin.z * ANG_TO_BOHR;
    writeln!(w, "{:5}  {:12.6}  {:12.6}  {:12.6}", n_atoms, ox, oy, oz)?;

    // 3行体素定义：n_i  vx  vy  vz（Å → Bohr）
    let ns = [nx, ny, nz];
    for i in 0..3 {
        let vx = cube.spacing[(i, 0)] * ANG_TO_BOHR;
        let vy = cube.spacing[(i, 1)] * ANG_TO_BOHR;
        let vz = cube.spacing[(i, 2)] * ANG_TO_BOHR;
        writeln!(w, "{:5}  {:12.6}  {:12.6}  {:12.6}", ns[i], vx, vy, vz)?;
    }

    // 原子行：Z  0.000000  x  y  z（Å → Bohr）
    for atom in &frame.atoms {
        let z = by_symbol(&atom.element)
            .map(|e| e.atomic_number as i32)
            .unwrap_or(0);
        let x = atom.position.x * ANG_TO_BOHR;
        let y = atom.position.y * ANG_TO_BOHR;
        let z_pos = atom.position.z * ANG_TO_BOHR;
        writeln!(w, "{:5}  {:12.6}  {:12.6}  {:12.6}  {:12.6}", z, 0.0, x, y, z_pos)?;
    }

    // 体积数据：每行 6 个值，科学计数法
    let data_flat = cube.data.as_slice().context("data not contiguous")?;
    for (i, val) in data_flat.iter().enumerate() {
        if i > 0 && i % 6 == 0 { writeln!(w)?; }
        write!(w, " {:12.5e}", val)?;
    }
    writeln!(w)?;

    Ok(())
}

// ─── 测试 ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::cube::read_cube;
    use ferro_core::{Atom, Cell, CubeData, Frame};
    use nalgebra::{Matrix3, Vector3};
    use ndarray::Array3;

    /// 构建一个简单的 2×2×2 CubeData（C 立方，密度递增）
    fn make_test_cube() -> CubeData {
        let a = 2.0_f64; // Å
        let cell = Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("C", Vector3::new(0.5, 0.5, 0.5)));

        // 体素 = a/n = 1.0 Å；spacing 对角矩阵
        let voxel = a / 2.0;
        let spacing = Matrix3::from_diagonal(&Vector3::new(voxel, voxel, voxel));

        let data = Array3::from_shape_vec(
            (2, 2, 2),
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
        ).unwrap();

        CubeData { frame, data, origin: Vector3::zeros(), spacing }
    }

    #[test]
    fn test_roundtrip() {
        let original = make_test_cube();
        let path = std::env::temp_dir().join("roundtrip.cube");
        let path_str = path.to_str().unwrap();

        write_cube(path_str, &original).expect("write_cube failed");
        let loaded = read_cube(path_str).expect("read_cube failed");

        // 原子数
        assert_eq!(loaded.frame.n_atoms(), original.frame.n_atoms());
        // 元素
        assert_eq!(loaded.frame.atom(0).element, "C");
        // 原子位置（Å，允许小数值误差）
        let dp = (loaded.frame.atom(0).position - original.frame.atom(0).position).norm();
        assert!(dp < 1e-4, "position drift = {:.2e} Å", dp);
        // cell 边长
        let [a0, ..] = original.frame.cell.as_ref().unwrap().lengths();
        let [a1, ..] = loaded.frame.cell.as_ref().unwrap().lengths();
        assert!((a0 - a1).abs() < 1e-4, "cell drift = {:.2e}", (a0-a1).abs());
        // 网格形状
        assert_eq!(loaded.shape(), original.shape());
        // 体积数据值
        for idx in ndarray::indices(original.data.raw_dim()) {
            let diff = (loaded.data[idx] - original.data[idx]).abs();
            assert!(diff < 1e-4, "data[{:?}] diff = {:.2e}", idx, diff);
        }
    }

    #[test]
    fn test_write_produces_valid_header() {
        let cube = make_test_cube();
        let path = std::env::temp_dir().join("header_check.cube");
        write_cube(path.to_str().unwrap(), &cube).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        let mut lines = content.lines();
        assert!(lines.next().is_some(), "line 1 comment");
        assert!(lines.next().is_some(), "line 2 comment");
        // 第3行：原子数=1
        let line3 = lines.next().unwrap();
        assert!(line3.trim_start().starts_with('1'), "n_atoms=1: {}", line3);
    }

    #[test]
    fn test_write_atomic_number() {
        let cube = make_test_cube();
        let path = std::env::temp_dir().join("z_check.cube");
        write_cube(path.to_str().unwrap(), &cube).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        // 原子行在第 7 行（2注释 + 1原点 + 3体素 + 1原子）
        let atom_line = content.lines().nth(6).unwrap();
        // C 的原子序数 = 6
        assert!(atom_line.trim_start().starts_with('6'), "Z=6 for C: {}", atom_line);
    }
}
