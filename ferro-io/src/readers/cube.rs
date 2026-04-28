//! Gaussian cube 格式读取器
//!
//! 格式规范：https://paulbourke.net/dataformats/cube/
//!
//! 文件结构：
//!   行1   注释
//!   行2   注释（可含轴顺序信息）
//!   行3   n_atoms  ox  oy  oz   （Bohr → Å）
//!   行4-6 n_i  vx  vy  vz      （体素数 + 体素向量，Bohr → Å）
//!   原子行 Z  charge  x  y  z  （Bohr → Å）
//!   体积数据 平铺浮点，XYZ 顺序

use ndarray::Array3;
use nalgebra::{Matrix3, Vector3};
use anyhow::{bail, Context, Result};
use ferro_core::{Atom, Cell, CubeData, Frame};
use ferro_core::units::BOHR_TO_ANG;
use ferro_core::data::elements::by_number;

pub fn read_cube(path: &str) -> Result<CubeData> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path}"))?;
    parse_cube(&content).with_context(|| format!("parsing {path}"))
}

fn parse_cube(content: &str) -> Result<CubeData> {
    let mut lines = content.lines();

    // 2行注释，直接跳过
    lines.next().context("missing comment line 1")?;
    lines.next().context("missing comment line 2")?;

    // 第3行：n_atoms  ox  oy  oz
    let line3 = lines.next().context("missing atom count line")?;
    let fields3: Vec<&str> = line3.split_whitespace().collect();
    if fields3.len() < 4 { bail!("atom count line too short"); }
    let n_atoms: usize = fields3[0].parse::<i64>().context("n_atoms")?.unsigned_abs() as usize;
    let origin = Vector3::new(
        fields3[1].parse::<f64>().context("ox")? * BOHR_TO_ANG,
        fields3[2].parse::<f64>().context("oy")? * BOHR_TO_ANG,
        fields3[3].parse::<f64>().context("oz")? * BOHR_TO_ANG,
    );

    // 第4-6行：n_i  vx  vy  vz（体素定义）
    let mut shape = [0usize; 3];
    let mut spacing = Matrix3::<f64>::zeros();
    let mut cell_mat = Matrix3::<f64>::zeros();

    for i in 0..3 {
        let line = lines.next().with_context(|| format!("missing voxel line {i}"))?;
        let fs: Vec<f64> = line.split_whitespace()
            .map(|s| s.parse::<f64>())
            .collect::<std::result::Result<_, _>>()
            .with_context(|| format!("parsing voxel line {i}"))?;
        if fs.len() < 4 { bail!("voxel line {i} too short"); }
        let n = fs[0] as usize;
        let vx = fs[1] * BOHR_TO_ANG;
        let vy = fs[2] * BOHR_TO_ANG;
        let vz = fs[3] * BOHR_TO_ANG;
        shape[i] = n;
        // spacing 行 i = 体素向量 i
        spacing[(i, 0)] = vx;
        spacing[(i, 1)] = vy;
        spacing[(i, 2)] = vz;
        // cell 向量 i = n_i × 体素向量 i
        cell_mat[(i, 0)] = n as f64 * vx;
        cell_mat[(i, 1)] = n as f64 * vy;
        cell_mat[(i, 2)] = n as f64 * vz;
    }

    let cell = Cell::from_matrix(cell_mat);
    // pbc：cell 向量非零则视为周期
    let pbc = [
        cell_mat.row(0).norm() > 1e-10,
        cell_mat.row(1).norm() > 1e-10,
        cell_mat.row(2).norm() > 1e-10,
    ];

    // 原子行：Z  charge  x  y  z（位置单位 Bohr）
    let mut frame = Frame::with_cell(cell, pbc);
    for idx in 0..n_atoms {
        let line = lines.next()
            .with_context(|| format!("missing atom line {idx}"))?;
        let fs: Vec<f64> = line.split_whitespace()
            .map(|s| s.parse::<f64>())
            .collect::<std::result::Result<_, _>>()
            .with_context(|| format!("parsing atom line {idx}"))?;
        if fs.len() < 5 { bail!("atom line {idx} too short"); }
        let z = fs[0] as u8;
        let pos = Vector3::new(
            fs[2] * BOHR_TO_ANG,
            fs[3] * BOHR_TO_ANG,
            fs[4] * BOHR_TO_ANG,
        );
        // Z → 元素符号；未知 Z 退回 "X"
        let symbol = by_number(z)
            .map(|e| e.symbol.to_string())
            .unwrap_or_else(|| format!("X{z}"));
        frame.add_atom(Atom::new(symbol, pos));
    }

    // 剩余全部内容为体积数据（空白分隔的浮点数）
    let raw: Vec<f64> = lines
        .flat_map(|l| l.split_whitespace().map(|s| s.parse::<f64>().unwrap_or(0.0)))
        .collect();

    let expected = shape[0] * shape[1] * shape[2];
    if raw.len() < expected {
        bail!("volumetric data too short: got {}, expected {}", raw.len(), expected);
    }

    // C-order (X outer, Z inner) = Array3 默认行主序
    let data = Array3::from_shape_vec((shape[0], shape[1], shape[2]), raw[..expected].to_vec())
        .context("reshaping volumetric data")?;

    Ok(CubeData { frame, data, origin, spacing })
}

// ─── 测试 ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // 最小 cube：2×2×2 网格，2个原子（H2 分子）
    const CUBE_H2: &str = "\
H2 molecule test cube
OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z
  2   0.000000   0.000000   0.000000
  2   1.889726   0.000000   0.000000
  2   0.000000   1.889726   0.000000
  2   0.000000   0.000000   1.889726
  1   0.000000   0.000000   0.000000   0.000000
  1   0.000000   1.889726   0.000000   0.000000
 1.0e+00
 2.0e+00
 3.0e+00
 4.0e+00
 5.0e+00
 6.0e+00
 7.0e+00
 8.0e+00
";

    fn write_tmp(name: &str, content: &str) -> String {
        let p = std::env::temp_dir().join(name);
        std::fs::write(&p, content).unwrap();
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn test_read_atom_count() {
        let cd = parse_cube(CUBE_H2).unwrap();
        assert_eq!(cd.frame.n_atoms(), 2);
    }

    #[test]
    fn test_read_element_symbol() {
        let cd = parse_cube(CUBE_H2).unwrap();
        assert_eq!(cd.frame.atom(0).element, "H");
        assert_eq!(cd.frame.atom(1).element, "H");
    }

    #[test]
    fn test_read_positions_in_angstrom() {
        let cd = parse_cube(CUBE_H2).unwrap();
        // 原子 0 在原点
        let p0 = cd.frame.atom(0).position;
        assert!(p0.norm() < 1e-10, "atom0 should be at origin");
        // 原子 1：x = 1.889726 Bohr = 1.0 Å（H-H 键长）
        let p1 = cd.frame.atom(1).position;
        assert!((p1.x - 1.0).abs() < 1e-4, "atom1 x = {:.4} Å", p1.x);
    }

    #[test]
    fn test_read_cell_in_angstrom() {
        let cd = parse_cube(CUBE_H2).unwrap();
        // 体素 1.889726 Bohr × 2 = 2.0 Å
        let [a, b, c] = cd.frame.cell.as_ref().unwrap().lengths();
        assert!((a - 2.0).abs() < 1e-4, "a = {:.4}", a);
        assert!((b - 2.0).abs() < 1e-4, "b = {:.4}", b);
        assert!((c - 2.0).abs() < 1e-4, "c = {:.4}", c);
    }

    #[test]
    fn test_read_grid_shape() {
        let cd = parse_cube(CUBE_H2).unwrap();
        assert_eq!(cd.shape(), (2, 2, 2));
    }

    #[test]
    fn test_read_grid_values() {
        let cd = parse_cube(CUBE_H2).unwrap();
        // 值按 XYZ 顺序：1..8
        assert!((cd.data[[0,0,0]] - 1.0).abs() < 1e-10);
        assert!((cd.data[[0,0,1]] - 2.0).abs() < 1e-10);
        assert!((cd.data[[1,1,1]] - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_read_origin() {
        let cd = parse_cube(CUBE_H2).unwrap();
        assert!(cd.origin.norm() < 1e-10, "origin should be zero");
    }

    #[test]
    fn test_read_from_file() {
        let path = write_tmp("test_h2.cube", CUBE_H2);
        let cd = read_cube(&path).unwrap();
        assert_eq!(cd.frame.n_atoms(), 2);
    }
}
