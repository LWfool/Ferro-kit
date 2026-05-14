//! Gaussian cube format reader.
//!
//! Specification: <https://paulbourke.net/dataformats/cube/>
//!
//! File structure:
//!   line 1   comment
//!   line 2   comment
//!   line 3   n_atoms  ox  oy  oz        (origin, Bohr)
//!   line 4-6 n_i  vx  vy  vz            (voxel count + step vector, Bohr)
//!   atom rows  Z  charge  x  y  z       (Bohr)
//!   volumetric data (i1 outer/slowest, i3 inner/fastest)
//!
//! Two public entry points:
//!   - `read_cube`       → `CubeData`         (visualisation / density maps)
//!   - `read_cube_as_chg` → `(Frame, ChargeGrid)` (Bader charge analysis)

use nalgebra::{Matrix3, Vector3};
use anyhow::{bail, Context, Result};
use ferro_core::{Atom, Cell, ChargeGrid, CubeData, Frame};
use ferro_core::units::BOHR_TO_ANG;
use ferro_core::data::elements::by_number;

// ─── 内部中间结构 ─────────────────────────────────────────────────────────────

struct ParsedCube {
    /// 原点（Bohr，绝对 Cartesian）
    origin_bohr: Vector3<f64>,
    /// 网格尺寸 [N1, N2, N3]
    shape: [usize; 3],
    /// 体素步向量（Bohr），step_bohr[i] = v_i
    step_bohr: [[f64; 3]; 3],
    /// 体素步向量（Å），step_ang[i] = v_i × BOHR_TO_ANG
    step_ang: [[f64; 3]; 3],
    /// 原子列表（位置已转换为 Å）
    frame: Frame,
    /// 原始密度值，cube 存储顺序：i1 最慢，i3 最快（i3-fastest）
    /// flat index = i1*(N2*N3) + i2*N3 + i3
    raw: Vec<f64>,
}

fn parse_header(content: &str) -> Result<ParsedCube> {
    let mut lines = content.lines();
    let mut next = |what: &str| -> Result<&str> {
        lines.next().with_context(|| format!("unexpected EOF before {what}"))
    };

    // 两行注释
    next("comment 1")?;
    next("comment 2")?;

    // 第 3 行：n_atoms  ox  oy  oz（Bohr）
    let f3 = parse_floats(next("atom count line")?, 4)?;
    let n_atoms = f3[0].abs() as usize;
    // 原点（Bohr）
    let origin_bohr = Vector3::new(f3[1], f3[2], f3[3]);

    // 第 4-6 行：N_i  v_ix  v_iy  v_iz（Bohr）
    let mut shape = [0usize; 3];
    let mut step_bohr = [[0.0f64; 3]; 3];
    let mut step_ang  = [[0.0f64; 3]; 3];

    for i in 0..3 {
        let fs = parse_floats(next(&format!("grid line {i}"))?, 4)?;
        shape[i] = fs[0] as usize;
        step_bohr[i] = [fs[1], fs[2], fs[3]];
        step_ang[i]  = [fs[1] * BOHR_TO_ANG, fs[2] * BOHR_TO_ANG, fs[3] * BOHR_TO_ANG];
    }

    // 原子行：Z  charge  x  y  z（Bohr）
    // 单元格用 Å 的步向量构建
    let cell_mat = Matrix3::new(
        shape[0] as f64 * step_ang[0][0], shape[0] as f64 * step_ang[0][1], shape[0] as f64 * step_ang[0][2],
        shape[1] as f64 * step_ang[1][0], shape[1] as f64 * step_ang[1][1], shape[1] as f64 * step_ang[1][2],
        shape[2] as f64 * step_ang[2][0], shape[2] as f64 * step_ang[2][1], shape[2] as f64 * step_ang[2][2],
    );
    let cell = Cell::from_matrix(cell_mat);
    let pbc = [
        cell_mat.row(0).norm() > 1e-10,
        cell_mat.row(1).norm() > 1e-10,
        cell_mat.row(2).norm() > 1e-10,
    ];
    let mut frame = Frame::with_cell(cell, pbc);

    // 原子位置：绝对 Bohr → Å，减去原点使坐标相对于格点原点
    let origin_ang = origin_bohr * BOHR_TO_ANG;
    for idx in 0..n_atoms {
        let fs = parse_floats(next(&format!("atom line {idx}"))?, 5)?;
        let z = fs[0] as u8;
        let pos = Vector3::new(
            fs[2] * BOHR_TO_ANG - origin_ang.x,
            fs[3] * BOHR_TO_ANG - origin_ang.y,
            fs[4] * BOHR_TO_ANG - origin_ang.z,
        );
        let symbol = by_number(z)
            .map(|e| e.symbol.to_string())
            .unwrap_or_else(|| format!("X{z}"));
        frame.add_atom(Atom::new(symbol, pos));
    }

    // 读取剩余所有密度值
    let nrho = shape[0] * shape[1] * shape[2];
    let raw: Vec<f64> = lines
        .flat_map(|l| l.split_whitespace().map(|s| s.parse::<f64>().unwrap_or(0.0)))
        .collect();
    if raw.len() < nrho {
        bail!("volumetric data too short: got {}, expected {}", raw.len(), nrho);
    }
    let raw = raw[..nrho].to_vec();

    Ok(ParsedCube { origin_bohr, shape, step_bohr, step_ang, frame, raw })
}

// ─── 公共入口：CubeData（可视化 / 密度图）────────────────────────────────────

/// Read a Gaussian cube file into `CubeData` (for visualisation and density maps).
pub fn read_cube(path: &str) -> Result<CubeData> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path}"))?;
    parse_cube(&content).with_context(|| format!("parsing {path}"))
}

fn parse_cube(content: &str) -> Result<CubeData> {
    let pc = parse_header(content)?;
    let spacing = Matrix3::new(
        pc.step_ang[0][0], pc.step_ang[0][1], pc.step_ang[0][2],
        pc.step_ang[1][0], pc.step_ang[1][1], pc.step_ang[1][2],
        pc.step_ang[2][0], pc.step_ang[2][1], pc.step_ang[2][2],
    );
    // CubeData 保留 cube 原始数据顺序（i3 最快）和原始数值（不缩放）
    Ok(CubeData {
        frame: pc.frame,
        data: pc.raw,
        shape: pc.shape,
        origin: pc.origin_bohr * BOHR_TO_ANG,
        spacing,
    })
}

// ─── 公共入口：(Frame, ChargeGrid)（Bader 电荷分析）──────────────────────────

/// Read a Gaussian cube file into `(Frame, ChargeGrid)` for Bader charge analysis.
///
/// Handles:
/// - Unit conversion: positions Bohr → Å; density e/Bohr³ → ChargeGrid convention.
/// - Index transposition: cube (i1 slowest, i3 fastest) → ChargeGrid (n1 fastest, n3 slowest).
/// - Non-zero origin: atom positions are shifted to be relative to the grid origin.
///
/// Density convention: the cube file must store physical charge density in e/Bohr³,
/// which is the standard used by Quantum ESPRESSO `pp.x` and Gaussian.
/// ChargeGrid internally stores `rho_stored = ρ_phys × V_cell`, so the conversion is:
/// `rho_stored = cube_value × V_cell_Bohr` (unit-independent invariant).
pub fn read_cube_as_chg(path: &str) -> Result<(Frame, ChargeGrid)> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path}"))?;
    parse_cube_as_chg(&content).with_context(|| format!("parsing {path}"))
}

fn parse_cube_as_chg(content: &str) -> Result<(Frame, ChargeGrid)> {
    let pc = parse_header(content)?;
    let [n1, n2, n3] = pc.shape;
    let nrho = n1 * n2 * n3;

    // 步向量矩阵（Bohr）：行 i = v_i
    let vm = Matrix3::new(
        pc.step_bohr[0][0], pc.step_bohr[0][1], pc.step_bohr[0][2],
        pc.step_bohr[1][0], pc.step_bohr[1][1], pc.step_bohr[1][2],
        pc.step_bohr[2][0], pc.step_bohr[2][1], pc.step_bohr[2][2],
    );
    // V_cell_Bohr = N1*N2*N3 × |det(v1, v2, v3)|
    let v_voxel_bohr = vm.determinant().abs();
    let v_cell_bohr = nrho as f64 * v_voxel_bohr;

    // 建 Å 单元格（ChargeGrid 要求）
    let cell_mat = Matrix3::new(
        n1 as f64 * pc.step_ang[0][0], n1 as f64 * pc.step_ang[0][1], n1 as f64 * pc.step_ang[0][2],
        n2 as f64 * pc.step_ang[1][0], n2 as f64 * pc.step_ang[1][1], n2 as f64 * pc.step_ang[1][2],
        n3 as f64 * pc.step_ang[2][0], n3 as f64 * pc.step_ang[2][1], n3 as f64 * pc.step_ang[2][2],
    );
    let cell = Cell::from_matrix(cell_mat);

    // 索引转置 + 密度缩放
    //
    // cube 存储顺序：flat = i1*(N2*N3) + i2*N3 + i3  （i3 最快）
    // ChargeGrid 顺序：flat = n1 + N1*n2 + N1*N2*n3  （n1 最快，n3 最慢）
    //
    // 映射：cube i1 → ChargeGrid n1（同为第一轴方向）
    //       cube i2 → ChargeGrid n2
    //       cube i3 → ChargeGrid n3
    //
    // 缩放：rho_stored = ρ_cube × V_cell_Bohr
    let mut rho = vec![0.0_f64; nrho];
    for i1 in 0..n1 {
        for i2 in 0..n2 {
            for i3 in 0..n3 {
                let cube_flat = i1 * n2 * n3 + i2 * n3 + i3;
                let chg_flat  = i1 + n1 * (i2 + n2 * i3);
                rho[chg_flat] = pc.raw[cube_flat] * v_cell_bohr;
            }
        }
    }

    let chg = ChargeGrid::new(rho, pc.shape, &cell);
    Ok((pc.frame, chg))
}

// ─── 内部辅助 ─────────────────────────────────────────────────────────────────

fn parse_floats(line: &str, min: usize) -> Result<Vec<f64>> {
    let v: Vec<f64> = line.split_whitespace()
        .map_while(|s| s.parse::<f64>().ok())
        .collect();
    if v.len() < min {
        bail!("expected ≥{min} numbers on line '{line}', got {}", v.len());
    }
    Ok(v)
}

// ─── 测试 ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // 2×2×2 格，2 个 H 原子，原点在 (0,0,0)
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

    // 非零原点（1×1×1 Bohr 偏移）
    const CUBE_OFFSET: &str = "\
Offset origin test
comment
  1   1.000000   0.000000   0.000000
  2   1.889726   0.000000   0.000000
  2   0.000000   1.889726   0.000000
  2   0.000000   0.000000   1.889726
  1   0.000000   2.889726   0.000000   0.000000
 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0
";

    fn write_tmp(name: &str, content: &str) -> String {
        let p = std::env::temp_dir().join(name);
        std::fs::write(&p, content).unwrap();
        p.to_str().unwrap().to_string()
    }

    // ── CubeData tests ────────────────────────────────────────────────────────

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
    fn test_read_positions_angstrom() {
        let cd = parse_cube(CUBE_H2).unwrap();
        let p0 = cd.frame.atom(0).position;
        assert!(p0.norm() < 1e-10, "atom0 should be at origin");
        // atom1: x = 1.889726 Bohr = 1.0 Å
        let p1 = cd.frame.atom(1).position;
        assert!((p1.x - 1.0).abs() < 1e-4, "atom1 x = {:.4} Å", p1.x);
    }

    #[test]
    fn test_read_cell_angstrom() {
        let cd = parse_cube(CUBE_H2).unwrap();
        let [a, b, c] = cd.frame.cell.as_ref().unwrap().lengths();
        // 2 × 1.889726 Bohr = 2.0 Å
        assert!((a - 2.0).abs() < 1e-4);
        assert!((b - 2.0).abs() < 1e-4);
        assert!((c - 2.0).abs() < 1e-4);
    }

    #[test]
    fn test_read_grid_shape() {
        let cd = parse_cube(CUBE_H2).unwrap();
        assert_eq!(cd.shape(), (2, 2, 2));
    }

    #[test]
    fn test_read_grid_values_cube_order() {
        // CubeData 保持原始 cube 顺序：i3 最快
        // get(ix, iy, iz) → ix*(ny*nz)+iy*nz+iz
        let cd = parse_cube(CUBE_H2).unwrap();
        // cube file: 1 2 3 4 5 6 7 8
        // (0,0,0)=1, (0,0,1)=2, (0,1,0)=3, ...
        assert!((cd.get(0, 0, 0) - 1.0).abs() < 1e-10);
        assert!((cd.get(0, 0, 1) - 2.0).abs() < 1e-10);
        assert!((cd.get(1, 1, 1) - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_read_from_file() {
        let path = write_tmp("test_h2.cube", CUBE_H2);
        let cd = read_cube(&path).unwrap();
        assert_eq!(cd.frame.n_atoms(), 2);
    }

    // ── ChargeGrid tests ──────────────────────────────────────────────────────

    #[test]
    fn test_chg_grid_shape() {
        let (_, chg) = parse_cube_as_chg(CUBE_H2).unwrap();
        assert_eq!(chg.shape, [2, 2, 2]);
        assert_eq!(chg.nrho, 8);
    }

    #[test]
    fn test_chg_grid_cell_angstrom() {
        let (frame, _) = parse_cube_as_chg(CUBE_H2).unwrap();
        let [a, b, c] = frame.cell.as_ref().unwrap().lengths();
        assert!((a - 2.0).abs() < 1e-4);
        assert!((b - 2.0).abs() < 1e-4);
        assert!((c - 2.0).abs() < 1e-4);
    }

    #[test]
    fn test_chg_grid_index_transposition() {
        // cube 顺序：i3 最快 → ChargeGrid 顺序：n1 最快
        // cube 2×2×2，值 1..8：cube(0,0,0)=1, cube(0,0,1)=2, cube(0,1,0)=3...
        // ChargeGrid: rho[n1 + 2*n2 + 4*n3]
        // n1=0,n2=0,n3=0 → cube(0,0,0)=1 → chg_flat=0
        // n1=1,n2=0,n3=0 → cube(1,0,0)=5 → chg_flat=1
        // n1=0,n2=0,n3=1 → cube(0,0,1)=2 → chg_flat=4  (flat = 0+0+4)
        let (_, chg) = parse_cube_as_chg(CUBE_H2).unwrap();
        // V_cell_bohr 用 2×1.889726 步向量立方体计算
        let v_vox_bohr = 1.889726_f64.powi(3); // 只适用于正交情形
        let v_cell_bohr = 8.0 * v_vox_bohr;
        let scale = v_cell_bohr;

        // ChargeGrid flat 0 = (n1=0,n2=0,n3=0) → cube(i1=0,i2=0,i3=0) = raw[0] = 1.0
        assert!((chg.rho[0] - 1.0 * scale).abs() < 1e-6, "rho[0] = {}", chg.rho[0]);
        // ChargeGrid flat 1 = (n1=1,n2=0,n3=0) → cube(i1=1,i2=0,i3=0) = raw[4] = 5.0
        assert!((chg.rho[1] - 5.0 * scale).abs() < 1e-6, "rho[1] = {}", chg.rho[1]);
        // ChargeGrid flat 4 = (n1=0,n2=0,n3=1) → cube(i1=0,i2=0,i3=1) = raw[1] = 2.0
        assert!((chg.rho[4] - 2.0 * scale).abs() < 1e-6, "rho[4] = {}", chg.rho[4]);
    }

    #[test]
    fn test_chg_charge_conservation() {
        // Σ rho[i] / nrho 应等于 Σ density × V_voxel_bohr
        // = mean_density × V_cell_bohr = mean_density × nrho × V_vox
        // 原始值 1..8，均值 4.5
        let (_, chg) = parse_cube_as_chg(CUBE_H2).unwrap();
        let total_chg: f64 = chg.rho.iter().sum::<f64>() / chg.nrho as f64;
        let v_vox_bohr = 1.889726_f64.powi(3);
        let v_cell_bohr = 8.0 * v_vox_bohr;
        let expected = 4.5 * v_cell_bohr; // mean density × V_cell
        assert!((total_chg - expected).abs() / expected < 1e-6,
            "total_chg = {total_chg:.6}, expected = {expected:.6}");
    }

    #[test]
    fn test_chg_offset_origin() {
        // 非零原点时，原子位置应相对于格点原点
        let (frame, _) = parse_cube_as_chg(CUBE_OFFSET).unwrap();
        // atom at x=2.889726 Bohr (absolute) − origin.x=1.0 Bohr = 1.889726 Bohr = 1.0 Å
        let p = frame.atom(0).position;
        assert!((p.x - 1.0).abs() < 1e-4, "p.x = {:.4}", p.x);
        assert!((p.y).abs() < 1e-4, "p.y = {:.4}", p.y);
    }

    #[test]
    fn test_chg_from_file() {
        let path = write_tmp("test_h2_chg.cube", CUBE_H2);
        let (frame, chg) = read_cube_as_chg(&path).unwrap();
        assert_eq!(frame.n_atoms(), 2);
        assert_eq!(chg.shape, [2, 2, 2]);
    }
}
