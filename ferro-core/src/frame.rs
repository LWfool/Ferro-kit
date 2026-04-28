//! Single-frame structure (analogous to ASE's `Atoms` object).

use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

use crate::atom::Atom;
use crate::cell::Cell;

/// 单帧结构，是 ferro 中的核心结构单元。
///
/// - 非周期性体系（分子）：`cell = None`，`pbc = [false; 3]`
/// - 三维周期性体系（晶体）：`cell = Some(…)`，`pbc = [true; 3]`
/// - 二维周期性体系（表面）：`cell = Some(…)`，`pbc = [true, true, false]`
///
/// NPT 模拟轨迹中每帧盒子不同，由 [`crate::trajectory::Trajectory`] 中各帧
/// 各自持有自己的 `cell` 自然处理，无需特殊设计。
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Frame {
    /// 原子列表；原子的唯一标识是其在此 Vec 中的下标
    pub atoms: Vec<Atom>,
    /// 晶格；`None` 表示非周期性体系
    pub cell: Option<Cell>,
    /// 各方向是否周期性
    pub pbc: [bool; 3],

    // ── 体系属性 ────────────────────────────────────────────────────────────
    /// 体系总电荷（单位 e）
    pub charge: i32,
    /// 自旋多重度 2S+1；未成对电子数 = multiplicity − 1
    pub multiplicity: u32,

    // ── 可选键连接（需要时再填充）───────────────────────────────────────────
    /// 键列表 (i, j)，i/j 为 atoms 中的下标
    pub bonds: Option<Vec<(usize, usize)>>,

    // ── 计算结果（后处理写回）───────────────────────────────────────────────
    /// 体系总能量（eV）
    pub energy: Option<f64>,
    /// 每个原子的受力（eV/Å），顺序与 atoms 一致
    pub forces: Option<Vec<Vector3<f64>>>,
    /// 应力张量（eV/Å³），Voigt 顺序由调用方自行约定
    pub stress: Option<Matrix3<f64>>,
    /// 每个原子的速度（Å/fs），顺序与 atoms 一致
    pub velocities: Option<Vec<Vector3<f64>>>,
}

impl Frame {
    /// 创建空帧（非周期性，中性单重态）。
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            cell: None,
            pbc: [false; 3],
            charge: 0,
            multiplicity: 1,
            bonds: None,
            energy: None,
            forces: None,
            stress: None,
            velocities: None,
        }
    }

    /// 创建周期性帧（晶体/表面常用入口）。
    pub fn with_cell(cell: Cell, pbc: [bool; 3]) -> Self {
        Self {
            cell: Some(cell),
            pbc,
            ..Self::new()
        }
    }

    // ── 原子访问 ─────────────────────────────────────────────────────────────

    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn atom(&self, index: usize) -> &Atom {
        &self.atoms[index]
    }

    pub fn atom_mut(&mut self, index: usize) -> &mut Atom {
        &mut self.atoms[index]
    }

    pub fn add_atom(&mut self, atom: Atom) {
        self.atoms.push(atom);
    }

    /// 返回所有原子的元素符号切片引用。
    pub fn symbols(&self) -> Vec<&str> {
        self.atoms.iter().map(|a| a.element.as_str()).collect()
    }

    /// 返回某元素在此帧中的原子数量。
    pub fn count_element(&self, element: &str) -> usize {
        self.atoms.iter().filter(|a| a.element == element).count()
    }

    // ── 几何量 ───────────────────────────────────────────────────────────────

    /// 总质量（amu）。
    pub fn total_mass(&self) -> f64 {
        self.atoms.iter().map(|a| a.effective_mass()).sum()
    }

    /// 质心（Å），按质量加权。
    pub fn center_of_mass(&self) -> Vector3<f64> {
        if self.atoms.is_empty() {
            return Vector3::zeros();
        }
        let total = self.total_mass();
        self.atoms
            .iter()
            .map(|a| a.position * a.effective_mass())
            .fold(Vector3::zeros(), |acc, v| acc + v)
            / total
    }

    /// 几何中心（Å），不加权。
    pub fn geometric_center(&self) -> Vector3<f64> {
        if self.atoms.is_empty() {
            return Vector3::zeros();
        }
        let sum = self.atoms
            .iter()
            .map(|a| a.position)
            .fold(Vector3::zeros(), |acc, v| acc + v);
        sum / self.atoms.len() as f64
    }

    // ── 变换 ─────────────────────────────────────────────────────────────────

    /// 平移所有原子。
    pub fn translate(&mut self, disp: Vector3<f64>) {
        for atom in &mut self.atoms {
            atom.position += disp;
        }
    }

    /// 将质心移到原点。
    pub fn center(&mut self) {
        let com = self.center_of_mass();
        self.translate(-com);
    }

    // ── 周期性 ───────────────────────────────────────────────────────────────

    /// 任意方向是否周期性。
    pub fn is_periodic(&self) -> bool {
        self.pbc.iter().any(|&p| p)
    }

    /// 将所有原子位置折叠回盒子内（需要 cell 存在且对应方向周期性）。
    pub fn wrap_all(&mut self) {
        if let Some(cell) = &self.cell {
            let cell = cell.clone();
            for atom in &mut self.atoms {
                atom.position = cell.wrap_position(atom.position);
            }
        }
    }
}

impl Default for Frame {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::Cell;

    #[test]
    fn test_new_frame() {
        let f = Frame::new();
        assert_eq!(f.n_atoms(), 0);
        assert!(!f.is_periodic());
        assert_eq!(f.charge, 0);
        assert_eq!(f.multiplicity, 1);
    }

    #[test]
    fn test_add_and_count() {
        let mut f = Frame::new();
        f.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        f.add_atom(Atom::new("O",  Vector3::new(2.0, 0.0, 0.0)));
        f.add_atom(Atom::new("O",  Vector3::new(0.0, 2.0, 0.0)));
        assert_eq!(f.n_atoms(), 3);
        assert_eq!(f.count_element("O"), 2);
    }

    #[test]
    fn test_center_of_mass() {
        let mut f = Frame::new();
        f.add_atom(Atom::new("C", Vector3::new(0.0, 0.0, 0.0)));
        f.add_atom(Atom::new("C", Vector3::new(2.0, 0.0, 0.0)));
        let com = f.center_of_mass();
        assert!((com.x - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_periodic_frame() {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let f = Frame::with_cell(cell, [true; 3]);
        assert!(f.is_periodic());
        assert!(f.cell.is_some());
    }

    #[test]
    fn test_wrap_all() {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut f = Frame::with_cell(cell, [true; 3]);
        f.add_atom(Atom::new("Fe", Vector3::new(11.0, -1.0, 5.0)));
        f.wrap_all();
        let pos = f.atom(0).position;
        assert!(pos.x >= 0.0 && pos.x < 10.0);
        assert!(pos.y >= 0.0 && pos.y < 10.0);
    }
}
