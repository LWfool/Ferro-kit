//! Single-frame structure (analogous to ASE's `Atoms` object).

use nalgebra::{Matrix3, Vector3};

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
#[derive(Debug, Clone)]
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
    /// Stress tensor (eV/Å³), row-major, **positive = compression**.
    ///
    /// This is the sign CP2K, VASP and QE print, and the opposite of ASE's
    /// (positive = tension); readers of ASE-flavoured formats negate on the way
    /// in. The virial follows from it without a sign change: `virial = stress *
    /// V` (eV), which is what DeePMD, QUIP and GPUMD call `virial`.
    ///
    /// Voigt order is left to the caller — nothing in this crate flattens it.
    pub stress: Option<Matrix3<f64>>,
    /// 每个原子的速度（Å/fs），顺序与 atoms 一致
    pub velocities: Option<Vec<Vector3<f64>>>,
    /// Instantaneous ionic temperature (K), as an MD engine reports it per step.
    ///
    /// `None` for everything that is not an MD frame — a relaxation, a POSCAR, a
    /// CIF. It lives here rather than in a parallel array because it is a
    /// per-frame quantity: `select`/`stride` must carry it along, and a side
    /// channel would silently shift out of step with the frames.
    ///
    /// CP2K and VASP OUTCAR print it directly; vasprun.xml does not, so that
    /// reader derives it from the ionic kinetic energy and says so.
    pub temperature: Option<f64>,
    /// The MD step number the engine gave this frame, when it prints one.
    ///
    /// Kept because `collect` stitches restart segments without de-duplicating
    /// them: overlapping frames are genuinely identical, and the step number is
    /// what makes that overlap visible instead of silently doubling a run.
    /// CP2K and VASP OUTCAR number their ionic steps; vasprun.xml does not.
    pub step: Option<i64>,
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
            temperature: None,
            step: None,
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

    /// 按首次出现顺序返回去重后的元素符号（用于 KIND/类型 ID/POSCAR 元素行等）。
    pub fn unique_elements(&self) -> Vec<String> {
        let mut seen = std::collections::HashSet::new();
        let mut out = Vec::new();
        for a in &self.atoms {
            if seen.insert(a.element.as_str()) {
                out.push(a.element.clone());
            }
        }
        out
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

    /// 把周期方向上越出盒子的原子折回 `[0, 1)`；非周期方向、盒内原子都不动。
    ///
    /// 越界的原子平移整数个晶格矢量，不做「分数 → 笛卡尔」往返：往返会给每个原子
    /// 带来 1e-15 级舍入，盒内原子的坐标末位也跟着变。判越界留 `WRAP_TOL` 的余量，
    /// 坐标正好在 0 而算出分数 -1e-17 的原子不该被搬到盒子另一头。
    /// 无 cell 或晶胞奇异时不动。
    pub fn wrap_all(&mut self) {
        const WRAP_TOL: f64 = 1e-10;
        let Some(cell) = self.cell.as_ref() else { return };
        let pbc = self.pbc;
        for atom in &mut self.atoms {
            let Ok(f) = cell.cartesian_to_fractional(atom.position) else { return };
            let shift = Vector3::from_fn(|i, _| {
                let out = f[i] < -WRAP_TOL || f[i] >= 1.0 + WRAP_TOL;
                if pbc[i] && out { -f[i].floor() } else { 0.0 }
            });
            if shift != Vector3::zeros() {
                atom.position += cell.fractional_to_cartesian(shift);
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
    fn test_unique_elements_first_appearance() {
        let mut f = Frame::new();
        for e in ["O", "H", "O", "Fe", "H", "O"] {
            f.add_atom(Atom::new(e, Vector3::zeros()));
        }
        // 去重并保持首次出现顺序
        assert_eq!(f.unique_elements(), vec!["O", "H", "Fe"]);
        assert!(Frame::new().unique_elements().is_empty());
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

    #[test]
    fn test_wrap_all_respects_pbc_per_axis_and_leaves_inside_atoms_untouched() {
        // 三斜，z 非周期（slab）
        let cell = Cell::from_lengths_angles(9.0, 11.0, 8.0, 80.0, 105.0, 70.0).unwrap();
        let mut f = Frame::with_cell(cell.clone(), [true, true, false]);
        let inside = cell.fractional_to_cartesian(Vector3::new(0.3, 0.7, 0.2)) + Vector3::new(1e-13, 0.0, 0.0);
        let at_zero = Vector3::new(0.0, 0.0, 0.0);
        let out = cell.fractional_to_cartesian(Vector3::new(1.25, -0.4, 1.6));
        for p in [inside, at_zero, out] { f.add_atom(Atom::new("Si", p)); }
        f.wrap_all();
        // 盒内原子逐位不动（不做分数往返）
        assert_eq!(f.atom(0).position, inside);
        assert_eq!(f.atom(1).position, at_zero);
        // 越界原子：x、y 折回，z 非周期保持 1.6
        let g = cell.cartesian_to_fractional(f.atom(2).position).unwrap();
        for (k, w) in [0.25, 0.6, 1.6].iter().enumerate() {
            assert!((g[k] - w).abs() < 1e-12, "第 {k} 轴分数坐标 {}，应为 {w}", g[k]);
        }
        // 位移恰为整数个晶格矢量
        let d = cell.cartesian_to_fractional(f.atom(2).position - out).unwrap();
        assert!((d - Vector3::new(-1.0, 1.0, 0.0)).norm() < 1e-12, "{d:?}");
    }
}
