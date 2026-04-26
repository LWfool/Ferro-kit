//! Lattice (periodic simulation box).

use nalgebra::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};
use crate::error::{ChemError, Result};

/// 周期性晶格，用行向量存储三个晶格矢量 a、b、c。
///
/// ```text
/// matrix = [a_x  a_y  a_z]   ← 行 0 = 晶格矢量 a
///          [b_x  b_y  b_z]   ← 行 1 = 晶格矢量 b
///          [c_x  c_y  c_z]   ← 行 2 = 晶格矢量 c
/// ```
///
/// Cartesian ↔ 分数坐标转换关系：
/// `cart = matrix.transpose() * frac`
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cell {
    pub matrix: Matrix3<f64>,
}

impl Cell {
    /// 直接由 3×3 矩阵构造（行 = 晶格矢量）。
    pub fn from_matrix(matrix: Matrix3<f64>) -> Self {
        Self { matrix }
    }

    /// 由晶格常数 (a, b, c, α, β, γ) 构造。
    ///
    /// 约定：a 沿 x 轴，b 在 xy 平面，c 由右手定则确定。
    /// 角度单位：度。
    pub fn from_lengths_angles(
        a: f64, b: f64, c: f64,
        alpha: f64, beta: f64, gamma: f64,
    ) -> Result<Self> {
        if a <= 0.0 || b <= 0.0 || c <= 0.0 {
            return Err(ChemError::ValidationError(
                "Lattice constants must be positive".into(),
            ));
        }
        let (al, be, ga) = (
            alpha.to_radians(),
            beta.to_radians(),
            gamma.to_radians(),
        );
        let cos_a = al.cos();
        let cos_b = be.cos();
        let cos_g = ga.cos();
        let sin_g = ga.sin();

        let cx = cos_b;
        let cy = (cos_a - cos_b * cos_g) / sin_g;
        let cz_sq = 1.0 - cx * cx - cy * cy;
        if cz_sq < 0.0 {
            return Err(ChemError::ValidationError(
                "Invalid lattice angles: cannot form a valid cell".into(),
            ));
        }
        let cz = cz_sq.sqrt();

        // 行向量：a、b、c
        let matrix = Matrix3::new(
            a,              0.0,        0.0,
            b * cos_g,      b * sin_g,  0.0,
            c * cx,         c * cy,     c * cz,
        );
        Ok(Self { matrix })
    }

    // ── 查询 ──────────────────────────────────────────────────────────────────

    /// 三个晶格矢量的长度 [a, b, c]（Å）。
    pub fn lengths(&self) -> [f64; 3] {
        [
            self.matrix.row(0).norm(),
            self.matrix.row(1).norm(),
            self.matrix.row(2).norm(),
        ]
    }

    /// 三个晶格角 [α, β, γ]（度）：α = ∠(b,c)，β = ∠(a,c)，γ = ∠(a,b)。
    pub fn angles(&self) -> [f64; 3] {
        let row = |i: usize| Vector3::new(
            self.matrix[(i, 0)],
            self.matrix[(i, 1)],
            self.matrix[(i, 2)],
        );
        let a = row(0);
        let b = row(1);
        let c = row(2);
        let angle_between = |u: Vector3<f64>, v: Vector3<f64>| {
            (u.dot(&v) / (u.norm() * v.norm())).clamp(-1.0, 1.0).acos().to_degrees()
        };
        [angle_between(b, c), angle_between(a, c), angle_between(a, b)]
    }

    /// 晶胞体积（Å³）。
    pub fn volume(&self) -> f64 {
        self.matrix.determinant().abs()
    }

    // ── 坐标转换 ──────────────────────────────────────────────────────────────

    /// 分数坐标 → Cartesian 坐标（Å）。
    pub fn fractional_to_cartesian(&self, frac: Vector3<f64>) -> Vector3<f64> {
        self.matrix.transpose() * frac
    }

    /// Cartesian 坐标（Å）→ 分数坐标。
    pub fn cartesian_to_fractional(&self, cart: Vector3<f64>) -> Vector3<f64> {
        self.matrix
            .transpose()
            .try_inverse()
            .expect("Cell matrix is singular")
            * cart
    }

    // ── 周期性处理 ────────────────────────────────────────────────────────────

    /// 将 Cartesian 坐标折叠回 [0, 1) 的分数坐标范围内（周期性包裹）。
    pub fn wrap_position(&self, cart: Vector3<f64>) -> Vector3<f64> {
        let frac = self.cartesian_to_fractional(cart);
        let wrapped = Vector3::new(
            frac.x - frac.x.floor(),
            frac.y - frac.y.floor(),
            frac.z - frac.z.floor(),
        );
        self.fractional_to_cartesian(wrapped)
    }

    /// 对位移向量应用最小镜像约定，使每个分量落在 [-0.5, 0.5)。
    ///
    /// 用法示例：计算两原子的周期性最短距离
    /// ```ignore
    /// let diff = atom_b.position - atom_a.position;
    /// let mic_diff = cell.minimum_image(diff);
    /// let dist = mic_diff.norm();
    /// ```
    pub fn minimum_image(&self, diff: Vector3<f64>) -> Vector3<f64> {
        let frac = self.cartesian_to_fractional(diff);
        let mic_frac = Vector3::new(
            frac.x - frac.x.round(),
            frac.y - frac.y.round(),
            frac.z - frac.z.round(),
        );
        self.fractional_to_cartesian(mic_frac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cubic(a: f64) -> Cell {
        Cell::from_lengths_angles(a, a, a, 90.0, 90.0, 90.0).unwrap()
    }

    #[test]
    fn test_cubic_volume() {
        let cell = cubic(10.0);
        assert!((cell.volume() - 1000.0).abs() < 1e-8);
    }

    #[test]
    fn test_cubic_lengths_angles() {
        let cell = cubic(5.0);
        let [a, b, c] = cell.lengths();
        assert!((a - 5.0).abs() < 1e-10 && (b - 5.0).abs() < 1e-10 && (c - 5.0).abs() < 1e-10);
        let [al, be, ga] = cell.angles();
        assert!((al - 90.0).abs() < 1e-8 && (be - 90.0).abs() < 1e-8 && (ga - 90.0).abs() < 1e-8);
    }

    #[test]
    fn test_frac_cart_roundtrip() {
        let cell = Cell::from_lengths_angles(4.0, 5.0, 6.0, 80.0, 90.0, 100.0).unwrap();
        let frac = Vector3::new(0.3, 0.5, 0.7);
        let cart = cell.fractional_to_cartesian(frac);
        let back = cell.cartesian_to_fractional(cart);
        assert!((back - frac).norm() < 1e-10);
    }

    #[test]
    fn test_wrap_position() {
        let cell = cubic(10.0);
        let pos = Vector3::new(12.0, -1.0, 5.0);
        let wrapped = cell.wrap_position(pos);
        assert!(wrapped.x >= 0.0 && wrapped.x < 10.0);
        assert!(wrapped.y >= 0.0 && wrapped.y < 10.0);
    }

    #[test]
    fn test_minimum_image() {
        let cell = cubic(10.0);
        // 距离 9 Å，MIC 应返回 -1 Å
        let diff = Vector3::new(9.0, 0.0, 0.0);
        let mic = cell.minimum_image(diff);
        assert!((mic.x + 1.0).abs() < 1e-10);
    }
}
