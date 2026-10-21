//! Lattice (periodic simulation box).

use nalgebra::{Matrix3, Vector3};
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
#[derive(Debug, Clone)]
pub struct Cell {
    pub matrix: Matrix3<f64>,
    /// 盒子原点（LAMMPS 的 lo），构造器一律填 0。只有认原点的格式（LAMMPS dump / data）
    /// 读写它；分数坐标、wrap、最小镜像与全部分析都**不用**它（同 ASE 的 `celldisp`）。
    /// 读入时把坐标平移到原点 0 会扰动 NPT 的 MSD，已否决（审查 M5），故原点作为数据带着走
    pub origin: Vector3<f64>,
}

impl Cell {
    /// 直接由 3×3 矩阵构造（行 = 晶格矢量）。
    pub fn from_matrix(matrix: Matrix3<f64>) -> Self {
        Self { matrix, origin: Vector3::zeros() }
    }

    /// 由晶格常数 (a, b, c, α, β, γ) 构造。
    ///
    /// 约定：a 沿 x 轴，b 在 xy 平面，c 由右手定则确定。
    /// 角度单位：度。
    pub fn from_lengths_angles(
        a: f64, b: f64, c: f64,
        alpha: f64, beta: f64, gamma: f64,
    ) -> Result<Self> {
        // Written as "valid or reject" so NaN, which fails every comparison, is rejected
        let length_ok = |v: f64| v > 0.0 && v.is_finite();
        if !(length_ok(a) && length_ok(b) && length_ok(c)) {
            return Err(ChemError::ValidationError(
                "Lattice constants must be finite and positive".into(),
            ));
        }
        let angle_ok = |v: f64| v > 0.0 && v < 180.0;
        if !(angle_ok(alpha) && angle_ok(beta) && angle_ok(gamma)) {
            return Err(ChemError::ValidationError(
                "Lattice angles must lie strictly between 0 and 180 degrees".into(),
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
        Ok(Self { matrix, origin: Vector3::zeros() })
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

    /// The volume of a cell that can serve as a 3-D box: all nine numbers finite and
    /// a volume above 1e-12 Å³; `what` names the cell in the error.
    ///
    /// The one test behind [`crate::PeriodicCell::new`], [`Self::interplanar_spacings`],
    /// the stress / virial conversion and `ferro-structure`'s vacuum, merge and
    /// supercell (review A-13: five tests before, and the one in `interplanar_spacings`
    /// let a NaN cell through as infinite spacings).
    pub fn check(&self, what: &str) -> Result<f64> {
        if !self.matrix.iter().all(|v| v.is_finite()) {
            return Err(ChemError::ValidationError(format!("{what} is not finite")));
        }
        let v = self.volume();
        if v > 1e-12 {
            Ok(v)
        } else {
            Err(ChemError::ValidationError(format!("{what} is singular (volume {v:e} Å³)")))
        }
    }

    /// The volume, or an error for a zero-volume cell — the guard of the stress /
    /// virial conversion, which divides or multiplies by it.
    pub(crate) fn nonzero_volume(&self) -> Result<f64> {
        self.check("the cell (needed to convert between stress and virial)")
    }

    /// The stress (eV/Å³) a virial (eV) stands for in this cell: `virial / V`, no
    /// sign change (see [`crate::Frame::stress`]).
    ///
    /// With [`crate::Frame::virial`] the one place the two are converted (review
    /// B13); a zero-volume cell is an error rather than an infinity.
    pub fn stress_from_virial(&self, virial: &Matrix3<f64>) -> Result<Matrix3<f64>> {
        Ok(virial / self.nonzero_volume()?)
    }

    /// Interplanar spacings d₁, d₂, d₃ \[Å\] — the perpendicular distances between
    /// the pairs of cell faces normal to each reciprocal-lattice direction.
    ///
    /// `dᵢ = 1 / ‖(Mᵀ)⁻¹.row(i)‖`, where the rows of `(Mᵀ)⁻¹` are the reciprocal
    /// lattice vectors.
    ///
    /// This — not `lengths()` — is the correct geometric bound for the minimum-image
    /// convention: MIC is valid only for `r < min(d)/2`. For non-orthogonal cells
    /// `d < L`, so using `lengths()` overestimates the usable cutoff.
    ///
    /// Returns `Err` when the cell is refused by [`Self::check`].
    pub fn interplanar_spacings(&self) -> Result<[f64; 3]> {
        self.check("cell matrix")?;
        let inv = self.matrix.transpose().try_inverse()
            .ok_or_else(|| ChemError::ValidationError("cell matrix is singular".into()))?;
        let d = |i: usize| {
            let n = inv.row(i).norm();
            if n > 0.0 { 1.0 / n } else { f64::INFINITY }
        };
        Ok([d(0), d(1), d(2)])
    }

    /// Upper bound on a cutoff radius compatible with the minimum-image convention:
    /// half the smallest interplanar spacing. Same formula as `PeriodicCell::mic_bound`
    /// (slim scan A-11).
    pub fn minimum_image_cutoff(&self) -> Result<f64> {
        let d = self.interplanar_spacings()?;
        Ok(d[0].min(d[1]).min(d[2]) * 0.5)
    }

    // ── 坐标转换 ──────────────────────────────────────────────────────────────

    /// 分数坐标 → Cartesian 坐标（Å）。
    pub fn fractional_to_cartesian(&self, frac: Vector3<f64>) -> Vector3<f64> {
        self.matrix.transpose() * frac
    }

    /// Cartesian 坐标（Å）→ 分数坐标。晶胞矩阵奇异时返回 `Err`。
    pub fn cartesian_to_fractional(&self, cart: Vector3<f64>) -> Result<Vector3<f64>> {
        self.matrix
            .transpose()
            .try_inverse()
            .map(|inv| inv * cart)
            .ok_or_else(|| ChemError::ValidationError("cell matrix is singular".into()))
    }

    // ── 周期性处理 ────────────────────────────────────────────────────────────

    /// Folds a Cartesian position into the cell, fractional coordinates in [0, 1) along
    /// all three axes.
    ///
    /// `rem_euclid(1.0)` keeps the output strictly in [0, 1), so the ~1e-15 error of the
    /// matrix inverse cannot make `floor` jump at a boundary.
    /// Returns `Err` when the cell matrix is singular.
    ///
    /// One of three folds, kept apart (review A-16): this one ignores `pbc` and
    /// `origin` and round-trips every position through fractions; [`crate::Frame::wrap_all`]
    /// follows `pbc`, moves only atoms outside by whole lattice vectors and folds about
    /// 0 (the cube grid it serves starts there); `ferro-structure`'s `merge_frames`
    /// folds about `origin` and refuses an atom outside along a non-periodic axis.
    /// Three contracts, three callers that need each.
    pub fn wrap_position(&self, cart: Vector3<f64>) -> Result<Vector3<f64>> {
        let frac = self.cartesian_to_fractional(cart)?;
        let wrapped = Vector3::new(
            frac.x.rem_euclid(1.0),
            frac.y.rem_euclid(1.0),
            frac.z.rem_euclid(1.0),
        );
        Ok(self.fractional_to_cartesian(wrapped))
    }

    /// Minimum-image form of a displacement: each fractional component in [-0.5, 0.5],
    /// both ends included (`c - c.round()`, and `round` takes ±0.5 away from zero), as
    /// [`crate::PeriodicCell::minimum_image`].
    ///
    /// Example: the shortest periodic distance between two atoms
    /// ```ignore
    /// let diff = atom_b.position - atom_a.position;
    /// let mic_diff = cell.minimum_image(diff)?;
    /// let dist = mic_diff.norm();
    /// ```
    /// Returns `Err` when the cell matrix is singular.
    ///
    /// Hot loops use `PeriodicCell::minimum_image` (inverse cached); the callers left here
    /// are in the suspended Bader / cube_sdf groups (slim scan A-11).
    pub fn minimum_image(&self, diff: Vector3<f64>) -> Result<Vector3<f64>> {
        let frac = self.cartesian_to_fractional(diff)?;
        let mic_frac = Vector3::new(
            frac.x - frac.x.round(),
            frac.y - frac.y.round(),
            frac.z - frac.z.round(),
        );
        Ok(self.fractional_to_cartesian(mic_frac))
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
    fn test_interplanar_spacings_cubic_equal_lengths() {
        // 正交晶胞：面间距 = 边长
        let cell = cubic(10.0);
        let d = cell.interplanar_spacings().unwrap();
        for di in d {
            assert!((di - 10.0).abs() < 1e-10, "d = {di}, expected 10");
        }
        assert!((cell.minimum_image_cutoff().unwrap() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_interplanar_spacings_smaller_than_lengths_when_skewed() {
        // 三斜晶胞：面间距严格小于边长 —— 这正是 MIC 上界不能用 lengths() 的原因
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 60.0, 70.0, 80.0).unwrap();
        let l = cell.lengths();
        let d = cell.interplanar_spacings().unwrap();
        for i in 0..3 {
            assert!(
                d[i] < l[i] - 1e-6,
                "axis {i}: d = {:.6} should be < L = {:.6}", d[i], l[i]
            );
        }
        // MIC 上界随之小于 min(L)/2
        let l_min_half = l.iter().cloned().fold(f64::INFINITY, f64::min) * 0.5;
        assert!(cell.minimum_image_cutoff().unwrap() < l_min_half);
    }

    #[test]
    fn test_interplanar_spacings_matches_volume_over_area() {
        // d_i = V / |a_j × a_k|，与倒格矢定义一致
        let cell = Cell::from_lengths_angles(7.0, 9.0, 11.0, 75.0, 85.0, 95.0).unwrap();
        let row = |i: usize| Vector3::new(cell.matrix[(i, 0)], cell.matrix[(i, 1)], cell.matrix[(i, 2)]);
        let v = cell.volume();
        let d = cell.interplanar_spacings().unwrap();
        let expect = [
            v / row(1).cross(&row(2)).norm(),
            v / row(2).cross(&row(0)).norm(),
            v / row(0).cross(&row(1)).norm(),
        ];
        for i in 0..3 {
            assert!(
                (d[i] - expect[i]).abs() < 1e-9,
                "axis {i}: got {:.9}, expected {:.9}", d[i], expect[i]
            );
        }
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
        // cartesian_to_fractional 内部做矩阵求逆，残差约 1e-15，容差取 1e-10
        let cell = Cell::from_lengths_angles(4.0, 5.0, 6.0, 80.0, 90.0, 100.0).unwrap();
        let frac = Vector3::new(0.3, 0.5, 0.7);
        let cart = cell.fractional_to_cartesian(frac);
        let back = cell.cartesian_to_fractional(cart).unwrap();
        assert!((back - frac).norm() < 1e-10);
    }

    #[test]
    fn test_wrap_position() {
        let cell = cubic(10.0);
        let pos = Vector3::new(12.0, -1.0, 5.0);
        let wrapped = cell.wrap_position(pos).unwrap();
        assert!(wrapped.x >= 0.0 && wrapped.x < 10.0);
        assert!(wrapped.y >= 0.0 && wrapped.y < 10.0);
    }

    #[test]
    fn test_wrap_position_triclinic_boundary() {
        // 非正交晶胞：矩阵求逆误差可能让分数坐标略为负（如 -1e-15），
        // rem_euclid 应把这类值折叠到 [0,1) 而不是触发 floor 边界跳变
        let cell = Cell::from_lengths_angles(5.0, 5.0, 5.0, 80.0, 85.0, 95.0).unwrap();
        // 原点处的 Cartesian 坐标经过来回转换后应仍在盒内
        let origin = cell.fractional_to_cartesian(Vector3::new(0.0, 0.0, 0.0));
        let wrapped = cell.wrap_position(origin).unwrap();
        let frac = cell.cartesian_to_fractional(wrapped).unwrap();
        assert!(frac.x >= -1e-10 && frac.x < 1.0 + 1e-10);
        assert!(frac.y >= -1e-10 && frac.y < 1.0 + 1e-10);
        assert!(frac.z >= -1e-10 && frac.z < 1.0 + 1e-10);
        // 盒外坐标必须被折叠回合理范围
        let out = cell.fractional_to_cartesian(Vector3::new(1.3, -0.5, 2.1));
        let wrapped2 = cell.wrap_position(out).unwrap();
        let frac2 = cell.cartesian_to_fractional(wrapped2).unwrap();
        assert!(frac2.x >= -1e-10 && frac2.x < 1.0 + 1e-10);
        assert!(frac2.y >= -1e-10 && frac2.y < 1.0 + 1e-10);
        assert!(frac2.z >= -1e-10 && frac2.z < 1.0 + 1e-10);
    }

    #[test]
    fn nan_and_out_of_range_lattice_parameters_are_rejected() {
        let ok = |a: f64, al: f64| Cell::from_lengths_angles(a, 5.0, 5.0, al, 90.0, 90.0).is_ok();
        assert!(ok(5.0, 90.0));
        assert!(!ok(f64::NAN, 90.0));
        assert!(!ok(f64::INFINITY, 90.0));
        assert!(!ok(5.0, f64::NAN));
        assert!(!ok(5.0, 0.0));
        assert!(!ok(5.0, 180.0));
    }

    #[test]
    fn test_minimum_image() {
        let cell = cubic(10.0);
        // 距离 9 Å，MIC 应返回 -1 Å
        let diff = Vector3::new(9.0, 0.0, 0.0);
        let mic = cell.minimum_image(diff).unwrap();
        assert!((mic.x + 1.0).abs() < 1e-10);
    }

    // Review A-1 / A-13: a NaN cell used to pass `interplanar_spacings` as [inf; 3]
    #[test]
    fn check_refuses_non_finite_and_singular_cells() {
        let cubic = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        assert!((cubic.check("cell").unwrap() - 1000.0).abs() < 1e-9);
        for bad in [f64::NAN, f64::INFINITY] {
            let mut c = cubic.clone();
            c.matrix[(0, 1)] = bad;
            let e = c.check("cell").unwrap_err().to_string();
            assert!(e.contains("cell is not finite"), "{e}");
            assert!(c.interplanar_spacings().is_err(), "{bad}");
            assert!(c.minimum_image_cutoff().is_err(), "{bad}");
        }
        let mut flat = cubic.clone();
        flat.matrix[(2, 2)] = 0.0;
        assert!(flat.check("cell").unwrap_err().to_string().contains("cell is singular"));
        assert!(flat.interplanar_spacings().is_err());
    }
}
