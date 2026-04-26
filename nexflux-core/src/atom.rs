//! Atom data structure.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// 单个原子实例。
///
/// 原子在帧中的唯一标识是其在 `Frame::atoms` 中的下标，不单独存储 index 字段。
/// `mass` 为 `None` 时，调用 [`Atom::effective_mass`] 将自动从元素表中查找标准值。
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    /// 元素符号，如 "Fe"、"O"
    pub element: String,
    /// Cartesian 坐标，单位 Å
    pub position: Vector3<f64>,
    /// 可选的人类可读标签，如 "Fe1"、"Fe2"，用于区分同种元素的不同位点
    pub label: Option<String>,
    /// 质量覆盖值，单位 amu；`None` 表示使用元素表标准值
    pub mass: Option<f64>,
    /// 初始磁矩（用于 DFT 自旋极化计算的初猜），单位 μ_B
    pub magmom: Option<f64>,
    /// 原子电荷，单位 e；通常由 Bader/DDEC 后处理写回
    pub charge: Option<f64>,
}

impl Atom {
    /// 最简构造器，仅需元素符号和位置。
    pub fn new(element: impl Into<String>, position: Vector3<f64>) -> Self {
        Self {
            element: element.into(),
            position,
            label: None,
            mass: None,
            magmom: None,
            charge: None,
        }
    }

    /// 返回有效质量（amu）：优先使用 `mass` 字段，否则查元素表，查不到返回 1.0。
    pub fn effective_mass(&self) -> f64 {
        self.mass.unwrap_or_else(|| {
            crate::data::elements::by_symbol(&self.element)
                .map(|e| e.atomic_mass)
                .unwrap_or(1.0)
        })
    }

    /// 计算与另一原子的 Cartesian 距离（Å），不考虑周期性。
    /// 周期性距离请通过 [`crate::cell::Cell::minimum_image`] 处理。
    pub fn distance_to(&self, other: &Atom) -> f64 {
        (self.position - other.position).norm()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let atom = Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0));
        assert_eq!(atom.element, "Fe");
        assert!(atom.label.is_none());
        assert!(atom.mass.is_none());
    }

    #[test]
    fn test_effective_mass_from_table() {
        let atom = Atom::new("Fe", Vector3::zeros());
        assert!((atom.effective_mass() - 55.845).abs() < 1e-3);
    }

    #[test]
    fn test_effective_mass_override() {
        let mut atom = Atom::new("Fe", Vector3::zeros());
        atom.mass = Some(56.0);
        assert!((atom.effective_mass() - 56.0).abs() < 1e-10);
    }

    #[test]
    fn test_distance() {
        let a = Atom::new("O", Vector3::new(0.0, 0.0, 0.0));
        let b = Atom::new("H", Vector3::new(1.0, 0.0, 0.0));
        assert!((a.distance_to(&b) - 1.0).abs() < 1e-10);
    }
}
