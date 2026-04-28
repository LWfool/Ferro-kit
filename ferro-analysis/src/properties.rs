//! 分子性质计算（占位实现）

use ferro_core::Frame;

/// 偶极矩（需要原子电荷；目前返回 0）
pub fn dipole_moment(_frame: &Frame) -> f64 {
    // TODO: Σ q_i * r_i
    0.0
}

/// 体系总能量（从 frame.energy 读取）
pub fn total_energy(frame: &Frame) -> Option<f64> {
    frame.energy
}
