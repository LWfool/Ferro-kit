//! 轨迹分析工具（旧接口，保持编译兼容）
//!
//! MSD、RMSD、Rg 等完整分析已迁移至 `md/` 子模块。
//! 此文件仅保留基础实现供直接导入使用。

use ferro_core::Trajectory;
use nalgebra::Vector3;
use rayon::prelude::*;

/// 计算均方位移 (MSD)
pub fn mean_squared_displacement(traj: &Trajectory) -> Vec<f64> {
    if traj.frames.is_empty() {
        return vec![];
    }

    let n_frames = traj.n_frames();
    let n_atoms = traj.n_atoms().unwrap_or(0);
    let mut msd = vec![0.0; n_frames];

    for dt in 0..n_frames {
        let mut sum = 0.0;
        let mut count = 0usize;
        for t0 in 0..(n_frames - dt) {
            let frame0 = &traj.frames[t0];
            let frame1 = &traj.frames[t0 + dt];
            for i in 0..n_atoms {
                if let (Some(a0), Some(a1)) = (frame0.atoms.get(i), frame1.atoms.get(i)) {
                    sum += (a1.position - a0.position).norm_squared();
                    count += 1;
                }
            }
        }
        if count > 0 {
            msd[dt] = sum / count as f64;
        }
    }
    msd
}

/// 计算相对于第一帧的 RMSD
pub fn rmsd_to_first(traj: &Trajectory) -> Vec<f64> {
    if traj.frames.is_empty() {
        return vec![];
    }
    let reference = &traj.frames[0];
    let n_atoms = reference.n_atoms();

    traj.frames
        .par_iter()
        .map(|frame| {
            let sum: f64 = (0..n_atoms)
                .filter_map(|i| {
                    let r = reference.atoms.get(i)?;
                    let f = frame.atoms.get(i)?;
                    Some((f.position - r.position).norm_squared())
                })
                .sum();
            if n_atoms > 0 { (sum / n_atoms as f64).sqrt() } else { 0.0 }
        })
        .collect()
}

/// 质心随时间的轨迹
pub fn center_of_mass_trajectory(traj: &Trajectory) -> Vec<Vector3<f64>> {
    traj.frames.iter().map(|f| f.center_of_mass()).collect()
}

/// 回转半径随时间变化
pub fn radius_of_gyration_trajectory(traj: &Trajectory) -> Vec<f64> {
    traj.frames
        .par_iter()
        .map(|frame| {
            let com = frame.center_of_mass();
            let total_mass = frame.total_mass();
            if total_mass == 0.0 { return 0.0; }
            let sum: f64 = frame
                .atoms
                .iter()
                .map(|a| {
                    let r = (a.position - com).norm();
                    a.effective_mass() * r * r
                })
                .sum();
            (sum / total_mass).sqrt()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Frame, Trajectory};
    use nalgebra::Vector3;

    #[test]
    fn test_rmsd_same_frames() {
        let mut traj = Trajectory::new();
        let mut frame = Frame::new();
        frame.add_atom(Atom::new("C", Vector3::new(0.0, 0.0, 0.0)));
        traj.add_frame(frame.clone());
        traj.add_frame(frame);

        let rmsd = rmsd_to_first(&traj);
        assert_eq!(rmsd.len(), 2);
        assert!(rmsd[0].abs() < 1e-10);
    }

    #[test]
    fn test_msd_empty() {
        let traj = Trajectory::new();
        assert!(mean_squared_displacement(&traj).is_empty());
    }
}
