//! Trajectory data structure.

use serde::{Deserialize, Serialize};

use crate::frame::Frame;

/// 轨迹元数据。
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct TrajectoryMetadata {
    /// 帧间时间步长（fs）；`None` 表示静态结构或未知
    pub timestep: Option<f64>,
    /// 来源文件或软件名称，如 "VASP OUTCAR"、"LAMMPS dump"
    pub source: Option<String>,
}

/// 轨迹：一个或多个帧的时间序列。
///
/// 单帧结构文件也以 `Trajectory { frames: vec![frame] }` 形式存储，
/// 保证所有模块的 API 签名统一，无需区分单帧/多帧。
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Trajectory {
    pub frames: Vec<Frame>,
    pub metadata: TrajectoryMetadata,
}

impl Trajectory {
    /// 创建空轨迹。
    pub fn new() -> Self {
        Self {
            frames: Vec::new(),
            metadata: TrajectoryMetadata::default(),
        }
    }

    /// 由单帧构造轨迹（结构文件读取的常用路径）。
    pub fn from_frame(frame: Frame) -> Self {
        Self {
            frames: vec![frame],
            metadata: TrajectoryMetadata::default(),
        }
    }

    // ── 帧访问 ───────────────────────────────────────────────────────────────

    pub fn n_frames(&self) -> usize {
        self.frames.len()
    }

    /// 原子数（取自第一帧；假设各帧原子数一致）。
    pub fn n_atoms(&self) -> Option<usize> {
        self.frames.first().map(|f| f.n_atoms())
    }

    pub fn frame(&self, index: usize) -> Option<&Frame> {
        self.frames.get(index)
    }

    pub fn frame_mut(&mut self, index: usize) -> Option<&mut Frame> {
        self.frames.get_mut(index)
    }

    pub fn first(&self) -> Option<&Frame> {
        self.frames.first()
    }

    pub fn last(&self) -> Option<&Frame> {
        self.frames.last()
    }

    pub fn add_frame(&mut self, frame: Frame) {
        self.frames.push(frame);
    }

    /// Return a new trajectory containing only the last `n` frames.
    ///
    /// If `n` ≥ the trajectory length, all frames are returned (no panic).
    /// Metadata is preserved unchanged.
    pub fn tail(&self, n: usize) -> Trajectory {
        let start = self.n_frames().saturating_sub(n);
        Trajectory {
            frames: self.frames[start..].to_vec(),
            metadata: self.metadata.clone(),
        }
    }

    pub fn iter_frames(&self) -> impl Iterator<Item = &Frame> {
        self.frames.iter()
    }

    // ── 时间 ─────────────────────────────────────────────────────────────────

    /// 返回第 `index` 帧对应的时间（fs）；需要 `metadata.timestep` 不为 `None`。
    pub fn time_at(&self, index: usize) -> Option<f64> {
        self.metadata.timestep.map(|dt| dt * index as f64)
    }
}

impl Default for Trajectory {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use nalgebra::Vector3;

    fn make_frame(x: f64) -> Frame {
        let mut f = Frame::new();
        f.add_atom(Atom::new("Fe", Vector3::new(x, 0.0, 0.0)));
        f
    }

    #[test]
    fn test_from_frame() {
        let traj = Trajectory::from_frame(make_frame(0.0));
        assert_eq!(traj.n_frames(), 1);
        assert_eq!(traj.n_atoms(), Some(1));
    }

    #[test]
    fn test_add_frames() {
        let mut traj = Trajectory::new();
        traj.add_frame(make_frame(0.0));
        traj.add_frame(make_frame(1.0));
        assert_eq!(traj.n_frames(), 2);
    }

    #[test]
    fn test_time_at() {
        let mut traj = Trajectory::new();
        traj.metadata.timestep = Some(2.0);
        traj.add_frame(make_frame(0.0));
        traj.add_frame(make_frame(1.0));
        assert_eq!(traj.time_at(0), Some(0.0));
        assert_eq!(traj.time_at(1), Some(2.0));
    }

    #[test]
    fn test_first_last() {
        let mut traj = Trajectory::new();
        traj.add_frame(make_frame(0.0));
        traj.add_frame(make_frame(5.0));
        assert!((traj.first().unwrap().atom(0).position.x).abs() < 1e-10);
        assert!((traj.last().unwrap().atom(0).position.x - 5.0).abs() < 1e-10);
    }
}
