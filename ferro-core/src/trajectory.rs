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
        self.select(start, None, 1)
    }

    /// Return the frames in `[start, end]` — both ends INCLUSIVE, 0-based — taking
    /// every `stride`-th one.
    ///
    /// The inclusive upper bound matches the frame numbering `ferro info` prints,
    /// so `--end 4` selects the frame shown as `Frame 4`. `end` of `None` means the
    /// last frame. Out-of-range bounds clamp rather than panic; a `start` past the
    /// end yields an empty trajectory. A `stride` of 0 is treated as 1.
    ///
    /// Metadata is preserved unchanged — note that `timestep` then no longer
    /// describes the spacing of the returned frames when `stride > 1`.
    pub fn select(&self, start: usize, end: Option<usize>, stride: usize) -> Trajectory {
        let last = self.n_frames().saturating_sub(1);
        let end = end.unwrap_or(last).min(last);
        let stride = stride.max(1);

        let frames = if start > end || self.frames.is_empty() {
            Vec::new()
        } else {
            self.frames[start..=end].iter().step_by(stride).cloned().collect()
        };
        Trajectory { frames, metadata: self.metadata.clone() }
    }

    /// The frame indices [`Trajectory::select`] would pick.
    ///
    /// Callers that name one output file per frame need the ORIGINAL index, not the
    /// position within the selection: frames taken with `stride = 10` should be
    /// named 0, 10, 20 so the products can be traced back to the trajectory.
    pub fn select_indices(&self, start: usize, end: Option<usize>, stride: usize) -> Vec<usize> {
        let last = self.n_frames().saturating_sub(1);
        let end = end.unwrap_or(last).min(last);
        if start > end || self.frames.is_empty() {
            return Vec::new();
        }
        (start..=end).step_by(stride.max(1)).collect()
    }

    /// Indices of `count` frames spread evenly over `[start, end]`, both ends
    /// INCLUSIVE and both endpoints always taken when `count >= 2`.
    ///
    /// This is the `--number` form of selection: the caller asks for a total, not a
    /// spacing. Endpoints are included because the last frame of a run is usually
    /// the most equilibrated one, and a `floor(n/count)` walk would systematically
    /// drop it. Returns at most `end - start + 1` indices — asking for more frames
    /// than exist yields every frame once, never a duplicate.
    pub fn spread_indices(&self, start: usize, end: Option<usize>, count: usize) -> Vec<usize> {
        let last = self.n_frames().saturating_sub(1);
        let end = end.unwrap_or(last).min(last);
        if count == 0 || start > end || self.frames.is_empty() {
            return Vec::new();
        }

        let available = end - start + 1;
        if count >= available {
            return (start..=end).collect();
        }
        if count == 1 {
            return vec![start];
        }

        // linspace 含两端：i=0 给 start，i=count-1 给 end
        let span = (end - start) as f64;
        (0..count)
            .map(|i| start + (span * i as f64 / (count - 1) as f64).round() as usize)
            .collect()
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

    /// 用 x 坐标标记帧序号，便于断言选出来的到底是哪几帧
    fn numbered_traj(n: usize) -> Trajectory {
        let mut traj = Trajectory::new();
        for i in 0..n {
            traj.add_frame(make_frame(i as f64));
        }
        traj
    }

    fn picked(traj: &Trajectory) -> Vec<usize> {
        traj.frames.iter().map(|f| f.atoms[0].position.x as usize).collect()
    }

    #[test]
    fn test_select_bounds_are_inclusive_on_both_ends() {
        // 闭区间：`--end 4` 取到 `ferro info` 显示的那个 Frame 4
        let traj = numbered_traj(10);
        assert_eq!(picked(&traj.select(0, Some(4), 1)), vec![0, 1, 2, 3, 4]);
        assert_eq!(picked(&traj.select(2, Some(2), 1)), vec![2]);
    }

    #[test]
    fn test_select_stride_counts_from_start_not_from_zero() {
        let traj = numbered_traj(10);
        assert_eq!(picked(&traj.select(1, Some(9), 3)), vec![1, 4, 7]);
        assert_eq!(picked(&traj.select(0, None, 4)), vec![0, 4, 8]);
    }

    #[test]
    fn test_select_indices_agree_with_select() {
        // 两者必须选出同一批帧，否则「文件名里的索引」会指向别的帧
        let traj = numbered_traj(10);
        for (start, end, stride) in
            [(0, None, 1), (1, Some(9), 3), (0, Some(4), 2), (7, None, 5), (4, Some(2), 1)]
        {
            assert_eq!(
                traj.select_indices(start, end, stride),
                picked(&traj.select(start, end, stride)),
                "start={start} end={end:?} stride={stride}"
            );
        }
    }

    #[test]
    fn test_select_clamps_instead_of_panicking() {
        let traj = numbered_traj(5);
        // end 越界 → 夹到最后一帧
        assert_eq!(picked(&traj.select(3, Some(999), 1)), vec![3, 4]);
        // start 越界 → 空轨迹，不是 panic
        assert!(traj.select(99, None, 1).frames.is_empty());
        // start > end → 空
        assert!(traj.select(4, Some(2), 1).frames.is_empty());
        // stride 0 当作 1，不是除零
        assert_eq!(picked(&traj.select(0, Some(2), 0)), vec![0, 1, 2]);
        // 空轨迹不 panic
        assert!(Trajectory::new().select(0, None, 1).frames.is_empty());
    }

    #[test]
    fn test_tail_still_takes_the_last_n() {
        let traj = numbered_traj(10);
        assert_eq!(picked(&traj.tail(3)), vec![7, 8, 9]);
        // n 超长度取全部，不 panic
        assert_eq!(picked(&traj.tail(99)).len(), 10);
    }

    #[test]
    fn test_spread_includes_both_endpoints() {
        let traj = numbered_traj(100);
        // 100 帧取 3 个：含首含尾
        assert_eq!(traj.spread_indices(0, None, 3), vec![0, 50, 99]);
        assert_eq!(traj.spread_indices(0, None, 2), vec![0, 99]);
        // 在子区间里同样含两端
        assert_eq!(traj.spread_indices(10, Some(20), 3), vec![10, 15, 20]);
    }

    #[test]
    fn test_spread_never_repeats_a_frame() {
        let traj = numbered_traj(5);
        // 要的比有的多 → 每帧一次，不是补重复
        assert_eq!(traj.spread_indices(0, None, 99), vec![0, 1, 2, 3, 4]);
        assert_eq!(traj.spread_indices(0, None, 5), vec![0, 1, 2, 3, 4]);
        // 要 1 个给起点；要 0 个给空
        assert_eq!(traj.spread_indices(2, None, 1), vec![2]);
        assert!(traj.spread_indices(0, None, 0).is_empty());
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
