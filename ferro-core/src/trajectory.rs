//! Trajectory data structure.


use crate::data::elements::canonical_order;
use crate::error::{ChemError, Result};
use crate::frame::Frame;

/// Indices `select` would pick out of a sequence of `n` items.
///
/// The index-level form of [`Trajectory::select_indices`], for callers whose
/// sequence is not the trajectory itself — dataset filtering applies the range
/// to the frames that SURVIVED earlier criteria, so the bound is the survivor
/// count, not `n_frames()`. Keeping one implementation is what stops the two
/// from drifting apart on the inclusive-end and clamping rules.
pub fn select_range(n: usize, start: usize, end: Option<usize>, stride: usize) -> Vec<usize> {
    let last = n.saturating_sub(1);
    let end = end.unwrap_or(last).min(last);
    if n == 0 || start > end {
        return Vec::new();
    }
    (start..=end).step_by(stride.max(1)).collect()
}

/// Indices of `count` items spread evenly over `[start, end]`, both ends taken.
///
/// The index-level form of [`Trajectory::spread_indices`]; see [`select_range`]
/// for why both forms exist.
pub fn spread_range(n: usize, start: usize, end: Option<usize>, count: usize) -> Vec<usize> {
    let last = n.saturating_sub(1);
    let end = end.unwrap_or(last).min(last);
    if n == 0 || count == 0 || start > end {
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

/// Trajectory metadata.
///
/// Two fields, two meanings (review B4): a structure file's own title, and the
/// program that produced an output.  Writers put back only the title — a program
/// name is not a title for the next file.
#[derive(Debug, Clone, Default)]
pub struct TrajectoryMetadata {
    /// The title the file gives itself: a POSCAR / xyz / LAMMPS data comment line,
    /// a CIF block name, a PDB TITLE, an extxyz `config_type` or `comment`
    pub title: Option<String>,
    /// The program, and its version when the output prints one, that wrote an
    /// AIMD output, e.g. "CP2K 2024.1 out"
    pub program: Option<String>,
}

/// 轨迹：一个或多个帧的时间序列。
///
/// 单帧结构文件也以 `Trajectory { frames: vec![frame] }` 形式存储，
/// 保证所有模块的 API 签名统一，无需区分单帧/多帧。
#[derive(Debug, Clone)]
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

    /// Appends a frame, sharing the previous frame's site table when the two are equal
    /// — a trajectory of one system then stores its elements and labels once.
    pub fn add_frame(&mut self, mut frame: Frame) {
        if let Some(last) = self.frames.last() {
            frame.share_sites_with(last);
        }
        self.frames.push(frame);
    }

    /// Links every frame to its predecessor's site table where the two are equal, for
    /// a trajectory whose frames were pushed onto `frames` directly.
    pub fn share_sites(&mut self) {
        for k in 1..self.frames.len() {
            let (done, rest) = self.frames.split_at_mut(k);
            rest[0].share_sites_with(&done[k - 1]);
        }
    }

    /// Check that every frame has the same atoms as frame 0: same count, same element
    /// at every index. The error states which frame differs; the caller adds why that
    /// matters to it.
    ///
    /// Time-correlation analyses (MSD, VACF, rotcorr, bond lifetimes, van Hove, jump
    /// maps) follow atom `i` across frames by its index, so a frame with a different
    /// atom list makes them index out of bounds or silently pair different atoms; a
    /// DeePMD system writes one `type.raw` for all its frames. g(r) and S(q) need a
    /// fixed count of each type but not a fixed order, and check that themselves.
    /// Angles, networks and density maps need neither: a varying atom count is
    /// legitimate for them (GCMC, reactive MD).
    pub fn check_same_atoms(&self) -> Result<()> {
        let Some(first) = self.frames.first() else { return Ok(()) };
        for (k, frame) in self.frames.iter().enumerate().skip(1) {
            if frame.n_atoms() != first.n_atoms() {
                return Err(ChemError::ValidationError(format!(
                    "frame {k} has {} atoms but frame 0 has {}",
                    frame.n_atoms(), first.n_atoms()
                )));
            }
            if frame.shares_sites_with(first) {
                continue;
            }
            // Elements only, not labels: in a trajectory `net` exported, labels change
            // from frame to frame while the atoms stay the same
            if let Some(j) = (0..first.n_atoms())
                .find(|&j| frame.element(j) != first.element(j))
            {
                return Err(ChemError::ValidationError(format!(
                    "atom {j} is {} in frame {k} but {} in frame 0",
                    frame.element(j), first.element(j)
                )));
            }
        }
        Ok(())
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
    /// Metadata is preserved unchanged.
    pub fn select(&self, start: usize, end: Option<usize>, stride: usize) -> Trajectory {
        self.subset(&self.select_indices(start, end, stride))
    }

    /// The frame indices [`Trajectory::select`] would pick.
    ///
    /// Callers that name one output file per frame need the ORIGINAL index, not the
    /// position within the selection: frames taken with `stride = 10` should be
    /// named 0, 10, 20 so the products can be traced back to the trajectory.
    pub fn select_indices(&self, start: usize, end: Option<usize>, stride: usize) -> Vec<usize> {
        select_range(self.n_frames(), start, end, stride)
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
        spread_range(self.n_frames(), start, end, count)
    }

    /// The frames at `indices`, in the order given.
    ///
    /// [`Trajectory::select`] covers ranges; this covers an arbitrary verdict —
    /// dataset filtering yields a list of surviving frame numbers that is not a
    /// range at all. Indices past the end are skipped rather than panicking, so
    /// a stale index list cannot take down a long run.
    pub fn subset(&self, indices: &[usize]) -> Trajectory {
        let frames = indices
            .iter()
            .filter_map(|&i| self.frames.get(i).cloned())
            .collect();
        Trajectory { frames, metadata: self.metadata.clone() }
    }

    pub fn iter_frames(&self) -> impl Iterator<Item = &Frame> {
        self.frames.iter()
    }

    /// The trajectory with every frame's atoms — and the per-atom arrays — in
    /// canonical (Z, symbol) order ([`crate::data::elements::canonical_cmp`]).
    ///
    /// Each frame takes its own permutation: xyz and other formats do not promise a
    /// constant atom order, and files written independently disagree. Consecutive
    /// frames with the same element sequence reuse the previous permutation.
    pub fn canonically_sorted(&self) -> Trajectory {
        let mut out = self.clone();
        // (element sequence, permutation), reused while consecutive frames keep the sequence
        let mut cached: Option<(Vec<String>, Vec<usize>)> = None;
        for f in &mut out.frames {
            let elements: Vec<String> = f.symbols().into_iter().map(str::to_string).collect();
            if cached.as_ref().is_none_or(|(e, _)| *e != elements) {
                let order = canonical_order(&elements);
                cached = Some((elements, order));
            }
            let (_, order) = cached.as_ref().expect("cached above");
            if order.iter().enumerate().all(|(i, &j)| i == j) {
                continue;
            }
            // Frame::permuted moves forces, velocities and charges with the atoms;
            // energy, cell and stress do not depend on atom numbering
            *f = f.permuted(order);
        }
        // Each permuted frame got a site table of its own; equal ones are shared again
        out.share_sites();
        out
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

    /// One frame of `elements` at x = 0, 1, 2…, with forces 0.0, 0.1, 0.2… along x.
    fn sort_input(elements: &[&str]) -> Trajectory {
        let mut f = Frame::new();
        f.set_atoms(elements.iter().enumerate().map(|(i, e)| Atom::new(*e, Vector3::new(i as f64, 0.0, 0.0))));
        f.set_forces(Some((0..elements.len()).map(|i| Vector3::new(i as f64 * 0.1, 0.0, 0.0)).collect()))
            .unwrap();
        Trajectory::from_frame(f)
    }

    #[test]
    fn sorting_carries_the_per_atom_arrays_along() {
        let s = sort_input(&["Zn", "O"]).canonically_sorted();
        let f = &s.frames[0];
        assert_eq!(f.element(0), "O");
        // O was atom 1: its position and its force must come along
        assert_eq!(f.positions()[0].x, 1.0);
        assert!((f.forces().unwrap()[0].x - 0.1).abs() < 1e-12);
        assert_eq!(f.element(1), "Zn");
        assert_eq!(f.positions()[1].x, 0.0);
    }

    /// Velocities are per-atom too, and were once left behind.
    ///
    /// A DeePMD system carries none, so `merge` never showed the bug; the
    /// diagnostics `collect` exports (a dump with velocities) do.
    #[test]
    fn sorting_carries_velocities_too() {
        let mut t = sort_input(&["Zn", "O", "Al"]);
        t.frames[0].set_velocities(Some(
            (0..3).map(|i| Vector3::new(i as f64 * 10.0, 0.0, 0.0)).collect(),
        )).unwrap();
        let s = t.canonically_sorted();
        let f = &s.frames[0];
        let v = f.velocities().unwrap();
        assert_eq!(f.element(0), "O");
        assert!((v[0].x - 10.0).abs() < 1e-12, "the velocity of O must move with it");
        assert!((v[2].x - 0.0).abs() < 1e-12, "and so must that of Zn");
    }

    /// Each frame takes its own permutation: one file may list the atoms in a
    /// different order from frame to frame.
    #[test]
    fn a_frame_with_its_own_atom_order_gets_its_own_permutation() {
        let mut t = sort_input(&["Zn", "O"]);
        let mut f1 = sort_input(&["O", "Zn"]).frames.remove(0);
        // frame 1: the same configuration with the atoms reversed, forces following them
        f1.positions_mut()[0].x = 1.0;
        f1.positions_mut()[1].x = 0.0;
        f1.set_forces(Some(vec![Vector3::new(0.1, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0)])).unwrap();
        t.frames.push(f1);

        let s = t.canonically_sorted();
        for (k, f) in s.frames.iter().enumerate() {
            assert_eq!(f.symbols(), vec!["O", "Zn"], "frame {k} is not in canonical order");
            assert_eq!(f.positions()[0].x, 1.0, "frame {k}: the position of O did not come along");
            assert!((f.forces().unwrap()[0].x - 0.1).abs() < 1e-12, "frame {k}: the force of O did not come along");
        }
    }

    #[test]
    fn an_already_sorted_system_is_untouched() {
        let t = sort_input(&["O", "O", "Al"]);
        let s = t.canonically_sorted();
        assert_eq!(t.frames[0].positions(), s.frames[0].positions());
    }

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

    #[test]
    fn test_check_same_atoms() {
        // 原子数相同、只有 label 变：net 标注轨迹的常态，应放行
        let mut traj = numbered_traj(3);
        traj.frames[2].sites_mut()[0].label = Some("Fe_x".into());
        assert!(traj.check_same_atoms().is_ok(), "只有 label 不同应放行");
        assert!(Trajectory::new().check_same_atoms().is_ok(), "空轨迹应放行");

        // 末帧少原子：截断 dump 的样子
        let mut short = numbered_traj(3);
        short.frames[2].set_atoms([]);
        let msg = short.check_same_atoms().unwrap_err().to_string();
        assert!(msg.contains("frame 2 has 0 atoms but frame 0 has 1"), "应点名帧与原子数：{msg}");

        // 原子数相同但元素换了位置
        let mut swapped = numbered_traj(3);
        swapped.frames[1].sites_mut()[0].element = "O".into();
        let msg = swapped.check_same_atoms().unwrap_err().to_string();
        assert!(msg.contains("atom 0 is O in frame 1 but Fe in frame 0"), "应点名原子与元素：{msg}");
    }

    fn picked(traj: &Trajectory) -> Vec<usize> {
        traj.frames.iter().map(|f| f.positions()[0].x as usize).collect()
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

    /// `select` and `tail` carry the metadata over unchanged, as `subset` does
    #[test]
    fn test_select_keeps_metadata() {
        let mut traj = numbered_traj(5);
        traj.metadata.title = Some("run".into());
        assert_eq!(traj.select(1, Some(3), 2).metadata.title.as_deref(), Some("run"));
        assert_eq!(traj.tail(2).metadata.title.as_deref(), Some("run"));
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
    fn test_first_last() {
        let mut traj = Trajectory::new();
        traj.add_frame(make_frame(0.0));
        traj.add_frame(make_frame(5.0));
        assert!((traj.first().unwrap().atom(0).position.x).abs() < 1e-10);
        assert!((traj.last().unwrap().atom(0).position.x - 5.0).abs() < 1e-10);
    }
}
