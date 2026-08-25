//! Quality selection over a labelled trajectory destined for training.
//!
//! Pure computation: it decides WHICH frames to keep and never touches the
//! filesystem. The caller applies the verdict (`Trajectory::select` with the
//! kept indices) and writes the result.
//!
//! # The funnel
//!
//! Criteria run in a fixed order and each one narrows the survivors:
//!
//! ```text
//! all frames -> |F|max -> |sigma|max -> [start:end:stride|number]
//! ```
//!
//! **The range applies to the survivors, not to original frame numbers.**
//! `--start 10` means "the 10th frame that passed the quality criteria", which
//! is what makes a range meaningful after an unknown number of frames were
//! dropped. Original indices are what [`FilterResult::keep`] carries, so the
//! kept frames can always be traced back.
//!
//! # Why a cross-tabulation
//!
//! The funnel alone cannot tell whether a criterion earns its place: each step
//! reports on what the previous step left, so a criterion that only ever flags
//! frames another criterion already caught looks productive. [`FilterResult`]
//! therefore also records, per criterion, how many frames it flags on its own
//! and how many of those NO other criterion flagged — the exclusive count. A
//! criterion whose exclusive count is ~0 is redundant and can be switched off.

use std::collections::BTreeMap;

use ferro_core::error::ChemError;
use ferro_core::{select_range, spread_range, Table, Trajectory, TypeParams};
use rayon::prelude::*;

use super::geometry::{count_with_coordination, min_pair_distance};
use super::merge::shuffle_order;

/// A quality criterion; the bit position is its slot in [`FrameVerdict::flags`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Criterion {
    /// Largest force VECTOR magnitude over the atoms of the frame
    Force,
    /// Largest absolute value among the 9 stress components
    Stress,
    /// Smallest O-O distance in the frame, minimum image
    OoMin,
    /// Frames holding no 6-coordinated Al at all
    Al6,
}

impl Criterion {
    pub const ALL: [Criterion; 4] =
        [Criterion::Force, Criterion::Stress, Criterion::OoMin, Criterion::Al6];

    pub fn name(self) -> &'static str {
        match self {
            Criterion::Force => "force",
            Criterion::Stress => "stress",
            Criterion::OoMin => "oo_min",
            Criterion::Al6 => "al6",
        }
    }

    fn bit(self) -> u32 {
        match self {
            Criterion::Force => 1,
            Criterion::Stress => 2,
            Criterion::OoMin => 4,
            Criterion::Al6 => 8,
        }
    }
}

/// Thresholds and range, all in ferro's internal units.
///
/// A threshold of `0.0` switches that criterion off, matching the convention of
/// the reference Python implementation: an explicit zero says "do not judge",
/// which a small positive number could not express.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterParams {
    /// Drop frames whose largest force magnitude exceeds this (eV/Å); 0 = off
    pub f_max: f64,
    /// Drop frames whose largest |stress component| exceeds this (eV/Å³); 0 = off
    pub s_max: f64,
    /// First survivor to take (0-based, inclusive)
    pub start: usize,
    /// Last survivor to take (0-based, INCLUSIVE); `None` = to the end
    pub end: Option<usize>,
    /// Take every Nth survivor within the range
    pub stride: usize,
    /// Take this many survivors, spread evenly; mutually exclusive with stride
    pub number: Option<usize>,
    /// Drop frames whose smallest O-O distance is below this (Å); 0 = off
    pub oo_min: f64,
    /// Drop frames holding no 6-coordinated Al; the Al-O cutoff (Å). `None` = off
    pub al6_rcut: Option<f64>,
    /// Shuffle the kept frames with this seed once every criterion has run.
    /// `None` keeps them in trajectory order.
    pub shuffle: Option<u64>,
}

impl Default for FilterParams {
    fn default() -> Self {
        Self {
            f_max: 0.0,
            s_max: 0.0,
            start: 0,
            end: None,
            stride: 1,
            number: None,
            oo_min: 0.0,
            al6_rcut: None,
            shuffle: None,
        }
    }
}

/// What every criterion said about one frame.
#[derive(Debug, Clone, PartialEq)]
pub struct FrameVerdict {
    /// Index in the input trajectory
    pub index: usize,
    pub max_force: Option<f64>,
    pub max_stress: Option<f64>,
    /// Smallest O-O distance, only computed when that criterion is on
    pub min_oo: Option<f64>,
    /// Number of 6-coordinated Al, only computed when that criterion is on
    pub n_al6: Option<usize>,
    /// Bit set of the criteria that flagged this frame
    pub flags: u32,
}

impl FrameVerdict {
    pub fn flagged_by(&self, c: Criterion) -> bool {
        self.flags & c.bit() != 0
    }
    pub fn is_clean(&self) -> bool {
        self.flags == 0
    }
}

#[derive(Debug, Clone)]
pub struct FilterResult {
    pub n_input: usize,
    /// Indices INTO THE INPUT trajectory that survived every step, in order
    pub keep: Vec<usize>,
    pub verdicts: Vec<FrameVerdict>,
    /// `(step name, frames remaining after it)`, in execution order
    pub funnel: Vec<(String, usize)>,
    pub params: FilterParams,
}

/// Decides which frames of `traj` to keep.
///
/// Fails before judging anything when a threshold is set but the trajectory
/// lacks that label — a missing force array is a broken dataset, not a frame
/// that happens to pass.
pub fn filter_frames(traj: &Trajectory, params: &FilterParams) -> Result<FilterResult, ChemError> {
    let n = traj.n_frames();
    if params.stride != 1 && params.number.is_some() {
        return Err(ChemError::ValidationError(
            "stride and number are two ways of saying the same thing; give only one".into(),
        ));
    }
    if params.f_max > 0.0 && traj.frames.iter().any(|f| f.forces.is_none()) {
        return Err(ChemError::ValidationError(
            "a force threshold was given but some frames carry no forces".into(),
        ));
    }
    if params.s_max > 0.0 && traj.frames.iter().any(|f| f.stress.is_none()) {
        return Err(ChemError::ValidationError(
            "a stress threshold was given but some frames carry no stress".into(),
        ));
    }
    let geometric = params.oo_min > 0.0 || params.al6_rcut.is_some();
    if geometric && traj.frames.iter().any(|f| f.cell.is_none()) {
        return Err(ChemError::ValidationError(
            "a geometric criterion was given but some frames carry no cell".into(),
        ));
    }
    // 最小镜像只在关心的距离小于最小面间距的一半时严格成立
    if let Some(f) = traj.frames.first() {
        if let (true, Some(cell)) = (geometric, f.cell.as_ref()) {
            let bound = cell.minimum_image_cutoff()?;
            let want = params.oo_min.max(params.al6_rcut.unwrap_or(0.0));
            if want > bound {
                return Err(ChemError::ValidationError(format!(
                    "cutoff {want:.3} A exceeds the minimum-image bound {bound:.3} A of the cell"
                )));
            }
        }
    }

    // Al6 走 network 的分类器，故只需给出 Al-O 一个截断；
    // Al 不在默认 Qn 名单 {B,P,Si} 里，落在「非 Qn 形成子」分支，cn 即配位数
    let type_params = params.al6_rcut.map(|r| {
        let mut cut = BTreeMap::new();
        cut.insert(("Al".to_string(), "O".to_string()), r);
        TypeParams::new(cut, Default::default())
    });

    // 逐帧算出各判据的量与判定，全部帧都算 —— 交叉表要的是「每个判据单独
    // 判坏多少」，只在存活帧上算就永远看不出判据之间的重叠
    let verdicts: Vec<FrameVerdict> = traj
        .frames
        .par_iter()
        .enumerate()
        .map(|(i, f)| {
            let max_force = f.forces.as_ref().map(|v| {
                v.iter().map(|x| x.norm()).fold(0.0_f64, f64::max)
            });
            let max_stress = f
                .stress
                .map(|s| s.iter().map(|x| x.abs()).fold(0.0_f64, f64::max));
            let mut flags = 0;
            if params.f_max > 0.0 {
                if let Some(m) = max_force {
                    if m > params.f_max {
                        flags |= Criterion::Force.bit();
                    }
                }
            }
            if params.s_max > 0.0 {
                if let Some(m) = max_stress {
                    if m > params.s_max {
                        flags |= Criterion::Stress.bit();
                    }
                }
            }
            let mut min_oo = None;
            let mut n_al6 = None;
            if let Some(cell) = f.cell.as_ref() {
                if params.oo_min > 0.0 {
                    let d = min_pair_distance(f, cell, "O", "O");
                    if let Some(d) = d {
                        if d < params.oo_min {
                            flags |= Criterion::OoMin.bit();
                        }
                    }
                    min_oo = d;
                }
                if let Some(tp) = &type_params {
                    let n = count_with_coordination(f, cell, tp, "Al", 6);
                    // 「保留含 Al6 的帧」与「删除不含 Al6 的帧」是同一件事，
                    // 表达成后者，判据语义就和其余三条一致，交叉表也不必分裂
                    if n == 0 {
                        flags |= Criterion::Al6.bit();
                    }
                    n_al6 = Some(n);
                }
            }
            FrameVerdict { index: i, max_force, max_stress, min_oo, n_al6, flags }
        })
        .collect();

    let mut funnel = vec![("input".to_string(), n)];

    let mut survivors: Vec<usize> = (0..n).collect();
    for c in Criterion::ALL {
        let enabled = match c {
            Criterion::Force => params.f_max > 0.0,
            Criterion::Stress => params.s_max > 0.0,
            Criterion::OoMin => params.oo_min > 0.0,
            Criterion::Al6 => params.al6_rcut.is_some(),
        };
        if !enabled {
            continue;
        }
        survivors.retain(|&i| !verdicts[i].flagged_by(c));
        funnel.push((c.name().to_string(), survivors.len()));
    }

    // 区间与抽样作用在**存活序列的位置**上，不是原始帧号
    let picked: Vec<usize> = match params.number {
        Some(k) => spread_range(survivors.len(), params.start, params.end, k),
        None => select_range(survivors.len(), params.start, params.end, params.stride),
    };
    let mut keep: Vec<usize> = picked.iter().map(|&p| survivors[p]).collect();
    funnel.push(("range".to_string(), keep.len()));

    // 打乱**在所有判据与抽帧之后**：--stride / --number 是「每隔多久取一帧」，
    // 在乱序上说不通；打乱一旦发生，时间序就再也取不回来了
    if let Some(seed) = params.shuffle {
        let order = shuffle_order(keep.len(), seed);
        keep = order.into_iter().map(|i| keep[i]).collect();
    }

    Ok(FilterResult { n_input: n, keep, verdicts, funnel, params: params.clone() })
}

impl FilterResult {
    /// The criteria that were switched on for this run, in funnel order.
    ///
    /// Reporting the disabled ones too would fill the tables with rows of zeros
    /// and bury the counts that mean something.
    pub fn enabled(&self) -> Vec<Criterion> {
        let p = &self.params;
        Criterion::ALL
            .into_iter()
            .filter(|c| match c {
                Criterion::Force => p.f_max > 0.0,
                Criterion::Stress => p.s_max > 0.0,
                Criterion::OoMin => p.oo_min > 0.0,
                Criterion::Al6 => p.al6_rcut.is_some(),
            })
            .collect()
    }

    /// Frames each criterion flagged, and how many of those no other criterion flagged.
    pub fn cross_tab(&self) -> Vec<(Criterion, usize, usize)> {
        self.enabled()
            .iter()
            .map(|&c| {
                let flagged = self.verdicts.iter().filter(|v| v.flagged_by(c)).count();
                let exclusive = self
                    .verdicts
                    .iter()
                    .filter(|v| v.flags == c.bit())
                    .count();
                (c, flagged, exclusive)
            })
            .collect()
    }

    /// Frames flagged by both members of each criterion pair.
    pub fn overlaps(&self) -> Vec<(Criterion, Criterion, usize)> {
        let mut out = Vec::new();
        let on = self.enabled();
        for (i, &a) in on.iter().enumerate() {
            for &b in &on[i + 1..] {
                let n = self
                    .verdicts
                    .iter()
                    .filter(|v| v.flagged_by(a) && v.flagged_by(b))
                    .count();
                out.push((a, b, n));
            }
        }
        out
    }

    pub fn to_tables(&self) -> Vec<(String, Table)> {
        let mut funnel = Table::new();
        funnel.meta_line("frames remaining after each step, in execution order");
        funnel
            .push_text("step", self.funnel.iter().map(|(s, _)| s.clone()).collect())
            .push_text("kept", self.funnel.iter().map(|(_, k)| k.to_string()).collect());

        let ct = self.cross_tab();
        let mut cross = Table::new();
        cross.meta_line("exclusive = flagged by this criterion and by no other");
        cross
            .push_text("criterion", ct.iter().map(|(c, _, _)| c.name().to_string()).collect())
            .push_text("flagged", ct.iter().map(|(_, f, _)| f.to_string()).collect())
            .push_text("exclusive", ct.iter().map(|(_, _, e)| e.to_string()).collect());

        let ov = self.overlaps();
        let mut overlap = Table::new();
        overlap
            .push_text("a", ov.iter().map(|(a, _, _)| a.name().to_string()).collect())
            .push_text("b", ov.iter().map(|(_, b, _)| b.name().to_string()).collect())
            .push_text("both", ov.iter().map(|(_, _, n)| n.to_string()).collect());

        vec![
            ("funnel".to_string(), funnel),
            ("criteria".to_string(), cross),
            ("overlap".to_string(), overlap),
        ]
    }

    pub fn meta_lines(&self) -> Vec<String> {
        let p = &self.params;
        let off = |v: f64| if v > 0.0 { format!("{v}") } else { "off".to_string() };
        let mut v = vec![
            format!("frames in  = {}", self.n_input),
            format!("frames out = {}", self.keep.len()),
            format!("f_max      = {} eV/Ang", off(p.f_max)),
            match p.s_max > 0.0 {
                true => format!(
                    "s_max      = {:.6} eV/Ang^3  ({:.4} GPa)",
                    p.s_max,
                    ferro_core::units::convert_pressure(
                        p.s_max,
                        ferro_core::units::PressureUnit::EVPerAng3,
                        ferro_core::units::PressureUnit::GPa,
                    )
                ),
                false => "s_max      = off".to_string(),
            },
            format!("range      = [{}:{}]", p.start, p.end.map(|e| e.to_string()).unwrap_or_default()),
        ];
        match p.number {
            Some(k) => v.push(format!("number     = {k}")),
            None => v.push(format!("stride     = {}", p.stride)),
        }
        v.push(format!("oo_min     = {} Ang", off(p.oo_min)));
        v.push(match p.al6_rcut {
            Some(r) => format!("al6_rcut   = {r:.3} Ang"),
            None => "al6_rcut   = off".to_string(),
        });
        v.push(match p.shuffle {
            Some(seed) => format!("shuffle    = seed {seed}"),
            None => "shuffle    = off (trajectory order kept)".to_string(),
        });
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame};
    use nalgebra::{Matrix3, Vector3};

    /// `forces[i]` and `stresses[i]` give frame i its peak values.
    fn traj_of(forces: &[f64], stresses: &[f64]) -> Trajectory {
        let cell = Cell::from_matrix(Matrix3::identity() * 5.0);
        let frames = forces
            .iter()
            .zip(stresses)
            .map(|(&f, &s)| {
                let mut fr = Frame::with_cell(cell.clone(), [true; 3]);
                fr.atoms = vec![
                    Atom::new("O", Vector3::zeros()),
                    Atom::new("O", Vector3::new(1.0, 0.0, 0.0)),
                ];
                // 第二个原子拿峰值，确保取的是逐原子最大而不是第一个
                fr.forces = Some(vec![Vector3::zeros(), Vector3::new(f, 0.0, 0.0)]);
                fr.stress = Some(Matrix3::from_row_slice(&[
                    0.0, 0.0, 0.0, 0.0, 0.0, s, 0.0, s, 0.0,
                ]));
                fr
            })
            .collect();
        Trajectory { frames, metadata: Default::default() }
    }


    /// A frame with two O too close, and one where the same pair is fine.
    fn oo_traj(gaps: &[f64]) -> Trajectory {
        let cell = Cell::from_matrix(Matrix3::identity() * 12.0);
        let frames = gaps
            .iter()
            .map(|&g| {
                let mut f = Frame::with_cell(cell.clone(), [true; 3]);
                f.atoms = vec![
                    Atom::new("O", Vector3::new(1.0, 1.0, 1.0)),
                    Atom::new("O", Vector3::new(1.0 + g, 1.0, 1.0)),
                ];
                f
            })
            .collect();
        Trajectory { frames, metadata: Default::default() }
    }

    #[test]
    fn oo_min_drops_the_close_contact_frames() {
        let t = oo_traj(&[2.5, 1.8, 2.1]);
        let p = FilterParams { oo_min: 2.0, ..Default::default() };
        let r = filter_frames(&t, &p).unwrap();
        assert_eq!(r.keep, vec![0, 2]);
        assert!(r.verdicts[1].flagged_by(Criterion::OoMin));
        // 判据打开时该量被记录下来，供只读模式画分布
        assert!((r.verdicts[1].min_oo.unwrap() - 1.8).abs() < 1e-12);
        assert!(r.verdicts[0].n_al6.is_none(), "al6 未开启就不该计算");
    }

    /// "Keep frames containing Al6" and "drop frames containing none" are the
    /// same rule; expressing it as the latter keeps one sense for all criteria.
    #[test]
    fn al6_drops_frames_without_any_six_coordinated_al() {
        let cell = Cell::from_matrix(Matrix3::identity() * 20.0);
        let mut frames = Vec::new();
        for n_o in [6usize, 4] {
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            let mut atoms = vec![Atom::new("Al", Vector3::new(5.0, 5.0, 5.0))];
            let dirs = [
                [1.9, 0.0, 0.0], [-1.9, 0.0, 0.0], [0.0, 1.9, 0.0],
                [0.0, -1.9, 0.0], [0.0, 0.0, 1.9], [0.0, 0.0, -1.9],
            ];
            for d in dirs.iter().take(n_o) {
                atoms.push(Atom::new("O", Vector3::new(5.0 + d[0], 5.0 + d[1], 5.0 + d[2])));
            }
            f.atoms = atoms;
            frames.push(f);
        }
        let t = Trajectory { frames, metadata: Default::default() };
        let p = FilterParams { al6_rcut: Some(2.4), ..Default::default() };
        let r = filter_frames(&t, &p).unwrap();
        assert_eq!(r.keep, vec![0]);
        assert_eq!(r.verdicts[0].n_al6, Some(1));
        assert_eq!(r.verdicts[1].n_al6, Some(0));
        assert!(r.verdicts[1].flagged_by(Criterion::Al6));
    }

    #[test]
    fn a_cutoff_past_the_minimum_image_bound_is_an_error() {
        let t = oo_traj(&[2.5]);
        // 盒子 12 Å，最小镜像上界 6 Å
        let p = FilterParams { oo_min: 7.0, ..Default::default() };
        let err = filter_frames(&t, &p).unwrap_err().to_string();
        assert!(err.contains("minimum-image"), "{err}");
    }

    #[test]
    fn disabled_criteria_stay_out_of_the_tables() {
        let t = oo_traj(&[2.5, 1.8]);
        let p = FilterParams { oo_min: 2.0, ..Default::default() };
        let r = filter_frames(&t, &p).unwrap();
        assert_eq!(r.enabled(), vec![Criterion::OoMin]);
        assert_eq!(r.cross_tab().len(), 1);
        assert!(r.overlaps().is_empty(), "只有一条判据时没有两两重叠");
    }

    #[test]
    fn force_and_stress_thresholds_drop_their_own_frames() {
        //           帧:   0    1     2    3
        let t = traj_of(&[1.0, 99.0, 1.0, 1.0], &[0.1, 0.1, 9.9, 0.1]);
        let p = FilterParams { f_max: 20.0, s_max: 1.0, ..Default::default() };
        let r = filter_frames(&t, &p).unwrap();
        assert_eq!(r.keep, vec![0, 3]);
        assert!(r.verdicts[1].flagged_by(Criterion::Force));
        assert!(r.verdicts[2].flagged_by(Criterion::Stress));
        assert!(r.verdicts[0].is_clean());
    }

    #[test]
    fn a_zero_threshold_switches_the_criterion_off() {
        let t = traj_of(&[1.0, 99.0], &[0.1, 0.1]);
        let r = filter_frames(&t, &FilterParams::default()).unwrap();
        assert_eq!(r.keep, vec![0, 1]);
        assert!(r.verdicts[1].is_clean());
    }

    /// The exclusive count is what tells a redundant criterion from a useful one.
    #[test]
    fn cross_tab_separates_exclusive_from_overlapping_hits() {
        //  帧 0 干净; 帧 1 只超力; 帧 2 只超应力; 帧 3 两者都超
        let t = traj_of(&[1.0, 99.0, 1.0, 99.0], &[0.1, 0.1, 9.9, 9.9]);
        let p = FilterParams { f_max: 20.0, s_max: 1.0, ..Default::default() };
        let r = filter_frames(&t, &p).unwrap();
        assert_eq!(r.keep, vec![0]);
        let ct = r.cross_tab();
        assert_eq!(ct[0], (Criterion::Force, 2, 1));
        assert_eq!(ct[1], (Criterion::Stress, 2, 1));
        assert_eq!(r.overlaps(), vec![(Criterion::Force, Criterion::Stress, 1)]);
    }

    /// The range counts survivors, not original frame numbers.
    #[test]
    fn the_range_applies_to_survivors_not_to_original_indices() {
        // 帧 1 被力判掉，存活序列是 [0, 2, 3, 4]
        let t = traj_of(&[1.0, 99.0, 1.0, 1.0, 1.0], &[0.1; 5]);
        let p = FilterParams {
            f_max: 20.0,
            start: 1,
            end: Some(2),
            ..Default::default()
        };
        let r = filter_frames(&t, &p).unwrap();
        // 存活序列的第 1..=2 个是原始帧 2 和 3
        assert_eq!(r.keep, vec![2, 3]);
    }

    #[test]
    fn number_spreads_over_survivors_and_takes_both_ends() {
        let t = traj_of(&[1.0; 10], &[0.1; 10]);
        let p = FilterParams { number: Some(3), ..Default::default() };
        let r = filter_frames(&t, &p).unwrap();
        assert_eq!(r.keep, vec![0, 5, 9]);
    }

    #[test]
    fn shuffle_permutes_the_kept_frames_after_everything_else() {
        let t = traj_of(&[1.0, 99.0, 1.0, 1.0, 1.0, 1.0], &[0.1; 6]);
        let base = FilterParams { f_max: 20.0, ..Default::default() };
        let ordered = filter_frames(&t, &base).unwrap();
        assert_eq!(ordered.keep, vec![0, 2, 3, 4, 5]);

        let p = FilterParams { shuffle: Some(7), ..base.clone() };
        let mixed = filter_frames(&t, &p).unwrap();
        // 同一批帧，顺序不同
        let mut sorted = mixed.keep.clone();
        sorted.sort();
        assert_eq!(sorted, ordered.keep);
        assert_ne!(mixed.keep, ordered.keep);
        // 同一 seed 可复现
        assert_eq!(filter_frames(&t, &p).unwrap().keep, mixed.keep);
    }

    /// Sampling has to see time order, so the shuffle runs last.
    #[test]
    fn stride_still_samples_in_time_before_the_shuffle() {
        let t = traj_of(&[1.0; 10], &[0.1; 10]);
        let p = FilterParams { stride: 3, shuffle: Some(1), ..Default::default() };
        let r = filter_frames(&t, &p).unwrap();
        let mut got = r.keep.clone();
        got.sort();
        // 等间隔取到的仍是 0/3/6/9，只是写出顺序被打乱
        assert_eq!(got, vec![0, 3, 6, 9]);
    }

    #[test]
    fn stride_and_number_together_are_an_error() {
        let t = traj_of(&[1.0; 4], &[0.1; 4]);
        let p = FilterParams { stride: 2, number: Some(2), ..Default::default() };
        assert!(filter_frames(&t, &p).is_err());
    }

    #[test]
    fn a_threshold_without_the_label_fails_before_judging() {
        let mut t = traj_of(&[1.0, 1.0], &[0.1, 0.1]);
        t.frames[1].forces = None;
        let p = FilterParams { f_max: 20.0, ..Default::default() };
        let err = filter_frames(&t, &p).unwrap_err().to_string();
        assert!(err.contains("forces"), "{err}");
    }

    #[test]
    fn funnel_records_every_enabled_step() {
        let t = traj_of(&[1.0, 99.0, 1.0], &[0.1, 0.1, 9.9]);
        let p = FilterParams { f_max: 20.0, s_max: 1.0, ..Default::default() };
        let r = filter_frames(&t, &p).unwrap();
        let steps: Vec<&str> = r.funnel.iter().map(|(s, _)| s.as_str()).collect();
        assert_eq!(steps, vec!["input", "force", "stress", "range"]);
        assert_eq!(r.funnel.iter().map(|(_, k)| *k).collect::<Vec<_>>(), vec![3, 2, 1, 1]);
        assert_eq!(r.to_tables().len(), 3);
    }
}
