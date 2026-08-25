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

use ferro_core::error::ChemError;
use ferro_core::{select_range, spread_range, Table, Trajectory};

/// A quality criterion; the bit position is its slot in [`FrameVerdict::flags`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Criterion {
    /// Largest force VECTOR magnitude over the atoms of the frame
    Force,
    /// Largest absolute value among the 9 stress components
    Stress,
}

impl Criterion {
    pub const ALL: [Criterion; 2] = [Criterion::Force, Criterion::Stress];

    pub fn name(self) -> &'static str {
        match self {
            Criterion::Force => "force",
            Criterion::Stress => "stress",
        }
    }

    fn bit(self) -> u32 {
        match self {
            Criterion::Force => 1,
            Criterion::Stress => 2,
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
}

impl Default for FilterParams {
    fn default() -> Self {
        Self { f_max: 0.0, s_max: 0.0, start: 0, end: None, stride: 1, number: None }
    }
}

/// What every criterion said about one frame.
#[derive(Debug, Clone, PartialEq)]
pub struct FrameVerdict {
    /// Index in the input trajectory
    pub index: usize,
    pub max_force: Option<f64>,
    pub max_stress: Option<f64>,
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

    // 逐帧算出各判据的量与判定，全部帧都算 —— 交叉表要的是「每个判据单独
    // 判坏多少」，只在存活帧上算就永远看不出判据之间的重叠
    let verdicts: Vec<FrameVerdict> = traj
        .frames
        .iter()
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
            FrameVerdict { index: i, max_force, max_stress, flags }
        })
        .collect();

    let mut funnel = vec![("input".to_string(), n)];

    let mut survivors: Vec<usize> = (0..n).collect();
    for c in Criterion::ALL {
        let enabled = match c {
            Criterion::Force => params.f_max > 0.0,
            Criterion::Stress => params.s_max > 0.0,
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
    let keep: Vec<usize> = picked.iter().map(|&p| survivors[p]).collect();
    funnel.push(("range".to_string(), keep.len()));

    Ok(FilterResult { n_input: n, keep, verdicts, funnel, params: params.clone() })
}

impl FilterResult {
    /// Frames each criterion flagged, and how many of those no other criterion flagged.
    pub fn cross_tab(&self) -> Vec<(Criterion, usize, usize)> {
        Criterion::ALL
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
        for (i, &a) in Criterion::ALL.iter().enumerate() {
            for &b in &Criterion::ALL[i + 1..] {
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
            .push_num("kept", self.funnel.iter().map(|(_, k)| *k as f64).collect());

        let ct = self.cross_tab();
        let mut cross = Table::new();
        cross.meta_line("exclusive = flagged by this criterion and by no other");
        cross
            .push_text("criterion", ct.iter().map(|(c, _, _)| c.name().to_string()).collect())
            .push_num("flagged", ct.iter().map(|(_, f, _)| *f as f64).collect())
            .push_num("exclusive", ct.iter().map(|(_, _, e)| *e as f64).collect());

        let ov = self.overlaps();
        let mut overlap = Table::new();
        overlap
            .push_text("a", ov.iter().map(|(a, _, _)| a.name().to_string()).collect())
            .push_text("b", ov.iter().map(|(_, b, _)| b.name().to_string()).collect())
            .push_num("both", ov.iter().map(|(_, _, n)| *n as f64).collect());

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
            format!("s_max      = {} eV/Ang^3", off(p.s_max)),
            format!("range      = [{}:{}]", p.start, p.end.map(|e| e.to_string()).unwrap_or_default()),
        ];
        match p.number {
            Some(k) => v.push(format!("number     = {k}")),
            None => v.push(format!("stride     = {}", p.stride)),
        }
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
