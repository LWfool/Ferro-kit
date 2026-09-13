//! Helpers shared by more than two `md` analyses.

use ferro_core::{Frame, Trajectory};
use nalgebra::Vector3;

/// Build a time-averaged representative frame for the cube file header.
///
/// Frames whose atom count differs from the first frame are skipped rather than
/// truncated — a changed count means the trajectory is not one system.
pub(crate) fn build_avg_frame(traj: &Trajectory) -> Frame {
    let ref_f = traj.frames.first().unwrap();
    let n = ref_f.atoms.len();
    let mut pos_sum = vec![Vector3::<f64>::zeros(); n];
    let mut valid = 0usize;
    for frame in &traj.frames {
        if frame.atoms.len() != n { continue; }
        for (s, a) in pos_sum.iter_mut().zip(frame.atoms.iter()) {
            *s += a.position;
        }
        valid += 1;
    }
    let mut out = ref_f.clone();
    if valid > 0 {
        for (a, s) in out.atoms.iter_mut().zip(pos_sum.iter()) {
            a.position = s / valid as f64;
        }
    }
    out
}
