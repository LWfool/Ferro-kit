//! Helpers shared by more than two `md` analyses.

use ferro_core::{Cell, Frame, Trajectory};
use nalgebra::Vector3;

/// Build a time-averaged representative frame for the cube file header.
///
/// Frames whose atom count differs from the first frame are skipped rather than
/// truncated — a changed count means the trajectory is not one system.
///
/// Averaging happens in **unwrapped fractional coordinates**: an atom that crosses a
/// periodic boundary has Cartesian components that jump between 0 and L, and averaging
/// those puts it in the middle of the box.
///
/// Every frame is converted through the **reference frame's** cell, not its own.  That
/// keeps this a pure bug fix: a fractional mean taken under one fixed matrix maps back to
/// exactly the Cartesian mean, so atoms that never cross a boundary come out bit-for-bit
/// unchanged and only the wrapped ones move.  Using each frame's own cell would also
/// change what the average *means* under NPT — that is the open question recorded for
/// `cube_density`'s reference-frame volume normalisation, and it is not settled here.
///
/// Without a cell there is no periodicity and nothing to unwrap, so that case keeps the
/// plain Cartesian mean.
pub(crate) fn build_avg_frame(traj: &Trajectory) -> Frame {
    let ref_f = traj.frames.first().unwrap();
    let n = ref_f.atoms.len();
    let mut out = ref_f.clone();

    let frames: Vec<&Frame> = traj.frames.iter().filter(|f| f.atoms.len() == n).collect();
    if frames.is_empty() {
        return out;
    }
    let inv_n = 1.0 / frames.len() as f64;

    if let (Some(cell), Some(sum)) = (
        ref_f.cell.as_ref(),
        ref_f.cell.as_ref().and_then(|c| unwrapped_frac_sum(&frames, n, c)),
    ) {
        for (a, s) in out.atoms.iter_mut().zip(sum.iter()) {
            a.position = cell.fractional_to_cartesian(s * inv_n);
        }
        return out;
    }

    let mut pos_sum = vec![Vector3::<f64>::zeros(); n];
    for frame in &frames {
        for (s, a) in pos_sum.iter_mut().zip(frame.atoms.iter()) {
            *s += a.position;
        }
    }
    for (a, s) in out.atoms.iter_mut().zip(pos_sum.iter()) {
        a.position = s * inv_n;
    }
    out
}

/// Per-atom sum of the unwrapped fractional coordinates under `cell`, or `None` if the
/// cell is singular (the caller then falls back to the Cartesian mean).
///
/// Unwrapping is relative to the *previous* frame, not the first, so an atom that drifts
/// across several periodic images stays continuous — the same scheme as
/// [`super::cube_jump`]'s `unwrap_single`.
fn unwrapped_frac_sum(frames: &[&Frame], n: usize, cell: &Cell) -> Option<Vec<Vector3<f64>>> {
    let mut sum = vec![Vector3::<f64>::zeros(); n];
    let mut prev: Vec<Vector3<f64>> = Vec::with_capacity(n);
    for (fi, frame) in frames.iter().enumerate() {
        for (ai, atom) in frame.atoms.iter().enumerate() {
            let mut f = cell.cartesian_to_fractional(atom.position).ok()?;
            if fi == 0 {
                prev.push(f);
            } else {
                for k in 0..3 {
                    f[k] -= (f[k] - prev[ai][k]).round();
                }
                prev[ai] = f;
            }
            sum[ai] += f;
        }
    }
    Some(sum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::Atom;
    use nalgebra::Matrix3;

    /// One atom sitting either side of the x = 0 face of a 10 Å cube.
    fn boundary_crosser(with_cell: bool) -> Trajectory {
        let cube = Cell::from_matrix(Matrix3::new(10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0));
        let mut traj = Trajectory::new();
        for x in [0.5, 9.5] {
            let mut f = Frame::new();
            f.add_atom(Atom::new("O", Vector3::new(x, 5.0, 5.0)));
            if with_cell {
                f.cell = Some(cube.clone());
            }
            traj.frames.push(f);
        }
        traj
    }

    #[test]
    fn an_atom_crossing_the_boundary_does_not_average_to_the_box_centre() {
        let avg = build_avg_frame(&boundary_crosser(true));
        // 0.5 与 9.5 隔着 x = 0 这个面，最短路径的中点是 0（等价于 10），不是 5。
        // 直接对笛卡尔分量取平均会给出 5.0 —— 原子被放到了盒子正中央。
        let x = avg.atoms[0].position.x;
        assert!(x.abs() < 1e-10 || (x - 10.0).abs() < 1e-10, "got x = {x}, expected 0 or 10");
        assert!((avg.atoms[0].position.y - 5.0).abs() < 1e-10);
    }

    #[test]
    fn without_a_cell_the_cartesian_mean_is_kept() {
        // 无周期性就没有「跨边界」可言，5.0 是这两个位置唯一讲得通的平均
        let avg = build_avg_frame(&boundary_crosser(false));
        assert!((avg.atoms[0].position.x - 5.0).abs() < 1e-10);
    }

    #[test]
    fn frames_with_a_different_atom_count_are_skipped() {
        let mut traj = boundary_crosser(true);
        let mut odd = Frame::new();
        odd.add_atom(Atom::new("O", Vector3::new(1.0, 1.0, 1.0)));
        odd.add_atom(Atom::new("O", Vector3::new(2.0, 2.0, 2.0)));
        traj.frames.push(odd);
        let avg = build_avg_frame(&traj);
        assert_eq!(avg.atoms.len(), 1);
        let x = avg.atoms[0].position.x;
        assert!(x.abs() < 1e-10 || (x - 10.0).abs() < 1e-10, "got x = {x}");
    }
}
