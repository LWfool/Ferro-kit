//! Helpers shared by more than two `md` analyses.

use ferro_core::error::{ChemError, Result};
use ferro_core::{Cell, CubeData, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use ndarray::Array3;

/// The binning rule of every histogram analysis (CLAUDE.md, "Histogram bins"): the bin
/// count is `ceil((hi - lo) / step)` and the width actually used is `(hi - lo) / count`,
/// so the last bin ends exactly at `hi` and no bin is wider than requested. A ratio
/// within 1e-9 (relative) of a whole number counts as that number, or `10.004 / 0.002`
/// would round up to one bin too many.
///
/// The three analyses used to differ: gr dropped a partial last bin (floor), vanhove and
/// angle kept it but normalised it as a full one, so with `--r-max 1.0 --dr 0.07` the
/// last vanhove bin read 3.5 times too low. Returns 0 bins for `hi <= lo`.
pub(crate) fn bins(lo: f64, hi: f64, step: f64) -> (usize, f64) {
    let ratio = (hi - lo) / step;
    if ratio.is_nan() || ratio <= 0.0 {
        return (0, step);
    }
    let n = (ratio - ratio * 1e-9).ceil().max(1.0) as usize;
    // A whole number of steps keeps the request as given: (hi - lo) / n would be one ulp
    // off (0.0020000000000000005)
    if (ratio - n as f64).abs() <= ratio * 1e-9 {
        return (n, step);
    }
    (n, (hi - lo) / n as f64)
}

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
/// Maps `items` in parallel and hands each result to `consume` in item order.
///
/// The map runs in batches of a few items per thread, so only one batch of results
/// is alive at a time: collecting every frame's grid before summing held hundreds of
/// MB for a `map density` (review D-P1). Consuming in item order keeps a float sum
/// bit-identical to the sequential one, which a parallel `reduce` would not — its
/// grouping follows the thread schedule.
pub(crate) fn par_map_in_order<T: Sync, R: Send>(
    items: &[T],
    map: impl Fn(&T) -> R + Sync,
    mut consume: impl FnMut(R),
) {
    use rayon::prelude::*;
    let batch = 2 * rayon::current_num_threads().max(1);
    for chunk in items.chunks(batch) {
        let out: Vec<R> = chunk.par_iter().map(&map).collect();
        out.into_iter().for_each(&mut consume);
    }
}

/// Indices of the atoms of `frame` whose element is in `elements` (`None` = every atom).
///
/// Selection is made once, on one frame, and applied to the whole trajectory; callers
/// check `same_atoms` first. An empty selection is an error: an analysis of no atoms
/// would write zeros that look like a measurement.
pub(crate) fn select_atoms(frame: &Frame, elements: Option<&[String]>) -> Result<Vec<usize>> {
    let picked: Vec<usize> = frame.atoms().enumerate()
        .filter(|(_, a)| elements.is_none_or(|els| els.iter().any(|e| e == a.element)))
        .map(|(i, _)| i)
        .collect();
    if picked.is_empty() {
        return Err(ChemError::ValidationError(format!(
            "no atom matches the element filter {:?}", elements.unwrap_or_default()
        )));
    }
    Ok(picked)
}

/// Packs a grid accumulated on `ref_cell` into a cube file: voxel `i` along an axis is
/// `cell row / n` wide, and the values sit at voxel centres, so the origin is half a
/// voxel in (the cube format places point `i` at `origin + i·spacing`; review B-2).
/// Shared by `cube_density`, `cube_radius` and `cube_jump`.
///
/// The reference structure is the time-averaged frame with atoms that left the cell
/// along a periodic axis folded back into `[0, 1)`. **This is a special case: only cube
/// output does it** (review M5, user decision 2026-10-02). Everywhere else Ferro keeps
/// the coordinates the reader produced, neither shifted nor folded — folding on read
/// perturbs an NPT MSD (a ΔL-sized residue under TOR unwrapping) and subtracting `lo`
/// per frame adds a spurious drift, see `dev/issues.md`. A cube differs: the grid always
/// spans `[0, L)` from origin 0 and voxels are assigned from folded fractional
/// coordinates, so an unfolded reference would draw atoms outside the grid whenever a
/// LAMMPS box has `lo ≠ 0` (659 of them in `tests/70Z30P00A_NVT_5`). Densities are
/// unaffected.
pub(crate) fn cube_from_grid(traj: &Trajectory, ref_cell: &Cell, data: Array3<f64>) -> CubeData {
    let [nx, ny, nz] = [data.dim().0, data.dim().1, data.dim().2];
    let m = ref_cell.matrix;
    let spacing = Matrix3::new(
        m[(0, 0)] / nx as f64, m[(0, 1)] / nx as f64, m[(0, 2)] / nx as f64,
        m[(1, 0)] / ny as f64, m[(1, 1)] / ny as f64, m[(1, 2)] / ny as f64,
        m[(2, 0)] / nz as f64, m[(2, 1)] / nz as f64, m[(2, 2)] / nz as f64,
    );
    let mut frame = build_avg_frame(traj);
    frame.wrap_all();
    CubeData {
        frame,
        data: data.into_iter().collect(),
        shape: [nx, ny, nz],
        origin: spacing.transpose() * Vector3::repeat(0.5),
        spacing,
    }
}

pub(crate) fn build_avg_frame(traj: &Trajectory) -> Frame {
    let ref_f = traj.frames.first().unwrap();
    let n = ref_f.n_atoms();
    let mut out = ref_f.clone();

    let frames: Vec<&Frame> = traj.frames.iter().filter(|f| f.n_atoms() == n).collect();
    if frames.is_empty() {
        return out;
    }
    let inv_n = 1.0 / frames.len() as f64;

    if let (Some(cell), Some(sum)) = (
        ref_f.cell.as_ref(),
        ref_f.cell.as_ref().and_then(|c| unwrapped_frac_sum(&frames, n, c)),
    ) {
        for (p, s) in out.positions_mut().iter_mut().zip(sum.iter()) {
            *p = cell.fractional_to_cartesian(s * inv_n);
        }
        return out;
    }

    let mut pos_sum = vec![Vector3::<f64>::zeros(); n];
    for frame in &frames {
        for (s, a) in pos_sum.iter_mut().zip(frame.atoms()) {
            *s += a.position;
        }
    }
    for (p, s) in out.positions_mut().iter_mut().zip(pos_sum.iter()) {
        *p = s * inv_n;
    }
    out
}

/// Lattice-view unwrap of one atom's fractional series, in place: each step is moved
/// by the whole lattice vector that brings it nearest the previous one (`|Δ| > 0.5` is a
/// boundary crossing). Shared by `vanhove` and `cube_jump`; `msd` uses the TOR scheme
/// instead, which is the correct one under NPT — moving these two to it is on
/// `dev/plan.md`.
///
/// [`unwrapped_frac_sum`] applies the same rule frame by frame while summing, so it never
/// holds a whole series; it keeps its own loop for that reason.
pub(crate) fn unwrap_lattice(frac: &mut [[f64; 3]]) {
    for i in 1..frac.len() {
        let (prev, curr_and_later) = frac.split_at_mut(i);
        for (c, p) in curr_and_later[0].iter_mut().zip(prev[i - 1].iter()) {
            *c -= (*c - *p).round();
        }
    }
}

/// Per-atom sum of the unwrapped fractional coordinates under `cell`, or `None` if the
/// cell is singular (the caller then falls back to the Cartesian mean).
///
/// Unwrapping is relative to the *previous* frame, not the first, so an atom that drifts
/// across several periodic images stays continuous — the same rule as
/// [`unwrap_lattice`], applied while summing.
fn unwrapped_frac_sum(frames: &[&Frame], n: usize, cell: &Cell) -> Option<Vec<Vector3<f64>>> {
    let mut sum = vec![Vector3::<f64>::zeros(); n];
    let mut prev: Vec<Vector3<f64>> = Vec::with_capacity(n);
    for (fi, frame) in frames.iter().enumerate() {
        for (ai, atom) in frame.atoms().enumerate() {
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

/// 时间相关分析的前提：帧等间隔、不重复。轨迹带步号（LAMMPS dump、CP2K、OUTCAR）时
/// 逐帧核对；没有步号的格式无从查起，只能信调用方给的 `dt`。
///
/// 重启拼接的 dump 常把重启点写两次（同一个 TIMESTEP），或中途改过 dump 间隔；
/// 全原点平均把每一对帧都当成 `lag × dt` 相隔，这两种情形都会静默算错。
pub(crate) fn check_frame_spacing(traj: &Trajectory) -> Result<()> {
    let steps: Option<Vec<i64>> = traj.frames.iter().map(|f| f.step).collect();
    let Some(steps) = steps else { return Ok(()) };
    let Some(first_gap) = steps.get(1).map(|s1| s1 - steps[0]) else { return Ok(()) };
    for (k, w) in steps.windows(2).enumerate() {
        let gap = w[1] - w[0];
        if gap <= 0 {
            return Err(ChemError::ValidationError(format!(
                "frames {k} and {} have steps {} and {}: duplicated or out-of-order frames \
                 (a restart written twice?); remove them before a time-correlation analysis",
                k + 1, w[0], w[1]
            )));
        }
        if gap != first_gap {
            return Err(ChemError::ValidationError(format!(
                "frames are not evenly spaced: step gap {first_gap} at the start, {gap} between \
                 frames {k} and {} (steps {} → {}); dt assumes one constant spacing",
                k + 1, w[0], w[1]
            )));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::Atom;
    use nalgebra::Matrix3;

    /// Spacing rows are the cell rows over the grid counts, and the origin sits half a
    /// voxel in along each — checked on a triclinic cell, where the two differ per axis
    #[test]
    fn cube_from_grid_puts_values_at_voxel_centres() {
        let cell = Cell::from_matrix(Matrix3::new(8.0, 0.0, 0.0, 2.0, 6.0, 0.0, 1.0, 1.0, 5.0));
        let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
        frame.add_atom(Atom::new("O", Vector3::new(9.0, 1.0, 1.0)));
        let mut traj = Trajectory::new();
        traj.frames.push(frame);
        let cube = cube_from_grid(&traj, &cell, Array3::zeros((2, 4, 5)));
        assert_eq!(cube.shape, [2, 4, 5]);
        let want = Matrix3::new(4.0, 0.0, 0.0, 0.5, 1.5, 0.0, 0.2, 0.2, 1.0);
        assert!((cube.spacing - want).norm() < 1e-12, "{}", cube.spacing);
        assert!((cube.origin - Vector3::new(2.35, 0.85, 0.5)).norm() < 1e-12, "{}", cube.origin);
        // the reference atom past the a face is folded back into the cell
        assert!((cube.frame.positions()[0] - Vector3::new(1.0, 1.0, 1.0)).norm() < 1e-12);
    }

    /// A walk across the x face, twice, comes back continuous; y and z stay put
    #[test]
    fn unwrap_lattice_follows_an_atom_across_faces() {
        let mut f = [[0.9, 0.5, 0.5], [0.1, 0.5, 0.5], [0.3, 0.5, 0.5], [0.1, 0.5, 0.5], [0.9, 0.5, 0.5]];
        unwrap_lattice(&mut f);
        let x: Vec<f64> = f.iter().map(|p| p[0]).collect();
        for (got, want) in x.iter().zip([0.9, 1.1, 1.3, 1.1, 0.9]) {
            assert!((got - want).abs() < 1e-12, "{x:?}");
        }
        assert!(f.iter().all(|p| p[1] == 0.5 && p[2] == 0.5));
    }

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
        let x = avg.positions()[0].x;
        assert!(x.abs() < 1e-10 || (x - 10.0).abs() < 1e-10, "got x = {x}, expected 0 or 10");
        assert!((avg.positions()[0].y - 5.0).abs() < 1e-10);
    }

    #[test]
    fn without_a_cell_the_cartesian_mean_is_kept() {
        // 无周期性就没有「跨边界」可言，5.0 是这两个位置唯一讲得通的平均
        let avg = build_avg_frame(&boundary_crosser(false));
        assert!((avg.positions()[0].x - 5.0).abs() < 1e-10);
    }

    #[test]
    fn frames_with_a_different_atom_count_are_skipped() {
        let mut traj = boundary_crosser(true);
        let mut odd = Frame::new();
        odd.add_atom(Atom::new("O", Vector3::new(1.0, 1.0, 1.0)));
        odd.add_atom(Atom::new("O", Vector3::new(2.0, 2.0, 2.0)));
        traj.frames.push(odd);
        let avg = build_avg_frame(&traj);
        assert_eq!(avg.n_atoms(), 1);
        let x = avg.positions()[0].x;
        assert!(x.abs() < 1e-10 || (x - 10.0).abs() < 1e-10, "got x = {x}");
    }



    fn frames_with_steps(steps: &[i64]) -> Trajectory {
        let mut t = Trajectory::new();
        for &s in steps {
            let mut f = Frame::new();
            f.step = Some(s);
            t.add_frame(f);
        }
        t
    }

    #[test]
    fn test_frame_spacing() {
        assert!(check_frame_spacing(&frames_with_steps(&[0, 100, 200, 300])).is_ok());
        assert!(check_frame_spacing(&frames_with_steps(&[0])).is_ok(), "单帧无间隔可查");
        assert!(check_frame_spacing(&Trajectory::from_frame(Frame::new())).is_ok(), "无步号的格式只能信调用方给的 dt");
        let dup = format!("{:#}", check_frame_spacing(&frames_with_steps(&[0, 100, 100, 200])).unwrap_err());
        assert!(dup.contains("frames 1 and 2") && dup.contains("duplicated"), "重启点写两次应点名帧，实际 {dup}");
        let uneven = format!("{:#}", check_frame_spacing(&frames_with_steps(&[0, 100, 200, 400])).unwrap_err());
        assert!(uneven.contains("not evenly spaced") && uneven.contains("frames 2 and 3"), "改过 dump 间隔应点名帧，实际 {uneven}");
    }

    #[test]
    fn bins_round_up_and_end_on_the_upper_bound() {
        // A whole number: no extra bin (the float quotient 10.004 / 0.002 is just above 5002)
        let (n, w) = bins(0.001, 10.005, 0.002);
        assert_eq!(n, 5002);
        assert_eq!(w, 0.002, "a whole number of steps keeps the request exactly");
        // Not a whole number: round up, narrowing the width to end exactly at the bound
        let (n, w) = bins(0.0, 1.0, 0.07);
        assert_eq!(n, 15);
        assert!((w * 15.0 - 1.0).abs() < 1e-15 && w < 0.07);
        // Less than one step still gives one bin; an empty range none
        assert_eq!(bins(0.0, 0.05, 0.07).0, 1);
        assert_eq!(bins(1.0, 1.0, 0.1).0, 0);
    }
}
