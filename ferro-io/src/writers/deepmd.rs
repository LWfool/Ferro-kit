//! DeePMD-kit system directory writer (the `set.*/*.npy` layout).
//!
//! A DeePMD *system* is a directory, not a file:
//!
//! ```text
//! <dir>/
//!   type.raw        one integer per atom, 0-based, indexing type_map.raw
//!   type_map.raw    one element symbol per line
//!   nopbc           empty marker file, when no frame is periodic (no cell, or
//!                   pbc off on all three axes; box.npy is still written, as dpdata does)
//!   set.000/
//!     coord.npy     (nframes, natoms*3)  Angstrom
//!     box.npy       (nframes, 9)         Angstrom, row-major lattice vectors
//!     energy.npy    (nframes, 1)         eV
//!     force.npy     (nframes, natoms*3)  eV/Angstrom
//!     virial.npy    (nframes, 9)         eV
//! ```
//!
//! On disk every array is **two-dimensional** — dpdata flattens with
//! `reshape([nframes, -1])` before saving — so the documented
//! `nframes x natoms x 3` is a logical shape, not the stored one.
//!
//! ferro writes `float64` where dpdata defaults to `float32`. This directory is
//! the head of the pipeline: filtering and merging read it back, and a lossy
//! head cannot be undone downstream. Narrowing to `float32` belongs at the step
//! that feeds the training framework.
//!
//! # Virial
//!
//! `virial = stress * V`, with no sign flip: [`ferro_core::Frame::stress`]
//! keeps the sign CP2K/VASP/QE print (positive = compression), and DeePMD's
//! virial uses that same orientation. In ASE terms both equal `-V * sigma_ase`.
//!
//! # One system, one composition
//!
//! Every frame in a system must share the atom count and the type sequence —
//! `type.raw` is written once for the whole directory. Trajectories that change
//! composition are rejected here rather than producing a directory DeePMD
//! silently misreads.

use std::fs;
use std::path::Path;

use anyhow::{bail, Context, Result};
use ferro_core::data::elements::canonical_cmp;
use ferro_core::{matrix3_row_major, Trajectory};
use ndarray::ArrayView2;
use ndarray_npy::write_npy;

/// Writes `traj` as one DeePMD system directory, all frames in `set.000`.
pub fn write_deepmd_npy(traj: &Trajectory, dir: &Path) -> Result<()> {
    write_deepmd_npy_sets(traj, dir, 0)
}

/// Writes `traj` as one DeePMD system directory, split into sets of `set_size`.
///
/// `set_size` of 0 means one set holding everything. Splitting exists to give a
/// later shuffle its boundaries, so a trailing set of a dozen frames is worse
/// than useless — the remainder is spread over the sets rather than left as a
/// stub.
pub fn write_deepmd_npy_sets(traj: &Trajectory, dir: &Path, set_size: usize) -> Result<()> {
    let bounds = set_bounds(traj.n_frames(), set_size);
    write_deepmd_npy_bounds(traj, dir, &bounds)
}

/// Writes `traj` with the set boundaries given explicitly as `[lo, hi)` ranges.
///
/// Merging by source needs boundaries that fall exactly where one source ends
/// and the next begins, so that every set holds frames from a single condition.
/// An even split cannot express that.
pub fn write_deepmd_npy_bounds(
    traj: &Trajectory,
    dir: &Path,
    bounds: &[(usize, usize)],
) -> Result<()> {
    let frames = &traj.frames;
    if frames.is_empty() {
        bail!("cannot write an empty trajectory as a DeePMD system");
    }
    let natoms = frames[0].n_atoms();
    if natoms == 0 {
        bail!("cannot write a DeePMD system with zero atoms");
    }

    // A system writes type.raw once, so every frame must have the same type sequence
    traj.check_same_atoms().map_err(|e| anyhow::anyhow!(
        "{e}; one DeePMD system holds one composition, so split the trajectory into one system per composition"
    ))?;
    let elements: Vec<&str> = frames[0].atoms().map(|a| a.element).collect();

    // type_map 按 (Z, 符号) 排序 —— 与 gr 的分组排序同一条规则，保证不同
    // system 之间的 type_map 可对齐（merge 依赖这一点）
    let mut species: Vec<&str> = elements.clone();
    species.sort_by(|a, b| canonical_cmp(a, b));
    species.dedup();
    let type_index: Vec<usize> = elements
        .iter()
        .map(|e| species.iter().position(|s| s == e).expect("element in species"))
        .collect();

    let has_cell = frames.iter().all(|f| f.cell.is_some());
    if !has_cell && frames.iter().any(|f| f.cell.is_some()) {
        bail!("some frames carry a cell and some do not; a DeePMD system is either periodic or not");
    }
    // DeePMD only knows periodic along all three axes or not at all (`nopbc`); a slab
    // (`pbc = T T F`) used to be written as fully periodic.
    if let Some(k) = frames.iter().position(|f| f.cell.is_some() && f.pbc != [true; 3] && f.pbc != [false; 3]) {
        bail!("frame {k} is periodic along some axes only (pbc = {:?}); a DeePMD system is \
               periodic along all three or none", frames[k].pbc);
    }
    let nopbc_frames = frames.iter().filter(|f| f.cell.is_none() || f.pbc == [false; 3]).count();
    if nopbc_frames != 0 && nopbc_frames != frames.len() {
        bail!("some frames are periodic and some are not; a DeePMD system is either periodic or not");
    }

    let nf = frames.len();
    let covered: usize = bounds.iter().map(|(lo, hi)| hi.saturating_sub(*lo)).sum();
    if covered != nf {
        bail!("set boundaries cover {covered} frame(s) but the trajectory holds {nf}");
    }

    // 逐项先在整条轨迹上摊平，再按 set 边界切片写出
    let mut arrays: Vec<(String, Vec<f64>, usize)> = Vec::new();

    let mut coord = Vec::with_capacity(nf * natoms * 3);
    for f in frames {
        for a in f.atoms() {
            coord.extend_from_slice(&[a.position.x, a.position.y, a.position.z]);
        }
    }
    arrays.push(("coord".into(), coord, natoms * 3));

    if has_cell {
        let mut boxes = Vec::with_capacity(nf * 9);
        for f in frames {
            boxes.extend_from_slice(&matrix3_row_major(&f.cell.as_ref().unwrap().matrix));
        }
        arrays.push(("box".into(), boxes, 9));
    }

    // 逐项要么全帧都有要么全帧都无：半有半无地补零会把「没算」写成「算出来是 0」
    if let Some(energies) = all_or_none(frames.iter().map(|f| f.energy), "energy")? {
        arrays.push(("energy".into(), energies, 1));
    }
    if let Some(force_frames) = all_or_none(frames.iter().map(|f| f.forces()), "forces")? {
        let mut force = Vec::with_capacity(nf * natoms * 3);
        for fs_ in force_frames {
            if fs_.len() != natoms {
                bail!("a frame carries {} force vectors for {natoms} atoms", fs_.len());
            }
            for v in fs_ {
                force.extend_from_slice(&[v.x, v.y, v.z]);
            }
        }
        arrays.push(("force".into(), force, natoms * 3));
    }
    if all_or_none(frames.iter().map(|f| f.stress), "stress")?.is_some() {
        if !has_cell {
            bail!("stress is present but the frames have no cell; the virial needs a volume");
        }
        let mut virial = Vec::with_capacity(nf * 9);
        for (i, f) in frames.iter().enumerate() {
            let v = f.virial()
                .with_context(|| format!("frame {i}"))?
                .expect("all_or_none: every frame has a stress");
            virial.extend(matrix3_row_major(&v));
        }
        arrays.push(("virial".into(), virial, 9));
    }

    // Every refusal above comes before the directory is touched: clearing the old set.* first used
    // to leave a rejected rewrite with nothing but type.raw
    fs::create_dir_all(dir).with_context(|| format!("cannot create {}", dir.display()))?;
    // A rewrite clears the old set.* first, or fewer frames would leave the previous tail behind
    for entry in fs::read_dir(dir)?.filter_map(|e| e.ok()) {
        let p = entry.path();
        if p.is_dir()
            && p.file_name().and_then(|n| n.to_str()).is_some_and(|n| n.starts_with("set."))
        {
            fs::remove_dir_all(&p)?;
        }
    }

    let type_raw: String = type_index.iter().map(|t| format!("{t}\n")).collect();
    fs::write(dir.join("type.raw"), type_raw)?;
    let map_raw: String = species.iter().map(|s| format!("{s}\n")).collect();
    fs::write(dir.join("type_map.raw"), map_raw)?;

    let nopbc = dir.join("nopbc");
    if nopbc_frames == 0 {
        let _ = fs::remove_file(&nopbc);
    } else {
        fs::write(&nopbc, "")?;
    }

    for (i, &(lo, hi)) in bounds.iter().enumerate() {
        let set_dir = dir.join(format!("set.{i:03}"));
        fs::create_dir_all(&set_dir)
            .with_context(|| format!("cannot create {}", set_dir.display()))?;
        for (name, data, per_frame) in &arrays {
            save2d(
                &set_dir.join(format!("{name}.npy")),
                &data[lo * per_frame..hi * per_frame],
                hi - lo,
                *per_frame,
            )?;
        }
    }

    Ok(())
}

/// `[lo, hi)` frame ranges of at most `set_size` frames, the remainder spread instead
/// of left as a stub set: 500 frames at 400 gives 250 + 250, not 400 + 100 — the
/// lopsided pair is worse for both training balance and for using a set as a
/// validation split. `set_size == 0` means one set.
///
/// `dataset merge` splits each system with it too, so the two paths cannot split
/// sets differently (review A-9: it once had a copy of its own).
pub fn set_bounds(nf: usize, set_size: usize) -> Vec<(usize, usize)> {
    if set_size == 0 || nf <= set_size {
        return vec![(0, nf)];
    }
    let n_sets = nf.div_ceil(set_size);
    let base = nf / n_sets;
    let extra = nf % n_sets;
    let mut out = Vec::with_capacity(n_sets);
    let mut lo = 0;
    for i in 0..n_sets {
        let take = base + usize::from(i < extra);
        out.push((lo, lo + take));
        lo += take;
    }
    out
}

/// `Some(values)` when every frame has the property, `None` when none has.
fn all_or_none<T>(it: impl Iterator<Item = Option<T>>, what: &str) -> Result<Option<Vec<T>>> {
    let mut out = Vec::new();
    let mut missing = 0usize;
    for v in it {
        match v {
            Some(v) => out.push(v),
            None => missing += 1,
        }
    }
    if missing == 0 {
        Ok(Some(out))
    } else if out.is_empty() {
        Ok(None)
    } else {
        bail!("{missing} frame(s) lack {what} while {} have it; a DeePMD set cannot hold a partially labelled property", out.len())
    }
}

// Vec -> 二维视图（零拷贝，from_shape 恒为 C 行序）-> .npy
fn save2d(path: &Path, data: &[f64], rows: usize, cols: usize) -> Result<()> {
    let view = ArrayView2::from_shape((rows, cols), data)
        .with_context(|| format!("{} : {rows}x{cols} does not match {} values", path.display(), data.len()))?;
    write_npy(path, &view).with_context(|| format!("cannot write {}", path.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame};
    use nalgebra::{Matrix3, Vector3};

    pub(super) fn demo_traj() -> Trajectory {
        let cell = Cell::from_matrix(Matrix3::from_row_slice(&[
            5.0, 0.0, 0.0, 0.1, 6.0, 0.0, 0.2, 0.3, 7.0,
        ]));
        let mut frames = Vec::new();
        for k in 0..2 {
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            f.set_atoms(vec![
                Atom::new("O", Vector3::new(k as f64, 0.0, 0.0)),
                Atom::new("Si", Vector3::new(1.0, 0.0, 0.0)),
            ]);
            f.energy = Some(-10.0 - k as f64);
            f.set_forces(Some(vec![Vector3::new(0.5, 0.0, 0.0), Vector3::new(-0.5, 0.0, 0.0)])).unwrap();
            f.stress = Some(Matrix3::from_row_slice(&[
                1.0, 2.0, 3.0, 2.0, 4.0, 5.0, 3.0, 5.0, 6.0,
            ]));
            frames.push(f);
        }
        Trajectory { frames, metadata: Default::default() }
    }

    pub(super) fn out_dir(name: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(name);
        let _ = fs::remove_dir_all(&d);
        d
    }

    #[test]
    fn writes_a_system_directory() {
        let d = out_dir("ferro_dp_basic");
        write_deepmd_npy(&demo_traj(), &d).unwrap();
        // type_map 按 (Z, 符号)：O(8) 在 Si(14) 之前
        assert_eq!(fs::read_to_string(d.join("type_map.raw")).unwrap(), "O\nSi\n");
        assert_eq!(fs::read_to_string(d.join("type.raw")).unwrap(), "0\n1\n");
        assert!(!d.join("nopbc").exists());
        for f in ["coord", "box", "energy", "force", "virial"] {
            assert!(d.join("set.000").join(format!("{f}.npy")).exists(), "{f}.npy missing");
        }
    }

    #[test]
    fn rejects_a_changing_composition() {
        let mut t = demo_traj();
        t.frames[1].sites_mut()[0].element = "N".to_string();
        let err = write_deepmd_npy(&t, &out_dir("ferro_dp_comp")).unwrap_err().to_string();
        assert!(err.contains("atom 0 is N in frame 1") && err.contains("one composition"), "{err}");
    }

    /// The virial is stress times volume; a flat cell has none (review B13: this
    /// writer used to multiply by zero and write it).
    #[test]
    fn rejects_a_stress_in_a_zero_volume_cell() {
        let mut t = demo_traj();
        t.frames[1].cell = Some(Cell::from_matrix(Matrix3::from_row_slice(&[
            5.0, 0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0,
        ])));
        let err = format!("{:#}", write_deepmd_npy(&t, &out_dir("ferro_dp_flat")).unwrap_err());
        assert!(err.contains("frame 1") && err.contains("singular (volume 0"), "{err}");
    }

    #[test]
    fn rejects_a_partially_labelled_property() {
        let mut t = demo_traj();
        t.frames[1].energy = None;
        let err = write_deepmd_npy(&t, &out_dir("ferro_dp_part")).unwrap_err().to_string();
        assert!(err.contains("partially labelled"), "{err}");
    }

    /// A refused rewrite used to clear the old set.* and write type.raw before the checks ran,
    /// leaving a directory with nothing but type.raw; now the old system survives byte for byte.
    #[test]
    fn a_refused_rewrite_leaves_the_old_system_intact() {
        let d = out_dir("ferro_dp_refused");
        write_deepmd_npy(&demo_traj(), &d).unwrap();
        let snapshot = |d: &Path| -> Vec<Vec<u8>> {
            ["type.raw", "type_map.raw", "set.000/coord.npy", "set.000/energy.npy", "set.000/virial.npy"]
                .iter().map(|f| fs::read(d.join(f)).unwrap()).collect()
        };
        let before = snapshot(&d);

        let mut t = demo_traj();
        t.frames[1].energy = None;
        t.frames[0].sites_mut()[0].element = "N".to_string();
        t.frames[1].sites_mut()[0].element = "N".to_string();
        assert!(write_deepmd_npy(&t, &d).is_err());
        assert_eq!(snapshot(&d), before);

        let err = write_deepmd_npy_bounds(&demo_traj(), &d, &[(0, 1)]).unwrap_err().to_string();
        assert!(err.contains("cover 1 frame(s)"), "{err}");
        assert_eq!(snapshot(&d), before);
    }

    #[test]
    fn no_cell_writes_nopbc_and_no_box() {
        let mut t = demo_traj();
        for f in &mut t.frames {
            f.cell = None;
            f.stress = None;
        }
        let d = out_dir("ferro_dp_nopbc");
        write_deepmd_npy(&t, &d).unwrap();
        assert!(d.join("nopbc").exists());
        assert!(!d.join("set.000/box.npy").exists());
    }
}

#[cfg(test)]
mod set_tests {
    use super::tests::{demo_traj, out_dir};
    use super::*;

    #[test]
    fn a_remainder_is_spread_not_left_as_a_stub() {
        // 410 帧、每 set 400：朴素切法给出 400 + 10，那个 10 帧的 set 当验证集没有意义
        assert_eq!(set_bounds(410, 400), vec![(0, 205), (205, 410)]);
        assert_eq!(set_bounds(1000, 400), vec![(0, 334), (334, 667), (667, 1000)]);
        assert_eq!(set_bounds(500, 400), vec![(0, 250), (250, 500)]);
        assert_eq!(
            set_bounds(2000, 400),
            vec![(0, 400), (400, 800), (800, 1200), (1200, 1600), (1600, 2000)]
        );
    }

    #[test]
    fn zero_or_small_input_stays_one_set() {
        assert_eq!(set_bounds(50, 0), vec![(0, 50)]);
        assert_eq!(set_bounds(50, 400), vec![(0, 50)]);
    }

    #[test]
    fn splitting_writes_every_set_and_drops_stale_ones() {
        let d = out_dir("ferro_dp_sets");
        let mut t = demo_traj();
        // demo 只有 2 帧，扩到 5 帧
        while t.frames.len() < 5 {
            t.frames.push(t.frames[0].clone());
        }
        write_deepmd_npy_sets(&t, &d, 2).unwrap();
        for i in 0..3 {
            assert!(d.join(format!("set.{i:03}/coord.npy")).exists(), "set {i}");
        }
        // 再写一次单 set，旧的 set.001/002 必须消失
        write_deepmd_npy(&t, &d).unwrap();
        assert!(d.join("set.000/coord.npy").exists());
        assert!(!d.join("set.001").exists());
    }

    /// A slab used to be written as fully periodic: DeePMD has no per-axis periodicity.
    #[test]
    fn partially_periodic_frames_are_refused() {
        let mut t = demo_traj();
        t.frames[1].pbc = [true, true, false];
        let err = write_deepmd_npy(&t, &out_dir("ferro_dp_slab")).unwrap_err().to_string();
        assert!(err.contains("frame 1") && err.contains("some axes only"), "{err}");
        // Non-periodic along all three axes still writes nopbc
        let mut t = demo_traj();
        for f in &mut t.frames { f.pbc = [false; 3]; }
        let d = out_dir("ferro_dp_nopbc_all");
        write_deepmd_npy(&t, &d).unwrap();
        assert!(d.join("nopbc").exists());
    }
}
