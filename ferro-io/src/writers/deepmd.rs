//! DeePMD-kit system directory writer (the `set.*/*.npy` layout).
//!
//! A DeePMD *system* is a directory, not a file:
//!
//! ```text
//! <dir>/
//!   type.raw        one integer per atom, 0-based, indexing type_map.raw
//!   type_map.raw    one element symbol per line
//!   nopbc           empty marker file, only when the frames carry no cell
//!   set.000/
//!     coord.npy     (nframes, natoms*3)  Angstrom
//!     box.npy       (nframes, 9)         Angstrom, row-major lattice vectors
//!     energy.npy    (nframes,)           eV
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
use ferro_core::data::elements::symbol_to_z;
use ferro_core::{matrix3_row_major, Trajectory};
use ndarray::ArrayView2;
use ndarray_npy::write_npy;

/// Writes `traj` as one DeePMD system directory.
pub fn write_deepmd_npy(traj: &Trajectory, dir: &Path) -> Result<()> {
    let frames = &traj.frames;
    if frames.is_empty() {
        bail!("cannot write an empty trajectory as a DeePMD system");
    }
    let natoms = frames[0].atoms.len();
    if natoms == 0 {
        bail!("cannot write a DeePMD system with zero atoms");
    }

    // 一个 system 内 type.raw 只写一次，故各帧的类型序列必须逐项相同
    let elements: Vec<&str> = frames[0].atoms.iter().map(|a| a.element.as_str()).collect();
    for (i, f) in frames.iter().enumerate() {
        if f.atoms.len() != natoms {
            bail!(
                "frame {i} has {} atoms but frame 0 has {natoms}; split the trajectory into one system per composition",
                f.atoms.len()
            );
        }
        if f.atoms.iter().zip(&elements).any(|(a, e)| a.element != *e) {
            bail!("frame {i} has a different type sequence than frame 0; one DeePMD system holds one composition");
        }
    }

    // type_map 按 (Z, 符号) 排序 —— 与 gr 的分组排序同一条规则，保证不同
    // system 之间的 type_map 可对齐（merge 依赖这一点）
    let mut species: Vec<&str> = elements.clone();
    species.sort_by_key(|s| (symbol_to_z(s), *s));
    species.dedup();
    let type_index: Vec<usize> = elements
        .iter()
        .map(|e| species.iter().position(|s| s == e).expect("element in species"))
        .collect();

    let has_cell = frames.iter().all(|f| f.cell.is_some());
    if !has_cell && frames.iter().any(|f| f.cell.is_some()) {
        bail!("some frames carry a cell and some do not; a DeePMD system is either periodic or not");
    }

    fs::create_dir_all(dir).with_context(|| format!("cannot create {}", dir.display()))?;
    let set_dir = dir.join("set.000");
    fs::create_dir_all(&set_dir).with_context(|| format!("cannot create {}", set_dir.display()))?;

    let type_raw: String = type_index.iter().map(|t| format!("{t}\n")).collect();
    fs::write(dir.join("type.raw"), type_raw)?;
    let map_raw: String = species.iter().map(|s| format!("{s}\n")).collect();
    fs::write(dir.join("type_map.raw"), map_raw)?;

    let nopbc = dir.join("nopbc");
    if has_cell {
        let _ = fs::remove_file(&nopbc);
    } else {
        fs::write(&nopbc, "")?;
    }

    let nf = frames.len();

    // coord: (nframes, natoms*3)
    let mut coord = Vec::with_capacity(nf * natoms * 3);
    for f in frames {
        for a in &f.atoms {
            coord.extend_from_slice(&[a.position.x, a.position.y, a.position.z]);
        }
    }
    save2d(&set_dir.join("coord.npy"), &coord, nf, natoms * 3)?;

    if has_cell {
        let mut boxes = Vec::with_capacity(nf * 9);
        for f in frames {
            boxes.extend_from_slice(&matrix3_row_major(&f.cell.as_ref().unwrap().matrix));
        }
        save2d(&set_dir.join("box.npy"), &boxes, nf, 9)?;
    }

    // 逐项要么全帧都有要么全帧都无：半有半无地补零会把「没算」写成「算出来是 0」
    if let Some(energies) = all_or_none(frames.iter().map(|f| f.energy), "energy")? {
        save2d(&set_dir.join("energy.npy"), &energies, nf, 1)?;
    }
    if let Some(force_frames) = all_or_none(
        frames.iter().map(|f| f.forces.as_ref()),
        "forces",
    )? {
        let mut force = Vec::with_capacity(nf * natoms * 3);
        for fs_ in force_frames {
            if fs_.len() != natoms {
                bail!("a frame carries {} force vectors for {natoms} atoms", fs_.len());
            }
            for v in fs_ {
                force.extend_from_slice(&[v.x, v.y, v.z]);
            }
        }
        save2d(&set_dir.join("force.npy"), &force, nf, natoms * 3)?;
    }
    if let Some(stresses) = all_or_none(frames.iter().map(|f| f.stress), "stress")? {
        if !has_cell {
            bail!("stress is present but the frames have no cell; the virial needs a volume");
        }
        let mut virial = Vec::with_capacity(nf * 9);
        for (f, s) in frames.iter().zip(stresses) {
            let v = f.cell.as_ref().unwrap().volume();
            for x in matrix3_row_major(&s) {
                virial.push(x * v);
            }
        }
        save2d(&set_dir.join("virial.npy"), &virial, nf, 9)?;
    }

    Ok(())
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

    fn demo_traj() -> Trajectory {
        let cell = Cell::from_matrix(Matrix3::from_row_slice(&[
            5.0, 0.0, 0.0, 0.1, 6.0, 0.0, 0.2, 0.3, 7.0,
        ]));
        let mut frames = Vec::new();
        for k in 0..2 {
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            f.atoms = vec![
                Atom::new("O", Vector3::new(k as f64, 0.0, 0.0)),
                Atom::new("Si", Vector3::new(1.0, 0.0, 0.0)),
            ];
            f.energy = Some(-10.0 - k as f64);
            f.forces = Some(vec![Vector3::new(0.5, 0.0, 0.0), Vector3::new(-0.5, 0.0, 0.0)]);
            f.stress = Some(Matrix3::from_row_slice(&[
                1.0, 2.0, 3.0, 2.0, 4.0, 5.0, 3.0, 5.0, 6.0,
            ]));
            frames.push(f);
        }
        Trajectory { frames, metadata: Default::default() }
    }

    fn out_dir(name: &str) -> std::path::PathBuf {
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
        t.frames[1].atoms[0].element = "N".to_string();
        assert!(write_deepmd_npy(&t, &out_dir("ferro_dp_comp")).is_err());
    }

    #[test]
    fn rejects_a_partially_labelled_property() {
        let mut t = demo_traj();
        t.frames[1].energy = None;
        let err = write_deepmd_npy(&t, &out_dir("ferro_dp_part")).unwrap_err().to_string();
        assert!(err.contains("partially labelled"), "{err}");
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
