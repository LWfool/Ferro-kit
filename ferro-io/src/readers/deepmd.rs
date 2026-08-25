//! DeePMD-kit system directory reader (the `set.*/*.npy` layout).
//!
//! The inverse of [`crate::writers::write_deepmd_npy`], returning a plain
//! [`Trajectory`] like every other reader — the filtering step downstream works
//! on `Trajectory`, so no dataset-specific type crosses a crate boundary.
//!
//! # Two dtypes in the wild
//!
//! dpdata writes `float32` by default, ferro writes `float64`. Both are read;
//! the `float32` path widens on load. A file whose dtype is neither is an
//! error rather than a guess.
//!
//! # Shapes
//!
//! Arrays are stored two-dimensional (`reshape([nframes, -1])`), but nothing
//! forbids a `(nframes,)` energy written by hand, so the loader takes the total
//! element count and the frame count and derives the rest.
//!
//! # Virial
//!
//! `virial.npy` holds eV; [`ferro_core::Frame::stress`] holds eV/Å³. The volume
//! comes from `box.npy` (`|det|`) because a DeePMD system carries no
//! `volume.npy`. Reading a virial without a box is an error — there is nothing
//! to divide by.

use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use ferro_core::{matrix3_from_row_major, Atom, Cell, Frame, Trajectory};
use nalgebra::Vector3;
use ndarray::ArrayD;
use ndarray_npy::read_npy;

/// Keys this reader maps onto `Frame`; anything else in a set is reported.
const KNOWN: &[&str] = &["coord", "box", "energy", "force", "virial"];

/// Reads a DeePMD system directory, discarding the warnings.
pub fn read_deepmd_npy(dir: &Path) -> Result<Trajectory> {
    Ok(read_deepmd_npy_with_warnings(dir)?.0)
}

/// Reads a DeePMD system directory and reports data this reader cannot carry.
///
/// Keys such as `atom_ener.npy` have no home in `Frame`; they are named in the
/// warnings rather than dropped in silence, because a lost label looks exactly
/// like a label that was never computed.
pub fn read_deepmd_npy_with_warnings(dir: &Path) -> Result<(Trajectory, Vec<String>)> {
    let types = read_type_indices(&dir.join("type.raw"))?;
    let type_map = read_lines_tokens(&dir.join("type_map.raw"))?;
    let natoms = types.len();
    if natoms == 0 {
        bail!("{}: type.raw is empty", dir.display());
    }
    for t in &types {
        if *t >= type_map.len() {
            bail!(
                "{}: type.raw holds index {t} but type_map.raw has {} entries",
                dir.display(), type_map.len()
            );
        }
    }
    let elements: Vec<&str> = types.iter().map(|t| type_map[*t].as_str()).collect();

    let sets = set_dirs(dir)?;
    let nopbc = dir.join("nopbc").exists();

    let mut warnings = Vec::new();
    let mut frames: Vec<Frame> = Vec::new();

    for set in &sets {
        let coord = load(set, "coord")?
            .with_context(|| format!("{}: coord.npy is required", set.display()))?;
        if coord.len() % (natoms * 3) != 0 {
            bail!(
                "{}: coord.npy holds {} values, not a multiple of natoms*3 = {}",
                set.display(), coord.len(), natoms * 3
            );
        }
        let nf = coord.len() / (natoms * 3);

        let boxes = load(set, "box")?;
        let energy = load(set, "energy")?;
        let force = load(set, "force")?;
        let virial = load(set, "virial")?;

        for (key, arr, per_frame) in [
            ("box", &boxes, 9usize),
            ("energy", &energy, 1),
            ("force", &force, natoms * 3),
            ("virial", &virial, 9),
        ] {
            if let Some(v) = arr {
                if v.len() != nf * per_frame {
                    bail!(
                        "{}: {key}.npy holds {} values, expected {} for {nf} frame(s)",
                        set.display(), v.len(), nf * per_frame
                    );
                }
            }
        }
        if virial.is_some() && boxes.is_none() {
            bail!("{}: virial.npy without box.npy — no volume to divide by", set.display());
        }

        for name in extra_keys(set)? {
            warnings.push(format!(
                "{}: {name}.npy is not carried by Frame and will be lost if this dataset is written back",
                set.display()
            ));
        }

        for k in 0..nf {
            let mut frame = Frame::new();
            frame.atoms = (0..natoms)
                .map(|i| {
                    let b = (k * natoms + i) * 3;
                    Atom::new(
                        elements[i],
                        Vector3::new(coord[b], coord[b + 1], coord[b + 2]),
                    )
                })
                .collect();

            if let Some(b) = &boxes {
                let nine: [f64; 9] = b[k * 9..k * 9 + 9].try_into().expect("9 values");
                let cell = Cell::from_matrix(matrix3_from_row_major(&nine));
                if let Some(v) = &virial {
                    let vol = cell.volume();
                    if vol <= 0.0 {
                        bail!("{}: frame {k} has a degenerate cell (volume {vol})", set.display());
                    }
                    let nine: [f64; 9] = v[k * 9..k * 9 + 9].try_into().expect("9 values");
                    let mut m = matrix3_from_row_major(&nine);
                    m /= vol;
                    frame.stress = Some(m);
                }
                frame.cell = Some(cell);
                frame.pbc = [!nopbc; 3];
            }
            if let Some(e) = &energy {
                frame.energy = Some(e[k]);
            }
            if let Some(f) = &force {
                frame.forces = Some(
                    (0..natoms)
                        .map(|i| {
                            let b = (k * natoms + i) * 3;
                            Vector3::new(f[b], f[b + 1], f[b + 2])
                        })
                        .collect(),
                );
            }
            frames.push(frame);
        }
    }

    if frames.is_empty() {
        bail!("{}: no frames found under set.*", dir.display());
    }
    let mut traj = Trajectory { frames, metadata: Default::default() };
    traj.metadata.source = Some("DeePMD npy".to_string());
    Ok((traj, warnings))
}

/// `set.*` subdirectories in name order, which is also frame order.
fn set_dirs(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut out: Vec<PathBuf> = std::fs::read_dir(dir)
        .with_context(|| format!("cannot list {}", dir.display()))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            p.is_dir()
                && p.file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.starts_with("set."))
        })
        .collect();
    out.sort();
    if out.is_empty() {
        bail!("{}: no set.* directory (is this a DeePMD system?)", dir.display());
    }
    Ok(out)
}

fn extra_keys(set: &Path) -> Result<Vec<String>> {
    let mut out: Vec<String> = std::fs::read_dir(set)?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter_map(|p| {
            let stem = p.file_stem()?.to_str()?.to_string();
            (p.extension().and_then(|e| e.to_str()) == Some("npy")
                && !KNOWN.contains(&stem.as_str()))
            .then_some(stem)
        })
        .collect();
    out.sort();
    Ok(out)
}

fn load(set: &Path, key: &str) -> Result<Option<Vec<f64>>> {
    let path = set.join(format!("{key}.npy"));
    if !path.exists() {
        return Ok(None);
    }
    // dpdata 默认 float32、ferro 写 float64，两种都收；别的 dtype 报错而不是猜
    match read_npy::<_, ArrayD<f64>>(&path) {
        Ok(a) => Ok(Some(a.iter().copied().collect())),
        Err(e64) => match read_npy::<_, ArrayD<f32>>(&path) {
            Ok(a) => Ok(Some(a.iter().map(|&x| x as f64).collect())),
            Err(_) => Err(anyhow::Error::new(e64))
                .with_context(|| format!("{}: not a float32/float64 .npy", path.display())),
        },
    }
}

fn read_type_indices(path: &Path) -> Result<Vec<usize>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {}", path.display()))?;
    // 一行一个还是一行全部都收：dpdata 的文档说是一行，np.savetxt 写的是一行一个
    text.split_whitespace()
        .map(|t| {
            t.parse::<usize>()
                .with_context(|| format!("{}: `{t}` is not an atom type index", path.display()))
        })
        .collect()
}

fn read_lines_tokens(path: &Path) -> Result<Vec<String>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {}", path.display()))?;
    Ok(text.split_whitespace().map(|s| s.to_string()).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::writers::write_deepmd_npy;
    use ferro_core::Trajectory;
    use nalgebra::Matrix3;
    use ndarray::Array2;
    use ndarray_npy::write_npy;

    fn demo() -> Trajectory {
        let cell = Cell::from_matrix(Matrix3::from_row_slice(&[
            5.0, 0.0, 0.0, 0.1, 6.0, 0.0, 0.2, 0.3, 7.0,
        ]));
        let mut frames = Vec::new();
        for k in 0..3 {
            let mut f = Frame::with_cell(cell.clone(), [true; 3]);
            f.atoms = vec![
                Atom::new("O", Vector3::new(k as f64, 0.5, 0.25)),
                Atom::new("Si", Vector3::new(1.0, 0.0, 0.0)),
            ];
            f.energy = Some(-10.0 - k as f64);
            f.forces = Some(vec![
                Vector3::new(0.5, -0.25, 0.125),
                Vector3::new(-0.5, 0.25, -0.125),
            ]);
            f.stress = Some(Matrix3::from_row_slice(&[
                1.0, 2.0, 3.0, 2.0, 4.0, 5.0, 3.0, 5.0, 6.0,
            ]));
            frames.push(f);
        }
        Trajectory { frames, metadata: Default::default() }
    }

    fn scratch(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(name);
        let _ = std::fs::remove_dir_all(&d);
        d
    }

    #[test]
    fn round_trip_preserves_every_label() {
        let d = scratch("ferro_dp_rt");
        let a = demo();
        write_deepmd_npy(&a, &d).unwrap();
        let (b, warn) = read_deepmd_npy_with_warnings(&d).unwrap();
        assert!(warn.is_empty(), "{warn:?}");
        assert_eq!(b.n_frames(), a.n_frames());
        for (x, y) in a.frames.iter().zip(&b.frames) {
            assert_eq!(x.atoms.len(), y.atoms.len());
            for (p, q) in x.atoms.iter().zip(&y.atoms) {
                assert_eq!(p.element, q.element);
                assert_eq!(p.position, q.position);
            }
            assert_eq!(x.energy, y.energy);
            assert_eq!(x.forces, y.forces);
            assert_eq!(x.cell.as_ref().unwrap().matrix, y.cell.as_ref().unwrap().matrix);
            // stress -> virial -> stress 经过一次乘除体积，不要求逐位相同
            let (s, t) = (x.stress.unwrap(), y.stress.unwrap());
            assert!((s - t).abs().max() < 1e-12, "{s} vs {t}");
        }
    }

    /// dpdata writes float32; the reader must widen rather than refuse.
    #[test]
    fn float32_arrays_are_accepted() {
        let d = scratch("ferro_dp_f32");
        write_deepmd_npy(&demo(), &d).unwrap();
        let set = d.join("set.000");
        // 用 float32 覆盖 coord.npy，模拟 dpdata 的产物
        let a: Array2<f32> = Array2::from_shape_vec(
            (3, 6),
            vec![
                0.0, 0.5, 0.25, 1.0, 0.0, 0.0,
                1.0, 0.5, 0.25, 1.0, 0.0, 0.0,
                2.0, 0.5, 0.25, 1.0, 0.0, 0.0,
            ],
        )
        .unwrap();
        write_npy(set.join("coord.npy"), &a).unwrap();
        let t = read_deepmd_npy(&d).unwrap();
        assert_eq!(t.n_frames(), 3);
        assert_eq!(t.frames[2].atoms[0].position, Vector3::new(2.0, 0.5, 0.25));
    }

    #[test]
    fn unknown_keys_are_reported_not_dropped_silently() {
        let d = scratch("ferro_dp_extra");
        write_deepmd_npy(&demo(), &d).unwrap();
        let a: Array2<f64> = Array2::zeros((3, 2));
        write_npy(d.join("set.000").join("atom_ener.npy"), &a).unwrap();
        let (_, warn) = read_deepmd_npy_with_warnings(&d).unwrap();
        assert_eq!(warn.len(), 1);
        assert!(warn[0].contains("atom_ener"), "{warn:?}");
    }

    #[test]
    fn a_wrong_length_array_is_an_error() {
        let d = scratch("ferro_dp_bad");
        write_deepmd_npy(&demo(), &d).unwrap();
        let a: Array2<f64> = Array2::zeros((2, 9)); // 只有 2 帧，coord 有 3 帧
        write_npy(d.join("set.000").join("box.npy"), &a).unwrap();
        let err = read_deepmd_npy(&d).unwrap_err().to_string();
        assert!(err.contains("box.npy"), "{err}");
    }
}
