//! `ferro dataset` — build and maintain machine-learning training sets.
//!
//! Three steps, deliberately separate commands rather than one pipeline flag,
//! because the first one is expensive and its output is what gets backed up:
//!
//! - `collect` — AIMD output → DeePMD system directories
//! - `filter`  — quality selection on an existing dataset
//! - `merge`   — combine same-composition datasets, resize sets
//!
//! `collect` writes one system directory per input DIRECTORY: the `.out` files
//! sitting together are the restart segments of one run, so putting them back
//! together is restoring a trajectory, not merging datasets. That is the line
//! between the two commands — `collect` reassembles the pieces of ONE run,
//! `merge` combines DIFFERENT runs of the same composition.

mod collect;
mod filter;
mod inspect;
mod merge;

pub use collect::CollectCmd;
pub use filter::FilterCmd;
pub use merge::MergeCmd;
use collect::run_collect;
use filter::run_filter;
use merge::run_merge;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use clap::{Args, Subcommand};

use crate::batch::expand_inputs;
use ferro_analysis::ml::diagnostics::{filter_diagnostics, quantiles};
use ferro_analysis::ml::merge::{
    composition_key, group_name, shuffle_order, shuffle_within, stitch_segments, DEFAULT_SEED,
};
use ferro_analysis::ml::Split;
use ferro_analysis::first_shell_cutoffs;
use ferro_analysis::ml::{filter_frames, FilterParams, FilterResult};
use ferro_core::units::{convert_pressure, PressureUnit};
use ferro_core::Trajectory;
use ferro_io::{
    read_aimd_as, read_deepmd_npy_with_warnings, set_bounds, write_deepmd_npy,
    write_deepmd_npy_bounds, write_deepmd_npy_sets, write_extxyz_with, AimdFormat, AimdStats,
    StressKey,
};

#[derive(Subcommand, Debug)]
pub enum DatasetCmd {
    /// Extract AIMD output into DeePMD system directories
    Collect(CollectCmd),
    /// Drop low-quality frames from existing datasets
    Filter(FilterCmd),
    /// Combine datasets of the same composition
    Merge(MergeCmd),
}

/// What `collect` produces.
///
/// A separate enum from [`OutType`] on purpose. `collect` cannot write nep/extxyz
/// today, so the two value sets do not overlap beyond `deepmd`; folding `Inspect`
/// into `OutType` would list it in `filter --help` and `merge --help`, where it is
/// rejected at run time — and the help/clap drift test checks option *names*, not
/// value sets, so that dishonesty would go unnoticed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, clap::ValueEnum)]
pub enum CollectType {
    /// DeePMD system directory (`type.raw` + `set.NNN/*.npy`)
    #[default]
    Deepmd,
    /// Diagnostics only: trajectory, last structure and per-frame scalars.
    /// Writes no dataset.
    Inspect,
}

/// What `filter` / `merge` write out.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, clap::ValueEnum)]
pub enum OutType {
    /// DeePMD system directory (`type.raw` + `set.NNN/*.npy`)
    #[default]
    Deepmd,
    /// One extxyz file per system, stress carried as `virial=` (GPUMD / NEP)
    Nep,
    /// One extxyz file per system, stress carried as `stress=` (ASE convention)
    Extxyz,
}

impl OutType {
    /// The stress key, or `None` for the DeePMD path.
    fn stress_key(self) -> Option<StressKey> {
        match self {
            Self::Deepmd => None,
            Self::Nep => Some(StressKey::Virial),
            Self::Extxyz => Some(StressKey::Stress),
        }
    }
}

/// The three split parts, in output order.
const SPLIT_PARTS: [(&str, &str); 3] =
    [("train", ".train"), ("valid", ".valid"), ("test", ".test")];

/// dpgen / dpdata split suffixes, recognised on input and written on output.
const SPLIT_SUFFIXES: [&str; 3] = [".train", ".test", ".valid"];

/// What an unsplit product is called: everything filter and merge write is training
/// data unless a split says otherwise.
const TRAIN_SUFFIX: &str = ".train";

/// What `collect` calls its products — raw, unfiltered data straight off an AIMD run.
///
/// Not a split suffix: it says where the data came from, not which part of a split
/// it is, so it never enters `SPLIT_SUFFIXES`. Putting it there would make the
/// "a `.valid` input must not be split again" guard reject `.db` inputs too.
/// `filter` and `merge` strip it before naming their own products, so a collected
/// `run1.db` filters to `run1.train` rather than `run1.db.train`.
const COLLECT_SUFFIX: &str = ".db";

/// Appends a split suffix to the last component of a path.
fn with_suffix(p: &Path, suffix: &str) -> PathBuf {
    if suffix.is_empty() {
        return p.to_path_buf();
    }
    let name = p.file_name().map(|s| s.to_string_lossy().to_string()).unwrap_or_default();
    p.with_file_name(format!("{name}{suffix}"))
}

/// Refuses to write over an existing dataset unless told to.
fn ensure_writable(dest: &Path, overwrite: bool) -> Result<()> {
    if overwrite {
        return Ok(());
    }
    let occupied = if dest.is_dir() {
        std::fs::read_dir(dest)?.next().is_some()
    } else {
        dest.exists()
    };
    if occupied {
        bail!("{} exists and is not empty (pass --overwrite)", dest.display());
    }
    Ok(())
}

/// Writes one trajectory in the requested form; returns the path written.
///
/// `set_size` only reaches the DeePMD path — an extxyz file has no sets.
fn write_as(
    traj: &Trajectory,
    base: &Path,
    ty: OutType,
    set_size: usize,
    overwrite: bool,
) -> Result<PathBuf> {
    match ty.stress_key() {
        None => {
            ensure_writable(base, overwrite)?;
            write_deepmd_npy_sets(traj, base, set_size)?;
            Ok(base.to_path_buf())
        }
        Some(key) => {
            // 不用 with_extension:base 常带 .train 这类后缀,那会被它当扩展名换掉
            let path = PathBuf::from(format!("{}.xyz", base.display()));
            ensure_writable(&path, overwrite)?;
            if let Some(parent) = path.parent() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("cannot create {}", parent.display()))?;
            }
            write_extxyz_with(traj, &path, key, false)?;
            Ok(path)
        }
    }
}

/// The split suffix `name` already carries, if any.
fn split_suffix_of(name: &str) -> Option<&'static str> {
    SPLIT_SUFFIXES.iter().copied().find(|s| name.ends_with(s))
}

/// Strips a trailing `.train` so the parts of a split can be named afresh.
///
/// Only `.train` is stripped, and only the split path reaches this: `.valid` / `.test`
/// inputs are refused before any file is read, because slicing a training set out of a
/// held-out set is a mistake, not a workflow.
/// Strips a trailing `.db` so a collected system can be renamed by filter / merge.
///
/// Applied before the split logic: `.db` is not one of the split suffixes, so
/// leaving it on would stack a second suffix onto the name and break `merge`'s
/// shared-suffix check, which reads only the trailing segment.
fn without_collect_suffix(base: &Path) -> PathBuf {
    let Some(name) = base.file_name().and_then(|n| n.to_str()) else {
        return base.to_path_buf();
    };
    match name.strip_suffix(COLLECT_SUFFIX) {
        Some(stem) => base.with_file_name(stem),
        None => base.to_path_buf(),
    }
}

fn without_train_suffix(base: &Path) -> PathBuf {
    let Some(name) = base.file_name().and_then(|n| n.to_str()) else {
        return base.to_path_buf();
    };
    match name.strip_suffix(TRAIN_SUFFIX) {
        Some(stem) => base.with_file_name(stem),
        None => base.to_path_buf(),
    }
}

/// Writes every non-empty part of a split and reports each line.
///
/// Unsplit output still gets `.train`: everything `filter` and `merge` write is meant
/// for training unless a split says otherwise, and an unlabelled directory downstream is
/// one nobody can place. A name that already ends in a split suffix keeps it — appending
/// would give `sysA.train.train`.
fn write_split(
    traj: &Trajectory,
    base: &Path,
    split: &Split,
    ty: OutType,
    set_size: usize,
    overwrite: bool,
) -> Result<()> {
    // collect 的 `.db` 先剥掉：它说的是「数据从哪来」，不是划分的哪一部分。
    // 留着就会叠成 run1.db.train，而 merge 的共享后缀检查只看末尾那一段
    let base = &without_collect_suffix(base);

    if split.is_off() {
        let named = match base.file_name().and_then(|n| n.to_str()) {
            Some(name) if split_suffix_of(name).is_some() => base.to_path_buf(),
            _ => with_suffix(base, TRAIN_SUFFIX),
        };
        let dest = write_as(traj, &named, ty, set_size, overwrite)?;
        println!("  -> {}", dest.display());
        return Ok(());
    }

    let stem = without_train_suffix(base);
    let parts = split.parts(traj.n_frames(), &base.display().to_string())?;
    for (idx, (label, suffix)) in SPLIT_PARTS.iter().enumerate() {
        if parts[idx].is_empty() {
            continue;
        }
        let sub = traj.subset(&parts[idx]);
        let dest = write_as(&sub, &with_suffix(&stem, suffix), ty, set_size, overwrite)?;
        println!("  {:<5} {:5} frames -> {}", label, parts[idx].len(), dest.display());
    }
    Ok(())
}

/// True when `ferro dataset collect` was typed with no input.
pub fn wants_help(cmd: &DatasetCmd) -> bool {
    match cmd {
        DatasetCmd::Collect(c) => c.input.is_empty(),
        DatasetCmd::Filter(c) => c.input.is_empty(),
        DatasetCmd::Merge(c) => c.input.is_empty(),
    }
}

pub fn print_help(cmd: &DatasetCmd) {
    match cmd {
        DatasetCmd::Collect(_) => crate::help::print_dataset_collect(),
        DatasetCmd::Filter(_) => crate::help::print_dataset_filter(),
        DatasetCmd::Merge(_) => crate::help::print_dataset_merge(),
    }
}

/// Returns the number of inputs that failed, for the process exit code.
pub fn run(cmd: &DatasetCmd) -> Result<usize> {
    match cmd {
        DatasetCmd::Collect(c) => run_collect(c),
        DatasetCmd::Filter(c) => run_filter(c),
        DatasetCmd::Merge(c) => run_merge(c),
    }
}

/// Directories holding a `type.raw`, searched depth-first and not descended into.
fn find_systems(root: &Path) -> Result<Vec<PathBuf>> {
    if !root.is_dir() {
        bail!("{} is not a directory", root.display());
    }
    if root.join("type.raw").exists() {
        return Ok(vec![root.to_path_buf()]);
    }
    let mut out = Vec::new();
    let mut stack = vec![root.to_path_buf()];
    while let Some(dir) = stack.pop() {
        for entry in std::fs::read_dir(&dir)
            .with_context(|| format!("cannot list {}", dir.display()))?
            .filter_map(|e| e.ok())
        {
            let p = entry.path();
            if !p.is_dir() {
                continue;
            }
            if p.join("type.raw").exists() {
                out.push(p); // 认作 system 就不再往下走
            } else {
                stack.push(p);
            }
        }
    }
    out.sort();
    Ok(out)
}

/// Refuses to split inputs that are already a held-out part.
///
/// `.train` is fine and common — everything unsplit now carries it, and re-splitting a
/// training set is a normal thing to want; `without_train_suffix` keeps the names from
/// growing into `X.train.test`. `.valid` / `.test` are refused: carving a training set
/// out of a held-out set silently destroys what it was held out for.
fn refuse_held_out_inputs(systems: &[PathBuf]) -> Result<()> {
    let held = systems.iter().find_map(|p| {
        let name = p.file_name()?.to_str()?;
        split_suffix_of(name).filter(|s| *s != TRAIN_SUFFIX).map(|s| (p.clone(), s))
    });
    if let Some((path, suffix)) = held {
        bail!(
            "{} is the `{suffix}` part of a split; splitting it again would carve a \
             training set out of held-out data. Split the `.train` part, or the \
             unsplit sources",
            path.display()
        );
    }
    Ok(())
}

// Frames whose energy, cell, stress, coordinates and forces name the frame and atom
// they belong to; shared by the merge and filter tests.
#[cfg(test)]
mod test_ids {
    use super::*;

    // 逐帧、逐原子可反查身份的 system：帧号 g 编进能量、盒子、应力与坐标，原子号 id
    // 编进力的 x 分量。原子按 `order` 给出的 id 顺序写，用来检验重排
    pub(super) const ID_ELEM: [&str; 3] = ["O", "Si", "O"];

    pub(super) fn id_base(id: usize) -> nalgebra::Vector3<f64> {
        nalgebra::Vector3::repeat(1.0 + 2.0 * id as f64)
    }

    pub(super) fn id_cell(g: usize) -> nalgebra::Matrix3<f64> {
        // 三斜且非对称，转置或错帧都会变样
        nalgebra::Matrix3::from_row_slice(&[
            10.0 + 0.01 * g as f64, 0.1, 0.2, 0.3, 11.0, 0.4, 0.5, 0.6, 12.0,
        ])
    }

    pub(super) fn id_stress(g: usize) -> nalgebra::Matrix3<f64> {
        let d = 1e-3 * g as f64;
        nalgebra::Matrix3::from_row_slice(&[d, 1e-3, 2e-3, 3e-3, d, 4e-3, 5e-3, 6e-3, 7e-3])
    }

    pub(super) fn id_system(keys: &[usize], order: &[usize]) -> Trajectory {
        use ferro_core::{Atom, Cell, Frame};
        let frames = keys
            .iter()
            .map(|&g| {
                let mut f = Frame::with_cell(Cell::from_matrix(id_cell(g)), [true; 3]);
                let shift = nalgebra::Vector3::new(0.001 * g as f64, 0.0, 0.0);
                f.set_atoms(order.iter().map(|&id| Atom::new(ID_ELEM[id], id_base(id) + shift)));
                f.set_forces(Some(
                    order.iter().map(|&id| nalgebra::Vector3::new(id as f64, 0.01 * g as f64, 0.5)).collect(),
                )).unwrap();
                f.energy = Some(-(g as f64) - 1.0);
                f.stress = Some(id_stress(g));
                f
            })
            .collect();
        Trajectory { frames, metadata: Default::default() }
    }

    /// 读回 `root` 下全部 system，逐帧逐原子核对身份，返回帧号 g 的出现顺序
    pub(super) fn check_ids(root: &Path) -> Vec<usize> {
        let mut seen = Vec::new();
        for sys in find_systems(root).unwrap() {
            let (traj, _) = read_deepmd_npy_with_warnings(&sys).unwrap();
            for f in &traj.frames {
                let g = (-f.energy.unwrap() - 1.0).round() as usize;
                assert_eq!(f.cell.as_ref().unwrap().matrix, id_cell(g), "帧 {g} 的盒子错配");
                let ds = (f.stress.unwrap() - id_stress(g)).abs().max();
                assert!(ds < 1e-12, "帧 {g} 的应力错配（差 {ds}）");
                let forces = f.forces().unwrap();
                let mut ids: Vec<usize> = Vec::new();
                for (atom, fv) in f.atoms().zip(forces) {
                    let id = fv.x.round() as usize;
                    ids.push(id);
                    assert_eq!(atom.element, ID_ELEM[id], "帧 {g} 原子 {id} 的元素错配");
                    let want = id_base(id) + nalgebra::Vector3::new(0.001 * g as f64, 0.0, 0.0);
                    assert!((atom.position - want).norm() < 1e-9, "帧 {g} 原子 {id} 的坐标错配");
                    assert!((fv.y - 0.01 * g as f64).abs() < 1e-12 && fv.z == 0.5, "帧 {g} 原子 {id} 的力错配");
                }
                ids.sort();
                assert_eq!(ids, vec![0, 1, 2], "帧 {g} 的原子不是原来那三个");
                seen.push(g);
            }
        }
        seen
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 已经划分过的输入:.train 可以再划(剥掉后缀重新命名),.valid / .test 拒绝。
    #[test]
    fn splitting_accepts_train_and_refuses_held_out_parts() {
        refuse_held_out_inputs(&[PathBuf::from("a/x.train"), PathBuf::from("b/y")]).unwrap();
        let err = refuse_held_out_inputs(&[PathBuf::from("a/x.valid")]).unwrap_err().to_string();
        assert!(err.contains("held-out"), "{err}");
        assert_eq!(without_train_suffix(Path::new("o/sysA.train")), PathBuf::from("o/sysA"));
        assert_eq!(without_train_suffix(Path::new("o/sysA")), PathBuf::from("o/sysA"));
    }

    /// 不叠加:名字已经以划分后缀结尾时不再追加。
    #[test]
    fn a_name_that_already_carries_a_suffix_keeps_it() {
        assert_eq!(split_suffix_of("sysA.train"), Some(".train"));
        assert_eq!(split_suffix_of("sysA.test"), Some(".test"));
        assert_eq!(split_suffix_of("sysA"), None);
        assert_eq!(split_suffix_of("init.011"), None);
    }

    #[test]
    fn a_split_suffix_lands_on_the_last_component_only() {
        // sys.train 这种名字里的点不是扩展名,不能被 with_extension 吃掉
        assert_eq!(with_suffix(Path::new("a/b/sysA"), ".train"),
                   PathBuf::from("a/b/sysA.train"));
        assert_eq!(with_suffix(Path::new("a/b/sysA"), ""), PathBuf::from("a/b/sysA"));
    }

    #[test]
    fn the_stress_key_follows_the_output_type() {
        assert_eq!(OutType::Deepmd.stress_key(), None);
        assert_eq!(OutType::Nep.stress_key(), Some(StressKey::Virial));
        assert_eq!(OutType::Extxyz.stress_key(), Some(StressKey::Stress));
    }
}
