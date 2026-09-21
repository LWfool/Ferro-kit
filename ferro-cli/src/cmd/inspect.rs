//! Diagnostic products for one collected AIMD run.
//!
//! `ferro dataset collect --type inspect` writes three files into a
//! `ferro_inspect/` folder **inside** the AIMD directory, so each run keeps its
//! own diagnostics next to the output that produced them — CP2K logs are very
//! often all called `total.out` and are told apart only by their directory.
//!
//! The three are a trajectory to look at, the last structure to carry on from,
//! and one table of the per-frame scalars that say whether the run was healthy.
//! Nothing here is part of a dataset: no npy is written on this path.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use ferro_core::table::Table;
use ferro_core::units::AMU_ANG3_TO_G_CM3;
use ferro_core::Trajectory;
use ferro_io::readers::lammps_dump::LammpsUnits;
use ferro_io::writers::table::{write_table, TableFormat};
use ferro_io::writers::{lammps_data::write_lammps_data, lammps_dump::write_lammps_dump};

/// The folder the three products go into, inside the AIMD directory.
///
/// Prefixed because it is created in the user's own data directory rather than
/// under `-o`: a bare `inspect/` gives nobody a way to tell, three months later,
/// what put it there or whether it is safe to delete.
pub const DIR_NAME: &str = "ferro_inspect";

/// Writes the trajectory, the last structure and the per-frame table.
///
/// `sources` names the file each frame came from, one entry per frame; `collect`
/// stitches restart segments, and that column plus the step numbers are what make
/// a seam — or a doubled run — visible.
pub fn write_all(
    traj: &Trajectory,
    dir: &Path,
    name: &str,
    sources: &[String],
    overwrite: bool,
) -> Result<Vec<PathBuf>> {
    if !overwrite && dir.exists() && std::fs::read_dir(dir)?.next().is_some() {
        anyhow::bail!("{} exists and is not empty (pass --overwrite)", dir.display());
    }
    std::fs::create_dir_all(dir)
        .with_context(|| format!("cannot create {}", dir.display()))?;

    let trj = dir.join(format!("{name}.lammpstrj"));
    write_lammps_dump(traj, &trj, LammpsUnits::Metal)?;

    // 末帧而不是第 0 帧：第 0 帧是喂进 CP2K 的初始构型，用户手上本来就有；
    // 末帧才是这次 AIMD 产出的新东西，接着跑或做下一步 DFT 用的是它
    let data = dir.join(format!("{name}.data"));
    let last = Trajectory::from_frame(
        traj.last().context("the trajectory has no frame")?.clone(),
    );
    write_lammps_data(&last, &data)?;

    let info = dir.join(format!("{name}_info.csv"));
    write_table(&build_table(traj, sources), &info, TableFormat::Csv)?;

    Ok(vec![trj, data, info])
}

/// One row per frame, with the run summary carried in the `#` metadata block.
///
/// The summary is an aggregate of the very columns below it, so it lives in the
/// comment block rather than in a second file: `pandas.read_csv(..., comment="#")`
/// gets the table, and a human reading the head of the file gets the summary.
fn build_table(traj: &Trajectory, sources: &[String]) -> Table {
    let n = traj.n_frames();
    let mut frame_idx = Vec::with_capacity(n);
    let mut steps = Vec::with_capacity(n);
    let mut temps = Vec::with_capacity(n);
    let mut energies = Vec::with_capacity(n);
    let mut volumes = Vec::with_capacity(n);
    let mut densities = Vec::with_capacity(n);

    for (i, f) in traj.frames.iter().enumerate() {
        frame_idx.push(i as f64);
        // 缺失一律 NaN，渲染成空字段。补 0 会把「这条路读不出步号」说成「第 0 步」
        steps.push(f.step.map_or(f64::NAN, |s| s as f64));
        temps.push(f.temperature.unwrap_or(f64::NAN));
        energies.push(f.energy.unwrap_or(f64::NAN));
        let (v, d) = match &f.cell {
            // 体积走 |det|，不是三个对角线相乘 —— 后者只对正交胞成立，
            // 三斜胞下静默给出错的密度
            Some(c) => {
                let v = c.volume();
                (v, f.total_mass() / v * AMU_ANG3_TO_G_CM3)
            }
            None => (f64::NAN, f64::NAN),
        };
        volumes.push(v);
        densities.push(d);
    }

    let mut t = Table::new();
    for line in summary_lines(traj, &temps, &densities) {
        t.meta_line(line);
    }
    t.push_num("frame", frame_idx);
    t.push_num("step", steps);
    t.push_num("temperature", temps);
    t.push_num("energy", energies);
    t.push_num("volume", volumes);
    t.push_num("density", densities);
    t.push_text("source", sources.to_vec());
    t
}

/// The `#` block: what this run is, and the two quantities worth a glance.
fn summary_lines(traj: &Trajectory, temps: &[f64], densities: &[f64]) -> Vec<String> {
    let mut out = Vec::new();
    let Some(first) = traj.first() else { return out };

    let mut composition: Vec<String> = Vec::new();
    for e in first.unique_elements() {
        composition.push(format!("{e}{}", first.count_element(&e)));
    }
    out.push(format!("atoms = {}  ({})", first.n_atoms(), composition.join(" ")));
    out.push(format!("frames = {}", traj.n_frames()));

    if let Some(c) = &first.cell {
        let [a, b, cc] = c.lengths();
        let [al, be, ga] = c.angles();
        out.push(format!(
            "cell (frame 0) = a {a:.4}  b {b:.4}  c {cc:.4} Ang   \
             alpha {al:.2}  beta {be:.2}  gamma {ga:.2} deg"
        ));
    }
    if let Some(l) = mean_sd_line("temperature", temps, "K") {
        out.push(l);
    }
    if let Some(l) = mean_sd_line("density", densities, "g/cm^3") {
        out.push(l);
    }
    // vasprun 没有温度字段，这一列是反算来的；不说清楚会让同一次运行换个来源
    // 时数字对不上却无人解释得了
    if traj.frames.iter().any(|f| f.temperature.is_some())
        && traj.metadata.source.as_deref() == Some("vasprun.xml")
    {
        out.push(
            "note: temperature is derived from the ionic kinetic energy \
             (T = 2*E_kin/(3N*k_B)); vasprun.xml does not print it"
                .to_string(),
        );
    }
    out
}

/// `mean +/- sd` over the finite values, or `None` when there are none.
///
/// Two-pass, as `ferro-analysis` does elsewhere: `<x^2> - <x>^2` loses most of its
/// significant digits when the spread is small next to the mean, which is exactly
/// the case for a thermostatted temperature.
fn mean_sd_line(label: &str, values: &[f64], unit: &str) -> Option<String> {
    let xs: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    if xs.is_empty() {
        return None;
    }
    let mean = xs.iter().sum::<f64>() / xs.len() as f64;
    let sd = if xs.len() > 1 {
        (xs.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (xs.len() - 1) as f64).sqrt()
    } else {
        0.0
    };
    Some(format!("{label} = {mean:.4} +/- {sd:.4} {unit}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::{Atom, Cell, Frame};
    use nalgebra::Vector3;

    fn traj_of(temps: &[Option<f64>], steps: &[Option<i64>]) -> Trajectory {
        let mut t = Trajectory::new();
        for (i, (temp, step)) in temps.iter().zip(steps).enumerate() {
            let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
            let mut f = Frame::with_cell(cell, [true; 3]);
            f.add_atom(Atom::new("Si", Vector3::new(0.0, 0.0, 0.0)));
            f.add_atom(Atom::new("O", Vector3::new(1.0, 1.0, 1.0)));
            f.energy = Some(-10.0 - i as f64);
            f.temperature = *temp;
            f.step = *step;
            t.add_frame(f);
        }
        t
    }

    /// Missing scalars render as empty fields, never as zero.
    #[test]
    fn a_missing_scalar_is_empty_not_zero() {
        let traj = traj_of(&[Some(300.0), None], &[Some(1), None]);
        let table = build_table(&traj, &["a.out".into(), "a.out".into()]);
        let step = table.column("step").unwrap();
        assert_eq!(step.cell(0), "1.000000e0");
        assert_eq!(step.cell(1), "", "缺步号必须是空字段,不是 0");
        let temp = table.column("temperature").unwrap();
        assert_eq!(temp.cell(1), "", "缺温度必须是空字段,不是 0");
    }

    /// The summary skips a quantity none of the frames carry.
    #[test]
    fn the_summary_omits_a_quantity_that_is_absent_everywhere() {
        let traj = traj_of(&[None, None], &[Some(1), Some(2)]);
        let table = build_table(&traj, &["a.out".into(), "a.out".into()]);
        let meta = table.meta.join("\n");
        assert!(meta.contains("atoms = 2"), "{meta}");
        assert!(meta.contains("frames = 2"), "{meta}");
        assert!(!meta.contains("temperature ="), "无温度时不该报一个假的均值: {meta}");
        assert!(meta.contains("density ="), "{meta}");
    }

    #[test]
    fn the_sd_is_zero_for_a_constant_series_and_not_noise() {
        let line = mean_sd_line("temperature", &[300.0, 300.0, 300.0], "K").unwrap();
        assert_eq!(line, "temperature = 300.0000 +/- 0.0000 K");
        assert!(mean_sd_line("x", &[f64::NAN], "K").is_none());
    }
}
