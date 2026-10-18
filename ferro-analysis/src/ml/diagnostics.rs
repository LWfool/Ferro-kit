//! Read-only diagnostics for choosing filter thresholds.
//!
//! These answer the question the criteria themselves cannot: *should* you be
//! filtering at all, and at what value. They exist because both geometric
//! criteria are far more sensitive to their parameter than to the data — an Al6
//! selection can go from flagging every frame to flagging none over 0.3 Å of
//! cutoff, and a minimum O–O distance is an extreme-value statistic whose
//! location shifts with the number of O atoms in the cell, so a threshold
//! carried over from another system means nothing.

use std::collections::BTreeMap;

use crate::network::TypeParams;
use ferro_core::units::AMU_ANG3_TO_G_CM3;
use ferro_core::{Cell, Table, Trajectory};

use super::filter::{FilterParams, FilterResult};

use super::geometry::coordination_histogram;
use crate::check;

/// How a coordination selection responds to its cutoff.
#[derive(Debug, Clone)]
pub struct CutoffScan {
    pub rcut: f64,
    /// Fraction of frames holding at least one atom with the target coordination
    pub frame_fraction: f64,
    /// Mean number of such atoms per frame
    pub per_frame: f64,
}

/// Scans `rcuts` and reports how many frames an "at least one" selection keeps.
///
/// Sampling at most `max_frames` frames: the shape of this curve is what
/// matters, and it does not need every frame to be visible. A selection sitting
/// on a steep part of the curve is chosen by its cutoff, not by the physics.
pub fn cutoff_scan(
    traj: &Trajectory,
    elem: &str,
    ligand: &str,
    target_cn: u32,
    rcuts: &[f64],
    max_frames: usize,
) -> ferro_core::Result<Vec<CutoffScan>> {
    // The scan reaches past the selection's own cutoff, which is all the entry of
    // the filter checked (review D-A5)
    let widest = rcuts.iter().copied().fold(0.0, f64::max);
    check::periodic_cells(traj, Some((&format!("{elem}-{ligand} cutoff scan"), widest)))?;
    let n = traj.n_frames();
    let step = (n / max_frames.max(1)).max(1);
    let sampled: Vec<usize> = (0..n).step_by(step).collect();

    rcuts
        .iter()
        .map(|&rcut| {
            let mut cut = BTreeMap::new();
            cut.insert((elem.to_string(), ligand.to_string()), rcut);
            let params = TypeParams::new(cut, Default::default());
            let mut hit = 0usize;
            let mut total = 0usize;
            for &i in &sampled {
                let f = &traj.frames[i];
                if f.cell.is_none() { continue }
                let n_target: usize = coordination_histogram(f, &params, elem)?
                    .into_iter()
                    .filter(|(cn, _)| *cn == target_cn)
                    .map(|(_, c)| c)
                    .sum();
                if n_target > 0 {
                    hit += 1;
                }
                total += n_target;
            }
            let m = sampled.len().max(1) as f64;
            Ok(CutoffScan {
                rcut,
                frame_fraction: hit as f64 / m,
                per_frame: total as f64 / m,
            })
        })
        .collect()
}

/// Pooled coordination-number histogram for one element over the trajectory.
pub fn pooled_coordination(
    traj: &Trajectory,
    params: &TypeParams,
    elem: &str,
) -> ferro_core::Result<Vec<(u32, usize)>> {
    let mut hist: BTreeMap<u32, usize> = BTreeMap::new();
    for f in &traj.frames {
        if f.cell.is_none() { continue }
        for (cn, n) in coordination_histogram(f, params, elem)? {
            *hist.entry(cn).or_default() += n;
        }
    }
    Ok(hist.into_iter().collect())
}

/// `[min, p1, p50, p99, max, mean]` over the finite values; `None` when there are none.
///
/// Reported instead of a single "how many are below the threshold" count
/// because that count cannot distinguish an outlier tail from a smooth
/// distribution — and only the first is worth filtering.
pub fn quantiles(values: &[f64]) -> Option<[f64; 6]> {
    let mut v: Vec<f64> = values.iter().copied().filter(|x| x.is_finite()).collect();
    // 先滤后判空：全是非有限值时下面的 `len() - 1` 会下溢
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |p: f64| -> f64 {
        let i = ((v.len() - 1) as f64 * p).round() as usize;
        v[i]
    };
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    Some([v[0], q(0.01), q(0.5), q(0.99), v[v.len() - 1], mean])
}

/// A coarse histogram of a per-frame quantity, over its finite values.
///
/// Carries no header line: the quantiles are per input, so they belong in the
/// `[inputs]` block (see [`quantiles`]) — a header line would survive stacking
/// only for the first input and pass for the whole batch.
pub fn distribution_table(values: &[f64], bins: usize) -> Table {
    let mut t = Table::new();
    let Some([lo, _, _, _, hi, _]) = quantiles(values) else { return t };
    let width = ((hi - lo) / bins as f64).max(f64::MIN_POSITIVE);
    let mut counts = vec![0usize; bins];
    for &x in values.iter().filter(|x| x.is_finite()) {
        let k = (((x - lo) / width) as usize).min(bins - 1);
        counts[k] += 1;
    }
    t.push_text(
        "bin_lo",
        (0..bins).map(|i| format!("{:.3}", lo + i as f64 * width)).collect(),
    )
    .push_text(
        "bin_hi",
        (0..bins).map(|i| format!("{:.3}", lo + (i + 1) as f64 * width)).collect(),
    )
    .push_text("count", counts.iter().map(|c| c.to_string()).collect());
    t
}

/// `value -> frame count` for an integer per-frame quantity.
pub fn count_histogram(name: &str, values: &[usize]) -> Table {
    let mut hist: BTreeMap<usize, usize> = BTreeMap::new();
    for &v in values {
        *hist.entry(v).or_default() += 1;
    }
    let mut t = Table::new();
    t.push_text(name.to_string(), hist.keys().map(|k| k.to_string()).collect())
        .push_text("frames", hist.values().map(|v| v.to_string()).collect());
    t
}

/// Builds the cutoff-scan table.
pub fn scan_table(scan: &[CutoffScan]) -> Table {
    let mut t = Table::new();
    t.meta_line("a steep column here means the selection is decided by the cutoff, not the data");
    t.push_text("rcut", scan.iter().map(|s| format!("{:.2}", s.rcut)).collect())
        .push_text(
            "frames_pct",
            scan.iter().map(|s| format!("{:.1}", 100.0 * s.frame_fraction)).collect(),
        )
        .push_text(
            "per_frame",
            scan.iter().map(|s| format!("{:.2}", s.per_frame)).collect(),
        );
    t
}

/// Coordination histogram as a table, with percentages.
pub fn coordination_table(hist: &[(u32, usize)]) -> Table {
    let total: usize = hist.iter().map(|(_, c)| *c).sum();
    let mut t = Table::new();
    t.push_text("cn", hist.iter().map(|(cn, _)| cn.to_string()).collect())
        .push_text("atoms", hist.iter().map(|(_, c)| c.to_string()).collect())
        .push_text(
            "pct",
            hist.iter()
                .map(|(_, c)| format!("{:.2}", 100.0 * *c as f64 / total.max(1) as f64))
                .collect(),
        );
    t
}

/// A [`Cell`]-carrying frame is required by every function here; this is the
/// check callers make before asking.
pub fn all_frames_have_cells(traj: &Trajectory) -> bool {
    traj.frames.iter().all(|f| f.cell.is_some())
}

/// Convenience for callers holding a cell they already validated.
pub fn cell_of(traj: &Trajectory, i: usize) -> Option<&Cell> {
    traj.frames.get(i)?.cell.as_ref()
}

/// The tables that answer "should I be filtering, and at what value".
///
/// A selection whose outcome swings with its cutoff is chosen by the cutoff
/// rather than by the structure, and a minimum distance drawn from a smooth
/// distribution has no outliers to remove — neither is visible from the funnel
/// alone. Only the criteria that are switched on contribute a table.
pub fn filter_diagnostics(
    traj: &Trajectory,
    r: &FilterResult,
    params: &FilterParams,
) -> ferro_core::Result<Vec<(String, Table)>> {
    let mut out = Vec::new();
    if params.oo_min > 0.0 {
        let v: Vec<f64> = r.verdicts.iter().filter_map(|x| x.min_oo).collect();
        out.push(("min_oo".to_string(), distribution_table(&v, 16)));
    }

    let Some(rcut) = params.al6_rcut else { return Ok(out) };

    let n6: Vec<usize> = r.verdicts.iter().filter_map(|x| x.n_al6).collect();
    if !n6.is_empty() {
        out.push(("al6".to_string(), count_histogram("n_al6", &n6)));
    }

    let mut cut = BTreeMap::new();
    cut.insert(("Al".to_string(), "O".to_string()), rcut);
    let tp = TypeParams::new(cut, Default::default());
    let hist = pooled_coordination(traj, &tp, "Al")?;
    if !hist.is_empty() {
        // The cutoff differs from system to system (--al6 without a value), so it
        // is recorded in the al6_rcut column of [inputs]
        out.push(("al_cn".to_string(), coordination_table(&hist)));
    }

    // A sweep around the current cutoff: how steep it is is what this table is for
    let rcuts: Vec<f64> = (-3..=3).map(|k| rcut + k as f64 * 0.1).filter(|v| *v > 0.0).collect();
    let scan = cutoff_scan(traj, "Al", "O", 6, &rcuts, 200)?;
    out.push(("rcut_scan".to_string(), scan_table(&scan)));
    Ok(out)
}

/// One row per frame of an AIMD run — step, temperature, energy, volume, density
/// and the file the frame came from — with the run summary in the `#` metadata block.
///
/// The summary is an aggregate of the very columns below it, so it lives in the
/// comment block rather than in a second file: `pandas.read_csv(..., comment="#")`
/// gets the table, and a human reading the head of the file gets the summary.
/// `sources` names the file of each frame, one entry per frame.
pub fn frame_table(traj: &Trajectory, sources: &[String]) -> Table {
    let n = traj.n_frames();
    let mut frame_idx = Vec::with_capacity(n);
    let mut steps = Vec::with_capacity(n);
    let mut temps = Vec::with_capacity(n);
    let mut energies = Vec::with_capacity(n);
    let mut volumes = Vec::with_capacity(n);
    let mut densities = Vec::with_capacity(n);

    for (i, f) in traj.frames.iter().enumerate() {
        frame_idx.push(i as f64);
        // A missing value is NaN, rendered as an empty field. A 0 would say "step 0"
        // where the truth is "this reader has no step number"
        steps.push(f.step.map_or(f64::NAN, |s| s as f64));
        temps.push(f.temperature.unwrap_or(f64::NAN));
        energies.push(f.energy.unwrap_or(f64::NAN));
        let (v, d) = match &f.cell {
            // The volume is |det|, not the product of the diagonal: that holds for
            // an orthogonal cell only and gives a silently wrong density otherwise
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

    fn run_of(temps: &[Option<f64>], steps: &[Option<i64>]) -> Trajectory {
        use ferro_core::{Atom, Frame};
        use nalgebra::Vector3;
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
        let traj = run_of(&[Some(300.0), None], &[Some(1), None]);
        let table = frame_table(&traj, &["a.out".into(), "a.out".into()]);
        let step = table.column("step").unwrap();
        assert_eq!(step.cell(0), "1.000000e0");
        assert_eq!(step.cell(1), "", "a missing step must be an empty field, not 0");
        let temp = table.column("temperature").unwrap();
        assert_eq!(temp.cell(1), "", "a missing temperature must be an empty field, not 0");
    }

    /// The summary skips a quantity none of the frames carry.
    #[test]
    fn the_summary_omits_a_quantity_that_is_absent_everywhere() {
        let traj = run_of(&[None, None], &[Some(1), Some(2)]);
        let table = frame_table(&traj, &["a.out".into(), "a.out".into()]);
        let meta = table.meta.join("\n");
        assert!(meta.contains("atoms = 2"), "{meta}");
        assert!(meta.contains("frames = 2"), "{meta}");
        assert!(!meta.contains("temperature ="), "no temperature, so no made-up mean: {meta}");
        assert!(meta.contains("density ="), "{meta}");
    }

    #[test]
    fn the_sd_is_zero_for_a_constant_series_and_not_noise() {
        let line = mean_sd_line("temperature", &[300.0, 300.0, 300.0], "K").unwrap();
        assert_eq!(line, "temperature = 300.0000 +/- 0.0000 K");
        assert!(mean_sd_line("x", &[f64::NAN], "K").is_none());
    }

    #[test]
    fn distribution_of_only_non_finite_values_is_empty_not_a_panic() {
        // 过滤掉非有限值后为空时，`len() - 1` 曾下溢
        let t = distribution_table(&[f64::NAN, f64::INFINITY], 4);
        assert_eq!(t.n_rows(), 0);
        assert!(quantiles(&[f64::NAN]).is_none());
    }

    #[test]
    fn distribution_reports_quantiles_and_bins() {
        let v: Vec<f64> = (0..100).map(|i| i as f64 * 0.01).collect();
        let q = quantiles(&v).unwrap();
        assert!((q[2] - 0.5).abs() < 1e-12, "p50 应为 0.5：{q:?}");
        let t = distribution_table(&v, 4);
        // 分位数逐输入，进 [inputs]；表头留空，否则堆叠后只剩第一个输入的
        assert!(t.meta.is_empty(), "{:?}", t.meta);
        assert_eq!(t.cols.len(), 3);
    }

    #[test]
    fn empty_input_gives_an_empty_table() {
        assert_eq!(distribution_table(&[], 4).cols.len(), 0);
    }

    #[test]
    fn count_histogram_groups_by_value() {
        let t = count_histogram("n", &[0, 0, 1, 2, 2, 2]);
        assert_eq!(t.cols.len(), 2);
    }
}
