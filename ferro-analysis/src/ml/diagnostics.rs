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

use ferro_core::{Cell, Table, Trajectory, TypeParams};

use super::geometry::coordination_histogram;

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
) -> Vec<CutoffScan> {
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
                let Some(cell) = f.cell.as_ref() else { continue };
                let n_target: usize = coordination_histogram(f, cell, &params, elem)
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
            CutoffScan {
                rcut,
                frame_fraction: hit as f64 / m,
                per_frame: total as f64 / m,
            }
        })
        .collect()
}

/// Pooled coordination-number histogram for one element over the trajectory.
pub fn pooled_coordination(
    traj: &Trajectory,
    params: &TypeParams,
    elem: &str,
) -> Vec<(u32, usize)> {
    let mut hist: BTreeMap<u32, usize> = BTreeMap::new();
    for f in &traj.frames {
        let Some(cell) = f.cell.as_ref() else { continue };
        for (cn, n) in coordination_histogram(f, cell, params, elem) {
            *hist.entry(cn).or_default() += n;
        }
    }
    hist.into_iter().collect()
}

/// Quantiles and a coarse histogram of a per-frame quantity.
///
/// Reported instead of a single "how many are below the threshold" count
/// because that count cannot distinguish an outlier tail from a smooth
/// distribution — and only the first is worth filtering.
pub fn distribution_table(name: &str, values: &[f64], bins: usize) -> Table {
    let mut t = Table::new();
    if values.is_empty() {
        return t;
    }
    let mut v: Vec<f64> = values.iter().copied().filter(|x| x.is_finite()).collect();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |p: f64| -> f64 {
        let i = ((v.len() - 1) as f64 * p).round() as usize;
        v[i]
    };
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    t.meta_line(format!(
        "{name}: min {:.3}  p1 {:.3}  p50 {:.3}  p99 {:.3}  max {:.3}  mean {:.3}",
        v[0], q(0.01), q(0.5), q(0.99), v[v.len() - 1], mean
    ));

    let (lo, hi) = (v[0], v[v.len() - 1]);
    let width = ((hi - lo) / bins as f64).max(f64::MIN_POSITIVE);
    let mut counts = vec![0usize; bins];
    for &x in &v {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distribution_reports_quantiles_and_bins() {
        let v: Vec<f64> = (0..100).map(|i| i as f64 * 0.01).collect();
        let t = distribution_table("d", &v, 4);
        assert!(t.meta.iter().any(|m| m.contains("p50 0.500")), "{:?}", t.meta);
        assert_eq!(t.cols.len(), 3);
    }

    #[test]
    fn empty_input_gives_an_empty_table() {
        assert_eq!(distribution_table("d", &[], 4).cols.len(), 0);
    }

    #[test]
    fn count_histogram_groups_by_value() {
        let t = count_histogram("n", &[0, 0, 1, 2, 2, 2]);
        assert_eq!(t.cols.len(), 2);
    }
}
