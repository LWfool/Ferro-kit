//! RDF-derived coordination cutoffs: the first minimum of g(r) past its first peak.
//!
//! Shared by `ferro net` (`--P-O=auto`) and `ferro dataset filter --al6`. Both need the
//! same answer for the same trajectory, so there is one definition of "the edge of the
//! first shell" in the project.

use crate::md::gr::{calc_gr, GrParams};

/// The first minimum of g(r) past its first peak, i.e. the edge of the first
/// coordination shell.
///
/// This is how a coordination cutoff is chosen in practice: the first peak is
/// the bonded shell, and the trough behind it is where "bonded" stops meaning
/// anything. Returning it lets a caller pick its own cutoff per system instead of
/// carrying a number copied from another composition.
///
/// **Only meaningful when the pair actually has a coordination shell.** For a
/// bonded pair (Al–O, P–O) the trough is deep and unambiguous; for a pair that
/// does not bond (O–O) the "first minimum" is a shallow feature several Å out
/// and means nothing. The caller is expected to know which case it is in — this
/// function cannot tell them apart, and a shallow trough is reported as
/// [`ShellCutoff::depth`] so the caller can judge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShellCutoff {
    /// Position of the first peak \[Å\]
    pub peak_r: f64,
    pub peak_g: f64,
    /// Position of the first minimum past the peak \[Å\] — the cutoff
    pub min_r: f64,
    /// g(r) at that minimum; near zero means a clean shell, ~1 means no shell
    pub depth: f64,
}

/// Locates the first coordination shell of each pair `(a, b)` in a trajectory.
///
/// One g(r) pass serves every pair: `calc_gr` computes all element pairs anyway, and
/// it is the expensive part. `None` for a pair the trajectory does not contain or whose
/// g(r) never rises above 1.
///
/// Uses a coarse 0.02 Å bin on purpose: the default 0.002 Å of `GrParams` gives
/// a noisy curve whose local minima are sampling artefacts rather than
/// structure.
pub fn first_shell_cutoffs(
    traj: &ferro_core::Trajectory,
    pairs: &[(&str, &str)],
) -> ferro_core::Result<Vec<Option<ShellCutoff>>> {
    let params = GrParams { r_min: 0.001, r_max: 6.0, dr: 0.02, ..Default::default() };
    let res = calc_gr(traj, &params)?;
    Ok(pairs
        .iter()
        .map(|(a, b)| res.gr.get(&format!("{a}-{b}")).and_then(|g| shell_from_curve(&res.r, g)))
        .collect())
}

/// The peak/trough search itself, split out so it can be tested on a synthetic curve.
fn shell_from_curve(r: &[f64], g: &[f64]) -> Option<ShellCutoff> {
    if r.len() != g.len() || r.len() < 3 {
        return None;
    }
    // 第一峰：g 首次越过 1 之后的最大值；越不过 1 说明这对根本没有近邻壳层
    let start = g.iter().position(|&v| v > 1.0)?;
    let mut ip = start;
    for i in start..g.len() {
        if g[i] > g[ip] {
            ip = i;
        }
        // 已经从峰上下来一大截就停，避免把远处更高的峰当第一峰
        if g[i] < g[ip] * 0.5 && i > ip {
            break;
        }
    }
    // 峰后 3 Å 窗口内的最小值
    let dr = r[1] - r[0];
    let window = ((3.0 / dr) as usize).max(1);
    let hi = (ip + window).min(g.len());
    if ip + 1 >= hi {
        return None;
    }
    let mut im = ip + 1;
    for i in ip + 1..hi {
        if g[i] < g[im] {
            im = i;
        }
    }
    // 极小是一段平台（壳层间隙里 g 恒为 0）时，取等于最小值的**最长**连续段的
    // 中点。取首格会把截断贴在第一峰的尾巴上：帧少时峰尾是零星采样，孤立的零格
    // 后面还跟着一个样本（审查 C-D5：tests/cp2k_md_3frames 的 Al-O 在 2.051 为 0、
    // 2.071 为 0.32、2.091 起才是真空隙；取首格得 cn6=46，空隙内任取一点都是 48）
    let (mut best, mut best_len) = (im, 0);
    let mut i = im;
    while i < hi {
        if g[i] != g[im] {
            i += 1;
            continue;
        }
        let s = i;
        while i < hi && g[i] == g[im] {
            i += 1;
        }
        if i - s > best_len {
            (best, best_len) = (s, i - s);
        }
    }
    let mid = best + (best_len - 1) / 2;
    Some(ShellCutoff { peak_r: r[ip], peak_g: g[ip], min_r: r[mid], depth: g[mid] })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A clean shell: peak at 1.8, trough at 2.4, then rising again.
    #[test]
    fn finds_the_trough_behind_the_first_peak() {
        let r: Vec<f64> = (0..300).map(|i| i as f64 * 0.02).collect();
        let g: Vec<f64> = r
            .iter()
            .map(|&x| {
                let peak = 5.0 * (-((x - 1.8) / 0.25).powi(2)).exp();
                let bulk = 1.0 / (1.0 + (-(x - 3.5) * 2.0).exp());
                peak + bulk
            })
            .collect();
        let s = shell_from_curve(&r, &g).unwrap();
        assert!((s.peak_r - 1.8).abs() < 0.05, "peak at {}", s.peak_r);
        assert!((2.3..2.9).contains(&s.min_r), "trough at {}", s.min_r);
        assert!(s.depth < 0.3, "depth {}", s.depth);
    }

    /// A curve that never rises above 1 has no shell to find.
    #[test]
    fn a_curve_without_a_peak_gives_none() {
        let r: Vec<f64> = (0..100).map(|i| i as f64 * 0.02).collect();
        let g = vec![0.4; 100];
        assert!(shell_from_curve(&r, &g).is_none());
    }

    /// 壳层间隙是 g = 0 的平台时，截断取平台中点而不是第一个零格
    #[test]
    fn a_zero_plateau_gives_its_midpoint() {
        // 峰 1.8–2.1，2.2–2.6 恒为 0，2.6 之后回到 1
        let r: Vec<f64> = (0..200).map(|i| i as f64 * 0.02).collect();
        let g: Vec<f64> = r
            .iter()
            .map(|&x| match x {
                x if x < 1.79 => 0.0,
                x if x < 2.19 => 6.0,
                x if x < 2.61 => 0.0,
                _ => 1.0,
            })
            .collect();
        let s = shell_from_curve(&r, &g).unwrap();
        assert!((s.min_r - 2.4).abs() < 0.03, "平台中点应在 2.4 附近，得到 {}", s.min_r);
        assert_eq!(s.depth, 0.0);
    }

    /// 峰尾有孤立零格、后面还有一个样本时，取的是真空隙而不是那个孤立零格
    #[test]
    fn an_isolated_zero_on_the_peak_tail_is_skipped() {
        let r: Vec<f64> = (0..200).map(|i| i as f64 * 0.02).collect();
        let g: Vec<f64> = r
            .iter()
            .map(|&x| match x {
                x if x < 1.79 => 0.0,
                x if x < 2.03 => 6.0,
                x if x < 2.05 => 0.0, // 孤立零格 2.04
                x if x < 2.07 => 0.3, // 峰尾最后一个样本 2.06
                x if x < 2.61 => 0.0,
                _ => 1.0,
            })
            .collect();
        let s = shell_from_curve(&r, &g).unwrap();
        assert!((2.3..2.4).contains(&s.min_r), "应取 2.08–2.60 空隙的中点，得到 {}", s.min_r);
    }
}
