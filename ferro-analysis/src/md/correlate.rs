//! FFT autocorrelation shared by the time-correlation analyses.
//!
//! `Σ_t x[t]·x[t+m]` over every time origin in `O(N log N)` (Wiener–Khinchin):
//! zero-pad to at least `2N` so the circular correlation of the FFT equals the
//! linear one, transform, take `|X|²`, transform back.
//!
//! Only `msd` uses it so far; `vacf` and `rotcorr` are to follow (see `dev/plan.md`).

use std::sync::Arc;

use rustfft::num_complex::Complex;
use rustfft::{Fft, FftPlanner};

/// A planned forward/inverse transform pair for series of one length.
///
/// Planning is the expensive part and the plans are `Send + Sync`, so plan once
/// with [`AutocorrPlan::new`] and give each worker thread its own
/// [`Autocorr`] (scratch buffers) via [`AutocorrPlan::worker`].
#[derive(Clone)]
pub struct AutocorrPlan {
    n: usize,
    len: usize,
    fwd: Arc<dyn Fft<f64>>,
    inv: Arc<dyn Fft<f64>>,
}

impl AutocorrPlan {
    /// Plans transforms for series of length `n`.
    pub fn new(n: usize) -> Self {
        // ≥ 2N 才能让循环相关不回卷；取 2 的幂让 rustfft 走最快的基数
        let len = (2 * n).next_power_of_two();
        let mut planner = FftPlanner::new();
        AutocorrPlan {
            n,
            len,
            fwd: planner.plan_fft_forward(len),
            inv: planner.plan_fft_inverse(len),
        }
    }

    /// Scratch buffers for one thread.
    pub fn worker(&self) -> Autocorr {
        let scratch = self.fwd.get_inplace_scratch_len().max(self.inv.get_inplace_scratch_len());
        Autocorr {
            plan: self.clone(),
            buf: vec![Complex::new(0.0, 0.0); self.len],
            scratch: vec![Complex::new(0.0, 0.0); scratch],
        }
    }
}

/// One thread's view of an [`AutocorrPlan`].
pub struct Autocorr {
    plan: AutocorrPlan,
    buf: Vec<Complex<f64>>,
    scratch: Vec<Complex<f64>>,
}

impl Autocorr {
    /// Writes `out[m] = Σ_{t=0}^{N-1-m} x[t]·x[t+m]` for `m < out.len()` —
    /// the **sum** over origins, not the mean; divide by `N − m` for that.
    pub fn sums(&mut self, x: &[f64], out: &mut [f64]) {
        let p = &self.plan;
        assert_eq!(x.len(), p.n, "series length differs from the planned length");
        assert!(out.len() <= p.n, "lag beyond the series length");
        for (b, &v) in self.buf.iter_mut().zip(x) {
            *b = Complex::new(v, 0.0);
        }
        for b in &mut self.buf[p.n..] {
            *b = Complex::new(0.0, 0.0);
        }
        p.fwd.process_with_scratch(&mut self.buf, &mut self.scratch);
        for b in &mut self.buf {
            *b = Complex::new(b.norm_sqr(), 0.0);
        }
        p.inv.process_with_scratch(&mut self.buf, &mut self.scratch);
        // rustfft 不归一化：正反变换一轮放大 len 倍
        let inv_len = 1.0 / p.len as f64;
        for (o, b) in out.iter_mut().zip(&self.buf) {
            *o = b.re * inv_len;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sums_match_direct_correlation() {
        // 非 2 的幂长度 + 非平凡数据，与 O(N²) 直接求和逐点比
        let n = 37;
        let x: Vec<f64> = (0..n).map(|i| ((i * 7 % 11) as f64 - 5.0) * 0.3 + (i as f64).sin()).collect();
        let mut out = vec![0.0; n];
        AutocorrPlan::new(n).worker().sums(&x, &mut out);
        for m in 0..n {
            let direct: f64 = (0..n - m).map(|t| x[t] * x[t + m]).sum();
            assert!((out[m] - direct).abs() < 1e-10, "lag {m}: fft {} vs direct {direct}", out[m]);
        }
    }

    #[test]
    fn test_worker_reusable_across_series() {
        // 同一 worker 连续处理两条序列，缓冲区残留不得串到第二条
        let plan = AutocorrPlan::new(4);
        let mut w = plan.worker();
        let mut a = vec![0.0; 4];
        w.sums(&[9.0, 9.0, 9.0, 9.0], &mut a);
        let mut b = vec![0.0; 4];
        w.sums(&[1.0, 0.0, 0.0, 0.0], &mut b);
        assert!((b[0] - 1.0).abs() < 1e-12 && b[1..].iter().all(|v| v.abs() < 1e-12), "{b:?}");
    }
}
