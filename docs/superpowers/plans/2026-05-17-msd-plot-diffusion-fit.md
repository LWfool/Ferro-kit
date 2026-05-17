# MSD plot_msd + Self-Diffusion Coefficient Fit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `plot_msd`, wire `--plot` into `fe-traj -m msd`, and add `--fit-range FMIN,FMAX` that fits the linear segment of the MSD curve and reports the self-diffusion coefficient.

**Architecture:** Pure fit logic (`fit_diffusion`, `MsdFit`) lives in `ferro-analysis/src/md/msd.rs`; `MsdResult` carries `Option<MsdFit>`. `ferro-cli/src/plot.rs` renders it; `ferro-cli/src/bin/traj.rs` wires the CLI flag; help + CLAUDE.md document it. Fraction range maps onto the lag-time axis of the MSD curve.

**Tech Stack:** Rust, cargo workspace, `nalgebra`, `rayon`, `plotters`, `clap`, `anyhow`, `thiserror`.

**Spec:** `docs/superpowers/specs/2026-05-17-msd-plot-diffusion-fit-design.md`

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `ferro-analysis/src/md/msd.rs` | MSD calc + new pure diffusion fit | Add `fit_range` to `MsdParams`, `MsdFit` struct, `fit` to `MsdResult`, `fit_diffusion`, header block in `write_msd`, tests |
| `ferro-cli/src/plot.rs` | CLI PNG rendering | New `plot_msd`; import `MsdResult` |
| `ferro-cli/src/bin/traj.rs` | CLI arg wiring | `--fit-range` arg; `run_msd` computes/prints D, plots |
| `ferro-cli/src/help.rs` | Mode help text | `print_msd()` adds `--fit-range`, `--plot`, example |
| `CLAUDE.md` | Project docs | msd flag block + plot-scope note |

Conversions used everywhere: `D[cm²/s] = D[Å²/fs]·0.1`, `D[m²/s] = D[Å²/fs]·1e-5`.

---

### Task 1: Plumb `fit_range` / `MsdFit` / `MsdResult.fit` (compile-green, no behavior)

**Files:**
- Modify: `ferro-analysis/src/md/msd.rs`

- [ ] **Step 1: Add `fit_range` to `MsdParams` and its `Default`**

In `ferro-analysis/src/md/msd.rs`, change the `MsdParams` struct (currently lines ~22-32) to:

```rust
/// Parameters for MSD calculation.
#[derive(Debug, Clone)]
pub struct MsdParams {
    /// Lag window size in frames (`None` = use all frames)
    pub tau: Option<usize>,
    /// Time shift between origins in frames (default: 1)
    pub shift: usize,
    /// Time step per frame \[fs\] (default: 1.0)
    pub dt: f64,
    /// Elements to include (`None` = all atoms)
    pub elements: Option<Vec<String>>,
    /// Linear-fit window as fractions of the MSD lag-time axis
    /// (`(fmin, fmax)`, `0 <= fmin < fmax <= 1`). `None` = no fit.
    pub fit_range: Option<(f64, f64)>,
}

impl Default for MsdParams {
    fn default() -> Self {
        MsdParams { tau: None, shift: 1, dt: 1.0, elements: None, fit_range: None }
    }
}
```

- [ ] **Step 2: Add `MsdFit` struct and `fit` field on `MsdResult`**

Immediately after the `MsdResult` struct definition block (after the closing `}` of `MsdResult`, ~line 66), add:

```rust
/// Linear-fit result for self-diffusion coefficient extraction.
///
/// `D = slope / 6` (Einstein relation, 3-D isotropic). Unit conversions:
/// `D[cm²/s] = d_ang2_per_fs · 0.1`, `D[m²/s] = d_ang2_per_fs · 1e-5`.
#[derive(Debug, Clone)]
pub struct MsdFit {
    /// Lower fraction of the lag-time axis used for the fit
    pub frac_lo: f64,
    /// Upper fraction of the lag-time axis used for the fit
    pub frac_hi: f64,
    /// First time point of the fit window \[fs\]
    pub t_lo: f64,
    /// Last time point of the fit window \[fs\]
    pub t_hi: f64,
    /// Fitted slope of MSD vs time \[Å²/fs\]
    pub slope: f64,
    /// Fitted intercept \[Å²\]
    pub intercept: f64,
    /// Self-diffusion coefficient `slope / 6` \[Å²/fs\]
    pub d_ang2_per_fs: f64,
    /// Coefficient of determination of the linear fit
    pub r2: f64,
    /// Number of points used in the fit
    pub n_points: usize,
}
```

Then add `fit` as the last field of `MsdResult` (inside its struct, after `pub elements: Vec<String>,`):

```rust
    /// Linear-fit / self-diffusion result (`None` unless `fit_range` was set)
    pub fit: Option<MsdFit>,
```

- [ ] **Step 3: Set `fit: None` in `build_result`**

In `build_result` (~lines 308-323), change the returned struct literal so it ends with `fit: None`:

```rust
    MsdResult {
        time, msd, msd_a, msd_b, msd_c,
        n_atoms, n_origins, params: params.clone(), elements,
        fit: None,
    }
```

(Leave the `build_result` signature `-> MsdResult` unchanged for now.)

- [ ] **Step 4: Add `fit_range: None` to the three explicit `MsdParams { .. }` test constructions**

In the `#[cfg(test)] mod tests`, three tests build `MsdParams` with explicit fields. Update each:

`test_msd_linear_motion` (~line 416):
```rust
        let result = calc_msd(&traj, &MsdParams {
            tau: Some(n), shift: 1, dt: 1.0, elements: None, fit_range: None,
        }).unwrap();
```

`test_msd_unwrap_across_boundary` (~line 451):
```rust
        let result = calc_msd(&traj, &MsdParams {
            tau: Some(n), shift: 1, dt: 1.0, elements: None, fit_range: None,
        }).unwrap();
```

`test_msd_element_filter` (~line 478):
```rust
        let result = calc_msd(&traj, &MsdParams {
            tau: Some(5), shift: 1, dt: 1.0,
            elements: Some(vec!["Li".to_string()]), fit_range: None,
        }).unwrap();
```

`test_msd_time_shift_averaging` (~line 497):
```rust
        let result = calc_msd(&traj, &MsdParams {
            tau: Some(3), shift: 1, dt: 2.0, elements: None, fit_range: None,
        }).unwrap();
```

(Tests using `MsdParams::default()` need no change — the new field defaults to `None`.)

- [ ] **Step 5: Build and run existing MSD tests**

Run: `cargo test --package ferro-analysis msd`
Expected: PASS (all existing MSD tests green; crate compiles with new fields).

- [ ] **Step 6: Commit**

```bash
git add ferro-analysis/src/md/msd.rs
git commit -m "Plumb MsdFit/fit_range/MsdResult.fit (no behavior yet)

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 2: TDD `fit_diffusion` (pure least-squares + validation)

**Files:**
- Modify: `ferro-analysis/src/md/msd.rs`

- [ ] **Step 1: Write the failing tests**

Inside `#[cfg(test)] mod tests`, add these three tests (place after `test_write_msd`):

```rust
    #[test]
    fn test_fit_diffusion_exact_line() {
        // msd = 6*D*t + c with known D → recovered slope/6 == D, R² == 1
        let d_true = 1.5e-4_f64; // Å²/fs
        let c = 0.7_f64;
        let dt = 2.0_f64;
        let n = 500;
        let time: Vec<f64> = (0..n).map(|i| i as f64 * dt).collect();
        let msd: Vec<f64> = time.iter().map(|&t| 6.0 * d_true * t + c).collect();

        let fit = fit_diffusion(&time, &msd, (0.2, 0.9)).unwrap();
        assert!((fit.d_ang2_per_fs - d_true).abs() < 1e-12,
            "D: expected {d_true}, got {}", fit.d_ang2_per_fs);
        assert!((fit.slope - 6.0 * d_true).abs() < 1e-12);
        assert!((fit.intercept - c).abs() < 1e-9);
        assert!((fit.r2 - 1.0).abs() < 1e-12, "R² = {}", fit.r2);
    }

    #[test]
    fn test_fit_range_index_mapping() {
        // n=11, dt=10 → time 0..100; frac (0.3,0.8) → indices 3..=8
        let dt = 10.0;
        let n = 11;
        let time: Vec<f64> = (0..n).map(|i| i as f64 * dt).collect();
        let msd: Vec<f64> = time.clone(); // slope 1
        let fit = fit_diffusion(&time, &msd, (0.3, 0.8)).unwrap();
        assert_eq!(fit.n_points, 6);
        assert!((fit.t_lo - 30.0).abs() < 1e-12);
        assert!((fit.t_hi - 80.0).abs() < 1e-12);
        assert!((fit.slope - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_fit_range_invalid() {
        let time: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let msd = time.clone();
        assert!(fit_diffusion(&time, &msd, (0.8, 0.3)).is_err()); // fmin>=fmax
        assert!(fit_diffusion(&time, &msd, (-0.1, 0.5)).is_err()); // out of range
        assert!(fit_diffusion(&time, &msd, (0.2, 1.5)).is_err()); // out of range
        let t2 = vec![0.0, 1.0];
        let m2 = vec![0.0, 1.0];
        assert!(fit_diffusion(&t2, &m2, (0.0, 0.001)).is_err()); // <2 points
    }
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --package ferro-analysis msd::tests::test_fit_diffusion_exact_line`
Expected: FAIL — compile error `cannot find function fit_diffusion in this scope`.

- [ ] **Step 3: Implement `fit_diffusion`**

In `ferro-analysis/src/md/msd.rs`, add this function in the `// ─── 计算 ───` region (e.g. immediately before `build_result`):

```rust
/// Ordinary least-squares fit of total MSD vs time over a fractional window
/// of the lag-time axis. `frac = (fmin, fmax)` with `0 <= fmin < fmax <= 1`
/// mapped to indices `i_lo = round(fmin·(n-1))`, `i_hi = round(fmax·(n-1))`.
///
/// Returns slope/intercept, `D = slope / 6` (Einstein, 3-D isotropic) and the
/// fit `R²`. Errors on invalid range, length mismatch, or a window with
/// fewer than 2 points / zero x-variance.
pub fn fit_diffusion(
    time: &[f64],
    msd: &[f64],
    frac: (f64, f64),
) -> ferro_core::Result<MsdFit> {
    let (fmin, fmax) = frac;
    if !(0.0..=1.0).contains(&fmin) || !(0.0..=1.0).contains(&fmax) || fmin >= fmax {
        return Err(ChemError::ValidationError(format!(
            "fit-range must satisfy 0 <= fmin < fmax <= 1, got [{fmin}, {fmax}]"
        )));
    }
    let n = time.len();
    if n != msd.len() {
        return Err(ChemError::ValidationError(
            "fit_diffusion: time/msd length mismatch".into(),
        ));
    }
    if n < 2 {
        return Err(ChemError::ValidationError(
            "MSD curve has fewer than 2 points; cannot fit".into(),
        ));
    }
    let last = (n - 1) as f64;
    let i_lo = (fmin * last).round() as usize;
    let i_hi = ((fmax * last).round() as usize).min(n - 1);
    if i_hi <= i_lo || (i_hi - i_lo + 1) < 2 {
        return Err(ChemError::ValidationError(format!(
            "fit-range [{fmin}, {fmax}] selects fewer than 2 points (i_lo={i_lo}, i_hi={i_hi})"
        )));
    }

    let xs = &time[i_lo..=i_hi];
    let ys = &msd[i_lo..=i_hi];
    let m = xs.len() as f64;
    let sx: f64 = xs.iter().sum();
    let sy: f64 = ys.iter().sum();
    let sxx: f64 = xs.iter().map(|x| x * x).sum();
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| x * y).sum();
    let denom = m * sxx - sx * sx;
    if denom.abs() < f64::EPSILON {
        return Err(ChemError::ValidationError(
            "degenerate fit window (zero x-variance)".into(),
        ));
    }
    let slope = (m * sxy - sx * sy) / denom;
    let intercept = (sy - slope * sx) / m;

    let mean_y = sy / m;
    let ss_tot: f64 = ys.iter().map(|y| (y - mean_y).powi(2)).sum();
    let ss_res: f64 = xs
        .iter()
        .zip(ys)
        .map(|(x, y)| (y - (slope * x + intercept)).powi(2))
        .sum();
    let r2 = if ss_tot.abs() < f64::EPSILON {
        1.0
    } else {
        1.0 - ss_res / ss_tot
    };

    Ok(MsdFit {
        frac_lo: fmin,
        frac_hi: fmax,
        t_lo: xs[0],
        t_hi: xs[xs.len() - 1],
        slope,
        intercept,
        d_ang2_per_fs: slope / 6.0,
        r2,
        n_points: xs.len(),
    })
}
```

- [ ] **Step 4: Run the new tests to verify they pass**

Run: `cargo test --package ferro-analysis msd::tests::test_fit_diffusion_exact_line msd::tests::test_fit_range_index_mapping msd::tests::test_fit_range_invalid`
Expected: PASS (3 passed).

- [ ] **Step 5: Run the full MSD test module**

Run: `cargo test --package ferro-analysis msd`
Expected: PASS (all green).

- [ ] **Step 6: Commit**

```bash
git add ferro-analysis/src/md/msd.rs
git commit -m "Add fit_diffusion: least-squares self-diffusion coefficient

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 3: Wire `fit_diffusion` into `calc_msd` / `build_result`

**Files:**
- Modify: `ferro-analysis/src/md/msd.rs`

- [ ] **Step 1: Write the failing test**

Add to `#[cfg(test)] mod tests`:

```rust
    #[test]
    fn test_calc_msd_populates_fit() {
        let traj = make_traj_static(10);

        let none = calc_msd(&traj, &MsdParams::default()).unwrap();
        assert!(none.fit.is_none());

        let with = calc_msd(&traj, &MsdParams {
            fit_range: Some((0.0, 1.0)),
            ..MsdParams::default()
        }).unwrap();
        assert!(with.fit.is_some());
        let f = with.fit.unwrap();
        assert!((f.d_ang2_per_fs).abs() < 1e-12); // static traj → D = 0
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --package ferro-analysis msd::tests::test_calc_msd_populates_fit`
Expected: FAIL — `assert!(with.fit.is_some())` fails (still `None`, `build_result` hard-codes `fit: None`).

- [ ] **Step 3: Add early fit-range validation in `calc_msd`**

In `calc_msd`, right after the `if n_steps < 2 { ... }` block (~line 125), insert:

```rust
    // Fail fast on an obviously bad fit-range before the heavy parallel loop.
    if let Some((fmin, fmax)) = params.fit_range {
        if !(0.0..=1.0).contains(&fmin) || !(0.0..=1.0).contains(&fmax) || fmin >= fmax {
            return Err(ChemError::ValidationError(format!(
                "fit-range must satisfy 0 <= fmin < fmax <= 1, got [{fmin}, {fmax}]"
            )));
        }
    }
```

- [ ] **Step 4: Change `build_result` to compute the fit and return `Result`**

Replace the entire `build_result` function with:

```rust
/// Build an `MsdResult` from the parallel-reduction accumulation array,
/// computing the diffusion fit when `params.fit_range` is set.
fn build_result(
    accum: Vec<[f64; 4]>,
    tau: usize,
    n_origins: usize,
    n_atoms: usize,
    elements: Vec<String>,
    params: &MsdParams,
) -> ferro_core::Result<MsdResult> {
    let inv = 1.0 / n_origins as f64;
    let time:  Vec<f64> = (0..tau).map(|i| i as f64 * params.dt).collect();
    let msd:   Vec<f64> = (0..tau).map(|i| accum[i][0] * inv).collect();
    let msd_a: Vec<f64> = (0..tau).map(|i| accum[i][1] * inv).collect();
    let msd_b: Vec<f64> = (0..tau).map(|i| accum[i][2] * inv).collect();
    let msd_c: Vec<f64> = (0..tau).map(|i| accum[i][3] * inv).collect();
    let fit = match params.fit_range {
        Some(fr) => Some(fit_diffusion(&time, &msd, fr)?),
        None => None,
    };
    Ok(MsdResult {
        time, msd, msd_a, msd_b, msd_c,
        n_atoms, n_origins, params: params.clone(), elements, fit,
    })
}
```

- [ ] **Step 5: Update the two `build_result` call sites**

In `calc_msd_periodic`, the final line is `Ok(build_result(accum, tau, n_origins, n_atoms, elements, params))`. Change it to:

```rust
    build_result(accum, tau, n_origins, n_atoms, elements, params)
```

In `calc_msd_nonperiodic`, the final line is `Ok(build_result(accum, tau, n_origins, n_atoms, elements, params))`. Change it to:

```rust
    build_result(accum, tau, n_origins, n_atoms, elements, params)
```

- [ ] **Step 6: Run the test and the module**

Run: `cargo test --package ferro-analysis msd`
Expected: PASS (`test_calc_msd_populates_fit` and all existing MSD tests green).

- [ ] **Step 7: Commit**

```bash
git add ferro-analysis/src/md/msd.rs
git commit -m "Compute diffusion fit in calc_msd when fit_range is set

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 4: `write_msd` header block for the fit

**Files:**
- Modify: `ferro-analysis/src/md/msd.rs`

- [ ] **Step 1: Write the failing test**

Add to `#[cfg(test)] mod tests`:

```rust
    #[test]
    fn test_write_msd_with_fit() {
        use std::io::Read;
        let traj = make_traj_linear(10.0, 0.2, 6);
        let result = calc_msd(&traj, &MsdParams {
            fit_range: Some((0.0, 1.0)),
            ..MsdParams::default()
        }).unwrap();
        let path = "/tmp/test_ferro_fit.msd";
        write_msd(&result, path).expect("write_msd failed");

        let mut content = String::new();
        std::fs::File::open(path).unwrap().read_to_string(&mut content).unwrap();
        assert!(content.contains("# D (total)"), "missing D line:\n{content}");
        assert!(content.contains("cm^2/s"), "missing cm^2/s conversion");
        assert!(content.contains("# R^2"), "missing R^2 line");
        assert!(content.contains("# time[fs]"), "column header lost");
    }
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test --package ferro-analysis msd::tests::test_write_msd_with_fit`
Expected: FAIL — `missing D line` (header not written yet).

- [ ] **Step 3: Add the fit header block to `write_msd`**

In `write_msd`, find the closing separator line written just before the column header:

```rust
    writeln!(w, "# {}", "-".repeat(60))?;
    writeln!(w, "# time[fs]\tmsd[Ang^2]\tmsd_a[Ang^2]\tmsd_b[Ang^2]\tmsd_c[Ang^2]")?;
```

Insert the fit block *between* those two lines so it becomes:

```rust
    writeln!(w, "# {}", "-".repeat(60))?;
    if let Some(f) = &result.fit {
        writeln!(w, "# fit range  = [{:.2}, {:.2}]  ->  t in [{:.1}, {:.1}] fs",
            f.frac_lo, f.frac_hi, f.t_lo, f.t_hi)?;
        writeln!(w, "# points     = {}", f.n_points)?;
        writeln!(w, "# slope      = {:.6e} Ang^2/fs", f.slope)?;
        writeln!(w, "# D (total)  = {:.6e} Ang^2/fs", f.d_ang2_per_fs)?;
        writeln!(w, "#            = {:.6e} cm^2/s = {:.6e} m^2/s",
            f.d_ang2_per_fs * 0.1, f.d_ang2_per_fs * 1e-5)?;
        writeln!(w, "# R^2        = {:.6}", f.r2)?;
        writeln!(w, "# {}", "-".repeat(60))?;
    }
    writeln!(w, "# time[fs]\tmsd[Ang^2]\tmsd_a[Ang^2]\tmsd_b[Ang^2]\tmsd_c[Ang^2]")?;
```

- [ ] **Step 4: Run tests**

Run: `cargo test --package ferro-analysis msd`
Expected: PASS (`test_write_msd_with_fit` plus all existing MSD tests green).

- [ ] **Step 5: Commit**

```bash
git add ferro-analysis/src/md/msd.rs
git commit -m "write_msd: emit fit/diffusion header block when fit present

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 5: `plot_msd` in `ferro-cli/src/plot.rs`

**Files:**
- Modify: `ferro-cli/src/plot.rs`

- [ ] **Step 1: Import `MsdResult`**

Change the import on line 5 from:

```rust
use ferro_analysis::{AngleResult, GrResult, SqResult};
```

to:

```rust
use ferro_analysis::{AngleResult, GrResult, MsdResult, SqResult};
```

- [ ] **Step 2: Add `plot_msd`**

Append this function at the end of `ferro-cli/src/plot.rs`, immediately *before* the `// ─── 系统查看器 ───` section (i.e. before `pub fn open_plot`):

```rust
// ─── MSD ─────────────────────────────────────────────────────────────────────

/// Plot total MSD plus a/b/c components vs time; overlay the linear fit and
/// D / R² when `result.fit` is present. Saved as PNG.
pub fn plot_msd(result: &MsdResult, dat_path: &str) -> Result<String> {
    let out = png_path(dat_path);

    let x_max = result.time.last().copied().unwrap_or(1.0).max(1e-9);
    let y_max = result.msd.iter()
        .chain(result.msd_a.iter())
        .chain(result.msd_b.iter())
        .chain(result.msd_c.iter())
        .copied()
        .fold(0.0f64, f64::max) * 1.1;
    let y_max = y_max.max(1e-6);

    {
        let root = BitMapBackend::new(&out, (900, 540)).into_drawing_area();
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .caption("Mean Squared Displacement", ("sans-serif", 18))
            .margin(15)
            .x_label_area_size(45)
            .y_label_area_size(70)
            .build_cartesian_2d(0.0..x_max, 0.0..y_max)?;

        chart.configure_mesh()
            .x_desc("t  [fs]")
            .y_desc("MSD  [Å²]")
            .x_labels(10).y_labels(8)
            .draw()?;

        // a/b/c 分量：淡色细线
        for (k, (label, data)) in [
            ("MSD_a", &result.msd_a),
            ("MSD_b", &result.msd_b),
            ("MSD_c", &result.msd_c),
        ].into_iter().enumerate() {
            let c = color(k + 1);
            let pts: Vec<(f64, f64)> =
                result.time.iter().copied().zip(data.iter().copied()).collect();
            let s = ShapeStyle { color: c.mix(0.55), filled: false, stroke_width: 1 };
            chart.draw_series(LineSeries::new(pts, s))?
                .label(label)
                .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], s));
        }

        // total MSD：主曲线
        {
            let c = color(0);
            let pts: Vec<(f64, f64)> =
                result.time.iter().copied().zip(result.msd.iter().copied()).collect();
            chart.draw_series(LineSeries::new(pts, style(c, 2)))?
                .label("MSD total")
                .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], style(c, 2)));
        }

        // 拟合直线 + D / R²（黑色实线，避开调色板）
        if let Some(f) = &result.fit {
            let s = ShapeStyle { color: BLACK.to_rgba(), filled: false, stroke_width: 2 };
            let line = vec![
                (f.t_lo, f.slope * f.t_lo + f.intercept),
                (f.t_hi, f.slope * f.t_hi + f.intercept),
            ];
            let label =
                format!("fit: D={:.3e} Å²/fs (R²={:.4})", f.d_ang2_per_fs, f.r2);
            chart.draw_series(LineSeries::new(line, s))?
                .label(label)
                .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], s));
        }

        chart.configure_series_labels()
            .background_style(WHITE.mix(0.85))
            .border_style(BLACK)
            .position(SeriesLabelPosition::UpperLeft)
            .draw()?;

        root.present()?;
    }
    Ok(out)
}
```

- [ ] **Step 3: Build to verify it compiles**

Run: `cargo build --package ferro`
Expected: builds with no errors (warnings about unused `plot_msd` are acceptable until Task 6 wires it).

- [ ] **Step 4: Commit**

```bash
git add ferro-cli/src/plot.rs
git commit -m "Add plot_msd: total + a/b/c curves with fit overlay

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 6: Wire `--fit-range` and `plot_msd` into `fe-traj`

**Files:**
- Modify: `ferro-cli/src/bin/traj.rs`

- [ ] **Step 1: Import `plot_msd`**

Change the `plot::{...}` import (lines 7) from:

```rust
    plot::{open_plot, plot_angle, plot_gr, plot_sq},
```

to:

```rust
    plot::{open_plot, plot_angle, plot_gr, plot_msd, plot_sq},
```

- [ ] **Step 2: Add the `--fit-range` CLI arg**

In the `// ── msd ───` block of `struct Cli` (after the `elements` field, ~line 88), add:

```rust
    /// [msd] Linear-fit window as trajectory fractions FMIN,FMAX (e.g. 0.3,0.8) -> self-diffusion D
    #[arg(long, value_delimiter = ',')]
    fit_range: Option<Vec<f64>>,
```

- [ ] **Step 3: Replace `run_msd`**

Replace the entire `run_msd` function (~lines 244-257) with:

```rust
fn run_msd(args: &Cli, traj: &ferro_core::Trajectory) -> Result<()> {
    let fit_range = match &args.fit_range {
        None => None,
        Some(v) if v.len() == 2 => Some((v[0], v[1])),
        Some(v) => return Err(anyhow!(
            "--fit-range expects exactly two comma-separated fractions, e.g. 0.3,0.8 (got {} value(s))",
            v.len()
        )),
    };

    let params = MsdParams {
        dt: args.dt,
        shift: args.shift,
        elements: args.elements.clone(),
        fit_range,
        ..MsdParams::default()
    };
    let result = calc_msd(traj, &params)?;

    let out = args.output.as_deref().unwrap_or(Path::new("msd.dat"));
    let out_str = out.to_str().unwrap_or("msd.dat");
    write_msd(&result, out_str)?;
    println!("MSD -> {out_str}");

    if let Some(f) = &result.fit {
        println!(
            "Fit  range [{:.2}, {:.2}]  ->  t in [{:.1}, {:.1}] fs  ({} pts)",
            f.frac_lo, f.frac_hi, f.t_lo, f.t_hi, f.n_points
        );
        println!(
            "D (total) = {:.6e} Ang^2/fs = {:.6e} cm^2/s = {:.6e} m^2/s  (R^2={:.4})",
            f.d_ang2_per_fs,
            f.d_ang2_per_fs * 0.1,
            f.d_ang2_per_fs * 1e-5,
            f.r2
        );
    }

    if args.plot {
        let png = plot_msd(&result, out_str)?;
        println!("Plot -> {png}");
        open_plot(&png);
    }
    Ok(())
}
```

- [ ] **Step 4: Build to verify it compiles**

Run: `cargo build --package ferro`
Expected: builds with no errors or warnings about `plot_msd`/`fit_range`.

- [ ] **Step 5: Commit**

```bash
git add ferro-cli/src/bin/traj.rs
git commit -m "fe-traj: --fit-range arg + wire plot_msd/D output into run_msd

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 7: Update `print_msd` help text

**Files:**
- Modify: `ferro-cli/src/help.rs`

- [ ] **Step 1: Replace `print_msd`**

Replace the entire `print_msd` function (~lines 345-363) with:

```rust
fn print_msd() {
    println!(
        r#"fe-traj -m msd — Mean Square Displacement
  Computes MSD(t) = <|r(t₀+t) − r(t₀)|²> averaged over time origins.
  Outputs total MSD and per-axis (a/b/c) components.

Parameters:
  --dt        FLOAT      Timestep between frames [fs]   default: 1.0
  --shift     INT        Time-origin stride             default: 1
  --elements  Fe,O,...   Track only these elements      default: all
  --fit-range FMIN,FMAX  Linear-fit window as fractions of the MSD
                         curve; reports self-diffusion D = slope/6
                         (Einstein, 3-D) and R²
  --last-n    INT        Use only the last N frames
  --ncore     INT        Parallel threads
  --plot                 Generate PNG and open in viewer
  -o PATH                Output file                    default: msd.dat

Example:
  fe-traj -m msd -i traj.xyz --dt 2.0
  fe-traj -m msd -i traj.dump --elements Li --dt 1.0 --last-n 2000
  fe-traj -m msd -i traj.dump --dt 1.0 --fit-range 0.3,0.8 --plot"#
    );
}
```

- [ ] **Step 2: Build to verify it compiles**

Run: `cargo build --package ferro`
Expected: builds with no errors.

- [ ] **Step 3: Sanity-check the help output**

Run: `cargo run --quiet --bin fe-traj -- -m msd`
Expected: prints the new help text including the `--fit-range` and `--plot` lines and the third example.

- [ ] **Step 4: Commit**

```bash
git add ferro-cli/src/help.rs
git commit -m "help: document --fit-range and --plot for fe-traj -m msd

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 8: Update `CLAUDE.md` docs

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add `--fit-range` to the msd-specific flag block**

In `CLAUDE.md`, locate the `# msd-specific` block:

```
# msd-specific
--dt    FLOAT   timestep [fs]              default 1.0
--shift INT     time-origin stride         default 1
--elements El,El,...  filter elements
```

Add one line after the `--elements` line so the block becomes:

```
# msd-specific
--dt    FLOAT   timestep [fs]              default 1.0
--shift INT     time-origin stride         default 1
--elements El,El,...  filter elements
--fit-range FMIN,FMAX  linear-fit window (traj. fractions) → self-diffusion D
```

- [ ] **Step 2: Widen the plot-scope note**

Locate:

```
# plot (gr / sq / angle only)
--plot   write PNG next to output file and open with system viewer
```

Replace the comment line so it reads:

```
# plot (gr / sq / angle / msd)
--plot   write PNG next to output file and open with system viewer
```

- [ ] **Step 3: Add an `msd` bullet to the "Plot behaviour:" list**

Locate the `**Plot behaviour:**` list (the bullets for `gr`, `sq`, `angle`). Add a fourth bullet after the `angle` bullet:

```
- `msd`: total MSD + a/b/c components; with `--fit-range` overlays the linear fit and prints D = slope/6 + R²
```

- [ ] **Step 4: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: document MSD --fit-range and --plot in CLAUDE.md

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 9: Full verification + manual smoke test

**Files:** none modified unless `cargo fmt` reformats.

- [ ] **Step 1: Format**

Run: `cargo fmt`
Then: `git diff --stat`
Expected: either no changes, or only whitespace reformatting in the touched files.

- [ ] **Step 2: Clippy on the touched crates**

Run: `cargo clippy --package ferro-analysis --package ferro --all-targets`
Expected: no *new* warnings pointing at the files changed in this plan
(`msd.rs`, `plot.rs`, `traj.rs`, `help.rs`). Pre-existing warnings elsewhere
are out of scope — fix only clippy findings in the files this plan touched,
then re-run.

- [ ] **Step 3: Full test suite**

Run: `cargo test`
Expected: all tests pass, including the four new MSD tests
(`test_fit_diffusion_exact_line`, `test_fit_range_index_mapping`,
`test_fit_range_invalid`, `test_calc_msd_populates_fit`, `test_write_msd_with_fit`).

- [ ] **Step 4: Create a tiny multi-frame trajectory fixture**

Run:

```bash
cat > /tmp/msd_smoke.xyz <<'EOF'
2
frame 0
Li 0.0 0.0 0.0
O  3.0 0.0 0.0
2
frame 1
Li 0.5 0.0 0.0
O  3.0 0.0 0.0
2
frame 2
Li 1.0 0.0 0.0
O  3.0 0.0 0.0
2
frame 3
Li 1.5 0.0 0.0
O  3.0 0.0 0.0
EOF
```

- [ ] **Step 5: Run `fe-traj -m msd` with `--fit-range` (no plot)**

Run:
```bash
cargo run --quiet --bin fe-traj -- -m msd -i /tmp/msd_smoke.xyz -o /tmp/msd_smoke.dat --dt 1.0 --fit-range 0.0,1.0
```
Expected: stdout shows `MSD -> /tmp/msd_smoke.dat`, a `Fit range [...]` line, and a `D (total) = ... Ang^2/fs = ... cm^2/s = ... m^2/s (R^2=...)` line. No plot line.

- [ ] **Step 6: Confirm the fit header is in the file**

Run: `grep -E "D \(total\)|R\^2|fit range" /tmp/msd_smoke.dat`
Expected: prints the `# fit range`, `# D (total)`, and `# R^2` header lines.

- [ ] **Step 7: Run with `--plot` and confirm the PNG**

Run:
```bash
cargo run --quiet --bin fe-traj -- -m msd -i /tmp/msd_smoke.xyz -o /tmp/msd_smoke.dat --dt 1.0 --fit-range 0.2,0.9 --plot
ls -l /tmp/msd_smoke.png
```
Expected: stdout includes `Plot -> /tmp/msd_smoke.png`; `ls` shows a non-empty PNG file.

- [ ] **Step 8: Negative check — invalid range fails cleanly**

Run:
```bash
cargo run --quiet --bin fe-traj -- -m msd -i /tmp/msd_smoke.xyz --dt 1.0 --fit-range 0.9,0.2
```
Expected: process exits non-zero with an error mentioning `0 <= fmin < fmax <= 1`.

- [ ] **Step 9: Final commit (only if fmt changed our files)**

```bash
git add ferro-analysis/src/md/msd.rs ferro-cli/src/plot.rs ferro-cli/src/bin/traj.rs ferro-cli/src/help.rs
git commit -m "Apply cargo fmt for MSD plot/fit feature

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

(If `git status --porcelain` shows none of those four files modified, skip this
commit. Do not stage unrelated files reformatted elsewhere by `cargo fmt`.)

---

## Notes for the implementer

- `ferro` is the package name of `ferro-cli` (binaries `fe-traj` etc. live there); `cargo build --package ferro` builds the CLI. If `--package ferro` errors with "package not found", run `cargo build` (whole workspace) instead.
- The early validation in `calc_msd` (Task 3 Step 3) intentionally duplicates the bounds check inside `fit_diffusion`. This is deliberate: it fails fast *before* the expensive parallel time-origin loop for the common user mistake (swapped/out-of-range fractions). Do not "DRY it away" by removing the early check.
- Do not bump any version number and do not touch dev-record / progress files — that is an explicit standing rule for this project and is out of scope here.
- `MsdParams` is `#[derive(Debug, Clone)]`; `MsdFit` must also derive `Debug, Clone` (already specified) because `MsdResult` derives them.
