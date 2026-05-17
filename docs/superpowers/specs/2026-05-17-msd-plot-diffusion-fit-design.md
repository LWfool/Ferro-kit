# MSD `plot_msd` + self-diffusion coefficient fit — Design

Date: 2026-05-17
Branch: `feature/msd-plot-diffusion-fit`

## Problem

`fe-traj -m msd` has no `--plot` support: the global `--plot` flag is parsed
but silently ignored in MSD mode (`run_msd` never reads `args.plot`, and there
is no `plot_msd` function — `plot.rs` only implements gr/sq/angle). Users also
have no way to extract a self-diffusion coefficient from an MSD curve.

## Goals

1. Add `plot_msd` and wire `--plot` into MSD mode.
2. Add `--fit-range FMIN,FMAX`: the user picks the linear diffusive segment of
   the MSD curve; ferro fits it and reports the self-diffusion coefficient `D`.

## Decisions (locked during brainstorming)

- **Fit-range unit:** trajectory fraction. `FMIN,FMAX ∈ [0,1]` map onto the
  **lag-time axis of the MSD curve** (the plotted x-axis, `0 .. (tau-1)·dt`),
  *not* raw frame count. `0.3,0.8` means "skip the early ballistic 30%, drop
  the noisy last 20%, fit the middle". Index window:
  `i_lo = round(fmin·(tau-1))`, `i_hi = round(fmax·(tau-1))`.
- **D reported:** total only. Ordinary least-squares of total MSD vs time over
  the window; `D = slope / 6` (Einstein relation, 3-D isotropic). Also report
  `R²`. No per-axis D.
- **When D is computed:** whenever `--fit-range` is given, independent of
  `--plot`. Result printed to console and appended to the `msd.dat` header.
  When `--plot` is also set, the fitted line and `D`/`R²` are overlaid on the
  plot.
- **Architecture:** fit logic lives in `ferro-analysis` (pure, tested,
  CLI-agnostic, reusable by `ferro-python`), not in the CLI entry layer.

## Unit conversions

Internal slope unit is `Å²/fs`. Reported as:

- `D [Å²/fs]` (internal)
- `D [cm²/s] = D[Å²/fs] × 0.1`
- `D [m²/s]  = D[Å²/fs] × 1e-5`

## Changes by file

### `ferro-analysis/src/md/msd.rs`

- `MsdParams` gains `fit_range: Option<(f64, f64)>` (default `None`).
- New public struct:

  ```rust
  pub struct MsdFit {
      pub frac_lo: f64,
      pub frac_hi: f64,
      pub t_lo: f64,            // fs
      pub t_hi: f64,            // fs
      pub slope: f64,           // Å²/fs
      pub intercept: f64,       // Å²
      pub d_ang2_per_fs: f64,   // slope / 6
      pub r2: f64,
      pub n_points: usize,
  }
  ```

- `MsdResult` gains `fit: Option<MsdFit>`.
- New pure fn
  `fit_diffusion(time: &[f64], msd: &[f64], frac: (f64, f64)) -> ferro_core::Result<MsdFit>`:
  ordinary least-squares (slope, intercept) of `msd` vs `time` over the index
  window, `d_ang2_per_fs = slope / 6.0`, plus `R²` (coefficient of
  determination). Validation (fail fast, *before* the heavy parallel loop in
  `calc_msd`): require `0.0 ≤ fmin < fmax ≤ 1.0` and a window of ≥2 points,
  else `ChemError::ValidationError`.
- `calc_msd` calls `fit_diffusion` after `build_result` when
  `params.fit_range` is `Some`, storing the result in `MsdResult.fit`.
- `write_msd`: when `fit` is `Some`, append a header block before the column
  header (data columns and numeric format unchanged):

  ```
  # ----------------------------------------------------------
  # fit range  = [0.30, 0.80]  ->  t in [2000.0, 10000.0] fs
  # points     = 6001
  # slope      = 1.234000e-04 Ang^2/fs
  # D (total)  = 2.056667e-05 Ang^2/fs
  #            = 2.056667e-06 cm^2/s = 2.056667e-10 m^2/s
  # R^2        = 0.998700
  # ----------------------------------------------------------
  ```

### `ferro-cli/src/plot.rs`

- Add `MsdResult` to the `use ferro_analysis::{...}` import.
- New `pub fn plot_msd(result: &MsdResult, dat_path: &str) -> Result<String>`:
  - x = time [fs], y = MSD [Å²].
  - Total MSD as the prominent curve; `msd_a/b/c` as thin lighter lines.
  - When `result.fit` is `Some`: overlay the fitted straight line across
    `[t_lo, t_hi]` and add a legend entry
    `fit: D=2.06e-05 Å²/fs (R²=0.999)`.
  - PNG path via existing `png_path`; reuses the tab10 `PALETTE`/`style`
    helpers.

### `ferro-cli/src/bin/traj.rs`

- New arg:

  ```rust
  /// [msd] Linear-fit range as trajectory fractions, e.g. 0.3,0.8 -> D
  #[arg(long, value_delimiter = ',', num_args = 2)]
  fit_range: Option<Vec<f64>>,
  ```

  Mapped to `Option<(f64, f64)>` (require exactly 2 values; `clap` enforces
  `num_args = 2`).
- `run_msd`: build `MsdParams.fit_range`; after `calc_msd`, if `result.fit` is
  `Some`, print the D summary to the console (always, regardless of `--plot`);
  if `args.plot`, call `plot_msd` + `open_plot` and print `Plot ->`.

### `ferro-cli/src/help.rs`

- `print_msd()` parameter list gains:
  - `--fit-range FMIN,FMAX  Linear-fit range (traj. fractions) -> self-diffusion D`
  - `--plot                 Generate PNG and open in viewer`
  - example: `fe-traj -m msd -i traj.dump --dt 1.0 --fit-range 0.3,0.8 --plot`

### `CLAUDE.md`

- msd-specific flag block: add `--fit-range` and `--plot`.
- Change the plot-scope note from "gr / sq / angle only" to include `msd`, and
  add an `msd` bullet to the "Plot behaviour" list (total + a/b/c curves; fitted
  line + D/R² when `--fit-range` given).
- No version bump, no dev-record/progress edits (per standing rule, those wait
  until explicitly requested).

## Tests (TDD, in `msd.rs` `#[cfg(test)]`)

1. `test_fit_diffusion_exact_line`: synthetic `msd = 6·D·t + c` with known
   `D`; assert recovered `slope/6 ≈ D` and `R² ≈ 1.0` to tight tolerance.
2. `test_fit_range_index_mapping`: given `tau`/`dt`, assert
   `(fmin,fmax)` → expected `i_lo,i_hi,t_lo,t_hi,n_points`.
3. `test_fit_range_invalid`: `fmin ≥ fmax`, out of `[0,1]`, and a window with
   <2 points each return `ChemError::ValidationError`.
4. `test_calc_msd_populates_fit`: `calc_msd` with `fit_range = Some` yields
   `result.fit.is_some()`; with `None` yields `result.fit.is_none()`.

Existing MSD tests must still pass (they construct `MsdParams` with explicit
fields — the new field is added there too).

## Out of scope (YAGNI)

Per-axis D, weighted/robust regression, automatic linear-regime detection.
Note: `--plot` is still allowed without `--fit-range` — it produces the MSD
plot with no fit overlay.
