//! Quick-look PNG plots.
//!
//! **These are self-check plots, not publication figures.** The data files are
//! long-format CSV, so a real figure is one line of seaborn
//! (`sns.lineplot(data=df, x="r", y="gr", hue="file")`). Chasing matplotlib here is an
//! arms race with no end; requests for log axes, error bars or themes belong in Python.
//!
//! Layout model: one PNG, split into panels — **one panel per quantity, one curve per
//! input file**. The file is the only variable dimension, so colours are assigned per
//! file and stay consistent across panels; only the first panel carries a legend.

use anyhow::{bail, Result};
use plotters::prelude::*;
use ferro_analysis::{AngleResult, GrResult, MsdResult, SqResult};

// matplotlib tab10 色板
const PALETTE: &[RGBColor] = &[
    RGBColor(31, 119, 180),
    RGBColor(255, 127, 14),
    RGBColor(44, 160, 44),
    RGBColor(214, 39, 40),
    RGBColor(148, 103, 189),
    RGBColor(140, 86, 75),
    RGBColor(227, 119, 194),
    RGBColor(127, 127, 127),
    RGBColor(188, 189, 34),
    RGBColor(23, 190, 207),
];

fn color(i: usize) -> RGBColor { PALETTE[i % PALETTE.len()] }

fn style(c: RGBColor, width: u32) -> ShapeStyle {
    ShapeStyle { color: c.to_rgba(), filled: false, stroke_width: width }
}

/// One curve: a file's values for one quantity.
pub struct Series {
    pub label: String,
    pub points: Vec<(f64, f64)>,
}

/// One sub-plot: a single quantity, with one [`Series`] per input file.
pub struct Panel {
    pub title: String,
    pub x_desc: String,
    pub y_desc: String,
    pub series: Vec<Series>,
}

/// PNG name derived from the data file: `gr_run1.csv` → `gr_run1.png`.
fn png_path(data_path: &str) -> String {
    let p = std::path::Path::new(data_path);
    let stem = p.file_stem().and_then(|s| s.to_str()).unwrap_or("out");
    match p.parent().and_then(|d| d.to_str()).filter(|s| !s.is_empty()) {
        Some(dir) => format!("{dir}/{stem}.png"),
        None      => format!("{stem}.png"),
    }
}

/// Data range with a small margin; falls back to a unit range for empty/flat data.
fn range(vals: impl Iterator<Item = f64>) -> (f64, f64) {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for v in vals.filter(|v| v.is_finite()) {
        lo = lo.min(v);
        hi = hi.max(v);
    }
    if !lo.is_finite() || !hi.is_finite() {
        return (0.0, 1.0);
    }
    let pad = ((hi - lo) * 0.05).max(1e-9);
    (lo - pad, hi + pad)
}

/// Renders `panels` into one PNG, `cols` panels per row.
fn render(out: &str, panels: &[Panel], cols: usize) -> Result<()> {
    if panels.is_empty() {
        bail!("nothing to plot");
    }
    let rows = panels.len().div_ceil(cols);
    let (w, h) = (520 * cols as u32, 400 * rows as u32);

    {
        let root = BitMapBackend::new(out, (w, h)).into_drawing_area();
        root.fill(&WHITE)?;
        let areas = root.split_evenly((rows, cols));

        for (pi, panel) in panels.iter().enumerate() {
            let area = &areas[pi];
            let (x_lo, x_hi) = range(panel.series.iter().flat_map(|s| s.points.iter().map(|p| p.0)));
            let (y_lo, y_hi) = range(panel.series.iter().flat_map(|s| s.points.iter().map(|p| p.1)));

            let mut chart = ChartBuilder::on(area)
                .caption(&panel.title, ("sans-serif", 16))
                .margin(12)
                .x_label_area_size(42)
                .y_label_area_size(58)
                .build_cartesian_2d(x_lo..x_hi, y_lo..y_hi)?;

            chart.configure_mesh()
                .x_desc(&panel.x_desc)
                .y_desc(&panel.y_desc)
                .x_labels(8)
                .y_labels(6)
                .draw()?;

            for (si, s) in panel.series.iter().enumerate() {
                let c = color(si);
                let st = style(c, 2);
                let pts: Vec<(f64, f64)> =
                    s.points.iter().copied().filter(|(x, y)| x.is_finite() && y.is_finite()).collect();
                let series = chart.draw_series(LineSeries::new(pts, st))?;
                // 图例只画在第一格:颜色按文件分配且各格一致,重复画只是占地方
                if pi == 0 {
                    series
                        .label(s.label.clone())
                        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 18, y)], st));
                }
            }

            if pi == 0 && panel.series.len() > 1 {
                chart.configure_series_labels()
                    .background_style(WHITE.mix(0.85))
                    .border_style(BLACK)
                    .position(SeriesLabelPosition::UpperRight)
                    .draw()?;
            }
        }
        root.present()?;
    }
    Ok(())
}

/// A panel definition: title, axis labels, and how to pull `(x, y)` out of one result.
type Quantity<T> = (&'static str, &'static str, &'static str, fn(&T) -> (Vec<f64>, Vec<f64>));

/// Builds one panel per quantity from a per-file extractor.
///
/// `quantities` names each panel and supplies `(x, y)` for one file.
fn panels_from<T>(
    results: &[(String, &T)],
    quantities: &[Quantity<T>],
) -> Vec<Panel> {
    quantities
        .iter()
        .map(|(title, x_desc, y_desc, extract)| Panel {
            title: (*title).to_string(),
            x_desc: (*x_desc).to_string(),
            y_desc: (*y_desc).to_string(),
            series: results
                .iter()
                .map(|(label, r)| {
                    let (xs, ys) = extract(r);
                    Series {
                        label: label.clone(),
                        points: xs.into_iter().zip(ys).collect(),
                    }
                })
                .collect(),
        })
        .collect()
}

// ─── GR ─────────────────────────────────────────────────────────────────────

/// g(r) and CN(r) as two panels, one curve per input file.
///
/// The old dual-Y-axis layout is gone: with several files there would be 2N curves
/// fighting over two scales. Separate panels also drop the axis-side ambiguity.
pub fn plot_gr(results: &[(String, &GrResult)], data_path: &str, a: &str, b: &str) -> Result<String> {
    let out = png_path(data_path);
    let key = format!("{a}-{b}");

    let mut gr_series = Vec::new();
    let mut cn_series = Vec::new();
    for (label, r) in results {
        let g = r.gr.get(&key)
            .ok_or_else(|| anyhow::anyhow!("pair '{key}' not present in {label}"))?;
        let c = r.cn.get(&key)
            .ok_or_else(|| anyhow::anyhow!("pair '{key}' not present in {label}"))?;
        gr_series.push(Series { label: label.clone(), points: r.r.iter().copied().zip(g.iter().copied()).collect() });
        cn_series.push(Series { label: label.clone(), points: r.r.iter().copied().zip(c.iter().copied()).collect() });
    }

    let panels = vec![
        Panel { title: format!("g(r)  {a}-{b}"), x_desc: "r  [Å]".into(), y_desc: "g(r)".into(), series: gr_series },
        Panel { title: format!("CN(r)  {a} -> {b}"), x_desc: "r  [Å]".into(), y_desc: "CN(r)".into(), series: cn_series },
    ];
    render(&out, &panels, 2)?;
    Ok(out)
}

// ─── SQ ─────────────────────────────────────────────────────────────────────

/// The two weighted totals as two panels — they are different quantities, not two
/// styles of the same one, so they do not share an axis.
pub fn plot_sq(results: &[(String, &SqResult)], data_path: &str) -> Result<String> {
    let out = png_path(data_path);

    let collect = |pick: fn(&SqResult) -> Option<&Vec<f64>>| -> Vec<Series> {
        results
            .iter()
            .filter_map(|(label, r)| {
                pick(r).map(|v| Series {
                    label: label.clone(),
                    points: r.q.iter().copied().zip(v.iter().copied()).collect(),
                })
            })
            .collect()
    };

    let mut panels = Vec::new();
    let xrd = collect(|r| r.total_xrd.as_ref());
    if !xrd.is_empty() {
        panels.push(Panel { title: "S(q)  XRD-weighted".into(), x_desc: "q  [Å⁻¹]".into(), y_desc: "S(q)".into(), series: xrd });
    }
    let neu = collect(|r| r.total_neutron.as_ref());
    if !neu.is_empty() {
        panels.push(Panel { title: "S(q)  neutron-weighted".into(), x_desc: "q  [Å⁻¹]".into(), y_desc: "S(q)".into(), series: neu });
    }
    if panels.is_empty() {
        bail!("no weighted total to plot (use --weighting xrd|neutron|both)");
    }

    render(&out, &panels, panels.len())?;
    Ok(out)
}

// ─── MSD ────────────────────────────────────────────────────────────────────

/// Total MSD plus the three lattice-direction components, 2×2.
pub fn plot_msd(results: &[(String, &MsdResult)], data_path: &str) -> Result<String> {
    let out = png_path(data_path);
    let quantities: &[Quantity<MsdResult>] = &[
        ("MSD total", "t  [fs]", "MSD  [Å²]", |r| (r.time.clone(), r.msd.clone())),
        ("MSD a",     "t  [fs]", "MSD  [Å²]", |r| (r.time.clone(), r.msd_a.clone())),
        ("MSD b",     "t  [fs]", "MSD  [Å²]", |r| (r.time.clone(), r.msd_b.clone())),
        ("MSD c",     "t  [fs]", "MSD  [Å²]", |r| (r.time.clone(), r.msd_c.clone())),
    ];
    let panels = panels_from(results, quantities);
    render(&out, &panels, 2)?;
    Ok(out)
}

// ─── Angle ──────────────────────────────────────────────────────────────────

/// One named triplet's angle distribution, one curve per input file.
///
/// Requires a named triplet for the same reason g(r) requires a named pair: without
/// one the figure would carry every triplet of every file at once.
pub fn plot_angle(results: &[(String, &AngleResult)], data_path: &str, triplet: &str) -> Result<String> {
    let out = png_path(data_path);

    let mut series = Vec::new();
    for (label, r) in results {
        let Some(h) = r.hist.get(triplet) else { continue };
        let total: f64 = h.iter().sum::<u64>() as f64;
        if total <= 0.0 { continue; }
        let stat = r.stats.get(triplet);
        let legend = match stat {
            Some(s) => format!("{label}  ({:.1}±{:.1}°)", s.mean, s.std),
            None => label.clone(),
        };
        series.push(Series {
            label: legend,
            points: r.angle.iter().copied().zip(h.iter().map(|&c| c as f64 / total)).collect(),
        });
    }
    if series.is_empty() {
        bail!("triplet '{triplet}' has no counts in any input");
    }

    let panels = vec![Panel {
        title: format!("Angle distribution  {triplet}"),
        x_desc: "θ  [°]".into(),
        y_desc: "P(θ)".into(),
        series,
    }];
    render(&out, &panels, 1)?;
    Ok(out)
}

// ─── 打开查看器 ──────────────────────────────────────────────────────────────

/// Opens a PNG with the system viewer. Only called for a single input — a batch would
/// otherwise fill the screen with windows.
pub fn open_plot(path: &str) {
    #[cfg(target_os = "macos")]
    let _ = std::process::Command::new("open").arg(path).spawn();
    #[cfg(target_os = "linux")]
    let _ = std::process::Command::new("xdg-open").arg(path).spawn();
    #[cfg(target_os = "windows")]
    let _ = std::process::Command::new("cmd").args(["/C", "start", "", path]).spawn();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_png_path_replaces_extension() {
        assert_eq!(png_path("gr_run1.csv"), "gr_run1.png");
        assert_eq!(png_path("out/sq.csv"), "out/sq.png");
    }

    #[test]
    fn test_range_pads_and_survives_degenerate_input() {
        let (lo, hi) = range([1.0, 3.0].into_iter());
        assert!(lo < 1.0 && hi > 3.0, "范围要留边距");
        let (lo2, hi2) = range([f64::NAN].into_iter());
        assert!(lo2 < hi2, "全 NaN 时回退到可用范围而不是 panic");
    }

    #[test]
    fn test_render_rejects_empty_panels() {
        assert!(render("/tmp/ferro_empty.png", &[], 2).is_err());
    }
}
