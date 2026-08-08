//! Simple PNG plots for GR, angle, and SQ results.

use anyhow::Result;
use plotters::prelude::*;
use ferro_analysis::{AngleResult, GrResult, MsdResult, SqResult};
use std::path::Path;

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

fn png_path(dat_path: &str) -> String {
    let p = Path::new(dat_path);
    let stem = p.file_stem().and_then(|s| s.to_str()).unwrap_or("out");
    match p.parent().and_then(|d| d.to_str()).filter(|s| !s.is_empty()) {
        Some(dir) => format!("{dir}/{stem}.png"),
        None      => format!("{stem}.png"),
    }
}

fn style(c: RGBColor, width: u32) -> ShapeStyle {
    ShapeStyle { color: c.to_rgba(), filled: false, stroke_width: width }
}

// ─── GR ─────────────────────────────────────────────────────────────────────

/// Plot g(r) on the left axis and CN(r) on the right axis, saved as PNG.
pub fn plot_gr(result: &GrResult, dat_path: &str, a: &str, b: &str) -> Result<String> {
    let out = png_path(dat_path);
    let key = format!("{a}-{b}");

    let g = result.gr.get(&key)
        .ok_or_else(|| anyhow::anyhow!("pair '{key}' not present in the trajectory"))?;
    let cn = result.cn.get(&key)
        .ok_or_else(|| anyhow::anyhow!("pair '{key}' not present in the trajectory"))?;

    let y1_max = (g.iter().copied().fold(0.0f64, f64::max) * 1.1).max(1.5);
    let y2_max = (cn.iter().copied().fold(0.0f64, f64::max) * 1.1).max(1.0);

    {
        let root = BitMapBackend::new(&out, (900, 540)).into_drawing_area();
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .caption(format!("g(r) and CN(r)   {a} -> {b}"), ("sans-serif", 18))
            .margin(15)
            .x_label_area_size(45)
            .y_label_area_size(60)
            .right_y_label_area_size(60)
            .build_cartesian_2d(result.params.r_min..result.params.r_max, 0.0..y1_max)?
            .set_secondary_coord(result.params.r_min..result.params.r_max, 0.0..y2_max);

        chart.configure_mesh()
            .x_desc("r  [Å]")
            .y_desc("g(r)")
            .x_labels(10).y_labels(8)
            .draw()?;

        chart.configure_secondary_axes()
            .y_desc("CN(r)")
            .draw()?;

        // g(r) 实线（左轴）
        let cg = color(0);
        let pts: Vec<(f64, f64)> = result.r.iter().copied().zip(g.iter().copied()).collect();
        chart.draw_series(LineSeries::new(pts, style(cg, 2)))?
            .label(format!("g(r)  {a}-{b}"))
            .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], style(cg, 2)));

        // CN(r) 淡线（右轴）。方向由 -a/-b 决定：每个 A 周围的 B 数
        let cc = color(1);
        let s = ShapeStyle { color: cc.mix(0.75), filled: false, stroke_width: 2 };
        let pts: Vec<(f64, f64)> = result.r.iter().copied().zip(cn.iter().copied()).collect();
        chart.draw_secondary_series(LineSeries::new(pts, s))?
            .label(format!("CN(r)  {a} -> {b}"))
            .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], s));

        chart.configure_series_labels()
            .background_style(WHITE.mix(0.85))
            .border_style(BLACK)
            .position(SeriesLabelPosition::UpperRight)
            .draw()?;

        root.present()?;
    }
    Ok(out)
}

// ─── Angle ──────────────────────────────────────────────────────────────────

/// Plot bond angle distributions (normalised to peak = 1) and save as PNG.
pub fn plot_angle(result: &AngleResult, dat_path: &str) -> Result<String> {
    let out = png_path(dat_path);

    let all_max = result.hist.values()
        .flat_map(|v| v.iter().copied())
        .max().unwrap_or(1) as f64;

    {
        let root = BitMapBackend::new(&out, (900, 540)).into_drawing_area();
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .caption("Bond Angle Distribution", ("sans-serif", 18))
            .margin(15)
            .x_label_area_size(45)
            .y_label_area_size(60)
            .build_cartesian_2d(0.0f64..180.0f64, 0.0f64..1.05f64)?;

        chart.configure_mesh()
            .x_desc("Angle  [°]")
            .y_desc("Normalised count")
            .x_labels(10).y_labels(8)
            .draw()?;

        let keys: Vec<&String> = result.hist.keys().collect();
        for (i, key) in keys.iter().enumerate() {
            if let Some(hist) = result.hist.get(*key) {
                let pts: Vec<(f64, f64)> = result.angle.iter().copied()
                    .zip(hist.iter().map(|&c| c as f64 / all_max))
                    .collect();
                let c = color(i);
                let label = match result.stats.get(*key) {
                    Some(s) => format!("{key}  ({:.1}° ± {:.1}°)", s.mean, s.std),
                    None    => key.to_string(),
                };
                chart.draw_series(LineSeries::new(pts, style(c, 2)))?
                    .label(label)
                    .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], style(c, 2)));
            }
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

// ─── SQ ─────────────────────────────────────────────────────────────────────

/// Plot S(q): only XRD-weighted and neutron-weighted total curves.
pub fn plot_sq(result: &SqResult, dat_path: &str) -> Result<String> {
    let out = png_path(dat_path);

    // 只绘制加权总曲线 —— 它们才是能和实验衍射谱对比的量
    let visible: Vec<(&str, &Vec<f64>)> = [
        ("XRD", result.total_xrd.as_ref()),
        ("Neutron", result.total_neutron.as_ref()),
    ].into_iter()
        .filter_map(|(label, v)| v.map(|v| (label, v)))
        .collect();

    let all_vals: Vec<f64> = visible.iter()
        .flat_map(|(_, v)| v.iter().copied())
        .collect();
    let y_min = all_vals.iter().copied().fold(f64::INFINITY, f64::min);
    let y_max = all_vals.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let pad = (y_max - y_min).max(0.5) * 0.1;

    {
        let root = BitMapBackend::new(&out, (900, 540)).into_drawing_area();
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .caption("Structure Factor  S(q)", ("sans-serif", 18))
            .margin(15)
            .x_label_area_size(45)
            .y_label_area_size(60)
            .build_cartesian_2d(
                result.params.q_min..result.params.q_max,
                (y_min - pad).min(0.0)..(y_max + pad),
            )?;

        chart.configure_mesh()
            .x_desc("q  [Å⁻¹]")
            .y_desc("S(q)")
            .x_labels(10).y_labels(8)
            .draw()?;

        for (i, (label, vals)) in visible.iter().enumerate() {
            let pts: Vec<(f64, f64)> = result.q.iter().copied().zip(vals.iter().copied()).collect();
            let c = color(i);
            chart.draw_series(LineSeries::new(pts, style(c, 2)))?
                .label(*label)
                .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], style(c, 2)));
        }

        chart.configure_series_labels()
            .background_style(WHITE.mix(0.85))
            .border_style(BLACK)
            .position(SeriesLabelPosition::UpperRight)
            .draw()?;

        root.present()?;
    }
    Ok(out)
}

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

// ─── 系统查看器 ──────────────────────────────────────────────────────────────

/// Open PNG with the OS default viewer (non-blocking).
pub fn open_plot(path: &str) {
    #[cfg(target_os = "macos")]
    { std::process::Command::new("open").arg(path).spawn().ok(); }
    #[cfg(target_os = "linux")]
    { std::process::Command::new("xdg-open").arg(path).spawn().ok(); }
}
