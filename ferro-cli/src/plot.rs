//! Simple PNG plots for GR, angle, and SQ results.

use anyhow::Result;
use plotters::prelude::*;
use ferro_analysis::{AngleResult, GrResult, SqResult};
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
pub fn plot_gr(result: &GrResult, dat_path: &str) -> Result<String> {
    let out = png_path(dat_path);

    let pair_keys: Vec<&String> = result.gr.keys()
        .filter(|k| k.as_str() != "total")
        .collect();

    let y1_max = pair_keys.iter()
        .filter_map(|k| result.gr.get(*k))
        .flat_map(|v| v.iter().copied())
        .fold(0.0f64, f64::max) * 1.1;
    let y1_max = y1_max.max(1.5);

    let y2_max = pair_keys.iter()
        .filter_map(|k| result.cn.get(*k))
        .flat_map(|v| v.iter().copied())
        .fold(0.0f64, f64::max) * 1.1;
    let y2_max = y2_max.max(1.0);

    {
        let root = BitMapBackend::new(&out, (900, 540)).into_drawing_area();
        root.fill(&WHITE)?;

        let mut chart = ChartBuilder::on(&root)
            .caption("g(r)  and  CN(r)", ("sans-serif", 18))
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

        // total g(r): 灰色细线垫底，不画 CN
        if let Some(g) = result.gr.get("total") {
            let pts: Vec<(f64, f64)> = result.r.iter().copied().zip(g.iter().copied()).collect();
            chart.draw_series(LineSeries::new(pts, style(PALETTE[7], 1)))?
                .label("total")
                .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], style(PALETTE[7], 1)));
        }

        // 每个 pair：g(r) 实线（左轴）+ CN 淡色细线（右轴）
        for (i, key) in pair_keys.iter().enumerate() {
            let c = color(i);

            if let Some(g) = result.gr.get(*key) {
                let pts: Vec<(f64, f64)> = result.r.iter().copied().zip(g.iter().copied()).collect();
                chart.draw_series(LineSeries::new(pts, style(c, 2)))?
                    .label(key.as_str())
                    .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], style(c, 2)));
            }

            if let Some(cn) = result.cn.get(*key) {
                let pts: Vec<(f64, f64)> = result.r.iter().copied().zip(cn.iter().copied()).collect();
                let s = ShapeStyle { color: c.mix(0.55), filled: false, stroke_width: 1 };
                let s2 = s.clone();
                chart.draw_secondary_series(LineSeries::new(pts, s))?
                    .label(format!("{key}  CN"))
                    .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], s2.clone()));
            }
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

    // 只绘制 total_xrd / total_neutron，附上可读标签
    const SHOW: &[(&str, &str)] = &[
        ("total_xrd",     "XRD"),
        ("total_neutron", "Neutron"),
    ];

    let visible: Vec<(&str, &str, &Vec<f64>)> = SHOW.iter()
        .filter_map(|&(key, label)| result.sq.get(key).map(|v| (key, label, v)))
        .collect();

    let all_vals: Vec<f64> = visible.iter()
        .flat_map(|(_, _, v)| v.iter().copied())
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

        for (i, (_key, label, vals)) in visible.iter().enumerate() {
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

// ─── 系统查看器 ──────────────────────────────────────────────────────────────

/// Open PNG with the OS default viewer (non-blocking).
pub fn open_plot(path: &str) {
    #[cfg(target_os = "macos")]
    { std::process::Command::new("open").arg(path).spawn().ok(); }
    #[cfg(target_os = "linux")]
    { std::process::Command::new("xdg-open").arg(path).spawn().ok(); }
}
