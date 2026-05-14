use anyhow::{bail, Result};
use clap::{CommandFactory, Parser, ValueEnum};
use ferro::io_dispatch::read_trajectory;
use ferro_analysis::{calc_network, NetworkResult};
use ferro_core::TypeParams;
use ferro_io::{write_xyz, LammpsUnits};
use ferro_structure::apply_type_labels;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, ValueEnum)]
enum NetworkMode {
    /// Qn 物种占比 + 氧类型占比 + CN 分布（多帧时间平均）
    Qn,
    /// 每个原子的类型标签（单帧，默认最后一帧）
    Type,
}

#[derive(Parser)]
#[command(
    name = "fe-network",
    about = "Glass network analysis: atom type labeling and Qn/CN/oxygen-type statistics",
    long_about = None,
    after_help = HELP_EXTRA,
)]
struct Cli {
    /// Input trajectory file (omit to show help with examples)
    #[arg(short, long)]
    input: Option<PathBuf>,

    /// Analysis mode
    #[arg(short, long = "mode", value_enum)]
    mode: Option<NetworkMode>,

    /// Output file.
    ///   -m Qn   → CSV/XLSX statistics (default: network.csv / network.xlsx)
    ///   -m type → structure file with type labels (omit to print to screen only)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output format for -m Qn: csv | xlsx  (default: csv)
    #[arg(long, default_value = "csv")]
    format: String,

    /// Use only the last N frames (default: all)
    #[arg(long)]
    last_n: Option<usize>,

    /// Frame index for -m type (0-based; default: last frame)
    #[arg(long)]
    frame: Option<usize>,

    /// Parallel threads (default: all cores)
    #[arg(long)]
    ncore: Option<usize>,

    /// Use LAMMPS metal units for dump files
    #[arg(long)]
    metal_units: bool,

    /// Modifier cation elements, comma-separated (e.g. Zn or Zn,Na).
    /// Supply their cutoffs via the same --Elem-O=cutoff syntax.
    #[arg(long)]
    modifier: Option<String>,
}

const HELP_EXTRA: &str = "\
PAIR ARGUMENTS:
  Specify cutoff radii with --Former-Ligand=<cutoff>:
    --P-O=2.3         P former, O ligand, cutoff 2.3 Å
    --Al-O=2.0 --Al-F=2.1

MODIFIER:
  --modifier Zn       classify Zn as modifier; supply cutoff via --Zn-O=3.5
  Modifier labels (by NBO count): Zn_f (0), Zn_t (1), Zn_b (2), X (≥3)

MODES:
  -m Qn   Time-averaged statistics over all frames:
            • Qn species distribution per former
            • Oxygen type distribution (Of / On_X / Ob_X_Y / X)
            • Total CN distribution per former
          Output: <stem>_qn.csv, <stem>_oxy.csv, <stem>_cn.csv  (or single xlsx)

  -m type  Per-atom type labels for one frame (default: last frame):
            Former → P0 P1 P2 …   (digit = bridging-O count)
            Free O → Of
            NBO    → On_P  On_Al …
            BO     → Ob_P_P  Ob_Al_P …
            Over-BO → X
            Modifier → Zn_f  Zn_t  Zn_b  X
           If -o is given, writes a structure file with labels as element names.
           If -o is omitted, prints the type table to stdout.

EXAMPLES:
  fe-network -i traj.dump --P-O=2.3 -m Qn
  fe-network -i traj.dump --P-O=2.3 -m Qn -o result.xlsx --format xlsx
  fe-network -i traj.dump --P-O=2.3 --Zn-O=3.5 --modifier Zn -m Qn
  fe-network -i traj.dump --P-O=2.3 -m type
  fe-network -i traj.dump --P-O=2.3 -m type -o typed.xyz";

fn main() -> Result<()> {
    let all_args: Vec<String> = std::env::args().collect();
    let (pair_args, clap_args) = split_pair_args(&all_args);
    let cli = Cli::parse_from(clap_args);

    if cli.input.is_none() || cli.mode.is_none() {
        println!("{}", Cli::command().render_long_help());
        return Ok(());
    }

    if pair_args.is_empty() {
        bail!("No pair cutoffs specified. Use --Former-Ligand=cutoff, e.g. --P-O=2.3");
    }

    // 解构 cli，避免后续部分移动问题
    let Cli { input, mode, output, format, last_n, frame, ncore, metal_units, modifier } = cli;
    let input  = input.unwrap();
    let mode   = mode.unwrap();

    let all_pairs = parse_pairs(&pair_args)?;
    let modifier_elems: std::collections::HashSet<String> = modifier
        .as_deref()
        .map(|s| s.split(',').map(|e| e.trim().to_string()).collect())
        .unwrap_or_default();

    let mut cutoffs = BTreeMap::new();
    let mut modifier_cutoffs = BTreeMap::new();
    for ((elem, ligand), cutoff) in all_pairs {
        if modifier_elems.contains(&elem) {
            modifier_cutoffs.insert((elem, ligand), cutoff);
        } else {
            cutoffs.insert((elem, ligand), cutoff);
        }
    }

    if let Some(n) = ncore {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
    }

    let units = if metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let mut traj = read_trajectory(&input, units)?;
    if let Some(n) = last_n {
        traj = traj.tail(n);
    }

    let params = TypeParams { cutoffs, modifier_cutoffs };

    match mode {
        NetworkMode::Qn   => run_qn(&traj, &params, output.as_deref(), &format)?,
        NetworkMode::Type => run_type(&traj, &params, output.as_deref(), frame)?,
    }
    Ok(())
}

// ─── -m Qn ───────────────────────────────────────────────────────────────────

fn run_qn(
    traj: &ferro_core::Trajectory,
    params: &TypeParams,
    output: Option<&Path>,
    format: &str,
) -> Result<()> {
    let result = calc_network(traj, params)
        .ok_or_else(|| anyhow::anyhow!(
            "Network analysis failed — trajectory empty or frames missing cell (PBC required)"
        ))?;

    let out_base = output.unwrap_or_else(|| Path::new("network"));
    match format.to_lowercase().as_str() {
        "csv"  => write_qn_csv(&result, params, out_base)?,
        "xlsx" => write_qn_xlsx(&result, params, out_base)?,
        other  => bail!("Unknown format '{other}'. Use csv or xlsx."),
    }
    Ok(())
}

// ─── -m type ─────────────────────────────────────────────────────────────────

fn run_type(
    traj: &ferro_core::Trajectory,
    params: &TypeParams,
    output: Option<&Path>,
    frame_arg: Option<usize>,
) -> Result<()> {
    if traj.frames.is_empty() {
        bail!("Trajectory is empty");
    }

    let frame_idx = frame_arg.unwrap_or(traj.frames.len() - 1);
    if frame_idx >= traj.frames.len() {
        bail!("Frame index {frame_idx} out of range (trajectory has {} frames)", traj.frames.len());
    }
    let frame = &traj.frames[frame_idx];
    let cell = frame.cell.as_ref()
        .ok_or_else(|| anyhow::anyhow!("Frame {frame_idx} has no cell (PBC required)"))?;

    let labels = ferro_core::classify_frame(frame, cell, params);
    print_type_table(&labels, frame_idx, traj.frames.len());

    if let Some(out) = output {
        let typed_frame = apply_type_labels(frame, &labels);
        let typed_traj = ferro_core::Trajectory {
            frames: vec![typed_frame],
            metadata: traj.metadata.clone(),
        };
        let out_str = out.to_str()
            .ok_or_else(|| anyhow::anyhow!("Output path is not valid UTF-8"))?;
        write_xyz(&typed_traj, out_str)?;
        println!("Structure → {}", out.display());
    }
    Ok(())
}

fn print_type_table(labels: &[String], frame_idx: usize, total_frames: usize) {
    use std::collections::HashMap;
    let mut counts: HashMap<&str, usize> = HashMap::new();
    for lbl in labels { *counts.entry(lbl.as_str()).or_insert(0) += 1; }

    let total = labels.len();
    println!("Frame {frame_idx}/{total_frames}  |  {total} atoms");
    println!("{:<18} {:>8} {:>10}", "Type", "Count", "Fraction");
    println!("{}", "-".repeat(38));

    // 排序：形成子先（Pn、Aln…），然后氧（Of→On→Ob→X），然后修饰子，然后其他
    let mut types: Vec<(&str, usize)> = counts.into_iter().collect();
    types.sort_by(|a, b| {
        type_sort_key(a.0).cmp(&type_sort_key(b.0)).then(a.0.cmp(b.0))
    });
    for (lbl, cnt) in &types {
        println!("{:<18} {:>8} {:>10.4}", lbl, cnt, *cnt as f64 / total as f64);
    }
}

/// 排序键：0=形成子(Qn), 1=Of, 2=On, 3=Ob, 4=修饰子(_f/_t/_b), 5=X, 6=其他
fn type_sort_key(label: &str) -> u8 {
    if label == "Of"               { return 1; }
    if label.starts_with("On_")    { return 2; }
    if label.starts_with("Ob_")    { return 3; }
    if label == "X"                { return 5; }
    if label.ends_with("_f") || label.ends_with("_t") || label.ends_with("_b") { return 4; }
    // 形成子：以字母开头，接数字
    if label.chars().next().map(|c| c.is_ascii_uppercase()).unwrap_or(false)
        && label.chars().last().map(|c| c.is_ascii_digit()).unwrap_or(false)
    { return 0; }
    6
}

// ─── CSV 输出（-m Qn）────────────────────────────────────────────────────────

fn write_qn_csv(result: &NetworkResult, params: &TypeParams, base: &Path) -> Result<()> {
    let stem = base.file_stem().and_then(|s| s.to_str()).unwrap_or("network");
    let dir  = base.parent().map(|p| p.to_str().unwrap_or("")).unwrap_or("");
    let prefix = if dir.is_empty() { stem.to_string() } else { format!("{dir}/{stem}") };

    write_qn_table(result, &format!("{prefix}_qn.csv"))?;
    write_oxy_table(result, &format!("{prefix}_oxy.csv"))?;
    write_cn_table(result, &format!("{prefix}_cn.csv"))?;
    println!("Qn      → {prefix}_qn.csv");
    println!("Oxygen  → {prefix}_oxy.csv");
    println!("CN      → {prefix}_cn.csv");

    if !result.modifier_dist.is_empty() {
        write_modifier_table(result, &format!("{prefix}_modifier.csv"))?;
        println!("Modifier→ {prefix}_modifier.csv");
    }
    let _ = params;
    Ok(())
}

fn write_qn_table(result: &NetworkResult, path: &str) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "Former,Qn,Count,Fraction,MeanQn")?;
    let mut formers: Vec<&String> = result.qn_dist.keys().collect();
    formers.sort();
    for former in formers {
        let rows = &result.qn_dist[former];
        let mean = result.mean_qn.get(former).copied().unwrap_or(0.0);
        for (i, &(qn, cnt, frac)) in rows.iter().enumerate() {
            if i == 0 {
                writeln!(f, "{former},{qn},{cnt},{frac:.6},{mean:.4}")?;
            } else {
                writeln!(f, "{former},{qn},{cnt},{frac:.6},")?;
            }
        }
    }
    Ok(())
}

fn write_oxy_table(result: &NetworkResult, path: &str) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "Type,Count,Fraction")?;
    for (lbl, cnt, frac) in &result.oxy_dist {
        writeln!(f, "{lbl},{cnt},{frac:.6}")?;
    }
    Ok(())
}

fn write_cn_table(result: &NetworkResult, path: &str) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "Former,CN,Count,Fraction,MeanCN")?;
    let mut formers: Vec<&String> = result.cn_dist.keys().collect();
    formers.sort();
    for former in formers {
        let rows = &result.cn_dist[former];
        let mean = result.mean_cn.get(former).copied().unwrap_or(0.0);
        for (i, &(cn, cnt, frac)) in rows.iter().enumerate() {
            if i == 0 {
                writeln!(f, "{former},{cn},{cnt},{frac:.6},{mean:.4}")?;
            } else {
                writeln!(f, "{former},{cn},{cnt},{frac:.6},")?;
            }
        }
    }
    Ok(())
}

fn write_modifier_table(result: &NetworkResult, path: &str) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "Modifier,Role,Count,Fraction")?;
    let mut mods: Vec<&String> = result.modifier_dist.keys().collect();
    mods.sort();
    for mod_elem in mods {
        for (lbl, cnt, frac) in &result.modifier_dist[mod_elem] {
            writeln!(f, "{mod_elem},{lbl},{cnt},{frac:.6}")?;
        }
    }
    Ok(())
}

// ─── XLSX 输出（-m Qn）───────────────────────────────────────────────────────

fn write_qn_xlsx(result: &NetworkResult, params: &TypeParams, base: &Path) -> Result<()> {
    use rust_xlsxwriter::*;

    let stem = base.file_stem().and_then(|s| s.to_str()).unwrap_or("network");
    let dir  = base.parent().map(|p| p.to_str().unwrap_or("")).unwrap_or("");
    let path = if dir.is_empty() { format!("{stem}.xlsx") } else { format!("{dir}/{stem}.xlsx") };

    let mut wb = Workbook::new();

    // Sheet: Qn
    {
        let ws = wb.add_worksheet(); ws.set_name("Qn")?;
        ws.write_row(0, 0, ["Former", "Qn", "Count", "Fraction", "MeanQn"])?;
        let mut row = 1u32;
        let mut formers: Vec<&String> = result.qn_dist.keys().collect(); formers.sort();
        for former in formers {
            let rows = &result.qn_dist[former];
            let mean = result.mean_qn.get(former).copied().unwrap_or(0.0);
            for (i, &(qn, cnt, frac)) in rows.iter().enumerate() {
                ws.write_row(row, 0, [former.as_str()])?;
                ws.write(row, 1, qn)?; ws.write(row, 2, cnt as u32)?; ws.write(row, 3, frac)?;
                if i == 0 { ws.write(row, 4, mean)?; }
                row += 1;
            }
        }
    }

    // Sheet: Oxygen
    {
        let ws = wb.add_worksheet(); ws.set_name("Oxygen")?;
        ws.write_row(0, 0, ["Type", "Count", "Fraction"])?;
        for (row, (lbl, cnt, frac)) in result.oxy_dist.iter().enumerate() {
            let row = (row + 1) as u32;
            ws.write(row, 0, lbl.as_str())?;
            ws.write(row, 1, *cnt as u32)?; ws.write(row, 2, *frac)?;
        }
    }

    // Sheet: CN
    {
        let ws = wb.add_worksheet(); ws.set_name("CN")?;
        ws.write_row(0, 0, ["Former", "CN", "Count", "Fraction", "MeanCN"])?;
        let mut row = 1u32;
        let mut formers: Vec<&String> = result.cn_dist.keys().collect(); formers.sort();
        for former in formers {
            let rows = &result.cn_dist[former];
            let mean = result.mean_cn.get(former).copied().unwrap_or(0.0);
            for (i, &(cn, cnt, frac)) in rows.iter().enumerate() {
                ws.write_row(row, 0, [former.as_str()])?;
                ws.write(row, 1, cn)?; ws.write(row, 2, cnt as u32)?; ws.write(row, 3, frac)?;
                if i == 0 { ws.write(row, 4, mean)?; }
                row += 1;
            }
        }
    }

    // Sheet: Modifier (optional)
    if !result.modifier_dist.is_empty() {
        let ws = wb.add_worksheet(); ws.set_name("Modifier")?;
        ws.write_row(0, 0, ["Modifier", "Role", "Count", "Fraction"])?;
        let mut row = 1u32;
        let mut mods: Vec<&String> = result.modifier_dist.keys().collect(); mods.sort();
        for mod_elem in mods {
            for (lbl, cnt, frac) in &result.modifier_dist[mod_elem] {
                ws.write_row(row, 0, [mod_elem.as_str(), lbl.as_str()])?;
                ws.write(row, 2, *cnt as u32)?; ws.write(row, 3, *frac)?;
                row += 1;
            }
        }
    }

    let sheets = if result.modifier_dist.is_empty() { "Qn, Oxygen, CN" }
                 else { "Qn, Oxygen, CN, Modifier" };
    wb.save(&path)?;
    println!("Network → {path}  (sheets: {sheets})");
    let _ = params;
    Ok(())
}

// ─── Pair 参数解析 ────────────────────────────────────────────────────────────

fn split_pair_args(all: &[String]) -> (Vec<String>, Vec<String>) {
    let mut pairs = Vec::new();
    let mut clap  = Vec::new();
    for arg in all {
        if is_pair_arg(arg) { pairs.push(arg.clone()); } else { clap.push(arg.clone()); }
    }
    (pairs, clap)
}

fn is_pair_arg(s: &str) -> bool {
    if !s.starts_with("--") { return false; }
    let inner = &s[2..];
    inner.starts_with(|c: char| c.is_ascii_uppercase()) && inner.contains('=')
}

fn parse_pairs(pair_args: &[String]) -> Result<BTreeMap<(String, String), f64>> {
    let mut map = BTreeMap::new();
    for arg in pair_args {
        let inner = arg.trim_start_matches('-');
        let (pair, cutoff_str) = inner.split_once('=')
            .ok_or_else(|| anyhow::anyhow!("Invalid pair argument (missing '='): {arg}"))?;
        let (former, ligand) = pair.split_once('-')
            .ok_or_else(|| anyhow::anyhow!("Invalid pair argument (missing '-'): {arg}"))?;
        let cutoff: f64 = cutoff_str.parse()
            .map_err(|_| anyhow::anyhow!("Invalid cutoff value in '{arg}'"))?;
        if cutoff <= 0.0 { bail!("Cutoff must be positive, got {cutoff} in '{arg}'"); }
        map.insert((former.to_string(), ligand.to_string()), cutoff);
    }
    Ok(map)
}
