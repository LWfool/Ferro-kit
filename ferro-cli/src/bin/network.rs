use anyhow::{bail, Result};
use clap::{CommandFactory, Parser};
use ferro::io_dispatch::read_trajectory;
use ferro_analysis::{calc_network, NetworkParams, NetworkResult};
use ferro_io::LammpsUnits;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(
    name = "fe-network",
    about = "Glass network analysis: CN, ligand classification (FO/NBO/BO/OBO), Qn speciation, modifier roles",
    long_about = None,
    after_help = HELP_EXTRA,
)]
struct Cli {
    /// Input trajectory file (omit to show help with examples)
    #[arg(short, long)]
    input: Option<PathBuf>,

    /// Output file (default: network.csv)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output format: csv | xlsx
    #[arg(long, default_value = "csv")]
    format: String,

    /// Use only the last N frames
    #[arg(long)]
    last_n: Option<usize>,

    /// Parallel threads (default: all cores)
    #[arg(long)]
    ncore: Option<usize>,

    /// Use LAMMPS metal units for dump files (velocities Å/ps, forces eV/Å)
    #[arg(long)]
    metal_units: bool,

    /// Modifier cation elements, comma-separated (e.g. Zn or Zn,Na,Ca).
    /// Their cutoffs are supplied via the same --Elem-O=cutoff syntax.
    #[arg(long)]
    modifier: Option<String>,
}

const HELP_EXTRA: &str = "\
PAIR ARGUMENTS:
  Specify cutoff radii using --Former-Ligand=<cutoff> flags:
    --P-O=2.3         P as network former, O as ligand, cutoff 2.3 Å
    --Si-O=1.8        Si as former, O as ligand, cutoff 1.8 Å
    --Al-O=2.0 --Al-F=2.1   multiple pairs for the same former

  Network formers (examples): Si, Al, P, B, Ge, Zn
  Ligands        (examples): O, F

MODIFIER ROLE ANALYSIS:
  Use --modifier to name modifier cation elements (comma-separated).
  Supply their cutoffs with the same --Elem-O=cutoff syntax.
  Modifier role labels (by NBO coordination count):
    Free  — 0 NBO neighbors
    T     — 1 NBO neighbor  (terminal)
    B     — 2 NBO neighbors (bridging)
    M     — ≥3 NBO neighbors (multi)

EXAMPLES:
  fe-network -i traj.dump --P-O=2.3
  fe-network -i traj.dump --P-O=2.3 --P-F=2.1 --Si-O=1.8 -o result.csv
  fe-network -i traj.dump --P-O=2.3 --format xlsx -o result.xlsx
  fe-network -i traj.dump --P-O=2.3 --last-n 500
  fe-network -i traj.dump --P-O=2.3 --Zn-O=3.5 --modifier Zn
  fe-network -i traj.dump --P-O=2.3 --Zn-O=3.5 --Na-O=3.2 --modifier Zn,Na

OUTPUT FILES:
  <stem>_cn.csv         per-pair and total CN distribution
  <stem>_ligand.csv     FO / NBO / BO / OBO classification
  <stem>_qn.csv         Qn species distribution
  <stem>_modifier.csv   modifier role distribution (Free/T/B/M)  [only with --modifier]";

fn main() -> Result<()> {
    // ── 预解析 --X-Y=N 形式的 pair 参数 ──────────────────────────────────────
    let all_args: Vec<String> = std::env::args().collect();
    let (pair_args, clap_args) = split_pair_args(&all_args);

    let cli = Cli::parse_from(clap_args);

    // 无输入文件时显示帮助
    let input = match cli.input {
        Some(ref p) => p.clone(),
        None => {
            println!("{}", Cli::command().render_long_help());
            return Ok(());
        }
    };

    if pair_args.is_empty() {
        bail!(
            "No pair cutoffs specified. Use --Former-Ligand=cutoff, e.g. --P-O=2.3\n\
             Run without -i for usage examples."
        );
    }

    let all_pairs = parse_pairs(&pair_args)?;

    // 修饰子集合（逗号分隔，如 "Zn,Na"）
    let modifier_elems: std::collections::HashSet<String> = cli.modifier
        .as_deref()
        .map(|s| s.split(',').map(|e| e.trim().to_string()).collect())
        .unwrap_or_default();

    // 按修饰子集合将 pair args 路由到 former cutoffs / modifier cutoffs
    let mut cutoffs = BTreeMap::new();
    let mut modifier_cutoffs = BTreeMap::new();
    for ((elem, ligand), cutoff) in all_pairs {
        if modifier_elems.contains(&elem) {
            modifier_cutoffs.insert((elem, ligand), cutoff);
        } else {
            cutoffs.insert((elem, ligand), cutoff);
        }
    }

    if let Some(n) = cli.ncore {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().ok();
    }

    let units = if cli.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let mut traj = read_trajectory(&input, units)?;
    if let Some(n) = cli.last_n {
        traj = traj.tail(n);
    }

    let params = NetworkParams { cutoffs, modifier_cutoffs };
    let result = calc_network(&traj, &params)
        .ok_or_else(|| anyhow::anyhow!(
            "Network analysis failed — trajectory empty or frames missing cell (PBC required)"
        ))?;

    let out_base = cli.output.clone().unwrap_or_else(|| PathBuf::from("network"));
    let fmt = cli.format.to_lowercase();

    match fmt.as_str() {
        "csv"  => write_csv(&result, &params, &out_base)?,
        "xlsx" => write_xlsx(&result, &params, &out_base)?,
        other  => bail!("Unknown format '{other}'. Use csv or xlsx."),
    }

    Ok(())
}

// ─── Pair 参数解析 ────────────────────────────────────────────────────────────

/// 将 CLI args 分为 pair 参数（--X-Y=N）和普通参数（交给 clap）。
///
/// 判别规则：去掉 `--` 前缀后，第一个字符是大写字母 → pair 参数。
fn split_pair_args(all: &[String]) -> (Vec<String>, Vec<String>) {
    let mut pairs = Vec::new();
    let mut clap = Vec::new();
    for arg in all {
        if is_pair_arg(arg) {
            pairs.push(arg.clone());
        } else {
            clap.push(arg.clone());
        }
    }
    (pairs, clap)
}

fn is_pair_arg(s: &str) -> bool {
    if !s.starts_with("--") { return false; }
    let inner = &s[2..];
    // 首字符大写 + 包含 '=' → 是 pair 参数
    inner.starts_with(|c: char| c.is_ascii_uppercase()) && inner.contains('=')
}

/// 解析 ["--P-O=2.3", "--Si-O=1.8"] → CutoffTable
fn parse_pairs(pair_args: &[String]) -> Result<BTreeMap<(String, String), f64>> {
    let mut map = BTreeMap::new();
    for arg in pair_args {
        let inner = arg.trim_start_matches('-');
        // inner = "P-O=2.3"
        let (pair, cutoff_str) = inner.split_once('=')
            .ok_or_else(|| anyhow::anyhow!("Invalid pair argument (missing '='): {arg}"))?;
        // pair = "P-O" → split on first '-'
        let (former, ligand) = pair.split_once('-')
            .ok_or_else(|| anyhow::anyhow!("Invalid pair argument (missing '-' between elements): {arg}"))?;
        let cutoff: f64 = cutoff_str.parse()
            .map_err(|_| anyhow::anyhow!("Invalid cutoff value in '{arg}'"))?;
        if cutoff <= 0.0 {
            bail!("Cutoff must be positive, got {cutoff} in '{arg}'");
        }
        map.insert((former.to_string(), ligand.to_string()), cutoff);
    }
    Ok(map)
}

// ─── CSV 输出 ─────────────────────────────────────────────────────────────────

fn write_csv(result: &NetworkResult, params: &NetworkParams, base: &Path) -> Result<()> {
    let stem = base.file_stem().and_then(|s| s.to_str()).unwrap_or("network");
    let dir  = base.parent().map(|p| p.to_str().unwrap_or("")).unwrap_or("");
    let prefix = if dir.is_empty() { stem.to_string() } else { format!("{dir}/{stem}") };

    write_cn_csv(result, params, &format!("{prefix}_cn.csv"))?;
    write_ligand_csv(result, &format!("{prefix}_ligand.csv"))?;
    write_qn_csv(result, &format!("{prefix}_qn.csv"))?;

    println!("CN       -> {prefix}_cn.csv");
    println!("Ligand   -> {prefix}_ligand.csv");
    println!("Qn       -> {prefix}_qn.csv");

    if !result.modifier_roles.is_empty() {
        let path = format!("{prefix}_modifier.csv");
        write_modifier_csv(result, &path)?;
        println!("Modifier -> {path}");
    }
    Ok(())
}

fn write_cn_csv(result: &NetworkResult, params: &NetworkParams, path: &str) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);

    writeln!(f, "Former,Ligand,CN,Count,Fraction")?;

    let mut pairs: Vec<&(String, String)> = result.cn_dist.keys().collect();
    pairs.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    for pair in &pairs {
        if let Some(rows) = result.cn_dist.get(*pair) {
            for &(cn, count, frac) in rows {
                writeln!(f, "{},{},{},{},{:.6}", pair.0, pair.1, cn, count, frac)?;
            }
        }
        // mean line
        if let Some(&mean) = result.mean_cn.get(*pair) {
            writeln!(f, "{},{},mean,{:.4},", pair.0, pair.1, mean)?;
        }
    }

    // Total CN per former
    let mut formers: Vec<&String> = result.cn_total.keys().collect();
    formers.sort();
    for former in formers {
        if let Some(rows) = result.cn_total.get(former) {
            for &(cn, count, frac) in rows {
                writeln!(f, "{},total,{},{},{:.6}", former, cn, count, frac)?;
            }
        }
        if let Some(&mean) = result.mean_cn_total.get(former) {
            writeln!(f, "{},total,mean,{:.4},", former, mean)?;
        }
    }

    // Summary header
    writeln!(f)?;
    writeln!(f, "# Mean CN summary")?;
    for pair in &pairs {
        if let Some(&mean) = result.mean_cn.get(*pair) {
            writeln!(f, "# {}-{}: {:.4}", pair.0, pair.1, mean)?;
        }
    }
    let _ = params;
    Ok(())
}

fn write_ligand_csv(result: &NetworkResult, path: &str) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "Ligand,Class,Count,Fraction")?;

    let mut ligs: Vec<&String> = result.ligand_classes.keys().collect();
    ligs.sort();
    for lig in ligs {
        if let Some(rows) = result.ligand_classes.get(lig) {
            for (label, count, frac) in rows {
                writeln!(f, "{lig},{label},{count},{frac:.6}")?;
            }
        }
    }
    Ok(())
}

fn write_qn_csv(result: &NetworkResult, path: &str) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "Former,Species,Count,Fraction")?;

    let mut formers: Vec<&String> = result.qn_species.keys().collect();
    formers.sort();
    for former in formers {
        if let Some(rows) = result.qn_species.get(former) {
            for (label, count, frac) in rows {
                writeln!(f, "{former},{label},{count},{frac:.6}")?;
            }
        }
    }
    Ok(())
}

fn write_modifier_csv(result: &NetworkResult, path: &str) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "Modifier,Role,Count,Fraction")?;

    let mut mods: Vec<&String> = result.modifier_roles.keys().collect();
    mods.sort();
    for mod_elem in mods {
        if let Some(rows) = result.modifier_roles.get(mod_elem) {
            for (label, count, frac) in rows {
                writeln!(f, "{mod_elem},{label},{count},{frac:.6}")?;
            }
        }
    }
    Ok(())
}

// ─── XLSX 输出 ────────────────────────────────────────────────────────────────

fn write_xlsx(result: &NetworkResult, params: &NetworkParams, base: &Path) -> Result<()> {
    use rust_xlsxwriter::*;

    let stem = base.file_stem().and_then(|s| s.to_str()).unwrap_or("network");
    let dir  = base.parent().map(|p| p.to_str().unwrap_or("")).unwrap_or("");
    let path = if dir.is_empty() {
        format!("{stem}.xlsx")
    } else {
        format!("{dir}/{stem}.xlsx")
    };

    let mut wb = Workbook::new();

    // ── Sheet 1: CN ──────────────────────────────────────────────────────────
    {
        let ws = wb.add_worksheet();
        ws.set_name("CN")?;
        ws.write_row(0, 0, ["Former", "Ligand", "CN", "Count", "Fraction"])?;
        let mut row = 1u32;

        let mut pairs: Vec<&(String, String)> = result.cn_dist.keys().collect();
        pairs.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        for pair in &pairs {
            if let Some(rows) = result.cn_dist.get(*pair) {
                for &(cn, count, frac) in rows {
                    ws.write_row(row, 0, [pair.0.as_str(), pair.1.as_str()])?;
                    ws.write(row, 2, cn)?;
                    ws.write(row, 3, count as u32)?;
                    ws.write(row, 4, frac)?;
                    row += 1;
                }
            }
            if let Some(&mean) = result.mean_cn.get(*pair) {
                ws.write_row(row, 0, [pair.0.as_str(), pair.1.as_str(), "mean"])?;
                ws.write(row, 4, mean)?;
                row += 1;
            }
        }
        // Total CN
        let mut formers: Vec<&String> = result.cn_total.keys().collect();
        formers.sort();
        for former in formers {
            if let Some(rows) = result.cn_total.get(former) {
                for &(cn, count, frac) in rows {
                    ws.write_row(row, 0, [former.as_str(), "total"])?;
                    ws.write(row, 2, cn)?;
                    ws.write(row, 3, count as u32)?;
                    ws.write(row, 4, frac)?;
                    row += 1;
                }
            }
            if let Some(&mean) = result.mean_cn_total.get(former) {
                ws.write_row(row, 0, [former.as_str(), "total", "mean"])?;
                ws.write(row, 4, mean)?;
                row += 1;
            }
        }
        let _ = params;
    }

    // ── Sheet 2: Ligand Classification ───────────────────────────────────────
    {
        let ws = wb.add_worksheet();
        ws.set_name("Ligand")?;
        ws.write_row(0, 0, ["Ligand", "Class", "Count", "Fraction"])?;
        let mut row = 1u32;

        let mut ligs: Vec<&String> = result.ligand_classes.keys().collect();
        ligs.sort();
        for lig in ligs {
            if let Some(rows) = result.ligand_classes.get(lig) {
                for (label, count, frac) in rows {
                    ws.write_row(row, 0, [lig.as_str(), label.as_str()])?;
                    ws.write(row, 2, *count as u32)?;
                    ws.write(row, 3, *frac)?;
                    row += 1;
                }
            }
        }
    }

    // ── Sheet 3: Qn Species ──────────────────────────────────────────────────
    {
        let ws = wb.add_worksheet();
        ws.set_name("Qn")?;
        ws.write_row(0, 0, ["Former", "Species", "Count", "Fraction"])?;
        let mut row = 1u32;

        let mut formers: Vec<&String> = result.qn_species.keys().collect();
        formers.sort();
        for former in formers {
            if let Some(rows) = result.qn_species.get(former) {
                for (label, count, frac) in rows {
                    ws.write_row(row, 0, [former.as_str(), label.as_str()])?;
                    ws.write(row, 2, *count as u32)?;
                    ws.write(row, 3, *frac)?;
                    row += 1;
                }
            }
        }
    }

    // ── Sheet 4: Modifier Roles (optional) ───────────────────────────────────
    if !result.modifier_roles.is_empty() {
        let ws = wb.add_worksheet();
        ws.set_name("Modifier")?;
        ws.write_row(0, 0, ["Modifier", "Role", "Count", "Fraction"])?;
        let mut row = 1u32;

        let mut mods: Vec<&String> = result.modifier_roles.keys().collect();
        mods.sort();
        for mod_elem in mods {
            if let Some(rows) = result.modifier_roles.get(mod_elem) {
                for (label, count, frac) in rows {
                    ws.write_row(row, 0, [mod_elem.as_str(), label.as_str()])?;
                    ws.write(row, 2, *count as u32)?;
                    ws.write(row, 3, *frac)?;
                    row += 1;
                }
            }
        }
    }

    let sheet_list = if result.modifier_roles.is_empty() {
        "CN, Ligand, Qn"
    } else {
        "CN, Ligand, Qn, Modifier"
    };
    wb.save(&path)?;
    println!("Network  -> {path}  (sheets: {sheet_list})");
    Ok(())
}
