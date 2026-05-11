use anyhow::{bail, Result};
use clap::Parser;
use ferro_io::read_chgcar;
use ferro_analysis::{BaderAnalyzer, BaderMethod};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "fe-bader", about = "Bader charge analysis from VASP CHGCAR")]
struct Cli {
    /// Input CHGCAR file
    #[arg(short, long)]
    input: PathBuf,

    /// Bader method: ongrid | neargrid | offgrid | weight
    #[arg(short, long, default_value = "neargrid")]
    method: String,

    /// Edge refinement: -1 = auto, -2 = single pass, N = N passes
    #[arg(short, long, default_value_t = -1)]
    refine: i32,

    /// Vacuum density threshold (e/Å³)
    #[arg(short, long, default_value_t = 1e-3)]
    vacval: f64,
}

fn main() -> Result<()> {
    let args = Cli::parse();
    let path = args.input.to_str().unwrap_or_default();

    let (frame, chg) = read_chgcar(path)?;

    let method = match args.method.to_lowercase().as_str() {
        "ongrid"   => BaderMethod::OnGrid,
        "neargrid" => BaderMethod::NearGrid,
        "offgrid"  => BaderMethod::OffGrid,
        "weight"   => BaderMethod::Weight,
        other => bail!("Unknown method: {other}.  Use ongrid|neargrid|offgrid|weight"),
    };

    println!("Bader analysis: method={:?}, refine={}, vacval={:.1e}", method, args.refine, args.vacval);
    println!("Grid: {} × {} × {} ({} points)", chg.shape[0], chg.shape[1], chg.shape[2], chg.nrho);
    println!("Atoms: {}", frame.n_atoms());

    let result = BaderAnalyzer::new(chg, frame)
        .method(method)
        .refine(args.refine)
        .vacval(args.vacval)
        .run();

    println!("Bader volumes found: {}", result.nvols);

    // Write output files
    let stem = args.input.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("bader");
    let acf_path = format!("{stem}_ACF.dat");
    let bcf_path = format!("{stem}_BCF.dat");
    let avf_path = format!("{stem}_AVF.dat");

    result.write_acf(&acf_path)?;
    result.write_bcf(&bcf_path)?;
    result.write_avf(&avf_path)?;

    println!("Output: {acf_path}, {bcf_path}, {avf_path}");

    // Summary
    let total_e: f64 = result.ionchg.iter().sum();
    println!("\nTotal ionic charge: {:.4} e", total_e);
    println!("Vacuum charge:      {:.4} e", result.vacchg);
    println!("Total:              {:.4} e", total_e + result.vacchg);

    Ok(())
}
