use anyhow::{bail, Result};
use clap::Args;
use ferro_io::{read_chgcar, read_cube_as_chg};
use ferro_analysis::{BaderAnalyzer, BaderMethod};
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct BaderCmd {
    /// Input file: CHGCAR (VASP) or .cube (Gaussian / QE pp.x)
    /// (omit to show the methods and output files)
    #[arg(short, long)]
    pub input: Option<PathBuf>,

    /// Bader method: ongrid | neargrid | offgrid | weight
    #[arg(short, long, default_value = "neargrid")]
    pub method: String,

    /// Edge refinement: -1 = auto, -2 = single pass, N = N passes
    #[arg(short, long, default_value_t = -1)]
    pub refine: i32,

    /// Vacuum density threshold (e/Å³)
    #[arg(short, long, default_value_t = 1e-3)]
    pub vacval: f64,
}

/// True when `ferro bader` was typed with no input: print the help page rather
/// than clap's bare "required argument" error.
pub fn wants_help(args: &BaderCmd) -> bool {
    args.input.is_none()
}

pub fn run(args: &BaderCmd) -> Result<()> {
    // input 为空由 main 分派到帮助页，这里的 bail 只是防御
    let Some(input) = &args.input else {
        bail!("bader needs an input file: -i <CHGCAR|FILE.cube>");
    };
    let path = input.to_str().unwrap_or_default();

    // 参数校验在读文件之前：CHGCAR 动辄几百 MB，读完再报「方法名打错了」是白等
    let method = match args.method.to_lowercase().as_str() {
        "ongrid"   => BaderMethod::OnGrid,
        "neargrid" => BaderMethod::NearGrid,
        "offgrid"  => BaderMethod::OffGrid,
        "weight"   => BaderMethod::Weight,
        other => bail!("Unknown method: {other}.  Use ongrid|neargrid|offgrid|weight"),
    };

    let is_cube = input
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("cube"))
        .unwrap_or(false);

    let (frame, chg) = if is_cube {
        read_cube_as_chg(path)?
    } else {
        read_chgcar(path)?
    };

    println!("Bader analysis: method={:?}, refine={}, vacval={:.1e}", method, args.refine, args.vacval);
    println!("Grid: {} × {} × {} ({} points)", chg.shape[0], chg.shape[1], chg.shape[2], chg.nrho);
    println!("Atoms: {}", frame.n_atoms());

    let result = BaderAnalyzer::new(chg, frame.clone())
        .method(method)
        .refine(args.refine)
        .vacval(args.vacval)
        .run();

    println!("Bader volumes found: {}", result.nvols);

    let stem = input.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("bader");
    let acf_path = format!("{stem}_ACF.dat");
    let bcf_path = format!("{stem}_BCF.dat");
    let avf_path = format!("{stem}_AVF.dat");

    result.write_acf(&acf_path, &frame)?;
    result.write_bcf(&bcf_path)?;
    result.write_avf(&avf_path)?;

    println!("Output: {acf_path}, {bcf_path}, {avf_path}");

    let total_e: f64 = result.ionchg.iter().sum();
    println!("\nTotal ionic charge: {:.4} e", total_e);
    println!("Vacuum charge:      {:.4} e", result.vacchg);
    println!("Total:              {:.4} e", total_e + result.vacchg);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cmd(input: Option<&str>) -> BaderCmd {
        BaderCmd {
            input: input.map(PathBuf::from),
            method: "neargrid".into(),
            refine: -1,
            vacval: 1e-3,
        }
    }

    #[test]
    fn test_wants_help_only_without_input() {
        assert!(wants_help(&cmd(None)));
        assert!(!wants_help(&cmd(Some("CHGCAR"))));
    }

    /// 方法名在读文件**之前**校验：拿错方法名跑一个 GB 级 CHGCAR 再报错是白等
    #[test]
    fn test_unknown_method_is_rejected() {
        let mut c = cmd(Some("no_such_CHGCAR"));
        c.method = "nosuchmethod".into();
        let err = run(&c).unwrap_err();
        assert!(err.to_string().contains("Unknown method"), "{err}");
    }
}
