use anyhow::{bail, Result};
use clap::Args;
use ferro_io::{read_chgcar, read_cube_as_chg};
use ferro_analysis::{BaderAnalyzer, BaderMethod};
use std::path::{Path, PathBuf};

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

    /// Directory for the three reports (default: next to the input file)
    #[arg(short, long, value_name = "DIR")]
    pub output: Option<PathBuf>,

    /// Tag appended to the report names: <stem>_ACF_<suffix>.dat
    #[arg(short, long, value_name = "SUFFIX")]
    pub suffix: Option<String>,

    /// Create the output directory without asking (required when there is no terminal)
    #[arg(long)]
    pub mkdir: bool,
}

/// True when `ferro bader` was typed with no input: print the help page rather
/// than clap's bare "required argument" error.
pub fn wants_help(args: &BaderCmd) -> bool {
    args.input.is_none()
}

/// Where the three reports go: `<-o dir>/<input stem>_ACF[_<suffix>].dat`.
///
/// `-o` defaults to the **input's own directory**, the one exception to "default is the
/// current directory". VASP names every charge density `CHGCAR`, so running `run1/CHGCAR`
/// and `run2/CHGCAR` from one shell used to write `CHGCAR_ACF.dat` twice, the second
/// silently over the first. Next to the input, each run keeps its own reports.
fn report_paths(args: &BaderCmd, input: &Path) -> [PathBuf; 3] {
    let stem = input.file_stem().and_then(|s| s.to_str()).unwrap_or("bader");
    let tag = match args.suffix.as_deref().filter(|s| !s.is_empty()) {
        Some(s) => format!("_{s}"),
        None => String::new(),
    };
    let dir = args
        .output
        .clone()
        .unwrap_or_else(|| input.parent().unwrap_or(Path::new("")).to_path_buf());
    ["ACF", "BCF", "AVF"].map(|kind| dir.join(format!("{stem}_{kind}{tag}.dat")))
}

pub fn run(args: &BaderCmd) -> Result<()> {
    // input 为空由 main 分派到帮助页，这里的 bail 只是防御
    let Some(input) = &args.input else {
        bail!("bader needs an input file: -i <CHGCAR|FILE.cube>");
    };
    // 参数校验在读文件之前：CHGCAR 动辄几百 MB，读完再报「方法名打错了」是白等
    let method = match args.method.to_lowercase().as_str() {
        "ongrid"   => BaderMethod::OnGrid,
        "neargrid" => BaderMethod::NearGrid,
        "offgrid"  => BaderMethod::OffGrid,
        "weight"   => BaderMethod::Weight,
        other => bail!("Unknown method: {other}.  Use ongrid|neargrid|offgrid|weight"),
    };

    // 目录在读 CHGCAR 之前建好:几百 MB 读完才发现路径打不开是白等
    if let Some(dir) = &args.output {
        crate::outpath::ensure_dir(dir, args.mkdir)?;
    }

    let is_cube = input
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("cube"))
        .unwrap_or(false);

    let (frame, chg) = if is_cube {
        read_cube_as_chg(input)?
    } else {
        read_chgcar(input)?
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

    let [acf_path, bcf_path, avf_path] = report_paths(args, input);
    std::fs::write(&acf_path, result.acf_text(&frame))?;
    std::fs::write(&bcf_path, result.bcf_text())?;
    std::fs::write(&avf_path, result.avf_text())?;

    println!(
        "Output: {}, {}, {}",
        acf_path.display(), bcf_path.display(), avf_path.display()
    );

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
            output: None,
            suffix: None,
            mkdir: false,
        }
    }

    #[test]
    fn test_wants_help_only_without_input() {
        assert!(wants_help(&cmd(None)));
        assert!(!wants_help(&cmd(Some("CHGCAR"))));
    }

    /// 产物落在输入旁边:两个体系的 CHGCAR 同名,写当前目录会互相覆盖。
    #[test]
    fn test_reports_land_next_to_the_input() {
        let c = cmd(Some("run1/CHGCAR"));
        let [acf, bcf, avf] = report_paths(&c, Path::new("run1/CHGCAR"));
        assert_eq!(acf, PathBuf::from("run1/CHGCAR_ACF.dat"));
        assert_eq!(bcf, PathBuf::from("run1/CHGCAR_BCF.dat"));
        assert_eq!(avf, PathBuf::from("run1/CHGCAR_AVF.dat"));
    }

    /// -o 覆盖目录,-s 追加标记 —— 同一个输入跑两套参数时用得上。
    #[test]
    fn test_output_dir_and_suffix_apply() {
        let mut c = cmd(Some("run1/CHGCAR"));
        c.output = Some(PathBuf::from("reports"));
        c.suffix = Some("weight".into());
        let [acf, ..] = report_paths(&c, Path::new("run1/CHGCAR"));
        assert_eq!(acf, PathBuf::from("reports/CHGCAR_ACF_weight.dat"));
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
