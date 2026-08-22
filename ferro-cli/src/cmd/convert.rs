use anyhow::{bail, Result};
use clap::Args;
use crate::io_dispatch::{read_trajectory, write_trajectory};
use ferro_io::LammpsUnits;
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct ConvertCmd {
    /// Input file (format auto-detected; omit to list the supported formats)
    #[arg(short, long)]
    pub input: Option<PathBuf>,

    /// Output file (format auto-detected)
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Use LAMMPS metal units for dump files (velocities Å/ps, forces eV/Å)
    #[arg(long)]
    pub metal_units: bool,
}

/// True when `ferro convert` was typed with no input: print the format matrix
/// rather than clap's bare "required argument" error.
pub fn wants_help(args: &ConvertCmd) -> bool {
    args.input.is_none()
}

pub fn run(args: &ConvertCmd) -> Result<()> {
    // input 为空由 main 分派到帮助页，这里的 bail 只是防御
    let Some(input) = &args.input else {
        bail!("convert needs an input file: -i <FILE>");
    };
    let Some(output) = &args.output else {
        bail!("convert needs an output file: -o <FILE>  (run `ferro convert` for the format list)");
    };

    if input == output {
        bail!("Input and output paths are identical");
    }

    let units = if args.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let traj = read_trajectory(input, units)?;
    let n = traj.frames.len();
    write_trajectory(&traj, output, units)?;

    println!(
        "Converted {} ({} frame{}) -> {}",
        input.display(),
        n,
        if n == 1 { "" } else { "s" },
        output.display()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cmd(input: Option<&str>, output: Option<&str>) -> ConvertCmd {
        ConvertCmd {
            input: input.map(PathBuf::from),
            output: output.map(PathBuf::from),
            metal_units: false,
        }
    }

    #[test]
    fn test_wants_help_only_without_input() {
        assert!(wants_help(&cmd(None, None)));
        assert!(wants_help(&cmd(None, Some("out.xyz"))));
        assert!(!wants_help(&cmd(Some("in.xyz"), None)));
    }

    /// 有输入却没有输出是**用法错误**，不是求助：报错而不是把帮助页当成功输出，
    /// 否则 shell 里 `ferro convert -i a.xyz && next` 会把漏参数当成转换完成
    #[test]
    fn test_missing_output_is_an_error_not_a_help_page() {
        let err = run(&cmd(Some("in.xyz"), None)).unwrap_err();
        assert!(err.to_string().contains("-o"), "{err}");
    }

    #[test]
    fn test_identical_paths_are_rejected() {
        let err = run(&cmd(Some("same.xyz"), Some("same.xyz"))).unwrap_err();
        assert!(err.to_string().contains("identical"), "{err}");
    }
}
