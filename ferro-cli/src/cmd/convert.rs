use anyhow::{bail, Context, Result};
use clap::Args;
use crate::io_dispatch::{holds_multiple_frames, read_trajectory, write_trajectory};
use ferro_core::Trajectory;
use ferro_io::LammpsUnits;
use std::path::{Path, PathBuf};

#[derive(Args, Debug)]
pub struct ConvertCmd {
    /// Input file (format auto-detected; omit to list the supported formats)
    #[arg(short, long)]
    pub input: Option<PathBuf>,

    /// Output file (format auto-detected)
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// First frame to take (0-based, inclusive)               [default: 0]
    #[arg(long, value_name = "N")]
    pub start: Option<usize>,

    /// Last frame to take (0-based, INCLUSIVE)     [default: the last frame]
    #[arg(long, value_name = "N")]
    pub end: Option<usize>,

    /// Take every Nth frame within [start, end]
    #[arg(long, value_name = "N")]
    pub stride: Option<usize>,

    /// Take this many frames, spread evenly over [start, end], both ends included
    #[arg(long, value_name = "N", conflicts_with = "stride")]
    pub number: Option<usize>,

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

    // 参数级错误在读第一个文件之前失败：轨迹可能是几个 GB，读完再报「--start 比
    // --end 大」是白等。跨参数的一致性能在这里查的都查掉
    let start = args.start.unwrap_or(0);
    if let (Some(s), Some(e)) = (args.start, args.end) {
        if s > e {
            bail!("--start {s} is past --end {e} (both are 0-based and inclusive)");
        }
    }
    if args.stride == Some(0) {
        bail!("--stride must be at least 1");
    }
    if args.number == Some(0) {
        bail!("--number must be at least 1");
    }

    let units = if args.metal_units { LammpsUnits::Metal } else { LammpsUnits::Real };
    let traj = read_trajectory(input, units)?;
    let n_read = traj.n_frames();

    let indices = match args.number {
        Some(count) => traj.spread_indices(start, args.end, count),
        None => traj.select_indices(start, args.end, args.stride.unwrap_or(1)),
    };
    if indices.is_empty() {
        bail!(
            "no frames selected: {} holds {n_read} frame(s) (0..={}), \
             but --start {start}{} selects none",
            input.display(),
            n_read.saturating_sub(1),
            args.end.map(|e| format!(" --end {e}")).unwrap_or_default()
        );
    }

    let picked = Trajectory {
        frames: indices.iter().map(|&i| traj.frames[i].clone()).collect(),
        metadata: traj.metadata.clone(),
    };

    // 目标格式装得下多帧就写一个文件，装不下就一帧一个 —— 往 POSCAR 写 20 帧
    // 本来就只能是 20 个文件，不必再要用户记一个开关
    if holds_multiple_frames(output) || picked.n_frames() == 1 {
        write_trajectory(&picked, output, units)?;
        println!(
            "Converted {} -> {}  ({} of {n_read} frame{} in 1 file)",
            input.display(),
            output.display(),
            picked.n_frames(),
            if n_read == 1 { "" } else { "s" },
        );
    } else {
        let width = index_width(*indices.last().unwrap_or(&0));
        for (&frame_idx, frame) in indices.iter().zip(&picked.frames) {
            let path = indexed_path(output, frame_idx, width);
            let single = Trajectory::from_frame(frame.clone());
            write_trajectory(&single, &path, units)
                .with_context(|| format!("writing frame {frame_idx} to {}", path.display()))?;
        }
        println!(
            "Converted {} -> {} files  ({} of {n_read} frames, {} holds one structure each)",
            input.display(),
            indices.len(),
            indices.len(),
            output
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("this format"),
        );
        println!(
            "  {}  ...  {}",
            indexed_path(output, indices[0], width).display(),
            indexed_path(output, *indices.last().unwrap(), width).display(),
        );
    }
    Ok(())
}

/// Zero-padded width for frame indices, at least 4 so `ls` sorts the products in
/// frame order.
fn index_width(max_index: usize) -> usize {
    max_index.to_string().len().max(4)
}

/// Inserts `_<frame index>` before the extension: `out.vasp` -> `out_0010.vasp`,
/// and `POSCAR` -> `POSCAR_0010`.
///
/// The number goes before the extension so the product keeps the extension its
/// format is detected by; extension-less VASP names keep their POSCAR prefix, which
/// is what detects them.
fn indexed_path(base: &Path, frame_idx: usize, width: usize) -> PathBuf {
    let stem = base.file_stem().and_then(|s| s.to_str()).unwrap_or("frame");
    let mut name = format!("{stem}_{frame_idx:0width$}");
    if let Some(ext) = base.extension().and_then(|e| e.to_str()) {
        name.push('.');
        name.push_str(ext);
    }
    base.with_file_name(name)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cmd(input: Option<&str>, output: Option<&str>) -> ConvertCmd {
        ConvertCmd {
            input: input.map(PathBuf::from),
            output: output.map(PathBuf::from),
            start: None,
            end: None,
            stride: None,
            number: None,
            metal_units: false,
        }
    }

    fn fixture() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../tests/43Z43P15A_NPT_5.lammpstrj")
    }

    #[test]
    fn test_indexed_path_keeps_the_name_its_format_is_detected_by() {
        // 扩展名格式：序号插在扩展名之前，扩展名必须留住
        assert_eq!(indexed_path(Path::new("out.lmp"), 10, 4), PathBuf::from("out_0010.lmp"));
        // POSCAR 没有扩展名，靠前缀识别 —— 序号追加在后面，前缀仍在
        assert_eq!(indexed_path(Path::new("POSCAR"), 2, 4), PathBuf::from("POSCAR_0002"));
        // 目录部分原样保留
        assert_eq!(
            indexed_path(Path::new("runs/700K/conf.lmp"), 7, 4),
            PathBuf::from("runs/700K/conf_0007.lmp")
        );
    }

    /// 序号补零到至少 4 位，否则 `ls` 会把 frame 10 排在 frame 2 前面
    #[test]
    fn test_index_width_pads_for_ls_ordering() {
        assert_eq!(index_width(9), 4);
        assert_eq!(index_width(9999), 4);
        assert_eq!(index_width(10000), 5);
        let w = index_width(12345);
        assert_eq!(indexed_path(Path::new("a.xyz"), 7, w), PathBuf::from("a_00007.xyz"));
    }

    #[test]
    fn test_frame_range_is_validated_before_the_file_is_read() {
        // 输入根本不存在：若先读文件，报的会是「打不开」而不是参数错误
        let mut c = cmd(Some("no_such_file.dump"), Some("out.xyz"));
        c.start = Some(9);
        c.end = Some(2);
        let err = run(&c).unwrap_err();
        assert!(err.to_string().contains("--start 9 is past --end 2"), "{err}");

        let mut c = cmd(Some("no_such_file.dump"), Some("out.xyz"));
        c.stride = Some(0);
        assert!(run(&c).unwrap_err().to_string().contains("--stride"));

        let mut c = cmd(Some("no_such_file.dump"), Some("out.xyz"));
        c.number = Some(0);
        assert!(run(&c).unwrap_err().to_string().contains("--number"));
    }

    #[test]
    fn test_selecting_nothing_is_an_error_that_names_the_available_range() {
        let out = std::env::temp_dir().join("ferro_convert_empty.xyz");
        let mut c = cmd(fixture().to_str(), out.to_str());
        c.start = Some(99);
        let err = run(&c).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("no frames selected"), "{msg}");
        assert!(msg.contains("0..=4"), "错误信息要说清可用范围: {msg}");
        assert!(!out.exists(), "选不到帧时不该留下文件");
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
