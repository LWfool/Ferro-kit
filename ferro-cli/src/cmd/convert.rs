use anyhow::{bail, Context, Result};
use clap::Args;
use crate::args::common::ReadArgs;
use crate::io_dispatch::{holds_multiple_frames, read_trajectory, write_trajectory};
use ferro_analysis::ml::merge::sort_atoms;
use ferro_core::Trajectory;
use std::path::{Path, PathBuf};

#[derive(Args, Debug)]
pub struct ConvertCmd {
    /// Input file (format auto-detected; omit to list the supported formats)
    #[arg(short, long)]
    pub input: Option<PathBuf>,

    /// Output file (format from its name), or a directory when --format is given;
    /// missing directories are created
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Treat -o as a directory and write this format into it (an extension: cif, vasp, extxyz…)
    #[arg(long, value_name = "EXT")]
    pub format: Option<String>,

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

    #[command(flatten)]
    pub read: ReadArgs,
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

    // 给了 --format 则 -o 是目录，产物名为 frame.<EXT>（逐帧时 frame_0000.<EXT>）；
    // 否则 -o 是文件，目标格式从文件名推断。目录形状的 -o 不带 --format 时无格式
    // 可判（Rust 与 pathlib 一样丢掉末尾分隔符，只有原文能看出意图），报错并给出
    // 可直接抄的写法
    let output = &match &args.format {
        Some(ext) => output.join(format!("frame.{}", ext.trim_start_matches('.'))),
        None => {
            if ends_with_separator(output) || output.is_dir() {
                let shown = output.display();
                bail!(
                    "-o '{shown}' names a directory, but no format was given. \
                     Add --format (e.g. --format cif), or name a file: -o '{}'",
                    output.join("conf.extxyz").display()
                );
            }
            output.clone()
        }
    };
    // 参数级错误在读第一个文件之前失败：轨迹可能是几个 GB，读完再报「--start 比
    // --end 大」是白等。跨参数的一致性能在这里查的都查掉。
    // ferro-python 的 `Trajectory.frame_indices` 有同义的一份（报错用 Python 参数名，按 R6 不合）
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
    // 目标格式也是参数：认不出就别建目录、别读输入（审查 D-M15）
    crate::io_dispatch::out_format(output).with_context(|| match &args.format {
        Some(ext) => format!("--format {ext}"),
        None => "-o names the output FILE; to write into a directory, add --format EXT".into(),
    })?;
    // convert 只写这一个目标，路径是用户逐字写出的，不必再问是否建目录
    crate::outpath::ensure_parent(output, true)?;

    let traj = read_trajectory(input, &args.read)?;
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

    // 所有产物一律按 (Z, 符号) 重排原子，像 POSCAR 那样同元素成块、块按周期表序，
    // 同一体系无论来自哪个程序，写出的原子顺序都一致。sort_atoms 逐帧各取置换 ——
    // xyz 等格式不保证帧间原子顺序不变
    let picked = sort_atoms(&Trajectory {
        frames: indices.iter().map(|&i| traj.frames[i].clone()).collect(),
        metadata: traj.metadata.clone(),
    });

    // 目标格式装得下多帧就写一个文件，装不下就一帧一个 —— 往 POSCAR 写 20 帧
    // 本来就只能是 20 个文件，不必再要用户记一个开关
    if holds_multiple_frames(output) || picked.n_frames() == 1 {
        write_trajectory(&picked, output, args.read.units)?;
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
            write_trajectory(&single, &path, args.read.units)
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

/// True when the user wrote a trailing separator (`-o out/`).
///
/// `Path` drops it the way `pathlib` does — `Path::new("out/").file_name()` is `Some("out")`
/// — so the raw text is the only place that intent survives. Without `--format`,
/// `convert` has no format to write into such a directory.
fn ends_with_separator(p: &Path) -> bool {
    let s = p.to_string_lossy();
    s.ends_with(std::path::MAIN_SEPARATOR) || s.ends_with('/')
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
            read: ReadArgs::default(),
            format: None,
        }
    }

    /// `-o out/` 没有文件名可判格式：不给 --format 就报错，并给出两种可抄的写法
    #[test]
    fn a_directory_without_format_is_refused_with_a_usable_hint() {
        let existing = std::env::temp_dir();
        for out in ["out/", existing.to_str().unwrap()] {
            let err = run(&cmd(Some("in.xyz"), Some(out))).unwrap_err().to_string();
            assert!(err.contains("names a directory"), "{out}: {err}");
            assert!(err.contains("--format"), "{out}: 要提示 --format: {err}");
            assert!(err.contains("conf.extxyz"), "{out}: 要给出可直接抄的文件写法: {err}");
        }
    }

    /// --format 让 -o 成为目录：缺失的目录直接建（不问、无需 --mkdir），单结构格式
    /// 逐帧 frame_NNNN.EXT，多帧格式一个 frame.EXT
    #[test]
    fn format_turns_output_into_a_directory_created_on_demand() {
        let root = std::env::temp_dir().join(format!("ferro_convert_dir_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        for (ext, want) in [
            ("cif", vec!["frame_0000.cif", "frame_0004.cif"]),
            (".vasp", vec!["frame_0000.vasp", "frame_0004.vasp"]),
            ("extxyz", vec!["frame.extxyz"]),
        ] {
            let dir = root.join(ext.trim_start_matches('.')).join("nested");
            let mut c = cmd(fixture().to_str(), dir.to_str());
            c.read.units = Some(ferro_io::LammpsUnits::Metal);
            c.format = Some(ext.into());
            run(&c).unwrap_or_else(|e| panic!("--format {ext}: {e}"));
            for name in want {
                assert!(dir.join(name).is_file(), "--format {ext}: 缺 {name}");
            }
        }
        let _ = std::fs::remove_dir_all(&root);
    }

    /// 审查 D-M15：参数与目标格式错了，应在建 -o 的目录、读输入之前报错，不留空目录
    #[test]
    fn bad_values_fail_before_the_directory_is_made() {
        let dir = std::env::temp_dir().join("ferro_d_m15_convert_must_not_exist");
        let _ = std::fs::remove_dir_all(&dir);
        let good = dir.join("out.extxyz");
        let bad_fmt = dir.join("out.nosuchfmt");
        type Case<'a> = (&'a Path, fn(&mut ConvertCmd), &'a str);
        let cases: [Case; 5] = [
            (&good, |c| c.stride = Some(0), "--stride"),
            (&good, |c| c.number = Some(0), "--number"),
            (&good, |c| { c.start = Some(5); c.end = Some(2); }, "--start 5 is past --end 2"),
            (&bad_fmt, |_| {}, "Unsupported output format"),
            (&dir, |c| c.format = Some("nosuchfmt".into()), "Unsupported output format"),
        ];
        for (out, tweak, want) in cases {
            let mut c = cmd(Some("no_such_input.xyz"), Some(out.to_str().unwrap()));
            tweak(&mut c);
            let err = format!("{:#}", run(&c).unwrap_err());
            assert!(err.contains(want), "应报「{want}」，实际 {err}");
            assert!(!dir.exists(), "「{want}」：参数错误时不该建 -o 的目录");
        }
    }

    #[test]
    fn a_path_with_directories_is_taken_as_written() {
        assert!(!ends_with_separator(Path::new("out/conf.extxyz")));
        assert!(ends_with_separator(Path::new("out/")));
        assert!(!ends_with_separator(Path::new("out")));
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

    /// 多帧写 CIF：每帧一个文件、每个文件一个 data_ 块。VESTA 等只认第一个块，
    /// 写成一个多块文件时用户看到的只有第 0 帧
    #[test]
    fn test_cif_is_written_one_frame_per_file() {
        let dir = std::env::temp_dir().join(format!("ferro_convert_cif_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let mut c = cmd(fixture().to_str(), dir.join("conf.cif").to_str());
        c.read.units = Some(ferro_io::LammpsUnits::Metal);
        run(&c).unwrap();
        for i in 0..5 {
            let p = dir.join(format!("conf_{i:04}.cif"));
            let text = std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
            assert_eq!(text.matches("data_").count(), 1, "{} 应只有一个 data_ 块", p.display());
        }
        assert!(!dir.join("conf.cif").exists(), "不该再写出多块的 conf.cif");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// 输出的原子按 (Z, 符号) 排：Zn Al P O → O Al P Zn。每个原子的坐标与力
    /// 跟着自己走；逐帧取置换，第 1 帧的原子顺序与第 0 帧不同也照样排对
    #[test]
    fn atoms_are_written_in_atomic_number_order() {
        use ferro_core::{Atom, Cell, Frame};
        use nalgebra::{Matrix3, Vector3};
        let frame_of = |elements: &[&str]| {
            let mut f = Frame::with_cell(Cell::from_matrix(Matrix3::identity() * 10.0), [true; 3]);
            // 坐标与力都由元素决定，排序后能逐原子核对是否跟对了
            let x = |e: &str| match e { "O" => 1.0, "Al" => 2.0, "P" => 3.0, _ => 4.0 };
            f.atoms = elements.iter().map(|e| Atom::new(*e, Vector3::new(x(e), 0.5, 0.5))).collect();
            f.forces = Some(elements.iter().map(|e| Vector3::new(-x(e), 0.0, 0.0)).collect());
            f
        };
        let dir = std::env::temp_dir().join(format!("ferro_convert_sort_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let input = dir.join("in.extxyz");
        let t = Trajectory {
            frames: vec![frame_of(&["Zn", "Al", "P", "O"]), frame_of(&["O", "P", "Zn", "Al"])],
            metadata: Default::default(),
        };
        write_trajectory(&t, &input, None).unwrap();

        let output = dir.join("out.extxyz");
        run(&cmd(input.to_str(), output.to_str())).unwrap();
        let back = read_trajectory(&output, &ReadArgs::default()).unwrap();
        for (k, f) in back.frames.iter().enumerate() {
            let order: Vec<&str> = f.atoms.iter().map(|a| a.element.as_str()).collect();
            assert_eq!(order, ["O", "Al", "P", "Zn"], "第 {k} 帧应按原子序数排");
            let forces = f.forces.as_ref().expect("力应写出");
            for (a, fv) in f.atoms.iter().zip(forces) {
                assert_eq!(fv.x, -a.position.x, "第 {k} 帧 {} 的力没跟着原子走", a.element);
            }
        }
        let _ = std::fs::remove_dir_all(&dir);
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
