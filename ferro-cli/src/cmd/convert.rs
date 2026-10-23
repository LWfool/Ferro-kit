use anyhow::{bail, Context, Result};
use clap::Args;
use crate::args::common::ReadArgs;
use crate::io_dispatch::read_trajectory;
use ferro_core::Trajectory;
use std::path::{Path, PathBuf};

#[derive(Args, Debug)]
pub struct ConvertCmd {
    /// Input file (format auto-detected; omit to list the supported formats)
    #[arg(short, long)]
    pub input: Option<PathBuf>,

    /// Output directory, created if missing; products are named after the input
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Format to write (an extension: cif, vasp, extxyz…)
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
    #[arg(short = 'N', long, value_name = "N", conflicts_with = "stride")]
    pub number: Option<usize>,

    /// Leave out a charge / magnetic-moment column only some atoms carry (refused otherwise)
    #[arg(long)]
    pub drop_partial: bool,

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
    let Some(dir) = &args.output else {
        bail!("convert needs an output directory: -o <DIR>  (run `ferro convert` for the format list)");
    };
    // Not clap's `required`: bare `ferro convert` must still reach the format list
    let Some(ext) = &args.format else {
        bail!(
            "convert needs a target format: --format <EXT>, e.g. `-o {} --format extxyz`",
            dir.display()
        );
    };

    // The product is named after the input: `traj.dump` -> DIR/traj.<EXT>, and a
    // format holding one structure per file numbers it, DIR/traj_0000.<EXT>.
    // A name without an extension (POSCAR) or a DeePMD directory is the stem as a whole
    let stem = input.file_stem().unwrap_or(input.as_os_str()).to_string_lossy();
    let output = &dir.join(format!("{stem}.{}", ext.trim_start_matches('.')));

    // Compared as resolved paths: `-o .` for `-i x.extxyz --format extxyz` is the
    // same file, and writing it truncates the input. The output need not exist yet,
    // so its parent is resolved.
    let resolved = |p: &Path| -> PathBuf {
        p.canonicalize().unwrap_or_else(|_| {
            let parent = p.parent().filter(|d| !d.as_os_str().is_empty()).unwrap_or(Path::new("."));
            let name = p.file_name().map(PathBuf::from).unwrap_or_default();
            parent.canonicalize().map(|d| d.join(&name)).unwrap_or_else(|_| p.to_path_buf())
        })
    };
    if resolved(input) == resolved(output) {
        bail!("Input and output paths are identical");
    }

    // 参数级错误在读第一个文件之前失败：轨迹可能是几个 GB，读完再报「--start 比
    // --end 大」是白等。跨参数的一致性能在这里查的都查掉。
    // ferro-python 的 `Trajectory.frame_indices` 有同义的一份（报错用 Python 参数名，按 R6 不合）
    // `dataset filter` checks its frame selection the same way (`FilterParams::validate`)
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
    ferro_io::Format::for_output(output).with_context(|| format!("--format {ext}"))?;
    // convert writes into this one directory, typed out by the user: no need to ask
    crate::outpath::ensure_dir(dir, true)?;

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

    // Every product has its atoms in (Z, symbol) order, one block per element in
    // periodic-table order as in a POSCAR, so one system is written in the same atom
    // order whichever program it came from. The sort takes a permutation per frame:
    // xyz and the like do not promise a constant atom order
    let picked = Trajectory {
        frames: indices.iter().map(|&i| traj.frames[i].clone()).collect(),
        metadata: traj.metadata.clone(),
    }
    .canonically_sorted();

    // One file when the format holds several frames, one per frame when it does not:
    // 20 frames to POSCAR can only be 20 files, so there is no switch to remember.
    // The per-frame files are numbered by the frame they came from in the input
    let opts = ferro_io::WriteOptions { units: args.read.units, drop_partial: args.drop_partial };
    let written = ferro_io::write_numbered(&picked, output, &opts, &indices)?;
    if written.len() == 1 {
        println!(
            "Converted {} -> {}  ({} of {n_read} frame{} in 1 file)",
            input.display(),
            output.display(),
            picked.n_frames(),
            if n_read == 1 { "" } else { "s" },
        );
    } else {
        println!(
            "Converted {} -> {} files  ({} of {n_read} frames, {} holds one structure each)",
            input.display(),
            indices.len(),
            indices.len(),
            ext.trim_start_matches('.'),
        );
        println!(
            "  {}  ...  {}",
            written[0].display(),
            written[written.len() - 1].display(),
        );
    }
    Ok(())
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
            drop_partial: false,
            read: ReadArgs::default(),
            format: Some("extxyz".into()),
        }
    }

    /// `-o .` / `-o sub/..` for `-i same.extxyz --format extxyz` lands on the input itself:
    /// comparing the text would let it through and the write truncate the input.
    #[test]
    fn the_same_file_spelt_differently_is_refused() {
        let dir = std::env::temp_dir().join(format!("ferro_convert_same_{}", std::process::id()));
        std::fs::create_dir_all(dir.join("sub")).unwrap();
        let input = dir.join("same.extxyz");
        std::fs::write(&input, "1\nProperties=species:S:1:pos:R:3\nO 0 0 0\n").unwrap();
        for out in [dir.join("."), dir.join("sub/..")] {
            let err = run(&cmd(Some(input.to_str().unwrap()), Some(out.to_str().unwrap())))
                .unwrap_err().to_string();
            assert!(err.contains("identical"), "{}: {err}", out.display());
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Without --format there is nothing to name the target format: refused before
    /// anything is read, with a line to copy
    #[test]
    fn a_missing_format_is_refused_with_a_usable_hint() {
        let mut c = cmd(Some("no_such_input.xyz"), Some("out"));
        c.format = None;
        let err = run(&c).unwrap_err().to_string();
        assert!(err.contains("-o out --format extxyz"), "should give a line to copy: {err}");
    }

    /// -o is a directory, created when missing (no question, no --mkdir). Products are
    /// named after the input: one STEM_NNNN.EXT per frame for a single-structure format,
    /// one STEM.EXT for a trajectory format
    #[test]
    fn output_is_a_directory_created_on_demand_and_named_after_the_input() {
        let root = std::env::temp_dir().join(format!("ferro_convert_dir_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        for (ext, want) in [
            ("cif", vec!["43Z43P15A_NPT_5_0000.cif", "43Z43P15A_NPT_5_0004.cif"]),
            (".vasp", vec!["43Z43P15A_NPT_5_0000.vasp", "43Z43P15A_NPT_5_0004.vasp"]),
            ("extxyz", vec!["43Z43P15A_NPT_5.extxyz"]),
        ] {
            let dir = root.join(ext.trim_start_matches('.')).join("nested");
            let mut c = cmd(fixture().to_str(), dir.to_str());
            c.read.units = Some(ferro_io::LammpsUnits::Metal);
            c.format = Some(ext.into());
            run(&c).unwrap_or_else(|e| panic!("--format {ext}: {e}"));
            for name in want {
                assert!(dir.join(name).is_file(), "--format {ext}: missing {name}");
            }
        }
        let _ = std::fs::remove_dir_all(&root);
    }

    /// `-N` is `--number`, as in `dataset filter` (review D-15)
    #[test]
    fn short_n_is_number() {
        #[derive(clap::Parser, Debug)]
        struct Wrap {
            #[command(flatten)]
            cmd: ConvertCmd,
        }
        use clap::Parser;
        let w = Wrap::try_parse_from(["convert", "-i", "a.xyz", "-o", "b.xyz", "-N", "3"]).unwrap();
        assert_eq!(w.cmd.number, Some(3));
        assert!(Wrap::try_parse_from(["convert", "-i", "a.xyz", "-o", "b.xyz", "-N", "3", "--stride", "2"]).is_err());
    }

    /// Review D-M15: bad values and an unknown target format fail before -o is created
    /// or the input read, leaving no empty directory behind
    #[test]
    fn bad_values_fail_before_the_directory_is_made() {
        let dir = std::env::temp_dir().join("ferro_d_m15_convert_must_not_exist");
        let _ = std::fs::remove_dir_all(&dir);
        type Case<'a> = (fn(&mut ConvertCmd), &'a str);
        let cases: [Case; 5] = [
            (|c| c.stride = Some(0), "--stride"),
            (|c| c.number = Some(0), "--number"),
            (|c| { c.start = Some(5); c.end = Some(2); }, "--start 5 is past --end 2"),
            (|c| c.format = Some("nosuchfmt".into()), "Unsupported output format"),
            (|c| c.format = None, "--format"),
        ];
        for (tweak, want) in cases {
            let mut c = cmd(Some("no_such_input.xyz"), dir.to_str());
            tweak(&mut c);
            let err = format!("{:#}", run(&c).unwrap_err());
            assert!(err.contains(want), "should report '{want}', got {err}");
            assert!(!dir.exists(), "'{want}': a bad value must not create the -o directory");
        }
    }

    /// A name without an extension (POSCAR) is the stem as a whole
    #[test]
    fn an_input_without_extension_names_the_product_whole() {
        let dir = std::env::temp_dir().join(format!("ferro_convert_noext_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let read = ReadArgs { units: Some(ferro_io::LammpsUnits::Metal), ..ReadArgs::default() };
        let mut t = read_trajectory(&fixture(), &read).unwrap();
        t.frames.truncate(1);
        let poscar = dir.join("POSCAR");
        ferro_io::write(&t, &poscar, &ferro_io::WriteOptions::default()).unwrap();
        run(&cmd(poscar.to_str(), dir.join("out").to_str())).unwrap();
        assert!(dir.join("out/POSCAR.extxyz").is_file());
        let _ = std::fs::remove_dir_all(&dir);
    }

    fn fixture() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../tests/43Z43P15A_NPT_5.lammpstrj")
    }

    #[test]
    fn test_frame_range_is_validated_before_the_file_is_read() {
        // 输入根本不存在：若先读文件，报的会是「打不开」而不是参数错误
        let mut c = cmd(Some("no_such_file.dump"), Some("out"));
        c.start = Some(9);
        c.end = Some(2);
        let err = run(&c).unwrap_err();
        assert!(err.to_string().contains("--start 9 is past --end 2"), "{err}");

        let mut c = cmd(Some("no_such_file.dump"), Some("out"));
        c.stride = Some(0);
        assert!(run(&c).unwrap_err().to_string().contains("--stride"));

        let mut c = cmd(Some("no_such_file.dump"), Some("out"));
        c.number = Some(0);
        assert!(run(&c).unwrap_err().to_string().contains("--number"));
    }

    #[test]
    fn test_selecting_nothing_is_an_error_that_names_the_available_range() {
        let out = std::env::temp_dir().join(format!("ferro_convert_empty_{}", std::process::id()));
        let mut c = cmd(fixture().to_str(), out.to_str());
        c.start = Some(99);
        let err = run(&c).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("no frames selected"), "{msg}");
        assert!(msg.contains("0..=4"), "错误信息要说清可用范围: {msg}");
        assert!(!out.join("43Z43P15A_NPT_5.extxyz").exists(), "selecting nothing must leave no file");
        let _ = std::fs::remove_dir_all(&out);
    }

    /// 多帧写 CIF：每帧一个文件、每个文件一个 data_ 块。VESTA 等只认第一个块，
    /// 写成一个多块文件时用户看到的只有第 0 帧
    #[test]
    fn test_cif_is_written_one_frame_per_file() {
        let dir = std::env::temp_dir().join(format!("ferro_convert_cif_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let mut c = cmd(fixture().to_str(), dir.to_str());
        c.format = Some("cif".into());
        c.read.units = Some(ferro_io::LammpsUnits::Metal);
        run(&c).unwrap();
        for i in 0..5 {
            let p = dir.join(format!("43Z43P15A_NPT_5_{i:04}.cif"));
            let text = std::fs::read_to_string(&p).unwrap_or_else(|e| panic!("{}: {e}", p.display()));
            assert_eq!(text.matches("data_").count(), 1, "{} should hold one data_ block", p.display());
        }
        assert!(!dir.join("43Z43P15A_NPT_5.cif").exists(), "no multi-block CIF should be written");
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
            f.set_atoms(elements.iter().map(|e| Atom::new(*e, Vector3::new(x(e), 0.5, 0.5))));
            f.set_forces(Some(elements.iter().map(|e| Vector3::new(-x(e), 0.0, 0.0)).collect())).unwrap();
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
        ferro_io::write(&t, &input, &ferro_io::WriteOptions::default()).unwrap();

        let output = dir.join("out");
        run(&cmd(input.to_str(), output.to_str())).unwrap();
        let output = output.join("in.extxyz");
        let back = read_trajectory(&output, &ReadArgs::default()).unwrap();
        for (k, f) in back.frames.iter().enumerate() {
            let order: Vec<&str> = f.atoms().map(|a| a.element).collect();
            assert_eq!(order, ["O", "Al", "P", "Zn"], "第 {k} 帧应按原子序数排");
            let forces = f.forces().expect("力应写出");
            for (a, fv) in f.atoms().zip(forces) {
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
        let mut c = cmd(Some("same.xyz"), Some("."));
        c.format = Some("xyz".into());
        let err = run(&c).unwrap_err();
        assert!(err.to_string().contains("identical"), "{err}");
    }
}
