//! Shared plumbing for AIMD output readers (CP2K, VASP).
//!
//! The three formats print different things, but a caller collecting frames
//! into a dataset asks the same questions of all of them: how many steps did
//! the file claim, how many survived, and why did the rest go. [`AimdStats`]
//! answers those; [`sniff`] decides which reader to hand the file to.

use std::path::Path;

use anyhow::{bail, Context, Result};
use ferro_core::Trajectory;

/// Which AIMD program wrote a file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AimdFormat {
    Cp2kMd,
    Cp2kSp,
    VaspOutcar,
    VaspXml,
}

impl AimdFormat {
    /// Human-readable name, for messages.
    pub fn name(self) -> &'static str {
        match self {
            Self::Cp2kMd => "CP2K MD output",
            Self::Cp2kSp => "CP2K single-point output",
            Self::VaspOutcar => "VASP OUTCAR",
            Self::VaspXml => "VASP vasprun.xml",
        }
    }

    /// How this reader decides a frame's SCF converged.
    ///
    /// The two VASP paths do not use the same criterion, and a caller comparing
    /// drop counts between them deserves to be told which one applied.
    pub fn convergence_rule(self) -> &'static str {
        match self {
            Self::Cp2kMd => "SCF run converged",
            Self::Cp2kSp => "SCF run converged",
            Self::VaspOutcar => "EDIFF reached",
            Self::VaspXml => "SCF steps < NELM",
        }
    }
}

/// What a reader dropped, and what the file claimed to hold.
///
/// Frame dropping happens inside the readers because every criterion needs the
/// surrounding text (the SCF marker above the anchor, the block line count, the
/// first frame's composition); the counts travel out so the caller can report
/// them instead of a reader printing behind its back.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AimdStats {
    pub format: AimdFormat,
    /// Frame anchors seen, i.e. steps the file claims
    pub n_steps: usize,
    pub n_kept: usize,
    pub n_scf_failed: usize,
    pub n_incomplete: usize,
    pub n_bad_composition: usize,
    /// Frames whose block offsets differ from the first frame's.
    ///
    /// CP2K only — the VASP readers anchor every block by token and never
    /// count on an offset, so they have nothing to drift.
    pub n_layout_drift: usize,
    /// Restart markers in the file; > 1 means the run was restarted in place.
    ///
    /// CP2K only: a restarted VASP run writes a separate OUTCAR rather than
    /// appending, so this stays 0 there.
    pub n_restarts: usize,
    /// First and last step number the file claims.
    ///
    /// The span a file covers, not the frames that survived. Restarting from a
    /// checkpoint makes two files overlap here, and a caller concatenating them
    /// can only show that overlap if it knows the spans.
    pub steps: Option<(i64, i64)>,
    /// Program version the file reports, when it prints one.
    pub version: Option<String>,
    /// Why that version deserves a second look, or `None` when it does not.
    ///
    /// CP2K reprints its log layout between releases, and ferro has fixtures
    /// for only some of them. A reader that recognises the blocks anyway still
    /// says so here rather than printing behind the caller's back — same rule
    /// as the drop counts above.
    pub version_note: Option<String>,
}

impl AimdStats {
    pub fn new(format: AimdFormat) -> Self {
        Self {
            format,
            n_steps: 0,
            n_kept: 0,
            n_scf_failed: 0,
            n_incomplete: 0,
            n_bad_composition: 0,
            n_layout_drift: 0,
            n_restarts: 0,
            steps: None,
            version: None,
            version_note: None,
        }
    }

    pub fn n_dropped(&self) -> usize {
        self.n_scf_failed + self.n_incomplete + self.n_bad_composition
    }
}

/// The value of CP2K's `GLOBAL| Run type` line, when the head holds one.
fn run_type(head: &str) -> Option<String> {
    head.lines()
        .find(|l| {
            let mut it = l.split_whitespace();
            ["GLOBAL|", "Run", "type"].iter().all(|t| it.next() == Some(*t))
        })
        .and_then(|l| l.split_whitespace().last())
        .map(|v| v.to_string())
}

/// Identifies an AIMD output file by its content, not by its name.
///
/// Naming cannot carry this: VASP writes `OUTCAR` with no extension at all,
/// people rename it to `run.outcar`, and `.out` is too generic to belong to
/// any one program. The banners are unambiguous and cost one read of the head
/// of the file.
pub fn sniff(path: &Path) -> Result<AimdFormat> {
    use std::io::{BufRead, BufReader};

    let file = std::fs::File::open(path)
        .with_context(|| format!("cannot open {}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut head = String::new();
    // CP2K 的 `GLOBAL| Run type` 在第 40~50 行之间,比横幅还靠后一点
    for _ in 0..96 {
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        head.push_str(&line);
    }

    // vasprun 的 <?xml 恒在第一行;OUTCAR 的 " vasp.6.4.2 ..." 也在最前面几行。
    // CP2K 的横幅稍靠后,但仍在头 64 行内
    if head.contains("<?xml") || head.contains("<modeling>") {
        return Ok(AimdFormat::VaspXml);
    }
    if head.contains(" vasp.") || head.contains("vasp.5") || head.contains("vasp.6") {
        return Ok(AimdFormat::VaspOutcar);
    }
    if head.contains("CP2K|") || head.contains("**** **** ******  **  PROGRAM STARTED") {
        // 同一个 .out 扩展名下 CP2K 写两种完全不同的布局,由它自己的 run type
        // 区分。读不到那一行就说出来:此前默认按 MD 走,于是一份被裁过头的单点
        // 日志会静默变成「MD 找不到 Step number 锚点」,报的错与真正的原因无关
        return match run_type(&head).as_deref() {
            Some("ENERGY") | Some("ENERGY_FORCE") => Ok(AimdFormat::Cp2kSp),
            Some(_) => Ok(AimdFormat::Cp2kMd),
            None => bail!(
                "{}: a CP2K banner but no `GLOBAL| Run type` line, so the layout \
                 cannot be told apart. Pass --format cp2k/md or --format cp2k/sp",
                path.display()
            ),
        };
    }
    bail!(
        "{}: cannot tell which program wrote this. Recognised: VASP OUTCAR \
         (a `vasp.X.Y` banner), VASP vasprun.xml (an `<?xml` declaration) and \
         CP2K output (a `CP2K|` banner)",
        path.display()
    )
}

/// Reads an AIMD output as the given format, reporting what was dropped.
///
/// The format is a parameter rather than something decided here, so a caller
/// holding a better answer than the banner — the user, typically — can say so.
pub fn read_aimd_as(path: &Path, fmt: AimdFormat) -> Result<(Trajectory, AimdStats)> {
    match fmt {
        AimdFormat::Cp2kMd => super::cp2k_md::read_cp2k_md_with_stats(path),
        AimdFormat::Cp2kSp => super::cp2k_sp::read_cp2k_sp_with_stats(path),
        AimdFormat::VaspOutcar => super::vasp_outcar::read_vasp_outcar_with_stats(path),
        AimdFormat::VaspXml => super::vasprun::read_vasprun_with_stats(path),
    }
}

/// Reads any recognised AIMD output, identifying the format by its banner.
pub fn read_aimd_with_stats(path: &Path) -> Result<(Trajectory, AimdStats)> {
    read_aimd_as(path, sniff(path)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn p(name: &str) -> PathBuf {
        PathBuf::from("../tests").join(name)
    }

    #[test]
    fn sniffs_by_content_not_by_name() {
        // 这两份 fixture 的名字都不是 VASP 自己写出来的名字
        assert_eq!(sniff(&p("vasp_OUTCAR_2frames")).unwrap(), AimdFormat::VaspOutcar);
        assert_eq!(sniff(&p("vasp_vasprun_2frames.xml")).unwrap(), AimdFormat::VaspXml);
    }

    /// Two CP2K layouts share the `.out` extension; the run type tells them apart.
    #[test]
    fn a_cp2k_out_is_split_by_its_run_type_not_its_extension() {
        assert_eq!(sniff(&p("cp2k_md_3frames.out")).unwrap(), AimdFormat::Cp2kMd);
        assert_eq!(sniff(&p("5Al_0003_1500K_f394.out")).unwrap(), AimdFormat::Cp2kSp);
    }

    #[test]
    fn a_run_type_line_is_found_wherever_cp2k_puts_it() {
        let head = " GLOBAL| Run type                                          ENERGY_FORCE\n";
        assert_eq!(run_type(head).as_deref(), Some("ENERGY_FORCE"));
        // MD_PAR| 之类的别行不该被当成 run type
        assert_eq!(run_type(" MD_PAR| Ensemble type   NVT\n"), None);
    }

    /// A CP2K file with no run type is ambiguous, and says so.
    ///
    /// Defaulting to MD used to hide this: a truncated single-point log became
    /// "MD| Step number not found", an error about the wrong thing entirely.
    #[test]
    fn a_cp2k_file_without_a_run_type_asks_for_format() {
        let f = std::env::temp_dir().join("ferro_no_run_type.out");
        std::fs::write(&f, " CP2K| version string: CP2K version 2025.2\n hello\n").unwrap();
        let msg = format!("{:#}", sniff(&f).unwrap_err());
        assert!(msg.contains("--format cp2k/md"), "{msg}");
        assert!(msg.contains("Run type"), "{msg}");
        let _ = std::fs::remove_file(&f);
    }

    #[test]
    fn an_unrecognised_file_names_what_is_recognised() {
        let f = std::env::temp_dir().join("nothing.out");
        std::fs::write(&f, "hello\nworld\n").unwrap();
        let e = sniff(&f).unwrap_err();
        let msg = format!("{e:#}");
        for expect in ["OUTCAR", "vasprun.xml", "CP2K"] {
            assert!(msg.contains(expect), "{msg}");
        }
    }

    #[test]
    fn the_convergence_rule_is_named_per_format() {
        // 两条 VASP 路径判据不同,丢帧数对不上时得说得清为什么
        assert_ne!(
            AimdFormat::VaspOutcar.convergence_rule(),
            AimdFormat::VaspXml.convergence_rule()
        );
    }

    #[test]
    #[ignore = "reads the multi-hundred-MB files under examples/"]
    fn the_full_files_agree_with_dpdata() {
        let (o, os) = read_aimd_with_stats(
            std::path::Path::new("../examples/50Z50P_0.970_3000K.outcar")).unwrap();
        assert_eq!(o.n_frames(), 2000, "dpdata reads 2000 frames");
        assert_eq!(os.n_restarts, 1, "1..1575 then 1..425");
        assert_eq!(os.n_incomplete, 1, "step 1576 was cut off mid-run");

        let (x, xs) = read_aimd_with_stats(
            std::path::Path::new("../examples/vasprun.xml")).unwrap();
        assert_eq!(x.n_frames(), 425);
        assert_eq!(xs.n_dropped(), 0);

        // vasprun 恰好是那份 OUTCAR 的第二段:末帧必须是同一个构型
        let a = o.frames.last().unwrap();
        let b = x.frames.last().unwrap();
        assert!((a.energy.unwrap() - b.energy.unwrap()).abs() < 1e-8);
        for i in [0usize, 296] {
            let (pa, pb) = (a.atom(i).position, b.atom(i).position);
            // OUTCAR 打 5 位小数,vasprun 打分数坐标全精度
            assert!((pa - pb).norm() < 1e-4, "atom {i}: {pa:?} vs {pb:?}");
        }
    }
}
