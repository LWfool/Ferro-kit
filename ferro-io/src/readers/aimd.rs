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
    Cp2kOut,
    VaspOutcar,
    VaspXml,
}

impl AimdFormat {
    /// Human-readable name, for messages.
    pub fn name(self) -> &'static str {
        match self {
            Self::Cp2kOut => "CP2K output",
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
            Self::Cp2kOut => "SCF run converged",
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
        }
    }

    pub fn n_dropped(&self) -> usize {
        self.n_scf_failed + self.n_incomplete + self.n_bad_composition
    }
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
    for _ in 0..64 {
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
        return Ok(AimdFormat::Cp2kOut);
    }
    bail!(
        "{}: cannot tell which program wrote this. Recognised: VASP OUTCAR \
         (a `vasp.X.Y` banner), VASP vasprun.xml (an `<?xml` declaration) and \
         CP2K output (a `CP2K|` banner)",
        path.display()
    )
}

/// Reads any recognised AIMD output, reporting what was dropped.
pub fn read_aimd_with_stats(path: &Path) -> Result<(Trajectory, AimdStats)> {
    let fmt = sniff(path)?;
    let name = path.to_string_lossy();
    match fmt {
        AimdFormat::Cp2kOut => super::cp2k_out::read_cp2k_out_with_stats(&name),
        AimdFormat::VaspOutcar => super::vasp_outcar::read_vasp_outcar_with_stats(&name),
        AimdFormat::VaspXml => super::vasprun::read_vasprun_with_stats(&name),
    }
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
