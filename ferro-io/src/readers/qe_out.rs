//! Quantum ESPRESSO `pw.x` output reader: `scf`, `md` and `vc-md`.
//!
//! Every ionic step runs one SCF, and pw.x opens each with the same header:
//!
//! ```text
//!      Self-consistent Calculation                       <- the step anchor
//!  ...
//! !    total energy              =    -262.21150703 Ry   <- `!` only when converged
//!      Forces acting on atoms (cartesian axes, Ry/au):   <- total forces: the FIRST group
//!      atom    1 type  1   force =    -0.00284143   -0.00306774    0.00000000
//!  ...
//!      The non-local contrib.  to forces                 <- verbosity 'high': six more
//!      atom    1 type  1   force =     0.30871553    0.33329197    0.00000000
//!  ...
//!           total   stress  (Ry/bohr**3)                   (kbar)     P=      -15.51
//!   -0.00012136   0.00004600  -0.00000000          -17.85        6.77       -0.00
//!  ...
//!      Entering Dynamics:    iteration =     1           <- md / vc-md
//! CELL_PARAMETERS (angstrom)                             <- vc-md: the NEXT step's cell
//! ATOMIC_POSITIONS (angstrom)                            <- the NEXT step's positions
//! ```
//!
//! # Frame layout
//!
//! Step k runs from its anchor to the next one. Its energy, forces and stress
//! are inside that span; its structure is NOT: step 0 uses the cell and
//! positions of the output header, and every later step the cards printed at
//! the end of the step before. The cards printed after the last step belong to
//! a step that never ran and are ignored, so an N-step run gives N frames.
//!
//! The force components printed under `verbosity = 'high'` (non-local, ionic,
//! local, core correction, Hubbard, SCF correction, DFT-D3) have exactly the
//! layout of the total forces; only the first group after the header is read.
//!
//! # What is dropped
//!
//! Steps whose SCF did not converge (`convergence NOT achieved`), and steps
//! that stop short: no `!` energy, no forces, or no stress although other steps
//! of the file print it (a run cut off mid-step). Forces are required, as in
//! the CP2K readers: a frame without them cannot enter a DeePMD system, so an
//! scf run with `tprnfor` off gives no frame and the error says so. Stress is
//! optional: a file that never prints it (`tstress` off) gives frames without.
//!
//! Optimisations (`relax`, `vc-relax`) are refused: they print a final SCF at
//! the relaxed structure whose step layout ferro has no fixture for.
//!
//! # Units
//!
//! Ry, Ry/bohr and Ry/bohr³ are converted to eV, eV/Å and eV/Å³. The stress
//! keeps pw.x's sign (positive = compression), matching
//! [`ferro_core::Frame::stress`]; ASE flips it.
//!
//! The header cell takes its scale from `celldm(1)` (six decimals) rather than
//! `lattice parameter (alat)` (four), and the header positions are printed with
//! seven decimals in alat units: frame 0 is good to about 1e-6 Å, later frames
//! to the nine or ten decimals of the cards.

use std::path::Path;

use anyhow::{bail, Context, Result};
use ferro_core::units::{BOHR_TO_ANG, HARTREE_TO_EV};
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};

use super::aimd::{AimdFormat, AimdStats};
use super::util::element_from_label;

const RY_TO_EV: f64 = HARTREE_TO_EV / 2.0;

/// (species label, Cartesian position in Å), in the file's atom order.
type Sites = Vec<(String, Vector3<f64>)>;

/// Major releases whose output layout ferro has a fixture for (7.5, plus a
/// 7.3.1 scf without forces under `examples/`).
const VERIFIED_MAJORS: &[&str] = &["7"];

/// True when the line's leading tokens match `tokens`.
///
/// Also in `cp2k_sp.rs`, `cp2k_md.rs` and `aimd.rs`, unshared on purpose; the
/// reason is on `cp2k_sp.rs::starts_with_tokens` (review D-2).
fn starts_with_tokens(line: &str, tokens: &[&str]) -> bool {
    let mut it = line.split_whitespace();
    tokens.iter().all(|t| it.next() == Some(*t))
}

/// The three numbers inside `= ( x y z )`, as in `a(1) = ( ... )` and `tau( 1) = ( ... )`.
fn vector_in_parens(line: &str) -> Option<Vector3<f64>> {
    let open = line.find("= (")? + 3;
    let close = line[open..].find(')')? + open;
    let v: Vec<f64> = line[open..close].split_whitespace().map(|t| t.parse().ok()).collect::<Option<_>>()?;
    (v.len() == 3).then(|| Vector3::new(v[0], v[1], v[2]))
}

/// The first three numbers of a line, other columns ignored.
fn three_numbers<'a>(mut it: impl Iterator<Item = &'a str>) -> Option<Vector3<f64>> {
    let x = it.next()?.parse().ok()?;
    let y = it.next()?.parse().ok()?;
    let z = it.next()?.parse().ok()?;
    Some(Vector3::new(x, y, z))
}

/// The number after `key` on its line, e.g. `celldm(1)=   8.104972`.
fn number_after(line: &str, key: &str) -> Option<f64> {
    let at = line.find(key)? + key.len();
    line[at..].split_whitespace().next()?.parse().ok()
}

/// Card option `(angstrom)` / `{bohr}` / `(alat= 8.10497200)` → (unit, alat if given).
fn card_unit(line: &str) -> (String, Option<f64>) {
    let opt = line.split_whitespace().skip(1).collect::<Vec<_>>().join(" ");
    let opt = opt.trim_matches(|c| matches!(c, '(' | ')' | '{' | '}')).to_lowercase();
    match opt.split_once('=') {
        Some((unit, value)) => (unit.trim().to_string(), value.trim().parse().ok()),
        None => (opt.trim().to_string(), None),
    }
}

/// A `CELL_PARAMETERS` card at line `i`: rows are lattice vectors, in Å. `base` is
/// the file line number of `lines[0]`, for messages.
fn read_cell_card(lines: &[&str], i: usize, base: usize, alat_ang: f64) -> Result<Matrix3<f64>> {
    let (unit, alat) = card_unit(lines[i]);
    let scale = match unit.as_str() {
        "angstrom" => 1.0,
        "bohr" => BOHR_TO_ANG,
        "alat" => alat.map_or(alat_ang, |a| a * BOHR_TO_ANG),
        other => bail!("line {}: CELL_PARAMETERS in unit `{other}` is not supported", base + i + 1),
    };
    let mut m = Matrix3::zeros();
    for row in 0..3 {
        let v = lines.get(i + 1 + row).and_then(|l| three_numbers(l.split_whitespace()))
            .with_context(|| format!("line {}: CELL_PARAMETERS needs three rows of three numbers", base + i + 1))?;
        m.set_row(row, &(v * scale).transpose());
    }
    Ok(m)
}

/// An `ATOMIC_POSITIONS` card at line `i`: (label, Cartesian Å) for `nat` atoms.
/// `base` is the file line number of `lines[0]`, for messages.
fn read_positions_card(
    lines: &[&str], i: usize, base: usize, nat: usize, cell: &Matrix3<f64>, alat_ang: f64,
) -> Result<Sites> {
    let (unit, _) = card_unit(lines[i]);
    let mut out = Vec::with_capacity(nat);
    for k in 0..nat {
        let line = lines.get(i + 1 + k).copied().unwrap_or("");
        let mut it = line.split_whitespace();
        let (Some(label), Some(v)) = (it.next(), three_numbers(&mut it)) else {
            bail!("line {}: ATOMIC_POSITIONS has fewer than {nat} atoms", base + i + 2 + k);
        };
        let pos = match unit.as_str() {
            "angstrom" => v,
            "bohr" => v * BOHR_TO_ANG,
            "alat" => v * alat_ang,
            // 行向量是晶格矢量：r = s · M
            "crystal" => (v.transpose() * cell).transpose(),
            other => bail!("line {}: ATOMIC_POSITIONS in unit `{other}` is not supported", base + i + 1),
        };
        out.push((label.to_string(), pos));
    }
    Ok(out)
}

/// Header cell and positions: `crystal axes` and `positions (alat units)`, scaled by alat.
fn read_header(lines: &[&str]) -> Result<(usize, f64, Matrix3<f64>, Sites)> {
    let nat = lines.iter()
        .find(|l| starts_with_tokens(l, &["number", "of", "atoms/cell"]))
        .and_then(|l| l.split('=').nth(1)?.trim().parse::<usize>().ok())
        .context("no `number of atoms/cell` line in the header")?;
    // celldm(1) 有六位小数，`lattice parameter (alat)` 只有四位
    let alat_bohr = lines.iter().find_map(|l| number_after(l, "celldm(1)="))
        .or_else(|| lines.iter()
            .find(|l| starts_with_tokens(l, &["lattice", "parameter", "(alat)"]))
            .and_then(|l| l.split('=').nth(1)?.split_whitespace().next()?.parse().ok()))
        .context("no `celldm(1)` or `lattice parameter (alat)` line in the header")?;
    let alat = alat_bohr * BOHR_TO_ANG;

    let axes = lines.iter().position(|l| l.trim_start().starts_with("crystal axes:"))
        .context("no `crystal axes` block in the header")?;
    let mut cell = Matrix3::zeros();
    for row in 0..3 {
        let v = lines.get(axes + 1 + row).and_then(|l| vector_in_parens(l))
            .with_context(|| format!("line {}: crystal axes need three `a(i) = ( x y z )` rows", axes + 1))?;
        cell.set_row(row, &(v * alat).transpose());
    }

    let site = lines.iter()
        .position(|l| starts_with_tokens(l, &["site", "n."]) && l.contains("positions (alat units)"))
        .context("no `positions (alat units)` block in the header")?;
    let mut atoms = Vec::with_capacity(nat);
    for k in 0..nat {
        let line = lines.get(site + 1 + k).copied().unwrap_or("");
        let (Some(label), Some(v)) = (line.split_whitespace().nth(1), vector_in_parens(line)) else {
            bail!("line {}: the header lists fewer than {nat} atom positions", site + 2 + k);
        };
        atoms.push((label.to_string(), v * alat));
    }
    Ok((nat, alat, cell, atoms))
}

/// The first force group after `Forces acting on atoms`: `nat` rows, Ry/bohr → eV/Å.
fn read_forces(lines: &[&str], head: usize, end: usize, nat: usize) -> Option<Vec<Vector3<f64>>> {
    let mut forces = Vec::with_capacity(nat);
    for line in lines[head + 1..end].iter().skip_while(|l| l.trim().is_empty()) {
        if forces.len() == nat || !starts_with_tokens(line, &["atom"]) {
            break;
        }
        let v = three_numbers(line.split('=').nth(1)?.split_whitespace())?;
        forces.push(v * (RY_TO_EV / BOHR_TO_ANG));
    }
    (forces.len() == nat).then_some(forces)
}

/// The `total stress` block: first three columns (Ry/bohr³; the last three are kbar) → eV/Å³.
fn read_stress(lines: &[&str], head: usize) -> Option<Matrix3<f64>> {
    let mut s = Matrix3::zeros();
    for row in 0..3 {
        let v = three_numbers(lines.get(head + 1 + row)?.split_whitespace())?;
        s.set_row(row, &(v * (RY_TO_EV / BOHR_TO_ANG.powi(3))).transpose());
    }
    Some(s)
}

// Judges QE releases; the two CP2K readers' `version_note`s judge CP2K ones, so
// the three change for different reasons and stay apart (review D-4).
fn version_note(version: Option<&str>) -> Option<String> {
    let Some(v) = version else {
        return Some("no `Program PWSCF v.X` banner, so the release could not be checked".into());
    };
    let major = v.split('.').next().unwrap_or(v);
    (!VERIFIED_MAJORS.contains(&major)).then(|| format!(
        "pw.x {v}: ferro has fixtures only from QE 7.x; spot-check one frame's energy \
         and forces before training on it"
    ))
}

/// Reads a pw.x output file, discarding the statistics.
pub fn read_qe_out(path: &Path) -> Result<Trajectory> {
    Ok(read_qe_out_with_stats(path)?.0)
}

/// Reads a pw.x output file and reports what was dropped.
pub fn read_qe_out_with_stats(path: &Path) -> Result<(Trajectory, AimdStats)> {
    let path_ = path.display();
    parse_qe_out_lines(|| super::util::text_lines(path)).with_context(|| format!("parsing {path_}"))
}

#[cfg(test)]
fn parse_qe_out(content: &str) -> Result<(Trajectory, AimdStats)> {
    parse_qe_out_lines(|| Ok(content.lines().map(|l| Ok(l.to_string()))))
}

/// Parses a pw.x output in two passes over its lines, never holding the whole file
/// (review C-1).
///
/// The first pass keeps no text: it finds the facts every step's verdict depends on
/// and that may be printed after the step — one `Program PWSCF` banner, no BFGS,
/// md or not, whether any step printed forces / stress. The second holds two steps'
/// worth of lines: a step's structure is printed at the end of the step before it,
/// so a step is parsed when the next `Self-consistent Calculation` arrives.
fn parse_qe_out_lines<I: Iterator<Item = Result<String>>>(
    open: impl Fn() -> Result<I>,
) -> Result<(Trajectory, AimdStats)> {
    let mut stats = AimdStats::new(AimdFormat::QeOut);

    let mut banner: Option<String> = None;
    let mut n_banners = 0usize;
    let (mut bfgs, mut is_md, mut has_forces, mut has_stress) = (false, false, false, false);
    for l in open()? {
        let l = l?;
        if l.contains("Program PWSCF") {
            n_banners += 1;
            banner.get_or_insert(l.clone());
        }
        bfgs |= l.contains("BFGS Geometry Optimization");
        // md 是 `Entering Dynamics:`，vc-md 是 `Entering Dynamics;`
        is_md |= l.trim_start().starts_with("Entering Dynamics");
        // 应力：别的步打印了而这一步没有，就是这一步没写完；整份文件都没有则是没开 tstress。
        // 力恒为必需（同 cp2k_sp），has_forces 只用于报错时点名 tprnfor
        has_forces |= starts_with_tokens(&l, &["Forces", "acting", "on", "atoms"]);
        has_stress |= starts_with_tokens(&l, &["total", "stress"]);
    }
    let banner = match (n_banners, banner) {
        (1, Some(b)) => b,
        (0, _) | (_, None) => bail!("no `Program PWSCF` banner; this is not a pw.x output"),
        (n, _) => bail!("{n} `Program PWSCF` banners: several pw.x runs in one file. Split it, one run per file"),
    };
    if bfgs {
        bail!("a relax / vc-relax output; ferro reads pw.x scf, md and vc-md");
    }
    stats.version = banner.split_whitespace()
        .find_map(|t| t.strip_prefix("v.")).map(|v| v.to_string());
    stats.version_note = version_note(stats.version.as_deref());

    let is_anchor = |l: &str| l.trim() == "Self-consistent Calculation";
    let mut steps = Steps { frames: Vec::new(), header: None, k: 0, is_md, has_stress };
    // Lines from the previous step's anchor (or the start of the file) on; `base` is
    // the file line number of window[0], `pending` the anchor whose step waits
    let mut window: Vec<String> = Vec::new();
    let mut base = 0usize;
    let mut pending: Option<usize> = None;
    for l in open()? {
        let l = l?;
        if is_anchor(&l) {
            let view: Vec<&str> = window.iter().map(|l| l.as_str()).collect();
            if steps.header.is_none() {
                steps.header = Some(read_header(&view)?);
            }
            if let Some(a) = pending.take() {
                steps.parse(&view, a, base, &mut stats)?;
                window.drain(..a);
                base += a;
            }
            pending = Some(window.len());
        }
        window.push(l);
    }
    if steps.header.is_none() {
        let view: Vec<&str> = window.iter().map(|l| l.as_str()).collect();
        read_header(&view)?;
        bail!("no `Self-consistent Calculation`: the run stopped before its first SCF");
    }
    if let Some(a) = pending {
        let view: Vec<&str> = window.iter().map(|l| l.as_str()).collect();
        steps.parse(&view, a, base, &mut stats)?;
    }
    stats.n_steps = steps.k;
    if is_md {
        stats.steps = Some((0, steps.k as i64 - 1));
    }
    let frames = steps.frames;

    // 同 cp2k_sp：一帧不剩就报错，collect 据此 SKIP 这份文件并说出原因，
    // 而不是拿一条空轨迹去和同目录的其他文件比组成
    if frames.is_empty() {
        bail!(
            "no usable frame: {} SCF step(s), {} dropped ({} SCF not converged, {} incomplete){}",
            stats.n_steps, stats.n_dropped(), stats.n_scf_failed, stats.n_incomplete,
            // 没收敛就停的那份也一行力都没有，那里提 tprnfor 是误导
            if !has_forces && stats.n_incomplete > 0 { "; no forces printed at all: set tprnfor = .true." } else { "" }
        );
    }
    let mut traj = Trajectory { frames, metadata: Default::default() };
    traj.metadata.program = Some(match &stats.version {
        Some(v) => format!("Quantum ESPRESSO {v} pw.x"),
        None => "Quantum ESPRESSO pw.x".to_string(),
    });
    Ok((traj, stats))
}

/// The steps parsed so far, and what the next one inherits from them.
struct Steps {
    frames: Vec<Frame>,
    /// `(nat, alat, cell, sites)`: the header's, then each step's own cell and sites
    header: Option<(usize, f64, Matrix3<f64>, Sites)>,
    /// Index of the next step
    k: usize,
    is_md: bool,
    has_stress: bool,
}

impl Steps {
    /// Parses step `self.k` from `lines`, the window from the previous anchor (or the
    /// start of the file) to the next one (or the end), with this step's anchor at `a`.
    /// `base` is the window's first line in the file, for messages.
    fn parse(&mut self, lines: &[&str], a: usize, base: usize, stats: &mut AimdStats) -> Result<()> {
        let k = self.k;
        self.k += 1;
        let (nat, alat, cell, atoms) = self.header.as_mut().expect("the header is read at the first anchor");
        let (nat, alat) = (*nat, *alat);
        let hi = lines.len();
        // 这一步的结构印在上一步的末尾（第 0 步用输出头）；先更新，后面丢帧也不影响下一步
        if k > 0 {
            let prev = 0..a;
            if let Some(c) = prev.clone().rev().find(|&i| lines[i].starts_with("CELL_PARAMETERS")) {
                *cell = read_cell_card(lines, c, base, alat)?;
            }
            match prev.rev().find(|&i| lines[i].starts_with("ATOMIC_POSITIONS")) {
                Some(p) => *atoms = read_positions_card(lines, p, base, nat, cell, alat)?,
                None => {
                    stats.n_incomplete += 1;
                    return Ok(());
                }
            }
        }
        let window = a..hi;
        if window.clone().any(|i| lines[i].contains("convergence NOT achieved")) {
            stats.n_scf_failed += 1;
            return Ok(());
        }
        let Some(energy) = window.clone()
            .find(|&i| lines[i].starts_with('!') && lines[i].contains("total energy"))
            .and_then(|i| lines[i].split('=').nth(1)?.split_whitespace().next()?.parse::<f64>().ok())
        else {
            stats.n_incomplete += 1;
            return Ok(());
        };
        let forces = match window.clone().find(|&i| starts_with_tokens(lines[i], &["Forces", "acting", "on", "atoms"])) {
            Some(f) => read_forces(lines, f, hi, nat),
            None => None,
        };
        let stress = window.clone()
            .find(|&i| starts_with_tokens(lines[i], &["total", "stress"]))
            .and_then(|s| read_stress(lines, s));
        if forces.is_none() || (self.has_stress && stress.is_none()) {
            stats.n_incomplete += 1;
            return Ok(());
        }

        let mut frame = Frame::with_cell(Cell::from_matrix(*cell), [true; 3]);
        frame.set_atoms(atoms.iter().map(|(label, pos)| {
            let mut atom = Atom::new(element_from_label(label), *pos);
            // 物种名与元素相同时不填 label（同 cp2k_sp）
            if atom.element != *label {
                atom.label = Some(label.clone());
            }
            atom
        }));
        frame.energy = Some(energy * RY_TO_EV);
        frame.set_forces(forces)?;
        frame.stress = stress;
        if self.is_md {
            frame.step = Some(k as i64);
        }
        self.frames.push(frame);
        stats.n_kept += 1;
        Ok(())
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::write_tmp;

    fn fixture(name: &str) -> String {
        std::fs::read_to_string(format!("../tests/{name}")).unwrap()
    }

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    /// Energy and forces agree with dpdata 1.0.2 (CODATA 2018, as pw.x 7.5 itself);
    /// the force must be the first group, not one of the six component groups after it.
    #[test]
    fn a_single_point_reads_energy_total_forces_and_stress() {
        let (traj, st) = parse_qe_out(&fixture("qe75_sp.out")).unwrap();
        assert_eq!((traj.frames.len(), st.n_steps, st.n_kept), (1, 1, 1));
        assert_eq!(st.version.as_deref(), Some("7.5"));
        assert!(st.version_note.is_none() && st.steps.is_none());
        let f = &traj.frames[0];
        assert!(close(f.energy.unwrap(), -3567.569297968025, 1e-8), "dpdata 的能量：{:?}", f.energy);
        let force = f.forces().unwrap()[0] / (RY_TO_EV / BOHR_TO_ANG);
        assert!(close(force.x, -0.00284143, 1e-12) && close(force.y, -0.00306774, 1e-12),
            "取总力（第一组），不是非局域分量 0.30871553：{force:?}");
        // 正 = 受压；输出里 P = -15.51 kbar（受拉），由 Ry/bohr³ 列算出的压力须与之相符
        let s = f.stress.unwrap();
        let p_kbar = s.trace() / 3.0 * 1602.1766208;
        assert!(close(p_kbar, -15.51, 0.01), "应力不翻号、单位对：P = {p_kbar}");
        assert!(close(s[(0, 1)], 0.00004574 * RY_TO_EV / BOHR_TO_ANG.powi(3), 1e-12));
        assert_eq!(f.symbols(), ["O", "O", "O", "O", "Al", "P"]);
        assert!(f.atoms().all(|a| a.label.is_none()), "物种名就是元素名时不填 label");
        assert_eq!(f.pbc, [true; 3]);
    }

    /// Frame 0 comes from the header (alat units, celldm(1) to six decimals);
    /// pw.in holds the same structure to ten.
    #[test]
    fn the_header_structure_matches_the_input_to_a_few_microangstroms() {
        let (traj, _) = parse_qe_out(&fixture("qe75_sp.out")).unwrap();
        let f = &traj.frames[0];
        let cell = f.cell.as_ref().unwrap().matrix;
        let input = Matrix3::new(4.2889663600, 0.0, 0.0, -0.3278899130, 4.2764144610, 0.0, 0.0, 0.0, 2.7079830000);
        assert!((cell - input).abs().max() < 5e-6, "晶胞行 = 晶格矢量，与 pw.in 一致：{cell}");
        let p = f.positions()[5];
        assert!((p - Vector3::new(1.9805382235, 2.1382072305, 1.3539915000)).abs().max() < 5e-6, "{p}");
    }

    /// Step k's structure is printed at the end of step k-1; the cards after the
    /// last step belong to a step that never ran.
    #[test]
    fn md_pairs_each_energy_with_the_structure_it_was_computed_for() {
        let (traj, st) = parse_qe_out(&fixture("qe75_md.out")).unwrap();
        assert_eq!((traj.frames.len(), st.n_steps, st.steps), (5, 5, Some((0, 4))));
        let e: Vec<f64> = traj.frames.iter().map(|f| f.energy.unwrap() / RY_TO_EV).collect();
        for (got, want) in e.iter().zip([-262.21150703, -262.21091702, -262.20919745, -262.20649405, -262.20288005]) {
            assert!(close(*got, want, 1e-9), "{e:?}");
        }
        let o1 = |k: usize| traj.frames[k].positions()[0];
        assert!((o1(0) - Vector3::new(0.8373002615, 0.9039570416, 1.3539915)).abs().max() < 5e-6, "第 0 帧取输出头");
        assert_eq!(o1(1), Vector3::new(0.8433996476, 0.9045494384, 1.3606785696),
            "第 1 帧取第 1 次 Entering Dynamics 之后的卡片");
        assert_eq!(o1(2), Vector3::new(0.8492833299, 0.9050318646, 1.3671255795));
        assert_eq!(traj.frames.iter().map(|f| f.step).collect::<Vec<_>>(), [0, 1, 2, 3, 4].map(Some));
        let m = |k: usize| traj.frames[k].cell.as_ref().unwrap().matrix;
        assert!((1..5).all(|k| m(k) == m(0)), "md 不变胞");
    }

    #[test]
    fn vc_md_takes_each_steps_cell_from_the_card_before_it() {
        let (traj, st) = parse_qe_out(&fixture("qe75_vcmd.out")).unwrap();
        assert_eq!((traj.frames.len(), st.n_kept), (5, 5));
        let row0 = |k: usize| traj.frames[k].cell.as_ref().unwrap().matrix.row(0).transpose();
        assert_eq!(row0(1), Vector3::new(4.288879115, 0.000027097, -0.000000005));
        assert!((row0(0) - Vector3::new(4.28896636, 0.0, 0.0)).abs().max() < 5e-6);
        let vol: Vec<f64> = traj.frames.iter().map(|f| f.cell.as_ref().unwrap().matrix.determinant()).collect();
        assert!(vol.windows(2).all(|w| w[1] < w[0]), "这条 vc-md 的体积逐步缩小：{vol:?}");
        assert!(close(traj.frames[4].energy.unwrap() / RY_TO_EV, -262.21155589, 1e-9));
    }

    #[test]
    fn spin_polarised_output_parses_the_same_way() {
        let (traj, st) = parse_qe_out(&fixture("qe75_sp_spin.out")).unwrap();
        assert_eq!(st.n_kept, 1);
        assert!(close(traj.frames[0].energy.unwrap(), -3561.1819324030753, 1e-8), "dpdata 的能量");
        assert_eq!(traj.frames[0].forces().unwrap().len(), 6);
    }

    #[test]
    fn an_unconverged_scf_is_counted_and_dropped() {
        // 唯一的一步没收敛：一帧不剩，报错并说出原因（同 cp2k_sp）
        let err = format!("{:#}", parse_qe_out(&fixture("qe75_sp_noconv.out")).unwrap_err());
        assert!(err.contains("no usable frame") && err.contains("1 SCF not converged"), "{err}");
        assert!(!err.contains("tprnfor"), "没收敛不是没开 tprnfor：{err}");
    }

    #[test]
    fn a_file_without_forces_gives_no_frame_and_names_tprnfor() {
        // 力是必需的（同 cp2k_sp）：没力的帧进不了 DeePMD system
        let text = fixture("qe75_sp.out").replace("Forces acting on atoms", "Forces gone");
        let err = format!("{:#}", parse_qe_out(&text).unwrap_err());
        assert!(err.contains("1 incomplete") && err.contains("tprnfor = .true."), "{err}");
    }

    /// The file is cut inside the last step's force block: that step is incomplete,
    /// the four before it are intact.
    #[test]
    fn a_run_cut_off_mid_step_keeps_the_finished_steps() {
        let text = fixture("qe75_md.out");
        let cut: String = text.lines().take(3257).map(|l| format!("{l}\n")).collect();
        let (traj, st) = parse_qe_out(&cut).unwrap();
        assert_eq!((traj.frames.len(), st.n_steps, st.n_incomplete), (4, 5, 1));
    }

    #[test]
    fn a_file_without_stress_gives_frames_without_stress() {
        // tstress 关掉时整份文件没有 total stress：帧照收，只是没有应力
        let text = fixture("qe75_sp.out").replace("total   stress", "no      stress");
        let (traj, st) = parse_qe_out(&text).unwrap();
        assert_eq!(st.n_kept, 1);
        assert!(traj.frames[0].stress.is_none() && traj.frames[0].forces().is_some());
    }

    #[test]
    fn relaxations_and_concatenated_runs_are_refused() {
        let relax = fixture("qe75_sp.out").replace("     Self-consistent Calculation",
            "     BFGS Geometry Optimization\n     Self-consistent Calculation");
        assert!(format!("{:#}", parse_qe_out(&relax).unwrap_err()).contains("relax"));
        let two = fixture("qe75_sp.out").repeat(2);
        assert!(format!("{:#}", parse_qe_out(&two).unwrap_err()).contains("2 `Program PWSCF` banners"));
    }

    #[test]
    fn cards_in_crystal_bohr_and_alat_units() {
        let cell = Matrix3::new(4.0, 0.0, 0.0, 1.0, 4.0, 0.0, 0.0, 0.0, 5.0);
        let card = ["ATOMIC_POSITIONS (crystal)", "O 0.5 0.5 0.5", "P 0 0 1.0"];
        let at = read_positions_card(&card, 0, 0, 2, &cell, 1.0).unwrap();
        assert_eq!(at[0].1, Vector3::new(2.5, 2.0, 2.5), "r = s·M（行 = 晶格矢量）");
        let card = ["ATOMIC_POSITIONS (bohr)", "O 1 0 0"];
        assert!(close(read_positions_card(&card, 0, 0, 1, &cell, 1.0).unwrap()[0].1.x, BOHR_TO_ANG, 1e-15));
        let card = ["CELL_PARAMETERS (alat= 2.00000000)", "1 0 0", "0 1 0", "0 0 1"];
        assert!(close(read_cell_card(&card, 0, 0, 99.0).unwrap()[(2, 2)], 2.0 * BOHR_TO_ANG, 1e-15),
            "卡片自带的 alat 优先于输出头的");
    }

    #[test]
    fn sniff_tells_pw_x_from_other_quantum_espresso_programs() {
        use super::super::aimd::sniff;
        assert_eq!(sniff(Path::new("../tests/qe75_md.out")).unwrap(), AimdFormat::QeOut);
        let gipaw = fixture("qe75_sp.out").replacen("Program PWSCF", "Program GIPAW", 1);
        let err = sniff(&write_tmp("ferro_qe_gipaw.out", &gipaw)).unwrap_err().to_string();
        assert!(err.contains("`GIPAW` output") && err.contains("pw.x"), "{err}");
    }

    /// A 396-atom scf with tprnfor / tstress off (QE 7.3.1): no frame, and the
    /// header still parses (alat, axes, 396 positions) up to that point.
    #[test]
    #[ignore = "reads examples/qe/, which is not in git"]
    fn an_energy_only_scf_from_another_release() {
        let text = std::fs::read_to_string("../examples/qe/5Zn-J17_scf.out").unwrap();
        let err = format!("{:#}", parse_qe_out(&text).unwrap_err());
        assert!(err.contains("tprnfor = .true."), "{err}");
        let lines: Vec<&str> = text.lines().collect();
        let (nat, _, cell, atoms) = read_header(&lines).unwrap();
        assert_eq!((nat, atoms.len()), (396, 396));
        assert!(close(cell[(0, 0)], 17.623106396720, 2e-5), "{cell}");
        for name in ["5Zn-J17_nmr.out", "pp.out"] {
            let err = super::super::aimd::sniff(Path::new(&format!("../examples/qe/{name}"))).unwrap_err();
            assert!(err.to_string().contains("ferro reads pw.x"), "{name}: {err}");
        }
    }
}
