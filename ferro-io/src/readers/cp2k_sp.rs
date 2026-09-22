//! CP2K single-point (`ENERGY` / `ENERGY_FORCE`) output reader.
//!
//! A single-point run prints no `MD| Step number`, so it shares no anchor with
//! [`super::cp2k_md`] — and none of its blocks either. Where the MD log dumps
//! coordinates and forces as two xyz blocks, a single point prints tables:
//!
//! ```text
//!  MODULE QUICKSTEP: ATOMIC COORDINATES IN ANGSTROM      <- coordinates
//!    Atom Kind Element         X             Y             Z       Z(eff)  Mass
//!       1    1 O     8      0.010480     15.694100     10.158000  6.0000  15.9994
//!  ...
//!  *** SCF run converged in     4 steps ***
//!  ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]     -8399.118581052152877  <- the anchor
//!  FORCES| Atomic forces [hartree/bohr]                 <- forces
//!  FORCES|   Atom     x               y               z               |f|
//!  FORCES|      1  1.22541206E-02 -1.80351825E-02  7.92236219E-03   2.31990324E-02
//!  ...
//!  STRESS| Analytical stress tensor [bar]               <- stress
//! ```
//!
//! # Frame layout
//!
//! The anchor is the `ENERGY|` line. Coordinates and cell are looked for
//! BACKWARDS from it, forces and stress FORWARDS, each bounded by the
//! neighbouring anchors. That ordering is CP2K's own, and anchoring this way
//! means a file made by concatenating N single-point outputs — which is how
//! people hand a batch of FP calculations to a dataset builder — reads as N
//! frames with no extra work.
//!
//! A file carrying two CP2K headers but one energy (a run that died and was
//! redone in place) therefore takes the SECOND header's coordinates: the ones
//! belonging to the energy that was actually reached.
//!
//! # What is dropped
//!
//! Frames whose SCF did not converge, whose coordinates / cell / forces are
//! missing or truncated, whose force count differs from the atom count, or
//! whose composition differs from the first frame. A dropped frame is counted
//! and skipped; it never aborts the file, because a batch of FP jobs always
//! contains a few that did not finish.
//!
//! Stress is the one optional quantity: a run without `STRESS_TENSOR` prints no
//! such block, and the frame is still perfectly good training data — it just
//! carries no virial.
//!
//! # Units
//!
//! Energy and stress carry their unit in the text and it is read from there.
//! Forces are printed in atomic units (Hartree/Bohr) in both block formats, the
//! new one saying so in its header and the old one not.
//!
//! The stress tensor keeps CP2K's own sign (positive = compression), matching
//! [`ferro_core::Frame::stress`].

use std::collections::HashMap;
use std::path::Path;

use anyhow::{bail, Context, Result};
use ferro_core::units::{convert_pressure, PressureUnit, BOHR_TO_ANG, HARTREE_TO_EV};
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};

use super::aimd::{AimdFormat, AimdStats};

/// Releases whose single-point log layout ferro has been run against.
///
/// 2025 is the fixture (`tests/5Al_0003_1500K_f394.out`); 2026 carried the 2025
/// layout forward unchanged. Older releases print a different force block and a
/// different stress block; both are handled below, but only ever checked
/// against cp2kdata's reading of them, so they earn a note.
const VERIFIED_MAJORS: &[&str] = &["2025", "2026"];

/// Says whether `version` is one ferro has a fixture for, and why not if not.
///
/// Travels out through [`AimdStats`] rather than being printed here: `collect`
/// reads dozens of files in a row and a library that prints cannot be silenced.
fn version_note(version: Option<&str>) -> Option<String> {
    let Some(v) = version else {
        return Some(
            "no `CP2K| version string` line, so the release could not be checked; \
             ferro's CP2K layout is verified for 2025 and 2026"
                .to_string(),
        );
    };
    let major = v.split('.').next().unwrap_or(v);
    if VERIFIED_MAJORS.contains(&major) {
        return None;
    }
    Some(format!(
        "CP2K {v}: ferro's single-point layout is verified for 2025 and 2026 only. \
         This file uses the older `ATOMIC FORCES` block, read with rules ported \
         from cp2kdata — spot-check one frame's forces before training on it"
    ))
}

/// True when the line's leading tokens match `tokens`.
fn starts_with_tokens(line: &str, tokens: &[&str]) -> bool {
    let mut it = line.split_whitespace();
    tokens.iter().all(|t| it.next() == Some(*t))
}

// 单位标注在方括号里:"STRESS| Analytical stress tensor [bar]" -> "bar"
fn unit_in_brackets(line: &str) -> Option<&str> {
    let start = line.find('[')?;
    let end = line[start..].find(']')? + start;
    Some(line[start + 1..end].trim())
}

// 能量行的单位标注在版本间换过两次括号:
//   <=7.1    ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):   -1766.2
//   8.1-2024 ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]:    -551.5
//   2025+    ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]   -75.6
// 取"数值前一个 token",不能取"第一个括号" —— 行里更早还有一个 `( QS )`
fn energy_unit(line: &str) -> Option<&str> {
    let toks: Vec<&str> = line.split_whitespace().collect();
    let t = toks.get(toks.len().checked_sub(2)?)?;
    Some(
        t.trim_end_matches(':')
            .trim_start_matches(['[', '('])
            .trim_end_matches([']', ')']),
    )
}

fn energy_to_ev(line: &str) -> Result<f64> {
    let factor = match energy_unit(line).unwrap_or("").to_ascii_lowercase().as_str() {
        "hartree" | "a.u." | "au" => HARTREE_TO_EV,
        "ev" => 1.0,
        other => bail!("unknown energy unit `{other}` in: {}", line.trim()),
    };
    let v: f64 = line
        .split_whitespace()
        .last()
        .context("energy line has no value")?
        .parse()
        .with_context(|| format!("cannot parse energy in: {}", line.trim()))?;
    Ok(v * factor)
}

fn stress_unit(line: &str) -> Result<PressureUnit> {
    match unit_in_brackets(line).unwrap_or("").to_ascii_lowercase().as_str() {
        "bar" => Ok(PressureUnit::Bar),
        "gpa" => Ok(PressureUnit::GPa),
        "kbar" => Ok(PressureUnit::Kbar),
        other => bail!(
            "unknown stress unit `[{other}]`; STRESS_UNIT is a CP2K input keyword, \
             so ferro reads it from the text rather than guessing"
        ),
    }
}

/// One row of the coordinate table: `idx kind symbol Z x y z Z(eff) mass`.
///
/// Nine fields exactly. The element column is really two — the symbol and the
/// atomic number — which is why `Kind` cannot simply be read as "the column
/// before the symbol" in the force table below.
fn coord_row(line: &str) -> Option<(usize, String, Vector3<f64>)> {
    let f: Vec<&str> = line.split_whitespace().collect();
    if f.len() != 9 {
        return None;
    }
    f[0].parse::<usize>().ok()?;
    let kind: usize = f[1].parse().ok()?;
    if !f[2].chars().all(|c| c.is_ascii_alphabetic()) {
        return None;
    }
    f[3].parse::<u32>().ok()?;
    let x: f64 = f[4].parse().ok()?;
    let y: f64 = f[5].parse().ok()?;
    let z: f64 = f[6].parse().ok()?;
    Some((kind, f[2].to_string(), Vector3::new(x, y, z)))
}

/// Reads the coordinate table whose header is at `head`.
///
/// Returns the kind index, element symbol and position of every atom. The rows
/// do not start at a fixed offset — 6.1 puts a blank line between the column
/// header and the first row and 2025 does not — so the search is for the first
/// line that parses as a row, then rows until one does not.
fn read_coord_table(lines: &[&str], head: usize, hi: usize) -> Option<Vec<(usize, String, Vector3<f64>)>> {
    let mut out = Vec::new();
    for l in lines.iter().take(hi).skip(head + 1) {
        match coord_row(l) {
            Some(r) => out.push(r),
            // 还没开始收就继续找表头下面的空行与列名行;已经收过就说明表结束了
            None if out.is_empty() => continue,
            None => break,
        }
    }
    (!out.is_empty()).then_some(out)
}

/// `1. Atomic kind: Fe1   Number of atoms:  6` -> `{1: "Fe1"}`.
///
/// The kind NAME is what the CP2K input called the species and it need not be
/// the element: `Fe1`/`Fe2` for two spin guesses, or an outright different
/// symbol when a site was substituted. It is kept as [`Atom::label`] so that
/// both readings survive; see `docs/src/data-model.md`.
fn read_kind_names(lines: &[&str], lo: usize, hi: usize) -> HashMap<usize, String> {
    let mut out = HashMap::new();
    for l in lines.iter().take(hi).skip(lo) {
        let f: Vec<&str> = l.split_whitespace().collect();
        if f.len() < 4 || f[1] != "Atomic" || f[2] != "kind:" {
            continue;
        }
        if let Ok(idx) = f[0].trim_end_matches('.').parse::<usize>() {
            out.insert(idx, f[3].to_string());
        }
    }
    out
}

/// Reads the three `CELL| Vector` rows into a row-major matrix.
///
/// Matches `CELL|` exactly: a 2025 output also prints `CELL_TOP|` and
/// `CELL_REF|` blocks with the same shape, and `starts_with("CELL")` would take
/// whichever came last.
///
/// The vector components are printed to three decimals while `|a|` gets six.
/// Rebuilding from lengths and angles would be more precise but would force the
/// a axis onto x, and releases up to 7.1 print only three decimals for the
/// lengths too, so there is nothing to rebuild from there. The resulting volume
/// error is ~1e-4 relative, i.e. <4e-3 eV of virial at 20 GPa — well under the
/// error a model trained on it would carry anyway.
fn read_cell(lines: &[&str], lo: usize, hi: usize) -> Option<(Matrix3<f64>, [bool; 3])> {
    let mut m = [[0.0_f64; 3]; 3];
    let mut seen = [false; 3];
    let mut pbc = [true; 3];
    for l in lines.iter().take(hi).skip(lo) {
        let f: Vec<&str> = l.split_whitespace().collect();
        if f.first() != Some(&"CELL|") {
            continue;
        }
        if f.get(1) == Some(&"Periodicity") {
            let p = f.get(2).copied().unwrap_or("XYZ").to_ascii_uppercase();
            pbc = if p == "NONE" {
                [false; 3]
            } else {
                [p.contains('X'), p.contains('Y'), p.contains('Z')]
            };
            continue;
        }
        if f.get(1) != Some(&"Vector") {
            continue;
        }
        let row = match f.get(2) {
            Some(&"a") => 0,
            Some(&"b") => 1,
            Some(&"c") => 2,
            _ => continue,
        };
        // 分量在 `[angstrom]:` 之后。2025 的 CELL_TOP| 那行少一个 `]:`,
        // 按 token 找而不是按固定列号,两种写法都对
        let at = f.iter().position(|t| t.starts_with("[angstrom"))? + 1;
        for (c, slot) in m[row].iter_mut().enumerate() {
            *slot = f.get(at + c)?.parse().ok()?;
        }
        seen[row] = true;
    }
    seen.iter().all(|s| *s).then(|| {
        let flat: Vec<f64> = m.iter().flatten().copied().collect();
        (Matrix3::from_row_slice(&flat), pbc)
    })
}

/// Reads the force block starting at `lo`, in whichever of the two formats.
///
/// Returns forces in eV/Å. Both formats print Hartree/Bohr; only the 2025 one
/// says so in its header.
fn read_forces(lines: &[&str], lo: usize, hi: usize) -> Option<Vec<Vector3<f64>>> {
    let factor = HARTREE_TO_EV / BOHR_TO_ANG;
    let mut out: Vec<Vector3<f64>> = Vec::new();
    for l in lines.iter().take(hi).skip(lo) {
        let f: Vec<&str> = l.split_whitespace().collect();
        // 2025+: `FORCES|  <i>  <x> <y> <z> <|f|>`,块尾是 `FORCES| Sum`
        let v = if f.first() == Some(&"FORCES|") {
            if f.len() != 6 || f[1].parse::<usize>().is_err() {
                if out.is_empty() { continue } else { break }
            }
            Some((f[2], f[3], f[4]))
        // <=2024: `<i> <kind> <symbol> <x> <y> <z>`,块尾是 `SUM OF ATOMIC FORCES`
        } else if f.len() == 6
            && f[0].parse::<usize>().is_ok()
            && f[1].parse::<usize>().is_ok()
            && f[2].chars().all(|c| c.is_ascii_alphabetic())
        {
            Some((f[3], f[4], f[5]))
        } else {
            None
        };
        match v {
            Some((x, y, z)) => {
                let (x, y, z) = (x.parse::<f64>().ok()?, y.parse::<f64>().ok()?, z.parse::<f64>().ok()?);
                out.push(Vector3::new(x, y, z) * factor);
            }
            None if out.is_empty() => continue,
            None => break,
        }
    }
    (!out.is_empty()).then_some(out)
}

/// Reads the stress tensor whose header is at `head`, in eV/Å³.
///
/// Rows are `STRESS| x  <xx> <xy> <xz>` from 8.1 on and `X  <xx> <xy> <xz>`
/// before it. Both are accepted, but the row must begin with `STRESS|` or with
/// the row label itself — the eigenvector rows that follow the block are three
/// bare floats and would otherwise be read as the matrix.
fn read_stress(lines: &[&str], head: usize, hi: usize) -> Result<Option<Matrix3<f64>>> {
    let unit = stress_unit(lines[head])?;
    let mut t = Matrix3::zeros();
    let mut row = 0usize;
    for l in lines.iter().take(hi.min(head + 12)).skip(head + 1) {
        if row == 3 {
            break;
        }
        let f: Vec<&str> = l.split_whitespace().collect();
        let labelled = f
            .first()
            .is_some_and(|t| matches!(*t, "X" | "Y" | "Z" | "x" | "y" | "z"));
        if f.first() != Some(&"STRESS|") && !labelled {
            continue;
        }
        if f.len() < 4 {
            continue;
        }
        // 取末尾三个:前面是 `STRESS|` 加行标,或只有行标,列数变了都不影响。
        // 表头行与 `1/3 Trace` 摘要行解析不出三个浮点,于是自动跳过
        let vals: Option<Vec<f64>> = f[f.len() - 3..].iter().map(|x| x.parse::<f64>().ok()).collect();
        if let Some(v) = vals {
            for (col, x) in v.iter().enumerate() {
                t[(row, col)] = convert_pressure(*x, unit, PressureUnit::EVPerAng3);
            }
            row += 1;
        }
    }
    Ok((row == 3).then_some(t))
}

/// Reads a CP2K single-point output file, discarding the statistics.
pub fn read_cp2k_sp(path: &Path) -> Result<Trajectory> {
    Ok(read_cp2k_sp_with_stats(path)?.0)
}

/// Reads a CP2K single-point output file and reports what was dropped.
pub fn read_cp2k_sp_with_stats(path: &Path) -> Result<(Trajectory, AimdStats)> {
    let path_ = path.display();
    let content = std::fs::read_to_string(path).with_context(|| format!("cannot open {path_}"))?;
    parse_cp2k_sp(&content).with_context(|| format!("parsing {path_}"))
}

fn parse_cp2k_sp(content: &str) -> Result<(Trajectory, AimdStats)> {
    let lines: Vec<&str> = content.lines().collect();
    let mut stats = AimdStats::new(AimdFormat::Cp2kSp);

    let mut anchors: Vec<usize> = Vec::new();
    let mut converged: Vec<usize> = Vec::new();
    let mut version: Option<String> = None;
    let mut run_type: Option<String> = None;
    for (i, l) in lines.iter().enumerate() {
        if starts_with_tokens(l, &["ENERGY|", "Total", "FORCE_EVAL"])
            || starts_with_tokens(l, &["ENERGY|", "Total", "force_eval"])
        {
            anchors.push(i);
        } else if l.contains("SCF run converged") {
            converged.push(i);
        } else if version.is_none() && starts_with_tokens(l, &["CP2K|", "version", "string:"]) {
            version = l.split_whitespace().last().map(|v| v.to_string());
        } else if run_type.is_none() && starts_with_tokens(l, &["GLOBAL|", "Run", "type"]) {
            run_type = l.split_whitespace().last().map(|v| v.to_string());
        }
    }

    match run_type.as_deref() {
        Some("ENERGY") | Some("ENERGY_FORCE") | None => {}
        Some(other) => bail!(
            "`GLOBAL| Run type {other}` is not a single point. ferro reads CP2K \
             ENERGY and ENERGY_FORCE here and MD through the MD reader; {other} \
             is not supported"
        ),
    }
    if anchors.is_empty() {
        bail!("no `ENERGY| Total FORCE_EVAL` line found; this is not a CP2K single-point output");
    }

    stats.version_note = version_note(version.as_deref());
    stats.version = version.clone();
    stats.n_steps = anchors.len();

    let mut frames: Vec<Frame> = Vec::with_capacity(anchors.len());
    let mut ref_comp: Option<Vec<String>> = None;

    for (k, &a) in anchors.iter().enumerate() {
        // 帧的前半在锚点之前,下界是上一个锚点;后半在锚点之后,上界是下一个
        let lo = if k == 0 { 0 } else { anchors[k - 1] };
        let hi = anchors.get(k + 1).copied().unwrap_or(lines.len());

        // SCF 未收敛时 CP2K 照样打能量和力,所以锚点在、数据不能要
        if !converged.iter().any(|&c| c > lo && c < a) {
            stats.n_scf_failed += 1;
            continue;
        }

        // 坐标表与 kind 名表都在锚点之前,取最近的那一份
        let Some(head) = (lo..a)
            .rev()
            .find(|&i| starts_with_tokens(lines[i], &["MODULE", "QUICKSTEP:", "ATOMIC", "COORDINATES"]))
        else {
            stats.n_incomplete += 1;
            continue;
        };
        let Some(rows) = read_coord_table(&lines, head, a) else {
            stats.n_incomplete += 1;
            continue;
        };
        let Some((matrix, pbc)) = read_cell(&lines, lo, a) else {
            stats.n_incomplete += 1;
            continue;
        };

        // 力块紧跟在能量之后;缺力的帧进不了数据集,force.npy 是必需项
        let Some(fhead) = (a..hi).find(|&i| {
            starts_with_tokens(lines[i], &["FORCES|", "Atomic", "forces"])
                || starts_with_tokens(lines[i], &["ATOMIC", "FORCES", "in"])
        }) else {
            stats.n_incomplete += 1;
            continue;
        };
        let Some(forces) = read_forces(&lines, fhead + 1, hi) else {
            stats.n_incomplete += 1;
            continue;
        };
        if forces.len() != rows.len() {
            // 截断通常发生在力表中途:原子数对不上就是这一帧没写完
            stats.n_incomplete += 1;
            continue;
        }

        let mut comp: Vec<String> = rows.iter().map(|(_, s, _)| s.clone()).collect();
        comp.sort();
        match &ref_comp {
            None => ref_comp = Some(comp),
            Some(r) if &comp != r => {
                stats.n_bad_composition += 1;
                continue;
            }
            Some(_) => {}
        }

        // 应力可缺:没开 STRESS_TENSOR 的单点照样是好数据,只是没有 virial
        let stress = match (a..hi).find(|&i| {
            starts_with_tokens(lines[i], &["STRESS|", "Analytical", "stress", "tensor"])
                || starts_with_tokens(lines[i], &["STRESS", "TENSOR"])
        }) {
            Some(s) => read_stress(&lines, s, hi)?,
            None => None,
        };

        let kinds = read_kind_names(&lines, lo, a);
        let atoms: Vec<Atom> = rows
            .iter()
            .map(|(kind, sym, pos)| {
                let mut atom = Atom::new(sym.as_str(), *pos);
                // kind 名与元素相同时不填 label —— 绝大多数体系如此,
                // 填了只会让每个原子都带一个复述元素符号的字段
                if let Some(name) = kinds.get(kind).filter(|n| n.as_str() != sym.as_str()) {
                    atom.label = Some(name.clone());
                }
                atom
            })
            .collect();

        let mut frame = Frame::with_cell(Cell::from_matrix(matrix), pbc);
        frame.atoms = atoms;
        frame.energy = Some(energy_to_ev(lines[a])?);
        frame.forces = Some(forces);
        frame.stress = stress;
        frames.push(frame);
        stats.n_kept += 1;
    }

    if frames.is_empty() {
        bail!(
            "no usable frame: {} energy block(s), {} dropped ({} SCF, {} incomplete, {} composition)",
            stats.n_steps, stats.n_dropped(), stats.n_scf_failed,
            stats.n_incomplete, stats.n_bad_composition
        );
    }

    let mut traj = Trajectory { frames, metadata: Default::default() };
    traj.metadata.source = Some(match version {
        Some(v) => format!("CP2K {v} single point"),
        None => "CP2K single point".to_string(),
    });
    Ok((traj, stats))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 2025-layout single point of two atoms, one kind renamed.
    const MINI: &str = concat!(
        " CP2K| version string:                                       CP2K version 2025.2\n",
        " GLOBAL| Run type                                                   ENERGY_FORCE\n",
        "  1. Atomic kind: Si                                   Number of atoms:       1\n",
        "  2. Atomic kind: O_b                                  Number of atoms:       1\n",
        " CELL| Volume [angstrom^3]:                                           125.000000\n",
        " CELL| Vector a [angstrom]:       5.000     0.000     0.000   |a| =     5.000000\n",
        " CELL| Vector b [angstrom]:       0.000     5.000     0.000   |b| =     5.000000\n",
        " CELL| Vector c [angstrom]:       0.000     0.000     5.000   |c| =     5.000000\n",
        " CELL| Periodicity                                                           XYZ\n",
        "\n",
        " MODULE QUICKSTEP: ATOMIC COORDINATES IN ANGSTROM\n",
        "\n",
        "   Atom Kind Element         X             Y             Z       Z(eff)     Mass\n",
        "      1    1 Si   14      0.000000      0.000000      0.000000   4.0000  28.0855\n",
        "      2    2 O     8      1.000000      0.000000      0.000000   6.0000  15.9994\n",
        "\n",
        " *** SCF run converged in     4 steps ***\n",
        "\n",
        " ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]            -10.000000000000000\n",
        "\n",
        " FORCES| Atomic forces [hartree/bohr]\n",
        " FORCES|   Atom     x               y               z               |f|\n",
        " FORCES|      1  1.00000000E-02  0.00000000E+00  0.00000000E+00   1.00000000E-02\n",
        " FORCES|      2 -1.00000000E-02  0.00000000E+00  0.00000000E+00   1.00000000E-02\n",
        " FORCES| Sum     0.00000000E+00  0.00000000E+00  0.00000000E+00\n",
        "\n",
        " STRESS| Analytical stress tensor [bar]\n",
        " STRESS|                        x                   y                   z\n",
        " STRESS|      x        1.00000000000E+04   0.00000000000E+00   0.00000000000E+00\n",
        " STRESS|      y        0.00000000000E+00   2.00000000000E+04   0.00000000000E+00\n",
        " STRESS|      z        0.00000000000E+00   0.00000000000E+00   3.00000000000E+04\n",
        " STRESS| 1/3 Trace                                            2.00000000000E+04\n",
    );

    #[test]
    fn reads_energy_forces_cell_and_stress() {
        let (traj, st) = parse_cp2k_sp(MINI).unwrap();
        assert_eq!(traj.n_frames(), 1);
        assert_eq!(st.n_steps, 1);
        assert_eq!(st.n_dropped(), 0);
        assert_eq!(st.version.as_deref(), Some("2025.2"));
        assert_eq!(st.version_note, None, "2025 已验证,不该有提示");

        let f = &traj.frames[0];
        assert_eq!(f.n_atoms(), 2);
        assert!((f.energy.unwrap() - -10.0 * HARTREE_TO_EV).abs() < 1e-9);
        assert_eq!(f.cell.as_ref().unwrap().lengths(), [5.0, 5.0, 5.0]);
        assert_eq!(f.pbc, [true; 3]);

        let fx = f.forces.as_ref().unwrap()[0].x;
        assert!((fx - 1e-2 * HARTREE_TO_EV / BOHR_TO_ANG).abs() < 1e-12, "{fx}");

        // 1e4 bar = 1 GPa,正 = 压缩,不变号
        let s = f.stress.unwrap();
        assert!((s[(0, 0)] * 160.217_663_4 - 1.0).abs() < 1e-6, "{s:?}");
        assert!((s[(2, 2)] * 160.217_663_4 - 3.0).abs() < 1e-6, "{s:?}");
        assert_eq!(s[(0, 1)], 0.0);
    }

    /// The kind name is kept as a label only when it is not the element.
    #[test]
    fn a_renamed_kind_becomes_a_label_and_never_replaces_the_element() {
        let (traj, _) = parse_cp2k_sp(MINI).unwrap();
        let atoms = &traj.frames[0].atoms;
        assert_eq!(atoms[0].element, "Si");
        assert_eq!(atoms[0].label, None, "kind 名与元素相同时不该填 label");
        assert_eq!(atoms[1].element, "O", "element 必须是真元素,不是 kind 名");
        assert_eq!(atoms[1].label.as_deref(), Some("O_b"));
    }

    /// Concatenating N outputs is how a batch of FP jobs is handed over.
    #[test]
    fn a_concatenated_file_reads_as_one_frame_per_energy_block() {
        let two = format!("{MINI}{}", MINI.replace("-10.000000", "-11.000000"));
        let (traj, st) = parse_cp2k_sp(&two).unwrap();
        assert_eq!(traj.n_frames(), 2);
        assert_eq!(st.n_steps, 2);
        let (a, b) = (traj.frames[0].energy.unwrap(), traj.frames[1].energy.unwrap());
        assert!((a - -10.0 * HARTREE_TO_EV).abs() < 1e-9, "第一帧该是 -10 Ha: {a}");
        assert!((b - -11.0 * HARTREE_TO_EV).abs() < 1e-9, "第二帧该是 -11 Ha: {b}");
    }

    /// An unconverged frame is dropped and counted, never silently kept.
    #[test]
    fn an_unconverged_frame_is_dropped_not_kept() {
        let bad = MINI.replace("*** SCF run converged in     4 steps ***", "*** SCF run NOT converged ***");
        let err = parse_cp2k_sp(&bad).unwrap_err().to_string();
        assert!(err.contains("1 SCF"), "{err}");

        // 但同一批里收敛的那份照常留下,不因为有坏帧就整份作废
        let mixed = format!("{bad}{MINI}");
        let (traj, st) = parse_cp2k_sp(&mixed).unwrap();
        assert_eq!(traj.n_frames(), 1);
        assert_eq!(st.n_scf_failed, 1);
    }

    /// A run without STRESS_TENSOR is good training data, minus the virial.
    #[test]
    fn a_missing_stress_block_is_not_a_missing_frame() {
        let cut = MINI.split(" STRESS| Analytical").next().unwrap();
        let (traj, st) = parse_cp2k_sp(cut).unwrap();
        assert_eq!(traj.n_frames(), 1);
        assert_eq!(st.n_dropped(), 0);
        assert!(traj.frames[0].stress.is_none());
        assert!(traj.frames[0].forces.is_some(), "缺应力不该连带丢掉力");
    }

    /// Truncation lands in the force table; the atom count is what catches it.
    #[test]
    fn a_force_table_shorter_than_the_atom_count_drops_the_frame() {
        let cut = MINI.replace(
            " FORCES|      2 -1.00000000E-02  0.00000000E+00  0.00000000E+00   1.00000000E-02\n",
            "",
        );
        let err = parse_cp2k_sp(&cut).unwrap_err().to_string();
        assert!(err.contains("1 incomplete"), "{err}");
    }

    /// `CELL_TOP|` and `CELL_REF|` print the same shape and must not be read.
    #[test]
    fn the_cell_top_and_cell_ref_blocks_are_not_mistaken_for_the_cell() {
        let noisy = MINI.replace(
            " CELL| Volume [angstrom^3]:                                           125.000000\n",
            concat!(
                " CELL_TOP| Vector a [angstrom     9.000     0.000     0.000   |a| =     9.000000\n",
                " CELL_TOP| Vector b [angstrom     0.000     9.000     0.000   |b| =     9.000000\n",
                " CELL_TOP| Vector c [angstrom     0.000     0.000     9.000   |c| =     9.000000\n",
            ),
        );
        let (traj, _) = parse_cp2k_sp(&noisy).unwrap();
        assert_eq!(
            traj.frames[0].cell.as_ref().unwrap().lengths(),
            [5.0, 5.0, 5.0],
            "读到 9.0 说明抓了 CELL_TOP|"
        );
    }

    /// The old `ATOMIC FORCES` block of releases up to 2024.
    #[test]
    fn the_old_force_table_is_read_and_the_release_is_flagged() {
        let old = MINI
            .replace("CP2K version 2025.2", "CP2K version 7.1")
            .replace(
                "energy [hartree]            -10.000000000000000",
                "energy (a.u.):              -10.000000000000000",
            )
            .replace(
                concat!(
                    " FORCES| Atomic forces [hartree/bohr]\n",
                    " FORCES|   Atom     x               y               z               |f|\n",
                    " FORCES|      1  1.00000000E-02  0.00000000E+00  0.00000000E+00   1.00000000E-02\n",
                    " FORCES|      2 -1.00000000E-02  0.00000000E+00  0.00000000E+00   1.00000000E-02\n",
                    " FORCES| Sum     0.00000000E+00  0.00000000E+00  0.00000000E+00\n",
                ),
                concat!(
                    " ATOMIC FORCES in [a.u.]\n",
                    "\n",
                    " # Atom   Kind   Element          X              Y              Z\n",
                    "      1      1      Si          0.01000000     0.00000000     0.00000000\n",
                    "      2      2      O          -0.01000000     0.00000000     0.00000000\n",
                    " SUM OF ATOMIC FORCES           0.00000000     0.00000000     0.00000000\n",
                ),
            );
        let (traj, st) = parse_cp2k_sp(&old).unwrap();
        assert_eq!(traj.n_frames(), 1);
        let note = st.version_note.expect("7.1 未验证,该有提示");
        assert!(note.contains("7.1"), "{note}");
        // 老格式的力与新格式一字不差,单位换算两条路必须一致
        let (new, _) = parse_cp2k_sp(MINI).unwrap();
        for (a, b) in traj.frames[0].forces.as_ref().unwrap().iter()
            .zip(new.frames[0].forces.as_ref().unwrap())
        {
            assert!((a - b).norm() < 1e-12, "{a:?} vs {b:?}");
        }
        // 圆括号的能量也得读出来
        assert!((traj.frames[0].energy.unwrap() - new.frames[0].energy.unwrap()).abs() < 1e-12);
    }

    /// An MD output must not be read here — it is a different layout entirely.
    #[test]
    fn an_md_run_type_is_refused_by_name() {
        let md = MINI.replace("ENERGY_FORCE", "MD");
        let err = parse_cp2k_sp(&md).unwrap_err().to_string();
        assert!(err.contains("MD") && err.contains("single point"), "{err}");
    }

    #[test]
    fn reads_the_reference_single_point() {
        let (traj, st) =
            read_cp2k_sp_with_stats(Path::new("../tests/5Al_0003_1500K_f394.out")).unwrap();
        assert_eq!(traj.n_frames(), 1);
        assert_eq!(st.n_dropped(), 0);
        assert_eq!(st.version.as_deref(), Some("2025.2"));
        assert_eq!(st.version_note, None);

        let f = &traj.frames[0];
        assert_eq!(f.n_atoms(), 457, "303 O + 96 P + 48 Zn + 10 Al");
        assert_eq!(f.count_element("O"), 303);
        assert_eq!(f.count_element("Al"), 10);
        // 这份体系的 kind 名就是元素,没有一个原子该带 label
        assert!(f.atoms.iter().all(|a| a.label.is_none()));

        assert!((f.energy.unwrap() - -8_399.118_581_052_153 * HARTREE_TO_EV).abs() < 1e-6);
        assert_eq!(f.cell.as_ref().unwrap().lengths().map(|v| (v * 1e3).round()), [18742.0; 3]);
        assert_eq!(f.forces.as_ref().unwrap().len(), 457);
        // STRESS| ... [bar] 对角三个都是负的(拉伸),符号照抄 CP2K 不变
        let s = f.stress.unwrap();
        assert!(s[(0, 0)] < 0.0 && s[(1, 1)] < 0.0 && s[(2, 2)] < 0.0, "{s:?}");
        assert!((s[(0, 1)] - s[(1, 0)]).abs() < 1e-12, "应力张量必须对称");
    }
}
