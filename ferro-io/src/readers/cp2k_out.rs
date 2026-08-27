//! CP2K molecular-dynamics output (`.out`) reader.
//!
//! Reads the run log CP2K writes to stdout when the trajectory, the forces and
//! the stress tensor are all printed to `__STD_OUT__` instead of to sibling
//! files. Such a file is self-contained — it carries coordinates, forces, cell
//! and energy for every MD step — which is exactly what building a
//! machine-learning training set needs.
//!
//! # Frame layout
//!
//! Every MD step writes, in this order:
//!
//! ```text
//!  ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]   -2059.593202975
//!  STRESS| Analytical stress tensor [bar]
//!  STRESS|              x            y            z
//!  STRESS|    x   2.5718E+03  -1.5646E+03   8.6671E+02
//!  STRESS|    y   ...
//!  STRESS|    z   ...
//!  MD| Step number                       1          <- the anchor
//!  MD| Time [fs]                  2.000000
//!      112                                          <- coordinates, Angstrom
//!  i =  1, time = 2.000, E = -2059.5932029753
//!  Al  -4.196  -1.527   6.015
//!  ...
//!      112                                          <- forces, atomic units
//!  i =  1, time = 2.000, E = -2059.5932029753
//!  Al  -0.0021  -0.0013  -0.0006
//!  ...
//!  1  2.000  11.41 0 0  0.11 11.41 0  -5.64 -5.83 8.02   1046.109
//! ```
//!
//! The two xyz blocks are byte-for-byte indistinguishable — same atom count,
//! same `i = …, time = …, E = …` comment, same element column. Only their
//! ORDER separates them, so the scan is bounded by the next anchor: a frame
//! whose blocks are missing is dropped rather than silently borrowing the next
//! frame's data.
//!
//! # What is dropped
//!
//! Frames whose SCF did not converge, whose blocks are truncated (a job killed
//! mid-step), or whose composition differs from the first frame. The last one
//! is a guard against block misalignment rather than a real change of system:
//! if a warning is printed inside an xyz block, the "element" column stops
//! holding element symbols and the composition check catches it.
//!
//! # Units
//!
//! Energy and stress carry their unit in the text (`[hartree]`, `[bar]`) and it
//! is read from there — CP2K's `STRESS_UNIT` is an input keyword, so the unit
//! is NOT a function of the code version and a version table cannot answer it.
//! An unknown unit is an error, never a default. Forces are the one quantity
//! printed with no unit at all; they are atomic units (Hartree/Bohr).
//!
//! The stress tensor is stored with CP2K's own sign (positive = compression),
//! matching [`ferro_core::Frame::stress`]. Note it is the POTENTIAL part only:
//! `MD| Pressure` additionally contains the kinetic term and must not be used
//! for a training set.

use anyhow::{bail, Context, Result};
use ferro_core::units::{convert_pressure, PressureUnit, BOHR_TO_ANG, HARTREE_TO_EV};
use super::aimd::{AimdFormat, AimdStats};
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};

/// Text anchors, as token sequences rather than literal substrings.
///
/// CP2K reformats its log between releases — column alignment shifts, a field
/// widens, an extra space appears. Matching `" MD| Step number"` literally ties
/// the parser to one release's whitespace; matching the token sequence
/// `["MD|", "Step", "number"]` survives any amount of padding, because
/// `split_whitespace` has already thrown the padding away.
///
/// Each tag is a LIST of acceptable token sequences. Supporting a release that
/// renamed something is then one more line here, not a new branch in the
/// scanner — which is the whole point of keeping them in one table.
mod tag {
    /// Start of one MD step; the frame anchor everything else hangs off.
    pub const FRAME: &[&[&str]] = &[&["MD|", "Step", "number"]];
    /// Total potential energy of the step.
    pub const ENERGY: &[&[&str]] = &[
        &["ENERGY|", "Total", "FORCE_EVAL"],
        &["ENERGY|", "Total", "force_eval"],
    ];
    /// Header of the analytical stress tensor block; carries the unit.
    pub const STRESS: &[&[&str]] = &[&["STRESS|", "Analytical", "stress", "tensor"]];
    /// Any `STRESS|` line — the numeric rows of the block are a subset.
    pub const STRESS_ROW: &[&[&str]] = &[&["STRESS|"]];
    /// SCF convergence verdict; the wording after this differs by method.
    pub const SCF: &[&[&str]] = &[&["SCF", "run"]];
    /// One per MD initialisation; more than one means the run was restarted.
    pub const RESTART: &[&[&str]] = &[&["MD_INI|", "MD", "initialization"]];
    /// Banner line carrying the CP2K release, kept for the trajectory metadata.
    pub const VERSION: &[&[&str]] = &[&["CP2K|", "version", "string:"]];
}

/// True when the line's leading tokens match any of `candidates`.
fn line_matches(line: &str, candidates: &[&[&str]]) -> bool {
    candidates.iter().any(|tokens| {
        let mut it = line.split_whitespace();
        tokens.iter().all(|t| it.next() == Some(*t))
    })
}

/// Per-file account of what was parsed and what was thrown away.
///
/// Frame dropping happens inside the reader because every criterion needs the
/// surrounding text (the SCF line above the anchor, the block line count, the
/// first frame's composition); the counts travel out so the caller can report
/// them instead of the reader printing behind its back.
/// Reads a CP2K MD output file, discarding the statistics.
pub fn read_cp2k_out(path: &str) -> Result<Trajectory> {
    Ok(read_cp2k_out_with_stats(path)?.0)
}

/// Reads a CP2K MD output file and reports what was dropped.
pub fn read_cp2k_out_with_stats(path: &str) -> Result<(Trajectory, AimdStats)> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path}"))?;
    parse_cp2k_out(&content).with_context(|| format!("parsing {path}"))
}

// 方括号里的单位标注: "ENERGY| ... [hartree]  -2059.5" -> "hartree"
fn unit_in_brackets(line: &str) -> Option<&str> {
    let start = line.find('[')?;
    let end = line[start..].find(']')? + start;
    Some(line[start + 1..end].trim())
}

fn energy_to_ev(line: &str) -> Result<f64> {
    let unit = unit_in_brackets(line).unwrap_or("");
    let factor = match unit.to_ascii_lowercase().as_str() {
        "hartree" | "a.u." | "au" => HARTREE_TO_EV,
        "ev" => 1.0,
        other => bail!("unknown energy unit `[{other}]` in: {}", line.trim()),
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
        other => bail!("unknown stress unit `[{other}]`; STRESS_UNIT is a CP2K input keyword, so ferro reads it from the text rather than guessing"),
    }
}

// 一行 xyz: "Al  -4.196  -1.527   6.015" -> ("Al", 三个数)
fn xyz_line(line: &str) -> Option<(&str, Vector3<f64>)> {
    let mut it = line.split_whitespace();
    let sym = it.next()?;
    let x: f64 = it.next()?.parse().ok()?;
    let y: f64 = it.next()?.parse().ok()?;
    let z: f64 = it.next()?.parse().ok()?;
    Some((sym, Vector3::new(x, y, z)))
}

/// Recognises an xyz block head structurally, without reading the comment line.
///
/// The comment CP2K writes (`i = 1, time = 2.000, E = -2059.59`) is not part of
/// the xyz format's contract and its wording has no guarantee across releases,
/// so it is only required NOT to parse as an atom row — which is what makes it
/// a comment. What must hold is the shape: an atom count, then a line that is
/// not data, then that many parsable atom rows.
fn block_head(lines: &[&str], i: usize) -> Option<usize> {
    let n: usize = lines.get(i)?.trim().parse().ok()?;
    if n == 0 {
        return None;
    }
    // 第二行必须是注释（即解析不成原子行），否则那个"数字"是数据的一部分
    if xyz_line(lines.get(i + 1)?).is_some() {
        return None;
    }
    let first = i + 2;
    if first + n > lines.len() {
        return None;
    }
    // 只验首尾两行；全部 n 行随后真正读取时才逐行解析，不在这里做两遍
    xyz_line(lines[first])?;
    xyz_line(lines[first + n - 1])?;
    Some(n)
}

/// Reads `n` atom rows starting at the block-head line, returning symbols and vectors.
fn read_xyz_block(lines: &[&str], head: usize, n: usize) -> Option<(Vec<String>, Vec<Vector3<f64>>)> {
    let first = head + 2;
    if first + n > lines.len() {
        return None;
    }
    let mut syms = Vec::with_capacity(n);
    let mut vecs = Vec::with_capacity(n);
    for line in &lines[first..first + n] {
        let (s, v) = xyz_line(line)?;
        syms.push(s.to_string());
        vecs.push(v);
    }
    Some((syms, vecs))
}

fn parse_cp2k_out(content: &str) -> Result<(Trajectory, AimdStats)> {
    let lines: Vec<&str> = content.lines().collect();

    let mut stats = AimdStats::new(AimdFormat::Cp2kOut);
    let mut anchors: Vec<usize> = Vec::new();
    let mut steps: Vec<i64> = Vec::new();
    let mut scf: Vec<(usize, bool)> = Vec::new();
    let mut version: Option<String> = None;
    for (i, l) in lines.iter().enumerate() {
        if line_matches(l, tag::FRAME) {
            anchors.push(i);
            // 行末是步号本身；解析不出就不记，锚点本身仍然有效
            if let Some(n) = l.split_whitespace().last().and_then(|t| t.parse::<i64>().ok()) {
                steps.push(n);
            }
        } else if line_matches(l, tag::SCF) {
            // 收敛与否看有没有 NOT 这个词，而不是整句措辞
            let not_converged = l.split_whitespace().any(|t| t == "NOT");
            scf.push((i, !not_converged));
        } else if line_matches(l, tag::RESTART) {
            stats.n_restarts += 1;
        } else if version.is_none() && line_matches(l, tag::VERSION) {
            // 行尾是版本号本身（`CP2K| version string: CP2K version 2024.1`）
            version = l.split_whitespace().last().map(|v| v.to_string());
        }
    }
    stats.n_steps = anchors.len();
    stats.steps = match (steps.first(), steps.last()) {
        (Some(&a), Some(&b)) => Some((a, b)),
        _ => None,
    };
    if anchors.is_empty() {
        bail!("no `MD| Step number` line found; this is not a CP2K MD output");
    }
    // 重启段数比初始化次数少一（第一次不是重启）
    stats.n_restarts = stats.n_restarts.saturating_sub(1);

    let mut frames: Vec<Frame> = Vec::with_capacity(anchors.len());
    let mut ref_comp: Option<Vec<String>> = None;
    let mut ref_offsets: Option<(usize, usize)> = None;

    for (k, &a) in anchors.iter().enumerate() {
        // SCF 收敛与否: 锚点上方最近的那一行
        let converged = match scf.partition_point(|&(pos, _)| pos < a) {
            0 => true, // 没有 SCF 行可依，不因此丢帧
            j => scf[j - 1].1,
        };
        if !converged {
            stats.n_scf_failed += 1;
            continue;
        }

        // 帧的前半在锚点之前，下界是上一个锚点，避免越界抓到别帧的量
        let lo = if k == 0 { 0 } else { anchors[k - 1] };
        let mut energy = None;
        let mut stress_at = None;
        for i in (lo..a).rev() {
            let l = lines[i];
            if energy.is_none() && line_matches(l, tag::ENERGY) {
                energy = Some(energy_to_ev(l)?);
            } else if stress_at.is_none() && line_matches(l, tag::STRESS) {
                stress_at = Some(i);
            }
            if energy.is_some() && stress_at.is_some() {
                break;
            }
        }

        // 帧的后半在锚点之后，上界是下一个锚点
        let hi = anchors.get(k + 1).copied().unwrap_or(lines.len());
        let mut heads: Vec<(usize, usize)> = Vec::new(); // (行号, 原子数)
        let mut i = a;
        while i < hi && heads.len() < 2 {
            if let Some(n) = block_head(&lines, i) {
                heads.push((i, n));
                i += n + 2;
            } else {
                i += 1;
            }
        }
        if heads.len() < 2 || heads[0].1 != heads[1].1 {
            stats.n_incomplete += 1;
            continue;
        }
        let natoms = heads[0].1;

        let Some((c_syms, coords)) = read_xyz_block(&lines, heads[0].0, natoms) else {
            stats.n_incomplete += 1;
            continue;
        };
        let Some((f_syms, forces)) = read_xyz_block(&lines, heads[1].0, natoms) else {
            stats.n_incomplete += 1;
            continue;
        };

        // 力块之后第一行非空行是 cell 行: step time + 9 分量 + 体积
        let mut cell_line = None;
        for l in lines.iter().take(hi).skip(heads[1].0 + 2 + natoms) {
            if !l.trim().is_empty() {
                cell_line = Some(*l);
                break;
            }
        }
        let Some(cell_fields) = cell_line.map(|l| l.split_whitespace().collect::<Vec<_>>()) else {
            stats.n_incomplete += 1;
            continue;
        };
        // 末位是体积，它之前的九个是晶胞；从尾部取，前缀多一列少一列都无所谓
        if cell_fields.len() < 10 {
            stats.n_incomplete += 1;
            continue;
        }
        let nine = &cell_fields[cell_fields.len() - 10..cell_fields.len() - 1];
        let mut m = [0.0_f64; 9];
        let mut ok = true;
        for (slot, s) in m.iter_mut().zip(nine) {
            match s.parse::<f64>() {
                Ok(v) => *slot = v,
                Err(_) => {
                    ok = false;
                    break;
                }
            }
        }
        if !ok {
            stats.n_incomplete += 1;
            continue;
        }

        // 组成必须与首帧一致；不一致多半是块错位而非体系真变了
        let mut comp: Vec<String> = c_syms.clone();
        comp.sort();
        let mut fcomp: Vec<String> = f_syms.clone();
        fcomp.sort();
        match &ref_comp {
            None => ref_comp = Some(comp.clone()),
            Some(r) => {
                if &comp != r || &fcomp != r {
                    stats.n_bad_composition += 1;
                    continue;
                }
            }
        }

        // 块相对锚点的偏移在同一 ensemble 下是固定的；漂移说明有额外输出插入
        let offsets = (heads[0].0 - a, heads[1].0 - a);
        match ref_offsets {
            None => ref_offsets = Some(offsets),
            Some(r) if r != offsets => stats.n_layout_drift += 1,
            Some(_) => {}
        }

        let stress = match stress_at {
            Some(s) => {
                let unit = stress_unit(lines[s])?;
                // 表头行（`STRESS| x y z`）与摘要行（`1/3 Trace`、`Determinant`）
                // 都解析不出三个浮点，于是自动被跳过 —— 不必知道块里有几行表头
                let mut t = Matrix3::zeros();
                let mut row = 0usize;
                for l in &lines[s + 1..(s + 12).min(lines.len())] {
                    if row == 3 {
                        break;
                    }
                    if !line_matches(l, tag::STRESS_ROW) {
                        continue;
                    }
                    let f: Vec<&str> = l.split_whitespace().collect();
                    if f.len() < 4 {
                        continue;
                    }
                    // 取末尾三个：前面是 `STRESS|` 加行标，列数变了也不影响
                    let vals: Option<Vec<f64>> = f[f.len() - 3..]
                        .iter()
                        .map(|x| x.parse::<f64>().ok())
                        .collect();
                    if let Some(v) = vals {
                        for (col, x) in v.iter().enumerate() {
                            t[(row, col)] = convert_pressure(*x, unit, PressureUnit::EVPerAng3);
                        }
                        row += 1;
                    }
                }
                if row == 3 { Some(t) } else { None }
            }
            None => None,
        };

        let force_factor = HARTREE_TO_EV / BOHR_TO_ANG;
        let atoms: Vec<Atom> = c_syms
            .iter()
            .zip(&coords)
            .map(|(s, p)| Atom::new(s.as_str(), *p))
            .collect();
        let mut frame = Frame::with_cell(
            Cell::from_matrix(Matrix3::from_row_slice(&m)),
            [true; 3],
        );
        frame.atoms = atoms;
        frame.energy = energy;
        frame.forces = Some(forces.iter().map(|f| f * force_factor).collect());
        frame.stress = stress;
        frames.push(frame);
        stats.n_kept += 1;
    }

    if frames.is_empty() {
        bail!(
            "no usable frame: {} step(s), {} dropped ({} SCF, {} incomplete, {} composition)",
            stats.n_steps, stats.n_dropped(), stats.n_scf_failed,
            stats.n_incomplete, stats.n_bad_composition
        );
    }

    let mut traj = Trajectory { frames, metadata: Default::default() };
    traj.metadata.source = Some(match version {
        Some(v) => format!("CP2K {v} out"),
        None => "CP2K out".to_string(),
    });
    Ok((traj, stats))
}

#[cfg(test)]
mod tests {
    use super::*;

    // 两帧的迷你 out：第二帧的 SCF 故意不收敛
    const MINI: &str = concat!(
        " SCF run converged in     5 steps\n",
        " ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]          -10.500000000000000\n",
        " STRESS| Analytical stress tensor [bar]\n",
        " STRESS|                        x                   y                   z\n",
        " STRESS|      x        1.00000000000E+04   0.00000000000E+00   0.00000000000E+00\n",
        " STRESS|      y        0.00000000000E+00   2.00000000000E+04   0.00000000000E+00\n",
        " STRESS|      z        0.00000000000E+00   0.00000000000E+00   3.00000000000E+04\n",
        " MD| Step number                                                               1\n",
        " MD| Time [fs]                                                          1.000000\n",
        "      2\n",
        " i =        1, time =        1.000, E =       -10.5000000000\n",
        " Si         0.0000000000        0.0000000000        0.0000000000\n",
        "  O         1.0000000000        0.0000000000        0.0000000000\n",
        "      2\n",
        " i =        1, time =        1.000, E =       -10.5000000000\n",
        " Si         0.0100000000        0.0000000000        0.0000000000\n",
        "  O        -0.0100000000        0.0000000000        0.0000000000\n",
        "       1       1.000       5.0000000000        0.0000000000        0.0000000000",
        "        0.0000000000       5.0000000000        0.0000000000",
        "        0.0000000000        0.0000000000       5.0000000000          125.0000000000\n",
        " SCF run NOT converged\n",
        " ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]          -10.400000000000000\n",
        " STRESS| Analytical stress tensor [bar]\n",
        " STRESS|                        x                   y                   z\n",
        " STRESS|      x        1.00000000000E+04   0.00000000000E+00   0.00000000000E+00\n",
        " STRESS|      y        0.00000000000E+00   2.00000000000E+04   0.00000000000E+00\n",
        " STRESS|      z        0.00000000000E+00   0.00000000000E+00   3.00000000000E+04\n",
        " MD| Step number                                                               2\n",
        "      2\n",
        " i =        2, time =        2.000, E =       -10.4000000000\n",
        " Si         0.0000000000        0.0000000000        0.0000000000\n",
        "  O         1.0000000000        0.0000000000        0.0000000000\n",
        "      2\n",
        " i =        2, time =        2.000, E =       -10.4000000000\n",
        " Si         0.0100000000        0.0000000000        0.0000000000\n",
        "  O        -0.0100000000        0.0000000000        0.0000000000\n",
        "       2       2.000       5.0000000000        0.0000000000        0.0000000000",
        "        0.0000000000       5.0000000000        0.0000000000",
        "        0.0000000000        0.0000000000       5.0000000000          125.0000000000\n",
    );


    // 同一份数据的几种「换了排版」的写法；解析结果必须逐位相同。
    // 这是 token 匹配相对字面匹配的全部理由，所以它得有测试兜着。
    fn variants() -> Vec<(&'static str, String)> {
        let relayout = |f: fn(&str) -> String| -> String {
            MINI.lines().map(f).collect::<Vec<_>>().join("\n") + "\n"
        };
        vec![
            // 所有空白变三倍，模拟列宽调整
            ("wider columns", MINI.replace(' ', "   ")),
            // 行首缩进消失
            ("no indent", relayout(|l| l.trim_start().to_string())),
            // 缩进换成制表符
            ("tab indent", relayout(|l| format!("\t{}", l.trim_start()))),
            // cell 行前面多出一列（某个版本加了个计数器）
            (
                "extra cell column",
                MINI.replace("\n       1       1.000       5.0000000000",
                             "\n     42       1       1.000       5.0000000000"),
            ),
        ]
    }

    #[test]
    fn relayouts_parse_identically() {
        let (base, base_st) = parse_cp2k_out(MINI).unwrap();
        let b = &base.frames[0];
        for (name, text) in variants() {
            let (t, st) = parse_cp2k_out(&text)
                .unwrap_or_else(|e| panic!("{name}: {e:#}"));
            assert_eq!(st, base_st, "{name}: stats differ");
            assert_eq!(t.n_frames(), 1, "{name}");
            let f = &t.frames[0];
            assert_eq!(f.energy, b.energy, "{name}: energy");
            assert_eq!(f.forces, b.forces, "{name}: forces");
            assert_eq!(f.stress, b.stress, "{name}: stress");
            assert_eq!(
                f.cell.as_ref().unwrap().matrix,
                b.cell.as_ref().unwrap().matrix,
                "{name}: cell"
            );
            assert_eq!(
                f.atoms.iter().map(|a| a.position).collect::<Vec<_>>(),
                b.atoms.iter().map(|a| a.position).collect::<Vec<_>>(),
                "{name}: positions"
            );
        }
    }

    #[test]
    fn stress_block_tolerates_extra_header_rows() {
        // 在表头与数值之间插一行；靠"解析得出三个浮点"筛选，多几行表头无所谓
        let text = MINI.replace(
            " STRESS|                        x                   y                   z\n",
            " STRESS|                        x                   y                   z\n STRESS|  (in the cell frame)\n",
        );
        let (t, _) = parse_cp2k_out(&text).unwrap();
        let (base, _) = parse_cp2k_out(MINI).unwrap();
        assert_eq!(t.frames[0].stress, base.frames[0].stress);
    }

    #[test]
    fn version_string_reaches_the_metadata() {
        let text = format!(" CP2K| version string:                 CP2K version 2024.1\n{MINI}");
        let (t, _) = parse_cp2k_out(&text).unwrap();
        assert_eq!(t.metadata.source.as_deref(), Some("CP2K 2024.1 out"));
    }

    #[test]
    fn parses_units_and_drops_unconverged_frames() {
        let (traj, st) = parse_cp2k_out(MINI).unwrap();
        assert_eq!(st.n_steps, 2);
        assert_eq!(st.n_kept, 1);
        assert_eq!(st.n_scf_failed, 1);
        assert_eq!(st.n_incomplete, 0);
        assert_eq!(traj.n_frames(), 1);

        let f = &traj.frames[0];
        assert_eq!(f.atoms.len(), 2);
        assert_eq!(f.atoms[0].element, "Si");
        // -10.5 hartree
        assert!((f.energy.unwrap() - -10.5 * HARTREE_TO_EV).abs() < 1e-9);
        // 0.01 Hartree/Bohr
        let fx = f.forces.as_ref().unwrap()[0].x;
        assert!((fx - 0.01 * HARTREE_TO_EV / BOHR_TO_ANG).abs() < 1e-9);
        // 1e4 bar = 1 GPa；符号原样保留（正 = 压缩）
        let s = f.stress.unwrap();
        assert!((s[(0, 0)] - 1.0 / 160.217_663_4).abs() < 1e-9);
        assert!(s[(1, 1)] > s[(0, 0)]);
        assert_eq!(s[(0, 1)], 0.0);
        // cell 行优先：对角 5 Å
        let m = f.cell.as_ref().unwrap().matrix;
        assert_eq!((m[(0, 0)], m[(1, 1)], m[(2, 2)]), (5.0, 5.0, 5.0));
    }

    #[test]
    fn unknown_stress_unit_is_an_error() {
        let bad = MINI.replace("[bar]", "[atm]");
        let err = parse_cp2k_out(&bad).unwrap_err().to_string();
        assert!(err.contains("atm"), "{err}");
    }

    #[test]
    fn truncated_last_frame_is_dropped_not_fatal() {
        // 砍掉最后一帧的力块之后的一切，模拟作业被 kill
        let cut = MINI.find("       2       2.000").unwrap();
        let (traj, st) = parse_cp2k_out(&MINI[..cut]).unwrap();
        assert_eq!(traj.n_frames(), 1);
        assert_eq!(st.n_scf_failed, 1);
    }

    /// Runs against a real 40 MB restarted run; `examples/` is gitignored, so
    /// this is opt-in: `cargo test -p ferro-io -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn reads_the_reference_output() {
        let (traj, st) = read_cp2k_out_with_stats("../examples/total.out").unwrap();
        println!("{st:?}");
        println!("frames={} atoms={:?}", traj.n_frames(), traj.n_atoms());
        let f = &traj.frames[0];
        println!("energy={:?}", f.energy);
        println!("stress[0][0]={:?}", f.stress.map(|s| s[(0, 0)]));
        println!("force[0]={:?}", f.forces.as_ref().map(|v| v[0]));
        println!("cell={:?}", f.cell.as_ref().map(|c| c.matrix));
    }
}
