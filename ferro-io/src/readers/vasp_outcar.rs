//! VASP `OUTCAR` reader for AIMD trajectories.
//!
//! # What is taken, and why that one
//!
//! * **Energy** — `free  energy   TOTEN`, not `energy(sigma->0)`. The forces
//!   printed alongside are the derivatives of *that* free energy, so pairing
//!   the extrapolated energy with them would hand a model two halves of
//!   different functionals. dpdata makes the same choice, which also keeps
//!   datasets converted by either tool comparable.
//! * **Stress** — the `in kB` line, six Voigt components in VASP's own order
//!   `XX YY ZZ XY YZ ZX` (note: *not* the `XX YY ZZ YZ XZ XY` of the extxyz
//!   spec), scaled to eV/Å³ and stored **without a sign change**:
//!   [`ferro_core::Frame::stress`] is positive-for-compression, which is what
//!   VASP prints. ASE negates the same line on its way in because its own
//!   convention is the opposite one.
//! * **Cell** — read per frame from the block under `VOLUME and BASIS-vectors`.
//!   A frame without one is dropped rather than inheriting its predecessor's:
//!   in a fixed-cell run the two are identical and the bug would be invisible,
//!   which is exactly how a variable-cell run would come out silently wrong.
//! * **Elements** — `VRHFIN` symbols zipped with `ions per type`. Both appear
//!   more than once in a real OUTCAR; the repeats must agree.
//!
//! Blocks are located by token sequence and bounded by the next frame anchor,
//! never by a fixed line offset — dpdata reads the stress 14 lines below its
//! anchor, and that offset is a property of one VASP version's layout.

use std::io::{BufRead, BufReader};

use anyhow::{bail, Context, Result};
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};

use super::aimd::{AimdFormat, AimdStats};
use ferro_core::units::{convert_pressure, PressureUnit};

/// Reads a VASP OUTCAR, discarding the statistics.
pub fn read_vasp_outcar(path: &str) -> Result<Trajectory> {
    Ok(read_vasp_outcar_with_stats(path)?.0)
}

/// True when `needle` appears as a contiguous run of whitespace-separated
/// tokens in `line`.
///
/// Matching on tokens rather than on the raw line makes every anchor immune to
/// column widths, alignment and indentation, which differ between VASP builds
/// and between ensembles.
fn has_tokens(line: &str, needle: &[&str]) -> bool {
    let toks: Vec<&str> = line.split_whitespace().collect();
    if needle.len() > toks.len() {
        return false;
    }
    toks.windows(needle.len()).any(|w| w == needle)
}

/// The floats of a line, in order.
fn floats(line: &str) -> Vec<f64> {
    line.split_whitespace().filter_map(|t| t.parse::<f64>().ok()).collect()
}

/// `Iteration  <ionic>( <scf>)` → the two integers.
fn iteration_numbers(line: &str) -> Option<(i64, i64)> {
    let rest = line.split("Iteration").nth(1)?;
    let mut nums = Vec::new();
    let mut cur = String::new();
    for c in rest.chars() {
        if c.is_ascii_digit() {
            cur.push(c);
        } else if !cur.is_empty() {
            nums.push(cur.parse::<i64>().ok()?);
            cur.clear();
            if nums.len() == 2 {
                break;
            }
        }
    }
    if !cur.is_empty() && nums.len() < 2 {
        nums.push(cur.parse::<i64>().ok()?);
    }
    (nums.len() == 2).then(|| (nums[0], nums[1]))
}

/// Everything a frame needs, filled as its lines go by.
#[derive(Default)]
struct Pending {
    step: i64,
    converged: Option<bool>,
    cell: Option<Matrix3<f64>>,
    stress_kb: Option<[f64; 6]>,
    positions: Vec<Vector3<f64>>,
    forces: Vec<Vector3<f64>>,
    energy: Option<f64>,
}

/// Reads a VASP OUTCAR and reports what was dropped.
pub fn read_vasp_outcar_with_stats(path: &str) -> Result<(Trajectory, AimdStats)> {
    let file = std::fs::File::open(path).with_context(|| format!("cannot open {path}"))?;
    let reader = BufReader::new(file);
    let mut stats = AimdStats::new(AimdFormat::VaspOutcar);

    let mut symbols: Vec<String> = Vec::new();
    let mut symbols_seen: Vec<Vec<String>> = Vec::new();
    let mut counts: Option<Vec<usize>> = None;
    let mut elements: Vec<String> = Vec::new();

    let mut traj = Trajectory::new();
    let mut pending: Option<Pending> = None;
    // 需要跨行读取的块:剩余行数与去向
    let mut want_cell = 0usize;
    let mut cell_rows: Vec<[f64; 3]> = Vec::new();
    let mut in_cell_block = false;
    let mut want_atoms = 0usize;
    let mut skip_before_atoms = 0usize;

    for line in reader.lines() {
        let line = line.with_context(|| format!("reading {path}"))?;

        // ── 头部:元素与计数 ──────────────────────────────────────────────
        if line.contains("VRHFIN") {
            if let Some(sym) = line.split("VRHFIN").nth(1).and_then(|r| {
                r.trim_start().trim_start_matches('=').split(':').next()
            }) {
                let sym = sym.trim().to_string();
                if !sym.is_empty() {
                    symbols.push(sym);
                }
            }
            continue;
        }
        if has_tokens(&line, &["ions", "per", "type"]) {
            let here: Vec<usize> = line
                .split('=')
                .nth(1)
                .map(|r| r.split_whitespace().filter_map(|t| t.parse().ok()).collect())
                .unwrap_or_default();
            match &counts {
                None => counts = Some(here),
                Some(prev) if *prev != here => bail!(
                    "{path}: two `ions per type` lines disagree ({prev:?} vs {here:?})"
                ),
                Some(_) => {}
            }
            continue;
        }

        // ── 帧锚点 ───────────────────────────────────────────────────────
        if line.contains("Iteration") {
            if let Some((ionic, scf)) = iteration_numbers(&line) {
                if scf == 1 {
                    // 第一帧开始前把头部定下来
                    if elements.is_empty() {
                        elements = resolve_elements(&symbols, &counts, &mut symbols_seen, path)?;
                    }
                    finish(&mut pending, &mut traj, &mut stats, &elements);
                    // 离子步倒退 = 又一段运行接在后面。VASP 重启通常另写一个
                    // OUTCAR,但被 cat 到一起的文件在实际数据里很常见,而拼接处
                    // 那一步往往是被中断的半帧
                    if let Some((_, last)) = stats.steps {
                        if ionic < last {
                            stats.n_restarts += 1;
                        }
                    }
                    stats.n_steps += 1;
                    stats.steps = Some(match stats.steps {
                        None => (ionic, ionic),
                        Some((lo, _)) => (lo, ionic),
                    });
                    pending = Some(Pending { step: ionic, ..Default::default() });
                }
            }
            continue;
        }
        let Some(p) = pending.as_mut() else { continue };

        // ── 跨行块的续读 ─────────────────────────────────────────────────
        if want_cell > 0 {
            let v = floats(&line);
            if v.len() >= 3 {
                cell_rows.push([v[0], v[1], v[2]]);
                want_cell -= 1;
                if want_cell == 0 && cell_rows.len() == 3 {
                    p.cell = Some(Matrix3::new(
                        cell_rows[0][0], cell_rows[0][1], cell_rows[0][2],
                        cell_rows[1][0], cell_rows[1][1], cell_rows[1][2],
                        cell_rows[2][0], cell_rows[2][1], cell_rows[2][2],
                    ));
                }
            }
            continue;
        }
        if skip_before_atoms > 0 {
            skip_before_atoms -= 1;
            continue;
        }
        if want_atoms > 0 {
            let v = floats(&line);
            if v.len() < 6 {
                // 块被截断:交给完整性检查处理
                want_atoms = 0;
                continue;
            }
            p.positions.push(Vector3::new(v[0], v[1], v[2]));
            p.forces.push(Vector3::new(v[3], v[4], v[5]));
            want_atoms -= 1;
            continue;
        }

        // ── 帧内锚点 ─────────────────────────────────────────────────────
        if has_tokens(&line, &["aborting", "loop"]) {
            p.converged = Some(line.contains("EDIFF is reached"));
            continue;
        }
        if line.trim_start().starts_with("in kB") {
            let v = floats(&line);
            if v.len() >= 6 {
                p.stress_kb = Some([v[0], v[1], v[2], v[3], v[4], v[5]]);
            }
            continue;
        }
        if has_tokens(&line, &["VOLUME", "and", "BASIS-vectors"]) {
            in_cell_block = true;
            continue;
        }
        if in_cell_block && has_tokens(&line, &["direct", "lattice", "vectors"]) {
            in_cell_block = false;
            want_cell = 3;
            cell_rows.clear();
            continue;
        }
        if has_tokens(&line, &["POSITION", "TOTAL-FORCE"]) {
            want_atoms = elements.len();
            skip_before_atoms = 1; // 分隔线
            p.positions.clear();
            p.forces.clear();
            continue;
        }
        if has_tokens(&line, &["free", "energy", "TOTEN"]) {
            let v = floats(&line);
            p.energy = v.last().copied();
            continue;
        }
    }
    finish(&mut pending, &mut traj, &mut stats, &elements);

    if elements.is_empty() {
        bail!("{path}: no ionic step found (no `Iteration N( 1)` anchor)");
    }
    Ok((traj, stats))
}

/// Zips `VRHFIN` symbols with `ions per type`, checking the repeats agree.
fn resolve_elements(
    symbols: &[String],
    counts: &Option<Vec<usize>>,
    seen: &mut Vec<Vec<String>>,
    path: &str,
) -> Result<Vec<String>> {
    let counts = counts
        .as_ref()
        .with_context(|| format!("{path}: no `ions per type` line before the first step"))?;
    if counts.is_empty() {
        bail!("{path}: `ions per type` names no species");
    }
    if symbols.len() < counts.len() {
        bail!(
            "{path}: {} VRHFIN symbol(s) for {} species in `ions per type`. \
             A VASP <= 6.3 bug drops species names from OUTCAR; vasprun.xml \
             carries them properly",
            symbols.len(), counts.len()
        );
    }
    // 头部会重复打印几遍;每一遍必须逐字相同,否则说明这不是同一个体系
    for chunk in symbols.chunks(counts.len()) {
        if chunk.len() == counts.len() {
            seen.push(chunk.to_vec());
        }
    }
    if let Some(first) = seen.first() {
        for other in seen.iter().skip(1) {
            if other != first {
                bail!(
                    "{path}: the species list changes within the file ({first:?} \
                     then {other:?})"
                );
            }
        }
    }
    let names = &seen[0];
    let mut out = Vec::new();
    for (sym, n) in names.iter().zip(counts) {
        for _ in 0..*n {
            out.push(sym.clone());
        }
    }
    Ok(out)
}

/// Turns a completed [`Pending`] into a frame, or counts why it cannot.
fn finish(
    pending: &mut Option<Pending>,
    traj: &mut Trajectory,
    stats: &mut AimdStats,
    elements: &[String],
) {
    let Some(p) = pending.take() else { return };
    if p.converged == Some(false) {
        stats.n_scf_failed += 1;
        return;
    }
    let n = elements.len();
    let (Some(cell), Some(energy)) = (p.cell, p.energy) else {
        stats.n_incomplete += 1;
        return;
    };
    if p.positions.len() != n || p.forces.len() != n || p.converged.is_none() {
        stats.n_incomplete += 1;
        return;
    }
    let mut frame = Frame::with_cell(Cell::from_matrix(cell), [true; 3]);
    for (i, sym) in elements.iter().enumerate() {
        frame.add_atom(Atom::new(sym.clone(), p.positions[i]));
    }
    frame.forces = Some(p.forces.clone());
    frame.energy = Some(energy);
    // in kB 的六个数是 VASP 自己的 Voigt 顺序 XX YY ZZ XY YZ ZX,与 extxyz 规格的
    // XX YY ZZ YZ XZ XY 不同。符号不变:VASP 与 Frame::stress 都是正 = 压缩
    if let Some(kb) = p.stress_kb {
        let e = |v: f64| convert_pressure(v, PressureUnit::Kbar, PressureUnit::EVPerAng3);
        let (xx, yy, zz) = (e(kb[0]), e(kb[1]), e(kb[2]));
        let (xy, yz, zx) = (e(kb[3]), e(kb[4]), e(kb[5]));
        frame.stress = Some(Matrix3::new(xx, xy, zx, xy, yy, yz, zx, yz, zz));
    }
    let _ = p.step;
    traj.add_frame(frame);
    stats.n_kept += 1;
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIXTURE: &str = "../tests/vasp_OUTCAR_2frames";

    /// Two frames with DIFFERENT cells, so a reader that reads the cell once
    /// and reuses it fails here. Real fixed-cell runs cannot catch that bug —
    /// the wrong answer is bit-identical to the right one.
    const VARIABLE_CELL: &str = r#" vasp.6.4.2 20Jul23 (build Aug 01 2023 12:00:00) complex
   VRHFIN =H: 1s1
   VRHFIN =O: 2s2 2p4
   ions per type =               1   1
--------------------------------------- Iteration      1(   1)  ---------------------------------------
------------------------ aborting loop because EDIFF is reached ----------------------------------------
  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Total       1.00000     2.00000     3.00000     0.10000     0.20000     0.30000
  in kB      10.00000    20.00000    30.00000     1.00000     2.00000     3.00000
 VOLUME and BASIS-vectors are now :
 -----------------------------------------------------------------------------
      direct lattice vectors                 reciprocal lattice vectors
    10.000000000  0.000000000  0.000000000     0.100000000  0.000000000  0.000000000
     0.000000000 10.000000000  0.000000000     0.000000000  0.100000000  0.000000000
     0.000000000  0.000000000 10.000000000     0.000000000  0.000000000  0.100000000

 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      0.00000      0.00000      0.00000         0.100000      0.200000      0.300000
      1.00000      1.00000      1.00000        -0.100000     -0.200000     -0.300000
 -----------------------------------------------------------------------------------
    total drift:                                0.000000      0.000000      0.000000

  free  energy   TOTEN  =       -10.00000000 eV

--------------------------------------- Iteration      2(   1)  ---------------------------------------
------------------------ aborting loop because EDIFF is reached ----------------------------------------
  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Total       1.00000     2.00000     3.00000     0.10000     0.20000     0.30000
  in kB      10.00000    20.00000    30.00000     1.00000     2.00000     3.00000
 VOLUME and BASIS-vectors are now :
 -----------------------------------------------------------------------------
      direct lattice vectors                 reciprocal lattice vectors
    11.000000000  0.000000000  0.000000000     0.090909091  0.000000000  0.000000000
     0.000000000 12.000000000  0.000000000     0.000000000  0.083333333  0.000000000
     0.000000000  0.000000000 13.000000000     0.000000000  0.000000000  0.076923077

 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
      0.00000      0.00000      0.00000         0.100000      0.200000      0.300000
      1.00000      1.00000      1.00000        -0.100000     -0.200000     -0.300000
 -----------------------------------------------------------------------------------
  free  energy   TOTEN  =       -20.00000000 eV
"#;

    fn tmp(name: &str, body: &str) -> String {
        let p = std::env::temp_dir().join(name);
        std::fs::write(&p, body).unwrap();
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn reads_two_real_frames() {
        let (traj, stats) = read_vasp_outcar_with_stats(FIXTURE).unwrap();
        assert_eq!(traj.n_frames(), 2);
        assert_eq!(stats.n_kept, 2);
        assert_eq!(stats.n_dropped(), 0);
        assert_eq!(stats.steps, Some((1, 2)));
        let f = traj.first().unwrap();
        assert_eq!(f.n_atoms(), 297);
        // VRHFIN 的顺序 zip ions per type: O 198, P 66, Zn 33
        assert_eq!(f.atom(0).element, "O");
        assert_eq!(f.atom(197).element, "O");
        assert_eq!(f.atom(198).element, "P");
        assert_eq!(f.atom(263).element, "P");
        assert_eq!(f.atom(264).element, "Zn");
        assert_eq!(f.atom(296).element, "Zn");
    }

    /// Reference values from dpdata 1.0.2 reading the same file — an outside
    /// implementation, not this reader's own round trip.
    #[test]
    fn matches_dpdata_on_the_real_fixture() {
        let (traj, _) = read_vasp_outcar_with_stats(FIXTURE).unwrap();
        let f0 = &traj.frames[0];
        assert!((f0.energy.unwrap() - -2041.07506914).abs() < 1e-8);
        assert!((traj.frames[1].energy.unwrap() - -2045.47272826).abs() < 1e-8);

        let p0 = f0.atom(0).position;
        for (got, want) in p0.iter().zip(&[2.15183, 5.81516, 1.20982]) {
            assert!((got - want).abs() < 1e-9, "{p0:?}");
        }
        let fl = f0.forces.as_ref().unwrap().last().unwrap();
        for (got, want) in fl.iter().zip(&[-1.609282, 0.303173, -0.607801]) {
            assert!((got - want).abs() < 1e-9, "{fl:?}");
        }
        // sigma = virial/V from dpdata; 容差放到 1e-7 相对量,因为 dpdata 用的
        // eV/Å³→GPa 常数是 160.2176621 而 units.rs 用 CODATA 2018 的 160.2176634
        let s = f0.stress.unwrap();
        let want = [
            [0.10444031438653728, -0.000673727219491115, 0.005413310796275813],
            [-0.000673727219491115, 0.1040530100207972, -0.003939890220130731],
            [0.005413310796275813, -0.003939890220130731, 0.10906719503304999],
        ];
        for i in 0..3 {
            for j in 0..3 {
                let (g, w) = (s[(i, j)], want[i][j]);
                assert!((g - w).abs() <= 1e-7 * w.abs().max(1e-3), "({i},{j}) {g} vs {w}");
            }
        }
    }

    #[test]
    fn the_stress_is_symmetric_and_positive_for_compression() {
        let (traj, _) = read_vasp_outcar_with_stats(FIXTURE).unwrap();
        let s = traj.frames[0].stress.unwrap();
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            assert!((s[(i, j)] - s[(j, i)]).abs() < 1e-15);
        }
        // 3000 K 的致密玻璃处于受压状态,VASP 打的 in kB 为正,Ferro 同号
        assert!(s[(0, 0)] > 0.0);
    }

    #[test]
    fn every_frame_reads_its_own_cell() {
        let (traj, _) =
            read_vasp_outcar_with_stats(&tmp("vc.outcar", VARIABLE_CELL)).unwrap();
        assert_eq!(traj.n_frames(), 2);
        let a = traj.frames[0].cell.as_ref().unwrap().matrix;
        let b = traj.frames[1].cell.as_ref().unwrap().matrix;
        assert!((a[(0, 0)] - 10.0).abs() < 1e-12);
        assert!((b[(0, 0)] - 11.0).abs() < 1e-12, "frame 1 reused frame 0's cell");
        assert!((b[(1, 1)] - 12.0).abs() < 1e-12);
        assert!((b[(2, 2)] - 13.0).abs() < 1e-12);
    }

    #[test]
    fn the_voigt_order_is_vasps_own() {
        // in kB = 10 20 30 1 2 3 → XX YY ZZ XY YZ ZX（不是 extxyz 规格的
        // XX YY ZZ YZ XZ XY）。搞混的话 xy 与 yz 会互换,而对角线看不出来
        let (traj, _) =
            read_vasp_outcar_with_stats(&tmp("voigt.outcar", VARIABLE_CELL)).unwrap();
        let s = traj.frames[0].stress.unwrap();
        let k = |kb: f64| convert_pressure(kb, PressureUnit::Kbar, PressureUnit::EVPerAng3);
        assert!((s[(0, 1)] - k(1.0)).abs() < 1e-15, "xy");
        assert!((s[(1, 2)] - k(2.0)).abs() < 1e-15, "yz");
        assert!((s[(0, 2)] - k(3.0)).abs() < 1e-15, "zx");
    }

    #[test]
    fn an_unconverged_frame_is_dropped_and_counted() {
        let body = VARIABLE_CELL.replace(
            "------------------------ aborting loop because EDIFF is reached ---",
            "------------------------ aborting loop EDIFF was not reached ---",
        );
        let (traj, stats) =
            read_vasp_outcar_with_stats(&tmp("unconv.outcar", &body)).unwrap();
        assert_eq!(traj.n_frames(), 0);
        assert_eq!(stats.n_scf_failed, 2);
        assert_eq!(stats.n_kept, 0);
    }

    #[test]
    fn a_truncated_last_frame_is_dropped_not_half_read() {
        let cut = VARIABLE_CELL.find("  free  energy   TOTEN  =       -20").unwrap();
        let (traj, stats) =
            read_vasp_outcar_with_stats(&tmp("trunc.outcar", &VARIABLE_CELL[..cut])).unwrap();
        assert_eq!(traj.n_frames(), 1);
        assert_eq!(stats.n_incomplete, 1);
    }

    #[test]
    fn fewer_species_than_ions_per_type_names_vasprun_as_the_way_out() {
        let body = VARIABLE_CELL.replace("   VRHFIN =O: 2s2 2p4\n", "");
        let e = read_vasp_outcar_with_stats(&tmp("nospec.outcar", &body)).unwrap_err();
        let msg = format!("{e:#}");
        assert!(msg.contains("vasprun.xml"), "{msg}");
    }

    #[test]
    fn a_restart_appended_to_the_same_file_is_counted() {
        // 离子步倒退 = 又一段运行。真实数据里就是这样:examples 那份 OUTCAR
        // 是 1..1575 接 1..425
        let body = format!("{VARIABLE_CELL}{}", VARIABLE_CELL
            .split("--------------------------------------- Iteration      1(   1)")
            .nth(1)
            .map(|rest| format!(
                "--------------------------------------- Iteration      1(   1){rest}"))
            .unwrap());
        let (traj, stats) =
            read_vasp_outcar_with_stats(&tmp("restart.outcar", &body)).unwrap();
        assert_eq!(traj.n_frames(), 4);
        assert_eq!(stats.n_restarts, 1);
    }
}
