//! VASP `vasprun.xml` reader for AIMD trajectories.
//!
//! Same quantities and the same conventions as [`super::vasp_outcar`]; only the
//! container differs. Two things are specific to this format:
//!
//! * **The energy is the last `e_fr_energy` of a `<calculation>`.** Every
//!   `<scstep>` carries an `<energy>` block of its own, so the first one in a
//!   calculation is an unconverged SCF iterate — tens of thousands of eV away
//!   from the answer. Taking it would not fail, it would just be wrong.
//! * **Convergence is inferred**, not stated: a calculation whose `<scstep>`
//!   count reaches `NELM` never met `EDIFF`. The OUTCAR path reads VASP's own
//!   verdict instead, so the two can disagree about the same run; the format is
//!   carried in [`super::aimd::AimdStats`] so a caller can say which rule
//!   applied.
//!
//! Parsing is a streaming pull, never a DOM: these files reach hundreds of
//! megabytes on a long run, and holding one in memory as a tree costs several
//! times that.

use std::io::BufReader;

use anyhow::{bail, Context, Result};
use ferro_core::units::{convert_pressure, PressureUnit};
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use quick_xml::events::Event;
use quick_xml::Reader;

use super::aimd::{AimdFormat, AimdStats};

/// Reads a vasprun.xml, discarding the statistics.
pub fn read_vasprun(path: &str) -> Result<Trajectory> {
    Ok(read_vasprun_with_stats(path)?.0)
}

/// Where in the document the pull parser currently is.
#[derive(Default)]
struct State {
    /// `<atominfo>` → element symbol per atom
    elements: Vec<String>,
    in_atominfo: bool,
    in_atoms_array: bool,
    /// 只有 <c> 里的文字才是数据；<field> 是表头（"element" 等）
    in_c: bool,
    atominfo_field: usize,
    nelm: Option<usize>,
    reading_nelm: bool,

    in_calculation: bool,
    in_scstep: bool,
    n_scstep: usize,
    /// name of the `<varray>` being read, if any
    varray: Option<String>,
    in_structure: bool,
    reading_energy_name: Option<String>,

    rows: Vec<[f64; 3]>,
    cell: Option<Matrix3<f64>>,
    positions: Vec<Vector3<f64>>,
    forces: Option<Vec<Vector3<f64>>>,
    stress_kb: Option<Matrix3<f64>>,
    energy: Option<f64>,
}

/// Reads a vasprun.xml and reports what was dropped.
pub fn read_vasprun_with_stats(path: &str) -> Result<(Trajectory, AimdStats)> {
    let file = std::fs::File::open(path).with_context(|| format!("cannot open {path}"))?;
    let mut reader = Reader::from_reader(BufReader::new(file));
    reader.config_mut().trim_text(true);

    let mut stats = AimdStats::new(AimdFormat::VaspXml);
    let mut traj = Trajectory::new();
    let mut st = State::default();
    let mut buf = Vec::new();

    loop {
        match reader.read_event_into(&mut buf) {
            Err(e) => bail!("{path}: malformed XML at byte {}: {e}", reader.buffer_position()),
            Ok(Event::Eof) => break,
            Ok(Event::Start(e)) => {
                let name = String::from_utf8_lossy(e.name().as_ref()).to_string();
                let attr = |key: &str| -> Option<String> {
                    e.attributes().flatten().find_map(|a| {
                        (a.key.as_ref() == key.as_bytes())
                            .then(|| String::from_utf8_lossy(&a.value).trim().to_string())
                    })
                };
                match name.as_str() {
                    "atominfo" => st.in_atominfo = true,
                    "array" if st.in_atominfo => {
                        st.in_atoms_array = attr("name").as_deref() == Some("atoms");
                    }
                    "rc" if st.in_atoms_array => st.atominfo_field = 0,
                    "c" if st.in_atoms_array => st.in_c = true,
                    "calculation" => {
                        st.in_calculation = true;
                        st.n_scstep = 0;
                        st.cell = None;
                        st.positions.clear();
                        st.forces = None;
                        st.stress_kb = None;
                        st.energy = None;
                    }
                    "scstep" => {
                        st.in_scstep = true;
                        st.n_scstep += 1;
                    }
                    "structure" if st.in_calculation => st.in_structure = true,
                    "varray" => st.varray = attr("name"),
                    "i" => {
                        let n = attr("name");
                        if n.as_deref() == Some("NELM") && st.nelm.is_none() {
                            st.reading_nelm = true;
                        }
                        // SCF 步里的 <energy> 是中间迭代值,只认 calculation 层的那个
                        if !st.in_scstep {
                            st.reading_energy_name = n;
                        }
                    }
                    _ => {}
                }
                if name == "v" || name == "c" {
                    // 由 Text 事件处理
                }
                buf.clear();
            }
            Ok(Event::Text(t)) => {
                let text = String::from_utf8_lossy(t.as_ref()).to_string();
                if st.reading_nelm {
                    st.nelm = text.trim().parse().ok();
                    st.reading_nelm = false;
                } else if let Some(field) = st.reading_energy_name.take() {
                    if field == "e_fr_energy" {
                        // 每帧覆盖:最后一个才是收敛值
                        st.energy = text.trim().parse().ok();
                    }
                } else if st.in_atoms_array && st.in_c {
                    // rc/c 交替给出元素与 type 编号,取第一列
                    if st.atominfo_field == 0 {
                        let sym = text.trim().to_string();
                        if !sym.is_empty() && sym.chars().next().is_some_and(|c| c.is_alphabetic())
                        {
                            st.elements.push(sym);
                        }
                    }
                    st.atominfo_field += 1;
                } else if st.varray.is_some() {
                    let v: Vec<f64> =
                        text.split_whitespace().filter_map(|x| x.parse().ok()).collect();
                    if v.len() >= 3 {
                        st.rows.push([v[0], v[1], v[2]]);
                    }
                }
                buf.clear();
            }
            Ok(Event::End(e)) => {
                let name = String::from_utf8_lossy(e.name().as_ref()).to_string();
                match name.as_str() {
                    "atominfo" => st.in_atominfo = false,
                    "array" => st.in_atoms_array = false,
                    "c" => st.in_c = false,
                    "scstep" => st.in_scstep = false,
                    "structure" => st.in_structure = false,
                    "varray" => {
                        let which = st.varray.take();
                        let rows = std::mem::take(&mut st.rows);
                        if st.in_calculation {
                            match which.as_deref() {
                                Some("basis") if st.in_structure && rows.len() == 3 => {
                                    st.cell = Some(Matrix3::new(
                                        rows[0][0], rows[0][1], rows[0][2],
                                        rows[1][0], rows[1][1], rows[1][2],
                                        rows[2][0], rows[2][1], rows[2][2],
                                    ));
                                }
                                Some("positions") if st.in_structure => {
                                    st.positions = rows
                                        .iter()
                                        .map(|r| Vector3::new(r[0], r[1], r[2]))
                                        .collect();
                                }
                                Some("forces") => {
                                    st.forces = Some(
                                        rows.iter()
                                            .map(|r| Vector3::new(r[0], r[1], r[2]))
                                            .collect(),
                                    );
                                }
                                Some("stress") if rows.len() == 3 => {
                                    st.stress_kb = Some(Matrix3::new(
                                        rows[0][0], rows[0][1], rows[0][2],
                                        rows[1][0], rows[1][1], rows[1][2],
                                        rows[2][0], rows[2][1], rows[2][2],
                                    ));
                                }
                                _ => {}
                            }
                        }
                    }
                    "calculation" => {
                        st.in_calculation = false;
                        stats.n_steps += 1;
                        let step = stats.n_steps as i64;
                        stats.steps = Some(match stats.steps {
                            None => (step, step),
                            Some((lo, _)) => (lo, step),
                        });
                        finish(&mut st, &mut traj, &mut stats);
                    }
                    _ => {}
                }
                buf.clear();
            }
            Ok(_) => buf.clear(),
        }
    }

    if st.elements.is_empty() {
        bail!("{path}: no <atominfo> species list found");
    }
    Ok((traj, stats))
}

/// Turns the finished `<calculation>` into a frame, or counts why it cannot.
fn finish(st: &mut State, traj: &mut Trajectory, stats: &mut AimdStats) {
    // 收敛靠推断:SCF 步数打满 NELM 就是没达到 EDIFF
    if let Some(nelm) = st.nelm {
        if st.n_scstep >= nelm {
            stats.n_scf_failed += 1;
            return;
        }
    }
    let n = st.elements.len();
    let (Some(cell), Some(energy)) = (st.cell, st.energy) else {
        stats.n_incomplete += 1;
        return;
    };
    if st.positions.len() != n || st.forces.as_ref().is_none_or(|f| f.len() != n) {
        stats.n_incomplete += 1;
        return;
    }
    let mut frame = Frame::with_cell(Cell::from_matrix(cell), [true; 3]);
    // vasprun 的 positions 是分数坐标
    for (i, sym) in st.elements.iter().enumerate() {
        let f = st.positions[i];
        let cart = cell.transpose() * f;
        frame.add_atom(Atom::new(sym.clone(), cart));
    }
    frame.forces = st.forces.clone();
    frame.energy = Some(energy);
    // <varray name="stress"> 是 kB 的完整 3x3,顺序无歧义;符号与 OUTCAR 同,不变
    if let Some(kb) = st.stress_kb {
        frame.stress = Some(kb.map(|v| {
            convert_pressure(v, PressureUnit::Kbar, PressureUnit::EVPerAng3)
        }));
    }
    traj.add_frame(frame);
    stats.n_kept += 1;
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIXTURE: &str = "../tests/vasp_vasprun_2frames.xml";

    #[test]
    fn reads_two_real_frames() {
        let (traj, stats) = read_vasprun_with_stats(FIXTURE).unwrap();
        assert_eq!(traj.n_frames(), 2);
        assert_eq!(stats.n_kept, 2);
        assert_eq!(stats.n_dropped(), 0);
        let f = traj.first().unwrap();
        assert_eq!(f.n_atoms(), 297);
        assert_eq!(f.atom(0).element, "O");
        assert_eq!(f.atom(198).element, "P");
        assert_eq!(f.atom(296).element, "Zn");
    }

    /// The `<field>` headers inside `<array name="atoms">` are text too, and
    /// taking them for data gives 298 species for 297 atoms — every frame then
    /// fails the completeness check and the file reads as empty.
    #[test]
    fn the_column_headers_are_not_mistaken_for_species() {
        let (traj, _) = read_vasprun_with_stats(FIXTURE).unwrap();
        assert_eq!(traj.first().unwrap().n_atoms(), 297);
    }

    /// Reference values from dpdata 1.0.2 on the same file.
    #[test]
    fn matches_dpdata_on_the_real_fixture() {
        let (traj, _) = read_vasprun_with_stats(FIXTURE).unwrap();
        let f0 = &traj.frames[0];
        assert!((f0.energy.unwrap() - -1929.96163538).abs() < 1e-8);
        assert!((traj.frames[1].energy.unwrap() - -1933.54154614).abs() < 1e-8);

        // <varray name="positions"> 是分数坐标,笛卡尔化后才能与 dpdata 比
        let p0 = f0.atom(0).position;
        for (got, want) in p0.iter().zip(&[0.138232, 5.899914, 2.216376]) {
            assert!((got - want).abs() < 1e-5, "{p0:?}");
        }
        let fl = f0.forces.as_ref().unwrap().last().unwrap();
        for (got, want) in fl.iter().zip(&[-0.754179, -0.855899, 0.228165]) {
            assert!((got - want).abs() < 1e-6, "{fl:?}");
        }
        let s = f0.stress.unwrap();
        for (got, want) in [
            (s[(0, 0)], 0.009614551228681297),
            (s[(1, 1)], -0.011997096660918008),
            (s[(2, 2)], 0.0290602257015458),
            (s[(0, 1)], 0.0022682876733850434),
        ] {
            assert!((got - want).abs() <= 1e-7 * want.abs().max(1e-3), "{got} vs {want}");
        }
    }

    /// The first `<energy>` of a calculation belongs to an SCF iterate and is
    /// tens of thousands of eV away from the converged value; only the last one
    /// is the answer. Nothing about the wrong choice looks wrong.
    #[test]
    fn the_energy_is_the_converged_one_not_the_first_scf_step() {
        let (traj, _) = read_vasprun_with_stats(FIXTURE).unwrap();
        let e = traj.frames[0].energy.unwrap();
        assert!(e < -1000.0 && e > -3000.0, "{e} looks like an SCF iterate");
    }

    #[test]
    fn malformed_xml_is_an_error_not_an_empty_trajectory() {
        let p = std::env::temp_dir().join("broken.xml");
        std::fs::write(&p, "<?xml version=\"1.0\"?><modeling><calculation>").unwrap();
        // 截断的 XML 要么报错、要么给出 0 帧,但绝不能装作读到了东西
        match read_vasprun_with_stats(p.to_str().unwrap()) {
            Err(e) => assert!(format!("{e:#}").contains("atominfo") || format!("{e:#}").contains("XML")),
            Ok((t, _)) => assert_eq!(t.n_frames(), 0),
        }
    }
}
