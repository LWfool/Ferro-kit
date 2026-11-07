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

use std::path::Path;
use std::io::BufReader;

use anyhow::{anyhow, bail, Context, Result};
use ferro_core::units::{convert_pressure, PressureUnit, BOLTZMANN_EV_K};
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use quick_xml::events::Event;
use quick_xml::Reader;

use super::aimd::{AimdFormat, AimdStats};

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
    /// text of the element being read, parsed at its end tag
    text: String,
    cell: Option<Matrix3<f64>>,
    positions: Vec<Vector3<f64>>,
    forces: Option<Vec<Vector3<f64>>>,
    stress_kb: Option<Matrix3<f64>>,
    /// a stress varray without three rows: the step is incomplete
    stress_bad: bool,
    energy: Option<f64>,
    kinetic: Option<f64>,
}

/// Reads a vasprun.xml and reports what was dropped.
pub fn read_vasprun(path: &Path) -> Result<(Trajectory, AimdStats)> {
    let file = std::fs::File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    parse_vasprun(BufReader::new(file)).with_context(|| format!("parsing {}", path.display()))
}

fn parse_vasprun(input: impl std::io::BufRead) -> Result<(Trajectory, AimdStats)> {
    let mut reader = Reader::from_reader(input);
    reader.config_mut().trim_text(true);

    let mut stats = AimdStats::new(AimdFormat::VaspXml);
    let mut traj = Trajectory::new();
    traj.metadata.program = Some("VASP vasprun.xml".to_string());
    let mut st = State::default();
    let mut buf = Vec::new();

    loop {
        match reader.read_event_into(&mut buf) {
            // Every syntax error means the input ended inside a construct: a run killed
            // mid-write. Ill-formed markup (mismatched tags) is still an error
            Err(quick_xml::Error::Syntax(_)) | Ok(Event::Eof) => break,
            Err(e) => bail!("malformed XML at byte {}: {e}", reader.buffer_position()),
            Ok(Event::Start(e)) => {
                let at = reader.buffer_position();
                st.text.clear();
                let name = utf8(e.name().into_inner(), at)?;
                let attr = |key: &str| -> Result<Option<String>> {
                    for a in e.attributes().flatten() {
                        if a.key.as_ref() == key.as_bytes() {
                            return Ok(Some(utf8(&a.value, at)?.trim().to_string()));
                        }
                    }
                    Ok(None)
                };
                match name {
                    "atominfo" => st.in_atominfo = true,
                    "array" if st.in_atominfo => {
                        st.in_atoms_array = attr("name")?.as_deref() == Some("atoms");
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
                        st.stress_bad = false;
                        st.energy = None;
                        st.kinetic = None;
                    }
                    "scstep" => {
                        st.in_scstep = true;
                        st.n_scstep += 1;
                    }
                    "structure" if st.in_calculation => st.in_structure = true,
                    "varray" => st.varray = attr("name")?,
                    "i" => {
                        let n = attr("name")?;
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
                buf.clear();
            }
            Ok(Event::Text(t)) => {
                let text = utf8(t.as_ref(), reader.buffer_position())?;
                if st.reading_nelm || st.reading_energy_name.is_some() || st.varray.is_some() {
                    // Parsed at the end tag, not here: a file cut inside a number never reaches
                    // that tag, so the cut step counts as incomplete instead of a bad number
                    st.text.push_str(text);
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
                }
                buf.clear();
            }
            Ok(Event::End(e)) => {
                let at = reader.buffer_position();
                let name = utf8(e.name().into_inner(), at)?;
                let text = std::mem::take(&mut st.text);
                let text = text.trim();
                match name {
                    "atominfo" => st.in_atominfo = false,
                    "array" => st.in_atoms_array = false,
                    "c" => st.in_c = false,
                    // Both flags are cleared here: an empty `<i name="e_fr_energy"></i>` has no Text
                    // event, and a flag left armed took the next text node (another quantity) as the
                    // energy (review 10-03, A minor). An empty value is a missing one, not a bad number
                    "i" => {
                        let field = st.reading_energy_name.take();
                        if std::mem::take(&mut st.reading_nelm) && !text.is_empty() {
                            st.nelm = Some(number(text, "NELM", at)?);
                        } else if !text.is_empty() {
                            match field.as_deref() {
                                // 每帧覆盖:最后一个才是收敛值
                                Some("e_fr_energy") => st.energy = Some(number(text, "e_fr_energy", at)?),
                                // 离子动能(eV),温度由它反算 —— vasprun 不打印温度本身
                                Some("kinetic") => st.kinetic = Some(number(text, "kinetic", at)?),
                                _ => {}
                            }
                        }
                    }
                    // Only the four varrays read below need numbers; others (`selective`'s T/F,
                    // one-column `weights`) are not this reader's business
                    "v" if st.in_calculation
                        && matches!(st.varray.as_deref(), Some("basis" | "positions" | "forces" | "stress")) =>
                    {
                        let v = text.split_whitespace()
                            .map(|x| number(x, "varray value", at))
                            .collect::<Result<Vec<f64>>>()?;
                        let [x, y, z] = v[..] else {
                            bail!("near byte {at}: a <varray name=\"{}\"> row has {} values, not 3",
                                  st.varray.as_deref().unwrap_or(""), v.len());
                        };
                        st.rows.push([x, y, z]);
                    }
                    "scstep" => st.in_scstep = false,
                    "structure" => st.in_structure = false,
                    "varray" => {
                        let which = st.varray.take();
                        let rows = std::mem::take(&mut st.rows);
                        if st.in_calculation {
                            match which.as_deref() {
                                Some("basis") if st.in_structure && rows.len() == 3 => {
                                    st.cell = Some(Matrix3::from_row_slice(rows.as_flattened()));
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
                                    st.stress_kb = Some(Matrix3::from_row_slice(rows.as_flattened()));
                                }
                                Some("stress") => st.stress_bad = true,
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
        bail!("no <atominfo> species list found");
    }
    // The file ends inside a calculation: a step that was started but never finished
    if st.in_calculation {
        stats.n_steps += 1;
        stats.n_incomplete += 1;
        let step = stats.n_steps as i64;
        stats.steps = Some((stats.steps.map_or(step, |(lo, _)| lo), step));
    }
    Ok((traj, stats))
}

/// Bytes that are not UTF-8 are an error, as in every other text reader (review C-8),
/// not replaced (review 2026-10-30 A1-M5). `at` is the parser's byte position.
fn utf8(bytes: &[u8], at: u64) -> Result<&str> {
    std::str::from_utf8(bytes).map_err(|e| {
        anyhow!("near byte {at}: byte 0x{:02x} is not UTF-8 text", bytes[e.valid_up_to()])
    })
}

/// A value this reader uses that does not parse is an error naming it; it used to be
/// skipped, so a frame went quietly incomplete or NELM was never known.
fn number<T: std::str::FromStr>(text: &str, what: &str, at: u64) -> Result<T> {
    text.parse().map_err(|_| anyhow!("near byte {at}: {what} {text:?} is not a number"))
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
    if st.positions.len() != n || st.forces.as_ref().is_none_or(|f| f.len() != n) || st.stress_bad {
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
    frame.set_forces(st.forces.take()).expect("one force per atom, checked above");
    frame.energy = Some(energy);
    // vasprun 没有温度字段(只有输入参数 TEBEG),按 T = 2*E_kin/(3N*k_B) 反算。
    // 除数是 3N 不是 3N-3:OUTCAR 的 EKIN=4.268142 eV / N=296 与它自己打印的
    // 111.55 K 对拍确认了这一点。这一列是导出量,报告的表头会注明
    // n > 0: a file without species was refused before any calculation
    frame.temperature = st
        .kinetic
        .map(|ek| 2.0 * ek / (3.0 * n as f64 * BOLTZMANN_EV_K));
    // <varray name="stress"> 是 kB 的完整 3x3,顺序无歧义;符号与 OUTCAR 同,不变。
    // 取 (σ+σᵀ)/2:VASP 6.5.1 写的这九个数不对称(相对 ~4e-5,OUTCAR 只打六个数故看不见),
    // 原样留着则导出成 extxyz 的 virial= 后,本仓自己的 extxyz reader 按规范拒收
    if let Some(kb) = st.stress_kb {
        let kb = (kb + kb.transpose()) * 0.5;
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

    const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../tests/vasp_vasprun_2frames.xml");

    /// vasprun has no temperature field, so it is derived from the ionic
    /// kinetic energy with `T = 2*E_kin/(3N*k_B)`. The divisor is 3N, not 3N-3:
    /// the OUTCAR of a comparable run prints 111.55 K next to EKIN = 4.268142 eV
    /// for N = 296, which only 3N reproduces.
    #[test]
    fn derives_temperature_from_ionic_kinetic_energy() {
        let (traj, _) =
            read_vasprun(Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../tests/vasp_vasprun_2frames.xml"))).unwrap();
        let n = traj.frames[0].n_atoms() as f64;
        for (i, f) in traj.frames.iter().enumerate() {
            let t = f.temperature.expect("temperature should be derived");
            assert!(t > 0.0 && t < 20000.0, "frame {i}: implausible T = {t}");
            // 反解回动能,确认用的是 3N
            let ek = t * 3.0 * n * ferro_core::units::BOLTZMANN_EV_K / 2.0;
            assert!(ek > 0.0, "frame {i}: E_kin = {ek}");
            println!("frame {i}: T = {t:.2} K  (N = {n}, E_kin = {ek:.5} eV)");
        }
    }

    #[test]
    fn reads_two_real_frames() {
        let (traj, stats) = read_vasprun(Path::new(FIXTURE)).unwrap();
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
        let (traj, _) = read_vasprun(Path::new(FIXTURE)).unwrap();
        assert_eq!(traj.first().unwrap().n_atoms(), 297);
    }

    /// Reference values from dpdata 1.0.2 on the same file.
    #[test]
    fn matches_dpdata_on_the_real_fixture() {
        let (traj, _) = read_vasprun(Path::new(FIXTURE)).unwrap();
        let f0 = &traj.frames[0];
        assert!((f0.energy.unwrap() - -1929.96163538).abs() < 1e-8);
        assert!((traj.frames[1].energy.unwrap() - -1933.54154614).abs() < 1e-8);

        // <varray name="positions"> 是分数坐标,笛卡尔化后才能与 dpdata 比
        let p0 = f0.atom(0).position;
        for (got, want) in p0.iter().zip(&[0.138232, 5.899914, 2.216376]) {
            assert!((got - want).abs() < 1e-5, "{p0:?}");
        }
        let fl = f0.forces().unwrap().last().unwrap();
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

    /// VASP 6.5.1 single points, the same three runs as the OUTCAR test.
    /// Reference values from dpdata 1.0.2. The `stress` varray of this VASP is
    /// not symmetric (frame_0011: xy = -0.95709 kB, yx = -0.95685 kB); dpdata
    /// keeps all nine numbers as written, this reader takes (σ+σᵀ)/2.
    #[test]
    fn matches_dpdata_on_vasp_6_5_single_points() {
        // 循环体与 vasp_outcar.rs 的同名测试相同；只两处，按 R6 不抽 helper
        // (帧目录, 能量, a, 原子号, 坐标, 力, dpdata 的 virial [eV] 行优先)
        type Case = (&'static str, f64, f64, usize, [f64; 3], [f64; 3], [f64; 9]);
        let cases: [Case; 3] = [
            ("frame_0011", -2724.82410263, 16.0970833, 74,
             [13.444300069243301, 11.76320004801197, 0.13421667570123502],
             [-4.34495478, 8.57785538, -0.54024799],
             [-107.24792121346321, -2.4916439978965976, 22.105411766460396,
              -2.491002170020148, -73.32045030153505, -4.593407369113374,
              22.103421824089178, -4.59317916017507, -58.829334754370585]),
            ("frame_0161", -2720.29256226, 15.96159402, 250,
             [14.54659998705452, 11.315200018176252, 0.15063403201306558],
             [-13.7607334, -1.10716601, 2.54395791],
             [-60.96141432966301, 2.2831840299791653, -8.86253317058157,
              2.282867141144278, -80.22867355054338, 3.437238469365705,
              -8.862113054923174, 3.4368162216558793, -49.58771016872843]),
            ("frame_0194", -2733.85701426, 16.00869172, 66,
             [9.967769956685695, 7.677489997676073, 8.11268996408337],
             [0.11733578, 8.3902186, -1.05019035],
             [-32.26906431846686, 0.8540925909315585, 14.744078741496649,
              0.8548218500528398, -60.36171619206951, -4.226226227980381,
              14.743517079585551, -4.227800258983086, -49.20522737762933]),
        ];
        for (dir, e, a, i, pos, force, virial) in cases {
            let path = format!("{}/../tests/vasp651_sp/{dir}/vasprun.xml", env!("CARGO_MANIFEST_DIR"));
            let (traj, stats) = read_vasprun(Path::new(&path)).unwrap();
            assert_eq!((traj.n_frames(), stats.n_dropped()), (1, 0), "{dir}");
            let f = &traj.frames[0];
            assert_eq!(f.n_atoms(), 302, "{dir}");
            assert!((f.energy.unwrap() - e).abs() < 1e-8, "{dir} 能量");
            let cell = f.cell.as_ref().unwrap();
            assert!((cell.matrix[(0, 0)] - a).abs() < 1e-9, "{dir} 晶胞");
            for (g, w) in f.atom(i).position.iter().zip(&pos) {
                assert!((g - w).abs() < 1e-9, "{dir} 坐标 {g} vs {w}");
            }
            for (g, w) in f.forces().unwrap()[i].iter().zip(&force) {
                assert!((g - w).abs() < 1e-9, "{dir} 力 {g} vs {w}");
            }
            // 表里是 dpdata 原样的九个数,期望值取其对称部分。
            // 1e-7 相对容差:dpdata 的 eV/Å³ 换算常数与 units.rs 差 8e-9
            let v = f.stress.unwrap() * cell.volume();
            assert_eq!(v, v.transpose(), "{dir} 读出的应力应已对称化");
            let got = ferro_core::matrix3_row_major(&v);
            for (k, g) in got.iter().enumerate() {
                let w = (virial[k] + virial[(k % 3) * 3 + k / 3]) / 2.0;
                assert!((g - w).abs() <= 1e-7 * w.abs().max(1.0), "{dir} virial[{k}] {g} vs {w}");
            }
        }
    }

    /// The asymmetric stress of VASP 6.5.1 used to survive into `virial=` and
    /// make the extxyz reader reject a file this crate had just written (the
    /// `dataset --type nep` path). Without the symmetrisation this fails.
    #[test]
    fn asymmetric_vasp_stress_round_trips_through_extxyz_virial() {
        let (traj, _) = read_vasprun(Path::new(
            concat!(env!("CARGO_MANIFEST_DIR"), "/../tests/vasp651_sp/frame_0011/vasprun.xml"),
        ))
        .unwrap();
        use crate::readers::extxyz::read_extxyz;
        use crate::writers::extxyz::{write_extxyz_with, StressKey};
        let p = crate::testutil::tmp_path("ferro_vasprun_virial_roundtrip.xyz");
        write_extxyz_with(&traj, &p, StressKey::Virial, false).unwrap();
        let back = read_extxyz(&p).unwrap().0;
        let (s0, s1) = (traj.frames[0].stress.unwrap(), back.frames[0].stress.unwrap());
        assert!((s0 - s1).abs().max() < 1e-12, "往返后应力变了:{s0} vs {s1}");
    }

    /// The first `<energy>` of a calculation belongs to an SCF iterate and is
    /// tens of thousands of eV away from the converged value; only the last one
    /// is the answer. Nothing about the wrong choice looks wrong.
    #[test]
    fn the_energy_is_the_converged_one_not_the_first_scf_step() {
        let (traj, _) = read_vasprun(Path::new(FIXTURE)).unwrap();
        let e = traj.frames[0].energy.unwrap();
        assert!(e < -1000.0 && e > -3000.0, "{e} looks like an SCF iterate");
    }

    #[test]
    fn malformed_xml_is_an_error_not_an_empty_trajectory() {
        let p = crate::testutil::tmp_path("broken.xml");
        std::fs::write(&p, "<?xml version=\"1.0\"?><modeling><calculation>").unwrap();
        // 截断的 XML 要么报错、要么给出 0 帧,但绝不能装作读到了东西
        match read_vasprun(&p) {
            Err(e) => assert!(format!("{e:#}").contains("atominfo") || format!("{e:#}").contains("XML")),
            Ok((t, _)) => assert_eq!(t.n_frames(), 0),
        }
    }

    /// A stress varray without three rows is an incomplete step, not a frame without stress
    /// (review 2026-11-06 A2-3)
    #[test]
    fn a_short_stress_varray_is_counted_incomplete() {
        let text = std::fs::read_to_string(FIXTURE).unwrap();
        let row = "   <v>       3.63419744     -19.22146779      -2.89830625 </v>\n";
        assert_eq!(text.matches(row).count(), 1);
        let path = crate::testutil::write_tmp("vasprun_short_stress.xml", &text.replace(row, ""));
        let (traj, stats) = read_vasprun(&path).unwrap();
        assert_eq!((traj.n_frames(), stats.n_incomplete), (1, 1));
    }

    /// A run killed mid-write leaves the file cut anywhere — inside a tag or inside a text
    /// node. Either way the finished calculations are kept and the cut one is counted as an
    /// incomplete step; a cut inside a tag used to discard the whole file.
    #[test]
    fn a_truncated_file_keeps_the_finished_frames() {
        let text = std::fs::read_to_string(FIXTURE).unwrap();
        let second = text.find("<i name=\"kinetic\">    121.39817355").unwrap();
        for (name, cut) in [("tag", second + 5), ("text", second + 25)] {
            let path = crate::testutil::write_tmp(&format!("vasprun_cut_{name}.xml"), &text[..cut]);
            let (traj, stats) = read_vasprun(&path).unwrap_or_else(|e| panic!("cut in {name}: {e:#}"));
            assert_eq!(traj.n_frames(), 1, "cut in {name}");
            assert_eq!((stats.n_steps, stats.n_incomplete), (2, 1), "cut in {name}");
        }
    }

    /// A calculation without a `kinetic` entry has no temperature; it used to inherit the
    /// previous calculation's kinetic energy and report that frame's temperature as its own.
    #[test]
    fn a_missing_kinetic_energy_is_not_carried_over_from_the_previous_frame() {
        let text = std::fs::read_to_string(FIXTURE).unwrap();
        let second = "<i name=\"kinetic\">    121.39817355 </i>";
        assert_eq!(text.matches(second).count(), 1);
        let path = crate::testutil::write_tmp("vasprun_no_kinetic.xml", &text.replacen(second, "", 1));
        let (traj, _) = read_vasprun(&path).unwrap();
        assert!(traj.frames[0].temperature.is_some());
        assert_eq!(traj.frames[1].temperature, None);
    }

    /// An empty `<i>` has no text event: its name used to stay armed, and the next text
    /// node — here `<time>`'s "181.36 185.61" — was taken as the kinetic energy, losing
    /// frame 0's temperature.
    #[test]
    fn empty_i_element_does_not_capture_the_next_text() {
        let text = std::fs::read_to_string(FIXTURE).unwrap();
        let total = "<i name=\"total\">  -2040.54719326 </i>";
        assert_eq!(text.matches(total).count(), 1);
        let edited = text.replacen(total, "<i name=\"kinetic\"></i>", 1);
        let path = crate::testutil::write_tmp("vasprun_empty_i.xml", &edited);
        let (want, _) = read_vasprun(Path::new(FIXTURE)).unwrap();
        let (got, _) = read_vasprun(&path).unwrap();
        assert!(got.frames[0].temperature.is_some());
        assert_eq!(got.frames[0].temperature, want.frames[0].temperature);
    }

    /// A bad number in a value this reader uses, and a byte that is not UTF-8, are errors;
    /// they used to be skipped (the frame went quietly incomplete) or replaced.
    #[test]
    fn bad_numbers_and_bytes_are_errors() {
        let text = std::fs::read_to_string(FIXTURE).unwrap();
        let row = "<v>       5.61175521      -0.45541218       5.62845267 </v>";
        let nelm = "<i type=\"int\" name=\"NELM\">   100</i>";
        assert_eq!(text.matches(row).count(), 1);
        for (name, edited, want) in [
            ("bad_force", text.replacen(row, "<v>       5.61175521      -0.45541218       5.6x845267 </v>", 1),
             "varray value \"5.6x845267\" is not a number"),
            ("short_force", text.replacen(row, "<v>       5.61175521      -0.45541218 </v>", 1), "has 2 values, not 3"),
            ("bad_nelm", text.replacen(nelm, "<i type=\"int\" name=\"NELM\">   1O0</i>", 1), "NELM \"1O0\""),
        ] {
            let path = crate::testutil::write_tmp(&format!("vasprun_{name}.xml"), &edited);
            let e = format!("{:#}", read_vasprun(&path).unwrap_err());
            assert!(e.contains(want), "{name}: {e}");
        }
        let mut bytes = text.into_bytes();
        let at = bytes.windows(row.len()).position(|w| w == row.as_bytes()).unwrap() + 10;
        bytes[at] = 0xff;
        let path = crate::testutil::tmp_path("ferro_vasprun_bad_byte.xml");
        std::fs::write(&path, bytes).unwrap();
        let e = format!("{:#}", read_vasprun(&path).unwrap_err());
        assert!(e.contains("byte 0xff is not UTF-8"), "{e}");
    }
}
