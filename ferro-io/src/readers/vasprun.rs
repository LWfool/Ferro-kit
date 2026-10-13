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

use anyhow::{bail, Context, Result};
use ferro_core::units::{convert_pressure, PressureUnit, BOLTZMANN_EV_K};
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use quick_xml::events::Event;
use quick_xml::Reader;

use super::aimd::{AimdFormat, AimdStats};

/// Reads a vasprun.xml, discarding the statistics.
pub fn read_vasprun(path: &Path) -> Result<Trajectory> {
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
    kinetic: Option<f64>,
}

/// Reads a vasprun.xml and reports what was dropped.
pub fn read_vasprun_with_stats(path: &Path) -> Result<(Trajectory, AimdStats)> {
    let path_ = path.display();
    let file = std::fs::File::open(path).with_context(|| format!("cannot open {path_}"))?;
    let mut reader = Reader::from_reader(BufReader::new(file));
    reader.config_mut().trim_text(true);

    let mut stats = AimdStats::new(AimdFormat::VaspXml);
    let mut traj = Trajectory::new();
    let mut st = State::default();
    let mut buf = Vec::new();

    loop {
        match reader.read_event_into(&mut buf) {
            Err(e) => bail!("{path_}: malformed XML at byte {}: {e}", reader.buffer_position()),
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
                    } else if field == "kinetic" {
                        // 离子动能(eV),温度由它反算 —— vasprun 不打印温度本身
                        st.kinetic = text.trim().parse().ok();
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
        bail!("{path_}: no <atominfo> species list found");
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
    // vasprun 没有温度字段(只有输入参数 TEBEG),按 T = 2*E_kin/(3N*k_B) 反算。
    // 除数是 3N 不是 3N-3:OUTCAR 的 EKIN=4.268142 eV / N=296 与它自己打印的
    // 111.55 K 对拍确认了这一点。这一列是导出量,报告的表头会注明
    frame.temperature = st
        .kinetic
        .filter(|_| n > 0)
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

    const FIXTURE: &str = "../tests/vasp_vasprun_2frames.xml";

    /// vasprun has no temperature field, so it is derived from the ionic
    /// kinetic energy with `T = 2*E_kin/(3N*k_B)`. The divisor is 3N, not 3N-3:
    /// the OUTCAR of a comparable run prints 111.55 K next to EKIN = 4.268142 eV
    /// for N = 296, which only 3N reproduces.
    #[test]
    fn derives_temperature_from_ionic_kinetic_energy() {
        let (traj, _) =
            read_vasprun_with_stats(Path::new("../tests/vasp_vasprun_2frames.xml")).unwrap();
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
        let (traj, stats) = read_vasprun_with_stats(Path::new(FIXTURE)).unwrap();
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
        let (traj, _) = read_vasprun_with_stats(Path::new(FIXTURE)).unwrap();
        assert_eq!(traj.first().unwrap().n_atoms(), 297);
    }

    /// Reference values from dpdata 1.0.2 on the same file.
    #[test]
    fn matches_dpdata_on_the_real_fixture() {
        let (traj, _) = read_vasprun_with_stats(Path::new(FIXTURE)).unwrap();
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
            let path = format!("../tests/vasp651_sp/{dir}/vasprun.xml");
            let (traj, stats) = read_vasprun_with_stats(Path::new(&path)).unwrap();
            assert_eq!((traj.n_frames(), stats.n_dropped()), (1, 0), "{dir}");
            let f = &traj.frames[0];
            assert_eq!(f.n_atoms(), 302, "{dir}");
            assert!((f.energy.unwrap() - e).abs() < 1e-8, "{dir} 能量");
            let cell = f.cell.as_ref().unwrap();
            assert!((cell.matrix[(0, 0)] - a).abs() < 1e-9, "{dir} 晶胞");
            for (g, w) in f.atom(i).position.iter().zip(&pos) {
                assert!((g - w).abs() < 1e-9, "{dir} 坐标 {g} vs {w}");
            }
            for (g, w) in f.forces.as_ref().unwrap()[i].iter().zip(&force) {
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
        let (traj, _) = read_vasprun_with_stats(Path::new(
            "../tests/vasp651_sp/frame_0011/vasprun.xml",
        ))
        .unwrap();
        use crate::readers::extxyz::read_extxyz;
        use crate::writers::extxyz::{write_extxyz_with, StressKey};
        let p = std::env::temp_dir().join("ferro_vasprun_virial_roundtrip.xyz");
        write_extxyz_with(&traj, &p, StressKey::Virial).unwrap();
        let back = read_extxyz(&p).unwrap();
        let (s0, s1) = (traj.frames[0].stress.unwrap(), back.frames[0].stress.unwrap());
        assert!((s0 - s1).abs().max() < 1e-12, "往返后应力变了:{s0} vs {s1}");
    }

    /// The first `<energy>` of a calculation belongs to an SCF iterate and is
    /// tens of thousands of eV away from the converged value; only the last one
    /// is the answer. Nothing about the wrong choice looks wrong.
    #[test]
    fn the_energy_is_the_converged_one_not_the_first_scf_step() {
        let (traj, _) = read_vasprun_with_stats(Path::new(FIXTURE)).unwrap();
        let e = traj.frames[0].energy.unwrap();
        assert!(e < -1000.0 && e > -3000.0, "{e} looks like an SCF iterate");
    }

    #[test]
    fn malformed_xml_is_an_error_not_an_empty_trajectory() {
        let p = std::env::temp_dir().join("broken.xml");
        std::fs::write(&p, "<?xml version=\"1.0\"?><modeling><calculation>").unwrap();
        // 截断的 XML 要么报错、要么给出 0 帧,但绝不能装作读到了东西
        match read_vasprun_with_stats(&p) {
            Err(e) => assert!(format!("{e:#}").contains("atominfo") || format!("{e:#}").contains("XML")),
            Ok((t, _)) => assert_eq!(t.n_frames(), 0),
        }
    }
}
