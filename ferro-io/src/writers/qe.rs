use std::fmt::Write as _;
use ferro_core::Frame;
use anyhow::{Context, Result};

/// 单帧的三张结构卡片 `[(卡片名, 原文)]`：ATOMIC_SPECIES（质量 + 赝势文件名）、CELL_PARAMETERS 与
/// ATOMIC_POSITIONS（均为 angstrom），原文从卡片首行起、以换行结尾。`pseudos` 为 元素 → 赝势文件名，
/// 须覆盖帧内全部元素。`ferro.jobs` 生成 QE 任务用它，QE 卡片格式只此一份。
/// Each card is returned on its own, so the caller never has to split a joined text again.
pub fn qe_structure_cards(frame: &Frame, pseudos: &[(String, String)]) -> Result<Vec<(String, String)>> {
    // 无胞时 pw.x 拒收 ibrav = 0 的输入（审查 D-M9），宁可在写之前报错
    let cell = frame.cell.as_ref().context("QE input needs a periodic cell (CELL_PARAMETERS)")?;

    let mut w = String::new();
    writeln!(w, "ATOMIC_SPECIES")?;
    for elem in frame.unique_elements() {
        // The atoms' own mass (explicit, else the table; neither is an error, never 1 amu,
        // review A-3). QE has one species per element, so one element with two masses
        // cannot be written
        let mut masses = frame.atoms().filter(|a| a.element == elem).map(|a| a.effective_mass());
        let mass = masses.next().expect("an element of the frame has an atom")?;
        for m in masses {
            if m? != mass {
                anyhow::bail!("{elem} atoms carry different masses; QE writes one species (one mass) per element");
            }
        }
        let pseudo = pseudos.iter().find(|(e, _)| *e == elem)
            .map(|(_, p)| p)
            .with_context(|| format!("no pseudopotential given for {elem}"))?;
        writeln!(w, "  {elem}  {mass:.4}  {pseudo}")?;
    }
    let species = std::mem::take(&mut w);

    writeln!(w, "CELL_PARAMETERS {{angstrom}}")?;
    for i in 0..3 {
        let r = cell.matrix.row(i);
        writeln!(w, "  {:.10}  {:.10}  {:.10}", r[0], r[1], r[2])?;
    }
    let cell_parameters = std::mem::take(&mut w);

    writeln!(w, "ATOMIC_POSITIONS {{angstrom}}")?;
    for atom in frame.atoms() {
        writeln!(w, "  {}  {:.10}  {:.10}  {:.10}",
            atom.element,
            atom.position.x, atom.position.y, atom.position.z)?;
    }
    Ok(vec![
        ("ATOMIC_SPECIES".into(), species),
        ("CELL_PARAMETERS".into(), cell_parameters),
        ("ATOMIC_POSITIONS".into(), w),
    ])
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::qe::parse_qe_input_parts;
    use crate::readers::read_qe_input;
    use crate::testutil::write_tmp as tmp;
    use ferro_core::Atom;
    use nalgebra::Vector3;

    #[test]
    fn test_cards_round_trip() {
        let cards = qe_structure_cards(&crate::testutil::bcc_fe(), &[("Fe".into(), "Fe.pbe.UPF".into())]).unwrap();
        let names: Vec<&str> = cards.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(names, ["ATOMIC_SPECIES", "CELL_PARAMETERS", "ATOMIC_POSITIONS"]);
        for (name, text) in &cards {
            assert!(text.starts_with(name.as_str()) && text.ends_with('\n'), "{text}");
        }
        let cards: String = cards.iter().map(|(_, t)| format!("{t}\n")).collect();
        assert!(cards.contains("  Fe  55.8450  Fe.pbe.UPF"), "{cards}");
        let text = format!("&SYSTEM\n  ibrav = 0, nat = 2, ntyp = 1\n/\n{cards}");
        let f = read_qe_input(&tmp("bcc_cards.qe", &text)).unwrap().frames.remove(0);
        assert_eq!(f.n_atoms(), 2);
        assert!((f.atom(1).position - Vector3::repeat(1.435)).norm() < 1e-9);
        let [a, ..] = f.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-9);
        // 卡片部分不含 namelist，拆分后原样落在 cards 里
        assert_eq!(parse_qe_input_parts(&text).unwrap().cards.trim_end(), cards.trim_end());
    }

    /// ATOMIC_SPECIES fell back to 1 amu for an unknown element and ignored explicit
    /// masses (review R2-1; review A-3 removed that fallback from the other writers).
    #[test]
    fn test_species_mass_is_the_atoms_mass() {
        let with_masses = |masses: [Option<f64>; 2], element: &str| {
            let mut f = crate::testutil::bcc_fe();
            let atoms: Vec<Atom> = f.atoms().zip(masses).map(|(a, m)| {
                Atom { element: element.into(), mass: m, ..a.to_atom() }
            }).collect();
            f.set_atoms(atoms);
            f
        };
        let pseudo = |e: &str| [(e.to_string(), "X.UPF".to_string())];
        let species = |f: &Frame, e: &str| qe_structure_cards(f, &pseudo(e)).map(|c| c[0].1.clone());

        // An explicit mass (an isotope) is written, not the table's
        let s = species(&with_masses([Some(57.9353), Some(57.9353)], "Fe"), "Fe").unwrap();
        assert!(s.contains("  Fe  57.9353  X.UPF"), "{s}");
        // No mass anywhere is an error, never 1 amu
        let e = species(&with_masses([None, None], "Fe1"), "Fe1").unwrap_err();
        assert!(format!("{e:#}").contains("no mass for element 'Fe1'"), "{e:#}");
        // One species per element: two masses for one element cannot be written
        let e = species(&with_masses([Some(55.845), Some(57.9353)], "Fe"), "Fe").unwrap_err();
        assert!(format!("{e:#}").contains("different masses"), "{e:#}");
    }

    #[test]
    fn test_structure_cards_refuse_incomplete_input() {
        let e = qe_structure_cards(&crate::testutil::bcc_fe(), &[]).unwrap_err();
        assert!(format!("{e:#}").contains("no pseudopotential given for Fe"), "{e:#}");
        let mut mol = Frame::new();
        mol.add_atom(Atom::new("O", Vector3::zeros()));
        let e = qe_structure_cards(&mol, &[("O".into(), "O.UPF".into())]).unwrap_err();
        assert!(format!("{e:#}").contains("periodic cell"), "无胞应报错（D-M9）：{e:#}");
    }
}
