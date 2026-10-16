use std::fmt::Write as _;
use ferro_core::Frame;
use anyhow::{Context, Result};

/// 单帧的三张结构卡片：ATOMIC_SPECIES（质量 + 赝势文件名）、CELL_PARAMETERS 与
/// ATOMIC_POSITIONS（均为 angstrom），每张后跟一空行。`pseudos` 为 元素 → 赝势文件名，
/// 须覆盖帧内全部元素。`ferro.jobs` 生成 QE 任务用它，QE 卡片格式只此一份
pub fn qe_structure_cards(frame: &Frame, pseudos: &[(String, String)]) -> Result<String> {
    // 无胞时 pw.x 拒收 ibrav = 0 的输入（审查 D-M9），宁可在写之前报错
    let cell = frame.cell.as_ref().context("QE input needs a periodic cell (CELL_PARAMETERS)")?;
    let mut w = String::new();

    writeln!(w, "ATOMIC_SPECIES")?;
    for elem in frame.unique_elements() {
        let mass = ferro_core::data::elements::by_symbol(&elem)
            .map(|e| e.atomic_mass)
            .unwrap_or(1.0);
        let pseudo = pseudos.iter().find(|(e, _)| *e == elem)
            .map(|(_, p)| p)
            .with_context(|| format!("no pseudopotential given for {elem}"))?;
        writeln!(w, "  {elem}  {mass:.4}  {pseudo}")?;
    }
    writeln!(w)?;

    writeln!(w, "CELL_PARAMETERS {{angstrom}}")?;
    for i in 0..3 {
        let r = cell.matrix.row(i);
        writeln!(w, "  {:.10}  {:.10}  {:.10}", r[0], r[1], r[2])?;
    }
    writeln!(w)?;

    writeln!(w, "ATOMIC_POSITIONS {{angstrom}}")?;
    for atom in &frame.atoms {
        writeln!(w, "  {}  {:.10}  {:.10}  {:.10}",
            atom.element,
            atom.position.x, atom.position.y, atom.position.z)?;
    }
    writeln!(w)?;
    Ok(w)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::qe::parse_qe_input_parts;
    use crate::readers::read_qe_input;
    use crate::testutil::write_tmp as tmp;
    use ferro_core::{Atom, Cell};
    use nalgebra::Vector3;

    fn bcc() -> Frame {
        let cell = Cell::from_lengths_angles(2.87, 2.87, 2.87, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(1.435, 1.435, 1.435)));
        frame
    }

    #[test]
    fn test_cards_round_trip() {
        let cards = qe_structure_cards(&bcc(), &[("Fe".into(), "Fe.pbe.UPF".into())]).unwrap();
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

    #[test]
    fn test_structure_cards_refuse_incomplete_input() {
        let e = qe_structure_cards(&bcc(), &[]).unwrap_err();
        assert!(format!("{e:#}").contains("no pseudopotential given for Fe"), "{e:#}");
        let mut mol = Frame::new();
        mol.add_atom(Atom::new("O", Vector3::zeros()));
        let e = qe_structure_cards(&mol, &[("O".into(), "O.UPF".into())]).unwrap_err();
        assert!(format!("{e:#}").contains("periodic cell"), "无胞应报错（D-M9）：{e:#}");
    }
}
