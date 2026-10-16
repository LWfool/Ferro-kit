use std::path::Path;
use ferro_core::{Frame, Trajectory};
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Context, Result};

/// 写 Quantum ESPRESSO pw.x 输入文件（ibrav=0，ATOMIC_POSITIONS angstrom）。
///
/// 这是 **格式转换** 用途（`fe-convert`，与 [`crate::readers::qe::read_qe_input`]
/// 往返）：只写结构骨架 + 占位 `&CONTROL/&SYSTEM/&ELECTRONS`。需要完整计算
/// 设置（任务/泛函/k 点/自旋等）请用 `ferro_workflow::QeJobBuilder`（`fe-job -s qe`）。
pub fn write_qe_input(trajectory: &Trajectory, path: &Path) -> Result<()> {
    let path_ = path.display();
    let frame = trajectory.first().context("trajectory is empty")?;

    let file = File::create(path).with_context(|| format!("cannot create {path_}"))?;
    let mut w = BufWriter::new(file);

    let prefix = trajectory.metadata.source.as_deref().unwrap_or("ferro");

    // Unique element labels in first-appearance order
    let elem_order = frame.unique_elements();
    let ntyp = elem_order.len();
    let nat = frame.n_atoms();

    // &CONTROL
    writeln!(w, "&CONTROL")?;
    writeln!(w, "  calculation = 'scf',")?;
    writeln!(w, "  prefix = '{prefix}',")?;
    writeln!(w, "/")?;
    writeln!(w)?;

    // &SYSTEM
    writeln!(w, "&SYSTEM")?;
    writeln!(w, "  ibrav = 0,")?;
    writeln!(w, "  nat = {nat},")?;
    writeln!(w, "  ntyp = {ntyp},")?;
    let charge = frame.charge;
    if charge != 0 {
        writeln!(w, "  tot_charge = {charge},")?;
    }
    let mult = frame.multiplicity;
    if mult != 1 {
        writeln!(w, "  tot_magnetization = {},", mult - 1)?;
    }
    writeln!(w, "/")?;
    writeln!(w)?;

    // &ELECTRONS (empty placeholder)
    writeln!(w, "&ELECTRONS")?;
    writeln!(w, "/")?;
    writeln!(w)?;

    let pseudos: Vec<(String, String)> = elem_order.iter()
        .map(|e| (e.clone(), format!("{e}.UPF")))
        .collect();
    write!(w, "{}", qe_structure_cards(frame, &pseudos)?)?;

    // K_POINTS gamma (placeholder)
    writeln!(w, "K_POINTS {{gamma}}")?;

    w.flush()?;
    Ok(())
}

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
    use crate::readers::qe::read_qe_input;
    use ferro_core::{Atom, Cell, Frame, Trajectory};
    use nalgebra::Vector3;

    fn bcc_traj() -> Trajectory {
        let cell = Cell::from_lengths_angles(2.87, 2.87, 2.87, 90.0, 90.0, 90.0).unwrap();
        let mut frame = Frame::with_cell(cell, [true; 3]);
        frame.add_atom(Atom::new("Fe", Vector3::new(0.0, 0.0, 0.0)));
        frame.add_atom(Atom::new("Fe", Vector3::new(1.435, 1.435, 1.435)));
        Trajectory::from_frame(frame)
    }

    #[test]
    fn test_roundtrip() {
        let path = std::env::temp_dir().join("bcc_rt.qe");
        let p = &path;
        let orig = bcc_traj();
        write_qe_input(&orig, p).unwrap();

        let loaded = read_qe_input(p).unwrap();
        let f = loaded.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
        assert!((f.atom(1).position.x - 1.435).abs() < 1e-6);
        let [a, ..] = f.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-4);
    }

    #[test]
    fn test_structure_cards_refuse_incomplete_input() {
        let f = bcc_traj().frames.remove(0);
        let e = qe_structure_cards(&f, &[]).unwrap_err();
        assert!(format!("{e:#}").contains("no pseudopotential given for Fe"), "{e:#}");
        let mut mol = Frame::new();
        mol.add_atom(Atom::new("O", Vector3::zeros()));
        let e = qe_structure_cards(&mol, &[("O".into(), "O.UPF".into())]).unwrap_err();
        assert!(format!("{e:#}").contains("periodic cell"), "无胞应报错（D-M9）：{e:#}");
    }
}
