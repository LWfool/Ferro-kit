use std::fmt::Write as _;
use ferro_core::Frame;
use anyhow::{Context, Result};

/// 单帧的 CP2K `&CELL` 与 `&COORD` 两段（Å，kind 名取元素符号），供 `ferro.jobs`
/// 嵌进 `&FORCE_EVAL/&SUBSYS`；`&KIND`（基组、赝势）由调用方按元素另写。
/// 与 [`crate::readers::read_cp2k_inp`] 往返
pub fn cp2k_structure_sections(frame: &Frame) -> Result<String> {
    let cell = frame.cell.as_ref().context("CP2K input needs a cell (&CELL)")?;
    // PERIODIC 取值见 input_cp2k_subsys.F：X Y Z XY XZ YZ XYZ NONE
    let axes: String = ['X', 'Y', 'Z'].iter().zip(frame.pbc).filter(|(_, p)| *p).map(|(c, _)| *c).collect();
    let periodic = if axes.is_empty() { "NONE".to_string() } else { axes };

    let mut w = String::new();
    writeln!(w, "&CELL")?;
    for (name, i) in [("A", 0), ("B", 1), ("C", 2)] {
        let r = cell.matrix.row(i);
        writeln!(w, "  {name}  {:.10}  {:.10}  {:.10}", r[0], r[1], r[2])?;
    }
    writeln!(w, "  PERIODIC {periodic}")?;
    writeln!(w, "&END CELL")?;
    writeln!(w, "&COORD")?;
    for atom in frame.atoms() {
        writeln!(w, "  {}  {:.10}  {:.10}  {:.10}",
            atom.element, atom.position.x, atom.position.y, atom.position.z)?;
    }
    writeln!(w, "&END COORD")?;
    Ok(w)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::readers::read_cp2k_inp;
    use crate::testutil::write_tmp as tmp;
    use ferro_core::{Atom, Cell};
    use nalgebra::{Matrix3, Vector3};

    #[test]
    fn test_round_trip_triclinic() {
        let m = Matrix3::new(10.0, 0.0, 0.0, 2.0, 12.0, 0.0, -3.0, 1.0, 14.0);
        let mut frame = Frame::with_cell(Cell::from_matrix(m), [true, true, false]);
        frame.add_atom(Atom::new("Zr", Vector3::new(1.0, 2.0, 3.0)));
        frame.add_atom(Atom::new("O", Vector3::new(4.5, 5.25, 6.125)));
        let text = format!("&FORCE_EVAL\n&SUBSYS\n{}&END SUBSYS\n&END FORCE_EVAL\n",
                           cp2k_structure_sections(&frame).unwrap());
        let back = read_cp2k_inp(&tmp("rt.inp", &text)).unwrap();
        let f = back.first().unwrap();
        assert!((f.cell.as_ref().unwrap().matrix - m).norm() < 1e-9, "胞应原样读回（行 = 晶格矢量）");
        assert_eq!(f.pbc, [true, true, false], "PERIODIC XY 应读回");
        assert_eq!((f.atom(0).element, f.atom(1).element), ("Zr", "O"));
        assert!((f.atom(1).position - Vector3::new(4.5, 5.25, 6.125)).norm() < 1e-9);
    }

    #[test]
    fn test_no_cell_is_an_error() {
        let mut frame = Frame::new();
        frame.add_atom(Atom::new("O", Vector3::zeros()));
        assert!(cp2k_structure_sections(&frame).is_err(), "无胞应报错，不写不完整的 &CELL");
    }
}
