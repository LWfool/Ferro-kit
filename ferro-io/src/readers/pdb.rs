use std::path::Path;
use crate::Notice;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::Vector3;
use anyhow::{bail, Context, Result};
use ferro_core::data::elements::by_symbol;
use super::util::capitalize_symbol;

/// Reads a PDB file, multi-model (MODEL / ENDMDL records), returning what the reader
/// noticed: a `CRYST1` placeholder taken as no
/// cell, and a last MODEL with no `ENDMDL` or last END-split block with no `END`
/// (kept when it has the atom count of the frame before it, dropped otherwise).
///
/// A frame ends at `ENDMDL` or at `END`: CP2K and VMD write trajectories as blocks
/// closed by `END` with no MODEL records (as ASE reads them). Those used to run
/// together into one frame holding every block's atoms.
pub fn read_pdb(path: &Path) -> Result<(Trajectory, Vec<Notice>)> {
    let mut traj = Trajectory::new();
    let mut warnings = Vec::new();
    let mut current = Frame::new();
    let mut open_model: Option<usize> = None; // line of an open MODEL, until its ENDMDL
    let mut cell: Option<Cell> = None;
    let close = |frame: &mut Frame, cell: &Option<Cell>| {
        if let Some(c) = cell {
            frame.cell = Some(c.clone());
            frame.pbc = [true; 3];
        }
    };

    let mut lines = super::util::text_lines(path)?.enumerate().peekable();
    while let Some((ln, line)) = lines.next() {
        let line = line?;
        let tag = col(&line, 0, 6);

        match tag.trim_end() {
            "HEADER" => {
                // 记录名 6 字节是 ASCII，6 处必是字符边界；以前从 10 切，多字节标题会 panic
                let source = line[6..].trim();
                if !source.is_empty() {
                    traj.metadata.title = Some(source.to_string());
                }
            }
            "CRYST1" => {
                cell = Some(parse_cryst1(&line).with_context(|| format!("line {}: {line:?}", ln + 1))?);
                // PDB format: a structure not determined crystallographically carries the placeholder
                // CRYST1 1 1 1 90 90 90 (MDAnalysis reads it as no cell too)
                if cell.as_ref().is_some_and(is_placeholder_cell) {
                    if !warnings.iter().any(|w: &Notice| w.message.contains("CRYST1")) {
                        warnings.push(Notice::warning(Some(ln + 1),
                            "CRYST1 1 1 1 90 90 90 is the PDB placeholder for \
                             a structure with no unit cell; read as non-periodic"));
                    }
                    cell = None;
                }
            }
            "MODEL" => open_model = Some(ln + 1),
            "ENDMDL" | "END" => {
                open_model = None;
                // An END after MODEL … ENDMDL makes no extra frame
                if current.n_atoms() > 0 {
                    let mut frame = std::mem::replace(&mut current, Frame::new());
                    close(&mut frame, &cell);
                    traj.add_frame(frame);
                }
            }
            "ATOM" | "HETATM" => {
                let last = lines.peek().is_none();
                match parse_atom_record(&line, last).with_context(|| format!("line {}: {line:?}", ln + 1))? {
                    Some(atom) => current.add_atom(atom),
                    None => {
                        let why = "the last ATOM line is cut short";
                        super::util::drop_truncated(&mut warnings, traj.n_frames(), ln + 1, why)?;
                        return Ok((traj, warnings));
                    }
                }
            }
            _ => {}
        }
    }

    if current.n_atoms() > 0 {
        let k = traj.n_frames();
        match open_model {
            // A one-frame file with no MODEL / END record: fine
            None if k == 0 => {}
            // Frames split by END (CP2K / VMD write no MODEL): a last block with no END is held to the
            // frame before, as a last MODEL with no ENDMDL is below
            None => {
                let n = current.n_atoms();
                let prev = traj.frames[k - 1].n_atoms();
                if prev != n {
                    warnings.push(Notice::warning(None, format!(
                        "last block (frame {k}) has no END and {n} atoms where the frame \
                         before has {prev}; dropped as incomplete")));
                    return Ok((traj, warnings));
                }
                warnings.push(Notice::warning(None, format!("last block (frame {k}) has no END; kept")));
            }
            // A last MODEL with no ENDMDL: trusted only with the previous frame's atom count, otherwise
            // most likely cut off
            Some(model_line) => {
                let n = current.n_atoms();
                match traj.frames.last().map(|f| f.n_atoms()) {
                    Some(prev) if prev != n => {
                        warnings.push(Notice::warning(Some(model_line), format!(
                            "last MODEL (frame {k}) has no ENDMDL and \
                             {n} atoms where the frame before has {prev}; dropped as incomplete")));
                        return Ok((traj, warnings));
                    }
                    _ => warnings.push(Notice::warning(Some(model_line), format!(
                        "last MODEL (frame {k}) has no ENDMDL; kept"))),
                }
            }
        }
        close(&mut current, &cell);
        traj.add_frame(current);
    }

    Ok((traj, warnings))
}

fn is_placeholder_cell(c: &Cell) -> bool {
    let near = |x: f64, y: f64| (x - y).abs() < 1e-6;
    c.lengths().iter().all(|&l| near(l, 1.0)) && c.angles().iter().all(|&g| near(g, 90.0))
}

/// 第 a..b 字节列（0 起、左闭右开），行短或落在多字节字符中间时为空串。PDB 是按列的
/// ASCII 格式，非 ASCII 只会出现在不该出现的地方，空串让各字段自己报错或回落
fn col(line: &str, a: usize, b: usize) -> &str {
    line.get(a..b.min(line.len())).unwrap_or("")
}

/// CRYST1 format: CRYST1  a  b  c  α  β  γ  sgroup  z
///
/// A record that does not give a cell is an error (as in ASE): it used to read as no cell, so
/// the file went on as non-periodic without a word (scan B-6).
fn parse_cryst1(line: &str) -> Result<Cell> {
    let field = |name: &str, a: usize, b: usize| -> Result<f64> {
        let s = col(line, a, b).trim();
        s.parse().with_context(|| format!("CRYST1 {name} in columns {}-{b} is {s:?}, not a number", a + 1))
    };
    let (a, b, c) = (field("a", 6, 15)?, field("b", 15, 24)?, field("c", 24, 33)?);
    let (al, be, ga) = (field("alpha", 33, 40)?, field("beta", 40, 47)?, field("gamma", 47, 54)?);
    Ok(Cell::from_lengths_angles(a, b, c, al, be, ga)?)
}

/// One ATOM / HETATM line → an atom. Bad coordinates or a line shorter than 54 columns are
/// an error (as in ASE) — atoms used to be dropped silently, fewer atoms and exit code 0.
///
/// `None` when this is the file's last line and its coordinates are missing or cut inside a
/// number: an MD run cut off, not a bad file (scan B-2, as the dump reader judges it).
fn parse_atom_record(line: &str, last: bool) -> Result<Option<Atom>> {
    let xyz: Option<Vec<f64>> = [(30, 38), (38, 46), (46, 54)].iter()
        .map(|&(a, b)| col(line, a, b).trim().parse().ok())
        .collect();
    let Some(&[x, y, z]) = xyz.as_deref() else {
        if last { return Ok(None); }
        bail!("invalid or missing coordinates in columns 31-54");
    };
    // Columns 77-78 are the standard element field; RCSB writes it upper-case (`FE`) — made `Fe`,
    // not checked further (the element table stops at Rn, a check would refuse U and the like).
    // When blank, the atom name gives it
    let sym = col(line, 76, 78).trim();
    let element = if sym.is_empty() { element_from_atom_name(col(line, 12, 16))? } else { capitalize_symbol(sym) };
    Ok(Some(Atom::new(element, Vector3::new(x, y, z))))
}

/// 原子名（13–16 列）→ 元素，按 PDB 的列对齐约定：元素符号右对齐在 13–14 列，
/// 所以 13 列为空格或数字时元素是 14 列的单字母（` CA ` 是 α 碳），否则取 13–14 列
/// （`CA  ` 是钙）；两字母不是元素时退到 13 列单字母（`HD21` 是氢）。ASE 不看对齐、
/// 先试两字母，会把 ` CA ` 读成钙。仍有歧义：写在 13 列的 `HG21` 读成汞 —— 规范文件
/// 写了 77–78 列就走不到这里
fn element_from_atom_name(name: &str) -> Result<String> {
    let b = name.as_bytes();
    let candidates: Vec<String> = match b.first() {
        Some(c) if c.is_ascii_whitespace() || c.is_ascii_digit() => vec![col(name, 1, 2).to_string()],
        Some(_) => vec![col(name, 0, 2).to_string(), col(name, 0, 1).to_string()],
        None => vec![],
    };
    candidates.into_iter()
        .map(|c| capitalize_symbol(&c))
        .find(|c| by_symbol(c).is_some())
        .with_context(|| format!(
            "no element in columns 77-78 and atom name {name:?} does not give one; \
             write the element symbol in columns 77-78"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER_PDB: &str = "\
HEADER    Water molecule
ATOM      1  O   UNK A   1       0.000   0.000   0.119  1.00  0.00           O
ATOM      2  H   UNK A   1       0.000   0.763  -0.477  1.00  0.00           H
ATOM      3  H   UNK A   1       0.000  -0.763  -0.477  1.00  0.00           H
END
";

    const MULTI_PDB: &str = "\
MODEL        1
ATOM      1  C   UNK A   1       0.000   0.000   0.000  1.00  0.00           C
ENDMDL
MODEL        2
ATOM      1  C   UNK A   1       1.400   0.000   0.000  1.00  0.00           C
ENDMDL
";

    use crate::testutil::write_tmp;

    #[test]
    fn test_single_frame() {
        let path = write_tmp("test_water.pdb", WATER_PDB);
        let traj = read_pdb(&path).unwrap().0;
        assert_eq!(traj.n_frames(), 1);
        let frame = traj.first().unwrap();
        assert_eq!(frame.n_atoms(), 3);
        assert_eq!(frame.atom(0).element, "O");
        assert_eq!(traj.metadata.title.as_deref(), Some("Water molecule"));
    }

    /// 原子名 + 77–78 列 → 一行 ATOM（列位照 PDB 规范）
    fn atom_line(name: &str, elem: &str) -> String {
        format!("ATOM      1 {name:<4} UNK A   1       0.000   0.000   0.119  1.00  0.00          {elem:>2}")
    }
    fn element_of(line: &str) -> String {
        read_pdb(&write_tmp("test_elem.pdb", &format!("{line}\n"))).unwrap().0.frames[0].element(0).to_string()
    }

    #[test]
    fn test_element_from_columns_or_atom_name() {
        assert_eq!(element_of(&atom_line("FE", "FE")), "Fe", "RCSB 的大写元素");
        assert_eq!(element_of(&atom_line(" O", "")), "O", "77–78 列空，从原子名推");
        assert_eq!(element_of(&atom_line(" CA", "")), "C", "13 列空格：α 碳");
        assert_eq!(element_of(&atom_line("CA", "")), "Ca", "13 列起：钙");
        assert_eq!(element_of(&atom_line("HD21", "")), "H", "两字母非元素，退单字母");
        assert_eq!(element_of(&atom_line("1HB", "")), "H", "13 列数字");
        // 行短于 77 列（packmol 等）同样从原子名推
        assert_eq!(element_of(&atom_line(" CA", "")[..66]), "C");
        let e = read_pdb(&write_tmp("test_elem_bad.pdb", &format!("{}\n", atom_line(" 1", "")))).unwrap_err();
        assert!(format!("{e:#}").contains("columns 77-78"), "{e:#}");
    }

    #[test]
    fn test_bad_or_short_atom_line_is_an_error_naming_the_line() {
        for bad in [WATER_PDB.replace("0.763  -0.477", "abc    -0.477"),
                    WATER_PDB.replace("END\n", "ATOM      4  H   UNK A   1       0.000\nEND\n")] {
            let e = read_pdb(&write_tmp("test_bad.pdb", &bad)).unwrap_err();
            let msg = format!("{e:#}");
            assert!(msg.contains("line ") && msg.contains("coordinates"), "应点名行号：{msg}");
        }
        assert!(format!("{:#}", read_pdb(&write_tmp("test_bad.pdb",
            &WATER_PDB.replace("0.763  -0.477", "abc    -0.477"))).unwrap_err()).contains("line 3"));
    }

    /// The file's last ATOM line cut short is an MD run cut off: the frame it belongs to is
    /// dropped, and when that was the only frame it is an error (scan B-2).
    #[test]
    fn test_cut_last_atom_line() {
        let cut = format!("{WATER_PDB}ATOM      1  O   UNK A   1       0.000   0.0");
        let (traj, w) = read_pdb(&write_tmp("test_cut.pdb", &cut)).unwrap();
        assert_eq!(traj.n_frames(), 1);
        assert!(w.len() == 1 && w[0].line == Some(6) && w[0].message.contains("frame 1"), "{w:?}");
        let only = &WATER_PDB[..WATER_PDB.find("-0.763").unwrap() + 3]; // cut inside atom 3's y
        let err = format!("{:#}", read_pdb(&write_tmp("test_cut_only.pdb", only)).unwrap_err());
        assert!(err.contains("only frame"), "{err}");
    }

    #[test]
    fn test_non_ascii_does_not_panic() {
        let text = WATER_PDB.replace("Water molecule", "水分子测试");
        let traj = read_pdb(&write_tmp("test_utf8.pdb", &text)).unwrap().0;
        assert_eq!(traj.metadata.title.as_deref(), Some("水分子测试"));
        // 坐标列里混进多字节字符：报错，不 panic
        let bad = WATER_PDB.replace("0.000   0.763", "0.000  水763");
        assert!(read_pdb(&write_tmp("test_utf8_bad.pdb", &bad)).is_err());
    }

    #[test]
    fn test_multi_model() {
        let path = write_tmp("test_multi.pdb", MULTI_PDB);
        let traj = read_pdb(&path).unwrap().0;
        assert_eq!(traj.n_frames(), 2);
    }

    const C1: &str = "ATOM      1  C   UNK A   1       0.000   0.000   0.000  1.00  0.00           C\n";
    const C2: &str = "ATOM      2  C   UNK A   1       1.400   0.000   0.000  1.00  0.00           C\n";

    /// A last MODEL with no ENDMDL used to be dropped without a word.
    #[test]
    fn test_last_model_without_endmdl() {
        // Same atom count as the frame before: kept, with a warning
        let text = format!("MODEL 1\n{C1}{C2}ENDMDL\nMODEL 2\n{C1}{C2}");
        let (traj, w) = read_pdb(&write_tmp("test_noend.pdb", &text)).unwrap();
        assert_eq!(traj.n_frames(), 2);
        assert!(w.len() == 1 && w[0].line == Some(5) && w[0].message.contains("frame 1")
            && w[0].message.contains("kept"), "{w:?}");
        // Fewer atoms: cut off, dropped with a warning
        let text = format!("MODEL 1\n{C1}{C2}ENDMDL\nMODEL 2\n{C1}");
        let (traj, w) = read_pdb(&write_tmp("test_noend_cut.pdb", &text)).unwrap();
        assert_eq!(traj.n_frames(), 1);
        assert!(w.len() == 1 && w[0].message.contains("dropped"), "{w:?}");
        // A single MODEL with no ENDMDL: kept
        let (traj, _) = read_pdb(&write_tmp("test_noend_one.pdb", &format!("MODEL 1\n{C1}"))).unwrap();
        assert_eq!(traj.n_frames(), 1);
    }

    /// CP2K / VMD trajectories close each block with END and have no MODEL records;
    /// they used to run together into one frame of every block's atoms.
    #[test]
    fn test_end_separated_blocks_are_frames() {
        let text = format!("{C1}{C2}END\n{C1}{C2}END\n");
        let (traj, w) = read_pdb(&write_tmp("test_end_blocks.pdb", &text)).unwrap();
        assert_eq!(traj.n_frames(), 2);
        assert_eq!(traj.frames[1].n_atoms(), 2);
        assert!(w.is_empty(), "{w:?}");
        // MODEL … ENDMDL followed by END gives no extra frame
        let (traj, _) = read_pdb(&write_tmp("test_model_end.pdb", MULTI_PDB_END)).unwrap();
        assert_eq!(traj.n_frames(), 2);
    }

    /// The last END-split block cut short used to be kept as a frame of fewer atoms,
    /// with no warning (scan B-1).
    #[test]
    fn test_last_end_block_without_end() {
        let text = format!("{C1}{C2}END\n{C1}{C2}END\n{C1}");
        let (traj, w) = read_pdb(&write_tmp("test_end_cut.pdb", &text)).unwrap();
        assert_eq!(traj.n_frames(), 2);
        assert!(w.len() == 1 && w[0].message.contains("frame 2") && w[0].message.contains("dropped"), "{w:?}");
        // Same atom count as the frame before: kept, with a warning
        let text = format!("{C1}{C2}END\n{C1}{C2}");
        let (traj, w) = read_pdb(&write_tmp("test_end_noend.pdb", &text)).unwrap();
        assert_eq!(traj.n_frames(), 2);
        assert!(w.len() == 1 && w[0].message.contains("frame 1") && w[0].message.contains("kept"), "{w:?}");
    }

    const MULTI_PDB_END: &str = "\
MODEL        1
ATOM      1  C   UNK A   1       0.000   0.000   0.000  1.00  0.00           C
ENDMDL
MODEL        2
ATOM      1  C   UNK A   1       1.400   0.000   0.000  1.00  0.00           C
ENDMDL
END
";

    /// `CRYST1 1 1 1 90 90 90` is the PDB placeholder for "no unit cell" (as MDAnalysis
    /// reads it); it used to become a 1 Å periodic box.
    #[test]
    fn test_cryst1_placeholder_is_no_cell() {
        let cryst = |abc: &str| format!("CRYST1{abc}  90.00  90.00  90.00 P 1           1\n{C1}END\n");
        let (traj, w) = read_pdb(&write_tmp("test_cryst_1.pdb",
            &cryst("    1.000    1.000    1.000"))).unwrap();
        assert!(traj.frames[0].cell.is_none());
        assert_eq!(traj.frames[0].pbc, [false; 3]);
        assert!(w.len() == 1 && w[0].message.contains("placeholder"), "{w:?}");
        let (traj, w) = read_pdb(&write_tmp("test_cryst_10.pdb",
            &cryst("   10.000   10.000   10.000"))).unwrap();
        assert!(traj.frames[0].cell.is_some());
        assert!(w.is_empty());
    }

    /// A CRYST1 that gives no cell used to read as non-periodic, silently (scan B-6).
    #[test]
    fn test_bad_cryst1_is_an_error_naming_the_line() {
        let text = format!("HEADER    x\nCRYST1   1o.000   10.000   10.000  90.00  90.00  90.00 P 1           1\n{C1}END\n");
        let err = format!("{:#}", read_pdb(&write_tmp("test_cryst_bad.pdb", &text)).unwrap_err());
        assert!(err.contains("line 2") && err.contains("CRYST1 a in columns 7-15"), "{err}");
        // Cut before gamma
        let text = format!("CRYST1   10.000   10.000   10.000  90.00  90.00\n{C1}END\n");
        let err = format!("{:#}", read_pdb(&write_tmp("test_cryst_short.pdb", &text)).unwrap_err());
        assert!(err.contains("gamma"), "{err}");
    }
}
