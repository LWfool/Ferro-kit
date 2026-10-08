use std::path::Path;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::Vector3;
use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::{bail, Context, Result};
use ferro_core::data::elements::by_symbol;

/// 读取 PDB 文件，支持多模型（MODEL/ENDMDL 记录）。
pub fn read_pdb(path: &Path) -> Result<Trajectory> {
    let path_ = path.display();
    let file = File::open(path).context(format!("cannot open {path_}"))?;
    let reader = BufReader::new(file);

    let mut traj = Trajectory::new();
    let mut current = Frame::new();
    let mut has_model = false;
    let mut cell: Option<Cell> = None;

    for (ln, line) in reader.lines().enumerate() {
        let line = line.context("read error")?;
        let tag = col(&line, 0, 6);

        match tag {
            "HEADER" => {
                // 记录名 6 字节是 ASCII，6 处必是字符边界；以前从 10 切，多字节标题会 panic
                let source = line[6..].trim();
                if !source.is_empty() {
                    traj.metadata.source = Some(source.to_string());
                }
            }
            "CRYST1" => {
                cell = parse_cryst1(&line);
            }
            "MODEL " => {
                has_model = true;
            }
            "ENDMDL" => {
                let mut frame = std::mem::replace(&mut current, Frame::new());
                if let Some(c) = &cell {
                    frame.cell = Some(c.clone());
                    frame.pbc = [true; 3];
                }
                traj.add_frame(frame);
            }
            "ATOM  " | "HETATM" => {
                let atom = parse_atom_record(&line)
                    .with_context(|| format!("line {}: {line:?}", ln + 1))?;
                current.add_atom(atom);
            }
            _ => {}
        }
    }

    // 没有 MODEL 记录时，视为单帧
    if !has_model && current.n_atoms() > 0 {
        if let Some(c) = cell {
            current.cell = Some(c);
            current.pbc = [true; 3];
        }
        traj.add_frame(current);
    }

    Ok(traj)
}

/// 第 a..b 字节列（0 起、左闭右开），行短或落在多字节字符中间时为空串。PDB 是按列的
/// ASCII 格式，非 ASCII 只会出现在不该出现的地方，空串让各字段自己报错或回落
fn col(line: &str, a: usize, b: usize) -> &str {
    line.get(a..b.min(line.len())).unwrap_or("")
}

/// CRYST1 format: CRYST1  a  b  c  α  β  γ  sgroup  z
fn parse_cryst1(line: &str) -> Option<Cell> {
    if line.len() < 54 { return None; }
    let a:  f64 = col(line, 6, 15).trim().parse().ok()?;
    let b:  f64 = col(line, 15, 24).trim().parse().ok()?;
    let c:  f64 = col(line, 24, 33).trim().parse().ok()?;
    let al: f64 = col(line, 33, 40).trim().parse().ok()?;
    let be: f64 = col(line, 40, 47).trim().parse().ok()?;
    let ga: f64 = col(line, 47, 54).trim().parse().ok()?;
    Cell::from_lengths_angles(a, b, c, al, be, ga).ok()
}

/// ATOM / HETATM 一行 → 原子。坐标坏或行短于 54 列报错（同 ASE）—— 以前静默丢原子，
/// 原子数少了却退出码 0
fn parse_atom_record(line: &str) -> Result<Atom> {
    let xyz: Option<Vec<f64>> = [(30, 38), (38, 46), (46, 54)].iter()
        .map(|&(a, b)| col(line, a, b).trim().parse().ok())
        .collect();
    let Some(&[x, y, z]) = xyz.as_deref() else {
        bail!("invalid or missing coordinates in columns 31-54");
    };
    // 第 77–78 列是规范的元素位，RCSB 写大写（`FE`）—— 规范成 `Fe`，不另校验
    // （元素表只到 Rn，校验会误拒 U 之类）。空着时从原子名推
    let sym = col(line, 76, 78).trim();
    let element = if sym.is_empty() { element_from_atom_name(col(line, 12, 16))? } else { capitalize(sym) };
    Ok(Atom::new(element, Vector3::new(x, y, z)))
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
        .map(|c| capitalize(&c))
        .find(|c| by_symbol(c).is_some())
        .with_context(|| format!(
            "no element in columns 77-78 and atom name {name:?} does not give one; \
             write the element symbol in columns 77-78"))
}

/// `FE` / `fe` → `Fe`
fn capitalize(s: &str) -> String {
    let mut c = s.chars();
    c.next().map(|f| f.to_ascii_uppercase().to_string() + &c.as_str().to_ascii_lowercase())
        .unwrap_or_default()
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
        let traj = read_pdb(&path).unwrap();
        assert_eq!(traj.n_frames(), 1);
        let frame = traj.first().unwrap();
        assert_eq!(frame.n_atoms(), 3);
        assert_eq!(frame.atom(0).element, "O");
        assert_eq!(traj.metadata.source.as_deref(), Some("Water molecule"));
    }

    /// 原子名 + 77–78 列 → 一行 ATOM（列位照 PDB 规范）
    fn atom_line(name: &str, elem: &str) -> String {
        format!("ATOM      1 {name:<4} UNK A   1       0.000   0.000   0.119  1.00  0.00          {elem:>2}")
    }
    fn element_of(line: &str) -> String {
        read_pdb(&write_tmp("test_elem.pdb", &format!("{line}\n"))).unwrap().frames[0].atoms[0].element.clone()
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
                    format!("{}ATOM      4  H   UNK A   1       0.000\n", WATER_PDB)] {
            let e = read_pdb(&write_tmp("test_bad.pdb", &bad)).unwrap_err();
            let msg = format!("{e:#}");
            assert!(msg.contains("line ") && msg.contains("coordinates"), "应点名行号：{msg}");
        }
        assert!(format!("{:#}", read_pdb(&write_tmp("test_bad.pdb",
            &WATER_PDB.replace("0.763  -0.477", "abc    -0.477"))).unwrap_err()).contains("line 3"));
    }

    #[test]
    fn test_non_ascii_does_not_panic() {
        let text = WATER_PDB.replace("Water molecule", "水分子测试");
        let traj = read_pdb(&write_tmp("test_utf8.pdb", &text)).unwrap();
        assert_eq!(traj.metadata.source.as_deref(), Some("水分子测试"));
        // 坐标列里混进多字节字符：报错，不 panic
        let bad = WATER_PDB.replace("0.000   0.763", "0.000  水763");
        assert!(read_pdb(&write_tmp("test_utf8_bad.pdb", &bad)).is_err());
    }

    #[test]
    fn test_multi_model() {
        let path = write_tmp("test_multi.pdb", MULTI_PDB);
        let traj = read_pdb(&path).unwrap();
        assert_eq!(traj.n_frames(), 2);
    }
}
