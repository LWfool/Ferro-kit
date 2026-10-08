use std::path::Path;
use std::collections::HashMap;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{bail, Context, Result};
use super::util::floats;

/// LAMMPS `atom_style`，决定 `Atoms` 段每行的列布局。
///
/// 由调用者显式给出：`Atoms # full` 这类注释在 LAMMPS 里是可选的，列数也推不出
/// style（charge 与 molecular 都是 6 列），两条路都会静默读错（审查 M1，用户裁定
/// 2026-10-02：不看注释、不猜）。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AtomStyle { Atomic, Charge, Full }

impl AtomStyle {
    /// 不含 image flag 的列数；带 image flag 时再多 3 列。
    fn columns(self) -> usize {
        match self { AtomStyle::Atomic => 5, AtomStyle::Charge => 6, AtomStyle::Full => 7 }
    }
}

impl std::str::FromStr for AtomStyle {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, String> {
        match s {
            "atomic" => Ok(AtomStyle::Atomic),
            "charge" => Ok(AtomStyle::Charge),
            "full" => Ok(AtomStyle::Full),
            _ => Err(format!("unknown atom style '{s}' (expected atomic, charge or full)")),
        }
    }
}

impl std::fmt::Display for AtomStyle {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            AtomStyle::Atomic => "atomic",
            AtomStyle::Charge => "charge",
            AtomStyle::Full => "full",
        })
    }
}

pub fn read_lammps_data(path: &Path, style: AtomStyle) -> Result<Trajectory> {
    let path_ = path.display();
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path_}"))?;
    parse_lammps_data(&content, style).with_context(|| format!("parsing {path_} as atom style {style}"))
}

fn parse_lammps_data(content: &str, style: AtomStyle) -> Result<Trajectory> {
    let lines: Vec<&str> = content.lines().collect();
    let mut i = 0;

    let skip = |l: &str| -> bool { l.trim().is_empty() || l.trim().starts_with('#') };

    // ── Header ────────────────────────────────────────────────────────────────
    // First non-blank line is the comment
    while i < lines.len() && skip(lines[i]) { i += 1; }
    // 空文件 / 只有空行时以前 `lines[i]` 越界 panic
    let comment = lines.get(i).context("empty file: no header comment line")?.trim().to_string();
    i += 1;

    // Box bounds
    let mut xlo = 0.0_f64; let mut xhi = 0.0_f64;
    let mut ylo = 0.0_f64; let mut yhi = 0.0_f64;
    let mut zlo = 0.0_f64; let mut zhi = 0.0_f64;
    let mut xy = 0.0_f64; let mut xz = 0.0_f64; let mut yz = 0.0_f64;
    let mut is_triclinic = false;
    let mut n_atoms: Option<usize> = None;

    while i < lines.len() {
        let line = strip_comment(lines[i]);
        if line.is_empty() { i += 1; continue; }

        if line.ends_with("atoms") && !line.contains("atom types") {
            let n = line.split_whitespace().next().unwrap_or("");
            n_atoms = Some(n.parse().with_context(|| format!("bad atom count in '{line}'"))?);
        } else if line.ends_with("atom types") {
            // type count stored but not needed
        } else if line.contains("xlo xhi") {
            let v = floats(line, 2)?;
            xlo = v[0]; xhi = v[1];
        } else if line.contains("ylo yhi") {
            let v = floats(line, 2)?;
            ylo = v[0]; yhi = v[1];
        } else if line.contains("zlo zhi") {
            let v = floats(line, 2)?;
            zlo = v[0]; zhi = v[1];
        } else if line.contains("xy xz yz") {
            let v = floats(line, 3)?;
            xy = v[0]; xz = v[1]; yz = v[2];
            is_triclinic = true;
        } else if line == "Masses" || line == "Atoms" || line == "Atoms # full"
                || line == "Atoms # atomic" || line == "Atoms # charge"
                || line.starts_with("Masses") || line.starts_with("Atoms") {
            break;
        }
        i += 1;
    }

    // Build cell
    let lx = xhi - xlo;
    let ly = yhi - ylo;
    let lz = zhi - zlo;
    let cell = if is_triclinic {
        Cell::from_matrix(Matrix3::new(
            lx,  0.0, 0.0,
            xy,  ly,  0.0,
            xz,  yz,  lz,
        ))
    } else {
        Cell::from_matrix(Matrix3::new(
            lx,  0.0, 0.0,
            0.0, ly,  0.0,
            0.0, 0.0, lz,
        ))
    };

    // ── Masses section ────────────────────────────────────────────────────────
    let mut type_element: HashMap<usize, String> = HashMap::new();
    while i < lines.len() {
        let raw = lines[i];
        let line = strip_comment(raw).trim();
        if line == "Masses" { i += 1; continue; }
        if line.is_empty() { i += 1; continue; }

        // Check if this line is a section header (Atoms, Velocities, ...)
        if is_section_header(line) { break; }

        // Masses line: type_id mass  # element_comment
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            if let Ok(tid) = parts[0].parse::<usize>() {
                let mass: f64 = parts[1].parse().unwrap_or(1.0);
                // Try comment for element symbol: "# Fe" or "# Fe ..."
                let elem = raw.find('#')
                    .and_then(|p| {
                        raw[p+1..].split_whitespace().next().map(|s| s.to_string())
                    })
                    .unwrap_or_else(|| element_from_mass(mass).to_string());
                type_element.insert(tid, elem);
            }
        }
        i += 1;
    }

    // ── Atoms section ─────────────────────────────────────────────────────────
    // 段头后的 `# style` 注释不看，布局只由调用者给的 style 决定
    while i < lines.len() {
        let line = strip_comment(lines[i]);
        i += 1;
        if line.starts_with("Atoms") { break; }
    }

    let Some(n_atoms) = n_atoms else { bail!("header has no 'N atoms' line") };
    let base = style.columns();
    // (id, type, charge, 坐标)；image flag 已展开进坐标
    let mut atoms_raw: Vec<(usize, usize, Option<f64>, Vector3<f64>)> = Vec::new();
    // 行优先：行 = a、b、c 三条晶格矢量，image flag (ix,iy,iz) 平移 ix·a+iy·b+iz·c
    let (a, b, c) = (cell.matrix.row(0).transpose(), cell.matrix.row(1).transpose(),
                     cell.matrix.row(2).transpose());

    while i < lines.len() {
        let line = strip_comment(lines[i]);
        let lineno = i + 1;
        i += 1;
        if line.is_empty() { continue; }
        if is_section_header(line) { break; }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() != base && parts.len() != base + 3 {
            bail!("line {lineno}: {} columns, atom style {style} has {base} (or {} with image flags)",
                parts.len(), base + 3);
        }
        let int = |k: usize| parts[k].parse::<usize>()
            .with_context(|| format!("line {lineno}: column {} '{}' is not an integer", k + 1, parts[k]));
        let real = |k: usize| parts[k].parse::<f64>()
            .with_context(|| format!("line {lineno}: column {} '{}' is not a number", k + 1, parts[k]));
        // 各 style 的 type、charge、x 所在列
        let (tp_at, q_at, x_at) = match style {
            AtomStyle::Atomic => (1, None, 2),     // id type x y z
            AtomStyle::Charge => (1, Some(2), 3),  // id type q x y z
            AtomStyle::Full => (2, Some(3), 4),    // id mol type q x y z
        };
        let mut pos = Vector3::new(real(x_at)?, real(x_at + 1)?, real(x_at + 2)?);
        if parts.len() == base + 3 {
            let flag = |k: usize| parts[k].parse::<i64>()
                .with_context(|| format!("line {lineno}: image flag '{}' is not an integer", parts[k]));
            pos += a * flag(base)? as f64 + b * flag(base + 1)? as f64 + c * flag(base + 2)? as f64;
        }
        let q = q_at.map(real).transpose()?;
        atoms_raw.push((int(0)?, int(tp_at)?, q, pos));
    }
    if atoms_raw.len() != n_atoms {
        bail!("header declares {n_atoms} atoms, the Atoms section has {}", atoms_raw.len());
    }

    // Sort by atom id
    atoms_raw.sort_by_key(|(id, ..)| *id);

    let mut frame = Frame::with_cell(cell, [true; 3]);
    for (_, tp, q, pos) in atoms_raw {
        let elem = type_element.get(&tp)
            .cloned()
            .unwrap_or_else(|| format!("X{tp}"));
        let mut atom = Atom::new(elem, pos);
        if let Some(q) = q.filter(|&q| q != 0.0) { atom.charge = Some(q); }
        frame.add_atom(atom);
    }

    let mut traj = Trajectory::from_frame(frame);
    if !comment.is_empty() { traj.metadata.source = Some(comment); }
    Ok(traj)
}

fn strip_comment(l: &str) -> &str {
    l.split('#').next().unwrap_or("").trim()
}

fn is_section_header(line: &str) -> bool {
    matches!(line, "Masses" | "Atoms" | "Velocities" | "Bonds" | "Angles"
        | "Dihedrals" | "Impropers" | "Pair Coeffs" | "Bond Coeffs"
        | "Angle Coeffs")
    || line.starts_with("Atoms #")
}


fn element_from_mass(mass: f64) -> &'static str {
    ferro_core::data::elements::ELEMENTS
        .iter()
        .filter(|e| (e.atomic_mass - mass).abs() < 0.5)
        .min_by(|a, b| {
            (a.atomic_mass - mass).abs()
                .partial_cmp(&(b.atomic_mass - mass).abs())
                .unwrap()
        })
        .map(|e| e.symbol)
        .unwrap_or("X")
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER_FULL: &str = "LAMMPS data file

6 atoms
2 atom types

0.0 20.0 xlo xhi
0.0 20.0 ylo yhi
0.0 20.0 zlo zhi

Masses

1 15.999  # O
2 1.008   # H

Atoms # full

1 1 1 -0.834 0.000  0.000  0.000
2 1 2  0.417 0.758  0.587  0.000
3 1 2  0.417 -0.758 0.587  0.000
4 2 1 -0.834 10.000 0.000  0.000
5 2 2  0.417 10.758 0.587  0.000
6 2 2  0.417 9.242  0.587  0.000
";

    use crate::testutil::write_tmp as tmp;

    #[test]
    fn test_water_full() {
        let traj = read_lammps_data(&tmp("water.lammps", WATER_FULL), AtomStyle::Full).unwrap();
        let frame = traj.first().unwrap();
        assert_eq!(frame.n_atoms(), 6);
        assert_eq!(frame.atom(0).element, "O");
        assert_eq!(frame.atom(1).element, "H");
        assert!((frame.atom(0).charge.unwrap() - (-0.834)).abs() < 1e-6);
        let [a, ..] = frame.cell.as_ref().unwrap().lengths();
        assert!((a - 20.0).abs() < 1e-6);
    }

    // 以下三份由 ASE 3.29 `write(..., atom_style=st, write_image_flags=True)` 生成：
    // 三斜胞，第 2、3 原子在胞外，写出时折回并记 image flag
    const ASE_ATOMIC: &str = "(written by ASE)

3 atoms
2 atom types

0.0                       6  xlo xhi
0.0                       5  ylo yhi
0.0                       7  zlo zhi
                    1.5                      -1     0.80000000000000004  xy xz yz

Masses

1      1.0079999997406976 # H
2      15.998999995884349 # O

Atoms # atomic

     1   2     0.49999999999999989                     0.5                     0.5      0      0      0
     2   1     0.20000000000000029                     1.8      6.3999999999999995      1      0     -1
     3   1      4.1999999999999993      0.3000000000000001      1.3999999999999997     -1      1      1
";

    const ASE_CHARGE: &str = "(written by ASE)

3 atoms
2 atom types

0.0                       6  xlo xhi
0.0                       5  ylo yhi
0.0                       7  zlo zhi
                    1.5                      -1     0.80000000000000004  xy xz yz

Masses

1      1.0079999997406976 # H
2      15.998999995884349 # O

Atoms # charge

     1   2  -0.8     0.49999999999999989                     0.5                     0.5      0      0      0
     2   1   0.4     0.20000000000000029                     1.8      6.3999999999999995      1      0     -1
     3   1   0.4      4.1999999999999993      0.3000000000000001      1.3999999999999997     -1      1      1
";

    const ASE_FULL: &str = "(written by ASE)

3 atoms
2 atom types

0.0                       6  xlo xhi
0.0                       5  ylo yhi
0.0                       7  zlo zhi
                    1.5                      -1     0.80000000000000004  xy xz yz

Masses

1      1.0079999997406976 # H
2      15.998999995884349 # O

Atoms # full

     1   0   2  -0.8     0.49999999999999989                     0.5                     0.5      0      0      0
     2   0   1   0.4     0.20000000000000029                     1.8      6.3999999999999995      1      0     -1
     3   0   1   0.4      4.1999999999999993      0.3000000000000001      1.3999999999999997     -1      1      1
";

    #[test]
    fn test_each_style_matches_ase_with_image_flags_unwrapped() {
        // 期望值 = ASE read(..., atom_style=st) 的结果，即折回前的原始坐标
        let want = [[0.5, 0.5, 0.5], [7.2, 1.0, -0.6], [-1.3, 6.1, 8.4]];
        for (style, text, q) in [
            (AtomStyle::Atomic, ASE_ATOMIC, [None, None, None]),
            (AtomStyle::Charge, ASE_CHARGE, [Some(-0.8), Some(0.4), Some(0.4)]),
            (AtomStyle::Full, ASE_FULL, [Some(-0.8), Some(0.4), Some(0.4)]),
        ] {
            let traj = read_lammps_data(&tmp(&format!("ase_{style}.data"), text), style).unwrap();
            let f = traj.first().unwrap();
            let elems: Vec<&str> = f.atoms.iter().map(|a| a.element.as_str()).collect();
            assert_eq!(elems, ["O", "H", "H"], "{style}");
            for (k, w) in want.iter().enumerate() {
                let p = f.atom(k).position;
                for d in 0..3 {
                    assert!((p[d] - w[d]).abs() < 1e-9, "{style} 原子 {k}：{p:?}，ASE 给 {w:?}");
                }
                assert_eq!(f.atom(k).charge, q[k], "{style} 原子 {k} 电荷");
            }
        }
    }

    #[test]
    fn test_empty_file_is_an_error() {
        for text in ["", "\n  \n# 只有注释\n"] {
            let e = parse_lammps_data(text, AtomStyle::Atomic).unwrap_err();
            assert!(format!("{e:#}").contains("empty file"), "{e:#}");
        }
    }

    #[test]
    fn test_wrong_style_or_broken_lines_are_errors() {
        let err = |name: &str, text: &str, style| {
            format!("{:#}", read_lammps_data(&tmp(name, text), style).unwrap_err())
        };
        // style 与列数对不上：full 文件按 atomic、atomic 文件按 full
        assert!(err("w1.data", ASE_FULL, AtomStyle::Atomic).contains("10 columns"));
        assert!(err("w2.data", ASE_ATOMIC, AtomStyle::Full).contains("8 columns"));
        // charge（9 列）按 full 读：既不是 7 也不是 10
        assert!(err("w3.data", ASE_CHARGE, AtomStyle::Full).contains("9 columns"));
        // 少一行原子
        let short = WATER_FULL.replace("6 2 2  0.417 9.242  0.587  0.000\n", "");
        assert!(err("w4.data", &short, AtomStyle::Full).contains("declares 6 atoms"));
        // 字段坏了
        let bad = WATER_FULL.replace("10.758", "10.7x8");
        assert!(err("w5.data", &bad, AtomStyle::Full).contains("'10.7x8' is not a number"));
        let flag = ASE_ATOMIC.replacen("      1      0     -1", "      1      0     -1.5", 1);
        assert!(err("w6.data", &flag, AtomStyle::Atomic).contains("image flag '-1.5'"));
        // 没有 N atoms
        let no_n = WATER_FULL.replace("6 atoms\n", "");
        assert!(err("w7.data", &no_n, AtomStyle::Full).contains("no 'N atoms'"));
    }

    #[test]
    fn test_atoms_comment_is_ignored() {
        // 注释写 atomic，实际按调用者给的 full 读
        let text = WATER_FULL.replace("Atoms # full", "Atoms # atomic");
        let traj = read_lammps_data(&tmp("comment.data", &text), AtomStyle::Full).unwrap();
        assert_eq!(traj.first().unwrap().n_atoms(), 6);
    }

    #[test]
    fn test_style_names() {
        for (name, style) in [("atomic", AtomStyle::Atomic), ("charge", AtomStyle::Charge), ("full", AtomStyle::Full)] {
            assert_eq!(name.parse::<AtomStyle>(), Ok(style));
            assert_eq!(style.to_string(), name);
        }
        // 只收这三种：molecular 与 charge 同为 6 列，收进来就得靠调用者分清
        assert!("molecular".parse::<AtomStyle>().is_err());
    }
}
