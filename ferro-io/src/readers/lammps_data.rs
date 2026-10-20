use std::path::Path;
use crate::Notice;
use std::collections::HashMap;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use ferro_core::data::elements::by_symbol;
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
    Ok(read_lammps_data_with_warnings(path, style)?.0)
}

/// [`read_lammps_data`], returning what the reader noticed: a Masses comment that is
/// not an element symbol, for which the element came from the mass instead.
pub fn read_lammps_data_with_warnings(path: &Path, style: AtomStyle) -> Result<(Trajectory, Vec<Notice>)> {
    let path_ = path.display();
    let content = super::util::read_text(path)?;
    parse_lammps_data(&content, style).with_context(|| format!("parsing {path_} as atom style {style}"))
}

fn parse_lammps_data(content: &str, style: AtomStyle) -> Result<(Trajectory, Vec<Notice>)> {
    let lines: Vec<&str> = content.lines().collect();
    let mut i = 0;
    let mut warnings = Vec::new();

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
        } else if is_section_header(line) {
            break;
        }
        i += 1;
    }

    // Build cell
    let lx = xhi - xlo;
    let ly = yhi - ylo;
    let lz = zhi - zlo;
    let mut cell = if is_triclinic {
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
    // data 的 lo 就是真实原点（不同于 dump 三斜的 *_bound），writer 写回（审查 A-12）
    cell.origin = nalgebra::Vector3::new(xlo, ylo, zlo);

    // LAMMPS 的段可任意排序（read_data 文档），按段名定位；段内容到下一个段头为止。
    // 段头后的 `# style` 注释被 strip_comment 去掉，布局只由调用者给的 style 决定
    let section = |name: &str| (i..lines.len()).find(|&k| strip_comment(lines[k]) == name);

    // ── Masses section ────────────────────────────────────────────────────────
    let mut type_element: HashMap<usize, String> = HashMap::new();
    let mut type_mass: HashMap<usize, f64> = HashMap::new();
    if let Some(start) = section("Masses") {
        for (k, raw) in lines.iter().enumerate().skip(start + 1) {
            let line = strip_comment(raw);
            if line.is_empty() { continue; }
            if is_section_header(line) { break; }
            // 恰好 `type mass` 两项：多出的列多半是认不出段头的另一段混了进来（审查 A-10）
            let parts: Vec<&str> = line.split_whitespace().collect();
            let entry = match parts[..] {
                [t, m] => t.parse::<usize>().ok().zip(m.parse::<f64>().ok()),
                _ => None,
            };
            let Some((tid, mass)) = entry else {
                bail!("line {}: Masses entry '{line}' is not 'type mass'", k + 1);
            };
            // An element symbol in the comment wins: "# Fe" or "# Fe ...". A comment that is not one
            // (`# oxygen`, `# Ow`) used to become the element name; the mass decides instead, with a warning
            let comment = raw.find('#')
                .and_then(|p| raw[p+1..].split_whitespace().next());
            let from_comment = comment
                .and_then(|c| by_symbol(&super::util::capitalize_symbol(c)).map(|e| e.symbol.to_string()));
            let elem = match (from_comment, comment) {
                (Some(e), _) => e,
                (None, c) => {
                    let e = element_from_mass(mass).to_string();
                    if let Some(c) = c {
                        warnings.push(Notice::warning(Some(k + 1), format!(
                            "Masses comment {c:?} is not an element symbol; \
                             atom type {tid} taken as {e} from its mass {mass}")));
                    }
                    e
                }
            };
            type_mass.insert(tid, mass);
            if type_element.insert(tid, elem).is_some() {
                bail!("line {}: Masses gives atom type {tid} twice", k + 1);
            }
        }
    }

    // ── Atoms section ─────────────────────────────────────────────────────────
    i = section("Atoms").context("no Atoms section")? + 1;

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
        atom.mass = type_mass.get(&tp).copied();
        // A charge / full charge column holds values, and 0 is one; it used to become None, the same
        // as having no charge column
        atom.charge = q;
        frame.add_atom(atom);
    }

    let mut traj = Trajectory::from_frame(frame);
    if !comment.is_empty() { traj.metadata.title = Some(comment); }
    Ok((traj, warnings))
}

fn strip_comment(l: &str) -> &str {
    l.split('#').next().unwrap_or("").trim()
}

/// `line` 须已经过 `strip_comment`。名单取自 LAMMPS `read_data` 文档的全部段名
/// （审查 A-10：只认十个时，`PairIJ Coeffs` 等段的行被当成上一段的内容）
fn is_section_header(line: &str) -> bool {
    matches!(line,
        "Atoms" | "Velocities" | "Masses" | "Ellipsoids" | "Lines" | "Triangles" | "Bodies"
        | "Bonds" | "Angles" | "Dihedrals" | "Impropers"
        | "Pair Coeffs" | "PairIJ Coeffs" | "Bond Coeffs" | "Angle Coeffs" | "Dihedral Coeffs"
        | "Improper Coeffs" | "BondBond Coeffs" | "BondAngle Coeffs" | "MiddleBondTorsion Coeffs"
        | "EndBondTorsion Coeffs" | "AngleTorsion Coeffs" | "AngleAngleTorsion Coeffs"
        | "BondBond13 Coeffs" | "AngleAngle Coeffs"
        | "Atom Type Labels" | "Bond Type Labels" | "Angle Type Labels" | "Dihedral Type Labels"
        | "Improper Type Labels")
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

    // The Masses value is the atom's mass, not just a hint for guessing the element:
    // D2O must keep 2.014 for H, or every density and centre of mass is silently off.
    #[test]
    fn test_masses_are_kept_as_atom_masses() {
        let text = WATER_FULL.replace("2 1.008   # H", "2 2.014   # H");
        let path = tmp("d2o.data", &text);
        let traj = read_lammps_data(&path, AtomStyle::Full).unwrap();
        let f = traj.first().unwrap();
        assert_eq!(f.atom(1).element, "H");
        assert_eq!(f.atom(1).mass, Some(2.014));
        assert_eq!(f.atom(0).mass, Some(15.999));

        let out = path.with_file_name("d2o_roundtrip.data");
        crate::writers::lammps_data::write_lammps_data(&traj, &out, false).unwrap();
        let back = read_lammps_data(&out, AtomStyle::Full).unwrap();
        assert!((back.first().unwrap().atom(1).effective_mass().unwrap() - 2.014).abs() < 1e-9);

        // H and D in one frame share an element but not a mass: refuse, write nothing
        let mut mixed = traj.clone();
        mixed.frames[0].sites_mut()[2].mass = Some(1.008);
        let refused = path.with_file_name("d2o_mixed.data");
        let err = crate::writers::lammps_data::write_lammps_data(&mixed, &refused, false).unwrap_err();
        assert!(err.to_string().contains("different masses"), "{err}");
        assert!(!refused.exists(), "a refused write must not leave a file");
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
            let elems: Vec<&str> = f.atoms().map(|a| a.element).collect();
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

    // 审查 A-10：Masses 后接未登记的段时，其行被当成 Masses，O/H 两 type 都成 X 并合并
    #[test]
    fn test_sections_after_masses_are_not_masses() {
        let pairij = "\nPairIJ Coeffs # lj/cut\n\n1 1 0.15 3.15\n1 2 0.0 1.0\n2 2 0.0 1.0\n";
        let text = WATER_FULL.replace("2 1.008   # H\n", &format!("2 1.008   # H\n{pairij}"));
        let traj = read_lammps_data(&tmp("pairij.data", &text), AtomStyle::Full).unwrap();
        let f = traj.first().unwrap();
        assert_eq!((f.atom(0).element, f.atom(1).element), ("O", "H"));
    }

    // LAMMPS 的段可任意排序：Masses 写在 Atoms 之后也得读到
    #[test]
    fn test_masses_after_atoms() {
        let masses = "Masses\n\n1 15.999  # O\n2 1.008   # H\n\n";
        let text = WATER_FULL.replace(masses, "") + "\n" + masses;
        let traj = read_lammps_data(&tmp("late_masses.data", &text), AtomStyle::Full).unwrap();
        let f = traj.first().unwrap();
        assert_eq!((f.atom(0).element, f.atom(1).element), ("O", "H"));
    }

    #[test]
    fn test_bad_masses_lines_are_errors() {
        let err = |name: &str, text: &str| {
            format!("{:#}", read_lammps_data(&tmp(name, text), AtomStyle::Full).unwrap_err())
        };
        // 一行三个数：多半是认不出的段混了进来
        let three = WATER_FULL.replace("2 1.008   # H", "2 1.008 3.0");
        assert!(err("m1.data", &three).contains("Masses"), "三个数应报错");
        // 同一 type 写两次
        let dup = WATER_FULL.replace("2 1.008   # H", "1 1.008   # H");
        assert!(err("m2.data", &dup).contains("twice"), "重复 type 应报错");
        // mass 不是数
        let nan = WATER_FULL.replace("2 1.008   # H", "2 abc   # H");
        assert!(err("m3.data", &nan).contains("Masses"), "坏 mass 应报错");
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

    fn data_with(masses: &str, atoms: &str, style: &str) -> String {
        format!("t\n\n2 atoms\n2 atom types\n\n0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\n\
                 Masses\n\n{masses}\nAtoms # {style}\n\n{atoms}")
    }

    /// A comment that is not an element symbol used to become the element name.
    #[test]
    fn test_masses_comment_must_be_an_element() {
        let text = data_with("1 15.999  # oxygen\n2 1.008  # h\n",
                             "1 1 0 0 0\n2 2 1 0 0\n", "atomic");
        let (traj, w) = parse_lammps_data(&text, AtomStyle::Atomic).unwrap();
        let f = &traj.frames[0];
        assert_eq!(f.atom(0).element, "O", "falls back to the mass");
        assert_eq!(f.atom(1).element, "H", "case is normalised");
        assert!(w.len() == 1 && w[0].message.contains("\"oxygen\"") && w[0].message.contains("atom type 1"), "{w:?}");
    }

    /// A charge of 0 is a value: it used to read as no charge at all.
    #[test]
    fn test_zero_charge_is_kept() {
        let text = data_with("1 15.999  # O\n2 1.008  # H\n",
                             "1 1 0.0 0 0 0\n2 2 0.4 1 0 0\n", "charge");
        let f = parse_lammps_data(&text, AtomStyle::Charge).unwrap().0.frames.remove(0);
        assert_eq!(f.atom(0).charge, Some(0.0));
        assert_eq!(f.atom(1).charge, Some(0.4));
        // atomic has no charge column: still None
        let text = data_with("1 15.999  # O\n2 1.008  # H\n", "1 1 0 0 0\n2 2 1 0 0\n", "atomic");
        assert_eq!(parse_lammps_data(&text, AtomStyle::Atomic).unwrap().0.frames[0].atom(0).charge, None);
    }

    /// LAMMPS does not order the sections; Atoms before Masses still gets the elements.
    #[test]
    fn test_atoms_before_masses() {
        let text = "t\n\n2 atoms\n2 atom types\n\n0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\n\
                    Atoms # atomic\n\n1 1 0 0 0\n2 2 1 0 0\n\nMasses\n\n1 15.999  # O\n2 1.008  # H\n";
        let f = parse_lammps_data(text, AtomStyle::Atomic).unwrap().0.frames.remove(0);
        assert_eq!((f.atom(0).element, f.atom(1).element), ("O", "H"));
    }
}
