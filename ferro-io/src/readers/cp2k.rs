use std::path::Path;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{bail, Context, Result};

pub fn read_cp2k_inp(path: &Path) -> Result<Trajectory> {
    let path_ = path.display();
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path_}"))?;
    parse_cp2k(&content).with_context(|| format!("parsing {path_}"))
}

/// CP2K restart 与 inp 格式相同，共用同一解析器。
pub fn read_cp2k_restart(path: &Path) -> Result<Trajectory> {
    read_cp2k_inp(path)
}

const BOHR: f64 = 0.52917721; // Bohr → Å

fn parse_cp2k(content: &str) -> Result<Trajectory> {
    let lines: Vec<&str> = content.lines().collect();

    // Build a stripped-comment line iterator
    let stripped: Vec<&str> = lines.iter()
        .map(|l| {
            let l = l.trim();
            // CP2K comments start with ! or #
            if let Some(p) = l.find('!').or_else(|| l.find('#')) { &l[..p] } else { l }
        })
        .collect();

    // ── Find &SUBSYS ──────────────────────────────────────────────────────────
    let subsys_start = stripped.iter().position(|l| {
        let low = l.to_lowercase();
        low.starts_with("&subsys")
    }).context("&SUBSYS section not found")?;

    let subsys_end = find_end(&stripped, subsys_start)
        .context("&END SUBSYS not found")?;

    let subsys: &[&str] = &stripped[subsys_start..=subsys_end];

    // ── &KIND label → element map ─────────────────────────────────────────────
    let mut kind_map: std::collections::HashMap<String, String> = Default::default();
    {
        let mut j = 0;
        while j < subsys.len() {
            let low = subsys[j].to_lowercase();
            if low.starts_with("&kind") {
                let label = low.strip_prefix("&kind").unwrap_or("").trim().to_string();
                let end = find_end(subsys, j).unwrap_or(subsys.len() - 1);
                for l in &subsys[j..=end] {
                    let kv: Vec<&str> = l.split_whitespace().collect();
                    if kv.first().map(|s| s.to_lowercase()) == Some("element".to_string()) {
                        if let Some(elem) = kv.get(1) {
                            kind_map.insert(label.clone(), elem.to_string());
                        }
                    }
                }
            }
            j += 1;
        }
    }

    // ── &CELL ─────────────────────────────────────────────────────────────────
    let cell_start = subsys.iter().position(|l| l.to_lowercase().starts_with("&cell"))
        .context("&CELL not found")?;
    let cell_end = find_end(subsys, cell_start).context("&END CELL not found")?;
    let cell_section: &[&str] = &subsys[cell_start..=cell_end];

    let cell = parse_cell_section(cell_section)?;

    // ── &COORD ────────────────────────────────────────────────────────────────
    let coord_start = subsys.iter().position(|l| l.to_lowercase().starts_with("&coord"))
        .context("&COORD not found")?;
    let coord_end = find_end(subsys, coord_start).context("&END COORD not found")?;
    let coord_section: &[&str] = &subsys[coord_start..=coord_end];

    // &COORD 的两个关键字。UNIT 默认 angstrom；SCALED 默认假、单写即真
    // （input_cp2k_subsys.F：default_l_val=.FALSE., lone_keyword_l_val=.TRUE.）
    let mut coord_scale = 1.0;
    let mut scaled = false;
    for l in &coord_section[1..coord_section.len().saturating_sub(1)] {
        let parts: Vec<&str> = l.split_whitespace().collect();
        match parts.first().map(|s| s.to_lowercase()).as_deref() {
            Some("unit") => {
                let u = parts.get(1).context("&COORD UNIT has no value")?;
                coord_scale = length_unit(u)?;
            }
            Some("scaled") => scaled = cp2k_logical(parts.get(1).copied(), true)?,
            _ => {}
        }
    }

    let mut frame = Frame::with_cell(cell.clone(), [true; 3]);

    for l in &coord_section[1..coord_section.len().saturating_sub(1)] {
        let l = l.trim();
        let first = l.split_whitespace().next().map(|s| s.to_lowercase());
        if l.is_empty() || matches!(first.as_deref(), Some("unit" | "scaled")) { continue; }

        // 行格式 `KIND x y z [分子名 ...]`（第 5 列起是分子名，不是速度 ——
        // 速度在独立的 &VELOCITY 段，单位 bohr/au_time，本 reader 不读）
        let parts: Vec<&str> = l.split_whitespace().collect();
        let xyz: Option<Vec<f64>> = parts.get(1..4)
            .map(|p| p.iter().map_while(|s| s.parse().ok()).collect());
        let Some([x, y, z]) = xyz.as_deref().and_then(|v| <[f64; 3]>::try_from(v).ok()) else {
            bail!("invalid &COORD line {l:?}: expected KIND x y z");
        };

        let label = parts[0].to_string();
        let elem = kind_map.get(&label.to_lowercase())
            .cloned()
            .unwrap_or_else(|| extract_element(&label));

        let v = Vector3::new(x, y, z);
        let position = if scaled { cell.fractional_to_cartesian(v) } else { v * coord_scale };

        let mut atom = Atom::new(elem, position);
        atom.label = Some(label);
        frame.add_atom(atom);
    }

    Ok(Trajectory::from_frame(frame))
}

fn parse_cell_section(section: &[&str]) -> Result<Cell> {
    // Supports:
    //   A x y z / B x y z / C x y z  (explicit vectors)
    //   ABC a b c  (diagonal only, Angstrom)
    //   ALPHA_BETA_GAMMA α β γ (angles, paired with ABC)
    let mut vecs: [Option<[f64; 3]>; 3] = [None; 3];
    let mut abc: Option<[f64; 3]> = None;
    let mut angles: [f64; 3] = [90.0; 3];

    // CP2K 的 &CELL 没有 UNIT 关键字：单位逐关键字写在数值前，`A [bohr] 10 0 0`，
    // 缺省 angstrom / deg（input_cp2k_subsys.F 的 unit_str）
    for l in section {
        let parts: Vec<&str> = l.split_whitespace().collect();
        let Some(key) = parts.first().map(|s| s.to_lowercase()) else { continue };
        let (unit, vals) = match parts.get(1) {
            Some(u) if u.starts_with('[') && u.ends_with(']') => (Some(&u[1..u.len() - 1]), &parts[2..]),
            _ => (None, &parts[1..]),
        };
        let parsed = || -> Result<[f64; 3]> {
            parse_vec3(vals).with_context(|| format!("invalid &CELL line {l:?}"))
        };
        match key.as_str() {
            "a" | "b" | "c" | "abc" => {
                let scale = length_unit(unit.unwrap_or("angstrom"))?;
                let v = parsed()?.map(|x| x * scale);
                match key.as_str() {
                    "a" => vecs[0] = Some(v),
                    "b" => vecs[1] = Some(v),
                    "c" => vecs[2] = Some(v),
                    _ => abc = Some(v),
                }
            }
            "alpha_beta_gamma" => {
                let to_deg = match unit.map(str::to_lowercase).as_deref() {
                    None | Some("deg") => 1.0,
                    Some("rad") => 180.0 / std::f64::consts::PI,
                    Some(u) => bail!("unsupported angle unit [{u}] in &CELL (expected deg or rad)"),
                };
                angles = parsed()?.map(|x| x * to_deg);
            }
            "unit" => bail!(
                "&CELL has no UNIT keyword in CP2K; give the unit per keyword, e.g. A [bohr] 10 0 0"
            ),
            _ => {}
        }
    }

    if let (Some(a), Some(b), Some(c)) = (vecs[0], vecs[1], vecs[2]) {
        Ok(Cell::from_matrix(Matrix3::new(
            a[0], a[1], a[2],
            b[0], b[1], b[2],
            c[0], c[1], c[2],
        )))
    } else if let Some([a, b, c]) = abc {
        Cell::from_lengths_angles(a, b, c, angles[0], angles[1], angles[2])
            .map_err(|e| anyhow::anyhow!("{e}"))
    } else {
        bail!("cannot parse &CELL: need A/B/C vectors or ABC lengths")
    }
}

/// 长度单位 → Å 的倍数。CP2K 认任意单位，这里只认 angstrom 与 bohr，其余报错 ——
/// 此前非 bohr 一律当 Å，`UNIT nm` 会静默差 10 倍
fn length_unit(u: &str) -> Result<f64> {
    match u.to_lowercase().as_str() {
        "angstrom" => Ok(1.0),
        "bohr" => Ok(BOHR),
        other => bail!("unsupported length unit {other:?} (supported: angstrom, bohr)"),
    }
}

/// CP2K 的逻辑值，照 cp_parser_methods.F `parser_get_logical` 的那张表（不区分大小写）；
/// 关键字单写时取 `lone`。表外的值报错 —— CP2K 自己也会 abort
fn cp2k_logical(value: Option<&str>, lone: bool) -> Result<bool> {
    let Some(v) = value else { return Ok(lone) };
    match v.to_uppercase().as_str() {
        "1" | "T" | ".T." | "TRUE" | ".TRUE." | "Y" | "YES" | "ON" => Ok(true),
        "0" | "F" | ".F." | "FALSE" | ".FALSE." | "N" | "NO" | "OFF" => Ok(false),
        other => bail!("{other:?} is not a CP2K logical value"),
    }
}

fn find_end(lines: &[&str], start: usize) -> Option<usize> {
    let _section_name: String = lines[start].to_lowercase()
        .trim_start_matches('&')
        .split_whitespace().next().unwrap_or("").to_string();
    let mut depth = 0usize;
    for (j, l) in lines[start..].iter().enumerate() {
        let low = l.to_lowercase().trim().to_string();
        if low.starts_with('&') && !low.starts_with("&end") {
            depth += 1;
        } else if low.starts_with("&end") {
            depth = depth.saturating_sub(1);
            if depth == 0 { return Some(start + j); }
        }
    }
    None
}

fn parse_vec3(parts: &[&str]) -> Option<[f64; 3]> {
    if parts.len() < 3 { return None; }
    Some([
        parts[0].parse().ok()?,
        parts[1].parse().ok()?,
        parts[2].parse().ok()?,
    ])
}

fn extract_element(label: &str) -> String {
    let alpha: String = label.chars().take_while(|c| c.is_alphabetic()).collect();
    if alpha.is_empty() { return "X".to_string(); }
    let mut it = alpha.chars();
    let first = it.next().unwrap().to_uppercase().to_string();
    let rest: String = it.collect::<String>().to_lowercase();
    format!("{first}{rest}")
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER_INP: &str = "
&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    &CELL
      ABC 10.0 10.0 10.0
    &END CELL
    &COORD
      O  0.000  0.000  0.119
      H  0.000  0.763 -0.477
      H  0.000 -0.763 -0.477
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
";

    const FE_INP: &str = "
&FORCE_EVAL
  &SUBSYS
    &CELL
      A 2.87 0.0  0.0
      B 0.0  2.87 0.0
      C 0.0  0.0  2.87
    &END CELL
    &COORD
      Fe1  0.0   0.0   0.0
      Fe2  1.435 1.435 1.435
    &END COORD
    &KIND Fe1
      ELEMENT Fe
    &END KIND
    &KIND Fe2
      ELEMENT Fe
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
";

    use crate::testutil::write_tmp as tmp;

    #[test]
    fn test_water_abc() {
        let traj = read_cp2k_inp(&tmp("water.inp", WATER_INP)).unwrap();
        let f = traj.first().unwrap();
        assert_eq!(f.n_atoms(), 3);
        assert_eq!(f.atom(0).element, "O");
    }

    #[test]
    fn test_fe_bcc_explicit_vectors() {
        let traj = read_cp2k_inp(&tmp("fe.inp", FE_INP)).unwrap();
        let f = traj.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
        let [a, ..] = f.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-4);
    }

    /// Si 原胞：`cell` 是 &CELL 段内容，`coord_kw` 是 &COORD 里坐标前的关键字行
    fn si(cell: &str, coord_kw: &str, x2: &str) -> String {
        format!("&FORCE_EVAL\n&SUBSYS\n&CELL\n{cell}\n&END CELL\n&COORD\n{coord_kw}\n\
                 Si 0.0 0.0 0.0\nSi {x2} {x2} {x2}\n&END COORD\n&END SUBSYS\n&END FORCE_EVAL\n")
    }
    const FCC: &str = "A 0.0 2.715 2.715\nB 2.715 0.0 2.715\nC 2.715 2.715 0.0";

    /// 第 2 个 Si 的 x 坐标与 a1 的 y 分量
    fn x2_and_a1y(text: &str) -> (f64, f64) {
        let f = parse_cp2k(text).unwrap_or_else(|e| panic!("{e:#}\n{text}"));
        (f.first().unwrap().atom(1).position.x, f.first().unwrap().cell.as_ref().unwrap().matrix[(0, 1)])
    }

    #[test]
    fn test_scaled_follows_cp2k_logical_table() {
        // 期望值手算（ASE 没有 CP2K 输入的 reader）；真值表取自 CP2K 2026.2
        // cp_parser_methods.F parser_get_logical，单写取 lone_keyword_l_val=.TRUE.
        for kw in ["SCALED", "SCALED T", "scaled yes", "SCALED .TRUE.", "SCALED on", "SCALED 1"] {
            let (x, _) = x2_and_a1y(&si(FCC, kw, "0.25"));
            assert!((x - 1.3575).abs() < 1e-9, "{kw:?} 应为真：分数坐标 0.25 → 1.3575 Å，实际 {x}");
        }
        for kw in ["", "SCALED F", "SCALED .FALSE.", "scaled no", "SCALED OFF", "SCALED 0"] {
            let (x, _) = x2_and_a1y(&si(FCC, kw, "1.3575"));
            assert!((x - 1.3575).abs() < 1e-9, "{kw:?} 应为假：坐标按 Å 原样，实际 {x}");
        }
        let err = parse_cp2k(&si(FCC, "SCALED maybe", "0.25")).expect_err("表外的逻辑值应报错");
        assert!(format!("{err:#}").contains("logical"), "{err:#}");
    }

    #[test]
    fn test_units() {
        const B: f64 = 0.52917721;
        // &COORD 的 UNIT；&CELL 的单位逐关键字写在 [] 里
        let (x, _) = x2_and_a1y(&si(FCC, "UNIT bohr", "2.0"));
        assert!((x - 2.0 * B).abs() < 1e-9, "UNIT bohr：{x}");
        let (x, _) = x2_and_a1y(&si(FCC, "UNIT Angstrom", "1.5"));
        assert!((x - 1.5).abs() < 1e-9, "UNIT Angstrom：{x}");
        let (_, a1y) = x2_and_a1y(&si("A [bohr] 0 5.13 5.13\nB [bohr] 5.13 0 5.13\nC [bohr] 5.13 5.13 0", "", "1"));
        assert!((a1y - 5.13 * B).abs() < 1e-9, "A [bohr]：{a1y}");
        let (_, a1y) = x2_and_a1y(&si("ABC [bohr] 10 10 10\nALPHA_BETA_GAMMA [rad] 1.5707963267948966 \
                                       1.5707963267948966 1.5707963267948966", "", "1"));
        assert!(a1y.abs() < 1e-9, "ABC [bohr] + [rad]：a1 应沿 x");

        for (what, text, needle) in [
            ("&COORD UNIT nm", si(FCC, "UNIT nm", "0.1"), "unsupported length unit"),
            ("&CELL [nm]", si("ABC [nm] 1 1 1", "", "0.1"), "unsupported length unit"),
            ("&CELL UNIT", si(&format!("UNIT bohr\n{FCC}"), "", "0.1"), "no UNIT keyword"),
            ("坐标行解析失败", si(FCC, "", "x"), "invalid &COORD line"),
        ] {
            let err = parse_cp2k(&text).expect_err(what);
            assert!(format!("{err:#}").contains(needle), "{what}：报错应含 {needle:?}，实际 {err:#}");
        }
    }

    #[test]
    fn test_fifth_column_is_a_molecule_name_not_a_velocity() {
        // &COORD 第 5 列是分子名；速度在 &VELOCITY 段。此前第 5–7 列能解析成数就当速度
        let text = si(FCC, "", "1.0 1 2 3");
        let f = parse_cp2k(&text).unwrap();
        assert!(f.first().unwrap().velocities.is_none());
    }
}
