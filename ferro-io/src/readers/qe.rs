use std::path::Path;
use std::collections::HashMap;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{bail, Context, Result};

const BOHR: f64 = 0.52917721;

pub fn read_qe_input(path: &Path) -> Result<Trajectory> {
    let path_ = path.display();
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path_}"))?;
    parse_qe(&content).with_context(|| format!("parsing {path_}"))
}

fn parse_qe(content: &str) -> Result<Trajectory> {
    let lines: Vec<&str> = content.lines().collect();

    // Strip comments (`!` and `#`, as pw.x does)
    let stripped: Vec<&str> = lines.iter()
        .map(|l| {
            let t = l.trim();
            if let Some(p) = t.find(['!', '#']) { &t[..p] } else { t }
        })
        .collect();

    // ── &SYSTEM namelist ──────────────────────────────────────────────────────
    let namelists = parse_namelists(&stripped)?;
    let sys: HashMap<String, String> = namelists.iter()
        .find(|n| n.name == "system")
        .context("&SYSTEM namelist not found")?
        .entries.iter()
        .map(|(k, v)| (k.clone(), v.trim_matches('\'').trim_matches('"').to_string()))
        .collect();
    // 写了却解析不了要报错：以前回落 0，空白分隔切坏的值会让非零 ibrav 蒙混过关。
    // 缺省也报错：INPUT_PW 标 ibrav 为 REQUIRED，回落 0 会把丢掉的 &SYSTEM 内容藏起来
    let s = sys.get("ibrav").context("ibrav missing in &SYSTEM (pw.x requires it)")?;
    let ibrav: i32 = s.parse().with_context(|| format!("invalid ibrav = {s:?}"))?;
    if ibrav != 0 {
        bail!("ibrav={ibrav} is not supported; use ibrav=0 with explicit CELL_PARAMETERS");
    }
    let _ntyp: usize = sys.get("ntyp").and_then(|s| s.parse().ok()).unwrap_or(0);

    // starting_magnetization(n) → type_index → magmom
    let mut type_magmom: HashMap<usize, f64> = HashMap::new();
    for (k, v) in &sys {
        if k.starts_with("starting_magnetization(") {
            let idx: usize = k.trim_start_matches("starting_magnetization(")
                .trim_end_matches(')')
                .parse().unwrap_or(0);
            if let Ok(m) = v.parse::<f64>() { type_magmom.insert(idx, m); }
        }
    }

    // ── ATOMIC_SPECIES card ───────────────────────────────────────────────────
    // label → element (first line word treated as label, second as mass (ignored), third as pseudo)
    let mut species: HashMap<String, (String, usize)> = HashMap::new(); // label → (element, type_index)
    if let Some(start) = find_card(&stripped, "ATOMIC_SPECIES") {
        let mut tidx = 1usize;
        for l in stripped[start+1..].iter() {
            let l = l.trim();
            if l.is_empty() || is_card_or_namelist(l) { break; }
            let parts: Vec<&str> = l.split_whitespace().collect();
            if parts.len() >= 2 {
                let label = parts[0].to_string();
                // Try to extract element: first alphabetic chars of label
                let elem = extract_element(&label);
                species.insert(label, (elem, tidx));
                tidx += 1;
            }
        }
    }

    // ── 晶格参数 alat ─────────────────────────────────────────────────────────
    // 口径照 pw.x 源码（Modules/cell_base.f90 cell_base_init）：celldm(1) [Bohr] 优先，
    // 其次 A [Å]；两者为 0 视同未给
    let lat = |key: &str| -> Result<Option<f64>> {
        match sys.get(key) {
            None => Ok(None),
            Some(v) => {
                let x = ffloat(v).with_context(|| format!("invalid {key} = {v}"))?;
                Ok((x != 0.0).then_some(x))
            }
        }
    };
    let alat: Option<f64> = match lat("celldm(1)")? {
        Some(c) => Some(c * BOHR),
        None => lat("a")?,
    };

    // ── CELL_PARAMETERS card ──────────────────────────────────────────────────
    let cell_start = find_card(&stripped, "CELL_PARAMETERS").context("CELL_PARAMETERS card not found")?;
    let cell_opt = card_option(stripped[cell_start], "CELL_PARAMETERS");
    let cell_scale = match cell_opt.as_str() {
        "angstrom" | "bohr" if alat.is_some() => bail!(
            "CELL_PARAMETERS {cell_opt} together with celldm(1) or A in &SYSTEM: the lattice \
             parameter is given twice (pw.x rejects this too)"
        ),
        "angstrom" => 1.0,
        "bohr" => BOHR,
        "alat" => alat.context("CELL_PARAMETERS alat needs celldm(1) or A in &SYSTEM")?,
        // pw.x 的旧默认：有 celldm(1)/A 按 alat，否则按 bohr（ASE 同）
        "" => alat.unwrap_or(BOHR),
        other => bail!("unknown CELL_PARAMETERS option {other:?} (expected alat, bohr or angstrom)"),
    };

    let mut cell_vecs: Vec<Vector3<f64>> = Vec::with_capacity(3);
    for l in stripped[cell_start + 1..].iter() {
        if cell_vecs.len() == 3 { break; }
        let l = l.trim();
        if l.is_empty() { continue; }
        if is_card_or_namelist(l) { break; }
        let v = fortran_vec3(l).with_context(|| format!("invalid CELL_PARAMETERS line {l:?}"))?;
        cell_vecs.push(v * cell_scale);
    }
    anyhow::ensure!(cell_vecs.len() == 3, "CELL_PARAMETERS must have 3 vectors");
    let cell = Cell::from_matrix(Matrix3::from_rows(&[
        cell_vecs[0].transpose(), cell_vecs[1].transpose(), cell_vecs[2].transpose(),
    ]));

    // ── ATOMIC_POSITIONS card ─────────────────────────────────────────────────
    let pos_start = find_card(&stripped, "ATOMIC_POSITIONS").context("ATOMIC_POSITIONS card not found")?;
    let pos_opt = card_option(stripped[pos_start], "ATOMIC_POSITIONS");
    // 坐标的 alat 在没有 celldm(1)/A 时取 |a1|（pw.x 同上处；ASE 不支持这一情形）
    let pos_alat = alat.unwrap_or_else(|| cell.matrix.row(0).norm());
    let to_cart = |v: Vector3<f64>| -> Vector3<f64> {
        match pos_opt.as_str() {
            "angstrom" => v,
            "bohr" => v * BOHR,
            "crystal" => cell.fractional_to_cartesian(v),
            _ => v * pos_alat, // alat 或未写（pw.x 的默认）
        }
    };
    match pos_opt.as_str() {
        "angstrom" | "bohr" | "crystal" | "alat" | "" => {}
        "crystal_sg" => bail!("ATOMIC_POSITIONS crystal_sg (Wyckoff positions) is not supported"),
        other => bail!(
            "unknown ATOMIC_POSITIONS option {other:?} (expected alat, bohr, angstrom or crystal)"
        ),
    }

    let mut frame = Frame::with_cell(cell.clone(), [true; 3]);

    for l in stripped[pos_start+1..].iter() {
        let l = l.trim();
        if l.is_empty() { continue; }
        if is_card_or_namelist(l) { break; }
        let parts: Vec<&str> = l.split_whitespace().collect();
        let label = parts[0].to_string();
        let v = fortran_vec3(&parts[1..].join(" "))
            .with_context(|| format!("invalid ATOMIC_POSITIONS line {l:?}"))?;
        // 未声明的标签 pw.x 拒收；以前按标签字母猜元素，磁矩也随之静默丢失
        let (elem, tidx) = species.get(&label).cloned().with_context(|| {
            format!("ATOMIC_POSITIONS label {label:?} is not declared in ATOMIC_SPECIES")
        })?;

        let mut atom = Atom::new(elem, to_cart(v));
        atom.label = Some(label);
        atom.magmom = type_magmom.get(&tidx).copied();
        frame.add_atom(atom);
    }
    if let Some(nat) = sys.get("nat")
        .map(|s| s.parse::<usize>().with_context(|| format!("invalid nat = {s:?}")))
        .transpose()?
    {
        anyhow::ensure!(
            frame.n_atoms() == nat,
            "ATOMIC_POSITIONS has {} atoms but nat = {nat}", frame.n_atoms()
        );
    }

    Ok(Trajectory::from_frame(frame))
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

/// 一个 namelist：名字小写；条目按出现顺序，键小写，值保持原样（含引号）
struct Namelist {
    name: String,
    entries: Vec<(String, String)>,
}

/// 全部 namelist。按 token 流收：`&name` 之后（同行剩余也算）到引号外的 `/`
/// （或旧式 `&end`）为止，可跨行。以前按行收、要求 `&name` 独占一行，
/// `&SYSTEM ibrav=0, celldm(1)=10.26, … /` 写在一行（Fortran 合法，pw.x 照收）时
/// 整段被忽略，celldm(1) 丢失、无选项的 CELL_PARAMETERS 回落 Bohr（审查 A-13）
fn parse_namelists(lines: &[&str]) -> Result<Vec<Namelist>> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        let Some(rest) = lines[i].trim().strip_prefix('&') else { i += 1; continue };
        let (name, first) = rest.split_once(char::is_whitespace).unwrap_or((rest, ""));
        let first_line = i + 1;
        let mut body = first.to_string();
        let end = loop {
            if let Some(end) = namelist_end(&body) { break end; }
            i += 1;
            anyhow::ensure!(i < lines.len(), "namelist &{name} (line {first_line}) is not closed with /");
            body.push('\n');
            body.push_str(lines[i]);
        };
        anyhow::ensure!(
            body[end..].starts_with('/') || body[end..].to_lowercase().starts_with("&end"),
            "namelist &{name} (line {first_line}) is not closed with / before the next namelist"
        );
        let entries = namelist_items(&body[..end]).into_iter()
            .filter_map(|item| {
                let (k, v) = item.split_once('=')?;
                (!k.is_empty()).then(|| (k.to_lowercase(), v.to_string()))
            })
            .collect();
        out.push(Namelist { name: name.to_lowercase(), entries });
        i += 1;
    }
    Ok(out)
}

/// namelist 正文里引号外第一个 `/` 或 `&` 的字节位置（`pseudo_dir='./pp/'` 里的
/// `/` 不算）；`&` 是否为合法的 `&end` 由调用方判断
fn namelist_end(body: &str) -> Option<usize> {
    let mut quote = None;
    for (p, c) in body.char_indices() {
        match quote {
            Some(q) => if c == q { quote = None },
            None if c == '\'' || c == '"' => quote = Some(c),
            None if c == '/' || c == '&' => return Some(p),
            None => {}
        }
    }
    None
}

/// namelist 一行 → `key=value` 项。Fortran namelist 里逗号与空白都是分隔符
/// （`ibrav=0 nat=3 ntyp=1` 合法，pw.x 与 ASE 都认）；`=` 两侧、括号内的空白不算，
/// 引号内原样。以前只按逗号切，空白分隔的一整行成了 `ibrav` 一个键的值
fn namelist_items(line: &str) -> Vec<String> {
    let mut words: Vec<String> = Vec::new();
    let mut cur = String::new();
    let (mut quote, mut paren) = (None, 0usize);
    for c in line.chars() {
        if let Some(q) = quote {
            cur.push(c);
            if c == q { quote = None; }
            continue;
        }
        match c {
            '\'' | '"' => { quote = Some(c); cur.push(c); }
            '(' => { paren += 1; cur.push(c); }
            ')' => { paren = paren.saturating_sub(1); cur.push(c); }
            c if paren > 0 && c.is_whitespace() => {}
            c if c == ',' || c.is_whitespace() => {
                if !cur.is_empty() { words.push(std::mem::take(&mut cur)); }
            }
            _ => cur.push(c),
        }
    }
    if !cur.is_empty() { words.push(cur); }
    // `key = value` / `key= value` / `key =value` 拆成了几个词，粘回 `key=value`
    let mut items: Vec<String> = Vec::new();
    for w in words {
        match items.last_mut() {
            Some(last) if last.ends_with('=') || w.starts_with('=') => last.push_str(&w),
            _ => items.push(w),
        }
    }
    items
}

fn find_card(lines: &[&str], card: &str) -> Option<usize> {
    lines.iter().position(|l| {
        l.to_uppercase().split_whitespace().next() == Some(card)
    })
}

/// 卡片行上的选项：去掉卡片名后忽略 `{}()`、空白与大小写。`{bohr}`、`(bohr)`、
/// `bohr` 是同一个选项（pw.x 与 ASE 都按关键字匹配，不要求括号）；未写返回空串
fn card_option(line: &str, card: &str) -> String {
    line.trim()[card.len()..]
        .chars()
        .filter(|c| !matches!(c, '{' | '}' | '(' | ')') && !c.is_whitespace())
        .collect::<String>()
        .to_lowercase()
}

/// Fortran 实数：指数字母可写 d/D（`0.25d0`、`1.026D+01`），QE 输入里常见
fn ffloat(s: &str) -> Option<f64> {
    s.replace(['d', 'D'], "e").parse().ok()
}

/// 行首三个 Fortran 实数；不足三个或解析失败即 `Err`，不补零
fn fortran_vec3(line: &str) -> Result<Vector3<f64>> {
    let v: Vec<f64> = line.split_whitespace().take(3).map_while(ffloat).collect();
    anyhow::ensure!(v.len() == 3, "expected 3 numbers");
    Ok(Vector3::new(v[0], v[1], v[2]))
}

fn is_card_or_namelist(l: &str) -> bool {
    let up = l.to_uppercase();
    matches!(
        up.split_whitespace().next().unwrap_or(""),
        "ATOMIC_POSITIONS" | "ATOMIC_SPECIES" | "CELL_PARAMETERS" |
        "K_POINTS" | "CONSTRAINTS" | "ATOMIC_FORCES"
    ) || up.trim().starts_with('&')
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

    const WATER_QE: &str = "
&CONTROL
  calculation = 'scf',
  prefix = 'water',
/
&SYSTEM
  ibrav = 0,
  nat = 3,
  ntyp = 2,
  ecutwfc = 30.0,
/
&ELECTRONS
/
ATOMIC_SPECIES
O  15.999  O.pbesol.UPF
H  1.008   H.pbesol.UPF

ATOMIC_POSITIONS {angstrom}
O  0.000  0.000  0.119
H  0.000  0.763 -0.477
H  0.000 -0.763 -0.477

CELL_PARAMETERS {angstrom}
10.0  0.0  0.0
0.0  10.0  0.0
0.0   0.0 10.0

K_POINTS {gamma}
";

    const CRYSTAL_POS: &str = "
&SYSTEM
  ibrav = 0,
  nat = 2,
  ntyp = 1,
/
ATOMIC_SPECIES
Fe  55.845  Fe.UPF

ATOMIC_POSITIONS {crystal}
Fe  0.0  0.0  0.0
Fe  0.5  0.5  0.5

CELL_PARAMETERS {angstrom}
2.87  0.0   0.0
0.0   2.87  0.0
0.0   0.0   2.87
";

    use crate::testutil::write_tmp as tmp;

    #[test]
    fn test_water() {
        let traj = read_qe_input(&tmp("water.qe", WATER_QE)).unwrap();
        let f = traj.first().unwrap();
        assert_eq!(f.n_atoms(), 3);
        assert_eq!(f.atom(0).element, "O");
        assert_eq!(f.atom(1).element, "H");
    }

    /// Si 原胞输入：`sys` 追加进 &SYSTEM，`cell`/`pos` 是两张卡的表头行
    fn si(sys: &str, cell_hdr: &str, cell: &str, pos_hdr: &str, pos: &str) -> String {
        format!("&SYSTEM\n ibrav=0, nat=2, ntyp=1{sys}\n/\nATOMIC_SPECIES\nSi 28.086 Si.upf\n\
                 {cell_hdr}\n{cell}\n{pos_hdr}\n{pos}\n")
    }
    const FCC_BOHR: &str = "0.0 5.13 5.13\n5.13 0.0 5.13\n5.13 5.13 0.0";
    const FCC_ANG: &str = "0.0 2.715 2.715\n2.715 0.0 2.715\n2.715 2.715 0.0";
    const FCC_UNIT: &str = "0.0 0.5 0.5\n0.5 0.0 0.5\n0.5 0.5 0.0";
    const P: &str = "Si 0.00 0.00 0.00\nSi 0.25 0.25 0.25";

    #[test]
    fn test_card_options_and_units_match_pw_x() {
        // 期望值：除注明者外取自 ASE 3.29.0 read_espresso_in 读同一份输入
        // （ASE 的 Bohr 是 CODATA 2014，与本仓差 ~1e-9 相对，容差 1e-6 Å）
        let (b, c) = (2.714679078854, 1.357339539427); // 5.13 Bohr、0.25×10.26 Bohr
        let cases: [(&str, String, f64, f64); 6] = [
            ("bohr + crystal（不带括号）",
             si("", "CELL_PARAMETERS bohr", FCC_BOHR, "ATOMIC_POSITIONS crystal", P), b, c),
            ("angstrom + (crystal)",
             si("", "CELL_PARAMETERS angstrom", FCC_ANG, "ATOMIC_POSITIONS (crystal)", P), 2.715, 1.3575),
            ("两张卡都不写选项：有 celldm(1) 即 alat",
             si(", celldm(1)=10.26", "CELL_PARAMETERS", FCC_UNIT, "ATOMIC_POSITIONS", P), b, c),
            // 手算：ASE 查 'A' 而它的 namelist 键已转小写，A 在 ASE 里从不生效
            ("A + {alat} / alat",
             si(", A=5.43", "CELL_PARAMETERS {alat}", FCC_UNIT, "ATOMIC_POSITIONS alat", P), 2.715, 1.3575),
            // 与第三条数值等价；ASE 的 namelist 解析读不了 1.026D+01
            ("Fortran d 指数",
             si(", celldm(1)=1.026D+01", "CELL_PARAMETERS alat", &FCC_UNIT.replace("0.5", "0.5d0"),
                "ATOMIC_POSITIONS crystal", "Si 0.0d0 0.0d0 0.0d0\nSi 0.25d0 0.25D0 2.5d-1"), b, c),
            // 手算：无 celldm/A 时坐标的 alat = |a1|（pw.x cell_base_init；ASE 不支持）
            ("angstrom 胞 + alat 坐标",
             si("", "CELL_PARAMETERS angstrom", FCC_ANG, "ATOMIC_POSITIONS alat", P), 2.715, 0.959897455461),
        ];
        for (what, text, a12, p2) in cases {
            let traj = parse_qe(&text).unwrap_or_else(|e| panic!("{what}：{e:#}"));
            let f = traj.first().unwrap();
            let m = f.cell.as_ref().unwrap().matrix;
            assert!((m[(0, 1)] - a12).abs() < 1e-6 && m[(0, 0)].abs() < 1e-12, "{what}：a1 = {:?}", m.row(0));
            let x = f.atom(1).position;
            assert!((x - Vector3::repeat(p2)).norm() < 1e-6, "{what}：第 2 个 Si 在 {x:?}，应为 {p2}");
        }
    }

    #[test]
    fn test_namelist_items_split_on_commas_and_blanks() {
        assert_eq!(namelist_items("ibrav=0 nat=3 ntyp=1"), ["ibrav=0", "nat=3", "ntyp=1"]);
        assert_eq!(namelist_items("ibrav = 0, nat= 3 ,ntyp =1,"), ["ibrav=0", "nat=3", "ntyp=1"]);
        assert_eq!(namelist_items("starting_magnetization( 1 ) = 0.5 title='a b, c'"),
                   ["starting_magnetization(1)=0.5", "title='a b, c'"]);
    }

    #[test]
    fn test_blank_separated_namelist() {
        let fe = |sys: &str| read_qe_input(&tmp("fe_blank.qe", &CRYSTAL_POS
            .replace("  ibrav = 0,\n  nat = 2,\n  ntyp = 1,\n", sys)));
        // 空白分隔：nat 校验生效、磁矩读到
        let traj = fe("  ibrav=0 nat=2 ntyp=1 starting_magnetization(1)=0.5\n").unwrap();
        assert_eq!(traj.first().unwrap().atom(0).magmom, Some(0.5));
        let e = fe("  ibrav=0 nat=3 ntyp=1\n").unwrap_err();
        assert!(format!("{e:#}").contains("nat = 3"), "{e:#}");
        // 非零 ibrav 不再因切坏的值蒙混过关
        let e = fe("  ibrav=2 nat=2 ntyp=1\n").unwrap_err();
        assert!(format!("{e:#}").contains("ibrav=2"), "{e:#}");
        // 写坏的值报错
        let e = fe("  ibrav=0, nat=two, ntyp=1\n").unwrap_err();
        assert!(format!("{e:#}").contains("invalid nat"), "{e:#}");
    }

    #[test]
    fn test_one_line_namelists() {
        // 审查 A-13 的复现：&SYSTEM 写在一行，无选项的 CELL_PARAMETERS 须按 celldm(1) 取 alat。
        // 修前 celldm(1) 丢失、胞回落 Bohr：a = 0.3742 Å
        let text = "&CONTROL calculation='scf', pseudo_dir='./pp/' /\n\
                    &SYSTEM ibrav=0, celldm(1)=10.26, nat=2, ntyp=1 /\n\
                    &ELECTRONS /\n\
                    ATOMIC_SPECIES\nSi 28.086 Si.upf\n\
                    CELL_PARAMETERS\n{FCC_UNIT}\nATOMIC_POSITIONS crystal\n{P}\n"
            .replace("{FCC_UNIT}", FCC_UNIT).replace("{P}", P);
        let traj = parse_qe(&text).unwrap();
        let m = traj.first().unwrap().cell.as_ref().unwrap().matrix;
        assert!((m[(0, 1)] - 2.714679078854).abs() < 1e-6, "一行 namelist 的 celldm(1) 丢了：a1 = {:?}", m.row(0));

        // 与上同义的多行跨行写法、`/` 与末项同行
        let split = text.replace("&SYSTEM ibrav=0, celldm(1)=10.26, nat=2, ntyp=1 /",
                                 "&SYSTEM\n ibrav=0,\n celldm(1)=10.26, nat=2,\n ntyp=1 /");
        let m2 = parse_qe(&split).unwrap().first().unwrap().cell.as_ref().unwrap().matrix;
        assert!((m - m2).norm() < 1e-12, "跨行写法与一行写法结果应相同");

        let nl = parse_namelists(&["&control calculation='scf', pseudo_dir='./pp/' /"]).unwrap();
        assert_eq!(nl[0].entries, [("calculation".into(), "'scf'".into()),
                                   ("pseudo_dir".into(), "'./pp/'".into())],
                   "引号内的 / 不应结束 namelist，值保持原样");
        let old = parse_namelists(&["&SYSTEM", "ibrav=0", "&END"]).unwrap();
        assert_eq!(old[0].entries, [("ibrav".into(), "0".into())], "旧式 &END 结束");
    }

    #[test]
    fn test_card_options_that_must_fail() {
        let cases = [
            ("celldm 与 bohr 同时给", si(", celldm(1)=10.26", "CELL_PARAMETERS bohr", FCC_BOHR,
                                       "ATOMIC_POSITIONS crystal", P), "given twice"),
            ("crystal_sg", si("", "CELL_PARAMETERS bohr", FCC_BOHR, "ATOMIC_POSITIONS crystal_sg", P), "crystal_sg"),
            ("未知选项", si("", "CELL_PARAMETERS {bohrr}", FCC_BOHR, "ATOMIC_POSITIONS crystal", P), "unknown"),
            ("alat 胞缺 celldm", si("", "CELL_PARAMETERS alat", FCC_UNIT, "ATOMIC_POSITIONS crystal", P), "celldm"),
            ("坐标解析失败", si("", "CELL_PARAMETERS bohr", FCC_BOHR, "ATOMIC_POSITIONS crystal",
                               "Si 0 0 0\nSi 0.25 x 0.25"), "ATOMIC_POSITIONS line"),
            ("原子数与 nat 不符", si("", "CELL_PARAMETERS bohr", FCC_BOHR, "ATOMIC_POSITIONS crystal",
                                   "Si 0 0 0"), "nat = 2"),
            ("未声明的标签", si("", "CELL_PARAMETERS bohr", FCC_BOHR, "ATOMIC_POSITIONS crystal",
                               "Si 0 0 0\nSi1 0.25 0.25 0.25"), "not declared"),
            ("缺 ibrav", si("", "CELL_PARAMETERS bohr", FCC_BOHR, "ATOMIC_POSITIONS crystal", P)
                .replace("ibrav=0, ", ""), "ibrav missing"),
            ("缺 &SYSTEM", "ATOMIC_SPECIES\nSi 28.086 Si.upf\n".into(), "&SYSTEM namelist not found"),
            ("namelist 未闭合", si("", "CELL_PARAMETERS bohr", FCC_BOHR, "ATOMIC_POSITIONS crystal", P)
                .replace("\n/\n", "\n"), "not closed"),
        ];
        for (what, text, needle) in cases {
            let err = parse_qe(&text).expect_err(what);
            assert!(format!("{err:#}").contains(needle), "{what}：报错应含 {needle:?}，实际 {err:#}");
        }
    }

    #[test]
    fn test_crystal_coords() {
        let traj = read_qe_input(&tmp("bcc.qe", CRYSTAL_POS)).unwrap();
        let f = traj.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        // Fe at (0,0,0) should be at Cartesian (0,0,0)
        assert!(f.atom(0).position.norm() < 1e-6);
        // Fe at (0.5,0.5,0.5) fractional = (1.435,1.435,1.435) Cartesian
        assert!((f.atom(1).position.x - 1.435).abs() < 1e-3);
    }
}
