use std::collections::HashMap;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{bail, ensure, Context, Result};

// ─── Tokenizer ───────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
enum Token {
    DataBlock(String),
    Loop,
    Tag(String),
    Value(String),
}

fn tokenize(input: &str) -> Vec<Token> {
    let mut tokens = Vec::new();
    let lines: Vec<&str> = input.lines().collect();
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i];

        // Semicolon-delimited multiline string — MUST start at column 0
        if line.starts_with(';') {
            let mut text = String::new();
            i += 1;
            while i < lines.len() && !lines[i].starts_with(';') {
                if !text.is_empty() { text.push('\n'); }
                text.push_str(lines[i]);
                i += 1;
            }
            tokens.push(Token::Value(text.trim().to_string()));
            i += 1; // closing ';'
            continue;
        }

        let line = strip_comment(line);
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            tokenize_line(trimmed, &mut tokens);
        }
        i += 1;
    }

    tokens
}

fn strip_comment(line: &str) -> &str {
    let mut in_sq = false;
    let mut in_dq = false;
    for (idx, ch) in line.char_indices() {
        match ch {
            '\'' if !in_dq => in_sq = !in_sq,
            '"' if !in_sq => in_dq = !in_dq,
            '#' if !in_sq && !in_dq => return &line[..idx],
            _ => {}
        }
    }
    line
}

fn tokenize_line(mut s: &str, out: &mut Vec<Token>) {
    loop {
        s = s.trim_start();
        if s.is_empty() { break; }

        if let Some(quote) = s.chars().next().filter(|&c| c == '\'' || c == '"') {
            s = &s[1..];
            let end = s.find(quote).unwrap_or(s.len());
            out.push(Token::Value(s[..end].to_string()));
            s = if end < s.len() { &s[end + 1..] } else { "" };
        } else {
            let end = s.find(char::is_whitespace).unwrap_or(s.len());
            let tok = &s[..end];
            out.push(make_token(tok));
            s = &s[end..];
        }
    }
}

fn make_token(s: &str) -> Token {
    let low = s.to_lowercase();
    if low == "loop_" {
        Token::Loop
    } else if low.starts_with("data_") {
        Token::DataBlock(s[5..].to_string())
    } else if s.starts_with('_') {
        Token::Tag(low)
    } else {
        Token::Value(s.to_string())
    }
}

// ─── Block parser ─────────────────────────────────────────────────────────────

#[derive(Debug, Default)]
struct CifBlock {
    name: String,
    singles: HashMap<String, String>,
    loops: Vec<(Vec<String>, Vec<Vec<String>>)>,
}

fn parse_blocks(tokens: &[Token]) -> Vec<CifBlock> {
    let mut blocks = Vec::new();
    let mut current: Option<CifBlock> = None;
    let mut i = 0;

    while i < tokens.len() {
        match &tokens[i] {
            Token::DataBlock(name) => {
                if let Some(b) = current.take() { blocks.push(b); }
                current = Some(CifBlock { name: name.clone(), ..Default::default() });
                i += 1;
            }
            Token::Tag(tag) => {
                let tag = tag.clone();
                i += 1;
                if let Some(Token::Value(val)) = tokens.get(i) {
                    if let Some(b) = current.as_mut() {
                        b.singles.insert(tag, val.clone());
                    }
                    i += 1;
                }
            }
            Token::Loop => {
                i += 1;
                let mut headers = Vec::new();
                while let Some(Token::Tag(h)) = tokens.get(i) {
                    headers.push(h.clone());
                    i += 1;
                }
                let ncols = headers.len();
                if ncols == 0 { continue; }
                let mut rows: Vec<Vec<String>> = Vec::new();
                let mut row: Vec<String> = Vec::new();
                while let Some(Token::Value(v)) = tokens.get(i) {
                    row.push(v.clone());
                    if row.len() == ncols {
                        rows.push(std::mem::take(&mut row));
                    }
                    i += 1;
                }
                if let Some(b) = current.as_mut() {
                    b.loops.push((headers, rows));
                }
            }
            Token::Value(_) => { i += 1; }
        }
    }
    if let Some(b) = current { blocks.push(b); }
    blocks
}

// ─── Cell ─────────────────────────────────────────────────────────────────────

fn parse_cell(block: &CifBlock) -> Result<Cell> {
    let get = |tag: &str| -> Result<f64> {
        let s = block.singles.get(tag)
            .with_context(|| format!("missing {tag}"))?;
        parse_cif_float(s).with_context(|| format!("invalid value for {tag}: {s}"))
    };
    Cell::from_lengths_angles(
        get("_cell_length_a")?,
        get("_cell_length_b")?,
        get("_cell_length_c")?,
        get("_cell_angle_alpha")?,
        get("_cell_angle_beta")?,
        get("_cell_angle_gamma")?,
    ).map_err(|e| anyhow::anyhow!("{e}"))
}

// ─── Atom sites ───────────────────────────────────────────────────────────────

struct AtomSites {
    elements: Vec<String>,
    labels: Vec<Option<String>>,
    frac: Vec<Vector3<f64>>,
}

fn parse_atom_sites(block: &CifBlock, cell: &Cell) -> Result<AtomSites> {
    let atom_loop = block.loops.iter().find(|(hdrs, _)| {
        hdrs.iter().any(|h| h.starts_with("_atom_site_"))
    });

    let (headers, rows) = match atom_loop {
        Some(l) => l,
        None => bail!("no _atom_site_ loop found"),
    };

    let col = |names: &[&str]| -> Option<usize> {
        names.iter().find_map(|n| headers.iter().position(|h| h == n))
    };

    let idx_sym = col(&["_atom_site_type_symbol"]);
    let idx_lbl = col(&["_atom_site_label"]);
    let idx_fx  = col(&["_atom_site_fract_x"]);
    let idx_fy  = col(&["_atom_site_fract_y"]);
    let idx_fz  = col(&["_atom_site_fract_z"]);
    let idx_cx  = col(&["_atom_site_cartn_x"]);
    let idx_cy  = col(&["_atom_site_cartn_y"]);
    let idx_cz  = col(&["_atom_site_cartn_z"]);

    ensure!(idx_sym.is_some() || idx_lbl.is_some(), "no element/label column");
    ensure!(
        (idx_fx.is_some() && idx_fy.is_some() && idx_fz.is_some())
            || (idx_cx.is_some() && idx_cy.is_some() && idx_cz.is_some()),
        "no position columns (need fract_x/y/z or cartn_x/y/z)"
    );

    let mut sites = AtomSites { elements: vec![], labels: vec![], frac: vec![] };

    for row in rows {
        let get = |i: usize| row.get(i).map(|s| s.as_str()).unwrap_or("?");

        let element = idx_sym
            .map(|i| clean_element(get(i)))
            .or_else(|| idx_lbl.map(|i| element_from_label(get(i))))
            .unwrap_or_else(|| "X".to_string());

        let label = idx_lbl.map(|i| get(i).to_string());

        let frac = if let (Some(ix), Some(iy), Some(iz)) = (idx_fx, idx_fy, idx_fz) {
            let fx = parse_cif_float(get(ix)).unwrap_or(0.0);
            let fy = parse_cif_float(get(iy)).unwrap_or(0.0);
            let fz = parse_cif_float(get(iz)).unwrap_or(0.0);
            Vector3::new(fx, fy, fz)
        } else {
            // Cartesian → fractional
            let x = parse_cif_float(get(idx_cx.unwrap())).unwrap_or(0.0);
            let y = parse_cif_float(get(idx_cy.unwrap())).unwrap_or(0.0);
            let z = parse_cif_float(get(idx_cz.unwrap())).unwrap_or(0.0);
            cell.cartesian_to_fractional(Vector3::new(x, y, z))
        };

        sites.elements.push(element);
        sites.labels.push(label);
        sites.frac.push(frac);
    }

    Ok(sites)
}

// ─── Symmetry ─────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct SymOp {
    rot: Matrix3<f64>,
    trans: Vector3<f64>,
}

impl SymOp {
    fn identity() -> Self {
        Self { rot: Matrix3::identity(), trans: Vector3::zeros() }
    }

    fn apply(&self, frac: Vector3<f64>) -> Vector3<f64> {
        self.rot * frac + self.trans
    }
}

fn collect_symops(block: &CifBlock) -> Vec<SymOp> {
    // Tag names used by old and new CIF dictionaries
    const TAGS: &[&str] = &[
        "_space_group_symop_operation_xyz",
        "_symmetry_equiv_pos_as_xyz",
    ];

    // 1. Try loop
    for (headers, rows) in &block.loops {
        for tag in TAGS {
            if let Some(col) = headers.iter().position(|h| h == tag) {
                let ops: Vec<SymOp> = rows.iter()
                    .filter_map(|row| row.get(col))
                    .filter_map(|s| parse_symop(s).ok())
                    .collect();
                if !ops.is_empty() { return ops; }
            }
        }
    }

    // 2. Try single tag
    for tag in TAGS {
        if let Some(val) = block.singles.get(*tag) {
            if let Ok(op) = parse_symop(val) {
                return vec![op];
            }
        }
    }

    // 3. No symops: identity only
    vec![SymOp::identity()]
}

fn parse_symop(s: &str) -> Result<SymOp> {
    let s = s.trim().replace(' ', "").to_lowercase();
    // Strip surrounding quotes
    let s = s.trim_matches('\'').trim_matches('"');
    let parts: Vec<&str> = s.split(',').collect();
    ensure!(parts.len() == 3, "symop must have 3 components, got: {s}");

    let mut rot = Matrix3::zeros();
    let mut trans = Vector3::zeros();
    for (row, part) in parts.iter().enumerate() {
        parse_symop_component(part, row, &mut rot, &mut trans)?;
    }
    Ok(SymOp { rot, trans })
}

fn parse_symop_component(
    s: &str,
    row: usize,
    rot: &mut Matrix3<f64>,
    trans: &mut Vector3<f64>,
) -> Result<()> {
    let mut sign = 1.0_f64;
    let mut term = String::new();

    for ch in s.chars() {
        match ch {
            '+' | '-' => {
                apply_term(&term, sign, row, rot, trans);
                term.clear();
                sign = if ch == '+' { 1.0 } else { -1.0 };
            }
            c => term.push(c),
        }
    }
    apply_term(&term, sign, row, rot, trans);
    Ok(())
}

fn apply_term(term: &str, sign: f64, row: usize, rot: &mut Matrix3<f64>, trans: &mut Vector3<f64>) {
    match term {
        "x" => rot[(row, 0)] += sign,
        "y" => rot[(row, 1)] += sign,
        "z" => rot[(row, 2)] += sign,
        "" => {}
        s => { if let Some(v) = parse_fraction(s) { trans[row] += sign * v; } }
    }
}

fn parse_fraction(s: &str) -> Option<f64> {
    if let Some((n, d)) = s.split_once('/') {
        let n: f64 = n.parse().ok()?;
        let d: f64 = d.parse().ok()?;
        if d == 0.0 { return None; }
        Some(n / d)
    } else {
        s.parse().ok()
    }
}

// ─── Symmetry expansion ───────────────────────────────────────────────────────

fn expand_by_symmetry(sites: AtomSites, symops: &[SymOp]) -> AtomSites {
    const TOL_SQ: f64 = 1e-4; // (0.01)² in fractional space

    let mut out = AtomSites { elements: vec![], labels: vec![], frac: vec![] };

    for ((element, label), &frac) in sites.elements.iter()
        .zip(sites.labels.iter())
        .zip(sites.frac.iter())
    {
        for op in symops {
            let mut p = op.apply(frac);
            // Wrap to [0, 1)
            p.x = p.x.rem_euclid(1.0);
            p.y = p.y.rem_euclid(1.0);
            p.z = p.z.rem_euclid(1.0);

            // Duplicate check using minimum-image fractional distance
            let is_dup = out.frac.iter().zip(out.elements.iter()).any(|(&q, e)| {
                e == element && frac_dist_sq(p, q) < TOL_SQ
            });
            if !is_dup {
                out.elements.push(element.clone());
                out.labels.push(label.clone());
                out.frac.push(p);
            }
        }
    }
    out
}

fn frac_dist_sq(a: Vector3<f64>, b: Vector3<f64>) -> f64 {
    let mut d = a - b;
    d.x -= d.x.round();
    d.y -= d.y.round();
    d.z -= d.z.round();
    d.norm_squared()
}

// ─── Block → Frame ────────────────────────────────────────────────────────────

fn block_to_frame(block: &CifBlock) -> Result<Frame> {
    let cell = parse_cell(block)?;
    let sites = parse_atom_sites(block, &cell)?;

    let symops = collect_symops(block);
    let sites = if symops.len() > 1 {
        expand_by_symmetry(sites, &symops)
    } else {
        sites
    };

    let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
    for ((element, label), frac) in sites.elements.iter()
        .zip(sites.labels.iter())
        .zip(sites.frac.iter())
    {
        let cart = cell.fractional_to_cartesian(*frac);
        let mut atom = Atom::new(element.clone(), cart);
        atom.label = label.clone();
        frame.add_atom(atom);
    }
    Ok(frame)
}

// ─── Public API ───────────────────────────────────────────────────────────────

/// 读取 CIF 文件，每个 `data_` block 对应轨迹中的一帧。
pub fn read_cif(path: &str) -> Result<Trajectory> {
    let input = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path}"))?;
    let tokens = tokenize(&input);
    let blocks = parse_blocks(&tokens);

    ensure!(!blocks.is_empty(), "no data_ blocks found in {path}");

    let mut traj = Trajectory::new();
    if let Some(name) = blocks.first().map(|b| b.name.clone()) {
        if !name.is_empty() { traj.metadata.source = Some(name); }
    }

    for block in &blocks {
        let frame = block_to_frame(block)
            .with_context(|| format!("error in data_{}", block.name))?;
        traj.add_frame(frame);
    }
    Ok(traj)
}

// ─── Helpers ──────────────────────────────────────────────────────────────────

/// 解析 CIF 浮点值，剥离不确定度括号如 `5.431(2)` → `5.431`。
/// `?` 和 `.` 返回 None。
fn parse_cif_float(s: &str) -> Option<f64> {
    match s.trim() {
        "?" | "." => None,
        s => {
            let s = s.split('(').next().unwrap_or(s);
            s.parse().ok()
        }
    }
}

/// 从原子标签提取元素符号（如 "Fe1a" → "Fe", "H12" → "H"）。
fn element_from_label(label: &str) -> String {
    let alpha: String = label.chars().take_while(|c| c.is_alphabetic()).collect();
    capitalize_element(&alpha)
}

/// 清理 `_atom_site_type_symbol`：去除电荷标注（如 "Fe2+" → "Fe"），
/// 处理 "?" / "."。
fn clean_element(s: &str) -> String {
    match s.trim() {
        "?" | "." | "" => "X".to_string(),
        s => {
            let alpha: String = s.chars().take_while(|c| c.is_alphabetic()).collect();
            capitalize_element(&alpha)
        }
    }
}

fn capitalize_element(s: &str) -> String {
    if s.is_empty() { return "X".to_string(); }
    let mut it = s.chars();
    match it.next() {
        None => "X".to_string(),
        Some(first) => {
            let rest: String = it.collect::<String>().to_lowercase();
            format!("{}{rest}", first.to_uppercase())
        }
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // P1 结构：简单立方，两个原子，分数坐标
    const P1_CIF: &str = "\
data_test_p1
_cell_length_a   5.000000
_cell_length_b   5.000000
_cell_length_c   5.000000
_cell_angle_alpha   90.000
_cell_angle_beta    90.000
_cell_angle_gamma   90.000
_space_group_name_H-M  'P 1'
_space_group_IT_number  1
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Fe1  Fe  0.00000  0.00000  0.00000
Fe2  Fe  0.50000  0.50000  0.50000
";

    // 体心立方 Fe：1 个不等价原子 + 体心平移对称操作 → 扩展为 2 个原子
    const BCC_CIF: &str = "\
data_bcc_Fe
_cell_length_a   2.870000
_cell_length_b   2.870000
_cell_length_c   2.870000
_cell_angle_alpha   90.000
_cell_angle_beta    90.000
_cell_angle_gamma   90.000
loop_
_space_group_symop_operation_xyz
'x,y,z'
'x+1/2,y+1/2,z+1/2'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Fe1  Fe  0.00000  0.00000  0.00000
";

    // 带不确定度的晶格常数
    const UNCERTAINTY_CIF: &str = "\
data_unc
_cell_length_a   5.431(2)
_cell_length_b   5.431(2)
_cell_length_c   5.431(2)
_cell_angle_alpha   90.00(1)
_cell_angle_beta    90.00(1)
_cell_angle_gamma   90.00(1)
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1  Si  0.0  0.0  0.0
";

    fn write_tmp(name: &str, content: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(name);
        std::fs::write(&path, content).unwrap();
        path
    }

    #[test]
    fn test_p1_two_atoms() {
        let path = write_tmp("test_p1.cif", P1_CIF);
        let traj = read_cif(path.to_str().unwrap()).unwrap();
        assert_eq!(traj.n_frames(), 1);
        let frame = traj.first().unwrap();
        assert_eq!(frame.n_atoms(), 2);
        assert!(frame.is_periodic());
        let [a, b, c] = frame.cell.as_ref().unwrap().lengths();
        assert!((a - 5.0).abs() < 1e-4);
        assert!((b - 5.0).abs() < 1e-4);
        assert!((c - 5.0).abs() < 1e-4);
    }

    #[test]
    fn test_bcc_symmetry_expansion() {
        let path = write_tmp("test_bcc.cif", BCC_CIF);
        let traj = read_cif(path.to_str().unwrap()).unwrap();
        let frame = traj.first().unwrap();
        // 1 asymmetric atom × 2 symops = 2 Fe atoms
        assert_eq!(frame.n_atoms(), 2);
        assert_eq!(frame.atom(0).element, "Fe");
        assert_eq!(frame.atom(1).element, "Fe");
    }

    #[test]
    fn test_uncertainty_stripped() {
        let path = write_tmp("test_unc.cif", UNCERTAINTY_CIF);
        let traj = read_cif(path.to_str().unwrap()).unwrap();
        let [a, ..] = traj.first().unwrap().cell.as_ref().unwrap().lengths();
        assert!((a - 5.431).abs() < 1e-3);
    }

    #[test]
    fn test_parse_symop_identity() {
        let op = parse_symop("x,y,z").unwrap();
        let p = Vector3::new(0.3, 0.5, 0.7);
        assert!((op.apply(p) - p).norm() < 1e-10);
    }

    #[test]
    fn test_parse_symop_inversion() {
        let op = parse_symop("-x,-y,-z").unwrap();
        let p = Vector3::new(0.3, 0.5, 0.7);
        let q = op.apply(p);
        assert!((q - Vector3::new(-0.3, -0.5, -0.7)).norm() < 1e-10);
    }

    #[test]
    fn test_parse_symop_translation() {
        let op = parse_symop("x+1/2,y+1/2,z").unwrap();
        let p = Vector3::new(0.0, 0.0, 0.0);
        let q = op.apply(p);
        assert!((q.x - 0.5).abs() < 1e-10);
        assert!((q.y - 0.5).abs() < 1e-10);
        assert!((q.z - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_element_from_label() {
        assert_eq!(element_from_label("Fe1"), "Fe");
        assert_eq!(element_from_label("H12"), "H");
        assert_eq!(element_from_label("Ca2a"), "Ca");
        assert_eq!(element_from_label("O"), "O");
    }
}
