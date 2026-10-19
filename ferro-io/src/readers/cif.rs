use std::path::Path;
use std::collections::HashMap;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{bail, ensure, Context, Result};

use super::spacegroup;

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
    let idx_occ = col(&["_atom_site_occupancy"]);

    ensure!(idx_sym.is_some() || idx_lbl.is_some(), "no element/label column");
    ensure!(
        (idx_fx.is_some() && idx_fy.is_some() && idx_fz.is_some())
            || (idx_cx.is_some() && idx_cy.is_some() && idx_cz.is_some()),
        "no position columns (need fract_x/y/z or cartn_x/y/z)"
    );

    let mut sites = AtomSites { elements: vec![], labels: vec![], frac: vec![] };
    // 部分占位的位点：(行号, 标签, 占位度)。Frame 只能装整原子，混占位（同一位点
    // Mg 0.5 / Fe 0.5）照写会变成两个重叠原子，空位（O 0.9）会多出 10% 的原子 ——
    // 密度、配位、喂给 QC 的结构都错，且没有任何症状。有序化怎么做是用户的物理判断
    let mut partial: Vec<(usize, String, f64)> = Vec::new();

    for (k, row) in rows.iter().enumerate() {
        let get = |i: usize| row.get(i).map(|s| s.as_str()).unwrap_or("?");
        let name = idx_lbl.or(idx_sym).map(|i| get(i).to_string()).unwrap_or_default();
        // 坐标与占位度的坏值报错并点名位点 —— 以前坐标 `unwrap_or(0.0)`，`?` 或写坏的数
        // 把原子静默放到原点（与审查 L6 同一类）
        let num = |i: usize, tag: &str| -> Result<f64> {
            parse_cif_float(get(i)).with_context(|| {
                format!("site {} '{name}': {tag} = {:?} is not a number", k + 1, get(i))
            })
        };

        // 缺列或写 ?/. 时取 CIF 字典的默认值 1（写坏的数仍报错）；容差只吸收 0.9995 这类舍入
        let occ = match idx_occ {
            Some(i) if !matches!(get(i).trim(), "?" | ".") => num(i, "_atom_site_occupancy")?,
            _ => 1.0,
        };
        if occ < 1.0 - 1e-3 {
            partial.push((k + 1, name.clone(), occ));
        }

        let element = idx_sym
            .map(|i| clean_element(get(i)))
            .or_else(|| idx_lbl.map(|i| element_from_label(get(i))))
            .unwrap_or_else(|| "X".to_string());

        let label = idx_lbl.map(|i| get(i).to_string());

        let frac = if let (Some(ix), Some(iy), Some(iz)) = (idx_fx, idx_fy, idx_fz) {
            Vector3::new(
                num(ix, "_atom_site_fract_x")?,
                num(iy, "_atom_site_fract_y")?,
                num(iz, "_atom_site_fract_z")?,
            )
        } else if let (Some(ix), Some(iy), Some(iz)) = (idx_cx, idx_cy, idx_cz) {
            // Cartesian → fractional. Column presence is guaranteed by the
            // ensure! above; checked locally so the invariant is self-contained
            // and a future edit to ensure! cannot silently introduce a panic.
            let cart = Vector3::new(
                num(ix, "_atom_site_Cartn_x")?,
                num(iy, "_atom_site_Cartn_y")?,
                num(iz, "_atom_site_Cartn_z")?,
            );
            cell.cartesian_to_fractional(cart)?
        } else {
            bail!("atom site has neither fractional nor Cartesian coordinate columns");
        };

        sites.elements.push(element);
        sites.labels.push(label);
        sites.frac.push(frac);
    }

    if !partial.is_empty() {
        let list: Vec<String> = partial.iter()
            .map(|(k, name, occ)| format!("site {k} '{name}' occupancy {occ}"))
            .collect();
        bail!(
            "{} atom site(s) are partially occupied ({}); a structure holds whole atoms only, \
             so they would become overlapping atoms or extra atoms. Order the structure first \
             (pick one species per site, build a supercell for vacancies) and write occupancy 1",
            partial.len(),
            list.join(", ")
        );
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

/// 对称操作的来源按可靠程度排：CIF 显式列出的操作 → 空间群符号查表（Hall / H-M /
/// IT 号，见 `spacegroup.rs`）→ 都没写才当 P1。以前没有操作循环就一律当 P1，只写
/// `F m -3 m` 的 NaCl 读成 2 个原子、密度 0.54 g/cm³（2026-10-03 审查 A1）
fn collect_symops(block: &CifBlock, cell: &Cell) -> Result<Vec<SymOp>> {
    // Tag names used by old and new CIF dictionaries (CIF2 uses a dot)
    const TAGS: &[&str] = &[
        "_space_group_symop_operation_xyz",
        "_space_group_symop.operation_xyz",
        "_symmetry_equiv_pos_as_xyz",
    ];
    // 写坏的操作报错并点名 —— 以前 `filter_map(.ok())` 静默丢掉，全坏时退回 P1
    let parse = |s: &str| parse_symop(s).with_context(|| format!("bad symmetry operation '{s}'"));

    // 1. Try loop
    for (headers, rows) in &block.loops {
        for tag in TAGS {
            if let Some(col) = headers.iter().position(|h| h == tag) {
                let ops = rows.iter()
                    .filter_map(|row| row.get(col))
                    .map(|s| parse(s))
                    .collect::<Result<Vec<SymOp>>>()?;
                if !ops.is_empty() { return Ok(ops); }
            }
        }
    }

    // 2. Try single tag
    for tag in TAGS {
        if let Some(val) = block.singles.get(*tag) {
            return Ok(vec![parse(val)?]);
        }
    }

    // 3. Space-group symbol / number
    let single = |tags: &[&str]| -> Option<&str> {
        tags.iter()
            .filter_map(|t| block.singles.get(*t))
            .map(|v| v.trim())
            .find(|v| !matches!(*v, "" | "?" | "."))
    };
    let hall = single(&["_space_group_name_hall", "_space_group.name_hall",
                        "_symmetry_space_group_name_hall"]);
    let hm = single(&["_space_group_name_h-m_alt", "_space_group.name_h-m_alt",
                      "_symmetry_space_group_name_h-m"]);
    let number = single(&["_space_group_it_number", "_space_group.it_number",
                          "_symmetry_int_tables_number"])
        .map(|v| v.parse::<u32>().with_context(|| format!("bad space-group number '{v}'")))
        .transpose()?;
    if hall.is_none() && hm.is_none() && number.is_none() {
        // 4. No symmetry information at all: P1
        return Ok(vec![SymOp::identity()]);
    }
    let metric = cell.matrix * cell.matrix.transpose();
    let ops = spacegroup::lookup(&spacegroup::Query { hall, hm, number }, &metric)
        .context("no symmetry operations listed; resolving the space group from its symbol failed")?;
    Ok(ops.into_iter().map(|(rot, trans)| SymOp { rot, trans }).collect())
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
                apply_term(&term, sign, row, rot, trans)?;
                term.clear();
                sign = if ch == '+' { 1.0 } else { -1.0 };
            }
            c => term.push(c),
        }
    }
    apply_term(&term, sign, row, rot, trans)
}

/// 认不出的项报错 —— 以前静默丢掉，`2x`、`x*1` 之类让该行旋转部分少一项，
/// 展开出的原子落在错误位置
fn apply_term(term: &str, sign: f64, row: usize, rot: &mut Matrix3<f64>, trans: &mut Vector3<f64>) -> Result<()> {
    match term {
        "x" => rot[(row, 0)] += sign,
        "y" => rot[(row, 1)] += sign,
        "z" => rot[(row, 2)] += sign,
        "" => {}
        s => match parse_fraction(s) {
            Some(v) => trans[row] += sign * v,
            None => bail!("unrecognised term {s:?}"),
        },
    }
    Ok(())
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

fn expand_by_symmetry(sites: AtomSites, symops: &[SymOp]) -> Result<AtomSites> {
    // 判重口径同 ASE `Spacegroup.equivalent_sites`（审查 A-7，用户 2026-10-11 定）：分数坐标
    // 逐分量取最小镜像，|Δ| < 1e-3 即同一点。以前是欧氏距离 < 0.01，且对全部同元素原子判重，
    // 两个相距不到 0.01 的不等价位点后一个被整体吞掉
    const SYMPREC: f64 = 1e-3;
    let same = |a: Vector3<f64>, b: Vector3<f64>| {
        let d = a - b;
        [d.x, d.y, d.z].iter().all(|c| (c - c.round()).abs() < SYMPREC)
    };

    let mut out = AtomSites { elements: vec![], labels: vec![], frac: vec![] };

    for ((element, label), &frac) in sites.elements.iter()
        .zip(sites.labels.iter())
        .zip(sites.frac.iter())
    {
        // 只在本位点的像之间去重：特殊位置上的多个像是同一个原子
        let mut orbit: Vec<Vector3<f64>> = Vec::new();
        for op in symops {
            let p = op.apply(frac).map(|c| c.rem_euclid(1.0));
            if !orbit.iter().any(|&q| same(p, q)) {
                orbit.push(p);
            }
        }
        // 与先前位点的像重合 = 同一位置写了两个非对称位点：不替用户挑一个
        // （部分占位已在 parse_atom_sites 拒收，这里剩下的是 CIF 写错了）
        for &p in &orbit {
            if let Some(k) = out.frac.iter().position(|&q| same(p, q)) {
                bail!(
                    "CIF sites {} and {} land on the same position (within 1e-3 \
                     fractional) after symmetry expansion",
                    label.as_deref().unwrap_or(element),
                    out.labels[k].as_deref().unwrap_or(&out.elements[k])
                );
            }
        }
        for p in orbit {
            out.elements.push(element.clone());
            out.labels.push(label.clone());
            out.frac.push(p);
        }
    }
    Ok(out)
}

// ─── Block → Frame ────────────────────────────────────────────────────────────

fn block_to_frame(block: &CifBlock) -> Result<Frame> {
    let cell = parse_cell(block)?;
    let sites = parse_atom_sites(block, &cell)?;

    let symops = collect_symops(block, &cell)?;
    let sites = if symops.len() > 1 {
        expand_by_symmetry(sites, &symops)?
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
pub fn read_cif(path: &Path) -> Result<Trajectory> {
    let path_ = path.display();
    let input = super::util::read_text(path)?;
    let tokens = tokenize(&input);
    let blocks = parse_blocks(&tokens);

    ensure!(!blocks.is_empty(), "no data_ blocks found in {path_}");

    // 没有 `_atom_site_` 循环的 block 不是结构（Acta Cryst / COD 常见的 `data_global`
    // 只放作者与期刊信息）—— 跳过，同 ASE / pymatgen；以前整份文件因它报错
    let (structs, skipped): (Vec<&CifBlock>, Vec<&CifBlock>) = blocks.iter().partition(|b| {
        b.loops.iter().any(|(hdrs, _)| hdrs.iter().any(|h| h.starts_with("_atom_site_")))
    });
    if structs.is_empty() {
        let names: Vec<String> = skipped.iter().map(|b| format!("data_{}", b.name)).collect();
        bail!("no data block in {path_} has an _atom_site_ loop ({})", names.join(", "));
    }

    let mut traj = Trajectory::new();
    if let Some(name) = structs.first().map(|b| b.name.clone()) {
        if !name.is_empty() { traj.metadata.title = Some(name); }
    }

    for block in structs {
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

    use crate::testutil::write_tmp;

    // 混占位：同一位点 Mg 0.5 / Fe 0.5（ICSD 无序结构的典型写法）
    const MIXED_OCC_CIF: &str = "\
data_mixed
_cell_length_a   4.2
_cell_length_b   4.2
_cell_length_c   4.2
_cell_angle_alpha   90
_cell_angle_beta    90
_cell_angle_gamma   90
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Mg1  Mg  0.0  0.0  0.0  0.5
Fe1  Fe  0.0  0.0  0.0  0.5
O1   O   0.5  0.5  0.5  1.0(0)
";

    #[test]
    fn test_bad_coordinate_or_occupancy_is_an_error_naming_the_site() {
        let err = |text: String| format!("{:#}", read_cif(&write_tmp("test_badcoord.cif", &text)).unwrap_err());
        // 未知坐标 `?`、写坏的数：以前原子静默落在原点
        for (bad, tag) in [("0.0  ?  0.0  0.5\nFe1", "fract_y"), ("0.0  0.0  0.O  0.5\nFe1", "fract_z")] {
            let text = MIXED_OCC_CIF.replacen("0.0  0.0  0.0  0.5\nFe1", bad, 1);
            let msg = err(text);
            assert!(msg.contains("site 1 'Mg1'") && msg.contains(tag), "应点名位点与列，实际：{msg}");
        }
        // 写坏的占位度报错；`?` / `.` 取默认值 1
        let msg = err(MIXED_OCC_CIF.replace("1.0(0)", "1.O"));
        assert!(msg.contains("site 3 'O1'") && msg.contains("occupancy"), "实际：{msg}");
        let unknown = MIXED_OCC_CIF.replace("0.0  0.5\nFe1  Fe  0.0  0.0  0.0  0.5\n", "0.0  ?\n");
        assert_eq!(read_cif(&write_tmp("test_occ_unknown.cif", &unknown)).unwrap().first().unwrap().n_atoms(), 2);
    }

    #[test]
    fn test_partial_occupancy_is_an_error() {
        let err = read_cif(&write_tmp("test_mixed.cif", MIXED_OCC_CIF)).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("2 atom site(s)") && msg.contains("'Mg1' occupancy 0.5")
                && msg.contains("'Fe1'") && !msg.contains("'O1'"),
                "应点名两个半占位位点、不点名满占位的 O1，实际：{msg}");
        // 全满占位（含 1.0(0) 的不确定度写法）照常读
        let full = MIXED_OCC_CIF.replace("0.0  0.5\nFe1  Fe  0.0  0.0  0.0  0.5\n", "0.0  1\n");
        let traj = read_cif(&write_tmp("test_full_occ.cif", &full)).unwrap();
        assert_eq!(traj.first().unwrap().n_atoms(), 2);
    }

    #[test]
    fn test_blocks_without_atom_sites_are_skipped() {
        // Acta Cryst / COD 的 data_global：只有出版信息，含一个非 atom_site 的 loop
        const GLOBAL: &str = "data_global\n_journal_name_full 'Acta Cryst.'\n\
            loop_\n_publ_author_name\n'Doe, J.'\n'Roe, R.'\n\n";
        let traj = read_cif(&write_tmp("test_global.cif", &format!("{GLOBAL}{P1_CIF}"))).unwrap();
        assert_eq!(traj.n_frames(), 1);
        assert_eq!(traj.metadata.title.as_deref(), Some("test_p1"), "source 取第一个结构 block");

        let err = read_cif(&write_tmp("test_global_only.cif", GLOBAL)).unwrap_err();
        assert!(format!("{err:#}").contains("data_global"), "{err:#}");
    }

    #[test]
    fn test_p1_two_atoms() {
        let path = write_tmp("test_p1.cif", P1_CIF);
        let traj = read_cif(&path).unwrap();
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
        let traj = read_cif(&path).unwrap();
        let frame = traj.first().unwrap();
        // 1 asymmetric atom × 2 symops = 2 Fe atoms
        assert_eq!(frame.n_atoms(), 2);
        assert_eq!(frame.atom(0).element, "Fe");
        assert_eq!(frame.atom(1).element, "Fe");
    }

    // 审查 A-7：判重容差 0.01 → ASE 的 1e-3（分数坐标逐分量），且只在同一非对称位点的像之间判重
    fn inversion_cif(sites: &str) -> String {
        BCC_CIF.replace("'x+1/2,y+1/2,z+1/2'", "'-x,-y,-z'")
               .replace("Fe1  Fe  0.00000  0.00000  0.00000\n", sites)
    }

    #[test]
    fn test_images_merge_only_within_ase_symprec() {
        let n = |name: &str, sites: &str| {
            read_cif(&write_tmp(name, &inversion_cif(sites))).unwrap().frames[0].n_atoms()
        };
        // 像相距 0.0008 < 1e-3：同一个原子
        assert_eq!(n("a7_close.cif", "O1 O 0.0004 0.5 0.5\n"), 1);
        // 像相距 0.006：以前 0.01 的容差吞成一个，实为两个原子
        assert_eq!(n("a7_apart.cif", "O1 O 0.003 0.5 0.5\n"), 2);
        // 两个不等价位点相距 0.005：以前后一个整体被吞
        assert_eq!(n("a7_two_sites.cif", "O1 O 0.100 0.1 0.1\nO2 O 0.105 0.1 0.1\n"), 4);
    }

    #[test]
    fn test_two_sites_on_one_position_is_an_error() {
        let text = inversion_cif("O1 O 0.1 0.2 0.3\nO2 O -0.1 -0.2 -0.3\n");
        let err = format!("{:#}", read_cif(&write_tmp("a7_same.cif", &text)).unwrap_err());
        assert!(err.contains("O1") && err.contains("O2"), "应点名两个位点：{err}");
    }

    // 只写空间群符号、不列对称操作的 NaCl（审查 A1 的复现用例）
    const NACL_HM_ONLY_CIF: &str = "\
data_NaCl
_cell_length_a   5.64
_cell_length_b   5.64
_cell_length_c   5.64
_cell_angle_alpha   90
_cell_angle_beta    90
_cell_angle_gamma   90
_symmetry_space_group_name_H-M  'F m -3 m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na1  Na  0.0  0.0  0.0
Cl1  Cl  0.5  0.5  0.5
";

    fn n_atoms_of(name: &str, text: &str) -> Result<usize> {
        Ok(read_cif(&write_tmp(name, text))?.first().unwrap().n_atoms())
    }

    #[test]
    fn test_symbol_only_cif_is_expanded() {
        // 以前当 P1 读成 2 个原子；ASE 读出 Cl4Na4
        assert_eq!(n_atoms_of("test_nacl_hm.cif", NACL_HM_ONLY_CIF).unwrap(), 8);
        let by_number = NACL_HM_ONLY_CIF.replace(
            "_symmetry_space_group_name_H-M  'F m -3 m'", "_space_group_IT_number  225");
        assert_eq!(n_atoms_of("test_nacl_it.cif", &by_number).unwrap(), 8);
        let by_hall = NACL_HM_ONLY_CIF.replace(
            "_symmetry_space_group_name_H-M  'F m -3 m'", "_space_group_name_Hall  '-F 4 2 3'");
        assert_eq!(n_atoms_of("test_nacl_hall.cif", &by_hall).unwrap(), 8);
    }

    #[test]
    fn test_symbol_only_cif_errors_instead_of_guessing() {
        let err = |name: &str, text: String| format!("{:#}", read_cif(&write_tmp(name, &text)).unwrap_err());
        // 原点选择不定：报错并给修法，不默认取 1
        let diamond = NACL_HM_ONLY_CIF.replace("F m -3 m", "F d -3 m");
        let msg = err("test_fd3m.cif", diamond.clone());
        assert!(msg.contains("ambiguous") && msg.contains(":2"), "{msg}");
        // 同样的坐标，原点 1 下是 8a+8b、原点 2 下是 16c+16d —— 选错原点原子数都会变
        // （与 ASE Spacegroup(227, setting=1/2).equivalent_sites 一致）
        assert_eq!(n_atoms_of("test_fd3m_1.cif", &diamond.replace("F d -3 m", "F d -3 m :1")).unwrap(), 16);
        assert_eq!(n_atoms_of("test_fd3m_2.cif", &diamond.replace("F d -3 m", "F d -3 m :2")).unwrap(), 32);
        // 符号与编号矛盾
        let msg = err("test_conflict.cif", NACL_HM_ONLY_CIF.replace(
            "_atom_site_label", "_space_group_IT_number 221\nloop_\n_atom_site_label").replacen("loop_\n_space", "_space", 1));
        assert!(msg.contains("contradicts"), "{msg}");
        // 不认识的符号
        let msg = err("test_unknown_hm.cif", NACL_HM_ONLY_CIF.replace("F m -3 m", "F m -3 x"));
        assert!(msg.contains("unknown Hermann-Mauguin"), "{msg}");
    }

    #[test]
    fn test_cif2_dotted_symop_tag_and_bad_symop() {
        let cif2 = BCC_CIF.replace("_space_group_symop_operation_xyz", "_space_group_symop.operation_xyz");
        assert_eq!(n_atoms_of("test_cif2_symop.cif", &cif2).unwrap(), 2);
        // 写坏的操作报错并点名，以前静默丢弃
        let bad = BCC_CIF.replace("'x+1/2,y+1/2,z+1/2'", "'x+1/2,y+1/2'");
        let msg = format!("{:#}", read_cif(&write_tmp("test_bad_symop.cif", &bad)).unwrap_err());
        assert!(msg.contains("bad symmetry operation 'x+1/2,y+1/2'"), "{msg}");
    }

    #[test]
    fn test_uncertainty_stripped() {
        let path = write_tmp("test_unc.cif", UNCERTAINTY_CIF);
        let traj = read_cif(&path).unwrap();
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

    /// A term the parser does not know used to be dropped, leaving the row short.
    #[test]
    fn test_parse_symop_unknown_term_is_an_error() {
        for bad in ["2x,y,z", "x,y*1,z", "x,y,1/0", "x,y,w"] {
            let err = parse_symop(bad).unwrap_err().to_string();
            assert!(err.contains("unrecognised term"), "{bad}: {err}");
        }
        // 小数平移、平移写在前面仍认
        let q = parse_symop("1/2+x,y-0.25,-z").unwrap().apply(Vector3::new(0.1, 0.5, 0.2));
        assert!((q - Vector3::new(0.6, 0.25, -0.2)).norm() < 1e-12);
    }
}
