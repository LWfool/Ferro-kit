//! Mixed-compound box construction: estimate size, random placement, soft-core relaxation.

use std::collections::HashMap;

use nalgebra::Vector3;
use rand::RngExt;

use ferro_core::atom::Atom;
use ferro_core::cell::Cell;
use ferro_core::data::{compounds, elements};
use ferro_core::error::{ChemError, Result};
use ferro_core::frame::Frame;

/// A component in a mixed-compound system.
pub struct Component {
    /// Compound name or formula, looked up in the COMPOUNDS database.
    pub compound: String,
    /// Number of molecules of this compound.
    pub n_molecules: usize,
}

/// Read an optional multiplicity at `*i`, advancing past it. Absent digits mean 1.
///
/// Rejects an explicit zero: `Ca0` and `(PO4)0` are almost certainly typos, and silently
/// dropping the group would leave the caller with a box missing a whole component.
fn read_count(bytes: &[u8], i: &mut usize, formula: &str) -> Result<usize> {
    let start = *i;
    while *i < bytes.len() && bytes[*i].is_ascii_digit() {
        *i += 1;
    }
    if start == *i {
        return Ok(1);
    }
    let text = &formula[start..*i];
    let count: usize = text.parse().map_err(|_| {
        ChemError::ParseError(format!("invalid count '{}' in '{}'", text, formula))
    })?;
    if count == 0 {
        return Err(ChemError::ParseError(format!(
            "zero multiplicity at position {} in '{}'",
            start, formula
        )));
    }
    Ok(count)
}

/// Parse a chemical formula into element counts.
///
/// Supports nested groups in round or square brackets, which must be balanced and
/// correctly paired: `Ca3(PO4)2`, `(NH4)2SO4`, `K4[Fe(CN)6]`, `Mg[Al(OH)4]2`.
/// Elements are one uppercase letter optionally followed by lowercase letters; a count
/// following an element or a closing bracket defaults to 1 when omitted.
///
/// Hydrate dot notation (`CuSO4·5H2O`) is *not* accepted — write it out as `CuSO9H10`
/// only if that is genuinely what you mean, or list the water as a separate component,
/// which is what the box builder wants anyway.
fn parse_formula(formula: &str) -> Result<HashMap<String, usize>> {
    if formula.is_empty() {
        return Err(ChemError::ParseError("empty formula".into()));
    }
    let bytes = formula.as_bytes();
    let n = bytes.len();
    let mut i = 0;

    // 每层一个计数表；栈底是整个式子，每遇到左括号压入一层
    let mut stack: Vec<HashMap<String, usize>> = vec![HashMap::new()];
    // 与 stack 平行，记录每层是被哪个字符打开的，用来拒绝 `[..)` 这类错配
    let mut open: Vec<(u8, usize)> = Vec::new();

    while i < n {
        let b = bytes[i];
        if b == b'(' || b == b'[' {
            open.push((b, i));
            stack.push(HashMap::new());
            i += 1;
        } else if b == b')' || b == b']' {
            let Some((opener, open_pos)) = open.pop() else {
                return Err(ChemError::ParseError(format!(
                    "unmatched '{}' at position {} in '{}'",
                    b as char, i, formula
                )));
            };
            let expected = if opener == b'(' { b')' } else { b']' };
            if b != expected {
                return Err(ChemError::ParseError(format!(
                    "'{}' at position {} closed by '{}' at position {} in '{}'",
                    opener as char, open_pos, b as char, i, formula
                )));
            }
            let group = stack.pop().expect("stack depth tracks open brackets");
            if group.is_empty() {
                return Err(ChemError::ParseError(format!(
                    "empty group at position {} in '{}'",
                    open_pos, formula
                )));
            }
            i += 1;
            let mult = read_count(bytes, &mut i, formula)?;
            let outer = stack.last_mut().expect("bottom frame is never popped");
            for (elem, c) in group {
                *outer.entry(elem).or_insert(0) += c * mult;
            }
        } else if b.is_ascii_uppercase() {
            let start = i;
            i += 1;
            while i < n && bytes[i].is_ascii_lowercase() {
                i += 1;
            }
            let elem = formula[start..i].to_string();
            let count = read_count(bytes, &mut i, formula)?;
            *stack.last_mut().expect("bottom frame is never popped")
                .entry(elem).or_insert(0) += count;
        } else {
            return Err(ChemError::ParseError(format!(
                "unexpected character '{}' at position {} in '{}'",
                b as char, i, formula
            )));
        }
    }

    if let Some((opener, pos)) = open.pop() {
        return Err(ChemError::ParseError(format!(
            "unclosed '{}' at position {} in '{}'",
            opener as char, pos, formula
        )));
    }
    let counts = stack.pop().expect("bottom frame is never popped");
    if counts.is_empty() {
        return Err(ChemError::ParseError(format!("no elements in '{}'", formula)));
    }
    Ok(counts)
}

/// Resolve a component spec into its element counts and molar mass (g/mol).
///
/// The COMPOUNDS database is tried first, by name or formula, so the curated
/// `molecular_mass` values keep being used verbatim for known compounds. Anything not in
/// the database is parsed as a chemical formula and its mass summed from atomic weights —
/// that fallback is what lets callers ask for `Ca3(PO4)2` without the database knowing it.
fn resolve_component(spec: &str) -> Result<(HashMap<String, usize>, f64)> {
    if let Some(cd) = compounds::find(spec) {
        return Ok((parse_formula(cd.formula)?, cd.molecular_mass));
    }
    let counts = parse_formula(spec).map_err(|e| {
        ChemError::ValidationError(format!(
            "'{}' is neither a known compound nor a valid formula ({})",
            spec, e
        ))
    })?;
    let mut mass = 0.0_f64;
    for (elem, count) in &counts {
        let ed = elements::by_symbol(elem).ok_or_else(|| {
            ChemError::ValidationError(format!("unknown element '{}' in formula '{}'", elem, spec))
        })?;
        mass += ed.atomic_mass * *count as f64;
    }
    Ok((counts, mass))
}

/// Estimate the cubic box edge length (Å) for a mixed-compound system.
///
/// - `components`: list of (compound, n_molecules).
/// - `density`: target total density in g/cm³.
///
/// V = Σ(n_i × M_i) / (ρ × N_A), L = V^(1/3).
pub fn estimate_box_length(components: &[Component], density: f64) -> Result<f64> {
    if density <= 0.0 {
        return Err(ChemError::ValidationError(format!(
            "density must be > 0, got {}",
            density
        )));
    }

    let mut total_mass = 0.0_f64; // g/mol
    for comp in components {
        if comp.n_molecules == 0 {
            continue;
        }
        let (_, molar_mass) = resolve_component(&comp.compound)?;
        total_mass += comp.n_molecules as f64 * molar_mass;
    }

    if total_mass == 0.0 {
        return Err(ChemError::ValidationError(
            "no molecules specified".into(),
        ));
    }

    // V (cm³) = total_mass (g/mol) / (density (g/cm³) × N_A (1/mol))
    // V (Å³) = V (cm³) × 10^24
    const N_A: f64 = 6.022_140_76e23;
    let volume_cm3 = total_mass / (density * N_A);
    let volume_ang3 = volume_cm3 * 1e24;
    Ok(volume_ang3.cbrt())
}

/// Build a mixed-compound box.
///
/// 1. Parse compound formulas to get element stoichiometry.
/// 2. Estimate cubic box size from density.
/// 3. Place atoms randomly.
/// 4. Soft-core relaxation to remove overlaps.
///
/// - `min_dist`: minimum allowed distance between atoms (Å). Suggested: 1.5.
/// - `relax_steps`: number of relaxation steps. Suggested: 100.
pub fn build_box(
    components: &[Component],
    density: f64,
    min_dist: f64,
    relax_steps: usize,
) -> Result<Frame> {
    if min_dist <= 0.0 {
        return Err(ChemError::ValidationError(format!(
            "min_dist must be > 0, got {}",
            min_dist
        )));
    }

    // ── 步骤 1：解析化学式，统计各元素原子数 ─────────────────────────────────
    let mut element_counts: HashMap<String, usize> = HashMap::new();
    for comp in components {
        if comp.n_molecules == 0 {
            continue;
        }
        let (formula_counts, _) = resolve_component(&comp.compound)?;
        for (elem, count) in formula_counts {
            *element_counts.entry(elem).or_insert(0) += count * comp.n_molecules;
        }
    }

    if element_counts.is_empty() {
        return Err(ChemError::ValidationError(
            "no atoms to place".into(),
        ));
    }

    // ── 步骤 2：估算盒子尺寸 ─────────────────────────────────────────────────
    let box_len = estimate_box_length(components, density)?;
    let cell = Cell::from_lengths_angles(box_len, box_len, box_len, 90.0, 90.0, 90.0)?;

    // ── 步骤 3：随机放置原子 ─────────────────────────────────────────────────
    let mut rng = rand::rng();
    let mut atoms: Vec<Atom> = Vec::new();
    for (elem, count) in &element_counts {
        for _ in 0..*count {
            let pos = Vector3::new(
                rng.random_range(0.0..box_len),
                rng.random_range(0.0..box_len),
                rng.random_range(0.0..box_len),
            );
            atoms.push(Atom::new(elem, pos));
        }
    }

    // ── 步骤 4：软核弛豫 ─────────────────────────────────────────────────────
    relax(&mut atoms, box_len, min_dist, relax_steps);

    Ok(Frame {
        atoms,
        cell: Some(cell),
        pbc: [true; 3],
        charge: 0,
        multiplicity: 1,
        bonds: None,
        energy: None,
        forces: None,
        stress: None,
        velocities: None,
    })
}

/// Cell list 辅助结构：将空间划分为 bin，加速近邻搜索。
struct CellList {
    n_bins: usize,
    inv_bin_size: f64,
    bins: Vec<Vec<usize>>,
}

impl CellList {
    fn new(box_len: f64, cutoff: f64) -> Self {
        // 每个 bin 至少为 cutoff 大小
        let n_bins = (box_len / cutoff).floor().max(1.0) as usize;
        let inv_bin_size = n_bins as f64 / box_len;
        let total = n_bins * n_bins * n_bins;
        Self {
            n_bins,
            inv_bin_size,
            bins: vec![Vec::new(); total],
        }
    }

    fn bin_index(&self, pos: &Vector3<f64>) -> (usize, usize, usize) {
        let ix = ((pos.x * self.inv_bin_size).floor() as usize).min(self.n_bins - 1);
        let iy = ((pos.y * self.inv_bin_size).floor() as usize).min(self.n_bins - 1);
        let iz = ((pos.z * self.inv_bin_size).floor() as usize).min(self.n_bins - 1);
        (ix, iy, iz)
    }

    fn flat_index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix * self.n_bins * self.n_bins + iy * self.n_bins + iz
    }

    fn rebuild(&mut self, atoms: &[Atom]) {
        for bin in &mut self.bins {
            bin.clear();
        }
        for (i, atom) in atoms.iter().enumerate() {
            let (ix, iy, iz) = self.bin_index(&atom.position);
            let idx = self.flat_index(ix, iy, iz);
            self.bins[idx].push(i);
        }
    }

    /// 返回 atom_i 的近邻原子下标（排除自身）。
    fn neighbors(&self, atoms: &[Atom], i: usize) -> Vec<usize> {
        let (ix, iy, iz) = self.bin_index(&atoms[i].position);
        let mut result = Vec::new();
        // 检查 27 个相邻 bin；盒子过小（n_bins ≤ 2）时 rem 折叠会让同一 bin 被多次
        // 映射，用 seen 去重避免把近邻原子重复计入（否则弛豫斥力被成倍累加）。
        let mut seen: Vec<usize> = Vec::with_capacity(27);
        for dx in [self.n_bins - 1, 0, 1] {
            for dy in [self.n_bins - 1, 0, 1] {
                for dz in [self.n_bins - 1, 0, 1] {
                    let nx = (ix + dx) % self.n_bins;
                    let ny = (iy + dy) % self.n_bins;
                    let nz = (iz + dz) % self.n_bins;
                    let idx = self.flat_index(nx, ny, nz);
                    if seen.contains(&idx) { continue; }
                    seen.push(idx);
                    for &j in &self.bins[idx] {
                        if j != i {
                            result.push(j);
                        }
                    }
                }
            }
        }
        result
    }
}

/// 最小镜像位移（正交盒子）。
fn minimum_image_ortho(diff: &Vector3<f64>, box_len: f64) -> Vector3<f64> {
    Vector3::new(
        diff.x - box_len * (diff.x / box_len).round(),
        diff.y - box_len * (diff.y / box_len).round(),
        diff.z - box_len * (diff.z / box_len).round(),
    )
}

/// 软核弛豫：最速下降法消除原子重叠。
fn relax(atoms: &mut [Atom], box_len: f64, min_dist: f64, max_steps: usize) {
    if atoms.is_empty() {
        return;
    }
    let cutoff = min_dist * 1.5; // cell list 截断略大于 min_dist
    let step_size = min_dist * 0.05;
    let n = atoms.len();

    let mut cl = CellList::new(box_len, cutoff);

    for _step in 0..max_steps {
        cl.rebuild(atoms);

        // 计算每个原子的斥力
        let mut forces = vec![Vector3::zeros(); n];
        for i in 0..n {
            for j in cl.neighbors(atoms, i) {
                if j <= i {
                    continue; // 避免重复计算
                }
                let diff = minimum_image_ortho(&(atoms[j].position - atoms[i].position), box_len);
                let dist = diff.norm();
                if dist > 0.0 && dist < min_dist {
                    let mag = 1.0 - dist / min_dist;
                    let dir = diff / dist;
                    forces[i] -= dir * mag;
                    forces[j] += dir * mag;
                }
            }
        }

        // 更新位置
        let mut max_disp = 0.0_f64;
        for i in 0..n {
            let disp = forces[i] * step_size;
            let d = disp.norm();
            if d > max_disp {
                max_disp = d;
            }
            atoms[i].position += disp;
            // wrap 回盒子
            atoms[i].position.x = atoms[i].position.x.rem_euclid(box_len);
            atoms[i].position.y = atoms[i].position.y.rem_euclid(box_len);
            atoms[i].position.z = atoms[i].position.z.rem_euclid(box_len);
        }

        // 早停：位移足够小
        if max_disp < 1e-6 {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_formula_water() {
        let counts = parse_formula("H2O").unwrap();
        assert_eq!(counts["H"], 2);
        assert_eq!(counts["O"], 1);
    }

    #[test]
    fn test_parse_formula_p2o5() {
        let counts = parse_formula("P2O5").unwrap();
        assert_eq!(counts["P"], 2);
        assert_eq!(counts["O"], 5);
    }

    #[test]
    fn test_parse_formula_zno() {
        let counts = parse_formula("ZnO").unwrap();
        assert_eq!(counts["Zn"], 1);
        assert_eq!(counts["O"], 1);
    }

    #[test]
    fn test_parse_formula_ch3oh() {
        let counts = parse_formula("CH3OH").unwrap();
        assert_eq!(counts["C"], 1);
        assert_eq!(counts["H"], 4); // 3 + 1
        assert_eq!(counts["O"], 1);
    }

    #[test]
    fn test_parse_formula_ccl4() {
        let counts = parse_formula("CCl4").unwrap();
        assert_eq!(counts["C"], 1);
        assert_eq!(counts["Cl"], 4);
    }

    #[test]
    fn test_parse_formula_single_element() {
        let counts = parse_formula("Fe").unwrap();
        assert_eq!(counts["Fe"], 1);
    }

    #[test]
    fn test_parse_formula_empty() {
        // 空串不是合法化学式。旧实现返回 Ok(空表)，把报错推迟到 build_box 的
        // "no atoms to place"，错误信息离病因更远。
        assert!(parse_formula("").is_err());
    }

    // ── 括号 ────────────────────────────────────────────────────────────────

    #[test]
    fn test_parse_formula_parenthesised_group() {
        let counts = parse_formula("Ca3(PO4)2").unwrap();
        assert_eq!(counts["Ca"], 3);
        assert_eq!(counts["P"], 2);
        assert_eq!(counts["O"], 8);
        assert_eq!(counts.len(), 3);
    }

    #[test]
    fn test_parse_formula_leading_group() {
        let counts = parse_formula("(NH4)2SO4").unwrap();
        assert_eq!(counts["N"], 2);
        assert_eq!(counts["H"], 8);
        assert_eq!(counts["S"], 1);
        assert_eq!(counts["O"], 4);
    }

    #[test]
    fn test_parse_formula_nested_groups() {
        // 嵌套：内层 (CN)6 先归并进 [Fe...] 层，再整体乘 1
        let counts = parse_formula("K4[Fe(CN)6]").unwrap();
        assert_eq!(counts["K"], 4);
        assert_eq!(counts["Fe"], 1);
        assert_eq!(counts["C"], 6);
        assert_eq!(counts["N"], 6);
    }

    #[test]
    fn test_parse_formula_nested_with_outer_multiplier() {
        // 外层乘数必须传播到内层：Al 2、O 8、H 8
        let counts = parse_formula("Mg[Al(OH)4]2").unwrap();
        assert_eq!(counts["Mg"], 1);
        assert_eq!(counts["Al"], 2);
        assert_eq!(counts["O"], 8);
        assert_eq!(counts["H"], 8);
    }

    #[test]
    fn test_parse_formula_group_merges_with_outer_element() {
        // 同一元素既在组内又在组外时必须累加，而不是相互覆盖
        let counts = parse_formula("O(O2)3").unwrap();
        assert_eq!(counts["O"], 7);
        assert_eq!(counts.len(), 1);
    }

    #[test]
    fn test_parse_formula_bracket_errors() {
        for bad in ["Ca3(PO4", "Ca3PO4)2", "Ca3(PO4]2", "Ca3()2", "Ca0", "(PO4)0", "ca3"] {
            assert!(parse_formula(bad).is_err(), "'{bad}' should be rejected");
        }
    }

    // ── 数据库外公式回退 ─────────────────────────────────────────────────────

    #[test]
    fn test_resolve_known_compound_uses_curated_mass() {
        let (counts, mass) = resolve_component("water").unwrap();
        assert_eq!(counts["H"], 2);
        assert_eq!(counts["O"], 1);
        assert_eq!(mass, compounds::find("water").unwrap().molecular_mass);
    }

    #[test]
    fn test_resolve_unknown_formula_sums_atomic_masses() {
        // Ca3(PO4)2 不在化合物库中，走公式回退；分子量 310.18 g/mol
        let (counts, mass) = resolve_component("Ca3(PO4)2").unwrap();
        assert_eq!(counts["Ca"], 3);
        assert_eq!(counts["O"], 8);
        assert!((mass - 310.18).abs() < 0.5, "molar mass {mass} off expected 310.18");
    }

    #[test]
    fn test_resolve_rejects_unknown_element() {
        let err = resolve_component("Xx2O3").unwrap_err();
        assert!(format!("{err}").contains("unknown element"), "got: {err}");
    }

    #[test]
    fn test_build_box_accepts_formula_outside_database() {
        // 端到端：库外带括号公式一路走到建盒
        let comps = vec![Component { compound: "Ca3(PO4)2".to_string(), n_molecules: 8 }];
        let frame = build_box(&comps, 3.0, 1.5, 20).unwrap();
        assert_eq!(frame.atoms.iter().filter(|a| a.element == "Ca").count(), 24);
        assert_eq!(frame.atoms.iter().filter(|a| a.element == "P").count(), 16);
        assert_eq!(frame.atoms.iter().filter(|a| a.element == "O").count(), 64);
    }

    #[test]
    fn test_estimate_box_length_water() {
        // 1000 个水分子，密度 1.0 g/cm³
        let comps = vec![Component {
            compound: "water".into(),
            n_molecules: 1000,
        }];
        let box_len = estimate_box_length(&comps, 1.0).unwrap();
        // 1000 × 18.015 / (1.0 × 6.022e23) = 2.992e-20 cm³ = 29920 Å³
        // L = 29920^(1/3) ≈ 31.04 Å
        assert!(box_len > 30.0 && box_len < 32.0, "box_len = {}", box_len);
    }

    #[test]
    fn test_estimate_box_length_multi_component() {
        let comps = vec![
            Component {
                compound: "water".into(),
                n_molecules: 500,
            },
            Component {
                compound: "ethanol".into(),
                n_molecules: 500,
            },
        ];
        let box_len = estimate_box_length(&comps, 0.9).unwrap();
        // 总质量 = 500×18.015 + 500×46.069 = 32042 g/mol
        // V = 32042 / (0.9 × 6.022e23) = 5.91e-20 cm³ = 59100 Å³
        // L ≈ 38.97 Å
        assert!(box_len > 38.0 && box_len < 40.0, "box_len = {}", box_len);
    }

    #[test]
    fn test_estimate_box_length_zero_density() {
        let comps = vec![Component {
            compound: "water".into(),
            n_molecules: 100,
        }];
        assert!(estimate_box_length(&comps, 0.0).is_err());
    }

    #[test]
    fn test_estimate_box_length_unknown_compound() {
        let comps = vec![Component {
            compound: "unknown_xyz".into(),
            n_molecules: 100,
        }];
        assert!(estimate_box_length(&comps, 1.0).is_err());
    }

    #[test]
    fn test_build_box_atom_count() {
        let comps = vec![Component {
            compound: "water".into(),
            n_molecules: 10,
        }];
        let frame = build_box(&comps, 1.0, 1.5, 50).unwrap();
        // 10 个 H2O = 30 个原子
        assert_eq!(frame.atoms.len(), 30);
    }

    #[test]
    fn test_build_box_element_counts() {
        let comps = vec![Component {
            compound: "water".into(),
            n_molecules: 5,
        }];
        let frame = build_box(&comps, 1.0, 1.5, 50).unwrap();
        assert_eq!(frame.count_element("H"), 10);
        assert_eq!(frame.count_element("O"), 5);
    }

    #[test]
    fn test_build_box_periodic() {
        let comps = vec![Component {
            compound: "water".into(),
            n_molecules: 5,
        }];
        let frame = build_box(&comps, 1.0, 1.5, 50).unwrap();
        assert!(frame.is_periodic());
        assert_eq!(frame.pbc, [true; 3]);
    }

    #[test]
    fn test_build_box_atoms_within_cell() {
        let comps = vec![Component {
            compound: "water".into(),
            n_molecules: 20,
        }];
        let frame = build_box(&comps, 1.0, 1.5, 100).unwrap();
        let box_len = frame.cell.as_ref().unwrap().lengths()[0];
        for atom in &frame.atoms {
            assert!(atom.position.x >= 0.0 && atom.position.x < box_len);
            assert!(atom.position.y >= 0.0 && atom.position.y < box_len);
            assert!(atom.position.z >= 0.0 && atom.position.z < box_len);
        }
    }

    #[test]
    fn test_build_box_no_overlaps_after_relax() {
        let comps = vec![Component {
            compound: "water".into(),
            n_molecules: 20,
        }];
        let min_dist = 1.5;
        let frame = build_box(&comps, 0.5, min_dist, 200).unwrap();
        let box_len = frame.cell.as_ref().unwrap().lengths()[0];
        let n = frame.atoms.len();
        // 检查所有原子对的最小距离
        for i in 0..n {
            for j in (i + 1)..n {
                let diff = frame.atoms[j].position - frame.atoms[i].position;
                let mic = minimum_image_ortho(&diff, box_len);
                let dist = mic.norm();
                // 允许小量误差
                assert!(
                    dist > min_dist * 0.9,
                    "atoms {} and {} too close: {:.3} Å",
                    i,
                    j,
                    dist
                );
            }
        }
    }

    #[test]
    fn test_build_box_multi_component() {
        let comps = vec![
            Component {
                compound: "water".into(),
                n_molecules: 5,
            },
            Component {
                compound: "methanol".into(),
                n_molecules: 3,
            },
        ];
        let frame = build_box(&comps, 0.8, 1.5, 50).unwrap();
        // 5×3 + 3×6 = 33 个原子
        assert_eq!(frame.atoms.len(), 33);
    }

    #[test]
    fn test_celllist_neighbors_no_duplicates_small_box() {
        // 小盒子使每轴 cell 数退化为 2：旧实现的 27-bin 遍历缺去重，
        // 同一 bin 被多次访问会把近邻原子重复计入（弛豫斥力被成倍累加）。
        let box_len = 4.0;
        let cutoff = 2.0; // n_bins = floor(4 / 2) = 2
        let mut cl = CellList::new(box_len, cutoff);
        assert_eq!(cl.n_bins, 2, "test requires n_bins == 2 to exercise the wrap-around path");

        let atoms = vec![
            Atom::new("Ar", Vector3::new(0.5, 0.5, 0.5)), // bin (0,0,0)
            Atom::new("Ar", Vector3::new(3.5, 3.5, 3.5)), // bin (1,1,1)
            Atom::new("Ar", Vector3::new(1.0, 1.0, 1.0)), // bin (0,0,0)
        ];
        cl.rebuild(&atoms);

        let nbrs = cl.neighbors(&atoms, 0);
        let mut deduped = nbrs.clone();
        deduped.sort_unstable();
        deduped.dedup();
        assert_eq!(
            nbrs.len(), deduped.len(),
            "neighbors must not contain duplicate atom indices, got {nbrs:?}"
        );
    }

    #[test]
    fn test_relax_empty() {
        let mut atoms: Vec<Atom> = Vec::new();
        relax(&mut atoms, 10.0, 1.5, 100);
        assert!(atoms.is_empty());
    }
}
