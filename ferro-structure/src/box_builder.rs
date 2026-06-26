//! Mixed-compound box construction: estimate size, random placement, soft-core relaxation.

use std::collections::HashMap;

use nalgebra::Vector3;
use rand::Rng;

use ferro_core::atom::Atom;
use ferro_core::cell::Cell;
use ferro_core::data::compounds;
use ferro_core::error::{ChemError, Result};
use ferro_core::frame::Frame;

/// A component in a mixed-compound system.
pub struct Component {
    /// Compound name or formula, looked up in the COMPOUNDS database.
    pub compound: String,
    /// Number of molecules of this compound.
    pub n_molecules: usize,
}

/// Parse a simple chemical formula into element counts.
///
/// Supports formulas without parentheses: H2O, P2O5, ZnO, CH3OH, CCl4.
/// Elements are one uppercase letter optionally followed by one lowercase letter.
/// Count defaults to 1 if omitted.
fn parse_formula(formula: &str) -> Result<HashMap<String, usize>> {
    let mut counts: HashMap<String, usize> = HashMap::new();
    let bytes = formula.as_bytes();
    let n = bytes.len();
    let mut i = 0;

    while i < n {
        // 解析元素符号：大写字母 + 可选小写字母
        if !bytes[i].is_ascii_uppercase() {
            return Err(ChemError::ParseError(format!(
                "expected uppercase letter at position {} in '{}'",
                i, formula
            )));
        }
        let start = i;
        i += 1;
        while i < n && bytes[i].is_ascii_lowercase() {
            i += 1;
        }
        let elem = &formula[start..i];

        // 解析数字（可选，默认 1）
        let mut count_str = String::new();
        while i < n && bytes[i].is_ascii_digit() {
            count_str.push(bytes[i] as char);
            i += 1;
        }
        let count: usize = if count_str.is_empty() {
            1
        } else {
            count_str.parse().map_err(|_| {
                ChemError::ParseError(format!("invalid count '{}' in '{}'", count_str, formula))
            })?
        };

        *counts.entry(elem.to_string()).or_insert(0) += count;
    }

    Ok(counts)
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
        let cd = compounds::find(&comp.compound).ok_or_else(|| {
            ChemError::ValidationError(format!("compound '{}' not found", comp.compound))
        })?;
        total_mass += comp.n_molecules as f64 * cd.molecular_mass;
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
        let cd = compounds::find(&comp.compound).ok_or_else(|| {
            ChemError::ValidationError(format!("compound '{}' not found", comp.compound))
        })?;
        let formula_counts = parse_formula(cd.formula)?;
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
    let mut rng = rand::thread_rng();
    let mut atoms: Vec<Atom> = Vec::new();
    for (elem, count) in &element_counts {
        for _ in 0..*count {
            let pos = Vector3::new(
                rng.gen_range(0.0..box_len),
                rng.gen_range(0.0..box_len),
                rng.gen_range(0.0..box_len),
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
        let counts = parse_formula("").unwrap();
        assert!(counts.is_empty());
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
