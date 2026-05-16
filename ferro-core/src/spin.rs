//! 从结构推断未成对电子数 / 自旋多重度。
//!
//! 三条可靠性递减的路径，由 [`guess_spin`] 按优先级选用：
//!
//! 1. **磁矩求和**（最可靠）：若原子带 `magmom`，`n_unpaired = round(|Σ magmom|)`。
//! 2. **氧化态 + Hund 规则**：对离子型固体用电负性规则定氧化态，
//!    过渡金属按 `n_d = 族号 − 氧化态` 取高自旋未成对数，主族按 s/p 填充。
//! 3. **电子数奇偶下限**（兜底）：奇数电子 → 至少 1 个未成对（双重态）。
//!
//! 结果始终做电子数奇偶校验；矛盾时回退到奇偶下限并给出警告。
//!
//! ## 已知限制
//! - 多磁中心间的磁耦合（铁磁 / 反铁磁）无法从结构推断，按铁磁求和给出上限。
//! - 低自旋判定需配位场分析（暂未实现），过渡金属一律按**高自旋**估计。
//! - 共价过渡金属配合物、单质 / 纯共价分子不适用氧化态法，回退到奇偶下限。
//! - 镧系 f 区元素的 f 电子未计入。

use std::collections::BTreeMap;

use crate::data::elements::{
    by_symbol, by_number, group_number, is_transition_metal, symbol_to_z, valence_electrons,
    ElementData,
};
use crate::frame::Frame;

/// 自旋推断采用的方法。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpinMethod {
    /// 由原子 `magmom` 求和得出
    Magmom,
    /// 由氧化态 + Hund 规则得出
    OxidationState,
    /// 仅由电子数奇偶性给出下限
    Parity,
}

/// [`guess_spin`] 的返回结果。
#[derive(Debug, Clone)]
pub struct SpinGuess {
    /// 未成对电子数
    pub n_unpaired: u32,
    /// 自旋多重度 2S+1 = n_unpaired + 1
    pub multiplicity: u32,
    /// 采用的推断方法
    pub method: SpinMethod,
    /// 各 distinct 元素的氧化态（仅 `OxidationState` 方法时为 `Some`）
    pub oxidation_states: Option<Vec<(String, i8)>>,
    /// 推断过程中产生的警告（不可靠、矛盾、歧义等）
    pub warnings: Vec<String>,
}

// ─── 元素数据查找 ─────────────────────────────────────────────────────────────

/// 按符号查元素数据，带 `symbol_to_z` 兜底（兼容 "Fe1" 之类标签）。
fn elem_data(sym: &str) -> Option<&'static ElementData> {
    by_symbol(sym).or_else(|| by_number(symbol_to_z(sym)))
}

// ─── 电子数与奇偶 ─────────────────────────────────────────────────────────────

/// 体系总电子数 = Σ Z_i − 总电荷。
pub fn total_electron_count(frame: &Frame) -> i64 {
    let z_sum: i64 = frame
        .atoms
        .iter()
        .map(|a| symbol_to_z(&a.element) as i64)
        .sum();
    z_sum - frame.charge as i64
}

/// 由电子数奇偶性给出最低自旋多重度：奇数 → 2（双重态），偶数 → 1（单重态）。
pub fn parity_min_multiplicity(frame: &Frame) -> u32 {
    if total_electron_count(frame).rem_euclid(2) == 1 {
        2
    } else {
        1
    }
}

// ─── 氧化态推断（电负性规则）─────────────────────────────────────────────────

/// 用电负性规则 + 电荷守恒推断每个原子的形式氧化态。
///
/// 适用于离子型固体（一个阴离子元素 + 若干阳离子元素）。返回值与
/// `frame.atoms` 顺序对齐。无法判定（单质、纯共价、无明确阴离子、
/// 电荷无法配平）时返回 `None`。
pub fn assign_oxidation_states(frame: &Frame) -> Option<Vec<i8>> {
    // distinct 元素 → 数量，保持确定性顺序
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for a in &frame.atoms {
        *counts.entry(a.element.as_str()).or_insert(0) += 1;
    }
    if counts.len() < 2 {
        return None; // 单质：氧化态法不适用
    }

    // 选阴离子：电负性最高且存在负氧化态的 distinct 元素
    let anion = counts
        .keys()
        .filter_map(|&s| {
            let e = elem_data(s)?;
            let en = e.electronegativity?;
            let min_ox = *e.common_oxidation_states.iter().min()?;
            (min_ox < 0).then_some((s, en, min_ox))
        })
        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    let (anion_sym, _, anion_ox) = anion?;

    let anion_count = counts[anion_sym] as i32;
    let anion_total = anion_count * anion_ox as i32;

    // 阳离子需配平的总氧化态
    let cation_target = frame.charge - anion_total;

    // 每个阳离子 distinct 元素的可选正氧化态
    let mut cations: Vec<(&str, i32, Vec<i8>)> = Vec::new();
    for (&sym, &cnt) in &counts {
        if sym == anion_sym {
            continue;
        }
        let e = elem_data(sym)?;
        let pos: Vec<i8> = e
            .common_oxidation_states
            .iter()
            .copied()
            .filter(|&o| o > 0)
            .collect();
        if pos.is_empty() {
            return None; // 该元素无正氧化态，无法作阳离子
        }
        cations.push((sym, cnt as i32, pos));
    }
    if cations.is_empty() {
        return None;
    }

    // 枚举各阳离子正氧化态的笛卡尔积，找所有配平解
    let mut solutions: Vec<Vec<i8>> = Vec::new();
    let mut current = vec![0i8; cations.len()];
    enumerate_solutions(&cations, cation_target, 0, &mut current, &mut solutions);
    if solutions.is_empty() {
        return None;
    }

    // 多解时取氧化态总和最大者（磷酸盐 / 硅酸盐等形成子取最高价），
    // 由调用方负责给出歧义警告。
    let chosen = solutions
        .into_iter()
        .max_by_key(|sol| {
            sol.iter()
                .zip(&cations)
                .map(|(&o, c)| o as i32 * c.1)
                .sum::<i32>()
        })
        .unwrap();

    // 展开回每个原子
    let mut ox_of: BTreeMap<&str, i8> = BTreeMap::new();
    ox_of.insert(anion_sym, anion_ox);
    for ((sym, _, _), &o) in cations.iter().zip(&chosen) {
        ox_of.insert(sym, o);
    }
    Some(
        frame
            .atoms
            .iter()
            .map(|a| ox_of[a.element.as_str()])
            .collect(),
    )
}

/// 递归枚举笛卡尔积，收集使 Σ count_i·ox_i == target 的组合。
fn enumerate_solutions(
    cations: &[(&str, i32, Vec<i8>)],
    target: i32,
    idx: usize,
    current: &mut Vec<i8>,
    out: &mut Vec<Vec<i8>>,
) {
    if idx == cations.len() {
        let sum: i32 = current
            .iter()
            .zip(cations)
            .map(|(&o, c)| o as i32 * c.1)
            .sum();
        if sum == target {
            out.push(current.clone());
        }
        return;
    }
    for &o in &cations[idx].2 {
        current[idx] = o;
        enumerate_solutions(cations, target, idx + 1, current, out);
    }
}

// ─── 单离子未成对电子数 ───────────────────────────────────────────────────────

/// 给定元素与形式氧化态，估计该离子的未成对电子数（过渡金属取高自旋）。
fn ion_unpaired(z: u8, ox: i8) -> u32 {
    if is_transition_metal(z) {
        let g = group_number(z).unwrap() as i32;
        let nd = (g - ox as i32).clamp(0, 10) as u32;
        return hund(nd, 5); // 5 个 d 轨道，高自旋
    }
    // 主族：离子价电子数填入 s(1)+p(3) 八隅
    let v0 = valence_electrons(z).map(|v| v as i32).unwrap_or(0);
    let v_ion = (v0 - ox as i32).rem_euclid(8) as u32; // 0..8
    if v_ion <= 2 {
        v_ion % 2 // s 亚层：1 个电子 → 1，其余 0
    } else {
        hund(v_ion - 2, 3) // s 满，p 亚层 3 轨道
    }
}

/// Hund 规则：`e` 个电子填入 `n_orb` 个简并轨道的未成对数。
fn hund(e: u32, n_orb: u32) -> u32 {
    let cap = 2 * n_orb;
    let e = e.min(cap);
    if e <= n_orb {
        e
    } else {
        cap - e
    }
}

// ─── 顶层入口 ─────────────────────────────────────────────────────────────────

/// 综合三条路径推断自旋多重度，按可靠性优先：magmom → 氧化态 → 奇偶下限。
pub fn guess_spin(frame: &Frame) -> SpinGuess {
    let total_e = total_electron_count(frame);
    let parity_odd = total_e.rem_euclid(2) == 1;
    let mut warnings = Vec::new();

    // 路径 1：磁矩求和
    if frame.atoms.iter().any(|a| a.magmom.is_some()) {
        let sum: f64 = frame.atoms.iter().filter_map(|a| a.magmom).sum();
        let n = sum.abs().round() as u32;
        let (n, mult) = reconcile_parity(n, parity_odd, &mut warnings, "magmom");
        return SpinGuess {
            n_unpaired: n,
            multiplicity: mult,
            method: SpinMethod::Magmom,
            oxidation_states: None,
            warnings,
        };
    }

    // 路径 2：氧化态 + Hund
    if let Some(ox) = assign_oxidation_states(frame) {
        let n: u32 = frame
            .atoms
            .iter()
            .zip(&ox)
            .map(|(a, &o)| ion_unpaired(symbol_to_z(&a.element), o))
            .sum();

        // distinct 元素氧化态（按符号排序，便于展示）
        let mut distinct: BTreeMap<String, i8> = BTreeMap::new();
        for (a, &o) in frame.atoms.iter().zip(&ox) {
            distinct.insert(a.element.clone(), o);
        }
        let ox_pairs: Vec<(String, i8)> = distinct.into_iter().collect();

        if frame.atoms.iter().any(|a| is_transition_metal(symbol_to_z(&a.element))) {
            warnings.push(
                "过渡金属一律按高自旋估计；低自旋 / 多中心磁耦合需 DFT 验证。".into(),
            );
        }
        let (n, mult) = reconcile_parity(n, parity_odd, &mut warnings, "氧化态");
        return SpinGuess {
            n_unpaired: n,
            multiplicity: mult,
            method: SpinMethod::OxidationState,
            oxidation_states: Some(ox_pairs),
            warnings,
        };
    }

    // 路径 3：奇偶下限
    warnings.push(
        "无法判定氧化态（单质 / 纯共价 / 无明确阴离子），仅给出电子数奇偶下限。".into(),
    );
    let n = u32::from(parity_odd);
    SpinGuess {
        n_unpaired: n,
        multiplicity: n + 1,
        method: SpinMethod::Parity,
        oxidation_states: None,
        warnings,
    }
}

/// 校验未成对数与电子数奇偶是否一致，矛盾则回退到奇偶下限并警告。
fn reconcile_parity(
    n: u32,
    parity_odd: bool,
    warnings: &mut Vec<String>,
    src: &str,
) -> (u32, u32) {
    if (n % 2 == 1) != parity_odd {
        warnings.push(format!(
            "{src} 推断的未成对数 {n} 与电子数奇偶矛盾，回退到奇偶下限。"
        ));
        let n = u32::from(parity_odd);
        (n, n + 1)
    } else {
        (n, n + 1)
    }
}

// ─── 测试 ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::cell::Cell;
    use nalgebra::Vector3;

    fn frame_of(elems: &[(&str, usize)]) -> Frame {
        let cell = Cell::from_lengths_angles(10.0, 10.0, 10.0, 90.0, 90.0, 90.0).unwrap();
        let mut f = Frame::with_cell(cell, [true; 3]);
        let mut x = 0.0;
        for &(sym, n) in elems {
            for _ in 0..n {
                f.add_atom(Atom::new(sym, Vector3::new(x, 0.0, 0.0)));
                x += 1.5;
            }
        }
        f
    }

    // ── 用户给定的两个验证案例 ──────────────────────────────────────────────

    #[test]
    fn test_znp2o6_diamagnetic() {
        let f = frame_of(&[("Zn", 1), ("P", 2), ("O", 6)]);
        let ox = assign_oxidation_states(&f).expect("应能定氧化态");
        // 顺序：Zn, P, P, O×6
        assert_eq!(ox[0], 2, "Zn 应为 +2");
        assert_eq!(ox[1], 5, "P 应为 +5");
        assert_eq!(ox[3], -2, "O 应为 −2");

        let g = guess_spin(&f);
        assert_eq!(g.n_unpaired, 0, "ZnP2O6 应抗磁");
        assert_eq!(g.multiplicity, 1);
        assert_eq!(g.method, SpinMethod::OxidationState);
    }

    #[test]
    fn test_mns_high_spin_d5() {
        let f = frame_of(&[("Mn", 1), ("S", 1)]);
        let ox = assign_oxidation_states(&f).expect("应能定氧化态");
        assert_eq!(ox[0], 2, "Mn 应为 +2");
        assert_eq!(ox[1], -2, "S 应为 −2");

        let g = guess_spin(&f);
        assert_eq!(g.n_unpaired, 5, "Mn²⁺ d⁵ 高自旋应有 5 个未成对电子");
        assert_eq!(g.multiplicity, 6);
        assert_eq!(g.method, SpinMethod::OxidationState);
    }

    // ── 其它已知体系 ────────────────────────────────────────────────────────

    #[test]
    fn test_fe2o3_fe3plus() {
        let f = frame_of(&[("Fe", 2), ("O", 3)]);
        let ox = assign_oxidation_states(&f).unwrap();
        assert_eq!(ox[0], 3, "Fe 应为 +3");
        assert_eq!(ox[2], -2);
        // Fe³⁺ = d⁵ → 高自旋 5，两个 Fe → 10
        let g = guess_spin(&f);
        assert_eq!(g.n_unpaired, 10);
        assert_eq!(g.multiplicity, 11);
    }

    #[test]
    fn test_feo_fe2plus_d6() {
        let f = frame_of(&[("Fe", 1), ("O", 1)]);
        let ox = assign_oxidation_states(&f).unwrap();
        assert_eq!(ox[0], 2);
        // Fe²⁺ = d⁶ → 高自旋 4
        let g = guess_spin(&f);
        assert_eq!(g.n_unpaired, 4);
        assert_eq!(g.multiplicity, 5);
    }

    #[test]
    fn test_nacl_closed_shell() {
        let f = frame_of(&[("Na", 1), ("Cl", 1)]);
        let ox = assign_oxidation_states(&f).unwrap();
        assert_eq!(ox[0], 1, "Na +1");
        assert_eq!(ox[1], -1, "Cl −1");
        let g = guess_spin(&f);
        assert_eq!(g.n_unpaired, 0);
        assert_eq!(g.multiplicity, 1);
    }

    #[test]
    fn test_single_element_falls_back_to_parity() {
        // 单个 Na 原子：氧化态法不适用，回退奇偶。Na Z=11 奇数 → 双重态
        let f = frame_of(&[("Na", 1)]);
        assert!(assign_oxidation_states(&f).is_none());
        let g = guess_spin(&f);
        assert_eq!(g.method, SpinMethod::Parity);
        assert_eq!(g.multiplicity, 2);
    }

    #[test]
    fn test_magmom_path_takes_priority() {
        let mut f = frame_of(&[("Fe", 1), ("O", 1)]);
        // 显式磁矩覆盖：总和 ≈ 2 → 三重态
        f.atoms[0].magmom = Some(2.0);
        f.atoms[1].magmom = Some(0.0);
        let g = guess_spin(&f);
        assert_eq!(g.method, SpinMethod::Magmom);
        assert_eq!(g.n_unpaired, 2);
        assert_eq!(g.multiplicity, 3);
    }

    #[test]
    fn test_parity_helpers() {
        let f = frame_of(&[("Mn", 1), ("S", 1)]); // 25+16 = 41 奇
        assert_eq!(total_electron_count(&f), 41);
        assert_eq!(parity_min_multiplicity(&f), 2);

        let f2 = frame_of(&[("Zn", 1), ("P", 2), ("O", 6)]); // 30+30+48 = 108 偶
        assert_eq!(total_electron_count(&f2), 108);
        assert_eq!(parity_min_multiplicity(&f2), 1);
    }

    #[test]
    fn test_hund_rule() {
        assert_eq!(hund(0, 5), 0);
        assert_eq!(hund(3, 5), 3);  // d³
        assert_eq!(hund(5, 5), 5);  // d⁵
        assert_eq!(hund(6, 5), 4);  // d⁶ 高自旋
        assert_eq!(hund(7, 5), 3);  // d⁷
        assert_eq!(hund(10, 5), 0); // d¹⁰
        assert_eq!(hund(3, 3), 3);  // p³
        assert_eq!(hund(4, 3), 2);  // p⁴
    }

    #[test]
    fn test_charged_system_parity() {
        // NH4+ 形式：N + 4H, 电荷 +1 → 电子数 7+4-1 = 10 偶
        let mut f = frame_of(&[("N", 1), ("H", 4)]);
        f.charge = 1;
        assert_eq!(total_electron_count(&f), 10);
        assert_eq!(parity_min_multiplicity(&f), 1);
    }
}
