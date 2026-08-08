//! Integration tests for per-frame volume normalisation on real trajectories.
//!
//! These live in `ferro-cli` because they need both `ferro-io` (to read the dumps) and
//! `ferro-analysis` (to run the analysis); the middle-layer crates are not allowed to
//! depend on each other, and `ferro-cli` is the entry-point layer that may combine them.
//!
//! Fixtures are 5-frame slices of the production trajectories, sampled at even intervals
//! so the NPT one retains the full volume span of its parent (a contiguous slice would
//! have shown almost no fluctuation and tested nothing).

use ferro_analysis::md::gr::{calc_gr, GrParams};
use ferro_analysis::md::sq::{calc_sq_from_gr, SqParams};
use ferro_core::Trajectory;
use ferro_io::{read_lammps_dump, LammpsUnits};
use std::collections::BTreeMap;

const NPT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../tests/43Z43P15A_NPT_5.lammpstrj");
const NVT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../tests/70Z30P00A_NVT_5.lammpstrj");

fn gr_params() -> GrParams {
    // r_max 固定且远小于任一帧的最小镜像上界，使逐帧与全轨迹的 bin 数一致
    GrParams { r_min: 0.005, r_max: 8.0, dr: 0.02, ..Default::default() }
}

fn sq_params() -> SqParams {
    SqParams { q_min: 0.5, q_max: 20.0, dq: 0.05, ..Default::default() }
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()).fold(0.0_f64, f64::max)
}

/// Time-average of per-frame results — the literal reading of code1/code2, used as the
/// oracle the production (folded) path must reproduce.
fn literal_per_frame(traj: &Trajectory) -> (BTreeMap<String, Vec<f64>>, BTreeMap<String, Vec<f64>>) {
    let (gp, sp) = (gr_params(), sq_params());
    let mut gr_acc: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut sq_acc: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for frame in &traj.frames {
        let single = Trajectory::from_frame(frame.clone());
        let g = calc_gr(&single, &gp).unwrap();
        let s = calc_sq_from_gr(&g, &sp);
        for (k, v) in &g.gr {
            let e = gr_acc.entry(k.clone()).or_insert_with(|| vec![0.0; v.len()]);
            for (a, b) in e.iter_mut().zip(v.iter()) { *a += b; }
        }
        for (k, v) in &s.sq {
            let e = sq_acc.entry(k.clone()).or_insert_with(|| vec![0.0; v.len()]);
            for (a, b) in e.iter_mut().zip(v.iter()) { *a += b; }
        }
    }
    let nf = traj.frames.len() as f64;
    for v in gr_acc.values_mut() { for x in v.iter_mut() { *x /= nf; } }
    for v in sq_acc.values_mut() { for x in v.iter_mut() { *x /= nf; } }
    (gr_acc, sq_acc)
}

#[test]
fn npt_fixture_actually_fluctuates() {
    // 守住前提：若哪天有人换掉 fixture 而新文件体积恒定，下面的测试会全部空转通过。
    let traj = read_lammps_dump(NPT, LammpsUnits::Real).unwrap();
    let g = calc_gr(&traj, &gr_params()).unwrap();
    assert_eq!(traj.frames.len(), 5);
    let rel = g.volume_std / g.avg_volume;
    assert!(rel > 0.02, "NPT fixture should fluctuate by >2% (got {:.4}%)", 100.0 * rel);
}

#[test]
fn nvt_fixture_has_constant_volume() {
    let traj = read_lammps_dump(NVT, LammpsUnits::Real).unwrap();
    let g = calc_gr(&traj, &gr_params()).unwrap();
    assert_eq!(traj.frames.len(), 5);
    // 逐位相同的体积经 Σv/n 仍可能与单个 v 差 1 ulp（V≈6e4 时约 7e-12），
    // 故按相对量判定「恒容」，而不是要求标准差恰为 0
    let rel = g.volume_std / g.avg_volume;
    assert!(rel < 1e-14, "NVT fixture must have a constant box (relative spread {rel:.3e})");
}

#[test]
fn nvt_is_unaffected_by_per_frame_normalisation() {
    // 定容轨迹上 ⟨ρ_f·g_f⟩ 精确等于 ρ·g，故 S(q) 与旧口径逐位相同 ——
    // 这是「本次改动不动存量 NVT 结果」的可执行证明。
    let traj = read_lammps_dump(NVT, LammpsUnits::Real).unwrap();
    let g = calc_gr(&traj, &gr_params()).unwrap();
    for (k, v) in &g.gr {
        for (i, &gv) in v.iter().enumerate() {
            let expect = g.rho * gv;
            assert!((g.rho_g[k][i] - expect).abs() <= 1e-9 * (1.0 + expect.abs()),
                "pair {k} bin {i}: rho_g deviates from rho*gr under NVT");
        }
    }
    assert!((g.rho - g.element_counts.values().sum::<usize>() as f64 / g.avg_volume).abs() < 1e-12,
        "under constant volume N*<1/V> and N/<V> must coincide");
}

#[test]
fn folded_matches_literal_per_frame_on_real_npt_trajectory() {
    // 生产路径只做一次傅里叶变换，靠 ⟨ρ_f·g_f⟩ 重构逐帧结果；这里在真实 NPT 轨迹上
    // 与字面逐帧实现对拍，确认两者是同一个数而非近似。
    let traj = read_lammps_dump(NPT, LammpsUnits::Real).unwrap();
    let g = calc_gr(&traj, &gr_params()).unwrap();
    let s = calc_sq_from_gr(&g, &sq_params());
    let (gr_lit, sq_lit) = literal_per_frame(&traj);

    for (k, v) in &g.gr {
        let d = max_abs_diff(v, &gr_lit[k]);
        assert!(d < 1e-10, "g(r) pair {k}: folded vs literal differs by {d:.3e}");
    }
    for (k, v) in &s.sq {
        let d = max_abs_diff(v, &sq_lit[k]);
        assert!(d < 1e-10, "S(q) pair {k}: folded vs literal differs by {d:.3e}");
    }
}

#[test]
fn npt_differs_from_average_volume_shortcut() {
    // 量化旧口径（先平均 g，再用 ρ=N/⟨V⟩ 变换一次）与逐帧口径的差距，确认改动
    // 在 NPT 上确实改变了结果 —— 否则前一个测试可能只是在比较两条相同的死路径。
    let traj = read_lammps_dump(NPT, LammpsUnits::Real).unwrap();
    let g = calc_gr(&traj, &gr_params()).unwrap();
    let s = calc_sq_from_gr(&g, &sq_params());

    let pi4 = 4.0 * std::f64::consts::PI;
    let n_total = g.element_counts.values().sum::<usize>() as f64;
    let rho_old = n_total / g.avg_volume;
    let dr = g.params.dr;

    let mut worst = 0.0_f64;
    for (k, gr_vals) in &g.gr {
        let Some(new) = s.sq.get(k) else { continue };
        let old: Vec<f64> = s.q.iter().map(|&qi| {
            if qi.abs() < 1e-10 { return 1.0; }
            let integral: f64 = g.r.iter().zip(gr_vals.iter())
                .map(|(&ri, &gv)| ri * (gv - 1.0) * (qi * ri).sin() * dr)
                .sum();
            1.0 + pi4 * rho_old / qi * integral
        }).collect();
        worst = worst.max(max_abs_diff(new, &old));
    }
    assert!(worst > 1e-4,
        "per-frame S(q) must differ measurably from the average-volume shortcut under NPT \
         (max diff {worst:.3e})");
}
