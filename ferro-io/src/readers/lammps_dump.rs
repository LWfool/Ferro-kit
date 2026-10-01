use std::path::Path;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use ferro_core::data::elements::{LabelSplit, split_element_label};
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{bail, Context, Result};

type AtomRaw = (usize, Atom, Option<Vector3<f64>>, Option<Vector3<f64>>);

const KCAL_TO_EV: f64 = 0.04336410;

/// LAMMPS unit system for dump files.
///
/// Positions are always Å in both systems. Differences:
/// - `Real`:  velocities Å/fs, forces kcal/(mol·Å)
/// - `Metal`: velocities Å/ps, forces eV/Å
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum LammpsUnits {
    #[default]
    Real,
    Metal,
}

pub fn read_lammps_dump(path: &Path, units: LammpsUnits) -> Result<Trajectory> {
    let path_ = path.display();
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path_}"))?;
    parse_lammps_dump(&content, units).with_context(|| format!("parsing {path_}"))
}

fn parse_lammps_dump(content: &str, units: LammpsUnits) -> Result<Trajectory> {
    let lines: Vec<&str> = content.lines().collect();
    let mut i = 0;
    let mut traj = Trajectory::new();
    // element 列原始串 → 拆分后的元素符号；仅用于读取结束时打印一次映射表
    let mut site_map: BTreeMap<String, String> = BTreeMap::new();
    let mut unknown_prefixes: BTreeSet<String> = BTreeSet::new();

    while i < lines.len() {
        // Find ITEM: TIMESTEP
        if !lines[i].trim().starts_with("ITEM: TIMESTEP") { i += 1; continue; }
        let start = i;
        match parse_frame(&lines, &mut i, units, &mut site_map, &mut unknown_prefixes) {
            Ok(frame) => traj.add_frame(frame),
            Err(why) => {
                let k = traj.n_frames();
                // 不完整的帧后面还有帧 = 文件中间坏了，报错；只有末帧不完整才是
                // MD 被中断的常态（dump 写到一半），丢掉它，前面的帧照常可用
                if lines[start + 1..].iter().any(|l| l.trim().starts_with("ITEM: TIMESTEP")) {
                    bail!("frame {k} (line {}) is incomplete: {why}", start + 1);
                }
                eprintln!(
                    "[ferro] warning: last frame (frame {k}, line {}) is incomplete and was \
                     dropped: {why}",
                    start + 1
                );
                break;
            }
        }
    }

    report_site_map(&site_map, &unknown_prefixes);
    Ok(traj)
}

/// Parse the frame whose `ITEM: TIMESTEP` line is at `*i`, advancing `*i` past it.
///
/// `Err` means only "this frame is incomplete" and carries the reason; whether that is
/// fatal is the caller's call (it depends on whether more frames follow).
fn parse_frame(
    lines: &[&str],
    i: &mut usize,
    units: LammpsUnits,
    site_map: &mut BTreeMap<String, String>,
    unknown_prefixes: &mut BTreeSet<String>,
) -> std::result::Result<Frame, String> {
    // 找本帧的下一个段头；先撞上下一帧的 TIMESTEP 或到文件尾，都说明本帧缺这一段
    let seek = |from: usize, tag: &str| -> std::result::Result<usize, String> {
        for (k, l) in lines.iter().enumerate().skip(from) {
            if l.contains(tag) { return Ok(k); }
            if l.trim().starts_with("ITEM: TIMESTEP") { break; }
        }
        Err(format!("no ITEM: {tag} section"))
    };

    // NUMBER OF ATOMS
    let na = seek(*i + 1, "NUMBER OF ATOMS")?;
    let n: usize = lines.get(na + 1)
        .and_then(|l| l.trim().parse().ok())
        .ok_or("NUMBER OF ATOMS has no count")?;

    // BOX BOUNDS
    let bb = seek(na + 2, "BOX BOUNDS")?;
    let is_triclinic = lines[bb].contains("xy");

    let mut lo = [0.0_f64; 3];
    let mut hi = [0.0_f64; 3];
    let mut tilt = [0.0_f64; 3]; // xy, xz, yz

    for dim in 0..3 {
        let vals: Vec<f64> = lines.get(bb + 1 + dim)
            .map(|l| l.split_whitespace().map_while(|s| s.parse().ok()).collect())
            .unwrap_or_default();
        if vals.len() < 2 {
            return Err(format!("BOX BOUNDS line {} of 3 is missing or cut short", dim + 1));
        }
        lo[dim] = vals[0];
        hi[dim] = vals[1];
        if vals.len() >= 3 {
            tilt[dim] = vals[2]; // xy on dim 0, xz on dim 1, yz on dim 2
        }
    }

    // User request C: also detect triclinic by value count (≥3 values per line)
    // already handled above: vals.len() >= 3 sets tilt

    let (xy, xz, yz) = (tilt[0], tilt[1], tilt[2]);
    // 三斜时这六个数是 *_bound（倾斜后的外接盒），真实边长要把倾斜量减回去：
    //   lx = (xhi_b - xlo_b) - |xy| - |xz|,  ly = (yhi_b - ylo_b) - |yz|
    // 该式与规格的 MAX(0,xy,xz,xy+xz) - MIN(...) 在四个象限上恒等，
    // 与 ase/io/lammpsrun.py::construct_cell 一致。正交时三个倾斜量为 0，
    // 退化成 hi - lo，故这个分支对既有的正交轨迹逐位不变
    let lx = (hi[0] - lo[0]) - xy.abs() - xz.abs();
    let ly = (hi[1] - lo[1]) - yz.abs();
    let lz = hi[2] - lo[2];

    let cell = if is_triclinic || xy != 0.0 || xz != 0.0 || yz != 0.0 {
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

    // ITEM: ATOMS col1 col2 ...
    let ah = seek(bb + 4, "ITEM: ATOMS")?;

    let col_names: Vec<&str> = lines[ah]
        .trim_start_matches("ITEM: ATOMS")
        .split_whitespace()
        .collect();
    let col: HashMap<&str, usize> = col_names.iter()
        .enumerate()
        .map(|(idx, &name)| (name, idx))
        .collect();

    let get_col = |name: &str| col.get(name).copied();

    let mut atoms_raw: Vec<AtomRaw> = Vec::new();

    for k in 0..n {
        // 文件尾、下一个段头、空行都说明原子段没写完；列数不足是最后一行被截在中间
        let line = lines.get(ah + 1 + k)
            .filter(|l| !l.trim().is_empty() && !l.trim_start().starts_with("ITEM:"))
            .ok_or_else(|| format!("only {k} of {n} atom lines"))?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < col_names.len() {
            return Err(format!(
                "atom line {} of {n} has {} of {} columns",
                k + 1, parts.len(), col_names.len()
            ));
        }

        let atom_id: usize = get_col("id")
            .and_then(|c| parts.get(c))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        let tp: usize = get_col("type")
            .and_then(|c| parts.get(c))
            .and_then(|s| s.parse().ok())
            .unwrap_or(1);

        // element 列可能写的是位点类型标签（`O_b_P_P`、`Zn_f`），按第一个下划线
        // 拆成 element + label；无下划线的普通符号原样通过。
        let raw = get_col("element")
            .and_then(|c| parts.get(c))
            .map(|s| s.to_string())
            .unwrap_or_else(|| format!("X{tp}"));
        let (element, label) = match split_element_label(&raw) {
            LabelSplit::Plain(e) => (e.to_string(), None),
            LabelSplit::Split { element, label } => {
                (element.to_string(), Some(label.to_string()))
            }
            LabelSplit::Unknown(s) => {
                unknown_prefixes.insert(s.to_string());
                (s.to_string(), None)
            }
        };
        site_map.entry(raw).or_insert_with(|| element.clone());

        // Position — try x/y/z first, then xs/ys/zs (scaled), then xu/yu/zu (unwrapped)
        let pos = if let (Some(cx), Some(cy), Some(cz)) =
            (get_col("x"), get_col("y"), get_col("z"))
        {
            let x: f64 = parts.get(cx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let y: f64 = parts.get(cy).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let z: f64 = parts.get(cz).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            Vector3::new(x, y, z)
        } else if let (Some(cx), Some(cy), Some(cz)) =
            (get_col("xs"), get_col("ys"), get_col("zs"))
        {
            // Scaled [0,1) → Cartesian via cell
            let sx: f64 = parts.get(cx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let sy: f64 = parts.get(cy).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let sz: f64 = parts.get(cz).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            cell.fractional_to_cartesian(Vector3::new(sx, sy, sz))
        } else if let (Some(cx), Some(cy), Some(cz)) =
            (get_col("xu"), get_col("yu"), get_col("zu"))
        {
            let x: f64 = parts.get(cx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let y: f64 = parts.get(cy).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let z: f64 = parts.get(cz).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            Vector3::new(x, y, z)
        } else {
            Vector3::zeros()
        };

        let mut atom = Atom::new(element, pos);
        atom.label = label;
        // Charge
        if let Some(c) = get_col("q") {
            atom.charge = parts.get(c).and_then(|s| s.parse().ok());
        }

        // Velocity: real Å/fs (internal), metal Å/ps → ×1e-3
        let vel = if let (Some(vx), Some(vy), Some(vz)) =
            (get_col("vx"), get_col("vy"), get_col("vz"))
        {
            let x: f64 = parts.get(vx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let y: f64 = parts.get(vy).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let z: f64 = parts.get(vz).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let vscale = match units {
                LammpsUnits::Real  => 1.0,
                LammpsUnits::Metal => 1e-3,
            };
            Some(Vector3::new(x, y, z) * vscale)
        } else { None };

        // Force: real kcal/(mol·Å) → eV/Å, metal eV/Å already
        let fscale = match units {
            LammpsUnits::Real  => KCAL_TO_EV,
            LammpsUnits::Metal => 1.0,
        };
        let force = if let (Some(fx), Some(fy), Some(fz)) =
            (get_col("fx"), get_col("fy"), get_col("fz"))
        {
            let x: f64 = parts.get(fx).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let y: f64 = parts.get(fy).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let z: f64 = parts.get(fz).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            Some(Vector3::new(x * fscale, y * fscale, z * fscale))
        } else { None };

        atoms_raw.push((atom_id, atom, vel, force));
    }

    // Sort by atom id
    atoms_raw.sort_by_key(|(id, _, _, _)| *id);

    let mut frame = Frame::with_cell(cell, [true; 3]);
    let mut all_vels = Vec::new();
    let mut all_forces = Vec::new();
    let mut has_vel = false;
    let mut has_force = false;

    for (_, atom, vel, force) in atoms_raw {
        frame.add_atom(atom);
        if let Some(v) = vel { all_vels.push(v); has_vel = true; }
        else { all_vels.push(Vector3::zeros()); }
        if let Some(f) = force { all_forces.push(f); has_force = true; }
        else { all_forces.push(Vector3::zeros()); }
    }

    if has_vel { frame.velocities = Some(all_vels); }
    if has_force { frame.forces = Some(all_forces); }

    *i = ah + 1 + n;
    Ok(frame)
}

/// Print the raw-string → element mapping once per file, so a silently wrong split
/// is visible rather than buried in the results.
///
/// Stays quiet for ordinary trajectories where every element column entry is already
/// a plain chemical symbol.
fn report_site_map(site_map: &BTreeMap<String, String>, unknown: &BTreeSet<String>) {
    let split: Vec<(&String, &String)> = site_map.iter()
        .filter(|(raw, elem)| raw.as_str() != elem.as_str())
        .collect();
    if split.is_empty() && unknown.is_empty() { return; }

    if !split.is_empty() {
        eprintln!("[ferro] LAMMPS dump: element column carries site labels, split into element + label:");
        for (raw, elem) in &split {
            eprintln!("[ferro]   {raw:<12} -> element {elem:<4} label {raw}");
        }
        eprintln!("[ferro]   -a/-b/-c select by element, -x/-y/-z select by label");
    }
    if !unknown.is_empty() {
        let list: Vec<&str> = unknown.iter().map(|s| s.as_str()).collect();
        eprintln!(
            "[ferro] warning: underscore present but prefix is not a known element, \
             kept as element verbatim: {}",
            list.join(", ")
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DUMP_ORTHO: &str = "ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 2.87
0 2.87
0 2.87
ITEM: ATOMS id type element x y z
1 1 Fe 0.0   0.0   0.0
2 1 Fe 1.435 1.435 1.435
ITEM: TIMESTEP
10
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 2.87
0 2.87
0 2.87
ITEM: ATOMS id type element x y z
1 1 Fe 0.01  0.0   0.0
2 1 Fe 1.445 1.435 1.435
";

    use crate::testutil::write_tmp as tmp;

    const DUMP_VEL: &str = "ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0 5.0
0 5.0
0 5.0
ITEM: ATOMS id type element x y z vx vy vz
1 1 Fe 0.0 0.0 0.0 2.0 0.0 0.0
";

    const DUMP_LABELS: &str = "ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
5
ITEM: BOX BOUNDS pp pp pp
0 10.0
0 10.0
0 10.0
ITEM: ATOMS id type element x y z
1 1 P_0     0.0 0.0 0.0
2 2 O_b_P_P 1.5 0.0 0.0
3 2 O_f     3.0 0.0 0.0
4 3 Zn_f    4.5 0.0 0.0
5 2 O       6.0 0.0 0.0
";

    const DUMP_BAD_LABEL: &str = "ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 10.0
0 10.0
0 10.0
ITEM: ATOMS id type element x y z
1 1 foo_bar 0.0 0.0 0.0
2 2 O       1.5 0.0 0.0
";

    #[test]
    fn test_label_split_fills_element_and_label() {
        let traj = read_lammps_dump(&tmp("labels.dump", DUMP_LABELS), LammpsUnits::Real).unwrap();
        let f = traj.first().unwrap();
        let got: Vec<(&str, Option<&str>)> = f.atoms.iter()
            .map(|a| (a.element.as_str(), a.label.as_deref()))
            .collect();
        assert_eq!(got, vec![
            ("P",  Some("P_0")),
            ("O",  Some("O_b_P_P")),   // 后缀含下划线，只在第一个处拆
            ("O",  Some("O_f")),
            ("Zn", Some("Zn_f")),
            ("O",  None),              // 普通符号：label 保持 None
        ]);
    }

    #[test]
    fn test_plain_element_column_leaves_label_none() {
        // 普通轨迹不受拆分影响
        let traj = read_lammps_dump(&tmp("plain.dump", DUMP_ORTHO), LammpsUnits::Real).unwrap();
        for atom in &traj.first().unwrap().atoms {
            assert_eq!(atom.element, "Fe");
            assert!(atom.label.is_none(), "plain symbol should not produce a label");
        }
    }

    #[test]
    fn test_unknown_prefix_kept_verbatim_as_element() {
        let traj = read_lammps_dump(&tmp("bad.dump", DUMP_BAD_LABEL), LammpsUnits::Real).unwrap();
        let f = traj.first().unwrap();
        // 前缀非法 → 整串当元素、label 为 None（告警走 stderr）
        assert_eq!(f.atom(0).element, "foo_bar");
        assert!(f.atom(0).label.is_none());
        assert_eq!(f.atom(1).element, "O");
    }

    #[test]
    fn test_multiframe() {
        let traj = read_lammps_dump(&tmp("bcc.dump", DUMP_ORTHO), LammpsUnits::Real).unwrap();
        assert_eq!(traj.n_frames(), 2);
        assert_eq!(traj.first().unwrap().n_atoms(), 2);
        assert_eq!(traj.first().unwrap().atom(0).element, "Fe");
    }

    /// `DUMP_ORTHO` 的前 `n` 行；每帧 11 行，第 1 帧从第 12 行开始
    fn head(n: usize) -> String {
        DUMP_ORTHO.lines().take(n).map(|l| format!("{l}\n")).collect()
    }

    #[test]
    fn test_truncated_last_frame_is_dropped() {
        // MD 被中断时 dump 末帧写到一半：截在原子段、截在一行中间、截在段头之间
        let cases = [
            ("原子段只写了 1 行", head(21)),
            ("最后一行被截在中间", head(21) + "2 1 Fe 1.44"),
            ("截在 BOX BOUNDS 之前", head(15)),
            ("BOX BOUNDS 只写了一行", head(17)),
            ("只有 TIMESTEP 一行", head(12)),
            ("ATOMS 表头之后为空", head(20)),
        ];
        for (what, text) in cases {
            let traj = parse_lammps_dump(&text, LammpsUnits::Real)
                .unwrap_or_else(|e| panic!("{what}：应丢掉末帧而不是报错，实际 {e:#}"));
            assert_eq!(traj.n_frames(), 1, "{what}：应只剩第 0 帧");
            assert_eq!(traj.first().unwrap().n_atoms(), 2, "{what}：第 0 帧应完整");
        }
    }

    #[test]
    fn test_incomplete_middle_frame_is_an_error() {
        // 第 0 帧少一行原子、后面还有完整的第 1 帧 = 文件中间坏了，不能静默丢帧
        let mut lines: Vec<&str> = DUMP_ORTHO.lines().collect();
        lines.remove(10);
        let err = parse_lammps_dump(&lines.join("\n"), LammpsUnits::Real)
            .expect_err("中间帧不完整应报错");
        let msg = format!("{err:#}");
        assert!(msg.contains("frame 0") && msg.contains("1 of 2 atom lines"), "报错应点名帧与行数：{msg}");
    }

    #[test]
    fn test_box_bounds() {
        let traj = read_lammps_dump(&tmp("box.dump", DUMP_ORTHO), LammpsUnits::Real).unwrap();
        let [a, ..] = traj.first().unwrap().cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-6);
    }

    #[test]
    fn test_velocity_real_units() {
        let traj = read_lammps_dump(&tmp("vel_real.dump", DUMP_VEL), LammpsUnits::Real).unwrap();
        let vx = traj.first().unwrap().velocities.as_ref().unwrap()[0].x;
        assert!((vx - 2.0).abs() < 1e-10, "real units: vx should be 2.0 Å/fs, got {vx}");
    }

    #[test]
    fn test_velocity_metal_units() {
        let traj = read_lammps_dump(&tmp("vel_metal.dump", DUMP_VEL), LammpsUnits::Metal).unwrap();
        let vx = traj.first().unwrap().velocities.as_ref().unwrap()[0].x;
        // 2.0 Å/ps × 1e-3 = 0.002 Å/fs
        assert!((vx - 0.002).abs() < 1e-12, "metal units: vx should be 0.002 Å/fs, got {vx}");
    }

    /// 三斜盒子行是 `*_bound`，真实边长要把倾斜量减回去。
    ///
    /// The numbers come from `tests/triclinic_2frames.lammpstrj`, whose cell was
    /// cross-checked against `ase.io.read(..., format="lammps-dump-text")`:
    /// lx/ly/lz = 10/12/14 with xy/xz/yz = 2/-3/1. Reading the bounds as xlo/xhi
    /// would give lx = 15 and ly = 13 — plausible numbers, wrong cell.
    #[test]
    fn test_triclinic_bounds_are_not_xlo_xhi() {
        const TRI: &str = "ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS xy xz yz pp pp pp
-3.0 12.0 2.0
0.0 13.0 -3.0
0.0 14.0 1.0
ITEM: ATOMS id type element x y z
1 1 Si 0.5 0.5 0.5
";
        let traj = read_lammps_dump(&tmp("tri.dump", TRI), LammpsUnits::Metal).unwrap();
        let m = traj.first().unwrap().cell.as_ref().unwrap().matrix;
        // 行优先：行 = 晶格矢量
        let want = [[10.0, 0.0, 0.0], [2.0, 12.0, 0.0], [-3.0, 1.0, 14.0]];
        for (i, row) in want.iter().enumerate() {
            for (j, w) in row.iter().enumerate() {
                assert!((m[(i, j)] - w).abs() < 1e-10, "cell[{i}][{j}]: got {}, want {w}", m[(i, j)]);
            }
        }
    }
}
