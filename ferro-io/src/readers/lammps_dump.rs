use std::path::Path;
use crate::Notice;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use ferro_core::data::elements::by_symbol;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use ferro_core::units::{convert_time, TimeUnit, KCAL_MOL_TO_EV};
use nalgebra::{Matrix3, Vector3};
use anyhow::{bail, Context, Result};

type AtomRaw = (usize, Atom, Option<Vector3<f64>>, Option<Vector3<f64>>);

/// LAMMPS unit system for dump files.
///
/// Positions are always Å in both systems. Differences:
/// - `Real`:  velocities Å/fs, forces kcal/(mol·Å)
/// - `Metal`: velocities Å/ps, forces eV/Å
///
/// 没有默认值：dump 不记录 `units` 命令，两种体系下速度差 1000 倍、力差 23 倍，
/// 猜错是静默的（DeePMD 只能跑 metal，而旧默认是 real）。调用者传 `Option`，
/// `None` 只在文件不含速度/力列时读得通。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LammpsUnits {
    Real,
    Metal,
}

impl std::str::FromStr for LammpsUnits {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, String> {
        match s {
            "real" => Ok(LammpsUnits::Real),
            "metal" => Ok(LammpsUnits::Metal),
            _ => Err(format!("unknown LAMMPS units '{s}' (expected real or metal)")),
        }
    }
}

/// `units` 为 `None` 时，含 `vx`/`fx` 等速度或力列的 dump 报错 —— 这些列的量纲取决于
/// LAMMPS 的 `units` 命令，文件里看不出来。
pub fn read_lammps_dump(path: &Path, units: Option<LammpsUnits>) -> Result<Trajectory> {
    Ok(read_lammps_dump_with_warnings(path, units)?.0)
}

/// Reads a LAMMPS dump, returning what the reader noticed alongside the frames.
///
/// The warnings name a dropped incomplete last frame and the element-column split
/// (`P_3` -> element `P`, label `P_3`); the caller prints them, the library never does.
pub fn read_lammps_dump_with_warnings(
    path: &Path,
    units: Option<LammpsUnits>,
) -> Result<(Trajectory, Vec<Notice>)> {
    let path_ = path.display();
    parse_dump_lines(super::util::text_lines(path)?, units).with_context(|| format!("parsing {path_}"))
}

/// Replaces the placeholder `X<type>` that a dump with no element column reads as, by `map`
/// (type number → element symbol).
///
/// Only this reader's own `X<type>` is taken: a dump that already names elements in an element
/// or typelabel column is an error, so two sources never overwrite each other; so are a type the
/// map leaves out and a symbol that is no element. The messages name no option: the CLI's
/// `--type-map` and Python's `type_map=` reach this alike (`dev/issues.md`, scan B-7)
pub fn apply_type_map(traj: &mut Trajectory, map: &[(u32, String)]) -> Result<()> {
    for (_, sym) in map {
        anyhow::ensure!(ferro_core::data::elements::by_symbol(sym).is_some(),
                        "type map: unknown element symbol {sym:?}");
    }
    for frame in &mut traj.frames {
        for atom in frame.sites_mut() {
            let tp: u32 = atom.element.strip_prefix('X').and_then(|n| n.parse().ok())
                .with_context(|| format!(
                    "a type map was given, but the dump already names elements ({:?}) in its element or typelabel column",
                    atom.element))?;
            let (_, sym) = map.iter().find(|(t, _)| *t == tp)
                .with_context(|| format!("atom type {tp} is missing from the type map"))?;
            atom.element = sym.clone();
        }
    }
    Ok(())
}

#[cfg(test)]
fn parse_lammps_dump(content: &str, units: Option<LammpsUnits>) -> Result<(Trajectory, Vec<Notice>)> {
    parse_dump_lines(content.lines().map(|l| Ok(l.to_string())), units)
}

/// Parses a dump one frame at a time: only the current frame's lines are held, never
/// the whole file (review C-1 — a 327 MB dump used to cost its size twice over in
/// text before a single atom was built).
///
/// A frame runs from its `ITEM: TIMESTEP` line to the line before the next one, and
/// [`parse_frame`] sees exactly that slice, so "is another frame after this one"
/// becomes "did this frame end at a TIMESTEP line rather than at the end of the file".
fn parse_dump_lines(
    mut lines: impl Iterator<Item = Result<String>>,
    units: Option<LammpsUnits>,
) -> Result<(Trajectory, Vec<Notice>)> {
    let is_step = |l: &str| l.trim().starts_with("ITEM: TIMESTEP");
    let mut traj = Trajectory::new();
    let mut warnings = Vec::new();
    // element 列原始串 → 拆分后的元素符号；仅用于读取结束时打印一次映射表
    let mut site_map: BTreeMap<String, String> = BTreeMap::new();
    let mut unknown_prefixes: BTreeSet<String> = BTreeSet::new();

    // Skip to the first frame; `start` is the 0-based line number of its TIMESTEP line
    let mut start = 0usize;
    let mut next = loop {
        match lines.next().transpose()? {
            None => return Ok((traj, warnings)),
            Some(l) if is_step(&l) => break Some(l),
            Some(_) => start += 1,
        }
    };

    while let Some(head) = next.take() {
        let mut frame_lines = vec![head];
        for l in lines.by_ref() {
            let l = l?;
            if is_step(&l) {
                next = Some(l);
                break;
            }
            frame_lines.push(l);
        }
        let view: Vec<&str> = frame_lines.iter().map(|l| l.as_str()).collect();
        match parse_frame(&view, &mut 0, next.is_none(), units, &mut site_map, &mut unknown_prefixes) {
            Ok(frame) => traj.add_frame(frame),
            Err(FrameError::Invalid(why)) => bail!("frame {} (line {}): {why}", traj.n_frames(), start + 1),
            Err(FrameError::Incomplete(why)) => {
                let k = traj.n_frames();
                // 不完整的帧后面还有帧 = 文件中间坏了，报错；只有末帧不完整才是
                // MD 被中断的常态（dump 写到一半），丢掉它，前面的帧照常可用
                if next.is_some() {
                    bail!("frame {k} (line {}) is incomplete: {why}", start + 1);
                }
                super::util::drop_truncated(&mut warnings, k, start + 1, &why)?;
                break;
            }
        }
        start += frame_lines.len();
    }

    warnings.extend(site_map_warnings(&site_map, &unknown_prefixes));
    Ok((traj, warnings))
}

/// 一帧读不成的两种原因，调用者对它们的处置不同。
enum FrameError {
    /// 写到一半被截断：末帧是 MD 中断的常态（丢帧 + 告警），中间帧才报错
    Incomplete(String),
    /// 文件本身写得不对（缺坐标列、字段不是数）：截断解释不了，恒报错
    Invalid(String),
}

impl From<String> for FrameError {
    fn from(s: String) -> Self { FrameError::Incomplete(s) }
}

impl From<&str> for FrameError {
    fn from(s: &str) -> Self { FrameError::Incomplete(s.to_string()) }
}

/// Parse the frame whose `ITEM: TIMESTEP` line is at `*i`, advancing `*i` past it.
///
/// Whether an `Err` is fatal is the caller's call: see [`FrameError`].
fn parse_frame(
    lines: &[&str],
    i: &mut usize,
    last_in_file: bool,
    units: Option<LammpsUnits>,
    site_map: &mut BTreeMap<String, String>,
    unknown_prefixes: &mut BTreeSet<String>,
) -> std::result::Result<Frame, FrameError> {
    // 找本帧的下一个段头；先撞上下一帧的 TIMESTEP 或到文件尾，都说明本帧缺这一段
    let seek = |from: usize, tag: &str| -> std::result::Result<usize, String> {
        for (k, l) in lines.iter().enumerate().skip(from) {
            if l.contains(tag) { return Ok(k); }
            if l.trim().starts_with("ITEM: TIMESTEP") { break; }
        }
        Err(format!("no ITEM: {tag} section"))
    };

    // TIMESTEP 的值紧跟段头；时间相关分析靠它查重复帧与不等间隔
    let step: i64 = lines.get(*i + 1)
        .and_then(|l| l.trim().parse().ok())
        .ok_or("ITEM: TIMESTEP has no step number")?;

    // NUMBER OF ATOMS
    let na = seek(*i + 1, "NUMBER OF ATOMS")?;
    let n: usize = lines.get(na + 1)
        .and_then(|l| l.trim().parse().ok())
        .ok_or("NUMBER OF ATOMS has no count")?;

    // BOX BOUNDS
    let bb = seek(na + 2, "BOX BOUNDS")?;
    let is_triclinic = lines[bb].contains("xy");
    // 边界标志 `pp pp ff`：首字母 p 为周期，f/s/m 为非周期。旧版 LAMMPS 不写标志，
    // 按周期处理（与此前的行为一致）
    let flags: Vec<&str> = lines[bb].split_whitespace()
        .filter(|t| t.len() == 2 && t.chars().all(|c| "pfsm".contains(c)))
        .collect();
    let pbc = match flags.as_slice() {
        [x, y, z] => [x.starts_with('p'), y.starts_with('p'), z.starts_with('p')],
        _ => [true; 3],
    };

    // General triclinic (`dump_modify triclinic/general yes`): the header is `BOX BOUNDS abc origin`,
    // then `ax ay az ox` / `bx by bz oy` / `cx cy cz oz`. Coordinates and vector columns such as
    // velocities and forces are written in that frame (LAMMPS dump / dump_modify manuals) and are
    // taken as written; scaled coordinates are still relative to A B C
    let general = {
        let t: Vec<&str> = lines[bb].split_whitespace().collect();
        t.contains(&"abc") && t.contains(&"origin")
    };
    let (mut cell, origin) = if general {
        let mut rows = [[0.0_f64; 4]; 3];
        for (dim, row) in rows.iter_mut().enumerate() {
            let vals: Vec<f64> = lines.get(bb + 1 + dim)
                .map(|l| l.split_whitespace().map_while(|s| s.parse().ok()).collect())
                .unwrap_or_default();
            let Ok(v) = <[f64; 4]>::try_from(vals.get(..4).unwrap_or_default()) else {
                return Err(format!("BOX BOUNDS abc origin line {} of 3 is missing or cut short", dim + 1).into());
            };
            *row = v;
        }
        let [a, b, c] = rows;
        let m = Matrix3::new(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]);
        (Cell::from_matrix(m), Vector3::new(a[3], b[3], c[3]))
    } else {
        let mut lo = [0.0_f64; 3];
        let mut hi = [0.0_f64; 3];
        let mut tilt = [0.0_f64; 3]; // xy, xz, yz

        for dim in 0..3 {
            let vals: Vec<f64> = lines.get(bb + 1 + dim)
                .map(|l| l.split_whitespace().map_while(|s| s.parse().ok()).collect())
                .unwrap_or_default();
            if vals.len() < 2 {
                return Err(format!("BOX BOUNDS line {} of 3 is missing or cut short", dim + 1).into());
            }
            lo[dim] = vals[0];
            hi[dim] = vals[1];
            if vals.len() >= 3 {
                tilt[dim] = vals[2]; // xy on dim 0, xz on dim 1, yz on dim 2
            }
        }

        let (xy, xz, yz) = (tilt[0], tilt[1], tilt[2]);
        // For a triclinic box these six numbers are the *_bound (the box around the tilted cell); the
        // real edges take the tilts back off:
        //   lx = (xhi_b - xlo_b) - |xy| - |xz|,  ly = (yhi_b - ylo_b) - |yz|
        // this equals the specification's MAX(0,xy,xz,xy+xz) - MIN(...) in all four quadrants,
        // as in ase/io/lammpsrun.py::construct_cell. Orthogonal, the three tilts are 0 and it
        // reduces to hi - lo, so orthogonal trajectories read bit for bit as before
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

        // Scaled coordinates are restored by LAMMPS's definition x = lo + s·L, lo being the box's real
        // origin, not the *_bound, so one atom written as x or as xs reads the same (ASE 3.29 does not
        // add lo to xs, inconsistent with its own x)
        let origin = Vector3::new(
            lo[0] - 0.0_f64.min(xy).min(xz).min(xy + xz),
            lo[1] - 0.0_f64.min(yz),
            lo[2],
        );
        (cell, origin)
    };
    // The origin travels with the box and the writer writes it back as lo (review A-12: losing it put
    // atoms outside the box after a round trip)
    cell.origin = origin;

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

    let dynamic: Vec<&str> = ["vx", "vy", "vz", "fx", "fy", "fz"].into_iter()
        .filter(|c| col.contains_key(c))
        .collect();
    if units.is_none() && !dynamic.is_empty() {
        return Err(FrameError::Invalid(format!(
            "dump has velocity/force columns ({}) whose unit depends on the LAMMPS `units` \
             command, which a dump does not record; give it explicitly \
             (CLI: --units real|metal, Python: units=\"real\"|\"metal\")",
            dynamic.join(" ")
        )));
    }

    // 坐标列按 x > xs > xu > xsu 取第一组完整的（同 ASE 的优先级）；四组都不全就
    // 没有坐标可读，报错而不是让坐标静默为 0
    let (xyz, scaled) = [("x", "y", "z", false), ("xs", "ys", "zs", true),
                         ("xu", "yu", "zu", false), ("xsu", "ysu", "zsu", true)]
        .into_iter()
        .find_map(|(a, b, c, s)| Some(([get_col(a)?, get_col(b)?, get_col(c)?], s)))
        .ok_or_else(|| FrameError::Invalid(
            "ITEM: ATOMS has no complete coordinate columns (x y z, xs ys zs, xu yu zu or xsu ysu zsu)".into()))?;

    let mut atoms_raw: Vec<AtomRaw> = Vec::new();

    for k in 0..n {
        // 文件尾、下一个段头、空行都说明原子段没写完；列数不足是最后一行被截在中间
        let line = lines.get(ah + 1 + k)
            .filter(|l| !l.trim().is_empty() && !l.trim_start().starts_with("ITEM:"))
            .ok_or_else(|| format!("only {k} of {n} atom lines"))?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        // 截断只会切在文件最后一行（常切在数字中间）；别处的列数不足或坏字段是文件写错了
        let at_eof = last_in_file && ah + 1 + k == lines.len() - 1;
        let fail = |msg: String| if at_eof { FrameError::Incomplete(msg) } else { FrameError::Invalid(msg) };
        if parts.len() < col_names.len() {
            return Err(fail(format!(
                "atom line {} of {n} has {} of {} columns",
                k + 1, parts.len(), col_names.len()
            )));
        }
        // 读不出来就报错并点名，不把「解析失败」伪装成「测到了 0」
        let bad = |c: usize| fail(format!("atom line {} of {n}: {} = '{}' is not a number",
                                          k + 1, col_names[c], parts[c]));
        let num = |c: usize| parts[c].parse::<f64>().map_err(|_| bad(c));
        let int = |c: usize| parts[c].parse::<usize>().map_err(|_| bad(c));
        let vec3 = |[cx, cy, cz]: [usize; 3]| -> std::result::Result<Vector3<f64>, FrameError> {
            Ok(Vector3::new(num(cx)?, num(cy)?, num(cz)?))
        };
        let triple = |a, b, c| Some([get_col(a)?, get_col(b)?, get_col(c)?]);

        let atom_id = get_col("id").map(int).transpose()?.unwrap_or(0);
        let tp = get_col("type").map(int).transpose()?.unwrap_or(1);

        // element 列可能写的是位点类型标签（`O_b_P_P`、`Zn_f`），按第一个下划线
        // into element + label; a plain symbol without one passes as is. Without an element column the
        // typelabel column (LAMMPS type labels, set by labelmap) is used, split the same way
        let raw = get_col("element").or_else(|| get_col("typelabel"))
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

        let pos = if scaled {
            origin + cell.fractional_to_cartesian(vec3(xyz)?)
        } else {
            vec3(xyz)?
        };

        let mut atom = Atom::new(element, pos);
        atom.label = label;
        atom.charge = get_col("q").map(num).transpose()?;
        // LAMMPS 的 mass 在 real 与 metal 下都是 g/mol，即 amu
        atom.mass = get_col("mass").map(num).transpose()?;

        // Velocity: real Å/fs (internal), metal Å/ps. A velocity in Å/ps is per ps, so it becomes
        // Å/fs by the factor that turns 1 fs into ps (×1e-3). With units None a file with velocity or
        // force columns was refused above, so the None arm is never used here
        let vscale = match units {
            Some(LammpsUnits::Metal) => convert_time(1.0, TimeUnit::Femtosecond, TimeUnit::Picosecond),
            _ => 1.0,
        };
        let vel = triple("vx", "vy", "vz").map(vec3).transpose()?.map(|v| v * vscale);

        // Force: real kcal/(mol·Å) → eV/Å, metal eV/Å already
        let fscale = match units {
            Some(LammpsUnits::Metal) => 1.0,
            _ => KCAL_MOL_TO_EV,
        };
        let force = triple("fx", "fy", "fz").map(vec3).transpose()?.map(|f| f * fscale);

        atoms_raw.push((atom_id, atom, vel, force));
    }

    // Sort by atom id
    atoms_raw.sort_by_key(|(id, _, _, _)| *id);

    let mut frame = Frame::with_cell(cell, pbc);
    frame.step = Some(step);
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

    // One vector per atom line read above, so the lengths always match
    if has_vel { frame.set_velocities(Some(all_vels)).expect("one velocity per atom"); }
    if has_force { frame.set_forces(Some(all_forces)).expect("one force per atom"); }

    *i = ah + 1 + n;
    Ok(frame)
}

/// Result of splitting a site string into a chemical element and an optional type label.
///
/// See [`split_element_label`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LabelSplit<'a> {
    /// No underscore — the whole string is the element symbol.
    Plain(&'a str),
    /// Underscore present and the prefix is a real element symbol.
    Split { element: &'a str, label: &'a str },
    /// Underscore present but the prefix is not a known element symbol.
    ///
    /// The caller should warn and fall back to treating the whole string as the element.
    Unknown(&'a str),
}

/// Split a site string into `(element, label)` at the **first** underscore.
///
/// Convention for pseudo-element / site-type labels: `<Element>_<suffix>`, where the
/// suffix may itself contain underscores (`O_b_P_P` → element `O`, label `O_b_P_P`).
///
/// ```text
/// "O"          → Plain("O")                                   element O, no label
/// "Zn_f"       → Split { element: "Zn", label: "Zn_f" }        element Zn, label Zn_f
/// "O_b_P_P"    → Split { element: "O",  label: "O_b_P_P" }     element O,  label O_b_P_P
/// "foo_bar"    → Unknown("foo_bar")                            warn; treat whole as element
/// ```
///
/// Splitting on an explicit separator — rather than greedily matching the longest
/// element prefix — avoids silently reading `Pb` as lead(82) when "P bridging" was
/// meant. Plain trajectories (`O`, `P`, `Zn`) are unaffected.
///
/// Note this is deliberately *not* the rule used by the CIF / CP2K / QE readers, whose
/// native site names (`Fe1`, `O2`) carry no separator and are correctly handled by an
/// alphabetic-prefix rule instead.
fn split_element_label(s: &str) -> LabelSplit<'_> {
    match s.split_once('_') {
        None => LabelSplit::Plain(s),
        Some((prefix, _)) if by_symbol(prefix).is_some() => {
            LabelSplit::Split { element: prefix, label: s }
        }
        Some(_) => LabelSplit::Unknown(s),
    }
}

/// The raw-string → element mapping, once per file, so a silently wrong split
/// is visible rather than buried in the results.
///
/// Empty for ordinary trajectories where every element column entry is already
/// a plain chemical symbol.
fn site_map_warnings(site_map: &BTreeMap<String, String>, unknown: &BTreeSet<String>) -> Vec<Notice> {
    let mut out = Vec::new();
    let split: Vec<String> = site_map.iter()
        .filter(|(raw, elem)| raw.as_str() != elem.as_str())
        .map(|(raw, elem)| format!("{raw} -> {elem}"))
        .collect();
    if !split.is_empty() {
        // One line, as every message is: the label is the raw string, so `raw -> element` says it all
        out.push(Notice::note(None, format!(
            "element / typelabel column carries site labels, split into element + label (label = the whole \
             string): {}", split.join(", "))));
    }
    if !unknown.is_empty() {
        let list: Vec<&str> = unknown.iter().map(|s| s.as_str()).collect();
        out.push(Notice::warning(None, format!(
            "underscore present but prefix is not a known element, \
             kept as element verbatim: {}",
            list.join(", ")
        )));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_split_element_label_plain_passthrough() {
        // a plain trajectory (no underscore) is unaffected
        for s in ["O", "P", "Zn", "Fe", "Al"] {
            assert_eq!(split_element_label(s), LabelSplit::Plain(s), "{s}");
        }
    }

    #[test]
    fn test_split_element_label_underscore_forms() {
        assert_eq!(
            split_element_label("Zn_f"),
            LabelSplit::Split { element: "Zn", label: "Zn_f" }
        );
        assert_eq!(
            split_element_label("P_0"),
            LabelSplit::Split { element: "P", label: "P_0" }
        );
        // an underscore inside the suffix: split at the first one only
        assert_eq!(
            split_element_label("O_b_P_P"),
            LabelSplit::Split { element: "O", label: "O_b_P_P" }
        );
        assert_eq!(
            split_element_label("O_n_P"),
            LabelSplit::Split { element: "O", label: "O_n_P" }
        );
    }

    #[test]
    fn test_split_element_label_unknown_prefix() {
        // the prefix is no element: the caller warns and takes the whole string
        assert_eq!(split_element_label("foo_bar"), LabelSplit::Unknown("foo_bar"));
        assert_eq!(split_element_label("On_P"), LabelSplit::Unknown("On_P"));
        assert_eq!(split_element_label("_x"), LabelSplit::Unknown("_x"));
    }

    #[test]
    fn test_split_element_label_no_greedy_misread() {
        // the point of the separator rule: "P bridging" is never misread as a real element
        assert_eq!(split_element_label("Pb"), LabelSplit::Plain("Pb"));   // lead, not split
        assert_eq!(split_element_label("Po"), LabelSplit::Plain("Po"));   // polonium, not split
        assert_eq!(
            split_element_label("P_b"),
            LabelSplit::Split { element: "P", label: "P_b" }
        );
    }

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

    fn parse_traj(content: &str, units: Option<LammpsUnits>) -> Result<Trajectory> {
        Ok(parse_lammps_dump(content, units)?.0)
    }

    /// The reader reports, it does not print: a dropped last frame and the site-label
    /// split come back as warnings, and an ordinary dump yields none.
    #[test]
    fn warnings_are_returned_not_printed() {
        let (_, plain) = parse_lammps_dump(DUMP_ORTHO, None).unwrap();
        assert!(plain.is_empty(), "{plain:?}");

        let cut: String = DUMP_ORTHO.lines().take(21).map(|l| format!("{l}\n")).collect();
        let (traj, w) = parse_lammps_dump(&cut, None).unwrap();
        assert_eq!(traj.n_frames(), 1);
        assert!(w.len() == 1 && w[0].message.contains("last frame (frame 1"), "{w:?}");

        let labelled = DUMP_ORTHO.replace(" 1 Fe ", " 1 Fe_a ");
        let (_, w) = parse_lammps_dump(&labelled, None).unwrap();
        assert!(w.len() == 1 && w[0].level == crate::Level::Note
            && w[0].message.contains("Fe_a -> Fe"), "{w:?}");
        assert!(!w[0].message.contains("-a/-b/-c"), "CLI option names stay out of the library: {w:?}");
    }

    #[test]
    fn test_type_map_on_dump_without_element_column() {
        let no_elem = DUMP_ORTHO.replace("id type element x y z", "id type x y z")
            .replace(" 1 Fe ", " 1 ").replace("2 1 1.4", "2 2 1.4");
        let read = || parse_traj(&no_elem, None).unwrap();
        let t = read();
        assert_eq!((t.frames[0].element(0), t.frames[0].element(1)),
                   ("X1", "X2"), "无映射时仍是占位符");

        let mut t = read();
        apply_type_map(&mut t, &[(1, "Zr".into()), (2, "O".into()), (3, "P".into())]).unwrap();
        for f in &t.frames {
            assert_eq!((f.atom(0).element, f.atom(1).element), ("Zr", "O"), "每帧都要换");
        }

        let err = |map: &[(u32, String)], t: &mut Trajectory| format!("{:#}", apply_type_map(t, map).unwrap_err());
        assert!(err(&[(1, "Zr".into())], &mut read()).contains("type 2 is missing"));
        assert!(err(&[(1, "Zr".into()), (2, "Xx".into())], &mut read()).contains("unknown element"));
        let mut named = parse_traj(DUMP_ORTHO, None).unwrap();
        assert!(err(&[(1, "Zr".into())], &mut named).contains("already names elements"),
                "有 element 列时不应被映射覆盖");
    }

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
        let traj = read_lammps_dump(&tmp("labels.dump", DUMP_LABELS), Some(LammpsUnits::Real)).unwrap();
        let f = traj.first().unwrap();
        let got: Vec<(&str, Option<&str>)> = f.atoms()
            .map(|a| (a.element, a.label))
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
        let traj = read_lammps_dump(&tmp("plain.dump", DUMP_ORTHO), Some(LammpsUnits::Real)).unwrap();
        for atom in traj.first().unwrap().atoms() {
            assert_eq!(atom.element, "Fe");
            assert!(atom.label.is_none(), "plain symbol should not produce a label");
        }
    }

    #[test]
    fn test_unknown_prefix_kept_verbatim_as_element() {
        let traj = read_lammps_dump(&tmp("bad.dump", DUMP_BAD_LABEL), Some(LammpsUnits::Real)).unwrap();
        let f = traj.first().unwrap();
        // Unknown prefix: the whole string is the element, label is None (the reader warns)
        assert_eq!(f.atom(0).element, "foo_bar");
        assert!(f.atom(0).label.is_none());
        assert_eq!(f.atom(1).element, "O");
    }

    #[test]
    fn test_multiframe() {
        let traj = read_lammps_dump(&tmp("bcc.dump", DUMP_ORTHO), Some(LammpsUnits::Real)).unwrap();
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
            let traj = parse_traj(&text, Some(LammpsUnits::Real))
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
        let err = parse_traj(&lines.join("\n"), Some(LammpsUnits::Real))
            .expect_err("中间帧不完整应报错");
        let msg = format!("{err:#}");
        assert!(msg.contains("frame 0") && msg.contains("1 of 2 atom lines"), "报错应点名帧与行数：{msg}");
    }

    #[test]
    fn test_box_bounds() {
        let traj = read_lammps_dump(&tmp("box.dump", DUMP_ORTHO), Some(LammpsUnits::Real)).unwrap();
        let [a, ..] = traj.first().unwrap().cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-6);
    }

    #[test]
    fn test_velocity_real_units() {
        let traj = read_lammps_dump(&tmp("vel_real.dump", DUMP_VEL), Some(LammpsUnits::Real)).unwrap();
        let vx = traj.first().unwrap().velocities().unwrap()[0].x;
        assert!((vx - 2.0).abs() < 1e-10, "real units: vx should be 2.0 Å/fs, got {vx}");
    }

    #[test]
    fn test_timestep_is_kept() {
        let traj = read_lammps_dump(&tmp("steps.dump", DUMP_ORTHO), None).unwrap();
        let steps: Vec<_> = traj.frames.iter().map(|f| f.step).collect();
        assert!(steps.iter().all(|s| s.is_some()), "每帧都应带上 TIMESTEP，实际 {steps:?}");
    }

    #[test]
    fn test_dynamic_columns_require_units() {
        // 有速度列而不给单位：两种体系差 1000 倍，不能猜
        let err = read_lammps_dump(&tmp("vel_none.dump", DUMP_VEL), None).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("vx") && msg.contains("--units"), "错误信息应点名速度列与 --units，实际：{msg}");
        // 只有坐标的 dump 不受单位影响，None 照常读
        let traj = read_lammps_dump(&tmp("plain_none.dump", DUMP_ORTHO), None).unwrap();
        assert!(traj.first().unwrap().velocities().is_none());
    }

    #[test]
    fn test_velocity_metal_units() {
        let traj = read_lammps_dump(&tmp("vel_metal.dump", DUMP_VEL), Some(LammpsUnits::Metal)).unwrap();
        let vx = traj.first().unwrap().velocities().unwrap()[0].x;
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
        let traj = read_lammps_dump(&tmp("tri.dump", TRI), Some(LammpsUnits::Metal)).unwrap();
        let m = traj.first().unwrap().cell.as_ref().unwrap().matrix;
        // 行优先：行 = 晶格矢量
        let want = [[10.0, 0.0, 0.0], [2.0, 12.0, 0.0], [-3.0, 1.0, 14.0]];
        for (i, row) in want.iter().enumerate() {
            for (j, w) in row.iter().enumerate() {
                assert!((m[(i, j)] - w).abs() < 1e-10, "cell[{i}][{j}]: got {}, want {w}", m[(i, j)]);
            }
        }
    }

    /// 三斜、原点不在 0 的盒子：真实 lo = (1, 2, -1)，a = (10,0,0)、b = (2,8,0)、
    /// c = (-1,1.5,6)；BOX BOUNDS 写的是 *_bound（xlo_b = 1 + min(0,2,-1,1) = 0）
    fn tri_dump(cols: &str, rows: [&str; 2]) -> String {
        format!("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n\
                 ITEM: BOX BOUNDS xy xz yz pp pp pp\n0 13 2\n2 11.5 -1\n-1 5 1.5\n\
                 ITEM: ATOMS id type element {cols}\n1 1 O {}\n2 1 O {}\n", rows[0], rows[1])
    }

    #[test]
    fn test_four_coordinate_kinds_give_the_same_positions() {
        // 期望值 = lo + s·cell（numpy 算），即 LAMMPS 对 xs 的定义 xs = (x - lo)/L 的逆
        let want = [[2.1, 4.05, 0.8], [10.25, 7.125, 3.5]];
        let cart = ["2.1 4.05 0.8", "10.25 7.125 3.5"];
        let frac = ["0.1 0.2 0.3", "0.9 0.5 0.75"];
        for (cols, rows) in [("x y z", cart), ("xs ys zs", frac), ("xu yu zu", cart), ("xsu ysu zsu", frac)] {
            let traj = parse_traj(&tri_dump(cols, rows), Some(LammpsUnits::Real)).unwrap();
            let f = traj.first().unwrap();
            for (k, w) in want.iter().enumerate() {
                let p = f.atom(k).position;
                for d in 0..3 {
                    assert!((p[d] - w[d]).abs() < 1e-9, "{cols} 原子 {k}：{p:?}，应为 {w:?}");
                }
            }
        }
    }

    #[test]
    fn test_missing_coordinates_or_bad_fields_are_errors() {
        let err = |text: String| format!("{:#}", parse_traj(&text, Some(LammpsUnits::Real)).unwrap_err());
        // 没有坐标列、坐标列不全：以前坐标静默为 0
        assert!(err(tri_dump("q", ["0.1", "0.2"])).contains("no complete coordinate columns"));
        assert!(err(tri_dump("x y", ["1 2", "3 4"])).contains("no complete coordinate columns"));
        // 字段读不出来：以前当 0
        assert!(err(tri_dump("x y z", ["2.1 4.05 nan?", "1 2 3"])).contains("z = 'nan?'"));
        assert!(err(tri_dump("x y z vx vy vz", ["1 2 3 0.1 x 0.3", "1 2 3 0 0 0"])).contains("vy = 'x'"));
        assert!(err(tri_dump("x y z q", ["1 2 3 -", "1 2 3 0"])).contains("q = '-'"));
    }

    #[test]
    fn test_mass_column_is_kept() {
        let traj = parse_traj(&tri_dump("x y z mass", ["1 2 3 2.014", "1 2 3 15.999"]), Some(LammpsUnits::Real)).unwrap();
        let f = traj.first().unwrap();
        assert_eq!(f.atom(0).mass, Some(2.014));
        assert_eq!(f.atom(1).mass, Some(15.999));
        // 没有 mass 列就是 None，不补零
        let plain = parse_traj(&tri_dump("x y z", ["1 2 3", "1 2 3"]), Some(LammpsUnits::Real)).unwrap();
        assert_eq!(plain.first().unwrap().atom(0).mass, None);
    }

    #[test]
    fn test_bad_field_is_truncation_only_on_the_last_line_of_the_file() {
        // 截断切在最后一行的数字中间（`1.4e`）：末帧丢弃，前面的帧照常
        let cut = DUMP_ORTHO.replace("2 1 Fe 1.445 1.435 1.435\n", "2 1 Fe 1.445 1.435 1.4e");
        let traj = parse_traj(&cut, Some(LammpsUnits::Real)).unwrap();
        assert_eq!(traj.n_frames(), 1);
        // 同样的坏字段不在最后一行：截断解释不了，即使在末帧也报错
        let bad = DUMP_ORTHO.replace("1 1 Fe 0.01  0.0   0.0", "1 1 Fe 0.01  0.0   0.0e");
        let msg = format!("{:#}", parse_traj(&bad, Some(LammpsUnits::Real)).unwrap_err());
        assert!(msg.contains("frame 1") && msg.contains("z = '0.0e'"), "{msg}");
    }

    #[test]
    fn test_boundary_flags_set_pbc() {
        let with = |flags: &str| {
            let text = DUMP_ORTHO.replace("ITEM: BOX BOUNDS pp pp pp", &format!("ITEM: BOX BOUNDS {flags}"));
            parse_traj(&text, Some(LammpsUnits::Real)).unwrap().frames[0].pbc
        };
        assert_eq!(with("pp pp pp"), [true; 3]);
        assert_eq!(with("pp pp ff"), [true, true, false]);
        assert_eq!(with("pp ss fm"), [true, false, false]);
        assert_eq!(with("xy xz yz pp fs pp"), [true, false, true]);
        // 旧版 LAMMPS 不写标志：按周期
        assert_eq!(with(""), [true; 3]);
    }

    /// `dump_modify triclinic/general yes` (LAMMPS dump manual): the box is three edge
    /// vectors plus the origin, and coordinates are already in that frame. It used to be
    /// read as `lo hi tilt`, giving a nonsense box.
    #[test]
    fn test_general_triclinic_box() {
        // A along +y, B along −x: not lower-triangular, so reading it as restricted goes wrong
        let text = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\n\
                    ITEM: BOX BOUNDS abc origin pp pp ff\n0 4 0 1\n-4 0 0 2\n0 0 4 3\n\
                    ITEM: ATOMS id type element x y z\n1 1 O 1 4 3\n2 1 H 0 2 4\n";
        let (traj, _) = parse_lammps_dump(text, None).unwrap();
        let f = &traj.frames[0];
        let cell = f.cell.as_ref().unwrap();
        assert_eq!(cell.matrix, Matrix3::new(0.0, 4.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 4.0));
        assert_eq!(cell.origin, Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(f.pbc, [true, true, false]);
        assert_eq!(f.positions()[1], Vector3::new(0.0, 2.0, 4.0), "x y z taken as written");
        // Scaled coordinates are relative to A B C: origin + 0.5·A = (1, 4, 3)
        let scaled = text.replace("x y z\n1 1 O 1 4 3", "xs ys zs\n1 1 O 0.5 0 0")
                         .replace("2 1 H 0 2 4", "2 1 H 0.5 0.25 0.25");
        let (traj, _) = parse_lammps_dump(&scaled, None).unwrap();
        assert!((traj.frames[0].positions()[0] - Vector3::new(1.0, 4.0, 3.0)).norm() < 1e-12);
        // One number short is an error (it used to give an empty trajectory and a warning)
        let bad = text.replace("-4 0 0 2", "-4 0 0");
        let err = format!("{:#}", parse_lammps_dump(&bad, None).unwrap_err());
        assert!(err.contains("only frame") && err.contains("abc origin line 2"), "{err}");
    }

    /// A dump with a `typelabel` column and no `element` column used to give `X<type>`.
    #[test]
    fn test_typelabel_column() {
        let text = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n\
                    0 10\n0 10\n0 10\nITEM: ATOMS id typelabel x y z\n1 O_b 0 0 0\n2 P 1 0 0\n";
        let f = parse_lammps_dump(text, None).unwrap().0.frames.remove(0);
        assert_eq!((f.atom(0).element, f.atom(0).label), ("O", Some("O_b")));
        assert_eq!((f.atom(1).element, f.atom(1).label), ("P", None));
        // The element column wins
        let both = text.replace("id typelabel x y z\n1 O_b", "id element typelabel x y z\n1 Si O_b")
                       .replace("2 P 1", "2 P P 1");
        assert_eq!(parse_lammps_dump(&both, None).unwrap().0.frames[0].atom(0).element, "Si");
    }
}
