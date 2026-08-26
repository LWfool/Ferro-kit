use std::collections::HashMap;
use ferro_core::{matrix3_row_major, Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{bail, Context, Result};

pub fn read_extxyz(path: &str) -> Result<Trajectory> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path}"))?;
    parse_extxyz(&content).with_context(|| format!("parsing {path}"))
}

fn parse_extxyz(content: &str) -> Result<Trajectory> {
    let mut lines = content.lines().peekable();
    let mut traj = Trajectory::new();

    loop {
        // Skip blank lines between frames
        while lines.peek().map(|l| l.trim().is_empty()).unwrap_or(false) {
            lines.next();
        }
        let count_line = match lines.next() {
            None => break,
            Some(l) => l.trim(),
        };
        if count_line.is_empty() { break; }

        let n: usize = count_line.parse()
            .with_context(|| format!("expected atom count, got {count_line:?}"))?;

        let comment = lines.next().context("missing comment line")?;
        let kv = parse_comment(comment);

        // Cell and pbc
        let (cell, pbc) = match kv.get("lattice") {
            Some(lat) => {
                let c = parse_lattice(lat).context("invalid Lattice")?;
                let pbc = kv.get("pbc").map(|s| parse_pbc(s)).unwrap_or([true; 3]);
                (Some(c), pbc)
            }
            None => (None, [false; 3]),
        };

        // Scalar properties in comment
        let energy: Option<f64> = kv.get("energy").and_then(|s| s.parse().ok());
        let stress = read_stress(&kv, cell.as_ref())
            .with_context(|| format!("frame {}", traj.n_frames()))?;

        // Properties column spec
        let props = kv.get("properties").map(|s| parse_properties(s))
            .unwrap_or_else(|| vec![
                ("species".to_string(), 'S', 1),
                ("pos".to_string(), 'R', 3),
            ]);
        let ncols: usize = props.iter().map(|(_, _, c)| c).sum();

        // Build column → field map
        let mut col_idx: HashMap<&str, usize> = HashMap::new();
        let mut offset = 0;
        for (name, _, count) in &props {
            col_idx.insert(name.as_str(), offset);
            offset += count;
        }

        let find = |name: &str| col_idx.get(name).copied();

        let mut frame = match cell {
            Some(c) => Frame::with_cell(c, pbc),
            None => Frame::new(),
        };
        frame.energy = energy;
        frame.stress = stress;

        let mut all_forces: Vec<Vector3<f64>> = Vec::new();
        let mut all_vels: Vec<Vector3<f64>> = Vec::new();

        for i in 0..n {
            let line = lines.next()
                .with_context(|| format!("missing atom line {i}"))?;
            let cols: Vec<&str> = line.split_whitespace().collect();
            if cols.len() < ncols {
                bail!("atom line {i}: expected {ncols} columns, got {}", cols.len());
            }

            let element = find("species")
                .map(|c| cols[c].to_string())
                .unwrap_or_else(|| "X".to_string());

            let pos = if let Some(c) = find("pos") {
                let x: f64 = cols[c].parse().context("invalid pos.x")?;
                let y: f64 = cols[c+1].parse().context("invalid pos.y")?;
                let z: f64 = cols[c+2].parse().context("invalid pos.z")?;
                Vector3::new(x, y, z)
            } else {
                Vector3::zeros()
            };

            let mut atom = Atom::new(element, pos);
            if let Some(c) = find("charges") {
                atom.charge = cols[c].parse().ok();
            }
            if let Some(c) = find("masses") {
                atom.mass = cols[c].parse().ok();
            }
            if let Some(c) = find("magmoms") {
                atom.magmom = cols[c].parse().ok();
            }
            // 位点标签走自己的一列,species 保持纯元素 —— 不像 LAMMPS dump
            // 那样需要按下划线拆分,故这里不做任何猜测
            if let Some(c) = find("label") {
                atom.label = Some(cols[c].to_string());
            }
            frame.add_atom(atom);

            if let Some(c) = find("forces").or_else(|| find("force")) {
                let fx: f64 = cols[c].parse().unwrap_or(0.0);
                let fy: f64 = cols[c+1].parse().unwrap_or(0.0);
                let fz: f64 = cols[c+2].parse().unwrap_or(0.0);
                all_forces.push(Vector3::new(fx, fy, fz));
            }
            if let Some(c) = find("velocities").or_else(|| find("momenta")) {
                let vx: f64 = cols[c].parse().unwrap_or(0.0);
                let vy: f64 = cols[c+1].parse().unwrap_or(0.0);
                let vz: f64 = cols[c+2].parse().unwrap_or(0.0);
                all_vels.push(Vector3::new(vx, vy, vz));
            }
        }

        if all_forces.len() == n { frame.forces = Some(all_forces); }
        if all_vels.len() == n { frame.velocities = Some(all_vels); }

        if traj.n_frames() == 0 {
            if let Some(src) = kv.get("config_type").or_else(|| kv.get("comment")) {
                traj.metadata.source = Some(src.clone());
            }
        }
        traj.add_frame(frame);
    }

    Ok(traj)
}

// ─── Comment line parser ──────────────────────────────────────────────────────

fn parse_comment(line: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    let mut s = line.trim();
    while !s.is_empty() {
        // Find next '='
        let eq = match s.find('=') {
            Some(p) => p,
            None => break,
        };
        let key = s[..eq].trim().to_lowercase();
        s = &s[eq + 1..];
        let (val, rest) = read_value(s);
        if !key.is_empty() { map.insert(key, val); }
        s = rest.trim_start();
    }
    map
}

fn read_value(s: &str) -> (String, &str) {
    let s = s.trim_start();
    if let Some(q) = s.chars().next().filter(|&c| c == '"' || c == '\'') {
        let inner = &s[1..];
        // Find closing quote (not escaped)
        let end = inner.find(q).unwrap_or(inner.len());
        (inner[..end].to_string(), &inner[end + 1..])
    } else {
        let end = s.find(char::is_whitespace).unwrap_or(s.len());
        (s[..end].to_string(), &s[end..])
    }
}

fn parse_lattice(s: &str) -> Result<Cell> {
    let v: Vec<f64> = s.split_whitespace()
        .map(|x| x.parse::<f64>().context("float"))
        .collect::<Result<_>>()?;
    anyhow::ensure!(v.len() == 9, "Lattice must have 9 values");
    Ok(Cell::from_matrix(Matrix3::new(
        v[0], v[1], v[2],
        v[3], v[4], v[5],
        v[6], v[7], v[8],
    )))
}

fn parse_pbc(s: &str) -> [bool; 3] {
    let mut pbc = [false; 3];
    for (i, tok) in s.split_whitespace().take(3).enumerate() {
        pbc[i] = matches!(tok.to_uppercase().as_str(), "T" | "TRUE" | "1");
    }
    pbc
}

/// Relative tolerance for the symmetry check and the stress/virial cross-check.
/// Both compare numbers that were rounded on their way into text.
const TENSOR_TOL: f64 = 1e-6;

/// Reads `stress=` / `virial=` into the convention [`ferro_core::Frame::stress`]
/// uses (eV/Å³, positive = compression).
///
/// extxyz carries the tensor under two keys that mean different things:
///
/// | key | file holds | conversion |
/// |---|---|---|
/// | `stress` | eV/Å³, positive = tension | negate |
/// | `virial` | eV, positive = compression | divide by the cell volume |
///
/// Three independent sources agree on this: the extxyz specification states
/// that `virial -> stress` is a multiplication by `-1/cell_vol`; the GPUMD
/// manual documents `virial` as positive-for-compressed and `stress` as
/// positive-for-stretched; dpdata 1.0.2 writes `virials = -volume * stress`
/// against ASE's stress (`dpdata/plugins/ase.py`).
///
/// A file carrying both keys is cross-checked rather than resolved by
/// precedence — the two disagreeing is a defect in the file, and picking a
/// winner would hide it.
fn read_stress(
    kv: &HashMap<String, String>,
    cell: Option<&Cell>,
) -> Result<Option<Matrix3<f64>>> {
    let stress = kv.get("stress").map(|v| parse_tensor9("stress", v)).transpose()?;
    let virial = kv.get("virial").map(|v| parse_tensor9("virial", v)).transpose()?;

    // virial 是 eV,要除体积才是应力;没有 Lattice 就没有体积可除
    let from_virial = match virial {
        None => None,
        Some(v) => {
            let cell = cell.context(
                "virial= without a Lattice — no cell volume to divide by")?;
            let vol = cell.volume();
            anyhow::ensure!(vol.abs() > 1e-12, "virial= with a zero-volume cell");
            Some(v / vol)
        }
    };
    let from_stress = stress.map(|s| -s);

    match (from_stress, from_virial) {
        (Some(s), Some(v)) => {
            if !close(&s, &v) {
                let vol = cell.map(|c| c.volume()).unwrap_or(1.0);
                // 行优先渲染 —— nalgebra 的 as_slice() 是列优先,连报错都会打出转置
                bail!(
                    "stress= and virial= disagree: virial/V gives {:?} but -stress gives {:?} \
                     (eV/Å³, row-major, cell volume {vol}). One of the two keys was written \
                     under a different convention; fix the file rather than letting this \
                     reader pick a winner",
                    matrix3_row_major(&v), matrix3_row_major(&s)
                );
            }
            Ok(Some(s))
        }
        (Some(s), None) => Ok(Some(s)),
        (None, Some(v)) => Ok(Some(v)),
        (None, None) => Ok(None),
    }
}

/// Parses the nine numbers of a `stress=` / `virial=` value.
///
/// The extxyz specification requires this tensor to be symmetric ("fail if not
/// symmetric"), and that requirement is what makes the row-major/column-major
/// question moot: ASE documents these nine numbers as Fortran-ordered while the
/// GPUMD manual spells them out row-major (`vxx vxy vxz vyx ...`), and both
/// readings agree on every symmetric tensor. So this parser checks symmetry
/// instead of betting on one of the two descriptions.
///
/// A six-number Voigt value is rejected on purpose: the component order is not
/// universal (the spec and ASE use `xx yy zz yz xz xy`, while VASP's `in kB`
/// line and GPUMD's own `stress_*.out` use `xx yy zz xy yz zx`), and nothing in
/// the file says which one produced it. Guessing would silently transpose the
/// off-diagonal components.
fn parse_tensor9(key: &str, s: &str) -> Result<Matrix3<f64>> {
    let v: Vec<f64> = s
        .split_whitespace()
        .map(|x| x.parse::<f64>().with_context(|| format!("{key}=: {x:?} is not a number")))
        .collect::<Result<_>>()?;

    if v.len() == 6 {
        bail!(
            "{key}= holds 6 numbers (Voigt). Ferro does not read this form: the component \
             order is not universal — the extxyz spec and ASE use `xx yy zz yz xz xy`, \
             VASP's `in kB` line and GPUMD's `stress_*.out` use `xx yy zz xy yz zx` — and \
             the file does not say which one wrote it. Re-emit the tensor as 9 numbers"
        );
    }
    anyhow::ensure!(v.len() == 9, "{key}= holds {} numbers, expected 9", v.len());

    let m = Matrix3::new(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
    let scale = v.iter().fold(0.0f64, |a, b| a.max(b.abs())).max(1.0);
    for (i, j) in [(0, 1), (0, 2), (1, 2)] {
        if (m[(i, j)] - m[(j, i)]).abs() > TENSOR_TOL * scale {
            bail!(
                "{key}= is not symmetric ({},{}) = {} but ({},{}) = {}. The extxyz \
                 specification requires a symmetric tensor; an asymmetric one means the \
                 nine numbers are not what this reader takes them for",
                i, j, m[(i, j)], j, i, m[(j, i)]
            );
        }
    }
    Ok(m)
}

/// Component-wise comparison with [`TENSOR_TOL`], scaled by the larger operand.
fn close(a: &Matrix3<f64>, b: &Matrix3<f64>) -> bool {
    a.iter().zip(b.iter()).all(|(x, y)| {
        (x - y).abs() <= TENSOR_TOL * x.abs().max(y.abs()).max(1.0)
    })
}

fn parse_properties(spec: &str) -> Vec<(String, char, usize)> {
    let parts: Vec<&str> = spec.split(':').collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i + 2 < parts.len() {
        let name = parts[i].to_lowercase();
        let tc = parts[i+1].chars().next().unwrap_or('S').to_ascii_uppercase();
        let count: usize = parts[i+2].parse().unwrap_or(1);
        out.push((name, tc, count));
        i += 3;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const SINGLE: &str = r#"2
Lattice="2.87 0.0 0.0 0.0 2.87 0.0 0.0 0.0 2.87" Properties=species:S:1:pos:R:3 pbc="T T T" energy=-17.43
Fe 0.0    0.0    0.0
Fe 1.435  1.435  1.435
"#;

    const WITH_FORCES: &str = r#"2
Lattice="5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0" Properties=species:S:1:pos:R:3:forces:R:3 pbc="T T T"
Fe 0.0 0.0 0.0 0.1 0.0 0.0
Fe 2.5 2.5 2.5 -0.1 0.0 0.0
"#;

    // ── stress / virial fixtures ─────────────────────────────────────────────
    //
    // 下面两份文本是 ASE 3.29.0 亲手写出的,不是照约定手敲的 —— 符号约定的依据
    // 必须来自外部实物:ferro 写出→ferro 读回的回路里符号错两次会互相抵消,
    // 测试照样全绿。生成方式(~/.miniforge3/envs/deepmd,ase 3.29.0 + dpdata 1.0.2):
    //
    //   cell   = [[10,0,0],[1,11,0],[2,3,12]]        V = 1320 Å³
    //   σ_ase  = [[0.01,0.002,0.003],[0.002,0.02,0.004],[0.003,0.004,0.03]]
    //   at.calc = SinglePointCalculator(at, stress=σ_ase, ...); ase.io.write(...)
    //   virial = -V * σ_ase                          (dpdata 1.0.2 的换算式)
    //
    // 两份的期望值相同:σ_ferro = -σ_ase = virial/V,亦即 dpdata 读同一份 virial
    // 键再除体积得到的数。
    const ASE_STRESS: &str = r#"2
Lattice="10.0 0.0 0.0 1.0 11.0 0.0 2.0 3.0 12.0" Properties=species:S:1:pos:R:3:forces:R:3 energy=-1.5 stress="0.01 0.002 0.003 0.002 0.02 0.004 0.003 0.004 0.03" pbc="T T T"
H        0.00000000       0.00000000       0.00000000       0.10000000       0.20000000       0.30000000
H        1.00000000       1.00000000       1.00000000      -0.10000000      -0.20000000      -0.30000000
"#;

    const ASE_VIRIAL: &str = r#"2
Lattice="10.0 0.0 0.0 1.0 11.0 0.0 2.0 3.0 12.0" Properties=species:S:1:pos:R:3:forces:R:3 virial="-13.200000000000012 -2.6400000000000023 -3.9600000000000035 -2.6400000000000023 -26.400000000000023 -5.280000000000005 -3.9600000000000035 -5.280000000000005 -39.60000000000003" energy=-1.5 pbc="T T T"
H        0.00000000       0.00000000       0.00000000       0.10000000       0.20000000       0.30000000
H        1.00000000       1.00000000       1.00000000      -0.10000000      -0.20000000      -0.30000000
"#;

    /// σ_ferro（正 = 压缩）for both fixtures above.
    const EXPECT: [[f64; 3]; 3] = [
        [-0.01,  -0.002, -0.003],
        [-0.002, -0.02,  -0.004],
        [-0.003, -0.004, -0.03 ],
    ];

    fn assert_expect(f: &Frame) {
        let s = f.stress.expect("stress");
        for i in 0..3 {
            for j in 0..3 {
                assert!((s[(i, j)] - EXPECT[i][j]).abs() < 1e-9,
                    "({i},{j}): {} vs {}", s[(i, j)], EXPECT[i][j]);
            }
        }
    }

    fn err_of(name: &str, text: &str) -> String {
        format!("{:#}", read_extxyz(&tmp(name, text)).unwrap_err())
    }

    fn tmp(name: &str, c: &str) -> String {
        let p = std::env::temp_dir().join(name);
        std::fs::write(&p, c).unwrap();
        p.to_str().unwrap().to_string()
    }

    #[test]
    fn test_basic() {
        let traj = read_extxyz(&tmp("bcc.extxyz", SINGLE)).unwrap();
        assert_eq!(traj.n_frames(), 1);
        let f = traj.first().unwrap();
        assert_eq!(f.n_atoms(), 2);
        assert_eq!(f.atom(0).element, "Fe");
        assert!((f.energy.unwrap() - (-17.43)).abs() < 1e-10);
        let [a, ..] = f.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-6);
    }

    #[test]
    fn test_forces() {
        let traj = read_extxyz(&tmp("forces.extxyz", WITH_FORCES)).unwrap();
        let f = traj.first().unwrap();
        let forces = f.forces.as_ref().unwrap();
        assert!((forces[0].x - 0.1).abs() < 1e-10);
        assert!((forces[1].x - (-0.1)).abs() < 1e-10);
    }
    #[test]
    fn test_stress_sign_from_ase() {
        // ASE 写的 stress= 是正 = 拉伸;Ferro 存正 = 压缩,故对角必须转负
        let traj = read_extxyz(&tmp("ase_stress.extxyz", ASE_STRESS)).unwrap();
        let f = traj.first().unwrap();
        assert_expect(f);
        assert!(f.stress.unwrap()[(0, 0)] < 0.0);
    }

    #[test]
    fn test_virial_matches_dpdata() {
        // virial= 是 eV 且正 = 压缩:除以体积、不变号,结果应与 stress fixture 相同
        let traj = read_extxyz(&tmp("ase_virial.extxyz", ASE_VIRIAL)).unwrap();
        assert_expect(traj.first().unwrap());
    }

    #[test]
    fn test_stress_and_virial_agree() {
        let both = ASE_STRESS.replace(
            "energy=-1.5",
            "energy=-1.5 virial=\"-13.2 -2.64 -3.96 -2.64 -26.4 -5.28 -3.96 -5.28 -39.6\"");
        let traj = read_extxyz(&tmp("both_ok.extxyz", &both)).unwrap();
        assert_expect(traj.first().unwrap());
    }

    #[test]
    fn test_stress_and_virial_conflict() {
        // virial 少一个负号 —— 正是「两处约定不一致」的真实故障形态
        let bad = ASE_STRESS.replace(
            "energy=-1.5",
            "energy=-1.5 virial=\"13.2 2.64 3.96 2.64 26.4 5.28 3.96 5.28 39.6\"");
        let e = err_of("both_bad.extxyz", &bad);
        assert!(e.contains("disagree"), "{e}");
    }

    #[test]
    fn test_virial_without_lattice() {
        let no_cell = "1\nvirial=\"1 0 0 0 1 0 0 0 1\"\nH 0.0 0.0 0.0\n";
        let e = err_of("no_lattice.extxyz", no_cell);
        assert!(e.contains("no cell volume"), "{e}");
    }

    #[test]
    fn test_asymmetric_tensor_rejected() {
        // 非对称张量:规格要求对称,而这也是行/列优先唯一能显形的地方
        let asym = ASE_STRESS.replace(
            "stress=\"0.01 0.002 0.003 0.002 0.02 0.004 0.003 0.004 0.03\"",
            "stress=\"0.01 0.002 0.003 0.009 0.02 0.004 0.003 0.004 0.03\"");
        let e = err_of("asym.extxyz", &asym);
        assert!(e.contains("not symmetric"), "{e}");
    }

    #[test]
    fn test_voigt6_rejected() {
        let v6 = ASE_STRESS.replace(
            "stress=\"0.01 0.002 0.003 0.002 0.02 0.004 0.003 0.004 0.03\"",
            "stress=\"0.01 0.02 0.03 0.004 0.003 0.002\"");
        let e = err_of("voigt6.extxyz", &v6);
        assert!(e.contains("6 numbers"), "{e}");
    }

    #[test]
    fn test_nep_singular_force_column() {
        // GPUMD 的 train.xyz 用 force:R:3,ASE 用 forces:R:3,两种拼法都要收 ——
        // 只认复数时,读 NEP 数据集会一声不吭地把受力全丢了
        let nep = ASE_STRESS.replace("forces:R:3", "force:R:3");
        let traj = read_extxyz(&tmp("nep_force.extxyz", &nep)).unwrap();
        let f = traj.first().unwrap();
        let forces = f.forces.as_ref().expect("force column not read");
        assert!((forces[0].x - 0.1).abs() < 1e-12);
        assert!((forces[1].z + 0.3).abs() < 1e-12);
    }

}
