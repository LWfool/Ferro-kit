use std::path::Path;
use ferro_core::{Atom, Cell, Frame, Trajectory};
use nalgebra::{Matrix3, Vector3};
use anyhow::{ensure, Context, Result};
use super::util::floats;

pub fn read_poscar(path: &Path) -> Result<Trajectory> {
    let path_ = path.display();
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("cannot open {path_}"))?;
    parse_poscar(&content).with_context(|| format!("parsing {path_}"))
}

pub fn read_contcar(path: &Path) -> Result<Trajectory> {
    read_poscar(path)
}

fn parse_poscar(content: &str) -> Result<Trajectory> {
    let mut lines = content.lines();
    let mut next = |what: &str| -> Result<&str> {
        lines.next().with_context(|| format!("unexpected EOF before {what}"))
    };

    let comment = next("comment")?.trim().to_string();

    let scale: f64 = next("scaling factor")?
        .trim().parse().context("invalid scaling factor")?;
    ensure!(scale > 0.0, "negative scaling factor (volume-based) is not supported");

    // Lattice vectors — rows of the Cell matrix
    let mut m = [[0.0_f64; 3]; 3];
    for (i, row) in m.iter_mut().enumerate() {
        let v = floats(next(&format!("lattice vector {i}"))?, 3)
            .with_context(|| format!("invalid lattice vector {i}"))?;
        *row = [v[0] * scale, v[1] * scale, v[2] * scale];
    }
    let cell = Cell::from_matrix(Matrix3::new(
        m[0][0], m[0][1], m[0][2],
        m[1][0], m[1][1], m[1][2],
        m[2][0], m[2][1], m[2][2],
    ));

    // VASP4 vs VASP5 detection
    let line5 = next("element/count line")?.trim();
    // Parse line5 once as usize counts: Some ⇒ VASP4, None ⇒ VASP5.
    // (Parsing directly as usize avoids the prior u64-check / usize-unwrap
    // mismatch that could panic on 32-bit targets for large values.)
    let line5_counts: Option<Vec<usize>> =
        line5.split_whitespace().map(|t| t.parse::<usize>().ok()).collect();
    let (elements, counts): (Vec<String>, Vec<usize>) =
        if let Some(c) = line5_counts {
            // VASP4: counts only
            let e = (1..=c.len()).map(|i| format!("X{i}")).collect();
            (e, c)
        } else {
            // VASP5: element symbols + counts on next line
            let e: Vec<String> = line5.split_whitespace().map(|s| s.to_string()).collect();
            let c: Vec<usize> = next("atom counts")?
                .split_whitespace()
                .map(|s| s.parse().context("invalid count"))
                .collect::<Result<_>>()?;
            ensure!(e.len() == c.len(), "element/count length mismatch");
            (e, c)
        };

    // Optional "Selective dynamics"
    let mut coord_type = next("coordinate type")?;
    if coord_type.trim().to_lowercase().starts_with('s') {
        coord_type = next("coordinate type after Selective dynamics")?;
    }
    // VASP 的规则：首字符 C/c/K/k 为 Cartesian，其余（含 `Fractional`、空行）一律 Direct。
    // 以前反过来只认 d 开头为 Direct，`Fractional` 之类被当 Cartesian 读
    let is_direct = !coord_type.trim_start().starts_with(['C', 'c', 'K', 'k']);

    // Atom positions
    let mut frame = Frame::with_cell(cell.clone(), [true; 3]);
    let total: usize = counts.iter().sum();
    for (elem, &count) in elements.iter().zip(counts.iter()) {
        for _ in 0..count {
            let v = floats(next("coordinate")?, 3).context("invalid coordinate")?;
            let pos = if is_direct {
                cell.fractional_to_cartesian(Vector3::new(v[0], v[1], v[2]))
            } else {
                Vector3::new(v[0] * scale, v[1] * scale, v[2] * scale)
            };
            frame.add_atom(Atom::new(elem.as_str(), pos));
        }
    }

    // 速度块（可选），口径同 VASP 文档与 ASE 的 read_vasp：
    // - 坐标之后一行是**模式行**，无论内容都占一行：空行或首字符 C/c/K/k = Cartesian
    //   （Å/fs，不乘缩放因子）；其余 = Direct（晶格矢量/时间步，换算要 POSCAR 里没有
    //   的 POTIM）→ 报错不猜。此前过滤掉空行后把模式行也当速度行，ASE 写的
    //   `Cartesian` 一行让速度少一条、整块静默丢弃
    // - NPT 的 CONTCAR 在这里先有一段 `Lattice velocities and vectors`（首字符 L）：
    //   1 行初始化状态 + 3 行晶格速度 + 3 行晶格矢量，其后一行才是速度的模式行
    // - 模式行之后首行不足两个词 = 没有速度（文件尾，或 MD 的预测校正块）
    let rest: Vec<&str> = lines.collect();
    let k = if rest.first().is_some_and(|l| l.trim_start().starts_with(['L', 'l'])) { 8 } else { 0 };
    let has_vel = rest.get(k + 1).is_some_and(|l| l.split_whitespace().count() >= 2);
    if has_vel {
        let mode = rest[k].trim();
        ensure!(
            mode.is_empty() || mode.starts_with(['C', 'c', 'K', 'k']),
            "velocity block is in direct coordinates (mode line {mode:?}); converting it to \
             Å/fs needs POTIM, which a POSCAR does not carry"
        );
        let vels: Vec<Vector3<f64>> = rest[k + 1..].iter().take(total)
            .map_while(|l| floats(l, 3).ok())
            .map(|v| Vector3::new(v[0], v[1], v[2]))
            .collect();
        if vels.len() == total {
            frame.velocities = Some(vels);
        } else {
            eprintln!(
                "[ferro] warning: velocity block has {} of {total} lines; velocities dropped",
                vels.len()
            );
        }
    }

    let mut traj = Trajectory::from_frame(frame);
    if !comment.is_empty() { traj.metadata.source = Some(comment); }
    Ok(traj)
}


#[cfg(test)]
mod tests {
    use super::*;

    const BCC_FE: &str = "BCC Fe
  1.00000000
     2.87000000   0.00000000   0.00000000
     0.00000000   2.87000000   0.00000000
     0.00000000   0.00000000   2.87000000
   Fe
   2
Direct
  0.0000000  0.0000000  0.0000000
  0.5000000  0.5000000  0.5000000
";

    const SCALE_POSCAR: &str = "Scaled
  2.00000000
     1.435 0.0 0.0
     0.0   1.435 0.0
     0.0   0.0   1.435
Fe
1
Direct
  0.0 0.0 0.0
";

    use crate::testutil::write_tmp as tmp;

    #[test]
    fn test_coordinate_mode_follows_vasp_rule() {
        // 第 2 个原子在分数坐标 (0.5,0.5,0.5) = 笛卡尔 (1.435,…)；Cartesian 读则是 (0.5,…)
        let x1 = |mode: &str| {
            let text = BCC_FE.replace("Direct\n", &format!("{mode}\n"));
            parse_poscar(&text).unwrap().frames[0].atoms[1].position.x
        };
        for mode in ["Direct", "direct", "Fractional", ""] {
            assert!((x1(mode) - 1.435).abs() < 1e-9, "{mode:?} 应按 Direct 读");
        }
        for mode in ["Cartesian", "cart", "K", "k-space"] {
            assert!((x1(mode) - 0.5).abs() < 1e-9, "{mode:?} 应按 Cartesian 读");
        }
    }

    #[test]
    fn test_bcc_fe() {
        let traj = read_poscar(&tmp("bcc.poscar", BCC_FE)).unwrap();
        let frame = traj.first().unwrap();
        assert_eq!(frame.n_atoms(), 2);
        assert_eq!(frame.atom(0).element, "Fe");
        let [a, ..] = frame.cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-6);
    }

    /// ASE 3.29.0 `write(..., format="vasp", direct=True)` 对 O Si O、速度 0.1/0.2/0.3 Å/fs
    /// 的原样输出（元素行不重排、速度块前一行 `Cartesian`），后接各测试自己的速度块
    const ASE_HEAD: &str = "O  Si O
 1.0000000000000000
    10.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000   10.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   10.0000000000000000
 O   Si  O
   1   1   1
Direct
  0.1000000000000000  0.0000000000000000  0.0000000000000000
  0.2000000000000000  0.0000000000000000  0.0000000000000000
  0.3000000000000000  0.0000000000000000  0.0000000000000000
";
    const VELS: &str = "  0.1 0.0 0.0\n  0.2 0.0 0.0\n  0.3 0.0 0.0\n";

    fn vel_x(text: &str) -> Option<Vec<f64>> {
        let traj = parse_poscar(text).unwrap();
        traj.first().unwrap().velocities.as_ref().map(|v| v.iter().map(|v| v.x).collect())
    }

    #[test]
    fn test_velocity_mode_line() {
        // ASE 写 `Cartesian`、VASP 的 CONTCAR 与 ferro 写空行、`K` 也是 Cartesian
        for mode in ["Cartesian", "", "  ", "K"] {
            let text = format!("{ASE_HEAD}{mode}\n{VELS}");
            assert_eq!(vel_x(&text), Some(vec![0.1, 0.2, 0.3]), "模式行 {mode:?} 下速度应完整读出");
        }
        // NPT 的 CONTCAR：晶格速度块在前，其后一行空行才是速度的模式行
        let lattice = "Lattice velocities and vectors\n  1\n  0 0 0\n  0 0 0\n  0 0 0\n  \
                       10 0 0\n  0 10 0\n  0 0 10\n\n";
        assert_eq!(vel_x(&format!("{ASE_HEAD}{lattice}{VELS}")), Some(vec![0.1, 0.2, 0.3]),
                   "晶格速度块之后的原子速度应读出");
    }

    #[test]
    fn test_no_velocity_block() {
        // 文件尾、只有空行、MD 预测校正块（单个数开头）都不是速度
        for tail in ["", "\n", "\n\n", "\n 1\n"] {
            assert_eq!(vel_x(&format!("{ASE_HEAD}{tail}")), None, "尾部 {tail:?} 不应读出速度");
        }
    }

    #[test]
    fn test_direct_velocities_are_an_error() {
        // Direct 速度的单位是 晶格矢量/时间步，换算要 POTIM —— 不猜
        let err = parse_poscar(&format!("{ASE_HEAD}Direct\n{VELS}")).expect_err("Direct 速度应报错");
        assert!(format!("{err:#}").contains("POTIM"), "{err:#}");
    }

    #[test]
    fn test_short_velocity_block_is_dropped() {
        // 少一行：告警并丢弃，而不是把两条速度错配给三个原子
        assert_eq!(vel_x(&format!("{ASE_HEAD}Cartesian\n  0.1 0 0\n  0.2 0 0\n")), None);
    }

    #[test]
    fn test_scaling_factor() {
        let traj = read_poscar(&tmp("scale.poscar", SCALE_POSCAR)).unwrap();
        let [a, ..] = traj.first().unwrap().cell.as_ref().unwrap().lengths();
        assert!((a - 2.87).abs() < 1e-6);
    }
}
