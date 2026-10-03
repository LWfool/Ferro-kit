//! 空间群符号 → 对称操作。CIF 只写符号（Hall / H-M / IT 号）不列操作时由 `cif.rs` 调用。
//!
//! 数据是 spglib 的 530 条 Hall 设置（`spacegroups.dat`，由 `scripts/gen_spacegroups.py`
//! 生成），经 `include_str!` 编进二进制，首次使用时解析一次。
//!
//! 同一个符号可能对应多个设置（原点选择 1/2、菱方群 H/R 胞、单斜唯一轴）。判定顺序：
//! 符号匹配 → `:1`/`:2`/`:H`/`:R` 后缀 → 度规筛选（旋转须保持晶胞度规）→ 操作集合
//! 相同者合并。仍剩多个就报错 —— 选错原点时原子数对、位置全错，比少原子更难察觉。

use std::sync::OnceLock;

use anyhow::{bail, Result};
use nalgebra::{Matrix3, Vector3};

const DATA: &str = include_str!("spacegroups.dat");

/// 一个操作：旋转矩阵按行 9 个整数 + 平移 3 个（单位 1/24）。整数存放，合并时可精确比较
type IntOp = [i32; 12];

struct Setting {
    hall: u32,
    number: u32,
    choice: &'static str,
    hall_symbol: &'static str,
    hm_setting: &'static str,
    hm_full: &'static str,
    ops: Vec<IntOp>,
}

fn table() -> &'static [Setting] {
    static TABLE: OnceLock<Vec<Setting>> = OnceLock::new();
    TABLE.get_or_init(|| {
        // 数据文件随源码编译进来，格式错是构建产物的错，panic 即可（测试会先抓到）
        let mut out: Vec<Setting> = Vec::with_capacity(530);
        for line in DATA.lines().filter(|l| !l.starts_with('#') && !l.is_empty()) {
            if line.contains('\t') {
                let f: Vec<&'static str> = line.split('\t').collect();
                assert_eq!(f.len(), 8, "spacegroups.dat 表头列数不对: {line}");
                out.push(Setting {
                    hall: f[0].parse().unwrap(),
                    number: f[1].parse().unwrap(),
                    choice: if f[2] == "-" { "" } else { f[2] },
                    hall_symbol: f[3],
                    hm_setting: f[4],
                    hm_full: f[6],
                    ops: Vec::with_capacity(f[7].parse().unwrap()),
                });
            } else {
                let v: Vec<i32> = line.split_whitespace().map(|s| s.parse().unwrap()).collect();
                let op: IntOp = v.try_into().expect("spacegroups.dat 操作行应为 12 个整数");
                out.last_mut().expect("操作行出现在第一个表头之前").ops.push(op);
            }
        }
        out
    })
}

/// CIF 里读到的空间群信息，三项都可缺。
pub(super) struct Query<'a> {
    pub hall: Option<&'a str>,
    pub hm: Option<&'a str>,
    pub number: Option<u32>,
}

/// 查出唯一的设置并返回其全部对称操作（分数坐标下的 旋转, 平移）。
/// `metric` 是晶胞度规 $G = M M^T$（`M` 的行为晶格矢量）。
pub(super) fn lookup(q: &Query, metric: &Matrix3<f64>) -> Result<Vec<(Matrix3<f64>, Vector3<f64>)>> {
    let table = table();
    let mut cands: Vec<&Setting> = if let Some(hall) = q.hall {
        let key = norm_hall(hall);
        let found: Vec<&Setting> = table.iter().filter(|s| norm_hall(s.hall_symbol) == key).collect();
        if found.is_empty() {
            bail!("unknown Hall symbol '{hall}'");
        }
        found
    } else if let Some(hm) = q.hm {
        match_hm(table, hm)?
    } else if let Some(n) = q.number {
        let found: Vec<&Setting> = table.iter().filter(|s| s.number == n).collect();
        if found.is_empty() {
            bail!("space group number {n} is out of range 1-230");
        }
        found
    } else {
        bail!("no space-group symbol or number given");
    };

    if let Some(n) = q.number {
        cands.retain(|s| s.number == n);
        if cands.is_empty() {
            bail!("space-group symbol '{}' contradicts space-group number {n}",
                  q.hall.or(q.hm).unwrap_or_default());
        }
    }

    // 度规筛选：单斜唯一轴、菱方群 H/R 胞由晶胞几何决定。容差取得宽（1e-2），宁可留下
    // 多个候选走报错，也不因 CIF 晶胞参数的舍入把正确设置筛掉
    let scale = metric.diagonal().max();
    let n_before = cands.len();
    cands.retain(|s| s.ops.iter().all(|op| {
        let w = Matrix3::from_fn(|i, j| op[3 * i + j] as f64);
        (w.transpose() * metric * w - metric).amax() <= 1e-2 * scale
    }));
    if cands.is_empty() {
        bail!("the cell parameters are incompatible with all {n_before} setting(s) of this \
               space group; check the cell or give the symmetry operations explicitly");
    }

    // 操作集合完全相同的设置（如 Cmme 的两个 Hall 号）不算歧义
    let sorted_ops = |s: &Setting| { let mut v = s.ops.clone(); v.sort(); v };
    let first = sorted_ops(cands[0]);
    if cands.iter().skip(1).any(|s| sorted_ops(s) != first) {
        let list: Vec<String> = cands.iter()
            .map(|s| format!("Hall #{} '{}' (setting '{}', choice {})", s.hall,
                             s.hall_symbol, s.hm_setting, if s.choice.is_empty() { "-" } else { s.choice }))
            .collect();
        bail!("space-group symbol is ambiguous: matches {}. Append the setting to the H-M symbol \
               (e.g. 'F d -3 m :2'), give the Hall symbol, or list the symmetry operations \
               explicitly", list.join(", "));
    }

    Ok(cands[0].ops.iter().map(|op| {
        let rot = Matrix3::from_fn(|i, j| op[3 * i + j] as f64);
        let trans = Vector3::new(op[9] as f64, op[10] as f64, op[11] as f64) / 24.0;
        (rot, trans)
    }).collect())
}

/// H-M 符号匹配。先比本设置的符号与全符号（精确），都不中再比去掉占位 `1` 的简写
/// （`P 21/c` = `P 1 21/c 1`）。简写不能先比：`P 3 1 2` 与 `P 3 2 1` 去 `1` 后撞名。
fn match_hm<'t>(table: &'t [Setting], hm: &str) -> Result<Vec<&'t Setting>> {
    let (base, suffix) = match hm.split_once(':') {
        Some((b, s)) => (b, Some(s.trim())),
        None => (hm, None),
    };
    let key = norm_hm(base);
    // 旧符号（ITA 引入 e 滑移面前）→ (IT 号, 轴序)。不能映射到新符号：67、68 号的
    // 新符号不分轴序（Cmma 与 Cmmb 都成 Cmme，两者差一个原点平移），旧符号反倒是唯一的
    let old = match key.as_str() {
        "abm2" => Some((39, "")),
        "aba2" => Some((41, "")),
        "cmca" => Some((64, "")),
        "cmma" => Some((67, "")),
        "cmmb" => Some((67, "ba-c")),
        "abmm" => Some((67, "cab")),
        "acmm" => Some((67, "-cba")),
        "bmcm" => Some((67, "bca")),
        "bmam" => Some((67, "a-cb")),
        "ccca" => Some((68, "")),
        "cccb" => Some((68, "ba-c")),
        "abaa" => Some((68, "cab")),
        "acaa" => Some((68, "-cba")),
        "bbcb" => Some((68, "bca")),
        "bbab" => Some((68, "a-cb")),
        _ => None,
    };

    let mut found: Vec<&Setting> = match old {
        // 68 号的 choice 是「原点 + 轴序」（如 `2ba-c`），去掉原点数字再比
        Some((n, axes)) => table.iter()
            .filter(|s| s.number == n && s.choice.trim_start_matches(['1', '2']) == axes)
            .collect(),
        None => table.iter()
            .filter(|s| norm_hm(s.hm_setting) == key || norm_hm(s.hm_full) == key)
            .collect(),
    };
    if found.is_empty() {
        found = table.iter().filter(|s| norm_hm(&drop_ones(s.hm_setting)) == key).collect();
    }
    if found.is_empty() {
        bail!("unknown Hermann-Mauguin symbol '{hm}'");
    }

    if let Some(sfx) = suffix {
        // `:1` 也要匹配正交群的 `1cab` 这类组合 choice（轴序已由符号本身确定）
        let sfx = sfx.to_lowercase();
        found.retain(|s| {
            let c = s.choice.to_lowercase();
            c == sfx || (matches!(sfx.as_str(), "1" | "2") && c.starts_with(&sfx))
        });
        if found.is_empty() {
            bail!("setting ':{}' does not exist for '{}'", suffix.unwrap_or_default(), base.trim());
        }
    }
    Ok(found)
}

/// 去空白、去下划线、转小写。首字母恒为晶格字母，转小写不会与滑移面撞名
fn norm_hm(s: &str) -> String {
    s.chars().filter(|c| !c.is_whitespace() && *c != '_').flat_map(char::to_lowercase).collect()
}

/// `P 1 2_1/c 1` → `P 2_1/c`：去掉作占位的独立 `1`（`P 1`、`P -1` 本身保留）
fn drop_ones(s: &str) -> String {
    let toks: Vec<&str> = s.split_whitespace().collect();
    if toks.len() <= 2 { return s.to_string(); }
    let mut out = vec![toks[0]];
    out.extend(toks[1..].iter().filter(|t| **t != "1"));
    out.join(" ")
}

/// Hall 符号：空白折叠、不分大小写。大写只出现在首位的晶格字母，不会撞名
fn norm_hall(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ").to_lowercase()
}

#[cfg(test)]
mod tests {
    use super::*;
    use ferro_core::Cell;

    fn metric(a: f64, b: f64, c: f64, al: f64, be: f64, ga: f64) -> Matrix3<f64> {
        let m = Cell::from_lengths_angles(a, b, c, al, be, ga).unwrap().matrix;
        m * m.transpose()
    }

    fn q<'a>(hall: Option<&'a str>, hm: Option<&'a str>, number: Option<u32>) -> Query<'a> {
        Query { hall, hm, number }
    }

    #[test]
    fn test_table_complete() {
        let t = table();
        assert_eq!(t.len(), 530, "spglib 的 Hall 设置应为 530 条");
        assert_eq!(t.iter().map(|s| s.ops.len()).sum::<usize>(), 7388);
        assert!((1..=530).all(|h| t[h as usize - 1].hall == h), "Hall 号应连续");
        assert!((1..=230).all(|n| t.iter().any(|s| s.number == n)), "230 个空间群都应在表里");
    }

    #[test]
    fn test_fm3m_unique() {
        let cubic = metric(5.64, 5.64, 5.64, 90.0, 90.0, 90.0);
        let ops = lookup(&q(None, Some("F m -3 m"), None), &cubic).unwrap();
        assert_eq!(ops.len(), 192);
        let ops = lookup(&q(None, None, Some(225)), &cubic).unwrap();
        assert_eq!(ops.len(), 192, "225 只有一个设置，只给 IT 号也应唯一");
    }

    #[test]
    fn test_origin_choice_needs_suffix() {
        let cubic = metric(8.0, 8.0, 8.0, 90.0, 90.0, 90.0);
        let err = lookup(&q(None, Some("F d -3 m"), None), &cubic).unwrap_err().to_string();
        assert!(err.contains("ambiguous"), "Fd-3m 不带后缀应报歧义: {err}");
        let o1 = lookup(&q(None, Some("F d -3 m :1"), None), &cubic).unwrap();
        let o2 = lookup(&q(None, Some("Fd-3m:2"), None), &cubic).unwrap();
        assert_eq!((o1.len(), o2.len()), (192, 192));
        assert!(o2.iter().any(|(r, _)| *r == -Matrix3::identity()), "原点 2 在反演中心上");
        let h = lookup(&q(Some("-F 4vw 2vw 3"), None, None), &cubic).unwrap();
        assert_eq!(h.len(), 192, "Hall 符号应直接唯一");
    }

    #[test]
    fn test_rhombohedral_by_metric() {
        let hex = metric(5.0, 5.0, 14.0, 90.0, 90.0, 120.0);
        let rho = metric(5.5, 5.5, 5.5, 60.0, 60.0, 60.0);
        assert_eq!(lookup(&q(None, Some("R -3 m"), None), &hex).unwrap().len(), 36, "六方胞含 R 心");
        assert_eq!(lookup(&q(None, Some("R -3 m"), None), &rho).unwrap().len(), 12, "菱方原胞");
    }

    #[test]
    fn test_monoclinic_short_symbol() {
        let beta = metric(5.0, 6.0, 7.0, 90.0, 100.0, 90.0);
        let ops = lookup(&q(None, Some("P 21/c"), Some(14)), &beta).unwrap();
        assert_eq!(ops.len(), 4);
        // b 唯一轴：二重螺旋轴 (-x, y+1/2, -z+1/2)
        let screw = Matrix3::new(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0);
        assert!(ops.iter().any(|(r, t)| *r == screw && *t == Vector3::new(0.0, 0.5, 0.5)));
        // 只给 IT 号：P 1 21/c 1、P 1 21/n 1、P 1 21/a 1 度规相同、操作不同
        let err = lookup(&q(None, None, Some(14)), &beta).unwrap_err().to_string();
        assert!(err.contains("ambiguous"), "{err}");
    }

    #[test]
    fn test_conflicts_and_unknowns() {
        let cubic = metric(5.0, 5.0, 5.0, 90.0, 90.0, 90.0);
        assert!(lookup(&q(None, Some("F m -3 m"), Some(221)), &cubic).is_err(), "符号与编号矛盾");
        assert!(lookup(&q(None, Some("X 9"), None), &cubic).is_err());
        assert!(lookup(&q(None, None, Some(231)), &cubic).is_err());
        assert!(lookup(&q(None, Some("P 6/m m m"), None), &cubic).is_err(), "立方胞不容六方群");
    }

    #[test]
    fn test_old_symbols_pick_axis_order() {
        let ortho = metric(5.0, 6.0, 7.0, 90.0, 90.0, 90.0);
        // 新符号 Cmme 同时是 abc 与 ba-c 两种设置，差一个原点平移，必须报歧义
        assert!(lookup(&q(None, Some("C m m e"), None), &ortho).is_err());
        let cmma = lookup(&q(None, Some("C m m a"), None), &ortho).unwrap();
        let cmmb = lookup(&q(None, Some("C m m b"), None), &ortho).unwrap();
        assert_eq!((cmma.len(), cmmb.len()), (16, 16));
        assert_ne!(cmma, cmmb, "Cmma 与 Cmmb 是不同设置");
        assert!(lookup(&q(None, Some("C c c a"), None), &ortho).is_err(), "68 号还要原点选择");
        assert_eq!(lookup(&q(None, Some("Ccca:2"), None), &ortho).unwrap().len(), 16);
        // spglib 把 Hall 331 写成 'B b c b'，生成脚本已改回 'B b e b'：:2 应有两个轴序候选
        let err = lookup(&q(None, Some("B b e b :2"), None), &ortho).unwrap_err().to_string();
        assert!(err.contains("2bca") && err.contains("2a-cb"), "{err}");
    }
}
