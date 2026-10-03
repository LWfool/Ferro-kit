"""CIF 空间群展开的对拍：ferro vs spglib（530 条 Hall 设置）与 ASE（独立的 LGPL 数据表）。

每条设置写一个只含空间群符号、不列对称操作的 CIF（一般位置单原子、按晶系给相容晶胞），
经 `ferro convert` 转 extxyz，比较分数坐标集合：

    (a) Hall 符号        → 应与 spglib 的展开一致
    (b) 本设置的 H-M 符号 → 一致，或按设计报歧义（67/68 号同符号差原点平移的设置）
    (c) ASE 的 (IT 号, setting) → 与 ASE `equivalent_sites` 一致

运行（先 `cargo build`，需要 spglib 与 ASE，deepmd 环境里都有）：

    ~/.miniforge3/envs/deepmd/bin/python scripts/crosscheck_spacegroups.py [ferro 路径]

2026-10-03 的结果：(a) 530/530；(b) 518 一致 + 12 报歧义；(c) 273/274，唯一差异是
67 号 `C m m e`（ASE 静默取标准设置，ferro 报歧义）。出现这些以外的差异时退出码 1。
"""

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import spglib
from ase.io import read as ase_read
from ase.spacegroup import Spacegroup

ROOT = Path(__file__).resolve().parent.parent
FERRO = sys.argv[1] if len(sys.argv) > 1 else str(ROOT / "target/debug/ferro")
POS = np.array([0.1234, 0.2345, 0.3456])
RHOMBOHEDRAL = {146, 148, 155, 160, 161, 166, 167}

# 按设计报歧义的情形（见 dev/issues.md「CIF 空间群符号展开」）
EXPECTED_HM_AMBIGUOUS = {316, 317, 318, 319, 320, 321, 323, 325, 327, 329, 331, 333}
EXPECTED_ASE_AMBIGUOUS = {(67, 1)}


def cell_for(number, axis="b", rhombo=False):
    """与该晶系相容、但不带多余对称性的晶胞参数。axis 为单斜唯一轴。"""
    if number <= 2:
        return (5.1, 6.2, 7.3, 81.0, 95.0, 103.0)
    if number <= 15:
        ang = {"a": (100.0, 90.0, 90.0), "b": (90.0, 100.0, 90.0), "c": (90.0, 90.0, 100.0)}[axis]
        return (5.1, 6.2, 7.3, *ang)
    if number <= 74:
        return (5.1, 6.2, 7.3, 90.0, 90.0, 90.0)
    if number <= 142:
        return (5.1, 5.1, 7.3, 90.0, 90.0, 90.0)
    if number <= 167 and rhombo:
        return (5.5, 5.5, 5.5, 70.0, 70.0, 70.0)
    if number <= 194:
        return (5.1, 5.1, 7.3, 90.0, 90.0, 120.0)
    return (6.0, 6.0, 6.0, 90.0, 90.0, 90.0)


def monoclinic_axis(rotations):
    """单斜群的非恒等旋转是 ±1 对角阵，与另两个对角元不同的那一维即唯一轴。"""
    for r in rotations:
        d = np.diag(r)
        if not np.allclose(r, np.eye(3)) and np.allclose(r, np.diag(d)):
            for i in range(3):
                if d[i] != d[(i + 1) % 3] and d[i] != d[(i + 2) % 3]:
                    return "abc"[i]
    raise ValueError("not a monoclinic rotation set")


def write_cif(path, cell, tag, value):
    a, b, c, al, be, ga = cell
    path.write_text(
        f"data_x\n_cell_length_a {a}\n_cell_length_b {b}\n_cell_length_c {c}\n"
        f"_cell_angle_alpha {al}\n_cell_angle_beta {be}\n_cell_angle_gamma {ga}\n"
        f"{tag} '{value}'\nloop_\n_atom_site_label\n_atom_site_type_symbol\n"
        f"_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
        f"Si1 Si {POS[0]} {POS[1]} {POS[2]}\n")


def ferro_frac(cif):
    """返回 (分数坐标, 错误信息)，二者恰有一个为 None。"""
    out = cif.with_suffix(".extxyz")
    r = subprocess.run([FERRO, "convert", "-i", str(cif), "-o", str(out)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return None, (r.stderr.strip().splitlines() or ["?"])[-1]
    return np.mod(ase_read(out).get_scaled_positions(wrap=False), 1.0), None


def same_set(a, b, tol=1e-4):
    if len(a) != len(b):
        return False
    for p in a:
        d = b - p
        d -= np.round(d)
        if not np.any(np.all(np.abs(d) < tol, axis=1)):
            return False
    return True


def orbit(rots, trans):
    pts = []
    for r, t in zip(rots, trans):
        p = np.mod(r @ POS + t, 1.0)
        if not any(np.all(np.abs((p - q) - np.round(p - q)) < 1e-4) for q in pts):
            pts.append(p)
    return np.array(pts)


def check(tmp):
    counts = {"hall ok": 0, "hm ok": 0, "hm ambiguous": 0, "ase ok": 0, "ase ambiguous": 0}
    unexpected = []

    for h in range(1, 531):
        t = spglib.get_spacegroup_type(h)
        sym = spglib.get_symmetry_from_database(h)
        want = orbit(sym["rotations"], sym["translations"])
        axis = monoclinic_axis(sym["rotations"]) if 3 <= t.number <= 15 else "b"
        cell = cell_for(t.number, axis, rhombo=(t.choice == "R"))

        p = tmp / f"hall_{h}.cif"
        write_cif(p, cell, "_space_group_name_Hall", t.hall_symbol)
        got, err = ferro_frac(p)
        if got is not None and same_set(got, want):
            counts["hall ok"] += 1
        else:
            unexpected.append(("hall", h, t.hall_symbol, err or f"{len(got)} vs {len(want)} atoms"))

        hm = t.international.split("=")[-1].strip()
        if h == 331:  # spglib 的笔误，见 gen_spacegroups.py
            hm = "B b e b"
        if t.choice[:1] in ("1", "2"):
            hm += f" :{t.choice[0]}"
        p = tmp / f"hm_{h}.cif"
        write_cif(p, cell, "_symmetry_space_group_name_H-M", hm)
        got, err = ferro_frac(p)
        if got is not None and same_set(got, want):
            counts["hm ok"] += 1
        elif err and "ambiguous" in err and h in EXPECTED_HM_AMBIGUOUS:
            counts["hm ambiguous"] += 1
        else:
            unexpected.append(("hm", h, hm, err or f"{len(got)} vs {len(want)} atoms"))

    for n in range(1, 231):
        has2 = True
        try:
            Spacegroup(n, 2)
        except Exception:
            has2 = False
        for setting in (1, 2) if has2 else (1,):
            sg = Spacegroup(n, setting)
            rhombo = n in RHOMBOHEDRAL and setting == 2
            axis = monoclinic_axis(sg.rotations) if 3 <= n <= 15 else "b"
            symbol = sg.symbol
            # ASE 的 setting 2 在单斜群是另一种晶胞选择、在 R 群是菱方胞，都不是原点选择
            if has2 and n > 15 and n not in RHOMBOHEDRAL:
                symbol += f" :{setting}"
            want, _ = sg.equivalent_sites([POS], onduplicates="keep")
            p = tmp / f"ase_{n}_{setting}.cif"
            write_cif(p, cell_for(n, axis, rhombo), "_symmetry_space_group_name_H-M", symbol)
            got, err = ferro_frac(p)
            if got is not None and same_set(got, np.mod(want, 1.0)):
                counts["ase ok"] += 1
            elif err and "ambiguous" in err and (n, setting) in EXPECTED_ASE_AMBIGUOUS:
                counts["ase ambiguous"] += 1
            else:
                unexpected.append(("ase", n, setting, symbol, err or f"{len(got)} vs {len(want)} atoms"))

    return counts, unexpected


def main():
    with tempfile.TemporaryDirectory() as d:
        counts, unexpected = check(Path(d))
    print(counts)
    for u in unexpected:
        print("UNEXPECTED", *u)
    sys.exit(1 if unexpected else 0)


if __name__ == "__main__":
    main()
