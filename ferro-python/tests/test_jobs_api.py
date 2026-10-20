"""ferro.jobs 依赖的 Rust 侧接口：读入映射、选帧、排序、结构块字符串、namelist。

只用仓库 tests/ 下的 fixture（examples/ 不进 git，CI 上没有）。
"""

from pathlib import Path

import pytest

import ferro

FIXTURES = Path(__file__).resolve().parents[2] / "tests"
NVT = FIXTURES / "70Z30P00A_NVT_5.lammpstrj"


@pytest.fixture(scope="module")
def nvt():
    return ferro.read(str(NVT))


def dump_without_element_column(tmp_path):
    """把 NVT fixture 的 element 列删掉，记下原文的 类型号 → 元素。"""
    out, mapping = [], {}
    for line in NVT.read_text().splitlines():
        parts = line.split()
        if line.startswith("ITEM: ATOMS"):
            out.append("ITEM: ATOMS id type x y z")
        elif not line.startswith("ITEM") and len(parts) == 6:
            mapping[int(parts[1])] = parts[2]
            out.append(" ".join(parts[:2] + parts[3:]))
        else:
            out.append(line)
    path = tmp_path / "noelem.lammpstrj"
    path.write_text("\n".join(out) + "\n")
    return path, mapping


def test_type_map_matches_element_column(tmp_path, nvt):
    path, mapping = dump_without_element_column(tmp_path)
    assert set(ferro.read(str(path)).symbols(0)) == {f"X{t}" for t in mapping}, "无映射时是占位符"
    mapped = ferro.read(str(path), type_map=mapping)
    for i in range(len(nvt)):
        assert mapped.symbols(i) == nvt.symbols(i), f"第 {i} 帧映射结果应与 element 列相同"
    with pytest.raises(RuntimeError, match="missing from the type map"):
        ferro.read(str(path), type_map={1: "P"})
    with pytest.raises(RuntimeError, match="LAMMPS dump input only"):
        ferro.read(str(FIXTURES / "cp2k_sp_v61.out"), type_map={1: "P"})


def test_frame_indices_follow_convert(nvt):
    assert nvt.frame_indices() == [0, 1, 2, 3, 4], "默认全部帧"
    assert nvt.frame_indices(number=3) == [0, 2, 4], "number 均匀取、两端都取"
    assert nvt.frame_indices(start=1, stride=2) == [1, 3]
    assert nvt.frame_indices(start=1, end=3) == [1, 2, 3], "end 是闭区间"
    for bad in [dict(stride=2, number=2), dict(stride=0), dict(number=0), dict(start=3, end=1),
                dict(start=9)]:
        with pytest.raises(RuntimeError):
            nvt.frame_indices(**bad)


def test_subset_pbc_and_sort(nvt):
    sub = nvt.subset([4, 0])
    assert len(sub) == 2 and sub.positions(0) == nvt.positions(4), "按给定顺序取帧"
    with pytest.raises(IndexError):
        nvt.subset([5])
    assert nvt.pbc(0) == (True, True, True)
    s = ferro.sort_atoms(nvt.subset([0])).symbols(0)
    blocks = [e for i, e in enumerate(s) if i == 0 or s[i - 1] != e]
    assert blocks == ["O", "P", "Zn"], f"应按原子序数成块：{blocks}"


def test_poscar_string_equals_write_poscar(tmp_path, nvt):
    one = nvt.subset([0])
    ferro.write(one, str(tmp_path / "POSCAR"))
    written = (tmp_path / "POSCAR").read_text().splitlines()
    text = ferro.poscar_string(one, comment="test").splitlines()
    assert text[0] == "test"
    assert text[1:] == written[1:], "除首行注释外应与 write_poscar 逐行相同"


def test_qe_cards_and_namelists_round_trip(tmp_path, nvt):
    one = ferro.sort_atoms(nvt.subset([0]))
    with pytest.raises(RuntimeError, match="no pseudopotential given"):
        ferro.qe_structure_cards(one, {"O": "O.UPF"})
    named = ferro.qe_structure_cards(one, {e: f"{e}.pbe.UPF" for e in one.elements()})
    assert [n for n, _ in named] == ["ATOMIC_SPECIES", "CELL_PARAMETERS", "ATOMIC_POSITIONS"]
    cards = "\n".join(text for _, text in named)
    namelists = [("control", [("calculation", "'scf'")]),
                 ("system", [("ibrav", "0"), ("nat", str(one.n_atoms())), ("ntyp", "3")])]
    path = tmp_path / "pw.in"
    path.write_text(ferro.write_namelists(namelists) + cards)
    back = ferro.read(str(path))
    assert back.symbols(0) == one.symbols(0)
    assert max(abs(a - b) for p, q in zip(back.positions(0), one.positions(0)) for a, b in zip(p, q)) < 1e-9
    parsed, rest = ferro.parse_qe_input_parts(path.read_text())
    assert parsed == namelists, "写回再拆应得到同一份 namelist"
    assert rest.strip() == cards.strip()


def test_one_line_namelist_keeps_quoted_slash():
    parsed, rest = ferro.parse_qe_input_parts(
        "&CONTROL calculation='md', pseudo_dir='./pp/' /\nK_POINTS gamma\n")
    assert parsed == [("control", [("calculation", "'md'"), ("pseudo_dir", "'./pp/'")])]
    assert rest == "K_POINTS gamma"


def test_cp2k_sections_round_trip(tmp_path, nvt):
    one = nvt.subset([0])
    path = tmp_path / "cp2k.inp"
    path.write_text("&FORCE_EVAL\n&SUBSYS\n" + ferro.cp2k_structure_sections(one)
                    + "&END SUBSYS\n&END FORCE_EVAL\n")
    back = ferro.read(str(path))
    assert back.symbols(0) == one.symbols(0)
    cell = [abs(a - b) for r, s in zip(back.cell(0), one.cell(0)) for a, b in zip(r, s)]
    assert max(cell) < 1e-9
