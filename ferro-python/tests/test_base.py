"""生成器基类：选项校验、自旋奇偶、k 网格、模板查找、共有的结构检查、InputSet 落盘。"""

import pytest

import ferro
from ferro import config
from ferro.jobs.base import Generator, InputSet, JobError, JobOptions, kmesh, multiplicity, write_batch

from conftest import NVT


@pytest.mark.parametrize("n, charge, requested, expected", [
    (16, 0, None, 1), (15, 0, None, 2), (16, 1, None, 2), (16, 0, 3, 3), (15, 0, 4, 4),
])
def test_multiplicity(n, charge, requested, expected):
    assert multiplicity(n, charge, requested) == expected


@pytest.mark.parametrize("n, charge, requested, message", [
    (16, 0, 2, "multiplicity 2 is impossible with 16 electrons"),
    (15, 0, 1, "impossible with 15 electrons"),
    (15.5, 0, None, "not a whole number"),
    (2, 3, None, "removes more than"),
])
def test_multiplicity_errors(n, charge, requested, message):
    with pytest.raises(JobError, match=message):
        multiplicity(n, charge, requested)


def test_kmesh():
    assert kmesh([[4, 0, 0], [0, 8, 0], [0, 0, 20]], 0.5) == [4, 2, 1], "⌈2π/(a·kspacing)⌉"
    assert kmesh([[2 * 3.141592653589793, 0, 0], [0, 100, 0], [0, 0, 100]], 0.5) == [2, 1, 1], \
        "恰好整数时不因舍入误差多进一格"
    hexagonal = [[3, 0, 0], [-1.5, 2.598076211353316, 0], [0, 0, 5]]
    assert kmesh(hexagonal, 0.3) == [9, 9, 5], "斜胞按倒格矢长度，不按 a/b/c"


@pytest.mark.parametrize("options, message", [
    (JobOptions("md"), "unknown task 'md'"),
    (JobOptions("sp", temperature=300), "temperature only apply to an aimd task"),
    (JobOptions("aimd", ensemble="nve"), "unknown ensemble"),
    (JobOptions("sp", vary_cell=True), "vary_cell only applies to an opt task"),
    (JobOptions("aimd", timestep=0.0), "timestep must be a positive"),
    (JobOptions("aimd", steps=0), "steps must be at least 1"),
    (JobOptions("sp", kspacing=float("nan")), "kspacing must be a positive"),
])
def test_option_errors(options, message):
    with pytest.raises(JobError, match=message):
        options.check()


class Echo(Generator):
    """测试用：把模板与结构原样放进 InputSet。"""
    program = "echo"
    template_name = "in.txt"

    def build(self, one, options, template_text):
        return InputSet({"in.txt": template_text, "POSCAR": ferro.poscar_string(one)}, shared={"POSCAR"})


@pytest.fixture
def user_templates(tmp_path):
    (tmp_path / "tpl" / "echo" / "sp").mkdir(parents=True)
    (tmp_path / "tpl" / "echo" / "sp" / "in.txt").write_text("user\n")
    (tmp_path / "c.toml").write_text('[templates]\ndir = "tpl"\n')
    return config.load(tmp_path / "c.toml")


def test_template_lookup(tmp_path, user_templates):
    g = Echo(user_templates)
    assert g.template("sp") == "user\n", "用户模板目录里同名的覆盖内置"
    (tmp_path / "mine.txt").write_text("explicit\n")
    assert g.template("sp", tmp_path / "mine.txt") == "explicit\n", "显式给的模板优先"
    with pytest.raises(JobError, match="no built-in echo template"):
        g.template("opt")
    with pytest.raises(JobError, match="does not exist"):
        g.template("sp", tmp_path / "missing.txt")


def test_generate_sorts_and_checks(tmp_path, user_templates):
    traj = ferro.read(str(NVT))
    s = Echo(user_templates).generate(traj, 2, JobOptions("sp"))
    assert s.files["POSCAR"] == ferro.poscar_string(ferro.sort_atoms(traj.subset([2]))), "取指定帧并按 (Z, 符号) 排序"

    ferro.write(traj.tail(1), str(tmp_path / "mol.xyz"), format="xyz")
    with pytest.raises(JobError, match="only periodic systems"):
        Echo(user_templates).generate(ferro.read(str(tmp_path / "mol.xyz")), 0, JobOptions("sp"))

    lines = NVT.read_text().splitlines()
    out = []
    for line in lines:
        parts = line.split()
        if line.startswith("ITEM: ATOMS"):
            out.append("ITEM: ATOMS id type x y z")
        elif not line.startswith("ITEM") and len(parts) == 6:
            out.append(" ".join(parts[:2] + parts[3:]))
        else:
            out.append(line)
    (tmp_path / "noelem.lammpstrj").write_text("\n".join(out) + "\n")
    with pytest.raises(JobError, match="X1, X2.*--type-map"):
        Echo(user_templates).generate(ferro.read(str(tmp_path / "noelem.lammpstrj")), 0, JobOptions("sp"))
    # a CIF label with no leading letter reads as a bare `X`
    (tmp_path / "x.cif").write_text(
        "data_x\n_cell_length_a 5\n_cell_length_b 5\n_cell_length_c 5\n_cell_angle_alpha 90\n"
        "_cell_angle_beta 90\n_cell_angle_gamma 90\n_symmetry_space_group_name_H-M 'P 1'\nloop_\n"
        "_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
        "1 0 0 0\nO1 0.5 0.5 0.5\n")
    with pytest.raises(JobError, match="atoms X have no element"):
        Echo(user_templates).generate(ferro.read(str(tmp_path / "x.cif")), 0, JobOptions("sp"))


def test_input_set_write_with_links(tmp_path):
    s = InputSet({"INCAR": "A\n", "POTCAR": "P\n"}, shared={"POTCAR"})
    (tmp_path / "POTCAR").write_text("P\n")
    (tmp_path / "f0").mkdir()
    s.write(tmp_path / "f0", links={"POTCAR": "../POTCAR"})
    s.write(tmp_path / "f0", links={"POTCAR": "../POTCAR"})   # 重写已有的链接与文件
    assert (tmp_path / "f0" / "POTCAR").is_symlink()
    assert (tmp_path / "f0" / "POTCAR").read_text() == "P\n", "相对链接指向上层的那一份"
    assert (tmp_path / "f0" / "INCAR").read_text() == "A\n"


def test_write_batch_refuses_other_composition(tmp_path):
    def batch(potcar):
        s = InputSet({"INCAR": "A\n", "POTCAR": potcar}, shared={"POTCAR"})
        return [("f0", s), ("f1", s)]
    write_batch(batch("P\n"), tmp_path)
    write_batch(batch("P\n"), tmp_path)   # 同一组成重写不报错
    with pytest.raises(JobError, match="another composition"):
        write_batch(batch("Q\n"), tmp_path)
    assert (tmp_path / "POTCAR").read_text() == "P\n", "拒绝时不得改动已有的共用 POTCAR"


def test_kmesh_rejects_a_singular_cell():
    with pytest.raises(ferro.FerroError, match="singular"):
        kmesh([[1, 0, 0], [2, 0, 0], [0, 0, 1]], 0.5)


def test_pressure_conversion_is_codata_2018():
    assert ferro.EV_ANG3_TO_KBAR == pytest.approx(1602.176634, rel=1e-12)
