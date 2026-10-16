"""生成器基类：选项校验、自旋奇偶、k 网格、模板查找、共有的结构检查、InputSet 落盘。"""

from pathlib import Path

import pytest

import ferro
from ferro import config
from ferro.jobs.base import Generator, InputSet, JobError, JobOptions, kmesh, multiplicity

NVT = Path(__file__).resolve().parents[2] / "tests" / "70Z30P00A_NVT_5.lammpstrj"


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
    with pytest.raises(JobError, match="fine_grid does not apply to echo"):
        Echo(user_templates).generate(traj, 0, JobOptions("sp", fine_grid=2))
    s = Echo(user_templates).generate(traj, 2, JobOptions("sp"))
    assert s.files["POSCAR"] == ferro.poscar_string(ferro.sort_atoms(traj.subset([2]))), "取指定帧并按 (Z, 符号) 排序"

    ferro.write(traj.tail(1), str(tmp_path / "mol.xyz"))
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


def test_input_set_write_with_links(tmp_path):
    s = InputSet({"INCAR": "A\n", "POTCAR": "P\n"}, shared={"POTCAR"})
    (tmp_path / "POTCAR").write_text("P\n")
    (tmp_path / "f0").mkdir()
    s.write(tmp_path / "f0", links={"POTCAR": "../POTCAR"})
    s.write(tmp_path / "f0", links={"POTCAR": "../POTCAR"})   # 重写已有的链接与文件
    assert (tmp_path / "f0" / "POTCAR").is_symlink()
    assert (tmp_path / "f0" / "POTCAR").read_text() == "P\n", "相对链接指向上层的那一份"
    assert (tmp_path / "f0" / "INCAR").read_text() == "A\n"
