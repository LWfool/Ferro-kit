"""VASP 生成器：模板 + 结构 → INCAR / POSCAR / POTCAR；模板导出再读回。"""

from pathlib import Path

import pytest

import ferro
from ferro import cli, config
from ferro.jobs.base import JobError, JobOptions
from ferro.jobs.vasp import Incar, VaspGenerator

NVT = Path(__file__).resolve().parents[2] / "tests" / "70Z30P00A_NVT_5.lammpstrj"
# 假 POTCAR：头部只有生成器要读的两项（价电子数取真实 PBE 值：O 6、P 5、Zn 12）
ZVAL = {"O": 6, "P": 5, "Zn": 12}


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    pp = tmp_path / "pp"
    for name, z in [("O", 6), ("P", 5), ("Zn", 12), ("O_h", 6)]:
        (pp / name).mkdir(parents=True)
        (pp / name / "POTCAR").write_text(f"  PAW_PBE {name}\n   ZVAL   =    {z}.000\n   ENMAX  =  {400 + z}.0\nEnd of Dataset\n")
    (tmp_path / "c.toml").write_text(f'[vasp]\npotcar_dir = "{pp}"\n')
    monkeypatch.setenv("FERRO_CONFIG", str(tmp_path / "c.toml"))
    return config.load(tmp_path / "c.toml")


@pytest.fixture(scope="module")
def traj():
    return ferro.read(str(NVT))


def electrons(traj):
    symbols = ferro.sort_atoms(traj.subset([0])).symbols(0)
    return sum(ZVAL[s] for s in symbols)


@pytest.mark.parametrize("task", ["sp", "opt", "aimd"])
def test_builtin_templates_generate(cfg, traj, task):
    s = VaspGenerator(cfg).generate(traj, 0, JobOptions(task))
    assert set(s.files) == {"INCAR", "POSCAR", "POTCAR"}, "不生成 KPOINTS（k 点在 INCAR 的 KSPACING）"
    incar = Incar.from_str(s.files["INCAR"])
    assert incar["KSPACING"] == "0.5" and incar["GGA"] == "PE" and incar["IVDW"] == "12"
    assert "NELECT" not in incar and "NUPDOWN" not in incar, "中性、偶数电子时不写"
    one = ferro.sort_atoms(traj.subset([0]))
    assert s.files["POSCAR"] == ferro.poscar_string(one, 0, " ".join(one.elements()))
    order = [line.split()[1] for line in s.files["POTCAR"].splitlines() if "PAW_PBE" in line]
    assert order == one.elements(), "POTCAR 顺序与 POSCAR 元素行一致"
    assert any("parallel settings" in n for n in s.notes), "告诉用户并行参数要自己加"
    assert any("mesh" in n for n in s.notes)


def test_charge_spin_and_md_options(cfg, traj):
    n = electrons(traj)
    s = VaspGenerator(cfg).generate(traj, 0, JobOptions(
        "aimd", charge=1, temperature=1500, steps=10, timestep=2.0, kspacing=0.3))
    incar = Incar.from_str(s.files["INCAR"])
    assert incar["NELECT"] == str(n - 1)
    assert incar["ISPIN"] == "2" and incar["NUPDOWN"] == "1", "奇数电子默认多重度 2"
    assert (incar["TEBEG"], incar["NSW"], incar["POTIM"], incar["KSPACING"]) == ("1500", "10", "2", "0.3")
    with pytest.raises(JobError, match="multiplicity 1 is impossible"):
        VaspGenerator(cfg).generate(traj, 0, JobOptions("sp", charge=1, multiplicity=1))
    with pytest.raises(JobError, match="NVT only"):
        VaspGenerator(cfg).generate(traj, 0, JobOptions("aimd", ensemble="npt_i"))


def test_vary_cell(cfg, traj):
    incar = Incar.from_str(VaspGenerator(cfg).generate(traj, 0, JobOptions("opt", vary_cell=True)).files["INCAR"])
    assert incar["ISIF"] == "3"
    assert incar["ENCUT"] == "536", "未写 ENCUT 时取 ⌈1.3 × 最大 ENMAX⌉（Zn 412 → 535.6）"


def test_template_values_are_kept(cfg, traj, tmp_path):
    tpl = tmp_path / "INCAR"
    tpl.write_text("KSPACING = 0.25  ! mine\nNCORE = 8\nTEBEG = 900\nIBRION = 0\nNSW = 5\n")
    s = VaspGenerator(cfg).generate(traj, 0, JobOptions("aimd"), tpl)
    assert s.files["INCAR"].startswith("KSPACING = 0.25  ! mine\nNCORE = 8\nTEBEG = 900"), "未给的量沿用模板"
    assert not any("parallel settings" in n for n in s.notes)


def test_structure_dependent_arrays_must_match(cfg, traj, tmp_path):
    tpl = tmp_path / "INCAR"
    tpl.write_text("MAGMOM = 2*1.0\n")
    with pytest.raises(JobError, match="MAGMOM has 2 values, this structure has"):
        VaspGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl)
    n = len(traj.symbols(0))
    tpl.write_text(f"MAGMOM = {n - 1}*0.0 1.0\nLANGEVIN_GAMMA = 10 10 10\n")
    VaspGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl)


def test_potcar_choice(cfg, traj, tmp_path):
    (tmp_path / "c.toml").write_text(f'[vasp]\npotcar_dir = "{tmp_path / "pp"}"\n[vasp.potcar]\nO = "O_h"\n')
    s = VaspGenerator(config.load(tmp_path / "c.toml")).generate(traj, 0, JobOptions("sp"))
    assert "PAW_PBE O_h" in s.files["POTCAR"], "配置里逐元素指定的变体"
    (tmp_path / "pp" / "Zn" / "POTCAR").unlink()
    with pytest.raises(JobError, match="no POTCAR for Zn"):
        VaspGenerator(cfg).generate(traj, 0, JobOptions("sp"))
    with pytest.raises(config.ConfigError, match="potcar_dir is not set"):
        VaspGenerator(config.load(tmp_path / "none.toml")).generate(traj, 0, JobOptions("sp"))


def test_export_then_read_back(cfg, tmp_path, monkeypatch):
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    out = tmp_path / "tpl"
    assert cli.main(["job", "template", "-s", "vasp", "-t", "aimd", "-o", str(out), "--mkdir"]) == 0
    assert (out / "INCAR").read_text() == VaspGenerator.builtin_template("aimd"), "导出的就是内置模板"
    assert cli.main(["job", "template", "-s", "vasp", "-t", "aimd", "-o", str(out)]) == 1, "不覆盖"
    (out / "INCAR").write_text((out / "INCAR").read_text().replace("POTIM = 1.0", "POTIM = 0.5"))
    jobs_dir = tmp_path / "jobs"
    assert cli.main(["job", "-s", "vasp", "-t", "aimd", "-i", str(NVT), "--number", "2",
                     "--template", str(out / "INCAR"), "-o", str(jobs_dir), "--mkdir"]) == 0
    first = jobs_dir / "70Z30P00A_NVT_5_0000"
    assert "POTIM = 0.5" in (first / "INCAR").read_text(), "改过的模板读回"
    assert (first / "POTCAR").is_symlink() and (jobs_dir / "POTCAR").is_file()
