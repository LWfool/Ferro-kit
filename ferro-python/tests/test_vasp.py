"""VASP 生成器：模板 + 结构 → INCAR / POSCAR / POTCAR；模板导出再读回。"""

from pathlib import Path

import pytest

import ferro
from ferro import cli, config
from ferro.jobs.base import JobError, JobOptions
from ferro.jobs.vasp import Incar, VaspGenerator, fine_grid

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


def test_vary_cell(cfg, traj, tmp_path):
    with pytest.raises(JobError, match="--vary-cell needs ENCUT"):
        VaspGenerator(cfg).generate(traj, 0, JobOptions("opt", vary_cell=True))
    tpl = tmp_path / "INCAR"
    tpl.write_text(VaspGenerator.builtin_template("opt") + "ENCUT = 700\nMETAGGA = R2SCAN\n")
    s = VaspGenerator(cfg).generate(traj, 0, JobOptions("opt", vary_cell=True), tpl)
    incar = Incar.from_str(s.files["INCAR"])
    assert (incar["ISIF"], incar["ENCUT"]) == ("3", "700")
    note = next(n for n in s.notes if n.startswith("fine grid"))
    coarse, fine = fine_grid(traj.cell(0), 700, 2)
    assert str(fine) in note and incar["NGXF"] == str(fine[0]), "细网格按模板写的 ENCUT 算"


def test_template_values_are_kept(cfg, traj, tmp_path):
    tpl = tmp_path / "INCAR"
    tpl.write_text("KSPACING = 0.25  ! mine\nNCORE = 8\nTEBEG = 900\nIBRION = 0\nNSW = 5\n")
    s = VaspGenerator(cfg).generate(traj, 0, JobOptions("aimd"), tpl)
    assert s.files["INCAR"].startswith("KSPACING = 0.25  ! mine\nNCORE = 8\nTEBEG = 900"), "未给的量沿用模板"
    assert not any("parallel settings" in n for n in s.notes)


def test_metagga_switches_algo(cfg, traj, tmp_path):
    tpl = tmp_path / "INCAR"
    tpl.write_text(VaspGenerator.builtin_template("sp").replace("GGA = PE", "METAGGA = R2SCAN"))
    s = VaspGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl)
    incar = Incar.from_str(s.files["INCAR"])
    assert (incar["ALGO"], incar["ISEARCH"]) == ("All", "1"), "meta-GGA 用 ALGO = All + ISEARCH = 1"
    assert any("METAGGA = R2SCAN" in n for n in s.notes)
    plain = Incar.from_str(VaspGenerator(cfg).generate(traj, 0, JobOptions("sp")).files["INCAR"])
    assert plain["ALGO"] == "Normal" and "ISEARCH" not in plain, "GGA 不动 ALGO"
    assert "NGXF" not in plain, "GGA 不写细网格"


def test_fine_grid_rule():
    # 与 make_vasp_label.py 同一规则：ENCUT 800、10 Å 立方 → 4·G_cut·L/2π = 92.25 → 96，2 倍 192
    assert fine_grid([[10, 0, 0], [0, 10, 0], [0, 0, 10]], 800, 2) == ([96, 96, 96], [192, 192, 192])
    assert fine_grid([[10, 0, 0], [0, 10, 0], [0, 0, 10]], 800, 3)[1] == [288, 288, 288]
    assert fine_grid([[10, 0, 0], [0, 10, 0], [0, 0, 10]], 800, 2.5)[1] == [240, 240, 240], "非整数倍向上取到 good size"
    # PREC = Normal：3·G_cut·L/2π = 69.2 → 70
    assert fine_grid([[10, 0, 0], [0, 10, 0], [0, 0, 10]], 800, 2, "Normal") == ([70, 70, 70], [140, 140, 140])


def test_metagga_writes_fine_grid(cfg, traj, tmp_path):
    tpl = tmp_path / "INCAR"
    r2scan = VaspGenerator.builtin_template("sp").replace("GGA = PE", "METAGGA = R2SCAN")
    tpl.write_text(r2scan + "ENCUT = 800\n")
    cell = ferro.sort_atoms(traj.subset([0])).cell(0)
    for factor, expected in [(None, 2), (3, 3)]:
        s = VaspGenerator(cfg).generate(traj, 0, JobOptions("sp", fine_grid=factor), tpl)
        incar = Incar.from_str(s.files["INCAR"])
        assert [int(incar[k]) for k in ("NGXF", "NGYF", "NGZF")] == fine_grid(cell, 800, expected)[1]
    assert any("Zn present" in n for n in VaspGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl).notes)
    tpl.write_text(r2scan + "NGXF = 100\nNGYF = 100\nNGZF = 100\n")
    kept = Incar.from_str(VaspGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl).files["INCAR"])
    assert kept["NGXF"] == "100", "模板写了且未给 --fine-grid 就沿用"
    tpl.write_text(r2scan.replace("PREC = Accurate", "PREC = Normal") + "ENCUT = 800\n")
    normal = Incar.from_str(VaspGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl).files["INCAR"])
    assert [int(normal[k]) for k in ("NGXF", "NGYF", "NGZF")] == fine_grid(cell, 800, 2, "Normal")[1], \
        "PREC 不是 Accurate 也照算，粗网格换 3/2·G_cut"
    s = VaspGenerator(cfg).generate(traj, 0, JobOptions("sp", fine_grid=2))
    assert "NGXF" in Incar.from_str(s.files["INCAR"]), "GGA 显式给了 --fine-grid 也照写"
    assert any("rarely needs" in n for n in s.notes), "并提示多半用不着"


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
