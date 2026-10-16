"""CP2K 生成器：&KIND 从数据文件查、GPW 与 GAPW 全电子两条路、k 点与 OT、MD 系综。"""

from pathlib import Path

import pytest

import ferro
from ferro import config
from ferro.jobs.base import JobError, JobOptions
from ferro.jobs.cp2k import Cp2kGenerator, Section, data_entries

NVT = Path(__file__).resolve().parents[2] / "tests" / "70Z30P00A_NVT_5.lammpstrj"

# CP2K 2026.2 data/ 的条目摘录：只留条目头与下一行（生成器只读这两行）
BASIS_MOLOPT = """\
# excerpt
 O  DZVP-MOLOPT-SR-GTH DZVP-MOLOPT-SR-GTH-q6
 1
 P  DZVP-MOLOPT-SR-GTH DZVP-MOLOPT-SR-GTH-q5
 1
 Zn DZVP-MOLOPT-SR-GTH DZVP-MOLOPT-SR-GTH-q12
 1
"""
GTH_POTENTIALS = """\
O GTH-PBE-q6 GTH-PBE
    2    4
P GTH-PBE-q5 GTH-PBE
    2    3
Zn GTH-PBE-q12 GTH-PBE
    2    0   10
"""
BASIS_POB = """\
O  pob-TZVP-rev2
 9
P  pob-TZVP-rev2
 10
Zn pob-TZVP-rev2
 14
Sn pob-TZVP-rev2 pob-TZVP-rev2-q22
 12
"""
ALL_POTENTIALS = """\
O ALLELECTRON ALL
    4    4    0
P ALLELECTRON ALL
    6    9    0
Zn ALLELECTRON ALL
    8   12   10
"""


@pytest.fixture
def cfg(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    for name, text in [("BASIS_MOLOPT", BASIS_MOLOPT), ("GTH_POTENTIALS", GTH_POTENTIALS),
                       ("BASIS_pob", BASIS_POB), ("ALL_POTENTIALS", ALL_POTENTIALS)]:
        (data / name).write_text(text)
    (tmp_path / "c.toml").write_text(f'[cp2k]\ndata_dir = "{data}"\n')
    return config.load(tmp_path / "c.toml")


@pytest.fixture(scope="module")
def traj():
    return ferro.read(str(NVT))


def gapw(task):
    return (Cp2kGenerator.builtin_template(task)
            .replace("METHOD GPW", "METHOD GAPW")
            .replace("BASIS_SET_FILE_NAME BASIS_MOLOPT", "BASIS_SET_FILE_NAME BASIS_pob")
            .replace("POTENTIAL_FILE_NAME GTH_POTENTIALS", "POTENTIAL_FILE_NAME ALL_POTENTIALS"))


def test_gapw_xc_is_refused(cfg, traj, tmp_path):
    tpl = tmp_path / "t.inp"
    tpl.write_text(Cp2kGenerator.builtin_template("sp").replace("METHOD GPW", "METHOD GAPW_XC"))
    with pytest.raises(JobError, match="GAPW_XC is not supported"):
        Cp2kGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl)


def test_data_entries(tmp_path):
    (tmp_path / "P").write_text(GTH_POTENTIALS)
    assert data_entries(tmp_path / "P")["Zn"] == [(["GTH-PBE-q12", "GTH-PBE"], [2, 0, 10])]


@pytest.mark.parametrize("task", ["sp", "opt", "aimd"])
def test_builtin_templates_generate(cfg, traj, task, tmp_path):
    s = Cp2kGenerator(cfg).generate(traj, 0, JobOptions(task))
    root = Section.from_str(s.files["cp2k.inp"])
    kinds = {k.param: (k.get("BASIS_SET"), k.get("POTENTIAL")) for k in root.section("FORCE_EVAL/SUBSYS").sections("KIND")}
    assert kinds == {"O": ("DZVP-MOLOPT-SR-GTH-q6", "GTH-PBE-q6"), "P": ("DZVP-MOLOPT-SR-GTH-q5", "GTH-PBE-q5"),
                     "Zn": ("DZVP-MOLOPT-SR-GTH-q12", "GTH-PBE-q12")}, "基组 -qN 配同 q 的 GTH"
    dft = root.section("FORCE_EVAL/DFT")
    assert (dft.get("CHARGE"), dft.get("MULTIPLICITY"), dft.get("UKS")) == ("0", "1", None)
    assert root.section("FORCE_EVAL/DFT/KPOINTS") is None, "大胞 Γ 点不写 &KPOINTS"
    assert root.section("FORCE_EVAL/DFT/SCF/OT") is not None
    path = tmp_path / "cp2k.inp"
    path.write_text(s.files["cp2k.inp"])
    back = ferro.read(str(path))
    one = ferro.sort_atoms(traj.subset([0]))
    assert back.symbols(0) == one.symbols(0), "写出的 &COORD 可被读回"
    ignore = root.section("FORCE_EVAL/DFT/SCF").get("IGNORE_CONVERGENCE_FAILURE")
    assert (ignore is not None) == (task == "aimd"), "不收敛照跑只在 aimd"


def test_all_electron(cfg, traj, tmp_path):
    tpl = tmp_path / "t.inp"
    tpl.write_text(gapw("sp"))
    s = Cp2kGenerator(cfg).generate(traj, 0, JobOptions("sp"))  # GPW 对照
    ae = Cp2kGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl)
    root = Section.from_str(ae.files["cp2k.inp"])
    kinds = {k.param: (k.get("BASIS_SET"), k.get("POTENTIAL")) for k in root.section("FORCE_EVAL/SUBSYS").sections("KIND")}
    assert set(kinds.values()) == {("pob-TZVP-rev2", "ALL")}
    symbols = ferro.sort_atoms(traj.subset([0])).symbols(0)
    z = {"O": 8, "P": 15, "Zn": 30}
    assert f"{sum(z[x] for x in symbols)} electrons" in " ".join(ae.notes), "全电子按原子序数数电子"
    assert s.files != ae.files


def test_ecp_basis_is_not_all_electron(cfg, tmp_path):
    (tmp_path / "SnO.vasp").write_text(
        "SnO\n1.0\n5 0 0\n0 5 0\n0 0 5\nO Sn\n1 1\nCartesian\n0 0 0\n2.5 2.5 2.5\n")
    tpl = tmp_path / "t.inp"
    tpl.write_text(gapw("sp"))
    with pytest.raises(JobError, match="pob-TZVP-rev2 for Sn is an ECP basis"):
        Cp2kGenerator(cfg).generate(ferro.read(str(tmp_path / "SnO.vasp")), 0, JobOptions("sp"), tpl)


def test_spin_and_kinds_override(cfg, traj, tmp_path):
    s = Cp2kGenerator(cfg).generate(traj, 0, JobOptions("sp", charge=1))
    dft = Section.from_str(s.files["cp2k.inp"]).section("FORCE_EVAL/DFT")
    assert (dft.get("CHARGE"), dft.get("MULTIPLICITY"), dft.get("UKS")) == ("1", "2", "T")
    (tmp_path / "c.toml").write_text(
        f'[cp2k]\ndata_dir = "{tmp_path / "data"}"\n[cp2k.kinds]\nO = {{ basis = "nope" }}\n')
    with pytest.raises(JobError, match="basis nope is not in BASIS_MOLOPT"):
        Cp2kGenerator(config.load(tmp_path / "c.toml")).generate(traj, 0, JobOptions("sp"))
    (tmp_path / "c.toml").write_text(
        f'[cp2k]\ndata_dir = "{tmp_path / "data"}"\nbasis_family = "TZV2P-MOLOPT-GTH"\n')
    with pytest.raises(JobError, match="no TZV2P-MOLOPT-GTH basis for O.*available: DZVP-MOLOPT-SR-GTH"):
        Cp2kGenerator(config.load(tmp_path / "c.toml")).generate(traj, 0, JobOptions("sp"))


def test_kpoints_switch_ot_to_diagonalisation(cfg, tmp_path):
    (tmp_path / "small.vasp").write_text(
        "ZnO\n1.0\n3.2 0 0\n0 3.2 0\n0 0 5.2\nO Zn\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5\n")
    small = ferro.read(str(tmp_path / "small.vasp"))
    s = Cp2kGenerator(cfg).generate(small, 0, JobOptions("sp"))
    root = Section.from_str(s.files["cp2k.inp"])
    assert root.section("FORCE_EVAL/DFT/KPOINTS").get("SCHEME") == "MONKHORST-PACK 4 4 3"
    scf = root.section("FORCE_EVAL/DFT/SCF")
    assert scf.section("OT") is None and scf.section("OUTER_SCF") is None, "OT 只能 Γ 点"
    assert scf.section("DIAGONALIZATION") is not None and scf.section("MIXING").get("METHOD") == "BROYDEN_MIXING"


def test_md_ensembles(cfg, traj):
    s = Cp2kGenerator(cfg).generate(traj, 0, JobOptions("aimd", temperature=1500, steps=10,
                                                        timestep=0.5, ensemble="npt_f"))
    root = Section.from_str(s.files["cp2k.inp"])
    md = root.section("MOTION/MD")
    assert (md.get("ENSEMBLE"), md.get("TEMPERATURE"), md.get("STEPS"), md.get("TIMESTEP")) == \
        ("NPT_F", "1500", "10", "0.5")
    assert md.section("BAROSTAT").get("PRESSURE") == "1.01325"
    assert md.section("BAROSTAT").get("VIRIAL") is None, "VIRIAL 不写，默认全部方向"
    assert root.section("FORCE_EVAL").get("STRESS_TENSOR") == "ANALYTICAL", "NPT 必须解析应力（D-M4）"
    nvt = Section.from_str(Cp2kGenerator(cfg).generate(traj, 0, JobOptions("aimd")).files["cp2k.inp"])
    assert nvt.section("MOTION/MD/BAROSTAT") is None and nvt.section("MOTION/MD").get("ENSEMBLE") == "NVT"


def test_vary_cell_and_errors(cfg, traj, tmp_path):
    root = Section.from_str(Cp2kGenerator(cfg).generate(traj, 0, JobOptions("opt", vary_cell=True)).files["cp2k.inp"])
    assert root.section("GLOBAL").get("RUN_TYPE") == "CELL_OPT"
    assert root.section("MOTION/CELL_OPT").get("TYPE") is None, \
        "CP2K 2026.2 删了 CELL_OPT/TYPE，写了即中止（不变胞的 GEO_OPT 也一样，因为节总在）"
    smear = Cp2kGenerator.builtin_template("sp").replace("      &OT\n", "      &SMEAR ON\n      &END SMEAR\n      &OT\n")
    tpl = tmp_path / "t.inp"
    tpl.write_text(smear)
    with pytest.raises(JobError, match="OT"):
        Cp2kGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl)
    s = Cp2kGenerator(cfg).generate(traj, 0, JobOptions("sp", fine_grid=2))
    grid = Section.from_str(s.files["cp2k.inp"]).section("FORCE_EVAL/DFT/XC/XC_GRID")
    assert grid.get("USE_FINER_GRID") == "T", "CP2K 的 --fine-grid 打开 USE_FINER_GRID"
    assert any("probably unnecessary" in n for n in s.notes), "PBE 不是 meta-GGA，提示多半用不着"
