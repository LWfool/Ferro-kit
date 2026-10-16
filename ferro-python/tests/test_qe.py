"""QE 生成器：pseudo_dir 里的 UPF 建索引（不靠 json）、截断取建议值、namelist 与卡片。"""

from pathlib import Path

import pytest

import ferro
from ferro import config
from ferro.jobs.base import JobError, JobOptions
from ferro.jobs.qe import QeGenerator, QeInput, upf_header

NVT = Path(__file__).resolve().parents[2] / "tests" / "70Z30P00A_NVT_5.lammpstrj"


def upf_v2(element, z, wfc, rho):
    return (f'<UPF version="2.0.1">\n<PP_INFO>x</PP_INFO>\n<PP_HEADER generated="test"\n'
            f'   element="{element:2s}" pseudo_type="PAW" z_valence="{z:.15E}"\n'
            f'   wfc_cutoff="{wfc:.15E}" rho_cutoff="{rho:.15E}"/>\n<PP_MESH/>\n</UPF>\n')


UPF_V1 = """\
<PP_INFO>
</PP_INFO>
<PP_HEADER>
   0                   Version Number
  Zn                   Element
   US                  Ultrasoft pseudopotential
    T                  Nonlinear Core Correction
 SLA PW PBX PBC    PBE Exchange-Correlation functional
   12.00000000000      Z valence
  -47.0                Total energy
    0.0000000    0.0000000 Suggested cutoff for wfc and rho
</PP_HEADER>
"""


@pytest.fixture
def cfg(tmp_path):
    pp = tmp_path / "upf"
    pp.mkdir()
    (pp / "O.pbe-n-kjpaw_psl.0.1.UPF").write_text(upf_v2("O", 6, 50, 400))
    (pp / "P.pbe-n-rrkjus_psl.1.0.0.UPF").write_text(upf_v2("P", 5, 30, 240))
    (pp / "zn_pbe_v1.uspp.F.UPF").write_text(UPF_V1)
    (pp / "README").write_text("not a pseudopotential\n")
    (tmp_path / "c.toml").write_text(f'[qe]\npseudo_dir = "{pp}"\n')
    return config.load(tmp_path / "c.toml")


@pytest.fixture(scope="module")
def traj():
    return ferro.read(str(NVT))


def test_upf_headers(cfg, tmp_path):
    assert upf_header(tmp_path / "upf" / "O.pbe-n-kjpaw_psl.0.1.UPF") == ("O", 6.0, 50.0, 400.0)
    assert upf_header(tmp_path / "upf" / "zn_pbe_v1.uspp.F.UPF") == ("Zn", 12.0, 0.0, 0.0), "UPF v1"


@pytest.mark.parametrize("task", ["sp", "opt", "aimd"])
def test_builtin_templates_generate(cfg, traj, task, tmp_path):
    s = QeGenerator(cfg).generate(traj, 0, JobOptions(task))
    qe = QeInput.from_str(s.files["pw.in"])
    one = ferro.sort_atoms(traj.subset([0]))
    assert qe.get("system", "nat") == str(one.n_atoms()) and qe.get("system", "ntyp") == "3"
    assert (qe.get("system", "ecutwfc"), qe.get("system", "ecutrho")) == ("50.0", "400.0"), "取 UPF 建议的最大值"
    assert qe.get("control", "pseudo_dir") == f"'{tmp_path / 'upf'}/'"
    assert "zn_pbe_v1.uspp.F.UPF" in qe.card("ATOMIC_SPECIES"), "按头部的元素归类，不按文件名"
    assert qe.card("K_POINTS").split() == ["K_POINTS", "automatic", "1", "1", "1", "0", "0", "0"]
    path = tmp_path / "pw.in"
    path.write_text(s.files["pw.in"])
    assert ferro.read(str(path)).symbols(0) == one.symbols(0), "生成的 pw.in 能被 ferro 读回"


def test_charge_spin_and_md(cfg, traj):
    s = QeGenerator(cfg).generate(traj, 0, JobOptions("aimd", charge=1, temperature=1500, steps=10, timestep=2.0))
    qe = QeInput.from_str(s.files["pw.in"])
    assert (qe.get("system", "tot_charge"), qe.get("system", "nspin"), qe.get("system", "tot_magnetization")) == \
        ("1", "2", "1")
    assert (qe.get("ions", "tempw"), qe.get("control", "nstep"), qe.get("control", "dt")) == ("1500", "10", "41.3411")
    with pytest.raises(JobError, match="NVT only"):
        QeGenerator(cfg).generate(traj, 0, JobOptions("aimd", ensemble="npt_i"))


def test_vary_cell_and_kpoints(cfg, traj, tmp_path):
    qe = QeInput.from_str(QeGenerator(cfg).generate(traj, 0, JobOptions("opt", vary_cell=True)).files["pw.in"])
    assert qe.get("control", "calculation") == "'vc-relax'" and qe.get("cell", "cell_dynamics") == "'bfgs'"
    assert list(qe.namelists).index("cell") > list(qe.namelists).index("ions"), "&CELL 在 &IONS 之后"
    (tmp_path / "small.vasp").write_text(
        "ZnO\n1.0\n3.2 0 0\n0 3.2 0\n0 0 5.2\nO Zn\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5\n")
    small = QeInput.from_str(QeGenerator(cfg).generate(ferro.read(str(tmp_path / "small.vasp")), 0,
                                                       JobOptions("sp")).files["pw.in"])
    assert small.card("K_POINTS").split()[2:5] == ["4", "4", "3"]
    tpl = tmp_path / "pw.in"
    tpl.write_text(QeGenerator.builtin_template("sp") + "K_POINTS gamma\n")
    kept = QeGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl)
    assert QeInput.from_str(kept.files["pw.in"]).card("K_POINTS").strip() == "K_POINTS gamma", "模板写死的不动"
    s = QeGenerator(cfg).generate(traj, 0, JobOptions("sp", fine_grid=2))
    assert any("--fine-grid ignored" in n for n in s.notes)


def test_pseudo_choice(cfg, traj, tmp_path):
    pp = tmp_path / "upf"
    (pp / "O.pbe-rrkjus.UPF").write_text(upf_v2("O", 6, 40, 320))
    with pytest.raises(JobError, match=r"2 UPF files for O.*choose one in \[qe.pseudos\]"):
        QeGenerator(cfg).generate(traj, 0, JobOptions("sp"))
    (tmp_path / "c.toml").write_text(f'[qe]\npseudo_dir = "{pp}"\n[qe.pseudos]\nO = "O.pbe-rrkjus.UPF"\n')
    s = QeGenerator(config.load(tmp_path / "c.toml")).generate(traj, 0, JobOptions("sp"))
    assert "O.pbe-rrkjus.UPF" in QeInput.from_str(s.files["pw.in"]).card("ATOMIC_SPECIES")
    (pp / "P.pbe-n-rrkjus_psl.1.0.0.UPF").unlink()
    with pytest.raises(JobError, match="no UPF files for P"):
        QeGenerator(config.load(tmp_path / "c.toml")).generate(traj, 0, JobOptions("sp"))


def test_cutoffs_needed(cfg, traj, tmp_path):
    (tmp_path / "upf" / "O.pbe-n-kjpaw_psl.0.1.UPF").write_text(upf_v2("O", 6, 0, 0))
    (tmp_path / "upf" / "P.pbe-n-rrkjus_psl.1.0.0.UPF").write_text(upf_v2("P", 5, 0, 0))
    with pytest.raises(JobError, match="suggest no ecutwfc"):
        QeGenerator(cfg).generate(traj, 0, JobOptions("sp"))
    tpl = tmp_path / "pw.in"
    tpl.write_text(QeGenerator.builtin_template("sp").replace("  ibrav = 0\n", "  ibrav = 0\n  ecutwfc = 60\n  ecutrho = 480\n"))
    qe = QeInput.from_str(QeGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl).files["pw.in"])
    assert qe.get("system", "ecutwfc") == "60", "模板写了截断就沿用"
