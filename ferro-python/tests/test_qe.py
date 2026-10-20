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
    assert any("parallel settings" in n for n in s.notes), "并行参数模板不写，生成时提示用户自己加"
    assert "parallelisation" in QeGenerator.builtin_template(task), "模板头注释也写明"


def test_charge_spin_and_md(cfg, traj):
    s = QeGenerator(cfg).generate(traj, 0, JobOptions("aimd", charge=1, temperature=1500, steps=10, timestep=2.0))
    qe = QeInput.from_str(s.files["pw.in"])
    assert (qe.get("system", "tot_charge"), qe.get("system", "nspin"), qe.get("system", "tot_magnetization")) == \
        ("1", "2", "1")
    assert (qe.get("ions", "tempw"), qe.get("control", "nstep"), qe.get("control", "dt")) == ("1500", "10", "41.3414")
    assert qe.get("system", "nosym") == ".true.", "MD 不关对称性，热浴给的速度破坏对称后 checkallsym 中止"
    assert any("nraise = 150 steps = 300 fs" in n for n in s.notes), "nraise 以步计，dt 2 fs 下是 300 fs"
    with pytest.raises(JobError, match="NVT only"):
        QeGenerator(cfg).generate(traj, 0, JobOptions("aimd", ensemble="npt_i"))


def test_explicit_singlet_with_spin_polarised_template(cfg, traj, tmp_path):
    # 模板 nspin = 2：显式 --multiplicity 1 固定总磁矩为 0；不给则磁矩照模板放开
    tpl = tmp_path / "pw.in"
    tpl.write_text(QeGenerator.builtin_template("sp").replace("  ibrav = 0\n", "  ibrav = 0\n  nspin = 2\n"))
    s = QeGenerator(cfg).generate(traj, 0, JobOptions("sp", multiplicity=1), tpl)
    qe = QeInput.from_str(s.files["pw.in"])
    assert (qe.get("system", "nspin"), qe.get("system", "tot_magnetization")) == ("2", "0")
    free = QeInput.from_str(QeGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl).files["pw.in"])
    assert (free.get("system", "nspin"), free.get("system", "tot_magnetization")) == ("2", None)


def test_vary_cell_and_kpoints(cfg, traj, tmp_path):
    with pytest.raises(JobError, match="--vary-cell needs ecutwfc in the template"):
        QeGenerator(cfg).generate(traj, 0, JobOptions("opt", vary_cell=True))
    tpl = tmp_path / "opt.in"
    tpl.write_text(QeGenerator.builtin_template("opt").replace("&SYSTEM", "&SYSTEM\n  ecutwfc = 70"))
    qe = QeInput.from_str(QeGenerator(cfg).generate(traj, 0, JobOptions("opt", vary_cell=True), tpl).files["pw.in"])
    assert qe.get("system", "ecutwfc") == "70" and float(qe.get("system", "ecutrho")) >= 280, \
        "ecutwfc from the template; ecutrho still falls back (UPF suggestion, at least 4 x ecutwfc)"
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


def test_cell_comes_from_cell_parameters(cfg, traj, tmp_path):
    tpl = tmp_path / "pw.in"
    tpl.write_text(QeGenerator.builtin_template("sp").replace(
        "  ibrav = 0\n", "  ibrav = 2\n  celldm(1) = 10.2\n  A = 5.0\n  cosAB = 0.5\n"))
    qe = QeInput.from_str(QeGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl).files["pw.in"])
    system = qe.namelists["system"]
    assert system["ibrav"] == "0", "ibrav ≠ 0 时 pw.x 不用 CELL_PARAMETERS，胞会静默换成模板的"
    assert not {"celldm(1)", "a", "cosab"} & set(system), "胞参数键一律删掉"
    tpl.write_text(QeGenerator.builtin_template("sp").replace("  ibrav = 0\n", ""))
    qe = QeInput.from_str(QeGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl).files["pw.in"])
    assert qe.get("system", "ibrav") == "0", "模板没写 ibrav（INPUT_PW 必填）也补上"


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


def test_upf_without_element_is_skipped(cfg, traj, tmp_path):
    pp = tmp_path / "upf"
    (pp / "1_unknown.UPF").write_text('<UPF version="2.0.1">\n<PP_HEADER z_valence="1.0"/>\n</UPF>\n')
    QeGenerator(cfg).generate(traj, 0, JobOptions("sp"))
    (pp / "P.pbe-n-rrkjus_psl.1.0.0.UPF").unlink()
    with pytest.raises(JobError, match="no UPF files for P.*skipped, element unknown: 1_unknown.UPF"):
        QeGenerator(cfg).generate(traj, 0, JobOptions("sp"))


def test_nraise_must_be_whole(cfg, traj):
    template = QeGenerator.builtin_template("aimd").replace("nraise = 150", "nraise = 10.0")
    with pytest.raises(JobError, match="nraise = 10.0 in the template"):
        QeGenerator(cfg).generate(traj, 0, JobOptions("aimd", timestep=2.0), template_text=template)


def test_cutoffs_needed(cfg, traj, tmp_path):
    (tmp_path / "upf" / "O.pbe-n-kjpaw_psl.0.1.UPF").write_text(upf_v2("O", 6, 0, 0))
    (tmp_path / "upf" / "P.pbe-n-rrkjus_psl.1.0.0.UPF").write_text(upf_v2("P", 5, 0, 0))
    with pytest.raises(JobError, match="suggest no ecutwfc"):
        QeGenerator(cfg).generate(traj, 0, JobOptions("sp"))
    tpl = tmp_path / "pw.in"
    tpl.write_text(QeGenerator.builtin_template("sp").replace("  ibrav = 0\n", "  ibrav = 0\n  ecutwfc = 60\n  ecutrho = 480\n"))
    qe = QeInput.from_str(QeGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl).files["pw.in"])
    assert qe.get("system", "ecutwfc") == "60", "模板写了截断就沿用"
    # wfc suggested, rho not: ecutrho falls back to 4 x ecutwfc instead of an error
    (tmp_path / "upf" / "O.pbe-n-kjpaw_psl.0.1.UPF").write_text(upf_v2("O", 6, 50, 0))
    s = QeGenerator(cfg).generate(traj, 0, JobOptions("sp"))
    qe = QeInput.from_str(s.files["pw.in"])
    assert (qe.get("system", "ecutwfc"), qe.get("system", "ecutrho")) == ("50.0", "200.0")
    assert "ecutrho = 200 Ry (4 x ecutwfc; the UPF headers suggest none)" in s.notes


@pytest.mark.parametrize("written, rho, note", [
    ("  ecutwfc = 120\n", "480.0", "ecutrho = 480 Ry (4 x ecutwfc"),     # UPF 建议 400 < 4 × 120
    ("  ecutwfc = 120\n  ecutrho = 300\n", "300", "ecutrho 300 Ry is below 4 x ecutwfc"),
    ("  ecutwfc = 80\n", "400.0", None),                                     # 400 ≥ 4 × 80，不动
])
def test_ecutrho_against_four_times_ecutwfc(cfg, traj, tmp_path, written, rho, note):
    tpl = tmp_path / "pw.in"
    tpl.write_text(QeGenerator.builtin_template("sp").replace("  ibrav = 0\n", "  ibrav = 0\n" + written))
    s = QeGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl)
    assert QeInput.from_str(s.files["pw.in"]).get("system", "ecutrho") == rho
    if note is None:
        assert not any("4 x ecutwfc" in n for n in s.notes)
    else:
        assert any(note in n for n in s.notes), f"该有提示 {note!r}: {s.notes}"


def test_ecutrho_not_above_ecutwfc_is_refused(cfg, traj, tmp_path):
    # pw.x 的 set_cutoff：dual ≤ 1 才 errore；dual < 4 只 infomsg
    tpl = tmp_path / "pw.in"
    tpl.write_text(QeGenerator.builtin_template("sp").replace(
        "  ibrav = 0\n", "  ibrav = 0\n  ecutwfc = 60\n  ecutrho = 60\n"))
    with pytest.raises(JobError, match="ecutrho must exceed ecutwfc"):
        QeGenerator(cfg).generate(traj, 0, JobOptions("sp"), tpl)
