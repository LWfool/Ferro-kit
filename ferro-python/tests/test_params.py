"""参数块：INCAR、pw.x namelist + 卡片、CP2K 节树的解析、编辑与写回。"""

from pathlib import Path

import pytest

import ferro
from ferro.jobs.cp2k import Cp2kInputError, Section
from ferro.jobs.qe import QeInput, QeInputError, format_value, unquote
from ferro.jobs.vasp import Incar, IncarError

NVT = Path(__file__).resolve().parents[2] / "tests" / "70Z30P00A_NVT_5.lammpstrj"

INCAR = """\
# Parallel
NCORE = 10         ! sqrt(COREs)
ISTART = 0; ICHARG = 2   # restart
ENCUT = 800

SYSTEM = "two
lines"
EDIFF = 1.0E-05
"""


def test_incar_round_trip_is_byte_identical():
    assert Incar.from_str(INCAR).to_str() == INCAR, "不改动时应逐字节写回"


def test_incar_edit_keeps_layout():
    incar = Incar.from_str(INCAR)
    assert incar["encut"] == "800" and incar["ICHARG"] == "2", "键不分大小写"
    assert incar["SYSTEM"] == '"two\nlines"', "引号内的多行值是一个值"
    incar["NCORE"] = 4
    incar["ICHARG"] = 1
    incar["LWAVE"] = False
    del incar["ISTART"]
    del incar["EDIFF"]
    lines = incar.to_str().splitlines()
    assert lines[1] == "NCORE = 4         ! sqrt(COREs)", "只换值，行尾注释留着"
    assert lines[2].split("#")[0].split() == ["ICHARG", "=", "1"], "删去同一行的另一个赋值"
    assert lines[-1] == "LWAVE = .FALSE.", "新键追加在末尾"
    assert "EDIFF" not in incar and "EDIFF" not in incar.to_str()


@pytest.mark.parametrize("text, message", [
    ("ENCUT = 400\nencut = 500\n", "sets ENCUT twice"),
    ("ENCUT 400\n", "line 1: not a KEY = value"),
    ('SYSTEM = "open\n', "inside a quoted value"),
])
def test_incar_errors(text, message):
    with pytest.raises(IncarError, match=message):
        Incar.from_str(text)


PW = """\
&control
    calculation = 'scf'
    pseudo_dir = './pseudo/'       ! Must end with the '/'
/
&system
    ibrav = 0, ecutwfc = 50
/
&electrons
    conv_thr =  1e-10
/

K_POINTS automatic
2 2 2  0 0 0

ATOMIC_SPECIES
 O 15.999 O.UPF
"""


def test_qe_edit_and_cards():
    qe = QeInput.from_str(PW)
    assert qe.get("CONTROL", "Pseudo_Dir") == "'./pseudo/'", "名与键不分大小写，值是原文"
    assert unquote(qe.get("control", "calculation")) == "scf"
    qe.set("control", "calculation", "md")
    qe.set("system", "nosym", True)
    qe.set("ions", "ion_temperature", "rescaling")
    qe.delete("electrons", "conv_thr")
    assert qe.card("K_POINTS").splitlines()[:2] == ["K_POINTS automatic", "2 2 2  0 0 0"]
    qe.set_card("ATOMIC_SPECIES", "ATOMIC_SPECIES\n O 15.999 O.pbe.UPF\n")
    qe.set_card("CELL_PARAMETERS", "CELL_PARAMETERS angstrom\n 1 0 0\n 0 1 0\n 0 0 1")
    with pytest.raises(QeInputError, match="must start with OCCUPATIONS"):
        qe.set_card("OCCUPATIONS", "K_POINTS gamma")
    again = QeInput.from_str(qe.to_str())
    assert again.namelists == {
        "control": {"calculation": "'md'", "pseudo_dir": "'./pseudo/'"},
        "system": {"ibrav": "0", "ecutwfc": "50", "nosym": ".true."},
        "electrons": {},
        "ions": {"ion_temperature": "'rescaling'"},
    }, "写回再解析得到同一份 namelist"
    assert "O.pbe.UPF" in again.card("ATOMIC_SPECIES")
    assert again.card("K_POINTS") == qe.card("K_POINTS"), "没替换的卡片原样保留"
    assert again.to_str() == qe.to_str(), "写回是稳定的"


def test_qe_values():
    assert [format_value(v) for v in (False, 3, 1e-10, "x")] == [".false.", "3", "1e-10", "'x'"]
    with pytest.raises(QeInputError, match="single quote"):
        format_value("it's")
    with pytest.raises(QeInputError, match="not closed"):
        QeInput.from_str("&control\n calculation='scf'\n")


def test_qe_structure_cards_drop_in(tmp_path):
    t = ferro.sort_atoms(ferro.read(str(NVT)).subset([0]))
    qe = QeInput.from_str(PW)
    qe.set("system", "nat", t.n_atoms())
    qe.set("system", "ntyp", len(t.elements()))
    cards = ferro.qe_structure_cards(t, {e: f"{e}.UPF" for e in t.elements()})
    for block in cards.strip().split("\n\n"):
        qe.set_card(block.split()[0], block)
    path = tmp_path / "pw.in"
    path.write_text(qe.to_str())
    assert ferro.read(str(path)).symbols(0) == t.symbols(0), "替换进模板的结构卡片可被读回"


CP2K = """\
@SET PROJECT x
&GLOBAL
  PROJECT ${PROJECT}
  RUN_TYPE MD
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    # whole-line comment
    &QS
      EPS_DEFAULT 5.0E-11  #Set all EPS_xxx

    &END QS
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &KIND O
      BASIS_SET DZVP-MOLOPT-SR-GTH
    &END KIND
    &KIND Zr
      BASIS_SET DZVP-MOLOPT-SR-GTH
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
"""


def test_cp2k_round_trip_and_lookup():
    root = Section.from_str(CP2K)
    assert root.to_str() == CP2K, "已是两格缩进的输入应逐字节写回"
    messy = "\n".join("\t " + line.strip().lower() if line.startswith("  ") else line
                      for line in CP2K.splitlines())
    assert Section.from_str(Section.from_str(messy).to_str()).to_str() == Section.from_str(messy).to_str()
    assert root.section("force_eval/dft/QS").get("eps_default") == "5.0E-11", "不分大小写"
    assert root.section("FORCE_EVAL/DFT/XC/XC_FUNCTIONAL").param == "PBE"
    assert [k.param for k in root.section("FORCE_EVAL/SUBSYS").sections("KIND")] == ["O", "Zr"]
    assert root.section("MOTION") is None


def test_cp2k_edit():
    root = Section.from_str(CP2K)
    qs = root.section("FORCE_EVAL/DFT/QS")
    qs.set("EPS_DEFAULT", "1.0E-12")
    qs.set("EXTRAPOLATION", "ASPC")
    md = root.section("MOTION/MD", create=True)
    md.set("STEPS", 100)
    md.section("THERMOSTAT/NOSE", create=True).set("TIMECON", 100.0)
    md.set("ENSEMBLE", "NVT")
    root.section("FORCE_EVAL").set("STRESS_TENSOR", "ANALYTICAL")
    root.section("FORCE_EVAL/SUBSYS").remove("KIND")
    text = root.to_str()
    assert "      EPS_DEFAULT 1.0E-12  #Set all EPS_xxx" in text, "改值保留行尾注释"
    assert text.index("ENSEMBLE NVT") < text.index("&THERMOSTAT"), "新关键字放在子节之前"
    assert text.index("STRESS_TENSOR ANALYTICAL") < text.index("&DFT")
    assert "&KIND" not in text
    assert Section.from_str(text).to_str() == text


@pytest.mark.parametrize("text, message", [
    ("&GLOBAL\n&END MOTION\n", "closes &GLOBAL"),
    ("&GLOBAL\n", "not closed"),
    ("&END\n", "without an open section"),
])
def test_cp2k_errors(text, message):
    with pytest.raises(Cp2kInputError, match=message):
        Section.from_str(text)
