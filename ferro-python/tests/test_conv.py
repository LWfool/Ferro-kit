"""收敛测试：扫描参数的解析、生成的目录与清单、读回输出后的偏差计算。"""

import json
import math
import shutil
from pathlib import Path

import pytest

import ferro
from ferro import cli, config, jobs
from ferro.jobs import conv
from ferro.jobs.base import JobError, JobOptions
from ferro.jobs.cp2k import Cp2kGenerator, Section
from ferro.jobs.vasp import Incar, VaspGenerator

FIXTURES = Path(__file__).resolve().parents[2] / "tests"
NVT = FIXTURES / "70Z30P00A_NVT_5.lammpstrj"
SP = [FIXTURES / "vasp651_sp" / f"frame_{n}" / "OUTCAR" for n in ("0011", "0161", "0194")]


@pytest.mark.parametrize("software, text, expected", [
    ("vasp", "encut=400:600:100", ("encut", [400.0, 500.0, 600.0])),
    ("vasp", "encut=700,500", ("encut", [500.0, 700.0])),
    ("vasp", "fine-grid=2,3", ("fine-grid", [2.0, 3.0])),
    ("cp2k", "finer-grid", ("finer-grid", ["off", "on"])),
    ("cp2k", "rel-cutoff=40:80:20", ("rel-cutoff", [40.0, 60.0, 80.0])),
])
def test_parse_scan(software, text, expected):
    assert conv.parse_scan(software, text) == expected


@pytest.mark.parametrize("software, text, message", [
    ("qe", "ecutwfc=30,40", "offered for vasp, cp2k"),
    ("vasp", "cutoff=300,400", "vasp scans encut, fine-grid"),
    ("vasp", "encut=400", "at least two distinct"),
    ("vasp", "encut=a:b:c", "START:STOP:STEP"),
    ("cp2k", "finer-grid=on", "is a switch"),
])
def test_parse_scan_errors(software, text, message):
    with pytest.raises(JobError, match=message):
        conv.parse_scan(software, text)


@pytest.fixture
def vasp_cfg(tmp_path, monkeypatch):
    pp = tmp_path / "pp"
    for name, z in [("O", 6), ("P", 5), ("Zn", 12)]:
        (pp / name).mkdir(parents=True)
        (pp / name / "POTCAR").write_text(f" PAW_PBE {name}\n ZVAL = {z}.0\n ENMAX = 400.0\n")
    (tmp_path / "c.toml").write_text(f'[vasp]\npotcar_dir = "{pp}"\n')
    monkeypatch.setenv("FERRO_CONFIG", str(tmp_path / "c.toml"))
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    return config.load(tmp_path / "c.toml")


def test_vasp_encut_scan_layout(vasp_cfg, tmp_path):
    out = tmp_path / "encut"
    assert cli.main(["conv", "-s", "vasp", "--scan", "encut=500:700:100", "-i", str(NVT), "--number", "2",
                     "-o", str(out), "--mkdir"]) == 0
    manifest = json.loads((out / conv.MANIFEST).read_text())
    assert manifest["values"] == [500.0, 600.0, 700.0] and len(manifest["jobs"]) == 6
    job = out / "encut_600" / "70Z30P00A_NVT_5_0004"
    assert Incar.from_str((job / "INCAR").read_text())["ENCUT"] == "600"
    assert Incar.from_str((job / "INCAR").read_text())["IBRION"] == "-1", "收敛测试用 sp 模板"
    assert (job / "POTCAR").is_symlink() and str((job / "POTCAR").readlink()) == "../../POTCAR", \
        "嵌套两层时链接仍是相对路径"
    assert (job / "POTCAR").read_text() == (out / "POTCAR").read_text()


def test_cp2k_switch_scans(tmp_path):
    text = Cp2kGenerator.builtin_template("sp")
    apply = conv.SCANS["cp2k"]["finer-grid"][0]
    on, _ = apply(text, "on", JobOptions("sp"))
    off, _ = apply(on, "off", JobOptions("sp"))
    assert Section.from_str(on).section("FORCE_EVAL/DFT/XC/XC_GRID").get("USE_FINER_GRID") == "T"
    assert Section.from_str(off).section("FORCE_EVAL/DFT/XC/XC_GRID").get("USE_FINER_GRID") is None
    cut, _ = conv.SCANS["cp2k"]["cutoff"][0](text, 600.0, JobOptions("sp"))
    assert Section.from_str(cut).section("FORCE_EVAL/DFT/MGRID").get("CUTOFF") == "600"
    with pytest.raises(JobError, match="only acts with &QS METHOD GAPW"):
        conv.SCANS["cp2k"]["accurate-xcint"][0](text, "on", JobOptions("sp"))


def fake_test(tmp_path, values=(400.0, 500.0, 600.0)):
    """三份真实 OUTCAR 当作同一结构的三档（原子数相同，数值上可比）。"""
    manifest = {"software": "vasp", "parameter": "encut", "unit": "eV", "values": list(values), "jobs": []}
    for value, outcar in zip(values, SP):
        d = tmp_path / f"encut_{value:g}" / "s"
        d.mkdir(parents=True)
        shutil.copy(outcar, d / "OUTCAR")
        manifest["jobs"].append({"dir": f"encut_{value:g}/s", "value": value, "structure": "s", "frame": 0})
    (tmp_path / conv.MANIFEST).write_text(json.dumps(manifest))
    return manifest


def test_compare_matches_an_independent_calculation(tmp_path):
    manifest = fake_test(tmp_path)
    rows, ref = conv.compare(manifest, tmp_path)
    assert ref == 600.0, "默认参考最高档"
    t = [ferro.read_output(str(p))[0] for p in SP]
    n = t[0].n_atoms()
    for row, x in zip(rows, t):
        d = [math.dist(a, b) * 1000 for a, b in zip(x.forces(0), t[2].forces(0))]
        assert row["dE_meV_atom"] == pytest.approx((x.energy(0) - t[2].energy(0)) / n * 1000)
        assert row["dF_max_meV_A"] == pytest.approx(max(d))
        assert row["dF_rms_meV_A"] == pytest.approx(math.sqrt(sum(v * v for v in d) / n))
        p = [sum(s[i][i] for i in range(3)) / 3 * conv.EV_A3_TO_KBAR for s in (x.stress(0), t[2].stress(0))]
        assert row["dP_kbar"] == pytest.approx(p[0] - p[1])
    assert rows[2]["within"] == "yes" and rows[2]["dF_rms_meV_A"] == 0, "参考档与自己比偏差为 0"
    rows, ref = conv.compare(manifest, tmp_path, ref=400.0)
    assert rows[0]["dE_meV_atom"] == 0
    with pytest.raises(JobError, match="not one of the scanned values"):
        conv.compare(manifest, tmp_path, ref=450.0)


def test_read_result_finds_outputs(tmp_path):
    shutil.copy(FIXTURES / "cp2k_sp_v61.out", tmp_path / "cp2k.out")
    energy, forces, pressure, status = conv.read_result(tmp_path, "cp2k")
    assert status == "ok" and isinstance(energy, float) and forces and pressure is not None
    assert conv.read_result(tmp_path / "missing", "vasp")[3] == "no output"
    (tmp_path / "OUTCAR").write_text("not an outcar\n")
    assert conv.read_result(tmp_path, "vasp")[3].startswith("unreadable"), "读不了的输出说明原因"


def test_analyze_cli(tmp_path, capsys):
    fake_test(tmp_path)
    shutil.rmtree(tmp_path / "encut_500" / "s")
    (tmp_path / "encut_500" / "s").mkdir()
    assert cli.main(["conv", "analyze", str(tmp_path), "--df-rms", "1e9", "--df-max", "1e9",
                     "--de", "1e9", "--dp", "1e9"]) == 1, "有任务缺输出时退出码 1"
    printed = capsys.readouterr().out
    assert "500 s: no output" in printed
    assert "First value within every threshold: 400" in printed
    csv = (tmp_path / "conv_encut.csv").read_text().splitlines()
    assert csv[0].startswith("# ferro-py conv analyze: vasp encut [eV]")
    header = csv[3].split(",")
    assert header[:4] == ["value", "structure", "frame", "status"] and "dF_rms_meV_A" in header
    assert len(csv) == 4 + 3
