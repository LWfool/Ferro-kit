"""读计算输出（能量、力、应力）与写表：收敛测试 ferro-py conv 依赖的 Rust 接口。"""

import pytest

import ferro

from conftest import FIXTURES


@pytest.mark.parametrize("path, fmt", [
    (FIXTURES / "vasp651_sp" / "frame_0011" / "OUTCAR", "VASP OUTCAR"),
    (FIXTURES / "vasp651_sp" / "frame_0011" / "vasprun.xml", "VASP vasprun.xml"),
    (FIXTURES / "cp2k_sp_v61.out", "CP2K single-point output"),
])
def test_single_point_outputs(path, fmt):
    traj, stats = ferro.read_output(str(path))
    assert stats["format"] == fmt and stats["kept"] == len(traj) == 1
    assert isinstance(traj.energy(0), float)
    forces = traj.forces(0)
    assert len(forces) == traj.n_atoms(), "每个原子一组力"
    stress = traj.stress(0)
    assert len(stress) == 3 and all(len(r) == 3 for r in stress)
    assert max(abs(stress[i][j] - stress[j][i]) for i in range(3) for j in range(3)) < 1e-6, "应力张量对称"


def test_outcar_and_vasprun_agree():
    a, _ = ferro.read_output(str(FIXTURES / "vasp651_sp" / "frame_0011" / "OUTCAR"))
    b, _ = ferro.read_output(str(FIXTURES / "vasp651_sp" / "frame_0011" / "vasprun.xml"))
    assert abs(a.energy(0) - b.energy(0)) < 1e-5, "同一次计算的两种输出给同一个能量"
    assert max(abs(x - y) for f, g in zip(a.forces(0), b.forces(0)) for x, y in zip(f, g)) < 1e-4


def test_structure_files_have_no_results():
    t = ferro.read(str(FIXTURES / "70Z30P00A_NVT_5.lammpstrj"))
    assert t.energy(0) is None and t.stress(0) is None


def test_write_table(tmp_path):
    path = tmp_path / "t.csv"
    ferro.write_table(str(path), [("label", ["a", "b"]), ("de", [1.5, float("nan")])], meta=["ref: b"])
    assert path.read_text() == "# ref: b\nlabel,de\na,1.500000e0\nb,\n", "NaN 写成空字段，数值格式同 ferro 的其他产物"
    with pytest.raises(RuntimeError, match="refusing to write"):
        ferro.write_table(str(path), [("a", [1.0]), ("b", [1.0, 2.0])])
    with pytest.raises(RuntimeError, match="all numbers or all strings"):
        ferro.write_table(str(path), [("a", [1.0, "x"])])


def test_read_output_carries_every_stat():
    _, stats = ferro.read_output(str(FIXTURES / "cp2k_sp_v61.out"))
    for key in ("layout_drift", "restarts", "version_note", "step_range"):
        assert key in stats, key
    # `steps` stays the count of steps the file claims; the span has its own key
    assert isinstance(stats["steps"], int) and stats["steps"] >= stats["kept"]
