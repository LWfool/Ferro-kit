"""read / write 的扩展名分派与 CLI 的 io_dispatch 一致。"""

from pathlib import Path

import pytest

import ferro

NVT = Path(__file__).resolve().parents[2] / "tests" / "70Z30P00A_NVT_5.lammpstrj"


@pytest.mark.parametrize("ext", ["data", "lmp", "lammps"])
def test_lammps_extension_is_data_file(tmp_path, ext):
    t = ferro.read(str(NVT)).tail(1)
    path = tmp_path / f"one.{ext}"
    ferro.write(t, str(path))
    assert "atoms" in path.read_text(), f".{ext} 应写成 LAMMPS data"
    with pytest.raises(RuntimeError, match="atom_style"):
        ferro.read(str(path))
    back = ferro.read(str(path), atom_style="full")
    assert back.symbols(0) == t.symbols(0), f".{ext} 应按 data 读回原元素"
    with pytest.raises(RuntimeError, match="LAMMPS dump input only"):
        ferro.read(str(path), atom_style="full", type_map={1: "P"})


@pytest.mark.parametrize("name", ["one.vasp", "one.poscar", "one.pos", "POSCAR_1", "CONTCAR"])
def test_vasp_names_are_poscar(tmp_path, name):
    t = ferro.read(str(NVT)).tail(1)
    path = tmp_path / name
    ferro.write(t, str(path))
    assert path.read_text() == ferro.poscar_string(t), f"{name} 应写成 POSCAR"
    assert ferro.read(str(path)).symbols(0) == t.symbols(0), f"{name} 应按 POSCAR 读回"


@pytest.mark.parametrize("name, made", [
    ("conf.cif", ["conf_0000.cif", "conf_0002.cif"]),
    ("conf.vasp", ["conf_0000.vasp", "conf_0002.vasp"]),
    ("POSCAR", ["POSCAR_0000", "POSCAR_0002"]),
    ("conf.data", ["conf_0000.data", "conf_0002.data"]),
])
def test_single_structure_formats_write_one_file_per_frame(tmp_path, name, made):
    t = ferro.read(str(NVT)).tail(3)
    ferro.write(t, str(tmp_path / name))
    assert not (tmp_path / name).exists(), "the unindexed name must not be written"
    files = sorted(p.name for p in tmp_path.iterdir())
    assert len(files) == 3 and files[0] == made[0] and files[-1] == made[1], files
    if name.endswith(".cif"):
        assert (tmp_path / made[0]).read_text().count("data_") == 1


def test_trajectory_formats_write_one_file(tmp_path):
    t = ferro.read(str(NVT)).tail(3)
    ferro.write(t, str(tmp_path / "all.extxyz"))
    assert [p.name for p in tmp_path.iterdir()] == ["all.extxyz"]
    assert ferro.read(str(tmp_path / "all.extxyz")).n_frames() == 3


def test_reader_warnings_are_python_warnings(tmp_path):
    # A dump cut inside its last frame: the frame is dropped and the reader says so
    # through the warnings module, not by printing to stderr behind the caller's back
    lines = NVT.read_text().splitlines(keepends=True)
    starts = [i for i, l in enumerate(lines) if l.startswith("ITEM: TIMESTEP")]
    path = tmp_path / "cut.lammpstrj"
    path.write_text("".join(lines[: starts[-1] + 5]))
    with pytest.warns(UserWarning, match="last frame"):
        t = ferro.read(str(path))
    assert len(t) == len(starts) - 1
