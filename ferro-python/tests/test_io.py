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
