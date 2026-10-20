"""read / write 的扩展名分派与 CLI 的 io_dispatch 一致。"""

import warnings
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


def test_reader_notices_name_their_file(tmp_path):
    # The same notice from two files is two warnings: the path leads the text, so Python's
    # once-per-(text, line) rule no longer folds a loop over inputs into one (scan E-12)
    paths = []
    for name in ("a.xyz", "b.xyz"):
        p = tmp_path / name
        p.write_text("2\nc\nO 0 0 0\nH 1 0 0\n2\nc\nO 0 0 0\n")
        paths.append(str(p))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("default")
        for p in paths:
            ferro.read(p)
    got = [w.message for w in caught]
    assert [w.category for w in caught] == [ferro.FerroWarning, ferro.FerroWarning]
    assert [(m.source, m.line) for m in got] == [(paths[0], 5), (paths[1], 5)]
    assert str(got[0]) == f"{paths[0]}:5: {got[0].detail}" and "1 of 2 atoms" in got[0].detail
    assert issubclass(ferro.FerroNote, UserWarning) and issubclass(ferro.FerroWarning, UserWarning)


def test_cli_shows_notices_in_the_cli_layout(tmp_path, capsys):
    from ferro import cli
    w = ferro.FerroWarning("x.pdb:57: last MODEL (frame 2) has no ENDMDL; kept")
    w.source, w.line, w.detail = "x.pdb", 57, "last MODEL (frame 2) has no ENDMDL; kept"
    cli._show_warning(w, ferro.FerroWarning, "caller.py", 1)
    n = ferro.FerroNote("sets/a.db: no virial")
    n.source, n.line, n.detail = "sets/a.db", None, "no virial"
    cli._show_warning(n, ferro.FerroNote, "caller.py", 1)
    assert capsys.readouterr().err == ("x.pdb:57: FerroWarning: last MODEL (frame 2) has no ENDMDL; kept\n"
                                       "sets/a.db: FerroNote: no virial\n")


def test_read_uses_the_same_format_table_as_convert(tmp_path):
    # A program output is recognised by its banner, as `ferro convert -i OUTCAR` does
    fixtures = NVT.parent
    t = ferro.read(str(fixtures / "vasp651_sp" / "frame_0011" / "vasprun.xml"))
    assert len(t) == 1 and t.energy(0) is not None
    # An upper-case extension is the same format
    path = tmp_path / "ONE.XYZ"
    ferro.write(ferro.read(str(NVT)).tail(1), str(path))
    assert len(ferro.read(str(path))) == 1


def test_errors_are_ferro_errors(tmp_path):
    assert issubclass(ferro.FerroError, RuntimeError), "except RuntimeError must keep catching them"
    with pytest.raises(ferro.FerroError, match="missing.xyz"):
        ferro.read(str(tmp_path / "missing.xyz"))
    with pytest.raises(ferro.FerroError):
        ferro.read(str(NVT)).frame_indices(stride=0)
