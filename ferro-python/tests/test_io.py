"""read / write 的扩展名分派与 CLI 的 io_dispatch 一致。"""

import warnings

import pytest

import ferro

from conftest import NVT


@pytest.mark.parametrize("ext", ["data", "lmp", "lammps"])
def test_lammps_extension_is_data_file(tmp_path, ext):
    t = ferro.read(str(NVT)).tail(1)
    path = tmp_path / f"one.{ext}"
    ferro.write(t, str(path), format=ext)
    assert "atoms" in path.read_text(), f"format={ext} should write LAMMPS data"
    with pytest.raises(RuntimeError, match="atom_style"):
        ferro.read(str(path))
    back = ferro.read(str(path), atom_style="full")
    assert back.symbols(0) == t.symbols(0), f".{ext} should read back as data, elements intact"
    with pytest.raises(RuntimeError, match="LAMMPS dump input only"):
        ferro.read(str(path), atom_style="full", type_map={1: "P"})


@pytest.mark.parametrize("name", ["one.vasp", "one.poscar", "one.pos", "POSCAR_1", "CONTCAR"])
def test_vasp_names_are_poscar(tmp_path, name):
    t = ferro.read(str(NVT)).tail(1)
    path = tmp_path / name
    ferro.write(t, str(path), format="vasp")
    assert path.read_text() == ferro.poscar_string(t), f"{name} should hold a POSCAR"
    assert ferro.read(str(path)).symbols(0) == t.symbols(0), f"{name} should read back as POSCAR"


def test_format_is_explicit_and_never_read_off_the_path(tmp_path):
    t = ferro.read(str(NVT)).tail(1)
    # A POSCAR-prefixed name with another extension is written in the format asked for
    path = tmp_path / "POSCAR_x.cif"
    ferro.write(t, str(path), format="cif")
    assert path.read_text().count("data_") == 1, "format='cif' must write CIF"
    # A leading dot and upper case are the same format name
    ferro.write(t, str(tmp_path / "a.out"), format=".EXTXYZ")
    assert "Properties=" in (tmp_path / "a.out").read_text().splitlines()[1]
    with pytest.raises(ferro.FerroError, match="Unsupported output format"):
        ferro.write(t, str(tmp_path / "b.xyz"), format="nosuchfmt")
    with pytest.raises(TypeError):
        ferro.write(t, str(tmp_path / "c.xyz"))


@pytest.mark.parametrize("name, fmt, made", [
    ("conf.cif", "cif", ["conf_0000.cif", "conf_0002.cif"]),
    ("conf.vasp", "vasp", ["conf_0000.vasp", "conf_0002.vasp"]),
    ("POSCAR", "vasp", ["POSCAR_0000", "POSCAR_0002"]),
    ("conf.data", "data", ["conf_0000.data", "conf_0002.data"]),
])
def test_single_structure_formats_write_one_file_per_frame(tmp_path, name, fmt, made):
    t = ferro.read(str(NVT)).tail(3)
    ferro.write(t, str(tmp_path / name), format=fmt)
    assert not (tmp_path / name).exists(), "the unindexed name must not be written"
    files = sorted(p.name for p in tmp_path.iterdir())
    assert len(files) == 3 and files[0] == made[0] and files[-1] == made[1], files
    if name.endswith(".cif"):
        assert (tmp_path / made[0]).read_text().count("data_") == 1


def test_trajectory_formats_write_one_file(tmp_path):
    t = ferro.read(str(NVT)).tail(3)
    ferro.write(t, str(tmp_path / "all.extxyz"), format="extxyz")
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
    ferro.write(ferro.read(str(NVT)).tail(1), str(path), format="xyz")
    assert len(ferro.read(str(path))) == 1


def test_format_reads_a_file_whatever_its_name(tmp_path):
    t = ferro.read(str(NVT)).tail(1)
    path = tmp_path / "one.txt"
    ferro.write(t, str(path), format="vasp")
    with pytest.raises(ferro.FerroError, match="Unsupported input format"):
        ferro.read(str(path))
    assert ferro.read(str(path), format="vasp/poscar").symbols(0) == t.symbols(0)
    with pytest.raises(ferro.FerroError, match="vasp/poscar"):
        ferro.read(str(path), format="vasp")
    # The checks follow the given format, not the name
    with pytest.raises(ferro.FerroError, match="atom_style"):
        ferro.read(str(path), format="lammps/data")


def test_errors_are_ferro_errors(tmp_path):
    assert issubclass(ferro.FerroError, RuntimeError), "except RuntimeError must keep catching them"
    with pytest.raises(ferro.FerroError, match="missing.xyz"):
        ferro.read(str(tmp_path / "missing.xyz"))
    with pytest.raises(ferro.FerroError):
        ferro.read(str(NVT)).frame_indices(stride=0)


def test_partial_charge_column_needs_drop_partial(tmp_path):
    (tmp_path / "q.data").write_text(
        "x\n\n2 atoms\n2 atom types\n\n0 5 xlo xhi\n0 5 ylo yhi\n0 5 zlo zhi\n\n"
        "Masses\n\n1 65.38 # Zn\n2 15.999 # O\n\nAtoms # charge\n\n1 1 2.0 0 0 0\n2 2 -2.0 2.5 2.5 2.5\n")
    (tmp_path / "POSCAR").write_text("ZnO\n1.0\n5 0 0\n0 5 0\n0 0 5\nZn O\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5\n")
    merged = ferro.merge(ferro.read(str(tmp_path / "q.data"), atom_style="charge"),
                         ferro.read(str(tmp_path / "POSCAR")), "z", 2.0)
    for name in ("m.extxyz", "m.data"):
        with pytest.raises(ferro.FerroError, match="2 of 4 atoms have a charge"):
            ferro.write(merged, str(tmp_path / name), format=name.split(".")[1])
        assert not (tmp_path / name).exists(), "refused before the file is created"
    ferro.write(merged, str(tmp_path / "m.extxyz"), format="extxyz", drop_partial=True)
    header = (tmp_path / "m.extxyz").read_text().splitlines()[1]
    assert "Properties=" in header and "charge" not in header, "drop_partial leaves the column out"
    ferro.write(merged, str(tmp_path / "m.data"), format="data", drop_partial=True)
    atoms = (tmp_path / "m.data").read_text().split("Atoms # full\n\n")[1].splitlines()[:4]
    assert [line.split()[3] for line in atoms] == ["0.000000"] * 4, \
        "LAMMPS data must have a charge column: every charge written as 0"


ZNO_POSCAR = "ZnO\n1.0\n5 0 0\n0 5 0\n0 0 5\nZn O\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5\n"


def test_supercell_smoke(tmp_path):
    # The binding has no other caller in the package: a signature change would go unnoticed
    (tmp_path / "POSCAR").write_text(ZNO_POSCAR)
    sc = ferro.supercell(ferro.read(str(tmp_path / "POSCAR")), 2, 2, 2)
    assert sc.n_atoms() == 16
    assert sum(sc.cell(0), []) == pytest.approx([10, 0, 0, 0, 10, 0, 0, 0, 10])


def test_supercell_over_the_trajectory_total_is_refused(tmp_path):
    # Each frame stays under the 1e8 cap (2 atoms x 3e7), the two frames together do not
    frame = 'Lattice="5 0 0 0 5 0 0 0 5" Properties=species:S:1:pos:R:3 pbc="T T T"\nZn 0 0 0\nO 2.5 2.5 2.5\n'
    (tmp_path / "two.extxyz").write_text(("2\n" + frame) * 2)
    traj = ferro.read(str(tmp_path / "two.extxyz"))
    with pytest.raises(ferro.FerroError, match="120000000 atoms in total"):
        ferro.supercell(traj, 6000, 5000, 1)


def test_add_vacuum_layer_smoke(tmp_path):
    (tmp_path / "POSCAR").write_text(ZNO_POSCAR)
    out = ferro.add_vacuum_layer(ferro.read(str(tmp_path / "POSCAR")), "z", 7.5)
    assert out.n_atoms() == 2
    assert sum(out.cell(0), []) == pytest.approx([5, 0, 0, 0, 5, 0, 0, 0, 12.5])


def test_build_box_is_reproducible():
    box = ferro.build_box([("H2O", 30), ("ZnO", 5)], 1.2, seed=7)
    assert box.n_atoms() == 30 * 3 + 5 * 2
    assert box.positions(0) == ferro.build_box([("H2O", 30), ("ZnO", 5)], 1.2, seed=7).positions(0)
    with pytest.raises(ferro.FerroError, match="density"):
        ferro.build_box([("H2O", 30)], -1.0)


def test_gr_defaults_are_those_of_the_cli():
    # ferro traj gr: --r-min 0.001, --dr 0.002; the bin centres start at r_min + dr/2.
    # r_max - r_min is a whole number of bins, so the width is not narrowed (md/util.rs::bins)
    t = ferro.read(str(NVT)).subset([0])
    g = ferro.gr_pair(t, "P", "O", r_max=3.001)
    assert g["r"][0] == pytest.approx(0.002) and g["r"][1] - g["r"][0] == pytest.approx(0.002)
    assert g == ferro.gr_pair(t, "P", "O", r_max=3.001, dr=0.002, r_min=0.001)
    assert ferro.gr_all(t, r_max=3.001)["r"] == g["r"]
