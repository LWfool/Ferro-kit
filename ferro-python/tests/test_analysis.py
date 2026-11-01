"""The analysis bindings: ferro.msd and ferro.interplanar_spacings (review 2026-10-30 A6-L7).

Both had no test of their own: msd's required dt and frame-spacing check were tested on the
CLI side only, and interplanar_spacings only through the k-point mesh of ferro.jobs.
"""

import math

import pytest

import ferro

from conftest import NVT


def test_msd_needs_dt():
    t = ferro.read(str(NVT))
    with pytest.raises(TypeError):
        ferro.msd(t)  # dt is required, as --dt is: a guessed default rescales the time axis


def test_msd_refuses_unevenly_spaced_frames():
    t = ferro.read(str(NVT)).subset([0, 1, 3])
    with pytest.raises(ferro.FerroError, match="not evenly spaced"):
        ferro.msd(t, dt=1.0)


def test_msd_keys_and_components():
    # Frames 0-2 only: the fixture's step gap changes between frames 2 and 3
    d = ferro.msd(ferro.read(str(NVT)).subset([0, 1, 2]), dt=1.0)
    assert set(d) == {"time", "msd", "msd_x", "msd_y", "msd_z"}
    assert len(d["time"]) == len(d["msd"]) and d["time"][0] == 0.0
    for tot, x, y, z in zip(d["msd"], d["msd_x"], d["msd_y"], d["msd_z"]):
        assert math.isclose(tot, x + y + z, rel_tol=1e-12, abs_tol=1e-12)


def test_interplanar_spacings_on_a_tilted_cell():
    # a along x, b at 60 degrees in the plane: the spacing between the planes normal to the
    # reciprocal vector of a is a·sin(60°), shorter than the edge
    a = 4.0
    cell = [[a, 0.0, 0.0], [a / 2, a * math.sqrt(3) / 2, 0.0], [0.0, 0.0, 5.0]]
    da, db, dc = ferro.interplanar_spacings(cell)
    assert math.isclose(da, a * math.sqrt(3) / 2, rel_tol=1e-12)
    assert math.isclose(db, a * math.sqrt(3) / 2, rel_tol=1e-12)
    assert math.isclose(dc, 5.0, rel_tol=1e-12)
    with pytest.raises(ferro.FerroError):
        ferro.interplanar_spacings([[1, 0, 0], [2, 0, 0], [0, 0, 1]])
