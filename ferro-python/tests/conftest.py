"""Fixture locations and the shared NVT trajectory, used across the test modules."""

from pathlib import Path

import pytest

import ferro

FIXTURES = Path(__file__).resolve().parents[2] / "tests"
NVT = FIXTURES / "70Z30P00A_NVT_5.lammpstrj"


@pytest.fixture(scope="module")
def nvt():
    return ferro.read(str(NVT))
