"""The wheel's version is ferro-python/Cargo.toml's, kept by hand in step with the workspace's."""

import tomllib
from pathlib import Path

import ferro

ROOT = Path(__file__).resolve().parents[2]


def test_version_matches_the_workspace():
    workspace = tomllib.loads((ROOT / "Cargo.toml").read_text())["workspace"]["package"]["version"]
    assert ferro.__version__ == workspace, (
        f"ferro-python/Cargo.toml says {ferro.__version__}, the workspace Cargo.toml {workspace}: "
        "bump both together"
    )
