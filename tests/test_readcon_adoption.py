"""readcon-core is pinned past 0.14 and readcon-db holds campaign frames."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_readcon_pin_is_past_0_14_and_frames_use_the_corpus():
    wrap = (ROOT / "subprojects" / "readcon-core.wrap").read_text(encoding="utf-8")
    assert "revision = v0.16.0" in wrap
    meson = (ROOT / "client" / "meson.build").read_text(encoding="utf-8")
    assert "version: '>=0.16.0'" in meson
    project = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert "readcon>=0.15.0" in project
    pixi = (ROOT / "pixi.toml").read_text(encoding="utf-8")
    assert 'readcon = ">=0.15.0"' in pixi
    db = (ROOT / "subprojects" / "readcon-db.wrap").read_text(encoding="utf-8")
    assert "lode-org/readcon-db" in db
    movie = (ROOT / "eon" / "movie.py").read_text(encoding="utf-8")
    assert "frame_texts_for_paths(" in movie
    catalog = (ROOT / "eon" / "process_catalog.py").read_text(encoding="utf-8")
    assert "_saddle_text_with_mode(" in catalog


def test_installed_readcon_core_is_at_least_0_16():
    pkg = shutil.which("pkg-config")
    if pkg is None:
        pytest.skip("pkg-config is not on PATH")
    env = os.environ.copy()
    probe = subprocess.run(
        [pkg, "--exists", "readcon-core"],
        check=False,
        env=env,
    )
    if probe.returncode != 0:
        pytest.skip("readcon-core is not installed")
    version = subprocess.check_output(
        [pkg, "--modversion", "readcon-core"],
        text=True,
        env=env,
    ).strip()
    parts = tuple(int(piece) for piece in version.split(".")[:2])
    assert parts >= (0, 16), version
