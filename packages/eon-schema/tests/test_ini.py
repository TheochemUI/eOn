"""INI write / hydrate / validate helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from eon_schema.config import (
    MainConfig,
    DimerConfig,
    PotentialConfig,
    XtsciConfig,
    defaults_from_catalog,
    hydrate_ini,
    model_to_ini_section,
    models_to_ini,
    read_ini,
    unknown_ini_keys,
    write_ini,
    write_models_ini,
)


def test_write_and_read_roundtrip(tmp_path: Path):
    sections = {
        "Main": {"job": "saddle_search", "random_seed": 42},
        "Potential": {"potential": "lj"},
    }
    path = write_ini(tmp_path / "config.ini", sections)
    assert path.is_file()
    got = read_ini(path)
    assert got["Main"]["job"] == "saddle_search"
    assert got["Main"]["random_seed"] == "42"
    assert got["Potential"]["potential"] == "lj"


def test_write_ini_into_directory(tmp_path: Path):
    path = write_ini(tmp_path, {"Main": {"job": "point"}})
    assert path.name == "config.ini"
    assert read_ini(tmp_path)["Main"]["job"] == "point"


def test_defaults_from_catalog_main():
    d = defaults_from_catalog("Main")
    assert isinstance(d, dict)
    assert "job" in d
    assert "temperature" in d


def test_hydrate_user_overrides():
    user = {"Main": {"job": "minimization", "temperature": 500.0}}
    h = hydrate_ini(user, base_sections=["Main"])
    assert h["Main"]["job"] == "minimization"
    assert float(h["Main"]["temperature"]) == 500.0
    # defaults still present for unset keys
    assert "random_seed" in h["Main"] or "finite_difference" in h["Main"]


def test_unknown_ini_keys_flags_typos():
    bad = unknown_ini_keys(
        {"Main": {"job": "point", "not_a_real_key": 1}},
        covered_only=True,
    )
    assert "Main.not_a_real_key" in bad


def test_unknown_ini_keys_skips_uncovered_sections():
    # SocketNWChemPot is not an L0 covered section
    bad = unknown_ini_keys(
        {
            "Main": {"job": "point"},
            "SocketNWChemPot": {"unix_socket_path": "/tmp/x"},
        },
        covered_only=True,
    )
    assert not any(x.startswith("SocketNWChemPot") for x in bad)


def test_model_to_ini_section_xtsci_qn():
    opts = model_to_ini_section(
        XtsciConfig(method="lbfgs", qn_step="newton", precon="pair")
    )
    assert opts["method"] == "lbfgs"
    assert opts["qn_step"] == "newton"
    assert opts["precon"] == "pair"
    assert "xtsci_qn_step" not in opts


def test_model_to_ini_section_dimer_aliases():
    d = DimerConfig(dimer_improved=True, dimer_remove_rotation=False)
    opts = model_to_ini_section(d)
    assert opts["improved"] == "true"
    assert opts["remove_rotation"] == "false"
    assert "dimer_improved" not in opts


def test_write_models_ini(tmp_path: Path):
    path = write_models_ini(
        tmp_path / "job.ini",
        MainConfig(job="saddle_search", random_seed=7),
        PotentialConfig(potential="lj"),
        DimerConfig(),
        extra={"Debug": {"write_movies": True}},
        validate=True,
    )
    got = read_ini(path)
    assert got["Main"]["job"] == "saddle_search"
    assert got["Potential"]["potential"] == "lj"
    assert got["Dimer"]["improved"] == "true"
    assert got["Debug"]["write_movies"] == "true"


def test_models_to_ini_combined():
    secs = models_to_ini(MainConfig(job="point"), PotentialConfig(potential="emt"))
    assert set(secs) >= {"Main", "Potential"}
