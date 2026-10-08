"""Confidence stays defined for an uncounted process and a zero prefactor."""

from __future__ import annotations

import textwrap
import warnings

import numpy as np

from eon import fileio as io
from eon.akmcstate import AKMCState
from eon.config import ConfigClass
from eon.structure import Structure


class _List:
    def __init__(self):
        self.kT = 300.0 / 11604.5
        self.thermal_window = 20.0
        self.max_thermal_window = 40.0


def _config(directory, scheme):
    path = directory / "config.ini"
    path.write_text(
        textwrap.dedent(
            f"""
            [Main]
            job = akmc
            temperature = 300

            [AKMC]
            confidence = 0.9
            confidence_scheme = {scheme}

            [Paths]
            main_directory = {directory}
            results = {directory}
            states = {directory / "states"}
            """
        ).lstrip()
    )
    cfg = ConfigClass()
    cfg.init(str(path))
    cfg.akmc_confidence_scheme = scheme
    return cfg


def _state(directory, scheme, prefactor):
    cfg = _config(directory, scheme)
    reactant = Structure(2)
    reactant.names = ["Pt", "Pt"]
    reactant.mass[:] = 195.084
    reactant.box = np.diag([20.0, 20.0, 20.0])
    reactant.r[1] = [2.5, 0.0, 0.0]
    con = directory / "reactant.con"
    io.savecon(str(con), reactant)
    states = directory / "states"
    states.mkdir()
    state = AKMCState(
        str(states / "0"),
        0,
        _List(),
        reactant_path=str(con),
        config=cfg,
    )
    state.append_process_table(
        3, -0.5, prefactor, -1, -1.0, 1.0e12, 0.2, 1.0e6, 0
    )
    return state


def test_sampling_without_repeat_counts_is_zero(tmp_path):
    state = _state(tmp_path, "sampling", 1.0e12)
    assert state.get_confidence() == 0.0


def test_dynamics_confidence_ignores_a_zero_prefactor(tmp_path):
    state = _state(tmp_path, "dynamics", 0.0)
    state.increment_time(1.0e6, 600.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        confidence = state.get_confidence()
    assert np.isfinite(confidence)
    assert confidence == 0.0


def test_old_confidence_rises_with_repeats(tmp_path):
    state = _state(tmp_path, "old", 1.0e12)
    assert state.get_confidence() == 0.0
    state.info.set("MetaData", "repeats", 5)
    confidence = state.get_confidence()
    assert 0.0 < confidence <= 1.0
