"""A recycling saddle must be visible to sampling confidence."""

from __future__ import annotations

import shutil
from io import StringIO
from pathlib import Path

from eon.akmcstate import AKMCState
from eon.akmcstatelist import AKMCStateList


def _blob(text: str) -> StringIO:
    return StringIO(text)


def _result(saddle: str, mode: str) -> dict:
    energies = {
        "potential_energy_reactant": -1.0,
        "potential_energy_saddle": -0.5,
        "potential_energy_product": -0.9,
        "prefactor_reactant_to_product": 1.0e12,
        "prefactor_product_to_reactant": 1.0e12,
        "barrier_reactant_to_product": 0.5,
        "displacement_saddle_distance": 0.2,
        "force_calls_saddle": 1,
        "force_calls_minimization": 1,
        "force_calls_prefactors": 1,
    }
    return {
        "type": "recycling",
        "wuid": 1,
        "results": energies,
        "saddle.con": _blob(saddle),
        "reactant.con": _blob(saddle),
        "product.con": _blob(saddle),
        "mode.dat": _blob(mode),
        "results.dat": _blob("barrier_reactant_to_product 0.5\n"),
        "saddle": True,
        "reactant": True,
        "product": True,
        "mode": True,
    }


def test_recycling_search_is_visible_to_sampling_confidence():
    root = Path(__file__).resolve().parent / "recycling-confidence-run"
    if root.exists():
        shutil.rmtree(root)
    root.mkdir()
    try:
        reactant = root / "reactant.con"
        source = Path(__file__).resolve().parent / "data" / "client" / "mini_1Pt_oneLayer" / "pos.con"
        reactant.write_text(source.read_text())
        config = type("Config", (), {})()
        config.akmc_server_side_process_search = False
        config.akmc_confidence_scheme = "sampling"
        config.akmc_confidence_correction = False
        config.akmc_eq_rate = 0.0
        config.recycling_on = False
        config.disp_moved_only = False
        config.sb_superbasin_confidence = False
        config.kdb_on = False
        config.comp_eps_e = 0.01
        config.comp_eps_r = 0.1
        config.comp_neighbor_cutoff = 3.3
        config.comp_check_rotation = False
        config.comp_use_identical = True
        config.comp_remove_translation = True
        config.path_states = str(root / "states")
        states = AKMCStateList(0.025, 20.0, 20.0, str(reactant), config=config)
        state = states.get_state(0)
        assert isinstance(state, AKMCState)
        state.add_process(_result(reactant.read_text(), "0.0 0.0 1.0\n"))
        confidence = state.get_confidence()
    finally:
        shutil.rmtree(root)
    assert confidence == 0.0
