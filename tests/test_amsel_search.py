"""[amsel] displacement = direction_cover."""

import sys
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

from eon.amsel_search import AmselSearch
from tests.test_library_bodies import _states


class _Info(dict):
    def get(self, section, key, default=None):
        return dict.get(self, (section, key), default)

    def set(self, section, key, value):
        self[(section, key)] = value


def _amsel(monkeypatch, mask, calls):
    module = ModuleType("amsel")

    def select_active_atoms(positions, cell, **kwargs):
        calls.append(("select", kwargs["rule"]))
        return SimpleNamespace(
            mask=mask, status="localized_defects", cluster_method="decaf"
        )

    def plan_generic_directions(positions, *, active_mask, cell_matrix, seed, budget):
        calls.append(("plan", seed, budget))
        atoms = [i for i, m in enumerate(active_mask) if m]
        return SimpleNamespace(
            proposals=tuple(
                SimpleNamespace(
                    atom=atoms[k % len(atoms)], direction=(0.0, 0.6, 0.8), arm_index=k
                )
                for k in range(budget)
            )
        )

    module.select_active_atoms = select_active_atoms
    module.plan_generic_directions = plan_generic_directions
    monkeypatch.setitem(sys.modules, "amsel", module)


def _search(tmp_path):
    cfg, _all, state, _product, _proc, _reactant = _states(tmp_path)
    cfg.amsel_active_atoms = "decaf_clustering"
    cfg.amsel_direction_seed = 5
    cfg.disp_magnitude = 0.1
    reactant = state.get_reactant()
    holder = SimpleNamespace(number=3, info=_Info())
    return AmselSearch(holder, reactant, cfg), holder, reactant, cfg


def test_the_start_is_on_an_active_atom_along_the_planned_direction(
    tmp_path, monkeypatch
):
    calls = []
    _amsel(monkeypatch, [False, True], calls)
    search, holder, reactant, _cfg = _search(tmp_path)
    displaced, mode = search.make_displacement()
    assert np.allclose(mode[0], 0.0)
    assert np.allclose(mode[1], (0.0, 0.6, 0.8))
    assert np.allclose(displaced.r - reactant.r, 0.1 * mode)
    assert calls[0] == ("select", "decaf_clustering")
    assert calls[1][:2] == ("plan", 8)
    assert holder.info.get("amsel", "directions_used") == 1


def test_a_second_process_resumes_the_plan(tmp_path, monkeypatch):
    calls = []
    _amsel(monkeypatch, [True, True], calls)
    search, holder, reactant, cfg = _search(tmp_path)
    for _ in range(3):
        search.make_displacement()
    again = AmselSearch(holder, reactant, cfg)
    _displaced, mode = again.make_displacement()
    assert holder.info.get("amsel", "directions_used") == 4
    assert np.allclose(mode[1], (0.0, 0.6, 0.8))


def test_the_plan_grows_past_its_first_chunk(tmp_path, monkeypatch):
    calls = []
    _amsel(monkeypatch, [True, True], calls)
    search, holder, _reactant, _cfg = _search(tmp_path)
    holder.info.set("amsel", "directions_used", 40)
    search.make_displacement()
    assert calls[-1] == ("plan", 8, 82)


def test_no_active_atom_is_an_error(tmp_path, monkeypatch):
    _amsel(monkeypatch, [False, False], [])
    search, *_rest = _search(tmp_path)
    with pytest.raises(RuntimeError, match="no free active atom"):
        search.make_displacement()
