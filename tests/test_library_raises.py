"""Library helpers raise, and a catalog key round-trips through the bytes it stores."""

from types import SimpleNamespace

import numpy as np
import pytest

from eon.akmc import akmc, kmc_step
from eon.config import ConfigClass
from eon.process_catalog import env_hash, pack_frame_key, unpack_frame_key
from eon.recycling import SB_Recycling
from eon.structure import Structure
from eon.superbasin import Superbasin


def test_missing_akmc_root_raises(tmp_path):
    config = SimpleNamespace(path_root=str(tmp_path / "absent"))
    with pytest.raises(FileNotFoundError, match="root directory does not exist"):
        akmc(config)


def test_energy_above_the_stop_criterion_raises(tmp_path):
    path = tmp_path / "config.ini"
    path.write_text(
        "\n".join(
            [
                "[Main]",
                "job = akmc",
                "temperature = 300",
                "[AKMC]",
                "confidence = 0.9",
                "max_kmc_steps = 1",
                "[Paths]",
                f"main_directory = {tmp_path}",
                f"results = {tmp_path}",
                "[Debug]",
                "use_mean_time = true",
                "stop_criterion = -1",
                "",
            ]
        )
    )
    cfg = ConfigClass()
    cfg.init(str(path))

    class _State:
        def __init__(self, number, procs, energy=0.0):
            self.number = number
            self.procs = procs
            self.energy = energy

        def get_process(self, pid):
            return self.procs[pid]

        def get_ratetable(self):
            return [(pid, proc["rate"], proc["prefactor"]) for pid, proc in self.procs.items()]

        def get_confidence(self, superbasin=None):
            return 1.0

        def get_energy(self):
            return self.energy

    class _States:
        def __init__(self, mapping):
            self.mapping = mapping

        def get_state(self, number):
            return self.mapping[int(number)]

        def get_product_state(self, reactant, proc_id):
            product = self.mapping[int(reactant)].procs[proc_id]["product"]
            return self.mapping[int(product)]

    state0 = _State(
        0,
        {3: {"rate": 1.0e8, "product": 1, "barrier": 0.2, "prefactor": 1.0e12}},
        energy=0.0,
    )
    state1 = _State(1, {}, energy=0.0)
    with pytest.raises(RuntimeError, match="exceeds the stop criterion"):
        kmc_step(state0, _States({0: state0, 1: state1}), 0.0, 0.025, None, config=cfg)


def test_unknown_recycling_target_names_the_state(tmp_path):
    current = SimpleNamespace(number=9)
    recycler = SB_Recycling.__new__(SB_Recycling)
    recycler.in_progress = True
    recycler.current_state = current
    recycler.sb_state_nums = [[1, 2]]
    recycler.path = tmp_path
    with pytest.raises(ValueError, match="state 9 is not a recycling target"):
        recycler.make_suggestion()


def test_superbasin_without_an_exit_raises():
    def _state(number, product):
        return SimpleNamespace(
            number=number,
            get_process_table=lambda product=product: {
                1: {"product": product, "rate": 1.0}
            },
        )

    basin = Superbasin.__new__(Superbasin)
    basin.id = 4
    basin.config = SimpleNamespace(amsel_discover_decide=False, sb_amsel_discover_decide=False)
    basin.state_numbers = [0, 1]
    basin.state_dict = {0: _state(0, 1), 1: _state(1, 0)}
    with pytest.raises(ValueError, match="superbasin 4 has no exit process"):
        basin.step(basin.state_dict[0], lambda number, pid: None)


def test_frame_key_round_trip():
    blob = pack_frame_key(3, 9)
    assert unpack_frame_key(blob) == (3, 9)
    assert len(blob) == 12


def test_env_hash_ignores_a_tenth_of_a_milliangstrom():
    atoms = Structure(1)
    atoms.names = ["Pt"]
    atoms.r[0] = [1.0, 0.0, 0.0]
    nudged = atoms.copy()
    nudged.r[0, 0] += 1.0e-6
    assert env_hash(atoms) == env_hash(nudged)
    hopped = atoms.copy()
    hopped.r[0, 0] += 0.01
    assert env_hash(atoms) != env_hash(hopped)
    assert len(env_hash(atoms)) == 16
