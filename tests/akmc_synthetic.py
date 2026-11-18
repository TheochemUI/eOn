"""Deterministic Cu vacancy plus <100> split-interstitial AKMC fixture.

A lattice walker stands in for the eOn client. Given a reactant it lists the
interstitial and vacancy hops, builds index-preserving saddle and product
structures and the results dictionary the server reads from a finished
search. Everything above that is production code: ``AKMCStateList``,
``AKMCState``, the superbasin scheme, the amsel gate and ``kmc_step``.

The hop barriers are 0.098 eV for the dumbbell and 0.689 eV for the vacancy.
Every state shares one reactant energy, so every state is energetically
close to every product and the registry has to settle each product by
geometry. Atom indices follow the client: an atom keeps its index through a
hop, so a state is the same state only when the indexed positions agree.

Run as a script it times the driver, with the client stand-in excluded, on a cell of
``--cells`` cubed fcc cells::

    python tests/akmc_synthetic.py --cells 9 --cycles 12 --root DIR
"""

from __future__ import annotations

import argparse
import itertools
import os
import time
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

import numpy as np

A0 = 3.615
HALF = A0 / 2.0
DUMBBELL_OFFSET = 0.6
E_REACTANT = -10318.9297636
SIA_BARRIER = 0.098
VACANCY_BARRIER = 0.689
BARRIER_STEP = 1.0e-4
SADDLE_NOISE = 1.0e-3
BAD_FRACTION = 0.1

_NN = np.array(
    [
        v
        for v in itertools.product((-1, 0, 1), repeat=3)
        if sum(abs(c) for c in v) == 2
    ],
    dtype=int,
)


def build_cell(cells: int):
    """Perfect fcc Cu with one vacancy at the centre and a [001] dumbbell.

    The dumbbell sits on the origin site. Its two atoms are the last two
    atoms of the structure.
    """
    from eon.structure import Structure

    length = cells * A0
    sites = []
    for i in range(2 * cells):
        for j in range(2 * cells):
            for k in range(2 * cells):
                if (i + j + k) % 2 == 0:
                    sites.append((i, j, k))
    sites = np.array(sites, dtype=int)
    vacancy = (cells, cells, cells)
    if (sum(vacancy) % 2) != 0:
        vacancy = (cells, cells, cells + 1)
    origin = (0, 0, 0)
    keep = [
        tuple(s) for s in sites if tuple(s) != vacancy and tuple(s) != origin
    ]
    n = len(keep) + 2
    p = Structure(n)
    p.box = np.diag([length] * 3)
    p.names = ["Cu"] * n
    p.mass[:] = 63.546
    for index, site in enumerate(keep):
        p.r[index] = np.array(site, dtype=float) * HALF
    p.r[n - 2] = np.array([0.0, 0.0, DUMBBELL_OFFSET])
    p.r[n - 1] = np.array([0.0, 0.0, -DUMBBELL_OFFSET])
    p.r %= length
    return p


@dataclass
class Hop:
    kind: str
    barrier: float
    product: object
    saddle: object


class Lattice:
    def __init__(self, cells: int):
        self.cells = cells
        self.length = cells * A0
        self.nsite = 2 * cells

    def _site_of(self, r):
        idx = np.rint(np.asarray(r) / HALF).astype(int) % self.nsite
        return idx

    def _position(self, site):
        return (np.asarray(site, dtype=float) % self.nsite) * HALF

    def _wrap(self, d):
        return d - self.length * np.rint(d / self.length)

    def classify(self, structure):
        sites = self._site_of(structure.r)
        keys = sites[:, 0] * self.nsite**2 + sites[:, 1] * self.nsite + sites[:, 2]
        order = np.argsort(keys, kind="stable")
        occupancy = {}
        for index in order:
            occupancy.setdefault(int(keys[index]), []).append(int(index))
        return sites, keys, occupancy

    def _key(self, site):
        s = np.asarray(site, dtype=int) % self.nsite
        return int(s[0] * self.nsite**2 + s[1] * self.nsite + s[2])

    def hops(self, structure) -> list:
        sites, keys, occupancy = self.classify(structure)
        pair_key = next(k for k, v in occupancy.items() if len(v) == 2)
        d1, d2 = occupancy[pair_key]
        centre = sites[d1]
        axis = int(np.argmax(np.abs(self._wrap(structure.r[d1] - structure.r[d2]))))
        vacancy_site = None
        for site in self._all_sites():
            if self._key(site) not in occupancy:
                vacancy_site = site
                break
        out = []
        direction = 0
        for v in _NN:
            if v[axis] == 0:
                continue
            target = (centre + v) % self.nsite
            tkey = self._key(target)
            holder = occupancy.get(tkey)
            if holder is None or len(holder) != 1:
                continue
            others = [a for a in range(3) if v[a] != 0 and a != axis]
            new_axis = others[0]
            m = holder[0]
            product = structure.copy()
            offset = np.zeros(3)
            offset[new_axis] = DUMBBELL_OFFSET
            product.r[d1] = self._position(centre)
            product.r[d2] = self._position(target) + offset
            product.r[m] = self._position(target) - offset
            product.r %= self.length
            out.append(
                self._finish(
                    "sia",
                    SIA_BARRIER + BARRIER_STEP * direction,
                    structure,
                    product,
                    (d1, d2, m),
                )
            )
            direction += 1
        direction = 0
        for v in _NN:
            target = (np.asarray(vacancy_site) + v) % self.nsite
            holder = occupancy.get(self._key(target))
            if holder is None or len(holder) != 1:
                continue
            m = holder[0]
            product = structure.copy()
            product.r[m] = self._position(vacancy_site)
            product.r %= self.length
            out.append(
                self._finish(
                    "vacancy",
                    VACANCY_BARRIER + BARRIER_STEP * direction,
                    structure,
                    product,
                    (m,),
                )
            )
            direction += 1
        return out

    def _all_sites(self):
        for i in range(self.nsite):
            for j in range(self.nsite):
                for k in range(self.nsite):
                    if (i + j + k) % 2 == 0:
                        yield np.array([i, j, k], dtype=int)

    def _finish(self, kind, barrier, reactant, product, moved):
        saddle = reactant.copy()
        for atom in moved:
            step = self._wrap(product.r[atom] - reactant.r[atom])
            saddle.r[atom] = reactant.r[atom] + 0.5 * step
        saddle.r %= self.length
        return Hop(kind, barrier, product, saddle)


def _con_text(structure) -> StringIO:
    from eon import fileio as io

    buf = StringIO()
    io.savecon(buf, structure)
    buf.seek(0)
    return buf


def _mode_text(mode) -> StringIO:
    from eon import fileio as io

    buf = StringIO()
    io.save_mode(buf, mode)
    buf.seek(0)
    return buf


class SyntheticClient:
    """Stands in for a finished batch of eonclient process searches."""

    def __init__(self, cells: int, seed: int):
        self.lattice = Lattice(cells)
        self.seed = int(seed)
        self.counter = 0
        self.seconds = 0.0

    def results(self, state, count: int) -> list:
        """Finished searches for *state*. The time spent building them is
        kept apart in :attr:`seconds` so the driver can be timed alone."""
        start = time.perf_counter()
        reactant = state.get_reactant()
        hops = self.lattice.hops(reactant)
        out = []
        for _ in range(count):
            rng = np.random.RandomState(self.seed * 1000003 + self.counter)
            wuid = self.counter
            self.counter += 1
            out.append(self._one(state, reactant, hops, rng, wuid))
        self.seconds += time.perf_counter() - start
        return out

    def _one(self, state, reactant, hops, rng, wuid):
        bad = rng.random_sample() < BAD_FRACTION
        hop = hops[int(rng.randint(len(hops)))]
        saddle = hop.saddle.copy()
        moved = np.flatnonzero(
            np.linalg.norm(hop.product.r - reactant.r, axis=1) > 1.0e-9
        )
        saddle.r[moved] += rng.normal(scale=SADDLE_NOISE, size=(len(moved), 3))
        mode = np.zeros_like(reactant.r)
        mode[moved] = rng.normal(size=(len(moved), 3))
        mode /= np.linalg.norm(mode)
        results = {
            "termination_reason": 5 if bad else 0,
            "potential_energy_saddle": E_REACTANT + hop.barrier,
            "potential_energy_reactant": E_REACTANT,
            "potential_energy_product": E_REACTANT,
            "barrier_reactant_to_product": hop.barrier,
            "barrier_product_to_reactant": hop.barrier,
            "displacement_saddle_distance": float(
                np.linalg.norm(saddle.r - reactant.r)
            ),
            "force_calls_saddle": 100,
            "force_calls_minimization": 50,
            "force_calls_prefactors": 0,
            "prefactor_reactant_to_product": 1.0e13,
            "prefactor_product_to_reactant": 1.0e13,
        }
        from eon import fileio as io

        results_text = StringIO()
        io.save_results_dat(results_text, results)
        return {
            "name": "%d_%d" % (state.number, wuid),
            "wuid": wuid,
            "type": "random",
            "results": results,
            "results.dat": results_text,
            "reactant.con": _con_text(reactant),
            "saddle.con": _con_text(saddle),
            "product.con": _con_text(hop.product),
            "mode.dat": _mode_text(mode),
        }


def write_config(
    root: Path, cells: int, extra: str = "", transitions: int = 5
) -> Path:
    text = f"""[Main]
job = akmc
temperature = 300
random_seed = 20260101

[AKMC]
confidence = 0.9
confidence_scheme = new
thermally_accessible_window = 20.0
max_kmc_steps = 1

[Paths]
main_directory = {root}
results = {root}

[Structure Comparison]
distance_difference = 0.1
energy_difference = 0.05
neighbor_cutoff = 3.0

[Debug]
use_mean_time = false

[Coarse Graining]
use_mcamc = true
superbasin_scheme = transition_counting
number_of_transitions = {transitions}

[amsel]
discover_decide = true
"""
    path = root / "config.ini"
    path.write_text(text + extra)
    return path


def setup_run(root: Path, cells: int, extra: str = "", transitions: int = 5):
    """Write pos.con and config.ini under *root* and return the config."""
    from eon import fileio as io
    from eon.config import ConfigClass

    root.mkdir(parents=True, exist_ok=True)
    io.savecon(str(root / "pos.con"), build_cell(cells))
    path = write_config(root, cells, extra, transitions)
    os.chdir(root)
    cfg = ConfigClass()
    cfg.init(str(path))
    return cfg


def cycle(config, client: SyntheticClient, searches: int, steps: int = 0):
    """One akmc() pass with the client's results in place of the explorer.

    Same order as :func:`eon.akmc.akmc`: metadata, state list, superbasin
    scheme, register the batch, take a KMC step, write the metadata back.
    """
    from configparser import ConfigParser

    from eon import akmc
    from eon import fileio as io

    kT = config.main_temperature / 11604.5
    start, t, prev, first_run, prev_temperature = akmc.get_akmc_metadata(config)
    states = akmc.get_statelist(kT, config=config)
    current = states.get_state(start)
    if prev == -1:
        previous = current
    else:
        previous = states.get_state(prev)
    superbasining = akmc.get_superbasin_scheme(states, config)
    sb = superbasining.get_containing_superbasin(current)
    explore = sb.get_lowest_confidence_state() if sb else current

    registered = 0
    for result in client.results(explore, searches):
        if result["results"]["termination_reason"] == 0:
            explore.add_process(result, sb)
        else:
            explore.register_bad_saddle(result, False, superbasin=sb)
        registered += 1
        confidence = sb.get_confidence() if sb else explore.get_confidence()
        if confidence >= config.akmc_confidence:
            break

    current, previous, t, steps = akmc.kmc_step(
        current, states, t, kT, superbasining, steps, config=config
    )
    parser = ConfigParser()
    if previous.number != current.number:
        prev = previous.number
    akmc.write_akmc_metadata(parser, current.number, t, prev, config.main_temperature)
    io.write_info_txt(config, parser)
    io.save_prng_state(io.prng_state_path(config))
    return current.number, t, registered


def read_dynamics(config) -> list:
    """Dynamics rows as (reactant, process, product, step_time, time, rate)."""
    path = Path(config.path_results) / "dynamics.txt"
    rows = []
    if not path.is_file():
        return rows
    for line in path.read_text().splitlines()[2:]:
        f = line.split()
        if len(f) < 9:
            continue
        rows.append(
            (int(f[1]), int(f[2]), int(f[3]), float(f[4]), float(f[5]), float(f[7]))
        )
    return rows


def snapshot(config) -> dict:
    """The server's record of a run as text: what a trajectory is compared on.

    Dynamics rows, the state table, every state's process table and the
    superbasin files. Nothing in it names a path or a time of day.
    """
    root = Path(config.path_results)
    states = Path(config.path_states)
    out = {
        "dynamics": (root / "dynamics.txt").read_text().splitlines()
        if (root / "dynamics.txt").is_file()
        else [],
        "state_table": (states / "state_table").read_text().splitlines(),
        "processtables": {},
        "superbasins": {},
    }
    basins = Path(config.sb_path)
    for path in sorted(basins.rglob("*")):
        if path.is_file():
            out["superbasins"][str(path.relative_to(basins))] = path.read_text()
    numbers = sorted(int(p.name) for p in states.iterdir() if p.name.isdigit())
    for number in numbers:
        table = states / str(number) / "processtable"
        out["processtables"][str(number)] = table.read_text().splitlines()
    return out


def run(root: Path, cells: int, cycles: int, searches: int, transitions: int = 5):
    cfg = setup_run(root, cells, transitions=transitions)
    client = SyntheticClient(cells, 7)
    per_cycle = []
    for _ in range(cycles):
        t0 = time.perf_counter()
        client_before = client.seconds
        cycle(cfg, client, searches)
        per_cycle.append(
            time.perf_counter() - t0 - (client.seconds - client_before)
        )
    return cfg, per_cycle


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=int, default=9)
    ap.add_argument("--cycles", type=int, default=12)
    ap.add_argument("--searches", type=int, default=16)
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--transitions", type=int, default=5)
    ap.add_argument(
        "--record", type=Path, help="write the run's snapshot to this JSON file"
    )
    args = ap.parse_args()
    cfg, per_cycle = run(
        args.root, args.cells, args.cycles, args.searches, args.transitions
    )
    if args.record is not None:
        import json

        args.record.write_text(
            json.dumps(snapshot(cfg), indent=1, sort_keys=True) + "\n"
        )
    rows = read_dynamics(cfg)
    n = len(per_cycle)
    mean = float(np.mean(per_cycle))
    se = float(np.std(per_cycle, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    print("cycles", n, "kmc_steps", len(rows))
    print("driver_seconds_per_cycle_mean %.4f se %.4f" % (mean, se))
    for index, seconds in enumerate(per_cycle):
        print("cycle %d %.4f" % (index, seconds))


if __name__ == "__main__":
    main()
