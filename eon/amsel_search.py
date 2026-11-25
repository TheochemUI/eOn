"""Saddle-search starts from amsel for ``[amsel] displacement = direction_cover``.

``amsel.select_active_atoms`` picks the atoms and
``amsel.plan_generic_directions`` orders the directions on them. The plan is a
function of the state and the seed, so a state resumes it from the number of
starts it has already used, which is kept in the state's ``info`` file.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

logger = logging.getLogger("amsel")

_SECTION = "amsel"
_KEY = "directions_used"
_CHUNK = 32


class AmselSearch:
    """The ordered (atom, direction) starts of one state."""

    def __init__(self, state: Any, reactant: Any, config: Any):
        self.state = state
        self.reactant = reactant
        self.config = config
        self._plan: tuple = ()
        self._mask: list[bool] | None = None

    def _active_mask(self) -> list[bool]:
        if self._mask is not None:
            return self._mask
        import amsel

        rule = self.config.amsel_active_atoms
        box = np.asarray(self.reactant.box, dtype=float)
        names: dict[str, int] = {}
        species = [names.setdefault(str(n), len(names)) for n in self.reactant.names]
        selection = amsel.select_active_atoms(
            np.asarray(self.reactant.r, dtype=float).tolist(),
            np.linalg.norm(box, axis=1).tolist(),
            rule=rule,
            species=species,
            cell_matrix=box.tolist(),
        )
        free = self.reactant.atom_is_free()
        mask = [bool(m) and bool(f) for m, f in zip(selection.mask, free)]
        if not any(mask):
            raise RuntimeError(
                f"amsel {rule} found no free active atom in state "
                f"{self.state.number} (status {selection.status})"
            )
        logger.info(
            "amsel %s: %d active atoms in state %s (%s, %s)",
            rule,
            sum(mask),
            self.state.number,
            selection.status,
            selection.cluster_method,
        )
        self._mask = mask
        return mask

    def _start(self, index: int) -> Any:
        if index >= len(self._plan):
            import amsel

            self._plan = amsel.plan_generic_directions(
                np.asarray(self.reactant.r, dtype=float).tolist(),
                active_mask=self._active_mask(),
                cell_matrix=np.asarray(self.reactant.box, dtype=float).tolist(),
                seed=int(self.config.amsel_direction_seed) + int(self.state.number),
                budget=max(_CHUNK, 2 * (index + 1)),
            ).proposals
        return self._plan[index]

    def make_displacement(self) -> tuple[Any, np.ndarray]:
        """The displaced structure and the unit mode of the next start."""
        index = int(self.state.info.get(_SECTION, _KEY, 0))
        start = self._start(index)
        mode = np.zeros(np.asarray(self.reactant.r).shape)
        mode[start.atom] = start.direction
        mode /= np.linalg.norm(mode)
        displaced = self.reactant.copy()
        displaced.r += float(self.config.disp_magnitude) * mode
        self.state.info.set(_SECTION, _KEY, index + 1)
        return displaced, mode
