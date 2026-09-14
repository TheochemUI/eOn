"""eon-schema — shared eOn parameter SSoT + job config + API models.

Layers
------
* **L0** — Cap'n Proto field graph (vendored copy of monorepo ``schema/``):
  params (``eon_schema.ssot``) and job request/result (``eon_schema.jobs``).
* **L1** — full job config pydantic models (``eon_schema.config``); also
  re-exported as ``eon.schema`` from eon-akmc for Sphinx / legacy imports.
* **L2** — in-process API composition (``DimerSpec``, ``NebSpec``).

Import::

    import eon_schema
    from eon_schema.ssot import capnp_path, load_catalog
    from eon_schema.jobs import job_result_capnp_path, results_dat_to_dict
    from eon_schema.config import MainConfig, Config  # needs pydantic
    from eon_schema.api import DimerSpec, NebSpec
"""

from __future__ import annotations

__version__ = "0.2.3"

from eon_schema.ssot import capnp_path, catalog_path, load_catalog
from eon_schema.jobs import (
    dict_to_results_dat,
    job_result_capnp_path,
    job_result_scalars_from_results_dat,
    results_dat_to_dict,
)

__all__ = [
    "__version__",
    "capnp_path",
    "catalog_path",
    "load_catalog",
    "job_result_capnp_path",
    "results_dat_to_dict",
    "dict_to_results_dat",
    "job_result_scalars_from_results_dat",
]
