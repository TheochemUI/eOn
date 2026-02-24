"""AKMC simulation of Al with dictionary-style configuration.

Equivalent to the config.ini in this directory, but uses
``write_eon_config`` from rgpycrumbs to generate the INI file
programmatically. All option names match the fields documented in
``eon.schema`` (the Pydantic models are the single source of truth).
"""

from pathlib import Path

from rgpycrumbs.eon.helpers import write_eon_config

settings = {
    "Main": {
        "job": "akmc",
        "temperature": 300,
        "random_seed": 42,
    },
    "Potential": {
        "potential": "eam_al",
    },
    "Optimizer": {
        "opt_method": "lbfgs",
        "converged_force": 1e-4,
        "max_iterations": 1000,
    },
    "AKMC": {
        "confidence": 0.95,
    },
    "Process Search": {
        "minimize_first": True,
    },
    "Prefactor": {
        "default_value": 1e12,
    },
    "Communicator": {
        "type": "local",
        "number_of_cpus": 2,
        "num_jobs": 8,
    },
    "Lanczos": {
        "tolerance": 0.05,
    },
    "Saddle Search": {
        "min_mode_method": "lanczos",
        "displace_radius": 5.0,
        "displace_magnitude": 0.2,
        "max_energy": 10.0,
        "displace_listed_atom_weight": 1.0,
        "displace_atom_list": -1,
    },
}

if __name__ == "__main__":
    write_eon_config(Path("."), settings)
