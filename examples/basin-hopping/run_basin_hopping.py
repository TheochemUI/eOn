"""Basin hopping global optimization with dictionary-style configuration.

Equivalent to the config.ini in this directory. Searches for the global
minimum of a Lennard-Jones cluster at 2000 K using Gaussian-distributed
displacements.
"""

from pathlib import Path

from rgpycrumbs.eon.helpers import write_eon_config

settings = {
    "Main": {
        "job": "basin_hopping",
        "temperature": 2000,
    },
    "Communicator": {
        "type": "local",
    },
    "Potential": {
        "potential": "lj",
    },
    "Basin Hopping": {
        "steps": 100,
        "displacement": 0.25,
        "significant_structure": True,
        "displacement_distribution": "gaussian",
    },
    "Optimizer": {
        "opt_method": "cg",
        "converged_force": 0.01,
        "max_iterations": 10000,
    },
}

if __name__ == "__main__":
    write_eon_config(Path("."), settings)
