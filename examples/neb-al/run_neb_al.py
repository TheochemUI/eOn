"""Nudged elastic band calculation with dictionary-style configuration.

Equivalent to the config.ini in this directory. Runs an NEB calculation
on Al with 7 intermediate images and LBFGS optimization.
"""

from pathlib import Path

from rgpycrumbs.eon.helpers import write_eon_config

settings = {
    "Main": {
        "job": "nudged_elastic_band",
    },
    "Potential": {
        "potential": "eam_al",
    },
    "Nudged Elastic Band": {
        "images": 7,
        "spring": 5.0,
    },
    "Optimizer": {
        "max_iterations": 1000,
        "opt_method": "lbfgs",
        "max_move": 0.1,
        "converged_force": 0.001,
    },
}

if __name__ == "__main__":
    write_eon_config(Path("."), settings)
