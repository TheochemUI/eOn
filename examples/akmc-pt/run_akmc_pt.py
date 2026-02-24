"""AKMC simulation of Pt heptamer with dictionary-style configuration.

Equivalent to the config.ini in this directory. Uses the dimer method
for saddle searches with displacement biased towards the least
coordinated atom.
"""

from pathlib import Path

from rgpycrumbs.eon.helpers import write_eon_config

settings = {
    "Main": {
        "job": "akmc",
        "temperature": 300,
    },
    "Potential": {
        "potential": "morse_pt",
    },
    "Optimizer": {
        "converged_force": 0.001,
        "max_iterations": 1000,
    },
    "AKMC": {
        "confidence": 0.95,
    },
    "Process Search": {
        "minimize_first": True,
    },
    "Communicator": {
        "type": "local",
        "number_of_CPUs": 2,
        "num_jobs": 2,
    },
    "Saddle Search": {
        "displace_least_coordinated_weight": 1.0,
        "displace_radius": 3.3,
        "displace_magnitude": 0.1,
        "min_mode_method": "dimer",
        "max_energy": 10.0,
    },
}

if __name__ == "__main__":
    write_eon_config(Path("."), settings)
