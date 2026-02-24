"""AKMC simulation of Cu vacancy with dictionary-style configuration.

Equivalent to the config.ini in this directory. Demonstrates the use of
a displacement script (ptmdisp.py) for targeted saddle searches around
the vacancy site, identified via polyhedral template matching.
"""

from pathlib import Path

from rgpycrumbs.eon.helpers import write_eon_config

settings = {
    "Main": {
        "job": "akmc",
        "temperature": 300,
        "random_seed": 757783492,
        "finite_difference": 0.01,
    },
    "Communicator": {
        "type": "local",
        "num_jobs": 4,
        "number_of_cpus": 4,
    },
    "AKMC": {
        "confidence": 0.005,
        "confidence_scheme": "new",
        "thermally_accessible_window": 100.0,
    },
    "Optimizer": {
        "max_iterations": 1000,
        "opt_method": "lbfgs",
        "converged_force": 0.01,
        "lbfgs_memory": 25,
        "lbfgs_inverse_curvature": 0.01,
        "lbfgs_auto_scale": True,
        "convergence_metric": "norm",
        "max_move": 0.05,
    },
    "Dimer": {
        "opt_method": "cg",
        "improved": True,
        "rotations_max": 20,
        "converged_angle": 1.0,
        "remove_rotation": False,
    },
    "Potential": {
        "potential": "emt",
    },
    "Process Search": {
        "minimize_first": True,
    },
    "Saddle Search": {
        "method": "min_mode",
        "min_mode_method": "dimer",
        "max_energy": 300,
        "displace_listed_atom_weight": 1.0,
        "displace_radius": 3.0,
        "displace_magnitude": 0.01,
        "displace_atom_kmc_state_script": "ptmdisp.py",
        "displace_all_listed": True,
        "remove_rotation": False,
    },
    "Structure Comparison": {
        "distance_difference": 0.1,
        "energy_difference": 0.08,
        "neighbor_cutoff": 3.0,
    },
    "Debug": {
        "write_movies": True,
    },
    "Prefactor": {
        "default_value": 1e13,
    },
}

if __name__ == "__main__":
    write_eon_config(Path("."), settings)
