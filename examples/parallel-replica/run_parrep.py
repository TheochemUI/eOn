"""Parallel replica dynamics with dictionary-style configuration.

Equivalent to the config.ini in this directory. Runs parallel replica
dynamics on an Al(100) surface with an Andersen thermostat at 500 K.
"""

from pathlib import Path

from rgpycrumbs.eon.helpers import write_eon_config

settings = {
    "Main": {
        "job": "parallel_replica",
        "temperature": 500,
        "random_seed": 1042,
    },
    "Potential": {
        "potential": "eam_al",
    },
    "Communicator": {
        "type": "local",
        "number_of_cpus": 1,
        "num_jobs": 2,
    },
    "Dynamics": {
        "time_step": 1.0,
        "time": 5000.0,
        "thermostat": "andersen",
        "andersen_alpha": 0.2,
        "andersen_collision_period": 10.0,
    },
    "Parallel Replica": {
        "dephase_time": 1000.0,
        "state_check_interval": 3000.0,
        "state_save_interval": 1000.0,
        "post_transition_time": 200.0,
        "stop_after_transition": False,
    },
    "Optimizer": {
        "opt_method": "cg",
        "converged_force": 0.005,
    },
}

if __name__ == "__main__":
    write_eon_config(Path("."), settings)
