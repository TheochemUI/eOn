"""Advanced NEB with climbing image and energy-weighted springs.

Demonstrates the full set of NEB options available via dictionary
configuration. This extends the basic NEB example with:

- Climbing image method for accurate saddle point location
- Energy-weighted springs for better resolution near the barrier
- IDPP initialization for a smoother initial path
- Off-path climbing image (MMF) for dimer-like refinement at the saddle

See ``eon.schema.NudgedElasticBandConfig`` for all available fields.
"""

from pathlib import Path

from rgpycrumbs.eon.helpers import write_eon_config

N_IMAGES = 10

settings = {
    "Main": {
        "job": "nudged_elastic_band",
    },
    "Potential": {
        "potential": "eam_al",
    },
    "Nudged Elastic Band": {
        "images": N_IMAGES,
        "spring": 5.0,
        "climbing_image_method": True,
        "energy_weighted": True,
        "ew_ksp_min": 0.972,
        "ew_ksp_max": 9.72,
        "initializer": "idpp",
        # Off-path climbing image with dimer-like refinement
        "ci_mmf": True,
        "ci_mmf_after": 0.5,
        "ci_mmf_nsteps": 1000,
    },
    "Optimizer": {
        "max_iterations": 1000,
        "opt_method": "lbfgs",
        "max_move": 0.1,
        "converged_force": 0.01,
    },
}

if __name__ == "__main__":
    write_eon_config(Path("."), settings)
