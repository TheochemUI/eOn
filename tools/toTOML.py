#!/usr/bin/env python3
import configparser
import sys
from typing import Any, Dict
import argparse
from pydantic import ValidationError
from pathlib import Path

import toml

from eon.schema import *  # Import the schema
from eon.migrator.config import ConfigClass  # Import the old config class


# Function to update the default Pydantic model with values from the old config
def update_model_with_config(
    default_model: BaseModel, config: Dict[str, Any]
) -> BaseModel:
    model_dict = default_model.dict()
    model_dict.update({k: v for k, v in config.items() if v is not None})
    return default_model.__class__(**model_dict)


# Function to transform the validated old config class into a dictionary
def transform_old_config_to_dict(config: ConfigClass) -> Dict[str, Any]:
    return {
        "main": {
            "job": config.main_job,
            "temperature": config.main_temperature,
            "checkpoint": config.main_checkpoint,
            "random_seed": config.main_random_seed,
        },
        "structure_comparison": {
            "energy_difference": config.comp_eps_e,
            "distance_difference": config.comp_eps_r,
            "indistinguishable_atoms": config.comp_use_identical,
            "check_rotation": config.comp_check_rotation,
            "brute_neighbors": config.comp_brute_neighbors,
            "neighbor_cutoff": config.comp_neighbor_cutoff,
            "use_covalent": config.comp_use_covalent,
            "covalent_scale": config.comp_covalent_scale,
            "remove_translation": config.comp_remove_translation,
        },
        "akmc": {
            "confidence": config.akmc_confidence,
            "server_side_process_search": config.akmc_server_side_process_search,
            "thermally_accessible_window": config.akmc_thermal_window,
            "thermally_accessible_buffer": config.akmc_max_thermal_window,
            "max_kmc_steps": config.akmc_max_kmc_steps,
            "confidence_scheme": config.akmc_confidence_scheme,
            "confidence_correction": config.akmc_confidence_correction,
        },
        "basin_hopping": {
            "initial_random_structure_probability": config.bh_initial_random_structure_probability,
            "initial_state_pool_size": config.bh_initial_state_pool_size,
        },
        "paths": {
            "main_directory": config.path_root,
            "jobs_out": config.path_jobs_out,
            "jobs_in": config.path_jobs_in,
            "incomplete": config.path_incomplete,
            "states": config.path_states,
            "results": config.path_results,
            "potential_files": config.path_pot,
            "bh_minima": config.path_bh_minima,
        },
        "communicator": {
            "type": config.comm_type,
            "jobs_per_bundle": config.comm_job_bundle_size,
            "num_jobs": config.comm_job_buffer_size,
            "max_jobs": config.comm_job_max_size,
            "client_path": getattr(config, "comm_local_client", None),
            "number_of_CPUs": getattr(config, "comm_local_ncpus", None),
            "script_path": getattr(config, "comm_script_path", None),
            "name_prefix": getattr(config, "comm_script_name_prefix", None),
            "queued_jobs": getattr(config, "comm_script_queued_jobs_cmd", None),
            "cancel_job": getattr(config, "comm_script_cancel_job_cmd", None),
            "submit_job": getattr(config, "comm_script_submit_job_cmd", None),
        },
        "process_search": {
            "minimization_offset": config.process_search_minimization_offset,
            "minimize_first": config.process_search_minimize_first,
        },
        "prefactor": {
            "default_value": config.process_search_default_prefactor,
        },
        "saddle_search": {
            "method": config.saddle_method,
            "max_iterations": config.saddle_search_max_iterations,
            "dynamics_temperature": config.saddle_dynamics_temperature,
            "displace_random_weight": config.displace_random_weight,
            "displace_listed_atom_weight": config.displace_listed_atom_weight,
            "displace_listed_type_weight": config.displace_listed_type_weight,
            "displace_all_listed": config.displace_all_listed,
            "displace_under_coordinated_weight": config.displace_under_coordinated_weight,
            "displace_least_coordinated_weight": config.displace_least_coordinated_weight,
            "displace_not_FCC_HCP_weight": config.displace_not_FCC_HCP_weight,
            "displace_not_TCP_BCC_weight": config.displace_not_TCP_BCC_weight,
            "displace_not_TCP_weight": config.displace_not_TCP_weight,
            "displace_water_weight": config.displace_water_weight,
            "stdev_translation": config.stdev_translation,
            "stdev_rotation": config.stdev_rotation,
            "molecule_list": config.molecule_list,
            "disp_at_random": config.disp_at_random,
            "displace_magnitude": config.disp_magnitude,
            "displace_radius": config.disp_radius,
            "displace_min_norm": config.disp_min_norm,
            "void_bias_fraction": config.void_bias_fraction,
            "displace_max_coordination": config.disp_max_coord,
            "random_mode": config.random_mode,
            "displace_atom_list": getattr(config, "disp_listed_atoms", []),
            "displace_type_list": getattr(config, "disp_listed_types", []),
            "displace_1d": config.displace_1d,
            "dynamics_max_init_curvature": config.dynamics_max_init_curvature,
            "zero_mode_abort_curvature": config.zero_mode_abort_curvature,
        },
        "kdb": {
            "use_kdb": config.kdb_on,
            "kdb_only": config.kdb_only,
            "kdb_scratch_path": config.kdb_scratch_path,
            "kdb_path": config.kdb_path,
            "remove_duplicates": config.kdb_nodupes,
            "kdb_name": config.kdb_name,
            "kdb_nf": config.kdb_nf,
            "kdb_dc": config.kdb_dc,
            "kdb_mac": config.kdb_mac,
        },
        "recycling": {
            "use_recycling": config.recycling_on,
            "save_suggestions": config.recycling_save_sugg,
            "displace_moved_only": config.disp_moved_only,
            "move_distance": config.recycling_move_distance,
            "active_region": config.recycling_active_region,
            "mass_weight_factor": config.recycling_mass_weight_factor,
            "use_sb_recycling": config.sb_recycling_on,
            "superbasin_recycling": config.sb_recycling_path,
        },
        "coarse_graining": {
            "use_mcamc": config.sb_on,
            "state_file": config.sb_state_file,
            "superbasin_scheme": config.sb_scheme,
            "max_size": config.sb_max_size,
            "number_of_transitions": getattr(config, "sb_tc_ntrans", None),
            "energy_increment": getattr(config, "sb_el_energy_increment", None),
            "superbasin_confidence": config.sb_superbasin_confidence,
            "use_askmc": config.askmc_on,
            "askmc_confidence": getattr(config, "askmc_confidence", None),
            "askmc_barrier_raise_param": getattr(config, "askmc_alpha", None),
            "askmc_high_barrier_def": getattr(config, "askmc_gamma", None),
            "askmc_barrier_test_on": getattr(config, "askmc_barrier_test_on", None),
            "askmc_connections_test_on": getattr(
                config, "askmc_connections_test_on", None
            ),
        },
        "optimizer": {
            "max_iterations": config.optimizers_max_iterations,
        },
        "debug": {
            "interactive_shell": config.debug_interactive_shell,
            "keep_bad_saddles": config.debug_keep_bad_saddles,
            "keep_all_result_files": config.debug_keep_all_results,
            "result_files_path": config.debug_results_path,
            "register_extra_results": config.debug_register_extra_results,
            "use_mean_time": config.debug_use_mean_time,
            "target_trajectory": config.debug_target_trajectory,
            "stop_criterion": config.debug_stop_criterion,
        },
    }


# Function to convert INI configuration to Pydantic model and then to TOML
def ini_to_toml(ini_file: str, toml_file: str, exclude_defaults: bool):
    config = ConfigClass()
    config.init(ini_file)

    try:
        transformed_config = transform_old_config_to_dict(config)
        pydantic_defaults = {
            "main": MainConfig(),
            "structure_comparison": StructureComparisonConfig(),
            "akmc": AKMCConfig(),
            "basin_hopping": BasinHoppingConfig(),
            "paths": PathsConfig(),
            "communicator": CommunicatorConfig(),
            "process_search": ProcessSearchConfig(),
            "potential": PotentialConfig(),
            "refine": RefineConfig(),
            "prefactor": PrefactorConfig(),
            "saddle_search": SaddleSearchConfig(),
            "kdb": KDBConfig(),
            "recycling": RecyclingConfig(),
            "coarse_graining": CoarseGrainingConfig(),
            "optimizer": OptimizerConfig(),
            "debug": DebugConfig(),
        }

        updated_config = {
            key: update_model_with_config(
                pydantic_defaults[key], transformed_config.get(key, {})
            )
            for key in pydantic_defaults.keys()
        }

        pydantic_config = Config(**updated_config)
    except ValidationError as e:
        print("Configuration validation error:", e)
        return

    # Convert Pydantic model to dictionary
    config_dict = pydantic_config.dict(exclude_unset=exclude_defaults)

    print(config_dict)
    # Write to TOML file
    with open(toml_file, "w") as f:
        toml.dump(config_dict, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert INI configuration file to TOML format."
    )
    parser.add_argument("ini_file", help="The input INI file to convert.")
    parser.add_argument("toml_file", help="The output TOML file.")
    parser.add_argument(
        "--exclude-defaults",
        action="store_true",
        help="Exclude default values from the output TOML file.",
    )

    args = parser.parse_args()

    ini_to_toml(args.ini_file, args.toml_file, args.exclude_defaults)
    if Path(args.toml_file).exists():
        print(f"Converted {args.ini_file} to {args.toml_file} successfully.")
    else:
        print("Failed conversion.")
