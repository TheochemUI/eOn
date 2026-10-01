"""AUTO-GENERATED from schema/eon_params.capnp — do not edit.

Regenerate: python tools/params_ssot/codegen.py
"""
from __future__ import annotations

CATALOG = {
  "flat_aliases": [
    {
      "default": 20,
      "flat_key": "lbfgs_memory",
      "ssot_path": "Optimizer.LBFGS.memory",
      "yaml_section": "Optimizer"
    },
    {
      "default": 0.01,
      "flat_key": "lbfgs_inverse_curvature",
      "ssot_path": "Optimizer.LBFGS.inverse_curvature",
      "yaml_section": "Optimizer"
    },
    {
      "default": 0.0,
      "flat_key": "lbfgs_max_inverse_curvature",
      "ssot_path": "Optimizer.LBFGS.max_inverse_curvature",
      "yaml_section": "Optimizer"
    },
    {
      "default": True,
      "flat_key": "lbfgs_auto_scale",
      "ssot_path": "Optimizer.LBFGS.auto_scale",
      "yaml_section": "Optimizer"
    },
    {
      "default": True,
      "flat_key": "lbfgs_angle_reset",
      "ssot_path": "Optimizer.LBFGS.angle_reset",
      "yaml_section": "Optimizer"
    },
    {
      "default": True,
      "flat_key": "lbfgs_distance_reset",
      "ssot_path": "Optimizer.LBFGS.distance_reset",
      "yaml_section": "Optimizer"
    },
    {
      "default": "reset",
      "flat_key": "lbfgs_curvature",
      "ssot_path": "Optimizer.LBFGS.curvature",
      "yaml_section": "Optimizer"
    },
    {
      "default": False,
      "flat_key": "lbfgs_project_rigid",
      "ssot_path": "Optimizer.LBFGS.project_rigid",
      "yaml_section": "Optimizer"
    },
    {
      "default": "standard",
      "flat_key": "lbfgs_secant",
      "ssot_path": "Optimizer.LBFGS.secant",
      "yaml_section": "Optimizer"
    },
    {
      "default": "none",
      "flat_key": "lbfgs_precon",
      "ssot_path": "Optimizer.LBFGS.precon",
      "yaml_section": "Optimizer"
    },
    {
      "default": "sy_yy",
      "flat_key": "lbfgs_h0",
      "ssot_path": "Optimizer.LBFGS.h0",
      "yaml_section": "Optimizer"
    },
    {
      "default": "none",
      "flat_key": "lbfgs_accept",
      "ssot_path": "Optimizer.LBFGS.accept",
      "yaml_section": "Optimizer"
    },
    {
      "default": 0,
      "flat_key": "lbfgs_extra_updates",
      "ssot_path": "Optimizer.LBFGS.extra_updates",
      "yaml_section": "Optimizer"
    },
    {
      "default": 1e-06,
      "flat_key": "lbfgs_cautious_eps",
      "ssot_path": "Optimizer.LBFGS.cautious_eps",
      "yaml_section": "Optimizer"
    },
    {
      "default": 0.01,
      "flat_key": "lbfgs_cautious_alpha",
      "ssot_path": "Optimizer.LBFGS.cautious_alpha",
      "yaml_section": "Optimizer"
    },
    {
      "default": 3.0,
      "flat_key": "lbfgs_precon_A",
      "ssot_path": "Optimizer.LBFGS.precon_A",
      "yaml_section": "Optimizer"
    },
    {
      "default": 1.0,
      "flat_key": "lbfgs_precon_mu",
      "ssot_path": "Optimizer.LBFGS.precon_mu",
      "yaml_section": "Optimizer"
    },
    {
      "default": 0.0,
      "flat_key": "lbfgs_precon_rcut",
      "ssot_path": "Optimizer.LBFGS.precon_rcut",
      "yaml_section": "Optimizer"
    },
    {
      "default": "lbfgs",
      "flat_key": "lbfgs_step",
      "ssot_path": "Optimizer.LBFGS.step",
      "yaml_section": "Optimizer"
    },
    {
      "default": False,
      "flat_key": "cg_no_overshooting",
      "ssot_path": "Optimizer.CG.no_overshooting",
      "yaml_section": "Optimizer"
    },
    {
      "default": False,
      "flat_key": "cg_knock_out_max_move",
      "ssot_path": "Optimizer.CG.knock_out_max_move",
      "yaml_section": "Optimizer"
    },
    {
      "default": False,
      "flat_key": "cg_line_search",
      "ssot_path": "Optimizer.CG.line_search",
      "yaml_section": "Optimizer"
    },
    {
      "default": 0.1,
      "flat_key": "cg_line_converged",
      "ssot_path": "Optimizer.CG.line_converged",
      "yaml_section": "Optimizer"
    },
    {
      "default": 10,
      "flat_key": "cg_max_iter_line_search",
      "ssot_path": "Optimizer.CG.line_search_max_iter",
      "yaml_section": "Optimizer"
    },
    {
      "default": 0,
      "flat_key": "cg_max_iter_before_reset",
      "ssot_path": "Optimizer.CG.max_iter_before_reset",
      "yaml_section": "Optimizer"
    },
    {
      "default": False,
      "flat_key": "qm_steepest_descent",
      "ssot_path": "Optimizer.Quickmin.steepest_descent",
      "yaml_section": "Optimizer"
    },
    {
      "default": 0.1,
      "flat_key": "sd_alpha",
      "ssot_path": "Optimizer.SD.alpha",
      "yaml_section": "Optimizer"
    },
    {
      "default": False,
      "flat_key": "sd_two_point",
      "ssot_path": "Optimizer.SD.two_point",
      "yaml_section": "Optimizer"
    },
    {
      "default": "lbfgs",
      "flat_key": "xtsci_method",
      "ssot_path": "Optimizer.Xtsci.method",
      "yaml_section": "Optimizer"
    },
    {
      "default": "lbfgs",
      "flat_key": "xtsci_qn_step",
      "ssot_path": "Optimizer.Xtsci.qn_step",
      "yaml_section": "Optimizer"
    },
    {
      "default": "none",
      "flat_key": "xtsci_precon",
      "ssot_path": "Optimizer.Xtsci.precon",
      "yaml_section": "Optimizer"
    },
    {
      "default": "none",
      "flat_key": "xtsci_accept",
      "ssot_path": "Optimizer.Xtsci.accept",
      "yaml_section": "Optimizer"
    },
    {
      "default": False,
      "flat_key": "xtsci_highs",
      "ssot_path": "Optimizer.Xtsci.highs",
      "yaml_section": "Optimizer"
    },
    {
      "default": "euclidean",
      "flat_key": "xtsci_manifold",
      "ssot_path": "Optimizer.Xtsci.manifold",
      "yaml_section": "Optimizer"
    }
  ],
  "schema_version": 1,
  "sections": {
    "Main": {
      "fields": [
        {
          "capnp": "job",
          "default": "process_search",
          "ordinal": 0,
          "snake": "job",
          "type": "Text"
        },
        {
          "capnp": "randomSeed",
          "default": -1,
          "ordinal": 1,
          "snake": "random_seed",
          "type": "Int64"
        },
        {
          "capnp": "temperature",
          "default": 300.0,
          "ordinal": 2,
          "snake": "temperature",
          "type": "Float64"
        },
        {
          "capnp": "quiet",
          "default": False,
          "ordinal": 3,
          "snake": "quiet",
          "type": "Bool"
        },
        {
          "capnp": "writeLog",
          "default": True,
          "ordinal": 4,
          "snake": "write_log",
          "type": "Bool"
        },
        {
          "capnp": "checkpoint",
          "default": False,
          "ordinal": 5,
          "snake": "checkpoint",
          "type": "Bool"
        },
        {
          "capnp": "iniFilename",
          "default": "config.ini",
          "ordinal": 6,
          "snake": "ini_filename",
          "type": "Text"
        },
        {
          "capnp": "conFilename",
          "default": "pos.con",
          "ordinal": 7,
          "snake": "con_filename",
          "type": "Text"
        },
        {
          "capnp": "finiteDifference",
          "default": 0.01,
          "ordinal": 8,
          "snake": "finite_difference",
          "type": "Float64"
        },
        {
          "capnp": "maxForceCalls",
          "default": 0,
          "ordinal": 9,
          "snake": "max_force_calls",
          "type": "Int64"
        },
        {
          "capnp": "removeNetForce",
          "default": True,
          "ordinal": 10,
          "snake": "remove_net_force",
          "type": "Bool"
        },
        {
          "capnp": "writeConForces",
          "default": False,
          "ordinal": 11,
          "snake": "write_con_forces",
          "type": "Bool"
        },
        {
          "capnp": "parallel",
          "default": True,
          "ordinal": 12,
          "snake": "parallel",
          "type": "Bool"
        }
      ],
      "struct": "MainOptions"
    },
    "Optimizer": {
      "fields": [
        {
          "capnp": "optMethod",
          "default": "cg",
          "ordinal": 0,
          "snake": "opt_method",
          "type": "Text"
        },
        {
          "capnp": "convergenceMetric",
          "default": "norm",
          "ordinal": 1,
          "snake": "convergence_metric",
          "type": "Text"
        },
        {
          "capnp": "maxIterations",
          "default": 1000,
          "ordinal": 2,
          "snake": "max_iterations",
          "type": "UInt64"
        },
        {
          "capnp": "maxMove",
          "default": 0.2,
          "ordinal": 3,
          "snake": "max_move",
          "type": "Float64"
        },
        {
          "capnp": "convergedForce",
          "default": 0.01,
          "ordinal": 4,
          "snake": "converged_force",
          "type": "Float64"
        },
        {
          "capnp": "timeStep",
          "default": 1.0,
          "ordinal": 5,
          "snake": "time_step",
          "type": "Float64"
        },
        {
          "capnp": "maxTimeStep",
          "default": 2.5,
          "ordinal": 6,
          "snake": "max_time_step",
          "type": "Float64"
        },
        {
          "capnp": "lbfgs",
          "default": None,
          "nested": True,
          "ordinal": 7,
          "snake": "lbfgs",
          "type": "OptimizerLbfgsOptions"
        },
        {
          "capnp": "cg",
          "default": None,
          "nested": True,
          "ordinal": 8,
          "snake": "cg",
          "type": "OptimizerCgOptions"
        },
        {
          "capnp": "quickmin",
          "default": None,
          "nested": True,
          "ordinal": 9,
          "snake": "quickmin",
          "type": "OptimizerQuickminOptions"
        },
        {
          "capnp": "sd",
          "default": None,
          "nested": True,
          "ordinal": 10,
          "snake": "sd",
          "type": "OptimizerSdOptions"
        },
        {
          "capnp": "xtsci",
          "default": None,
          "nested": True,
          "ordinal": 11,
          "snake": "xtsci",
          "type": "OptimizerXtsciOptions"
        }
      ],
      "struct": "OptimizerOptions"
    },
    "Optimizer.CG": {
      "fields": [
        {
          "capnp": "noOvershooting",
          "default": False,
          "ordinal": 0,
          "snake": "no_overshooting",
          "type": "Bool"
        },
        {
          "capnp": "knockOutMaxMove",
          "default": False,
          "ordinal": 1,
          "snake": "knock_out_max_move",
          "type": "Bool"
        },
        {
          "capnp": "lineSearch",
          "default": False,
          "ordinal": 2,
          "snake": "line_search",
          "type": "Bool"
        },
        {
          "capnp": "lineConverged",
          "default": 0.1,
          "ordinal": 3,
          "snake": "line_converged",
          "type": "Float64"
        },
        {
          "capnp": "lineSearchMaxIter",
          "default": 10,
          "ordinal": 4,
          "snake": "line_search_max_iter",
          "type": "Int64"
        },
        {
          "capnp": "maxIterBeforeReset",
          "default": 0,
          "ordinal": 5,
          "snake": "max_iter_before_reset",
          "type": "Int64"
        }
      ],
      "struct": "OptimizerCgOptions"
    },
    "Optimizer.LBFGS": {
      "fields": [
        {
          "capnp": "memory",
          "default": 20,
          "ordinal": 0,
          "snake": "memory",
          "type": "Int64"
        },
        {
          "capnp": "inverseCurvature",
          "default": 0.01,
          "ordinal": 1,
          "snake": "inverse_curvature",
          "type": "Float64"
        },
        {
          "capnp": "maxInverseCurvature",
          "default": 0.0,
          "ordinal": 2,
          "snake": "max_inverse_curvature",
          "type": "Float64"
        },
        {
          "capnp": "autoScale",
          "default": True,
          "ordinal": 3,
          "snake": "auto_scale",
          "type": "Bool"
        },
        {
          "capnp": "angleReset",
          "default": True,
          "ordinal": 4,
          "snake": "angle_reset",
          "type": "Bool"
        },
        {
          "capnp": "distanceReset",
          "default": True,
          "ordinal": 5,
          "snake": "distance_reset",
          "type": "Bool"
        },
        {
          "capnp": "curvature",
          "default": "reset",
          "ordinal": 6,
          "snake": "curvature",
          "type": "Text"
        },
        {
          "capnp": "projectRigid",
          "default": False,
          "ordinal": 7,
          "snake": "project_rigid",
          "type": "Bool"
        },
        {
          "capnp": "secant",
          "default": "standard",
          "ordinal": 8,
          "snake": "secant",
          "type": "Text"
        },
        {
          "capnp": "precon",
          "default": "none",
          "ordinal": 9,
          "snake": "precon",
          "type": "Text"
        },
        {
          "capnp": "h0",
          "default": "sy_yy",
          "ordinal": 10,
          "snake": "h0",
          "type": "Text"
        },
        {
          "capnp": "accept",
          "default": "none",
          "ordinal": 11,
          "snake": "accept",
          "type": "Text"
        },
        {
          "capnp": "extraUpdates",
          "default": 0,
          "ordinal": 12,
          "snake": "extra_updates",
          "type": "Int64"
        },
        {
          "capnp": "cautiousEps",
          "default": 1e-06,
          "ordinal": 13,
          "snake": "cautious_eps",
          "type": "Float64"
        },
        {
          "capnp": "cautiousAlpha",
          "default": 0.01,
          "ordinal": 14,
          "snake": "cautious_alpha",
          "type": "Float64"
        },
        {
          "capnp": "preconA",
          "default": 3.0,
          "ordinal": 15,
          "snake": "precon_A",
          "type": "Float64"
        },
        {
          "capnp": "preconMu",
          "default": 1.0,
          "ordinal": 16,
          "snake": "precon_mu",
          "type": "Float64"
        },
        {
          "capnp": "preconRcut",
          "default": 0.0,
          "ordinal": 17,
          "snake": "precon_rcut",
          "type": "Float64"
        },
        {
          "capnp": "step",
          "default": "lbfgs",
          "ordinal": 18,
          "snake": "step",
          "type": "Text"
        }
      ],
      "struct": "OptimizerLbfgsOptions"
    },
    "Optimizer.Quickmin": {
      "fields": [
        {
          "capnp": "steepestDescent",
          "default": False,
          "ordinal": 0,
          "snake": "steepest_descent",
          "type": "Bool"
        }
      ],
      "struct": "OptimizerQuickminOptions"
    },
    "Optimizer.SD": {
      "fields": [
        {
          "capnp": "alpha",
          "default": 0.1,
          "ordinal": 0,
          "snake": "alpha",
          "type": "Float64"
        },
        {
          "capnp": "twoPoint",
          "default": False,
          "ordinal": 1,
          "snake": "two_point",
          "type": "Bool"
        }
      ],
      "struct": "OptimizerSdOptions"
    },
    "Optimizer.Xtsci": {
      "fields": [
        {
          "capnp": "method",
          "default": "lbfgs",
          "ordinal": 0,
          "snake": "method",
          "type": "Text"
        },
        {
          "capnp": "qnStep",
          "default": "lbfgs",
          "ordinal": 1,
          "snake": "qn_step",
          "type": "Text"
        },
        {
          "capnp": "precon",
          "default": "none",
          "ordinal": 2,
          "snake": "precon",
          "type": "Text"
        },
        {
          "capnp": "accept",
          "default": "none",
          "ordinal": 3,
          "snake": "accept",
          "type": "Text"
        },
        {
          "capnp": "highs",
          "default": False,
          "ordinal": 4,
          "snake": "highs",
          "type": "Bool"
        },
        {
          "capnp": "manifold",
          "default": "euclidean",
          "ordinal": 5,
          "snake": "manifold",
          "type": "Text"
        }
      ],
      "struct": "OptimizerXtsciOptions"
    },
    "Potential": {
      "fields": [
        {
          "capnp": "potential",
          "default": "lj",
          "ordinal": 0,
          "snake": "potential",
          "type": "Text"
        },
        {
          "capnp": "mpiPollPeriod",
          "default": 0.25,
          "ordinal": 1,
          "snake": "mpi_poll_period",
          "type": "Float64"
        },
        {
          "capnp": "lammpsLogging",
          "default": False,
          "ordinal": 2,
          "snake": "lammps_logging",
          "type": "Bool"
        },
        {
          "capnp": "lammpsThreads",
          "default": 0,
          "ordinal": 3,
          "snake": "lammps_threads",
          "type": "Int32"
        },
        {
          "capnp": "emtRasmussen",
          "default": False,
          "ordinal": 4,
          "snake": "emt_rasmussen",
          "type": "Bool"
        },
        {
          "capnp": "logPotential",
          "default": False,
          "ordinal": 5,
          "snake": "log_potential",
          "type": "Bool"
        },
        {
          "capnp": "extPotPath",
          "default": "./ext_pot",
          "ordinal": 6,
          "snake": "ext_pot_path",
          "type": "Text"
        },
        {
          "capnp": "potentialsPath",
          "default": "",
          "ordinal": 7,
          "snake": "potentials_path",
          "type": "Text"
        }
      ],
      "struct": "PotentialOptions"
    },
    "Process Search": {
      "fields": [
        {
          "capnp": "minimizeFirst",
          "default": True,
          "ordinal": 0,
          "snake": "minimize_first",
          "type": "Bool"
        },
        {
          "capnp": "minimizationOffset",
          "default": 0.2,
          "ordinal": 1,
          "snake": "minimization_offset",
          "type": "Float64"
        }
      ],
      "struct": "ProcessSearchOptions"
    },
    "RgpotPot": {
      "accessor": "rgpot_options",
      "fields": [
        {
          "capnp": "backend",
          "default": "nwchemc",
          "ordinal": 0,
          "snake": "backend",
          "type": "Text"
        },
        {
          "capnp": "basis",
          "default": "sto-3g",
          "ini_order": [
            "basis",
            "nwchem_basis"
          ],
          "ordinal": 1,
          "snake": "basis",
          "type": "Text"
        },
        {
          "capnp": "theory",
          "default": "scf",
          "ini_order": [
            "theory",
            "nwchem_theory"
          ],
          "ordinal": 2,
          "snake": "theory",
          "type": "Text"
        },
        {
          "capnp": "scf_type",
          "default": "rhf",
          "ini_order": [
            "scf_type",
            "nwchem_scf_type"
          ],
          "ordinal": 3,
          "snake": "scf_type",
          "type": "Text"
        },
        {
          "capnp": "functional",
          "default": "BLYP",
          "ini_order": [
            "functional",
            "cpmd_functional"
          ],
          "ordinal": 4,
          "snake": "functional",
          "type": "Text"
        },
        {
          "capnp": "cutoff_ry",
          "default": 70.0,
          "ini_order": [
            "cutOffRy",
            "cutoff_ry",
            "cpmd_cut_off_ry"
          ],
          "ordinal": 5,
          "snake": "cutoff_ry",
          "type": "Float64"
        },
        {
          "capnp": "charge",
          "default": 0,
          "ini_order": [
            "charge",
            "nwchem_charge"
          ],
          "ordinal": 6,
          "overlay_order": [
            "charge"
          ],
          "snake": "charge",
          "type": "Int32"
        },
        {
          "capnp": "multiplicity",
          "default": 1,
          "ini_order": [
            "multiplicity",
            "nwchem_multiplicity"
          ],
          "ordinal": 7,
          "overlay_order": [
            "multiplicity"
          ],
          "snake": "multiplicity",
          "type": "Int32"
        },
        {
          "capnp": "engine_path",
          "default": "",
          "ordinal": 8,
          "snake": "engine_path",
          "type": "Text"
        },
        {
          "capnp": "engine_library",
          "default": "",
          "ordinal": 9,
          "snake": "engine_library",
          "type": "Text"
        },
        {
          "capnp": "engine_root",
          "default": "",
          "ordinal": 10,
          "snake": "engine_root",
          "type": "Text"
        },
        {
          "capnp": "title",
          "default": "",
          "ordinal": 11,
          "snake": "title",
          "type": "Text"
        },
        {
          "capnp": "memory_mb",
          "default": 0,
          "ordinal": 12,
          "snake": "memory_mb",
          "type": "Int32"
        },
        {
          "capnp": "scratch_dir",
          "default": "",
          "ordinal": 13,
          "snake": "scratch_dir",
          "type": "Text"
        },
        {
          "capnp": "input_block",
          "default": "",
          "ordinal": 14,
          "snake": "input_block",
          "type": "Text"
        },
        {
          "capnp": "permanent_dir",
          "default": "",
          "ordinal": 15,
          "snake": "permanent_dir",
          "type": "Text"
        },
        {
          "capnp": "params_path",
          "default": "",
          "ordinal": 16,
          "snake": "params_path",
          "type": "Text"
        },
        {
          "capnp": "ranks_per_image",
          "default": 0,
          "ordinal": 17,
          "snake": "ranks_per_image",
          "type": "Int32"
        },
        {
          "capnp": "model_path",
          "default": "",
          "ordinal": 18,
          "snake": "model_path",
          "type": "Text"
        },
        {
          "capnp": "device",
          "default": "cpu",
          "ordinal": 19,
          "snake": "device",
          "type": "Text"
        },
        {
          "capnp": "length_unit",
          "default": "angstrom",
          "ordinal": 20,
          "snake": "length_unit",
          "type": "Text"
        },
        {
          "capnp": "extensions_directory",
          "default": "",
          "ordinal": 21,
          "snake": "extensions_directory",
          "type": "Text"
        },
        {
          "capnp": "check_consistency",
          "default": False,
          "ordinal": 22,
          "snake": "check_consistency",
          "type": "Bool"
        },
        {
          "capnp": "uncertainty_threshold",
          "default": -1.0,
          "ordinal": 23,
          "snake": "uncertainty_threshold",
          "type": "Float64"
        },
        {
          "capnp": "torch_determinism_strict",
          "default": False,
          "ordinal": 24,
          "snake": "torch_determinism_strict",
          "type": "Bool"
        },
        {
          "capnp": "xtb_paramset",
          "default": "GFN2xTB",
          "ini_order": [
            "paramset",
            "xtb_paramset"
          ],
          "ordinal": 25,
          "snake": "xtb_paramset",
          "type": "Text"
        },
        {
          "capnp": "xtb_accuracy",
          "default": 1.0,
          "ini_order": [
            "accuracy",
            "xtb_accuracy"
          ],
          "ordinal": 26,
          "snake": "xtb_accuracy",
          "type": "Float64"
        },
        {
          "capnp": "xtb_electronic_temperature",
          "default": 300.0,
          "ini_order": [
            "electronic_temperature",
            "xtb_electronic_temperature"
          ],
          "ordinal": 27,
          "snake": "xtb_electronic_temperature",
          "type": "Float64"
        },
        {
          "capnp": "xtb_max_iterations",
          "cxx_cast": "static_cast<int>",
          "default": 250,
          "ini_order": [
            "max_iterations",
            "xtb_max_iterations"
          ],
          "ordinal": 28,
          "snake": "xtb_max_iterations",
          "type": "Int32"
        },
        {
          "capnp": "xtb_charge",
          "default": 0.0,
          "ini_fallback_member": "charge",
          "ordinal": 29,
          "snake": "xtb_charge",
          "type": "Float64"
        },
        {
          "capnp": "xtb_uhf",
          "cxx_cast": "static_cast<int>",
          "default": 0,
          "ini_order": [
            "uhf",
            "xtb_uhf"
          ],
          "ordinal": 30,
          "snake": "xtb_uhf",
          "type": "Int32"
        }
      ],
      "ini_when": "potential=RGPOT",
      "json_name": "RgpotPot",
      "overlays": [
        {
          "map": [
            {
              "key": "paramset",
              "member": "xtb_paramset"
            },
            {
              "key": "accuracy",
              "member": "xtb_accuracy"
            },
            {
              "key": "electronic_temperature",
              "member": "xtb_electronic_temperature"
            },
            {
              "key": "max_iterations",
              "member": "xtb_max_iterations"
            },
            {
              "key": "uhf",
              "member": "xtb_uhf"
            },
            {
              "key": "charge",
              "member": "xtb_charge"
            }
          ],
          "section": "XTBPot",
          "when": "backend:xtb,xtbpot,gfn,gfnxtb"
        },
        {
          "fields": [
            "functional",
            "cutoff_ry",
            "charge",
            "multiplicity",
            "title",
            "memory_mb",
            "input_block"
          ],
          "section": "cpmd",
          "when": "backend:cpmd,cpmdc,cpmdpot"
        }
      ],
      "project": [
        "ini",
        "json"
      ],
      "struct": "RgpotPotOptions"
    },
    "Structure Comparison": {
      "fields": [
        {
          "capnp": "distanceDifference",
          "default": 0.1,
          "ordinal": 0,
          "snake": "distance_difference",
          "type": "Float64"
        },
        {
          "capnp": "neighborCutoff",
          "default": 3.3,
          "ordinal": 1,
          "snake": "neighbor_cutoff",
          "type": "Float64"
        },
        {
          "capnp": "checkRotation",
          "default": False,
          "ordinal": 2,
          "snake": "check_rotation",
          "type": "Bool"
        },
        {
          "capnp": "indistinguishableAtoms",
          "default": True,
          "ordinal": 3,
          "snake": "indistinguishable_atoms",
          "type": "Bool"
        },
        {
          "capnp": "energyDifference",
          "default": 0.01,
          "ordinal": 4,
          "snake": "energy_difference",
          "type": "Float64"
        },
        {
          "capnp": "removeTranslation",
          "default": True,
          "ordinal": 5,
          "snake": "remove_translation",
          "type": "Bool"
        },
        {
          "capnp": "useCovalent",
          "default": False,
          "ordinal": 6,
          "snake": "use_covalent",
          "type": "Bool"
        },
        {
          "capnp": "covalentScale",
          "default": 1.3,
          "ordinal": 7,
          "snake": "covalent_scale",
          "type": "Float64"
        },
        {
          "capnp": "bruteNeighbors",
          "default": False,
          "ordinal": 8,
          "snake": "brute_neighbors",
          "type": "Bool"
        }
      ],
      "struct": "StructureComparisonOptions"
    }
  },
  "server_yaml_default_overrides": {
    "Main.job": "akmc",
    "Potential.ext_pot_path": "ext_pot"
  },
  "source": "schema/eon_params.capnp"
}
