/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#include "ParametersINI.h"
#include "BaseStructures.h"
#include "EpiCenters.h"
#include "HelperFunctions.h"
#include "Parameters.h"
#include "magic_enum/magic_enum.hpp"

#include <INIReader.h>

#include <cerrno>
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>
#include <string>

#include "EonLogger.h"

namespace {
std::string toLowerCase(std::string s) {
  for (std::string::size_type i = 0; i < s.length(); ++i) {
    s[i] = tolower(s[i]);
  }
  return s;
}
} // namespace

namespace eonc::config {

int load_ini(INIReader &ini, Parameters &params) {
  int error = 0;

  // [Main] //

  params.main_options.job =
      magic_enum::enum_cast<JobType>(ini.Get("Main", "job", ""),
                                     magic_enum::case_insensitive)
          .value_or(JobType::Unknown);
  params.main_options.temperature =
      ini.GetReal("Main", "temperature", params.main_options.temperature);
  params.main_options.randomSeed =
      ini.GetInteger("Main", "random_seed", params.main_options.randomSeed);
  params.main_options.checkpoint =
      ini.GetBoolean("Main", "checkpoint", params.main_options.checkpoint);
  params.main_options.quiet =
      ini.GetBoolean("Main", "quiet", params.main_options.quiet);
  params.main_options.writeLog =
      ini.GetBoolean("Main", "write_log", params.main_options.writeLog);
  params.main_options.finiteDifference = ini.GetReal(
      "Main", "finite_difference", params.main_options.finiteDifference);
  // Initialize random generator
  if (params.main_options.randomSeed < 0) {
    unsigned i = static_cast<unsigned>(std::time(nullptr));
    params.main_options.randomSeed = i;
    eonc::helpers::random(i);
  } else {
    eonc::helpers::random(params.main_options.randomSeed);
  }
  params.main_options.maxForceCalls = ini.GetInteger(
      "Main", "max_force_calls", params.main_options.maxForceCalls);
  params.main_options.removeNetForce = ini.GetBoolean(
      "Main", "remove_net_force", params.main_options.removeNetForce);
  params.main_options.parallel =
      ini.GetBoolean("Main", "parallel", params.main_options.parallel);

  // [Potential] //

  params.potential_options.potential =
      magic_enum::enum_cast<PotType>(ini.Get("Potential", "potential", ""),
                                     magic_enum::case_insensitive)
          .value_or(PotType::UNKNOWN);
  params.potential_options.MPIPollPeriod = ini.GetReal(
      "Potential", "mpi_poll_period", params.potential_options.MPIPollPeriod);
  params.potential_options.LAMMPSLogging = ini.GetBoolean(
      "Potential", "lammps_logging", params.potential_options.LAMMPSLogging);
  params.potential_options.LAMMPSThreads = static_cast<int>(ini.GetInteger(
      "Potential", "lammps_threads", params.potential_options.LAMMPSThreads));
  params.potential_options.EMTRasmussen = ini.GetBoolean(
      "Potential", "emt_rasmussen", params.potential_options.EMTRasmussen);
  params.potential_options.extPotPath =
      ini.Get("Potential", "ext_pot_path", params.potential_options.extPotPath);
  params.potential_options.potentialsPath = ini.Get(
      "Potential", "potentials_path", params.potential_options.potentialsPath);

  if (params.potential_options.potential == PotType::MPI ||
      params.potential_options.potential == PotType::VASP) {
    params.potential_options.LogPotential = true;
  } else {
    params.potential_options.LogPotential = false;
  }
  params.potential_options.LogPotential = ini.GetBoolean(
      "Potential", "log_potential", params.potential_options.LogPotential);

  // [AMS]
  if (params.potential_options.potential == PotType::AMS) {
    params.ams_options.engine =
        ini.Get("AMS", "engine", params.ams_options.engine);
    params.ams_options.forcefield =
        ini.Get("AMS", "forcefield", params.ams_options.forcefield);
    params.ams_options.resources =
        ini.Get("AMS", "resources", params.ams_options.resources);
    params.ams_options.model =
        ini.Get("AMS", "model", params.ams_options.model);
    params.ams_options.xc = ini.Get("AMS", "xc", params.ams_options.xc);
    params.ams_options.basis =
        ini.Get("AMS", "basis", params.ams_options.basis);
  }
  // [AMS_IO]
  if (params.potential_options.potential == PotType::AMS_IO) {
    params.ams_options.engine =
        ini.Get("AMS_IO", "engine", params.ams_options.engine);
    params.ams_options.forcefield =
        ini.Get("AMS_IO", "forcefield", params.ams_options.forcefield);
    params.ams_options.model =
        ini.Get("AMS_IO", "model", params.ams_options.model);
    params.ams_options.xc = ini.Get("AMS_IO", "xc", params.ams_options.xc);
  }
  // [AMS_ENV]
  if (params.potential_options.potential == PotType::AMS_IO ||
      params.potential_options.potential == PotType::AMS) {
    params.ams_options.env.amshome =
        ini.Get("AMS_ENV", "amshome", params.ams_options.env.amshome);
    params.ams_options.env.scm_tmpdir =
        ini.Get("AMS_ENV", "scm_tmpdir", params.ams_options.env.scm_tmpdir);
    params.ams_options.env.scmlicense =
        ini.Get("AMS_ENV", "scmlicense", params.ams_options.env.scmlicense);
    params.ams_options.env.scm_pythondir = ini.Get(
        "AMS_ENV", "scm_pythondir", params.ams_options.env.scm_pythondir);
    params.ams_options.env.amsbin =
        ini.Get("AMS_ENV", "amsbin", params.ams_options.env.amsbin);
    params.ams_options.env.amsresources =
        ini.Get("AMS_ENV", "amsresources", params.ams_options.env.amsresources);
  }
  // [XTBPot]
  if (params.potential_options.potential == PotType::XTB) {
    params.xtb_options.paramset =
        ini.Get("XTBPot", "paramset", params.xtb_options.paramset);
    params.xtb_options.acc =
        ini.GetReal("XTBPot", "accuracy", params.xtb_options.acc);
    params.xtb_options.elec_temperature =
        ini.GetReal("XTBPot", "electronic_temperature",
                    params.xtb_options.elec_temperature);
    params.xtb_options.maxiter =
        ini.GetInteger("XTBPot", "max_iterations", params.xtb_options.maxiter);
    params.xtb_options.uhf =
        ini.GetInteger("XTBPot", "uhf", params.xtb_options.uhf);
    params.xtb_options.charge =
        ini.GetReal("XTBPot", "charge", params.xtb_options.charge);
  }
  // [ZBLPot]
  if (params.potential_options.potential == PotType::ZBL) {
    params.zbl_options.cut_inner =
        ini.GetReal("ZBLPot", "cut_inner", params.zbl_options.cut_inner);
    params.zbl_options.cut_global =
        ini.GetReal("ZBLPot", "cut_global", params.zbl_options.cut_global);
    if (params.zbl_options.cut_inner > params.zbl_options.cut_global) {
      throw std::runtime_error(
          "Switching function must begin before the global cutoff!");
    }
  }
  // [SocketNWChemPot]
  if (params.potential_options.potential == PotType::SocketNWChem) {
    params.socket_nwchem_options.host =
        ini.Get("SocketNWChemPot", "host", params.socket_nwchem_options.host);
    params.socket_nwchem_options.port = ini.GetInteger(
        "SocketNWChemPot", "port", params.socket_nwchem_options.port);
    params.socket_nwchem_options.mem_in_gb = ini.GetInteger(
        "SocketNWChemPot", "mem_in_gb", params.socket_nwchem_options.mem_in_gb);
    params.socket_nwchem_options.nwchem_settings =
        ini.Get("SocketNWChemPot", "nwchem_settings",
                params.socket_nwchem_options.nwchem_settings);
    params.socket_nwchem_options.unix_socket_path =
        ini.Get("SocketNWChemPot", "unix_socket_path",
                params.socket_nwchem_options.unix_socket_path);
    params.socket_nwchem_options.unix_socket_mode =
        ini.GetBoolean("SocketNWChemPot", "unix_socket_mode",
                       params.socket_nwchem_options.unix_socket_mode);
    params.socket_nwchem_options.make_template_input =
        ini.GetBoolean("SocketNWChemPot", "make_template_input",
                       params.socket_nwchem_options.make_template_input);
  }

  // [Debug] //

  params.debug_options.write_movies = ini.GetBoolean(
      "Debug", "write_movies", params.debug_options.write_movies);
  params.debug_options.write_movies_interval =
      ini.GetInteger("Debug", "write_movies_interval",
                     params.debug_options.write_movies_interval);
  params.debug_options.write_deprecated_outs =
      ini.GetBoolean("Debug", "write_deprecated_outs",
                     params.debug_options.write_deprecated_outs);
  params.debug_options.estimate_neb_eigenvalues =
      ini.GetBoolean("Debug", "estimate_neb_eigenvalues",
                     params.debug_options.estimate_neb_eigenvalues);
  params.debug_options.neb_mmf = toLowerCase(
      ini.Get("Debug", "neb_mmf_estimator", params.debug_options.neb_mmf));

  // [Structure Comparison] //

  params.structure_comparison_options.distance_difference =
      ini.GetReal("Structure Comparison", "distance_difference",
                  params.structure_comparison_options.distance_difference);
  params.structure_comparison_options.neighbor_cutoff =
      ini.GetReal("Structure Comparison", "neighbor_cutoff",
                  params.structure_comparison_options.neighbor_cutoff);
  params.structure_comparison_options.check_rotation =
      ini.GetBoolean("Structure Comparison", "check_rotation",
                     params.structure_comparison_options.check_rotation);
  params.structure_comparison_options.energy_difference =
      ini.GetReal("Structure Comparison", "energy_difference",
                  params.structure_comparison_options.energy_difference);
  params.structure_comparison_options.indistinguishable_atoms = ini.GetBoolean(
      "Structure Comparison", "indistinguishable_atoms",
      params.structure_comparison_options.indistinguishable_atoms);
  params.structure_comparison_options.remove_translation =
      ini.GetBoolean("Structure Comparison", "remove_translation",
                     params.structure_comparison_options.remove_translation);

  // [Process Search] //

  params.process_search_options.minimize_first =
      ini.GetBoolean("Process Search", "minimize_first",
                     params.process_search_options.minimize_first);
  params.process_search_options.minimization_offset =
      ini.GetReal("Process Search", "minimization_offset",
                  params.process_search_options.minimization_offset);

  // [Optimizers] //
  auto inp_optMethod =
      magic_enum::enum_cast<OptType>(ini.Get("Optimizer", "opt_method", "none"),
                                     magic_enum::case_insensitive)
          .value_or(OptType::Unknown);
  if (inp_optMethod != OptType::None) {
    params.optimizer_options.method = inp_optMethod;
  }

  params.optimizer_options.convergence_metric =
      toLowerCase(ini.Get("Optimizer", "convergence_metric",
                          params.optimizer_options.convergence_metric));
  if (params.optimizer_options.convergence_metric == "max_atom") {
    params.optimizer_options.convergence_metric_label = "Max atom force";
  } else if (params.optimizer_options.convergence_metric == "max_component") {
    params.optimizer_options.convergence_metric_label = "Max force comp";
  } else if (params.optimizer_options.convergence_metric == "norm") {
    params.optimizer_options.convergence_metric_label = "||Force||";
  } else {
    EONC_LOG_ERROR("unknown convergence_metric {}",
                   params.optimizer_options.convergence_metric);
    error = 1;
  }

  if (ini.HasSection("Refine")) {
    params.optimizer_options.refine.method =
        magic_enum::enum_cast<OptType>(ini.Get("Refine", "opt_method", ""),
                                       magic_enum::case_insensitive)
            .value_or(OptType::None);
    params.optimizer_options.refine.threshold = ini.GetReal(
        "Refine", "threshold", params.optimizer_options.refine.threshold);
  }

  params.optimizer_options.converged_force = ini.GetReal(
      "Optimizer", "converged_force", params.optimizer_options.converged_force);
  params.optimizer_options.max_iterations = static_cast<size_t>(ini.GetInteger(
      "Optimizer", "max_iterations", params.optimizer_options.max_iterations));
  params.optimizer_options.max_move =
      ini.GetReal("Optimizer", "max_move", params.optimizer_options.max_move);
  // Handle each optimizer separately
  if (ini.HasSection("QuickMin")) {
    params.optimizer_options.time_step_input = ini.GetReal(
        "QuickMin", "time_step", params.optimizer_options.time_step_input);
    params.optimizer_options.time_step =
        params.optimizer_options.time_step_input / params.constants.timeUnit;
    params.optimizer_options.quickmin.steepest_descent =
        ini.GetBoolean("Optimizer", "qm_steepest_descent",
                       params.optimizer_options.quickmin.steepest_descent);
  }
  if (ini.HasSection("FIRE")) {
    EONC_LOG_WARNING("Overwriting QuickMin timestep with Fire timestep!!");
    params.optimizer_options.time_step_input = ini.GetReal(
        "FIRE", "time_step", params.optimizer_options.time_step_input);
    params.optimizer_options.time_step =
        params.optimizer_options.time_step_input / params.constants.timeUnit;
    params.optimizer_options.max_time_step_input = ini.GetReal(
        "FIRE", "time_step_max", params.optimizer_options.max_time_step_input);
    params.optimizer_options.max_time_step =
        params.optimizer_options.max_time_step_input /
        params.constants.timeUnit;
  }
  if (ini.HasSection("LBFGS")) {
    params.optimizer_options.lbfgs.memory = ini.GetInteger(
        "LBFGS", "lbfgs_memory", params.optimizer_options.lbfgs.memory);
    params.optimizer_options.lbfgs.inverse_curvature =
        ini.GetReal("LBFGS", "lbfgs_inverse_curvature",
                    params.optimizer_options.lbfgs.inverse_curvature);
    params.optimizer_options.lbfgs.max_inverse_curvature =
        ini.GetReal("LBFGS", "lbfgs_max_inverse_curvature",
                    params.optimizer_options.lbfgs.max_inverse_curvature);
    params.optimizer_options.lbfgs.auto_scale = ini.GetBoolean(
        "LBFGS", "lbfgs_auto_scale", params.optimizer_options.lbfgs.auto_scale);
    params.optimizer_options.lbfgs.angle_reset =
        ini.GetBoolean("LBFGS", "lbfgs_angle_reset",
                       params.optimizer_options.lbfgs.angle_reset);
    params.optimizer_options.lbfgs.distance_reset =
        ini.GetBoolean("LBFGS", "lbfgs_distance_reset",
                       params.optimizer_options.lbfgs.distance_reset);
  }
  if (ini.HasSection("CG")) {
    params.optimizer_options.cg.no_overshooting =
        ini.GetBoolean("CG", "cg_no_overshooting",
                       params.optimizer_options.cg.no_overshooting);
    params.optimizer_options.cg.knock_out_max_move =
        ini.GetBoolean("CG", "cg_knock_out_max_move",
                       params.optimizer_options.cg.knock_out_max_move);
    params.optimizer_options.cg.line_search = ini.GetBoolean(
        "CG", "cg_line_search", params.optimizer_options.cg.line_search);
    params.optimizer_options.cg.line_converged = ini.GetReal(
        "CG", "cg_line_converged", params.optimizer_options.cg.line_converged);
    params.optimizer_options.cg.max_iter_before_reset =
        ini.GetInteger("CG", "cg_max_iter_before_reset",
                       params.optimizer_options.cg.max_iter_before_reset);
    params.optimizer_options.cg.line_search_max_iter =
        ini.GetInteger("CG", "cg_max_iter_line_search",
                       params.optimizer_options.cg.line_search_max_iter);
  }
  if (ini.HasSection("SD")) {
    params.optimizer_options.sd.alpha =
        ini.GetReal("SD", "sd_alpha", params.optimizer_options.sd.alpha);
    params.optimizer_options.sd.two_point = ini.GetBoolean(
        "SD", "sd_twopoint", params.optimizer_options.sd.two_point);
  }

  // [Dimer] //

  params.dimer_options.rotation_angle =
      ini.GetReal("Dimer", "finite_angle", params.dimer_options.rotation_angle);
  params.dimer_options.improved =
      ini.GetBoolean("Dimer", "improved", params.dimer_options.improved);
  params.dimer_options.converged_angle = ini.GetReal(
      "Dimer", "converged_angle", params.dimer_options.converged_angle);
  params.dimer_options.max_iterations = ini.GetInteger(
      "Dimer", "max_iterations", params.dimer_options.max_iterations);
  params.dimer_options.opt_method = toLowerCase(
      ini.Get("Dimer", "opt_method", params.dimer_options.opt_method));
  params.dimer_options.rotations_min = ini.GetInteger(
      "Dimer", "rotations_min", params.dimer_options.rotations_min);
  params.dimer_options.rotations_max = ini.GetInteger(
      "Dimer", "rotations_max", params.dimer_options.rotations_max);
  params.dimer_options.torque_min =
      ini.GetReal("Dimer", "torque_min", params.dimer_options.torque_min);
  params.dimer_options.torque_max =
      ini.GetReal("Dimer", "torque_max", params.dimer_options.torque_max);
  params.dimer_options.remove_rotation = ini.GetBoolean(
      "Dimer", "remove_rotation", params.dimer_options.remove_rotation);
  {
    const auto rotTok =
        toLowerCase(ini.Get("Dimer", "rotation_backend", "classical"));
    params.dimer_options.rotation_backend =
        magic_enum::enum_cast<DimerRotationBackend>(
            rotTok, magic_enum::case_insensitive)
            .value_or(DimerRotationBackend::Classical);
  }

  // GP Surrogate Parameters
  params.gp_surrogate_options.enabled =
      ini.GetBoolean("Surrogate", "use_surrogate", false);
  if (params.gp_surrogate_options.enabled) {
    params.gp_surrogate_options.sub_job = params.main_options.job;
    params.main_options.job = JobType::GP_Surrogate;
  }
  params.gp_surrogate_options.uncertainty = ini.GetReal(
      "Surrogate", "gp_uncertainty", params.gp_surrogate_options.uncertainty);
  if (ini.HasSection("Surrogate")) {
    params.gp_surrogate_options.potential =
        magic_enum::enum_cast<PotType>(ini.Get("Surrogate", "potential", ""),
                                       magic_enum::case_insensitive)
            .value_or(PotType::UNKNOWN);
    if (params.gp_surrogate_options.potential != PotType::CatLearn) {
      throw std::runtime_error("We only support catlearn for GP right now");
    }
  }
  // [CatLearn]
  if (ini.HasSection("CatLearn")) {
    params.catlearn_options.path = ini.Get("CatLearn", "catl_path", "");
    params.catlearn_options.model = ini.Get("CatLearn", "model", "catl_model");
    params.catlearn_options.prior = ini.Get("CatLearn", "prior", "catl_prior");
    params.catlearn_options.use_deriv =
        ini.GetBoolean("CatLearn", "use_derivatives", "catl_deriv");
    params.catlearn_options.use_fingerprint =
        ini.GetBoolean("CatLearn", "use_fingerprint", "catl_fingerprint");
    params.catlearn_options.parallel = ini.GetBoolean(
        "CatLearn", "parallel_hyperparameter_opt", "catl_parallel");
  }
  // [ASE_ORCA]
  if (ini.HasSection("ASE_ORCA")) {
    params.ase_orca_options.path = ini.Get("ASE_ORCA", "orca_path", "");
    params.ase_orca_options.nproc = ini.Get("ASE_ORCA", "nproc", "1");
    params.ase_orca_options.simpleinput =
        ini.Get("ASE_ORCA", "simpleinput", "");
  }
  // [ASE_NWCHEM]
  if (ini.HasSection("ASE_NWCHEM")) {
    params.ase_nwchem_options.path = ini.Get("ASE_NWCHEM", "nwchem_path", "");
    params.ase_nwchem_options.nproc = ini.Get("ASE_NWCHEM", "nproc", "1");
    params.ase_nwchem_options.mpi_launcher = ini.Get(
        "ASE_NWCHEM", "mpi_launcher", params.ase_nwchem_options.mpi_launcher);
    params.ase_nwchem_options.multiplicity =
        ini.Get("ASE_NWCHEM", "multiplicity", "");
    params.ase_nwchem_options.scf_thresh =
        ini.GetReal("ASE_NWCHEM", "scf_thresh", 1e-5);
    params.ase_nwchem_options.scf_maxiter =
        ini.GetInteger("ASE_NWCHEM", "scf_maxiter", 200);
  }
  // [Metatomic]
  if (ini.HasSection("Metatomic")) {
    params.metatomic_options.model_path =
        ini.Get("Metatomic", "model_path", "");
    params.metatomic_options.device = ini.Get("Metatomic", "device", "cpu");
    params.metatomic_options.length_unit =
        ini.Get("Metatomic", "length_unit", "angstrom");
    params.metatomic_options.extensions_directory =
        ini.Get("Metatomic", "extensions_directory", "");
    params.metatomic_options.check_consistency =
        ini.GetBoolean("Metatomic", "check_consistency", false);
    params.metatomic_options.uncertainty_threshold =
        ini.GetReal("Metatomic", "uncertainty_threshold", -1.0);
    params.metatomic_options.energy_output = ini.Get(
        "Metatomic", "energy_output", params.metatomic_options.energy_output);
    params.metatomic_options.energy_uncertainty_output =
        ini.Get("Metatomic", "energy_uncertainty_output",
                params.metatomic_options.energy_uncertainty_output);
    auto &_variant = params.metatomic_options.variant;
    _variant.base = ini.Get("Metatomic", "variant_base", "");
    _variant.energy = ini.Get("Metatomic", "variant_energy", "");
    _variant.energy_uncertainty =
        ini.Get("Metatomic", "variant_energy_uncertainty", "");
  }
  // [Serve]
  if (ini.HasSection("Serve")) {
    params.serve_options.host =
        ini.Get("Serve", "host", params.serve_options.host);
    params.serve_options.port = static_cast<uint16_t>(
        ini.GetInteger("Serve", "port", params.serve_options.port));
    params.serve_options.replicas = static_cast<size_t>(
        ini.GetInteger("Serve", "replicas", params.serve_options.replicas));
    params.serve_options.gateway_port = static_cast<uint16_t>(ini.GetInteger(
        "Serve", "gateway_port", params.serve_options.gateway_port));
    params.serve_options.endpoints =
        ini.Get("Serve", "endpoints", params.serve_options.endpoints);
  }

  // GP_NEB only
  params.gp_surrogate_options.linear_path_always =
      ini.GetBoolean("Surrogate", "gp_linear_path_always",
                     params.gp_surrogate_options.linear_path_always);
  // [Lanczos] //

  params.lanczos_options.tolerance =
      ini.GetReal("Lanczos", "tolerance", params.lanczos_options.tolerance);
  params.lanczos_options.max_iterations = ini.GetInteger(
      "Lanczos", "max_iterations", params.lanczos_options.max_iterations);
  params.lanczos_options.quit_early = ini.GetBoolean(
      "Lanczos", "quit_early", params.lanczos_options.quit_early);

  // [Davidson] //
  params.davidson_options.tolerance =
      ini.GetReal("Davidson", "tolerance", params.davidson_options.tolerance);
  params.davidson_options.max_iterations = ini.GetInteger(
      "Davidson", "max_iterations", params.davidson_options.max_iterations);
  params.davidson_options.diagonal_preconditioner =
      ini.GetBoolean("Davidson", "diagonal_preconditioner",
                     params.davidson_options.diagonal_preconditioner);

  // [ARTn] //
  params.artn_options.push_step_size =
      ini.GetReal("ARTn", "push_step_size", params.artn_options.push_step_size);
  params.artn_options.force_threshold = ini.GetReal(
      "ARTn", "force_threshold", params.artn_options.force_threshold);
  params.artn_options.max_iterations = ini.GetInteger(
      "ARTn", "max_iterations", params.artn_options.max_iterations);
  params.artn_options.ninit =
      ini.GetInteger("ARTn", "ninit", params.artn_options.ninit);
  params.artn_options.nperp_limitation =
      ini.Get("ARTn", "nperp_limitation", params.artn_options.nperp_limitation);
  params.artn_options.lanczos_min_size = ini.GetInteger(
      "ARTn", "lanczos_min_size", params.artn_options.lanczos_min_size);
  params.artn_options.nsmooth =
      ini.GetInteger("ARTn", "nsmooth", params.artn_options.nsmooth);
  params.artn_options.filin =
      ini.Get("ARTn", "filin", params.artn_options.filin);

  // [IRA] //
  params.ira_options.distance_threshold = ini.GetReal(
      "IRA", "distance_threshold", params.ira_options.distance_threshold);
  params.ira_options.symmetry_threshold = ini.GetReal(
      "IRA", "symmetry_threshold", params.ira_options.symmetry_threshold);
  params.ira_options.use_pbc =
      ini.GetBoolean("IRA", "use_pbc", params.ira_options.use_pbc);

  // [GPR Dimer] //
  params.gpr_dimer_options.rotation_angle = ini.GetReal(
      "GPR Dimer", "finite_angle", params.gpr_dimer_options.rotation_angle);
  params.gpr_dimer_options.converged_angle = ini.GetReal(
      "GPR Dimer", "converged_angle", params.gpr_dimer_options.converged_angle);
  params.gpr_dimer_options.relax_conv_angle =
      ini.GetReal("GPR Dimer", "relaxation_converged_angle",
                  params.gpr_dimer_options.relax_conv_angle);
  params.gpr_dimer_options.init_rotations_max = static_cast<long>(
      ini.GetInteger("GPR Dimer", "max_initial_rotation_iterations",
                     params.gpr_dimer_options.init_rotations_max));
  params.gpr_dimer_options.relax_rotations_max = static_cast<long>(
      ini.GetInteger("GPR Dimer", "max_relaxation_rotation_iterations",
                     params.gpr_dimer_options.relax_rotations_max));
  params.gpr_dimer_options.divisor_t_dimer_gp = static_cast<long>(
      ini.GetInteger("GPR Dimer", "divisor_t_dimer",
                     params.gpr_dimer_options.divisor_t_dimer_gp));
  params.gpr_dimer_options.max_outer_iterations = static_cast<long>(
      ini.GetInteger("GPR Dimer", "max_outer_iterations",
                     params.gpr_dimer_options.max_outer_iterations));
  params.gpr_dimer_options.max_inner_iterations = static_cast<long>(
      ini.GetInteger("GPR Dimer", "max_inner_iterations",
                     params.gpr_dimer_options.max_inner_iterations));
  params.gpr_dimer_options.midpoint_max_disp =
      ini.GetReal("GPR Dimer", "max_midpoint_displacement",
                  params.gpr_dimer_options.midpoint_max_disp);
  params.gpr_dimer_options.rot_opt_method =
      ini.Get("GPR Dimer", "rotation_opt_method",
              params.gpr_dimer_options.rot_opt_method);
  params.gpr_dimer_options.trans_opt_method =
      ini.Get("GPR Dimer", "translation_opt_method",
              params.gpr_dimer_options.trans_opt_method);
  params.gpr_dimer_options.active_radius = ini.GetReal(
      "GPR Dimer", "active_radius", params.gpr_dimer_options.active_radius);
  params.gpr_dimer_options.dimer_sep = ini.GetReal(
      "GPR Dimer", "dimer_separation", params.gpr_dimer_options.dimer_sep);
  params.gpr_dimer_options.conv_step =
      ini.GetReal("GPR Dimer", "convex_region_step_size",
                  params.gpr_dimer_options.conv_step);
  params.gpr_dimer_options.max_step = ini.GetReal(
      "GPR Dimer", "max_step_size", params.gpr_dimer_options.max_step);
  params.gpr_dimer_options.ratio_at_limit = ini.GetReal(
      "GPR Dimer", "ratio_at_limit", params.gpr_dimer_options.ratio_at_limit);
  params.gpr_dimer_options.init_rot_gp =
      ini.GetBoolean("GPR Dimer", "nogp_initial_rotations",
                     params.gpr_dimer_options.init_rot_gp);
  params.gpr_dimer_options.init_trans_gp =
      ini.GetBoolean("GPR Dimer", "nogp_init_translations",
                     params.gpr_dimer_options.init_trans_gp);
  params.gpr_dimer_options.many_iterations =
      ini.GetBoolean("GPR Dimer", "has_many_iterations",
                     params.gpr_dimer_options.many_iterations);
  // GPR Params
  params.gpr_dimer_options.gpr_params.hyper_opt_method =
      ini.Get("GPR Dimer", "hyperparameter_opt_method",
              params.gpr_dimer_options.gpr_params.hyper_opt_method);
  params.gpr_dimer_options.gpr_params.sigma2 = ini.GetReal(
      "GPR Dimer", "gpr_variance", params.gpr_dimer_options.gpr_params.sigma2);
  params.gpr_dimer_options.gpr_params.jitter_sigma2 =
      ini.GetReal("GPR Dimer", "gpr_jitter_variance",
                  params.gpr_dimer_options.gpr_params.jitter_sigma2);
  params.gpr_dimer_options.gpr_params.noise_sigma2 =
      ini.GetReal("GPR Dimer", "gpr_noise_variance",
                  params.gpr_dimer_options.gpr_params.noise_sigma2);
  params.gpr_dimer_options.gpr_params.prior_mu = ini.GetReal(
      "GPR Dimer", "prior_mean", params.gpr_dimer_options.gpr_params.prior_mu);
  params.gpr_dimer_options.gpr_params.prior_sigma2 =
      ini.GetReal("GPR Dimer", "prior_variance",
                  params.gpr_dimer_options.gpr_params.prior_sigma2);
  params.gpr_dimer_options.gpr_params.prior_nu =
      ini.GetReal("GPR Dimer", "prior_degrees_of_freedom",
                  params.gpr_dimer_options.gpr_params.prior_nu);
  // GPR Optimization Parameters
  params.gpr_dimer_options.opt_params.check_derivatives =
      ini.GetBoolean("GPR Dimer", "check_derivatives",
                     params.gpr_dimer_options.opt_params.check_derivatives);
  params.gpr_dimer_options.opt_params.max_iterations = static_cast<int>(
      ini.GetInteger("GPR Dimer", "opt_max_iterations",
                     params.gpr_dimer_options.opt_params.max_iterations));
  params.gpr_dimer_options.opt_params.tol_func =
      ini.GetReal("GPR Dimer", "opt_tol_func",
                  params.gpr_dimer_options.opt_params.tol_func);
  params.gpr_dimer_options.opt_params.tol_sol = ini.GetReal(
      "GPR Dimer", "opt_tol_sol", params.gpr_dimer_options.opt_params.tol_sol);
  params.gpr_dimer_options.opt_params.lambda_limit =
      ini.GetReal("GPR Dimer", "opt_lambda_limit",
                  params.gpr_dimer_options.opt_params.lambda_limit);
  params.gpr_dimer_options.opt_params.lambda_init =
      ini.GetReal("GPR Dimer", "opt_lambda_init",
                  params.gpr_dimer_options.opt_params.lambda_init);
  // GPR Debugging Parameters
  params.gpr_dimer_options.debug_params.report_level = static_cast<int>(
      ini.GetInteger("GPR Dimer", "report_level",
                     params.gpr_dimer_options.debug_params.report_level));
  params.gpr_dimer_options.debug_params.debug_level = static_cast<int>(
      ini.GetInteger("GPR Dimer", "debug_level",
                     params.gpr_dimer_options.debug_params.debug_level));
  params.gpr_dimer_options.debug_params.out_dir =
      ini.Get("GPR Dimer", "debug_output_directory",
              params.gpr_dimer_options.debug_params.out_dir);
  params.gpr_dimer_options.debug_params.pos_file =
      ini.Get("GPR Dimer", "debug_position_basename",
              params.gpr_dimer_options.debug_params.pos_file);
  params.gpr_dimer_options.debug_params.energy_file =
      ini.Get("GPR Dimer", "debug_energy_basename",
              params.gpr_dimer_options.debug_params.energy_file);
  params.gpr_dimer_options.debug_params.grad_file =
      ini.Get("GPR Dimer", "debug_gradient_basename",
              params.gpr_dimer_options.debug_params.grad_file);
  params.gpr_dimer_options.debug_params.offset_mid_point =
      ini.GetReal("GPR Dimer", "debug_midpoint_offset",
                  params.gpr_dimer_options.debug_params.offset_mid_point);
  params.gpr_dimer_options.debug_params.dy = ini.GetReal(
      "GPR Dimer", "debug_y_step", params.gpr_dimer_options.debug_params.dy);
  params.gpr_dimer_options.debug_params.dz = ini.GetReal(
      "GPR Dimer", "debug_z_step", params.gpr_dimer_options.debug_params.dz);
  // GPR Prune
  params.gpr_dimer_options.prune_params.use_prune =
      ini.GetBoolean("GPR Dimer", "use_prune",
                     params.gpr_dimer_options.prune_params.use_prune);
  params.gpr_dimer_options.prune_params.begin = static_cast<int>(
      ini.GetInteger("GPR Dimer", "start_prune_at",
                     params.gpr_dimer_options.prune_params.begin));
  params.gpr_dimer_options.prune_params.n_vals = static_cast<int>(
      ini.GetInteger("GPR Dimer", "nprune_vals",
                     params.gpr_dimer_options.prune_params.n_vals));
  params.gpr_dimer_options.prune_params.threshold =
      ini.GetReal("GPR Dimer", "prune_threshold",
                  params.gpr_dimer_options.prune_params.threshold);

  // [Prefactor] //

  params.prefactor_options.default_value = ini.GetReal(
      "Prefactor", "default_value", params.prefactor_options.default_value);
  params.prefactor_options.max_value =
      ini.GetReal("Prefactor", "max_value", params.prefactor_options.max_value);
  params.prefactor_options.min_value =
      ini.GetReal("Prefactor", "min_value", params.prefactor_options.min_value);
  params.prefactor_options.within_radius = ini.GetReal(
      "Prefactor", "within_radius", params.prefactor_options.within_radius);
  params.prefactor_options.min_displacement =
      ini.GetReal("Prefactor", "min_displacement",
                  params.prefactor_options.min_displacement);
  params.prefactor_options.rate = toLowerCase(
      ini.Get("Prefactor", "rate_estimation", params.prefactor_options.rate));
  params.prefactor_options.configuration = toLowerCase(ini.Get(
      "Prefactor", "configuration", params.prefactor_options.configuration));
  params.prefactor_options.all_free_atoms = ini.GetBoolean(
      "Prefactor", "all_free_atoms", params.prefactor_options.all_free_atoms);
  params.prefactor_options.filter_scheme = toLowerCase(ini.Get(
      "Prefactor", "filter_scheme", params.prefactor_options.filter_scheme));
  params.prefactor_options.filter_fraction = ini.GetReal(
      "Prefactor", "filter_fraction", params.prefactor_options.filter_fraction);

  // [Hessian] //

  params.hessian_options.atom_list = toLowerCase(
      ini.Get("Hessian", "atom_list", params.hessian_options.atom_list));
  params.hessian_options.zero_freq_value = ini.GetReal(
      "Hessian", "zero_freq_value", params.hessian_options.zero_freq_value);

  // [Nudged Elastic Band] //
  const std::string neb_section = "Nudged Elastic Band";

  params.neb_options.image_count =
      ini.GetInteger(neb_section, "images", params.neb_options.image_count);
  params.neb_options.max_iterations = ini.GetInteger(
      neb_section, "max_iterations", params.optimizer_options.max_iterations);
  params.neb_options.force_tolerance = ini.GetReal(
      neb_section, "converged_force", params.optimizer_options.converged_force);
  auto neb_optMethod =
      magic_enum::enum_cast<OptType>(ini.Get(neb_section, "opt_method", "none"),
                                     magic_enum::case_insensitive)
          .value_or(OptType::Unknown);
  if (neb_optMethod != OptType::None) {
    params.neb_options.opt_method = neb_optMethod;
  }
  params.neb_options.mmf_peaks.enabled = ini.GetBoolean(
      neb_section, "setup_mmf_peaks", params.neb_options.mmf_peaks.enabled);
  params.neb_options.mmf_peaks.tolerance =
      ini.GetReal(neb_section, "mmf_peak_tolerance",
                  params.neb_options.mmf_peaks.tolerance);

  params.neb_options.spring.constant =
      ini.GetReal(neb_section, "spring", params.neb_options.spring.constant);
  params.neb_options.spring.use_elastic_band = ini.GetBoolean(
      neb_section, "elastic_band", params.neb_options.spring.use_elastic_band);
  params.neb_options.spring.doubly_nudged = ini.GetBoolean(
      neb_section, "doubly_nudged", params.neb_options.spring.doubly_nudged);
  params.neb_options.spring.use_switching =
      ini.GetBoolean(neb_section, "doubly_nudged_switching",
                     params.neb_options.spring.use_switching);

  params.neb_options.spring.weighting.enabled =
      ini.GetBoolean(neb_section, "energy_weighted",
                     params.neb_options.spring.weighting.enabled);
  params.neb_options.spring.weighting.trigger = ini.GetReal(
      neb_section, "ew_trigger", params.neb_options.spring.weighting.trigger);
  params.neb_options.spring.weighting.k_min = ini.GetReal(
      neb_section, "ew_ksp_min", params.neb_options.spring.weighting.k_min);
  params.neb_options.spring.weighting.k_max = ini.GetReal(
      neb_section, "ew_ksp_max", params.neb_options.spring.weighting.k_max);

  params.neb_options.spring.om.enabled = ini.GetBoolean(
      neb_section, "onsager_machlup", params.neb_options.spring.om.enabled);
  params.neb_options.spring.om.optimize_k = ini.GetBoolean(
      neb_section, "om_optimize_k", params.neb_options.spring.om.optimize_k);
  params.neb_options.spring.om.k_scale = ini.GetReal(
      neb_section, "om_k_scale", params.neb_options.spring.om.k_scale);
  params.neb_options.spring.om.k_min =
      ini.GetReal(neb_section, "om_k_min", params.neb_options.spring.om.k_min);
  params.neb_options.spring.om.k_max =
      ini.GetReal(neb_section, "om_k_max", params.neb_options.spring.om.k_max);

  params.neb_options.climbing_image.enabled =
      ini.GetBoolean(neb_section, "climbing_image_method",
                     params.neb_options.climbing_image.enabled);
  params.neb_options.climbing_image.converged_only =
      ini.GetBoolean(neb_section, "climbing_image_converged_only",
                     params.neb_options.climbing_image.converged_only);
  params.neb_options.climbing_image.use_old_tangent =
      ini.GetBoolean(neb_section, "old_tangent",
                     params.neb_options.climbing_image.use_old_tangent);
  params.neb_options.climbing_image.trigger_force = ini.GetReal(
      neb_section, "ci_after", params.neb_options.climbing_image.trigger_force);
  params.neb_options.climbing_image.trigger_factor =
      ini.GetReal(neb_section, "ci_after_rel",
                  params.neb_options.climbing_image.trigger_factor);

  auto &oci = params.neb_options.climbing_image.ocineb;
  oci.use_mmf = ini.GetBoolean(neb_section, "ci_mmf", oci.use_mmf);
  oci.trigger_force =
      ini.GetReal(neb_section, "ci_mmf_after", oci.trigger_force);
  oci.trigger_factor =
      ini.GetReal(neb_section, "ci_mmf_after_rel", oci.trigger_factor);
  oci.max_steps = ini.GetInteger(neb_section, "ci_mmf_nsteps", oci.max_steps);
  oci.ci_stability_count = ini.GetInteger(
      neb_section, "ci_mmf_ci_stability_count", oci.ci_stability_count);
  oci.angle_tol = ini.GetReal(neb_section, "ci_mmf_angle", oci.angle_tol);

  auto &init = params.neb_options.initialization;
  init.method =
      magic_enum::enum_cast<NEBInit>(ini.Get(neb_section, "initializer", ""),
                                     magic_enum::case_insensitive)
          .value_or(NEBInit::LINEAR);
  init.input_path = ini.Get(neb_section, "initial_path_in", init.input_path);
  init.max_iterations =
      ini.GetInteger(neb_section, "init_max_iterations", init.max_iterations);
  init.nsteps = ini.GetInteger(neb_section, "init_nsteps", init.nsteps);
  init.max_move = ini.GetReal(neb_section, "init_max_move", init.max_move);
  init.force_tolerance =
      ini.GetReal(neb_section, "init_force_threshold", init.force_tolerance);
  init.sidpp_alpha =
      ini.GetReal(neb_section, "sidpp_growth_alpha", init.sidpp_alpha);
  init.sidpp_frontier_tol =
      ini.GetReal(neb_section, "sidpp_frontier_tol", init.sidpp_frontier_tol);
  init.sidpp_reparam =
      ini.GetBoolean(neb_section, "sidpp_reparameterize", init.sidpp_reparam);
  init.sidpp_ideal_ksp =
      ini.GetBoolean(neb_section, "sidpp_ideal_ksp", init.sidpp_ideal_ksp);
  auto neb_ipath_optMethod =
      magic_enum::enum_cast<OptType>(
          ini.Get(neb_section, "ipath_opt_method", "none"),
          magic_enum::case_insensitive)
          .value_or(OptType::Unknown);
  if (neb_ipath_optMethod != OptType::None) {
    params.neb_options.initialization.opt_method = neb_ipath_optMethod;
  }
  init.oversampling =
      ini.GetBoolean(neb_section, "oversampling", init.oversampling);
  init.oversampling_factor = ini.GetInteger(neb_section, "oversampling_factor",
                                            init.oversampling_factor);

  params.neb_options.endpoints.minimize = ini.GetBoolean(
      neb_section, "minimize_endpoints", params.neb_options.endpoints.minimize);
  params.neb_options.endpoints.use_path_file =
      ini.GetBoolean(neb_section, "minimize_endpoints_for_ipath",
                     params.neb_options.endpoints.use_path_file);

  // [Dynamics] //

  params.dynamics_options.time_step_input = ini.GetReal(
      "Dynamics", "time_step", params.dynamics_options.time_step_input);
  params.dynamics_options.time_step =
      params.dynamics_options.time_step_input / params.constants.timeUnit;
  params.dynamics_options.time_input =
      ini.GetReal("Dynamics", "time", params.dynamics_options.time_input);
  params.dynamics_options.time =
      params.dynamics_options.time_input / params.constants.timeUnit;
  params.dynamics_options.steps = static_cast<long>(std::floor(
      params.dynamics_options.time / params.dynamics_options.time_step + 0.5));
  params.thermostat_options.kind =
      toLowerCase(ini.Get("Dynamics", "thermostat", "andersen"));
  params.thermostat_options.andersen_alpha = ini.GetReal(
      "Dynamics", "andersen_alpha", params.thermostat_options.andersen_alpha);
  params.thermostat_options.andersen_tcol_input =
      ini.GetReal("Dynamics", "andersen_collision_period",
                  params.thermostat_options.andersen_tcol_input);
  params.thermostat_options.andersen_tcol =
      params.thermostat_options.andersen_tcol_input / params.constants.timeUnit;
  params.thermostat_options.nose_mass =
      ini.GetReal("Dynamics", "nose_mass", params.thermostat_options.nose_mass);
  params.thermostat_options.langevin_friction_input =
      ini.GetReal("Dynamics", "langevin_friction",
                  params.thermostat_options.langevin_friction_input);
  params.thermostat_options.langevin_friction =
      params.thermostat_options.langevin_friction_input *
      params.constants.timeUnit;

  // [Parallel Replica]

  params.parallel_replica_options.auto_stop =
      ini.GetBoolean("Parallel Replica", "stop_after_transition",
                     params.parallel_replica_options.auto_stop);
  params.parallel_replica_options.refine_transition =
      ini.GetBoolean("Parallel Replica", "refine_transition",
                     params.parallel_replica_options.refine_transition);
  params.parallel_replica_options.dephase_loop_stop =
      ini.GetBoolean("Parallel Replica", "dephase_loop_stop",
                     params.parallel_replica_options.dephase_loop_stop);
  params.parallel_replica_options.dephase_time_input =
      ini.GetReal("Parallel Replica", "dephase_time",
                  params.parallel_replica_options.dephase_time_input);
  params.parallel_replica_options.dephase_time =
      params.parallel_replica_options.dephase_time_input /
      params.constants.timeUnit;
  params.parallel_replica_options.dephase_loop_max =
      ini.GetInteger("Parallel Replica", "dephase_loop_max",
                     params.parallel_replica_options.dephase_loop_max);
  params.parallel_replica_options.state_check_interval_input =
      ini.GetReal("Parallel Replica", "state_check_interval",
                  params.parallel_replica_options.state_check_interval_input);
  params.parallel_replica_options.state_check_interval =
      params.parallel_replica_options.state_check_interval_input /
      params.constants.timeUnit;
  params.parallel_replica_options.record_interval_input = ini.GetReal(
      "Parallel Replica", "state_save_interval",
      0.1 * params.parallel_replica_options.state_check_interval_input);
  params.parallel_replica_options.record_interval =
      params.parallel_replica_options.record_interval_input /
      params.constants.timeUnit;
  params.parallel_replica_options.corr_time_input =
      ini.GetReal("Parallel Replica", "post_transition_time",
                  params.parallel_replica_options.corr_time_input);
  params.parallel_replica_options.corr_time =
      params.parallel_replica_options.corr_time_input /
      params.constants.timeUnit;

  // [Temperature Accelerated Dynamics] //

  params.tad_options.low_temperature =
      ini.GetReal("TAD", "low_temperature", params.tad_options.low_temperature);
  params.tad_options.min_prefactor =
      ini.GetReal("TAD", "min_prefactor", params.tad_options.min_prefactor);
  params.tad_options.confidence =
      ini.GetReal("TAD", "confidence", params.tad_options.confidence);

  // [Replica Exchange] //

  params.replica_exchange_options.temperature_distribution = toLowerCase(
      ini.Get("Replica Exchange", "temperature_distribution",
              params.replica_exchange_options.temperature_distribution));
  params.replica_exchange_options.replicas = ini.GetInteger(
      "Replica Exchange", "replicas", params.replica_exchange_options.replicas);
  params.replica_exchange_options.exchange_trials =
      ini.GetInteger("Replica Exchange", "exchange_trials",
                     params.replica_exchange_options.exchange_trials);
  params.replica_exchange_options.sampling_time_input =
      ini.GetReal("Replica Exchange", "sampling_time",
                  params.replica_exchange_options.sampling_time_input);
  params.replica_exchange_options.sampling_time =
      params.replica_exchange_options.sampling_time_input /
      params.constants.timeUnit;
  params.replica_exchange_options.temperature_low = ini.GetReal(
      "Replica Exchange", "temperature_low", params.main_options.temperature);
  params.replica_exchange_options.temperature_high =
      ini.GetReal("Replica Exchange", "temperature_high",
                  params.replica_exchange_options.temperature_high);
  params.replica_exchange_options.exchange_period_input =
      ini.GetReal("Replica Exchange", "exchange_period",
                  params.replica_exchange_options.exchange_period_input);
  params.replica_exchange_options.exchange_period =
      params.replica_exchange_options.exchange_period_input /
      params.constants.timeUnit;

  // [Hyperdynamics] //

  params.hyperdynamics_options.rmd_time_input =
      ini.GetReal("Hyperdynamics", "bb_rmd_time",
                  params.hyperdynamics_options.rmd_time_input);
  params.hyperdynamics_options.rmd_time =
      params.hyperdynamics_options.rmd_time_input / params.constants.timeUnit;
  params.hyperdynamics_options.boost_atom_list =
      toLowerCase(ini.Get("Hyperdynamics", "bb_boost_atomlist",
                          params.hyperdynamics_options.boost_atom_list));
  params.hyperdynamics_options.dvmax = ini.GetReal(
      "Hyperdynamics", "bb_dvmax", params.hyperdynamics_options.dvmax);
  params.hyperdynamics_options.qrr =
      ini.GetReal("Hyperdynamics", "bb_stretch_threshold",
                  params.hyperdynamics_options.qrr);
  params.hyperdynamics_options.prr = ini.GetReal(
      "Hyperdynamics", "bb_ds_curvature", params.hyperdynamics_options.prr);
  params.hyperdynamics_options.qcut = ini.GetReal(
      "Hyperdynamics", "bb_rcut", params.hyperdynamics_options.qcut);
  params.hyperdynamics_options.bias_potential =
      toLowerCase(ini.Get("Hyperdynamics", "bias_potential",
                          params.hyperdynamics_options.bias_potential));

  // [Saddle Search] //

  params.saddle_search_options.method = toLowerCase(
      ini.Get("Saddle Search", "method", params.saddle_search_options.method));
  params.saddle_search_options.minmode_method =
      toLowerCase(ini.Get("Saddle Search", "min_mode_method",
                          params.saddle_search_options.minmode_method));
  params.saddle_search_options.displace_magnitude =
      ini.GetReal("Saddle Search", "displace_magnitude",
                  params.saddle_search_options.displace_magnitude);
  params.saddle_search_options.displace_radius =
      ini.GetReal("Saddle Search", "displace_radius",
                  params.saddle_search_options.displace_radius);
  params.saddle_search_options.max_energy = ini.GetReal(
      "Saddle Search", "max_energy", params.saddle_search_options.max_energy);
  params.saddle_search_options.max_iterations =
      ini.GetInteger("Saddle Search", "max_iterations",
                     params.optimizer_options.max_iterations);
  params.saddle_search_options.nonnegative_displacement_abort = ini.GetBoolean(
      "Saddle Search", "nonnegative_displacement_abort",
      params.saddle_search_options.nonnegative_displacement_abort);
  params.saddle_search_options.max_single_displace =
      ini.GetReal("Saddle Search", "max_single_displace",
                  params.saddle_search_options.max_single_displace);
  params.saddle_search_options.converged_force =
      ini.GetReal("Saddle Search", "converged_force",
                  params.optimizer_options.converged_force);
  params.saddle_search_options.perp_force_ratio =
      ini.GetReal("Saddle Search", "perp_force_ratio",
                  params.saddle_search_options.perp_force_ratio);
  params.saddle_search_options.displace_type = toLowerCase(ini.Get(
      "Saddle Search", "client_displace_type", eonc::EpiCenters::DISP_LOAD));
  params.saddle_search_options.nonlocal_count_abort =
      ini.GetInteger("Saddle Search", "nonlocal_count_abort",
                     params.saddle_search_options.nonlocal_count_abort);
  params.saddle_search_options.nonlocal_distance_abort =
      ini.GetReal("Saddle Search", "nonlocal_distance_abort",
                  params.saddle_search_options.nonlocal_distance_abort);
  if (params.saddle_search_options.displace_type !=
          eonc::EpiCenters::DISP_NOT_FCC_OR_HCP &&
      params.saddle_search_options.displace_type !=
          eonc::EpiCenters::DISP_MIN_COORDINATED &&
      params.saddle_search_options.displace_type !=
          eonc::EpiCenters::DISP_LAST_ATOM &&
      params.saddle_search_options.displace_type !=
          eonc::EpiCenters::DISP_RANDOM &&
      params.saddle_search_options.displace_type !=
          eonc::EpiCenters::DISP_LISTED_ATOMS) {
    params.saddle_search_options.displace_type = eonc::EpiCenters::DISP_LOAD;
  }
  // Parse comma-separated atom list
  {
    std::string atomListStr =
        ini.Get("Saddle Search", "displace_atom_list", "");
    if (!atomListStr.empty()) {
      std::stringstream ss(atomListStr);
      std::string token;
      while (std::getline(ss, token, ',')) {
        size_t start = token.find_first_not_of(" \t");
        size_t end = token.find_last_not_of(" \t");
        if (start != std::string::npos) {
          params.saddle_search_options.displace_atom_list.push_back(
              std::stol(token.substr(start, end - start + 1)));
        }
      }
    }
  }
  params.saddle_search_options.confine_positive.enabled =
      ini.GetBoolean("Saddle Search", "confine_positive",
                     params.saddle_search_options.confine_positive.enabled);
  if (params.saddle_search_options.confine_positive.enabled) {
    params.saddle_search_options.confine_positive.bowl_breakout =
        ini.GetBoolean(
            "Saddle Search", "bowl_breakout",
            params.saddle_search_options.confine_positive.bowl_breakout);
    params.saddle_search_options.confine_positive.bowl_active = ini.GetInteger(
        "Saddle Search", "bowl_active_atoms",
        params.saddle_search_options.confine_positive.bowl_active);
    params.saddle_search_options.confine_positive.min_force =
        ini.GetReal("Saddle Search", "confine_positive_min_move",
                    params.saddle_search_options.confine_positive.min_force);
    params.saddle_search_options.confine_positive.scale_ratio =
        ini.GetReal("Saddle Search", "confine_positive_scale_ratio",
                    params.saddle_search_options.confine_positive.scale_ratio);
    params.saddle_search_options.confine_positive.boost =
        ini.GetReal("Saddle Search", "confine_positive_boost",
                    params.saddle_search_options.confine_positive.boost);
    params.saddle_search_options.confine_positive.min_active = ini.GetInteger(
        "Saddle Search", "confine_positive_min_active",
        params.saddle_search_options.confine_positive.min_active);
  }
  params.saddle_search_options.dynamics.temperature =
      params.main_options.temperature;
  params.saddle_search_options.dynamics.temperature =
      ini.GetReal("Saddle Search", "dynamics_temperature",
                  params.saddle_search_options.dynamics.temperature);
  params.saddle_search_options.dynamics.state_check_interval_input =
      ini.GetReal(
          "Saddle Search", "dynamics_state_check_interval",
          params.saddle_search_options.dynamics.state_check_interval_input);
  params.saddle_search_options.dynamics.state_check_interval =
      params.saddle_search_options.dynamics.state_check_interval_input /
      params.constants.timeUnit;
  params.saddle_search_options.dynamics.record_interval_input =
      ini.GetReal("Saddle Search", "dynamics_record_interval",
                  params.saddle_search_options.dynamics.record_interval_input);
  params.saddle_search_options.dynamics.record_interval =
      params.saddle_search_options.dynamics.record_interval_input /
      params.constants.timeUnit;
  params.saddle_search_options.dynamics.linear_interpolation = ini.GetBoolean(
      "Saddle Search", "dynamics_linear_interpolation",
      params.saddle_search_options.dynamics.linear_interpolation);
  params.saddle_search_options.remove_rotation =
      ini.GetBoolean("Saddle Search", "remove_rotation",
                     params.saddle_search_options.remove_rotation);
  params.saddle_search_options.dynamics.max_init_curvature =
      ini.GetReal("Saddle Search", "dynamics_max_init_curvature",
                  params.saddle_search_options.dynamics.max_init_curvature);
  params.saddle_search_options.zero_mode_abort_curvature =
      ini.GetReal("Saddle Search", "zero_mode_abort_curvature",
                  params.saddle_search_options.zero_mode_abort_curvature);

  // [Basin Hopping] //

  params.basin_hopping_options.displacement =
      ini.GetReal("Basin Hopping", "displacement",
                  params.basin_hopping_options.displacement);
  params.basin_hopping_options.push_apart_distance =
      ini.GetReal("Basin Hopping", "push_apart_distance",
                  params.basin_hopping_options.push_apart_distance);
  params.basin_hopping_options.initial_random_structure_probability =
      ini.GetReal(
          "Basin Hopping", "initial_random_structure_probability",
          params.basin_hopping_options.initial_random_structure_probability);
  params.basin_hopping_options.steps = ini.GetInteger(
      "Basin Hopping", "steps", params.basin_hopping_options.steps);
  params.basin_hopping_options.quenching_steps =
      ini.GetInteger("Basin Hopping", "quenching_steps",
                     params.basin_hopping_options.quenching_steps);
  params.basin_hopping_options.single_atom_displace =
      ini.GetBoolean("Basin Hopping", "single_atom_displace",
                     params.basin_hopping_options.single_atom_displace);
  params.basin_hopping_options.significant_structure =
      ini.GetBoolean("Basin Hopping", "significant_structure",
                     params.basin_hopping_options.significant_structure);
  params.basin_hopping_options.displacement_algorithm =
      toLowerCase(ini.Get("Basin Hopping", "displacement_algorithm",
                          params.basin_hopping_options.displacement_algorithm));
  params.basin_hopping_options.displacement_distribution = toLowerCase(
      ini.Get("Basin Hopping", "displacement_distribution",
              params.basin_hopping_options.displacement_distribution));
  params.basin_hopping_options.swap_probability =
      ini.GetReal("Basin Hopping", "swap_probability",
                  params.basin_hopping_options.swap_probability);
  params.basin_hopping_options.jump_max = ini.GetInteger(
      "Basin Hopping", "jump_max", params.basin_hopping_options.jump_max);
  params.basin_hopping_options.jump_steps = ini.GetInteger(
      "Basin Hopping", "jump_steps", params.basin_hopping_options.jump_steps);
  params.basin_hopping_options.adjust_displacement =
      ini.GetBoolean("Basin Hopping", "adjust_displacement",
                     params.basin_hopping_options.adjust_displacement);
  params.basin_hopping_options.adjust_period =
      ini.GetInteger("Basin Hopping", "adjust_period",
                     params.basin_hopping_options.adjust_period);
  params.basin_hopping_options.adjust_fraction =
      ini.GetReal("Basin Hopping", "adjust_fraction",
                  params.basin_hopping_options.adjust_fraction);
  params.basin_hopping_options.target_ratio =
      ini.GetReal("Basin Hopping", "target_ratio",
                  params.basin_hopping_options.target_ratio);
  params.basin_hopping_options.write_unique =
      ini.GetBoolean("Basin Hopping", "write_unique",
                     params.basin_hopping_options.write_unique);
  params.basin_hopping_options.stop_energy = ini.GetReal(
      "Basin Hopping", "stop_energy", params.basin_hopping_options.stop_energy);

  // [Global Optimization] //

  params.global_optimization_options.move_method =
      toLowerCase(ini.Get("Global Optimization", "move_method",
                          params.global_optimization_options.move_method));
  params.global_optimization_options.decision_method =
      toLowerCase(ini.Get("Global Optimization", "decision_method",
                          params.global_optimization_options.decision_method));
  params.global_optimization_options.steps = ini.GetInteger(
      "Global Optimization", "steps", params.global_optimization_options.steps);
  params.global_optimization_options.beta = ini.GetReal(
      "Global Optimization", "beta", params.global_optimization_options.beta);
  params.global_optimization_options.alpha = ini.GetReal(
      "Global Optimization", "alpha", params.global_optimization_options.alpha);
  params.global_optimization_options.mdmin = ini.GetInteger(
      "Global Optimization", "mdmin", params.global_optimization_options.mdmin);
  params.global_optimization_options.target_energy =
      ini.GetReal("Global Optimization", "target_energy",
                  params.global_optimization_options.target_energy);

  // [BGSD] //

  params.bgsd_options.alpha =
      ini.GetReal("BGSD", "alpha", params.bgsd_options.alpha);
  params.bgsd_options.beta =
      ini.GetReal("BGSD", "beta", params.bgsd_options.beta);
  params.bgsd_options.gradient_finite_difference =
      ini.GetReal("BGSD", "gradientfinitedifference",
                  params.bgsd_options.gradient_finite_difference);
  params.bgsd_options.grad2energy_convergence =
      ini.GetReal("BGSD", "grad2energyconvergence",
                  params.bgsd_options.grad2energy_convergence);
  params.bgsd_options.grad2force_convergence =
      ini.GetReal("BGSD", "grad2forceconvergence",
                  params.bgsd_options.grad2force_convergence);

  // [Monte Carlo] //

  params.monte_carlo_options.step_size = ini.GetReal(
      "Monte Carlo", "step_size", params.monte_carlo_options.step_size);
  params.monte_carlo_options.steps = static_cast<int>(
      ini.GetInteger("Monte Carlo", "steps", params.monte_carlo_options.steps));

  return error;
}

void validate_and_link(Parameters &params) {
  // Time unit conversions
  double tu = params.constants.timeUnit;
  params.optimizer_options.time_step =
      params.optimizer_options.time_step_input / tu;
  params.optimizer_options.max_time_step =
      params.optimizer_options.max_time_step_input / tu;

  params.dynamics_options.time_step =
      params.dynamics_options.time_step_input / tu;
  params.dynamics_options.time = params.dynamics_options.time_input / tu;
  params.dynamics_options.steps = static_cast<long>(std::floor(
      params.dynamics_options.time / params.dynamics_options.time_step + 0.5));

  params.thermostat_options.langevin_friction =
      params.thermostat_options.langevin_friction_input * tu;

  params.parallel_replica_options.dephase_time =
      params.parallel_replica_options.dephase_time_input / tu;
  params.parallel_replica_options.state_check_interval =
      params.parallel_replica_options.state_check_interval_input / tu;
  params.parallel_replica_options.record_interval =
      params.parallel_replica_options.record_interval_input / tu;
  params.parallel_replica_options.corr_time =
      params.parallel_replica_options.corr_time_input / tu;

  params.replica_exchange_options.sampling_time =
      params.replica_exchange_options.sampling_time_input / tu;
  params.replica_exchange_options.exchange_trials =
      params.replica_exchange_options.replicas;

  params.hyperdynamics_options.rmd_time =
      params.hyperdynamics_options.rmd_time_input / tu;

  params.saddle_search_options.dynamics.state_check_interval =
      params.saddle_search_options.dynamics.state_check_interval_input / tu;
  params.saddle_search_options.dynamics.record_interval =
      params.saddle_search_options.dynamics.record_interval_input / tu;

  // Cross-group defaults
  params.process_search_options.minimization_offset =
      params.optimizer_options.max_move;
  params.neb_options.force_tolerance = params.optimizer_options.converged_force;
}

} // namespace eonc::config
