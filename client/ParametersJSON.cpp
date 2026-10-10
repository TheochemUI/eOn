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
#include "eon/ParametersJSON.h"
#include "eon/HelperFunctions.h"
#include "eon/PIQTST.h"
#include "eon/Parameters.h"
#include "eon/ParametersINI.h"
#include "magic_enum/magic_enum.hpp"

#include <cctype>
#include <cstdint>
#include <format>
#include <nlohmann/json.hpp>
#include <sstream>
#include <stdexcept>
#include <vector>

using json = nlohmann::json;

// Helper: enum <-> string via magic_enum
template <typename E> static json enum_to_json(E val) {
  return std::string(magic_enum::enum_name(val));
}

template <typename E> static E enum_from_json(const json &j, E fallback) {
  if (j.is_string()) {
    auto result = magic_enum::enum_cast<E>(j.get<std::string>(),
                                           magic_enum::case_insensitive);
    return result.value_or(fallback);
  }
  return fallback;
}

// Macro for optional JSON field extraction with default
#define JSON_OPT(j, key, target)                                               \
  if ((j).contains(key))                                                       \
  (target) = (j).at(key).get<decltype(target)>()

namespace eonc::config {

#include "eon/generated/ParametersSSOTJson.inc"

json to_json(const Parameters &p) {
  json j;

  // [Main]
  j["Main"] = {
      {"job", enum_to_json(ParametersLoadAccess::main_options(p).job)},
      {"random_seed", ParametersLoadAccess::main_options(p).randomSeed},
      {"temperature", ParametersLoadAccess::main_options(p).temperature},
      {"quiet", ParametersLoadAccess::main_options(p).quiet},
      {"write_log", ParametersLoadAccess::main_options(p).writeLog},
      {"checkpoint", ParametersLoadAccess::main_options(p).checkpoint},
      {"ini_filename", ParametersLoadAccess::main_options(p).iniFilename},
      {"con_filename", ParametersLoadAccess::main_options(p).conFilename},
      {"finite_difference",
       ParametersLoadAccess::main_options(p).finiteDifference},
      {"max_force_calls", ParametersLoadAccess::main_options(p).maxForceCalls},
      {"remove_net_force",
       ParametersLoadAccess::main_options(p).removeNetForce},
      {"write_con_forces",
       ParametersLoadAccess::main_options(p).writeConForces},
  };

  // [Potential]
  j["Potential"] = {
      {"potential",
       enum_to_json(ParametersLoadAccess::potential_options(p).potential)},
      {"mpi_poll_period",
       ParametersLoadAccess::potential_options(p).MPIPollPeriod},
      {"lammps_logging",
       ParametersLoadAccess::potential_options(p).LAMMPSLogging},
      {"lammps_threads",
       ParametersLoadAccess::potential_options(p).LAMMPSThreads},
      {"emt_rasmussen",
       ParametersLoadAccess::potential_options(p).EMTRasmussen},
      {"log_potential",
       ParametersLoadAccess::potential_options(p).LogPotential},
      {"ext_pot_path", ParametersLoadAccess::potential_options(p).extPotPath},
      {"potentials_path",
       ParametersLoadAccess::potential_options(p).potentialsPath},
  };

  // [Structure Comparison]
  j["Structure Comparison"] = {
      {"distance_difference",
       ParametersLoadAccess::structure_comparison_options(p)
           .distance_difference},
      {"neighbor_cutoff",
       ParametersLoadAccess::structure_comparison_options(p).neighbor_cutoff},
      {"check_rotation",
       ParametersLoadAccess::structure_comparison_options(p).check_rotation},
      {"indistinguishable_atoms",
       ParametersLoadAccess::structure_comparison_options(p)
           .indistinguishable_atoms},
      {"energy_difference",
       ParametersLoadAccess::structure_comparison_options(p).energy_difference},
      {"remove_translation",
       ParametersLoadAccess::structure_comparison_options(p)
           .remove_translation},
  };

  // [Optimizer]
  j["Optimizer"] = {
      {"opt_method",
       enum_to_json(ParametersLoadAccess::optimizer_options(p).method)},
      {"convergence_metric",
       ParametersLoadAccess::optimizer_options(p).convergence_metric},
      {"max_iterations",
       ParametersLoadAccess::optimizer_options(p).max_iterations},
      {"max_move", ParametersLoadAccess::optimizer_options(p).max_move},
      {"converged_force",
       ParametersLoadAccess::optimizer_options(p).converged_force},
      {"time_step", ParametersLoadAccess::optimizer_options(p).time_step_input},
      {"max_time_step",
       ParametersLoadAccess::optimizer_options(p).max_time_step_input},
  };
  j["Optimizer"]["LBFGS"] = {
      {"memory", ParametersLoadAccess::optimizer_options(p).lbfgs.memory},
      {"inverse_curvature",
       ParametersLoadAccess::optimizer_options(p).lbfgs.inverse_curvature},
      {"max_inverse_curvature",
       ParametersLoadAccess::optimizer_options(p).lbfgs.max_inverse_curvature},
      {"auto_scale",
       ParametersLoadAccess::optimizer_options(p).lbfgs.auto_scale},
      {"angle_reset",
       ParametersLoadAccess::optimizer_options(p).lbfgs.angle_reset},
      {"distance_reset",
       ParametersLoadAccess::optimizer_options(p).lbfgs.distance_reset},
      {"curvature", ParametersLoadAccess::optimizer_options(p).lbfgs.curvature},
      {"project_rigid",
       ParametersLoadAccess::optimizer_options(p).lbfgs.project_rigid},
      {"secant", ParametersLoadAccess::optimizer_options(p).lbfgs.secant},
      {"precon", ParametersLoadAccess::optimizer_options(p).lbfgs.precon},
      {"step", ParametersLoadAccess::optimizer_options(p).lbfgs.step},
      {"h0", ParametersLoadAccess::optimizer_options(p).lbfgs.h0},
      {"accept", ParametersLoadAccess::optimizer_options(p).lbfgs.accept},
      {"extra_updates",
       ParametersLoadAccess::optimizer_options(p).lbfgs.extra_updates},
      {"cautious_eps",
       ParametersLoadAccess::optimizer_options(p).lbfgs.cautious_eps},
      {"cautious_alpha",
       ParametersLoadAccess::optimizer_options(p).lbfgs.cautious_alpha},
      {"precon_A", ParametersLoadAccess::optimizer_options(p).lbfgs.precon_A},
      {"precon_mu", ParametersLoadAccess::optimizer_options(p).lbfgs.precon_mu},
      {"precon_rcut",
       ParametersLoadAccess::optimizer_options(p).lbfgs.precon_rcut},
  };
  j["Optimizer"]["Xtsci"] = {
      {"method", ParametersLoadAccess::optimizer_options(p).xtsci.method},
      {"qn_step", ParametersLoadAccess::optimizer_options(p).xtsci.qn_step},
      {"precon", ParametersLoadAccess::optimizer_options(p).xtsci.precon},
      {"accept", ParametersLoadAccess::optimizer_options(p).xtsci.accept},
      {"highs", ParametersLoadAccess::optimizer_options(p).xtsci.highs},
      {"manifold", ParametersLoadAccess::optimizer_options(p).xtsci.manifold},
  };

  // [Dynamics]
  j["Dynamics"] = {
      {"time_step", ParametersLoadAccess::dynamics_options(p).time_step_input},
      {"time", ParametersLoadAccess::dynamics_options(p).time_input},
  };

  // [Thermostat]
  j["Thermostat"] = {
      {"kind", ParametersLoadAccess::thermostat_options(p).kind},
      {"andersen_alpha",
       ParametersLoadAccess::thermostat_options(p).andersen_alpha},
      {"andersen_collision_period",
       ParametersLoadAccess::thermostat_options(p).andersen_tcol_input},
      {"nose_mass", ParametersLoadAccess::thermostat_options(p).nose_mass},
      {"langevin_friction",
       ParametersLoadAccess::thermostat_options(p).langevin_friction_input},
      {"path_beads", ParametersLoadAccess::thermostat_options(p).path_beads},
      {"path_springs",
       ParametersLoadAccess::thermostat_options(p).path_springs},
      {"path_eco_omega_max",
       ParametersLoadAccess::thermostat_options(p).path_eco_omega_max},
      {"path_gle_file",
       ParametersLoadAccess::thermostat_options(p).path_gle_file},
      {"path_pile_tau",
       ParametersLoadAccess::thermostat_options(p).path_pile_tau_input},
      {"path_pile_scale",
       ParametersLoadAccess::thermostat_options(p).path_pile_scale},
      {"path_seed", ParametersLoadAccess::thermostat_options(p).path_seed},
  };

  // [Nudged Elastic Band]
  j["Nudged Elastic Band"] = {
      {"images", ParametersLoadAccess::neb_options(p).image_count},
      {"max_iterations", ParametersLoadAccess::neb_options(p).max_iterations},
      {"opt_method",
       enum_to_json(ParametersLoadAccess::neb_options(p).opt_method)},
      {"converged_force", ParametersLoadAccess::neb_options(p).force_tolerance},
      {"solid_state", ParametersLoadAccess::neb_options(p).solid_state.enabled},
      {"solid_state_weight",
       ParametersLoadAccess::neb_options(p).solid_state.weight},
      {"solid_state_pressure",
       ParametersLoadAccess::neb_options(p).solid_state.pressure},
      {"temperature", ParametersLoadAccess::neb_options(p).quantum_temperature},
  };
  j["Nudged Elastic Band"]["spring"] = {
      {"constant", ParametersLoadAccess::neb_options(p).spring.constant},
      {"elastic_band",
       ParametersLoadAccess::neb_options(p).spring.use_elastic_band},
      {"doubly_nudged",
       ParametersLoadAccess::neb_options(p).spring.doubly_nudged},
      {"geometric", ParametersLoadAccess::neb_options(p).spring.geometric},
  };
  j["Nudged Elastic Band"]["climbing_image"] = {
      {"enabled", ParametersLoadAccess::neb_options(p).climbing_image.enabled},
      {"converged_only",
       ParametersLoadAccess::neb_options(p).climbing_image.converged_only},
      {"band_slack",
       ParametersLoadAccess::neb_options(p).climbing_image.band_slack},
  };
  {
    const auto &zoom = ParametersLoadAccess::neb_options(p).zoom;
    j["Nudged Elastic Band"]["zoom"] = {
        {"enabled", zoom.enabled},
        {"alpha", zoom.alpha},
        {"offset", zoom.offset},
        {"mode", enum_to_json(zoom.mode)},
        {"activation_threshold", zoom.activation_threshold},
        {"interpolation", enum_to_json(zoom.interpolation)},
        {"stability_count", zoom.stability_count},
        {"max_iterations", zoom.max_iterations},
    };
  }

  // [Dimer]
  j["Dimer"] = {
      {"rotation_angle", ParametersLoadAccess::dimer_options(p).rotation_angle},
      {"improved", ParametersLoadAccess::dimer_options(p).improved},
      {"converged_angle",
       ParametersLoadAccess::dimer_options(p).converged_angle},
      {"max_iterations", ParametersLoadAccess::dimer_options(p).max_iterations},
      {"opt_method",
       enum_to_json(ParametersLoadAccess::dimer_options(p).opt_method)},
      {"rotation_backend",
       enum_to_json(ParametersLoadAccess::dimer_options(p).rotation_backend)},
      {"lor_residual_tol",
       ParametersLoadAccess::dimer_options(p).lor_residual_tol},
  };

  // [Saddle Search]
  j["Saddle Search"] = {
      {"method", ParametersLoadAccess::saddle_search_options(p).method},
      {"min_mode_method",
       ParametersLoadAccess::saddle_search_options(p).minmode_method},
      {"max_energy", ParametersLoadAccess::saddle_search_options(p).max_energy},
      {"max_iterations",
       ParametersLoadAccess::saddle_search_options(p).max_iterations},
      {"displace_magnitude",
       ParametersLoadAccess::saddle_search_options(p).displace_magnitude},
      {"displace_radius",
       ParametersLoadAccess::saddle_search_options(p).displace_radius},
  };

  // [Prefactor]
  j["Prefactor"] = {
      {"default_value",
       ParametersLoadAccess::prefactor_options(p).default_value},
      {"max_value", ParametersLoadAccess::prefactor_options(p).max_value},
      {"min_value", ParametersLoadAccess::prefactor_options(p).min_value},
  };

  // [Lanczos]
  j["Lanczos"] = {
      {"tolerance", ParametersLoadAccess::lanczos_options(p).tolerance},
      {"max_iterations",
       ParametersLoadAccess::lanczos_options(p).max_iterations},
      {"quit_early", ParametersLoadAccess::lanczos_options(p).quit_early},
      {"phva_atoms", ParametersLoadAccess::lanczos_options(p).phva_atoms},
  };

  // [Davidson]
  j["Davidson"] = {
      {"tolerance", ParametersLoadAccess::davidson_options(p).tolerance},
      {"max_iterations",
       ParametersLoadAccess::davidson_options(p).max_iterations},
      {"diagonal_preconditioner",
       ParametersLoadAccess::davidson_options(p).diagonal_preconditioner},
      {"phva_atoms", ParametersLoadAccess::davidson_options(p).phva_atoms},
  };

  // [Hessian]
  j["Hessian"] = {
      {"phva_atoms", ParametersLoadAccess::hessian_options(p).phva_atoms},
      {"zero_freq_value",
       ParametersLoadAccess::hessian_options(p).zero_freq_value},
      {"fd_scheme", ParametersLoadAccess::hessian_options(p).fd_scheme},
      {"resume", ParametersLoadAccess::hessian_options(p).resume},
      {"checkpoint_path",
       ParametersLoadAccess::hessian_options(p).checkpoint_path},
      {"write_modes", ParametersLoadAccess::hessian_options(p).write_modes},
  };

  // [Instanton]
  {
    const auto &o = ParametersLoadAccess::instanton_options(p);
    j["Instanton"] = {
        {"mode", o.mode},
        {"reactant_filename", o.reactant_filename},
        {"product_filename", o.product_filename},
        {"initial_path", o.initial_path},
        {"beads", o.beads},
        {"beta_hbar_omega", o.beta_hbar_omega},
        {"max_iterations", o.max_iterations},
        {"force_tolerance", o.force_tolerance},
        {"hessian_stride", o.hessian_stride},
        {"saddle_filename", o.saddle_filename},
        {"temperature", o.temperature},
        {"temperatures", o.temperatures},
        {"half_ring", o.half_ring},
        {"symmetrize", o.symmetrize},
        {"initial_hessians", o.initial_hessians},
        {"energy_shift", o.energy_shift},
        {"discretization", o.discretization},
        {"friction", o.friction},
        {"friction_eta", o.friction_eta},
        {"friction_eta_beads", o.friction_eta_beads},
        {"bead_ladder", o.bead_ladder},
        {"hessian_final", o.hessian_final},
        {"springs", o.springs},
        {"pi_planes", o.pi_planes},
        {"pi_beads", o.pi_beads},
        {"pi_equilibration_steps", o.pi_equilibration_steps},
        {"pi_sampling_steps", o.pi_sampling_steps},
        {"pi_time_step", o.pi_time_step},
        {"pi_thermostat", o.pi_thermostat},
        {"pi_gle_file", o.pi_gle_file},
        {"pi_pile_tau", o.pi_pile_tau},
        {"pi_pile_scale", o.pi_pile_scale},
        {"pi_seed", o.pi_seed},
        {"pi_direction", o.pi_direction},
        {"pi_reactant_extent", o.pi_reactant_extent},
        {"pi_recrossing_parents", o.pi_recrossing_parents},
        {"pi_recrossing_children", o.pi_recrossing_children},
        {"pi_recrossing_time", o.pi_recrossing_time},
        {"pi_recrossing_spacing", o.pi_recrossing_spacing},
    };
  }

  // [Debug]
  j["Debug"] = {
      {"write_movies", ParametersLoadAccess::debug_options(p).write_movies},
      {"write_movies_interval",
       ParametersLoadAccess::debug_options(p).write_movies_interval},
      {"write_deprecated_outs",
       ParametersLoadAccess::debug_options(p).write_deprecated_outs},
  };

  // [Serve]
  j["Serve"] = {
      {"host", ParametersLoadAccess::serve_options(p).host},
      {"port", ParametersLoadAccess::serve_options(p).port},
      {"replicas", ParametersLoadAccess::serve_options(p).replicas},
      {"gateway_port", ParametersLoadAccess::serve_options(p).gateway_port},
      {"endpoints", ParametersLoadAccess::serve_options(p).endpoints},
  };

  project_ssot_json_write(j, p);
  return j;
}

static std::string lowerCopy(std::string value) {
  for (char &ch : value) {
    ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  }
  return value;
}

static std::vector<double> temperaturesFromJson(const json &value) {
  if (value.is_array()) {
    return value.get<std::vector<double>>();
  }
  if (!value.is_string()) {
    throw std::invalid_argument(
        "[Instanton] temperatures must be a list or a comma-separated "
        "string");
  }
  std::vector<double> out;
  std::stringstream ss(value.get<std::string>());
  std::string token;
  while (std::getline(ss, token, ',')) {
    const size_t start = token.find_first_not_of(" \t");
    const size_t end = token.find_last_not_of(" \t");
    if (start == std::string::npos) {
      continue;
    }
    try {
      out.push_back(std::stod(token.substr(start, end - start + 1)));
    } catch (const std::exception &) {
      throw std::invalid_argument(
          "[Instanton] temperatures must be comma-separated kelvin "
          "values, not " +
          token);
    }
  }
  return out;
}

// Path-integral keys live on [Dynamics] in an ini file and under Thermostat
// in the JSON writer. Either object may carry them. The damping time is
// stored in femtoseconds; from_json converts it once every object is read.
static void readPathIntegralKeys(const json &s, Parameters &p) {
  auto &th = ParametersLoadAccess::thermostat_options(p);
  JSON_OPT(s, "path_beads", th.path_beads);
  if (s.contains("path_springs")) {
    th.path_springs = lowerCopy(s.at("path_springs").get<std::string>());
    if (th.path_springs != "trotter" && th.path_springs != "eco") {
      throw std::invalid_argument(
          "[Dynamics] path_springs must be trotter or eco, not " +
          th.path_springs);
    }
  }
  JSON_OPT(s, "path_eco_omega_max", th.path_eco_omega_max);
  JSON_OPT(s, "path_gle_file", th.path_gle_file);
  JSON_OPT(s, "path_pile_tau", th.path_pile_tau_input);
  JSON_OPT(s, "path_pile_scale", th.path_pile_scale);
  if (s.contains("path_seed")) {
    const long seed = s.at("path_seed").get<long>();
    if (seed < 0) {
      throw std::invalid_argument("[Dynamics] path_seed must be non-negative");
    }
    th.path_seed = static_cast<std::uint64_t>(seed);
  }
}

static void readBathKeys(const json &s, Parameters &p) {
  auto &th = ParametersLoadAccess::thermostat_options(p);
  if (s.contains("thermostat")) {
    th.kind = lowerCopy(s.at("thermostat").get<std::string>());
  }
  if (s.contains("kind")) {
    th.kind = lowerCopy(s.at("kind").get<std::string>());
  }
  JSON_OPT(s, "andersen_alpha", th.andersen_alpha);
  JSON_OPT(s, "andersen_collision_period", th.andersen_tcol_input);
  JSON_OPT(s, "nose_mass", th.nose_mass);
  JSON_OPT(s, "langevin_friction", th.langevin_friction_input);
  readPathIntegralKeys(s, p);
}

void from_json(const json &j, Parameters &p) {
  // [Main]
  if (j.contains("Main")) {
    auto &m = j.at("Main");
    if (m.contains("job"))
      ParametersLoadAccess::main_options(p).job = enum_from_json(
          m.at("job"), ParametersLoadAccess::main_options(p).job);
    JSON_OPT(m, "random_seed",
             ParametersLoadAccess::main_options(p).randomSeed);
    JSON_OPT(m, "temperature",
             ParametersLoadAccess::main_options(p).temperature);
    JSON_OPT(m, "quiet", ParametersLoadAccess::main_options(p).quiet);
    JSON_OPT(m, "write_log", ParametersLoadAccess::main_options(p).writeLog);
    JSON_OPT(m, "checkpoint", ParametersLoadAccess::main_options(p).checkpoint);
    JSON_OPT(m, "ini_filename",
             ParametersLoadAccess::main_options(p).iniFilename);
    JSON_OPT(m, "con_filename",
             ParametersLoadAccess::main_options(p).conFilename);
    JSON_OPT(m, "finite_difference",
             ParametersLoadAccess::main_options(p).finiteDifference);
    JSON_OPT(m, "max_force_calls",
             ParametersLoadAccess::main_options(p).maxForceCalls);
    JSON_OPT(m, "remove_net_force",
             ParametersLoadAccess::main_options(p).removeNetForce);
    JSON_OPT(m, "write_con_forces",
             ParametersLoadAccess::main_options(p).writeConForces);
  }

  // [Potential]
  if (j.contains("Potential")) {
    auto &s = j.at("Potential");
    if (s.contains("potential"))
      ParametersLoadAccess::potential_options(p).potential =
          enum_from_json(s.at("potential"),
                         ParametersLoadAccess::potential_options(p).potential);
    JSON_OPT(s, "mpi_poll_period",
             ParametersLoadAccess::potential_options(p).MPIPollPeriod);
    JSON_OPT(s, "lammps_logging",
             ParametersLoadAccess::potential_options(p).LAMMPSLogging);
    JSON_OPT(s, "lammps_threads",
             ParametersLoadAccess::potential_options(p).LAMMPSThreads);
    JSON_OPT(s, "emt_rasmussen",
             ParametersLoadAccess::potential_options(p).EMTRasmussen);
    JSON_OPT(s, "log_potential",
             ParametersLoadAccess::potential_options(p).LogPotential);
    JSON_OPT(s, "ext_pot_path",
             ParametersLoadAccess::potential_options(p).extPotPath);
    JSON_OPT(s, "potentials_path",
             ParametersLoadAccess::potential_options(p).potentialsPath);
  }

  // [Structure Comparison]
  if (j.contains("Structure Comparison")) {
    auto &s = j.at("Structure Comparison");
    JSON_OPT(s, "distance_difference",
             ParametersLoadAccess::structure_comparison_options(p)
                 .distance_difference);
    JSON_OPT(
        s, "neighbor_cutoff",
        ParametersLoadAccess::structure_comparison_options(p).neighbor_cutoff);
    JSON_OPT(
        s, "check_rotation",
        ParametersLoadAccess::structure_comparison_options(p).check_rotation);
    JSON_OPT(s, "indistinguishable_atoms",
             ParametersLoadAccess::structure_comparison_options(p)
                 .indistinguishable_atoms);
    JSON_OPT(s, "energy_difference",
             ParametersLoadAccess::structure_comparison_options(p)
                 .energy_difference);
    JSON_OPT(s, "remove_translation",
             ParametersLoadAccess::structure_comparison_options(p)
                 .remove_translation);
  }

  // [Optimizer]
  if (j.contains("Optimizer")) {
    auto &s = j.at("Optimizer");
    if (s.contains("opt_method"))
      ParametersLoadAccess::optimizer_options(p).method =
          enum_from_json(s.at("opt_method"),
                         ParametersLoadAccess::optimizer_options(p).method);
    JSON_OPT(s, "convergence_metric",
             ParametersLoadAccess::optimizer_options(p).convergence_metric);
    for (char &c :
         ParametersLoadAccess::optimizer_options(p).convergence_metric) {
      c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    if (auto label = eonc::helpers::convergenceMetricLabel(
            ParametersLoadAccess::optimizer_options(p).convergence_metric)) {
      ParametersLoadAccess::optimizer_options(p).convergence_metric_label =
          std::string(*label);
    } else {
      throw std::invalid_argument(std::format(
          "unknown convergence_metric: {}",
          ParametersLoadAccess::optimizer_options(p).convergence_metric));
    }
    JSON_OPT(s, "max_iterations",
             ParametersLoadAccess::optimizer_options(p).max_iterations);
    JSON_OPT(s, "max_move",
             ParametersLoadAccess::optimizer_options(p).max_move);
    JSON_OPT(s, "converged_force",
             ParametersLoadAccess::optimizer_options(p).converged_force);
    JSON_OPT(s, "time_step",
             ParametersLoadAccess::optimizer_options(p).time_step_input);
    JSON_OPT(s, "max_time_step",
             ParametersLoadAccess::optimizer_options(p).max_time_step_input);
    if (s.contains("LBFGS")) {
      auto &l = s.at("LBFGS");
      JSON_OPT(l, "memory",
               ParametersLoadAccess::optimizer_options(p).lbfgs.memory);
      JSON_OPT(
          l, "inverse_curvature",
          ParametersLoadAccess::optimizer_options(p).lbfgs.inverse_curvature);
      JSON_OPT(l, "max_inverse_curvature",
               ParametersLoadAccess::optimizer_options(p)
                   .lbfgs.max_inverse_curvature);
      JSON_OPT(l, "auto_scale",
               ParametersLoadAccess::optimizer_options(p).lbfgs.auto_scale);
      JSON_OPT(l, "angle_reset",
               ParametersLoadAccess::optimizer_options(p).lbfgs.angle_reset);
      JSON_OPT(l, "distance_reset",
               ParametersLoadAccess::optimizer_options(p).lbfgs.distance_reset);
      JSON_OPT(l, "curvature",
               ParametersLoadAccess::optimizer_options(p).lbfgs.curvature);
      JSON_OPT(l, "project_rigid",
               ParametersLoadAccess::optimizer_options(p).lbfgs.project_rigid);
      JSON_OPT(l, "secant",
               ParametersLoadAccess::optimizer_options(p).lbfgs.secant);
      JSON_OPT(l, "precon",
               ParametersLoadAccess::optimizer_options(p).lbfgs.precon);
      JSON_OPT(l, "step",
               ParametersLoadAccess::optimizer_options(p).lbfgs.step);
      JSON_OPT(l, "h0", ParametersLoadAccess::optimizer_options(p).lbfgs.h0);
      JSON_OPT(l, "accept",
               ParametersLoadAccess::optimizer_options(p).lbfgs.accept);
      JSON_OPT(l, "extra_updates",
               ParametersLoadAccess::optimizer_options(p).lbfgs.extra_updates);
      JSON_OPT(l, "cautious_eps",
               ParametersLoadAccess::optimizer_options(p).lbfgs.cautious_eps);
      JSON_OPT(l, "cautious_alpha",
               ParametersLoadAccess::optimizer_options(p).lbfgs.cautious_alpha);
      JSON_OPT(l, "precon_A",
               ParametersLoadAccess::optimizer_options(p).lbfgs.precon_A);
      JSON_OPT(l, "precon_mu",
               ParametersLoadAccess::optimizer_options(p).lbfgs.precon_mu);
      JSON_OPT(l, "precon_rcut",
               ParametersLoadAccess::optimizer_options(p).lbfgs.precon_rcut);
    }
    if (s.contains("Xtsci")) {
      auto &x = s.at("Xtsci");
      JSON_OPT(x, "method",
               ParametersLoadAccess::optimizer_options(p).xtsci.method);
    }
    JSON_OPT(s, "xtsci_method",
             ParametersLoadAccess::optimizer_options(p).xtsci.method);
  }

  // [Dynamics] Thermostat names on this object match the ini file.
  // A later Thermostat object overrides the same keys.
  if (j.contains("Dynamics")) {
    auto &s = j.at("Dynamics");
    JSON_OPT(s, "time_step",
             ParametersLoadAccess::dynamics_options(p).time_step_input);
    JSON_OPT(s, "time", ParametersLoadAccess::dynamics_options(p).time_input);
    readBathKeys(s, p);
  }

  // [Thermostat]
  if (j.contains("Thermostat")) {
    readBathKeys(j.at("Thermostat"), p);
  }

  // [Nudged Elastic Band]
  if (j.contains("Nudged Elastic Band")) {
    auto &s = j.at("Nudged Elastic Band");
    JSON_OPT(s, "images", ParametersLoadAccess::neb_options(p).image_count);
    JSON_OPT(s, "max_iterations",
             ParametersLoadAccess::neb_options(p).max_iterations);
    if (s.contains("opt_method"))
      ParametersLoadAccess::neb_options(p).opt_method = enum_from_json(
          s.at("opt_method"), ParametersLoadAccess::neb_options(p).opt_method);
    JSON_OPT(s, "converged_force",
             ParametersLoadAccess::neb_options(p).force_tolerance);
    JSON_OPT(s, "solid_state",
             ParametersLoadAccess::neb_options(p).solid_state.enabled);
    JSON_OPT(s, "solid_state_weight",
             ParametersLoadAccess::neb_options(p).solid_state.weight);
    JSON_OPT(s, "solid_state_pressure",
             ParametersLoadAccess::neb_options(p).solid_state.pressure);
    JSON_OPT(s, "temperature",
             ParametersLoadAccess::neb_options(p).quantum_temperature);
    if (s.contains("spring")) {
      auto &sp = s.at("spring");
      JSON_OPT(sp, "constant",
               ParametersLoadAccess::neb_options(p).spring.constant);
      JSON_OPT(sp, "elastic_band",
               ParametersLoadAccess::neb_options(p).spring.use_elastic_band);
      JSON_OPT(sp, "doubly_nudged",
               ParametersLoadAccess::neb_options(p).spring.doubly_nudged);
      JSON_OPT(sp, "geometric",
               ParametersLoadAccess::neb_options(p).spring.geometric);
    }
    if (s.contains("climbing_image")) {
      auto &ci = s.at("climbing_image");
      JSON_OPT(ci, "enabled",
               ParametersLoadAccess::neb_options(p).climbing_image.enabled);
      JSON_OPT(
          ci, "converged_only",
          ParametersLoadAccess::neb_options(p).climbing_image.converged_only);
      JSON_OPT(ci, "band_slack",
               ParametersLoadAccess::neb_options(p).climbing_image.band_slack);
    }
    if (s.contains("zoom")) {
      auto &z = s.at("zoom");
      auto &zoom = ParametersLoadAccess::neb_options(p).zoom;
      JSON_OPT(z, "enabled", zoom.enabled);
      JSON_OPT(z, "alpha", zoom.alpha);
      JSON_OPT(z, "offset", zoom.offset);
      if (z.contains("mode")) {
        zoom.mode = enum_from_json(z.at("mode"), zoom.mode);
      }
      JSON_OPT(z, "activation_threshold", zoom.activation_threshold);
      if (z.contains("interpolation")) {
        zoom.interpolation =
            enum_from_json(z.at("interpolation"), zoom.interpolation);
      }
      JSON_OPT(z, "stability_count", zoom.stability_count);
      JSON_OPT(z, "max_iterations", zoom.max_iterations);
    }
  }

  // [Dimer]
  if (j.contains("Dimer")) {
    auto &s = j.at("Dimer");
    JSON_OPT(s, "rotation_angle",
             ParametersLoadAccess::dimer_options(p).rotation_angle);
    JSON_OPT(s, "improved", ParametersLoadAccess::dimer_options(p).improved);
    JSON_OPT(s, "converged_angle",
             ParametersLoadAccess::dimer_options(p).converged_angle);
    JSON_OPT(s, "max_iterations",
             ParametersLoadAccess::dimer_options(p).max_iterations);
    if (s.contains("opt_method"))
      ParametersLoadAccess::dimer_options(p).opt_method =
          enum_from_json(s.at("opt_method"),
                         ParametersLoadAccess::dimer_options(p).opt_method);
    if (s.contains("rotation_backend"))
      ParametersLoadAccess::dimer_options(p).rotation_backend = enum_from_json(
          s.at("rotation_backend"),
          ParametersLoadAccess::dimer_options(p).rotation_backend);
    JSON_OPT(s, "lor_residual_tol",
             ParametersLoadAccess::dimer_options(p).lor_residual_tol);
  }

  // [Saddle Search]
  if (j.contains("Saddle Search")) {
    auto &s = j.at("Saddle Search");
    JSON_OPT(s, "method",
             ParametersLoadAccess::saddle_search_options(p).method);
    JSON_OPT(s, "min_mode_method",
             ParametersLoadAccess::saddle_search_options(p).minmode_method);
    JSON_OPT(s, "max_energy",
             ParametersLoadAccess::saddle_search_options(p).max_energy);
    JSON_OPT(s, "max_iterations",
             ParametersLoadAccess::saddle_search_options(p).max_iterations);
    JSON_OPT(s, "displace_magnitude",
             ParametersLoadAccess::saddle_search_options(p).displace_magnitude);
    JSON_OPT(s, "displace_radius",
             ParametersLoadAccess::saddle_search_options(p).displace_radius);
  }

  // [Serve]
  if (j.contains("Serve")) {
    auto &s = j.at("Serve");
    JSON_OPT(s, "host", ParametersLoadAccess::serve_options(p).host);
    if (s.contains("port"))
      ParametersLoadAccess::serve_options(p).port =
          s.at("port").get<uint16_t>();
    if (s.contains("replicas"))
      ParametersLoadAccess::serve_options(p).replicas =
          s.at("replicas").get<size_t>();
    if (s.contains("gateway_port"))
      ParametersLoadAccess::serve_options(p).gateway_port =
          s.at("gateway_port").get<uint16_t>();
    JSON_OPT(s, "endpoints", ParametersLoadAccess::serve_options(p).endpoints);
  }

  // [Debug]
  if (j.contains("Debug")) {
    auto &s = j.at("Debug");
    JSON_OPT(s, "write_movies",
             ParametersLoadAccess::debug_options(p).write_movies);
    JSON_OPT(s, "write_movies_interval",
             ParametersLoadAccess::debug_options(p).write_movies_interval);
    JSON_OPT(s, "write_deprecated_outs",
             ParametersLoadAccess::debug_options(p).write_deprecated_outs);
  }

  // [Hessian] phva_atoms wins over the legacy atom_list, as in the ini.
  if (j.contains("Hessian")) {
    auto &s = j.at("Hessian");
    auto &h = ParametersLoadAccess::hessian_options(p);
    if (s.contains("phva_atoms")) {
      h.phva_atoms = lowerCopy(s.at("phva_atoms").get<std::string>());
    } else if (s.contains("atom_list")) {
      h.phva_atoms = lowerCopy(s.at("atom_list").get<std::string>());
    }
    JSON_OPT(s, "zero_freq_value", h.zero_freq_value);
    if (s.contains("fd_scheme")) {
      h.fd_scheme = lowerCopy(s.at("fd_scheme").get<std::string>());
    }
    JSON_OPT(s, "resume", h.resume);
    JSON_OPT(s, "checkpoint_path", h.checkpoint_path);
    JSON_OPT(s, "write_modes", h.write_modes);
  }

  // [Instanton]
  if (j.contains("Instanton")) {
    auto &s = j.at("Instanton");
    auto &o = ParametersLoadAccess::instanton_options(p);
    JSON_OPT(s, "mode", o.mode);
    JSON_OPT(s, "reactant_filename", o.reactant_filename);
    JSON_OPT(s, "product_filename", o.product_filename);
    JSON_OPT(s, "initial_path", o.initial_path);
    JSON_OPT(s, "beads", o.beads);
    JSON_OPT(s, "beta_hbar_omega", o.beta_hbar_omega);
    JSON_OPT(s, "max_iterations", o.max_iterations);
    JSON_OPT(s, "force_tolerance", o.force_tolerance);
    JSON_OPT(s, "hessian_stride", o.hessian_stride);
    JSON_OPT(s, "saddle_filename", o.saddle_filename);
    JSON_OPT(s, "temperature", o.temperature);
    if (s.contains("temperatures")) {
      o.temperatures = temperaturesFromJson(s.at("temperatures"));
    }
    JSON_OPT(s, "half_ring", o.half_ring);
    JSON_OPT(s, "symmetrize", o.symmetrize);
    if (s.contains("initial_hessians")) {
      o.initial_hessians = s.at("initial_hessians").get<std::string>();
    }
    if (o.initial_hessians != "saddle" &&
        o.initial_hessians != "finite_difference") {
      throw std::invalid_argument(
          "[Instanton] initial_hessians must be saddle or "
          "finite_difference, not " +
          o.initial_hessians);
    }
    JSON_OPT(s, "energy_shift", o.energy_shift);
    if (s.contains("discretization")) {
      o.discretization.clear();
      for (const auto &item : s.at("discretization")) {
        o.discretization.push_back(item.get<double>());
      }
    }
    if (s.contains("friction")) {
      o.friction = lowerCopy(s.at("friction").get<std::string>());
    }
    if (o.friction != "none" && o.friction != "implicit" &&
        o.friction != "explicit") {
      throw std::invalid_argument(
          "[Instanton] friction must be none, implicit or explicit, not " +
          o.friction);
    }
    JSON_OPT(s, "friction_eta", o.friction_eta);
    if (s.contains("friction_eta_beads")) {
      o.friction_eta_beads.clear();
      for (const auto &item : s.at("friction_eta_beads")) {
        o.friction_eta_beads.push_back(item.get<double>());
      }
    }
    JSON_OPT(s, "bead_ladder", o.bead_ladder);
    JSON_OPT(s, "hessian_final", o.hessian_final);
    if (s.contains("springs")) {
      o.springs = lowerCopy(s.at("springs").get<std::string>());
    }
    if (o.springs != "trotter" && o.springs != "eco") {
      throw std::invalid_argument(
          "[Instanton] springs must be trotter or eco, not " + o.springs);
    }
    if (o.mode != "splitting" && o.mode != "rate") {
      throw std::invalid_argument("[Instanton] mode must be splitting or "
                                  "rate, not " +
                                  o.mode);
    }
    JSON_OPT(s, "pi_planes", o.pi_planes);
    JSON_OPT(s, "pi_beads", o.pi_beads);
    JSON_OPT(s, "pi_equilibration_steps", o.pi_equilibration_steps);
    JSON_OPT(s, "pi_sampling_steps", o.pi_sampling_steps);
    JSON_OPT(s, "pi_time_step", o.pi_time_step);
    if (s.contains("pi_thermostat")) {
      o.pi_thermostat = lowerCopy(s.at("pi_thermostat").get<std::string>());
    }
    JSON_OPT(s, "pi_gle_file", o.pi_gle_file);
    JSON_OPT(s, "pi_pile_tau", o.pi_pile_tau);
    JSON_OPT(s, "pi_pile_scale", o.pi_pile_scale);
    JSON_OPT(s, "pi_seed", o.pi_seed);
    if (s.contains("pi_direction")) {
      o.pi_direction = lowerCopy(s.at("pi_direction").get<std::string>());
    }
    JSON_OPT(s, "pi_reactant_extent", o.pi_reactant_extent);
    JSON_OPT(s, "pi_recrossing_parents", o.pi_recrossing_parents);
    JSON_OPT(s, "pi_recrossing_children", o.pi_recrossing_children);
    JSON_OPT(s, "pi_recrossing_time", o.pi_recrossing_time);
    JSON_OPT(s, "pi_recrossing_spacing", o.pi_recrossing_spacing);
    eonc::piqtst::validateOptions(o);
  }

  // The ini loader derives these from the femtosecond inputs whether or
  // not the keys are present; validate_and_link does not.
  {
    auto &th = ParametersLoadAccess::thermostat_options(p);
    const double timeUnit = ParametersLoadAccess::constants(p).timeUnit;
    th.andersen_tcol = timeUnit > 0.0 ? th.andersen_tcol_input / timeUnit : 0.0;
    th.path_pile_tau = timeUnit > 0.0 ? th.path_pile_tau_input / timeUnit : 0.0;
  }

  project_ssot_json_read(j, p);

  // Resolve computed fields
  validate_and_link(p);
}

int load_json(std::string_view json_str, Parameters &params) {
  try {
    auto j = json::parse(json_str);
    from_json(j, params);
    return 0;
  } catch (const json::exception &e) {
    return 1;
  } catch (const std::invalid_argument &) {
    return 1;
  }
}

} // namespace eonc::config
