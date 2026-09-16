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
#include "eon/Parameters.h"
#include "eon/ParametersINI.h"
#include "magic_enum/magic_enum.hpp"

#include <cctype>
#include <format>
#include <nlohmann/json.hpp>
#include <stdexcept>

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
  };

  // [Nudged Elastic Band]
  j["Nudged Elastic Band"] = {
      {"images", ParametersLoadAccess::neb_options(p).image_count},
      {"max_iterations", ParametersLoadAccess::neb_options(p).max_iterations},
      {"opt_method",
       enum_to_json(ParametersLoadAccess::neb_options(p).opt_method)},
      {"converged_force", ParametersLoadAccess::neb_options(p).force_tolerance},
  };
  j["Nudged Elastic Band"]["spring"] = {
      {"constant", ParametersLoadAccess::neb_options(p).spring.constant},
      {"elastic_band",
       ParametersLoadAccess::neb_options(p).spring.use_elastic_band},
      {"doubly_nudged",
       ParametersLoadAccess::neb_options(p).spring.doubly_nudged},
  };
  j["Nudged Elastic Band"]["climbing_image"] = {
      {"enabled", ParametersLoadAccess::neb_options(p).climbing_image.enabled},
      {"converged_only",
       ParametersLoadAccess::neb_options(p).climbing_image.converged_only},
      {"band_slack",
       ParametersLoadAccess::neb_options(p).climbing_image.band_slack},
  };

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
  };

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

  return j;
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

  // [Dynamics]
  if (j.contains("Dynamics")) {
    auto &s = j.at("Dynamics");
    JSON_OPT(s, "time_step",
             ParametersLoadAccess::dynamics_options(p).time_step_input);
    JSON_OPT(s, "time", ParametersLoadAccess::dynamics_options(p).time_input);
  }

  // [Thermostat]
  if (j.contains("Thermostat")) {
    auto &s = j.at("Thermostat");
    JSON_OPT(s, "kind", ParametersLoadAccess::thermostat_options(p).kind);
    JSON_OPT(s, "andersen_alpha",
             ParametersLoadAccess::thermostat_options(p).andersen_alpha);
    JSON_OPT(s, "andersen_collision_period",
             ParametersLoadAccess::thermostat_options(p).andersen_tcol_input);
    JSON_OPT(s, "nose_mass",
             ParametersLoadAccess::thermostat_options(p).nose_mass);
    JSON_OPT(
        s, "langevin_friction",
        ParametersLoadAccess::thermostat_options(p).langevin_friction_input);
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
    if (s.contains("spring")) {
      auto &sp = s.at("spring");
      JSON_OPT(sp, "constant",
               ParametersLoadAccess::neb_options(p).spring.constant);
      JSON_OPT(sp, "elastic_band",
               ParametersLoadAccess::neb_options(p).spring.use_elastic_band);
      JSON_OPT(sp, "doubly_nudged",
               ParametersLoadAccess::neb_options(p).spring.doubly_nudged);
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
