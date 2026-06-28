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
#include "ParametersJSON.h"
#include "Parameters.h"
#include "ParametersINI.h"
#include "magic_enum/magic_enum.hpp"

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
      {"job", enum_to_json(p.main_options.job)},
      {"random_seed", p.main_options.randomSeed},
      {"temperature", p.main_options.temperature},
      {"quiet", p.main_options.quiet},
      {"write_log", p.main_options.writeLog},
      {"checkpoint", p.main_options.checkpoint},
      {"ini_filename", p.main_options.iniFilename},
      {"con_filename", p.main_options.conFilename},
      {"finite_difference", p.main_options.finiteDifference},
      {"max_force_calls", p.main_options.maxForceCalls},
      {"remove_net_force", p.main_options.removeNetForce},
  };

  // [Potential]
  j["Potential"] = {
      {"potential", enum_to_json(p.potential_options.potential)},
      {"mpi_poll_period", p.potential_options.MPIPollPeriod},
      {"lammps_logging", p.potential_options.LAMMPSLogging},
      {"lammps_threads", p.potential_options.LAMMPSThreads},
      {"emt_rasmussen", p.potential_options.EMTRasmussen},
      {"log_potential", p.potential_options.LogPotential},
      {"ext_pot_path", p.potential_options.extPotPath},
      {"potentials_path", p.potential_options.potentialsPath},
  };

  // [Structure Comparison]
  j["Structure Comparison"] = {
      {"distance_difference",
       p.structure_comparison_options.distance_difference},
      {"neighbor_cutoff", p.structure_comparison_options.neighbor_cutoff},
      {"check_rotation", p.structure_comparison_options.check_rotation},
      {"indistinguishable_atoms",
       p.structure_comparison_options.indistinguishable_atoms},
      {"energy_difference", p.structure_comparison_options.energy_difference},
      {"remove_translation", p.structure_comparison_options.remove_translation},
  };

  // [Optimizer]
  j["Optimizer"] = {
      {"opt_method", enum_to_json(p.optimizer_options.method)},
      {"convergence_metric", p.optimizer_options.convergence_metric},
      {"max_iterations", p.optimizer_options.max_iterations},
      {"max_move", p.optimizer_options.max_move},
      {"converged_force", p.optimizer_options.converged_force},
      {"time_step", p.optimizer_options.time_step_input},
      {"max_time_step", p.optimizer_options.max_time_step_input},
  };
  j["Optimizer"]["LBFGS"] = {
      {"memory", p.optimizer_options.lbfgs.memory},
      {"inverse_curvature", p.optimizer_options.lbfgs.inverse_curvature},
      {"auto_scale", p.optimizer_options.lbfgs.auto_scale},
      {"angle_reset", p.optimizer_options.lbfgs.angle_reset},
      {"distance_reset", p.optimizer_options.lbfgs.distance_reset},
  };

  // [Dynamics]
  j["Dynamics"] = {
      {"time_step", p.dynamics_options.time_step_input},
      {"time", p.dynamics_options.time_input},
  };

  // [Thermostat]
  j["Thermostat"] = {
      {"kind", p.thermostat_options.kind},
      {"andersen_alpha", p.thermostat_options.andersen_alpha},
      {"andersen_collision_period", p.thermostat_options.andersen_tcol_input},
      {"nose_mass", p.thermostat_options.nose_mass},
      {"langevin_friction", p.thermostat_options.langevin_friction_input},
  };

  // [Nudged Elastic Band]
  j["Nudged Elastic Band"] = {
      {"images", p.neb_options.image_count},
      {"max_iterations", p.neb_options.max_iterations},
      {"opt_method", enum_to_json(p.neb_options.opt_method)},
      {"converged_force", p.neb_options.force_tolerance},
  };
  j["Nudged Elastic Band"]["spring"] = {
      {"constant", p.neb_options.spring.constant},
      {"elastic_band", p.neb_options.spring.use_elastic_band},
      {"doubly_nudged", p.neb_options.spring.doubly_nudged},
  };
  j["Nudged Elastic Band"]["climbing_image"] = {
      {"enabled", p.neb_options.climbing_image.enabled},
      {"converged_only", p.neb_options.climbing_image.converged_only},
  };

  // [Dimer]
  j["Dimer"] = {
      {"rotation_angle", p.dimer_options.rotation_angle},
      {"improved", p.dimer_options.improved},
      {"converged_angle", p.dimer_options.converged_angle},
      {"max_iterations", p.dimer_options.max_iterations},
      {"opt_method", p.dimer_options.opt_method},
      {"rotation_backend", enum_to_json(p.dimer_options.rotation_backend)},
  };

  // [Saddle Search]
  j["Saddle Search"] = {
      {"method", p.saddle_search_options.method},
      {"min_mode_method", p.saddle_search_options.minmode_method},
      {"max_energy", p.saddle_search_options.max_energy},
      {"max_iterations", p.saddle_search_options.max_iterations},
      {"displace_magnitude", p.saddle_search_options.displace_magnitude},
      {"displace_radius", p.saddle_search_options.displace_radius},
  };

  // [Prefactor]
  j["Prefactor"] = {
      {"default_value", p.prefactor_options.default_value},
      {"max_value", p.prefactor_options.max_value},
      {"min_value", p.prefactor_options.min_value},
  };

  // [Hessian]
  j["Hessian"] = {
      {"atom_list", p.hessian_options.atom_list},
      {"zero_freq_value", p.hessian_options.zero_freq_value},
  };

  // [Debug]
  j["Debug"] = {
      {"write_movies", p.debug_options.write_movies},
      {"write_movies_interval", p.debug_options.write_movies_interval},
      {"write_deprecated_outs", p.debug_options.write_deprecated_outs},
  };

  // [Serve]
  j["Serve"] = {
      {"host", p.serve_options.host},
      {"port", p.serve_options.port},
      {"replicas", p.serve_options.replicas},
      {"gateway_port", p.serve_options.gateway_port},
      {"endpoints", p.serve_options.endpoints},
  };

  return j;
}

void from_json(const json &j, Parameters &p) {
  // [Main]
  if (j.contains("Main")) {
    auto &m = j.at("Main");
    if (m.contains("job"))
      p.main_options.job = enum_from_json(m.at("job"), p.main_options.job);
    JSON_OPT(m, "random_seed", p.main_options.randomSeed);
    JSON_OPT(m, "temperature", p.main_options.temperature);
    JSON_OPT(m, "quiet", p.main_options.quiet);
    JSON_OPT(m, "write_log", p.main_options.writeLog);
    JSON_OPT(m, "checkpoint", p.main_options.checkpoint);
    JSON_OPT(m, "ini_filename", p.main_options.iniFilename);
    JSON_OPT(m, "con_filename", p.main_options.conFilename);
    JSON_OPT(m, "finite_difference", p.main_options.finiteDifference);
    JSON_OPT(m, "max_force_calls", p.main_options.maxForceCalls);
    JSON_OPT(m, "remove_net_force", p.main_options.removeNetForce);
  }

  // [Potential]
  if (j.contains("Potential")) {
    auto &s = j.at("Potential");
    if (s.contains("potential"))
      p.potential_options.potential =
          enum_from_json(s.at("potential"), p.potential_options.potential);
    JSON_OPT(s, "mpi_poll_period", p.potential_options.MPIPollPeriod);
    JSON_OPT(s, "lammps_logging", p.potential_options.LAMMPSLogging);
    JSON_OPT(s, "lammps_threads", p.potential_options.LAMMPSThreads);
    JSON_OPT(s, "emt_rasmussen", p.potential_options.EMTRasmussen);
    JSON_OPT(s, "log_potential", p.potential_options.LogPotential);
    JSON_OPT(s, "ext_pot_path", p.potential_options.extPotPath);
    JSON_OPT(s, "potentials_path", p.potential_options.potentialsPath);
  }

  // [Structure Comparison]
  if (j.contains("Structure Comparison")) {
    auto &s = j.at("Structure Comparison");
    JSON_OPT(s, "distance_difference",
             p.structure_comparison_options.distance_difference);
    JSON_OPT(s, "neighbor_cutoff",
             p.structure_comparison_options.neighbor_cutoff);
    JSON_OPT(s, "check_rotation",
             p.structure_comparison_options.check_rotation);
    JSON_OPT(s, "indistinguishable_atoms",
             p.structure_comparison_options.indistinguishable_atoms);
    JSON_OPT(s, "energy_difference",
             p.structure_comparison_options.energy_difference);
    JSON_OPT(s, "remove_translation",
             p.structure_comparison_options.remove_translation);
  }

  // [Optimizer]
  if (j.contains("Optimizer")) {
    auto &s = j.at("Optimizer");
    if (s.contains("opt_method"))
      p.optimizer_options.method =
          enum_from_json(s.at("opt_method"), p.optimizer_options.method);
    JSON_OPT(s, "convergence_metric", p.optimizer_options.convergence_metric);
    JSON_OPT(s, "max_iterations", p.optimizer_options.max_iterations);
    JSON_OPT(s, "max_move", p.optimizer_options.max_move);
    JSON_OPT(s, "converged_force", p.optimizer_options.converged_force);
    JSON_OPT(s, "time_step", p.optimizer_options.time_step_input);
    JSON_OPT(s, "max_time_step", p.optimizer_options.max_time_step_input);
    if (s.contains("LBFGS")) {
      auto &l = s.at("LBFGS");
      JSON_OPT(l, "memory", p.optimizer_options.lbfgs.memory);
      JSON_OPT(l, "inverse_curvature",
               p.optimizer_options.lbfgs.inverse_curvature);
      JSON_OPT(l, "auto_scale", p.optimizer_options.lbfgs.auto_scale);
      JSON_OPT(l, "angle_reset", p.optimizer_options.lbfgs.angle_reset);
      JSON_OPT(l, "distance_reset", p.optimizer_options.lbfgs.distance_reset);
    }
  }

  // [Dynamics]
  if (j.contains("Dynamics")) {
    auto &s = j.at("Dynamics");
    JSON_OPT(s, "time_step", p.dynamics_options.time_step_input);
    JSON_OPT(s, "time", p.dynamics_options.time_input);
  }

  // [Thermostat]
  if (j.contains("Thermostat")) {
    auto &s = j.at("Thermostat");
    JSON_OPT(s, "kind", p.thermostat_options.kind);
    JSON_OPT(s, "andersen_alpha", p.thermostat_options.andersen_alpha);
    JSON_OPT(s, "andersen_collision_period",
             p.thermostat_options.andersen_tcol_input);
    JSON_OPT(s, "nose_mass", p.thermostat_options.nose_mass);
    JSON_OPT(s, "langevin_friction",
             p.thermostat_options.langevin_friction_input);
  }

  // [Nudged Elastic Band]
  if (j.contains("Nudged Elastic Band")) {
    auto &s = j.at("Nudged Elastic Band");
    JSON_OPT(s, "images", p.neb_options.image_count);
    JSON_OPT(s, "max_iterations", p.neb_options.max_iterations);
    if (s.contains("opt_method"))
      p.neb_options.opt_method =
          enum_from_json(s.at("opt_method"), p.neb_options.opt_method);
    JSON_OPT(s, "converged_force", p.neb_options.force_tolerance);
    if (s.contains("spring")) {
      auto &sp = s.at("spring");
      JSON_OPT(sp, "constant", p.neb_options.spring.constant);
      JSON_OPT(sp, "elastic_band", p.neb_options.spring.use_elastic_band);
      JSON_OPT(sp, "doubly_nudged", p.neb_options.spring.doubly_nudged);
    }
    if (s.contains("climbing_image")) {
      auto &ci = s.at("climbing_image");
      JSON_OPT(ci, "enabled", p.neb_options.climbing_image.enabled);
      JSON_OPT(ci, "converged_only",
               p.neb_options.climbing_image.converged_only);
    }
  }

  // [Dimer]
  if (j.contains("Dimer")) {
    auto &s = j.at("Dimer");
    JSON_OPT(s, "rotation_angle", p.dimer_options.rotation_angle);
    JSON_OPT(s, "improved", p.dimer_options.improved);
    JSON_OPT(s, "converged_angle", p.dimer_options.converged_angle);
    JSON_OPT(s, "max_iterations", p.dimer_options.max_iterations);
    JSON_OPT(s, "opt_method", p.dimer_options.opt_method);
    if (s.contains("rotation_backend"))
      p.dimer_options.rotation_backend = enum_from_json(
          s.at("rotation_backend"), p.dimer_options.rotation_backend);
  }

  // [Saddle Search]
  if (j.contains("Saddle Search")) {
    auto &s = j.at("Saddle Search");
    JSON_OPT(s, "method", p.saddle_search_options.method);
    JSON_OPT(s, "min_mode_method", p.saddle_search_options.minmode_method);
    JSON_OPT(s, "max_energy", p.saddle_search_options.max_energy);
    JSON_OPT(s, "max_iterations", p.saddle_search_options.max_iterations);
    JSON_OPT(s, "displace_magnitude",
             p.saddle_search_options.displace_magnitude);
    JSON_OPT(s, "displace_radius", p.saddle_search_options.displace_radius);
  }

  // [Serve]
  if (j.contains("Serve")) {
    auto &s = j.at("Serve");
    JSON_OPT(s, "host", p.serve_options.host);
    if (s.contains("port"))
      p.serve_options.port = s.at("port").get<uint16_t>();
    if (s.contains("replicas"))
      p.serve_options.replicas = s.at("replicas").get<size_t>();
    if (s.contains("gateway_port"))
      p.serve_options.gateway_port = s.at("gateway_port").get<uint16_t>();
    JSON_OPT(s, "endpoints", p.serve_options.endpoints);
  }

  // [Debug]
  if (j.contains("Debug")) {
    auto &s = j.at("Debug");
    JSON_OPT(s, "write_movies", p.debug_options.write_movies);
    JSON_OPT(s, "write_movies_interval", p.debug_options.write_movies_interval);
    JSON_OPT(s, "write_deprecated_outs", p.debug_options.write_deprecated_outs);
  }

  // Resolve computed fields
  validate_and_link(p);
}

int load_json(const std::string &json_str, Parameters &params) {
  try {
    auto j = json::parse(json_str);
    from_json(j, params);
    return 0;
  } catch (const json::exception &e) {
    return 1;
  }
}

} // namespace eonc::config
