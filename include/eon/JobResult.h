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
#pragma once

#include "BaseStructures.h"

#include <cstdint>
#include <format>
#include <fstream>
#include <optional>
#include <readcon-core.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "magic_enum/magic_enum.hpp"

namespace eonc {

/// In-memory JobResult scalars. Matches schema/eon_job_result.capnp
/// field names as results.dat keys. Geometries stay on Matter until capnp
/// codegen lands.
struct JobResultEnvelope {
  std::string job_type;
  int status_code{0};
  std::string status_text;
  std::string potential_type;
  std::int64_t random_seed{-1};
  std::uint64_t force_calls{0};
  double potential_energy{0.0};
  bool has_energy{false};
  double potential_energy_saddle{0.0};
  double potential_energy_reactant{0.0};
  double potential_energy_product{0.0};
  bool has_saddle{false};
  bool has_reactant{false};
  bool has_product{false};
  bool process_search_layout{false};
  std::uint64_t force_calls_minimization{0};
  std::uint64_t force_calls_saddle{0};
  std::uint64_t force_calls_prefactors{0};
  double barrier_reactant_to_product{0.0};
  double barrier_product_to_reactant{0.0};
  double prefactor_reactant_to_product{0.0};
  double prefactor_product_to_reactant{0.0};
  double displacement_saddle_distance{0.0};
  double simulation_time{0.0};
  double md_temperature{0.0};
  bool has_dynamics{false};
  std::optional<readcon::ConFrame> reactant_frame;
  std::optional<readcon::ConFrame> saddle_frame;
  std::optional<readcon::ConFrame> product_frame;
  std::vector<std::pair<std::string, double>> extras;
  std::vector<std::pair<std::string, std::string>> tags;

  std::string processSearchString() const {
    std::ostringstream out;
    out << status_code << " termination_reason\n";
    if (!status_text.empty()) {
      out << status_text << " termination_reason_text\n";
    }
    out << random_seed << " random_seed\n";
    if (!potential_type.empty()) {
      out << potential_type << " potential_type\n";
    }
    out << force_calls << " total_force_calls\n";
    out << force_calls_minimization << " force_calls_minimization\n";
    out << force_calls_saddle << " force_calls_saddle\n";
    out << std::format("{:.12e} potential_energy_saddle\n",
                       potential_energy_saddle);
    out << std::format("{:.12e} potential_energy_reactant\n",
                       potential_energy_reactant);
    out << std::format("{:.12e} potential_energy_product\n",
                       potential_energy_product);
    out << std::format("{:.12e} barrier_reactant_to_product\n",
                       barrier_reactant_to_product);
    out << std::format("{:.12e} barrier_product_to_reactant\n",
                       barrier_product_to_reactant);
    out << std::format("{:.12e} displacement_saddle_distance\n",
                       displacement_saddle_distance);
    if (has_dynamics) {
      out << std::format("{:.12e} simulation_time\n", simulation_time);
      out << std::format("{:.12e} md_temperature\n", md_temperature);
    }
    out << force_calls_prefactors << " force_calls_prefactors\n";
    out << std::format("{:.12e} prefactor_reactant_to_product\n",
                       prefactor_reactant_to_product);
    out << std::format("{:.12e} prefactor_product_to_reactant\n",
                       prefactor_product_to_reactant);
    return out.str();
  }

  std::string toString() const {
    if (process_search_layout) {
      return processSearchString();
    }
    std::ostringstream out;
    out << status_code << " termination_reason\n";
    if (!status_text.empty()) {
      out << status_text << " termination_reason_text\n";
    }
    if (!job_type.empty()) {
      out << job_type << " job_type\n";
    }
    if (!potential_type.empty()) {
      out << potential_type << " potential_type\n";
    }
    if (random_seed >= 0) {
      out << random_seed << " random_seed\n";
    }
    out << force_calls << " total_force_calls\n";
    if (has_energy) {
      out << std::format("{:.12e} potential_energy\n", potential_energy);
    }
    if (has_saddle) {
      out << std::format("{:.12e} potential_energy_saddle\n",
                         potential_energy_saddle);
    }
    if (has_reactant) {
      out << std::format("{:.12e} potential_energy_reactant\n",
                         potential_energy_reactant);
    }
    if (has_product) {
      out << std::format("{:.12e} potential_energy_product\n",
                         potential_energy_product);
    }
    for (const auto &kv : extras) {
      out << std::format("{:.12e} {}\n", kv.second, kv.first);
    }
    for (const auto &kv : tags) {
      out << kv.second << " " << kv.first << "\n";
    }
    return out.str();
  }

  void writeResultsDat(const std::string &path) const {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
      throw std::runtime_error("JobResultEnvelope: cannot open " + path);
    }
    out << toString();
    if (!out) {
      throw std::runtime_error("JobResultEnvelope: write failed for " + path);
    }
  }

  static JobResultEnvelope fromMinimization(RunStatus status, PotType pot,
                                            std::uint64_t fcalls, bool hasE,
                                            double energy) {
    JobResultEnvelope e;
    e.job_type = "minimization";
    e.status_code = static_cast<int>(status);
    e.status_text = std::string(magic_enum::enum_name<RunStatus>(status));
    e.potential_type = std::string(magic_enum::enum_name<PotType>(pot));
    e.force_calls = fcalls;
    e.has_energy = hasE;
    e.potential_energy = energy;
    return e;
  }

  static JobResultEnvelope fromProcessSearch(
      int status, std::string statusText, PotType pot, std::int64_t seed,
      std::uint64_t fMin, std::uint64_t fSaddle, std::uint64_t fPref,
      double eSaddle, double eReactant, double eProduct, double barrierFwd,
      double barrierRev, double displacement, double prefFwd, double prefRev,
      bool dynamics, double simTime, double mdTemp) {
    JobResultEnvelope e;
    e.process_search_layout = true;
    e.job_type = "process_search";
    e.status_code = status;
    e.status_text = std::move(statusText);
    e.potential_type = std::string(magic_enum::enum_name<PotType>(pot));
    e.random_seed = seed;
    e.force_calls_minimization = fMin;
    e.force_calls_saddle = fSaddle;
    e.force_calls_prefactors = fPref;
    e.force_calls = fMin + fSaddle + fPref;
    e.has_saddle = true;
    e.has_reactant = true;
    e.has_product = true;
    e.potential_energy_saddle = eSaddle;
    e.potential_energy_reactant = eReactant;
    e.potential_energy_product = eProduct;
    e.barrier_reactant_to_product = barrierFwd;
    e.barrier_product_to_reactant = barrierRev;
    e.displacement_saddle_distance = displacement;
    e.prefactor_reactant_to_product = prefFwd;
    e.prefactor_product_to_reactant = prefRev;
    e.has_dynamics = dynamics;
    e.simulation_time = simTime;
    e.md_temperature = mdTemp;
    return e;
  }
};

} // namespace eonc
