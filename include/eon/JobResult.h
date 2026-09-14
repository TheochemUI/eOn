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
#include <stdexcept>
#include <string>

#include "magic_enum/magic_enum.hpp"

namespace eonc {

/// In-memory JobResult scalars (eOn-4mrf). Matches schema/eon_job_result.capnp
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

  void writeResultsDat(const std::string &path) const {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
      throw std::runtime_error("JobResultEnvelope: cannot open " + path);
    }
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
};

} // namespace eonc
