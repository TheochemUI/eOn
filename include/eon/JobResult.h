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
#include <string_view>
#include <utility>
#include <vector>

#include "magic_enum/magic_enum.hpp"

namespace eonc {

/// Optimizer, xtsci, eindir, and rgpot identity carried on a JobResult.
/// Empty backend means the record is absent. Legacy results.dat parsing
/// is unchanged: these lines are extra value-key records.
struct JobResultProvenance {
  std::string backend;
  std::string engine_version;
  std::string engine_build_identity;
  std::string rgpot_name;
  std::string rgpot_version;
  bool xtsci{false};
  std::uint16_t xts_abi_major{0};
  std::uint16_t xts_abi_minor{0};
  std::uint16_t xts_abi_layout{0};
  bool has_eindir{false};
  std::uint32_t eindir_abi_major{0};
  std::uint32_t eindir_abi_minor{0};
  std::uint32_t eindir_objective_layout{0};
  std::uint64_t eindir_objective_size{0};
  std::uint64_t eindir_objective_align{0};
  std::uint32_t eindir_dlpack_major{0};
  std::uint32_t eindir_dlpack_minor{0};
  std::uint64_t eindir_features{0};

  static constexpr std::uint16_t readcon_spec_version = 3;
  static constexpr std::string_view readcon_min_version = "0.14.9";
  static constexpr std::string_view eon_schema_min_version = "0.2.3";
  static constexpr std::string_view rgpycrumbs_min_version = "1.10.4";
  static constexpr std::string_view chemparseplot_min_version = "1.9.17";
  static constexpr std::string_view rgpot_pin = "3.2.0";

  std::string text() const {
    if (backend.empty()) {
      return {};
    }
    std::ostringstream out;
    out << backend << " optimizer_backend\n";
    out << "eon.optimizer.v1 optimizer_provenance_schema\n";
    out << "eon.compatibility.v1 compatibility_schema\n";
    out << readcon_spec_version << " compatibility_readcon_spec_version\n";
    out << readcon_min_version << " compatibility_readcon_min_version\n";
    out << eon_schema_min_version << " compatibility_eon_schema_min_version\n";
    out << rgpycrumbs_min_version << " compatibility_rgpycrumbs_min_version\n";
    out << chemparseplot_min_version
        << " compatibility_chemparseplot_min_version\n";
    out << "eon engine_id\n";
    if (!engine_version.empty()) {
      out << engine_version << " engine_version\n";
    }
    if (!engine_build_identity.empty()) {
      out << engine_build_identity << " engine_build_identity\n";
    }
    if (!rgpot_name.empty()) {
      out << rgpot_name << " rgpot_name\n";
    }
    if (!rgpot_version.empty()) {
      out << "eon.rgpot.v1 rgpot_schema\n";
      out << rgpot_version << " rgpot_version\n";
    }
    if (xtsci) {
      out << "eon.objective compatibility_engine_protocol_family\n";
      out << "1 compatibility_engine_protocol_major\n";
      out << "0 compatibility_engine_protocol_minor\n";
      out << xts_abi_major << " compatibility_engine_abi_major\n";
      out << xts_abi_minor << " compatibility_engine_abi_minor\n";
      out << xts_abi_layout << " compatibility_engine_layout_revision\n";
      out << xts_abi_major << " optimizer_xts_abi_major\n";
      out << xts_abi_minor << " optimizer_xts_abi_minor\n";
      out << xts_abi_layout << " optimizer_xts_abi_layout\n";
    }
    if (has_eindir) {
      out << eindir_abi_major << " optimizer_eindir_abi_major\n";
      out << eindir_abi_minor << " optimizer_eindir_abi_minor\n";
      out << eindir_objective_layout << " optimizer_eindir_objective_layout\n";
      out << eindir_objective_size << " optimizer_eindir_objective_size\n";
      out << eindir_objective_align << " optimizer_eindir_objective_align\n";
      out << eindir_dlpack_major << " optimizer_eindir_dlpack_major\n";
      out << eindir_dlpack_minor << " optimizer_eindir_dlpack_minor\n";
      out << eindir_features << " optimizer_eindir_features\n";
    }
    return out.str();
  }
};

class Parameters;

/// Optimizer backend, engine identity, rgpot pin, and xts ABI from a job.
/// eindir fields stay at zero unless a real stamp filled them.
JobResultProvenance provenanceForJob(const Parameters &params);

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
  JobResultProvenance provenance;

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
    out << provenance.text();
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
    out << provenance.text();
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
