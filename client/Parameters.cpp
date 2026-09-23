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
#include "eon/Parameters.h"
#include "eon/ParametersINI.h"
#include "eon/ParametersJSON.h"
#include "eon/ParametersSSOT.h"
#include "magic_enum/magic_enum.hpp"

#include <INIReader.h>

#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <string_view>
#include <utility>

#include "eon/EonLogger.h"

namespace {

bool ams_engine_name_is(std::string_view got, std::string_view name) {
  if (got.size() != name.size()) {
    return false;
  }
  for (std::size_t i = 0; i < got.size(); ++i) {
    const auto left = static_cast<unsigned char>(got[i]);
    const auto right = static_cast<unsigned char>(name[i]);
    if (std::toupper(left) != std::toupper(right)) {
      return false;
    }
  }
  return true;
}

} // namespace

namespace eonc {

struct Parameters::Impl {
  std::string last_source;
  int last_error{0};
};

Parameters::Parameters()
    : impl_(std::make_unique<Impl>()) {
  // Covered groups: defaults originate from schema/eon_params.capnp via
  // apply_ssot_defaults (codegen). Uncovered groups still use NSDMI.
  eonc::config::apply_ssot_defaults(*this);
  // Resolve computed fields via validate_and_link.
  eonc::config::validate_and_link(*this);
}

Parameters::~Parameters() = default;
Parameters::Parameters(Parameters &&) noexcept = default;
Parameters &Parameters::operator=(Parameters &&) noexcept = default;

Parameters::Parameters(const Parameters &other)
    : constants_(other.constants_),
      main_options_(other.main_options_),
      potential_options_(other.potential_options_),
      ams_options_(other.ams_options_),
      xtb_options_(other.xtb_options_),
      zbl_options_(other.zbl_options_),
      dftd_options_(other.dftd_options_),
      expr_options_(other.expr_options_),
      mopac_options_(other.mopac_options_),
      socket_nwchem_options_(other.socket_nwchem_options_),
      rgpot_options_(other.rgpot_options_),
      structure_comparison_options_(other.structure_comparison_options_),
      process_search_options_(other.process_search_options_),
      saddle_search_options_(other.saddle_search_options_),
      optimizer_options_(other.optimizer_options_),
      dimer_options_(other.dimer_options_),
      gpr_dimer_options_(other.gpr_dimer_options_),
      gp_surrogate_options_(other.gp_surrogate_options_),
      catlearn_options_(other.catlearn_options_),
      ase_orca_options_(other.ase_orca_options_),
      ase_nwchem_options_(other.ase_nwchem_options_),
      metatomic_options_(other.metatomic_options_),
      lanczos_options_(other.lanczos_options_),
      davidson_options_(other.davidson_options_),
      prefactor_options_(other.prefactor_options_),
      hessian_options_(other.hessian_options_),
      neb_options_(other.neb_options_),
      dynamics_options_(other.dynamics_options_),
      parallel_replica_options_(other.parallel_replica_options_),
      tad_options_(other.tad_options_),
      thermostat_options_(other.thermostat_options_),
      replica_exchange_options_(other.replica_exchange_options_),
      hyperdynamics_options_(other.hyperdynamics_options_),
      basin_hopping_options_(other.basin_hopping_options_),
      global_optimization_options_(other.global_optimization_options_),
      monte_carlo_options_(other.monte_carlo_options_),
      bgsd_options_(other.bgsd_options_),
      serve_options_(other.serve_options_),
      artn_options_(other.artn_options_),
      ira_options_(other.ira_options_),
      debug_options_(other.debug_options_),
      oh_tst_options_(other.oh_tst_options_),
      impl_(other.impl_ ? std::make_unique<Impl>(*other.impl_)
                        : std::make_unique<Impl>()) {}

Parameters &Parameters::operator=(const Parameters &other) {
  if (this == &other) {
    return *this;
  }
  Parameters tmp(other);
  *this = std::move(tmp);
  return *this;
}

std::string_view Parameters::last_load_source() const {
  return impl_ ? std::string_view{impl_->last_source} : std::string_view{};
}

int Parameters::last_load_error() const {
  return impl_ ? impl_->last_error : 0;
}

void Parameters::record_load(std::string_view source, int error) {
  if (!impl_) {
    impl_ = std::make_unique<Impl>();
  }
  impl_->last_source.assign(source);
  impl_->last_error = error;
}

int Parameters::load(std::string_view filename) {
  INIReader ini{std::string(filename)};
  if (ini.ParseError() < 0) {
    EONC_LOG_ERROR("Can't load INI file: {}", filename);
    record_load(filename, 1);
    return 1;
  }

  int error = eonc::config::load_ini(ini, *this);

  // Sanity Checks
  if (parallel_replica_options_.state_check_interval > dynamics_options_.time &&
      magic_enum::enum_name<JobType>(main_options_.job) == "parallel_replica") {
    EONC_LOG_ERROR("[Parallel Replica] state_check_interval must be <= time");
    error = 1;
  }

  if (!neb_options_.initialization.input_path.empty() &&
      neb_options_.initialization.method == NEBInit::LINEAR) {
    EONC_LOG_WARNING(
        "[Nudged Elastic Band] 'initial_path_in' is provided, but "
        "'initializer' defaults to linear. "
        "Ensure this is intentional, as the loaded path will not be "
        "used without initializer set to file.");
  }

  if (saddle_search_options_.dynamics.record_interval_input >
      saddle_search_options_.dynamics.state_check_interval_input) {
    EONC_LOG_ERROR("[Saddle Search] dynamics_record_interval must be <= "
                   "dynamics_state_check_interval");
    error = 1;
  }

  if (potential_options_.potential == PotType::AMS ||
      potential_options_.potential == PotType::AMS_IO) {
    // generate_run allows DFTB with resources and FORCEFIELD alone.
    const bool dftb_with_resources =
        ams_engine_name_is(ams_options_.engine, "DFTB") &&
        !ams_options_.resources.empty();
    const bool forcefield_engine =
        ams_engine_name_is(ams_options_.engine, "FORCEFIELD");
    if (ams_options_.forcefield.empty() && ams_options_.model.empty() &&
        ams_options_.xc.empty() && !dftb_with_resources && !forcefield_engine) {
      EONC_LOG_ERROR("[AMS] Must provide atleast forcefield or model or xc");
      error = 1;
    }

    if (!ams_options_.forcefield.empty() && !ams_options_.model.empty() &&
        !ams_options_.xc.empty()) {
      EONC_LOG_ERROR("[AMS] Must provide either forcefield or model");
      error = 1;
    }
  }

  record_load(filename, error);
  return error;
}

int Parameters::load(FILE *file) {
  constexpr std::string_view kFileSource{"<FILE*>"};
  if (!file) {
    EONC_LOG_ERROR("Can't load INI from a null FILE*");
    record_load(kFileSource, 1);
    return 1;
  }
  if (fseek(file, 0, SEEK_END) != 0) {
    EONC_LOG_ERROR("Can't seek INI FILE*");
    record_load(kFileSource, 1);
    return 1;
  }
  const long size = ftell(file);
  if (size < 0) {
    EONC_LOG_ERROR("Can't tell INI FILE* size");
    record_load(kFileSource, 1);
    return 1;
  }
  if (fseek(file, 0, SEEK_SET) != 0) {
    EONC_LOG_ERROR("Can't rewind INI FILE*");
    record_load(kFileSource, 1);
    return 1;
  }

  std::string buffer(static_cast<size_t>(size), '\0');
  if (fread(buffer.data(), 1, static_cast<size_t>(size), file) !=
      static_cast<size_t>(size)) {
    EONC_LOG_ERROR("Couldn't read the ini file from FILE*");
    record_load(kFileSource, 1);
    return 1;
  }

  INIReader ini(buffer.c_str(), buffer.size());
  if (ini.ParseError() < 0) {
    EONC_LOG_ERROR("Couldn't parse the ini file from FILE*");
    record_load(kFileSource, 1);
    return 1;
  }

  int error = eonc::config::load_ini(ini, *this);
  // Same validation as the filename overload
  // (duplicated intentionally for now; will be extracted to validate())
  if (parallel_replica_options_.state_check_interval > dynamics_options_.time &&
      magic_enum::enum_name<JobType>(main_options_.job) == "parallel_replica") {
    EONC_LOG_ERROR("[Parallel Replica] state_check_interval must be <= time");
    error = 1;
  }
  if (saddle_search_options_.dynamics.record_interval_input >
      saddle_search_options_.dynamics.state_check_interval_input) {
    EONC_LOG_ERROR("[Saddle Search] dynamics_record_interval must be <= "
                   "dynamics_state_check_interval");
    error = 1;
  }
  record_load(kFileSource, error);
  return error;
}

int Parameters::load_ini_text(std::string_view ini_text) {
  constexpr std::string_view kIniSource{"<ini>"};
  INIReader ini(ini_text.data(), ini_text.size());
  if (ini.ParseError() < 0) {
    EONC_LOG_ERROR("Couldn't parse INI from memory");
    record_load(kIniSource, 1);
    return 1;
  }
  const int error = eonc::config::load_ini(ini, *this);
  record_load(kIniSource, error);
  return error;
}

int Parameters::load_json(std::string_view json_str) {
  constexpr std::string_view kJsonSource{"<json>"};
  const int error = eonc::config::load_json(json_str, *this);
  record_load(kJsonSource, error);
  return error;
}

std::string Parameters::to_json() const {
  return eonc::config::to_json(*this).dump(2);
}

} // namespace eonc
