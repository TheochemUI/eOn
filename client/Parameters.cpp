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
#include "ParametersImpl.h"
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


Parameters::Parameters()
    : impl_(std::make_unique<Impl>()) {
  // Covered groups: defaults originate from schema/eon_params.capnp via
  // apply_ssot_defaults (codegen). Uncovered groups still use NSDMI.
  eonc::config::apply_ssot_defaults(*this);
  // Resolve computed fields via validate_and_link.
  eonc::config::validate_and_link(*this);
}

int Parameters::load(std::string_view filename) {
  INIReader ini{std::string(filename)};
  if (ini.ParseError() != 0) {
    EONC_LOG_ERROR("Can't load INI file: {}", filename);
    record_load(filename, 1);
    return 1;
  }

  int error = eonc::config::load_ini(ini, *this);

  // Sanity Checks
  if (impl_->parallel_replica_options_.state_check_interval >
          impl_->dynamics_options_.time &&
      magic_enum::enum_name<JobType>(impl_->main_options_.job) ==
          "parallel_replica") {
    EONC_LOG_ERROR("[Parallel Replica] state_check_interval must be <= time");
    error = 1;
  }

  if (!impl_->neb_options_.initialization.input_path.empty() &&
      impl_->neb_options_.initialization.method == NEBInit::LINEAR) {
    EONC_LOG_WARNING(
        "[Nudged Elastic Band] 'initial_path_in' is provided, but "
        "'initializer' defaults to linear. "
        "Ensure this is intentional, as the loaded path will not be "
        "used without initializer set to file.");
  }

  if (impl_->saddle_search_options_.dynamics.record_interval_input >
      impl_->saddle_search_options_.dynamics.state_check_interval_input) {
    EONC_LOG_ERROR("[Saddle Search] dynamics_record_interval must be <= "
                   "dynamics_state_check_interval");
    error = 1;
  }

  if (impl_->potential_options_.potential == PotType::AMS ||
      impl_->potential_options_.potential == PotType::AMS_IO) {
    // generate_run allows DFTB with resources and FORCEFIELD alone.
    const bool dftb_with_resources =
        ams_engine_name_is(impl_->ams_options_.engine, "DFTB") &&
        !impl_->ams_options_.resources.empty();
    const bool forcefield_engine =
        ams_engine_name_is(impl_->ams_options_.engine, "FORCEFIELD");
    if (impl_->ams_options_.forcefield.empty() &&
        impl_->ams_options_.model.empty() && impl_->ams_options_.xc.empty() &&
        !dftb_with_resources && !forcefield_engine) {
      EONC_LOG_ERROR("[AMS] Must provide atleast forcefield or model or xc");
      error = 1;
    }

    if (!impl_->ams_options_.forcefield.empty() &&
        !impl_->ams_options_.model.empty() && !impl_->ams_options_.xc.empty()) {
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
  if (ini.ParseError() != 0) {
    EONC_LOG_ERROR("Couldn't parse the ini file from FILE*");
    record_load(kFileSource, 1);
    return 1;
  }

  int error = eonc::config::load_ini(ini, *this);
  // Same validation as the filename overload
  // (duplicated intentionally for now; will be extracted to validate())
  if (impl_->parallel_replica_options_.state_check_interval >
          impl_->dynamics_options_.time &&
      magic_enum::enum_name<JobType>(impl_->main_options_.job) ==
          "parallel_replica") {
    EONC_LOG_ERROR("[Parallel Replica] state_check_interval must be <= time");
    error = 1;
  }
  if (impl_->saddle_search_options_.dynamics.record_interval_input >
      impl_->saddle_search_options_.dynamics.state_check_interval_input) {
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
  if (ini.ParseError() != 0) {
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
