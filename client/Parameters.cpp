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
#include "Parameters.h"
#include "ParametersINI.h"
#include "ParametersJSON.h"
#include "magic_enum/magic_enum.hpp"

#include <INIReader.h>

#include <cerrno>
#include <cstdio>
#include <cstring>

#include "EonLogger.h"

Parameters::Parameters() {
  // All simple defaults are NSDMI in Parameters.h.
  // Resolve computed fields via validate_and_link.
  eonc::config::validate_and_link(*this);
}

int Parameters::load(std::string filename) {
  INIReader ini(filename);
  if (ini.ParseError() < 0) {
    EONC_LOG_ERROR("Can't load INI file: {}", filename);
    return 1;
  }

  int error = eonc::config::load_ini(ini, *this);

  // Sanity Checks
  if (parallel_replica_options.state_check_interval > dynamics_options.time &&
      magic_enum::enum_name<JobType>(main_options.job) == "parallel_replica") {
    EONC_LOG_ERROR("[Parallel Replica] state_check_interval must be <= time");
    error = 1;
  }

  if (!neb_options.initialization.input_path.empty() &&
      neb_options.initialization.method == NEBInit::LINEAR) {
    EONC_LOG_WARNING(
        "[Nudged Elastic Band] 'initial_path_in' is provided, but "
        "'initializer' defaults to linear. "
        "Ensure this is intentional, as the loaded path will not be "
        "used without initializer set to file.");
  }

  if (saddle_search_options.dynamics.record_interval_input >
      saddle_search_options.dynamics.state_check_interval_input) {
    EONC_LOG_ERROR("[Saddle Search] dynamics_record_interval must be <= "
                   "dynamics_state_check_interval");
    error = 1;
  }

  if (potential_options.potential == PotType::AMS ||
      potential_options.potential == PotType::AMS_IO) {
    if (ams_options.forcefield.empty() && ams_options.model.empty() &&
        ams_options.xc.empty()) {
      EONC_LOG_ERROR("[AMS] Must provide atleast forcefield or model or xc");
      error = 1;
    }

    if (!ams_options.forcefield.empty() && !ams_options.model.empty() &&
        !ams_options.xc.empty()) {
      EONC_LOG_ERROR("[AMS] Must provide either forcefield or model");
      error = 1;
    }
  }

  return error;
}

int Parameters::load(FILE *file) {
  // Legacy FILE* overload: read into string buffer, use INIReader buffer ctor
  fseek(file, 0, SEEK_END);
  long size = ftell(file);
  fseek(file, 0, SEEK_SET);

  std::string buffer(size, '\0');
  fread(buffer.data(), 1, size, file);

  INIReader ini(buffer.c_str(), buffer.size());
  if (ini.ParseError() < 0) {
    EONC_LOG_ERROR("Couldn't parse the ini file from FILE*");
    return 1;
  }

  int error = eonc::config::load_ini(ini, *this);
  // Same validation as the filename overload
  // (duplicated intentionally for now; will be extracted to validate())
  if (parallel_replica_options.state_check_interval > dynamics_options.time &&
      magic_enum::enum_name<JobType>(main_options.job) == "parallel_replica") {
    EONC_LOG_ERROR("[Parallel Replica] state_check_interval must be <= time");
    error = 1;
  }
  if (saddle_search_options.dynamics.record_interval_input >
      saddle_search_options.dynamics.state_check_interval_input) {
    EONC_LOG_ERROR("[Saddle Search] dynamics_record_interval must be <= "
                   "dynamics_state_check_interval");
    error = 1;
  }
  return error;
}

int Parameters::load_json(const std::string &json_str) {
  return eonc::config::load_json(json_str, *this);
}

std::string Parameters::to_json() const {
  return eonc::config::to_json(*this).dump(2);
}
