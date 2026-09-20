/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
#include "eon/JobRegistry.h"

#include <map>
#include <stdexcept>

namespace eonc {

std::map<JobType, JobFactory> &jobTable() {
  static std::map<JobType, JobFactory> t;
  return t;
}

void registerJob(JobType type, JobFactory factory) {
  jobTable()[type] = std::move(factory);
}

} // namespace eonc

namespace eonc::helpers {

std::unique_ptr<Job> makeJob(std::unique_ptr<Parameters> params,
                             Runtime runtime) {
  eonc::forceJobRegistration();
  const JobType type = params->main_options().job;
  auto &t = eonc::jobTable();
  auto it = t.find(type);
  if (it == t.end()) {
    throw std::runtime_error("No known job could be constructed");
  }
  return it->second(std::move(params), std::move(runtime));
}

} // namespace eonc::helpers
