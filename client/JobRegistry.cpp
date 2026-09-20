/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
*/
#include "eon/JobRegistry.h"

#include <map>
#include <memory>
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

namespace {

std::unique_ptr<Job> makeJobFromFactory(std::unique_ptr<Parameters> params,
                                        Runtime &runtime) {
  eonc::forceJobRegistration();
  const JobType type = params->main_options().job;
  auto &t = eonc::jobTable();
  auto it = t.find(type);
  if (it == t.end()) {
    throw std::runtime_error("No known job could be constructed");
  }
  return it->second(std::move(params), runtime);
}

} // namespace

std::unique_ptr<Job> makeJob(std::unique_ptr<Parameters> params,
                             Runtime &runtime) {
  return makeJobFromFactory(std::move(params), runtime);
}

std::unique_ptr<Job> makeJob(std::unique_ptr<Parameters> params,
                             std::unique_ptr<Runtime> runtime) {
  if (!runtime) {
    throw std::logic_error("makeJob: null Runtime");
  }
  Runtime &ref = *runtime;
  auto job = makeJobFromFactory(std::move(params), ref);
  job->adoptRuntime(std::move(runtime));
  return job;
}

std::unique_ptr<Job> makeJob(std::unique_ptr<Parameters> params,
                             Runtime &&runtime) {
  return makeJob(std::move(params),
                 std::make_unique<Runtime>(std::move(runtime)));
}

std::unique_ptr<Job> makeJob(std::unique_ptr<Parameters> params) {
  return makeJob(std::move(params), std::make_unique<Runtime>());
}

} // namespace eonc::helpers
