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
#include "Parser.hpp"
#include "RelaxJob.hpp"
#include "StructureComparisonJob.h"
#include "client/PointJob.h"
#include "thirdparty/toml.hpp"
#include <stdexcept>
#include <variant>

namespace eonc {
/** @defgroup Jobs
 *
 * \brief ClientEON main procedures
 *
 * This page provides links to all of the available jobs that can be run by the
 * ClientEON, as well as documentation on the job class, and the overview
 * section relating the job structure to the rest of the program.
 *
 */

/**
 * @file
 * @ingroup Jobs
 *
 * \brief The Job template class is used to serve as a base class for all job
 * types, providing a common interface to execute jobs at runtime based on the
 * derived job type.
 *
 * The variant ensures compile-time polymorphism and type safety. Each derived
 * job class must implement the runImpl method, providing the specific behavior
 * for that job.
 *
 * Jobs can be executed based on runtime parameters, and their behavior can be
 * configured through the config.toml file. The job execution framework supports
 * both standalone jobs and jobs that are part of larger routines. Some jobs do
 * not involve optimizers and are documented in their own respective files.

 *
 * \note The run method must be implemented by each derived job class.
 */

using JobVariant = std::variant<PointJob, RelaxJob,
                                StructureComparisonJob /*, OtherJobTypes */>;

JobVariant mkJob(const toml::table &);

struct JobRunnerImpl {
  // Generic template that produces a compile-time error
  template <typename JobType, typename... Args>
  bool operator()(JobType &, Args &&...) const {
    static_assert(
        always_false<JobType>::value,
        "No explicit specialization for this job type in JobRunnerImpl.");
    return false; // This return will never be reached, but is required
                  // syntactically.
  }

  template <typename... Args>
  bool operator()(StructureComparisonJob &job, Args &&...args) const {
    if constexpr (sizeof...(Args) == 2) {
      return job.runImpl(std::forward<Args>(args)...);
    }
  }

  template <typename... Args>
  bool operator()(PointJob &job, Args &&...args) const {
    if constexpr (sizeof...(Args) == 1) {
      return job.runImpl(std::forward<Args>(args)...);
    }
  }

  template <typename... Args>
  bool operator()(RelaxJob &job, Args &&...args) const {
    if constexpr (sizeof...(Args) == 1) {
      return job.runImpl(std::forward<Args>(args)...);
    }
    // Shouldn't reach here..
    throw std::runtime_error("Job dispatched incorrectly");
    return false;
  }

private:
  // Helper template that is always false to trigger static_assert
  template <typename T> struct always_false : std::false_type {};
};

template <typename... Args>
bool JobRunner(JobVariant &jobVariant, Args &&...args) {
  return std::visit(
      [&](auto &job) -> bool {
        return JobRunnerImpl()(job, std::forward<Args>(args)...);
      },
      jobVariant);
}

bool runJob(JobVariant &, std::vector<Matter> &);

} // namespace eonc
