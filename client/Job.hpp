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
#include "client/PointJob.h"
#include "thirdparty/toml.hpp"
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

using JobVariant = std::variant<PointJob /*, OtherJobTypes */>;

JobVariant mkJob(const toml::table &);

// Templated struct to handle job execution
template <typename... Args> struct JobRunnerImpl {
  std::tuple<Args &&...> args; // Store arguments as references

  // Perfectly forward the arguments to the tuple
  JobRunnerImpl(Args &&...args)
      : args(std::forward<Args>(args)...) {}

  // Callable operator to handle the job type
  template <typename JobType> bool operator()(JobType &job) const {
    return applyImpl(job, std::index_sequence_for<Args...>{});
  }

private:
  // Helper function to unpack the tuple and forward arguments to the job's
  // runImpl method
  template <typename JobType, std::size_t... I>
  bool applyImpl(JobType &job, std::index_sequence<I...>) const {
    return job.runImpl(std::forward<Args>(std::get<I>(args))...);
  }
};

// Wrapper function to execute the JobRunnerImpl with the variant
template <typename... Args>
bool JobRunner(JobVariant &jobVariant, Args &&...args) {
  return std::visit(JobRunnerImpl<Args &&...>(std::forward<Args>(args)...),
                    jobVariant);
}

} // namespace eonc
