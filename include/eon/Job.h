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
#include "Parameters.h"
#include "Potential.h"
#include "Runtime.h"
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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
 * \brief The job class is used to serve as an abstract class for all jobs,
 *  as well as to call a job at runtime based off of the passed in parameters.
 *
 * The Static members are used to tell at runtime which job to run as set by the
 * parameters, and therefore as set by the config.init file. About half of the
 * jobs are standalone, while others are run from routines with the same name. A
 * certain subset of jobs do not run optimizers (SEE OVERVIEW) and are
 * documented in their own files accordingly.
 *
 */

/**
 * Declaration of job class
 */

class Job {
private:
protected:
  // make const
  JobType jtype;
  Parameters params;
  /// Non-null when this Job owns the composition root (one-shot makeJob).
  std::unique_ptr<Runtime> owned_runtime_;
  /// Always valid: either owned_runtime_.get() or a caller-owned Runtime.
  Runtime *runtime_;
  std::shared_ptr<Potential> pot;

public:
  /// Take ownership of a Runtime previously passed as Runtime&.
  void adoptRuntime(std::unique_ptr<Runtime> rt) {
    if (!rt || runtime_ != rt.get()) {
      throw std::logic_error(
          "Job::adoptRuntime: Runtime is not the borrowed instance");
    }
    owned_runtime_ = std::move(rt);
  }

  /// Borrow: caller keeps Runtime alive (CLI stack / Python Session).
  Job(std::unique_ptr<Parameters> parameters, Runtime &rt)
      : jtype{parameters->main_options().job},
        params{*std::move(parameters)},
        owned_runtime_{},
        runtime_{&rt},
        pot{helpers::makePotential(params.potential_options().potential, params,
                                   rt)} {}
  /// Own a Runtime (one-shot makeJob / rvalue).
  Job(std::unique_ptr<Parameters> parameters, std::unique_ptr<Runtime> rt)
      : Job(std::move(parameters), *rt) {
    adoptRuntime(std::move(rt));
  }
  /// Own a default-constructed Runtime.
  explicit Job(std::unique_ptr<Parameters> parameters)
      : Job(std::move(parameters), std::make_unique<Runtime>()) {}
  Job(std::shared_ptr<Potential> potPassed, const Parameters &parameters)
      : jtype{parameters.main_options().job},
        params{parameters},
        owned_runtime_{std::make_unique<Runtime>()},
        runtime_{owned_runtime_.get()},
        pot{potPassed} {}
  virtual ~Job() = default;
  //! Virtual run; used solely for dynamic dispatch
  virtual std::vector<std::string> run() = 0;
  [[nodiscard]] JobType getType() { return this->jtype; };
  [[nodiscard]] PotRegistry &pots() noexcept { return runtime_->pots(); }
  /// Drop the Potential so on_destroyed is recorded before Runtime dies.
  void releasePotential() { pot.reset(); }
};

namespace helpers {
/// One-shot: Job owns a default Runtime.
std::unique_ptr<Job> makeJob(std::unique_ptr<Parameters> params);
/// Borrow: caller keeps Runtime alive.
std::unique_ptr<Job> makeJob(std::unique_ptr<Parameters> params,
                             Runtime &runtime);
/// Own the rvalue Runtime.
std::unique_ptr<Job> makeJob(std::unique_ptr<Parameters> params,
                             Runtime &&runtime);
/// Own the unique_ptr Runtime.
std::unique_ptr<Job> makeJob(std::unique_ptr<Parameters> params,
                             std::unique_ptr<Runtime> runtime);
} // namespace helpers

} // namespace eonc
