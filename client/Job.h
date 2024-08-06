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
#include <string>
#include <vector>

#include "client/matter/Matter.h"
#include "thirdparty/toml.hpp"

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
 * The template parameter is used to specify the specific job type, ensuring
 * compile-time polymorphism and type safety. Each derived job class must
 * implement the run method, providing the specific behavior for that job.
 *
 * Jobs can be executed based on runtime parameters, and their behavior can be
 * configured through the config.init file. The job execution framework supports
 * both standalone jobs and jobs that are part of larger routines. Some jobs do
 * not involve optimizers and are documented in their own respective files.
 *
 * \tparam T The specific job type derived from the Job base class.
 *
 * \note The run method must be implemented by each derived job class.
 */

class JobBase {
public:
  virtual ~JobBase() = default;
  virtual std::vector<std::string> run() = 0;
};

template <typename T> class Job : public JobBase {
public:
  std::vector<std::string> run() override {
    return static_cast<T *>(this)->run();
  }
};

std::unique_ptr<JobBase>
makeJob(const toml::table &config,
        std::optional<std::reference_wrapper<Matter>> mat = std::nullopt);

} // namespace eonc
