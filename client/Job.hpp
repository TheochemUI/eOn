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

template <typename... Args>
JobVariant mkJob(const toml::table &config, Args &&...args) {
  config_section(config, "Main");
  auto jtype = get_enum_toml<JobType>(config["Main"]["job"]);

  switch (jtype) {
  case JobType::Point: {
    return PointJob(std::forward<Args>(args)...);
  }
  default: {
    throw std::runtime_error("No known job could be constructed");
  }
  }
}

struct JobRunner {
  template <typename T> bool operator()(T &job) const { return job.runImpl(); }
};

} // namespace eonc
