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

#include <memory>

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
 * configured through the config.toml file. The job execution framework supports
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
  // No need to track the output files, which can vary by user parameter anyway,
  // a boolean is sufficient
  // TODO(rg) :: Consider populating a struct, RunResults with the boolean and
  // the files generated if needed (YAGNI)
  virtual bool run() = 0;
  // virtual JobBase* clone() const = 0;
};

template <typename T> class Job : public JobBase {
public:
  bool run() override { return static_cast<T *>(this)->runImpl(); }
  virtual bool runImpl() = 0;
  // For cloning
  std::unique_ptr<T> clone() const {
    return std::unique_ptr<T>(new T(*static_cast<const T *>(this)));
  }
};

} // namespace eonc
