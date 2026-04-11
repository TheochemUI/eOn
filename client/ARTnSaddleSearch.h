/*
 * This file is part of eOn.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * Copyright (c) 2010--present, eOn Development Team
 * All rights reserved.
 *
 * Repo:
 * https://github.com/TheochemUI/eOn
 */

#pragma once

#include "EonLogger.h"
#include "Matter.h"
#include "MinModeSaddleSearch.h"
#include "Parameters.h"
#include "SaddleSearchMethod.h"
#include <limits>

// Include the ARTn resource for thread-local access
#ifdef WITH_ARTN
#include "libs/ARTn/ARTnResource.h"
#endif

namespace eonc {

/*
 * ARTn integration with thread-local dynamic loading
 */
#ifdef WITH_ARTN

// Constants
constexpr double ARTN_MODE_TOLERANCE = 1e-10;
constexpr double ARTN_SMALL_DISPLACEMENT = 1e-6;

// WARNING: All ARTn operations are serialized through a global mutex
// due to the non-thread-safe Fortran backend. This can be a performance
// bottleneck in multi-threaded contexts. Consider process-level parallelism
// if true concurrent ARTn operations are required.

#endif // WITH_ARTN

/// Saddle search method using the Activation-Relaxation Technique nouveau.
/// Wraps the pARTn Fortran library via its C API (artn.h).
class ARTnSaddleSearch : public SaddleSearchMethod {
public:
  static constexpr int STATUS_GOOD = 0;
  static constexpr int STATUS_BAD_MAX_ITERATIONS =
      MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS;
  static constexpr int STATUS_BAD_ARTN_ERROR = 22;

  ARTnSaddleSearch(std::shared_ptr<Matter> matterPassed,
                   std::shared_ptr<Potential> potPassed, AtomMatrix modeInitial,
                   const Parameters &paramsPassed);
  ~ARTnSaddleSearch() override;

  int run() override;
  double getEigenvalue() override;
  AtomMatrix getEigenvector() override;
  std::string_view describeStatus(int status) const override;
  int getStatus() const override { return status; }
  int getIterationCount() const override { return iteration; }
  int getForceCalls() const override { return forcecalls; }

private:
  std::shared_ptr<Matter> matter;
  double eigenvalue{std::numeric_limits<double>::quiet_NaN()};
  AtomMatrix eigenvector, mode;
  int status{0};
  int iteration{0};
  int forcecalls{0};
  eonc::log::Scoped log;
};

} // namespace eonc

using eonc::ARTnSaddleSearch;
