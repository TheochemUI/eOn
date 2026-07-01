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

#include "Eigen.h"
#include "EigenmodeStrategy.h"
#include "Matter.h"
#include "Optimizer.h"
#include "SaddleSearchMethod.h"

#include <string>

namespace eonc {

class MinModeSaddleSearch : public SaddleSearchMethod {

public:
  enum Status : int {
    // DO NOT CHANGE THE ORDER OF THIS LIST
    STATUS_GOOD,                           // 0
    STATUS_INIT,                           // 1
    STATUS_BAD_NO_CONVEX,                  // 2
    STATUS_BAD_HIGH_ENERGY,                // 3
    STATUS_BAD_MAX_CONCAVE_ITERATIONS,     // 4
    STATUS_BAD_MAX_ITERATIONS,             // 5
    STATUS_BAD_NOT_CONNECTED,              // 6
    STATUS_BAD_PREFACTOR,                  // 7
    STATUS_BAD_HIGH_BARRIER,               // 8
    STATUS_BAD_MINIMA,                     // 9
    STATUS_FAILED_PREFACTOR,               // 10
    STATUS_POTENTIAL_FAILED,               // 11
    STATUS_NONNEGATIVE_ABORT,              // 12
    STATUS_NONLOCAL_ABORT,                 // 13
    STATUS_NEGATIVE_BARRIER,               // 14
    STATUS_BAD_MD_TRAJECTORY_TOO_SHORT,    // 15
    STATUS_BAD_NO_NEGATIVE_MODE_AT_SADDLE, // 16
    STATUS_BAD_NO_BARRIER,                 // 17
    STATUS_ZEROMODE_ABORT,                 // 18
    STATUS_OPTIMIZER_ERROR,                // 19
    STATUS_DIMER_LOST_MODE,                // 20
    STATUS_DIMER_RESTORED_BEST             // 21
  };

  /// Human-readable message for a status code.
  static constexpr std::string_view statusMessage(int status) {
    constexpr std::string_view msgs[] = {
        "Success",                                                // 0
        "Initialized",                                            // 1
        "Initial displacement unable to reach convex region",     // 2
        "Barrier too high",                                       // 3
        "Too many iterations in concave region",                  // 4
        "Too many iterations",                                    // 5
        "Saddle is not connected to initial state",               // 6
        "Prefactors not within window",                           // 7
        "Energy barrier not within window",                       // 8
        "Minimizations from saddle did not converge",             // 9
        "Hessian calculation failed",                             // 10
        "Potential evaluation failed",                            // 11
        "Nonnegative initial mode, aborting",                     // 12
        "Nonlocal abort",                                         // 13
        "Negative barrier detected",                              // 14
        "No reaction found during MD trajectory",                 // 15
        "Converged to stationary point with zero negative modes", // 16
        "No forward barrier found along minimized band",          // 17
        "Zero mode abort",                                        // 18
        "Optimizer error",                                        // 19
        "Dimer lost mode",                                        // 20
        "Dimer restored best",                                    // 21
    };
    if (status >= 0 &&
        status < static_cast<int>(sizeof(msgs) / sizeof(msgs[0])))
      return msgs[status];
    return "Unknown status";
  }

  /// Issue #20 policy: never leave climb as STATUS_GOOD when the climb
  /// objective is still unconverged (unfeasible systems must not look like
  /// success). Used by run(); unit-tested so deleting this rule fails CI.
  [[nodiscard]] static constexpr int
  finalizeClimbStatus(int climbStatus, bool objectiveConverged) noexcept {
    if (climbStatus == STATUS_GOOD && !objectiveConverged) {
      return STATUS_BAD_MAX_ITERATIONS;
    }
    return climbStatus;
  }

  MinModeSaddleSearch(std::shared_ptr<Matter> matterPassed,
                      AtomMatrix modePassed, double reactantEnergyPassed,
                      const Parameters &parametersPassed,
                      std::shared_ptr<Potential> potPassed);
  ~MinModeSaddleSearch() = default;
  AtomMatrix getEigenvector(); // lowest eigenmode
  double getEigenvalue();      // estimate for the lowest eigenvalue
  std::string_view describeStatus(int status) const override {
    return statusMessage(status);
  }
  int getStatus() const override { return status; }
  int getIterationCount() const override { return iteration; }
  int getForceCalls() const override { return forcecalls; }

  int run() override;
  int run(long max_iterations_override);

  int forcecalls{0};
  int iteration{0};
  int status{0};

private:
  AtomMatrix mode;
  AtomMatrix initialTangent_;
  std::shared_ptr<Matter> matter;
  std::shared_ptr<EigenmodeStrategy>
      minModeMethod; // shared with the objective func
  double reactantEnergy;
  eonc::log::Scoped log;
};

} // namespace eonc

using eonc::MinModeSaddleSearch;
