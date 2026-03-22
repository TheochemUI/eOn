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
#include "EonLogger.h"
#include "Job.h"
#include "Matter.h"

#include <memory>
#include <string>
#include <vector>

namespace eonc {

/// Base class for replica dynamics jobs (TAD, SafeHyper).
///
/// Provides shared infrastructure: state checking via minimize-and-compare,
/// binary search refinement of transition times, dephasing, and results I/O.
/// Subclasses implement initExtra() for job-specific matter objects and
/// dynamics() for the accelerated MD loop.
class ReplicaDynamicsJob : public Job {
public:
  using Job::Job;
  ~ReplicaDynamicsJob() override = default;
  std::vector<std::string> run() override;

protected:
  // -- Shared matter objects --
  std::shared_ptr<Matter> current;
  std::shared_ptr<Matter> reactant;
  std::shared_ptr<Matter> saddle;
  std::shared_ptr<Matter> product;
  std::shared_ptr<Matter> finalState;    // state at transition detection
  std::shared_ptr<Matter> finalStateTmp; // candidate before confirmation

  // -- Flags --
  bool metaStateFlag{false};
  bool newStateFlag{false};

  // -- Force call counters --
  size_t minimizeFCalls{0};
  size_t mdFCalls{0};
  size_t dephaseFCalls{0};
  size_t refineFCalls{0};

  // -- Dynamics state --
  long transitionStep{0};
  double time{0.0};
  double minCorrectedTime{1.0e200};
  double transitionTime{0.0};

  std::vector<std::string> returnFiles;
  eonc::log::Scoped log;

  // -- Shared methods --

  /// Minimize a copy of current and compare to reactant.
  /// Returns true if a transition (new state) was detected.
  bool checkState(Matter *current, Matter *reactant);

  /// Binary search for the transition frame in a snapshot buffer.
  long refine(const std::vector<std::shared_ptr<Matter>> &buff,
              Matter *reactant);

  /// Dephase the trajectory to ensure thermal independence.
  void dephase();

  /// Write results.dat and structure con files.
  void saveData(int status);

  // -- Template method hooks --

  /// Create any extra matter objects needed by the subclass.
  virtual void initExtra() {}

  /// The accelerated dynamics loop. Returns status (1 = transition, 0 = none).
  virtual int dynamics() = 0;

  /// Post-dynamics logging (override for job-specific messages).
  virtual void reportResults() {}
};

} // namespace eonc

using eonc::ReplicaDynamicsJob;
