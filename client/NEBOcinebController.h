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

namespace eonc {
class NudgedElasticBand; // forward declaration
}

namespace eonc::neb {

/// Goswami (in prep).
/// Controller for OCINEB hybrid dimer refinement of the climbing image.
/// Manages MMF triggering, threshold adaptation, and backoff logic.
class OCINEBController {
public:
  struct Config {
    bool use_mmf;
    double trigger_force;
    double trigger_factor;
    long max_steps;
    long ci_stability_count;
    double angle_tol;
    double penalty_strength;
    double penalty_base;
    double force_tolerance;
  };

  static Config fromParams(const Parameters &params);
  explicit OCINEBController(const Config &cfg);

  void initBaseline(double baseline_force);

  bool shouldTrigger(double convForce, bool ci_active, long climbingImage,
                     long numImages, int ciStabilityCounter) const;

  struct MMFResult {
    double newForce;
    bool convergedAfterMMF;
    bool shouldResetOptimizer;
  };

  MMFResult run(eonc::NudgedElasticBand &neb, double convForce);

  void resetStability() { ciStabilityCounter_ = 0; }
  void updateStability(long climbingImage);
  int stabilityCount() const { return ciStabilityCounter_; }
  double threshold() const { return current_threshold_; }

private:
  Config cfg_;
  double current_threshold_{-1.0};
  double baseline_force_{-1.0};
  long previousClimbingImage_{-1};
  int ciStabilityCounter_{0};
  int mmf_iterations_used_{0};

  int runDimer(eonc::NudgedElasticBand &neb, double &alignment);
  void updateThresholdSuccess(double convForce, double newForce);
  void updateThresholdBackoff(double alignment);
};

} // namespace eonc::neb
