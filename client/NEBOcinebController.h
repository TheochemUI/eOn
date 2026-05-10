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
#include "Parameters.h"

namespace eonc {
class NudgedElasticBand; // forward declaration
}

namespace eonc::neb {

/// Internal terminal status for the OCINEB MMF dimer step. Distinct
/// from the broader SaddleStatus vocabulary because positive
/// curvature on the climbing image is an *expected* short-circuit
/// for the controller (the CI sat at a local minimum, not a saddle)
/// rather than a generic saddle-search failure.
enum class MMFStatus : int {
  /// MMF converged to a saddle on the climbing image, alignment
  /// passed the angle_tol gate; cache the mode for the next call.
  Helped = 0,
  /// MMF hit the iteration cap without converging; CI position is
  /// still useable but flagged as "MMF didn't help this round".
  MaxIterations = 1,
  /// MMF refused to act -- input was malformed (invalid climbing
  /// image index or vanishing tangent), or the converged mode
  /// drifted outside angle_tol after MMF. Treated by run() as "MMF
  /// didn't help this round".
  Skipped = -1,
  /// Positive curvature detected -- the CI sat at a local minimum,
  /// not a saddle. run() restores the pre-MMF position to prevent
  /// a force explosion in the next NEB step.
  PositiveCurvature = -2,
};

[[nodiscard]] constexpr int to_int(MMFStatus s) noexcept {
  return static_cast<int>(s);
}

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

  // Warm-start: cache converged eigenvector for next dimer call
  bool has_cached_mode_{false};
  AtomMatrix cached_mode_;

  MMFStatus runDimer(eonc::NudgedElasticBand &neb, double &alignment);
  void updateThresholdSuccess(double convForce, double newForce);
  void updateThresholdBackoff(double alignment);
};

} // namespace eonc::neb
