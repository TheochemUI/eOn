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
#include "NEBOcinebController.h"
#include "EonLogger.h"
#include "MinModeSaddleSearch.h"
#include "NudgedElasticBand.h"
#include "eonExceptions.hpp"
#include <algorithm>
#include <cmath>

namespace eonc::neb {

OCINEBController::Config
OCINEBController::fromParams(const Parameters &params) {
  auto &ci = params.neb_options.climbing_image;
  auto &r = ci.ocineb;
  return Config{
      r.use_mmf,
      r.trigger_force,
      r.trigger_factor,
      r.max_steps,
      r.ci_stability_count,
      r.angle_tol,
      r.penalty.strength,
      r.penalty.base,
      params.neb_options.force_tolerance,
  };
}

OCINEBController::OCINEBController(const Config &cfg)
    : cfg_{cfg} {}

void OCINEBController::initBaseline(double baseline_force) {
  baseline_force_ = baseline_force;
  current_threshold_ = baseline_force_ * cfg_.trigger_factor;
}

bool OCINEBController::shouldTrigger(double convForce, bool ci_active,
                                     long climbingImage, long numImages,
                                     int ciStabilityCounter) const {
  if (!cfg_.use_mmf || !ci_active)
    return false;
  if (climbingImage <= 0 || climbingImage > numImages)
    return false;
  if (ciStabilityCounter <= static_cast<int>(cfg_.ci_stability_count))
    return false;
  if (convForce <= cfg_.force_tolerance)
    return false;
  return (convForce < current_threshold_ || convForce < cfg_.trigger_force);
}

void OCINEBController::updateStability(long climbingImage) {
  if (climbingImage == previousClimbingImage_) {
    ciStabilityCounter_++;
  } else {
    ciStabilityCounter_ = 0;
    previousClimbingImage_ = climbingImage;
  }
}

OCINEBController::MMFResult OCINEBController::run(eonc::NudgedElasticBand &neb,
                                                  double convForce) {
  auto *log = eonc::log::get();

  QUILL_LOG_DEBUG(log,
                  "Triggering MMF.  Force: {:.4f}, Threshold: {:.4f} "
                  "({:.2f}x baseline)",
                  convForce, current_threshold_,
                  current_threshold_ / baseline_force_);

  // Save climbing image state before MMF
  AtomMatrix savedPositions = neb.path[neb.climbingImage]->getPositions();

  double alignment = 0.0;
  int mmfResult = runDimer(neb, alignment);

  // Always update forces after MMF
  neb.movedAfterForceCall = true;
  neb.updateForces();
  double newForce = neb.convergenceForce();

  // Check convergence
  if (newForce < cfg_.force_tolerance) {
    QUILL_LOG_DEBUG(log, "NEB converged after MMF. Force: {:.4f}", newForce);
    return {newForce, true, false};
  }

  bool mmfHelped = (newForce < convForce) && mmfResult != -2;

  if (mmfHelped) {
    updateThresholdSuccess(convForce, newForce);
    QUILL_LOG_DEBUG(log,
                    "MMF helped (status={}). Force: {:.4f} -> {:.4f} "
                    "({:.2f}x baseline). New threshold: {:.4f}",
                    mmfResult, convForce, newForce, newForce / baseline_force_,
                    current_threshold_);
  } else {
    updateThresholdBackoff(alignment);
    QUILL_LOG_DEBUG(
        log,
        "MMF backoff (status={}). Force: {:.4f} -> {:.4f}, "
        "Alignment:  {:.3f}. New threshold: {:.4f} ({:.2f}x baseline)",
        mmfResult, convForce, newForce, alignment, current_threshold_,
        current_threshold_ / baseline_force_);
  }

  bool shouldReset =
      (savedPositions - neb.path[neb.climbingImage]->getPositions()).norm() >
      neb.params.optimizer_options.max_move *
          neb.params.neb_options.image_count;

  if (shouldReset) {
    QUILL_LOG_DEBUG(log, "Resetting optimization history.");
  }

  ciStabilityCounter_ = 0;

  return {newForce, false, shouldReset};
}

int OCINEBController::runDimer(eonc::NudgedElasticBand &neb,
                               double &alignment) {
  auto *log = eonc::log::get();
  alignment = 0.0;

  if (neb.climbingImage <= 0 || neb.climbingImage > neb.numImages) {
    QUILL_LOG_WARNING(log, "Invalid climbing image for MMF:  {}",
                      neb.climbingImage);
    return -1;
  }

  AtomMatrix initialMode = *neb.tangent[neb.climbingImage];
  double tangentNorm = initialMode.norm();
  if (tangentNorm < 1e-8) {
    QUILL_LOG_WARNING(log, "Tangent too small for MMF initialization");
    return -1;
  }
  initialMode /= tangentNorm;

  auto tempMinModeSearch = std::make_shared<MinModeSaddleSearch>(
      neb.path[neb.climbingImage], initialMode,
      neb.path[neb.climbingImage]->getPotentialEnergy(), neb.params, neb.pot);

  int minModeStatus;
  try {
    minModeStatus = tempMinModeSearch->run(cfg_.max_steps);
  } catch (const eonc::DimerModeRestoredException &) {
    minModeStatus = MinModeSaddleSearch::STATUS_DIMER_RESTORED_BEST;
    QUILL_LOG_DEBUG(log, "MMF:  Dimer restored to best state");
  } catch (const eonc::DimerModeLostException &) {
    minModeStatus = MinModeSaddleSearch::STATUS_DIMER_LOST_MODE;
    QUILL_LOG_WARNING(log, "Dimer lost mode during MMF refinement");
  }

  mmf_iterations_used_ += tempMinModeSearch->iteration;

  double eigenvalue = tempMinModeSearch->getEigenvalue();
  if (eigenvalue > 0.0) {
    QUILL_LOG_WARNING(log,
                      "MMF skipped: Positive curvature detected (eig={:.4f}).",
                      eigenvalue);
    return -2;
  }

  AtomMatrix finalModeMatrix = tempMinModeSearch->getEigenvector();
  VectorXd finalMode = VectorXd::Map(finalModeMatrix.data(), 3 * neb.atoms);
  VectorXd currentTangent =
      VectorXd::Map(neb.tangent[neb.climbingImage]->data(), 3 * neb.atoms);
  alignment = std::abs(finalMode.normalized().dot(currentTangent.normalized()));

  if (minModeStatus == MinModeSaddleSearch::STATUS_GOOD ||
      minModeStatus == MinModeSaddleSearch::STATUS_DIMER_RESTORED_BEST) {
    if (alignment < cfg_.angle_tol) {
      QUILL_LOG_WARNING(
          log,
          "MMF converged/restored but mode drifted (alignment={:.3f} < {:.3f})",
          alignment, cfg_.angle_tol);
      return -1;
    }
    return 0;
  } else if (minModeStatus == MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS) {
    return 1;
  } else {
    QUILL_LOG_WARNING(log, "MMF failed.  Mode-tangent alignment: {:.3f}",
                      alignment);
    return -1;
  }
}

void OCINEBController::updateThresholdSuccess(double convForce,
                                              double newForce) {
  current_threshold_ = newForce * (0.5 + 0.4 * (newForce / convForce));
  double max_threshold = baseline_force_ * cfg_.trigger_factor;
  current_threshold_ = std::min(current_threshold_, max_threshold);
}

void OCINEBController::updateThresholdBackoff(double alignment) {
  double alpha = std::clamp(alignment, 0.0, 1.0);
  double penalty_factor =
      cfg_.penalty_base +
      (1.0 - cfg_.penalty_base) * std::pow(alpha, cfg_.penalty_strength);
  penalty_factor = std::clamp(penalty_factor, cfg_.penalty_base, 1.0);

  current_threshold_ = baseline_force_ * cfg_.trigger_factor * penalty_factor;
  // Lower bound on the MMF trigger threshold. Scaled by force_tolerance so
  // the MMF gate never collapses to zero when the NEB is already near
  // convergence, but also capped by the trigger_factor envelope so a
  // loose force_tolerance cannot push min_threshold above the cap and
  // starve MMF activation.
  double min_threshold = std::min(cfg_.force_tolerance * 2.0,
                                  baseline_force_ * cfg_.trigger_factor);
  current_threshold_ = std::max(current_threshold_, min_threshold);
}

} // namespace eonc::neb
