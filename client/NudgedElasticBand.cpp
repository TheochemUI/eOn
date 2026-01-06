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
#include "NudgedElasticBand.h"
#include "BaseStructures.h"
#include "IDPPObjectiveFunction.hpp"
#include "ImprovedDimer.h"
#include "Lanczos.h"
#include "MinModeSaddleSearch.h"
#include "NEBInitialPaths.hpp"
#include "Optimizer.h"
#include "eonExceptions.hpp"
#include "magic_enum/magic_enum.hpp"

#include <spdlog/spdlog.h>

using namespace helper_functions;
namespace fs = std::filesystem;

// NEBObjectiveFunction definitions
VectorXd NEBObjectiveFunction::getGradient(bool fdstep) {
  VectorXd forceV;
  forceV.resize(3 * neb->atoms * neb->numImages);
  if (neb->movedAfterForceCall)
    neb->updateForces();
  for (long i = 1; i <= neb->numImages; i++) {
    forceV.segment(3 * neb->atoms * (i - 1), 3 * neb->atoms) =
        VectorXd::Map(neb->projectedForce[i]->data(), 3 * neb->atoms);
  }
  return -forceV;
}

double NEBObjectiveFunction::getEnergy() {
  double Energy{0};
  for (long i = 1; i <= neb->numImages; i++) {
    Energy += neb->path[i]->getPotentialEnergy();
  }
  return Energy;
}

void NEBObjectiveFunction::setPositions(VectorXd x) {
  neb->movedAfterForceCall = true;
  for (long i = 1; i <= neb->numImages; i++) {
    neb->path[i]->setPositions(MatrixXd::Map(
        x.segment(3 * neb->atoms * (i - 1), 3 * neb->atoms).data(), neb->atoms,
        3));
  }
}

VectorXd NEBObjectiveFunction::getPositions() {
  VectorXd posV;
  posV.resize(3 * neb->atoms * neb->numImages);
  for (long i = 1; i <= neb->numImages; i++) {
    posV.segment(3 * neb->atoms * (i - 1), 3 * neb->atoms) =
        VectorXd::Map(neb->path[i]->getPositions().data(), 3 * neb->atoms);
  }
  return posV;
}

int NEBObjectiveFunction::degreesOfFreedom() {
  return 3 * neb->numImages * neb->atoms;
}

bool NEBObjectiveFunction::isUncertain() {
  double maxMaxUnc = std::numeric_limits<double>::lowest();
  double currentMaxUnc{0};
  for (long idx = 0; idx <= neb->numImages; idx++) {
    currentMaxUnc = neb->path[idx]->getEnergyVariance();
    if (currentMaxUnc > maxMaxUnc) {
      maxMaxUnc = currentMaxUnc;
    }
  }
  bool unc_conv{maxMaxUnc > params->gp_uncertainity};
  if (unc_conv) {
    this->status = NudgedElasticBand::NEBStatus::MAX_UNCERTAINITY;
  }
  return unc_conv;
}

bool NEBObjectiveFunction::isConverged() {
  bool force_conv = getConvergence() < params->neb_options.force_tolerance;
  return force_conv;
}

double NEBObjectiveFunction::getConvergence() {
  return neb->convergenceForce();
}

VectorXd NEBObjectiveFunction::difference(VectorXd a, VectorXd b) {
  VectorXd pbcDiff(3 * neb->numImages * neb->atoms);
  for (int i = 1; i <= neb->numImages; i++) {
    int n = (i - 1) * 3 * neb->atoms;
    int m = 3 * neb->atoms;
    pbcDiff.segment(n, m) =
        neb->path[i]->pbcV(a.segment(n, m) - b.segment(n, m));
  }
  return pbcDiff;
}

// Nudged Elastic Band definitions
// First constructor: Now delegates to the second constructor
NudgedElasticBand::NudgedElasticBand(
    std::shared_ptr<Matter> initialPassed, std::shared_ptr<Matter> finalPassed,
    std::shared_ptr<Parameters> parametersPassed,
    std::shared_ptr<Potential> potPassed)
    : NudgedElasticBand(
          [&]() {
            auto &init_opt = parametersPassed->neb_options.initialization;
            const size_t base_count = parametersPassed->neb_options.image_count;

            // Apply oversampling factor if flag exists
            const size_t relax_count =
                init_opt.oversampling
                    ? base_count * init_opt.oversampling_factor
                    : base_count;

            std::vector<Matter> path;
            switch (init_opt.method) {
            case NEBInit::FILE: {
              std::vector<fs::path> file_paths =
                  helper_functions::neb_paths::readFilePaths(
                      init_opt.input_path);
              path = helper_functions::neb_paths::filePathInit(
                  file_paths, *initialPassed, base_count);
              break;
            }
            case NEBInit::IDPP:
              path = helper_functions::neb_paths::idppPath(
                  *initialPassed, *finalPassed, relax_count, parametersPassed);
              break;
            case NEBInit::IDPP_COLLECTIVE:
              path = helper_functions::neb_paths::idppCollectivePath(
                  *initialPassed, *finalPassed, relax_count, parametersPassed);
              break;
            case NEBInit::SIDPP:
            case NEBInit::SIDPP_ZBL:
              path = helper_functions::neb_paths::sidppPath(
                  *initialPassed, *finalPassed, relax_count, parametersPassed,
                  (init_opt.method == NEBInit::SIDPP_ZBL));
              break;
            case NEBInit::LINEAR:
            default:
              path = helper_functions::neb_paths::linearPath(
                  *initialPassed, *finalPassed, base_count);
              break;
            }

            // Decimate path back to the target count using the cubic spline
            if (init_opt.oversampling && path.size() > (base_count + 2)) {
              auto log = spdlog::get("combi");
              SPDLOG_LOGGER_INFO(log,
                                 "Decimating oversampled path ({} images) to "
                                 "{} images via cubic spline.",
                                 path.size() - 2, base_count);

              // Perform the Spline Resampling
              path =
                  helper_functions::neb_paths::resamplePath(path, base_count);

              // POST-DECIMATION RE-RELAXATION
              // The spline might have placed atoms in high-energy positions.
              SPDLOG_LOGGER_INFO(
                  log, "Relaxing decimated path to restore IDPP surface...");

              // Collective IDPP objective for the new reduced path
              std::shared_ptr<ObjectiveFunction> post_decim_objf =
                  std::make_shared<CollectiveIDPPObjectiveFunction>(
                      path, parametersPassed);

              // Wrap with ZBL if the original method used ZBL options
              bool use_zbl = (init_opt.method == NEBInit::SIDPP_ZBL);
              if (use_zbl) {
                auto zbl_pot =
                    helper_functions::neb_paths::createZBLPotential();
                post_decim_objf = std::make_shared<ZBLRepulsiveIDPPObjective>(
                    post_decim_objf, zbl_pot, path, parametersPassed, 1.0);
              }

              auto optim = helpers::create::mkOptim(
                  post_decim_objf, parametersPassed->neb_options.opt_method,
                  parametersPassed);

              // Run a short optimization
              optim->run(
                  parametersPassed->neb_options.initialization.max_iterations,
                  parametersPassed->neb_options.initialization.max_move);
            }
            return path;
          }(),
          parametersPassed, potPassed) {}

// Second constructor: Contains all the common setup code
NudgedElasticBand::NudgedElasticBand(
    std::vector<Matter> initPath, std::shared_ptr<Parameters> parametersPassed,
    std::shared_ptr<Potential> potPassed)
    : params{parametersPassed},
      pot{potPassed},
      E_ref{0.0},
      mmf_active{false},
      mmf_iterations_used{0} {

  log = spdlog::get("combi");
  this->status = NEBStatus::INIT;
  numImages = params->neb_options.image_count;
  atoms = initPath.front().numberOfAtoms();

  // Common initialization logic
  path.resize(numImages + 2);
  tangent.resize(numImages + 2);
  projectedForce.resize(numImages + 2);
  extremumPosition.resize(2 * (numImages + 1));
  extremumEnergy.resize(2 * (numImages + 1));
  extremumCurvature.resize(2 * (numImages + 1));
  numExtrema = 0;

  for (long i = 0; i <= numImages + 1; i++) {
    path[i] = std::make_shared<Matter>(std::move(initPath[i]));
    tangent[i] = std::make_shared<AtomMatrix>();
    tangent[i]->resize(atoms, 3);
    tangent[i]->setZero();
    projectedForce[i] = std::make_shared<AtomMatrix>();
    projectedForce[i]->resize(atoms, 3);
    projectedForce[i]->setZero();
  }

  // Common final setup
  movedAfterForceCall = true;
  path[0]->getPotentialEnergy();
  path[numImages + 1]->getPotentialEnergy();
  climbingImage = 0;

  // Setup springs
  k_u = params->neb_options.spring.weighting.k_max;
  k_l = params->neb_options.spring.weighting.k_min;
  if (params->neb_options.spring.weighting.enabled) {
    ksp = k_l;
  } else {
    ksp = params->neb_options.spring.constant;
  }

  // Optional debugging setup
  if (params->estNEBeig) {
    eigenmode_solvers.resize(numImages + 2);
    for (long i = 0; i <= numImages + 1; i++) {
      if (params->nebMMF == LowestEigenmode::MINMODE_DIMER) {
        eigenmode_solvers[i] =
            std::make_shared<ImprovedDimer>(path[i], parametersPassed, pot);
      } else if (params->nebMMF == LowestEigenmode::MINMODE_LANCZOS) {
        eigenmode_solvers[i] =
            std::make_shared<Lanczos>(path[i], parametersPassed, pot);
      } else {
        log = spdlog::get("_traceback");
        SPDLOG_LOGGER_CRITICAL(log, "[Debug] unknown neb_mmf_estimator: {}",
                               params->nebMMF);
        std::exit(1);
      }
    }
  }
}

NudgedElasticBand::NEBStatus NudgedElasticBand::compute(void) {
  if (current_mmf_threshold < 0) {
    current_mmf_threshold =
        params->neb_options.climbing_image.roneb.trigger_force;
  }
  long iteration = 0;
  this->status = NEBStatus::RUNNING;
  mmf_active = false;
  mmf_iterations_used = 0;

  long previousClimbingImage = -1;
  int ciStabilityCounter = 0;
  // Number of stable iterations required
  const int stability_threshold =
      params->neb_options.climbing_image.roneb.ci_stability_count;

  SPDLOG_LOGGER_DEBUG(log, "Nudged elastic band calculation started.");

  // Initialize E_ref for energy weighting
  E_ref = std::min(path[0]->getPotentialEnergy(),
                   path[numImages + 1]->getPotentialEnergy());

  updateForces();

  auto objf = std::make_shared<NEBObjectiveFunction>(this, params);

  bool switched{false};
  // TODO(rg): clear up the refine stanza with it's own opt method, use pre_post
  auto optim =
      helpers::create::mkOptim(objf, params->neb_options.opt_method, params);
  std::unique_ptr<Optimizer> refine_optim{nullptr};
  if (params->refineOptMethod != OptType::None) {
    refine_optim =
        helpers::create::mkOptim(objf, params->refineOptMethod, params);
  }
  while (this->status != NEBStatus::GOOD) {
    if (params->writeMovies) {
      bool append = (iteration != 0);
      path[maxEnergyImage]->matter2con("neb_maximage.con", append);
      std::string nebFilename(fmt::format("neb_path_{:03d}.con", iteration));
      FILE *fileNEBPath = fopen(nebFilename.c_str(), "wb");
      for (long idx = 0; idx <= numImages + 1; idx++) {
        path[idx]->matter2con(fileNEBPath);
      }
      fclose(fileNEBPath);
      printImageData(true, iteration);
    }

    VectorXd pos = objf->getPositions();
    double convForce = convergenceForce();

    // convergenceForce() calls updateForces(), which updates 'climbingImage'
    if (climbingImage == previousClimbingImage) {
      ciStabilityCounter++;
    } else {
      ciStabilityCounter = 0;
      previousClimbingImage = climbingImage;
    }
    // -------------------------------------

    if (iteration == 0) {
      baseline_force = convForce;
      // Initialize MMF threshold based on baseline
      current_mmf_threshold =
          baseline_force *
          params->neb_options.climbing_image.roneb.trigger_factor;

      // --- Improved Logging Start ---
      auto &ci_opt = params->neb_options.climbing_image;
      auto &mmf_opt = ci_opt.roneb;

      // Helper to format infinity
      auto fmt_trigger = [](double val) -> std::string {
        if (val > 1e100)
          return "INF";
        return fmt::format("{:.4f}", val);
      };

      SPDLOG_LOGGER_INFO(
          log,
          "===============================================================");
      SPDLOG_LOGGER_INFO(log, " NEB Optimization Configuration");
      SPDLOG_LOGGER_INFO(
          log,
          "===============================================================");
      SPDLOG_LOGGER_INFO(log, " {:<25} : {:.4f}", "Baseline Force",
                         baseline_force);

      // Climbing Image Logs
      std::string ci_status = ci_opt.enabled ? "ENABLED" : "DISABLED";
      SPDLOG_LOGGER_INFO(log, " {:<25} : {}", "Climbing Image (CI)", ci_status);
      if (ci_opt.enabled) {
        double ci_rel_val = baseline_force * ci_opt.trigger_factor;
        SPDLOG_LOGGER_INFO(log, "   - {:<21} : {} (Factor: {:.2f})",
                           "Relative Trigger", fmt_trigger(ci_rel_val),
                           ci_opt.trigger_factor);
        SPDLOG_LOGGER_INFO(log, "   - {:<21} : {}", "Absolute Trigger",
                           fmt_trigger(ci_opt.trigger_force));
        SPDLOG_LOGGER_INFO(log, "   - {:<21} : {}", "Converged Only",
                           ci_opt.converged_only);
      }

      // RONEB / MMF Logs
      std::string mmf_status =
          (ci_opt.enabled && mmf_opt.use_mmf) ? "ENABLED" : "DISABLED";
      SPDLOG_LOGGER_INFO(log, " {:<25} : {}", "Hybrid MMF (RONEB)", mmf_status);
      if (ci_opt.enabled && mmf_opt.use_mmf) {
        SPDLOG_LOGGER_INFO(log, "   - {:<21} : {:.4f} (Factor: {:.2f})",
                           "Initial Threshold", current_mmf_threshold,
                           mmf_opt.trigger_factor);
        SPDLOG_LOGGER_INFO(log, "   - {:<21} : {:.4f}", "Absolute Floor",
                           mmf_opt.trigger_force);
        SPDLOG_LOGGER_INFO(log,
                           "   - {:<21} : {:.2f} (Base: {:.2f}, Str: {:.2f})",
                           "Penalty Scheme", mmf_opt.penalty.base,
                           mmf_opt.penalty.base, mmf_opt.penalty.strength);
        SPDLOG_LOGGER_INFO(log, "   - {:<21} : {:.4f}", "Angle Tolerance",
                           mmf_opt.angle_tol);
      }
      SPDLOG_LOGGER_INFO(
          log,
          "---------------------------------------------------------------");
      // --- Improved Logging End ---
    }

    if (iteration == 0) {

      SPDLOG_DEBUG("{:>10s} {:>12s} {:>14s} {:>11s} {:>12s}", "iteration",
                   "step size", params->optConvergenceMetricLabel, "max image",
                   "max energy");
      SPDLOG_DEBUG(
          "---------------------------------------------------------------\n");
    }

    // CI active when force drops below relative threshold
    bool ci_active =
        params->neb_options.climbing_image.enabled &&
        (convForce < baseline_force *
                         params->neb_options.climbing_image.trigger_factor ||
         convForce < params->neb_options.climbing_image.trigger_force);

    if (iteration) {
      // MMF triggering
      if (params->neb_options.climbing_image.roneb.use_mmf && ci_active &&
          climbingImage > 0 && climbingImage <= numImages &&
          (ciStabilityCounter > stability_threshold) &&
          (convForce < current_mmf_threshold ||
           convForce <
               params->neb_options.climbing_image.roneb.trigger_force) &&
          convForce > params->neb_options.force_tolerance) {

        SPDLOG_LOGGER_DEBUG(log,
                            "Triggering MMF.  Force: {:.4f}, Threshold: {:.4f} "
                            "({:.2f}x baseline)",
                            convForce, current_mmf_threshold,
                            current_mmf_threshold / baseline_force);

        // Save climbing image state before MMF
        AtomMatrix savedPositions = path[climbingImage]->getPositions();

        // Run MMF and get alignment info
        double alignment = 0.0;
        int mmfResult = runMMFRefinement(alignment);

        // Always update forces after MMF - it moved the climbing image
        movedAfterForceCall = true;
        updateForces();
        double newForce = convergenceForce();

        // Check if we're done
        if (newForce < params->neb_options.force_tolerance) {
          SPDLOG_LOGGER_DEBUG(log, "NEB converged after MMF. Force: {:.4f}",
                              newForce);
          status = NEBStatus::GOOD;
          break;
        }

        // Evaluate whether MMF helped, regardless of its internal status
        bool mmfHelped = (newForce < convForce) && mmfResult != -2;

        if (mmfHelped) {
          // MMF made progress - set threshold relative to baseline
          double progress_ratio = newForce / baseline_force;
          // Allow retry when force drops to ~50-90% of current level
          current_mmf_threshold =
              newForce * (0.5 + 0.4 * (newForce / convForce));

          // But never above the initial MMF trigger
          double max_threshold =
              baseline_force *
              params->neb_options.climbing_image.roneb.trigger_factor;
          current_mmf_threshold =
              std::min(current_mmf_threshold, max_threshold);

          SPDLOG_LOGGER_DEBUG(log,
                              "MMF helped (status={}). Force: {:.4f} -> {:.4f} "
                              "({:.2f}x baseline). New threshold: {:.4f}",
                              mmfResult, convForce, newForce,
                              newForce / baseline_force, current_mmf_threshold);
        } else {
          // MMF didn't help - apply backoff relative to baseline
          double strength =
              params->neb_options.climbing_image.roneb.penalty.strength;
          double base = params->neb_options.climbing_image.roneb.penalty.base;

          // Ensure alignment stays within a physically meaningful [0, 1] range
          double alpha = std::clamp(alignment, 0.0, 1.0);

          // penalty_factor:  how much of the baseline trigger to use
          // Good alignment (1.0) -> factor approaches 1.0 (use full trigger)
          // Poor alignment (0.0) -> factor approaches 'base' (reduced trigger)
          double penalty_factor =
              base + (1.0 - base) * std::pow(alpha, strength);

          // Clamp the factor to ensure the threshold never increases
          penalty_factor = std::clamp(penalty_factor, base, 1.0);

          // New threshold as fraction of baseline
          double trigger_factor =
              params->neb_options.climbing_image.roneb.trigger_factor;
          current_mmf_threshold =
              baseline_force * trigger_factor * penalty_factor;

          // Floor at 2x force tolerance
          double min_threshold = params->neb_options.force_tolerance * 2.0;
          current_mmf_threshold =
              std::max(current_mmf_threshold, min_threshold);

          SPDLOG_LOGGER_DEBUG(
              log,
              "MMF backoff (status={}). Force: {:.4f} -> {:.4f}, "
              "Alignment:  {:.3f}. New threshold: {:.4f} ({:.2f}x baseline)",
              mmfResult, convForce, newForce, alignment, current_mmf_threshold,
              current_mmf_threshold / baseline_force);

          // NOTE(rg): Doesn't really help anyway
          // Only revert if MMF made things much worse
          // if (newForce > convForce * 2) {
          //   SPDLOG_LOGGER_DEBUG(log, "Reverting to older position.");
          //   path[climbingImage]->setPositions(savedPositions);
          // }
        }
        if ((savedPositions - path[climbingImage]->getPositions()).norm() >
            params->optMaxMove * params->neb_options.image_count) {
          // Reset optimizer after MMF - history is stale
          SPDLOG_LOGGER_DEBUG(log, "Resetting optimization history.");
          optim = helpers::create::mkOptim(objf, params->neb_options.opt_method,
                                           params);
        }

        // Reset stability counter since MMF moved the image
        ciStabilityCounter = 0;
      }

      if (iteration >= params->neb_options.max_iterations) {
        status = NEBStatus::BAD_MAX_ITERATIONS;
        break;
      }

      // Take optimizer step with CI enabled/disabled based on threshold
      bool originalCIflag = params->neb_options.climbing_image.enabled;
      params->neb_options.climbing_image.enabled = ci_active;

      if (refine_optim) {
        if (optim && convForce > params->refineThreshold) {
          optim->step(params->optMaxMove);
        } else {
          if (!switched) {
            switched = true;
            SPDLOG_DEBUG("Switched to {}", magic_enum::enum_name<OptType>(
                                               params->refineOptMethod));
          }
          refine_optim->step(params->optMaxMove);
        }
      } else {
        optim->step(params->optMaxMove);
      }

      params->neb_options.climbing_image.enabled = originalCIflag;
    }

    iteration++;

    double dE = path[maxEnergyImage]->getPotentialEnergy() -
                path[0]->getPotentialEnergy();
    double stepSize = helper_functions::maxAtomMotionV(
        path[0]->pbcV(objf->getPositions() - pos));
    SPDLOG_LOGGER_DEBUG(log, "{:>10} {:>12.4e} {:>14.4e} {:>11} {:>12.4}",
                        iteration, stepSize, convergenceForce(), maxEnergyImage,
                        dE);

    if (pot->getType() == PotType::CatLearn) {
      if (objf->isUncertain()) {
        SPDLOG_LOGGER_DEBUG(log, "NEB failed due to high uncertainity");
        status = NEBStatus::MAX_UNCERTAINITY;
        break;
      } else if (objf->isConverged()) {
        SPDLOG_LOGGER_DEBUG(log, "NEB converged\n");
        status = NEBStatus::GOOD;
        break;
      }
    } else {
      if (objf->isConverged()) {
        SPDLOG_LOGGER_DEBUG(log, "NEB converged\n");
        status = NEBStatus::GOOD;
        break;
      }
    }
  }
  return status;
}

// Modified signature to return alignment
int NudgedElasticBand::runMMFRefinement(double &alignment) {
  alignment = 0.0;

  if (climbingImage <= 0 || climbingImage > numImages) {
    SPDLOG_LOGGER_WARN(log, "Invalid climbing image for MMF:  {}",
                       climbingImage);
    return -1;
  }

  // Get the current tangent as initial mode guess
  AtomMatrix initialMode = *tangent[climbingImage];

  // Verify tangent is valid
  double tangentNorm = initialMode.norm();
  if (tangentNorm < 1e-8) {
    SPDLOG_LOGGER_WARN(log, "Tangent too small for MMF initialization");
    return -1;
  }
  initialMode /= tangentNorm;

  auto tempMinModeSearch = std::make_shared<MinModeSaddleSearch>(
      path[climbingImage], initialMode,
      path[climbingImage]->getPotentialEnergy(), params, pot);

  // Save and override saddle search parameters
  auto originalSaddleMaxIterations = params->saddleMaxIterations;

  params->saddleMaxIterations =
      params->neb_options.climbing_image.roneb.max_steps;

  int minModeStatus;
  try {
    minModeStatus = tempMinModeSearch->run();
  } catch (const eonc::DimerModeRestoredException &e) {
    // Dimer restored to valid state - treat as partial success
    minModeStatus = MinModeSaddleSearch::STATUS_DIMER_RESTORED_BEST;
    SPDLOG_LOGGER_DEBUG(log, "MMF:  Dimer restored to best state");
  } catch (const eonc::DimerModeLostException &e) {
    minModeStatus = MinModeSaddleSearch::STATUS_DIMER_LOST_MODE;
    SPDLOG_LOGGER_WARN(log, "Dimer lost mode during MMF refinement");
  }

  // Restore original parameters
  params->saddleMaxIterations = originalSaddleMaxIterations;

  // Track iterations used
  mmf_iterations_used += tempMinModeSearch->iteration;

  double eigenvalue = tempMinModeSearch->getEigenvalue();

  // If the surface is convex (positive curvature), MMF will drag us to a
  // minimum. We only want to refine if we are in a saddle region (negative
  // curvature).
  if (eigenvalue > 0.0) {
    SPDLOG_LOGGER_WARN(log,
                       "MMF skipped: Positive curvature detected (eig={:.4f}).",
                       eigenvalue);
    return -2;
  }

  // Calculate alignment for all outcomes
  AtomMatrix finalModeMatrix = tempMinModeSearch->getEigenvector();
  VectorXd finalMode = VectorXd::Map(finalModeMatrix.data(), 3 * atoms);
  VectorXd currentTangent =
      VectorXd::Map(tangent[climbingImage]->data(), 3 * atoms);
  alignment = std::abs(finalMode.normalized().dot(currentTangent.normalized()));

  if (minModeStatus == MinModeSaddleSearch::STATUS_GOOD ||
      minModeStatus == MinModeSaddleSearch::STATUS_DIMER_RESTORED_BEST) {
    // Check alignment - if mode drifted too far, treat as failure
    if (alignment < params->neb_options.climbing_image.roneb.angle_tol) {
      SPDLOG_LOGGER_WARN(
          log,
          "MMF converged/restored but mode drifted (alignment={:.3f} < {:.3f})",
          alignment, params->neb_options.climbing_image.roneb.angle_tol);
      return -1;
    }
    return 0;
  } else if (minModeStatus == MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS) {
    return 1;
  } else {
    SPDLOG_LOGGER_WARN(log, "MMF failed.  Mode-tangent alignment: {:.3f}",
                       alignment);
    return -1;
  }
}

// generate the force value that is compared to the convergence criterion
double NudgedElasticBand::convergenceForce(void) {
  if (movedAfterForceCall)
    updateForces();
  double fmax = 0;

  for (long i = 1; i <= numImages; i++) {

    if (params->neb_options.climbing_image.converged_only == true &&
        params->neb_options.climbing_image.enabled && climbingImage != 0) {
      i = climbingImage;
    }
    if (params->optConvergenceMetric == "norm") {
      fmax = max(fmax, projectedForce[i]->norm());
    } else if (params->optConvergenceMetric == "max_atom") {
      for (int j = 0; j < path[0]->numberOfAtoms(); j++) {
        if (path[0]->getFixed(j))
          continue;
        fmax = max(fmax, projectedForce[i]->row(j).norm());
      }
    } else if (params->optConvergenceMetric == "max_component") {
      fmax = max(fmax, projectedForce[i]->maxCoeff());
    } else {
      log = spdlog::get("_traceback");
      SPDLOG_LOGGER_CRITICAL(
          log, "[Nudged Elastic Band] unknown opt_convergence_metric: {}",
          params->optConvergenceMetric);
      std::exit(1);
    }
    if (params->neb_options.climbing_image.converged_only == true &&
        params->neb_options.climbing_image.enabled && climbingImage != 0) {
      break;
    }
  }
  return fmax;
}

// Update the forces, do the projections, and add spring forces
void NudgedElasticBand::updateForces(void) {
  // variables for tangent
  double maxDiffEnergy, minDiffEnergy;
  double energyDiffPrev, energyDiffNext;
  double energy, energyPrev, energyNext;

  // variables for climbing image
  double maxEnergy;

  // variables for the energy weighted springs
  k_l = params->neb_options.spring.weighting.k_min;
  k_u = params->neb_options.spring.weighting.k_max;
  std::vector<double> springConstants(numImages + 2, k_l);

  // variables for force projections
  AtomMatrix force(atoms, 3), forcePerp(atoms, 3), forcePar(atoms, 3);
  AtomMatrix forceSpringPar(atoms, 3), forceSpring(atoms, 3),
      forceSpringPerp(atoms, 3);
  AtomMatrix forceDNEB(atoms, 3);
  AtomMatrix pos(atoms, 3), posNext(atoms, 3), posPrev(atoms, 3),
      posDiffNext(atoms, 3), posDiffPrev(atoms, 3);
  double distNext, distPrev;

  // Update forces for all intermediate images FIRST
  for (long i = 1; i <= numImages; i++) {
    path[i]->getForces();
  }

  // Find the highest energy non-endpoint image
  auto first = path.begin() + 1;
  auto last = path.begin() + numImages + 1;
  auto it = std::max_element(
      first, last,
      [](const std::shared_ptr<Matter> &a, const std::shared_ptr<Matter> &b) {
        return a->getPotentialEnergy() < b->getPotentialEnergy();
      });
  maxEnergyImage = std::distance(path.begin(), it);
  maxEnergy = (*it)->getPotentialEnergy();

  // Update E_ref for energy weighting (use current endpoint energies)
  if (params->neb_options.spring.weighting.enabled) {
    E_ref = std::min(path[0]->getPotentialEnergy(),
                     path[numImages + 1]->getPotentialEnergy());
  }

  // Reset climbing image if CI is disabled
  if (!params->neb_options.climbing_image.enabled) {
    climbingImage = 0;
  }

  // Onsager-Machlup (OM) Action Logic
  // Mandelli & Parrinello (2021)
  // L_i = (1 / 2k) * Force_i
  bool omActive = params->neb_options.spring.om.enabled;
  const VectorXd &atomMasses = path[0]->getMasses();
  std::vector<AtomMatrix> L_vecs;

  if (omActive) {
    // Reuse the standard spring constant as k_OM
    double base_k = params->neb_options.spring.constant;
    if (params->neb_options.spring.om.optimize_k) {
      // Method A: Physical estimation (assuming v*dt=1)
      // base_k = m / (2 * dt^2). With v*dt=1, v = 1/dt.
      // base_k = m * (1/dt) / (2 * dt) = m / (2 * dt^2).
      // This requires a valid timestep 'dt'.
      // If dt is not reliable, we use Method B (Force Balancing).

      // Method B: Force Balancing (Heuristic)
      // We want F_spring ~ F_potential on average.
      // F_spring_OM ~ base_k * <|R_{i+1} + R_{i-1} - 2R_i|>
      // F_potential ~ <|Force|>
      // base_k = <|Force|> / <|2nd_derivative_finite_diff|>

      double avgPotForce = 0.0;
      double avgPathCurvature = 0.0;
      int count = 0;

      for (long j = 1; j <= numImages; j++) {
        avgPotForce += path[j]->getForces().norm();

        // Calculate finite difference curvature (R_{i+1} + R_{i-1} - 2R_i)
        // Need to handle PBC carefully
        AtomMatrix next = path[j + 1]->getPositions();
        AtomMatrix prev = path[j - 1]->getPositions();
        AtomMatrix curr = path[j]->getPositions();

        AtomMatrix curvVec = path[j]->pbc(next + prev - 2.0 * curr);
        avgPathCurvature += curvVec.norm();
        count++;
      }

      if (count > 0 && avgPathCurvature > 1e-6) { // Avoid division by zero
        avgPotForce /= count;
        avgPathCurvature /= count;

        // Scale factor can be a tuning parameter, default 1.0
        double scale = params->neb_options.spring.om.k_scale;
        base_k = scale * (avgPotForce / avgPathCurvature);

        // Optional: Clamp base_k to reasonable bounds to prevent instability
        base_k = std::max(base_k, params->neb_options.spring.om.k_min);
        base_k = std::min(base_k, params->neb_options.spring.om.k_max);
      }

      // Log the optimized k for debugging
      // TODO(rg): refactor and log separately
      // SPDLOG_LOGGER_DEBUG(log, "Optimized OM k: {}", base_k);
    }
    // Pre-calculate L vectors for all images (including endpoints)
    L_vecs.resize(numImages + 2);

    for (long j = 0; j <= numImages + 1; j++) {
      L_vecs[j].resize(atoms, 3);
      // Endpoints might not have updated forces if fixed, assume 0 (minima)
      if (j == 0 || j == numImages + 1) {
        L_vecs[j].setZero();
      } else {
        AtomMatrix forces = path[j]->getForces();
        for (int k = 0; k < atoms; k++) {
          // Apply 1/m scaling as per Eq (12)
          // Note: base_k here represents (nu / 2 dt).
          // L = (1 / (2 * base_k * m)) * Force
          double m = atomMasses(k);
          if (m > 1e-8) {
            double alpha_k = 1.0 / (2.0 * base_k * m);
            L_vecs[j].row(k) = alpha_k * forces.row(k);
          } else {
            L_vecs[j].row(k).setZero();
          }
        }
      }
    }
  } else if (params->neb_options.spring.weighting.enabled) {
    // Energy weighted springs
    // Protect against division by zero
    double energyRange = maxEnergy - E_ref;
    if (energyRange < 1e-10) {
      // All images at same energy - use uniform springs
      std::fill(springConstants.begin(), springConstants.end(), k_l);
    } else {
      for (int idx = 1; idx <= numImages + 1; idx++) {
        double Ei = std::max(path[idx]->getPotentialEnergy(),
                             path[idx - 1]->getPotentialEnergy());
        if (Ei > E_ref) {
          double alpha_i = (maxEnergy - Ei) / energyRange;
          // Clamp alpha to [0, 1] for numerical safety
          alpha_i = std::max(0.0, std::min(1.0, alpha_i));
          springConstants[idx - 1] = (1.0 - alpha_i) * k_u + alpha_i * k_l;
        } else {
          springConstants[idx - 1] = k_l;
        }
      }
    }
  } else {
    // If CI is not yet active, use the standard uniform spring constant
    std::fill(springConstants.begin(), springConstants.end(),
              params->neb_options.spring.constant);
  }

  // Projection Loop
  for (long i = 1; i <= numImages; i++) {
    // Zero out temporary force matrices
    forceSpring.setZero();
    forceSpringPar.setZero();
    forceSpringPerp.setZero();
    forceDNEB.setZero();

    // set local variables
    force = path[i]->getForces();
    pos = path[i]->getPositions();
    posPrev = path[i - 1]->getPositions();
    posNext = path[i + 1]->getPositions();
    energy = path[i]->getPotentialEnergy();
    energyPrev = path[i - 1]->getPotentialEnergy();
    energyNext = path[i + 1]->getPotentialEnergy();
    posDiffNext = path[i]->pbc(posNext - pos); // R[i+1] - R[i]
    posDiffPrev = path[i]->pbc(pos - posPrev); // R[i] - R[i-1]
    distNext = posDiffNext.norm();             // Distance to next image
    distPrev = posDiffPrev.norm();             // Distance to previous image

    // determine the tangent
    if (params->neb_options.climbing_image.use_old_tangent) {
      // old tangent
      *tangent[i] = posDiffNext;
    } else {
      // new improved tangent
      // higherEnergyPrev = energyPrev > energyNext;
      // higherEnergyNext = energyNext > energyPrev;

      if (energyNext > energy && energy > energyPrev) {
        *tangent[i] = posDiffNext;
      } else if (energy > energyNext && energyPrev > energy) {
        *tangent[i] = posDiffPrev;
      } else {
        // we are at an extremum
        energyDiffPrev = energyPrev - energy;
        energyDiffNext = energyNext - energy;

        // calculate the energy difference to neighboring numImages
        minDiffEnergy = min(abs(energyDiffPrev), abs(energyDiffNext));
        maxDiffEnergy = max(abs(energyDiffPrev), abs(energyDiffNext));

        // use these energy differences to weight the tangent
        if (energyDiffPrev > energyDiffNext) {
          *tangent[i] = posDiffNext * minDiffEnergy;
          *tangent[i] += posDiffPrev * maxDiffEnergy;
        } else {
          *tangent[i] = posDiffNext * maxDiffEnergy;
          *tangent[i] += posDiffPrev * minDiffEnergy;
        }
      }
    }

    // Normalize tangent with safety check
    double tangentNorm = tangent[i]->norm();
    if (tangentNorm > 1e-10) {
      *tangent[i] /= tangentNorm;
    } else {
      // Fallback:  use direction to next image
      *tangent[i] = posDiffNext;
      tangentNorm = tangent[i]->norm();
      if (tangentNorm > 1e-10) {
        *tangent[i] /= tangentNorm;
      }
    }

    // Calculate the force perpendicular to the tangent
    forcePerp =
        force - (force.array() * (*tangent[i]).array()).sum() * *tangent[i];

    // Calculate spring forces
    if (omActive) {
      double base_k = params->neb_options.spring.constant;
      if (params->neb_options.spring.om.optimize_k) {
        // This is recalculated unnecessarily here, could be cached if
        // optimizing k But for clarity and following the logic:
        double avgPotForce = 0.0;
        double avgPathCurvature = 0.0;
        int count = 0;
        for (long j = 1; j <= numImages; j++) {
          avgPotForce += path[j]->getForces().norm();
          AtomMatrix next = path[j + 1]->getPositions();
          AtomMatrix prev = path[j - 1]->getPositions();
          AtomMatrix curr = path[j]->getPositions();
          AtomMatrix curvVec = path[j]->pbc(next + prev - 2.0 * curr);
          avgPathCurvature += curvVec.norm();
          count++;
        }
        if (count > 0 && avgPathCurvature > 1e-6) {
          double scale = params->neb_options.spring.om.k_scale;
          base_k = scale * (avgPotForce / avgPathCurvature);
          base_k = std::max(base_k, params->neb_options.spring.om.k_min);
          base_k = std::min(base_k, params->neb_options.spring.om.k_max);
        }
      }

      // Mandelli Eq. 13: k * ( R(i+1) + R(i-1) - 2R(i) + L(i-1) - L(i) )
      AtomMatrix diff = path[i]->pbc(posNext + posPrev - 2.0 * pos +
                                     L_vecs[i + 1] - L_vecs[i]);

      for (int k = 0; k < atoms; k++) {
        diff.row(k) *= atomMasses(k);
      }

      // Now multiply by the global tuning constant
      AtomMatrix f_om_vec = base_k * diff;

      // Mandelli Eq. 15: Project onto tangent
      forceSpringPar =
          (f_om_vec.array() * (*tangent[i]).array()).sum() * *tangent[i];

    } else if (params->neb_options.spring.weighting.enabled) {
      // Use energy-weighted constants if trigger is met
      // Spring for segment (i) -> (i+1)
      double kspNext = springConstants[i];
      // Spring for segment (i-1) -> (i)
      double kspPrev = springConstants[i - 1];
      forceSpringPar =
          ((kspNext * distNext) - (kspPrev * distPrev)) * *tangent[i];
    } else {
      // Use uniform spring constant during initial relaxation
      this->ksp = params->neb_options.spring.constant;
      forceSpringPar = this->ksp * (distNext - distPrev) * *tangent[i];
      forceSpring = this->ksp * path[i]->pbc((posNext - pos) - (pos - posPrev));
    }

    // DNEB forces
    if (params->neb_options.spring.doubly_nudged && !omActive &&
        !params->neb_options.spring.weighting.enabled) {
      forceSpringPerp =
          forceSpring -
          (forceSpring.array() * (*tangent[i]).array()).sum() * *tangent[i];

      double forceSpringPerpNorm = forceSpringPerp.norm();
      double forcePerpNorm = forcePerp.norm();

      if (forceSpringPerpNorm > 1e-10 && forcePerpNorm > 1e-10) {
        AtomMatrix forcePerpNormalized = forcePerp / forcePerpNorm;
        forceDNEB =
            forceSpringPerp -
            (forceSpringPerp.array() * forcePerpNormalized.array()).sum() *
                forcePerpNormalized;

        double switching =
            2.0 / M_PI *
            atan(pow(forcePerpNorm, 2) / pow(forceSpringPerpNorm, 2));
        forceDNEB *= switching;
      } else {
        // Force norms too small for stable DNEB calculation
        forceDNEB.setZero();
      }
    } else {
      forceDNEB.setZero();
    }

    // Apply the Climbing Image projection or standard NEB force
    // Note: CI only activates if both the parameter is enabled AND the trigger
    // is met
    if (params->neb_options.climbing_image.enabled && i == maxEnergyImage) {
      climbingImage = maxEnergyImage;
      // Climbing image:  invert force component along tangent
      *projectedForce[i] =
          force -
          (2.0 * (force.array() * (*tangent[i]).array()).sum() * *tangent[i]) +
          forceDNEB;
    } else { // non-climbing images
      // Regular image
      // sum the spring and potential forces for the neb force
      if (params->neb_options.spring.use_elastic_band && !omActive &&
          !params->neb_options.spring.weighting.enabled) {
        *projectedForce[i] = forceSpring + force;
      } else {
        *projectedForce[i] = forceSpringPar + forcePerp + forceDNEB;
      }
    }

    // Zero net translational force (if all atoms are free)
    if (path[i]->numberOfFreeAtoms() == path[i]->numberOfAtoms()) {
      for (int j = 0; j <= 2; j++) {
        double translationMag = projectedForce[i]->col(j).sum();
        int natoms = projectedForce[i]->col(j).size();
        projectedForce[i]->col(j).array() -= translationMag / ((double)natoms);
      }
    }
  }

  // Flag that forces are fresh
  movedAfterForceCall = false;
  return;
}

// Print NEB image data
void NudgedElasticBand::printImageData(bool writeToFile, size_t idx) {
  double dist, distTotal = 0;
  AtomMatrix tangentStart =
      path[0]->pbc(path[1]->getPositions() - path[0]->getPositions());
  AtomMatrix tangentEnd = path[numImages]->pbc(
      path[numImages + 1]->getPositions() - path[numImages]->getPositions());
  AtomMatrix tang;
  std::string header, fmttr;
  if (params->estNEBeig) {
    header = fmt::format("{:>3s} {:>12s} {:>12s} {:>12s} {:>12s}", "img",
                         "rxn_coord", "energy", "f_para", "eigval");
    fmttr = "{:>3} {:>12.6f} {:>12.6f} {:>12.6f} {:>12.6f}";
  } else {
    header = fmt::format("{:>3s} {:>12s} {:>12s} {:>12s}", "img", "rxn_coord",
                         "energy", "f_para");
    fmttr = "{:>3} {:>12.6f} {:>12.6f} {:>12.6f}";
  }

  // Endpoint tangents must be normalized for a correct projection
  tangentStart.normalize();
  tangentEnd.normalize();

  std::shared_ptr<spdlog::logger> fileLogger;
  if (writeToFile) {
    if (spdlog::get("file_logger")) {
      spdlog::drop("file_logger");
    }
    // If idx is SIZE_MAX, treat this as the dedicated final output file
    // "neb.dat". This avoids overwriting iteration files named neb_000.dat,
    // neb_001.dat, ...
    std::string neb_dat_fs;
    if (idx == std::numeric_limits<size_t>::max()) {
      neb_dat_fs = "neb.dat";
    } else {
      neb_dat_fs = fmt::format("neb_{:03}.dat", idx);
    }
    if (fs::exists(neb_dat_fs)) {
      fs::remove(neb_dat_fs);
    }
    fileLogger = spdlog::basic_logger_mt("file_logger", neb_dat_fs);
    fileLogger->set_pattern("%v");
    fileLogger->info(header);
  }
  const double energy_reactant = path[0]->getPotentialEnergy();

  for (long i = 0; i <= numImages + 1; i++) {
    if (i == 0) {
      tang = tangentStart;
    } else if (i == numImages + 1) {
      tang = tangentEnd;
    } else {
      tang = *tangent[i];
    }

    if (i > 0) {
      dist = path[i]->distanceTo(*path[i - 1]);
      distTotal += dist;
    }

    // Calculate standard values for logging
    double relative_energy = path[i]->getPotentialEnergy() - energy_reactant;
    double parallel_force = (path[i]->getForces().array() * tang.array()).sum();

    // Conditionally run an MMF and log the eigenvalue
    if (params->estNEBeig) {
      eigenmode_solvers[i]->compute(path[i], tang);
      double lowest_eigenvalue = eigenmode_solvers[i]->getEigenvalue();
      if (fileLogger) {
        fileLogger->info(fmttr, i, distTotal, relative_energy, parallel_force,
                         lowest_eigenvalue);
      } else {
        SPDLOG_LOGGER_DEBUG(log, fmttr, i, distTotal, relative_energy,
                            parallel_force, lowest_eigenvalue);
      }
    } else { // Standard output without the eigenvalue
      if (fileLogger) {
        fileLogger->info(fmttr, i, distTotal, relative_energy, parallel_force);
      } else {
        SPDLOG_LOGGER_DEBUG(log, fmttr, i, distTotal, relative_energy,
                            parallel_force);
      }
    }
  }

  if (fileLogger) {
    spdlog::drop("file_logger");
  }
}

// Estimate the barrier using a cubic spline
void NudgedElasticBand::findExtrema(void) {
  // calculate the cubic parameters for each interval (a,b,c,d)

  AtomMatrix tangentEndpoint;
  double a[numImages + 1], b[numImages + 1], c[numImages + 1], d[numImages + 1];
  double F1, F2, U1, U2, dist;

  for (long i = 0; i <= numImages; i++) {
    dist = path[i]->distanceTo(*path[i + 1]);
    if (i == 0) {
      tangentEndpoint =
          path[i]->pbc(path[1]->getPositions() - path[0]->getPositions());
      tangentEndpoint.normalize();
      F1 =
          (path[i]->getForces().array() * tangentEndpoint.array()).sum() * dist;
    } else {
      F1 = (path[i]->getForces().array() * (*tangent[i]).array()).sum() * dist;
    }
    if (i == numImages) {
      tangentEndpoint = path[i + 1]->pbc(path[numImages + 1]->getPositions() -
                                         path[numImages]->getPositions());
      tangentEndpoint.normalize();
      F2 = (path[i + 1]->getForces().array() * tangentEndpoint.array()).sum() *
           dist;
    } else {
      F2 =
          (path[i + 1]->getForces().array() * (*tangent[i + 1]).array()).sum() *
          dist;
    }
    U1 = path[i]->getPotentialEnergy();
    U2 = path[i + 1]->getPotentialEnergy();
    a[i] = U1;
    b[i] = -F1;
    c[i] = 3. * (U2 - U1) + 2. * F1 + F2;
    d[i] = -2. * (U2 - U1) - (F1 + F2);
  }

  // finding extrema along the MEP

  //    long numExtrema = 0;
  //    double extremaEnergy[2*(numImages+1)]; // the maximum number of
  //    extrema double extremaPosition[2*(numImages+1)];
  double discriminant, f;

  for (long i = 0; i <= numImages; i++) {
    discriminant = pow(c[i], 2) - 3. * b[i] * d[i];
    if (discriminant >= 0) {
      f = -1;

      // quadratic case
      if ((d[i] == 0) && (c[i] != 0)) {
        f = (-b[i] / (2. * c[i]));
      }
      // cubic case 1
      else if (d[i] != 0) {
        f = -(c[i] + sqrt(discriminant)) / (3. * d[i]);
      }
      if ((f >= 0) && (f <= 1)) {
        extremumPosition[numExtrema] = i + f;
        extremumEnergy[numExtrema] =
            d[i] * pow(f, 3) + c[i] * pow(f, 2) + b[i] * f + a[i];
        extremumCurvature[numExtrema] = 6.0 * d[i] * f + 2 * c[i];
        numExtrema++;
      }
      // cubic case 2
      if (d[i] != 0) {
        f = (-(c[i] - sqrt(discriminant)) / (3. * d[i]));
      }
      if ((f >= 0) && (f <= 1)) {
        extremumPosition[numExtrema] = i + f;
        extremumEnergy[numExtrema] =
            d[i] * pow(f, 3) + c[i] * pow(f, 2) + b[i] * f + a[i];
        extremumCurvature[numExtrema] = 6 * d[i] * f + 2 * c[i];
        numExtrema++;
      }
    }
  }

  SPDLOG_LOGGER_DEBUG(log, "Found {} extrema", numExtrema);
  SPDLOG_LOGGER_DEBUG(log, "Energy reference: {}",
                      path[0]->getPotentialEnergy());
  for (long i = 0; i < numExtrema; i++) {
    SPDLOG_LOGGER_DEBUG(
        log, "extrema #{} at image position {} with energy {} and curvature {}",
        i + 1, extremumPosition[i],
        extremumEnergy[i] - path[0]->getPotentialEnergy(),
        extremumCurvature[i]);
  }
}
