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
#include "eon/NudgedElasticBand.h"
#include "eon/BaseStructures.h"
#include "eon/EigenmodeStrategy.h"
#include "eon/IDPPObjectiveFunction.hpp"
#include "eon/IRACompare.h"
#include "eon/NEBForceProjection.h"
#include "eon/NEBInitialPaths.hpp"
#include "eon/NEBOcinebController.h"
#include "eon/NEBProjection.h"
#include "eon/NEBSplineExtrema.h"
#include "eon/NEBSpringForce.h"
#include "eon/NEBTangent.h"
#include "eon/NEBZoom.h"
#include "eon/Optimizer.h"
#include "eon/PotCapabilities.h"
#include "eon/SafeMath.h"
#include "eon/SolidStateNEB.h"
#ifdef WITH_RGSADDLE
#include "eon/XtsciBand.h"
#endif
#include "magic_enum/magic_enum.hpp"

#include "ForEachImage.h"
#include "eon/EonLogger.h"
#include <cmath>
#include <exception>
#include <format>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#ifdef EON_PARALLEL_NEB
#include <algorithm>
#include <execution>
#include <numeric>
#include <vector>
#endif

namespace eonc {

namespace fs = std::filesystem;

// Nudged Elastic Band definitions
// First constructor: Now delegates to the second constructor
NudgedElasticBand::NudgedElasticBand(std::shared_ptr<Matter> initialPassed,
                                     std::shared_ptr<Matter> finalPassed,
                                     const Parameters &parametersPassed,
                                     std::shared_ptr<Potential> potPassed)
    : NudgedElasticBand(
          [&]() {
            auto &init_opt = parametersPassed.neb_options().initialization;
            const size_t base_count =
                parametersPassed.neb_options().image_count;
            if (parametersPassed.neb_options().solid_state.enabled &&
                init_opt.method != NEBInit::LINEAR &&
                init_opt.method != NEBInit::FILE) {
              throw std::invalid_argument(
                  "solid_state accepts initializer linear or file");
            }
            if (parametersPassed.neb_options().match_endpoints) {
              auto aligned = eonc::IRACompare::alignReactantToProduct(
                  *initialPassed, *finalPassed, 1.0);
              auto *log = eonc::log::get();
              if (aligned.error != 0) {
                QUILL_LOG_WARNING(
                    log,
                    "match_endpoints: IRA align failed (error {}), "
                    "interpolating the input order",
                    aligned.error);
              } else {
                QUILL_LOG_INFO(log,
                               "match_endpoints: Hausdorff {:.4f} A after IRA "
                               "permute+rotate of the reactant",
                               aligned.hausdorffDistance);
              }
            }

            // Apply oversampling factor if flag exists
            const size_t relax_count =
                init_opt.oversampling
                    ? base_count * init_opt.oversampling_factor
                    : base_count;

            std::vector<Matter> path;
            switch (init_opt.method) {
            case NEBInit::FILE: {
              std::vector<fs::path> file_paths =
                  eonc::helpers::neb_paths::readFilePaths(init_opt.input_path);
              path = eonc::helpers::neb_paths::filePathInit(
                  file_paths, *initialPassed, base_count);
              // filePathInit reloads every frame, including ends the caller
              // may already have minimized.
              path.front() = Matter(*initialPassed);
              path.back() = Matter(*finalPassed);
              break;
            }
            case NEBInit::IDPP:
              path = eonc::helpers::neb_paths::idppPath(
                  *initialPassed, *finalPassed, relax_count, parametersPassed);
              break;
            case NEBInit::IDPP_COLLECTIVE:
              path = eonc::helpers::neb_paths::idppCollectivePath(
                  *initialPassed, *finalPassed, relax_count, parametersPassed);
              break;
            case NEBInit::SIDPP:
            case NEBInit::SIDPP_ZBL:
              path = eonc::helpers::neb_paths::sidppPath(
                  *initialPassed, *finalPassed, relax_count, parametersPassed,
                  (init_opt.method == NEBInit::SIDPP_ZBL));
              break;
            case NEBInit::LINEAR:
            default:
              path = eonc::helpers::neb_paths::linearPath(
                  *initialPassed, *finalPassed, base_count);
              break;
            }

            // Decimate path back to the target count using the cubic spline
            if (init_opt.oversampling && path.size() > (base_count + 2)) {
              auto *log = eonc::log::get();
              QUILL_LOG_INFO(log,
                             "Decimating oversampled path ({} images) to "
                             "{} images via cubic spline.",
                             path.size() - 2, base_count);

              // Perform the Spline Resampling
              path = eonc::helpers::neb_paths::resamplePath(path, base_count);

              // POST-DECIMATION RE-RELAXATION
              // The spline might have placed atoms in high-energy positions.
              QUILL_LOG_INFO(
                  log, "Relaxing decimated path to restore IDPP surface...");

              // Collective IDPP objective for the new reduced path
              std::shared_ptr<ObjectiveFunction> post_decim_objf =
                  std::make_shared<CollectiveIDPPObjectiveFunction>(
                      path, parametersPassed);

              // Wrap with ZBL if the original method used ZBL options
              bool use_zbl = (init_opt.method == NEBInit::SIDPP_ZBL);
              if (use_zbl) {
                auto zbl_pot = eonc::helpers::neb_paths::createZBLPotential();
                post_decim_objf = std::make_shared<ZBLRepulsiveIDPPObjective>(
                    post_decim_objf, zbl_pot, path, parametersPassed, 1.0);
              }

              auto optim = eonc::helpers::create::mkOptim(
                  post_decim_objf, parametersPassed.neb_options().opt_method,
                  parametersPassed);

              // Run a short optimization
              optim->run(
                  parametersPassed.neb_options().initialization.max_iterations,
                  parametersPassed.neb_options().initialization.max_move);
            }
            if (parametersPassed.neb_options().solid_state.enabled &&
                init_opt.method == NEBInit::LINEAR) {
              if (init_opt.oversampling) {
                throw std::invalid_argument(
                    "solid_state NEB does not oversample the initial path");
              }
              for (Matter &image : path) {
                eonc::neb::orientSolidStateMatter(image);
              }
              eonc::neb::interpolateSolidStateLinear(path);
            } else if (parametersPassed.neb_options().solid_state.enabled &&
                       init_opt.method == NEBInit::FILE) {
              for (Matter &image : path) {
                eonc::neb::orientSolidStateMatter(image);
              }
            }
            return path;
          }(),
          parametersPassed, potPassed) {}

// Second constructor: Contains all the common setup code
NudgedElasticBand::NudgedElasticBand(std::vector<Matter> initPath,
                                     const Parameters &parametersPassed,
                                     std::shared_ptr<Potential> potPassed)
    : ci_enabled_{parametersPassed.neb_options().climbing_image.enabled},
      params{parametersPassed},
      pot{potPassed},
      E_ref{0.0} {

  log = eonc::log::get();
  eonc::helpers::requireKnownConvergenceMetric(
      params.optimizer_options().convergence_metric, "[Nudged Elastic Band]");
  this->status = NEBStatus::INIT;
  numImages = params.neb_options().image_count;
  if (initPath.empty()) {
    throw std::invalid_argument("NEB: initPath is empty");
  }
  if (initPath.size() != static_cast<size_t>(numImages + 2)) {
    throw std::invalid_argument("NEB: initPath.size() must be image_count + 2");
  }
  atoms = initPath.front().numberOfAtoms();
  for (size_t i = 1; i < initPath.size(); ++i) {
    eonc::helpers::neb_paths::requireSameAtomCount(initPath.front(),
                                                   initPath[i], "path images");
  }

  // Common initialization logic
  path.resize(numImages + 2);
  tangent.resize(numImages + 2);
  projectedForce.resize(numImages + 2);
  extremumPosition.resize(2 * (numImages + 1));
  extremumEnergy.resize(2 * (numImages + 1));
  extremumCurvature.resize(2 * (numImages + 1));
  numExtrema = 0;

  // Create per-image potentials if needed for true parallel force evaluation
  perImagePotentials_ =
      pot->needsPerImageInstance() && params.main_options().parallel;
  if (perImagePotentials_) {
    QUILL_LOG_INFO(log,
                   "NEB: Creating per-image potential instances for "
                   "parallel force evaluation ({} images)",
                   numImages + 2);
  }

  for (long i = 0; i <= numImages + 1; i++) {
    path[i] = std::make_shared<Matter>(std::move(initPath[i]));

    // Give each INTERMEDIATE image its own potential for true parallelism.
    // Endpoints (i=0, i=numImages+1) keep the shared pot -- they are only
    // evaluated once during initialization and never in parallel.
    if (perImagePotentials_ && i > 0 && i <= numImages) {
      auto cloned = pot->clonePotential();
      path[i]->setPotential(cloned ? cloned
                                   : eonc::helpers::sharePotential(
                                         eonc::helpers::makePotential(params)));
    }

    tangent[i] = std::make_shared<AtomMatrix>();
    tangent[i]->resize(atoms, 3);
    tangent[i]->setZero();
    projectedForce[i] = std::make_shared<AtomMatrix>();
    projectedForce[i]->resize(atoms, 3);
    projectedForce[i]->setZero();
  }

  // Common final setup
  movedAfterForceCall = true;
  prepareSolidState();
  // Name images before the endpoint evaluation. The host matches a
  // system to a band index by the position pointer.
  if (pot) {
    pot->bindBand(path);
  }
  // Both endpoints in one call, so two calculator groups take one each.
  {
    Matter *const ends[] = {path[0].get(), path[numImages + 1].get()};
    eonc::evaluateTogether(*pot, ends);
  }
  reactantEnergy = path[0]->getPotentialEnergy();
  climbingImage = 0;

  // Setup springs
  k_u = params.neb_options().spring.weighting.k_max;
  k_l = params.neb_options().spring.weighting.k_min;
  if (params.neb_options().spring.weighting.enabled) {
    ksp = k_l;
  } else {
    ksp = params.neb_options().spring.constant;
  }

  // Cache strategies that are constant across iterations
  tangentStrat_ = eonc::neb::buildTangentStrategy(params);
  projectionStrat_ = eonc::neb::buildProjectionStrategy(params);

  // Optional debugging setup
  if (params.debug_options().estimate_neb_eigenvalues) {
    eigenmode_solvers.resize(numImages + 2);
    for (long i = 0; i <= numImages + 1; i++) {
      eigenmode_solvers[i] =
          eonc::buildEigenmodeStrategy(path[i], parametersPassed, pot);
    }
  }
}

NudgedElasticBand::NEBStatus NudgedElasticBand::compute() {
  iteration_ = 0;
  long &iteration = iteration_;
  this->status = NEBStatus::RUNNING;

  QUILL_LOG_DEBUG(log, "Nudged elastic band calculation started.");

  // Higher endpoint: springs stay at k_min until an image exceeds it.
  E_ref = std::max(path[0]->getPotentialEnergy(),
                   path[numImages + 1]->getPotentialEnergy());

  updateForces();

  auto objf = std::make_shared<NEBObjectiveFunction>(this, params);

  bool switched{false};
#ifdef WITH_RGSADDLE
  std::unique_ptr<XtsciBand> rustBand;
#endif
  std::unique_ptr<Optimizer> optim;
  std::unique_ptr<Optimizer> refine_optim;
#ifdef WITH_RGSADDLE
  if (params.neb_options().opt_method == OptType::XTSCI) {
    rustBand = std::make_unique<XtsciBand>(*this, params);
  } else
#endif
  {
    optim = eonc::helpers::create::mkOptim(
        objf, params.neb_options().opt_method, params);
    if (params.optimizer_options().refine.method != OptType::None) {
      refine_optim = eonc::helpers::create::mkOptim(
          objf, params.optimizer_options().refine.method, params);
    }
  }
  // OCINEB controller
  auto ocinebCfg = eonc::neb::OCINEBController::fromParams(params);
  eonc::neb::OCINEBController ocineb(ocinebCfg);
  bool zoomDone{false};
  long zoomStable{0};
  long zoomPrevCI{-1};
  long zoomAt{-1};

  while (this->status != NEBStatus::GOOD) {
    if (!path.empty()) {
      path.front()->pollCancel("neb");
    }
    if (params.debug_options().write_movies &&
        (iteration % params.debug_options().write_movies_interval == 0)) {
      bool append = (iteration != 0);
      if (!eonc::io::io_ok(eonc::neb::writePathCon(
              path, tangent, eigenmode_solvers, numImages,
              params.debug_options().estimate_neb_eigenvalues,
              std::format("neb_path_{:03d}.con", iteration), iteration,
              reactantEnergy))) {
        QUILL_LOG_ERROR(log, "Failed to write NEB path movie for iteration {}",
                        iteration);
      }

      AtomMatrix maxTang;
      if (maxEnergyImage == 0) {
        maxTang =
            path[0]->pbc(path[1]->getPositions() - path[0]->getPositions());
      } else if (maxEnergyImage == static_cast<size_t>(numImages + 1)) {
        maxTang = path[numImages]->pbc(path[numImages + 1]->getPositions() -
                                       path[numImages]->getPositions());
      } else {
        maxTang = *tangent[maxEnergyImage];
      }
      eonc::safemath::safe_normalize_inplace(maxTang);
      auto maxImageMetadata = eonc::io::ConFrameMetadata{};
      maxImageMetadata.frame_index = static_cast<uint64_t>(maxEnergyImage);
      maxImageMetadata.energy = path[maxEnergyImage]->getPotentialEnergy();
      maxImageMetadata.neb_bead = static_cast<uint64_t>(maxEnergyImage);
      maxImageMetadata.neb_band = static_cast<uint64_t>(iteration);
      maxImageMetadata.scalars.push_back(
          {"relative_energy",
           path[maxEnergyImage]->getPotentialEnergy() - reactantEnergy});
      maxImageMetadata.scalars.push_back(
          {"parallel_force",
           matDot(path[maxEnergyImage]->getForces(), maxTang)});
      maxImageMetadata.strings.push_back({"movie_kind", "neb_maximage"});
      if (!eonc::io::io_ok(path[maxEnergyImage]->matter2con(
              "neb_maximage.con", append, &maxImageMetadata))) {
        EONC_LOG_WARNING("Failed to write neb_maximage.con");
      }
      printImageData(true, iteration);
    }

    VectorXd pos = objf->getPositions();
    double convForce = convergenceForce();

    ocineb.updateStability(climbingImage);

    if (iteration == 0) {
      baseline_force = convForce;
      ocineb.initBaseline(convForce);

      // Log configuration banner
      auto &ci_opt = params.neb_options().climbing_image;
      auto &mmf_opt = ci_opt.ocineb;
      auto fmt_trigger = [](double val) -> std::string {
        if (val > 1e100)
          return "INF";
        return std::format("{:.4f}", val);
      };

      QUILL_LOG_INFO(
          log,
          "===============================================================");
      QUILL_LOG_INFO(log, " NEB Optimization Configuration");
      QUILL_LOG_INFO(
          log,
          "===============================================================");
      QUILL_LOG_INFO(log, " {:<25} : {:.4f}", "Baseline Force", baseline_force);

      std::string ci_status = ci_opt.enabled ? "ENABLED" : "DISABLED";
      QUILL_LOG_INFO(log, " {:<25} : {}", "Climbing Image (CI)", ci_status);
      if (ci_opt.enabled) {
        double ci_rel_val = baseline_force * ci_opt.trigger_factor;
        QUILL_LOG_INFO(log, "   - {:<21} : {} (Factor: {:.2f})",
                       "Relative Trigger", fmt_trigger(ci_rel_val),
                       ci_opt.trigger_factor);
        QUILL_LOG_INFO(log, "   - {:<21} : {}", "Absolute Trigger",
                       fmt_trigger(ci_opt.trigger_force));
        QUILL_LOG_INFO(log, "   - {:<21} : {}", "Converged Only",
                       ci_opt.converged_only);
      }

      std::string mmf_status =
          (ci_opt.enabled && mmf_opt.use_mmf) ? "ENABLED" : "DISABLED";
      QUILL_LOG_INFO(log, " {:<25} : {}", "Hybrid MMF (OCINEB)", mmf_status);
      if (ci_opt.enabled && mmf_opt.use_mmf) {
        QUILL_LOG_INFO(log, "   - {:<21} : {:.4f} (Factor: {:.2f})",
                       "Initial Threshold", ocineb.threshold(),
                       mmf_opt.trigger_factor);
        QUILL_LOG_INFO(log, "   - {:<21} : {:.4f}", "Absolute Floor",
                       mmf_opt.trigger_force);
        QUILL_LOG_INFO(log, "   - {:<21} : {:.4f}", "Angle Tolerance",
                       mmf_opt.angle_tol);
      }
      QUILL_LOG_INFO(
          log,
          "---------------------------------------------------------------");

      EONC_LOG_DEBUG("{:>10s} {:>12s} {:>14s} {:>11s} {:>12s}", "iteration",
                     "step size",
                     params.optimizer_options().convergence_metric_label,
                     "max image", "max energy");
      QUILL_LOG_DEBUG(
          eonc::log::get(),
          "---------------------------------------------------------------\n");
    }

    // CI active when force drops below relative threshold
    bool ci_active =
        params.neb_options().climbing_image.enabled &&
        (convForce < baseline_force *
                         params.neb_options().climbing_image.trigger_factor ||
         convForce < params.neb_options().climbing_image.trigger_force);

    bool zoomedThisStep = false;
    if (iteration && !zoomDone && params.neb_options().zoom.enabled &&
        climbingImage > 0 && climbingImage + 1 < path.size()) {
      if (static_cast<long>(climbingImage) == zoomPrevCI) {
        ++zoomStable;
      } else {
        zoomPrevCI = static_cast<long>(climbingImage);
        zoomStable = 1;
      }
      const auto &zoom = params.neb_options().zoom;
      const double zoomForce =
          zoom.activation_threshold > 0.0
              ? zoom.activation_threshold
              : 10.0 * params.neb_options().force_tolerance;
      if (zoomStable >= zoom.stability_count && convForce < zoomForce) {
        std::vector<double> energy;
        energy.reserve(path.size());
        for (const auto &image : path) {
          energy.push_back(image->getPotentialEnergy());
        }
        const auto window =
            eonc::neb::zoom::selectWindow(energy, climbingImage, zoom);
        if (eonc::neb::zoom::redistributePath(path, window,
                                              zoom.interpolation)) {
          movedAfterForceCall = true;
          optim = eonc::helpers::create::mkOptim(
              objf, params.neb_options().opt_method, params);
          zoomDone = true;
          zoomedThisStep = true;
          zoomAt = iteration;
          QUILL_LOG_INFO(log, "Zoom-NEB: packed the band onto images [{}, {}]",
                         window.lo, window.hi);
        }
      }
    }

    if (iteration) {
      // MMF triggering via controller. Skipped on the zoom step so the
      // dimer sees the redistributed band, not the pre-zoom geometry.
      if (!zoomedThisStep &&
          ocineb.shouldTrigger(convForce, ci_active, climbingImage, numImages,
                               ocineb.stabilityCount())) {
        auto result = ocineb.run(*this, convForce);

        if (result.convergedAfterMMF) {
          status = NEBStatus::GOOD;
          break;
        }
        // Post-MMF arc-length reparameterization: pass full path
        // (endpoints are fixed by resamplePathInPlace, only interior
        // images are redistributed). Zero force-call cost; next NEB
        // iteration recomputes all forces anyway.
        bool didResample = false;
        if (!result.convergedAfterMMF && result.newForce < convForce &&
            path[climbingImage]->getPeriodic() &&
            path[0]->numberOfAtoms() > 6) {
          eonc::helpers::neb_paths::resamplePathInPlace(
              std::span{path.data(), path.size()});
          movedAfterForceCall = true;
          didResample = true;
        }

        // Reset optimizer AFTER reparameterization so fresh L-BFGS
        // starts from the redistributed positions.
        if (result.shouldResetOptimizer || didResample) {
#ifdef WITH_RGSADDLE
          if (rustBand) {
            rustBand->syncFromPath();
            rustBand->reset();
          } else
#endif
          {
            optim = eonc::helpers::create::mkOptim(
                objf, params.neb_options().opt_method, params);
          }
        }
      }

      long iterLimit = params.neb_options().max_iterations;
      if (zoomDone && params.neb_options().zoom.max_iterations > 0 &&
          zoomAt >= 0) {
        iterLimit = zoomAt + params.neb_options().zoom.max_iterations;
      }
      if (iteration >= iterLimit) {
        status = NEBStatus::BAD_MAX_ITERATIONS;
        break;
      }

      if (zoomedThisStep) {
        iteration++;
        continue;
      }

      // Set CI state so updateForces() inside the optimizer step
      // applies the correct force projection.
      setCIEnabled(ci_active);

#ifdef WITH_RGSADDLE
      if (rustBand) {
        rustBand->step(params.optimizer_options().max_move);
      } else
#endif
      {
        auto &activeOptim =
            (refine_optim &&
             convForce <= params.optimizer_options().refine.threshold)
                ? refine_optim
                : optim;
        if (refine_optim &&
            convForce <= params.optimizer_options().refine.threshold &&
            !switched) {
          switched = true;
          EONC_LOG_DEBUG("Switched to {}",
                         magic_enum::enum_name<OptType>(
                             params.optimizer_options().refine.method));
        }
        activeOptim->step(params.optimizer_options().max_move);
      }

      setCIEnabled(params.neb_options().climbing_image.enabled);
    }

    iteration++;

    double dE = path[maxEnergyImage]->getPotentialEnergy() - reactantEnergy;
    double stepSize = 0.0;
    if (solidState_) {
      const VectorXd delta = objf->difference(objf->getPositions(), pos);
      const long seg = 3L * atoms + 9L;
      for (long image = 0; image < numImages; ++image) {
        stepSize = std::max(
            stepSize, delta.segment(image * seg, seg).cwiseAbs().maxCoeff());
      }
    } else {
      stepSize = eonc::geometry::maxAtomMotionV(
          path[0]->pbcV(objf->getPositions() - pos));
    }
    QUILL_LOG_DEBUG(log, "{:>10} {:>12.4e} {:>14.4e} {:>11} {:>12.4}",
                    iteration, stepSize, convergenceForce(), maxEnergyImage,
                    dE);

    if (pot->getType() == PotType::CatLearn) {
      if (objf->isUncertain()) {
        QUILL_LOG_DEBUG(log, "NEB failed due to high uncertainty");
        status = NEBStatus::MAX_UNCERTAINTY;
        break;
      } else if (objf->isConverged()) {
        QUILL_LOG_DEBUG(log, "NEB converged\n");
        status = NEBStatus::GOOD;
        break;
      }
    } else {
      if (objf->isConverged()) {
        QUILL_LOG_DEBUG(log, "NEB converged\n");
        status = NEBStatus::GOOD;
        break;
      }
    }
  }
  return status;
}

// generate the force value that is compared to the convergence criterion
double NudgedElasticBand::convergenceForce() {
  if (movedAfterForceCall)
    updateForces();

  auto imageForce = [&](long i) -> double {
    const double cellNorm =
        solidState_ ? projectedCellForce[static_cast<size_t>(i)].norm() : 0.0;
    if (params.optimizer_options().convergence_metric == "norm") {
      return std::hypot(projectedForce[i]->norm(), cellNorm);
    }
    if (params.optimizer_options().convergence_metric == "max_atom") {
      // Every image shares the reactant constraint mask.
      return std::max(path[0]->maxFreeAtomForce(*projectedForce[i]), cellNorm);
    }
    if (params.optimizer_options().convergence_metric == "max_component") {
      double component = projectedForce[i]->cwiseAbs().maxCoeff();
      if (solidState_) {
        component = std::max(
            component,
            projectedCellForce[static_cast<size_t>(i)].cwiseAbs().maxCoeff());
      }
      return component;
    }
    log = eonc::log::traceback();
    QUILL_LOG_CRITICAL(
        log, "[Nudged Elastic Band] unknown opt_convergence_metric: {}",
        params.optimizer_options().convergence_metric);
    throw std::invalid_argument(
        std::format("[Nudged Elastic Band] unknown convergence_metric: {}",
                    params.optimizer_options().convergence_metric));
  };

  double bandMax = 0;
  for (long i = 1; i <= numImages; ++i) {
    bandMax = std::max(bandMax, imageForce(i));
  }

  const bool ciOnly = params.neb_options().climbing_image.converged_only &&
                      ci_enabled_ && climbingImage != 0;
  if (!ciOnly) {
    return bandMax;
  }

  const double ciForce = imageForce(climbingImage);
  const double slack = params.neb_options().climbing_image.band_slack;
  const double tol = params.neb_options().force_tolerance;
  if (slack > 0.0 && bandMax > slack * tol) {
    return bandMax;
  }
  return ciForce;
}

// Update the forces, do the projections, and add spring forces
void NudgedElasticBand::updateForces(bool ci_active) {
  // Update forces for all intermediate images. Prefer batched evaluation
  // (single model.forward() over all dirty images, e.g. MetatomicPotential
  // on GPU). Else fall back to per-image evaluation, which is itself
  // thread-parallel when (a) the potential is thread-safe on the same
  // instance, or (b) per-image instances were created (separate models).
  if (pot->supportsBatchEvaluation() && numImages > 1) {
    // Collect only images that need recomputation (positions changed).
    // Materialize atomic numbers and cells first, then build the raw-pointer
    // arrays after storage is stable. Otherwise vector growth can invalidate
    // earlier .data() pointers and hand garbage cells/types to forceBatch().
    std::vector<long> dirty; // indices into path[] (1-based)
    dirty.reserve(numImages);
    for (long i = 1; i <= numImages; i++) {
      if (path[i]->needsForceUpdate()) {
        dirty.push_back(i);
      }
    }

    if (!dirty.empty()) {
      if (!path.empty()) {
        path.front()->pollCancel("neb");
      }
      auto nDirty = static_cast<long>(dirty.size());
      std::vector<VectorXi> nrsStore;
      std::vector<Matrix3d> boxStore;
      std::vector<const double *> posVec, boxVec;
      std::vector<const int *> nrsVec;
      std::vector<double *> frcVec;
      nrsStore.reserve(static_cast<size_t>(nDirty));
      boxStore.reserve(static_cast<size_t>(nDirty));
      posVec.reserve(static_cast<size_t>(nDirty));
      boxVec.reserve(static_cast<size_t>(nDirty));
      nrsVec.reserve(static_cast<size_t>(nDirty));
      frcVec.reserve(static_cast<size_t>(nDirty));

      for (long idx : dirty) {
        nrsStore.push_back(path[idx]->getAtomicNrs());
        // Isolated molecules still store a box for I/O. Pots that infer
        // PBC from a non-zero cell must see a zero box, as
        // Matter::computePotential does on the endpoints.
        boxStore.push_back(
            (path[idx]->getPeriodic() || pot->forwardsStoredCell())
                ? path[idx]->getCell()
                : Matrix3d::Zero());
      }
      for (long j = 0; j < nDirty; j++) {
        auto idx = dirty[static_cast<size_t>(j)];
        posVec.push_back(path[idx]->getPositions().data());
        nrsVec.push_back(nrsStore[static_cast<size_t>(j)].data());
        frcVec.push_back(path[idx]->forcesData());
        boxVec.push_back(boxStore[static_cast<size_t>(j)].data());
      }

      std::vector<double> energies(nDirty), variances(nDirty);
      // Image i is system i - 1 to the potential's router, whether or not
      // the images before it are dirty.
      std::vector<long> owners(static_cast<size_t>(nDirty));
      for (long j = 0; j < nDirty; j++) {
        owners[static_cast<size_t>(j)] = dirty[static_cast<size_t>(j)] - 1;
      }
      pot->forceBatchOwned(nDirty, atoms, posVec.data(), nrsVec.data(),
                           frcVec.data(), energies.data(), variances.data(),
                           boxVec.data(), owners.data());
      for (long j = 0; j < nDirty; j++) {
        path[dirty[j]]->setComputedPotential(energies[j], variances[j]);
      }
    }
  } else {
    // Per-image evaluation (sequential or parallel threads)
    bool canParallel =
        eonc::potAllowsSharedInstance(*pot) || perImagePotentials_;
    if (numImages > 1 && params.main_options().parallel && canParallel) {
#ifdef EON_PARALLEL_NEB
      // nvc++ -stdpar=multicore|gpu (meson -Dstdpar=cpu|gpu).
      std::vector<long> beads(static_cast<size_t>(numImages));
      std::iota(beads.begin(), beads.end(), 1);
      std::exception_ptr neb_cancel;
      std::mutex neb_cancel_mu;
      std::for_each(std::execution::par, beads.begin(), beads.end(),
                    [this, &neb_cancel, &neb_cancel_mu](long i) {
                      try {
                        path[i]->getForcesRaw();
                      } catch (const JobCancelled &) {
                        std::lock_guard<std::mutex> lock(neb_cancel_mu);
                        if (!neb_cancel) {
                          neb_cancel = std::current_exception();
                        }
                      }
                    });
      if (neb_cancel) {
        std::rethrow_exception(neb_cancel);
      }
#else
      eonc::forEachImage(numImages,
                         [this](long i) { path[i]->getForcesRaw(); });
#endif
    } else {
      for (long i = 1; i <= numImages; i++) {
        path[i]->getForcesRaw();
      }
    }
  }

  if (solidState_) {
    projectSolidState(ci_active);
    movedAfterForceCall = false;
    return;
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
  double maxEnergy = (*it)->getPotentialEnergy();

  // Update E_ref for energy weighting. The higher endpoint keeps the
  // soft spring on the side below that minimum.
  if (params.neb_options().spring.weighting.enabled) {
    E_ref = std::max(path[0]->getPotentialEnergy(),
                     path[numImages + 1]->getPotentialEnergy());
  }

  // Climbing requires an interior peak above both fixed endpoints.
  // A monotonic path retains the spring force on every interior image.
  const double endpointEnergy = std::max(path.front()->getPotentialEnergy(),
                                         path.back()->getPotentialEnergy());
  const bool climb = ci_active && maxEnergy > endpointEnergy;
  const long previousCI = static_cast<long>(climbingImage);
  climbingImage = 0;

  // Spring strategy must be rebuilt each iteration (depends on maxEnergy,
  // E_ref). Tangent and projection strategies are cached as members.
  auto spring = eonc::neb::buildSpringStrategy(params, path, numImages, atoms,
                                               maxEnergy, E_ref);

  if (const auto *uniform = std::get_if<eonc::neb::UniformSpring>(&spring)) {
    ksp = uniform->ksp;
  }
  // The climbing image, chosen before the per-image work so that work
  // writes nothing outside its own image.
  // Re-picking the highest image every call hops the climber onto a
  // shoulder, so the current one stays while it is a strict interior
  // local maximum.
  long ciTarget = static_cast<long>(maxEnergyImage);
  if (climb && previousCI > 0 && previousCI <= numImages) {
    const double eCI = path[previousCI]->getPotentialEnergy();
    if (eCI > path[previousCI - 1]->getPotentialEnergy() &&
        eCI > path[previousCI + 1]->getPotentialEnergy()) {
      ciTarget = previousCI;
    }
  }
  if (climb) {
    climbingImage = static_cast<size_t>(ciTarget);
  }

  // Each image reads its neighbours' positions and energies and writes
  // only its own tangent and projected force, so with [Main] parallel the
  // images are projected on the image pool like their forces.
  auto projectImage = [&](long i) {
    AtomMatrix posDiffNext(atoms, 3), posDiffPrev(atoms, 3);
    const AtomMatrix &force = path[i]->getForces();
    const AtomMatrix &pos = path[i]->getPositions();
    const AtomMatrix &posPrev = path[i - 1]->getPositions();
    const AtomMatrix &posNext = path[i + 1]->getPositions();
    double energy = path[i]->getPotentialEnergy();
    double energyPrev = path[i - 1]->getPotentialEnergy();
    double energyNext = path[i + 1]->getPotentialEnergy();
    posDiffNext.noalias() = posNext - pos;
    posDiffNext = path[i]->pbc(posDiffNext);
    posDiffPrev.noalias() = pos - posPrev;
    posDiffPrev = path[i]->pbc(posDiffPrev);
    double distNext = posDiffNext.norm();
    double distPrev = posDiffPrev.norm();

    // Tangent via strategy dispatch
    *tangent[i] = std::visit(
        [&](auto &t) {
          return t.compute(posDiffNext, posDiffPrev, energy, energyPrev,
                           energyNext);
        },
        tangentStrat_);

    // Spring forces via strategy dispatch
    eonc::neb::SpringResult springResult = std::visit(
        [&](auto &s) -> eonc::neb::SpringResult {
          using T = std::decay_t<decltype(s)>;
          if constexpr (std::is_same_v<T, eonc::neb::UniformSpring>) {
            return s.compute(i, *tangent[i], distNext, distPrev, posDiffNext,
                             posDiffPrev, path[i]);
          } else if constexpr (std::is_same_v<T, eonc::neb::WeightedSpring>) {
            return s.compute(i, *tangent[i], distNext, distPrev);
          } else {
            return s.compute(i, *tangent[i], posNext, posPrev, pos, path[i]);
          }
        },
        spring);

    // Climbing image or projected force
    if (climb && i == ciTarget) {
      // CI force: F - 2*(F.t)*t, plus DNEB correction if active
      AtomMatrix forceDNEB = AtomMatrix::Zero(atoms, 3);
      if (const auto *dnebProj =
              std::get_if<eonc::neb::DNEB_Projection>(&projectionStrat_)) {
        AtomMatrix fPerp = eonc::neb::forcePerp(force, *tangent[i]);
        forceDNEB = eonc::neb::computeDNEBComponent(springResult.forceSpring,
                                                    *tangent[i], fPerp,
                                                    dnebProj->use_switching);
      }
      *projectedForce[i] =
          eonc::neb::climbingImageForce(force, *tangent[i], forceDNEB);
    } else {
      eonc::neb::ImageForceData data{force, *tangent[i], springResult,
                                     path[i]->numberOfFreeAtoms(),
                                     path[i]->numberOfAtoms()};
      *projectedForce[i] = std::visit([&](auto &p) { return p.project(data); },
                                      projectionStrat_);
    }

    // The spring force acts along the tangent, which has a component on a
    // fixed atom whenever the endpoints place it differently; a fixed atom
    // carries no force on the band.
    if (path[i]->numberOfFreeAtoms() < path[i]->numberOfAtoms()) {
      for (long j = 0; j < atoms; j++) {
        if (path[i]->getFixed(j)) {
          projectedForce[i]->row(j).setZero();
        }
      }
    }

    eonc::neb::zeroTranslation(*projectedForce[i], path[i]->numberOfFreeAtoms(),
                               path[i]->numberOfAtoms());
  };
  if (numImages > 1 && params.main_options().parallel) {
    eonc::forEachImage(numImages, projectImage);
  } else {
    for (long i = 1; i <= numImages; i++) {
      projectImage(i);
    }
  }

  movedAfterForceCall = false;
}

// Thin wrappers delegating to eonc::neb:: free functions

void NudgedElasticBand::printImageData(bool writeToFile, size_t idx) {
  eonc::neb::printImageData(path, tangent, eigenmode_solvers, numImages,
                            params.debug_options().estimate_neb_eigenvalues,
                            writeToFile, idx, log, reactantEnergy);
}

void NudgedElasticBand::findExtrema() {
  auto result = eonc::neb::findSplineExtrema(path, tangent, numImages);
  numExtrema = result.numExtrema;
  extremumPosition = std::move(result.positions);
  extremumEnergy = std::move(result.energies);
  extremumCurvature = std::move(result.curvatures);
}

std::vector<readcon::ConFrame>
NudgedElasticBand::pathFrames(std::optional<size_t> bandIndex) {
  return eonc::neb::pathToConFrames(
      path, tangent, eigenmode_solvers, numImages,
      params.debug_options().estimate_neb_eigenvalues, bandIndex,
      reactantEnergy);
}

namespace {

void refuseSolidStateCombination(const Parameters &params) {
  const auto &neb = params.neb_options();
  if (!neb.solid_state.enabled) {
    return;
  }
  if (neb.climbing_image.ocineb.use_mmf) {
    throw std::invalid_argument(
        "solid_state is set and ci_mmf is true. The min-mode walk does not "
        "move the cell");
  }
  if (neb.zoom.enabled) {
    throw std::invalid_argument(
        "solid_state is set and zoom_neb is true. Zoom does not move the cell");
  }
  if (neb.spring.om.enabled) {
    throw std::invalid_argument(
        "solid_state is set and onsager_machlup is true");
  }
  if (neb.spring.doubly_nudged) {
    throw std::invalid_argument(
        "solid_state is set and neb_doubly_nudged is true");
  }
  if (neb.spring.use_elastic_band) {
    throw std::invalid_argument(
        "solid_state is set and neb_elastic_band is true");
  }
  const auto method = neb.initialization.method;
  if (method != NEBInit::LINEAR && method != NEBInit::FILE) {
    throw std::invalid_argument(
        "solid_state accepts initializer linear or file");
  }
  if (!(neb.solid_state.weight > 0.0)) {
    throw std::invalid_argument("solid_state_weight must be positive");
  }
}

AtomMatrix packJoint(const AtomMatrix &atomic, const Matrix3d &cell) {
  AtomMatrix packed(atomic.rows() + 3, 3);
  packed.topRows(atomic.rows()) = atomic;
  packed.bottomRows(3) = cell;
  return packed;
}

double springScale(const eonc::neb::SpringStrategy &spring, long image,
                   double distNext, double distPrev) {
  if (const auto *uniform = std::get_if<eonc::neb::UniformSpring>(&spring)) {
    return uniform->ksp * (distNext - distPrev);
  }
  if (const auto *weighted = std::get_if<eonc::neb::WeightedSpring>(&spring)) {
    return weighted->springConstants[static_cast<size_t>(image)] * distNext -
           weighted->springConstants[static_cast<size_t>(image - 1)] * distPrev;
  }
  throw std::invalid_argument(
      "solid_state NEB does not use this spring strategy");
}

} // namespace

void NudgedElasticBand::prepareSolidState() {
  if (!params.neb_options().solid_state.enabled) {
    return;
  }
  refuseSolidStateCombination(params);
  for (const auto &image : path) {
    if (!image->getPeriodic()) {
      throw std::invalid_argument(
          "solid_state NEB requires periodic boundaries on every image");
    }
    eonc::neb::orientSolidStateMatter(*image);
  }
  const double volume0 = std::abs(path.front()->getCell().determinant());
  const double volume1 = std::abs(path.back()->getCell().determinant());
  solidJacobian_ =
      eonc::neb::solidStateJacobian(0.5 * (volume0 + volume1), atoms,
                                    params.neb_options().solid_state.weight);
  projectedCellForce.assign(static_cast<size_t>(numImages + 2),
                            Matrix3d::Zero());
  solidState_ = true;
  QUILL_LOG_INFO(log, "Solid-state NEB Jacobian {:.6f} Angstrom",
                 solidJacobian_);
}

void NudgedElasticBand::projectSolidState(bool ci_active) {
  const double pressure = params.neb_options().solid_state.pressure;
  const Matrix3d external = Matrix3d::Identity() * pressure;
  std::vector<double> enthalpy(static_cast<size_t>(numImages + 2));
  for (long i = 0; i <= numImages + 1; ++i) {
    enthalpy[static_cast<size_t>(i)] =
        eonc::neb::solidStateEnthalpy(*path[i], *path[0], pressure);
  }

  long highest = 1;
  for (long i = 2; i <= numImages; ++i) {
    if (enthalpy[static_cast<size_t>(i)] >
        enthalpy[static_cast<size_t>(highest)]) {
      highest = i;
    }
  }
  maxEnergyImage = static_cast<size_t>(highest);
  const double endpointEnthalpy = std::max(enthalpy.front(), enthalpy.back());
  const bool climb =
      ci_active && enthalpy[static_cast<size_t>(highest)] > endpointEnthalpy;
  climbingImage = 0;

  if (params.neb_options().spring.weighting.enabled) {
    E_ref = std::max(path[0]->getPotentialEnergy(),
                     path[numImages + 1]->getPotentialEnergy());
  }
  const double maxEnergy = path[highest]->getPotentialEnergy();
  auto spring = eonc::neb::buildSpringStrategy(params, path, numImages, atoms,
                                               maxEnergy, E_ref);
  if (const auto *uniform = std::get_if<eonc::neb::UniformSpring>(&spring)) {
    ksp = uniform->ksp;
  }

  if (!pot->computesStress()) {
    throw std::invalid_argument(
        "solid_state NEB requires a potential that reports the Cauchy stress");
  }

  for (long i = 1; i <= numImages; ++i) {
    const double volume = std::abs(path[i]->getCell().determinant());
    Matrix3d stress = path[i]->cauchyStress();
    Matrix3d cellTrue =
        eonc::neb::cellNebForce(stress, volume, solidJacobian_, external);
    AtomMatrix atomicTrue = path[i]->getForces();
    for (long atom = 0; atom < atoms; ++atom) {
      if (path[i]->getFixed(atom)) {
        atomicTrue.row(atom).setZero();
      }
    }

    const eonc::neb::JointBlock toNext =
        eonc::neb::jointDisplacement(*path[i], *path[i + 1], solidJacobian_);
    const eonc::neb::JointBlock toPrev =
        eonc::neb::jointDisplacement(*path[i - 1], *path[i], solidJacobian_);
    const AtomMatrix packedNext = packJoint(toNext.atomic, toNext.cell);
    const AtomMatrix packedPrev = packJoint(toPrev.atomic, toPrev.cell);
    const double energy = enthalpy[static_cast<size_t>(i)];
    const double energyPrev = enthalpy[static_cast<size_t>(i - 1)];
    const double energyNext = enthalpy[static_cast<size_t>(i + 1)];
    const AtomMatrix packedTangent = std::visit(
        [&](auto &tangentStrategy) {
          return tangentStrategy.compute(packedNext, packedPrev, energy,
                                         energyPrev, energyNext);
        },
        tangentStrat_);
    *tangent[i] = packedTangent.topRows(atoms);

    AtomMatrix packedForce = packJoint(atomicTrue, cellTrue);
    AtomMatrix total;
    if (climb && i == highest) {
      climbingImage = static_cast<size_t>(highest);
      const double parallel = matDot(packedForce, packedTangent);
      total = packedForce - 2.0 * parallel * packedTangent;
    } else {
      const double parallel = matDot(packedForce, packedTangent);
      const AtomMatrix perpendicular = packedForce - parallel * packedTangent;
      const double scale = springScale(spring, i, eonc::neb::jointNorm(toNext),
                                       eonc::neb::jointNorm(toPrev));
      total = perpendicular + scale * packedTangent;
    }
    *projectedForce[i] = total.topRows(atoms);
    projectedCellForce[static_cast<size_t>(i)] = total.bottomRows(3);
    projectedCellForce[static_cast<size_t>(i)](0, 1) = 0.0;
    projectedCellForce[static_cast<size_t>(i)](0, 2) = 0.0;
    projectedCellForce[static_cast<size_t>(i)](1, 2) = 0.0;
    for (long atom = 0; atom < atoms; ++atom) {
      if (path[i]->getFixed(atom)) {
        projectedForce[i]->row(atom).setZero();
      }
    }
    eonc::neb::zeroTranslation(*projectedForce[i], path[i]->numberOfFreeAtoms(),
                               path[i]->numberOfAtoms());
  }
}

} // namespace eonc
