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
#include "EigenmodeStrategy.h"
#include "IDPPObjectiveFunction.hpp"
#include "NEBForceProjection.h"
#include "NEBInitialPaths.hpp"
#include "NEBOcinebController.h"
#include "NEBProjection.h"
#include "NEBSplineExtrema.h"
#include "NEBSpringForce.h"
#include "NEBTangent.h"
#include "Optimizer.h"
#include "magic_enum/magic_enum.hpp"

#include "EonLogger.h"
#include <thread>
using namespace eonc::helpers;
namespace fs = std::filesystem;

// Nudged Elastic Band definitions
// First constructor: Now delegates to the second constructor
NudgedElasticBand::NudgedElasticBand(std::shared_ptr<Matter> initialPassed,
                                     std::shared_ptr<Matter> finalPassed,
                                     const Parameters &parametersPassed,
                                     std::shared_ptr<Potential> potPassed)
    : NudgedElasticBand(
          [&]() {
            auto &init_opt = parametersPassed.neb_options.initialization;
            const size_t base_count = parametersPassed.neb_options.image_count;

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
                  post_decim_objf, parametersPassed.neb_options.opt_method,
                  parametersPassed);

              // Run a short optimization
              optim->run(
                  parametersPassed.neb_options.initialization.max_iterations,
                  parametersPassed.neb_options.initialization.max_move);
            }
            return path;
          }(),
          parametersPassed, potPassed) {}

// Second constructor: Contains all the common setup code
NudgedElasticBand::NudgedElasticBand(std::vector<Matter> initPath,
                                     const Parameters &parametersPassed,
                                     std::shared_ptr<Potential> potPassed)
    : ci_enabled_{parametersPassed.neb_options.climbing_image.enabled},
      params{parametersPassed},
      pot{potPassed},
      E_ref{0.0} {

  log = eonc::log::get();
  this->status = NEBStatus::INIT;
  numImages = params.neb_options.image_count;
  atoms = initPath.front().numberOfAtoms();

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
      pot->needsPerImageInstance() && params.main_options.parallel;
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
      path[i]->setPotential(eonc::helpers::makePotential(params));
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
  path[0]->getPotentialEnergy();
  path[numImages + 1]->getPotentialEnergy();
  climbingImage = 0;

  // Setup springs
  k_u = params.neb_options.spring.weighting.k_max;
  k_l = params.neb_options.spring.weighting.k_min;
  if (params.neb_options.spring.weighting.enabled) {
    ksp = k_l;
  } else {
    ksp = params.neb_options.spring.constant;
  }

  // Cache strategies that are constant across iterations
  tangentStrat_ = eonc::neb::buildTangentStrategy(params);
  projectionStrat_ = eonc::neb::buildProjectionStrategy(params);

  // Optional debugging setup
  if (params.debug_options.estimate_neb_eigenvalues) {
    eigenmode_solvers.resize(numImages + 2);
    for (long i = 0; i <= numImages + 1; i++) {
      eigenmode_solvers[i] =
          eonc::buildEigenmodeStrategy(path[i], parametersPassed, pot);
    }
  }
}

NudgedElasticBand::NEBStatus NudgedElasticBand::compute() {
  long iteration = 0;
  this->status = NEBStatus::RUNNING;

  QUILL_LOG_DEBUG(log, "Nudged elastic band calculation started.");

  // Initialize E_ref for energy weighting
  E_ref = std::min(path[0]->getPotentialEnergy(),
                   path[numImages + 1]->getPotentialEnergy());

  updateForces();

  auto objf = std::make_shared<NEBObjectiveFunction>(this, params);

  bool switched{false};
  auto optim = eonc::helpers::create::mkOptim(
      objf, params.neb_options.opt_method, params);
  std::unique_ptr<Optimizer> refine_optim{nullptr};
  if (params.optimizer_options.refine.method != OptType::None) {
    refine_optim = eonc::helpers::create::mkOptim(
        objf, params.optimizer_options.refine.method, params);
  }

  // OCINEB controller
  auto ocinebCfg = eonc::neb::OCINEBController::fromParams(params);
  eonc::neb::OCINEBController ocineb(ocinebCfg);

  while (this->status != NEBStatus::GOOD) {
    if (params.debug_options.write_movies &&
        (iteration % params.debug_options.write_movies_interval == 0)) {
      bool append = (iteration != 0);
      path[maxEnergyImage]->matter2con("neb_maximage.con", append);
      std::string nebFilename(std::format("neb_path_{:03d}.con", iteration));
      for (long idx = 0; idx <= numImages + 1; idx++) {
        path[idx]->matter2con(nebFilename, /*append=*/idx > 0);
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
      auto &ci_opt = params.neb_options.climbing_image;
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
        QUILL_LOG_INFO(log, "   - {:<21} : {:.2f} (Base: {:.2f}, Str: {:.2f})",
                       "Penalty Scheme", mmf_opt.penalty.base,
                       mmf_opt.penalty.base, mmf_opt.penalty.strength);
        QUILL_LOG_INFO(log, "   - {:<21} : {:.4f}", "Angle Tolerance",
                       mmf_opt.angle_tol);
      }
      QUILL_LOG_INFO(
          log,
          "---------------------------------------------------------------");

      EONC_LOG_DEBUG("{:>10s} {:>12s} {:>14s} {:>11s} {:>12s}", "iteration",
                     "step size",
                     params.optimizer_options.convergence_metric_label,
                     "max image", "max energy");
      QUILL_LOG_DEBUG(
          eonc::log::get(),
          "---------------------------------------------------------------\n");
    }

    // CI active when force drops below relative threshold
    bool ci_active =
        params.neb_options.climbing_image.enabled &&
        (convForce < baseline_force *
                         params.neb_options.climbing_image.trigger_factor ||
         convForce < params.neb_options.climbing_image.trigger_force);

    if (iteration) {
      // MMF triggering via controller
      if (ocineb.shouldTrigger(convForce, ci_active, climbingImage, numImages,
                               ocineb.stabilityCount())) {
        auto result = ocineb.run(*this, convForce);

        if (result.convergedAfterMMF) {
          status = NEBStatus::GOOD;
          break;
        }
        if (result.shouldResetOptimizer) {
          optim = eonc::helpers::create::mkOptim(
              objf, params.neb_options.opt_method, params);
        }
      }

      if (iteration >= params.neb_options.max_iterations) {
        status = NEBStatus::BAD_MAX_ITERATIONS;
        break;
      }

      // Set CI state so updateForces() inside the optimizer step
      // applies the correct force projection.
      setCIEnabled(ci_active);

      auto &activeOptim =
          (refine_optim &&
           convForce <= params.optimizer_options.refine.threshold)
              ? refine_optim
              : optim;
      if (refine_optim &&
          convForce <= params.optimizer_options.refine.threshold && !switched) {
        switched = true;
        EONC_LOG_DEBUG("Switched to {}",
                       magic_enum::enum_name<OptType>(
                           params.optimizer_options.refine.method));
      }
      activeOptim->step(params.optimizer_options.max_move);

      setCIEnabled(params.neb_options.climbing_image.enabled);
    }

    iteration++;

    double dE = path[maxEnergyImage]->getPotentialEnergy() -
                path[0]->getPotentialEnergy();
    double stepSize = eonc::helpers::maxAtomMotionV(
        path[0]->pbcV(objf->getPositions() - pos));
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
  double fmax = 0;

  // Determine which images to check for convergence
  bool ciOnly = params.neb_options.climbing_image.converged_only &&
                ci_enabled_ && climbingImage != 0;
  long iStart = ciOnly ? climbingImage : 1;
  long iEnd = ciOnly ? climbingImage : numImages;

  for (long i = iStart; i <= iEnd; i++) {
    if (params.optimizer_options.convergence_metric == "norm") {
      fmax = std::max(fmax, projectedForce[i]->norm());
    } else if (params.optimizer_options.convergence_metric == "max_atom") {
      for (int j = 0; j < path[0]->numberOfAtoms(); j++) {
        if (path[0]->getFixed(j))
          continue;
        fmax = std::max(fmax, projectedForce[i]->row(j).norm());
      }
    } else if (params.optimizer_options.convergence_metric == "max_component") {
      fmax = std::max(fmax, projectedForce[i]->maxCoeff());
    } else {
      log = eonc::log::traceback();
      QUILL_LOG_CRITICAL(
          log, "[Nudged Elastic Band] unknown opt_convergence_metric: {}",
          params.optimizer_options.convergence_metric);
      std::exit(1);
    }
  }
  return fmax;
}

// Update the forces, do the projections, and add spring forces
void NudgedElasticBand::updateForces(bool ci_active) {
  // Update forces for all intermediate images.
  // Each image has its own Matter+Potential, so force evaluations are
  // independent. Parallelize when: (a) potential is thread-safe on same
  // instance, OR (b) per-image instances were created (separate models).
  bool canParallel = pot->isThreadSafe() || perImagePotentials_;
  if (numImages > 1 && params.main_options.parallel && canParallel) {
    // std::thread rather than std::jthread -- Apple Clang libc++ lacks the
    // latter. Wrap launch + join so a throw from any lambda still joins the
    // remaining threads before we rethrow; otherwise the unjoined std::thread
    // destructors call std::terminate().
    std::vector<std::thread> threads;
    threads.reserve(static_cast<size_t>(numImages));
    try {
      for (long i = 1; i <= numImages; i++) {
        threads.emplace_back([this, i] { path[i]->getForcesRaw(); });
      }
      for (auto &t : threads)
        t.join();
    } catch (...) {
      for (auto &t : threads)
        if (t.joinable())
          t.join();
      throw;
    }
  } else {
    for (long i = 1; i <= numImages; i++) {
      path[i]->getForcesRaw();
    }
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

  // Update E_ref for energy weighting
  if (params.neb_options.spring.weighting.enabled) {
    E_ref = std::min(path[0]->getPotentialEnergy(),
                     path[numImages + 1]->getPotentialEnergy());
  }

  if (!ci_active) {
    climbingImage = 0;
  }

  // Spring strategy must be rebuilt each iteration (depends on maxEnergy,
  // E_ref). Tangent and projection strategies are cached as members.
  auto spring = eonc::neb::buildSpringStrategy(params, path, numImages, atoms,
                                               maxEnergy, E_ref);

  // Pre-allocate temporaries outside the loop to avoid repeated heap
  // allocation of Nx3 matrices (each ~8KB for 337 atoms).
  AtomMatrix posDiffNext(atoms, 3), posDiffPrev(atoms, 3);

  for (long i = 1; i <= numImages; i++) {
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
            this->ksp = s.ksp;
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
    if (ci_active && i == static_cast<long>(maxEnergyImage)) {
      climbingImage = maxEnergyImage;
      // CI force: F - 2*(F.t)*t, plus DNEB correction if active
      AtomMatrix forceDNEB = AtomMatrix::Zero(atoms, 3);
      if (std::holds_alternative<eonc::neb::DNEB_Projection>(
              projectionStrat_)) {
        AtomMatrix fPerp = eonc::neb::forcePerp(force, *tangent[i]);
        forceDNEB = eonc::neb::computeDNEBComponent(springResult.forceSpring,
                                                    *tangent[i], fPerp);
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

    eonc::neb::zeroTranslation(*projectedForce[i], path[i]->numberOfFreeAtoms(),
                               path[i]->numberOfAtoms());
  }

  movedAfterForceCall = false;
}

// Thin wrappers delegating to eonc::neb:: free functions

void NudgedElasticBand::printImageData(bool writeToFile, size_t idx) {
  eonc::neb::printImageData(path, tangent, eigenmode_solvers, numImages,
                            params.debug_options.estimate_neb_eigenvalues,
                            writeToFile, idx, log);
}

void NudgedElasticBand::findExtrema() {
  auto result = eonc::neb::findSplineExtrema(path, tangent, numImages);
  numExtrema = result.numExtrema;
  extremumPosition = std::move(result.positions);
  extremumEnergy = std::move(result.energies);
  extremumCurvature = std::move(result.curvatures);
}
