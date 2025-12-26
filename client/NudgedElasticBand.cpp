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
          [&, initialPassed]() {
            switch (parametersPassed->neb_options.initialization.method) {
            case NEBInit::FILE: {
              if (parametersPassed->neb_options.initialization.input_path
                      .empty()) {
                throw std::runtime_error(
                    "NEB initializer set to FILE, but neb_ipath is empty.");
              }
              std::vector<fs::path> file_paths =
                  helper_functions::neb_paths::readFilePaths(
                      parametersPassed->neb_options.initialization.input_path);
              return helper_functions::neb_paths::filePathInit(
                  file_paths, *initialPassed,
                  parametersPassed->neb_options.image_count);
            }
            case NEBInit::IDPP: {
              return helper_functions::neb_paths::idppPath(
                  *initialPassed, *finalPassed,
                  parametersPassed->neb_options.image_count, parametersPassed);
            }
            case NEBInit::IDPP_COLLECTIVE: {
              return helper_functions::neb_paths::idppCollectivePath(
                  *initialPassed, *finalPassed,
                  parametersPassed->neb_options.image_count, parametersPassed);
            }
            case NEBInit::SIDPP: {
              return helper_functions::neb_paths::sidppPath(
                  *initialPassed, *finalPassed,
                  parametersPassed->neb_options.image_count, parametersPassed);
            }
            case NEBInit::LINEAR:
            default: {
              return helper_functions::neb_paths::linearPath(
                  *initialPassed, *finalPassed,
                  parametersPassed->neb_options.image_count);
            }
            }
          }(),
          parametersPassed, potPassed) {}

// Second constructor: Contains all the common setup code
NudgedElasticBand::NudgedElasticBand(
    std::vector<Matter> initPath, std::shared_ptr<Parameters> parametersPassed,
    std::shared_ptr<Potential> potPassed)
    : params{parametersPassed},
      pot{potPassed} {

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

  SPDLOG_LOGGER_DEBUG(log, "Nudged elastic band calculation started.");
  if (params->neb_options.spring.weighting.enabled) {
    E_ref = std::min(path[0]->getPotentialEnergy(),
                     path[numImages + 1]->getPotentialEnergy());
  }
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
  SPDLOG_DEBUG("{:>10s} {:>12s} {:>14s} {:>11s} {:>12s}", "iteration",
               "step size", params->optConvergenceMetricLabel, "max image",
               "max energy");
  SPDLOG_DEBUG(
      "---------------------------------------------------------------\n");

  while (this->status != NEBStatus::GOOD) {
    if (params->writeMovies) {
      bool append = true;
      if (iteration == 0) {
        append = false;
      }
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
    double convForce{convergenceForce()};
    bool originalCIflag = params->neb_options.climbing_image.enabled;
    if (convForce >= params->neb_options.climbing_image.trigger_force) {
      params->neb_options.climbing_image.enabled = false;
    }
    if (iteration) {
      // so that we print forces before taking an optimizer step
      if (climbingImage > 0 && climbingImage <= numImages &&
          params->neb_options.climbing_image.roneb.use_mmf &&
          params->neb_options.climbing_image.enabled &&
          convForce < current_mmf_threshold && iteration > 1) {

        SPDLOG_LOGGER_DEBUG(
            log, "Triggering MMF. Current Force: {:.4f}, Threshold: {:.4f}",
            convForce, current_mmf_threshold);

        auto tempMinModeSearch = std::make_shared<MinModeSaddleSearch>(
            path[climbingImage], *tangent[climbingImage],
            path[climbingImage]->getPotentialEnergy(), params, pot);

        // Save the original value of saddleMaxIterations
        auto originalSaddleMaxIterations = params->saddleMaxIterations;
        params->saddleMaxIterations =
            params->neb_options.climbing_image.roneb.max_steps;
        int minModeStatus;
        try {
          minModeStatus = tempMinModeSearch->run();
        } catch (const DimerModeLostException &e) {
          // Determine we crashed because of mode loss
          minModeStatus = MinModeSaddleSearch::STATUS_DIMER_LOST_MODE;
        }

        // Restore the original value of saddleMaxIterations
        params->saddleMaxIterations = originalSaddleMaxIterations;
        if (minModeStatus != MinModeSaddleSearch::STATUS_GOOD &&
            minModeStatus != MinModeSaddleSearch::STATUS_BAD_MAX_ITERATIONS) {
          SPDLOG_LOGGER_WARN(log,
                             "MMF Failed (Status {}). Adjusting threshold "
                             "based on Mode-Tangent alignment.",
                             static_cast<int>(minModeStatus));

          // Retrieve the calculated mode from the failed search
          AtomMatrix finalModeMatrix = tempMinModeSearch->getEigenvector();
          VectorXd finalMode = VectorXd::Map(finalModeMatrix.data(), 3 * atoms);

          // Retrieve the current NEB tangent
          VectorXd currentTangent =
              VectorXd::Map(tangent[climbingImage]->data(), 3 * atoms);

          // Calculate alignment (absolute dot product)
          // We use absolute value because the mode direction sign is arbitrary
          double alignment =
              std::abs(finalMode.normalized().dot(currentTangent.normalized()));

          // Define a penalty scaling factor based on alignment.
          // High alignment (1.0) -> Mild penalty (0.8x)
          // Low alignment (0.0) -> Harsh penalty (0.1x)
          // This uses a sigmoid-like or linear blend.
          double penalty_factor =
              params->neb_options.climbing_image.roneb.penalty.base +
              (params->neb_options.climbing_image.roneb.penalty.strength *
               alignment);

          double old_threshold = current_mmf_threshold;

          // Set new threshold based on the CURRENT force and the physical
          // alignment We discard the 'std::min' history because the physics of
          // the current configuration dictates the trust radius.
          current_mmf_threshold = convForce * penalty_factor;

          SPDLOG_LOGGER_WARN(
              log, "Alignment: {:.3f}. New Threshold: {:.4f} (Factor {:.2f})",
              alignment, current_mmf_threshold, penalty_factor);

        } else {
          // If it worked (or just ran out of steps without crashing),
          // we can relax the threshold back to the original parameter
          // just in case we need it again later (hysteresis).
          current_mmf_threshold =
              params->neb_options.climbing_image.roneb.trigger_force;
        }

        // Now, the climbing image (path[climbingImage]) has been optimized for
        // a few steps by the MinModeSaddleSearch.
        movedAfterForceCall = true;
        updateForces();
      }

      if (iteration >= params->neb_options.max_iterations) {
        status = NEBStatus::BAD_MAX_ITERATIONS;
        break;
      }
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
    }
    iteration++;
    params->neb_options.climbing_image.enabled = originalCIflag;

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
  // bool higherEnergyPrev, higherEnergyNext;

  // variables for climbing image
  double maxEnergy;

  // variables for the energy weighted springs
  k_l = params->neb_options.spring.weighting.k_min;
  k_u = params->neb_options.spring.weighting.k_max;
  std::vector<double> springConstants(numImages + 1, k_l);

  // variables for force projections
  AtomMatrix force(atoms, 3), forcePerp(atoms, 3), forcePar(atoms, 3);
  AtomMatrix forceSpringPar(atoms, 3), forceSpring(atoms, 3),
      forceSpringPerp(atoms, 3);
  AtomMatrix forceDNEB(atoms, 3);
  AtomMatrix pos(atoms, 3), posNext(atoms, 3), posPrev(atoms, 3),
      posDiffNext(atoms, 3), posDiffPrev(atoms, 3);
  double distNext, distPrev;

  // Update forces for all intermediate images
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
  // Determine if energy weighting applies
  double rawMaxForce = 0.0;
  for (long k = 1; k <= numImages; k++) {
    // Use the raw forces (already updated above) to check convergence proximity
    if (params->optConvergenceMetric == "norm") {
      rawMaxForce = std::max(rawMaxForce, path[k]->getForces().norm());
    } else {
      // Default to max component/atom proxy for safety
      rawMaxForce = std::max(rawMaxForce, path[k]->getForces().maxCoeff());
    }
  }

  double ciTrigger = params->neb_options.climbing_image.trigger_force;
  bool triggerMet = (rawMaxForce < ciTrigger);
  bool weightingActive =
      params->neb_options.spring.weighting.enabled && triggerMet;

  // Energy weighted springs, calculated here since all the springs are used
  // internally
  if (weightingActive) {
    for (int idx = 1; idx <= numImages + 1; idx++) {
      double Ei = std::max(path[idx]->getPotentialEnergy(),
                           path[idx - 1]->getPotentialEnergy());
      if (Ei > E_ref) {
        double alpha_i = (maxEnergy - Ei) / (maxEnergy - E_ref);
        springConstants[idx - 1] =
            (1 - alpha_i) * k_u + alpha_i * k_l; // Equation (3) and (4)
      } // else always k_l
    }
  } else {
    // If CI is not yet active, use the standard uniform spring constant
    std::fill(springConstants.begin(), springConstants.end(),
              params->neb_options.spring.constant);
  }

  for (long i = 1; i <= numImages; i++) {
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
    tangent[i]->normalize();

    // Calculate the force perpendicular to the tangent
    forcePerp =
        force - (force.array() * (*tangent[i]).array()).sum() * *tangent[i];

    // Calculate spring forces with the conditional trigger
    if (weightingActive) {
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

    // Doubly Nudged Elastic Band logic
    if (params->neb_options.spring.doubly_nudged) {
      if (!weightingActive) {
        forceSpringPerp =
            forceSpring -
            (forceSpring.array() * (*tangent[i]).array()).sum() * *tangent[i];
        forceDNEB =
            forceSpringPerp -
            (forceSpringPerp.array() * forcePerp.normalized().array()).sum() *
                forcePerp.normalized();

        double switching =
            2.0 / M_PI *
            atan(pow(forcePerp.norm(), 2) / pow(forceSpringPerp.norm(), 2));
        forceDNEB *= switching;
      } else {
        SPDLOG_WARN("Not using doubly nudged since energy_weighted is set");
      }
    } else {
      forceDNEB.setZero();
    }

    // Apply the Climbing Image projection or standard NEB force
    // Note: CI only activates if both the parameter is enabled AND the trigger
    // is met
    if (params->neb_options.climbing_image.enabled && triggerMet &&
        i == maxEnergyImage) {
      climbingImage = maxEnergyImage;
      *projectedForce[i] =
          force -
          (2.0 * (force.array() * (*tangent[i]).array()).sum() * *tangent[i]) +
          forceDNEB;
    } else // all non-climbing numImages
    {
      // sum the spring and potential forces for the neb force
      if (params->neb_options.spring.use_elastic_band && !weightingActive) {
        *projectedForce[i] = forceSpring + force;
      } else {
        *projectedForce[i] = forceSpringPar + forcePerp + forceDNEB;
      }
      //*projectedForce[i] = forceSpring + forcePerp;

      // if (params->nebFullSpring) {

      movedAfterForceCall = false; // so that we don't repeat a force call
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
