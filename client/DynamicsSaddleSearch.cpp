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

#include "DynamicsSaddleSearch.h"

using eonc::SaddleStatus;

#include "BondBoost.h"

#include "Dynamics.h"

#include "EigenmodeStrategy.h"

#include "MinModeSaddleSearch.h"

#include "NudgedElasticBand.h"

#include <cmath>
#include <filesystem>
#include <limits>

SaddleStatus DynamicsSaddleSearch::run() {
  std::vector<std::shared_ptr<Matter>> mdSnapshots;
  std::vector<double> mdTimes;
  QUILL_LOG_DEBUG(log, "Starting dynamics NEB saddle search");

  if (std::filesystem::exists("masses.dat")) {
    QUILL_LOG_DEBUG(log, "Found mass weights file");
    Eigen::VectorXd masses =
        eonc::helpers::loadMasses("masses.dat", saddle->numberOfAtoms());
    saddle->setMasses(masses);
    QUILL_LOG_DEBUG(log, "Applied mass weights");
  } else {
    QUILL_LOG_DEBUG(log, "No mass weights file found");
  }

  Dynamics dyn(saddle.get(), DynamicsConfig::fromParams(params));
  QUILL_LOG_DEBUG(
      log, "Initializing velocities from Maxwell-Boltzmann distribution");
  dyn.setTemperature(params.saddle_search_options.dynamics.temperature);
  dyn.setThermalVelocity();

  int dephaseSteps =
      static_cast<int>(std::floor(params.parallel_replica_options.dephase_time /
                                      params.dynamics_options.time_step +
                                  0.5));

  while (true) {

    QUILL_LOG_DEBUG(log, "Dephasing: {} steps", dephaseSteps);
    // always start from the initial configuration
    *saddle = *reactant;
    dyn.setThermalVelocity();

    // Dephase MD trajectory
    for (int step = 1; step <= dephaseSteps; step++) {
      dyn.oneStep(step);
    }

    // Check to see if a transition occured
    Matter min(pot, params);
    min = *saddle;
    min.relax();

    if (min.compare(*reactant)) {
      QUILL_LOG_DEBUG(log, "Dephasing successful");
      break;
    } else {
      QUILL_LOG_DEBUG(log, "Transition occured during dephasing; Restarting");
      dephaseSteps /= 2;
      if (dephaseSteps < 1)
        dephaseSteps = 1;
    }
  }

  BondBoost bondBoost(saddle.get(), params);
  if (params.hyperdynamics_options.bias_potential ==
      Hyperdynamics::BOND_BOOST) {
    QUILL_LOG_DEBUG(log, "Initializing Bond Boost");
    bondBoost.initialize();
  }

  // Round-to-nearest with std::lround instead of `(int)(x + 0.5)`,
  // which is broken for negative arguments and exact representations
  // (bugprone-incorrect-roundings).
  const auto checkInterval = static_cast<int>(
      std::lround(params.saddle_search_options.dynamics.state_check_interval /
                  params.dynamics_options.time_step));
  const auto recordInterval = static_cast<int>(
      std::lround(params.saddle_search_options.dynamics.record_interval /
                  params.dynamics_options.time_step));

  if (params.debug_options.write_movies) {
    saddle->matter2con("dynamics", false);
  }

  for (int step = 1; step <= params.dynamics_options.steps; step++) {
    dyn.oneStep(step);

    if (recordInterval != 0 && step % recordInterval == 0) {
      QUILL_LOG_DEBUG(
          log, "recording configuration at step {} time {:.3f}", step,
          step * params.dynamics_options.time_step * params.constants.timeUnit);
      // BUG FIX: was sharing ownership with saddle instead of copying
      auto snapshot = std::make_shared<Matter>(*saddle);
      mdSnapshots.push_back(snapshot);
      mdTimes.push_back(step * params.dynamics_options.time_step);
    }

    if (params.debug_options.write_movies) {
      saddle->matter2con("dynamics", true);
    }

    if (step % checkInterval == 0) {
      QUILL_LOG_DEBUG(log, "Minimizing trajectory, step {}", step);

      product = std::make_shared<Matter>(*saddle);
      product->relax(false, false);

      if (!product->compare(*reactant)) {
        QUILL_LOG_DEBUG(log, "Found new state");
        int image = refineTransition(mdSnapshots, product);
        *saddle = *mdSnapshots[image];
        QUILL_LOG_DEBUG(log, "Found transition at snapshot image {}", image);
        for (int ii = 0; ii < static_cast<int>(mdTimes.size()); ii++) {
          QUILL_LOG_DEBUG(log, "MDTimes[{}] = {:.3f}", ii,
                          mdTimes[ii] * params.constants.timeUnit);
        }
        // Subtract half the record interval to avoid systematic bias
        time = mdTimes[image] -
               params.saddle_search_options.dynamics.record_interval / 2.0;
        QUILL_LOG_DEBUG(log, "Transition time {:.2f} fs",
                        time * params.constants.timeUnit);

        NudgedElasticBand neb(reactant, product, params, pot);

        if (!params.saddle_search_options.dynamics.linear_interpolation) {
          QUILL_LOG_DEBUG(
              log, "Interpolating initial band through MD transition state");
          AtomMatrix reactantToSaddle =
              saddle->pbc(saddle->getPositions() - reactant->getPositions());
          AtomMatrix saddleToProduct =
              saddle->pbc(product->getPositions() - saddle->getPositions());
          QUILL_LOG_DEBUG(log, "Initial band saved to neb_initial_band.con");
          neb.path[0]->matter2con("neb_initial_band.con", false);
          int mid = neb.numImages / 2 + 1;
          for (int img = 1; img <= neb.numImages; img++) {
            if (img < mid) {
              double frac = static_cast<double>(img) / static_cast<double>(mid);
              neb.path[img]->setPositions(reactant->getPositions() +
                                          frac * reactantToSaddle);
            } else if (img > mid) {
              double frac = static_cast<double>(img - mid) /
                            static_cast<double>(neb.numImages - mid + 1);
              neb.path[img]->setPositions(saddle->getPositions() +
                                          frac * saddleToProduct);
            } else {
              neb.path[img]->setPositions(saddle->getPositions());
            }
            neb.path[img]->matter2con("neb_initial_band.con", true);
          }
          neb.path[neb.numImages + 1]->matter2con("neb_initial_band.con", true);
        } else {
          QUILL_LOG_DEBUG(
              log, "Linear interpolation between minima used for initial band");
          neb.path[0]->matter2con("neb_initial_band.con", false);
          for (int j = 1; j <= neb.numImages + 1; j++) {
            neb.path[j]->matter2con("neb_initial_band.con", true);
          }
        }

        AtomMatrix mode;
        if (params.neb_options.max_iterations > 0) {
          auto minModeMethod =
              eonc::buildEigenmodeStrategy(saddle, params, pot);

          neb.compute();
          neb.printImageData(true);
          int extremumImage = -1;
          int jExt = 0;
          for (jExt = 0; jExt < neb.numExtrema; jExt++) {
            if (neb.extremumCurvature[jExt] <
                params.saddle_search_options.dynamics.max_init_curvature) {
              extremumImage =
                  static_cast<int>(std::floor(neb.extremumPosition[jExt]));
              *saddle = *neb.path[extremumImage];
              double interpDist = neb.extremumPosition[jExt] -
                                  static_cast<double>(extremumImage);
              AtomMatrix bandDir =
                  saddle->pbc(neb.path[extremumImage + 1]->getPositions() -
                              neb.path[extremumImage]->getPositions());
              saddle->setPositions(interpDist * bandDir +
                                   saddle->getPositions());
              mode = saddle->pbc(neb.path[extremumImage + 1]->getPositions() -
                                 saddle->getPositions());
              mode.normalize();
              eonc::eigenmodeCompute(*minModeMethod, saddle, mode);
              double ev = eonc::eigenmodeGetEigenvalue(*minModeMethod);
              QUILL_LOG_DEBUG(log, "extrema #{} has eigenvalue {:.8f}",
                              jExt + 1, ev);

              if (ev < 0) {
                QUILL_LOG_DEBUG(
                    log, "chose image {} (extrema #{}) as extremum image",
                    extremumImage, jExt + 1);
                break;
              } else {
                extremumImage = -1;
              }
            }
          }

          if (extremumImage != -1) {
            *saddle = *neb.path[extremumImage];
            double interpDist =
                neb.extremumPosition[jExt] - static_cast<double>(extremumImage);
            QUILL_LOG_DEBUG(log, "interpDistance {}", interpDist);
            AtomMatrix bandDir =
                saddle->pbc(neb.path[extremumImage + 1]->getPositions() -
                            neb.path[extremumImage]->getPositions());
            saddle->setPositions(interpDist * bandDir + saddle->getPositions());
            mode = saddle->pbc(neb.path[extremumImage + 1]->getPositions() -
                               saddle->getPositions());
            mode.normalize();
          } else {
            QUILL_LOG_DEBUG(
                log, "no maxima found, using max energy non-endpoint image");
            double maxEnergy = -std::numeric_limits<double>::infinity();
            for (int img = 1; img <= neb.numImages; img++) {
              double U = neb.path[img]->getPotentialEnergy();
              if (U > maxEnergy) {
                maxEnergy = U;
                *saddle = *neb.path[img];
                mode = saddle->pbc(neb.path[img + 1]->getPositions() -
                                   saddle->getPositions());
                mode.normalize();
              }
            }
            if (maxEnergy <= reactant->getPotentialEnergy()) {
              QUILL_LOG_DEBUG(log, "warning: no barrier found");
              return SaddleStatus::BadNoBarrier;
            }
          }
        } else {
          neb.maxEnergyImage = neb.numImages / 2 + 1;
        }

        QUILL_LOG_DEBUG(
            log, "Initial saddle guess saved to saddle_initial_guess.con");
        saddle->matter2con("saddle_initial_guess.con");
        MinModeSaddleSearch search = MinModeSaddleSearch(
            saddle, mode, reactant->getPotentialEnergy(), params, pot);
        SaddleStatus minModeStatus = search.run();

        if (minModeStatus != SaddleStatus::Good) {
          QUILL_LOG_DEBUG(log, "error in min mode saddle search");
          return minModeStatus;
        }

        eigenvalue = search.getEigenvalue();
        eigenvector = search.getEigenvector();
        QUILL_LOG_DEBUG(log, "eigenvalue: {:.3f}", eigenvalue);

        double barrier =
            saddle->getPotentialEnergy() - reactant->getPotentialEnergy();
        QUILL_LOG_DEBUG(log, "found barrier of {:.3f}", barrier);
        mdSnapshots.clear();
        mdTimes.clear();
        return SaddleStatus::Good;
      } else {
        QUILL_LOG_DEBUG(log, "Still in original state");
        mdTimes.clear();
        mdSnapshots.clear();
      }
    }
  }

  mdSnapshots.clear();
  time = params.dynamics_options.steps * params.dynamics_options.time_step;
  return SaddleStatus::BadMdTrajectoryTooShort;
}

/// Binary search through MD snapshots to find the transition point.
int DynamicsSaddleSearch::refineTransition(
    const std::vector<std::shared_ptr<Matter>> &snapshots,
    const std::shared_ptr<Matter> &prod) {
  int lo = 0;
  int hi = static_cast<int>(snapshots.size()) - 1;
  if (hi == 0) {
    return 0;
  }

  QUILL_LOG_DEBUG(log, "refining transition time");

  while ((hi - lo) > 1) {
    int mid = lo + (hi - lo) / 2;
    QUILL_LOG_DEBUG(log, "minimizing image {}", mid);
    Matter snap(pot, params);
    snap = *snapshots[mid];
    snap.relax(false);

    if (snap.compare(*reactant)) {
      QUILL_LOG_DEBUG(log, "image {} minimizes to reactant", mid);
      lo = mid;
    } else {
      QUILL_LOG_DEBUG(log, "image {} minimizes to product", mid);
      *prod = snap;
      hi = mid;
    }
  }

  return (lo + hi) / 2;
}

double DynamicsSaddleSearch::getEigenvalue() { return eigenvalue; }

AtomMatrix DynamicsSaddleSearch::getEigenvector() { return eigenvector; }
