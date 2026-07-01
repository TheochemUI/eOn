/*
#include <stdexcept>
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
#include "ParallelReplicaJob.h"
#include "BaseStructures.h"
#include "BondBoost.h"
#include "Dynamics.h"
#include "ForceCallTimer.h"
#include "HelperFunctions.h"
#include "Matter.h"

#include <cmath>
#include <format>
#include <fstream>

std::vector<std::string> ParallelReplicaJob::run() {
  reactant = std::make_shared<Matter>(pot, params);
  {
    const auto posIn =
        eonc::helpers::getRelevantFile(params.main_options.conFilename);
    if (!eonc::io::io_ok(reactant->con2matter(posIn))) {
      QUILL_LOG_CRITICAL(log, "Failed to load {}", posIn);
      throw std::runtime_error("failed to load " + posIn);
    }
  }

  QUILL_LOG_DEBUG(log, "[ParallelReplica] Minimizing initial position");
  reactant->relax();
  if (!eonc::io::io_ok(reactant->matter2con("reactant.con"))) {
    QUILL_LOG_ERROR(log, "Failed to write reactant.con");
  }

  auto trajectory = std::make_shared<Matter>(pot, params);
  *trajectory = *reactant;
  Dynamics dynamics(trajectory.get(), params);
  BondBoost bondBoost(trajectory.get(), params);

  if (params.hyperdynamics_options.bias_potential ==
      Hyperdynamics::BOND_BOOST) {
    bondBoost.initialize();
    trajectory->setBiasPotential(&bondBoost);
  }

  dephase(*trajectory);

  int stateCheckInterval = static_cast<int>(
      std::floor(params.parallel_replica_options.state_check_interval /
                     params.dynamics_options.time_step +
                 0.5));
  int recordInterval = static_cast<int>(
      std::floor(params.parallel_replica_options.record_interval /
                     params.dynamics_options.time_step +
                 0.5));

  std::vector<std::shared_ptr<Matter>> mdSnapshots;
  std::vector<double> mdTimes;
  double transitionTime = 0;
  Matter transitionStructure(pot, params);
  size_t refineForceCalls = 0;

  double simulationTime = 0.0;
  if (params.hyperdynamics_options.bias_potential == Hyperdynamics::NONE) {
    QUILL_LOG_DEBUG(
        log, "[ParallelReplica] {:>8} {:>12} {:>10} {:>12} {:>12} {:>10}",
        "Step", "Time (s)", "KE", "PE", "TE", "KinT");
  } else {
    QUILL_LOG_DEBUG(
        log,
        "[ParallelReplica] {:>8} {:>12} {:>10} {:>10} {:>12} {:>12} {:>10}",
        "Step", "Time (s)", "Boost", "KE", "PE", "TE", "KinT");
  }

  for (int step = 1; step <= params.dynamics_options.steps; step++) {
    dynamics.oneStep();
    double boost = 1.0;
    if (params.hyperdynamics_options.bias_potential ==
        Hyperdynamics::BOND_BOOST) {
      double boostPotential = bondBoost.boost();
      double kB = params.constants.kB;
      boost = std::exp(boostPotential / kB / params.main_options.temperature);
      simulationTime += params.dynamics_options.time_step * boost;
    } else {
      simulationTime += params.dynamics_options.time_step;
    }

    double kinE = trajectory->getKineticEnergy();
    double potE = trajectory->getPotentialEnergy();
    double kinT = (2.0 * kinE / (trajectory->numberOfFreeAtoms() * 3) /
                   params.constants.kB);

    if (step % params.debug_options.write_movies_interval == 0) {
      if (params.hyperdynamics_options.bias_potential == Hyperdynamics::NONE) {
        QUILL_LOG_DEBUG(log,
                        "[ParallelReplica] {:>8} {:>12.4e} {:>10.4f} "
                        "{:>12.4f} {:>12.4f} {:>10.2f}",
                        step,
                        simulationTime * params.constants.timeUnit * 1e-15,
                        kinE, potE, kinE + potE, kinT);
      } else {
        double boostPotential = bondBoost.boost();
        QUILL_LOG_DEBUG(
            log,
            "[ParallelReplica] {:>8} {:>12.4e} {:>10.3e} "
            "{:>10.4f} {:>12.4f} {:>12.4f} {:>10.2f}",
            step, simulationTime * params.constants.timeUnit * 1e-15, boost,
            kinE, potE + boostPotential, kinE + potE + boostPotential, kinT);
      }
    }

    // Snapshots for refinement
    if (step % recordInterval == 0 &&
        params.parallel_replica_options.refine_transition) {
      auto snap = std::make_shared<Matter>(pot, params);
      *snap = *trajectory;
      mdSnapshots.push_back(std::move(snap));
      mdTimes.push_back(simulationTime);
    }

    // Check for transition
    if (step % stateCheckInterval == 0 ||
        step == params.dynamics_options.steps) {
      QUILL_LOG_DEBUG(log, "[ParallelReplica] Checking for transition");

      Matter minimized(pot, params);
      minimized = *trajectory;
      minimized.relax();

      if (!minimized.compare(*reactant) && transitionTime == 0) {
        QUILL_LOG_DEBUG(log, "[ParallelReplica] Transition occurred");

        if (params.parallel_replica_options.refine_transition) {
          QUILL_LOG_DEBUG(log, "[ParallelReplica] Refining transition time");
          int snapshotIndex;
          {
            eonc::ForceCallTimer timer(refineForceCalls);
            snapshotIndex = refineTransition(mdSnapshots);
          }

          transitionTime = mdTimes[snapshotIndex];
          transitionStructure = *mdSnapshots[snapshotIndex];
        } else {
          transitionStructure = *trajectory;
          transitionTime = simulationTime;
        }
        QUILL_LOG_DEBUG(log, "[ParallelReplica] Transition time: {:.3e} s",
                        transitionTime * params.constants.timeUnit * 1e-15);

      } else if (step + 1 == params.dynamics_options.steps &&
                 transitionTime == 0) {
        // Fake refinement to prevent force-call bias
        if (params.parallel_replica_options.refine_transition) {
          QUILL_LOG_DEBUG(
              log,
              "[ParallelReplica] Simulation ended without seeing a transition");
          QUILL_LOG_DEBUG(
              log, "[ParallelReplica] Refining anyways to prevent bias...");
          {
            eonc::ForceCallTimer timer(refineForceCalls);
            refineTransition(mdSnapshots, true);
          }
        }
        transitionStructure = *trajectory;
      }

      mdSnapshots.clear();
      mdTimes.clear();
    }
  }

  // Decorrelation dynamics
  int decorrelationSteps =
      static_cast<int>(std::floor(params.parallel_replica_options.corr_time /
                                      params.dynamics_options.time_step +
                                  0.5));
  QUILL_LOG_DEBUG(log, "[ParallelReplica] Decorrelating: {} steps",
                  decorrelationSteps);
  for (int step = 1; step <= decorrelationSteps; step++) {
    dynamics.oneStep(step);
  }
  QUILL_LOG_DEBUG(log, "[ParallelReplica] Decorrelation complete");

  // Minimize final structure
  Matter product(pot, params);
  product = *trajectory;
  product.relax();
  if (!eonc::io::io_ok(product.matter2con("product.con"))) {
    QUILL_LOG_ERROR(log, "Failed to write product.con");
  }

  // Write results
  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  {
    std::ofstream out(resultsFilename, std::ios::binary);
    if (out) {
      out << std::format(
          "{} potential_type\n",
          magic_enum::enum_name<PotType>(params.potential_options.potential));
      out << std::format("{} random_seed\n", params.main_options.randomSeed);
      out << std::format("{:f} potential_energy_reactant\n",
                         reactant->getPotentialEnergy());
      out << std::format("{} force_calls_refine\n", refineForceCalls);
      out << std::format("{} total_force_calls\n",
                         PotRegistry::get().total_force_calls());

      if (transitionTime == 0) {
        out << "0 transition_found\n";
        out << std::format("{:e} simulation_time_s\n",
                           simulationTime * params.constants.timeUnit *
                               1.0e-15);
      } else {
        out << "1 transition_found\n";
        out << std::format("{:e} transition_time_s\n",
                           transitionTime * params.constants.timeUnit *
                               1.0e-15);
        out << std::format("{:e} correlation_time_s\n",
                           params.parallel_replica_options.corr_time *
                               params.constants.timeUnit * 1.0e-15);
        out << std::format("{:f} potential_energy_product\n",
                           product.getPotentialEnergy());
      }
      out << std::format("{:f} speedup\n",
                         simulationTime / (params.dynamics_options.steps *
                                           params.dynamics_options.time_step));
    }
  }

  return returnFiles;
}

void ParallelReplicaJob::dephase(Matter &trajectory) {
  Dynamics dynamics(&trajectory, params);

  int dephaseSteps =
      static_cast<int>(std::floor(params.parallel_replica_options.dephase_time /
                                      params.dynamics_options.time_step +
                                  0.5));
  QUILL_LOG_DEBUG(log, "[ParallelReplica] Dephasing: {} steps", dephaseSteps);

  Matter initial(pot, params);
  initial = trajectory;

  while (true) {
    trajectory = initial;
    dynamics.setThermalVelocity();

    for (int step = 1; step <= dephaseSteps; step++) {
      dynamics.oneStep(step);
    }

    Matter minimized(pot, params);
    minimized = trajectory;
    minimized.relax();

    if (minimized.compare(*reactant)) {
      QUILL_LOG_DEBUG(log, "[ParallelReplica] Dephasing successful");
      break;
    } else {
      QUILL_LOG_DEBUG(
          log,
          "[ParallelReplica] Transition occured during dephasing; Restarting");
    }
  }
}

int ParallelReplicaJob::refineTransition(
    const std::vector<std::shared_ptr<Matter>> &snapshots, bool fake) {
  int lo = 0;
  int hi = static_cast<int>(snapshots.size()) - 1;

  while ((hi - lo) > 1) {
    int mid = lo + (hi - lo) / 2;
    Matter snapshot(*snapshots[mid]);
    snapshot.relax(true);

    bool stillReactant;
    if (!fake) {
      stillReactant = snapshot.compare(*reactant);
    } else {
      stillReactant = static_cast<bool>(eonc::helpers::randomInt(0, 1));
    }

    if (stillReactant) {
      lo = mid;
    } else {
      hi = mid;
    }
  }

  return (lo + hi) / 2 + 1;
}
