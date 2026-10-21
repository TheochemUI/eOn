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
#include "eon/ParallelReplicaJob.h"
#include "eon/BaseStructures.h"
#include "eon/BondBoost.h"
#include "eon/ConjugateGradients.h"
#include "eon/Dynamics.h"
#include "eon/ForceCallTimer.h"
#include "eon/HelperFunctions.h"
#include "eon/Matter.h"
#include <functional>
#include <stdexcept>

namespace eonc::helpers {
bool relaxMatterAccepted(Matter &matter, const Parameters &params, bool quiet,
                         bool writeMovie, bool checkpoint,
                         std::string prefixMovie,
                         std::string prefixCheckpoint,
                         std::vector<readcon::ConFrame> *outFrames,
                         std::function<bool(Matter &, double)> accept);
}

#include <cmath>
#include <format>
#include <fstream>
#include <memory>

namespace eonc {

std::vector<std::string> ParallelReplicaJob::run() {
  reactant = std::make_shared<Matter>(pot, params);
  {
    const auto posIn =
        eonc::helpers::getRelevantFile(params.main_options().conFilename);
    if (!eonc::io::io_ok(reactant->con2matter(posIn))) {
      QUILL_LOG_CRITICAL(log, "Failed to load {}", posIn);
      throw std::runtime_error("failed to load " + posIn);
    }
  }
  static_cast<void>(runFromMatter(reactant));
  return returnFiles;
}

std::shared_ptr<Matter>
ParallelReplicaJob::runFromMatter(std::shared_ptr<Matter> initial) {
  if (!initial) {
    throw std::runtime_error("ParallelReplicaJob::runFromMatter: null Matter");
  }
  reactant = initial;
  reactant->setPotential(pot);
  profile_ = {};

  QUILL_LOG_DEBUG(log, "[ParallelReplica] Minimizing initial position");
  {
    const size_t at = pot->forceCallCounter.load();
    reactant->relax();
    if (!eonc::io::io_ok(reactant->matter2con("reactant.con"))) {
      QUILL_LOG_ERROR(log, "Failed to write reactant.con");
    }
    profile_.initial = pot->forceCallCounter.load() - at;
  }

  auto trajectory = std::make_shared<Matter>(pot, params);
  *trajectory = *reactant;
  Dynamics dynamics(trajectory.get(), params);
  BondBoost bondBoost(trajectory.get(), params);

  BondBoost *bias = nullptr;
  if (params.hyperdynamics_options().bias_potential ==
      Hyperdynamics::BOND_BOOST) {
    bondBoost.initialize();
    trajectory->setBiasPotential(&bondBoost);
    bias = &bondBoost;
  }

  dephase(*trajectory, bias);

  int stateCheckInterval = static_cast<int>(
      std::floor(params.parallel_replica_options().state_check_interval /
                     params.dynamics_options().time_step +
                 0.5));
  int recordInterval = static_cast<int>(
      std::floor(params.parallel_replica_options().record_interval /
                     params.dynamics_options().time_step +
                 0.5));
  if (stateCheckInterval < 1)
    stateCheckInterval = 1;
  if (recordInterval < 1)
    recordInterval = 1;

  std::vector<std::shared_ptr<Matter>> mdSnapshots;
  std::vector<double> mdTimes;
  double transitionTime = 0;
  Matter transitionStructure(pot, params);
  size_t refineForceCalls = 0;

  double simulationTime = 0.0;
  if (params.hyperdynamics_options().bias_potential == Hyperdynamics::NONE) {
    QUILL_LOG_DEBUG(
        log, "[ParallelReplica] {:>8} {:>12} {:>10} {:>12} {:>12} {:>10}",
        "Step", "Time (s)", "KE", "PE", "TE", "KinT");
  } else {
    QUILL_LOG_DEBUG(
        log,
        "[ParallelReplica] {:>8} {:>12} {:>10} {:>10} {:>12} {:>12} {:>10}",
        "Step", "Time (s)", "Boost", "KE", "PE", "TE", "KinT");
  }

  for (int step = 1; step <= params.dynamics_options().steps; step++) {
    const size_t mdAt = pot->forceCallCounter.load();
    if (params.hyperdynamics_options().bias_potential ==
        Hyperdynamics::BOND_BOOST) {
      // oneStep() evaluates accelerations (and therefore boost()) more than
      // once. Advance the equilibration counter here, once per MD step, so
      // ParallelReplica and SafeHyper share the same rmd_time schedule.
      bondBoost.advance();
    }
    dynamics.oneStep();
    double boost = 1.0;
    if (params.hyperdynamics_options().bias_potential ==
        Hyperdynamics::BOND_BOOST) {
      double boostPotential = bondBoost.boost();
      double kB = params.constants().kB;
      boost = std::exp(boostPotential / kB / params.main_options().temperature);
      simulationTime += params.dynamics_options().time_step * boost;
    } else {
      simulationTime += params.dynamics_options().time_step;
    }

    double kinE = trajectory->getKineticEnergy();
    double potE = trajectory->getPotentialEnergy();
    double kinT = (2.0 * kinE / (trajectory->numberOfFreeAtoms() * 3) /
                   params.constants().kB);

    if (step % params.debug_options().write_movies_interval == 0) {
      if (params.hyperdynamics_options().bias_potential ==
          Hyperdynamics::NONE) {
        QUILL_LOG_DEBUG(log,
                        "[ParallelReplica] {:>8} {:>12.4e} {:>10.4f} "
                        "{:>12.4f} {:>12.4f} {:>10.2f}",
                        step,
                        simulationTime * params.constants().timeUnit * 1e-15,
                        kinE, potE, kinE + potE, kinT);
      } else {
        double boostPotential = bondBoost.boost();
        QUILL_LOG_DEBUG(
            log,
            "[ParallelReplica] {:>8} {:>12.4e} {:>10.3e} "
            "{:>10.4f} {:>12.4f} {:>12.4f} {:>10.2f}",
            step, simulationTime * params.constants().timeUnit * 1e-15, boost,
            kinE, potE + boostPotential, kinE + potE + boostPotential, kinT);
      }
    }
    profile_.md += pot->forceCallCounter.load() - mdAt;

    // Snapshots for refinement
    if (step % recordInterval == 0 &&
        params.parallel_replica_options().refine_transition) {
      auto snap = std::make_shared<Matter>(pot, params);
      *snap = *trajectory;
      mdSnapshots.push_back(std::move(snap));
      mdTimes.push_back(simulationTime);
    }

    // Check for transition
    if (step % stateCheckInterval == 0 ||
        step == params.dynamics_options().steps) {
      QUILL_LOG_DEBUG(log, "[ParallelReplica] Checking for transition");

      Matter minimized(pot, params);
      minimized = *trajectory;
      ++profile_.stateChecks;
      const double distTol =
          params.structure_comparison_options().distance_difference;
      if (distTol > 0.0 && minimized.perAtomNorm(*reactant) < distTol) {
        ++profile_.hotInside;
      }
      {
        const size_t quenchAt = pot->forceCallCounter.load();
        basinRelax(minimized, false);
        profile_.stateCheck += pot->forceCallCounter.load() - quenchAt;
      }

      if (!minimized.compare(*reactant) && transitionTime == 0) {
        QUILL_LOG_DEBUG(log, "[ParallelReplica] Transition occurred");

        if (params.parallel_replica_options().refine_transition &&
            !mdSnapshots.empty() && !mdTimes.empty()) {
          QUILL_LOG_DEBUG(log, "[ParallelReplica] Refining transition time");
          int snapshotIndex;
          {
            eonc::ForceCallTimer timer(refineForceCalls, *pot);
            snapshotIndex = refineTransition(mdSnapshots);
          }
          if (snapshotIndex < 0)
            snapshotIndex = 0;
          if (snapshotIndex >= static_cast<int>(mdSnapshots.size()))
            snapshotIndex = static_cast<int>(mdSnapshots.size()) - 1;

          transitionTime = mdTimes[static_cast<size_t>(snapshotIndex)];
          transitionStructure =
              *mdSnapshots[static_cast<size_t>(snapshotIndex)];
        } else {
          transitionStructure = *trajectory;
          transitionTime = simulationTime;
        }
        QUILL_LOG_DEBUG(log, "[ParallelReplica] Transition time: {:.3e} s",
                        transitionTime * params.constants().timeUnit * 1e-15);
        // A run that keeps going after the transition keeps the bond boost
        // it started with; plain copy assignment would clear it.
        trajectory->assignKeepingBias(transitionStructure);
        // A false stop_after_transition keeps the remaining dynamics steps.
        if (params.parallel_replica_options().auto_stop) {
          break;
        }

      } else if (step + 1 == params.dynamics_options().steps &&
                 transitionTime == 0) {
        // Fake refinement to prevent force-call bias (only if snapshots exist)
        if (params.parallel_replica_options().refine_transition &&
            !mdSnapshots.empty()) {
          QUILL_LOG_DEBUG(
              log,
              "[ParallelReplica] Simulation ended without seeing a transition");
          QUILL_LOG_DEBUG(
              log, "[ParallelReplica] Refining anyways to prevent bias...");
          {
            eonc::ForceCallTimer timer(refineForceCalls, *pot);
            refineTransition(mdSnapshots, true);
          }
        }
        transitionStructure = *trajectory;
      }

      mdSnapshots.clear();
      mdTimes.clear();
    }
  }

  std::unique_ptr<Matter> product;
  if (transitionTime != 0) {
    int decorrelationSteps = static_cast<int>(
        std::floor(params.parallel_replica_options().corr_time /
                       params.dynamics_options().time_step +
                   0.5));
    if (decorrelationSteps < 0) {
      decorrelationSteps = 0;
    }
    QUILL_LOG_DEBUG(log, "[ParallelReplica] Decorrelating: {} steps",
                    decorrelationSteps);
    {
      const size_t decorAt = pot->forceCallCounter.load();
      for (int dstep = 1; dstep <= decorrelationSteps; ++dstep) {
        dynamics.oneStep(dstep);
      }
      profile_.decorrelate = pot->forceCallCounter.load() - decorAt;
    }
    product = std::make_unique<Matter>(pot, params);
    *product = *trajectory;
    {
      const size_t productAt = pot->forceCallCounter.load();
      product->relax();
      if (!eonc::io::io_ok(product->matter2con("product.con"))) {
        QUILL_LOG_ERROR(log, "Failed to write product.con");
      }
      profile_.product = pot->forceCallCounter.load() - productAt;
    }
  }

  // Write results
  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  {
    std::ofstream out(resultsFilename, std::ios::binary);
    if (out) {
      out << std::format(
          "{} potential_type\n",
          magic_enum::enum_name<PotType>(params.potential_options().potential));
      out << std::format("{} random_seed\n", params.main_options().randomSeed);
      out << std::format("{:f} potential_energy_reactant\n",
                         reactant->getPotentialEnergy());
      out << std::format("{} force_calls_refine\n", refineForceCalls);
      const size_t totalCalls = pot->forceCallCounter.load();
      out << std::format("{} total_force_calls\n", totalCalls);
      out << std::format("{} force_calls_initial_minimize\n", profile_.initial);
      out << std::format("{} force_calls_dephase_md\n", profile_.dephaseMd);
      out << std::format("{} force_calls_dephase_quench\n",
                         profile_.dephaseQuench);
      out << std::format("{} force_calls_md\n", profile_.md);
      out << std::format("{} force_calls_state_check\n", profile_.stateCheck);
      out << std::format("{} force_calls_decorrelate\n", profile_.decorrelate);
      out << std::format("{} force_calls_product_minimize\n", profile_.product);
      const size_t accounted = profile_.initial + profile_.dephaseMd +
                               profile_.dephaseQuench + profile_.md +
                               profile_.stateCheck + refineForceCalls +
                               profile_.decorrelate + profile_.product;
      const long gap =
          static_cast<long>(totalCalls) - static_cast<long>(accounted);
      out << std::format("{} force_calls_other\n", gap);
      out << std::format("{} state_checks\n", profile_.stateChecks);
      out << std::format("{} state_checks_hot_inside\n", profile_.hotInside);
      out << std::format("{} basin_early_accepts\n", profile_.earlyAccepts);

      if (transitionTime == 0) {
        out << "0 transition_found\n";
        out << std::format("{:e} simulation_time_s\n",
                           simulationTime * params.constants().timeUnit *
                               1.0e-15);
      } else {
        out << "1 transition_found\n";
        out << std::format("{:e} transition_time_s\n",
                           transitionTime * params.constants().timeUnit *
                               1.0e-15);
        out << std::format("{:e} correlation_time_s\n",
                           params.parallel_replica_options().corr_time *
                               params.constants().timeUnit * 1.0e-15);
        out << std::format("{:f} potential_energy_product\n",
                           product->getPotentialEnergy());
      }
      out << std::format("{:f} speedup\n",
                         simulationTime /
                             (params.dynamics_options().steps *
                              params.dynamics_options().time_step));
    }
  }

  // bondBoost is a stack local of this function and setBiasPotential does not
  // own it, so the trajectory must not carry the pointer past the return.
  // getBiasForces reads a null bias potential as zero bias, which is what a
  // Matter that is no longer being boosted means.
  trajectory->setBiasPotential(nullptr);

  return trajectory;
}

void ParallelReplicaJob::dephase(Matter &trajectory, BondBoost *bias) {
  Dynamics dynamics(&trajectory, params);
  // The boost object names this Matter. Attach it here, and reset the
  // trajectory with assignKeepingBias so each dephase loop keeps it.
  trajectory.setBiasPotential(bias);

  const double dt = params.dynamics_options().time_step;
  if (!(dt > 0.0)) {
    throw std::invalid_argument(
        "ParallelReplicaJob::dephase: time_step must be positive");
  }
  int dephaseSteps = static_cast<int>(
      std::floor(params.parallel_replica_options().dephase_time / dt + 0.5));
  if (dephaseSteps < 1)
    dephaseSteps = 1;
  const long maxLoops =
      std::max(1L, params.parallel_replica_options().dephase_loop_max);
  QUILL_LOG_DEBUG(log, "[ParallelReplica] Dephasing: {} steps (max {} loops)",
                  dephaseSteps, maxLoops);

  Matter initial(pot, params);
  initial = trajectory;

  for (long loop = 0; loop < maxLoops; ++loop) {
    trajectory.assignKeepingBias(initial);
    dynamics.setThermalVelocity();

    {
      const size_t mdAt = pot->forceCallCounter.load();
      for (int step = 1; step <= dephaseSteps; step++) {
        dynamics.oneStep(step);
      }
      profile_.dephaseMd += pot->forceCallCounter.load() - mdAt;
    }

    Matter minimized(pot, params);
    minimized = trajectory;
    {
      const size_t quenchAt = pot->forceCallCounter.load();
      basinRelax(minimized, false);
      profile_.dephaseQuench += pot->forceCallCounter.load() - quenchAt;
    }

    if (minimized.compare(*reactant)) {
      QUILL_LOG_DEBUG(log, "[ParallelReplica] Dephasing successful");
      return;
    }
    QUILL_LOG_DEBUG(
        log,
        "[ParallelReplica] Transition occured during dephasing; Restarting");
  }
  // Exhausted retries: keep last dephased trajectory rather than hang.
  QUILL_LOG_DEBUG(log,
                  "[ParallelReplica] Dephase loop max reached; continuing");
}

int ParallelReplicaJob::refineTransition(
    const std::vector<std::shared_ptr<Matter>> &snapshots, bool fake) {
  int lo = 0;
  int hi = static_cast<int>(snapshots.size()) - 1;

  while ((hi - lo) > 1) {
    int mid = lo + (hi - lo) / 2;
    Matter snapshot(*snapshots[mid]);
    basinRelax(snapshot, true);

    bool stillReactant;
    if (!fake) {
      stillReactant = snapshot.compare(*reactant);
    } else {
      stillReactant = static_cast<bool>(eonc::rng::randomInt(0, 1));
    }

    if (stillReactant) {
      lo = mid;
    } else {
      hi = mid;
    }
  }

  return (lo + hi) / 2 + 1;
}

void ParallelReplicaJob::basinRelax(Matter &matter, bool quiet) {
  if (params.optimizer_options().method != OptType::CG) {
    matter.relax(quiet);
    return;
  }
  // The finite-difference curvature probe is one force call per CG step.
  // Basin identity does not need that probe once a previous step has a
  // positive curvature, and it does not need the norm driven to the usual
  // tolerance once the trial is inside the reactant ball with a small force
  // and a step that fits in the remaining margin.
  BasinQuenchSecant secant;
  const double distTol =
      params.structure_comparison_options().distance_difference;
  const double forceGate = params.optimizer_options().converged_force;
  auto accept = [&](Matter &trial, double stepSize) -> bool {
    if (!(distTol > 0.0)) {
      return false;
    }
    const double dist = trial.perAtomNorm(*reactant);
    if (!(dist < distTol)) {
      return false;
    }
    const double margin = distTol - dist;
    if (!(stepSize <= 0.5 * margin)) {
      return false;
    }
    if (!(trial.maxForce() <= forceGate)) {
      return false;
    }
    ++profile_.earlyAccepts;
    return true;
  };
  eonc::helpers::relaxMatterAccepted(matter, params, quiet, false, false,
                                     std::string(), std::string(), nullptr,
                                     accept);
}

} // namespace eonc
