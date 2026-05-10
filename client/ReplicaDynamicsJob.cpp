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
#include "ReplicaDynamicsJob.h"
#include "BaseStructures.h"
#include "Dynamics.h"
#include "ForceCallTimer.h"
#include "HelperFunctions.h"

#include <format>
#include <fstream>

std::vector<std::string> ReplicaDynamicsJob::run() {
  current = std::make_shared<Matter>(pot, params);
  reactant = std::make_shared<Matter>(pot, params);
  saddle = std::make_shared<Matter>(pot, params);
  product = std::make_shared<Matter>(pot, params);
  finalState = std::make_shared<Matter>(pot, params);
  finalStateTmp = std::make_shared<Matter>(pot, params);

  minimizeFCalls = mdFCalls = refineFCalls = dephaseFCalls = 0;
  time = 0.0;

  std::string reactantFilename =
      eonc::helpers::getRelevantFile(params.main_options.conFilename);
  current->con2matter(reactantFilename);

  QUILL_LOG_DEBUG(log, "Minimizing initial reactant");
  {
    eonc::ForceCallTimer timer(minimizeFCalls);
    *reactant = *current;
    reactant->relax();
  }

  initExtra();

  int status = dynamics();

  saveData(status);
  reportResults();

  return returnFiles;
}

bool ReplicaDynamicsJob::checkState(Matter *curr, Matter *react) {
  Matter tmp(pot, params);
  tmp = *curr;
  tmp.relax(true);
  return !tmp.compare(*react);
}

long ReplicaDynamicsJob::refine(
    const std::vector<std::shared_ptr<Matter>> &buff, Matter *react) {
  QUILL_LOG_TRACE_L1(log, "Refining transition time.");

  long lo = 0;
  long hi = static_cast<long>(buff.size()) - 1;

  while ((hi - lo) > 1) {
    long mid = lo + (hi - lo) / 2;
    if (!checkState(buff[mid].get(), react)) {
      lo = mid;
    } else {
      hi = mid;
    }
  }

  return (lo + hi) / 2 + 1;
}

void ReplicaDynamicsJob::dephase() {
  long DephaseSteps =
      static_cast<long>(params.parallel_replica_options.dephase_time /
                        params.dynamics_options.time_step);
  Dynamics dephaseDynamics(current.get(), DynamicsConfig::fromParams(params));
  QUILL_LOG_DEBUG(log, "Dephasing for {:.2f} fs",
                  params.parallel_replica_options.dephase_time *
                      params.constants.timeUnit);

  long step = 0, loop = 0;

  while (step < DephaseSteps) {
    long dephaseBufferLength = DephaseSteps - step;
    loop++;
    std::vector<std::shared_ptr<Matter>> dephaseBuffer(dephaseBufferLength);

    for (long i = 0; i < dephaseBufferLength; i++) {
      dephaseBuffer[i] = std::make_shared<Matter>(pot, params);
      dephaseDynamics.oneStep();
      *dephaseBuffer[i] = *current;
    }

    bool transitionFlag = checkState(current.get(), reactant.get());

    if (transitionFlag) {
      long dephaseRefineStep = refine(dephaseBuffer, reactant.get());
      QUILL_LOG_DEBUG(log, "loop = {}; dephase refine step = {}", loop,
                      dephaseRefineStep);
      long ts = dephaseRefineStep - 1;
      ts = (ts > 0) ? ts : 0;
      QUILL_LOG_DEBUG(
          log,
          "Dephasing warning: in a new state, inverse the momentum and restart "
          "from step {}",
          step + ts);
      *current = *dephaseBuffer[ts];
      AtomMatrix velocity = current->getVelocities();
      velocity = velocity * (-1);
      current->setVelocities(velocity);
      step = step + ts;
    } else {
      step = step + dephaseBufferLength;
      QUILL_LOG_TRACE_L1(log, "Successful dephasing for {} steps", step);
    }

    if ((params.parallel_replica_options.dephase_loop_stop) &&
        (loop > params.parallel_replica_options.dephase_loop_max)) {
      QUILL_LOG_DEBUG(
          log,
          "Reach dephase loop maximum, stop dephasing! Dephased for {} steps",
          step);
      break;
    }
    QUILL_LOG_DEBUG(log, "Successfully Dephased for {:.2f} fs",
                    step * params.dynamics_options.time_step *
                        params.constants.timeUnit);
  }
}

void ReplicaDynamicsJob::saveData(int status) {
  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  size_t totalFCalls = minimizeFCalls + mdFCalls + dephaseFCalls + refineFCalls;

  {
    std::ofstream out(resultsFilename, std::ios::binary);
    if (out) {
      out << std::format(
          "{} potential_type\n",
          magic_enum::enum_name<PotType>(params.potential_options.potential));
      out << std::format("{} random_seed\n", params.main_options.randomSeed);
      out << std::format("{:f} potential_energy_reactant\n",
                         reactant->getPotentialEnergy());
      out << std::format("{} total_force_calls\n", totalFCalls);
      out << std::format("{} force_calls_dephase\n", dephaseFCalls);
      out << std::format("{} force_calls_dynamics\n", mdFCalls);
      out << std::format("{} force_calls_minimize\n", minimizeFCalls);
      out << std::format("{} force_calls_refine\n", refineFCalls);
      out << std::format("{} transition_found\n", (newStateFlag) ? 1 : 0);

      if (newStateFlag) {
        out << std::format("{:e} transition_time_s\n",
                           minCorrectedTime * 1.0e-15 *
                               params.constants.timeUnit);
        out << std::format("{:f} potential_energy_product\n",
                           product->getPotentialEnergy());
        out << std::format("{:f} moved_distance\n",
                           product->distanceTo(*reactant));
      }

      out << std::format("{:e} simulation_time_s\n",
                         time * 1.0e-15 * params.constants.timeUnit);
      out << std::format("{:f} speedup\n",
                         time / params.dynamics_options.steps /
                             params.dynamics_options.time_step);
    }
  }

  std::string reactantFilename("reactant.con");
  returnFiles.push_back(reactantFilename);
  reactant->matter2con(reactantFilename);

  if (newStateFlag) {
    std::string productFilename("product.con");
    returnFiles.push_back(productFilename);
    product->matter2con(productFilename);

    if (params.parallel_replica_options.refine_transition) {
      std::string saddleFilename("saddle.con");
      returnFiles.push_back(saddleFilename);
      saddle->matter2con(saddleFilename);
    }
  }
}
