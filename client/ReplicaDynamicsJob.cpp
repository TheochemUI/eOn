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
#include "eon/ReplicaDynamicsJob.h"
#include "eon/BaseStructures.h"
#include "eon/Dynamics.h"
#include "eon/ForceCallTimer.h"
#include "eon/HelperFunctions.h"
#include <stdexcept>

#include <format>
#include <fstream>

namespace eonc {

std::vector<std::string> ReplicaDynamicsJob::run() {
  auto seed = std::make_shared<Matter>(pot, params);
  std::string reactantFilename =
      eonc::helpers::getRelevantFile(params.main_options().conFilename);
  if (!eonc::io::io_ok(seed->con2matter(reactantFilename))) {
    QUILL_LOG_CRITICAL(log, "Failed to load {}", reactantFilename);
    throw std::runtime_error("failed to load " + reactantFilename);
  }
  (void)runFromMatter(std::move(seed));
  return returnFiles;
}

std::shared_ptr<Matter>
ReplicaDynamicsJob::runFromMatter(std::shared_ptr<Matter> initial) {
  if (!initial) {
    throw std::runtime_error("runFromMatter: initial Matter is null");
  }
  current = initial;
  current->setPotential(pot);
  reactant = std::make_shared<Matter>(pot, params);
  saddle = std::make_shared<Matter>(pot, params);
  product = std::make_shared<Matter>(pot, params);
  finalState = std::make_shared<Matter>(pot, params);
  finalStateTmp = std::make_shared<Matter>(pot, params);

  minimizeFCalls = mdFCalls = refineFCalls = dephaseFCalls = 0;
  time = 0.0;

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

  return current;
}

bool ReplicaDynamicsJob::checkState(Matter *curr, Matter *react) {
  Matter tmp(pot, params);
  tmp = *curr;
  tmp.relax(true);
  return !tmp.compare(*react);
}

ReplicaDynamicsJob::PrdClock ReplicaDynamicsJob::prdClock() const {
  const double dt = params.dynamics_options().time_step;
  auto to_steps = [&](double interval) -> long {
    if (!(dt > 0.0) || !(interval > 0.0)) {
      return 1;
    }
    const long n = static_cast<long>(interval / dt);
    return n < 1 ? 1 : n;
  };
  PrdClock c;
  c.state_check =
      to_steps(params.parallel_replica_options().state_check_interval);
  c.record = to_steps(params.parallel_replica_options().record_interval);
  if (c.record > c.state_check) {
    c.record = c.state_check;
  }
  c.buffer = c.state_check / c.record;
  if (c.buffer < 1) {
    c.buffer = 1;
  }
  return c;
}

long ReplicaDynamicsJob::refine(
    const std::vector<std::shared_ptr<Matter>> &buff, Matter *react) {
  QUILL_LOG_TRACE_L1(log, "Refining transition time.");
  const long n = static_cast<long>(buff.size());
  if (n <= 1) {
    throw std::runtime_error(
        "ReplicaDynamics refine: need at least two snapshots");
  }

  long lo = 0;
  long hi = n - 1;

  while ((hi - lo) > 1) {
    long mid = lo + (hi - lo) / 2;
    if (!checkState(buff[static_cast<size_t>(mid)].get(), react)) {
      lo = mid;
    } else {
      hi = mid;
    }
  }

  long idx = (lo + hi) / 2 + 1;
  if (idx < 1) {
    idx = 1;
  }
  if (idx >= n) {
    idx = n - 1;
  }
  return idx;
}

void ReplicaDynamicsJob::dephase() {
  const double dt = params.dynamics_options().time_step;
  if (!(dt > 0.0)) {
    throw std::invalid_argument(
        "ReplicaDynamicsJob::dephase: time_step must be positive");
  }
  long DephaseSteps =
      static_cast<long>(params.parallel_replica_options().dephase_time / dt);
  Dynamics dephaseDynamics(current.get(), params);
  QUILL_LOG_DEBUG(log, "Dephasing for {:.2f} fs",
                  params.parallel_replica_options().dephase_time *
                      params.constants().timeUnit);

  long step = 0, loop = 0;

  while (step < DephaseSteps) {
    long dephaseBufferLength = DephaseSteps - step;
    if (dephaseBufferLength < 1) {
      break;
    }
    loop++;
    std::vector<std::shared_ptr<Matter>> dephaseBuffer(dephaseBufferLength);

    for (long i = 0; i < dephaseBufferLength; i++) {
      dephaseBuffer[i] = std::make_shared<Matter>(pot, params);
      dephaseDynamics.oneStep();
      *dephaseBuffer[i] = *current;
    }

    bool transitionFlag = checkState(current.get(), reactant.get());

    if (transitionFlag) {
      if (dephaseBuffer.size() < 2) {
        AtomMatrix velocity = current->getVelocities();
        velocity = velocity * (-1);
        current->setVelocities(velocity);
        continue;
      }
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

    const long loop_max = params.parallel_replica_options().dephase_loop_max;
    if (loop_max > 0 && loop >= loop_max) {
      QUILL_LOG_DEBUG(
          log,
          "Reach dephase loop maximum, stop dephasing! Dephased for {} steps",
          step);
      break;
    }
    QUILL_LOG_DEBUG(log, "Successfully Dephased for {:.2f} fs",
                    step * params.dynamics_options().time_step *
                        params.constants().timeUnit);
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
          magic_enum::enum_name<PotType>(params.potential_options().potential));
      out << std::format("{} random_seed\n", params.main_options().randomSeed);
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
                               params.constants().timeUnit);
        out << std::format("{:f} potential_energy_product\n",
                           product->getPotentialEnergy());
        out << std::format("{:f} moved_distance\n",
                           product->distanceTo(*reactant));
      }

      out << std::format("{:e} simulation_time_s\n",
                         time * 1.0e-15 * params.constants().timeUnit);
      out << std::format("{:f} speedup\n",
                         time / params.dynamics_options().steps /
                             params.dynamics_options().time_step);
    }
  }

  std::string reactantFilename("reactant.con");
  returnFiles.push_back(reactantFilename);
  if (!eonc::io::io_ok(reactant->matter2con(reactantFilename))) {
    QUILL_LOG_ERROR(log, "Failed to write {}", reactantFilename);
  }

  if (newStateFlag) {
    std::string productFilename("product.con");
    returnFiles.push_back(productFilename);
    if (!eonc::io::io_ok(product->matter2con(productFilename))) {
      QUILL_LOG_ERROR(log, "Failed to write {}", productFilename);
    }

    if (params.parallel_replica_options().refine_transition) {
      std::string saddleFilename("saddle.con");
      returnFiles.push_back(saddleFilename);
      if (!eonc::io::io_ok(saddle->matter2con(saddleFilename))) {
        QUILL_LOG_ERROR(log, "Failed to write {}", saddleFilename);
      }
    }
  }
}

} // namespace eonc
