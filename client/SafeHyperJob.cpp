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
#include "SafeHyperJob.h"
#include "BondBoost.h"
#include "Dynamics.h"
#include "ForceCallTimer.h"

#include <cmath>

void SafeHyperJob::reportResults() {
  if (newStateFlag) {
    QUILL_LOG_DEBUG(log, "Transition time: {:.2e} s",
                    minCorrectedTime * 1.0e-15 * params.constants.timeUnit);
  } else {
    QUILL_LOG_DEBUG(log,
                    "No new state was found in {} dynamics steps ({:.3e} s)",
                    params.dynamics_options.steps,
                    time * 1.0e-15 * params.constants.timeUnit);
  }
}

int SafeHyperJob::dynamics() {
  bool transitionFlag = false, recordFlag = true, stopFlag = false,
       firstTransitFlag = false;
  long nFreeCoord = reactant->numberOfFreeAtoms() * 3;
  long mdBufferLength;
  long step = 0, refineStep, newStateStep = 0;
  long nCheck = 0, nRecord = 0, nBoost = 0, nState = 0;
  long StateCheckInterval, RecordInterval;
  double kinE, kinT, avgT, varT;
  double kB = params.constants.kB;
  double correctedTime = 0.0, sumCorrectedTime = 0.0, firstTransitionTime = 0.0;
  double Temp = 0.0, sumT = 0.0, sumT2 = 0.0;
  double sumboost = 0.0, boost = 1.0, boostPotential = 0.0;
  double transitionTime_current = 0.0, transitionTime_pre = 0.0;
  AtomMatrix velocity;

  minCorrectedTime = 1.0e200;
  StateCheckInterval =
      static_cast<long>(params.parallel_replica_options.state_check_interval /
                        params.dynamics_options.time_step);
  RecordInterval =
      static_cast<long>(params.parallel_replica_options.record_interval /
                        params.dynamics_options.time_step);
  Temp = params.main_options.temperature;
  newStateFlag = metaStateFlag = false;

  mdBufferLength = static_cast<long>(StateCheckInterval / RecordInterval);
  std::vector<std::shared_ptr<Matter>> mdBuffer(mdBufferLength);
  for (long i = 0; i < mdBufferLength; i++) {
    mdBuffer[i] = std::make_shared<Matter>(pot, params);
  }
  timeBuffer.resize(mdBufferLength);
  biasBuffer.resize(mdBufferLength);

  Dynamics safeHyper(current.get(), DynamicsConfig::fromParams(params));
  BondBoost bondBoost(current.get(), params);

  if (params.hyperdynamics_options.bias_potential ==
      Hyperdynamics::BOND_BOOST) {
    bondBoost.initialize();
  }

  safeHyper.setThermalVelocity();

  {
    eonc::ForceCallTimer timer(dephaseFCalls);
    dephase();
  }

  QUILL_LOG_DEBUG(
      log,
      "Starting MD run\nTemperature: {:.2f} Kelvin\n"
      "Total Simulation Time: {:.2f} fs\nTime Step: {:.2f} fs\nTotal Steps: {}",
      Temp,
      params.dynamics_options.steps * params.dynamics_options.time_step *
          params.constants.timeUnit,
      params.dynamics_options.time_step * params.constants.timeUnit,
      params.dynamics_options.steps);
  QUILL_LOG_DEBUG(log, "MD buffer length: {}", mdBufferLength);

  long tenthSteps = params.dynamics_options.steps / 10;
  if (tenthSteps == 0) {
    tenthSteps = params.dynamics_options.steps;
  }

  while (!stopFlag) {
    if ((params.hyperdynamics_options.bias_potential ==
         Hyperdynamics::BOND_BOOST) &&
        !newStateFlag) {
      boostPotential = bondBoost.boost();
      QUILL_LOG_TRACE_L1(log, "step= {} , boost = {:.5f}", step,
                         boostPotential);
      boost = std::exp(boostPotential / kB / Temp);
      time += params.dynamics_options.time_step * boost;
      if (boost > 1.0) {
        sumboost += boost;
        nBoost++;
      }
    }

    kinE = current->getKineticEnergy();
    kinT = (2.0 * kinE / nFreeCoord / kB);
    sumT += kinT;
    sumT2 += kinT * kinT;
    QUILL_LOG_TRACE_L1(log, "steps = {:10} temp = {:10.5f}", step, kinT);

    safeHyper.oneStep();
    mdFCalls++;

    nCheck++;
    step++;
    QUILL_LOG_TRACE_L1(log, "step = {:4}, time = {:10.4f}", step, time);

    if (params.parallel_replica_options.refine_transition && recordFlag &&
        !newStateFlag) {
      if (nCheck % RecordInterval == 0) {
        *mdBuffer[nRecord] = *current;
        timeBuffer[nRecord] = time;
        biasBuffer[nRecord] = boostPotential;
        nRecord++;
      }
    }

    if ((nCheck == StateCheckInterval) && !newStateFlag) {
      nCheck = 0;
      nRecord = 0;
      {
        eonc::ForceCallTimer timer(minimizeFCalls);
        transitionFlag = checkState(current.get(), reactant.get());
      }
      if (transitionFlag) {
        nState++;
        QUILL_LOG_DEBUG(log, "New State {}: ", nState);
        *finalStateTmp = *current;
        transitionTime = time;
        newStateStep = step;
        transitionStep = newStateStep;
        firstTransitFlag = 1;
      }
    }

    if (transitionFlag) {
      QUILL_LOG_TRACE_L1(log, "Refining transition time.");
      {
        eonc::ForceCallTimer timer(refineFCalls);
        refineStep = refine(mdBuffer, reactant.get());
      }

      transitionStep =
          newStateStep - StateCheckInterval + refineStep * RecordInterval;
      transitionTime_current = timeBuffer[refineStep];
      transitionTime = transitionTime_current - transitionTime_pre;
      transitionTime_pre = transitionTime_current;
      transitionPot = biasBuffer[refineStep];
      correctedTime =
          transitionTime * std::exp((-1) * transitionPot / kB / Temp);
      sumCorrectedTime += correctedTime;
      if (nState == 1) {
        firstTransitionTime = transitionTime;
      }

      *current = *mdBuffer[refineStep - 1];
      velocity = current->getVelocities();
      velocity = velocity * (-1);
      current->setVelocities(velocity);

      if (correctedTime < minCorrectedTime) {
        minCorrectedTime = correctedTime;
        *saddle = *mdBuffer[refineStep];
        *finalState = *finalStateTmp;
      }
      QUILL_LOG_DEBUG(log,
                      "tranisitonTime= {:.3e} s, biasPot= {:.3f} eV, "
                      "correctedTime= {:.3e} s, "
                      "sumCorrectedTime= {:.3e} s, minCorTime= {:.3e} s",
                      transitionTime * 1e-15 * params.constants.timeUnit,
                      transitionPot,
                      correctedTime * 1e-15 * params.constants.timeUnit,
                      sumCorrectedTime * 1e-15 * params.constants.timeUnit,
                      minCorrectedTime * 1.0e-15 * params.constants.timeUnit);

      transitionFlag = false;
    }

    if (firstTransitFlag && sumCorrectedTime > firstTransitionTime) {
      stopFlag = true;
      newStateFlag = true;
    }

    if ((step % tenthSteps == 0) || (step == params.dynamics_options.steps)) {
      double maxAtomDistance = current->perAtomNorm(*reactant);
      QUILL_LOG_DEBUG(
          log, "progress: {:.0f}%, max displacement: {:.3f}, step {} / {}",
          static_cast<double>(100.0 * step / params.dynamics_options.steps),
          maxAtomDistance, step, params.dynamics_options.steps);
    }
  }

  avgT = sumT / step;
  varT = sumT2 / step - avgT * avgT;

  if (nBoost > 0) {
    QUILL_LOG_DEBUG(log,
                    "Temperature : Average = {:.6f} ; Stddev = {:.6f} ; "
                    "Factor = {:.6f}; Boost = {:.6f}",
                    avgT, std::sqrt(varT), varT / avgT / avgT * nFreeCoord / 2,
                    sumboost / nBoost);
  } else {
    QUILL_LOG_DEBUG(
        log,
        "Temperature : Average = {:.6f} ; Stddev = {:.6f} ; Factor = {:.6f}",
        avgT, std::sqrt(varT), varT / avgT / avgT * nFreeCoord / 2);
  }
  if (std::isfinite(avgT) == 0) {
    QUILL_LOG_DEBUG(log, "Infinite average temperature, something went wrong!");
    newStateFlag = false;
  }

  *product = *finalState;

  if (newStateFlag) {
    return 1;
  } else {
    return 0;
  }
}
