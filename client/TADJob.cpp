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
#include "TADJob.h"
#include "Dynamics.h"
#include "ForceCallTimer.h"

#include <cmath>

void TADJob::initExtra() {
  crossing = std::make_shared<Matter>(pot, params);

  QUILL_LOG_DEBUG(log, "Temperature Accelerated Dynamics, running");
  QUILL_LOG_DEBUG(log,
                  "High temperature MD simulation running at {:.2f} K to "
                  "simulate dynamics at {:.2f} K",
                  params.main_options.temperature,
                  params.tad_options.low_temperature);
}

void TADJob::reportResults() {
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

int TADJob::dynamics() {
  bool transitionFlag = false, recordFlag = true, stopFlag = false,
       firstTransitFlag = false;
  long nFreeCoord = reactant->numberOfFreeAtoms() * 3;
  long mdBufferLength;
  long step = 0, refineStep, newStateStep = 0;
  long nCheck = 0, nRecord = 0, nState = 0;
  long StateCheckInterval, RecordInterval;
  double kinE, kinT, avgT, varT;
  double kB = params.constants.kB;
  double correctedTime = 0.0;
  double stopTime = 0.0, sumSimulatedTime = 0.0;
  double Temp = 0.0, sumT = 0.0, sumT2 = 0.0;
  double correctionFactor = 1.0;
  double transitionTime_current = 0.0, transitionTime_previous = 0.0;
  double delta, minmu, factor, highT, lowT;

  AtomMatrix velocity;

  minCorrectedTime = 1.0e200;
  lowT = params.tad_options.low_temperature;
  highT = params.main_options.temperature;
  delta = params.tad_options.confidence;
  minmu = params.tad_options.min_prefactor;
  factor = std::log(1.0 / delta) / minmu;
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

  Dynamics TAD(current.get(), DynamicsConfig::fromParams(params));
  TAD.setThermalVelocity();

  {
    eonc::ForceCallTimer timer(dephaseFCalls);
    dephase();
  }

  QUILL_LOG_DEBUG(
      log,
      "Starting MD run\nTemperature: {:.2f} Kelvin"
      "Total Simulation Time: {:.2f} fs\nTime Step: {:.2f} fs\nTotal Steps: "
      "{}\n",
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
    kinE = current->getKineticEnergy();
    kinT = (2.0 * kinE / nFreeCoord / kB);
    sumT += kinT;
    sumT2 += kinT * kinT;
    QUILL_LOG_TRACE_L1(log, "steps = {:10d} temp = {:10.5f} ", step, kinT);

    TAD.oneStep();
    mdFCalls++;

    time += params.dynamics_options.time_step;
    nCheck++;
    step++;
    QUILL_LOG_TRACE_L1(log, "step = {:4d}, time= {:10.4f}", step, time);

    if (params.parallel_replica_options.refine_transition && recordFlag &&
        !newStateFlag) {
      if (nCheck % RecordInterval == 0) {
        *mdBuffer[nRecord] = *current;
        timeBuffer[nRecord] = time;
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
      *crossing = *mdBuffer[refineStep];
      transitionTime = transitionTime_current - transitionTime_previous;
      transitionTime_previous = transitionTime_current;
      barrier = crossing->getPotentialEnergy() - reactant->getPotentialEnergy();
      QUILL_LOG_DEBUG(log, "barrier= {:.3f}", barrier);
      correctionFactor = std::exp(barrier / kB * (1.0 / lowT - 1.0 / highT));
      correctedTime = transitionTime * correctionFactor;
      sumSimulatedTime += transitionTime;

      *current = *mdBuffer[refineStep - 1];
      velocity = current->getVelocities();
      velocity = velocity * (-1);
      current->setVelocities(velocity);

      if (correctedTime < minCorrectedTime) {
        minCorrectedTime = correctedTime;
        *saddle = *crossing;
        *finalState = *finalStateTmp;
      }
      stopTime = factor * std::pow(minCorrectedTime / factor, lowT / highT);
      QUILL_LOG_DEBUG(
          log,
          "tranisitonTime= {:.3e} s, Barrier= {:.3f} eV, correctedTime= {:.3e} "
          "s, "
          "SimulatedTime= {:.3e} s, minCorTime= {:.3e} s, stopTime= {:.3e} s",
          transitionTime * 1e-15 * params.constants.timeUnit, barrier,
          correctedTime * 1e-15 * params.constants.timeUnit,
          sumSimulatedTime * 1e-15 * params.constants.timeUnit,
          minCorrectedTime * 1.0e-15 * params.constants.timeUnit,
          stopTime * 1.0e-15 * params.constants.timeUnit);

      transitionFlag = false;
    }

    if (firstTransitFlag && sumSimulatedTime >= stopTime) {
      stopFlag = true;
      newStateFlag = true;
    }

    if ((step % tenthSteps == 0) || (step == params.dynamics_options.steps)) {
      double maxAtomDistance = current->perAtomNorm(*reactant);
      QUILL_LOG_DEBUG(
          log, "progress: {:.0f}%, max displacement: {:.3f}, step {}/{}",
          static_cast<double>(100.0 * step) / params.dynamics_options.steps,
          maxAtomDistance, step, params.dynamics_options.steps);
    }

    if (step == params.dynamics_options.steps) {
      stopFlag = true;
      if (firstTransitFlag) {
        QUILL_LOG_DEBUG(log, "Detected one transition");
      } else {
        QUILL_LOG_DEBUG(log, "Failed to detect any transition");
      }
    }
  }

  avgT = sumT / step;
  varT = sumT2 / step - avgT * avgT;

  QUILL_LOG_DEBUG(log,
                  "Temperature : Average = {} ; Stddev = {} ; Factor = {}; "
                  "Average_Boost = {}",
                  avgT, std::sqrt(varT), varT / avgT / avgT * nFreeCoord / 2,
                  minCorrectedTime / step / params.dynamics_options.time_step);
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
