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
#include "ReplicaExchangeJob.h"
#include "BaseStructures.h"
#include "Dynamics.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "SafeMath.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <string>
#include <thread>

std::vector<std::string> ReplicaExchangeJob::run() {
  long samplingSteps =
      static_cast<long>(params.replica_exchange_options.sampling_time /
                            params.dynamics_options.time_step +
                        0.5);
  long exchangePeriodSteps =
      static_cast<long>(params.replica_exchange_options.exchange_period /
                            params.dynamics_options.time_step +
                        0.5);
  const double kB = params.constants.kB;

  std::string posFilename =
      eonc::helpers::getRelevantFile(params.main_options.conFilename);
  pos = std::make_shared<Matter>(pot, params);
  pos->con2matter(posFilename);

  QUILL_LOG_DEBUG(log, "Running Replica Exchange");

  long refForceCalls = PotRegistry::get().total_force_calls();

  const long nReplicas = params.replica_exchange_options.replicas;
  std::vector<std::shared_ptr<Matter>> replica(nReplicas);
  std::vector<std::unique_ptr<Dynamics>> replicaDynamics(nReplicas);

  // Per-image potentials: each replica gets its own potential instance when
  // the potential requires it (e.g. ML potentials with internal state)
  const bool perImage = pot->needsPerImageInstance();
  for (long i = 0; i < nReplicas; i++) {
    auto replicaPot = perImage ? eonc::helpers::makePotential(params) : pot;
    replica[i] = std::make_shared<Matter>(replicaPot, params);
    *replica[i] = *pos;
    replicaDynamics[i] = std::make_unique<Dynamics>(replica[i].get(), params);
  }

  std::vector<double> replicaTemperature(nReplicas);

  QUILL_LOG_DEBUG(log, "Temperature distribution:");
  if (params.replica_exchange_options.temperature_distribution == "linear") {
    for (long i = 0; i < nReplicas; i++) {
      replicaTemperature[i] =
          params.replica_exchange_options.temperature_low +
          static_cast<double>(i) / static_cast<double>(nReplicas - 1) *
              (params.replica_exchange_options.temperature_high -
               params.replica_exchange_options.temperature_low);
      replicaDynamics[i]->setTemperature(replicaTemperature[i]);
    }
  } else if (params.replica_exchange_options.temperature_distribution ==
             "exponential") {
    double kTemp = std::log(params.replica_exchange_options.temperature_high /
                            params.replica_exchange_options.temperature_low) /
                   static_cast<double>(nReplicas - 1);
    for (long i = 0; i < nReplicas; i++) {
      replicaTemperature[i] = params.replica_exchange_options.temperature_low *
                              std::exp(kTemp * static_cast<double>(i));
      replicaDynamics[i]->setTemperature(replicaTemperature[i]);
      QUILL_LOG_DEBUG(log, "replica: {} temperature {:.0f}", i + 1,
                      replicaTemperature[i]);
    }
  }

  QUILL_LOG_DEBUG(
      log, "Replica Exchange sampling for {:.0f} fs; {} steps; {} replicas.",
      params.replica_exchange_options.sampling_time * 10.18, samplingSteps,
      params.replica_exchange_options.replicas);

  // Parallel replica dynamics when enabled and potential supports it
  const bool canParallel =
      params.main_options.parallel && (pot->isThreadSafe() || perImage);

  for (long step = 1; step <= samplingSteps; step++) {
    if (canParallel && nReplicas > 1) {
      // Parallel: each replica runs its MD step in a separate thread
      std::vector<std::thread> threads;
      threads.reserve(nReplicas);
      for (long i = 0; i < nReplicas; i++) {
        threads.emplace_back([&, i] { replicaDynamics[i]->oneStep(); });
      }
      for (auto &t : threads) {
        t.join();
      }
    } else {
      for (long i = 0; i < nReplicas; i++) {
        replicaDynamics[i]->oneStep();
      }
    }

    // Metropolis replica exchange
    if ((step % exchangePeriodSteps) == 0) {
      for (long trial = 0;
           trial < params.replica_exchange_options.exchange_trials; trial++) {
        long i = eonc::helpers::randomInt(0, nReplicas - 2);
        double energyLow = replica[i]->getPotentialEnergy();
        double energyHigh = replica[i + 1]->getPotentialEnergy();
        double kbTLow = kB * replicaTemperature[i];
        double kbTHigh = kB * replicaTemperature[i + 1];
        double pAcc =
            std::min(1.0, std::exp((energyHigh - energyLow) *
                                   (eonc::safemath::safe_recip(kbTHigh, 0.0) -
                                    eonc::safemath::safe_recip(kbTLow, 0.0))));
        double rnd = eonc::helpers::randomDouble();
        QUILL_LOG_INFO(log,
                       "step: {} trial swap, i {}, elow: {:.5f}, ehigh: "
                       "{:.5f}, pAcc: {:.5f}, rand: {}",
                       step, i, energyLow, energyHigh, pAcc, rnd);
        if (rnd < pAcc) {
          QUILL_LOG_INFO(log, "swap");
          std::swap(replica[i], replica[i + 1]);
          replicaDynamics[i]->setThermalVelocity();
          replicaDynamics[i + 1]->setThermalVelocity();
        } else {
          QUILL_LOG_INFO(log, "no swap");
        }
      }
    }
  }

  forceCalls = PotRegistry::get().total_force_calls() - refForceCalls;
  saveData();

  return returnFiles;
}

void ReplicaExchangeJob::saveData() {
  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);

  std::ofstream out(resultsFilename, std::ios::binary);
  if (out) {
    out << std::format("{} random_seed\n", params.main_options.randomSeed);
    out << std::format(
        "{} potential_type\n",
        magic_enum::enum_name<PotType>(params.potential_options.potential));
    out << std::format("{} force_calls_sampling\n", forceCalls);
  }

  std::string posFilename("pos_out.con");
  returnFiles.push_back(posFilename);
  // Note: original code opened the file but wrote nothing to it
  std::ofstream posOut(posFilename, std::ios::binary);
}
