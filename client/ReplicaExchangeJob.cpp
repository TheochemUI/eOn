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

std::vector<std::string> ReplicaExchangeJob::run(void) {
  long i, step,
      samplingSteps =
          long(params->repexc.samplingTime / params->md.timeStep + 0.5);
  long exchangePeriodSteps =
      long(params->repexc.exchangePeriod / params->md.timeStep + 0.5);
  double energyLow, energyHigh;
  double kbTLow, kbTHigh;
  double kB = params->kB;
  double pAcc;
  std::shared_ptr<Matter> tmpMatter;

  std::string posFilename =
      helper_functions::getRelevantFile(params->main.conFilename);
  pos = std::make_shared<Matter>(pot, params);
  pos->con2matter(posFilename);

  SPDLOG_LOGGER_DEBUG(log, "Running Replica Exchange");

  long refForceCalls = Potential::fcalls;

  // allocate a Matter and Dynamics object for each replica
  std::vector<std::shared_ptr<Matter>> replica;
  replica.resize(params->repexc.replicas);
  Dynamics *replicaDynamics[params->repexc.replicas];
  for (i = 0; i < params->repexc.replicas; i++) {
    replica[i] = std::make_shared<Matter>(pot, params);
    *replica[i] = *pos;
    replicaDynamics[i] = new Dynamics(replica[i].get(), params.get());
  }

  // assign temperatures
  double replicaTemperature[params->repexc.replicas];

  SPDLOG_LOGGER_DEBUG(log, "Temperature distribution:");
  if (params->repexc.temperatureDistribution == "linear") {
    for (i = 0; i < params->repexc.replicas; i++) {
      replicaTemperature[i] =
          params->repexc.temperatureLow +
          double(i) / double(params->repexc.replicas - 1.0) *
              (params->repexc.temperatureHigh - params->repexc.temperatureLow);
      replicaDynamics[i]->setTemperature(replicaTemperature[i]);
    }
  } else if (params->repexc.temperatureDistribution == "exponential") {
    double kTemp = std::log(params->repexc.temperatureHigh /
                            params->repexc.temperatureLow) /
                   (params->repexc.replicas - 1.0);
    // cout <<"kTemp: "<<kTemp<<endl;
    for (i = 0; i < params->repexc.replicas; i++) {
      replicaTemperature[i] = params->repexc.temperatureLow * exp(kTemp * i);
      replicaDynamics[i]->setTemperature(replicaTemperature[i]);
      SPDLOG_LOGGER_DEBUG(log, "replica: {} temperature {:.0f}", i + 1,
                          replicaTemperature[i]);
    }
  }

  SPDLOG_LOGGER_DEBUG(
      log, "Replica Exchange sampling for {:.0f} fs; {} steps; {} replicas.",
      params->repexc.samplingTime * 10.18, samplingSteps,
      params->repexc.replicas);

  for (step = 1; step <= samplingSteps; step++) {
    for (i = 0; i < params->repexc.replicas; i++) {
      replicaDynamics[i]->oneStep();
      SPDLOG_LOGGER_DEBUG(log, "step: {} i {} energy: {}", step, i,
                          replica[i]->getPotentialEnergy());
    }
    if ((step % exchangePeriodSteps) == 0) {
      for (long trial = 0; trial < params->repexc.exchangeTrials; trial++) {
        i = helper_functions::randomInt(0, params->repexc.replicas - 2);
        energyLow = replica[i]->getPotentialEnergy();
        energyHigh = replica[i + 1]->getPotentialEnergy();
        kbTLow = kB * replicaTemperature[i];
        kbTHigh = kB * replicaTemperature[i + 1];
        pAcc = std::min(1.0, exp((energyHigh - energyLow) *
                                 (1.0 / kbTHigh - 1.0 / kbTLow)));
        double tmp = helper_functions::randomDouble();
        SPDLOG_LOGGER_INFO(log,
                           "step: {} trial swap, i {}, elow: {:.5f}, ehigh: "
                           "{:.5f}, pAcc: {:.5f}, rand: {}",
                           step, i, energyLow, energyHigh, pAcc, tmp);
        // if(helper_functions::randomDouble()<pAcc)
        if (tmp < pAcc) {
          SPDLOG_LOGGER_INFO(log, "swap");
          // swap configurations
          tmpMatter = replica[i];
          replica[i] = replica[i + 1];
          replica[i + 1] = tmpMatter;
          // reset velocities
          replicaDynamics[i]->setThermalVelocity();
          replicaDynamics[i + 1]->setThermalVelocity();
        } else {
          SPDLOG_LOGGER_INFO(log, "no swap");
        }
      }
    }
  }

  // forceCalls = Potential::fcalls - refForceCalls;
  saveData();

  // delete Matter and Dynamics objects

  return returnFiles;
}

void ReplicaExchangeJob::saveData(void) {

  FILE *fileResults, *filePos;

  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  fileResults = fopen(resultsFilename.c_str(), "wb");

  fprintf(fileResults, "%ld random_seed\n", params->main.randomSeed);
  fprintf(fileResults, "%s potential_type\n",
          std::string{magic_enum::enum_name<PotType>(params->pot.potential)}
              .c_str());
  // fprintf(fileResults, "%ld force_calls_sampling\n", forceCalls);
  fclose(fileResults);

  std::string posFilename("pos_out.con");
  returnFiles.push_back(posFilename);
  filePos = fopen(posFilename.c_str(), "wb");
  fclose(filePos);

  return;
}
