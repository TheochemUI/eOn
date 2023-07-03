
#include <cstdlib>

#include "BaseStructures.h"
#include "Dynamics.h"
#include "HelperFunctions.h"
#include "Log.h"
#include "Matter.h"
#include "ReplicaExchangeJob.h"

std::vector<std::string> ReplicaExchangeJob::run(void) {
  long i, step,
      samplingSteps =
          long(params->repexcSamplingTime / params->mdTimeStep + 0.5);
  long exchangePeriodSteps =
      long(params->repexcExchangePeriod / params->mdTimeStep + 0.5);
  double energyLow, energyHigh;
  double kbTLow, kbTHigh;
  double kB = params->kB;
  double pAcc;
  std::shared_ptr<Matter> tmpMatter;

  string posFilename = helper_functions::getRelevantFile(params->conFilename);
  pos = std::make_shared<Matter>(pot, params);
  pos->con2matter(posFilename);

  log("\nRunning Replica Exchange\n\n");

  long refForceCalls = Potential::fcalls;

  // allocate a Matter and Dynamics object for each replica
  std::vector<std::shared_ptr<Matter>> replica;
  replica.resize(params->repexcReplicas);
  Dynamics *replicaDynamics[params->repexcReplicas];
  for (i = 0; i < params->repexcReplicas; i++) {
    replica[i] = std::make_shared<Matter>(pot, params);
    *replica[i] = *pos;
    replicaDynamics[i] = new Dynamics(replica[i].get(), params.get());
  }

  // assign temperatures
  double replicaTemperature[params->repexcReplicas];

  log("Temperature distribution:\n");
  if (params->repexcTemperatureDistribution == "linear") {
    for (i = 0; i < params->repexcReplicas; i++) {
      replicaTemperature[i] =
          params->repexcTemperatureLow +
          double(i) / double(params->repexcReplicas - 1.0) *
              (params->repexcTemperatureHigh - params->repexcTemperatureLow);
      replicaDynamics[i]->setTemperature(replicaTemperature[i]);
    }
  } else if (params->repexcTemperatureDistribution == "exponential") {
    double kTemp =
        log(params->repexcTemperatureHigh / params->repexcTemperatureLow) /
        (params->repexcReplicas - 1.0);
    // cout <<"kTemp: "<<kTemp<<endl;
    for (i = 0; i < params->repexcReplicas; i++) {
      replicaTemperature[i] = params->repexcTemperatureLow * exp(kTemp * i);
      replicaDynamics[i]->setTemperature(replicaTemperature[i]);
      log("replica: %ld temperature %.0f \n", i + 1, replicaTemperature[i]);
    }
  }

  log("\nReplica Exchange sampling for %.0f fs; %ld steps; %ld replicas.\n",
      params->repexcSamplingTime * 10.18, samplingSteps,
      params->repexcReplicas);

  for (step = 1; step <= samplingSteps; step++) {
    for (i = 0; i < params->repexcReplicas; i++) {
      replicaDynamics[i]->oneStep();
      cout << "step: " << step << " i " << i
           << " energy: " << replica[i]->getPotentialEnergy() << endl;
    }
    if ((step % exchangePeriodSteps) == 0) {
      for (long trial = 0; trial < params->repexcExchangeTrials; trial++) {
        i = helper_functions::randomInt(0, params->repexcReplicas - 2);
        energyLow = replica[i]->getPotentialEnergy();
        energyHigh = replica[i + 1]->getPotentialEnergy();
        kbTLow = kB * replicaTemperature[i];
        kbTHigh = kB * replicaTemperature[i + 1];
        pAcc = min(1.0, exp((energyHigh - energyLow) *
                            (1.0 / kbTHigh - 1.0 / kbTLow)));
        double tmp = helper_functions::randomDouble();
        cout << "step: " << step << " trial swap, i " << i
             << " elow: " << energyLow << " ehigh: " << energyHigh
             << " pAcc: " << pAcc << " rand: " << tmp;
        // if(helper_functions::randomDouble()<pAcc)
        if (tmp < pAcc) {
          cout << " swap\n";
          // swap configurations
          tmpMatter = replica[i];
          replica[i] = replica[i + 1];
          replica[i + 1] = tmpMatter;
          // reset velocities
          replicaDynamics[i]->setThermalVelocity();
          replicaDynamics[i + 1]->setThermalVelocity();
        } else {
          cout << " no_swap\n";
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

  fprintf(fileResults, "%ld random_seed\n", params->randomSeed);
  fprintf(fileResults, "%s potential_type\n",
          helper_functions::getPotentialName(params->potential).c_str());
  // fprintf(fileResults, "%ld force_calls_sampling\n", forceCalls);
  fclose(fileResults);

  std::string posFilename("pos_out.con");
  returnFiles.push_back(posFilename);
  filePos = fopen(posFilename.c_str(), "wb");
  fclose(filePos);

  return;
}
