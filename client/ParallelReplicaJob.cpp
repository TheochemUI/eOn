#include <algorithm>
#include <cstdlib>

#include "BondBoost.h"
#include "Dynamics.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "ParallelReplicaJob.h"

std::vector<std::string> ParallelReplicaJob::run(void) {
  // load pos.con
  reactant = new Matter(pot, params);
  reactant->con2matter(helper_functions::getRelevantFile(params->conFilename));

  // minimize the initial reactant
  SPDLOG_LOGGER_DEBUG(log, "[ParallelReplica] Minimizing initial position");
  reactant->relax();
  reactant->matter2con("reactant.con");

  // trajectory is the matter object for the current MD configuration
  Matter *trajectory = new Matter(pot, params);
  *trajectory = *reactant;
  Dynamics dynamics(trajectory, params.get());
  BondBoost bondBoost(trajectory, params.get());

  if (params->biasPotential == Hyperdynamics::BOND_BOOST) {
    bondBoost.initialize();
    trajectory->setBiasPotential(&bondBoost);
  }

  dephase(trajectory);

  // convert from simulation times to number of steps
  int stateCheckInterval =
      int(floor(params->parrepStateCheckInterval / params->mdTimeStep + 0.5));
  int recordInterval =
      int(floor(params->parrepRecordInterval / params->mdTimeStep + 0.5));

  std::vector<Matter *> MDSnapshots;
  std::vector<double> MDTimes;
  double transitionTime = 0;
  Matter transitionStructure(pot, params);
  int refineForceCalls = 0;

  // Main MD loop
  double simulationTime = 0.0;
  if (params->biasPotential == Hyperdynamics::NONE) {
    SPDLOG_LOGGER_DEBUG(
        log, "[ParallelReplica] {:>8} {:>12} {:>10} {:>12} {:>12} {:>10}",
        "Step", "Time (s)", "KE", "PE", "TE", "KinT");
  } else {
    SPDLOG_LOGGER_DEBUG(
        log,
        "[ParallelReplica] {:>8} {:>12} {:>10} {:>10} {:>12} {:>12} {:>10}",
        "Step", "Time (s)", "Boost", "KE", "PE", "TE", "KinT");
  }
  for (int step = 1; step <= params->mdSteps; step++) {
    dynamics.oneStep();
    double boost = 1.0;
    if (params->biasPotential == Hyperdynamics::BOND_BOOST) {
      double boostPotential = bondBoost.boost();
      double kB = params->kB;
      boost = exp(boostPotential / kB / params->temperature);

      simulationTime += params->mdTimeStep * boost;
    } else {
      simulationTime += params->mdTimeStep;
    }

    double kinE = trajectory->getKineticEnergy();
    double potE = trajectory->getPotentialEnergy();
    double kinT =
        (2.0 * kinE / (trajectory->numberOfFreeAtoms() * 3) / params->kB);

    if (step % params->writeMoviesInterval == 0) {
      if (params->biasPotential == Hyperdynamics::NONE) {
        SPDLOG_LOGGER_DEBUG(log,
                            "[ParallelReplica] {:>8} {:>12.4e} {:>10.4f} "
                            "{:>12.4f} {:>12.4f} {:>10.2f}",
                            step, simulationTime * params->timeUnit * 1e-15,
                            kinE, potE, kinE + potE, kinT);
      } else {
        double boostPotential = bondBoost.boost();
        SPDLOG_LOGGER_DEBUG(log,
                            "[ParallelReplica] {:>8} {:>12.4e} {:>10.3e} "
                            "{:>10.4f} {:>12.4f} {:>12.4f} {:>10.2f}",
                            step, simulationTime * params->timeUnit * 1e-15,
                            boost, kinE, potE + boostPotential,
                            kinE + potE + boostPotential, kinT);
      }
    }

    // Snapshots of the trajectory used for the refinement
    if (step % recordInterval == 0 && params->parrepRefineTransition) {
      Matter *tmp = new Matter(pot, params);
      *tmp = *trajectory;
      MDSnapshots.push_back(tmp);
      MDTimes.push_back(simulationTime);
    }

    // check for a transition every stateCheckInterval or at the end of the
    // simulation
    if (step % stateCheckInterval == 0 || step == params->mdSteps) {
      SPDLOG_LOGGER_DEBUG(log, "[ParallelReplica] Checking for transition");

      Matter min(pot, params);
      min = *trajectory;
      min.relax();

      // only check for a transition if one has yet to occur
      if (!min.compare(*reactant) && transitionTime == 0) {
        SPDLOG_LOGGER_DEBUG(log, "[ParallelReplica] Transition occurred");

        // perform the binary search for the transition structure
        if (params->parrepRefineTransition) {
          SPDLOG_LOGGER_DEBUG(log,
                              "[ParallelReplica] Refining transition time");
          // int tmpFcalls = Potential::fcalls;
          int snapshotIndex = refineTransition(MDSnapshots);

          // refineForceCalls += Potential::fcalls - tmpFcalls;

          transitionTime = MDTimes[snapshotIndex];
          transitionStructure = *MDSnapshots[snapshotIndex];

          // if not using refinement use the current configuration as the
          // transition structure
        } else {
          transitionStructure = *trajectory;
          transitionTime = simulationTime;
        }
        SPDLOG_LOGGER_DEBUG(log, "[ParallelReplica] Transition time: {:.3e} s",
                            transitionTime * params->timeUnit * 1e-15);

        // at the end of the simulation perform the refinement if it hasn't
        // happened yet this ensures that if a transition isn't seen that the
        // same number of force calls will be performed on average
      } else if (step + 1 == params->mdSteps && transitionTime == 0) {

        // fake refinement
        if (params->parrepRefineTransition) {
          SPDLOG_LOGGER_DEBUG(
              log,
              "[ParallelReplica] Simulation ended without seeing a transition");
          SPDLOG_LOGGER_DEBUG(
              log, "[ParallelReplica] Refining anyways to prevent bias...");
          int tmpFcalls = Potential::fcalls;
          refineTransition(MDSnapshots, true);
          // refineForceCalls += Potential::fcalls - tmpFcalls;
        }
        transitionStructure = *trajectory;
      }

      for (unsigned int i = 0; i < MDSnapshots.size(); i++) {
        delete MDSnapshots[i];
      }
      MDSnapshots.clear();
      MDTimes.clear();
    }
  }

  // start the decorrelation dynamics from the transition structure
  int decorrelationSteps =
      int(floor(params->parrepCorrTime / params->mdTimeStep + 0.5));
  SPDLOG_LOGGER_DEBUG(log, "[ParallelReplica] Decorrelating: {} steps",
                      decorrelationSteps);
  for (int step = 1; step <= decorrelationSteps; step++) {
    dynamics.oneStep(step);
  }
  SPDLOG_LOGGER_DEBUG(log, "[ParallelReplica] Decorrelation complete");

  // minimize the final structure
  Matter product(pot, params);
  product = *trajectory;
  product.relax();
  product.matter2con("product.con");

  // report the results
  FILE *fileResults;
  std::string resultsFilename("results.dat");
  returnFiles.push_back(resultsFilename);
  fileResults = fopen(resultsFilename.c_str(), "wb");
  fprintf(
      fileResults, "%s potential_type\n",
      std::string{magic_enum::enum_name<PotType>(params->potential)}.c_str());
  fprintf(fileResults, "%ld random_seed\n", params->randomSeed);
  fprintf(fileResults, "%f potential_energy_reactant\n",
          reactant->getPotentialEnergy());
  // fprintf(fileResults, "%i force_calls_refine\n", refineForceCalls);
  // fprintf(fileResults, "%d total_force_calls\n", Potential::fcalls);

  if (transitionTime == 0) {
    fprintf(fileResults, "0 transition_found\n");
    fprintf(fileResults, "%e simulation_time_s\n",
            simulationTime * params->timeUnit * 1.0e-15);
  } else {
    fprintf(fileResults, "1 transition_found\n");
    fprintf(fileResults, "%e transition_time_s\n",
            transitionTime * params->timeUnit * 1.0e-15);
    fprintf(fileResults, "%e correlation_time_s\n",
            params->parrepCorrTime * params->timeUnit * 1.0e-15);
    fprintf(fileResults, "%lf potential_energy_product\n",
            product.getPotentialEnergy());
  }
  fprintf(fileResults, "%lf speedup\n",
          simulationTime / (params->mdSteps * params->mdTimeStep));

  fclose(fileResults);

  MDSnapshots.clear();
  MDTimes.clear();
  delete trajectory;
  delete reactant;

  return returnFiles;
}

void ParallelReplicaJob::dephase(Matter *trajectory) {
  Dynamics dynamics(trajectory, params.get());

  int dephaseSteps =
      int(floor(params->parrepDephaseTime / params->mdTimeStep + 0.5));
  SPDLOG_LOGGER_DEBUG(log, "[ParallelReplica] Dephasing: {} steps",
                      dephaseSteps);

  Matter initial(pot, params);
  initial = *trajectory;

  while (true) {
    // always start from the initial configuration
    *trajectory = initial;
    dynamics.setThermalVelocity();

    // Dephase MD trajectory
    for (int step = 1; step <= dephaseSteps; step++) {
      dynamics.oneStep(step);
    }

    // Check to see if a transition occured
    Matter min(pot, params);
    min = *trajectory;
    min.relax();

    if (min.compare(*reactant)) {
      SPDLOG_LOGGER_DEBUG(log, "[ParallelReplica] Dephasing successful");
      break;
    } else {
      SPDLOG_LOGGER_DEBUG(
          log,
          "[ParallelReplica] Transition occured during dephasing; Restarting");
    }
  }
}

int ParallelReplicaJob::refineTransition(std::vector<Matter *> MDSnapshots,
                                         bool fake) {
  int min, max, mid;
  bool midTest;
  min = 0;
  max = MDSnapshots.size() - 1;

  while ((max - min) > 1) {
    mid = min + (max - min) / 2;
    Matter *snapshot = MDSnapshots[mid];
    snapshot->relax(true);

    if (fake == false) {
      midTest = snapshot->compare(*reactant);
    } else {
      // if we are faking the refinement just generate a random answer
      // for the comparison test
      midTest = bool(helper_functions::randomInt(0, 1));
    }

    if (midTest) {
      min = mid;
    } else {
      max = mid;
    }
  }

  return (min + max) / 2 + 1;
}
